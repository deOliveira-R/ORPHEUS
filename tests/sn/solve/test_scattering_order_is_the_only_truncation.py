r"""``scattering_order`` is the ONLY truncation — the data layer is lossless in ℓ (#426 step 1).

Two claims, both on the real 421-group library (`[M]` 2026-09-03):

* **The solver's order clamp reads the SCATTERING stack alone.** The (n,2n) stack of an
  isotope whose tape carries no MT=16 (H-1, B-10) is the zero P0 block — length 1 — and a
  clamp that read it would drop every mixture containing such an isotope to P0, deleting
  the elastic P1/P2 (`[M]` worth +5787 pcm-relative on the Be-reflected fixture, 14× the
  anisotropy #426 restores). Ruling O-1: a channel that stores fewer orders than requested
  is ZERO there and pads; it never lowers the solve's order. Witness: the same slab with
  B-10's (n,2n) stack as stored (length 1) and zero-padded to seven orders gives the SAME k
  to the bit — under a two-list clamp the first solve would be P0 and the second P2.
* **The library now serves the orders it stores.** Until step 1 the ingest cut every
  scattering channel to P2 and the clamp silently served P2 to a ``scattering_order = 3``
  request; the Be-reflected fast slab now reads a different k at P3 than at P2. The P2 value
  is a RECORD: at step 1 it was the pre-carve ``1.0953221881419453`` (`main` ``1e02f6b1``;
  lossless ingest moved no stored P0..P2 value), and since #426 step 2 (2026-09-04) the
  (n,2n) stack enters the solve at the scattering order, so the SAME fixture reads
  ``k(L=2) = 1.0911996566537725`` (`[M]` the flagship's ℓ ≤ 2 arm; the pre-step-2 value
  survives as that gate's P0-only control, where the (n,2n) ℓ ≥ 1 moments are zeroed).

Cost: four 421-group slab solves, `[M]` ≈ 3 + 3 + 7 + 8 s. Not ``slow`` (anti-pattern #36).
"""
from __future__ import annotations

import contextlib
import io
import warnings
from dataclasses import replace

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import compute_macro_xs
from orpheus.data.micro_xs import load_isotope
from orpheus.data.micro_xs.isotope import NG
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn

pytestmark = pytest.mark.l1

_TOL = dict(keff_tol=1e-9, flux_tol=1e-8, inner_tol=1e-10, max_outer=3000)
#: #426 step-2 record (2026-09-04), same fixture, same tolerances: the (n,2n) anisotropy
#: in the solve. Step 1's record was 1.0953221881419453 (now the flagship's P0-only control).
_K_P2_RECORD = 1.0911996566537725


def _mixture(name: str, density: float):
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # sigma0 = 0 clip on a pure isotope
        return compute_macro_xs([load_isotope(name, 294)], np.array([density]))


def _slab(refl_cm: float, core_cm: float, refl_cells: int, core_cells: int) -> Mesh1D:
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=1, outer_thickness_cm=refl_cm),
                 Region(mat_id=0, outer_thickness_cm=core_cm),
                 Region(mat_id=1, outer_thickness_cm=refl_cm)),
        bcs=(BC.vacuum, BC.vacuum),
    )
    return Mesh1D.from_geometry(geom, region_meshes=(
        RegionMesh(n_cells=refl_cells), RegionMesh(n_cells=core_cells), RegionMesh(n_cells=refl_cells)))


def _solve(materials, mesh, L: int, **tol):
    with contextlib.redirect_stdout(io.StringIO()):
        sol = solve_sn(materials, mesh, Quadrature.gauss_legendre(n_ordinates=8), scattering_order=L, **tol)
    assert sol.history is not None and sol.history.fully_converged
    assert sol.keff is not None
    return float(sol.keff)


class TestTheClampReadsTheScatteringStackAlone:
    def test_a_short_n2n_stack_does_not_lower_the_order(self):
        fuel = _mixture("U_235", 0.04894)
        b10 = _mixture("B_010", 0.05)
        assert len(b10.Sig2) == 1 and b10.Sig2[0].nnz == 0, "B-10 must be the (n,2n)-free witness"
        assert len(b10.SigS) == 7 and len(fuel.SigS) == 7
        padded = replace(b10, Sig2=[b10.Sig2[0]] + [csr_matrix((NG, NG)) for _ in range(6)])
        mesh = _slab(2.0, 4.0, 8, 16)
        tol = dict(keff_tol=1e-7, flux_tol=1e-6, inner_tol=1e-8, max_outer=3000)
        k_stored = _solve({0: fuel, 1: b10}, mesh, 2, **tol)
        k_padded = _solve({0: fuel, 1: padded}, mesh, 2, **tol)
        assert k_stored == k_padded, (
            f"k differs between a length-1 and a zero-padded (n,2n) stack "
            f"({k_stored!r} vs {k_padded!r}): the order clamp is reading Sig2"
        )


class TestTheLibraryServesTheOrdersItStores:
    @pytest.fixture(scope="class")
    def fixture(self):
        return {0: _mixture("U_235", 0.04894), 1: _mixture("BE009", 0.1236)}, _slab(3.0, 4.0, 12, 16)

    def test_p2_reproduces_the_record(self, fixture):
        materials, mesh = fixture
        assert _solve(materials, mesh, 2, **_TOL) == _K_P2_RECORD

    def test_p3_is_a_different_solve_than_p2(self, fixture):
        """Before step 1 a ``scattering_order = 3`` request was silently served P2."""
        materials, mesh = fixture
        assert all(len(m.SigS) == 7 for m in materials.values())
        k3 = _solve(materials, mesh, 3, **_TOL)
        assert abs(k3 - _K_P2_RECORD) > 100 * _TOL["keff_tol"], (
            f"P3 and P2 agree to {abs(k3 - _K_P2_RECORD):.2e}: the third elastic order is not reaching the solve"
        )
