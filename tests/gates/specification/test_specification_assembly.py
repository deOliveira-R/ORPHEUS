r"""The assembly rung of the specification (#405 P1 step 8, S8.4).

Specified by the test-architect (2026-10-02, ``.claude/plans/reference_p1_spec.md``
§1.8). An integration rung: green on arrival (no first red of its own), it
shows that a specification carries what a method needs, through the test's own
lift, since posing a method from a specification is P4's (no solver change in
P1). Its teeth are the reference side's activation legs and the battery arms
named in the spec.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.data.materials import Materials
from orpheus.geometry import BC, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.question import Eigen
from orpheus.specification import GeometrySpecification, InfiniteMediumSpecification
from tests.gates._content_identity_helpers import require
from tests.gates.specification._fixtures import fuel, moderator, slab2

pytestmark = pytest.mark.l1

F = Channel.FISSION_EMISSION


def _k() -> CellCoefficient:
    return CellCoefficient.every(F)


def _two_materials() -> Materials:
    return Materials({0: fuel(), 1: moderator()})


# ═════════════════════════════════════════════════════════════════════════════
# S8.4 — the assembly rung
# ═════════════════════════════════════════════════════════════════════════════


def _dense_k_inf(mixture) -> float:
    """k-infinity from the raw ``Mixture`` arrays, independent of every operator:
    A = diag(Sigma_t) - Sigma_s0^T - 2 Sigma_2^T, F = chi (x) nu Sigma_f, k = rho(A^-1 F)."""
    a = np.diag(mixture.SigT) - mixture.SigS[0].toarray().T - 2.0 * mixture.Sig2[0].toarray().T
    f = np.outer(mixture.chi, mixture.SigP)
    return float(max(abs(np.linalg.eigvals(np.linalg.solve(a, f)))))


def test_s8_4_a_two_interval_specification_assembles_a_material_mesh() -> None:
    from orpheus.transport.mesh import MaterialMesh

    spec = GeometrySpecification(Materials({**_two_materials(), 9: moderator()}), slab2(), Eigen(_k()))
    require(sorted(spec.materials.ids) == [0, 1], f"the spectator 9 was kept: {sorted(spec.materials.ids)}")
    mesh = Mesher(spec.geometry).partition(CellsByCount.uniform_width(3)).mesh
    mm = MaterialMesh(mesh, spec.materials)
    require(mm.ng == spec.materials.uniform_group_count() == 2, f"ng {mm.ng}")
    require(mm.materials == spec.materials, "the mesh holds another declaration")
    occupied = sorted(m for m, idx in mm.cells_by_material.items() if idx[0].size)
    require(occupied == [0, 1], f"occupied materials {occupied}")


@pytest.mark.rests_on(
    "tests/gates/sn/solve/test_d3_admission.py::test_kinf_3d_equals_2d_equals_1d_homogeneous_reflective",
    "tests/gates/derivations/test_adjoint_spectrum_reference.py::TestAdjointSpectrumReference",
)
@pytest.mark.parametrize("lift", ["homogeneous", "sn-reflective-slab"])
def test_s8_4_the_infinite_medium_k_is_the_dense_pencil_k(lift: str) -> None:
    """``InfiniteMediumSpecification(0, fuel, Eigen(every(F)))`` asks k-infinity of its one material.

    The fixture activates what a 2-group infinite medium can: upscatter and a
    non-zero (n,2n) emission (asserted, with the reference's sensitivity to the
    multiplicity 2). Tolerances: the 0-D solve is direct (``[M]`` 1 ULP from the
    dense value; 50 ULP admitted, never a LAPACK bit pin, vv anti-#38); the SN
    solve is iterative, 10 x its ``keff_tol`` (``[M]`` 1.5e-13 at 1e-10)."""
    from orpheus.homogeneous.solver import solve_homogeneous_infinite
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.solver import solve_sn

    spec = InfiniteMediumSpecification(0, fuel(), Eigen(_k()))
    require(isinstance(spec.question, Eigen) and spec.question.parameter == CellCoefficient({(0, F)}), "the rung does not ask k")
    mixture = spec.mixture
    require(mixture.Sig2[0].count_nonzero() > 0 and mixture.SigS[0].toarray()[1, 0] > 0, "activation: (n,2n), upscatter")
    k_ref = _dense_k_inf(mixture)
    a1 = np.diag(mixture.SigT) - mixture.SigS[0].toarray().T - mixture.Sig2[0].toarray().T
    k_mult1 = float(max(abs(np.linalg.eigvals(np.linalg.solve(a1, np.outer(mixture.chi, mixture.SigP))))))
    require(abs(k_ref - k_mult1) > 1e-2, "activation: the reference does not read the (n,2n) multiplicity")
    if lift == "homogeneous":
        k = solve_homogeneous_infinite(mixture).k_inf
        tol = 50 * np.spacing(k_ref)
    else:
        keff_tol = 1e-10
        mesh = Mesher(StructuredGeometry.from_homogeneous(0.5, BC.reflective)).partition(CellsByCount.uniform_width(2)).mesh
        k = solve_sn(dict(spec.materials.items()), mesh, Quadrature.gauss_legendre(n_ordinates=8), keff_tol=keff_tol, inner_tol=1e-11).outcome.keff
        tol = 10 * keff_tol
    require(abs(k - k_ref) <= tol, f"{lift}: k {k!r} vs the dense pencil {k_ref!r} (|diff| {abs(k - k_ref):.3e} > {tol:.1e})")
