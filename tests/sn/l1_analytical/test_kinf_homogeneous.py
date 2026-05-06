"""L1 analytical-reference test: SN reproduces k_inf and flux spectrum.

Verifies the SN solver against the closed-form continuous reference at
``orpheus/derivations/continuous/analytical/homogeneous.py``. The
reference is **structurally independent** of the SN sweep — it is the
dominant eigenvector + eigenvalue of :math:`A^{-1}F` solved by pure
linear algebra on the cross-section matrices, with no transport
discretisation involved. Cross-checks an SN solver run against this
reference exercise the discretisation chain end-to-end against a
mathematically-orthogonal ground.

Quadrature choices (claim-driven per `docs/testing/sn_verification_matrix.rst`):

* **Slab / Sphere**: GL-1D N=8. Exact for polynomials of degree ≤ 15
  in :math:`\\mu`. Homogeneous flux is degree 0 in :math:`\\mu`, so
  N=8 is far above what the claim demands.
* **Cylinder**: ProductQuadrature(n_mu=2, n_phi=4) — 8 ordinates total.
  The cylindrical sweep requires per-level azimuthal structure, so a
  ``level_indices``-bearing quadrature is mandatory.

Each test is a single SN run on a 10-cell mesh — peak ψ-array footprint
is ``(N_ord × n_cells × 1 × n_groups × 8 bytes) ≤ 8 × 10 × 1 × 4 × 8 =
2.5 KB`` per array. With ~5 intermediate arrays held simultaneously the
peak memory is well under 100 MB including pytest/numpy overhead.

Catches: any sign-flip / variable-swap / convention-drift bug that
disturbs ``A^{-1}F`` assembly or its dominant eigenvector — including
the eigenvector-flip class invisible to a 1G k_inf test (which is
flux-shape independent per Cardinal Rule 6).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations import continuous_get
from orpheus.geometry import (
    BC,
    Mesh1D,
    Region,
    RegionMesh,
    StructuredGeometry,
)
from orpheus.sn.quadrature import GaussLegendre1D, ProductQuadrature
from orpheus.sn.solver import solve_sn


pytestmark = pytest.mark.l1


# ─── helpers ─────────────────────────────────────────────────────────


_GEOMETRY_TAG = {"slab": "SLB", "sphere": "SPH", "cylinder": "CYL"}


def _homogeneous_mesh(coord: str, n_cells: int, length: float, mat_id: int) -> Mesh1D:
    """Single-region homogeneous mesh in the requested 1D coordinate.

    Reflective BCs on both ends (slab) or symmetry+reflective (sphere /
    cylinder) — guarantees no neutron leakage so the discrete
    eigenvalue converges to k_inf.
    """
    geom = StructuredGeometry(
        geometry=_GEOMETRY_TAG[coord],
        regions=(Region(mat_id=mat_id, outer_thickness_cm=length),),
        bcs=(BC.reflective, BC.reflective),
    )
    return Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=n_cells),))


def _quadrature_for(coord: str):
    if coord == "cylinder":
        return ProductQuadrature.create(n_mu=2, n_phi=4)  # 8 ordinates, level_indices populated
    return GaussLegendre1D.create(n_ordinates=8)


# ─── eigenvalue + spectrum tests ─────────────────────────────────────


@pytest.mark.verifies("matrix-eigenvalue", "fission-matrix", "removal-matrix")
@pytest.mark.parametrize("ng_key", ["1eg", "2eg", "4eg"])
@pytest.mark.parametrize("coord", ["slab", "sphere", "cylinder"])
def test_kinf_homogeneous(ng_key: str, coord: str) -> None:
    """SN reproduces analytical k_inf on homogeneous medium, every coord × ng."""
    case = continuous_get(f"homo_{ng_key}")
    mat_id = next(iter(case.problem.materials.keys()))

    mesh = _homogeneous_mesh(coord=coord, n_cells=10, length=2.0, mat_id=mat_id)
    quadrature = _quadrature_for(coord)

    result = solve_sn(
        case.problem.materials, mesh, quadrature,
        max_outer=300, keff_tol=1e-12, flux_tol=1e-10,
        max_inner=200, inner_tol=1e-10,
    )

    np.testing.assert_allclose(
        result.keff, case.k_eff, rtol=1e-9,
        err_msg=(
            f"SN k_inf disagrees with analytical reference: "
            f"got {result.keff!r}, expected {case.k_eff!r} "
            f"(coord={coord}, ng={ng_key})"
        ),
    )


@pytest.mark.verifies("matrix-eigenvalue", "fission-matrix", "removal-matrix")
@pytest.mark.parametrize("ng_key", ["2eg", "4eg"])
@pytest.mark.parametrize("coord", ["slab", "sphere", "cylinder"])
def test_kinf_homogeneous_spectrum(ng_key: str, coord: str) -> None:
    """SN reproduces the dominant ``A^{-1}F`` eigenvector (multi-group only).

    A 1G test cannot pin the spectrum (only one component), so this test
    is parametrised only over multi-group cases. This is the test that
    catches scattering-transpose bugs that preserve trace (and thus
    k_inf) but flip the eigenvector — failure mode #2 (variable swap)
    in `vv-principles`.
    """
    case = continuous_get(f"homo_{ng_key}")
    ng = case.problem.n_groups
    mat_id = next(iter(case.problem.materials.keys()))

    mesh = _homogeneous_mesh(coord=coord, n_cells=10, length=2.0, mat_id=mat_id)
    quadrature = _quadrature_for(coord)
    result = solve_sn(
        case.problem.materials, mesh, quadrature,
        max_outer=300, keff_tol=1e-12, flux_tol=1e-10,
        max_inner=200, inner_tol=1e-10,
    )

    # spatial average per group → spectrum vector (homogeneous → spatially flat)
    phi_per_group = result.scalar_flux.mean(axis=(0, 1))
    phi_solver = phi_per_group / np.linalg.norm(phi_per_group)

    phi_ref = np.array([case.phi(0.0, g) for g in range(ng)], dtype=float)
    phi_ref = phi_ref / np.linalg.norm(phi_ref)

    # eigenvector sign degree of freedom — pin first component positive
    if phi_solver[0] < 0:
        phi_solver = -phi_solver
    if phi_ref[0] < 0:
        phi_ref = -phi_ref

    np.testing.assert_allclose(
        phi_solver, phi_ref, rtol=1e-8,
        err_msg=(
            f"SN flux spectrum disagrees with analytical reference: "
            f"got {phi_solver!r}, expected {phi_ref!r} "
            f"(coord={coord}, ng={ng_key})"
        ),
    )
