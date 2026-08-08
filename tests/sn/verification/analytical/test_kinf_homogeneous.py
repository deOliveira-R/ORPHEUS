"""L1 analytical-reference test: SN reproduces k_inf and flux spectrum.

Verifies the SN solver against the closed-form continuous reference at
``orpheus/derivations/continuous/analytical/homogeneous.py``. The
reference is **structurally independent** of the SN sweep — it is the
dominant eigenvector + eigenvalue of :math:`A^{-1}F` solved by pure
linear algebra on the cross-section matrices, with no transport
discretisation involved. Cross-checks an SN solver run against this
reference exercise the discretisation chain end-to-end against a
mathematically-orthogonal ground.

Quadrature choices (claim-driven — each choice is exact for the claim
it verifies):

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
flux-shape independent per the 1-group-degeneracy rule).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.analytical.homogeneous import (
    derive_1g_continuous,
    derive_2g_continuous,
    derive_4g_continuous,
)
from orpheus.geometry import (
    BC,
    Mesh1D,
    Region,
    RegionMesh,
    StructuredGeometry,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn


pytestmark = pytest.mark.l1


# ─── helpers ─────────────────────────────────────────────────────────


_GEOMETRY_TAG = {"slab": "SLB", "sphere": "SPH", "cylinder": "CYL"}


_CASE_BUILDERS = {
    "1eg": derive_1g_continuous,
    "2eg": derive_2g_continuous,
    "4eg": derive_4g_continuous,
}


def _get_continuous_case(ng_key: str):
    """Bypass the ``continuous_get`` registry to skip its slow auto-discovery
    walk over every ``orpheus.derivations.*`` module. We point at the
    specific analytical homogeneous derivation module the test consumes.

    This is intentional V&V hygiene: name the reference path explicitly.
    """
    return _CASE_BUILDERS[ng_key]()


def _homogeneous_mesh(coord: str, n_cells: int, length: float, mat_id: int) -> Mesh1D:
    """Single-region homogeneous mesh in the requested 1D coordinate.

    Slab takes (left, right) BCs; sphere/cylinder take only the outer
    BC (the inner r=0 symmetry is implicit). Reflective everywhere
    guarantees no neutron leakage so the discrete eigenvalue converges
    to k_inf.
    """
    bcs = (BC.reflective, BC.reflective) if coord == "slab" else (BC.reflective,)
    geom = StructuredGeometry(
        geometry=_GEOMETRY_TAG[coord],
        regions=(Region(mat_id=mat_id, outer_thickness_cm=length),),
        bcs=bcs,
    )
    return Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=n_cells),))


def _quadrature_for(coord: str):
    if coord == "cylinder":
        return Quadrature.folded_product(n_mu=2, n_phi=4)  # 4 folded ordinates, level_indices populated
    return Quadrature.gauss_legendre(n_ordinates=8)


# ─── eigenvalue + spectrum tests ─────────────────────────────────────


# Tight tolerances so the analytical-reference comparison is solver-floor
# correct on EVERY coord × ng × inner_solver combination.  Per the
# 2026-05-19 numerics-investigator memo (curvilinear MG inner-tol
# amplification): default ``inner_tol=1e-9`` lets the within-group SI
# residual bound ``‖φ_n − φ_(n−1)‖ / ‖φ_n‖ < 1e-9`` propagate as a gap
# of ``ρ/(1−ρ) × 1e-9 ≈ 9 × 1e-9 ≈ 1e-7`` to ``k_inf`` for the c=0.9
# thermal-group case on curvilinear geometries (slab SI converges fast
# enough to mask this).  Tightening ``inner_tol`` to 1e-12 collapses the
# drift below 1e-11; running with ``inner_solver="krylov"`` removes the
# SI amplification entirely.  Both inner-solver paths are now gated by
# this L1 reference.
_TIGHT_KW = dict(
    max_outer=1000, keff_tol=1e-14, flux_tol=1e-12,
    # max_inner bumped from 300 to 1000 (2026-05-26) to accommodate the
    # ERR-052 production-rate-normalised iteration trajectory.  Empirically,
    # sphere-2eg-krylov needs ~600 GMRES iterations per outer step under
    # the new trajectory (the normalisation alters the initial guess each
    # outer iter, so GMRES re-builds its Krylov subspace from a non-warmed
    # state).  Below 1000, the inner solver hits the cap, returns an
    # under-converged result, and the outer iteration accumulates ~2.4e-7
    # drift in keff.  With max_inner=1000 every (coord, ng, inner_solver)
    # variant reaches FP-precision; the rtol=1e-10 keff and rtol=1e-9
    # spectrum gates hold for all 28 cases.  Issue #200 (block-inverse
    # preconditioner) will eventually let us reduce this back to 300.
    max_inner=1000, inner_tol=1e-12,
)


@pytest.mark.verifies("matrix-eigenvalue", "fission-matrix", "removal-matrix")
@pytest.mark.parametrize("inner_solver", ["source_iteration", "krylov"])
@pytest.mark.parametrize("ng_key", ["1eg", "2eg", "4eg"])
@pytest.mark.parametrize("coord", ["slab", "sphere", "cylinder"])
def test_kinf_homogeneous(ng_key: str, coord: str, inner_solver: str) -> None:
    """SN reproduces analytical k_inf on homogeneous medium, every coord × ng.

    Parametrised over BOTH ``inner_solver`` paths (Option B+D per the
    numerics-investigator memo): the source-iteration path is the
    default production path and MUST converge to the analytical
    reference at tight tolerance; the Krylov path is the
    symmetric-closure path and provides structural independence — the
    two paths must reach the SAME k_inf to better than ``rtol=1e-10``.
    """
    # R-1 Step D (2026-05-19) — the carved ``_solve_krylov`` runs GMRES
    # UNPRECONDITIONED on curvilinear meshes (sphere / cylinder) per
    # the user's R-1 direction "no preconditioner; consolidating the
    # foundational architecture first".  The default sweep-as-
    # preconditioner cannot be used because GMRES feeds residual
    # vectors that have no ``rhs(1)`` history, so the curvilinear
    # Carlson coupled-pole seed (Hébert §3.9.4) falls back to the
    # in-iteration-source default — not the algebraic inverse of
    # ``(L+C).apply``.  Unpreconditioned GMRES converges on
    # ``sphere-1g/-2g`` and ``cylinder-1g/-2g/-4g`` within the
    # ``max_inner=300`` budget, but ``sphere-4g`` is too ill-
    # conditioned (the 4-group cross-section coupling on the sphere
    # pole inflates the spectral radius beyond what 300 GMRES iters
    # can resolve).  Tracked by issue #200 — the block-inverse face
    # preconditioner will re-enable this case once it lands.
    if coord == "sphere" and ng_key == "4eg" and inner_solver == "krylov":
        pytest.xfail(
            "R-1 — unpreconditioned GMRES on sphere-4g exceeds the "
            "max_inner=300 budget without converging.  Issue #200 "
            "tracks the block-inverse face preconditioner that "
            "re-enables Krylov on curvilinear multi-group at high "
            "scattering ratio."
        )

    case = _get_continuous_case(ng_key)
    mat_id = next(iter(case.problem.materials.keys()))

    mesh = _homogeneous_mesh(coord=coord, n_cells=10, length=2.0, mat_id=mat_id)
    quadrature = _quadrature_for(coord)

    result = solve_sn(
        case.problem.materials, mesh, quadrature,
        inner_solver=inner_solver, **_TIGHT_KW,
    )

    np.testing.assert_allclose(
        result.keff, case.k_eff, rtol=1e-10,
        err_msg=(
            f"SN k_inf disagrees with analytical reference: "
            f"got {result.keff!r}, expected {case.k_eff!r} "
            f"(coord={coord}, ng={ng_key}, inner_solver={inner_solver})"
        ),
    )


# ── Sentinel (canary) — closed-form eigenvalue pillar ───────────────
# One sentinel for the ``verification_analytical`` capability node: the
# cheapest NON-DEGENERATE closed-form k_inf cross-check (2G slab krylov).
# This is the EIGENVALUE-claim pillar (closed-form analytical reference,
# NOT MMS — MMS does not prove eigenvalues per vv-principles). A flip
# here localizes to the eigenvalue + analytical-reference cluster.
# See .claude/plans/sn_sentinel_harness.md.


@pytest.mark.sentinel
@pytest.mark.verifies("matrix-eigenvalue", "fission-matrix", "removal-matrix")
def test_sentinel_kinf_slab_2g_krylov() -> None:
    """Sentinel: closed-form 2G slab k_inf via the krylov path."""
    test_kinf_homogeneous(ng_key="2eg", coord="slab", inner_solver="krylov")


@pytest.mark.verifies("matrix-eigenvalue", "fission-matrix", "removal-matrix")
@pytest.mark.parametrize("inner_solver", ["source_iteration", "krylov"])
@pytest.mark.parametrize("ng_key", ["2eg", "4eg"])
@pytest.mark.parametrize("coord", ["slab", "sphere", "cylinder"])
def test_kinf_homogeneous_spectrum(ng_key: str, coord: str, inner_solver: str) -> None:
    """SN reproduces the dominant ``A^{-1}F`` eigenvector (multi-group only).

    A 1G test cannot pin the spectrum (only one component), so this test
    is parametrised only over multi-group cases. This is the test that
    catches scattering-transpose bugs that preserve trace (and thus
    k_inf) but flip the eigenvector — failure mode #2 (variable swap)
    in `vv-principles`.

    Parametrised over both inner solvers per the 2026-05-19 memo.  The
    spectrum is solved at TIGHT tolerances so the SI inner-residual
    amplification does not bias the dominant eigenvector recovery.
    """
    # R-1 Step D — see ``test_kinf_homogeneous`` for the structural
    # rationale.  Sphere-4g unpreconditioned GMRES does not converge
    # within ``max_inner=300``.  Issue #200 tracks the fix.
    if coord == "sphere" and ng_key == "4eg" and inner_solver == "krylov":
        pytest.xfail(
            "R-1 — unpreconditioned GMRES on sphere-4g exceeds the "
            "max_inner=300 budget without converging.  Issue #200 "
            "tracks the block-inverse face preconditioner."
        )

    case = _get_continuous_case(ng_key)
    ng = case.problem.n_groups
    mat_id = next(iter(case.problem.materials.keys()))

    mesh = _homogeneous_mesh(coord=coord, n_cells=10, length=2.0, mat_id=mat_id)
    quadrature = _quadrature_for(coord)
    result = solve_sn(
        case.problem.materials, mesh, quadrature,
        inner_solver=inner_solver, **_TIGHT_KW,
    )

    # Spatial average per group → spectrum vector (homogeneous → spatially flat).
    # ``scalar_flux.values`` carries the canonical ``(ng, *spatial)`` layout
    # (Issue #197 PR-TYPED-5); averaging over every spatial axis (all axes
    # except the leading group axis 0) leaves the group axis ``ng``.
    flux_vals = result.scalar_flux.values
    spatial_axes = tuple(range(1, flux_vals.ndim))
    phi_per_group = flux_vals.mean(axis=spatial_axes)
    phi_solver = phi_per_group / np.linalg.norm(phi_per_group)

    phi_ref = np.array([case.phi(0.0, g) for g in range(ng)], dtype=float)
    phi_ref = phi_ref / np.linalg.norm(phi_ref)

    # eigenvector sign degree of freedom — pin first component positive
    if phi_solver[0] < 0:
        phi_solver = -phi_solver
    if phi_ref[0] < 0:
        phi_ref = -phi_ref

    np.testing.assert_allclose(
        phi_solver, phi_ref, rtol=1e-9,
        err_msg=(
            f"SN flux spectrum disagrees with analytical reference: "
            f"got {phi_solver!r}, expected {phi_ref!r} "
            f"(coord={coord}, ng={ng_key}, inner_solver={inner_solver})"
        ),
    )
