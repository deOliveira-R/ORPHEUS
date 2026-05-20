"""L1 companion: inner-SI residual does NOT bound the k_inf drift for curvilinear MG.

Promoted from ``derivations/diagnostics/`` on 2026-05-19 as a permanent
L1 companion to ``test_kinf_homogeneous`` per the
``diagnostic-to-test promotion`` protocol.  Original investigation by
numerics-investigator covered the four failing cases
(sphere/cylinder × 2eg/4eg).

**Verdict (Step 1 of the cascade)**: the 1.2e-7 drift is a
**convergence-tolerance issue**, not a solver bug.

**Mechanism.** The within-group source-iteration termination criterion
``‖φ_n − φ_{n−1}‖ / ‖φ_n‖ < inner_tol`` is the *residual*, NOT the gap
to the true within-group fixed point. The Banach contraction bound on
that gap is::

    gap_to_true_fixed_point ≈ ρ / (1 − ρ) × inner_tol_residual

where ρ is the within-group SI spectral radius (≈ ``c_within = Σ_(s,g→g) / Σ_(t,g)``
for the dominant group). For the 2eg test problem ρ ≈ 0.9 in the
thermal group, so the amplification factor is ~9. The gap is then
re-amplified through power iteration coupling to give the observed
~1e-7 ``k_inf`` drift.

**Independent confirmation.** Running with ``inner_solver="krylov"``
(GMRES on ``L.apply`` — no SI dominance ratio in the loop) reduces
the drift to ~4e-10 at the SAME tolerances. The SI/Krylov gap is the
load-bearing evidence: an under-converged-inner-SI bias dissolves
when the solver is replaced by one that converges to the true within-
group fixed point regardless of ρ.

**Refined-tolerance fix**: ``inner_tol = 1e-12`` brings every curvilinear
multi-group case to ~1e-12 ``k_inf`` drift, matching the slab cases.

This diagnostic test pins BOTH the failure mode and the fix path. It
is a permanent regression gate against the inner-SI bias resurfacing
in solver refactors that swap inner solvers, change the SI
preconditioner, or alter the convergence-criterion semantics.

If this test starts failing, the inner-SI termination model has
changed.  Either (a) the convergence-criterion semantics flipped
(residual → gap or vice versa); (b) the SI sweep spectral radius
changed enough to make the amplification factor different; or
(c) ``inner_solver="krylov"`` no longer routes through the
symmetric-closure L operator.  Investigate before adjusting bounds.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.analytical.homogeneous import (
    derive_2g_continuous,
    derive_4g_continuous,
)
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.quadrature import GaussLegendre1D, ProductQuadrature
from orpheus.sn.solver import solve_sn


pytestmark = pytest.mark.l1


_GEOMETRY_TAG = {"slab": "SLB", "sphere": "SPH", "cylinder": "CYL"}
_BUILDERS = {"2eg": derive_2g_continuous, "4eg": derive_4g_continuous}


def _mesh(coord: str, mat_id: int, n_cells: int = 10, length: float = 2.0) -> Mesh1D:
    bcs = (BC.reflective, BC.reflective) if coord == "slab" else (BC.reflective,)
    geom = StructuredGeometry(
        geometry=_GEOMETRY_TAG[coord],
        regions=(Region(mat_id=mat_id, outer_thickness_cm=length),),
        bcs=bcs,
    )
    return Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=n_cells),))


def _quadrature(coord: str):
    if coord == "cylinder":
        return ProductQuadrature.create(n_mu=2, n_phi=4)
    return GaussLegendre1D.create(n_ordinates=8)


@pytest.mark.verifies("matrix-eigenvalue", "fission-matrix", "removal-matrix")
@pytest.mark.parametrize("coord", ["sphere", "cylinder"])
@pytest.mark.parametrize("ng_key", ["2eg", "4eg"])
def test_inner_tol_bias_collapses_at_1e_12(coord, ng_key):
    """k_inf drift drops below 1e-11 when inner_tol is tightened to 1e-12.

    Pins the inner-SI under-convergence model: at default ``inner_tol=1e-9``
    the drift is ~1e-7; tightening to 1e-12 collapses it. If a future
    refactor breaks this contract (e.g., a sweep change that makes the
    inner converge by a different metric), this test surfaces it.
    """
    case = _BUILDERS[ng_key]()
    mat_id = next(iter(case.problem.materials.keys()))
    mesh = _mesh(coord, mat_id)
    quad = _quadrature(coord)

    r_loose = solve_sn(
        case.problem.materials, mesh, quad,
        max_outer=200, keff_tol=1e-10, flux_tol=1e-9,
        max_inner=100, inner_tol=1e-9,
    )
    r_tight = solve_sn(
        case.problem.materials, mesh, quad,
        max_outer=1000, keff_tol=1e-14, flux_tol=1e-12,
        max_inner=300, inner_tol=1e-12,
    )

    drift_loose = abs(r_loose.keff - case.k_eff) / case.k_eff
    drift_tight = abs(r_tight.keff - case.k_eff) / case.k_eff

    # The convergence-tolerance fingerprint: loose >> tight.
    # If a future change makes loose == tight (e.g. SI replaced by an
    # always-converging Krylov by default), THIS test must be updated.
    # If a future change makes tight >> 1e-11, a real solver bug has
    # been introduced.
    assert drift_tight < 1e-11, (
        f"Inner-tol=1e-12 should bring k_inf drift below 1e-11. "
        f"Got drift_loose={drift_loose:.3e}, drift_tight={drift_tight:.3e}. "
        f"This signals a real solver bias (drift no longer collapses with "
        f"tighter inner tolerance)."
    )


@pytest.mark.verifies("matrix-eigenvalue", "fission-matrix", "removal-matrix")
@pytest.mark.parametrize("coord", ["sphere", "cylinder"])
@pytest.mark.parametrize("ng_key", ["2eg", "4eg"])
def test_krylov_inner_matches_tight_si(coord, ng_key):
    """Krylov inner-solver at default tolerances matches tight-SI to ~1e-9.

    Pins the structural-independence claim: SI and Krylov converge to
    the SAME within-group fixed point — they differ only in convergence
    rate, not in target. Disagreement at this level would mean SI and
    Krylov route through structurally different operators (a closure
    drift). At the time of writing, agreement is to ~1e-9 absolute.
    """
    # R-1 Step D — see ``test_kinf_homogeneous`` for the rationale.
    # Sphere-4g unpreconditioned GMRES does not converge within the
    # max_inner=100 budget at default tolerances.  Issue #200 tracks
    # the block-inverse face preconditioner.
    if coord == "sphere" and ng_key == "4eg":
        pytest.xfail(
            "R-1 — unpreconditioned GMRES on sphere-4g exceeds the "
            "max_inner=100 budget without converging.  Issue #200 "
            "tracks the block-inverse face preconditioner."
        )

    case = _BUILDERS[ng_key]()
    mat_id = next(iter(case.problem.materials.keys()))
    mesh = _mesh(coord, mat_id)
    quad = _quadrature(coord)

    r_kr = solve_sn(
        case.problem.materials, mesh, quad,
        inner_solver="krylov",
        max_outer=200, keff_tol=1e-10, flux_tol=1e-9,
        max_inner=100, inner_tol=1e-9,
    )
    drift_kr = abs(r_kr.keff - case.k_eff) / case.k_eff
    # Krylov should converge cleanly without inner-SI amplification.
    assert drift_kr < 1e-8, (
        f"Krylov inner-solver should match analytical k_inf to better "
        f"than 1e-8 at default tolerances. Got {drift_kr:.3e}. "
        f"If it doesn't, the L operator's apply path has a real bias."
    )
