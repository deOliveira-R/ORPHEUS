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

History — the retired ``sphere-4eg`` exclusion
----------------------------------------------
From 2026-05-19 (R-1 Step D) until 2026-08-10,
``test_krylov_inner_matches_tight_si`` imperatively xfailed the
``sphere × 4eg`` row on the grounds that unpreconditioned GMRES "exceeds
the ``max_inner=100`` budget without converging", pending issue #200 (the
block-inverse face preconditioner).

Retired as healed 2026-08-10, and #200 was not what healed it — #200 is
still open and GMRES here still runs with an explicit identity
preconditioner.  The cure lives in the GMRES ``restart``-sizing lineage:
**ERR-053** (2026-05-28) removed the ``restart=min(50, full_size)`` clamp,
and **#282 / #280 route (a)** (2026-07-04) sized ``restart`` from the full
augmented ravel instead of the bulk alone.  See the "History" section of
``test_kinf_homogeneous.py`` for the full account.  Measured with the
imperative xfail neutralised, at
this test's OWN budget (``max_inner=100``, ``inner_tol=1e-9`` — the exact
configuration the retired reason string named): the row passes with
``rel = 4.462e-13`` against the closed-form reference, four orders inside
its own ``1e-8`` bound, and with no
:class:`~orpheus.numerics.convergence.ConvergenceWarning`.

The exclusion outlived its cause by eleven weeks because
``pytest.xfail()`` called imperatively reports ``xfail`` unconditionally and
can never report ``XPASS``.  Any future exclusion here uses
``@pytest.mark.xfail(strict=True)``, which retires itself (issue #340, R5).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.analytical.homogeneous import (
    derive_2g_continuous,
    derive_4g_continuous,
)
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.numerics.quadrature import Quadrature
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
        return Quadrature.folded_product(n_mu=2, n_phi=4)
    return Quadrature.gauss_legendre(n_ordinates=8)


@pytest.mark.filterwarnings(
    "ignore::orpheus.numerics.convergence.ConvergenceWarning"
)
@pytest.mark.verifies("matrix-eigenvalue", "fission-matrix", "removal-matrix")
@pytest.mark.parametrize("coord", ["sphere", "cylinder"])
@pytest.mark.parametrize("ng_key", ["2eg", "4eg"])
def test_inner_tol_bias_collapses_at_1e_12(coord, ng_key):
    """k_inf drift drops below 1e-11 when inner_tol is tightened to 1e-12.

    Pins the inner-SI under-convergence model: at ``inner_tol=1e-9`` the
    drift is `[M]` 3.9e-11 … 7.5e-11 across the four parametrisations;
    tightening to 1e-12 collapses it to 6.5e-14 … 8.4e-14. If a future
    refactor breaks this contract (e.g., a sweep change that makes the
    inner converge by a different metric), this test surfaces it.

    ⛔ Until 2026-08-10 the paragraph above said the loose drift is
    ``~1e-7`` — **four orders too large**, and never measured.  The
    contrast the test rests on is real but is ~3 orders, not the ~4+ the
    old prose implied.  A stale magnitude in a bias-model docstring is
    exactly the number a future session would design a tolerance against.

    ⭐ **The LOOSE leg's truncation is the FIXTURE, and it is declared**
    (#340 R2).  ``max_inner=100`` at ``inner_tol=1e-9`` is deliberately
    starved: the under-converged inner IS the bias this test exhibits, so
    the loose solve exits truncated and says so.  `[M]` 2026-08-10, 3–4 of
    its 6–8 inners hit the cap at ρ ≈ 0.89–0.93; converging it would need
    ``max_inner ≈ 190–320`` and would collapse ``drift_loose`` onto
    ``drift_tight``, deleting the contrast the test is built on.  The
    suppression above therefore names the ONE category this row expects
    rather than silencing the module.

    ⚠ **The suppression is per-TEST and this test runs TWO solves**, so it
    covers the TIGHT leg as well — which is why the tight leg is budgeted
    to actually converge.  `[M]` until 2026-08-10 ``max_inner=300`` left
    the tight inner truncated in three of the four parametrisations
    (``2eg-sphere`` 2/6, ``2eg-cylinder`` 2/5, ``4eg-sphere`` 2/5; only
    ``4eg-cylinder`` was clean), so the *reference* half of a bias
    measurement was itself biased — and a bare marker would have silenced
    that too.  ``max_inner=420`` covers the measured worst case (409) with
    headroom, and `[M]` at 420 all four parametrisations run 0 truncated
    inners.  ``drift_tight`` is unmoved in magnitude — 6.3e-14 … 8.2e-14
    at 300 against 6.5e-14 … 8.4e-14 at 420, both three orders inside the
    1e-11 bar — and ``4eg-cylinder``, the one that never truncated, is
    bit-identical at 6.253e-14 either way.
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
        # 420, not 300: the tight leg is the REFERENCE this test measures
        # the loose leg against, so it must actually converge.  Measured
        # worst case across the four params is 409 (#340 R2).
        max_inner=420, inner_tol=1e-12,
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
    # ``sphere-4eg`` is gated here too — see the module docstring's
    # "History" section for the retired exclusion.
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
