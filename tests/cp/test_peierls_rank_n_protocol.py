"""L1 protocol test — mandatory >=2-quadrature signed-error stability
for any rank-N closure claim against F.4.

GitHub Issue #123 (Direction N). See
``.claude/plans/rank-n-closure-research-log.md`` Session 2026-04-22 (late),
lessons L17-L19.

Background
----------
L16 (Issue #120 retraction, 2026-04-22): compare rank-N closures at matched
quadrature.

L17-L19 (Session 2026-04-22 late, Issues #121 falsification and #123 codification):
*even F.4's own RICH reference is quadrature-limited at sigma_t * R >= 10*.

- F.4 anchor (sigma_t*R=10, rho=0.3): RICH err 0.0033 %, ULTRA err 0.0010 %
  -- ~3x lower.
- F.4 (sigma_t*R=20, rho=0.5): RICH err +0.0054 %, RICH+panels err -0.0039 %
  -- SIGN FLIP.

Consequence: a closure C evaluated at quadrature Q whose residual is not
separated from |err(F.4, Q)| can show an apparent win that is just the
cancellation ``|eps_struct - eps_quad|`` evaluated in the lucky direction.
To rule this out the claim must be verified at >= 2 quadratures.

This module ships:

1. ``assert_rank_n_structural_win(closure_fn, f4_fn, point, quads)`` --
   the protocol helper. Any future rank-N closure claim imports this.
2. Unit tests of the helper: trivial-pass (F.4 vs itself) and
   falsification (synthetic signed-flip case).
3. A parametrized baseline test that pins L17-L19's F.4 sign-stability
   claim at a documented per-point reference quadrature on the
   canonical 6-point grid.

Verifies Sphinx label ``peierls-rank-n-stability`` (to be added to
``docs/theory/peierls_unified.rst`` when Issue #123 is closed; see
file-level TODO below).

Runtime
-------
The baseline test runs F.4 at two quadratures per point in child
processes, capped at 120 s/run. Whole file is ``@pytest.mark.slow``:
total wall ~6-15 min depending on machine. The helper unit tests
themselves are fast (no F.4 evaluation).
"""

# TODO(Issue #123 close-out): add the ``:label: peierls-rank-n-stability``
# block to ``docs/theory/peierls_unified.rst`` next to the existing
# ``peierls-rank-n-bc-closure`` label. The policy text should cite L17-L19.
from __future__ import annotations

import json
import subprocess
import sys
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Sequence

import pytest

K_INF = 1.5  # nu_sig_f / (sig_t - sig_s) = 1.0 / (1.0 - 1/3) = 1.5


# =========================================================================
# Protocol helper -- the contract any rank-N claim must pass.
# =========================================================================


@dataclass(frozen=True)
class StabilityReport:
    """Outcome of a 2+ quadrature signed-error stability check."""

    closure_signed: tuple[float, ...]       # per-quad signed residuals
    f4_signed: tuple[float, ...]            # per-quad signed residuals
    quads: tuple[str, ...]                  # labels (for diagnostic printing)
    point: tuple[float, float]              # (sigma_t*R, rho)
    # Derived flags:
    closure_monotone_down: bool
    closure_sign_stable: bool
    closure_beats_f4_every_quad: bool
    f4_monotone_down: bool     # false => F.4 is in its own quadrature-noise regime
    f4_sign_stable: bool       # false => F.4 itself has a sign flip (L17 (20,0.5))


def assert_rank_n_structural_win(
    closure_fn: Callable[[float, float, dict], float],
    f4_fn: Callable[[float, float, dict], float],
    point: tuple[float, float],
    quads: Sequence[tuple[str, dict]],
) -> StabilityReport:
    """Protocol gate: closure beats F.4 at ``point`` across ALL ``quads``.

    Parameters
    ----------
    closure_fn : callable
        ``closure_fn(tau, rho, quad_dict) -> k_eff``. ``quad_dict`` is a
        ``{"n_panels", "p_order", "n_ang"}`` triple.
    f4_fn : callable
        Same signature, running the production F.4 scalar closure.
    point : tuple
        (sigma_t * R, rho = r_0 / R).
    quads : sequence of (label, quad_dict)
        At least two quadratures in order of INCREASING refinement.
        Typical: RICH, RICH+panels, RICH+pp, ULTRA.

    Returns
    -------
    StabilityReport

    Raises
    ------
    AssertionError
        If any of the five L19 stability conditions fail:
          (i) len(quads) >= 2
          (ii) closure |err| < F.4 |err| at every quadrature
          (iii) closure signed err does not flip sign across quads
          (iv) closure |err| is monotone non-increasing across quads
          (v) F.4 |err| itself is monotone non-increasing (otherwise the
              reference is in F.4's own quadrature-noise regime and the
              claim is unverifiable at ``quads``)

    Notes
    -----
    Both ``closure_fn`` and ``f4_fn`` must return k_eff so the helper can
    compute the SIGNED residual (k - K_INF) / K_INF. Returning |err|
    directly would lose the sign flip signature.
    """
    if len(quads) < 2:
        raise AssertionError(
            "assert_rank_n_structural_win requires >= 2 quadratures "
            "(L17/L19 -- one-quadrature RICH-wins have proven to be "
            "crossing artifacts). Got "
            f"{len(quads)}: {[q[0] for q in quads]}."
        )

    tau, rho = point
    labels = tuple(lbl for (lbl, _) in quads)

    closure_signed = []
    f4_signed = []
    for (_lbl, qd) in quads:
        k_c = closure_fn(tau, rho, qd)
        k_f = f4_fn(tau, rho, qd)
        closure_signed.append((k_c - K_INF) / K_INF)
        f4_signed.append((k_f - K_INF) / K_INF)
    closure_signed = tuple(closure_signed)
    f4_signed = tuple(f4_signed)

    # Derived flags
    closure_abs = [abs(x) for x in closure_signed]
    f4_abs = [abs(x) for x in f4_signed]
    closure_monotone_down = all(
        closure_abs[i + 1] <= closure_abs[i] for i in range(len(closure_abs) - 1)
    )
    f4_monotone_down = all(
        f4_abs[i + 1] <= f4_abs[i] for i in range(len(f4_abs) - 1)
    )

    closure_sign_stable = (
        all(x >= 0 for x in closure_signed)
        or all(x <= 0 for x in closure_signed)
    )
    f4_sign_stable = (
        all(x >= 0 for x in f4_signed)
        or all(x <= 0 for x in f4_signed)
    )

    closure_beats_f4_every_quad = all(
        closure_abs[i] < f4_abs[i] for i in range(len(quads))
    )

    report = StabilityReport(
        closure_signed=closure_signed,
        f4_signed=f4_signed,
        quads=labels,
        point=(tau, rho),
        closure_monotone_down=closure_monotone_down,
        closure_sign_stable=closure_sign_stable,
        closure_beats_f4_every_quad=closure_beats_f4_every_quad,
        f4_monotone_down=f4_monotone_down,
        f4_sign_stable=f4_sign_stable,
    )

    # L19 assertions -- raise on failure with a diagnostic message.
    diag = (
        f"\n  point = (sigma_t*R = {tau}, rho = {rho})\n"
        f"  quads = {labels}\n"
        f"  closure signed = {tuple(f'{x*100:+.6f}%' for x in closure_signed)}\n"
        f"  F.4     signed = {tuple(f'{x*100:+.6f}%' for x in f4_signed)}"
    )
    if not closure_beats_f4_every_quad:
        raise AssertionError(
            "L19 FAIL: closure does NOT beat F.4 at every quadrature"
            + diag
            + f"\n  per-quad |closure| < |F.4|: "
            + str([closure_abs[i] < f4_abs[i] for i in range(len(quads))])
        )
    if not closure_sign_stable:
        raise AssertionError(
            "L19 FAIL: closure signed error FLIPS SIGN across quadratures "
            "-- crossing artifact, not a structural win"
            + diag
        )
    if not closure_monotone_down:
        raise AssertionError(
            "L19 FAIL: closure |err| is NOT monotone non-increasing across "
            "quadratures -- unconverged, not a structural win"
            + diag
        )
    if not f4_monotone_down:
        raise AssertionError(
            "L19 FAIL: F.4 |err| itself is NOT monotone non-increasing across "
            "quadratures -- reference is in F.4's own quadrature-noise "
            "regime, claim is unverifiable at these quadratures. "
            "Include a richer reference quadrature."
            + diag
        )
    if not f4_sign_stable:
        raise AssertionError(
            "L19 FAIL: F.4 signed err itself FLIPS SIGN across quadratures "
            "-- reference is in F.4's own quadrature-noise regime (L17 "
            "(20,0.5) signature), claim is unverifiable at these quadratures. "
            "Include a richer reference quadrature."
            + diag
        )
    return report


# =========================================================================
# Helper unit tests -- verify the gate passes a trivial case and catches
# synthetic failure modes.
# =========================================================================


@pytest.mark.l1
@pytest.mark.verifies("peierls-rank-n-stability")
def test_protocol_trivial_pass_f4_vs_itself_must_fail_beat_check():
    """F.4 compared to itself: closure never STRICTLY beats F.4, so the
    beat-every-quad check must fire. This pins the helper's strictness --
    a closure that ties F.4 is NOT a win."""
    series = [1e-3, 5e-4, 2e-4, 5e-5]

    def mkfn(vals):
        i = [0]
        def fn(tau, rho, qd):
            v = vals[i[0]]
            i[0] += 1
            return K_INF * (1.0 + v)
        return fn

    with pytest.raises(AssertionError, match="does NOT beat F.4"):
        assert_rank_n_structural_win(
            closure_fn=mkfn(series),
            f4_fn=mkfn(series),
            point=(10.0, 0.3),
            quads=[("A", {}), ("B", {}), ("C", {}), ("D", {})],
        )


@pytest.mark.l1
@pytest.mark.verifies("peierls-rank-n-stability")
def test_protocol_trivial_pass_closure_strictly_better():
    """Synthetic closure strictly better than synthetic F.4: every L19 flag
    should be True. The helper must return a report, not raise."""
    f4_series = [1e-3, 8e-4, 5e-4, 3e-4]
    closure_series = [5e-4, 3e-4, 1e-4, 5e-5]

    def mkfn(vals):
        i = [0]
        def fn(tau, rho, qd):
            v = vals[i[0]]
            i[0] += 1
            return K_INF * (1.0 + v)
        return fn

    rep = assert_rank_n_structural_win(
        closure_fn=mkfn(closure_series),
        f4_fn=mkfn(f4_series),
        point=(10.0, 0.3),
        quads=[("A", {}), ("B", {}), ("C", {}), ("D", {})],
    )
    assert rep.closure_beats_f4_every_quad
    assert rep.closure_sign_stable
    assert rep.closure_monotone_down
    assert rep.f4_monotone_down


@pytest.mark.l1
@pytest.mark.verifies("peierls-rank-n-stability")
def test_protocol_catches_sign_flip():
    """A closure that wins at RICH but flips sign at RICH+panels (the
    (5.0, 0.5) PCA case from Direction C) MUST raise AssertionError."""
    # F.4: small positive, monotone down in magnitude.
    f4_series = [+3.0e-3, +2.5e-3]
    # Closure: wins at RICH in |err|, but signed err flips sign at quad 2.
    closure_series = [+2.8e-4, -7.6e-4]  # PCA (5.0, 0.5) @ alpha*=1.0671

    def mkfn(vals):
        i = [0]
        def fn(tau, rho, qd):
            v = vals[i[0]]
            i[0] += 1
            return K_INF * (1.0 + v)
        return fn

    with pytest.raises(AssertionError, match="FLIPS SIGN"):
        assert_rank_n_structural_win(
            closure_fn=mkfn(closure_series),
            f4_fn=mkfn(f4_series),
            point=(5.0, 0.5),
            quads=[("RICH", {}), ("RICH+panels", {})],
        )


@pytest.mark.l1
@pytest.mark.verifies("peierls-rank-n-stability")
def test_protocol_catches_non_monotone():
    """A closure whose |err| grows at a richer quadrature must fail."""
    f4_series = [+3.0e-3, +2.5e-3, +2.0e-3]
    # Closure wins at each quad individually, but |err| grows between
    # quads 1 and 2 -- textbook unconverged.
    closure_series = [+5.0e-4, +1.0e-3, +1.5e-3]

    def mkfn(vals):
        i = [0]
        def fn(tau, rho, qd):
            v = vals[i[0]]
            i[0] += 1
            return K_INF * (1.0 + v)
        return fn

    with pytest.raises(AssertionError, match="monotone"):
        assert_rank_n_structural_win(
            closure_fn=mkfn(closure_series),
            f4_fn=mkfn(f4_series),
            point=(5.0, 0.3),
            quads=[("A", {}), ("B", {}), ("C", {})],
        )


@pytest.mark.l1
@pytest.mark.verifies("peierls-rank-n-stability")
def test_protocol_catches_non_beat_any_quad():
    """A closure that loses at ANY quadrature must fail."""
    f4_series = [+3.0e-3, +2.5e-3]
    # Closure loses at quad 1 even though it wins at quad 0.
    closure_series = [+1.0e-3, +3.0e-3]

    def mkfn(vals):
        i = [0]
        def fn(tau, rho, qd):
            v = vals[i[0]]
            i[0] += 1
            return K_INF * (1.0 + v)
        return fn

    with pytest.raises(AssertionError, match="does NOT beat F.4"):
        assert_rank_n_structural_win(
            closure_fn=mkfn(closure_series),
            f4_fn=mkfn(f4_series),
            point=(5.0, 0.3),
            quads=[("A", {}), ("B", {})],
        )


@pytest.mark.l1
@pytest.mark.verifies("peierls-rank-n-stability")
def test_protocol_catches_unresolved_f4_reference_sign_flip():
    """F.4 itself sign-flips between RICH and RICH+panels (L17, (20,0.5))
    -- the reference is in F.4's quadrature-noise regime, no claim is
    verifiable here."""
    # Signed trajectory mirrors L17 (20, 0.5): +0.0054% -> -0.0039%.
    # In magnitude this IS monotone down, so the magnitude check does not
    # fire -- only the SIGN check catches it.
    f4_series = [+5.4e-5, -3.9e-5]
    # Closure is carefully constructed to look good -- monotone down,
    # sign-stable, strictly smaller than F.4 in magnitude -- so ONLY the
    # F.4-sign-flip check can raise.
    closure_series = [+1.0e-5, +5.0e-6]

    def mkfn(vals):
        i = [0]
        def fn(tau, rho, qd):
            v = vals[i[0]]
            i[0] += 1
            return K_INF * (1.0 + v)
        return fn

    with pytest.raises(AssertionError, match="F.4 signed err itself FLIPS SIGN"):
        assert_rank_n_structural_win(
            closure_fn=mkfn(closure_series),
            f4_fn=mkfn(f4_series),
            point=(20.0, 0.5),
            quads=[("RICH", {}), ("RICH+panels", {})],
        )


@pytest.mark.l1
@pytest.mark.verifies("peierls-rank-n-stability")
def test_protocol_catches_unresolved_f4_reference_magnitude_grows():
    """F.4 magnitude grows under refinement -- unresolved."""
    f4_series = [+1.0e-5, +3.0e-5]   # magnitude grew -> not monotone down
    closure_series = [+1.0e-6, +5.0e-7]  # looks great

    def mkfn(vals):
        i = [0]
        def fn(tau, rho, qd):
            v = vals[i[0]]
            i[0] += 1
            return K_INF * (1.0 + v)
        return fn

    with pytest.raises(AssertionError, match="F.4 .* NOT monotone"):
        assert_rank_n_structural_win(
            closure_fn=mkfn(closure_series),
            f4_fn=mkfn(f4_series),
            point=(20.0, 0.3),
            quads=[("RICH", {}), ("RICH+panels", {})],
        )


@pytest.mark.l1
@pytest.mark.verifies("peierls-rank-n-stability")
def test_protocol_requires_at_least_two_quads():
    """One-quadrature claims are not even evaluable by L19."""
    with pytest.raises(AssertionError, match=">= 2 quadratures"):
        assert_rank_n_structural_win(
            closure_fn=lambda t, r, q: K_INF,
            f4_fn=lambda t, r, q: K_INF,
            point=(5.0, 0.3),
            quads=[("RICH", {})],
        )


# =========================================================================
# F.4 sign-stability baseline -- anchors L17/L19 on the canonical grid.
# =========================================================================
#
# Per-point reference quadratures are populated below from the
# ``diag_f4_structural_floor_baseline.py`` scan output. The pinning rule
# (shared with the scanner):
#     - The reference quadrature Q_ref is the earliest quadrature in the
#       ladder (RICH, RICH+panels, RICH+pp, ULTRA) such that
#         |signed(Q_ref)|         < 0.005%,
#         |signed(Q_ref)|         < |signed(Q_prev)|,  and
#         sign(signed(Q_ref))     == sign(signed(Q_prev)).
#     - If no such quadrature exists within a 120 s/run budget on the
#       container CI machine, the point is flagged "unresolved" and the
#       baseline test for it is SKIPPED with an informative reason.
#
# Each entry: (tau, rho): (ref_quad_label, prev_quad_label, status_note).
# ``status_note`` is "resolved" or a skip reason string.
#
# Updated 2026-04-22 from the baseline scan on the devcontainer.
# See .claude/agent-memory/numerics-investigator/direction_n_quadrature_baseline.md
# for the full signed-err trajectory at each point.
#

# Quadrature definitions (must match the scanner).
QUAD = {
    "RICH":        dict(n_panels=4, p_order=8,  n_ang=64),
    "RICH+panels": dict(n_panels=5, p_order=8,  n_ang=64),
    "RICH+pp":     dict(n_panels=5, p_order=10, n_ang=64),
    "ULTRA":       dict(n_panels=5, p_order=10, n_ang=96),
}


# Per-point reference state from the 2026-04-22 scan on the devcontainer
# (diag_f4_structural_floor_baseline.py). All six points REMAIN UNRESOLVED at
# the 120 s/run budget on this hardware: RICH+pp and ULTRA both time out at
# every point. The accessible pair RICH vs RICH+panels is documented here
# with the observed signed-err trajectory. No point meets the L19 resolution
# criteria (|err| < 0.005% AND monotone DOWN AND sign-stable), so the
# parametrized baseline test SKIPS every point with an informative reason.
#
# Trajectory legend for the `trajectory` field:
#   "monotone_down"    -- magnitude shrinks, sign stable (approaching floor)
#   "magnitude_grew"   -- magnitude larger at RICH+panels (quadrature-noise regime)
#   "sign_flip"        -- signed err flips sign (L17 canonical signature)
#
# See ``direction_n_quadrature_baseline.md`` for the full signed-err trajectory
# and the hardware cost to resolve each point (extrapolated from the
# observed scaling: RICH+pp ~ 180-240 s/run on this box; ULTRA ~ 300-400 s).
F4_REFERENCE_BASELINE: dict[tuple[float, float], dict] = {
    (5.0, 0.3): {
        "ref": None, "prev": None,
        "signed_rich": +5.7825e-4,
        "signed_rich_panels": +3.2331e-4,
        "trajectory": "monotone_down",
        "status": (
            "unresolved at 120s/run devcontainer budget: RICH+pp timed out; "
            "RICH/RICH+panels are +0.0578%/+0.0323% (monotone down, sign-stable) "
            "but still above 0.005% -- need p_order=10 or ULTRA to resolve"
        ),
    },
    (5.0, 0.5): {
        "ref": None, "prev": None,
        "signed_rich": +2.97309e-3,
        "signed_rich_panels": +3.22013e-3,
        "trajectory": "magnitude_grew",
        "status": (
            "unresolved: RICH +0.2973% -> RICH+panels +0.3220% (|err| GREW). "
            "One-panel refinement in the wrong direction means this point is "
            "not even on the trailing edge of F.4 convergence at RICH -- "
            "needs p_order=10 + more panels"
        ),
    },
    (10.0, 0.3): {
        "ref": None, "prev": None,
        "signed_rich": +3.285e-5,
        "signed_rich_panels": -8.221e-5,
        "trajectory": "sign_flip",
        "status": (
            "unresolved: RICH +0.00329% -> RICH+panels -0.00822% -- L17 canonical "
            "SIGN FLIP. F.4 structural floor here is below BOTH; historically "
            "F.4 ULTRA=0.0010% here (3x below RICH). Need ULTRA (5,10,96) as "
            "reference quadrature; requires > 120s/run on this hardware"
        ),
    },
    (10.0, 0.5): {
        "ref": None, "prev": None,
        "signed_rich": +1.6683e-4,
        "signed_rich_panels": +7.990e-5,
        "trajectory": "monotone_down",
        "status": (
            "unresolved: RICH +0.0167% -> RICH+panels +0.0080%. Monotone down "
            "and sign-stable but |err| still > 0.005%. Trend suggests one more "
            "refinement (RICH+pp or ULTRA) would cross the threshold; need the "
            "RICH+pp run to pin it, which did not fit in the 120s budget"
        ),
    },
    (20.0, 0.3): {
        "ref": None, "prev": None,
        "signed_rich": -5.982e-5,
        "signed_rich_panels": -7.481e-5,
        "trajectory": "magnitude_grew",
        "status": (
            "unresolved: RICH -0.00598% -> RICH+panels -0.00748% (|err| GREW "
            "marginally, sign stable). Both just above 0.005% threshold. Needs "
            "p_order=10 to separate structural from quadrature residual"
        ),
    },
    (20.0, 0.5): {
        "ref": None, "prev": None,
        "signed_rich": +5.430e-5,
        "signed_rich_panels": -3.942e-5,
        "trajectory": "sign_flip",
        "status": (
            "unresolved: RICH +0.00543% -> RICH+panels -0.00394% -- the L17 "
            "canonical (20,0.5) SIGN FLIP, reproduced exactly. F.4 structural "
            "floor here is BELOW RICH+panels. Needs ULTRA (5,10,96) at least; "
            "requires > 120s/run"
        ),
    },
}


# =========================================================================
# Subprocess worker to evaluate F.4 without binding the test process.
# Mirrors diag_f4_structural_floor_baseline.py.
# =========================================================================


_WORKER_CODE = textwrap.dedent(
    r"""
    import json, sys, time
    sys.path.insert(0, '/workspaces/ORPHEUS/scratch/derivations/diagnostics')
    from diag_cin_aware_split_basis_keff import run_scalar_f4
    tau = float(sys.argv[1]); rho = float(sys.argv[2])
    n_panels = int(sys.argv[3]); p_order = int(sys.argv[4]); n_ang = int(sys.argv[5])
    R = tau; r_0 = rho * R
    t0 = time.time()
    k = run_scalar_f4(r_0, R, 1.0, 1.0/3.0, 1.0,
                      n_panels=n_panels, p_order=p_order, n_ang=n_ang)
    dt = time.time() - t0
    print(json.dumps({"k": k, "wall": dt}))
    """
).strip()


def _run_f4_subprocess(tau: float, rho: float, quad: dict, timeout: float = 120.0) -> float:
    args = [
        sys.executable, "-c", _WORKER_CODE,
        str(tau), str(rho),
        str(quad["n_panels"]), str(quad["p_order"]), str(quad["n_ang"]),
    ]
    proc = subprocess.run(
        args, capture_output=True, text=True, timeout=timeout, check=True,
    )
    out = json.loads(proc.stdout.strip().splitlines()[-1])
    return out["k"]


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("peierls-rank-n-stability")
@pytest.mark.parametrize("tau,rho", list(F4_REFERENCE_BASELINE.keys()))
def test_f4_is_sign_stable_at_its_reference_quadrature(tau: float, rho: float):
    """L17/L19 baseline: at each point of the canonical 6-point grid,
    F.4's signed error at the documented reference quadrature must be

      (a) below 0.005 % in magnitude,
      (b) smaller in magnitude than at the immediately-coarser quadrature, and
      (c) the same sign as at that coarser quadrature.

    This is the ground truth that any future rank-N closure competes
    against. If this test fails, either the reference quadrature pinned
    by ``diag_f4_structural_floor_baseline.py`` is stale or F.4 itself
    has drifted.

    Points that remain unresolved at the 120 s/run budget on the scan
    machine are SKIPPED with an informative reason. On the 2026-04-22
    devcontainer scan ALL SIX POINTS skip -- RICH+pp and ULTRA both
    time out. See ``direction_n_quadrature_baseline.md`` for the full
    signed-err trajectory and the extrapolated hardware cost to resolve
    them.
    """
    rec = F4_REFERENCE_BASELINE[(tau, rho)]
    if rec["status"] != "resolved":
        pytest.skip(
            f"F.4 baseline at (sigma_t*R = {tau}, rho = {rho}): "
            f"{rec['status']}"
        )
    ref_label = rec["ref"]
    prev_label = rec["prev"]
    ref_quad = QUAD[ref_label]
    prev_quad = QUAD[prev_label]

    k_prev = _run_f4_subprocess(tau, rho, prev_quad)
    k_ref = _run_f4_subprocess(tau, rho, ref_quad)

    signed_prev = (k_prev - K_INF) / K_INF
    signed_ref = (k_ref - K_INF) / K_INF

    assert abs(signed_ref) < 5.0e-5, (
        f"F.4 at ({tau}, {rho}) {ref_label}: |signed err| = "
        f"{abs(signed_ref)*100:.6f}% NOT below 0.005%. Reference quadrature "
        f"pin may be stale."
    )
    assert abs(signed_ref) < abs(signed_prev), (
        f"F.4 at ({tau}, {rho}): magnitude GREW from "
        f"{abs(signed_prev)*100:.6f}% ({prev_label}) to "
        f"{abs(signed_ref)*100:.6f}% ({ref_label}). Not monotone -- this "
        "point is still in F.4's quadrature-noise regime at the documented "
        "reference; re-pin to a richer quadrature."
    )
    assert (signed_prev > 0) == (signed_ref > 0) or signed_prev == 0, (
        f"F.4 at ({tau}, {rho}): signed err FLIPPED SIGN from "
        f"{signed_prev*100:+.6f}% ({prev_label}) to "
        f"{signed_ref*100:+.6f}% ({ref_label}). This is the L17 "
        "crossing signature -- the documented reference quadrature is below "
        "F.4's own structural floor at this point. Re-pin richer."
    )


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("peierls-rank-n-stability")
@pytest.mark.parametrize("tau,rho", list(F4_REFERENCE_BASELINE.keys()))
def test_f4_rich_vs_rich_panels_matches_pinned_baseline(tau: float, rho: float):
    """Pin the F.4 RICH vs RICH+panels trajectory observed on 2026-04-22.

    Until richer reference quadratures are affordable on CI, this test
    IS the L17/L19 evidence -- it captures the per-point signed error
    at the two quadratures that do fit in the 120 s/run budget, and
    asserts that F.4 still reproduces them within a loose tolerance.

    Tolerance: 1e-6 (absolute, on the signed fractional residual).
    This is ~2 orders of magnitude below the smallest measured signed
    err (0.00329% = 3.3e-5), so a drift that flips the sign flip
    signature or dislodges the magnitude ordering will trip the test.

    Points whose trajectory is ``sign_flip`` encode L17 directly:
    (10.0, 0.3) and (20.0, 0.5). If a future change removes either
    sign flip by shrinking the F.4 RICH residual to zero, the closure
    has moved its structural floor below RICH -- in which case the
    whole baseline must be re-scanned at richer quadratures.
    """
    rec = F4_REFERENCE_BASELINE[(tau, rho)]
    expected_rich = rec["signed_rich"]
    expected_rich_panels = rec["signed_rich_panels"]

    k_rich = _run_f4_subprocess(tau, rho, QUAD["RICH"])
    k_rich_panels = _run_f4_subprocess(tau, rho, QUAD["RICH+panels"])

    signed_rich = (k_rich - K_INF) / K_INF
    signed_rich_panels = (k_rich_panels - K_INF) / K_INF

    tol = 1.0e-6  # absolute, on the signed fractional residual
    assert abs(signed_rich - expected_rich) < tol, (
        f"F.4 RICH at ({tau}, {rho}): observed {signed_rich*100:+.6f}% vs "
        f"pinned {expected_rich*100:+.6f}% (delta {(signed_rich - expected_rich)*1e6:+.3f} ppm). "
        "F.4 has drifted -- re-scan the baseline."
    )
    assert abs(signed_rich_panels - expected_rich_panels) < tol, (
        f"F.4 RICH+panels at ({tau}, {rho}): observed {signed_rich_panels*100:+.6f}% "
        f"vs pinned {expected_rich_panels*100:+.6f}% "
        f"(delta {(signed_rich_panels - expected_rich_panels)*1e6:+.3f} ppm). "
        "F.4 has drifted -- re-scan the baseline."
    )
