r"""Bit-identity gate for the #208 affine flux-algebra carve.

The affine typing (Piece 1: ``FluxDisplacement`` + the ``flux+flux`` gate),
the SI-loop displacement retype (Piece 2), and the typed equation-residual
evaluation (Piece 3) add **ZERO numerical change** to the converged flux:

* The SI stopping norm STAYS the flat ``_l2_norm(displacement.to_flat())``
  (the displacement composite's ``to_flat()`` reads the SAME ``.values`` as the
  pre-carve ``(psi-psi_prev).to_flat()`` — see the displacement mechanism in
  ``.claude/plans/issue_208_residual_typing_closeout.md`` §3a).
* The residual evaluation is additive / diagnostic — never in the convergence
  path.

That claim was verified at the pre-carve commit ``63719a2`` by freezing a
``sha256`` of the converged ``angular_flux.interior.values`` and
``scalar_flux.values``.

⛔ #333 re-pose (2026-08-12) — the sha256 instrument is RETIRED
===============================================================

**#208's zero-numerical-change claim above is now HISTORICAL.** It is
past-tense record, not a live assertion: it was true and verified at
``63719a2``, and this module no longer re-verifies it.

The ``sha256`` was the sharpest instrument available for a zero-change claim,
and the choice was defensible when it was made (see "Why a dedicated golden"
below).  Its defect is structural rather than a mistake: a zero-change claim
pinned against a FROZEN PAST has a shelf life.  It is falsified by the first
legitimate change anywhere upstream, and at that moment the instrument gives
no way to decide whether the new value is fine — **the magnitude cannot be
computed, because the old VALUES were never stored.  1 ULP and a catastrophic
error are the same red.**

That is exactly what happened.  Four quadrature commits moved these values,
every one of them a verified improvement:

* ``d6f53afe`` — LS derived nodes, 16 of 24 ordinates, **1 ULP** (rules now
  IMPOSE their symmetry, so derived ordinates are bit-copies of the seed
  octant).  Reddens the two 2-D arms.
* ``579d5eaf`` — ``gauss_legendre`` nodes AND weights, 3 ULP / 17 ULP.
  Reddens the slab arm.  [M] the new GL8 integrates the exact rational
  moments **5.8× better** than ``numpy.leggauss``, winning 7 of 7
  non-trivial even moments.
* ``df33913d`` (#327) / ``59bb38a0`` (#337) — the LS weights and the μ₁ seed,
  replacing a self-declared house convention (``1/√6``) with the
  Lewis & Miller / LA-3186 moment-matched value.

The re-pose stores VALUES and budgets against the solver's own stopping
criterion.  **It IS a weakening** — see the test's own docstring, which states
what is and is not asserted now.

⚠ The predecessor was ALSO re-baselined three times before this (the
regeneration history below), each with a written justification.  That history
is kept deliberately: it records why each move was legitimate, and that
reasoning outlives the instrument it was written for.

Why a dedicated golden, not the DD regression snapshots (the ORIGINAL rationale)
================================================================================

``tests/sn/regression/test_dd_regression.py`` compares iterative results at
``SAFETY × conv_tol ≈ 1e-11`` (escalatable to strict only via
``-W error::DriftWarning``), and the ``2d_2g_p1_aniso_dd_8x4_het_si`` snapshot
ALREADY pre-drifts ~6920 ULP / ~9.8e-13 (Phase-5b/5c inheritance) — so the
snapshot suite CANNOT verify this carve's stronger zero-change claim. A
``sha256`` of the raw bytes was the sharpest (sub-ULP) bit-identity assertion.

⚠ Read honestly, that rationale now argues the re-posed gate overlaps the DD
suite's tolerance. What it still adds: three DIFFERENT drivers (2-D SI, 2-D
Krylov, 1-D slab SI) pinned with stored values and an any-movement
``DriftWarning`` tripwire, which the DD suite does not cover for these configs.

Coverage (Cardinal-6 ≥2G + heterogeneous; vv-principles §H1/H2)
===============================================================

* ``si_2d_p1_aniso_het`` — 2-D Cartesian, 2G, fuel|moderator, P1-anisotropic,
  source-iteration: the **windowed SI** path whose iterate ``bulk`` is a
  :class:`HarmonicMomentFlux` (so ``psi-psi_prev`` becomes a
  ``MomentDisplacement``). The mirror of the frozen snapshot config.
* ``krylov_2d_p1_aniso_het`` — same config, full-angular Krylov (scipy's flat
  ``b-Ax`` — exercises the typing without the SI displacement path).
* ``si_slab_2g_het`` — 1-D slab, 2G, fuel|moderator, P1, source-iteration: the
  ``AngularFlux``-bulk SI path (1-D never windows → ``AngularDisplacement``).

Marks: ``foundation`` — this is a software-invariant gate (the converged value
is inherited from the pre-carve reference), not a physics ``:label:`` claim.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D
from orpheus.geometry.mesh import Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import solve_sn_fixed_source
from tests.sn._test_helpers import SN_TESTS_ROOT
from tests.sn.regression._regression_assert import assert_regression


pytestmark = pytest.mark.foundation


# Generated at pre-carve HEAD 63719a2 (commit with Piece 1). sha256 of the
# C-contiguous float64 bytes of the converged fields.
#
# REGENERATION HISTORY:
# * 2026-06-11 (S6.9 Fork-B2, #222) — the FOUR 2-D hashes regenerated in the
#   default-flip commit: the multi-D Cartesian production default changed
#   MovingFrontierWindow → ScanMarch, a SCHEDULE change (row-march vs
#   anti-diagonal) that shifts the converged bytes at FP-association level
#   (principled-equivalent, NOT a numerics change — vv §bit-identity-vs-
#   principled).  Output-identity evidence for the flip: the G4.a/G4.b
#   Mode-9 FP-invariance gates (test_scan_march_end_to_end.py, scanmarch
#   default ≡ forced window at solver tol) + the ScanMarch G2.c nulp oracle.
#   The 1-D slab hashes are UNTOUCHED (CumprodScan stays the 1-D default —
#   the flip's blast radius pin).
# * 2026-06-15 (#240 Phase 2 Step B) — the TWO ``krylov_2d_p1_aniso_het`` hashes
#   regenerated. The Step-B carve made ``StreamingCollisionOperator.apply`` OWN its
#   matvec (``loss_action(self.sigma)`` direct) instead of the inherited leaf
#   sum ``(loss_action(σ_t) − σ_t·ψ) + σ_t·ψ`` — the override DROPS that
#   ``−σ_t·ψ + σ_t·ψ`` round-trip, so the 2-D Cartesian matvec re-associates
#   ≤5 ULP per application; accumulated through GMRES (inner_tol=1e-12) the
#   converged Krylov bytes shift at the ~1e-12 GMRES-tolerance level. Verified
#   structurally-independent: the Krylov converged φ agrees with the SI 2-D φ
#   (which does NOT use the apply override — bit-identical to ITS golden) to
#   3.9e-12 rel (vv-principles 3-criteria: named ``loss_action`` / SI cross-
#   check / drift = GMRES-tol × ULP). The SI 2-D + slab hashes are UNTOUCHED
#   (SI rides ``solve``, not ``apply`` — the carve's apply-only blast-radius
#   pin).
# * 2026-06-16 (#240 Phase 2 Step D5a) — ALL FOUR 2-D hashes regenerated (BOTH
#   ``si_2d_p1_aniso_het`` AND ``krylov_2d_p1_aniso_het``). D5a folded the 2-D
#   row-march ``ScanMarch._sweep_interior`` (SOLVE) + ``_loss_action_interior``
#   (matvec) off inline DiamondDifference onto the scheme's coefficient model
#   (#239): the SOLVE now consumes ``scheme.cartesian_scan_coefficients`` →
#   ``(a, inverse_denom, w)`` + ``source_emission``/``cell_average`` (the
#   ``×inverse_denom`` reciprocal form, replacing the ONE remaining inline ``÷D``
#   DIVISION — so this is the same byte re-association the 1-D CumprodScan
#   already rides, NOT a numerics change); the matvec rides
#   ``scheme.residual_kernel_batch`` + ``reflect_scan_coefficients``. Both folds
#   re-associate the 2-D cell-balance reduction ~1 ULP. The SI shift accumulates
#   through source iteration (iter_count × ULP); the Krylov shift through GMRES
#   (GMRES-tol × ULP). Verified structurally-independent: the post-fold SI 2-D φ
#   (``.solve`` fold) and Krylov 2-D φ (``.apply`` fold) — DIFFERENT code paths —
#   agree to 4.2e-12 rel, AND the ScanMarch≡FullFieldWavefront oracle
#   (test_scan_march_equivalence.py) pins the value to the analytical
#   ``k_inf``/``φ=Q/Σ_t`` grounds (vv-principles §"Bit-identity vs principled-
#   equivalence" 3-criteria: named coefficient-model ops / two-path + oracle
#   cross-check / drift = iter-or-GMRES-tol × ULP). The 1-D slab hashes are
#   UNTOUCHED — D5a's blast radius is the 2-D row-march ONLY (the 1-D CumprodScan
#   scan path is byte-identical), which is exactly the D5a.3 negative-control
#   proof that the fold did not leak into the 1-D solve.
# * 2026-07-12 (step 5 #41, R-5.2/R-5.3) — the THREE SI hash pairs
#   regenerated (si_2d_p1_aniso_het + si_slab_2g_het). The SourceIteration
#   stop re-posed onto the ρ-honest equation residual ‖Σg·Δψ‖/‖q_ext‖
#   (was the iterate increment ‖Δψ‖/‖ψ‖) — a DELIBERATE tol
#   re-interpretation, so the SI converged bytes shift at the
#   iteration-count level (the same fixed point, reached at a different
#   stopping surface; vv §bit-identity-vs-principled). Verified
#   structurally-independent: the post-change SI 2-D φ agrees with the
#   Krylov 2-D φ (a DIFFERENT driver whose GMRES stop is UNTOUCHED) to
#   9.8e-13 rel at inner_tol=1e-12, and the closed-form kinf/Q-over-Σt
#   anchor walls held green through the change. The KRYLOV hashes are
#   UNTOUCHED — the stop change's blast radius is SI-only, which is
#   exactly the pin (GMRES was already residual-stopped).
# * 2026-08-12 (#333) — THE sha256 ERA ENDS HERE. The six digests above were
#   deleted and the gate re-posed onto stored VALUES; see the "#333 re-pose"
#   section of the module docstring. The history above is kept because it
#   records WHY each regeneration was legitimate — that reasoning outlives the
#   instrument it was written for.

# ─────────────────────────────────────────────────────────────────────
# #333 — stored VALUES, not digests.  Reuses the ROOT conftest's
# ``--capture-baseline`` flag (the same switch
# ``tests/sn/sweep/core/test_affine_carve_baseline.py`` uses) rather than
# minting a second capture mechanism.
# ─────────────────────────────────────────────────────────────────────
_BASELINE_DIR = SN_TESTS_ROOT / "_data" / "affine_carve_converged"

#: The solver's OWN stopping criterion for these runs.  ``assert_regression``
#: requires ``conv_tol`` to be read off the run config, never hardcoded at the
#: call site — so the solves below and the gate share this one name.
_INNER_TOL = 1e-12
_MAX_INNER = 3000


def _capturing(request) -> bool:
    return bool(request.config.getoption("--capture-baseline", default=False))


def _baseline_path(case: str, quantity: str) -> Path:
    return _BASELINE_DIR / f"{case}_{quantity}.npy"


def _capture_or_assert(request, case: str, quantity: str, actual) -> bool:
    """WRITE under ``--capture-baseline``, else READ + assert.

    Returns ``True`` if it WROTE.  The caller skips ONCE at the end when
    any leg captured — calling ``pytest.skip`` in here would short-circuit
    the remaining legs (the bug a multi-array gate must avoid).
    """
    arr = np.ascontiguousarray(np.asarray(actual, dtype=np.float64))
    path = _baseline_path(case, quantity)
    if _capturing(request):
        _BASELINE_DIR.mkdir(parents=True, exist_ok=True)
        np.save(path, arr)
        return True
    if not path.exists():
        pytest.fail(
            f"missing baseline {path}; run with --capture-baseline to write it."
        )
    assert_regression(
        arr,
        np.load(path),
        conv_tol=_INNER_TOL,
        case_name=f"{case}:{quantity}",
        kind="iterative",
        quantity=quantity,
    )
    return False


def _build_2d():
    """2-D Cartesian 2G het-aniso (mirror of ``2d_2g_p1_aniso_dd_8x4_het_si``)."""
    fuel = get_mixture("A", "2g")
    mod = get_mixture("B", "2g")
    nx, ny = 8, 4
    mat = np.zeros((nx, ny), dtype=int)
    mat[:4, :] = 2
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 1.0, ny + 1),
        mat_map=mat,
        bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    sum_w = float(quad.weights.sum())
    q_ext = np.full((quad.N, 2, nx, ny), 1.0 / sum_w)
    return {2: fuel, 0: mod}, mesh, quad, q_ext


def _build_slab():
    """1-D slab 2G het (fuel|moderator split, vacuum BC)."""
    fuel = get_mixture("A", "2g")
    mod = get_mixture("B", "2g")
    nx = 16
    mat_ids = np.zeros(nx, dtype=int)
    mat_ids[: nx // 2] = 2
    mesh = Mesh1D(
        edges=np.linspace(0.0, 4.0, nx + 1),
        mat_ids=mat_ids,
        bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    sum_w = float(quad.weights.sum())
    q_ext = np.full((quad.N, 2, nx), 1.0 / sum_w)
    return {2: fuel, 0: mod}, mesh, quad, q_ext


def _solve_case(case: str):
    if case == "si_2d_p1_aniso_het":
        mats, mesh, quad, q = _build_2d()
        return solve_sn_fixed_source(
            mats, mesh, quad, q, scattering_order=1,
            inner_solver="source_iteration", max_inner=_MAX_INNER, inner_tol=_INNER_TOL,
        )
    if case == "krylov_2d_p1_aniso_het":
        mats, mesh, quad, q = _build_2d()
        return solve_sn_fixed_source(
            mats, mesh, quad, q, scattering_order=1,
            inner_solver="krylov", max_inner=_MAX_INNER, inner_tol=_INNER_TOL,
        )
    if case == "si_slab_2g_het":
        mats, mesh, quad, q = _build_slab()
        return solve_sn_fixed_source(
            mats, mesh, quad, q, scattering_order=1,
            inner_solver="source_iteration", max_inner=_MAX_INNER, inner_tol=_INNER_TOL,
        )
    raise ValueError(f"unknown case {case!r}")


@pytest.mark.parametrize(
    "case",
    ["si_2d_p1_aniso_het", "krylov_2d_p1_aniso_het", "si_slab_2g_het"],
)
def test_converged_flux_matches_stored_reference(request, case: str) -> None:
    r"""The converged ψ / φ reproduce the stored reference arrays.

    ⛔ **This is deliberately a WEAKER claim than the gate it replaces**, and
    the weakening is irreversible (#333).  The predecessor asserted
    byte-equality against a ``sha256`` frozen at pre-carve ``63719a2``; the
    pre-carve VALUES were never stored, so when a legitimate upstream change
    moved the bytes there was no way to ask *by how much*.  1 ULP and a
    catastrophic error were the same red.  The information needed to do
    better was discarded in 2026-06; no amount of care here recovers it.

    **#208's zero-numerical-change claim is therefore HISTORICAL.** It was
    verified at ``63719a2`` and is recorded in the module docstring and the
    regeneration history above.  It is *not* re-verifiable from this gate and
    this gate no longer asserts it.  What this gate asserts now:

    * the converged values reproduce a **stored reference** — so a red
      reports a magnitude and is diagnosable;
    * hard-failing above the solver's OWN stopping criterion
      (``SAFETY × conv_tol``, ``conv_tol = _INNER_TOL`` read off the run
      config).  A budget tighter than that would assert precision the solver
      never promised; a budget fitted to today's observed drift would assert
      nothing at all (the refuted candidate in #333);
    * a :class:`~tests.sn.regression._regression_assert.DriftWarning`
      tripwire on ANY movement below that floor, so a sub-tolerance drift
      stays AUDIBLE.  ``pytest -W error::DriftWarning`` restores a strict
      bit-identity gate for anyone who wants one.

    The three configs are kept unchanged — the fixture was never the problem
    (2-D 2G P1-aniso het under SI *and* Krylov, plus the 1-D slab whose
    ``AngularFlux``-bulk SI path is the multi-D flip's blast-radius pin).

    Regenerate with ``--capture-baseline``.
    """
    sol = _solve_case(case)
    if not sol.history.converged:
        raise AssertionError(f"{case}: solve did not converge")
    captured = [
        _capture_or_assert(
            request, case, "angular_flux", sol.angular_flux.interior.values,
        ),
        _capture_or_assert(
            request, case, "scalar_flux", sol.scalar_flux.values,
        ),
    ]
    if any(captured):
        pytest.skip(f"{case}: baseline captured (--capture-baseline)")
