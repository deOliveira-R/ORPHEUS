"""DD regression test against frozen reference snapshots.

Each snapshot at ``snapshots/<name>.npz`` carries ``(scalar_flux[,
k_eff])`` captured by ``_generate_snapshots.py``. This test re-runs each
case and asserts the reproduction stays within a **principled** tolerance
— NOT a magic floor and NOT "bit-for-bit".

Why no magic floor.  The legacy gate used ``rtol=1e-12`` / ``atol=1e-13``
on slab and ``rtol=5e-6`` on curvilinear. Those numbers were chosen by
hand to "just clear" the observed drift; they encode no claim about what
the solver promised. The slab cases drift ~1e-11 (FP-non-associativity
past a 1e-12 floor); the curvilinear floor was a wide band papering over
the (now-fixed) curvilinear-k_inf bug. Both are replaced here by a single
defensible rule (see :mod:`._regression_assert`):

* **Correctness gate** = ``SAFETY × conv_tol``, where ``conv_tol`` is the
  solver's OWN convergence stopping criterion for the quantity being
  pinned, read off the run config — k_eff converges on ``keff_tol``, the
  eigen flux on ``flux_tol``, the fixed-source flux on ``inner_tol``.
  Asserting tighter than the solver converged is unphysical.
* **Drift tripwire** = a :class:`._regression_assert.DriftWarning` emitted
  when a value passes the correctness gate but moved beyond bit-identity.
  Informational by default; ``pytest -W error::DriftWarning`` escalates it
  to a hard failure — the strict bit-identity gate a "pure-refactor,
  zero-numerical-change" PR should run.

**Independent corroboration (NOT tautology).**  Regenerating a snapshot
and asserting the solver matches it is circular. Every regenerated value
in this suite was independently corroborated before the snapshot was
trusted (see the commit narrative + ``_generate_snapshots.py``):
homogeneous-reflective k_eff against the closed-form ``k_inf = νΣ_f/Σ_a``
(2G mixture A → 1.875; 1G → 1.5); heterogeneous k_eff against the
reflective-BC particle balance ``production/absorption = k_eff``; the
P1-anisotropic fixed-source flux against the closed-form flat
infinite-medium solution ``(diag(Σ_t) − Σ_{s0}^T)^{-1} Q``; the slab
vacuum fixed-source against global ``source = absorption + leakage``.

**Truncated inners are part of the frozen state** (#340 R2).  Five cases —
``2d_2g_LS4_dd_8x4_het_si``, ``cyl_2g_3reg_folded_4x8_dd_n40``,
``slab_2g_3reg_dd_n40``, ``sphere_2g_3reg_dd_n40``,
``sphere_2g_homogeneous_dd_n20`` — reach ``max_inner`` on at least one
inner before ``inner_tol`` is met, at ρ ≈ 0.93–0.99.  That is not a defect
in the snapshot: the generator and the test share one ``run_case`` /
``run_config`` path, so the *truncated* trajectory is deterministic and is
exactly what the ``.npz`` froze.  `[M]` 2026-08-10 all five reproduce their
snapshot **bit-identically** (``|Δk| = 0.0``, ``max|Δφ| = 0.0``).  Their
``ConvergenceWarning`` is suppressed at the parametrize site, per case —
never at function level, so a SIXTH truncating case announces itself.

⛔ **Raising ``max_inner`` re-baselines this suite, and for two rows it is a
hard failure rather than a drift.**  `[M]` at ``max_inner=2400``:
``sphere_2g_homogeneous_dd_n20`` moves ``|Δk| = 2.45e-11`` and
``slab_2g_3reg_dd_n40`` moves ``|Δk| = 1.81e-11``, both at or past the
correctness gate ``SAFETY × keff_tol = 1e-11``.  The other three stay
inside the correctness gate but every one of the five loses bit-identity,
so ``-W error::DriftWarning`` — the strict gate documented above — would
fail on all five.  Any budget change here is a snapshot regeneration with
the independent-corroboration discipline re-run, never a tolerance nudge.
(``2d_2g_LS4_dd_8x4_het_si`` is the stubborn one: ρ ≈ 0.994 leaves it still
truncated at 2400.)

Tests are skipped if the snapshot file is not yet present — this lets the
regression infrastructure land before the snapshots themselves.

Marker scheme (per ``docs/theory/verification/harness.rst``):

* ``@pytest.mark.regression`` — regression-tier opt-in (CI gate).
* ``@pytest.mark.foundation`` — software-invariant V&V level so the audit
  does not flag the snapshot comparison as an orphan equation.
"""
from __future__ import annotations

import numpy as np
import pytest

from ._generate_snapshots import (
    CASES,
    EIGEN_FLUX_TOL_KEY,
    EIGEN_KEFF_TOL_KEY,
    FIXED_SOURCE_FLUX_TOL_KEY,
    SNAPSHOT_DIR,
    SnapshotCase,
    run_case,
    run_config,
)
from ._regression_assert import assert_regression


pytestmark = [pytest.mark.regression, pytest.mark.foundation]


#: The cases whose snapshot was CAPTURED with a truncated inner (#340 R2) —
#: see the module docstring's "Truncated inners are part of the frozen
#: state".  Not a suppression of convenience: the truncation is IN the
#: frozen artefact, and raising the budget is a re-baseline (a hard
#: correctness failure for two of these five).  Scoped per case rather than
#: per function so that a sixth truncating case starts warning.
_TRUNCATED_INNER_CASES = frozenset({
    "2d_2g_LS4_dd_8x4_het_si",
    "cyl_2g_3reg_folded_4x8_dd_n40",
    "slab_2g_3reg_dd_n40",
    "sphere_2g_3reg_dd_n40",
    "sphere_2g_homogeneous_dd_n20",
})


@pytest.mark.parametrize(
    "case",
    [
        pytest.param(c, marks=pytest.mark.filterwarnings(
            "ignore::orpheus.numerics.convergence.ConvergenceWarning"))
        if c.name in _TRUNCATED_INNER_CASES else c
        for c in CASES
    ],
    ids=lambda c: c.name,
)
def test_dd_regression(case: SnapshotCase) -> None:
    """Re-run case and assert principled agreement with the frozen snapshot.

    Every snapshot case is an iterated solve (eigen power iteration or
    fixed-source source/Krylov iteration), so each quantity is pinned with
    ``kind="iterative"`` at ``SAFETY × conv_tol``. The ``conv_tol`` for
    each quantity is read off the case's run config (the single source of
    truth shared with the generator) — never hardcoded here.

    Eigenvalue cases (``case_kind="eigen"``) pin both ``k_eff``
    (``conv_tol = keff_tol``) and ``scalar_flux`` (``conv_tol = flux_tol``).
    Fixed-source cases (``case_kind="fixed_source"``) pin only
    ``scalar_flux`` (``conv_tol = inner_tol``); there is no eigenvalue for
    the pure transport operator.
    """
    snapshot_file = SNAPSHOT_DIR / f"{case.name}.npz"
    if not snapshot_file.exists():
        pytest.skip(
            f"snapshot {snapshot_file.name} not yet generated; run "
            "`python -m tests.sn.regression._generate_snapshots`"
        )

    snap = np.load(snapshot_file)
    case_kind = str(snap["case_kind"]) if "case_kind" in snap.files else "eigen"
    expected_flux = np.asarray(snap["scalar_flux"], dtype=np.float64)

    cfg = case.builder()
    rc = run_config(cfg)
    result = run_case(cfg)

    if case_kind == "eigen":
        expected_keff = float(snap["keff"])
        assert_regression(
            result.keff, expected_keff,
            conv_tol=rc[EIGEN_KEFF_TOL_KEY],
            case_name=case.name, kind="iterative", quantity="k_eff",
        )
        flux_conv_tol = rc[EIGEN_FLUX_TOL_KEY]
    else:
        flux_conv_tol = rc[FIXED_SOURCE_FLUX_TOL_KEY]

    assert_regression(
        np.asarray(result.scalar_flux.values, dtype=np.float64),
        expected_flux,
        conv_tol=flux_conv_tol,
        case_name=case.name, kind="iterative", quantity="scalar_flux",
    )
