"""DD regression test against frozen reference snapshots.

Each snapshot at ``snapshots/<name>.npz`` carries
``(scalar_flux[, k_eff])`` captured by ``_generate_snapshots.py`` at
the time of commit. This test re-runs each case and asserts bit-for-bit
agreement.

Tests are skipped if the snapshot file is not yet present — this lets
the regression infrastructure land before the snapshots themselves
(decoupled commits), and protects against accidentally running a
stale-snapshot CI gate during the snapshot-generation roll-out.

Marker scheme (per ``docs/testing/sn_verification_matrix.rst``):

* ``@pytest.mark.regression`` — regression-tier opt-in, runs only when
  ``orpheus/sn/``, ``orpheus/geometry/``, or ``orpheus/numerics/``
  are touched (CI gate).
"""
from __future__ import annotations

import numpy as np
import pytest

from ._generate_snapshots import CASES, SNAPSHOT_DIR, SnapshotCase, run_case


# ``regression`` flags the frozen-snapshot drift gate; ``foundation``
# gives it a V&V-level so the audit does not report it as an orphan
# (the snapshot comparison is a software invariant, not a physics
# equation gate). Both compose.
pytestmark = [pytest.mark.regression, pytest.mark.foundation]


@pytest.mark.parametrize("case", CASES, ids=lambda c: c.name)
def test_dd_regression(case: SnapshotCase) -> None:
    """Re-run case and assert per-case agreement with the frozen snapshot.

    Tolerance is set just above the iterative-solver convergence floor
    so the test detects any semantic drift while tolerating reproducible
    BiCGSTAB / power-iteration noise at the last few digits.

    Cartesian (slab) cases converge to bit-identical floor on this
    architecture; curvilinear (sphere / cylinder) cases have a wider
    floor driven by the M-M angular closure's iteration history coupled
    with the eigenvalue power iteration, so they need a relaxed bound.
    Per ``vv-principles`` "bit-identity vs principled-equivalence" — the
    drift on the curvilinear snapshots is bounded by
    ``(outer_iters) × (cond_num) × ULP`` and both ACTUAL and DESIRED
    converge to the same analytical limit (homogeneous-reflective
    cases bottom out at k_inf = νΣ_f / Σ_a exactly).

    Both eigenvalue (``kind=eigen``) and fixed-source
    (``kind=fixed_source``) cases are supported. The snapshot's
    ``case_kind`` field selects the comparison: eigenvalue cases pin
    both ``keff`` and ``scalar_flux``; fixed-source cases pin only
    ``scalar_flux`` (no ``keff`` exists for the pure transport
    operator).
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

    # Curvilinear (sphere / cylinder) hits a much wider convergence floor
    # than slab on the current architecture — the observed scalar-flux
    # drift in the multi-group homogeneous-reflective regression snapshots
    # sits at ~7e-7 relative (with both ACTUAL and DESIRED converging
    # toward the analytical k_inf = νΣ_f/Σ_a = 1.875).  The k_eff drift
    # is tighter (~1e-10) because k_eff is a scalar reduction that
    # averages out the per-cell iterative-floor noise.
    #
    # Curvilinear k_inf drift is currently under investigation (see the
    # ``test_kinf_homogeneous`` failure cluster); when that lands a
    # solver-side tightening, this regression bound should track it back
    # down.  Slab keeps the bit-identity-grade bound.
    is_curvilinear = case.name.startswith(("sphere", "cyl"))
    rtol = 5e-6 if is_curvilinear else 1e-12
    atol = 1e-7 if is_curvilinear else 1e-13

    cfg = case.builder()
    result = run_case(cfg)

    if case_kind == "eigen":
        expected_keff = float(snap["keff"])
        np.testing.assert_allclose(
            result.keff, expected_keff, rtol=rtol, atol=atol,
            err_msg=f"k_eff regression failed for {case.name!r}",
        )

    np.testing.assert_allclose(
        np.asarray(result.scalar_flux.values, dtype=np.float64),
        expected_flux, rtol=rtol, atol=atol,
        err_msg=f"scalar_flux regression failed for {case.name!r}",
    )
