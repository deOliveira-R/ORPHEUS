"""DD regression test against frozen reference snapshots.

Each snapshot at ``snapshots/<name>.npz`` carries (k_eff, scalar_flux)
captured by ``_generate_snapshots.py`` at the time of commit. This
test re-runs each case and asserts bit-for-bit agreement.

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

from orpheus.sn.solver import solve_sn

from ._generate_snapshots import CASES, SNAPSHOT_DIR, SnapshotCase


pytestmark = pytest.mark.regression


@pytest.mark.parametrize("case", CASES, ids=lambda c: c.name)
def test_dd_regression(case: SnapshotCase) -> None:
    """Re-run case and assert bit-identical match to frozen snapshot.

    Tolerance is set just above the iterative-solver convergence
    floor (``rtol=1e-12``, ``atol=1e-13``) so the test detects any
    semantic drift while tolerating reproducible BiCGSTAB / power-
    iteration noise at the last few digits.
    """
    snapshot_file = SNAPSHOT_DIR / f"{case.name}.npz"
    if not snapshot_file.exists():
        pytest.skip(
            f"snapshot {snapshot_file.name} not yet generated; run "
            "`python -m tests.sn.regression._generate_snapshots`"
        )

    snap = np.load(snapshot_file)
    expected_keff = float(snap["keff"])
    expected_flux = np.asarray(snap["scalar_flux"], dtype=np.float64)

    cfg = case.builder()
    result = solve_sn(
        cfg["materials"], cfg["mesh"], cfg["quadrature"],
        scattering_order=cfg.get("scattering_order", 0),
        max_outer=500, keff_tol=1e-12, flux_tol=1e-10,
        max_inner=300, inner_tol=1e-10,
    )

    np.testing.assert_allclose(
        result.keff, expected_keff, rtol=1e-12, atol=1e-13,
        err_msg=f"k_eff regression failed for {case.name!r}",
    )
    np.testing.assert_allclose(
        np.asarray(result.scalar_flux, dtype=np.float64),
        expected_flux, rtol=1e-12, atol=1e-13,
        err_msg=f"scalar_flux regression failed for {case.name!r}",
    )
