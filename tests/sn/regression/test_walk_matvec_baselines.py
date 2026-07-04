r"""Frozen pre-carve 1-D walk matvec baselines — the 2.5a relocation anchor.

Re-runs each :mod:`._generate_walk_baselines` case and asserts BOTH
orientations of the ``(L+C)`` matvec reproduce the frozen snapshot at
``kind="direct"``, ``reduction_depth=1`` (a single walk pass has no
iteration; ``nulp=1`` is the tightest principled gate, with the
:class:`~tests.sn.regression._regression_assert.DriftWarning` tripwire
reporting ANY bit movement — escalate with ``-W error::DriftWarning`` for
the strict bit-identity contract the 2.5a relocation claims).

Why this is not the self-referential canary: the snapshot was captured
BEFORE the relocation, so it cannot move with a carve that relocates both
the SUT path and its in-process reference together (gate-spec §10.2; the
removal-form ``array_equal`` gates keep the override-not-leak role).

Lifetime: the curvilinear rows freeze the LAGGED Carlson-seed formulation
and RE-CAPTURE at 2.5d (ruling R10 — principled re-baseline); the slab row
stays bitwise through the phase. Disposition decided explicitly at 2.5e.
"""
from __future__ import annotations

import numpy as np
import pytest

from ._generate_walk_baselines import CASES, WalkMatvecCase, run_case
from ._generate_snapshots import SNAPSHOT_DIR
from ._regression_assert import assert_regression


pytestmark = [pytest.mark.regression, pytest.mark.foundation]


@pytest.mark.parametrize("case", CASES, ids=lambda c: c.name)
def test_walk_matvec_baseline(case: WalkMatvecCase) -> None:
    snapshot_file = SNAPSHOT_DIR / case.snapshot_name
    if not snapshot_file.exists():
        pytest.skip(
            f"snapshot {snapshot_file.name} not yet generated; run "
            "`python -m tests.sn.regression._generate_walk_baselines`"
        )

    snap = np.load(snapshot_file)
    blocks = run_case(case)
    for quantity, actual in blocks.items():
        assert_regression(
            actual,
            np.asarray(snap[quantity], dtype=np.float64),
            conv_tol=0.0,  # ignored for kind="direct"
            case_name=f"walk_matvec_{case.name}",
            kind="direct",
            reduction_depth=1,
            quantity=quantity,
        )
