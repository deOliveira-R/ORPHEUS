r"""G04p: the three byte verbs on the dual of every axis-built harness space keep today's bytes (N04, spec 7.5, decision ``[ON04-1]``).

The record ``dual_probes/dual_probes.npz`` (147 arrays: ``apply_metric`` and ``apply_inverse_metric`` of the harness's seeded vector and ``metric_diagonal`` where it exists, on ``.dual()`` of
the 51 axis-built spaces the harness reaches, 18 of the arrays on the 6 diagonal-head moment spaces) was taken on the UNCHANGED base ``cf924aaa`` (the tip plus N03b) by the capture script of
the unit's opening over the tracked population builder ``_dual_probe_population.py``, before any production edit; it is re-taken by no one before N10. This gate rebuilds the population on the
tree under test and compares bytes (dtype, shape and ``.tobytes()``), key set first. The record sits beside ``snapshots/`` and ``axis_forms/`` (the precedent of ON03-12), so the harness' key-set
and directory-budget assertions and the 336 manifest are untouched.

Host-keyed: the arrays pass through LAPACK, so the gate compares stored bytes on the stamped host and FAILS on a foreign stamp with both stamps in the message (a skip would be a silent green);
it is not in the CI step. GREEN at the base tip by construction (the record is the base's). The input that reddens it is in the unit's own tree: ``[M]`` 18 of 147 arrays differ under the
N03 ``Factor.form`` alone, and dropping the density-end ``DualMetric`` branch of the head factor's form reddens 18 of 147 (with the branch 0 of 147). Mutation: that branch dropped.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tests.gates.cross_method.w1 import _dual_probe_population as pop
from tests.gates.cross_method.w1 import _snapshot as snap

pytestmark = pytest.mark.foundation

RECORD_DIR = Path(__file__).resolve().parent / "dual_probes"


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def test_the_dual_of_every_axis_built_space_has_the_bases_bytes() -> None:
    stamp = tuple((RECORD_DIR / "STAMP.txt").read_text().splitlines()[0].split("\t"))
    live_stamp = snap.host_stamp()
    _require(stamp == live_stamp, f"the dual record cannot be read on this host: it was taken on {stamp} and this host is {live_stamp}. Re-take it on this host at the last green commit before comparing.")
    record = np.load(RECORD_DIR / "dual_probes.npz")
    live = {f"{f}::{n}": a for (f, n), a in pop.capture_dual_probes().items()}
    _require(len(record.files) == 147, f"the record holds {len(record.files)} arrays, expected 147")
    _require(sorted(live) == sorted(record.files), f"the population's keys differ from the record's: {sorted(set(live) ^ set(record.files))[:6]}")
    bad = [
        k
        for k in sorted(live)
        if np.ascontiguousarray(live[k]).dtype != record[k].dtype or live[k].shape != record[k].shape or np.ascontiguousarray(live[k]).tobytes() != record[k].tobytes()
    ]
    _require(not bad, f"{len(bad)} of {len(live)} dual arrays differ from the base's bytes: {bad[:6]}")
