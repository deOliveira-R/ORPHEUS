"""Pyright error ratchet (ORPHEUS#226) — counts must only go down.

A software invariant about the codebase artifact itself, not a
physics equation: ``foundation`` per the V&V taxonomy. ``slow``
because a whole-package pyright pass is ~a minute. Skips when the
host pyright is absent (npm tool, see README setup) — the ratchet
enforces on machines that CAN measure, and silence-by-skip is
visible in the pytest summary, not a false green.

Both directions are failures by design:

- a module count ABOVE baseline is a type-error regression — fix the
  new errors (or, for a deliberate trade, justify and re-baseline);
- a module count BELOW baseline means someone burned errors down but
  left the ratchet slack — tighten it so the improvement can't
  silently erode: ``python -m tests._harness.pyright_ratchet --update``.
"""

from __future__ import annotations

import pytest

from tests._harness.pyright_ratchet import (
    UPDATE_CMD,
    collect_module_counts,
    find_pyright,
    read_baseline,
)

_PYRIGHT = find_pyright()

pytestmark = [
    pytest.mark.foundation,
    pytest.mark.slow,
    pytest.mark.skipif(
        _PYRIGHT is None,
        reason="host pyright not installed (npm install -g pyright)",
    ),
]


def test_pyright_error_counts_never_increase():
    baseline, baseline_version = read_baseline()
    assert _PYRIGHT is not None
    live, live_version = collect_module_counts(_PYRIGHT)

    version_hint = (
        "" if live_version == baseline_version else
        f" [pyright version drift: baseline measured with "
        f"{baseline_version}, this run is {live_version} — counts can "
        f"move without code changes; re-baseline if the drift explains "
        f"the delta]"
    )

    # One keyed pass over the union: every module compared exactly
    # once, in both directions — including burned-to-zero (absent
    # from live) and brand-new (absent from baseline) modules.
    deltas = {
        m: (baseline.get(m, 0), live.get(m, 0))
        for m in baseline.keys() | live.keys()
    }
    regressions = {m: bn for m, bn in deltas.items() if bn[1] > bn[0]}
    improvements = {m: bn for m, bn in deltas.items() if bn[1] < bn[0]}

    if regressions:
        pytest.fail(
            "pyright error count INCREASED (module: baseline -> now): "
            + ", ".join(f"{m}: {b} -> {n}" for m, (b, n) in sorted(regressions.items()))
            + ". Fix the new type errors; only re-baseline for a "
              "deliberate, justified trade." + version_hint
        )
    if improvements:
        pytest.fail(
            "pyright error count DECREASED (module: baseline -> now): "
            + ", ".join(f"{m}: {b} -> {n}" for m, (b, n) in sorted(improvements.items()))
            + f". Nice — now tighten the ratchet so it sticks: {UPDATE_CMD}"
            + version_hint
        )
