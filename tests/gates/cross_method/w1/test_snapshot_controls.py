r"""G00b, G00c and G00e: the controls of the bit-identity instrument. Each compares a capture with a capture in one process.

They are platform-free: no stored byte is read except the manifest's array count and the size of the snapshot directory. They
prove that the instrument can read a difference (a one-ulp step, a reversed slot order) and that its population is alive.

Every assertion is here and none in ``_snapshot.py``: a helper's ``assert`` is removed by ``python -O``. The gate bodies
(``check_*``) take the thing a mutation replaces as a parameter.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from pathlib import Path

import numpy as np
import pytest

from orpheus.numerics import metric as metric_module
from orpheus.numerics.metric import FactoredMetric
from tests.gates.cross_method.w1 import _fixtures, _snapshot as snap

MEGABYTE = 1_000_000


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


@pytest.fixture(scope="module")
def baseline():
    return snap.capture()


# ─────────────────────────── G00b ───────────────────────────


def check_one_ulp_control(captured: Mapping, compare: Callable = snap.differing_keys) -> None:
    """Every array, with entry 0 moved by one ``nextafter`` step, compares different from the capture."""
    moved = {}
    for key, array in captured.items():
        copy = np.array(array, copy=True)
        _require(copy.size and np.isfinite(copy.flat[0]), f"{key}: entry 0 cannot carry the control (empty or not finite)")
        copy.flat[0] = np.nextafter(copy.flat[0], np.inf)
        moved[key] = copy
    unseen = sorted(set(captured) - set(compare(captured, moved)))
    _require(not unseen, f"{len(unseen)} of {len(captured)} arrays read equal after a one-ulp step: {unseen[:10]}")


def test_one_ulp_control(baseline):
    check_one_ulp_control(baseline)


# ─────────────────────────── G00c ───────────────────────────

_ANGULAR_SPACES = ("angular_bulk_space", "angular_trial_space", "moment_space(2)")


def _reversed_walk(self, x, *, inverse):
    """``FactoredMetric._walk`` with the entries applied in reverse order."""
    verb = "apply_inverse" if inverse else "apply"
    out = metric_module.leading_coordinates(np.asarray(x), self._block_shape, owner=f"FactoredMetric.{verb}")
    starts, start = [], 0
    for block_shape, _ in self.entries:
        starts.append(start)
        start += len(block_shape)
    for (block_shape, factor), s in reversed(list(zip(self.entries, starts))):
        if factor is not None:
            out = factor.apply_block(out, start=s, block_shape=block_shape, inverse=inverse)
    return out


def check_slot_order_control(baseline: Mapping, walk: Callable = _reversed_walk) -> None:
    """With ``FactoredMetric._walk`` replaced by ``walk``, re-capture; the (d) arrays of the two populations must all differ."""
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(FactoredMetric, "_walk", walk)
        patched = snap.capture()
    differing = set(snap.differing_keys(baseline, patched))
    first = [k for k in baseline if any(k[1].startswith(f"d.{s}.") for s in _ANGULAR_SPACES)]
    second = [k for k in baseline if k[0].startswith("sn_") and k[1].startswith("d.full_field_space.")]
    _require(len(first) == 72 and len(second) == 24, f"the control populations are {len(first)} and {len(second)}, not 72 and 24")
    unseen1 = sorted(set(first) - differing)
    unseen2 = sorted(set(second) - differing)
    _require(not unseen1, f"{len(unseen1)} of 72 angular and moment-space arrays do not see the slot order: {unseen1[:6]}")
    _require(not unseen2, f"{len(unseen2)} of 24 SN full-field arrays do not see the slot order: {unseen2[:6]}")


def test_slot_order_control(baseline):
    check_slot_order_control(baseline)


# ─────────────────────────── G00e ───────────────────────────


def check_population_is_active(captured: Mapping, manifest: Path = snap.SNAPSHOT_DIR / "MANIFEST.txt",
                               directory: Path = snap.SNAPSHOT_DIR) -> None:
    """Four assertions on one capture: distinct Grams, the manifest's count, no all-zero array, the directory budget."""
    grams = {snap.sha1(a) for (_, name), a in captured.items() if name.endswith(".input_gram")}
    _require(len(grams) >= 21, f"{len(grams)} distinct Grams over the frame families, at least 21 expected")
    _, lines = snap.read_manifest(manifest)
    _require(len(captured) == len(lines), f"{len(captured)} arrays captured against {len(lines)} array lines in the manifest")
    zero = sorted(k for k, a in captured.items() if not np.any(a))
    _require(not zero, f"{len(zero)} captured arrays are all zero: {zero[:6]}")
    size = snap.directory_bytes(directory)
    _require(size <= MEGABYTE, f"{directory} holds {size} bytes, over the budget of {MEGABYTE}")


def test_population_is_active(baseline):
    check_population_is_active(baseline)


def test_the_families_are_the_21_of_the_specification():
    """The Gram check above counts what the capture builds; this pins how many families it builds."""
    _require(len(_fixtures.frame_families()) == 21, "the harness holds 21 frame families")
