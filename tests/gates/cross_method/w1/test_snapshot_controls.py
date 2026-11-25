r"""G00b, G00c and G00e: the controls of the bit-identity instrument. Each compares a capture with a capture in one process.

They are platform-free: no stored byte is read except the manifest's array count and the size of the snapshot directory. They
prove that the instrument can read a difference (a one-ulp step, a reversed slot order) and that its population is alive.

Every assertion is here and none in ``_snapshot.py``: a helper's ``assert`` is removed by ``python -O``. The gate bodies
(``check_*``) take the thing a mutation replaces as a parameter.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from pathlib import Path

import numpy as np
import pytest

from orpheus.numerics.coordinate_space import leading_coordinates
from orpheus.numerics.metric import FactoredMetric, NoForm
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


def _walk_in(order: Callable[[list], Iterable]) -> Callable:
    """A copy of ``FactoredMetric._walk`` (``orpheus/numerics/metric.py``) that applies the entries in the order ``order`` gives.

    ``order`` is the identity for the forward copy and ``reversed`` for the control. The two copies differ in that argument only,
    and ``check_slot_order_control`` asserts that the forward copy reproduces the production capture, so the copy cannot drift
    from production unseen.
    """
    def walk(self, x, *, inverse):
        verb = "apply_inverse" if inverse else "apply"
        out = leading_coordinates(np.asarray(x), self._block_shape, owner=f"FactoredMetric.{verb}")
        starts, start = [], 0
        for block_shape, _ in self.entries:
            starts.append(start)
            start += len(block_shape)
        for (block_shape, factor), s in order(list(zip(self.entries, starts))):
            if factor is not None and not isinstance(factor, NoForm):  # a NoForm entry applies nothing, as a None entry did
                out = factor.apply_block(out, start=s, block_shape=block_shape, inverse=inverse)
        return out
    return walk


def _capture_under(walk: Callable) -> dict:
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(FactoredMetric, "_walk", walk)
        return snap.capture()


def check_slot_order_control(baseline: Mapping, order: Callable[[list], Iterable] = reversed) -> None:
    """Re-capture under a copy of ``FactoredMetric._walk`` in the slot order ``order`` gives; three standing assertions.

    (1) The forward copy equals the production capture on every array (the copy is faithful). (2) Under ``order``, the (d)
    arrays of the two populations of 72 and 24 all differ. (3) Under ``order`` nothing else differs: the 72 order-blind (d)
    arrays and every other array read equal, so the copy is a pure slot permutation. Both population counts are reported.
    """
    forward = snap.differing_keys(baseline, _capture_under(_walk_in(lambda entries: entries)))
    _require(not forward, f"the forward copy of the walk differs from production on {len(forward)} of {len(baseline)} arrays: {forward[:6]}")
    differing = set(snap.differing_keys(baseline, _capture_under(_walk_in(order))))
    first = [k for k in baseline if any(k[1].startswith(f"d.{s}.") for s in _ANGULAR_SPACES)]
    second = [k for k in baseline if k[0].startswith("sn_") and k[1].startswith("d.full_field_space.")]
    _require(len(first) == 72 and len(second) == 24, f"the control populations are {len(first)} and {len(second)}, not 72 and 24")
    unseen1, unseen2 = sorted(set(first) - differing), sorted(set(second) - differing)
    stray = sorted(differing - set(first) - set(second))
    _require(
        not (unseen1 or unseen2),
        f"{len(unseen1)} of 72 angular and moment-space arrays and {len(unseen2)} of 24 SN full-field arrays do not see the slot order: "
        f"{unseen1[:4]} {unseen2[:4]}",
    )
    _require(not stray, f"{len(stray)} arrays outside the two populations differ under the reordered walk: {stray[:6]}")


def test_slot_order_control(baseline):
    check_slot_order_control(baseline)


# ─────────────────────────── G00e ───────────────────────────


def check_population_is_active(captured: Mapping, manifest: Path = snap.SNAPSHOT_DIR / "MANIFEST.txt",
                               directory: Path = snap.SNAPSHOT_DIR) -> None:
    """Four assertions on one capture: distinct Grams, the manifest's count and key set, no all-zero array, the directory budget (subdirectories included)."""
    grams = {snap.sha1(a) for (_, name), a in captured.items() if name.endswith(".input_gram")}
    _require(len(grams) >= 21, f"{len(grams)} distinct Grams over the frame families, at least 21 expected")
    _, lines = snap.read_manifest(manifest)
    _require(len(captured) == len(lines), f"{len(captured)} arrays captured against {len(lines)} array lines in the manifest")
    listed = snap.manifest_keys(lines)
    _require(set(captured) == listed, f"the captured keys differ from the manifest's: only captured {sorted(set(captured) - listed)[:4]}, "
                                      f"only in the manifest {sorted(listed - set(captured))[:4]}")
    zero = sorted(k for k, a in captured.items() if not np.any(a))
    _require(not zero, f"{len(zero)} captured arrays are all zero: {zero[:6]}")
    size = snap.directory_bytes(directory)
    _require(size <= MEGABYTE, f"{directory} holds {size} bytes, over the budget of {MEGABYTE}")


def test_population_is_active(baseline):
    check_population_is_active(baseline)


def test_the_families_are_the_21_of_the_specification():
    """The Gram check above counts what the capture builds; this pins how many families it builds."""
    _require(len(_fixtures.frame_families()) == 21, "the harness holds 21 frame families")
