r"""G04q: ``canonical_reals`` is the entrywise canonical form of an array of any rank, rank 0 included (N04, ruling W5b).

The law: for a ``bool``, integer or real array of rank 0 to 3 the result is the ``float64`` copy of its entries with ``-0.0`` made ``+0.0``; a NaN entry is a ``ValueError`` naming its index
(``()`` at rank 0), and the content digest of an array is a function of its canonical entries (equal arrays of rank 0 digest equal, and ``-0.0`` digests as ``+0.0``).

The rank-0 leg has a production consumer, measured: ``AngularTraceSpace._current_spaces`` declares an axis over the face's index shape after the ordinate axis is integrated out, which is ``()`` on a face layout
with no group axis, and its counting form is a rank-0 weight array that the space's name digest reads through ``canonical_reals`` (13 tests of ``tests/gates/sn/operators/test_sn_boundary_realizer.py`` build one). The widening
that admits it is ``np.argwhere`` for ``np.nonzero`` in the NaN scan, because ``np.nonzero`` refuses a rank-0 array (numpy 2.4).

Mutation: the NaN scan as ``zip(*np.nonzero(np.isnan(doubles)))`` reddens the rank-0 rows (a ``ValueError`` from numpy, not the keyed refusal, and the digest rows); ranks 1 to 3 stay green. Exact on any host: values and a raise.
"""

from __future__ import annotations

import re

import numpy as np
import pytest

from orpheus.values.content import content_digest
from orpheus.values.scalars import canonical_reals

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def hand_canonical(values: np.ndarray) -> np.ndarray:
    """The law spelled by hand: the float64 copy, with every zero made positive."""
    out = np.array(values, dtype=np.float64)
    out[out == 0.0] = 0.0
    return out


SAMPLES = {
    "rank 0": np.array(-0.0),
    "rank 0 value": np.array(2.5),
    "rank 1": np.array([-0.0, 1.5, -2.0]),
    "rank 2": np.array([[0.0, -0.0], [3.0, -4.5]]),
    "rank 3": np.arange(-4.0, 4.0).reshape(2, 2, 2),
    "rank 1 integers": np.array([0, 7, -3]),
    "rank 0 integer": np.array(5),
}


@pytest.mark.parametrize("name", sorted(SAMPLES))
def test_the_canonical_form_is_the_float_copy_with_positive_zeros(name: str) -> None:
    values = SAMPLES[name]
    got = canonical_reals(values, "sample")
    expected = hand_canonical(values)
    _require(got.shape == values.shape and got.dtype == np.float64, f"[{name}] shape {got.shape} or dtype {got.dtype} is not the input's shape and float64")
    _require(got.tobytes() == expected.tobytes(), f"[{name}] the canonical bytes differ from the hand form (a negative zero kept, or a value changed)")


@pytest.mark.parametrize("shape", [(), (3,), (2, 2), (2, 1, 2)])
def test_a_nan_entry_is_refused_naming_its_index(shape: tuple[int, ...]) -> None:
    values = np.zeros(shape)
    values[(0,) * len(shape)] = np.nan
    with pytest.raises(ValueError, match=re.escape(f"the entry {(0,) * len(shape)!r} is NaN, which is not a number")) as caught:
        canonical_reals(values, "sample")
    _require(type(caught.value) is ValueError and "nonzero" not in str(caught.value), f"[rank {len(shape)}] the refusal is not the keyed one: {caught.value}")


def test_the_digest_of_a_rank_zero_array_is_stable_and_reads_the_canonical_entry() -> None:
    first, again = content_digest(np.array(1.25)), content_digest(np.array(1.25))
    _require(first == again, "the digest of one rank-0 array differs between two calls")
    _require(content_digest(np.array(-0.0)) == content_digest(np.array(0.0)), "-0.0 and +0.0 digest differently at rank 0")
    _require(content_digest(np.array(1.25)) != content_digest(np.array(1.5)), "two different rank-0 values digest equal")
    _require(content_digest(np.array(1.25)) != content_digest(np.array([1.25])), "a rank-0 array and a one-entry rank-1 array digest equal")
