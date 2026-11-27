r"""G05u (N05, ruling X4 of the second post-verify pass): ``TensorSpace.axis_index_of_coordinates`` reads the factor's array position from ``positions_spanned_by``.

A synthetic tensor space with a factor of two array positions BEFORE the energy factor (and another AFTER it): the index of every flat coordinate along ``energy`` is
``np.indices(space.shape)[position]``. On every production fixture the factors before energy span one array position each, so no fixture could tell the preceding-factor
arithmetic from ``position = k``; this row can. First red: the mutation ``position = k`` (the slot index in place of the array position). Exact on any host: integers.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.gates.numerics.w1 import _tensor_probe as tp


def _wide_axis(shape: tuple[int, ...], label: str):
    from orpheus.numerics.axis import Axis, BasisKind, Declared
    from orpheus.numerics.metric import NoForm

    return Axis(label, shape, origin=Declared(BasisKind.NODAL, NoForm()))


def _wide_factor(shape: tuple[int, ...], label: str):
    from orpheus.numerics import factor as F

    return F.Factor(_wide_axis(shape, label), F.FormlessVariance.MEASURE, F.Subquotient.whole(int(np.prod(shape))))


def _space(before: tuple[int, ...], ng: int, after: tuple[int, ...]):
    factors = []
    if before:
        factors.append(_wide_factor(before, "before"))
    factors.append(tp.formless_factor(tp.formless_axis(ng, "energy")))
    if after:
        factors.append(_wide_factor(after, "after"))
    return tp.tensor_space().of_axes(*factors)


@pytest.mark.parametrize("before,ng,after", [((2, 3), 4, ()), ((2, 3), 4, (3, 2)), ((), 4, (2, 3)), ((3,), 2, (2, 2)), ((2, 2, 3), 3, (2,))])
def test_axis_index_reads_the_array_position_of_the_factor(before: tuple[int, ...], ng: int, after: tuple[int, ...]) -> None:
    space = _space(before, ng, after)
    position = len(before)
    expected = np.indices(space.shape)[position].ravel()
    got = space.axis_index_of_coordinates("energy")
    if got.shape != expected.shape or not np.array_equal(got, expected):
        pytest.fail(f"before={before} ng={ng} after={after}: axis_index_of_coordinates('energy') differs from np.indices(shape)[{position}]")
    if got.dtype != np.intp:
        pytest.fail(f"the index is {got.dtype}, expected intp")


def test_a_factor_of_two_positions_has_no_single_index() -> None:
    space = _space((2, 3), 4, ())
    with pytest.raises(ValueError, match="no single index"):
        space.axis_index_of_coordinates("before")
