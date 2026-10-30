r"""Membership is the space's verb, a parse with two legs, both on every run.

Plan ``.claude/plans/space_membership_covector_gauge.md``, units U1 and U1b:

* the IDENTITY leg: a typed member (:class:`~orpheus.numerics.space.Member`; every :class:`~orpheus.numerics.field.Field` is one) handed to a space verb is compared with the space by content
  identity (``x.space == V``). This is the only check that separates two spaces of equal shape and different metric.
* the SHAPE leg: a bare array must be an element of the coordinate space
  (:mod:`orpheus.numerics.coordinate_space`): EXACT shape for the scalar verbs (``inner_product``, ``norm``), LEADING
  shape for the actions (``apply_metric``, ``apply_inverse_metric``, the Riesz legs), where a trailing batch rides
  through as ``A (x) I_P``.

Until U1b (2026-10-30) the shape leg ran under ``__debug__`` only, and a companion row pinned that under ``-O`` the raw
path was unchecked; that row pinned the degradation as the contract and was deleted when the check became the parse.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.axis import Axis, BasisKind, Declared, EnergyAxis
from orpheus.numerics.space import FunctionSpace
from dataclasses import dataclass

from orpheus.numerics.field import Field
from orpheus.numerics.units import SCALAR_FLUX_UNITS

pytestmark = pytest.mark.l0


_E = EnergyAxis.synthetic(2)
#: Two spaces of EQUAL shape and DIFFERENT metric: no shape check can tell them apart (explorer map §1, the structural fact).
V = FunctionSpace.of_axes(_E, Axis("spatial", (3,), weights=np.array([1.0, 2.0, 3.0]), origin=Declared(BasisKind.NODAL)))
W = FunctionSpace.of_axes(_E, Axis("spatial", (3,), weights=np.array([3.0, 2.0, 1.0]), origin=Declared(BasisKind.NODAL)))
X = np.arange(1.0, 7.0).reshape(2, 3)


@dataclass(frozen=True, eq=False, kw_only=True)
class _Member(Field):
    """A minimal concrete member: a typed element of whatever space it is given, built in numerics."""

    UNITS = SCALAR_FLUX_UNITS


# ── rung 1: the formal tier (runs under -O) ────────────────────────────────────────────────────────────

class TestATypedMemberIsComparedByContent:
    def test_a_member_is_admitted_and_reads_the_same_number_as_its_values(self):
        """RAW-vs-FORMAL equivalence (R1's pair): the typed call and the raw call are one arithmetic; array_equal."""
        f = _Member(values=X, space=V)
        assert V.inner_product(f, f) == V.inner_product(X, X) == 210.0   # [M] 210.0 today on the raw path
        assert V.norm(f) == V.norm(X)
        np.testing.assert_array_equal(np.asarray(V.apply_metric(f)), np.asarray(V.apply_metric(X)))

    def test_a_member_of_an_equal_shape_space_with_another_metric_is_REFUSED(self):
        """Before membership (`[M]` e42fd498) this raised a TypeError with no membership reason, and the raw spelling
        reads 210 on V and 154 on W, undetectably: only the typed tier tells the two spaces apart."""
        fW = _Member(values=X, space=W)
        with pytest.raises(ValueError, match="not a member of"):
            V.inner_product(fW, fW)
        with pytest.raises(ValueError, match="not a member of"):
            V.apply_metric(fW)

    def test_coordinates_of_returns_the_coordinates_parse_dont_validate(self):
        f = _Member(values=X, space=V)
        assert V.coordinates_of(f) is f.values


# ── rung 2: the shape leg of a bare array — scalar verbs exact, actions leading ────────────────────────────

_COLUMN_SPACE = FunctionSpace("g", (2, 1), inner_product_weights=np.array([[1.0], [2.0]]))
_COL = np.array([[0.6], [0.4]])
_FLAT = np.array([0.6, 0.4])


class TestTheRawScalarVerbsTakeExactMembers:
    def test_a_flat_carrier_on_a_column_space_is_REFUSED(self):
        """Before membership (`[M]` e42fd498) a two-sided broadcast admitted it (qa measured 1.4 for 0.68 with a
        1-D [1, 2] weight)."""
        with pytest.raises(ValueError, match="exact shape"):
            _COLUMN_SPACE.inner_product(_COL, _FLAT)

    def test_a_payload_tail_is_REFUSED_by_a_scalar_verb(self):
        """A batch has no inner product (P1 refined): V x P is spelled, never inferred."""
        with pytest.raises(ValueError, match="exact shape"):
            V.inner_product(np.ones(V.shape + (5,)), np.ones(V.shape + (5,)))


class TestTheRawActionsTakeALeadingMember:
    def test_a_trailing_batch_rides_through_an_action(self):
        """Admitted on every run (A (x) I_P): the legitimate payload callers that are actions (the Riesz legs, ``.H``)
        keep their trailing batch; the two scalar payload callers were re-posed on ``V * P``."""
        y = np.random.default_rng(0).standard_normal(V.shape + (4, 2))
        np.testing.assert_array_equal(V.apply_metric(y)[..., 1, 0], V.apply_metric(y[..., 1, 0]))

    def test_a_lower_rank_or_misaligned_element_is_REFUSED_by_an_action(self):
        """Leading-shape law: ``x.shape[:V.ndim] == V.shape``. Before membership (`[M]` e42fd498), (2, 1) on V was
        admitted by broadcast."""
        with pytest.raises(ValueError, match="leading shape"):
            V.apply_metric(np.ones((2, 1)))


# ── rung 3: the explicit tensor product replaces the inferred payload (re-pose witnesses) ───────────────

def test_the_payload_pairing_is_the_tensor_product_space_pairing():
    """The re-pose premise of the two direct payload callers [M] probes/u1_repose.py: on V (x) P with a counting P the
    pairing equals today's payload reading bit-for-bit (8 of 8 heads, and the derived-object file's synthetic axes)."""
    P = FunctionSpace.of_axes(Axis("payload", (5,), origin=Declared(BasisKind.NODAL)))
    x = np.random.default_rng(1).standard_normal(V.shape + (5,))
    expected = float(np.sum(V.apply_metric(x) * x))  # the action applied, then the counting contraction
    assert (V * P).inner_product(x, x) == expected


def test_two_axes_less_spaces_of_one_name_and_different_weights_are_two_spaces():
    """The identity of an axes-less space carries its metric (the U1 qa review, F2): before, two spaces of one name and
    shape with different weights compared equal, so a member of one was admitted by the other."""
    a = FunctionSpace("tier", (3,), inner_product_weights=np.array([1.0, 2.0, 3.0]))
    b = FunctionSpace("tier", (3,), inner_product_weights=np.array([3.0, 2.0, 1.0]))
    assert a != b
    assert a == FunctionSpace("tier", (3,), inner_product_weights=np.array([1.0, 2.0, 3.0]))
    with pytest.raises(ValueError, match="not a member of"):
        a.inner_product(_Member(values=np.ones(3), space=b), _Member(values=np.ones(3), space=b))
