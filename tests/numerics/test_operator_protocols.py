r"""O.1 (Issue #208 / Wave O) — block-role classification mechanism.

Foundation tests for the bulk/full/boundary operator classification
introduced in Wave O: the :class:`~orpheus.numerics.operator.BlockRole`
enum, the value-based ``isinstance``-marker metaclass
(:class:`~orpheus.numerics.operator.BulkOperator` /
:class:`~orpheus.numerics.operator.FullOperator` /
:class:`~orpheus.numerics.operator.BoundaryOperator`), and the partition +
``None``-default semantics.

These pin the MECHANISM on toy / generic operators (no SN fixtures). The
SN leaf tagging (C/S/F → BULK, L / InvertibleOperator → FULL, the
realized boundary laws → BOUNDARY) lives in
``tests/sn/operators/test_operator_block_role.py``.

All three markers ship as of Wave O step O.4a.1-γ: the BULK / FULL
markers landed in O.1; the BOUNDARY marker landed in O.4a.1-γ once the
boundary laws ``B`` (``SNBoundaryRealizer.realize`` outputs) became
first-class operators carrying :attr:`BlockRole.BOUNDARY` and the
deprecated ``geometry.boundary.BoundaryOperator`` alias was retired
(O.4a.1-β), freeing the name for the marker.
"""
from __future__ import annotations

import pytest

from orpheus.numerics.operator import (
    BlockRole,
    BoundaryOperator,
    BulkOperator,
    FullOperator,
    IdentityOperator,
    OperatorSum,
    ZeroOperator,
)

pytestmark = [pytest.mark.foundation]


class _Tagged:
    """Minimal object carrying only a ``block_role``.

    Exercises the value-based ``isinstance`` metaclass without any
    operator machinery — proving the classification is an instance
    property read off ``block_role``, not class inheritance.
    """

    def __init__(self, role: "BlockRole | None") -> None:
        self.block_role = role


class TestBlockRoleEnum:
    def test_enum_has_the_three_roles(self) -> None:
        assert {r.name for r in BlockRole} == {"BULK", "FULL", "BOUNDARY"}

    def test_roles_are_distinct(self) -> None:
        assert len({BlockRole.BULK, BlockRole.FULL, BlockRole.BOUNDARY}) == 3


class TestIsinstanceMarkers:
    def test_bulk_role_matches_bulk_marker_only(self) -> None:
        op = _Tagged(BlockRole.BULK)
        assert isinstance(op, BulkOperator)
        assert not isinstance(op, FullOperator)
        assert not isinstance(op, BoundaryOperator)

    def test_full_role_matches_full_marker_only(self) -> None:
        op = _Tagged(BlockRole.FULL)
        assert isinstance(op, FullOperator)
        assert not isinstance(op, BulkOperator)
        assert not isinstance(op, BoundaryOperator)

    def test_boundary_role_matches_boundary_marker_only(self) -> None:
        op = _Tagged(BlockRole.BOUNDARY)
        assert isinstance(op, BoundaryOperator)
        assert not isinstance(op, BulkOperator)
        assert not isinstance(op, FullOperator)

    def test_none_role_matches_no_marker(self) -> None:
        op = _Tagged(None)
        assert not isinstance(op, BulkOperator)
        assert not isinstance(op, FullOperator)
        assert not isinstance(op, BoundaryOperator)

    def test_object_without_block_role_matches_no_marker(self) -> None:
        # getattr-default: a bare object is never a role marker.
        assert not isinstance(object(), BulkOperator)
        assert not isinstance(object(), FullOperator)
        assert not isinstance(object(), BoundaryOperator)

    def test_isinstance_is_value_based_not_inheritance(self) -> None:
        # _Tagged inherits nothing from the markers; membership is decided
        # solely by the block_role VALUE (the _BlockRoleMeta __instancecheck__).
        op = _Tagged(BlockRole.BOUNDARY)
        assert BoundaryOperator not in type(op).__mro__
        assert isinstance(op, BoundaryOperator)

    def test_classification_is_a_partition(self) -> None:
        """No single role satisfies more than one marker (exclusivity, vv L11)."""
        markers = (BulkOperator, FullOperator, BoundaryOperator)
        for role in (BlockRole.BULK, BlockRole.FULL, BlockRole.BOUNDARY):
            op = _Tagged(role)
            matched = [m for m in markers if isinstance(op, m)]
            assert len(matched) == 1, f"{role} matched {[m.__name__ for m in matched]}"


class TestGenericOperatorsAreUnclassified:
    """O.1 classifies only the fixed-role leaves; generic composition
    operators stay unclassified (``block_role is None``). Their role is
    derived from operands at O.2 (dispatch), not at O.1."""

    def test_identity_operator_is_unclassified(self) -> None:
        assert IdentityOperator().block_role is None
        assert not isinstance(IdentityOperator(), BulkOperator)
        assert not isinstance(IdentityOperator(), FullOperator)
        assert not isinstance(IdentityOperator(), BoundaryOperator)

    def test_zero_operator_is_unclassified(self) -> None:
        assert ZeroOperator().block_role is None

    def test_operator_sum_is_unclassified_at_O1(self) -> None:
        s = IdentityOperator() + ZeroOperator()
        assert isinstance(s, OperatorSum)
        assert s.block_role is None
        assert not isinstance(s, BulkOperator)
        assert not isinstance(s, FullOperator)
        assert not isinstance(s, BoundaryOperator)
