"""Faithfulness gate for the structural capability predicates (#226).

The inverse-as-operator carve replaces the stringly-typed
``capabilities: frozenset[str]`` advertisement with the per-axis, typed
:attr:`~orpheus.numerics.operator.LinearOperator.is_invertible` /
:attr:`~orpheus.numerics.operator.LinearOperator.is_adjointable` properties
(the runtime, instance-accurate successors to ``CAP_SOLVE`` /
``CAP_APPLY_TRANSPOSE``).

THIS gate is the keystone that licenses deleting the frozenset in Phase 4:
during coexistence it proves the derived predicates mirror the frozenset
EXACTLY for every operator. Its **permanent** successors (once the frozenset
is gone) are the recursive-composition pins below, which never reference
``capabilities`` — they assert the structural law directly
(``(a+b).is_adjointable == a.is_adjointable and b.is_adjointable``).

The two **asymmetry fixtures** — ``ZeroOperator`` (adjointable-but-NOT-
invertible) and ``_ApplyOnly`` (neither) — are load-bearing: without an
operator that distinguishes the two axes, a buggy ``is_adjointable`` that
merely returned ``is_invertible`` would pass on every both-True / both-False
leaf and ship silently (test-architect §2.3).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.operator import (
    CAP_APPLY,
    DiagonalOperator,
    IdentityOperator,
    IncomingOrdinateMaskTensor,
    LinearOperator,
    OperatorProduct,
    PeriodicWrapOperator,
    PermutationOperator,
    ScaledOperator,
    ZeroOperator,
)
from tests._harness.predicates import assert_capability_faithful as _assert_faithful

pytestmark = pytest.mark.foundation


class _ApplyOnly(LinearOperator):
    """Synthetic apply-only operator — the (¬invertible, ¬adjointable)
    corner. Inherits the base ``False`` predicates; advertises only
    ``CAP_APPLY``. The non-adjointable summand that makes a sum
    half-adjointable."""

    capabilities = frozenset({CAP_APPLY})

    def apply(self, x, /):
        return x


_C = np.array([1.0, 2.0, 3.0])
_CZ = np.array([1.0, 0.0, 3.0])  # a zero entry → singular

# Representative leaves spanning the (invertible × adjointable) quadrants.
_LEAVES = [
    IdentityOperator(),                            # both
    ZeroOperator(),                                # adjointable, NOT invertible
    DiagonalOperator(_C),                          # both
    DiagonalOperator(_CZ),                         # adjointable, NOT invertible (min|f|=0)
    PermutationOperator(np.array([1, 0, 2])),      # both
    IncomingOrdinateMaskTensor(np.array([0]), 3),  # adjointable, NOT invertible (rank-deficient)
    PeriodicWrapOperator(),                        # adjointable, NOT invertible
    _ApplyOnly(),                                  # neither
]


@pytest.mark.parametrize("op", _LEAVES, ids=lambda o: type(o).__name__)
def test_leaf_predicates_are_faithful_to_the_frozenset(op):
    """Every leaf's derived predicates match its ``capabilities`` exactly."""
    _assert_faithful(op)


def test_asymmetry_fixtures_break_the_axis_coincidence():
    """The two axes must be INDEPENDENTLY correct, not coincidentally equal."""
    z = ZeroOperator()
    assert z.is_adjointable is True and z.is_invertible is False
    dz = DiagonalOperator(_CZ)
    assert dz.is_adjointable is True and dz.is_invertible is False
    a = _ApplyOnly()
    assert a.is_adjointable is False and a.is_invertible is False


def test_sum_predicates_recursive_and_faithful():
    """``(A+B)`` adjointable iff BOTH; invertible iff the LEADING term is.

    The invertibility recursion is deliberately ORDER-SENSITIVE (#226
    taxonomy §12 step 4): the left-spine head is the Green splitting's
    designated preconditioner (canonical ordering — invertible term
    first), so the SAME operand pair flips the predicate when reversed.
    """
    ident = IdentityOperator()
    s_both = ident + DiagonalOperator(_C)
    # Leading Identity is invertible → the sum produces a GreenOperator
    # inverse (step 4): is_invertible True, CAP_SOLVE rides along.
    assert s_both.is_adjointable is True and s_both.is_invertible is True
    assert s_both.is_adjointable == (
        ident.is_adjointable and DiagonalOperator(_C).is_adjointable
    )
    assert s_both.is_invertible == ident.is_invertible  # the leading-term rule
    _assert_faithful(s_both)

    # The SAME summands reversed: leading _ApplyOnly is not invertible —
    # the ordering contract refuses to designate a preconditioner.
    s_reversed = _ApplyOnly() + ident
    assert s_reversed.is_invertible is False
    _assert_faithful(s_reversed)

    s_half = ident + _ApplyOnly()  # one summand non-adjointable
    assert s_half.is_adjointable is False
    assert s_half.is_invertible is True  # adjoint axis ≠ inverse axis
    _assert_faithful(s_half)


def test_product_predicates_recursive_and_faithful():
    """``(AB)`` invertible/adjointable iff BOTH factors."""
    p_both = DiagonalOperator(_C) @ PermutationOperator(np.array([1, 0, 2]))
    assert p_both.is_invertible is True and p_both.is_adjointable is True
    _assert_faithful(p_both)

    p_singular = DiagonalOperator(_CZ) @ IdentityOperator()  # one factor singular
    assert p_singular.is_invertible is False
    _assert_faithful(p_singular)


def test_scaled_and_adjoint_predicates_faithful():
    """Scaling (α≠0) passes the predicates through; the adjoint wrapper
    advertises neither solve nor transpose (the ``A.H.inverse()`` /
    ``A.H.H`` directions are deferred to #280)."""
    sc = ScaledOperator(2.0, DiagonalOperator(_C))
    assert sc.is_invertible is True and sc.is_adjointable is True
    _assert_faithful(sc)

    adj = DiagonalOperator(_C).H  # _AdjointOperator
    assert adj.is_invertible is False and adj.is_adjointable is False
    _assert_faithful(adj)


def test_inverse_objects_are_faithful():
    """I3 (#226 §12.3): the NEW inverse OBJECTS satisfy the keystone too.

    ``is_invertible ≡ CAP_SOLVE`` / ``is_adjointable ≡ CAP_APPLY_TRANSPOSE``
    must hold on what ``.inverse()`` RETURNS, not just on the advertisers:
    an ``InverseOperator`` advertises ``{APPLY, SOLVE}`` with
    ``is_invertible=True`` (its solve is the leaf's forward) and stays
    honest on the adjoint axis (False until #280).
    """
    D = DiagonalOperator(_C)
    _assert_faithful(D.inverse())
    _assert_faithful(OperatorProduct(D, PermutationOperator(np.roll(np.arange(3), 1))).inverse())
    _assert_faithful(ScaledOperator(2.0, D).inverse())
    assert D.inverse().inverse() is D  # (A⁻¹)⁻¹ — the leaf, by identity
    # G-I3 (#226 step 4): the GreenOperator a Green-invertible SUM returns
    # is faithful too — {APPLY, SOLVE} via the mixin back-half,
    # is_invertible=True, adjoint axis honestly deferred (#280).
    s = IdentityOperator() + DiagonalOperator(_C)
    green = s.inverse()
    _assert_faithful(green)
    assert green.inverse() is s  # involution by object identity (mixin)
    # M-I3 (#226 step 5): the 4th sibling — DIRECT construction (the
    # ``.inverse()`` factory routing is #138/CP; no dispatch pin here).
    from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator

    minv = MatrixInverseOperator(D, basis_shape=(3,))
    _assert_faithful(minv)
    assert minv.inverse() is D  # the inner, by object identity (mixin)
