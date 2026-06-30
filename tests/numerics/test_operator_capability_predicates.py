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
    CAP_APPLY_TRANSPOSE,
    CAP_SOLVE,
    DiagonalOperator,
    IdentityOperator,
    IncomingOrdinateMaskTensor,
    LinearOperator,
    PeriodicWrapOperator,
    PermutationOperator,
    ScaledOperator,
    ZeroOperator,
)

pytestmark = pytest.mark.foundation


class _ApplyOnly(LinearOperator):
    """Synthetic apply-only operator — the (¬invertible, ¬adjointable)
    corner. Inherits the base ``False`` predicates; advertises only
    ``CAP_APPLY``. The non-adjointable summand that makes a sum
    half-adjointable."""

    capabilities = frozenset({CAP_APPLY})

    def apply(self, x, /):
        return x


def _assert_faithful(op: LinearOperator) -> None:
    """``is_invertible ≡ CAP_SOLVE∈caps`` AND ``is_adjointable ≡ CAP_APPLY_TRANSPOSE∈caps``."""
    assert op.is_invertible == (CAP_SOLVE in op.capabilities), (
        f"{type(op).__name__}: is_invertible={op.is_invertible} but "
        f"(CAP_SOLVE in caps)={CAP_SOLVE in op.capabilities}"
    )
    assert op.is_adjointable == (CAP_APPLY_TRANSPOSE in op.capabilities), (
        f"{type(op).__name__}: is_adjointable={op.is_adjointable} but "
        f"(CAP_APPLY_TRANSPOSE in caps)={CAP_APPLY_TRANSPOSE in op.capabilities}"
    )


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
    """``(A+B)`` adjointable iff BOTH; never invertible (no general sum inverse)."""
    ident = IdentityOperator()
    s_both = ident + DiagonalOperator(_C)
    assert s_both.is_adjointable is True and s_both.is_invertible is False
    assert s_both.is_adjointable == (
        ident.is_adjointable and DiagonalOperator(_C).is_adjointable
    )
    _assert_faithful(s_both)

    s_half = ident + _ApplyOnly()  # one summand non-adjointable
    assert s_half.is_adjointable is False
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
