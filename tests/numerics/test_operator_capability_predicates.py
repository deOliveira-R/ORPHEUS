"""Keystone v2 — the permanent two-axis inverse/adjoint contract (#226 carve P4).

The frozenset-coexistence scaffold that used to live here (predicates ≡
``capabilities`` membership) DELETED with the frozenset at W2, its
licensing job done. What remains is the STANDING net (spec §36): the
recursive-composition pins (the closure laws asserted directly), the
keystone-v2 enumeration (predicate ⟺ method behaviour, with the inverse
axis's structural/value split pinned as CONTRACT), the eager-``.H``
raise pin, and the TypeGuard-bridge consistency leg.

The asymmetry fixtures — ``ZeroOperator`` (adjointable-but-NOT-invertible,
STRUCTURAL arm), the zero-entry ``DiagonalOperator`` (VALUE arm), and
``_ApplyOnly`` (neither) — are load-bearing: without operators that
distinguish the two axes, a buggy ``is_adjointable`` that merely returned
``is_invertible`` would pass on every both-True / both-False leaf and
ship silently (test-architect §2.3/§36).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.operator import (
    _AdjointOperator,
    DiagonalOperator,
    IdentityOperator,
    LinearOperator,
    MissingAdjoint,
    PermutationOperator,
    ScaledOperator,
    TensorProductOperator,
    TraceRestrictionOperator,
    ZeroOperator,
)
from orpheus.numerics.space import FunctionSpace
from tests._harness.predicates import (
    INVERTIBLE,
    STRUCTURAL_ABSENT,
    VALUE_RAISE,
    assert_inverse_adjoint_contract as _assert_contract,
)

pytestmark = pytest.mark.foundation


class _ApplyOnly(LinearOperator):
    """Synthetic apply-only operator — the (¬invertible, ¬adjointable)
    corner. Inherits the base ``False`` predicates. The non-adjointable
    summand that makes a sum half-adjointable."""
    # S4-amendment: the base DEMANDS an answer from every subclass; this
    # double is a deliberately-unbound probe, so it DECLARES the unbound
    # state instead of inheriting a silent default (which no longer exists).
    domain = None
    codomain = None

    def apply(self, x, /):
        return x


class _BoundDiagonalLeaf(LinearOperator):
    """A BOUND non-pointwise leaf for the _AdjointOperator rows: diagonal
    action with declared spaces, adjointable, invertible via the generic
    InverseOperator (whose adjoint axis is #280-deferred — the branch the
    wrapper rows exist to pin). Since the S4-amendment a pointwise member
    never builds the wrapper (its .H is itself) and an UNBOUND inner is
    REFUSED at the wrapper's constructor, so this bound stand-in is the
    honest fixture."""

    def __init__(self, c: np.ndarray) -> None:
        self._c = np.asarray(c, dtype=float)
        self._space = FunctionSpace("bound_leaf", self._c.shape)

    @property
    def domain(self):
        return self._space

    @property
    def codomain(self):
        return self._space

    def apply(self, x, /):
        return self._c * x

    def apply_transpose(self, x, /):
        return self._c * x

    @property
    def is_adjointable(self) -> bool:
        return True

    @property
    def is_invertible(self) -> bool:
        return bool(np.min(np.abs(self._c)) > 0.0)

    def inverse(self):
        from orpheus.numerics.operator import InverseOperator

        return InverseOperator(self)

    def solve(self, b, /):
        return b / self._c


_C = np.array([1.0, 2.0, 3.0])
_CZ = np.array([1.0, 0.0, 3.0])  # a zero entry → singular
#: One bound 3-space for the non-pointwise leaf fixtures (S4-amendment:
#: their .H demands declared metrics; a bare space = Euclidean).
_SP3 = FunctionSpace("p3", (3,))

# NOTE (B3.3): a ``_LEAVES`` inventory used to sit here, spanning the same
# (invertible × adjointable) quadrants. It was deleted with
# ``IncomingOrdinateMaskTensor``, and the deletion cost nothing: `[M]` it had
# ZERO consumers — nothing in this module or any other read it — while
# ``_CONTRACT_ROWS`` below covers the same leaves AND is actually
# parametrized over. It was a Pattern-2 duplicate that happened to be the
# unexercised copy.


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
    # inverse (step 4): is_invertible True.
    assert s_both.is_adjointable is True and s_both.is_invertible is True
    assert s_both.is_adjointable == (
        ident.is_adjointable and DiagonalOperator(_C).is_adjointable
    )
    assert s_both.is_invertible == ident.is_invertible  # the leading-term rule

    # The SAME summands reversed: leading _ApplyOnly is not invertible —
    # the ordering contract refuses to designate a preconditioner.
    s_reversed = _ApplyOnly() + ident
    assert s_reversed.is_invertible is False

    s_half = ident + _ApplyOnly()  # one summand non-adjointable
    assert s_half.is_adjointable is False
    assert s_half.is_invertible is True  # adjoint axis ≠ inverse axis


def test_product_predicates_recursive_and_faithful():
    """``(AB)`` invertible/adjointable iff BOTH factors."""
    p_both = DiagonalOperator(_C) @ PermutationOperator(np.array([1, 0, 2]))
    assert p_both.is_invertible is True and p_both.is_adjointable is True

    p_singular = DiagonalOperator(_CZ) @ IdentityOperator()  # one factor singular
    assert p_singular.is_invertible is False


def test_scaled_and_adjoint_predicates_faithful():
    """Scaling (α≠0) passes the predicates through. Two adjoint regimes
    since the S4-amendment's pointwise family: a POINTWISE member's ``.H``
    is itself (metric-free self-adjointness — the wrapper is never built,
    so the #280-deferred wrapper branches dissolve for the family), while
    the generic ``_AdjointOperator`` wrapper — pinned here by direct
    construction — keeps its CONDITIONAL swap law (#280 2.5c): False for
    this inner because ``DiagonalOperator``'s inverse
    (``InverseOperator``) is non-adjointable (#280-deferred), so
    ``(A⁻¹).H`` does not exist; the ``A.H.H`` direction stays deferred."""
    from orpheus.numerics.operator import _AdjointOperator

    sc = ScaledOperator(2.0, DiagonalOperator(_C))
    assert sc.is_invertible is True and sc.is_adjointable is True

    d = DiagonalOperator(_C)
    assert d.H is d  # the pointwise family law: no wrapper, full predicates
    assert d.H.is_invertible is True and d.H.is_adjointable is True

    adj = _AdjointOperator(_BoundDiagonalLeaf(_C))  # the generic wrapper
    assert adj.is_invertible is False and adj.is_adjointable is False


def test_inverse_objects_are_faithful():
    """I3 (#226 §12.3): the NEW inverse OBJECTS satisfy the contract too.

    The two-axis truth must hold on what ``.inverse()`` RETURNS, not just
    on the advertisers: every family member is ``is_invertible=True``
    (its ``solve`` is the wrapped forward) with the involution closing by
    object identity, and stays honest on the adjoint axis (False until
    #280). The full contract legs run in the keystone-v2 rows below;
    here the involution identities are pinned.
    """
    D = DiagonalOperator(_C)
    assert D.inverse().inverse() is D  # (A⁻¹)⁻¹ — the leaf, by identity
    # G-I3 (#226 step 4): the GreenOperator a Green-invertible SUM returns
    # is faithful too — {APPLY, SOLVE} via the mixin back-half,
    # is_invertible=True, adjoint axis honestly deferred (#280).
    s = IdentityOperator() + DiagonalOperator(_C)
    green = s.inverse()
    assert green.inverse() is s  # involution by object identity (mixin)
    # M-I3 (#226 step 5): the 4th sibling — DIRECT construction (the
    # ``.inverse()`` factory routing is #138/CP; no dispatch pin here).
    from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator

    minv = MatrixInverseOperator(D, basis_shape=(3,))
    assert minv.inverse() is D  # the inner, by object identity (mixin)


# ───────────────────────────────────────────────────────────────────────
# KEYSTONE v2 (carve P4, spec §36) — the PERMANENT two-axis contract:
# predicate ⟺ method behaviour, with the inverse axis's structural/value
# split pinned as CONTRACT (Design C), the adjoint axis's EAGER ``.H``
# raise, and the TypeGuard-bridge consistency. References NO frozenset.
# ───────────────────────────────────────────────────────────────────────

def _mixin_family():
    D = DiagonalOperator(_C)
    from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator

    # (name, op, is_adjointable) — the adjoint axis lands PER SIBLING
    # (#280 family): FREE for the materialized matrix (step 5 — the
    # ``trans=1`` backsolve on the SAME LU factors, the arm the swap law
    # rides), still deferred for the pointwise and splitting siblings.
    return [
        ("InverseOperator", D.inverse(), False),
        ("GreenOperator", (IdentityOperator() + D).inverse(), False),
        (
            # Over the BOUND leaf (S4-amendment): its adjointable arm takes
            # .H, which demands the inner's declared spaces.
            "MatrixInverseOperator",
            MatrixInverseOperator(_BoundDiagonalLeaf(_C), basis_shape=(3,)),
            True,
        ),
    ]


_CONTRACT_ROWS = [
    # (id, op, is_invertible, is_adjointable, inverse_contract)
    ("Identity", IdentityOperator(), True, True, INVERTIBLE),
    ("Zero", ZeroOperator(), False, True, STRUCTURAL_ABSENT),
    ("Diagonal", DiagonalOperator(_C), True, True, INVERTIBLE),
    ("Diagonal-singular", DiagonalOperator(_CZ), False, True, VALUE_RAISE),
    ("Permutation", PermutationOperator(np.array([1, 0, 2]), domain=_SP3, codomain=_SP3), True, True, INVERTIBLE),
    # The rank-deficient-but-adjointable representative. Was
    # ``IncomingOrdinateMaskTensor`` until B3.3 retired it; the successor
    # occupies the same quadrant for the same structural reason (a
    # restriction discards rows, so it is never invertible, and its
    # transpose is the scatter, so it is always adjointable). The row still
    # MOVES: the two are independent implementations, and the successor
    # declares both predicates on its own.
    ("TraceRestriction", TraceRestrictionOperator(np.array([0]), 3, domain=_SP3, codomain=FunctionSpace("p1", (1,))), False, True, STRUCTURAL_ABSENT),
    ("ApplyOnly", _ApplyOnly(), False, False, STRUCTURAL_ABSENT),
    # Composites — every arm of the recursive laws:
    ("Sum-leading-invertible", IdentityOperator() + DiagonalOperator(_C), True, True, INVERTIBLE),
    ("Sum-leading-not", _ApplyOnly() + IdentityOperator(), False, False, VALUE_RAISE),
    ("Sum-half-adjointable", IdentityOperator() + _ApplyOnly(), True, False, INVERTIBLE),
    ("Product-both", DiagonalOperator(_C) @ PermutationOperator(np.array([1, 0, 2]), domain=_SP3, codomain=_SP3), True, True, INVERTIBLE),
    ("Product-singular-factor", DiagonalOperator(_CZ) @ IdentityOperator(), False, True, VALUE_RAISE),
    ("Scaled", ScaledOperator(2.0, DiagonalOperator(_C)), True, True, INVERTIBLE),
    ("Scaled-singular", ScaledOperator(2.0, DiagonalOperator(_CZ)), False, True, VALUE_RAISE),
    ("TensorProduct-both", TensorProductOperator((DiagonalOperator(_C), IdentityOperator())), True, True, INVERTIBLE),
    ("TensorProduct-singular", TensorProductOperator((DiagonalOperator(_CZ), IdentityOperator())), False, True, VALUE_RAISE),
    # The adjoint wrapper: the swap law (A.H)⁻¹ = (A⁻¹).H is now DECLARED
    # (#280 2.5c), but THIS instance is not invertible — DiagonalOperator's
    # inverse (InverseOperator) is #280-deferred on the adjoint axis, so
    # is_invertible=False and inverse() RAISES NotInvertible (VALUE_RAISE, not
    # STRUCTURAL_ABSENT — the method now exists). A.H.H still raises
    # MissingAdjoint EAGERLY (adjoint-of-adjoint deferred). The inner is
    # the BOUND leaf since the S4-amendment: Diagonal is a
    # PointwiseOperator (its .H is itself, no wrapper), and the wrapper's
    # constructor now REFUSES an unbound inner — the wrapper's own
    # predicate derivation is what this row pins.
    ("AdjointWrapper", _AdjointOperator(_BoundDiagonalLeaf(_C)), False, False, VALUE_RAISE),
    # The wrap-delegate inverse family (SweepOperator = SN side, pinned
    # in tests/sn/operators/test_capability_survival.py); the adjoint
    # column is per-sibling — see _mixin_family:
    *[
        (name, op, True, adj, INVERTIBLE)
        for name, op, adj in _mixin_family()
    ],
]


@pytest.mark.parametrize(
    ("op", "inv", "adj", "contract"),
    [row[1:] for row in _CONTRACT_ROWS],
    ids=[row[0] for row in _CONTRACT_ROWS],
)
def test_inverse_adjoint_contract_keystone_v2(op, inv, adj, contract):
    """Spec §36 legs (a)/(b)/(c): the permanent two-axis faithfulness net.

    The three mandatory config-blindness fixtures (§36/§0.6) are in the
    rows: Zero (adjointable-not-invertible, STRUCTURAL), Diagonal with a
    TRUE zero entry (adjointable-not-invertible, VALUE), and the
    half-adjointable sum (invertible-not-adjointable) — without them a
    buggy predicate that mirrors the OTHER axis passes silently.
    """
    _assert_contract(op, invertible=inv, adjointable=adj, inverse_contract=contract)


def test_half_adjointable_sum_H_raises_eagerly():
    """Spec §38 (the M-ADJ-EAGER target): the raise site is ``.H`` itself.

    Pre-carve, ``A.H`` on a non-adjointable ``A`` SUCCEEDED and the
    failure was lazy at the first ``.apply`` — the broken-stub pattern.
    Design C raises ``MissingAdjoint`` at wrapper CONSTRUCTION.
    """
    s = DiagonalOperator(_C) + _ApplyOnly()
    assert s.is_adjointable is False
    with pytest.raises(MissingAdjoint, match="adjoint"):
        s.H
