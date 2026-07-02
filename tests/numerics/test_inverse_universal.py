r"""Universal inverse-operator gates — taxonomy §12 step 1 (#226).

The keystone threading the whole inverse family (``operator_machinery_
taxonomy.md`` §9/§12): ``inv.apply(A.apply(x)) ≈ x`` for EVERY
``(A, A.inverse())`` pair, paired with a **closed-form value anchor** per
type (round-trip alone is self-consistency, not correctness — both arms
share the coefficient; the closed form is the structurally-independent
"was right"). Spec: ``issue_226_inverse_operator_verification.md`` §12.

Tolerances are PER TYPE = the reduction depth of the FP path (§12.0):
zero-arithmetic gathers (Identity, Permutation) are exact ``array_equal``;
one multiply + one divide (Diagonal) is ``nulp(2)``; composites inherit the
sum. Config-blindness (§0.6/§12.6) drives the fixtures: non-unit non-uniform
``c`` (a unit ``c`` nulls M-INV-APPLY), a NON-commuting ``Diagonal @
Permutation`` (a commuting product cannot red the order-swap M-INV-ORDER),
a NON-involution 7-cycle permutation (an involution is blind to M-INV-PERM),
and ``α ≠ ±1`` (α=1 is blind to M-INV-SCALE).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.boundary._bound_compat import _BoundBoundaryOperator
from orpheus.numerics.operator import (
    DiagonalOperator,
    IdentityOperator,
    InverseOperator,
    MissingCapability,
    OperatorProduct,
    PermutationOperator,
    ScaledOperator,
    TensorProductOperator,
)

pytestmark = pytest.mark.foundation

_RNG = np.random.default_rng(226)

# Non-unit, non-uniform coefficient (activates M-INV-APPLY: uniform c=1
# would make A⁻¹(Ax) = c²x = x even under the apply-for-solve mutation).
_C7 = _RNG.uniform(0.5, 3.0, 7)
_C5 = _RNG.uniform(0.5, 3.0, 5)
# Non-involution permutations (a 7-cycle / 4-cycle: perm[perm] != identity),
# activating M-INV-PERM (an involution equals its own inverse).
_P7 = np.roll(np.arange(7), 1)
_P4 = np.roll(np.arange(4), 1)
_ALPHA = -2.5  # α ∉ {0, ±1} — activates M-INV-SCALE

_X7 = _RNG.standard_normal(7)
_X54 = _RNG.standard_normal((5, 4))


def _rows():
    """(name, A, closed_form_inverse, x, nulp) registry — §12.1."""
    diag = DiagonalOperator(_C7)
    perm = PermutationOperator(_P7)
    return [
        ("diagonal", diag, lambda q: q / _C7, _X7, 2),
        ("identity", IdentityOperator(), lambda q: q, _X7, 0),
        ("permutation", perm, lambda q: np.take(q, np.argsort(_P7)), _X7, 0),
        (
            "scaled_diagonal",
            ScaledOperator(_ALPHA, diag),
            lambda q: q / (_ALPHA * _C7),
            _X7,
            4,
        ),
        (
            # NON-COMMUTING product: D@P applies P first, then D — so the
            # inverse MUST reverse (P⁻¹ then... i.e. P⁻¹(D⁻¹ q)); a
            # commuting pair would hide the M-INV-ORDER swap (§12.6).
            "diag_at_perm",
            OperatorProduct(diag, perm),
            lambda q: np.take(q / _C7, np.argsort(_P7)),
            _X7,
            2,
        ),
        (
            "tensor_diag_perm",
            TensorProductOperator(
                (DiagonalOperator(_C5, axis=0), PermutationOperator(_P4, axis=1))
            ),
            lambda q: np.take(q / _C5[:, None], np.argsort(_P4), axis=1),
            _X54,
            2,
        ),
    ]


def _assert_close(actual, expected, nulp, msg):
    if nulp == 0:
        np.testing.assert_array_equal(actual, expected, err_msg=msg)
    else:
        np.testing.assert_array_almost_equal_nulp(actual, expected, nulp=nulp)


@pytest.mark.parametrize("name,A,closed,x,nulp", _rows(), ids=lambda r: str(r))
def test_i1_roundtrip_and_closed_form(name, A, closed, x, nulp):
    """I1 both ways + the structurally-independent closed-form anchor (§12.1)."""
    inv = A.inverse()
    _assert_close(inv.apply(A.apply(x)), x, nulp, f"{name}: A⁻¹(Ax) ≠ x")
    _assert_close(A.apply(inv.apply(x)), x, nulp, f"{name}: A(A⁻¹x) ≠ x")
    _assert_close(inv.apply(x), closed(x), nulp, f"{name}: A⁻¹ ≠ closed form")


@pytest.mark.parametrize("name,A,closed,x,nulp", _rows(), ids=lambda r: str(r))
def test_i1_inverse_apply_is_solve_bit_identical(name, A, closed, x, nulp):
    """``inverse().apply ≡ solve`` BIT-IDENTICALLY, all types (§12.0).

    The delegation design's defining fact: one realization, no reciprocal
    twin — exact for every advertiser, not a tolerance claim.
    """
    np.testing.assert_array_equal(
        A.inverse().apply(x), A.solve(x), err_msg=f"{name}: inverse ≢ solve"
    )


def test_i2_scaled():
    """(αA)⁻¹ = (1/α)·A⁻¹ as an action equality (§12.2)."""
    D = DiagonalOperator(_C7)
    lhs = ScaledOperator(_ALPHA, D).inverse().apply(_X7)
    rhs = (1.0 / _ALPHA) * D.inverse().apply(_X7)
    np.testing.assert_array_almost_equal_nulp(lhs, rhs, nulp=2)


def test_i2_product_reversed():
    """(AB)⁻¹ = B⁻¹A⁻¹ on a NON-commuting pair (§12.2 — order-swap teeth)."""
    D, P = DiagonalOperator(_C7), PermutationOperator(_P7)
    lhs = OperatorProduct(D, P).inverse().apply(_X7)
    rhs = P.inverse().apply(D.inverse().apply(_X7))
    np.testing.assert_array_almost_equal_nulp(lhs, rhs, nulp=2)


def test_i2_involution_per_type():
    """(A⁻¹)⁻¹ = A, in each type's OWN strongest form (§12.0/§12.2)."""
    D = DiagonalOperator(_C7)
    assert D.inverse().inverse() is D  # InverseOperator: OBJECT identity
    P = PermutationOperator(_P7)
    np.testing.assert_array_equal(P.inverse().inverse().perm, P.perm)  # exact
    sc = ScaledOperator(_ALPHA, D)
    np.testing.assert_array_almost_equal_nulp(  # Scaled: action (1/(1/α))
        sc.inverse().inverse().apply(_X7), sc.apply(_X7), nulp=4
    )
    I = IdentityOperator()
    assert I.inverse() is I


def test_singular_diagonal_inverse_raises():
    """NEGATIVE control — TRUE zero (1e-300 would pass the != 0 gate) §12.4."""
    with pytest.raises(MissingCapability, match="zero entry"):
        DiagonalOperator(np.array([2.0, 0.0, 0.5])).inverse()


def test_product_with_singular_factor_inverse_raises():
    with pytest.raises(MissingCapability, match="invertible"):
        OperatorProduct(
            DiagonalOperator(np.array([1.0, 0.0, 2.0])), IdentityOperator()
        ).inverse()


def test_inverse_operator_ctor_guards_singular_inner():
    with pytest.raises(MissingCapability, match="invertible"):
        InverseOperator(DiagonalOperator(np.array([1.0, 0.0])))


def test_bound_shim_forwards_inverse():
    """The NINTH advertiser: the shim forwards ``.inverse()`` (§12.5)."""
    inner = PermutationOperator(_P7)
    shim = _BoundBoundaryOperator(inner, kind="reflective")
    np.testing.assert_array_equal(
        shim.inverse().apply(_X7), inner.inverse().apply(_X7)
    )
    np.testing.assert_array_equal(  # every forwarded cap has a backing method
        shim.solve(_X7), inner.solve(_X7)
    )
    masked = _BoundBoundaryOperator(  # non-invertible inner → propagate raise
        DiagonalOperator(np.array([1.0, 0.0])), kind="vacuum"
    )
    with pytest.raises(MissingCapability):
        masked.inverse()
