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

from orpheus.geometry.boundary import ReflectiveBoundary, VacuumInflow
from orpheus.geometry.boundary._bound_compat import _BoundBoundaryOperator
from orpheus.numerics.elementary_operators import DiagonalOperator, PermutationOperator
from orpheus.numerics.operator import (
    IdentityOperator,
    InverseOperator,
    NotInvertible,
    OperatorProduct,
    OperatorSum,
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
            _ALPHA * diag,
            lambda q: q / (_ALPHA * _C7),
            _X7,
            4,
        ),
        (
            # NON-COMMUTING product: D@P applies P first, then D — so the
            # inverse MUST reverse (P⁻¹ then... i.e. P⁻¹(D⁻¹ q)); a
            # commuting pair would hide the M-INV-ORDER swap (§12.6).
            "diag_at_perm",
            diag @ perm,
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
    """``inverse().apply ≡ solve`` BIT-IDENTICALLY where solve EXISTS (§12.0).

    The delegation design's defining fact: one realization, no reciprocal
    twin — exact, not a tolerance claim. REWIRED at carve P4 (Design B):
    the algebra-closed advertisers (Identity/Permutation/Scaled/TP)
    retired their ``solve`` — their inverse IS a first-class forward,
    whose VALUE the closed-form anchor in ``test_i1_roundtrip`` pins —
    so on those rows this gate pins the RETIREMENT instead.
    """
    if hasattr(A, "solve"):
        np.testing.assert_array_equal(
            A.inverse().apply(x), A.solve(x), err_msg=f"{name}: inverse ≢ solve"
        )
    else:
        # S2: the scaled operator is the one-term ``OperatorSum`` (ruling 51).
        assert type(A) in (
            IdentityOperator, PermutationOperator, OperatorSum,
            TensorProductOperator,
        ), f"{name}: unexpected solve-less advertiser"



def test_i2_scaled():
    """(αA)⁻¹ = (1/α)·A⁻¹ as an action equality (§12.2)."""
    D = DiagonalOperator(_C7)
    lhs = (_ALPHA * D).inverse().apply(_X7)
    rhs = (1.0 / _ALPHA) * D.inverse().apply(_X7)
    np.testing.assert_array_almost_equal_nulp(lhs, rhs, nulp=2)


def test_i2_product_reversed():
    """(AB)⁻¹ = B⁻¹A⁻¹ on a NON-commuting pair (§12.2 — order-swap teeth)."""
    D, P = DiagonalOperator(_C7), PermutationOperator(_P7)
    lhs = (D @ P).inverse().apply(_X7)
    rhs = P.inverse().apply(D.inverse().apply(_X7))
    np.testing.assert_array_almost_equal_nulp(lhs, rhs, nulp=2)


def test_i2_involution_per_type():
    """(A⁻¹)⁻¹ = A, in each type's OWN strongest form (§12.0/§12.2)."""
    D = DiagonalOperator(_C7)
    assert D.inverse().inverse() is D  # InverseOperator: OBJECT identity
    P = PermutationOperator(_P7)
    np.testing.assert_array_equal(P.inverse().inverse().perm, P.perm)  # exact
    sc = _ALPHA * D
    np.testing.assert_array_almost_equal_nulp(  # Scaled: action (1/(1/α))
        sc.inverse().inverse().apply(_X7), sc.apply(_X7), nulp=4
    )
    # S2: the scaled inverse is the one-term combination (1/α)·D⁻¹, whose
    # inverse is (1/(1/α))·D = α·D exactly for α = −2.5 (1/(1/α) rounds back),
    # equal as combinations: the inverse arrow is keyed by its inner (Q5).
    assert sc.inverse().inverse() == sc
    I = IdentityOperator()
    assert I.inverse() is I


def test_singular_diagonal_inverse_raises():
    """NEGATIVE control — TRUE zero (1e-300 would pass the != 0 gate) §12.4."""
    with pytest.raises(NotInvertible, match="zero entry"):
        DiagonalOperator(np.array([2.0, 0.0, 0.5])).inverse()


def test_product_with_singular_factor_inverse_raises():
    with pytest.raises(NotInvertible, match="invertible"):
        (DiagonalOperator(np.array([1.0, 0.0, 2.0])) @ IdentityOperator()).inverse()


def test_inverse_operator_ctor_guards_singular_inner():
    with pytest.raises(NotInvertible, match="invertible"):
        InverseOperator(DiagonalOperator(np.array([1.0, 0.0])))


def test_a_sum_has_no_inverse_its_split_has_a_green():
    """The step-4 factory arm, RE-POSED at S2 (ruling 50): a sum of two terms
    has no inverse of its own (the left-spine head that designated the
    preconditioner is gone), and the Green of its NAMED split ``A = M − N``
    inverts it — the value anchored on the dense LU, not on the Green's
    own forward."""
    from orpheus.numerics.green_operator import GreenOperator

    D, N = DiagonalOperator(_C7), 0.25 * PermutationOperator(_P7)
    s = D - N
    with pytest.raises(NotInvertible, match="GreenOperator\\(M, N\\)"):
        s.inverse()
    green = GreenOperator(D, N)
    assert green.inverse() == s  # the involution: the forward IS M − N
    dense = np.diag(_C7) - 0.25 * np.eye(7)[_P7]
    ref = np.linalg.solve(dense, _X7)
    # The Green's promise is a NORM residual, so the band is normalised by
    # the solution's scale, not componentwise (a small entry carries the
    # whole vector's absolute error).
    np.testing.assert_allclose(
        green.apply(_X7), ref, rtol=0.0, atol=10 * green.tol * float(np.max(np.abs(ref))),
    )


def test_matrix_inverse_operator_universal_roundtrip():
    """§30.10 (step 5): the 4th sibling joins the universal net as a
    DIRECT-construction registry row — the ``.inverse()``-factory routing
    (``type(A.inverse()) is MatrixInverseOperator``) lands with #138/CP;
    a dispatch pin now would be a phantom."""
    from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator

    A = DiagonalOperator(_C7)
    minv = MatrixInverseOperator(A, basis_shape=(7,))
    np.testing.assert_allclose(minv.apply(A.apply(_X7)), _X7, atol=1e-12)
    np.testing.assert_allclose(minv.apply(_X7), _X7 / _C7, atol=1e-12)  # closed form
    assert minv.inverse() is A  # involution, object identity (mixin)


def test_product_inverse_accepts_seeded_apply():
    """§31.2 PROD-285-REPRO — the #285 witness: a composed inverse now
    carries the canonical seeded-apply kwarg (TypeError pre-#285: the raw
    reversed product's ``apply(x, /)`` was positional-only; accepted post:
    ``InverseOperator``'s mixin ``apply(x, /, *, initial_guess=None)``).
    The equality is the designed-green accept-and-ignore proof — the seed
    changes nothing, exactly like the exact leaves."""
    D, P = DiagonalOperator(_C7), PermutationOperator(_P7)
    inv = (D @ P).inverse()
    q = _RNG.standard_normal(7)
    out = inv.apply(q, initial_guess=q)  # <-- TypeError pre-#285
    np.testing.assert_array_equal(out, inv.apply(q))


def test_product_inverse_value_and_involution():
    """§31.3: the family wrapper is BIT-identical to the composition path
    the composed inverse action :math:`P^{-1}(D^{-1}q)` (the value the
    raw reversed product realized), and the involution STRENGTHENS to
    object identity (the raw product rebuilt fresh objects).
    NON-commuting factors (§0.6) so the M-PROD-FACTORORDER a/b swap
    reddens. Reference spelled through the factors' canonical inverses
    (carve P4: Permutation.solve retired; ``P.inverse().apply`` is the
    SAME integer gather, bit-identical)."""
    D, P = DiagonalOperator(_C7), PermutationOperator(_P7)
    prod = D @ P
    inv = prod.inverse()
    assert type(inv) is InverseOperator  # the generic member, not a raw product
    q = _RNG.standard_normal(7)
    np.testing.assert_array_equal(
        inv.apply(q), P.inverse().apply(D.inverse().apply(q)),
        err_msg="#285 product-inverse value ≠ P⁻¹(D⁻¹ q)",
    )
    assert inv.inverse() is prod  # strengthened involution


def test_algebra_closed_inverses_unchanged():
    """§31.4: the OTHER kind of inverse — algebra-closed structures whose
    inverse is a first-class FORWARD in the same closed family (a perm's
    inverse IS a perm). #285 does NOT route these through InverseOperator
    (documented in ``OperatorProduct.inverse``; not gated for seed)."""
    assert type(PermutationOperator(_P7).inverse()) is PermutationOperator
    assert type(IdentityOperator().inverse()) is IdentityOperator
    # S2: (2D)⁻¹ is the one-term combination ½·D⁻¹ (ruling 51), equal to the
    # same combination spelled directly (the inverse arrow keyed by its inner).
    scaled = 2.0 * DiagonalOperator(_C7)
    assert type(scaled.inverse()) is OperatorSum
    assert scaled.inverse() == 0.5 * scaled.terms[0][1].inverse()


def test_bound_shim_forwards_inverse():
    """The NINTH advertiser: the shim forwards ``.inverse()`` (§12.5)."""
    inner = PermutationOperator(_P7)
    shim = _BoundBoundaryOperator(inner, ReflectiveBoundary())
    np.testing.assert_array_equal(
        shim.inverse().apply(_X7), inner.inverse().apply(_X7)
    )
    # carve P4: the shim's solve forward retired with the surface (the
    # promised deletion) — solving through the shim IS .inverse().apply.
    assert not hasattr(_BoundBoundaryOperator, "solve")
    masked = _BoundBoundaryOperator(  # non-invertible inner → propagate raise
        DiagonalOperator(np.array([1.0, 0.0])), VacuumInflow()
    )
    with pytest.raises(NotInvertible):
        masked.inverse()
