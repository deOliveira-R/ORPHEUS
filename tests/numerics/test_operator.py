"""Foundation tests for ``orpheus.numerics.operator``.

These are software invariants — composition algebra closure, capability
propagation, scipy interop — not physics-equation verification. They
exercise the operator-algebra layer on numpy-matrix-backed fixtures so
that any regression in the protocol or its composers is caught
independently of any transport solver.

All tests use ``@pytest.mark.foundation`` per the V&V harness rules
(software invariants, not L0/L1/L2 claims about a solver).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.operator import (
    IdentityOperator,
    InverseOperator,
    LinearOperator,
    NotInvertible,
    OperatorProduct,
    OperatorSum,
    PermutationOperator,
    ScaledOperator,
    ZeroOperator,
)

pytestmark = pytest.mark.foundation


# ───────────────────────────────────────────────────────────────────────
# Fixtures: numpy-matrix-backed LinearOperator implementations
# ───────────────────────────────────────────────────────────────────────


class MatrixOperator(LinearOperator):
    """Test operator backed by a dense numpy matrix.

    Carries honest structural predicates per constructor flags. Used to
    drive composition + capability tests with a fully analytic ground
    truth (the matrix is stored, so any test can compare against
    direct ``M @ x``). Inherits :class:`LinearOperator` so that
    ``+`` / ``-`` / ``*`` / ``@`` invoke the composers.
    """

    def __init__(
        self,
        matrix: np.ndarray,
        *,
        can_solve: bool = False,
        can_transpose: bool = False,
    ) -> None:
        self.matrix = np.asarray(matrix, dtype=float)
        self._can_solve = bool(can_solve)
        self._can_transpose = bool(can_transpose)

    # Faithful test double: the predicates mirror the ctor flags exactly
    # (carve P4 §38 — the eager ``.H`` gate reads is_adjointable, so a
    # fixture advertising the cap without the predicate is a stale lie).
    @property
    def is_invertible(self) -> bool:
        return self._can_solve

    @property
    def is_adjointable(self) -> bool:
        return self._can_transpose

    def apply(self, x: np.ndarray) -> np.ndarray:
        return self.matrix @ x

    def solve(self, b: np.ndarray) -> np.ndarray:
        return np.linalg.solve(self.matrix, b)

    def inverse(self) -> InverseOperator:
        # is_invertible=True PROMISES a working inverse() (keystone v2) —
        # the solve-backed leaf serves it through the generic family
        # member, exactly like the production value leaves.
        if not self._can_solve:
            raise NotInvertible(
                "MatrixOperator(can_solve=False) produces no inverse."
            )
        return InverseOperator(self)

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        return self.matrix.T @ x


@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.default_rng(20260506)


@pytest.fixture
def n_dim() -> int:
    return 5


@pytest.fixture
def matrix_full(rng, n_dim):
    """A well-conditioned non-singular matrix advertising every capability."""
    M = rng.standard_normal((n_dim, n_dim))
    # Diagonal-shifted to be non-singular.
    M += n_dim * np.eye(n_dim)
    return MatrixOperator(M, can_solve=True, can_transpose=True)


@pytest.fixture
def matrix_apply_only(rng, n_dim):
    """An operator advertising only apply (e.g. a scattering source)."""
    M = rng.standard_normal((n_dim, n_dim))
    return MatrixOperator(M, can_solve=False, can_transpose=False)


@pytest.fixture
def vector(rng, n_dim):
    return rng.standard_normal(n_dim)


# ───────────────────────────────────────────────────────────────────────
# Protocol conformance
# ───────────────────────────────────────────────────────────────────────


def test_matrix_operator_is_a_linear_operator(matrix_full):
    """The runtime-checkable Protocol accepts conforming objects."""
    assert isinstance(matrix_full, LinearOperator)


def test_identity_and_zero_are_linear_operators():
    assert isinstance(IdentityOperator(), LinearOperator)
    assert isinstance(ZeroOperator(), LinearOperator)


# ───────────────────────────────────────────────────────────────────────
# Structural-predicate propagation — OperatorSum
# ───────────────────────────────────────────────────────────────────────


def test_sum_composes_with_both_applying(matrix_full, matrix_apply_only):
    s = OperatorSum(matrix_full, matrix_apply_only)  # eager guard accepts
    assert callable(getattr(s, "apply", None))


def test_sum_solve_rides_the_leading_term(matrix_full, matrix_apply_only):
    """Solve on a sum rides its LEADING term — order-sensitively.

    There is still no OPERAND-WISE algorithm for ``(A+B)^{-1}`` from
    ``A^{-1}``/``B^{-1}`` — but since #226 step 4 a sum whose left-spine
    head is invertible carries the preconditioned-splitting
    :class:`GreenOperator` inverse — invertibility rides the leading
    term, order-sensitively. REWIRED at carve P4: the old
    "solve never propagates" pin predated step 4 and survived only
    because the fixture's ``is_invertible`` was a stale ``False``.
    """
    from orpheus.numerics.green_operator import GreenOperator

    s = OperatorSum(matrix_full, matrix_apply_only)
    assert s.is_invertible is True
    assert type(s.inverse()) is GreenOperator
    # The SAME operands reversed: a non-invertible head has no
    # designated preconditioner — the ordering contract refuses.
    s_reversed = OperatorSum(matrix_apply_only, matrix_full)
    assert s_reversed.is_invertible is False


def test_sum_transpose_propagates_with_both(matrix_full):
    s = OperatorSum(matrix_full, matrix_full)
    assert s.is_adjointable is True


def test_sum_transpose_drops_when_one_lacks(matrix_full, matrix_apply_only):
    s = OperatorSum(matrix_full, matrix_apply_only)
    assert s.is_adjointable is False


# ───────────────────────────────────────────────────────────────────────
# Structural-predicate propagation — OperatorProduct
# ───────────────────────────────────────────────────────────────────────


def test_product_composes_with_both_applying(matrix_full, matrix_apply_only):
    p = OperatorProduct(matrix_full, matrix_apply_only)  # eager guard accepts
    assert callable(getattr(p, "apply", None))


def test_product_invertible_with_both(matrix_full):
    p = OperatorProduct(matrix_full, matrix_full)
    assert p.is_invertible is True


def test_product_not_invertible_when_one_factor_is_not(
    matrix_full, matrix_apply_only, vector
):
    p = OperatorProduct(matrix_full, matrix_apply_only)
    assert p.is_invertible is False
    with pytest.raises(NotInvertible, match="both factors"):
        p.solve(vector)  # the kept verb refuses eagerly (value arm)
    with pytest.raises(NotInvertible, match="both factors"):
        p.inverse()


def test_product_transpose_propagates_with_both(matrix_full):
    p = OperatorProduct(matrix_full, matrix_full)
    assert p.is_adjointable is True


# ───────────────────────────────────────────────────────────────────────
# Structural-predicate propagation — ScaledOperator
# ───────────────────────────────────────────────────────────────────────


def test_scaled_preserves_both_axes(matrix_full):
    s = ScaledOperator(2.5, matrix_full)
    # Scalar multiplication (α≠0) is unitary on both structural axes.
    assert s.is_invertible == matrix_full.is_invertible
    assert s.is_adjointable == matrix_full.is_adjointable


def test_scaled_apply_only_stays_apply_only(matrix_apply_only):
    s = ScaledOperator(0.5, matrix_apply_only)
    assert not s.is_invertible and not s.is_adjointable


def test_scaled_zero_scalar_rejected(matrix_full):
    with pytest.raises(ValueError, match="ZeroOperator"):
        ScaledOperator(0.0, matrix_full)


# ───────────────────────────────────────────────────────────────────────
# Structural predicates — Identity / Zero
# ───────────────────────────────────────────────────────────────────────


def test_identity_both_axes_and_self_inverse():
    i = IdentityOperator()
    assert i.is_invertible and i.is_adjointable
    assert i.inverse() is i           # algebra-closed: I^{-1} = I
    assert not hasattr(i, "solve")    # carve P4: solving IS .inverse().apply


def test_zero_is_structurally_non_invertible():
    """Critical: the zero map is the singular map par excellence.

    STRUCTURAL arm (Design C): no ``inverse()`` is declared at all —
    misuse is a static error — while the adjoint axis stays honest
    (0^T = 0).
    """
    z = ZeroOperator()
    assert z.is_invertible is False
    assert z.is_adjointable is True
    assert not hasattr(z, "inverse")


def test_apply_only_operator_is_not_solvable(matrix_apply_only):
    """GAP-3 (B4 carve safety net) — the advertised capability set is the
    single source of truth for invertibility, NOT the mere presence of a
    ``solve`` surface.

    The inverse-as-operator carve gave the operator family an
    ``inverse()`` surface; this pin guards that an apply-only leaf does
    NOT thereby become solvable: ``is_invertible`` stays ``False`` and
    ``inverse()`` refuses.

    GAP-3′ (carve P4, spec §37.4): the runtime truth is the predicate;
    this fixture DECLARES ``inverse()`` (a solve-backed test double), so
    it exercises the VALUE arm — ``inverse()`` refuses eagerly with
    ``NotInvertible``. The STRUCTURAL arm (no method at all) is pinned by
    ``ZeroOperator``/``_ApplyOnly`` in the keystone-v2 enumeration. The
    invertibility obligation's home, the inverse BUILDER, is pinned in
    ``test_iteration.py``.
    """
    assert not matrix_apply_only.is_invertible
    with pytest.raises(NotInvertible, match="no inverse|produces no inverse"):
        matrix_apply_only.inverse()


# ───────────────────────────────────────────────────────────────────────
# Negative tests: composition-time apply-guard (eager TypeError)
# ───────────────────────────────────────────────────────────────────────


class NoApplyOperator:
    """Pathological operand: genuinely has NO ``apply`` — the eager
    composition guard (carve P4: a callable check, no registry) must
    reject it at construction time, never mid-iteration."""


def test_sum_rejects_missing_apply_at_composition(matrix_full):
    bad = NoApplyOperator()
    with pytest.raises(TypeError, match="apply"):
        OperatorSum(matrix_full, bad)
    with pytest.raises(TypeError, match="apply"):
        OperatorSum(bad, matrix_full)


def test_product_rejects_missing_apply_at_composition(matrix_full):
    bad = NoApplyOperator()
    with pytest.raises(TypeError, match="apply"):
        OperatorProduct(matrix_full, bad)
    with pytest.raises(TypeError, match="apply"):
        OperatorProduct(bad, matrix_full)


def test_scaled_rejects_missing_apply_at_composition():
    bad = NoApplyOperator()
    with pytest.raises(TypeError, match="apply"):
        ScaledOperator(2.0, bad)


# ───────────────────────────────────────────────────────────────────────
# Numerical correctness — apply / solve / transpose
# ───────────────────────────────────────────────────────────────────────


def test_sum_apply_matches_dense(matrix_full, matrix_apply_only, vector):
    s = OperatorSum(matrix_full, matrix_apply_only)
    expected = matrix_full.matrix @ vector + matrix_apply_only.matrix @ vector
    np.testing.assert_allclose(s.apply(vector), expected, rtol=1e-12)


def test_product_apply_matches_dense(matrix_full, matrix_apply_only, vector):
    p = OperatorProduct(matrix_full, matrix_apply_only)
    expected = matrix_full.matrix @ (matrix_apply_only.matrix @ vector)
    np.testing.assert_allclose(p.apply(vector), expected, rtol=1e-12)


def test_product_solve_reverses_order(matrix_full, vector):
    """Verify (AB)^{-1} = B^{-1} A^{-1}."""
    p = OperatorProduct(matrix_full, matrix_full)
    composed = matrix_full.matrix @ matrix_full.matrix
    expected = np.linalg.solve(composed, vector)
    np.testing.assert_allclose(p.solve(vector), expected, rtol=1e-10)


def test_scaled_apply_matches_dense(matrix_full, vector):
    s = ScaledOperator(2.5, matrix_full)
    expected = 2.5 * (matrix_full.matrix @ vector)
    np.testing.assert_allclose(s.apply(vector), expected, rtol=1e-12)


def test_scaled_inverse_divides(matrix_full, vector):
    """(αL)^{-1} = (1/α) L^{-1} — rewired at carve P4: ScaledOperator.solve
    retired (the inverse is ALGEBRA-CLOSED), solving is .inverse().apply."""
    s = ScaledOperator(2.5, matrix_full)
    expected = (1.0 / 2.5) * np.linalg.solve(matrix_full.matrix, vector)
    np.testing.assert_allclose(s.inverse().apply(vector), expected, rtol=1e-10)
    assert not hasattr(ScaledOperator, "solve")  # the retirement pin


def test_identity_apply(vector):
    np.testing.assert_array_equal(IdentityOperator().apply(vector), vector)


def test_zero_apply(vector):
    np.testing.assert_array_equal(
        ZeroOperator().apply(vector), np.zeros_like(vector)
    )


def test_transpose_distributes_over_sum(matrix_full, vector):
    """(A + B)^T x = A^T x + B^T x."""
    s = OperatorSum(matrix_full, matrix_full)
    expected = (matrix_full.matrix.T @ vector) + (matrix_full.matrix.T @ vector)
    np.testing.assert_allclose(s.apply_transpose(vector), expected, rtol=1e-12)


def test_transpose_reverses_product(matrix_full, vector):
    """(A B)^T x = B^T A^T x."""
    p = OperatorProduct(matrix_full, matrix_full)
    expected = matrix_full.matrix.T @ (matrix_full.matrix.T @ vector)
    np.testing.assert_allclose(p.apply_transpose(vector), expected, rtol=1e-12)


# ───────────────────────────────────────────────────────────────────────
# Operator-algebra dunder smoke tests
# ───────────────────────────────────────────────────────────────────────


def test_dunder_add_returns_sum(matrix_full):
    s = matrix_full + matrix_full
    assert isinstance(s, OperatorSum)


def test_dunder_sub_returns_sum_with_negation(matrix_full, vector):
    """A - B == A + (-1) * B."""
    s = matrix_full - matrix_full
    expected = matrix_full.matrix @ vector - matrix_full.matrix @ vector
    np.testing.assert_allclose(s.apply(vector), expected, atol=1e-12)


def test_dunder_mul_scalar_returns_scaled(matrix_full):
    s = 3.0 * matrix_full
    assert isinstance(s, ScaledOperator)
    assert s.scalar == 3.0


def test_dunder_rmul_scalar_returns_scaled(matrix_full):
    s = matrix_full * 3.0
    assert isinstance(s, ScaledOperator)
    assert s.scalar == 3.0


def test_dunder_matmul_returns_product(matrix_full):
    p = matrix_full @ matrix_full
    assert isinstance(p, OperatorProduct)


def test_dunder_neg_returns_scaled_minus_one(matrix_full, vector):
    """``-A`` returns ``ScaledOperator(-1.0, A)``. Issue #185."""
    neg = -matrix_full
    assert isinstance(neg, ScaledOperator)
    assert neg.scalar == -1.0
    # Value check: ``(-A) x == -(A x)``.
    np.testing.assert_allclose(
        neg.apply(vector), -matrix_full.apply(vector), atol=1e-12,
    )


def test_dunder_neg_composes_with_subtraction(matrix_full, vector):
    """``A + (-B)`` is value-equivalent to ``A - B``. Issue #185."""
    a_plus_neg_b = matrix_full + (-matrix_full)
    a_minus_b = matrix_full - matrix_full
    np.testing.assert_allclose(
        a_plus_neg_b.apply(vector),
        a_minus_b.apply(vector),
        atol=1e-12,
    )


def test_dunder_truediv_returns_scaled_reciprocal(matrix_full, vector):
    """``A / α`` returns ``ScaledOperator(1/α, A)``. Issue #185."""
    halved = matrix_full / 2.0
    assert isinstance(halved, ScaledOperator)
    assert halved.scalar == 0.5
    # Value check: ``(A / α) x == (A x) / α``.
    np.testing.assert_allclose(
        halved.apply(vector), matrix_full.apply(vector) / 2.0, atol=1e-12,
    )


def test_dunder_truediv_by_zero_raises(matrix_full):
    """``A / 0`` raises ``ZeroDivisionError`` per Python convention."""
    with pytest.raises(ZeroDivisionError):
        matrix_full / 0.0


def test_dunder_truediv_rejects_non_numeric(matrix_full):
    """``A / 'string'`` returns ``NotImplemented`` → ``TypeError``."""
    with pytest.raises(TypeError):
        matrix_full / "two"  # type: ignore[operator]


def test_transport_eigenvalue_algebra_smoke(matrix_full, matrix_apply_only):
    """Smoke test the (L - S - F) ψ algebra without dispatching it.

    Build the natural transport operator expression with the algebra
    dunders and verify the capability set + numerical action match
    the explicit composition. This is the integration the SN/MoC/CP
    consumers will perform downstream.
    """
    L = matrix_full
    S = matrix_apply_only
    F = matrix_apply_only
    op = L - S - F
    # The invertibility RIDES the invertible leading loss term: the
    # transport normal form (L - S - F)^{-1} IS the preconditioned
    # splitting (#226 step 4; rewired at carve P4).
    assert op.is_invertible is True
    x = np.ones(L.matrix.shape[0])
    expected = L.matrix @ x - S.matrix @ x - F.matrix @ x
    np.testing.assert_allclose(op.apply(x), expected, rtol=1e-12)


# ═══════════════════════════════════════════════════════════════════════
# Issue 9.6 dunder + adjoint additions
# ═══════════════════════════════════════════════════════════════════════
#
# These tests cover the new affordances shipped in Phase B:
#   * ``A(x)`` aliases ``A.apply(x)``.
#   * ``A**n`` for non-negative ``n`` builds repeated composition.
#   * ``A.H`` aliases ``A.adjoint()``.
#   * Hilbert-adjoint identity ``<A x, y> = <x, A.H y>`` for both
#     Euclidean and quadrature-weight-aware inner products.
#   * ``IncompatibleOperatorComposition`` raised on
#     domain/codomain mismatch in the composers.
#   * ``__repr__`` smoke tests.

from orpheus.numerics.operator import IncompatibleOperatorComposition
from orpheus.numerics.space import FunctionSpace


def test_call_aliases_apply(matrix_full, vector):
    """``A(x)`` reads as math; should equal ``A.apply(x)``."""
    np.testing.assert_array_equal(
        matrix_full(vector), matrix_full.apply(vector),
    )


def test_pow_zero_returns_identity(matrix_full, vector):
    """``A**0`` is the identity operator."""
    I = matrix_full ** 0
    np.testing.assert_allclose(I.apply(vector), vector)


def test_pow_one_returns_self(matrix_full):
    """``A**1`` is ``self`` exactly."""
    assert matrix_full ** 1 is matrix_full


def test_pow_three_matches_repeated_apply(matrix_full, vector):
    """``(A**3)(x)`` equals ``A(A(A(x)))`` to round-off."""
    Acubed = matrix_full ** 3
    expected = matrix_full.apply(
        matrix_full.apply(matrix_full.apply(vector))
    )
    np.testing.assert_allclose(Acubed.apply(vector), expected, rtol=1e-12)


def test_pow_negative_rejected(matrix_full):
    """Negative powers raise — operator inverse is not auto-constructed."""
    with pytest.raises(ValueError, match="inverse"):
        _ = matrix_full ** -1


def test_pow_non_integer_returns_notimplemented(matrix_full):
    """Non-integer exponent is rejected by the dunder protocol."""
    with pytest.raises(TypeError):
        _ = matrix_full ** 1.5


def test_H_aliases_adjoint(matrix_full):
    """``A.H`` returns an adjoint wrapper just like ``A.adjoint()``."""
    H1 = matrix_full.H
    H2 = matrix_full.adjoint()
    # Both wrap the same inner.
    assert getattr(H1, "inner", None) is matrix_full
    assert getattr(H2, "inner", None) is matrix_full


def test_adjoint_euclidean_identity(matrix_full, rng, n_dim):
    """For Euclidean inner products, ``<A u, v> == <u, A.H v>``."""
    u = rng.standard_normal(n_dim)
    v = rng.standard_normal(n_dim)
    Au = matrix_full.apply(u)
    Hv = matrix_full.H.apply(v)
    lhs = float(np.dot(Au, v))
    rhs = float(np.dot(u, Hv))
    assert np.isclose(lhs, rhs, rtol=1e-12)


class _SpacedMatrixOperator(LinearOperator):
    """MatrixOperator carrying explicit domain / codomain function spaces.

    Used for the weight-aware Hilbert-adjoint identity tests, where the
    inner product on the codomain carries non-trivial weights and the
    adjoint must invert them. Mirrors ``MatrixOperator`` but exposes
    the function-space metadata shipped in Phase B.
    """

    def __init__(
        self,
        matrix: np.ndarray,
        domain_space: FunctionSpace,
        codomain_space: FunctionSpace,
    ) -> None:
        self.matrix = np.asarray(matrix, dtype=float)
        self._domain = domain_space
        self._codomain = codomain_space

    @property
    def is_adjointable(self) -> bool:
        # The eager .H gate reads the predicate (carve P4 §38).
        return True

    @property
    def domain(self) -> FunctionSpace:
        return self._domain

    @property
    def codomain(self) -> FunctionSpace:
        return self._codomain

    def apply(self, x: np.ndarray) -> np.ndarray:
        return self.matrix @ x

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        return self.matrix.T @ x


def test_adjoint_weighted_identity(rng):
    """Hilbert-adjoint identity holds with quadrature-weight metadata.

    Codomain inner product carries non-trivial diagonal weights ``w_W``;
    the adjoint must produce an output that, evaluated under the
    domain inner product, balances ``<A u, v>_W``.
    """
    n_dom = 4
    n_cod = 4
    M = rng.standard_normal((n_cod, n_dom))
    w_W = np.array([1.0, 2.0, 3.0, 4.0])  # codomain weights (non-trivial)
    domain = FunctionSpace(name="V", shape=(n_dom,))  # Euclidean
    codomain = FunctionSpace(
        name="W", shape=(n_cod,), inner_product_weights=w_W,
    )
    A = _SpacedMatrixOperator(M, domain, codomain)
    u = rng.standard_normal(n_dom)
    v = rng.standard_normal(n_cod)

    Au = A.apply(u)                     # in W
    Hv = A.H.apply(v)                   # in V
    lhs = codomain.inner_product(Au, v)
    rhs = domain.inner_product(u, Hv)
    assert np.isclose(lhs, rhs, rtol=1e-12)


def test_adjoint_swaps_domain_and_codomain():
    """``A.H.domain`` is ``A.codomain`` and vice versa."""
    n_dom = 3
    n_cod = 4
    M = np.random.default_rng(0).standard_normal((n_cod, n_dom))
    domain = FunctionSpace(name="V", shape=(n_dom,))
    codomain = FunctionSpace(name="W", shape=(n_cod,))
    A = _SpacedMatrixOperator(M, domain, codomain)
    assert A.H.domain == codomain
    assert A.H.codomain == domain


def test_sum_rejects_incompatible_domains():
    """OperatorSum with mismatched domains raises IncompatibleOperatorComposition."""
    n = 3
    M = np.eye(n)
    V1 = FunctionSpace(name="V1", shape=(n,))
    V2 = FunctionSpace(name="V2", shape=(n,))
    W = FunctionSpace(name="W", shape=(n,))
    A = _SpacedMatrixOperator(M, V1, W)
    B = _SpacedMatrixOperator(M, V2, W)
    with pytest.raises(IncompatibleOperatorComposition, match="domain"):
        _ = A + B


def test_sum_rejects_incompatible_codomains():
    n = 3
    M = np.eye(n)
    V = FunctionSpace(name="V", shape=(n,))
    W1 = FunctionSpace(name="W1", shape=(n,))
    W2 = FunctionSpace(name="W2", shape=(n,))
    A = _SpacedMatrixOperator(M, V, W1)
    B = _SpacedMatrixOperator(M, V, W2)
    with pytest.raises(IncompatibleOperatorComposition, match="codomain"):
        _ = A + B


def test_product_rejects_incompatible_inner_space():
    """``A @ B`` requires ``A.domain == B.codomain``."""
    n = 3
    M = np.eye(n)
    V = FunctionSpace(name="V", shape=(n,))
    W = FunctionSpace(name="W", shape=(n,))
    Z = FunctionSpace(name="Z", shape=(n,))
    A = _SpacedMatrixOperator(M, W, V)  # A: W → V
    B = _SpacedMatrixOperator(M, V, Z)  # B: V → Z; codomain Z != A.domain W
    with pytest.raises(
        IncompatibleOperatorComposition, match="A @ B|domain|codomain",
    ):
        _ = A @ B


def test_compatible_composition_succeeds_with_spaces():
    """Sum / product where the spaces match is accepted."""
    n = 3
    M = np.eye(n)
    V = FunctionSpace(name="V", shape=(n,))
    W = FunctionSpace(name="W", shape=(n,))
    A = _SpacedMatrixOperator(M, V, W)
    B = _SpacedMatrixOperator(M, V, W)
    s = A + B  # same domain V, same codomain W: OK
    assert s.domain == V
    assert s.codomain == W


def test_composition_skips_check_when_either_is_none(matrix_full):
    """Backward-compat: operators without function spaces pass freely.

    ``MatrixOperator`` (the existing fixture) returns ``None`` for
    domain and codomain. Sums and products against another ``None``
    operator must not raise — the composition compatibility check
    only fires when BOTH operators carry function spaces.
    """
    sum_ok = matrix_full + matrix_full
    assert sum_ok is not None
    prod_ok = matrix_full @ matrix_full
    assert prod_ok is not None


def test_repr_format(matrix_full, matrix_apply_only):
    """``repr(op)`` carries the class name, domain/codomain, and the
    two-axis predicate tokens (present iff True — carve P4 §44.F)."""
    r = repr(matrix_full)
    assert "MatrixOperator" in r
    # MatrixOperator has no function-space tagging; expect '?' placeholders.
    assert "domain='?'" in r
    assert "codomain='?'" in r
    assert "invertible" in r and "adjointable" in r
    r_bare = repr(matrix_apply_only)
    assert "invertible" not in r_bare and "adjointable" not in r_bare


def test_repr_with_function_spaces_shows_names():
    n = 3
    M = np.eye(n)
    V = FunctionSpace(name="phi", shape=(n,))
    W = FunctionSpace(name="psi", shape=(n,))
    A = _SpacedMatrixOperator(M, V, W)
    r = repr(A)
    assert "'phi'" in r
    assert "'psi'" in r


# ───────────────────────────────────────────────────────────────────────
# The OperatorProduct.solve RE-ROUTE (carve P4, spec §40) — sentinel +
# dense anchor + pre-carve baseline inheritance + the factor-kind matrix
# ───────────────────────────────────────────────────────────────────────

_REROUTE_BASELINE = "tests/numerics/data/step6_product_reroute_baseline.npz"


def _reroute_fixtures():
    """The EXACT pre-carve capture fixtures (seed 2860 — spec §40 leg iv).

    Mirrors the capture script that ran the OLD ``b.solve(a.solve(q))``
    spelling on the pre-carve tree; the snapshot is the INHERITANCE
    reference the re-routed spelling must reproduce.
    """
    from orpheus.numerics.operator import DiagonalOperator, TensorProductOperator  # noqa: F401

    rng = np.random.default_rng(2860)
    q = rng.standard_normal(7)
    D1 = DiagonalOperator(np.array([2.0, 0.5, 4.0, 1.25, 8.0, 0.25, 5.0]))
    D2 = DiagonalOperator(np.array([3.0, 0.75, 1.5, 6.0, 0.375, 2.25, 9.0]))
    P3 = PermutationOperator(np.array([2, 0, 1, 4, 5, 3, 6]))  # NON-involution
    DS = DiagonalOperator(np.full(7, 0.3))
    return q, D1, D2, P3, DS


def test_product_solve_reroute_executes_and_matches_dense(monkeypatch):
    """Spec §40.1 — the four legs on the NON-COMMUTING headline product.

    (ii) SENTINEL: ``(D@P).inverse().apply(q)`` routes InverseOperator →
    the re-routed ``OperatorProduct.solve`` — proven EXECUTED (Mode 11),
    because a value gate spelled ``.inverse().apply`` both sides would be
    tautological. (iii) INDEPENDENCE: the dense ``np.linalg.solve`` on
    the materialized product is structurally independent of the
    re-route. (i)/(iv) live in the factor-kind matrix below.
    """
    q, D1, _, P3, _ = _reroute_fixtures()
    prod = OperatorProduct(D1, P3)  # NON-COMMUTING (config-blindness §0.6)

    calls: list[int] = []
    real = OperatorProduct.solve
    monkeypatch.setattr(
        OperatorProduct, "solve",
        lambda self, b: (calls.append(1), real(self, b))[1],
    )
    out = prod.inverse().apply(q)
    assert calls, "OperatorProduct.solve NEVER executed — Mode-11 vacuous gate"

    dense = np.linalg.solve(prod.as_matrix(basis_shape=(7,)), q)
    np.testing.assert_allclose(out, dense, rtol=1e-12)


def test_product_solve_reroute_inherits_the_precarve_baseline():
    """Spec §40.2 — the factor-kind matrix vs the PRE-carve snapshot.

    Each row exercises one factor KIND through the re-routed recursion
    (``b.inverse().apply(a.inverse().apply(q))``) against the value the
    OLD spelling (``b.solve(a.solve(q))``) produced on the pre-carve
    tree. The four DIRECT rows are bit-identical by construction
    (identical delegated realizations) → ``assert_array_equal``; the
    Green-factor row runs the ITERATIVE splitting driver, where
    cross-run FP jitter forbids a bit-identity claim (lessons L7) →
    driver-tolerance ``assert_allclose``. Do NOT tighten the sum row to
    array_equal; do NOT loosen the direct rows.
    """
    q, D1, D2, P3, DS = _reroute_fixtures()
    baseline = np.load(_REROUTE_BASELINE)
    np.testing.assert_array_equal(q, baseline["q"])  # fixture drift guard

    direct_rows = {
        "diag": OperatorProduct(D2, D1),
        "permutation": OperatorProduct(D1, P3),
        "scaled": OperatorProduct(ScaledOperator(2.5, D2), D1),
        "identity": OperatorProduct(IdentityOperator(), D1),
    }
    for name, prod in direct_rows.items():
        np.testing.assert_array_equal(
            prod.solve(q), baseline[name],
            err_msg=f"factor-kind row {name!r} drifted from the pre-carve baseline",
        )
    # The Green-factor row: an Identity-headed sum — the factor whose own
    # solve retired; its inverse is the GreenOperator (this row also
    # end-to-end pins the algebra-closed-head seeding fix in
    # _seeded_inverse, which this exact fixture exposed).
    green_prod = OperatorProduct(IdentityOperator() + DS, D1)
    np.testing.assert_allclose(
        green_prod.solve(q), baseline["sum_green"], rtol=1e-8, atol=0.0,
    )
