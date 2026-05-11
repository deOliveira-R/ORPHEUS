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
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    CAP_SOLVE,
    IdentityOperator,
    LinearOperator,
    LinearOperatorMixin,
    MissingCapability,
    OperatorProduct,
    OperatorSum,
    ScaledOperator,
    ZeroOperator,
    as_scipy_linop,
)

pytestmark = pytest.mark.foundation


# ───────────────────────────────────────────────────────────────────────
# Fixtures: numpy-matrix-backed LinearOperator implementations
# ───────────────────────────────────────────────────────────────────────


class MatrixOperator(LinearOperatorMixin):
    """Test operator backed by a dense numpy matrix.

    Advertises capabilities according to constructor flags. Used to
    drive composition + capability tests with a fully analytic ground
    truth (the matrix is stored, so any test can compare against
    direct ``M @ x``). Inherits :class:`LinearOperatorMixin` so that
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
        caps = {CAP_APPLY}
        if can_solve:
            caps.add(CAP_SOLVE)
        if can_transpose:
            caps.add(CAP_APPLY_TRANSPOSE)
        self.capabilities = frozenset(caps)

    def apply(self, x: np.ndarray) -> np.ndarray:
        return self.matrix @ x

    def solve(self, b: np.ndarray) -> np.ndarray:
        return np.linalg.solve(self.matrix, b)

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
# Capability propagation — OperatorSum
# ───────────────────────────────────────────────────────────────────────


def test_sum_apply_propagates_with_both(matrix_full, matrix_apply_only):
    s = OperatorSum(matrix_full, matrix_apply_only)
    assert CAP_APPLY in s.capabilities


def test_sum_solve_does_not_propagate(matrix_full):
    """No general algorithm for (A+B)^{-1} from A^{-1}, B^{-1}."""
    s = OperatorSum(matrix_full, matrix_full)
    assert CAP_SOLVE not in s.capabilities


def test_sum_transpose_propagates_with_both(matrix_full):
    s = OperatorSum(matrix_full, matrix_full)
    assert CAP_APPLY_TRANSPOSE in s.capabilities


def test_sum_transpose_drops_when_one_lacks(matrix_full, matrix_apply_only):
    s = OperatorSum(matrix_full, matrix_apply_only)
    assert CAP_APPLY_TRANSPOSE not in s.capabilities


# ───────────────────────────────────────────────────────────────────────
# Capability propagation — OperatorProduct
# ───────────────────────────────────────────────────────────────────────


def test_product_apply_propagates(matrix_full, matrix_apply_only):
    p = OperatorProduct(matrix_full, matrix_apply_only)
    assert CAP_APPLY in p.capabilities


def test_product_solve_propagates_with_both(matrix_full):
    p = OperatorProduct(matrix_full, matrix_full)
    assert CAP_SOLVE in p.capabilities


def test_product_solve_drops_when_one_lacks(matrix_full, matrix_apply_only):
    p = OperatorProduct(matrix_full, matrix_apply_only)
    assert CAP_SOLVE not in p.capabilities


def test_product_transpose_propagates_with_both(matrix_full):
    p = OperatorProduct(matrix_full, matrix_full)
    assert CAP_APPLY_TRANSPOSE in p.capabilities


# ───────────────────────────────────────────────────────────────────────
# Capability propagation — ScaledOperator
# ───────────────────────────────────────────────────────────────────────


def test_scaled_preserves_all_capabilities(matrix_full):
    s = ScaledOperator(2.5, matrix_full)
    # Scalar multiplication is unitary on the cap set.
    assert s.capabilities == matrix_full.capabilities


def test_scaled_apply_only_stays_apply_only(matrix_apply_only):
    s = ScaledOperator(0.5, matrix_apply_only)
    assert s.capabilities == frozenset({CAP_APPLY})


def test_scaled_zero_scalar_rejected(matrix_full):
    with pytest.raises(ValueError, match="ZeroOperator"):
        ScaledOperator(0.0, matrix_full)


# ───────────────────────────────────────────────────────────────────────
# Capability propagation — Identity / Zero
# ───────────────────────────────────────────────────────────────────────


def test_identity_full_capabilities():
    cap = IdentityOperator().capabilities
    assert cap == frozenset({CAP_APPLY, CAP_SOLVE, CAP_APPLY_TRANSPOSE})


def test_zero_lacks_solve():
    """Critical: a ZeroOperator MUST NOT advertise solve.

    Inverting zero is meaningless; forcing a stub here would corrupt
    composer capability propagation downstream.
    """
    cap = ZeroOperator().capabilities
    assert CAP_SOLVE not in cap
    assert CAP_APPLY in cap
    assert CAP_APPLY_TRANSPOSE in cap


# ───────────────────────────────────────────────────────────────────────
# Negative tests: composition-time MissingCapability
# ───────────────────────────────────────────────────────────────────────


class NoApplyOperator:
    """Pathological operator: declares an empty capability set."""

    capabilities: frozenset[str] = frozenset()

    def apply(self, x: np.ndarray) -> np.ndarray:  # pragma: no cover
        # Never called — composers reject it at construction time.
        return x


def test_sum_rejects_missing_apply_at_composition(matrix_full):
    bad = NoApplyOperator()
    with pytest.raises(MissingCapability, match="apply"):
        OperatorSum(matrix_full, bad)
    with pytest.raises(MissingCapability, match="apply"):
        OperatorSum(bad, matrix_full)


def test_product_rejects_missing_apply_at_composition(matrix_full):
    bad = NoApplyOperator()
    with pytest.raises(MissingCapability, match="apply"):
        OperatorProduct(matrix_full, bad)
    with pytest.raises(MissingCapability, match="apply"):
        OperatorProduct(bad, matrix_full)


def test_scaled_rejects_missing_apply_at_composition():
    bad = NoApplyOperator()
    with pytest.raises(MissingCapability, match="apply"):
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


def test_scaled_solve_divides(matrix_full, vector):
    """(αL)^{-1} = (1/α) L^{-1}."""
    s = ScaledOperator(2.5, matrix_full)
    expected = (1.0 / 2.5) * np.linalg.solve(matrix_full.matrix, vector)
    np.testing.assert_allclose(s.solve(vector), expected, rtol=1e-10)


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
    # apply propagates (everyone has apply).
    assert CAP_APPLY in op.capabilities
    # solve does NOT propagate (sum kills it).
    assert CAP_SOLVE not in op.capabilities
    x = np.ones(L.matrix.shape[0])
    expected = L.matrix @ x - S.matrix @ x - F.matrix @ x
    np.testing.assert_allclose(op.apply(x), expected, rtol=1e-12)


# ───────────────────────────────────────────────────────────────────────
# scipy interop
# ───────────────────────────────────────────────────────────────────────


def test_as_scipy_linop_matvec_matches_apply(matrix_full, vector):
    n = matrix_full.matrix.shape[0]
    sp = as_scipy_linop(matrix_full, shape=(n, n))
    np.testing.assert_allclose(
        sp.matvec(vector), matrix_full.apply(vector), rtol=1e-12
    )


def test_as_scipy_linop_rmatvec_when_transpose_capable(matrix_full, vector):
    n = matrix_full.matrix.shape[0]
    sp = as_scipy_linop(matrix_full, shape=(n, n))
    np.testing.assert_allclose(
        sp.rmatvec(vector), matrix_full.apply_transpose(vector), rtol=1e-12
    )


def test_as_scipy_linop_no_rmatvec_when_not_transpose_capable(
    matrix_apply_only, vector
):
    n = matrix_apply_only.matrix.shape[0]
    sp = as_scipy_linop(matrix_apply_only, shape=(n, n))
    # scipy raises when rmatvec is not provided. We don't pin scipy's
    # exact exception type — only that calling rmatvec fails.
    with pytest.raises(Exception):  # noqa: BLE001 - scipy version-dependent
        sp.rmatvec(vector)


def test_as_scipy_linop_rejects_missing_apply():
    bad = NoApplyOperator()
    with pytest.raises(MissingCapability, match="apply"):
        as_scipy_linop(bad, shape=(3, 3))


def test_as_scipy_linop_works_with_composite(matrix_full, matrix_apply_only):
    """End-to-end: build (L - S) via algebra, hand to scipy."""
    composite = matrix_full - matrix_apply_only
    n = matrix_full.matrix.shape[0]
    sp = as_scipy_linop(composite, shape=(n, n))
    x = np.ones(n)
    expected = matrix_full.matrix @ x - matrix_apply_only.matrix @ x
    np.testing.assert_allclose(sp.matvec(x), expected, rtol=1e-12)


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
#     domain/range mismatch in the composers.
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


class _SpacedMatrixOperator(LinearOperatorMixin):
    """MatrixOperator carrying explicit domain / range function spaces.

    Used for the weight-aware Hilbert-adjoint identity tests, where the
    inner product on the range carries non-trivial weights and the
    adjoint must invert them. Mirrors ``MatrixOperator`` but exposes
    the function-space metadata shipped in Phase B.
    """

    def __init__(
        self,
        matrix: np.ndarray,
        domain_space: FunctionSpace,
        range_space: FunctionSpace,
    ) -> None:
        self.matrix = np.asarray(matrix, dtype=float)
        self._domain = domain_space
        self._range = range_space
        self.capabilities = frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})

    @property
    def domain(self) -> FunctionSpace:
        return self._domain

    @property
    def range(self) -> FunctionSpace:
        return self._range

    def apply(self, x: np.ndarray) -> np.ndarray:
        return self.matrix @ x

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        return self.matrix.T @ x


def test_adjoint_weighted_identity(rng):
    """Hilbert-adjoint identity holds with quadrature-weight metadata.

    Range inner product carries non-trivial diagonal weights ``w_W``;
    the adjoint must produce an output that, evaluated under the
    domain inner product, balances ``<A u, v>_W``.
    """
    n_dom = 4
    n_rng = 4
    M = rng.standard_normal((n_rng, n_dom))
    w_W = np.array([1.0, 2.0, 3.0, 4.0])  # range weights (non-trivial)
    domain = FunctionSpace(name="V", shape=(n_dom,))  # Euclidean
    range_ = FunctionSpace(
        name="W", shape=(n_rng,), inner_product_weights=w_W,
    )
    A = _SpacedMatrixOperator(M, domain, range_)
    u = rng.standard_normal(n_dom)
    v = rng.standard_normal(n_rng)

    Au = A.apply(u)                     # in W
    Hv = A.H.apply(v)                   # in V
    lhs = range_.inner_product(Au, v)
    rhs = domain.inner_product(u, Hv)
    assert np.isclose(lhs, rhs, rtol=1e-12)


def test_adjoint_swaps_domain_and_range():
    """``A.H.domain`` is ``A.range`` and vice versa."""
    n_dom = 3
    n_rng = 4
    M = np.random.default_rng(0).standard_normal((n_rng, n_dom))
    domain = FunctionSpace(name="V", shape=(n_dom,))
    range_ = FunctionSpace(name="W", shape=(n_rng,))
    A = _SpacedMatrixOperator(M, domain, range_)
    assert A.H.domain == range_
    assert A.H.range == domain


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


def test_sum_rejects_incompatible_ranges():
    n = 3
    M = np.eye(n)
    V = FunctionSpace(name="V", shape=(n,))
    W1 = FunctionSpace(name="W1", shape=(n,))
    W2 = FunctionSpace(name="W2", shape=(n,))
    A = _SpacedMatrixOperator(M, V, W1)
    B = _SpacedMatrixOperator(M, V, W2)
    with pytest.raises(IncompatibleOperatorComposition, match="range"):
        _ = A + B


def test_product_rejects_incompatible_inner_space():
    """``A @ B`` requires ``A.domain == B.range``."""
    n = 3
    M = np.eye(n)
    V = FunctionSpace(name="V", shape=(n,))
    W = FunctionSpace(name="W", shape=(n,))
    Z = FunctionSpace(name="Z", shape=(n,))
    A = _SpacedMatrixOperator(M, W, V)  # A: W → V
    B = _SpacedMatrixOperator(M, V, Z)  # B: V → Z; range Z != A.domain W
    with pytest.raises(
        IncompatibleOperatorComposition, match="A @ B|domain|range",
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
    s = A + B  # same domain V, same range W: OK
    assert s.domain == V
    assert s.range == W


def test_composition_skips_check_when_either_is_none(matrix_full):
    """Backward-compat: operators without function spaces pass freely.

    ``MatrixOperator`` (the existing fixture) returns ``None`` for
    domain and range. Sums and products against another ``None``
    operator must not raise — the composition compatibility check
    only fires when BOTH operators carry function spaces.
    """
    sum_ok = matrix_full + matrix_full
    assert sum_ok is not None
    prod_ok = matrix_full @ matrix_full
    assert prod_ok is not None


def test_repr_format(matrix_full):
    """``repr(op)`` carries the class name, domain/range, and caps."""
    r = repr(matrix_full)
    assert "MatrixOperator" in r
    # MatrixOperator has no function-space tagging; expect '?' placeholders.
    assert "domain='?'" in r
    assert "range='?'" in r
    # Capabilities list.
    assert "apply" in r


def test_repr_with_function_spaces_shows_names():
    n = 3
    M = np.eye(n)
    V = FunctionSpace(name="phi", shape=(n,))
    W = FunctionSpace(name="psi", shape=(n,))
    A = _SpacedMatrixOperator(M, V, W)
    r = repr(A)
    assert "'phi'" in r
    assert "'psi'" in r
