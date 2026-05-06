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
