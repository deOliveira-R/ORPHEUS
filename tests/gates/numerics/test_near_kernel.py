r"""U3b rows W1-W6 (spec ``scratch/membership_covector/u3b/gate_spec.md`` §2.1): the NearKernel witness on domain-free matrices.

The primitive is pinned against matrices with CHOSEN singular structure (X4: a verified primitive is trusted against a
domain-free closed form, never only an in-domain solver); its SN and diffusion rows are in their own packages."""
from __future__ import annotations

import inspect

import numpy as np
import pytest

from scipy.sparse import csr_array

from orpheus.numerics.assembled_operator import SparseAssembledOperator
from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator
from orpheus.numerics.near_kernel import DENSE_LIMIT, NearKernel

from orpheus.numerics.space import FunctionSpace

pytestmark = [pytest.mark.l0]
_EPS = np.finfo(float).eps


def DenseOperator(matrix, space):
    """A bound operator on ``space`` from a dense matrix (the witness reads only its action and its ends)."""
    return SparseAssembledOperator(csr_array(matrix), domain=space, codomain=space)


class _Invertible(SparseAssembledOperator):
    """An implicit part that inverts by its LU (``MatrixInverseOperator``), as the SN sweep inverts ``L + C``."""

    @property
    def is_invertible(self) -> bool:
        return True

    def inverse(self):
        return MatrixInverseOperator(self)


def Implicit(matrix, space):
    return _Invertible(csr_array(matrix), domain=space, codomain=space)


def dimension(w): return w.dimension
def right(w): return w.right
def left(w): return w.left
def sigma_max(w): return w.sigma_max
def singular_gap(w): return w.singular_gap
def splitting_gap(w): return w.splitting_gap


def _orthogonal(n, seed):
    return np.linalg.qr(np.random.default_rng(seed).standard_normal((n, n)))[0]


def _chosen(n=40, k=3, seed=5):
    U, V = _orthogonal(n, seed), _orthogonal(n, seed + 1)
    s = np.concatenate([np.geomspace(10.0, 0.05, n - k), np.zeros(k)])
    return U @ np.diag(s) @ V.T, U, s, V


def _space(n, g=None):
    return FunctionSpace("w", (n,)) if g is None else FunctionSpace("w", (n,), inner_product_weights=g)


def _span_residual(B, Z):
    if Z.shape[1] == 0:
        return 0.0
    Q = np.linalg.qr(B)[0]
    return float(np.linalg.norm(Z - Q @ (Q.T @ Z)) / np.linalg.norm(Z))


@pytest.mark.parametrize("k", [0, 1, 3, 7])
def test_W1_dense_route_on_chosen_singular_values(k):
    """THEOREM: dimension k, right = V[:, -k:], left = U[:, -k:] (G = I), gap >= 1e10. Arm rank_rtol_1e-2 reds (k+1)."""
    A, U, s, V = _chosen(k=k)
    w = NearKernel.dense(DenseOperator(A, _space(A.shape[0])))
    assert dimension(w) == k
    if k:
        assert _span_residual(V[:, -k:], right(w)) <= 64 * _EPS * 40
        assert _span_residual(U[:, -k:], left(w)) <= 64 * _EPS * 40
        gap = singular_gap(w)
        assert gap is not None and gap >= 1e10
    assert abs(sigma_max(w) / s[0] - 1.0) <= 64 * _EPS


def test_W2_left_vectors_in_H_coordinates_and_the_ker_G_quotient():
    """G diagonal 0.5-4 with three exact zeros (ker G): left = G^+ U[:, -k:]; a direction supported only on ker G is not
    counted. Arms left_euclidean, no_quotient."""
    n, k = 40, 3
    rng = np.random.default_rng(9)
    g = rng.uniform(0.5, 4.0, n)
    g[[4, 17, 31]] = 0.0
    A, U, s, V = _chosen(n=n, k=k)
    A[:, [4, 17, 31]] = 0.0  # ker G columns: e_4, e_17, e_31 are raw null vectors, outside the Hilbert space
    A[[4, 17, 31], :] = 0.0
    w = NearKernel.dense(DenseOperator(A, _space(n, g)))
    Ud, sd, Vtd = np.linalg.svd(A)
    raw = int(np.sum(sd < 1e-10 * sd[0]))
    pos = g > 0
    Rd = Vtd[n - raw:].T
    k_mod = int(np.sum(np.linalg.svd(Rd[pos], compute_uv=False) > 1e-8))
    assert raw > k_mod and dimension(w) == k_mod
    gp = np.divide(1.0, g, out=np.zeros_like(g), where=pos)
    LH = gp[:, None] * Ud[:, n - raw:]
    LH = LH[:, np.linalg.norm(LH, axis=0) > 1e-12]
    assert _span_residual(LH, left(w)) <= 1e-12


def _splitting(n=40, k=3, minus_one=3, seed=13):
    """M lower triangular dominant, N = M - A with M^-1 N spectrum {1 (k), -1 (minus_one), |lambda| <= 0.9}."""
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((n, n)) + n * np.eye(n)
    lam = np.concatenate([np.ones(k), -np.ones(minus_one), rng.uniform(-0.9, 0.9, n - k - minus_one)])
    K = X @ np.diag(lam) @ np.linalg.inv(X)          # the iteration map M^-1 N
    M = np.tril(rng.standard_normal((n, n))) + 4 * n * np.eye(n)
    N = M @ K
    return M, N, M - N, lam, X


@pytest.mark.parametrize("k", [1, 3, 14])
def test_W3_block_route_averages_away_the_minus_one_modes(k):
    """The block route on P = (I + M^-1 N)/2 returns exactly k (not k + 3), right span = X[:, :k], left span = ker A^T,
    and splitting_gap = 1 - |(1 + lambda_{k+1})/2|. Arms no_averaging (k+3), adjoint_order (left O(1)),
    fixed_block_1 (the positive control: 1 found for k = 3, 14)."""
    M, N, A, lam, X = _splitting(k=k, n=60)
    space = _space(60)
    w = NearKernel.of_splitting(Implicit(M, space), DenseOperator(N, space))
    assert dimension(w) == k
    assert _span_residual(X[:, :k], right(w)) <= 1e-10
    Ul = np.linalg.svd(A)[0][:, -k:]
    assert _span_residual(Ul, left(w)) <= 1e-10
    rest = np.sort(np.abs((1.0 + lam[k:]) / 2.0))[::-1]
    measured_gap = splitting_gap(w)
    assert measured_gap is not None
    np.testing.assert_allclose(measured_gap, 1.0 - rest[0], rtol=1e-6)


def test_W4_the_block_route_is_predictor_blind_by_signature():
    """No parameter names a dimension, projector or kernel (a signature-tautological gate otherwise, vv mode 8(3))."""
    names = set(inspect.signature(NearKernel.of_splitting).parameters)
    assert not names & {"k", "dimension", "kernel", "projector", "predicted", "block_size"}, names


def test_W5_own_sigma_max_is_the_largest_singular_value_not_the_spectral_radius():
    """Power on A^T A within 1 % (`[M]` p9 4.8e-4 / 6.3e-5 on SN slabs); the rho(A) spelling reads 2-5x low (arm)."""
    M, N, A, *_ = _splitting(k=3, n=60)
    space = _space(60)
    w = NearKernel.of_splitting(Implicit(M, space), DenseOperator(N, space))
    assert abs(sigma_max(w) / np.linalg.svd(A, compute_uv=False)[0] - 1.0) <= 1e-2


def test_W6_measure_dispatches_on_the_dense_limit(monkeypatch):
    """Both sides of the branch (vv #26): a counting spy per route; DENSE_LIMIT read from the module."""
    calls = {"dense": 0, "split": 0}
    real_dense, real_split = NearKernel.dense, NearKernel.of_splitting

    def counting(route, real):
        def spy(cls, *args, **kwargs):
            calls[route] += 1
            return real(*args, **kwargs)
        return classmethod(spy)

    monkeypatch.setattr(NearKernel, "dense", counting("dense", real_dense))
    monkeypatch.setattr(NearKernel, "of_splitting", counting("split", real_split))
    small = 40
    assert small <= DENSE_LIMIT
    M, N, *_ = _splitting(n=small)
    NearKernel.measure(Implicit(M, _space(small)), DenseOperator(N, _space(small)))
    assert calls == {"dense": 1, "split": 0}


def test_W6b_above_the_limit_the_block_route_runs(monkeypatch):
    import orpheus.numerics.near_kernel as nk
    monkeypatch.setattr(nk, "DENSE_LIMIT", 10)
    M, N, *_ = _splitting(n=40)
    w = nk.NearKernel.measure(Implicit(M, _space(40)), DenseOperator(N, _space(40)))
    assert dimension(w) == 3 and splitting_gap(w) is not None
