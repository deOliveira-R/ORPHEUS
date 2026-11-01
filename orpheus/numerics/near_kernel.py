r"""The near-kernel WITNESS: an operator's spectrum at zero, measured.

A posing's Fredholm verdict is decided from the kernel its method PREDICTS
from declared structure (:class:`~orpheus.numerics.fredholm.PredictedKernel`),
never from a measurement: a production verdict must not depend on a threshold
chosen against round-off. This module is the measurement that CHECKS the
prediction. It returns the right and left near-null vectors of an operator
:math:`A`, the count of them, and the gap that separates them from the rest of
the spectrum, so a gate can assert "measured dimension = predicted dimension,
and the predicted spans contain the measured vectors". It decides nothing in
production.

Two routes, one result
----------------------

* **Dense** (:meth:`NearKernel.dense`), for :math:`n \le` :data:`DENSE_LIMIT`:
  the singular value decomposition of the assembled matrix. The near-null
  vectors are the singular vectors whose singular value is at most
  :data:`RANK_RTOL` :math:`\sigma_{\max}`, and the gap is the ratio of the
  first singular value above that threshold to the largest below it.
* **Block** (:meth:`NearKernel.of_splitting`), at scale: a splitting
  :math:`A = M - N` whose implicit part inverts cheaply (the SN sweep) gives
  the map :math:`P = (I + M^{-1}N)/2`. Its eigenvalue-1 space is exactly
  :math:`\ker A` (:math:`Px = x \iff M^{-1}Nx = x \iff Ax = 0`), and every
  other eigenvalue of :math:`P` lies strictly inside the unit disc when the
  source iteration :math:`M^{-1}N` is a contraction off its kernel, because
  the averaging sends the iteration's :math:`\lambda = -1` modes (the Jacobi
  boundary's mirror-odd modes) to :math:`0`. Block subspace iteration with a
  Rayleigh–Ritz step resolves a REPEATED eigenvalue 1 (a single Krylov vector
  holds one copy of it: `[M]` ARPACK found 7 of 14 on the LS4 box), and the
  block grows until the converged count stops filling it. The left vectors
  are the fixed points of the DUAL split's map
  :math:`(I + M^{-\mathsf T}N^{\mathsf T})/2`, which are :math:`\ker A^{\mathsf T}`
  itself (:math:`A^{\mathsf T}y = 0 \iff M^{\mathsf T}y = N^{\mathsf T}y`);
  the transpose of the forward map has different fixed points.

Why not a matrix-free singular value route: `[M]` 2026-10-31 on the LS4
:math:`ng = 2` box (14-dimensional kernel), Lanczos on :math:`A^{\mathsf T}A`
did not converge in 20 000 iterations, LOBPCG stalled at the
:math:`\sqrt\varepsilon` floor (a gap of about :math:`10^5`, not
:math:`10^{10}`), and PROPACK returned 7 of the 14 clustered vectors
(``scratch/membership_covector/u3b/iterative_probe.py``).

The witness is PREDICTOR-BLIND: no parameter names a dimension, a projector or
a kernel, the block size is its own, and :math:`\sigma_{\max}` comes from its
own power iteration on :math:`A^{\mathsf T}A` (X4: a witness sized by the
prediction could only confirm it).

The two quotients
-----------------

* **Modulo** :math:`\ker G`. A coordinate whose metric weight is zero (the
  tangential trace slots of a closed face) lies outside the Hilbert space; a
  near-null vector supported only there is not an element of it, and the
  counts are taken on the :math:`G > 0` coordinates.
* **The left coordinates.** A left null vector :math:`y` of :math:`A` is a
  co-vector; its Riesz representative :math:`G^{+}y` is what
  :class:`~orpheus.numerics.fredholm.PredictedKernel.H` holds, so :attr:`NearKernel.left`
  is returned in those coordinates and a containment check is a plain
  comparison.

The metric must be diagonal: the one declared door,
:func:`~orpheus.numerics.kernel_projection_gauge.diagonal_metric`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.assembled_operator import dense as assembled_dense
from orpheus.numerics.kernel_projection_gauge import diagonal_metric
from orpheus.numerics.operator import invertible

if TYPE_CHECKING:
    from orpheus.numerics.operator import LinearOperator
    from orpheus.numerics.space import FunctionSpace

__all__ = ["DENSE_LIMIT", "RANK_RTOL", "NearKernel"]

#: The largest flat dimension the dense route takes. `[M]` 2026-10-31, the SN
#: thick boxes: assembly plus SVD took 14 s at n = 2304 and 217 s at n = 7680.
DENSE_LIMIT = 3000

#: The relative threshold below which a singular value counts as zero.
RANK_RTOL = 1e-10

#: Relative singular value below which a basis direction has no component on
#: the G > 0 coordinates (it lives in ker G only).
_QUOTIENT_RTOL = 1e-8

#: The block route's defaults: the stopping residual of a converged Ritz pair,
#: how close to 1 a Ritz value must be to count, the starting block, the slack
#: the block keeps above the count, and the sweep budget.
_STOP_RTOL = 1e-12
_ONE_TOL = 1e-8
_BLOCK = 8
_SLACK = 4
_MAX_ITERATIONS = 20_000
_CHECK_EVERY = 10
_SEED = 20261031


@dataclass(frozen=True)
class NearKernel:
    r"""The measured near-kernel of an operator: its right and left near-null vectors and the gap behind them.

    ``right`` is ``(n, k)``, an orthonormal basis of the measured
    :math:`\ker A` modulo :math:`\ker G`; ``left`` is ``(n, k_left)``, an
    orthonormal basis of :math:`G^{+}\ker A^{\mathsf T}` (the coordinates of
    :class:`~orpheus.numerics.fredholm.PredictedKernel.H`). ``singular_gap`` is
    the dense route's :math:`\sigma_{k+1}/\sigma_k` (``inf`` when the null
    singular values are exactly zero; ``None`` with no kernel or on the block
    route); ``splitting_gap`` is the block route's :math:`1 - |w_{k+1}|`, the
    averaged map's distance from 1 to the rest of its spectrum (a property of
    the splitting as well as of :math:`A`); ``stop_residual`` is the route's
    own reading of how null the returned vectors are,
    :math:`\max\lVert A z\rVert/\sigma_{\max}` over the right vectors.
    ``sigma_min`` is the dense route's smallest singular value ABOVE the
    kernel, :math:`\sigma_{k+1}` (the smallest one when the kernel is
    trivial): the operator's conditioning on the kernel's complement is
    :math:`\sigma_{k+1}/\sigma_{\max}`. The block route does not measure it
    (``None``).
    """

    space: "FunctionSpace"
    right: NDArray
    left: NDArray
    sigma_max: float
    sigma_min: float | None
    singular_gap: float | None
    splitting_gap: float | None
    stop_residual: float
    route: str

    @property
    def dimension(self) -> int:
        r"""The measured kernel's dimension modulo :math:`\ker G`."""
        return int(self.right.shape[1])

    # ── the dense route ────────────────────────────────────────────────

    @classmethod
    def dense(cls, op: "LinearOperator") -> "NearKernel":
        r"""The dense route: the SVD of ``op`` assembled on its flat chart."""
        space = _domain(op)
        weights = diagonal_metric(space)
        n = weights.size
        A = assembled_dense(op, basis_shape=(n,), max_dimension=max(DENSE_LIMIT, n))
        U, s, Vt = np.linalg.svd(A)
        sigma_max = float(s[0]) if s.size else 0.0
        null = s <= RANK_RTOL * sigma_max
        k = int(np.count_nonzero(null))
        right = _modulo_ker_g(Vt[n - k:].T, weights)
        left = _left_coordinates(U[:, n - k:], weights)
        if k == 0 or k == n:
            gap = None
        else:
            smallest_null = float(s[n - k])
            gap = float("inf") if smallest_null == 0.0 else float(s[n - k - 1]) / smallest_null
        residual = _residual(lambda x: A @ x, right, sigma_max)
        sigma_min = float(s[n - k - 1]) if k < n else None
        return cls(space, right, left, sigma_max, sigma_min, gap, None, residual, "dense")

    # ── the block route ────────────────────────────────────────────────

    @classmethod
    def of_splitting(cls, implicit: "LinearOperator", lagged: "LinearOperator") -> "NearKernel":
        r"""The block route on the splitting :math:`A = M - N` (``implicit`` :math:`M`, ``lagged`` :math:`N`).

        ``implicit`` must invert (its ``inverse()`` is the sweep) and both parts
        must transpose (``dual()``), so the dual split's map exists: on SN the
        multi-D transposed sweep, #616. The operators act on flat coordinates
        of their domain's chart (wrap a carrier-valued operator in
        :class:`~orpheus.numerics.flat_operator.FlattenedOperator`).
        """
        space = _domain(implicit)
        weights = diagonal_metric(space)
        if not invertible(implicit):
            raise ValueError(
                f"NearKernel.of_splitting: the implicit part {type(implicit).__name__} does not invert, so the "
                f"averaged splitting map does not exist."
            )
        inverse = implicit.inverse()

        def forward(x: NDArray) -> NDArray:
            return 0.5 * (x + _flat(inverse.apply(_flat(lagged.apply(x)))))

        inverse_dual, lagged_dual = inverse.dual(), lagged.dual()

        def dual(y: NDArray) -> NDArray:
            return 0.5 * (y + _flat(inverse_dual.apply(_flat(lagged_dual.apply(y)))))

        def loss(x: NDArray) -> NDArray:
            return _flat(implicit.apply(x)) - _flat(lagged.apply(x))

        def loss_transpose(y: NDArray) -> NDArray:
            return _flat(implicit.dual().apply(y)) - _flat(lagged_dual.apply(y))

        n = weights.size
        right_raw, gap, _ = _eigenvalue_one_space(forward, n)
        left_raw, _, _ = _eigenvalue_one_space(dual, n)
        sigma_max = _largest_singular_value(loss, loss_transpose, n)
        right = _modulo_ker_g(right_raw, weights)
        left = _left_coordinates(left_raw, weights)
        residual = _residual(loss, right, sigma_max)
        return cls(space, right, left, sigma_max, None, None, gap, residual, "block")

    # ── the dispatch ───────────────────────────────────────────────────

    @classmethod
    def measure(cls, implicit: "LinearOperator", lagged: "LinearOperator") -> "NearKernel":
        r"""The route by size: dense on :math:`A = M - N` up to :data:`DENSE_LIMIT`, the block route above it."""
        n = diagonal_metric(_domain(implicit)).size
        if n <= DENSE_LIMIT:
            return cls.dense(implicit - lagged)
        return cls.of_splitting(implicit, lagged)


# ═══════════════════════════════════════════════════════════════════════
# The pieces
# ═══════════════════════════════════════════════════════════════════════


def _domain(op: "LinearOperator") -> "FunctionSpace":
    domain = op.domain
    if domain is None:
        raise ValueError(
            f"NearKernel: {type(op).__name__} has no bound domain, so its metric (the quotient by ker G and the "
            f"left coordinates) is unknown."
        )
    return domain


def _flat(x: object) -> NDArray:
    return np.asarray(x, dtype=float).reshape(-1)


def _raise(left: NDArray, weights: NDArray) -> NDArray:
    r"""The Riesz representatives :math:`G^{+}y` of left null vectors (zero on :math:`\ker G`)."""
    inverse = np.divide(1.0, weights, out=np.zeros_like(weights), where=weights > 0.0)
    return inverse[:, None] * left


def _modulo_ker_g(basis: NDArray, weights: NDArray) -> NDArray:
    r"""The directions of an ORTHONORMAL ``basis`` with a component on the :math:`G > 0` coordinates, orthonormal.

    Judged on the orthonormal basis, so the tolerance is absolute: a direction
    whose :math:`G > 0` part is below :data:`_QUOTIENT_RTOL` lies in
    :math:`\ker G` (a near-null vector of round-off size there is not an
    element of the Hilbert space).
    """
    if basis.shape[1] == 0:
        return basis
    _, sv, vt = np.linalg.svd(basis[weights > 0.0], full_matrices=False)
    keep = sv > _QUOTIENT_RTOL
    return np.linalg.qr(basis @ vt[keep].T)[0] if np.any(keep) else basis[:, :0]


def _left_coordinates(left: NDArray, weights: NDArray) -> NDArray:
    r"""The Riesz representatives :math:`G^{+}y` of the left null vectors modulo :math:`\ker G`, orthonormal."""
    kept = _modulo_ker_g(left, weights)
    return np.linalg.qr(_raise(kept, weights))[0] if kept.shape[1] else kept


def _residual(action, basis: NDArray, sigma_max: float) -> float:
    if basis.shape[1] == 0 or sigma_max == 0.0:
        return 0.0
    return max(float(np.linalg.norm(action(z))) for z in basis.T) / sigma_max


def _apply_columns(action, X: NDArray) -> NDArray:
    return np.column_stack([action(x) for x in X.T])


def _eigenvalue_one_space(action, n: int) -> tuple[NDArray, float, float]:
    r"""The eigenvalue-1 space of a map whose other eigenvalues lie strictly inside the unit disc.

    Block subspace iteration with a Rayleigh–Ritz step every
    :data:`_CHECK_EVERY` applications. A Ritz pair counts when its value is
    within :data:`_ONE_TOL` of 1 and its residual is below :data:`_STOP_RTOL`;
    the block doubles whenever the count reaches the block minus
    :data:`_SLACK` (the block, not the caller, sizes itself). It stops when
    the count is unchanged over two checks with every counted pair converged.
    Returns the basis, the gap :math:`1 - |w_{k+1}|` to the largest Ritz value
    not counted, and the stopping residual.
    """
    rng = np.random.default_rng(_SEED)
    block = min(_BLOCK, n)
    X = np.linalg.qr(rng.standard_normal((n, block)))[0]
    previous = -1
    for iteration in range(1, _MAX_ITERATIONS + 1):
        X = np.linalg.qr(_apply_columns(action, X))[0]
        if iteration % _CHECK_EVERY:
            continue
        PX = _apply_columns(action, X)
        values, vectors = np.linalg.eig(X.T @ PX)
        order = np.argsort(-np.abs(values))
        values, vectors = values[order], vectors[:, order]
        ritz = X @ vectors
        residuals = np.linalg.norm(PX @ vectors - ritz * values, axis=0) / np.linalg.norm(ritz, axis=0)
        near_one = np.abs(values - 1.0) <= _ONE_TOL
        count = int(np.count_nonzero(near_one))
        if count > block - _SLACK and block < n:
            grown = min(2 * block, n)
            X = np.linalg.qr(np.column_stack([X, rng.standard_normal((n, grown - block))]))[0]
            block, previous = grown, -1
            continue
        converged = bool(np.all(residuals[near_one] <= _STOP_RTOL))
        # The certificate that no eigenvalue 1 is still hiding among the
        # uncounted Ritz values: the largest of them is resolved well enough
        # (its residual below half its distance from 1) that the eigenvalue it
        # approximates cannot be 1. Without it, "nothing counted, twice" would
        # stop the iteration before anything had converged.
        rest = ~near_one
        separated = True
        if np.any(rest):
            j = int(np.flatnonzero(rest)[0])
            separated = bool(residuals[j] <= 0.5 * (1.0 - abs(values[j])))
        if converged and separated and count == previous:
            # Eigenvalue 1 is real, so its eigenspace has a real basis; a
            # complex Ritz vector of it spans it through its two parts, and
            # the rank of those parts (never a truncation) is the basis.
            parts = np.column_stack([np.real(ritz[:, near_one]), np.imag(ritz[:, near_one])])
            u, sv, _ = np.linalg.svd(parts, full_matrices=False)
            basis = u[:, sv > _QUOTIENT_RTOL * (float(sv[0]) if sv.size else 0.0)] if count else parts[:, :0]
            rest = np.abs(values[~near_one])
            gap = 1.0 - float(rest.max()) if rest.size else 1.0
            stop = float(residuals[near_one].max()) if count else 0.0
            return basis, gap, stop
        previous = count
    raise ValueError(
        f"NearKernel: the eigenvalue-1 space did not converge in {_MAX_ITERATIONS} applications of the averaged "
        f"splitting map (block {block}); the splitting's next eigenvalue is too close to 1 for this budget."
    )


def _largest_singular_value(action, transpose, n: int) -> float:
    r""":math:`\sigma_{\max}` by power iteration on :math:`A^{\mathsf T}A` (fast at the top of the spectrum)."""
    x = np.random.default_rng(_SEED + 1).standard_normal(n)
    x /= np.linalg.norm(x)
    estimate = 0.0
    for _ in range(500):
        y = transpose(action(x))
        norm = float(np.linalg.norm(y))
        if norm == 0.0:
            return 0.0
        x = y / norm
        if abs(norm - estimate) <= 1e-10 * norm:
            return float(np.sqrt(norm))
        estimate = norm
    return float(np.sqrt(estimate))
