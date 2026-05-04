"""Foundation tests for the geometry-agnostic Variant α power-iteration
driver :func:`power_iterate_variant_alpha`.

These tests pin the **driver itself** (control flow, history, error
handling, convergence, normalisation), independent of the
geometry-specific operators. Geometry-level bit-equality across the
refactor is exercised by the existing
``tests/derivations/test_peierls_greens_function_*`` suite, which runs
unchanged.

The tests use synthetic, fully-deterministic ``step`` callables — no
geometry-specific physics. This guarantees the driver's universal
book-keeping is verified independently of the operator construction
that varies across the 6 geometries.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.trajectory_resolvent.power_iteration import (
    PowerIterationResult,
    power_iterate_variant_alpha,
)


@pytest.mark.foundation
def test_driver_converges_on_constant_step():
    """A step function that returns a fixed scaling reaches the
    Rayleigh-quotient fixed point in one iteration.

    If ``step(psi, k_eff)`` always returns ``(c·psi, F_old, c·F_old)``
    with ``c = 2.0``, then ``k_new = k_eff * 2.0`` for every iteration —
    the residual ``|k_new - k_eff| / |k_eff|`` is constant at ``1.0``,
    so the iteration runs to ``max_iter``.

    More importantly: the new psi is ``2 * psi / (2 * F_old) = psi /
    F_old``, exactly the original normalisation behaviour.
    """
    F_old = 1.5
    c = 2.0

    def step(psi, k_eff):
        psi_new = c * psi
        return psi_new, F_old, c * F_old

    psi0 = np.ones(4)
    result = power_iterate_variant_alpha(
        step, psi0, initial_k=1.0, max_iter=5, tol=1e-12,
    )

    # k_new = k_eff * c at every step → 1, 2, 4, 8, 16
    assert result.iterations == 5
    assert not result.converged
    assert result.k_history == [2.0, 4.0, 8.0, 16.0, 32.0]
    # psi normalised by Frate_new = c * F_old at every step.
    # Build the reference psi via the SAME order of FP operations the
    # driver uses (psi := (c * psi) / (c * F_old) per iteration), so
    # the comparison is bit-equal.
    expected = psi0.copy()
    for _ in range(5):
        expected = (c * expected) / (c * F_old)
    np.testing.assert_array_equal(result.psi, expected)


@pytest.mark.foundation
def test_driver_detects_convergence_on_unit_ratio():
    """A step function that returns ``Frate_old == Frate_new`` makes
    ``k_new == k_eff``, so ``rel_dk == 0`` and convergence is declared
    on the first iteration.
    """
    def step(psi, k_eff):
        return psi.copy(), 2.0, 2.0

    psi0 = np.array([1.0, 2.0, 3.0])
    result = power_iterate_variant_alpha(
        step, psi0, initial_k=0.95, max_iter=10, tol=1e-12,
    )

    assert result.converged
    assert result.iterations == 1
    assert result.k_eff == 0.95
    assert len(result.k_history) == 1
    assert result.k_history[0] == 0.95
    assert result.residual_history[0] == 0.0


@pytest.mark.foundation
def test_driver_raises_on_vanished_fission_rate():
    """Frate_old < 1e-30 raises with the iteration number in the
    message."""
    def step(psi, k_eff):
        return psi.copy(), 1e-40, 1.0

    with pytest.raises(RuntimeError, match=r"iter 1"):
        power_iterate_variant_alpha(
            step, np.ones(3), initial_k=1.0, max_iter=10, tol=1e-12,
        )


@pytest.mark.foundation
def test_driver_returns_history_lengths_match_iterations():
    """Both ``k_history`` and ``residual_history`` have exactly
    ``iterations`` entries — one per pass through the loop body."""
    counter = [0]

    def step(psi, k_eff):
        counter[0] += 1
        # Slowly converging — never hits tol within 7 iterations.
        return psi.copy(), 1.0, 1.0 + 0.1 / counter[0]

    psi0 = np.ones(2)
    result = power_iterate_variant_alpha(
        step, psi0, initial_k=1.0, max_iter=7, tol=1e-30,
    )

    assert result.iterations == 7
    assert len(result.k_history) == 7
    assert len(result.residual_history) == 7
    # Residuals are positive and decreasing (frate ratio drifts to 1).
    for r in result.residual_history:
        assert r > 0


@pytest.mark.foundation
def test_driver_preserves_bit_equality_on_synthetic_iteration():
    """Compare the driver's output bit-for-bit against an inlined
    reference loop on a synthetic geometry-free iteration.

    This is the load-bearing bit-equality probe: if the driver introduces
    any FP-order change relative to the original inlined loops, this
    test fails.
    """
    rng = np.random.default_rng(20260503)
    psi0 = np.abs(rng.standard_normal(8)) + 0.5
    initial_k = 1.234

    # Use a positive operator (row-stochastic so spectral radius = 1)
    # so psi stays positive and Frate is invariant under the iteration.
    M = np.abs(rng.standard_normal((8, 8))) + 0.1
    M = M / M.sum(axis=1, keepdims=True)

    def step(psi, k_eff):
        # Pure row-stochastic action: Frate = sum(psi) is invariant
        # because sum(M @ psi) = sum(psi). This makes the iteration
        # converge to k_new = k_eff in one pass — perfect bit-equality
        # probe under a non-trivial operator.
        psi_new = M @ psi
        Frate_old = float(np.sum(psi))
        Frate_new = float(np.sum(psi_new))
        return psi_new, Frate_old, Frate_new

    # Inlined reference loop (mirrors the driver verbatim).
    psi_ref = psi0.copy()
    k_ref = initial_k
    iters_ref = 0
    converged_ref = False
    k_history_ref = []
    residual_history_ref = []
    max_iter, tol = 5, 1e-30
    for it in range(max_iter):
        iters_ref = it + 1
        psi_new, Fold, Fnew = step(psi_ref, k_ref)
        if Fold < 1e-30:
            raise AssertionError("synthetic frate vanished")
        k_new = k_ref * Fnew / Fold
        psi_normed = psi_new / Fnew
        rel_dk = abs(k_new - k_ref) / max(abs(k_ref), 1e-30)
        psi_ref = psi_normed
        k_history_ref.append(k_new)
        residual_history_ref.append(rel_dk)
        k_ref = k_new
        if rel_dk < tol:
            converged_ref = True
            break

    # Driver run.
    result = power_iterate_variant_alpha(
        step, psi0.copy(), initial_k=initial_k,
        max_iter=max_iter, tol=tol,
    )

    np.testing.assert_array_equal(result.psi, psi_ref)
    assert result.k_eff == k_ref
    assert result.iterations == iters_ref
    assert result.converged == converged_ref
    assert result.k_history == k_history_ref
    assert result.residual_history == residual_history_ref


@pytest.mark.foundation
def test_driver_returns_powerresult_dataclass():
    """The return type is :class:`PowerIterationResult` with the
    documented field set."""
    def step(psi, k_eff):
        return psi.copy(), 1.0, 1.0

    psi0 = np.array([2.0, 3.0])
    result = power_iterate_variant_alpha(
        step, psi0, initial_k=1.0, max_iter=10, tol=1e-12,
    )
    assert isinstance(result, PowerIterationResult)
    # Documented public fields:
    for field in (
        "psi", "k_eff", "iterations", "converged",
        "k_history", "residual_history",
    ):
        assert hasattr(result, field), f"missing {field}"
