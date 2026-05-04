r"""Geometry-agnostic Variant α power-iteration driver.

This module collapses the 11 byte-identical power-iteration outer loops
that previously lived inside the per-geometry solvers
(:mod:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function`,
:mod:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function_cylinder`,
:mod:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function_slab_asymmetric`,
:mod:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function_hollow_sphere`,
:mod:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function_annulus`)
into a single function.

The fixed-source iteration in
:func:`solve_greens_function_sphere_mr_fixed_source` is structurally
distinct (no Rayleigh quotient, convergence on scalar flux instead of
on :math:`k_{\rm eff}`) and is NOT collapsed by this driver — it is a
fixed-source iteration, not a power iteration.

Design — preservation of bit-equality
--------------------------------------

Every per-geometry power-iteration loop has the same control-flow
shape (Rayleigh-quotient update of :math:`k_{\rm eff}`, fission-rate
normalisation, relative-:math:`k` convergence test). The
geometry-specific work — computing the source profile, applying the
Variant α operator, evaluating the volume-integrated fission rate —
varies in array-axis structure and call signature.

To collapse the 11 loops without changing IEEE-754 evaluation order at
any call site, the driver takes a **single per-step callable** that
encapsulates the geometry-specific arithmetic in exactly the order it
ran pre-refactor:

.. code-block:: python

    def step(psi, k_eff) -> (psi_new, Frate_old, Frate_new):
        '''One Variant α iteration body — geometry-specific arithmetic.

        Computes the new angular flux from the current iterate, plus
        the OLD and NEW volume-integrated fission rates needed for the
        Rayleigh-quotient update. The driver multiplies by
        ``k_eff * Frate_new / Frate_old`` and normalises by Frate_new.
        '''
        ...
        return psi_new, Frate_old, Frate_new

The driver itself only carries the universal book-keeping
(:math:`k_{\rm new} = k_{\rm eff} \cdot \mathrm{Frate}_{\rm new} /
\mathrm{Frate}_{\rm old}`, normalisation, rel-:math:`k` test, history
collection). Because the per-step callable does exactly the same
arithmetic in the same order as the inlined original, IEEE-754
reproducibility guarantees bit-equality of every numerical output.

Diagnostics
-----------

The driver collects ``k_history`` and ``residual_history`` (one entry
per iteration). These were computed implicitly by every original loop
(``rel_dk`` is the loop-final residual; ``k_new`` is the loop-final
k_eff) but were thrown away. They are now retained for free in
:class:`PowerIterationResult` and are used by the test-architect's
cross-method regression net for convergence-quality diagnostics.

The driver does NOT change the existing per-geometry public API. Each
``solve_greens_function_*`` function still returns its existing
geometry-specific result dataclass; the diagnostics dataclass is
unwrapped at the call site to populate the existing fields. R2 will
refactor the return shapes separately.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


@dataclass(frozen=True)
class PowerIterationResult:
    """Result of a single Variant α power-iteration solve.

    Attributes
    ----------
    psi : ndarray
        Converged (or final) angular flux iterate. Shape is
        geometry-specific — the driver is shape-agnostic.
    k_eff : float
        Final eigenvalue estimate.
    iterations : int
        Number of outer iterations performed (1-indexed: same convention
        as the pre-refactor loops, where the loop body sets
        ``iterations = it + 1`` on every pass).
    converged : bool
        ``True`` iff the relative-:math:`k` residual fell below ``tol``
        at the final iteration. If ``max_iter`` is exhausted without
        meeting the tolerance, ``converged = False``.
    k_history : list[float]
        ``k_eff`` value AFTER each iteration's Rayleigh-quotient update.
        Length equals ``iterations``.
    residual_history : list[float]
        Relative-:math:`k` residual ``|Δk|/|k|`` after each iteration.
        Length equals ``iterations``.
    """

    psi: np.ndarray
    k_eff: float
    iterations: int
    converged: bool
    k_history: list[float]
    residual_history: list[float]


def power_iterate_variant_alpha(
    step: Callable[[np.ndarray, float], tuple[np.ndarray, float, float]],
    initial_psi: np.ndarray,
    initial_k: float,
    *,
    max_iter: int = 200,
    tol: float = 1e-10,
) -> PowerIterationResult:
    r"""Geometry-agnostic Variant α power-iteration driver.

    Implements the universal book-keeping shared by every per-geometry
    Variant α solver in :mod:`~orpheus.derivations.continuous.trajectory_resolvent`:
    Rayleigh-quotient update of :math:`k_{\rm eff}`, fission-rate
    normalisation, relative-:math:`k` convergence test, history
    collection.

    The geometry-specific arithmetic — computing the source profile,
    applying the operator, evaluating the volume-integrated fission
    rate — is encapsulated in the ``step`` callable supplied by the
    caller.

    The per-iteration sequence is:

    1. ``psi_new, Frate_old, Frate_new = step(psi, k_eff)``
    2. If ``Frate_old < 1e-30``: raise ``RuntimeError``.
    3. ``k_new = k_eff * Frate_new / Frate_old`` (Rayleigh quotient).
    4. ``psi = psi_new / Frate_new`` (normalise to keep magnitude
       :math:`O(1)`).
    5. ``rel_dk = |k_new - k_eff| / max(|k_eff|, 1e-30)``.
    6. Append ``k_new`` to ``k_history`` and ``rel_dk`` to
       ``residual_history``.
    7. ``k_eff = k_new``.
    8. If ``rel_dk < tol``: declare converged, break.

    The order of FP operations is identical to the inlined per-geometry
    loops pre-refactor. Bit-equality of ``k_eff`` and ``psi`` is
    preserved by IEEE-754 reproducibility.

    Parameters
    ----------
    step : callable
        ``(psi, k_eff) -> (psi_new, Frate_old, Frate_new)``. One
        Variant α iteration body. ``Frate_old`` is the volume-integrated
        fission rate from the CURRENT ``psi``; ``Frate_new`` is from
        ``psi_new``. The driver does NOT compute these — they are
        geometry-specific (volume-integral Jacobians differ between
        sphere/cylinder/slab) and are needed in the same order they
        appear in the pre-refactor loops to preserve bit-equality.
    initial_psi : ndarray
        Initial angular flux iterate. Geometry-specific shape.
    initial_k : float
        Initial :math:`k_{\rm eff}` estimate. Typically ``k_inf`` (or
        the largest eigenvalue of the homogenised infinite-medium
        transfer matrix in MG).
    max_iter : int, default 200
        Maximum number of outer iterations.
    tol : float, default 1e-10
        Relative-:math:`k` convergence tolerance. Iteration stops the
        first time ``|k_new - k_eff| / max(|k_eff|, 1e-30) < tol``.

    Returns
    -------
    PowerIterationResult
        Converged ``psi``, ``k_eff``, iteration count, convergence
        flag, and per-iteration histories.

    Raises
    ------
    RuntimeError
        If the volume-integrated fission rate vanishes below ``1e-30``
        at any iteration (signals a non-multiplying medium or a
        degenerate iteration). The error message names the iteration
        number (1-indexed) for diagnostic clarity.
    """
    psi = initial_psi
    k_eff = initial_k
    iterations = 0
    converged = False
    k_history: list[float] = []
    residual_history: list[float] = []

    for it in range(max_iter):
        iterations = it + 1

        psi_new, Frate_old, Frate_new = step(psi, k_eff)

        if Frate_old < 1e-30:
            raise RuntimeError(
                f"Fission rate vanished at iter {iterations}; "
                "non-multiplying medium."
            )

        k_new = k_eff * Frate_new / Frate_old
        psi_normed = psi_new / Frate_new

        rel_dk = abs(k_new - k_eff) / max(abs(k_eff), 1e-30)
        psi = psi_normed
        k_history.append(k_new)
        residual_history.append(rel_dk)
        k_eff = k_new

        if rel_dk < tol:
            converged = True
            break

    return PowerIterationResult(
        psi=psi,
        k_eff=float(k_eff),
        iterations=iterations,
        converged=converged,
        k_history=k_history,
        residual_history=residual_history,
    )
