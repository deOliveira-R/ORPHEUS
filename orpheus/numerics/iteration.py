r"""Iteration primitives for the operator-algebra :math:`(L - S - F)`.

The neutron transport equation, in its operator-algebra form, is

.. math::

    (L - S - F)\,\psi = q_{\rm ext}
    \qquad\text{(fixed source)}

.. math::

    (L - S)\,\psi = \tfrac{1}{k}\,F\,\psi
    \qquad\text{(eigenvalue)}

where :math:`L` is the streaming-collision operator, :math:`S` is
the scattering source operator, and :math:`F` is the fission source
operator (Lewis & Miller §6.4 frame this decomposition; Trefethen &
Bau 1997 §3.2 give the matrix-free Krylov view).

This module installs **two stand-alone iteration primitives** that
consume the Wave A :class:`~orpheus.numerics.operator.LinearOperator`
Protocol and operate on the operator triple :math:`(L, S, F)`:

* :class:`SourceIteration` solves the within-group fixed-source
  problem :math:`(L - S - F)\,\psi = q_{\rm ext}` by a fixed-point
  iteration

  .. math::

      \psi_{n+1} \;=\; L^{-1}\,(S\,\psi_n + F\,\psi_n + q_{\rm ext}).

  The convergence rate is bounded by the spectral radius
  :math:`\rho(L^{-1}(S+F)) \le \max_{\rm cell}\,\Sigma_s/\Sigma_t`
  for an SN sweep.  Trefethen & Bau §3.2.

* :class:`KEigenvalue` solves the eigenvalue problem
  :math:`(L - S)\,\psi = F\,\psi/k` by classical power iteration on
  the outer :math:`k`-update, with :class:`SourceIteration` driving
  the inner :math:`(L - S)\,\psi = F\psi_{\rm old}/k_{\rm old}`
  fixed-source solve at each outer step.  Convergence is governed
  by the dominance ratio :math:`|k_1/k_0|`.

The primitives are deliberately kept **shape-agnostic**.  They make
no assumption about the rank or layout of :math:`\psi` — only that
the operator triple acts linearly on it and that
:func:`numpy.linalg.norm` returns a scalar that orders by relative
size.  The L0 synthetic tests use 4×4 dense matrices acting on
``(4,)`` flat vectors; the L1 SN gate uses
:class:`~orpheus.sn.operator.SNStreamingOperator`,
:class:`~orpheus.sn.scattering.ScatteringOperator`, and
:class:`~orpheus.sn.fission.FissionOperator` acting on
:math:`(n_x, n_y, n_g)` scalar-flux arrays.

The ``inverter`` parameter — Wave E's load-bearing design choice
================================================================

Both primitives accept an optional ``inverter: Callable[[ndarray],
ndarray]`` that supplies :math:`L^{-1}`.  When ``None``, the primitive
uses :meth:`L.solve` directly.  When provided, the caller controls
how :math:`L^{-1}` is realised:

* ``inverter = None`` (default):  :math:`L^{-1}\,q = L.solve(q)` —
  for an SN sweep, this is the WDD asymmetric closure (Wave D
  Round 2 unified sweep), which has a closure-bias-driven self-
  consistent fixed point on curvilinear meshes (ERR-026).
* ``inverter = lambda q: gmres(as_scipy_linop(L), q, M=...)`` —
  Krylov-on-:meth:`apply` (the symmetric closure), with the sweep
  injected as a preconditioner :math:`M`.  This is the Wave E
  Round 2 reconciliation that closes ERR-026 for curvilinear SN.

By making :math:`L^{-1}` a caller-supplied hook, the iteration
primitives do not need to be re-implemented when the inversion
strategy changes.  The same :class:`SourceIteration` runs in the
synthetic L0 case (where ``L`` is a plain dense matrix and
``inverter`` defaults to a direct solve), in the L1 SN case (where
``inverter`` defaults to the WDD sweep), and in the future
Krylov-on-:meth:`apply` SN case (where ``inverter`` is supplied
explicitly by the caller).

Forward references
==================

* :func:`orpheus.numerics.eigenvalue.power_iteration` — the legacy
  primitive that this module supersedes for new SN code.  Stays
  alive (with a :class:`DeprecationWarning` on import) through the
  cross-solver migration sequence (CP, diffusion, MoC,
  homogeneous).  See :ref:`cross-solver-migration-sequence`.
* :class:`~orpheus.sn.solver.SNSolver` — Wave E Round 2 wires this
  class to consume :class:`SourceIteration` / :class:`KEigenvalue`,
  retiring the per-method delegators (Issue #163's downstream
  consumer).
"""

from __future__ import annotations

from typing import Callable, Protocol

import numpy as np

from .operator import (
    CAP_APPLY,
    CAP_SOLVE,
    LinearOperator,
    MissingCapability,
    ZeroOperator,
)


__all__ = [
    "SourceIteration",
    "KEigenvalue",
]


# Type alias for the L^{-1} hook: takes a right-hand-side, returns L^{-1}·rhs.
Inverter = Callable[[np.ndarray], np.ndarray]
KeffEstimator = Callable[
    [LinearOperator, LinearOperator, LinearOperator, np.ndarray], float,
]


def _has(op: object, cap: str) -> bool:
    """Return ``True`` iff ``op`` advertises the capability ``cap``."""
    caps = getattr(op, "capabilities", None)
    return bool(caps) and cap in caps


def _default_keff_estimator(
    L: LinearOperator,
    S: LinearOperator,
    F: LinearOperator,
    psi: np.ndarray,
) -> float:
    r"""Generic Rayleigh-quotient :math:`k`-estimator.

    Computes

    .. math::

        k \;=\; \frac{\sum (F\,\psi)}{\sum (L\,\psi) - \sum (S\,\psi)}

    which is the operator-form analogue of the production/absorption
    ratio for the homogeneous infinite-medium balance.  Each term is
    summed over every entry of the array; the ratio is dimensionless
    so any consistent reduction works.

    For SN problems with explicit volume weighting in the production /
    loss balance (see :meth:`orpheus.sn.solver.SNSolver.compute_keff`)
    the caller should supply a custom :pydata:`keff_estimator` that
    matches the existing math; otherwise the generic ratio above is
    correct for synthetic L0 cases and for any operator triple where
    the action already carries the volume measure.
    """
    num = np.sum(F.apply(psi))
    den = np.sum(L.apply(psi)) - np.sum(S.apply(psi))
    return float(num / den)


# ───────────────────────────────────────────────────────────────────────
# SourceIteration
# ───────────────────────────────────────────────────────────────────────


class SourceIteration:
    r"""Fixed-point iteration for :math:`(L - S - F)\,\psi = q_{\rm ext}`.

    Solves the within-group fixed-source equation in its operator-
    algebra form.  Each iteration applies :math:`L^{-1}` to a
    right-hand-side built from the current iterate plus the external
    source:

    .. math::

        \psi_{n+1} \;=\; L^{-1}\,(S\,\psi_n + F\,\psi_n + q_{\rm ext}).

    The :math:`L^{-1}` action is supplied by the ``inverter``
    parameter.  By default it routes through :meth:`L.solve`; the
    caller may override with a Krylov method on :meth:`L.apply` to
    decouple from a closure-biased :meth:`L.solve` (see ERR-026).

    Convergence test (matching :class:`SNSolver` for bit-identical
    Round 2 wiring):

    .. math::

        {\rm res}_n \;=\; \frac{\|\psi_n - \psi_{n-1}\|_2}
                                {\max(\|\psi_n\|_2,\,10^{-30})}

    and the iteration breaks when :math:`{\rm res}_n < {\rm tol}`.

    Parameters
    ----------
    L : LinearOperator
        Streaming-collision operator.  Must advertise
        :pydata:`CAP_APPLY`.  Must advertise :pydata:`CAP_SOLVE` *or*
        the caller must supply ``inverter`` — otherwise the
        iteration cannot evaluate :math:`L^{-1}`.
    S : LinearOperator
        Scattering source operator.  Must advertise
        :pydata:`CAP_APPLY`.  Pass
        :class:`~orpheus.numerics.operator.ZeroOperator` for the
        scattering-free case (the operator advertises ``apply`` so
        the capability check passes; its action is identically zero).
    F : LinearOperator
        Fission source operator.  Must advertise
        :pydata:`CAP_APPLY`.  Pass
        :class:`~orpheus.numerics.operator.ZeroOperator` for the
        fission-free case.
    inverter : callable or None, optional
        Function ``rhs -> L^{-1}·rhs``.  When ``None`` (default), the
        primitive uses ``L.solve`` directly — which requires ``L`` to
        advertise :pydata:`CAP_SOLVE`.  When provided, ``inverter``
        shadows ``L.solve`` and ``L`` is no longer required to
        advertise :pydata:`CAP_SOLVE`.  The Wave E Round 2 SN-on-Krylov
        path uses this hook to inject a ``gmres(as_scipy_linop(L),
        rhs, M=sweep_preconditioner)`` Krylov solve.
    max_iter : int, optional
        Maximum fixed-point iterations.  Default ``1000``.
    tol : float, optional
        Convergence tolerance on the relative residual norm.  Default
        ``1e-8``.

    Raises
    ------
    MissingCapability
        At construction time if ``L`` or ``S`` or ``F`` lacks
        :pydata:`CAP_APPLY`, or if ``L`` lacks :pydata:`CAP_SOLVE`
        and no ``inverter`` was supplied.

    Notes
    -----
    The primitive is shape-agnostic: ``q_ext`` is whatever shape the
    operators consume.  The convergence test uses
    :func:`numpy.linalg.norm` on the flattened arrays.  Both the L0
    synthetic case (flat ``(N,)`` vector) and the L1 SN case
    (structured ``(nx, ny, ng)`` array) are handled by the same
    :func:`np.linalg.norm` call by virtue of numpy's behaviour on
    higher-rank arrays (it returns the Frobenius norm).
    """

    def __init__(
        self,
        L: LinearOperator,
        S: LinearOperator,
        F: LinearOperator,
        *,
        inverter: Inverter | None = None,
        max_iter: int = 1000,
        tol: float = 1e-8,
    ) -> None:
        # Capability checks at construction so a downstream caller
        # never sees a stub failure mid-iteration (Wave A philosophy).
        if not _has(L, CAP_APPLY):
            raise MissingCapability(
                f"SourceIteration requires {CAP_APPLY!r} on L; "
                f"{type(L).__name__} advertises "
                f"{getattr(L, 'capabilities', frozenset())}."
            )
        if not _has(S, CAP_APPLY):
            raise MissingCapability(
                f"SourceIteration requires {CAP_APPLY!r} on S; "
                f"{type(S).__name__} advertises "
                f"{getattr(S, 'capabilities', frozenset())}."
            )
        if not _has(F, CAP_APPLY):
            raise MissingCapability(
                f"SourceIteration requires {CAP_APPLY!r} on F; "
                f"{type(F).__name__} advertises "
                f"{getattr(F, 'capabilities', frozenset())}."
            )
        if inverter is None and not _has(L, CAP_SOLVE):
            raise MissingCapability(
                f"SourceIteration requires either {CAP_SOLVE!r} on L "
                f"or a caller-supplied inverter; {type(L).__name__} "
                f"advertises "
                f"{getattr(L, 'capabilities', frozenset())} and no "
                f"inverter was provided."
            )

        self.L = L
        self.S = S
        self.F = F
        self.max_iter = int(max_iter)
        self.tol = float(tol)

        # Pin the inverter at construction.  Default: route through
        # L.solve.  Capability check above guarantees L.solve exists
        # when inverter is None.
        if inverter is None:
            self._inverter: Inverter = lambda q: self.L.solve(q)  # type: ignore[attr-defined]
        else:
            self._inverter = inverter

    def solve(
        self,
        q_ext: np.ndarray,
        initial_guess: np.ndarray | None = None,
    ) -> tuple[np.ndarray, list[float]]:
        r"""Run fixed-point iteration to convergence.

        Parameters
        ----------
        q_ext : np.ndarray
            External source.  Shape determined by what the operator
            triple consumes — for SN this is ``(nx, ny, ng)``; for the
            L0 synthetic case it is a flat ``(N,)`` vector.
        initial_guess : np.ndarray or None, optional
            Initial iterate.  When ``None`` (default), the iteration
            starts from :func:`np.zeros_like` of ``q_ext``.

        Returns
        -------
        psi : np.ndarray
            Converged iterate (or the final iterate if ``max_iter``
            was hit before tolerance was reached).
        residual_history : list[float]
            Relative residual at every iteration.  Useful for
            plotting convergence and for diagnosing stalled
            iterations.
        """
        psi = (
            np.zeros_like(q_ext)
            if initial_guess is None
            else np.asarray(initial_guess).copy()
        )
        residual_history: list[float] = []

        for _ in range(self.max_iter):
            psi_prev = psi

            # Build the RHS of the fixed-point step.  All three
            # operators are LinearOperators; their .apply contracts
            # are the only thing this loop touches.
            rhs = self.F.apply(psi) + self.S.apply(psi) + q_ext

            # Apply L^{-1} via the caller-supplied inverter (default
            # routes through L.solve).
            psi = self._inverter(rhs)

            # SNSolver-compatible convergence test (bit-identical loop
            # structure for Round 2 wiring): relative L2 residual.
            norm = float(np.linalg.norm(psi))
            if norm > 0.0:
                res = float(np.linalg.norm(psi - psi_prev) / max(norm, 1e-30))
            else:
                res = float(np.linalg.norm(psi - psi_prev))
            residual_history.append(res)

            if res < self.tol:
                break

        return psi, residual_history


# ───────────────────────────────────────────────────────────────────────
# KEigenvalue
# ───────────────────────────────────────────────────────────────────────


class KEigenvalue:
    r"""Power iteration for the eigenvalue :math:`(L - S)\,\psi = F\psi/k`.

    Outer loop: classical power iteration on :math:`k`.  Each outer
    step builds the fission source :math:`q_n = F\,\psi_n/k_n`, then
    drives an inner :class:`SourceIteration` with operator triple
    :math:`(L, S, 0)` and external source :math:`q_n`:

    .. math::

        \psi_{n+1} \;=\; (L - S)^{-1}\,F\,\psi_n / k_n

    .. math::

        k_{n+1} \;=\; {\rm keff\_estimator}(L, F, \psi_{n+1})

    Convergence test:  both :math:`|k_{n+1} - k_n| < {\rm keff\_tol}`
    AND the relative flux residual :math:`\|\psi_{n+1} -
    \psi_n\|_2 / \|\psi_{n+1}\|_2 < {\rm flux\_tol}`.

    The convergence rate is governed by the dominance ratio
    :math:`|k_1/k_0|` (Trefethen & Bau §27 power iteration analysis).

    Parameters
    ----------
    L, S, F : LinearOperator
        Operator triple, same constraints as
        :class:`SourceIteration`.  ``F`` must advertise
        :pydata:`CAP_APPLY` (no degenerate zero-fission case for an
        eigenvalue solve — without fission the spectrum is empty).
    inverter : callable or None, optional
        :math:`L^{-1}` hook, propagated to the inner
        :class:`SourceIteration`.  See :class:`SourceIteration` for
        the design rationale.
    max_outer : int, optional
        Maximum outer (power) iterations.  Default ``500``.
    keff_tol, flux_tol : float, optional
        Outer convergence tolerances.  Defaults match
        :class:`~orpheus.sn.solver.SNSolver` (``1e-7`` / ``1e-6``).
    max_inner : int, optional
        Inner :class:`SourceIteration` ``max_iter`` budget.  Default
        ``1000``.
    inner_tol : float, optional
        Inner :class:`SourceIteration` ``tol``.  Default ``1e-8``.
    eigenvalue_method : str, optional
        Forward hook for FEAST-style contour-integral methods (Issue
        #163 acceptance criterion).  Currently only ``"power"`` is
        implemented; other values raise :class:`NotImplementedError`
        at construction time.
    keff_estimator : callable or None, optional
        ``(L, S, F, psi) -> keff`` function used to update the outer
        eigenvalue.  When ``None`` (default), uses the generic
        Rayleigh-quotient :func:`_default_keff_estimator`.  SN
        consumers that need explicit volume weighting (matching
        :meth:`SNSolver.compute_keff`) should supply a custom
        ``keff_estimator``; the SN-aware estimator is the load-
        bearing way Round 2 preserves bit-identity with the legacy
        :func:`power_iteration` path.

    Raises
    ------
    MissingCapability
        Same conditions as :class:`SourceIteration`.
    NotImplementedError
        If ``eigenvalue_method`` is not ``"power"`` — the FEAST hook
        is reserved for a future wave.
    """

    def __init__(
        self,
        L: LinearOperator,
        S: LinearOperator,
        F: LinearOperator,
        *,
        inverter: Inverter | None = None,
        max_outer: int = 500,
        keff_tol: float = 1e-7,
        flux_tol: float = 1e-6,
        max_inner: int = 1000,
        inner_tol: float = 1e-8,
        eigenvalue_method: str = "power",
        keff_estimator: KeffEstimator | None = None,
    ) -> None:
        if eigenvalue_method != "power":
            raise NotImplementedError(
                f"KEigenvalue currently only supports "
                f"eigenvalue_method='power'; got "
                f"{eigenvalue_method!r}.  FEAST-style contour-integral "
                f"methods are a forward hook reserved for a future wave."
            )

        # Defer capability checks to the inner SourceIteration: the
        # same constraints apply, and we want one source of truth for
        # the messages.  The inner constructor will raise
        # MissingCapability for any deficient operand.
        # (Note: we still do the F-apply check above implicitly via
        # the inner SourceIteration receiving F as fission term.)

        self.L = L
        self.S = S
        self.F = F
        self.inverter = inverter
        self.max_outer = int(max_outer)
        self.keff_tol = float(keff_tol)
        self.flux_tol = float(flux_tol)
        self.max_inner = int(max_inner)
        self.inner_tol = float(inner_tol)
        self.eigenvalue_method = eigenvalue_method
        self._keff_estimator = (
            keff_estimator if keff_estimator is not None
            else _default_keff_estimator
        )

        # Validate operator capabilities up front by trial-construction
        # of the inner SourceIteration shell.  This catches any L/S/F
        # apply-capability violations and the L.solve / inverter
        # requirement at construction time, NEVER mid-iteration.
        SourceIteration(L, S, F, inverter=inverter, max_iter=1, tol=1.0)

    def solve(
        self,
        initial_guess: np.ndarray | None = None,
    ) -> tuple[float, list[float], np.ndarray]:
        r"""Run power iteration to convergence.

        Parameters
        ----------
        initial_guess : np.ndarray or None, optional
            Initial flux guess.  When ``None`` (default), the
            iteration starts from
            :func:`np.ones_like` of ``F.apply(<probe>)``.  When the
            caller cannot easily produce a probe for shape inference,
            they should supply ``initial_guess`` explicitly — this is
            the recommended path.

        Returns
        -------
        keff : float
            Converged dominant eigenvalue.
        keff_history : list[float]
            Eigenvalue at every outer iteration.
        psi : np.ndarray
            Converged fundamental-mode iterate (arbitrary
            normalisation).
        """
        if initial_guess is None:
            raise ValueError(
                "KEigenvalue.solve requires initial_guess for shape "
                "inference; the operator triple is not constrained to "
                "expose its action shape.  Pass np.ones(...) of the "
                "appropriate shape, or use the SNSolver wrapper that "
                "already builds the initial guess."
            )

        psi = np.asarray(initial_guess).copy()
        keff = 1.0
        keff_history: list[float] = []

        # Build a single inner SourceIteration shell.  Its operator
        # triple is (L, S, ZeroOperator) — the fission contribution at
        # the inner level is zero because F·psi/k is the *external
        # source* that drives the inner solve, NOT a within-group
        # fixed-point term.
        zero_F = ZeroOperator()
        inner = SourceIteration(
            self.L,
            self.S,
            zero_F,
            inverter=self.inverter,
            max_iter=self.max_inner,
            tol=self.inner_tol,
        )

        for n in range(1, self.max_outer + 1):
            psi_old = psi
            keff_old = keff

            # Outer fission source: q_ext = F·ψ / k  (the eigenvalue
            # division stays at the outer level — F.apply returns
            # F·ψ alone, the operator-algebra discipline).
            q_fission = self.F.apply(psi) / keff

            # Inner solve: (L - S)·ψ_new = q_fission.  The inner
            # SourceIteration warms up from psi_old to amortise the
            # inner cost across outer iterations (the same pattern
            # SNSolver._solve_source_iteration uses).
            psi, _inner_residuals = inner.solve(
                q_fission, initial_guess=psi_old,
            )

            # Outer keff update.
            keff = self._keff_estimator(self.L, self.S, self.F, psi)
            keff_history.append(keff)

            # Convergence test mirrors SNSolver.converged: at least 3
            # outer iterations before convergence is declared, then
            # both keff and flux must satisfy their tolerances.
            if n <= 2:
                continue
            dk = abs(keff - keff_old)
            norm_psi = float(np.linalg.norm(psi))
            dphi = (
                float(np.linalg.norm(psi - psi_old) / max(norm_psi, 1e-30))
                if norm_psi > 0.0
                else float(np.linalg.norm(psi - psi_old))
            )
            if dk < self.keff_tol and dphi < self.flux_tol:
                break

        return keff, keff_history, psi
