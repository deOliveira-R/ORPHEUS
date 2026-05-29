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

L⁻¹ comes from L, preconditioner is a different concept
========================================================

R-1 Step B (2026-05-19) — the legacy ``inverter`` Callable hook
retires.  The two primitives now read the inverse action directly off
the operator:

* :class:`SourceIteration` calls ``L.solve(rhs)`` directly each step.
  The caller controls the inverse step by passing an ``L`` whose
  ``.solve`` realises the desired action — typically a composite
  :class:`~orpheus.sn.operator.InvertibleOperator` whose ``.solve`` IS
  the WDD sweep (R-1 Step C).  ``L`` MUST advertise
  :pydata:`CAP_SOLVE`; if it cannot be inverted, the caller is asking
  the wrong primitive.

* :class:`KrylovAcceleration` takes an explicit ``preconditioner``
  Callable hook (the previous ``inverter`` name was a category
  mistake — it was a GMRES left preconditioner :math:`M \approx
  A^{-1}`, never the inverse step).  Default behaviour preserves the
  prior choice: if ``L`` advertises :pydata:`CAP_SOLVE`, use
  ``L.solve`` as the preconditioner; otherwise run unpreconditioned.

R-1 Step A integration — typed-flux Carlson seed threading
==========================================================

When :class:`SourceIteration` consumes a typed
:class:`~orpheus.sn.angular_flux.AngularFlux` rhs, the previous
iterate is threaded onto the rhs via
:meth:`~orpheus.sn.angular_flux.AngularFlux.stash` so the sweep can
read it back as ``rhs(1)``.  This is the load-bearing plumbing for
curvilinear sweeps where the previous iterate IS the Carlson seed.
Bare-ndarray rhs paths skip the stash (no protocol; the synthetic L0
tests have no seed dependency).

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

import inspect
import warnings
from typing import Callable, Protocol

import numpy as np
import scipy.sparse.linalg as spla

from .operator import (
    CAP_APPLY,
    CAP_SOLVE,
    LinearOperator,
    MissingCapability,
    ZeroOperator,
)


__all__ = [
    "SourceIteration",
    "KrylovAcceleration",
    "KEigenvalue",
]


# Type alias for the L^{-1} hook: takes a right-hand-side, returns L^{-1}·rhs.
Preconditioner = Callable[[np.ndarray], np.ndarray]
KeffEstimator = Callable[
    [LinearOperator, LinearOperator, LinearOperator, np.ndarray], float,
]
# Production-rate estimator: returns the volume-integrated fission
# production ``∫νΣ_f·φ·dV`` (scalar).  Used by :class:`KEigenvalue` to
# renormalise the iterate to unit production rate each outer step so
# subcritical eigenmodes do not underflow to denormalised FP — ERR-052.
ProductionEstimator = Callable[
    [LinearOperator, LinearOperator, LinearOperator, np.ndarray], float,
]


def _has(op: object, cap: str) -> bool:
    """Return ``True`` iff ``op`` advertises the capability ``cap``."""
    caps = getattr(op, "capabilities", None)
    return bool(caps) and cap in caps


# ───────────────────────────────────────────────────────────────────────
# Ravellable protocol — bridges typed flux containers to scipy's
# flat-vector requirement.  :class:`KrylovAcceleration` and
# :class:`SourceIteration` accept typed flux containers
# (:class:`~orpheus.transport.timed_full_field.TimedFullField`) via
# duck-typing on the pair of methods (``to_flat()`` instance method +
# class-level ``from_flat(flat, template)`` factory).
#
# Keeping the protocol duck-typed here (not an ABC import) preserves
# the deliberate decoupling of ``orpheus.numerics`` from
# ``orpheus.transport`` — the iteration primitives still know nothing
# about transport-specific shapes; they just consume any object that
# advertises the ravel pair.
#
# Bare ndarrays match neither protocol and fall through to the
# numpy reshape path in the helpers below.
# ───────────────────────────────────────────────────────────────────────


def _is_ravellable(x: object) -> bool:
    """Detect the ravellable protocol (template-based).

    Matches any object exposing ``to_flat()`` instance method +
    ``from_flat(flat, template)`` classmethod.  The canonical instance
    is :class:`~orpheus.transport.timed_full_field.TimedFullField`.
    """
    return (
        hasattr(x, "to_flat")
        and hasattr(type(x), "from_flat")
    )


def _ravel(x):
    """Ravel typed flux or bare ndarray to a 1-D ``float64`` ndarray."""
    if _is_ravellable(x):
        return np.asarray(x.to_flat(), dtype=float)
    return np.asarray(x, dtype=float).ravel()


def _unravel_like(template, flat: np.ndarray):
    """Reconstruct the typed flux (or reshape the bare ndarray) from ``flat``.

    Uses ``template`` only to recover the shape / mesh / factory —
    ``flat`` is the new numeric content.
    """
    if _is_ravellable(template):
        return type(template).from_flat(flat, template)
    return flat.reshape(template.shape)


def _zeros_like(template):
    """Zero typed flux or bare ndarray matching ``template``'s shape/mesh."""
    if _is_ravellable(template):
        flat_size = template.to_flat().size
        return type(template).from_flat(
            np.zeros(flat_size, dtype=float), template,
        )
    return np.zeros_like(template)


def _l2_norm(x) -> float:
    """L2 norm — ravels typed flux via the protocol before delegating to numpy."""
    if _is_ravellable(x):
        return float(np.linalg.norm(x.to_flat()))
    return float(np.linalg.norm(np.asarray(x)))


def _default_production_estimator(
    L: LinearOperator,
    S: LinearOperator,
    F: LinearOperator,
    psi: np.ndarray,
) -> float:
    r"""Generic production-rate estimator: :math:`P(\psi) = \sum (F\,\psi)`.

    The :math:`F` operator already carries any volume weights its
    domain advertises; the unweighted sum over array entries is the
    discrete analog of :math:`\int \nu\Sigma_f\,\phi\,dV` when the
    operator's action absorbs the measure (as ORPHEUS's typed
    operators do).

    SN consumers that need explicit volume weighting (matching
    :meth:`orpheus.sn.solver.SNSolver.compute_production_rate`) should
    supply a custom :pydata:`production_estimator`.
    """
    return float(np.sum(F.apply(psi)))


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

    The :math:`L^{-1}` action is read off ``L.solve`` directly.  The
    caller controls the inverse step by constructing ``L`` with the
    desired ``.solve`` behaviour — for SN within-group inner solves
    this is typically
    :class:`~orpheus.sn.operator.InvertibleOperator` whose ``.solve``
    IS the WDD sweep.

    Convergence test (matching :class:`SNSolver` for bit-identical
    Round 2 wiring):

    .. math::

        {\rm res}_n \;=\; \frac{\|\psi_n - \psi_{n-1}\|_2}
                                {\max(\|\psi_n\|_2,\,10^{-30})}

    and the iteration breaks when :math:`{\rm res}_n < {\rm tol}`.

    Parameters
    ----------
    L : LinearOperator
        Streaming-collision operator.  Must advertise BOTH
        :pydata:`CAP_APPLY` and :pydata:`CAP_SOLVE` — the iteration
        step is :math:`\psi_{n+1} = L^{-1}(S\psi_n + F\psi_n +
        q_{\rm ext})`, so ``L.solve`` is non-negotiable.
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
    max_iter : int, optional
        Maximum fixed-point iterations.  Default ``1000``.
    tol : float, optional
        Convergence tolerance on the relative residual norm.  Default
        ``1e-8``.

    Raises
    ------
    MissingCapability
        At construction time if ``L`` or ``S`` or ``F`` lacks
        :pydata:`CAP_APPLY`, or if ``L`` lacks :pydata:`CAP_SOLVE`.

    Notes
    -----
    The primitive is shape-agnostic: ``q_ext`` is whatever shape the
    operators consume.  The convergence test uses
    :func:`numpy.linalg.norm` on the flattened arrays.  Both the L0
    synthetic case (flat ``(N,)`` vector) and the L1 SN case
    (structured :class:`~orpheus.sn.angular_flux.AngularFlux`) are
    handled by the same :func:`_l2_norm` call routed through the
    ravellable protocol.

    R-1 Step A/B (2026-05-19) — when ``q_ext`` is a typed
    :class:`~orpheus.sn.angular_flux.AngularFlux`, the previous
    iterate ``psi`` is threaded onto the rhs as the lag-1 frame so
    ``L.solve`` can read it via ``rhs(1)`` for Carlson-seed plumbing
    on curvilinear sweeps.  See :func:`_attach_previous`.
    """

    def __init__(
        self,
        L: LinearOperator,
        S: LinearOperator,
        F: LinearOperator,
        *,
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
        if not _has(L, CAP_SOLVE):
            raise MissingCapability(
                f"SourceIteration requires {CAP_SOLVE!r} on L — the "
                f"iteration step is psi_(n+1) = L.solve(...); "
                f"{type(L).__name__} advertises "
                f"{getattr(L, 'capabilities', frozenset())}.  Pass an "
                f"L whose .solve realises the desired inverse action "
                f"(typically an InvertibleOperator)."
            )

        self.L = L
        self.S = S
        self.F = F
        self.max_iter = int(max_iter)
        self.tol = float(tol)
        # Detect once whether ``L.solve`` accepts ``initial_guess`` —
        # InvertibleOperator does (Phase 1.2, post Carlson-seed
        # unification); MatrixOperator and other generic LinearOperators
        # do not.  Avoids per-iteration introspection.
        try:
            self._solve_accepts_seed = (
                "initial_guess" in inspect.signature(self.L.solve).parameters
            )
        except (TypeError, ValueError):
            self._solve_accepts_seed = False

    def _solve_with_seed(self, rhs, psi_prev):
        """Dispatch ``L.solve`` with or without ``initial_guess`` per L's signature."""
        if self._solve_accepts_seed:
            return self.L.solve(rhs, initial_guess=psi_prev)
        return self.L.solve(rhs)

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
            triple consumes — for SN this is ``(ng, nx, ny)``
            (principled storage; see :ref:`theory-sn-index-convention`);
            for the L0 synthetic case it is a flat ``(N,)`` vector.
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
        # R-1 Step 4a — ravellable protocol: typed flux containers
        # (:class:`AngularFlux`) and bare ndarrays both work.  ``psi``
        # carries the same type as ``q_ext``; arithmetic, norm, and
        # zeros routing through the protocol helpers.
        if initial_guess is None:
            psi = _zeros_like(q_ext)
        elif _is_ravellable(initial_guess):
            psi = initial_guess  # typed: trust frozen-arithmetic contract
        else:
            psi = np.asarray(initial_guess).copy()
        residual_history: list[float] = []

        for _ in range(self.max_iter):
            psi_prev = psi

            # Build the RHS of the fixed-point step.  All three
            # operators are LinearOperators; their .apply contracts
            # are the only thing this loop touches.
            rhs = self.F.apply(psi) + self.S.apply(psi) + q_ext

            # Apply L^{-1} directly.  Phase-1.2 — the curvilinear
            # Carlson coupled-pole seed and the reflective-BC partner-
            # flux trace travel through ``initial_guess`` explicitly
            # (the M-M closure reads the level's ψ at μ = −1 from
            # this argument).  For bare-ndarray rhs the kwarg is
            # ignored downstream — the synthetic L0 tests have no
            # seed dependency.
            psi = self._solve_with_seed(rhs, psi_prev)

            # SNSolver-compatible convergence test (bit-identical loop
            # structure for Round 2 wiring): relative L2 residual via
            # the ravellable protocol — typed flux uses the flat
            # representation including boundary face state.
            norm = _l2_norm(psi)
            if norm > 0.0:
                res = _l2_norm(psi - psi_prev) / max(norm, 1e-30)
            else:
                res = _l2_norm(psi - psi_prev)
            residual_history.append(res)

            if res < self.tol:
                break

        return psi, residual_history


# ───────────────────────────────────────────────────────────────────────
# KrylovAcceleration
# ───────────────────────────────────────────────────────────────────────


class KrylovAcceleration:
    r"""GMRES on :math:`(L - S - F)\,\psi = q_{\rm ext}` — sibling of
    :class:`SourceIteration` for the same algebra.

    Both primitives solve the same fixed-source equation; they differ
    in algorithm.  :class:`SourceIteration` lags :math:`(S + F)\,\psi`
    as the right-hand side and inverts :math:`L` at every step
    (geometric convergence at rate :math:`\rho(L^{-1}(S+F)) \le
    \max\Sigma_s/\Sigma_t`).  :class:`KrylovAcceleration` builds the
    composed matvec :math:`(L - S - F)\cdot` as a single linear operator
    and solves it with GMRES, optionally preconditioned by :math:`L^{-1}`
    (the sweep).  When the scattering ratio :math:`c = \Sigma_s/\Sigma_t`
    approaches 1, GMRES converges in :math:`\mathcal{O}(\sqrt{\kappa})`
    matvecs vs source iteration's :math:`\mathcal{O}(1/(1-c))` — the
    standard transport-Krylov win documented in Adams & Larsen 2002
    (the SAILOR / preconditioned-Krylov framework).

    Algebra-of-record:

    .. math::

        (L - S - F)\,\psi \;=\; q_{\rm ext}.

    The composed matvec is realised as ``L.apply(psi) - S.apply(psi) -
    F.apply(psi)`` per call — no intermediate :class:`OperatorSum`
    allocation.  The right-hand side is whatever shape the operator
    triple consumes; scipy GMRES requires a flat 1-D view internally,
    so the primitive ravels at the boundary and reshapes the solution
    back to ``q_ext.shape`` on return.

    The ``preconditioner`` parameter (R-1 Step B rename)
    =====================================================

    The GMRES PRECONDITIONER :math:`M \approx A^{-1}` where :math:`A =
    L - S - F`.  The natural choice for transport problems is :math:`M
    = L^{-1}` (the sweep) — this is the "transport-corrected"
    preconditioner from Adams & Larsen 2002 §III.  When :math:`c` is
    small, the sweep is an excellent preconditioner; when :math:`c` is
    near unity, the sweep is diffusion-like and GMRES needs more
    iterations.

    Previously named ``inverter`` — but that was a category mistake:
    in :class:`SourceIteration` the ``inverter`` realised the
    iteration's INVERSE step :math:`L^{-1}`, whereas here ``M`` is the
    GMRES left preconditioner (an approximation to the FULL inverse
    :math:`A^{-1}`).  The rename surfaces the distinction and
    decouples the two surfaces (R-1 Step B, 2026-05-19).

    * ``preconditioner = None`` (default): if ``L`` advertises
      :pydata:`CAP_SOLVE`, use ``L.solve`` as the preconditioner;
      otherwise, no preconditioner (identity ``M = I``).
    * ``preconditioner = lambda q: sweep_preconditioner(q)``:
      caller-supplied preconditioner.  Typically wraps a sweep that
      consumes the same packed/structured layout the operators
      consume.

    Parameters
    ----------
    L, S, F : LinearOperator
        Operator triple.  Must each advertise :pydata:`CAP_APPLY`.
        Pass :class:`ZeroOperator` for absent terms (e.g. ``F = Zero``
        for within-group fixed-source: :class:`KEigenvalue` builds the
        fission source as an EXTERNAL :math:`q_{\rm ext}` and zeroes
        the within-group ``F``).
    preconditioner : callable or None, optional
        GMRES left preconditioner.  See above.  When ``None`` and
        ``L`` has no :pydata:`CAP_SOLVE`, runs GMRES without
        preconditioner.
    max_iter : int, optional
        Maximum GMRES iterations (``maxiter`` in scipy).  Default
        ``1000``.
    tol : float, optional
        GMRES relative residual tolerance (``rtol`` in scipy).
        Default ``1e-8``.
    restart : int, optional
        GMRES restart length.  Default ``50``.  Clamped to ``n`` at
        :meth:`solve` time.

    Raises
    ------
    MissingCapability
        At construction time if any of ``L``, ``S``, ``F`` lacks
        :pydata:`CAP_APPLY`.

    Notes
    -----
    The primitive is shape-agnostic at the operator-triple level — it
    only requires that the operators all consume and return arrays of
    the same shape as ``q_ext``.  Internally it ravels to 1-D for
    scipy's GMRES requirement and reshapes the solution to
    ``q_ext.shape`` on return.
    """

    def __init__(
        self,
        L: LinearOperator,
        S: LinearOperator,
        F: LinearOperator,
        *,
        preconditioner: Preconditioner | None = None,
        max_iter: int = 1000,
        tol: float = 1e-8,
        restart: int = 50,
    ) -> None:
        if not _has(L, CAP_APPLY):
            raise MissingCapability(
                f"KrylovAcceleration requires {CAP_APPLY!r} on L; "
                f"{type(L).__name__} advertises "
                f"{getattr(L, 'capabilities', frozenset())}."
            )
        if not _has(S, CAP_APPLY):
            raise MissingCapability(
                f"KrylovAcceleration requires {CAP_APPLY!r} on S; "
                f"{type(S).__name__} advertises "
                f"{getattr(S, 'capabilities', frozenset())}."
            )
        if not _has(F, CAP_APPLY):
            raise MissingCapability(
                f"KrylovAcceleration requires {CAP_APPLY!r} on F; "
                f"{type(F).__name__} advertises "
                f"{getattr(F, 'capabilities', frozenset())}."
            )

        self.L = L
        self.S = S
        self.F = F
        self.max_iter = int(max_iter)
        self.tol = float(tol)
        self.restart = int(restart)

        # Pin the preconditioner choice at construction.  If caller
        # supplied one, use it.  Otherwise, fall back to L.solve when
        # available; if not, run GMRES without preconditioner.
        if preconditioner is not None:
            self._preconditioner: Preconditioner | None = preconditioner
        elif _has(L, CAP_SOLVE):
            self._preconditioner = lambda q: self.L.solve(q)  # type: ignore[attr-defined]
        else:
            self._preconditioner = None

    def solve(
        self,
        q_ext: np.ndarray,
        initial_guess: np.ndarray | None = None,
    ) -> tuple[np.ndarray, list[float]]:
        r"""Run GMRES on :math:`(L - S - F)\,\psi = q_{\rm ext}` to convergence.

        Parameters
        ----------
        q_ext : np.ndarray
            External source.  Whatever shape the operator triple
            consumes — ravelled to 1-D internally for scipy GMRES,
            reshaped back to ``q_ext.shape`` on return.
        initial_guess : np.ndarray or None, optional
            GMRES initial iterate.  When ``None`` (default), starts
            from :func:`np.zeros_like` of ``q_ext``.

        Returns
        -------
        psi : np.ndarray
            Converged solution, shape ``q_ext.shape``.
        residual_history : list[float]
            Preconditioned residual norm at every GMRES inner iteration
            (scipy's ``callback_type='pr_norm'``).  Empty if GMRES
            returned in zero iterations.
        """
        # R-1 Step 4a — ravellable protocol: when ``q_ext`` is a typed
        # flux container (:class:`AngularFlux`), the ravel/unravel goes
        # through ``to_flat_with_traces`` / ``from_flat_with_traces`` so
        # the flat vector for scipy carries the FULL state including the
        # boundary face block.  Bare-ndarray inputs route through the
        # legacy reshape path unchanged.
        b = _ravel(q_ext)
        n = b.size

        def A_matvec(psi_flat: np.ndarray) -> np.ndarray:
            # Lift back to the typed (or shaped) carrier, compose
            # (L − S − F)·ψ, ravel.  Operator arithmetic propagates
            # via dunders to ``.boundary`` (AngularFlux) or just the
            # ndarray (bare).
            psi = _unravel_like(q_ext, psi_flat)
            out = self.L.apply(psi) - self.S.apply(psi) - self.F.apply(psi)
            return _ravel(out)

        A_scipy = spla.LinearOperator((n, n), matvec=A_matvec, dtype=float)

        if self._preconditioner is not None:
            precond_fn = self._preconditioner

            def M_matvec(q_flat: np.ndarray) -> np.ndarray:
                q = _unravel_like(q_ext, q_flat)
                out = precond_fn(q)
                return _ravel(out)

            M_scipy: spla.LinearOperator | None = spla.LinearOperator(
                (n, n), matvec=M_matvec, dtype=float,
            )
        else:
            M_scipy = None

        x0 = (
            _ravel(initial_guess)
            if initial_guess is not None
            else np.zeros_like(b)
        )

        residual_history: list[float] = []

        def callback(rk: object) -> None:
            # scipy GMRES with callback_type='pr_norm' passes the
            # preconditioned-residual norm (a scalar).  Older versions
            # may pass the residual vector — handle both defensively.
            r = np.asarray(rk)
            if r.ndim == 0:
                residual_history.append(float(r))
            else:
                residual_history.append(float(np.linalg.norm(r)))

        try:
            solution, info = spla.gmres(
                A_scipy, b, x0=x0, M=M_scipy,
                rtol=self.tol, atol=0.0,
                maxiter=self.max_iter,
                restart=min(self.restart, n),
                callback=callback,
                callback_type='pr_norm',
            )
        except TypeError:
            # Older scipy versions use ``tol`` instead of ``rtol`` and
            # may not honour ``callback_type``.  Drop both keywords.
            solution, info = spla.gmres(
                A_scipy, b, x0=x0, M=M_scipy,
                tol=self.tol,
                maxiter=self.max_iter,
                restart=min(self.restart, n),
                callback=callback,
            )

        # D-H.1e (2026-05-28) — surface GMRES non-convergence.  Pre-fix
        # the ``info`` flag was discarded; an unconverged ``solution``
        # would silently be consumed as if it were the true inverse.
        # Conjunction with the legacy ``restart=min(50, full_size)``
        # clamp at the caller produced the ERR-053 keff drift on
        # curvilinear meshes (the GMRES iteration was structurally
        # truncated, then the failure was hidden by this discard).
        # scipy convention: ``info > 0`` means "not converged within
        # ``maxiter``"; ``info < 0`` means illegal-input.  Both
        # surface as warnings — raising would break long-standing
        # callers that tolerate slow convergence and need the
        # best-effort iterate.  See ERR-053.
        if info != 0:
            warnings.warn(
                f"KrylovAcceleration.solve: scipy.sparse.linalg.gmres "
                f"returned info={info} (not converged within "
                f"maxiter={self.max_iter}; restart={min(self.restart, n)}; "
                f"rtol={self.tol}).  Returning best-effort iterate; "
                f"residual_history tail = "
                f"{residual_history[-3:] if residual_history else '[]'}.  "
                f"Tighten ``restart`` to ``n`` (full size) if the Krylov "
                f"subspace is being truncated; see ERR-053.",
                RuntimeWarning,
                stacklevel=2,
            )

        return _unravel_like(q_ext, solution), residual_history


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
        :class:`SourceIteration` — ``L`` MUST advertise BOTH
        :pydata:`CAP_APPLY` and :pydata:`CAP_SOLVE`; ``S`` and ``F``
        must advertise :pydata:`CAP_APPLY`.  ``F`` is non-trivial for
        an eigenvalue solve (no degenerate zero-fission case — without
        fission the spectrum is empty).
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
        max_outer: int = 500,
        keff_tol: float = 1e-7,
        flux_tol: float = 1e-6,
        max_inner: int = 1000,
        inner_tol: float = 1e-8,
        eigenvalue_method: str = "power",
        keff_estimator: KeffEstimator | None = None,
        production_estimator: ProductionEstimator | None = None,
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

        self.L = L
        self.S = S
        self.F = F
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
        self._production_estimator = (
            production_estimator if production_estimator is not None
            else _default_production_estimator
        )

        # Validate operator capabilities up front by trial-construction
        # of the inner SourceIteration shell.  This catches any L/S/F
        # apply-capability violations and the L.solve requirement at
        # construction time, NEVER mid-iteration.
        SourceIteration(L, S, F, max_iter=1, tol=1.0)

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

            # Renormalise to unit production rate so the iterate stays
            # at a physically natural O(1) magnitude across iterations
            # regardless of whether the operator is supercritical
            # (k>1, would grow) or subcritical (k<1, would decay to
            # denormalised FP within ~30-60 iterations and the keff
            # ratio would become 0/0 numerically meaningless —
            # ERR-052).  Production rate is scale-invariant in
            # ``keff`` so the converged eigenvalue is unchanged; the
            # converged ``ψ`` carries the canonical reactor-physics
            # output convention :math:`\\int \\nu\\Sigma_f\\,\\phi\\,dV = 1`,
            # which makes rescaling to absolute flux at a target power
            # a single multiplication by
            # :math:`P_{\\text{target}} / \\kappa`.
            production = self._production_estimator(self.L, self.S, self.F, psi)
            if production > 0.0:
                psi = psi / production

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
