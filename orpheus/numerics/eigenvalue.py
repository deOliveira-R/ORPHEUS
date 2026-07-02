r"""Generalized eigenvalue solvers for neutron transport and diffusion.

This module hosts the **canonical power-iteration algorithm** for the
deterministic criticality problem, posed as the **generalized eigenproblem**

.. math::

    A_{\rm loss}\,\psi \;=\; \lambda\,M\,\psi ,

whose power-method realization is the dominant eigenpair of the **resolvent**
:math:`A_{\rm loss}^{-1} M`.  By Krein--Rutman / Perron--Frobenius the
fundamental mode is the unique non-negative eigenvector and the dominant
eigenvalue is real and positive (the only physically meaningful steady state);
all higher harmonics change sign in space.

Layering
========

The eigenvalue machinery is layered so the algorithm is **method-agnostic**
(see the operator-algebra theory page for the full posing table):

* **Operator leaves** (method-specific): :math:`L, C, S, F, B`.
* **Problem posing** arranges the leaves into :math:`(A_{\rm loss}, M)` and a
  :math:`\mu \to` physical-eigenvalue map.  The **k-eigenvalue** row is
  :math:`A_{\rm loss} = L+C-S-B`, :math:`M = F`, :math:`k = \mu`.  (The
  :math:`\alpha`-eigenvalue row — :math:`A_{\rm loss} = L+C-S-F-B`,
  :math:`M = 1/v`, :math:`\alpha = -1/\mu` — and the adjoint row
  :math:`A_{\rm loss}^\dagger \psi^\dagger = \lambda M^\dagger \psi^\dagger`
  are documented future seams.)
* **Resolvent** :math:`A_{\rm loss}^{-1}` (method-specific): the fixed-source
  inner solve — SN sweep / Krylov, CP BiCGSTAB, diffusion FD, ...
* **Algorithm** (this module): :func:`power_iteration` over the
  method-agnostic :class:`EigenvalueSolver` boundary.  It sees ONLY a
  normalized-source fixed-point procedure — never the method's operators or
  sweeps.

:func:`power_iteration` is the **canonical engine**: any deterministic solver
that expresses its physics in the :class:`EigenvalueSolver` Protocol plugs in —
INCLUDING methods (CP, diffusion, homogeneous) whose loss operator is a
monolithic matrix with no :math:`(L, S, F)` factorization.  The operator-triple
entry point :class:`~orpheus.numerics.iteration.KEigenvalue` is one such
implementer: it realizes the boundary from an :math:`(L, S, F)` triple and
delegates the loop here (single source of truth — one loop).

The eigenvector is a **flux distribution** — its shape is determined but its
absolute scale is arbitrary.  Per-step renormalization to unit production rate
(:meth:`EigenvalueSolver.compute_production_rate`) keeps the iterate at
:math:`O(1)` (ERR-052); rescaling to an absolute flux at a target reactor power
is a single multiplication (a future ``target_power`` hook).

Future solution algorithms (full-spectrum Arnoldi / Krylov--Schur,
shift-invert / FEAST for interior eigenvalues) slot in at this same boundary,
dispatched via :attr:`~orpheus.numerics.iteration.KEigenvalue.eigenvalue_method`.
"""

from __future__ import annotations

import warnings
from typing import Protocol, runtime_checkable

import numpy as np


class EigenvalueSolver(Protocol):
    """The method-agnostic boundary that :func:`power_iteration` consumes.

    Any deterministic eigenvalue solver — SN sweep, CP collision matrix,
    diffusion FD, MoC ray-trace, homogeneous direct — satisfies this Protocol
    and plugs into the canonical :func:`power_iteration` algorithm.  The
    algorithm sees ONLY these methods (a normalized-source fixed-point
    procedure), never the method's operators or its fixed-source solve.

    Power iteration structure (each outer iteration):

    1. ``compute_fission_source`` — build the fission RHS from the
       current flux distribution estimate and eigenvalue.
    2. ``solve_fixed_source`` — apply the transport (or diffusion)
       operator to the fission source, returning an updated flux
       distribution.  Scattering and (n,2n) sources are assembled
       **inside** this method because they couple to the transport
       solve (e.g. inner scattering iterations in SN).
    3. ``compute_keff`` — update the eigenvalue from the neutron
       production / loss balance.
    """

    def initial_flux_distribution(self) -> np.ndarray:
        """Return an initial guess for the flux distribution."""
        ...

    def compute_fission_source(
        self,
        flux_distribution: np.ndarray,
        keff: float,
    ) -> np.ndarray:
        """Fission source: Q_f = χ · (νΣ_f · φ) / k_eff."""
        ...

    def solve_fixed_source(
        self,
        fission_source: np.ndarray,
        flux_distribution: np.ndarray,
    ) -> np.ndarray:
        """Apply the transport operator and return an updated flux distribution.

        This method encapsulates the model-specific physics:

        * **Collision probability** — direct matrix multiplication with P_inf.
        * **Discrete ordinates (SN)** — diamond-difference sweep with inner
          scattering iterations.
        * **Method of characteristics** — ray-tracing sweep.
        * **Diffusion** — implicit solve with BiCGSTAB.
        * **Homogeneous** — sparse direct solve of the removal matrix.

        Scattering and (n,2n) sources must be assembled inside this method
        so that inner iteration schemes (e.g. source iteration in SN) can
        update them between sweeps.

        Numerical conditioning (e.g. dividing by max(φ) to prevent overflow)
        is an implementation detail of this method, not physics normalization.
        """
        ...

    def compute_keff(self, flux_distribution: np.ndarray) -> float:
        """Compute the eigenvalue from the neutron balance.

        k_eff = production / (absorption + leakage)

        For lattice models with reflective boundary conditions the leakage
        term is zero.  For whole-core models with vacuum boundary conditions
        (e.g. diffusion) it is non-zero.
        """
        ...

    def converged(
        self,
        keff: float,
        keff_old: float,
        flux_distribution: np.ndarray,
        flux_old: np.ndarray,
        iteration: int,
    ) -> bool:
        """Return True when the outer iteration has converged."""
        ...


@runtime_checkable
class ProductionRateSolver(EigenvalueSolver, Protocol):
    """An :class:`EigenvalueSolver` that has adopted the production-rate
    normalisation contract (ERR-052).

    Modelling the optionality HONESTLY: ``compute_production_rate`` is the
    OPTIONAL extension of the base contract, NOT a required member.
    :func:`power_iteration` renormalises with it when the solver provides it
    and falls back to the un-normalised legacy trajectory otherwise (the
    deprecation window for CP / diffusion / MoC / homogeneous).  Solvers still
    in that window -- e.g. :class:`~orpheus.cp.solver.CPSolver` -- conform to
    :class:`EigenvalueSolver` WITHOUT this method and plug into
    ``power_iteration`` with no suppression.  ``runtime_checkable`` lets the
    driver narrow with ``isinstance`` (which fires under ``python -O``, unlike
    a stripped ``assert``).
    """

    def compute_production_rate(self, flux_distribution: np.ndarray) -> float:
        """Total volume-integrated neutron production rate (scalar).

        Power iteration renormalises :math:`\\phi` to unit production
        rate at each outer step (ERR-052) so the iterate stays at a
        physically natural O(1) magnitude regardless of whether the
        operator is supercritical or subcritical.  Callers rescale to
        absolute flux at a target reactor power via a single
        multiplication.

        Solvers that have not yet adopted the production-rate contract
        may omit this method; :func:`power_iteration` falls back to
        the un-normalised legacy trajectory in that case.
        """
        ...


def power_iteration(
    solver: EigenvalueSolver,
    max_iter: int = 500,
) -> tuple[float, list[float], np.ndarray]:
    """Converge to the dominant eigenvalue and fundamental mode.

    Power iteration converges to the largest eigenvalue k_0 (= k_eff)
    and its eigenvector φ_0 (the fundamental mode).  The convergence
    rate is governed by the dominance ratio :math:`|k_1 / k_0|`.

    Returns
    -------
    keff : float
        Dominant eigenvalue (k_eff).
    keff_history : list[float]
        Eigenvalue estimate at each outer iteration.
    flux_distribution : np.ndarray
        Fundamental mode (arbitrary normalization).
    """
    flux_distribution = solver.initial_flux_distribution()
    keff = 1.0
    keff_history: list[float] = []

    for n in range(1, max_iter + 1):
        flux_old = flux_distribution.copy()
        keff_old = keff

        fission_source = solver.compute_fission_source(flux_distribution, keff)
        flux_distribution = solver.solve_fixed_source(fission_source, flux_distribution)
        # Renormalise to unit production rate so the iterate stays at a
        # physically natural O(1) magnitude across iterations regardless
        # of whether the operator is supercritical (k>1, would grow) or
        # subcritical (k<1, would decay to denormalised FP within ~30-60
        # iterations and the keff ratio would become 0/0 numerically
        # meaningless — ERR-052).  Production rate is scale-invariant in
        # ``keff`` so the converged eigenvalue is unchanged; the
        # converged ``φ`` carries the canonical reactor-physics output
        # convention :math:`\\int \\nu\\Sigma_f\\,\\phi\\,dV = 1`, which
        # makes rescaling to absolute flux at a target power a single
        # multiplication by :math:`P_{\\text{target}} / \\kappa`.
        #
        # Solvers that have not adopted ``compute_production_rate``
        # retain the legacy un-normalised trajectory (the deprecation
        # window for CP / diffusion / MoC / homogeneous migration).
        if isinstance(solver, ProductionRateSolver):
            p = float(solver.compute_production_rate(flux_distribution))
            if p > 0.0:
                flux_distribution = flux_distribution / p
        keff = solver.compute_keff(flux_distribution)
        keff_history.append(keff)

        if solver.converged(keff, keff_old, flux_distribution, flux_old, n):
            break

    return keff, keff_history, flux_distribution


def _sign_normalised(v: np.ndarray) -> np.ndarray:
    r"""Sign-normalise an eigenvector to the physical, non-negative convention.

    A criticality eigenvector is a flux distribution — defined up to scale,
    and in particular up to SIGN.  The family-wide output convention is
    ``v.sum() >= 0`` (the physical, non-negative spectrum); this helper is
    its single source, shared by every engine that emits an eigenvector
    (:func:`dominant_eigenpair` — and through it :func:`direct_eigenvalue` —
    and :func:`rayleigh_quotient_iteration`).
    """
    return v if v.sum() >= 0.0 else -v


def dominant_eigenpair(
    M: np.ndarray,
    *,
    imag_tol: float = 1e-9,
) -> tuple[float, np.ndarray]:
    r"""Dominant (Perron--Frobenius) eigenpair of a materialized resolvent.

    The shared eigen-EXTRACTION primitive of the direct engines: given the
    dense resolvent :math:`M = A^{-1}F` — *however it was formed* — return
    the dominant eigenpair with the criticality contract enforced.  This is
    the ONE home of the Perron--Frobenius validation (taxonomy step 5b): by
    Krein--Rutman / Perron--Frobenius the fundamental mode of a well-posed
    criticality resolvent is the unique non-negative eigenvector and its
    eigenvalue is **real and positive**, so a complex dominant eigenvalue is
    rejected as a malformed problem (Cardinal Rule 1 — fail loud, never
    return a non-eigenvalue).

    Two callers, one validation:

    * :func:`direct_eigenvalue` — forms the resolvent from the posed
      ``(A, F)`` pair via :func:`numpy.linalg.solve`, then delegates here.
    * the homogeneous K-path
      (:func:`~orpheus.homogeneous.solver.solve_homogeneous_infinite`) —
      forms the resolvent through the operator algebra
      ``K = MatrixInverseOperator(loss) @ production`` and materializes it
      with :meth:`~orpheus.numerics.operator.LinearOperator.as_matrix`,
      bypassing the ``(A, F)`` posing boundary entirely.

    Parameters
    ----------
    M : np.ndarray
        The materialized resolvent (square, ``(n, n)``, real) whose
        dominant eigenpair is sought.
    imag_tol : float, optional
        Tolerance on the dominant eigenvalue's imaginary part before it is
        rejected as complex.  Default ``1e-9``.

    Returns
    -------
    k : float
        The dominant eigenvalue :math:`\lambda_{\max}(M)` (largest real
        part; real by the Perron--Frobenius contract).
    phi : np.ndarray
        The corresponding right eigenvector (real-dtype, ``(n,)``),
        sign-normalised so ``phi.sum() >= 0``; its absolute scale is
        arbitrary (rescale to a target production rate / power downstream).

    Raises
    ------
    ValueError
        If ``M`` is not a square 2-D matrix, or the dominant eigenvalue is
        complex (``|Im λ| > imag_tol``).
    """
    M = np.asarray(M, dtype=float)
    if M.ndim != 2 or M.shape[0] != M.shape[1]:
        raise ValueError(
            f"dominant_eigenpair: M must be a square 2-D matrix; got shape {M.shape}."
        )
    eigvals, eigvecs = np.linalg.eig(M)

    # The dominant (Perron–Frobenius) mode: the largest real part.
    dominant = int(np.argmax(np.real(eigvals)))
    lam = eigvals[dominant]
    if abs(float(np.imag(lam))) > imag_tol:
        raise ValueError(
            f"dominant_eigenpair: the dominant eigenvalue is complex "
            f"({lam:.6g}); the resolvent A⁻¹F of a well-posed criticality "
            f"problem has a real dominant eigenvalue (Perron–Frobenius). "
            f"A complex dominant signals a malformed problem."
        )
    k = float(np.real(lam))
    phi = _sign_normalised(np.real(eigvecs[:, dominant]))
    return k, phi


def direct_eigenvalue(
    A: np.ndarray,
    F: np.ndarray,
    *,
    imag_tol: float = 1e-9,
) -> tuple[float, np.ndarray]:
    r"""Exact dominant eigenpair of the generalized eigenproblem
    :math:`A\,\varphi = \tfrac{1}{k}\,F\,\varphi`, via the dense resolvent.

    The **direct (non-iterative) sibling of** :func:`power_iteration`.  Where
    ``power_iteration`` converges to the dominant eigenpair at the
    dominance-ratio rate (the realization for large, sparse, sweep-only loss
    operators that can only be *applied*), ``direct_eigenvalue`` forms the dense
    resolvent :math:`A^{-1}F` and returns the EXACT dominant eigenpair in one
    LAPACK call:

    .. math::

        k \;=\; \lambda_{\max}\!\bigl(A^{-1}F\bigr), \qquad
        A^{-1}F\,\varphi \;=\; k\,\varphi .

    It is the right realization for **small, densifiable** operators — 0-D
    infinite-medium spectra, few-group / few-region problems — where the
    iterative engine would only approximate (at a dominance-ratio-dependent
    rate) an answer the dense solve gives exactly.  (The homogeneous solver
    :func:`~orpheus.homogeneous.solver.solve_homogeneous_infinite` poses that
    problem class through the operator algebra instead and reaches the shared
    kernel via :func:`dominant_eigenpair` directly — taxonomy step 5b — so a
    change here does NOT move its :math:`k_\infty`.)
    Both engines solve the SAME posed problem
    :math:`(A_{\rm loss}, M) = (A, F)`; they differ only in exact-dense vs
    iterative.

    By Perron--Frobenius / Krein--Rutman the fundamental mode of a well-posed
    criticality resolvent is the unique non-negative eigenvector and its
    eigenvalue is **real and positive**; a complex dominant eigenvalue therefore
    signals a malformed ``(A, F)`` and is rejected (Cardinal Rule 1 — fail loud,
    never return a non-eigenvalue).

    Parameters
    ----------
    A : np.ndarray
        The loss matrix :math:`A` (square, ``(n, n)``).  Must be invertible — a
        singular loss matrix is a malformed problem and the underlying
        :func:`numpy.linalg.solve` raises :class:`numpy.linalg.LinAlgError`.
    F : np.ndarray
        The production matrix :math:`F` (same shape as ``A``).  For the
        homogeneous problem this is the rank-1 fission dyad
        :math:`\chi \otimes \nu\Sigma_f`.
    imag_tol : float, optional
        Tolerance on the dominant eigenvalue's imaginary part before it is
        rejected as complex.  Default ``1e-9``.

    Returns
    -------
    k : float
        The dominant eigenvalue :math:`\lambda_{\max}(A^{-1}F)` (real).
    phi : np.ndarray
        The corresponding right eigenvector (real-dtype, ``(n,)``),
        sign-normalised so ``phi.sum() >= 0`` (a physical, non-negative
        spectrum); its absolute scale is arbitrary (rescale to a target
        production rate / power downstream).

    Raises
    ------
    ValueError
        If ``A`` / ``F`` are not square matrices of equal shape, or the dominant
        eigenvalue is complex (``|Im λ| > imag_tol`` — raised by
        :func:`dominant_eigenpair`, the shared extraction home).
    numpy.linalg.LinAlgError
        If ``A`` is singular (propagated from :func:`numpy.linalg.solve`).

    See Also
    --------
    power_iteration : the iterative sibling (large sweep-only operators).
    dominant_eigenpair : the shared Perron--Frobenius extraction primitive
        this engine delegates to after forming the resolvent.

    Notes
    -----
    A third realization — **Rayleigh-Quotient Iteration** (the bordered /
    augmented-matrix form, in which a previous-iterate eigenvector enters as a
    normalisation row and the eigenvalue update falls out of one linear solve) —
    sits between these two: iterative like ``power_iteration`` but
    *superlinearly* convergent, for large operators where a fast eigenvalue
    refinement from a flux estimate is wanted.  It is a documented future seam
    (alongside the Arnoldi / FEAST hooks above), not yet built.
    """
    A = np.asarray(A, dtype=float)
    F = np.asarray(F, dtype=float)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError(
            f"direct_eigenvalue: A must be a square 2-D matrix; got shape {A.shape}."
        )
    if F.shape != A.shape:
        raise ValueError(
            f"direct_eigenvalue: A and F must have the same shape; "
            f"got A {A.shape}, F {F.shape}."
        )

    # The dense resolvent A⁻¹F.  np.linalg.solve raises LinAlgError on a
    # singular A — a malformed loss matrix, left to fail loud.  The
    # Perron–Frobenius extraction (argmax-real selection, complex-dominant
    # rejection, sign normalisation) lives in ONE home — the shared
    # :func:`dominant_eigenpair` primitive (taxonomy step 5b).
    return dominant_eigenpair(np.linalg.solve(A, F), imag_tol=imag_tol)


def rayleigh_quotient_iteration(
    A: np.ndarray,
    F: np.ndarray,
    *,
    v0: np.ndarray | None = None,
    tol: float = 1e-12,
    max_iter: int = 50,
) -> tuple[float, np.ndarray]:
    r"""Refine an eigenvector estimate to the NEAREST eigenpair of :math:`A^{-1}F`
    by bordered Rayleigh-Quotient Iteration — *superlinearly* convergent.

    The **third realization** of the eigenvalue boundary, sitting between the
    iterative :func:`power_iteration` (linear convergence, the dominant mode of
    large sweep-only operators) and the exact dense :func:`direct_eigenvalue`
    (one LAPACK shot, small operators).  RQI takes an eigenvector estimate
    :math:`v`, forms the Rayleigh-quotient shift

    .. math::

        k \;=\; \frac{v^{T} F\,v}{v^{T} A\,v}

    (the production/loss ratio — the eigenvector *as a row*), and refines the
    pair :math:`(v, k)` by a Newton step on the eigenpair residual
    :math:`[\,A^{-1}F v - k v;\; c^{T} v - 1\,]` written in the well-conditioned
    **bordered** form, in which the previous-iterate eigenvector :math:`c = v`
    enters as the normalisation **row**:

    .. math::

        \begin{bmatrix} F - kA & -Av \\ v^{T} & 0 \end{bmatrix}
        \begin{bmatrix} \Delta v \\ \Delta k \end{bmatrix}
        \;=\;
        \begin{bmatrix} -(F - kA)\,v \\ 0 \end{bmatrix} .

    The first block-row is the eigen-residual Newton equation premultiplied by
    :math:`A` (so :math:`M = A^{-1}F` is never formed); the bottom row is the
    normalisation constraint :math:`v^{T}\Delta v = 0` (the Newton update stays
    tangent to the unit sphere).  The border keeps the
    :math:`(n{+}1)\times(n{+}1)` system non-singular even as :math:`F - kA`
    becomes singular at convergence — precisely where raw shifted inverse
    iteration breaks down.  Convergence is locally **quadratic** for the
    non-symmetric resolvent :math:`A^{-1}F` (cubic for symmetric / normal).

    .. warning::

        RQI converges to the eigenpair **nearest the initial Rayleigh
        quotient**, NOT necessarily the dominant one.  For the dominant
        (k-eigenvalue) mode, warm-start ``v0`` near it — a few
        :func:`power_iteration` steps, or a prior outer iterate.  With the
        default ``v0`` (all-ones) the result is whichever eigenvalue is nearest
        :math:`(\mathbf{1}^{T} F \mathbf{1}) / (\mathbf{1}^{T} A \mathbf{1})`.

    Parameters
    ----------
    A, F : np.ndarray
        The generalized pair :math:`A\,\varphi = (1/k)\,F\,\varphi` (square,
        equal shape).
    v0 : np.ndarray, optional
        Initial eigenvector estimate.  Default all-ones.  Determines WHICH
        eigenpair is found (the one nearest its Rayleigh quotient).
    tol : float, optional
        Convergence tolerance on the Newton step :math:`\lVert\Delta v\rVert`.
        Default ``1e-12``.
    max_iter : int, optional
        Maximum RQI steps.  Default ``50``; RQI converges in a handful of steps
        from a warm start, so reaching this bound signals a poor ``v0`` or a
        defective pair (a :class:`RuntimeWarning` is emitted and the best
        iterate returned).

    Returns
    -------
    k : float
        The converged eigenvalue (real — the iteration is real arithmetic
        throughout).
    phi : np.ndarray
        The converged eigenvector (real-dtype), sign-normalised so
        ``phi.sum() >= 0``.

    Raises
    ------
    ValueError
        If ``A`` / ``F`` are not square matrices of equal shape.

    See Also
    --------
    power_iteration : iterative, linear, dominant mode (large operators).
    direct_eigenvalue : exact dense dominant eigenpair (small operators).
    """
    A = np.asarray(A, dtype=float)
    F = np.asarray(F, dtype=float)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError(
            f"rayleigh_quotient_iteration: A must be a square 2-D matrix; "
            f"got shape {A.shape}."
        )
    if F.shape != A.shape:
        raise ValueError(
            f"rayleigh_quotient_iteration: A and F must have the same shape; "
            f"got A {A.shape}, F {F.shape}."
        )
    n = A.shape[0]

    v = np.ones(n) if v0 is None else np.asarray(v0, dtype=float).ravel().copy()
    v = v / np.linalg.norm(v)
    # Initial Rayleigh-quotient shift k = vᵀMv (v unit) with M = A⁻¹F.
    k = float(v @ np.linalg.solve(A, F @ v))

    border = np.zeros((n + 1, n + 1))
    rhs = np.zeros(n + 1)
    converged = False
    step_norm = float("inf")
    for _ in range(max_iter):
        shifted = F - k * A
        border[:n, :n] = shifted
        border[:n, n] = -(A @ v)
        border[n, :n] = v
        rhs[:n] = -(shifted @ v)
        # rhs[n] = 0: with v unit, the normalisation residual vᵀv − 1 vanishes.
        step = np.linalg.solve(border, rhs)
        v = v + step[:n]
        k = k + float(step[n])
        v = v / np.linalg.norm(v)
        step_norm = float(np.linalg.norm(step[:n]))
        if step_norm <= tol:
            converged = True
            break

    if not converged:
        warnings.warn(
            f"rayleigh_quotient_iteration: no convergence in max_iter="
            f"{max_iter} (last ‖Δv‖={step_norm:.3e}, tol={tol:.1e}). "
            f"RQI converges in a handful of steps from a warm start; a poor v0 "
            f"(far from any eigenvector) or a defective pair can stall it. "
            f"Returning the best iterate.",
            RuntimeWarning,
            stacklevel=2,
        )

    phi = _sign_normalised(v)  # the family-wide non-negative-spectrum convention
    return k, phi
