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

from typing import Protocol

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
        production = getattr(solver, "compute_production_rate", None)
        if production is not None:
            p = float(production(flux_distribution))
            if p > 0.0:
                flux_distribution = flux_distribution / p
        keff = solver.compute_keff(flux_distribution)
        keff_history.append(keff)

        if solver.converged(keff, keff_old, flux_distribution, flux_old, n):
            break

    return keff, keff_history, flux_distribution
