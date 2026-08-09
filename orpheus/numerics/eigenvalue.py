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
  :math:`A_{\rm loss} = L+C-S-B`, :math:`M = F`, :math:`k = \mu`.  The
  **adjoint row is LIVE** (#276 A4/A5): the daggered posing
  :math:`A_{\rm loss}^\dagger \psi^* = \frac{1}{k} F^\dagger \psi^*` is
  realized purely by DAGGER-ing the triple —
  ``KEigenvalue((L+C).H, (S+B).H, F.H)`` (the daggered RESOLVENT, gain,
  and fission; the loss dagger is formed inside) — and runs through THIS SAME loop
  unchanged (:func:`~orpheus.sn.solver.solve_sn_adjoint` is the SN
  entry; :math:`k^\dagger = k` is an exact algebraic identity, and the
  daggered vector is certified by the defining-law residual rows —
  see the SN adjoint theory chapter).  (The :math:`\alpha`-eigenvalue
  row — :math:`A_{\rm loss} = L+C-S-F-B`, :math:`M = 1/v`,
  :math:`\alpha = -1/\mu` — remains a documented future seam.)
* **Resolvent** :math:`A_{\rm loss}^{-1}` (method-specific): the fixed-source
  inner solve — SN sweep / Krylov; CP's ``P_inf`` collision-probability kernel;
  diffusion's and homogeneous's eager LU
  (:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`).
* **Algorithm** (this module): :func:`power_iteration` over the
  method-agnostic :class:`EigenvalueSolver` boundary.  It sees ONLY a
  normalized-source fixed-point procedure — never the method's operators or
  sweeps.

:func:`power_iteration` is the **canonical engine**: any deterministic solver
that expresses its physics in the :class:`EigenvalueSolver` Protocol plugs in.
The Protocol binds the inner resolvent **late** because the inner-solve
*strategy* and its *carrier* vary — a sweep, an eager LU, a single Jacobi step,
or Gauss-Seidel with inners; a typed composite or a bare ``(N, ng)`` array —
**not** because the loss operator resists factorization.  Every deterministic
method ORPHEUS ships DOES factor: SN and diffusion as :math:`L+C-S-B`
(``DiffusionSolver.loss = leakage + collision - scattering - boundary``),
homogeneous as :math:`A = C - K_{\rm iso}`, and CP applies :math:`S`,
:math:`(n,2n)`, and :math:`F` explicitly *outside* its ``P_inf`` kernel — that
kernel IS its :math:`(L+C)^{-1}`.  The operator-triple entry point
:class:`~orpheus.numerics.iteration.KEigenvalue` is one such implementer: it
realizes the boundary from an :math:`(L, S, F)` triple and delegates the loop
here (single source of truth — one loop).

The eigenvector is a **flux distribution** — its shape is determined but its
absolute scale is arbitrary.  Per-step renormalization to unit production rate
(:meth:`EigenvalueSolver.compute_production_rate`) keeps the iterate at
:math:`O(1)` (ERR-052); rescaling to an absolute flux at a target reactor power
is a single multiplication (a future ``target_power`` hook).

Future solution algorithms (full-spectrum Arnoldi / Krylov--Schur,
shift-invert / FEAST for interior eigenvalues) slot in at this same boundary,
dispatched via the ``eigenvalue_method`` constructor selector on
:class:`~orpheus.numerics.iteration.KEigenvalue` (today only ``"power"`` is
implemented; any other value raises at construction).
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Generic, Protocol, cast, runtime_checkable

import numpy as np

from .convergence import IterationRecord, StoppingCriterion
from .vector import Carrier, Vector


MINIMUM_OUTER_ITERATIONS = 3
"""Outers below which power iteration refuses to claim convergence.

The stop test reads INCREMENTS (``|Δk|``, ``‖Δφ‖/‖φ‖``), and one increment is
a difference against an arbitrary initial guess — it can be small because the
iteration is converging or because the guess happened to sit near the first
iterate.  Two increments are the minimum that can show a trend, so the first
claimable outer is the third.

This is a property of the ALGORITHM, not of any solver, which is why it lives
here and reaches :attr:`~orpheus.numerics.convergence.IterationRecord.min_iterations`
from one place.  `[M]` 2026-08-09 it was transcribed as ``if iteration <= 2:
return False`` in all five realizers of :class:`EigenvalueSolver` — identical
in every one, and invisible to any single-file review.
"""


class EigenvalueSolver(Protocol[Carrier]):
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

    def initial_flux_distribution(self) -> Carrier:
        """Return an initial guess for the flux distribution."""
        ...

    def compute_fission_source(
        self,
        flux_distribution: Carrier,
        keff: float,
    ) -> Carrier:
        """Fission source: Q_f = χ · (νΣ_f · φ) / k_eff."""
        ...

    def solve_fixed_source(
        self,
        fission_source: Carrier,
        flux_distribution: Carrier,
    ) -> Carrier:
        r"""Apply the transport operator and return an updated flux distribution.

        This method encapsulates the model-specific physics:

        * **Collision probability** — direct matrix multiplication with P_inf.
        * **Discrete ordinates (SN)** — diamond-difference sweep with inner
          scattering iterations.
        * **Method of characteristics** — ray-tracing sweep.
        * **Diffusion** — eager LU of the composed loss :math:`L+C-S-B`.
        * **Homogeneous** — eager LU of :math:`A = C - K_{\rm iso}`.

        Scattering and (n,2n) sources must be assembled inside this method
        so that inner iteration schemes (e.g. source iteration in SN) can
        update them between sweeps.

        Numerical conditioning (e.g. dividing by max(φ) to prevent overflow)
        is an implementation detail of this method, not physics normalization.
        """
        ...

    def compute_keff(self, flux_distribution: Carrier) -> float:
        """Compute the eigenvalue from the neutron balance.

        k_eff = fission production / (absorption + leakage − (n,2n) emission)

        The estimator MUST be the eigenvalue of the problem
        ``solve_fixed_source`` poses (#259 P1 / R7): every implementer
        scales ONLY fission by 1/k — scattering and the (n,2n) emission
        enter the inner solve as plain gains — so the numerator is
        fission-only and the (n,2n) gain sits on the net-removal side.
        For lattice models with reflective boundary conditions the
        leakage term is a structural zero; for vacuum-bounded models
        (diffusion, vacuum-bounded SN) it is non-zero (#291).
        """
        ...

    def measure_stopping_criteria(
        self,
        keff: float,
        keff_old: float,
        flux_distribution: Carrier,
        flux_old: Carrier,
    ) -> tuple[StoppingCriterion, ...]:
        r"""Read every quantity this solver stops on, at this iterate.

        Returns one :meth:`~orpheus.numerics.convergence.StoppingCriterion.reading`
        per criterion — the measured magnitude AND the tolerance it is judged
        against, welded together so neither can travel without the other.
        :func:`power_iteration` accumulates the readings into trajectories and
        DERIVES the verdict; the solver states what it measured and does not
        decide.

        ⛔ This replaced ``converged(...) -> bool`` on 2026-08-09 (#340 N2b).
        The predicate computed ``|Δk|`` and ``‖Δφ‖/‖φ‖``, compared both, and
        returned one bit — so a stalled solve could not say WHICH criterion
        was lagging, how fast it was closing, or what budget would reach it.
        Same lossy-return-type defect the inner drivers carried until N2a
        (`[M]` #340 F2), one level up.

        Each solver reports in its own metric and that is the point: CP judges
        the flux change in :math:`\ell^\infty`, SN / MoC / diffusion in
        relative :math:`\ell^2`.  A loop that computed the readings itself
        would have to pick one, and would then be a twin of all five.

        ⚠ There is deliberately NO ``iteration`` parameter.  A reading is a
        function of the two iterates and nothing else; the "don't claim
        convergence before iteration 3" rule is a property of power iteration
        (you need at least two increments to see a trend), not of any solver,
        and it lives once in :data:`MINIMUM_OUTER_ITERATIONS`.  `[M]` it was
        transcribed identically in all FIVE realizers before this change.
        Removing the parameter is what makes the sixth transcription
        unspellable.
        """
        ...


@runtime_checkable
class RecordingSolver(EigenvalueSolver[Carrier], Protocol[Carrier]):
    """An :class:`EigenvalueSolver` that retains its inner solves' records.

    The OPTIONAL extension that lets the outer record carry a subtree, modelled
    exactly like :class:`ProductionRateSolver` below: solvers that keep their
    inner trajectories (SN via
    :class:`~orpheus.numerics.iteration.KEigenvalue`) expose them here and
    :func:`power_iteration` narrows with ``isinstance``; CP / MoC / diffusion
    conform to the base contract without it and plug in with no suppression.

    Without this the outer level could only ever report itself, and "the outer
    stalled because its inner starved" — the #340 failure that motivated the
    whole campaign — would stay unanswerable from the returned object.
    """

    @property
    def inner_records(self) -> Sequence[IterationRecord]:
        """Every inner solve this instance has run, in order, appended-to.

        ⚠ Deliberately NOT "the records of the current solve".  A realizer is
        free to accumulate across solves — :func:`power_iteration` slices off
        what was appended during ITS loop and never reads the earlier entries,
        so an instance reused for two solves cannot contribute a stale child.
        Requiring a reset instead would put the correctness of the tree in the
        hands of five separate realizers, each of which would have to find its
        own "start of solve" hook; `[M]` #340 F12 measured exactly that failure
        (``SNSolver`` had no such hook and its counter double-counted).
        """
        ...


@runtime_checkable
class ProductionRateSolver(EigenvalueSolver[Carrier], Protocol[Carrier]):
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

    def compute_production_rate(self, flux_distribution: Carrier) -> float:
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


@dataclass(frozen=True)
class PowerIterationOutcome(Generic[Carrier]):
    r"""What :func:`power_iteration` produced, INCLUDING why it stopped.

    The loop knows whether it broke on the convergence test or ran out of
    iterations.  Before this type existed the return was a bare
    ``(keff, keff_history, flux)`` triple, so that fact was discarded at the
    source and every one of the five consumers had to reconstruct it —
    ``solve_sn_adjoint`` inferred it from ``len(keff_history) < max_outer``,
    ``solve_sn`` hardcoded ``converged=True`` (issue #342), and CP / MoC /
    diffusion did not report it at all.  A primitive that throws away a fact
    its callers need does not have a documentation problem; it has a return
    type problem.

    Deliberately NOT tuple-unpackable.  Making it destructure like the old
    triple would preserve the very idiom that lost the flag, and would let a
    consumer keep ignoring it by accident.
    """

    keff: float
    """Dominant eigenvalue :math:`k_{\\rm eff}`."""

    keff_history: list[float]
    """Eigenvalue estimate at each outer iteration, in order.

    The eigenvalue's own trajectory — a physics output, NOT a stopping
    criterion.  What the stop test reads is its per-iteration INCREMENT, which
    lives in :attr:`record` under the name the solver gave it.  Keeping the two
    apart is why the loop no longer has to re-difference this list to say what
    it was judging (`[M]` #340: the warning path did exactly that, and could
    only ever recover ONE of the two criteria that way).
    """

    flux_distribution: Carrier
    """Fundamental mode (unit production rate where the solver supports it)."""

    record: IterationRecord
    """Everything the loop measured: per-criterion trajectories, the budget it
    was given, and — for a :class:`RecordingSolver` — one child record per
    inner solve.

    This is the object that answers "where did it stall, and what do I set".
    :attr:`converged` is one bit of it.
    """

    @property
    def converged(self) -> bool:
        """Did the OUTER loop meet all of its own criteria?

        DERIVED from :attr:`record`, never stored — the campaign's founding
        rule (#342).  A stored flag is a transcription, and a transcription of
        a convergence verdict is the exact defect this type was minted to fix:
        it was written by hand at five sites, one of them the literal
        ``converged=True``.  Deriving it makes the sixth unspellable.

        ``False`` means the iteration cap was reached with a criterion unmet
        and :attr:`flux_distribution` is a **best-effort** iterate — mid-descent,
        not the answer.  A caller that asserts physics against it is asserting
        an arbitrary point on the trajectory.

        ⚠ Scoped to THIS level.  A converged outer whose inners starved reads
        ``True`` here and ``False`` at
        :attr:`~orpheus.numerics.convergence.IterationRecord.fully_converged`,
        and that gap is not a wart — it is the #340 headline, since an
        increment-only outer stop CANNOT see an upstream throttle (a truncated
        inner suppresses the very increments the outer reads, so it stalls and
        calls the stall convergence).
        """
        return self.record.converged


def power_iteration(
    solver: EigenvalueSolver[Carrier],
    max_iter: int = 500,
) -> PowerIterationOutcome[Carrier]:
    """Converge to the dominant eigenvalue and fundamental mode.

    Power iteration converges to the largest eigenvalue k_0 (= k_eff)
    and its eigenvector φ_0 (the fundamental mode).  The convergence
    rate is governed by the dominance ratio :math:`|k_1 / k_0|`.

    Returns
    -------
    PowerIterationOutcome
        The eigenpair, the eigenvalue trajectory, and — the reason this is
        an object rather than a triple — whether the loop actually
        converged or merely ran out of iterations.
    """
    flux_distribution = solver.initial_flux_distribution()
    keff = 1.0
    keff_history: list[float] = []
    # Accumulated per criterion NAME, so a solver may report any number of
    # them and the loop never has to know which.  There is no `converged`
    # local any more: exhausting the budget is not a flag to default to False,
    # it is what the trajectories SAY when the last readings did not clear.
    criteria: dict[str, StoppingCriterion] = {}
    # How many inner records the solver had ALREADY accumulated before this
    # loop began.  Slicing from here is what makes a reused solver instance
    # unable to contribute a stale child, without asking any realizer to
    # implement a reset hook (see RecordingSolver.inner_records).
    inners_before = (
        len(solver.inner_records) if isinstance(solver, RecordingSolver) else 0
    )

    for n in range(1, max_iter + 1):
        # Stash the previous iterate for the convergence test.  Typed
        # carriers are frozen (immutable — an alias is a faithful stash);
        # a bare ndarray is mutable and a solver may write into it, so it
        # keeps the defensive copy (the SourceIteration guess idiom).
        flux_old = (
            cast(Carrier, flux_distribution.copy())
            if isinstance(flux_distribution, np.ndarray)
            else flux_distribution
        )
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
                # ``Carrier`` is deliberately unbounded (see
                # :data:`~orpheus.numerics.vector.Carrier` — numpy's stubs
                # cannot prove ndarray ⊨ Vector), so the scalar division —
                # part of the runtime Vector contract every carrier
                # honours — is stated through the cast pair.
                flux_distribution = cast(
                    Carrier, cast(Vector, flux_distribution) / p,
                )
        keff = solver.compute_keff(flux_distribution)
        keff_history.append(keff)

        readings = solver.measure_stopping_criteria(
            keff, keff_old, flux_distribution, flux_old,
        )
        if not readings:
            # An empty conjunction is vacuously true, so a solver that measures
            # nothing would "converge" at MINIMUM_OUTER_ITERATIONS with an
            # empty record to show for it — a #342-class lie assembled out of
            # nothing.  Refuse instead of certifying silence.
            raise ValueError(
                f"{type(solver).__name__}.measure_stopping_criteria returned "
                f"no criteria — a loop with nothing to measure cannot "
                f"converge, and an empty conjunction would claim it did"
            )
        for reading in readings:
            seen = criteria.get(reading.name)
            criteria[reading.name] = (
                reading if seen is None else seen.extended_with(reading)
            )

        # The verdict is the conjunction of THIS iteration's readings, and the
        # min-outer rule is applied here rather than inside any solver (see
        # MINIMUM_OUTER_ITERATIONS).  It is stated again as `min_iterations` on
        # the record below, where it is the same number serving the same rule —
        # `IterationRecord.converged` must reach the identical verdict off the
        # trajectories alone, since that derived reading is what every consumer
        # actually sees.
        if n >= MINIMUM_OUTER_ITERATIONS and all(r.cleared for r in readings):
            break

    return PowerIterationOutcome(
        keff=keff,
        keff_history=keff_history,
        flux_distribution=flux_distribution,
        record=IterationRecord(
            label="outer(power-iteration)",
            criteria=tuple(criteria.values()),
            budget=max_iter,
            # One reading per outer, so this equals every criterion's length —
            # stated anyway because only the producer knows that (the SI inner
            # measures DIFFERENCES and runs one more pass than it records;
            # #340 F10/F11).  The record's co-indexing invariant checks the two
            # against each other.
            iterations_run=len(keff_history),
            min_iterations=MINIMUM_OUTER_ITERATIONS,
            children=(
                tuple(solver.inner_records[inners_before:])
                if isinstance(solver, RecordingSolver)
                else ()
            ),
        ),
    )


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
