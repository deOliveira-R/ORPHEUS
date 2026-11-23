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
  realized purely by DAGGER-ing the split — ``KEigenvalue(posing.H(),
  GreenOperator((L+C).H, (S+B).H))`` (the daggered posing and the Green of
  the daggered split, whose relaxation is the inner) — and runs through THIS SAME loop
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
(:meth:`ProductionRateSolver.compute_production_rate` — the OPTIONAL extension
of the base contract, deliberately *not* a member of it) keeps the iterate at
:math:`O(1)` (ERR-052); a solver still outside that contract falls back to the
un-normalised trajectory.  Rescaling to an absolute flux at a target reactor
power is a single multiplication (a future ``target_power`` hook).

Future solution algorithms (full-spectrum Arnoldi / Krylov--Schur,
shift-invert / FEAST for interior eigenvalues) slot in at this same boundary,
dispatched via the ``eigenvalue_method`` constructor selector on
:class:`~orpheus.numerics.iteration.KEigenvalue` (today only ``"power"`` is
implemented; any other value raises at construction).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Generic, Protocol, cast, runtime_checkable

import numpy as np

from .convergence import IterationBudget, IterationRecord, StoppingCriterion
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


class EigenvalueSolver(Protocol[Carrier]):  # concept: driver | role: head | kind: morphism
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
class RecordingSolver(EigenvalueSolver[Carrier], Protocol[Carrier]):  # concept: driver | role: arm | kind: morphism
    """An :class:`EigenvalueSolver` that retains its inner solves' records.

    Structurally OPTIONAL, modelled exactly like :class:`ProductionRateSolver`
    below — :func:`power_iteration` narrows with ``isinstance`` and falls back
    to an empty subtree — but as of #340 N4 **every shipped family declares
    it**, because a family that does not is a family whose starved inner is
    invisible:

    * SN, via :class:`~orpheus.numerics.iteration.KEigenvalue`;
    * CP — ``ng`` children per outer under Gauss-Seidel (the per-group
      within-group scatter solves); NONE under the default Jacobi mode, which
      is the honest answer rather than an omission, since Jacobi deliberately
      LAGS the scattering source and has no inner level to converge;
    * MoC — one child per outer, declaring BOTH readings its sweep loop is
      driven by (scalar flux AND boundary angular flux);
    * diffusion — one ``DIRECT`` child per outer: an exact LU
      back-substitution did not iterate and DID converge, recorded so that
      *"no inner level"* stays distinguishable from *"never recorded"*.

    ⛔ Until 2026-08-11 this paragraph read "CP / MoC / diffusion conform to
    the base contract without it and plug in with no suppression". That was
    true, and it was the #340 F6 hole stated as a design: all three built a
    correct record and dropped it three lines later, so ``converged`` was
    unanswerable outside SN.

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
class ProductionRateSolver(EigenvalueSolver[Carrier], Protocol[Carrier]):  # concept: driver | role: arm | kind: morphism
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
class PowerIterationOutcome(Generic[Carrier]):  # concept: solution | role: twin | kind: record
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
    budget_name: str = "max_iter",
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
            # One unit of ``max_iter`` is one outer, so the identity
            # conversion is the honest statement here.
            budget=IterationBudget(max_iter, budget_name),
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
