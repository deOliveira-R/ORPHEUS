r"""The shared vocabulary for *"did this iteration reach its promise?"*.

An iterative solver can stop for two structurally different reasons, and a
result object that cannot tell them apart is a trap:

* it **converged** — the stopping criterion was met, and the answer is the
  one the caller asked for;
* its **budget ran out** — the loop hit ``max_iter`` with the criterion
  unmet, and the answer is a *best-effort iterate*, mid-descent.

Both come back through the same return type, so unless the distinction is
carried explicitly a caller cannot tell a converged answer from an arbitrary
point on the way to one.

Why this module exists — the measured failure
=============================================

`[M]` 2026-08-08.  ``test_d3_pure_absorber_per_ordinate_psi_exact`` asserted
a closed-form identity to ``rtol=1e-10`` on an all-reflective 3-D box.  That
problem needs **1631** source-iteration sweeps; ``max_inner`` defaulted to
**1000**.  The solve returned the 999th iterate — honestly flagged
``history.converged = False`` — and the gate, which never read the flag,
asserted physics against it.  The error was ``3.287e-10``, and it was
*bit-identical* at every ``inner_tol`` from ``1e-9`` to ``1e-15``, because
the running residual (``1.185e-09``) never fell below even the loosest of
them: all four runs hit the same cap and returned the same bytes.

Two lessons are worth carrying forward, both counter-intuitive:

#. **Tolerance-insensitivity is the signature of a budget truncation, not of
   a discretization floor.**  If an error does not move when you tighten the
   tolerance, the first thing to check is whether the tolerance ever *bound*
   — read the iteration count against the cap before concluding anything
   about the discretization.
#. **A value gate that does not assert convergence is asserting an arbitrary
   iterate.**  It will be green or red by luck.  That gate had ridden an
   unconverged exit since it was written; it passed for months because the
   truncated error happened to land inside the tolerance, and a *correct*
   quadrature change (#337) later moved it out.

The remedy has two halves.  The **structural** half is that a loop which
knows why it stopped must say so in its return value — see
:class:`~orpheus.numerics.eigenvalue.PowerIterationOutcome`, and the
``converged`` field every solver result now derives rather than asserts.
The **loudness** half is this module's warning.

Why a warning and not an exception
==================================

Raising would be louder, and it is the wrong default here.  The project has
already taken this decision once, for the same defect class: when
``KrylovAcceleration`` stopped discarding scipy's ``info`` flag
(``orpheus/numerics/iteration.py``, D-H.1e / ERR-053), the ruling was

    *"Both surface as warnings — raising would break long-standing callers
    that tolerate slow convergence and need the best-effort iterate."*

That reasoning holds at the public solver entries too: rate studies harvest
the residual history rather than the answer, spy tests discard the result
entirely, and a diverging configuration is sometimes exactly what is being
measured.  `[M]` an audit of the tree found **zero** production callers and
**zero** ``examples/`` callers that depend on a truncated answer, but six
legitimate test consumers that do.

The warning is nonetheless *escalatable*, which is what makes it a gate
rather than a decoration.  ``pyproject.toml`` sets no ``filterwarnings`` and
no ``-W error``, so emitting it costs nothing today, and CI (or any careful
caller) turns the whole category into a hard failure with one flag::

    python -O -m pytest -W error::orpheus.numerics.convergence.ConvergenceWarning

That is the same named-warning + opt-in-escalation recipe the regression
suite already uses for its bit-identity tripwire.

.. important::

   ⛔ **The category must be spelled DOTTED.**  This recipe shipped on
   2026-08-09 as ``-W error::ConvergenceWarning`` at four sites (including
   inside the emitted message itself) and **that string does not parse** —
   Python resolves an undotted ``-W`` category against ``builtins``, so the
   flag dies with ``AttributeError: module 'builtins' has no attribute
   'ConvergenceWarning'`` and pytest collects **zero** tests.  The CI
   contract was imaginary for as long as it was published.

   It survived review because the gate that "proved" it,
   ``test_it_is_escalatable_to_an_error``, installs the filter through
   ``warnings.simplefilter`` — the in-process API.  That is a true claim
   about the *category* and says nothing about the *spelling*, and the
   doc, the runtime message and the test all agreed on a spelling no
   interpreter accepts.

   Hence :data:`ESCALATION_FLAG` below: the spelling is now **derived from
   the class**, not retyped, so a module move or a rename cannot desync it,
   and one gate parses the string itself.

Two emission points, at two different levels
============================================

⛔ **This category is raised from exactly two places, and they must NOT be
consolidated** (#340 N4.7, 2026-08-11).  They answer different questions and
a pass that "unifies the ConvergenceWarning emission points" — a natural
tidying instinct, since the category is shared — would destroy one of them:

.. list-table::
   :header-rows: 1
   :widths: 12 30 58

   * - level
     - site
     - the proposition it asserts
   * - **ENTRY**
     - :func:`warn_if_unconverged`, called by every public solver entry
     - *"the iteration TREE this call returns is not fully converged"* —
       names the level that failed (``first_failure`` searches children
       first), its budget, its knob, and the budget its observed rate
       projects.  Fires **at most once per solve.**
   * - **LEAF**
     - :meth:`~orpheus.numerics.iteration.KrylovAcceleration.solve`
     - *"SciPy's ``gmres`` returned ``info != 0`` for THIS inner solve"* —
       mechanism-specific, and the only place that fact is observable.
       Fires per inner solve.

Merging them either drowns the tree verdict in per-inner noise (if the leaf
wins) or silently discards the only report of a SciPy-level failure (if the
entry wins).  The right relationship is the one that exists: the leaf reports
a mechanism, the entry reports a contract, and the entry's record is what
carries the leaf's truncation up to the caller.

What a stalled solve owes its caller
====================================

The warning above is the *loudness* half, and on its own it is thin: it says
one level of one solve ran out of budget.  A real convergence problem is
almost never legible at that resolution, because the object that failed is a
**tree** — an outer power iteration over an inner fixed-source solve — and
the question a user actually has is *which level* stalled, on *which
criterion*, and *what budget would have sufficed*.

The tree is **ragged, and its shape is a property of the ENTRY POINT rather
than of the solver**, which is why the record is recursive instead of a fixed
outer/inner pair.  `[M]` 2026-08-09, traced on real solves: SN is two levels
deep and never three — ``build_within_group_system`` builds ONE
:class:`~orpheus.numerics.iteration.SourceIteration` over the whole
``(N, ng, nx, ny)`` state, so "within-group" there means *fission-external*,
not *one group at a time*.  CP genuinely does loop per group.  MoC's inner is
a fixed sweep count with no tolerance and no residual — a level with no
criteria at all.  Diffusion's inner is one LU back-substitution, so its record
is a **one-level** tree.  A type that assumed any single one of those shapes
would have to grow an arm for each of the others.

:class:`IterationRecord` is that tree, and :class:`StoppingCriterion` is one
quantity-driven-below-one-tolerance within it.  The load-bearing property is
that **nothing about convergence is stored** — every verdict is derived from
the trajectory:

* a fact you cannot transcribe cannot be transcribed *wrongly*, which is the
  whole of #342 (five sites each spelling their own ``converged``, one of
  them a hardcoded ``True``) made unspellable rather than merely fixed;
* the derivation is the loop's own stopping predicate re-evaluated on the
  final state, never a transcription of which way control flow left the loop.

`[M]` 2026-08-09, the measured failure this shape exists to answer to: on a
2-D all-reflective 2-group eigenvalue solve (LS-S4, ``keff_tol=1e-7``,
``inner_tol=1e-10``) a throttled inner moved ``keff`` by **5.3x keff_tol**
while the returned object reported ``converged = True`` and emitted nothing.
The outer's stop test is *entirely increments* (``dk`` and ``dphi``), so a
starved inner suppresses the very quantities the outer reads: it does not
converge, it **stalls and reads its own stall as convergence**.  A flat
boolean has nowhere to put that fact; :attr:`IterationRecord.fully_converged`
is where it goes.

⭐ **And a level owes the reader one fact its trajectory cannot supply: the
NAME of the knob that capped it** (:attr:`IterationBudget.name`, #340
N6).  Everything else in a diagnostic is derivable — which level, which
criterion, how far, what rate, what budget would have sufficed — but *"set
what?"* is answerable only by the site that constructed the level, because
the same :class:`~orpheus.numerics.iteration.SourceIteration` is somebody's
``max_inner`` here and its own ``max_iter`` there.  So the producer states it
and every consumer reads it, rather than each consumer being handed a name
that describes the level IT was thinking about — which is the defect N6
removed: the SN entries passed their own top-level knob, so on a tree the
advice named a knob that could not help.

One state per question, not one predicate for all of them
---------------------------------------------------------

`[M]` the same campaign's audit reused a production predicate as a detector
and got **44 of 90 rows wrong**, all in the flattering direction: the
predicate was also false for an *empty* history — which for GMRES means
"converged on the initial guess", not "truncated".  A predicate is written to
make one decision; it does not classify states.

So the record exposes the states themselves — :attr:`~IterationRecord.iterated`,
:attr:`~IterationRecord.converged`, :attr:`~IterationRecord.exhausted_budget`,
:attr:`~IterationRecord.truncated` — and each consumer names the one it
means.  Their four combinations are exactly the four ways a loop can stop,
and :attr:`~IterationRecord.status` spells each one aloud.

Related, and deliberately NOT merged with this
==============================================

* :class:`~orpheus.sn.solver.ConvergenceCertificateError` asserts a
  *different* proposition — *"a convergence claim was made and it was
  FALSE"* (the in-M lag-death class).  It stays a hard error: a false claim
  is a bug, whereas an honest best-effort answer is a legitimate result.
* :class:`~orpheus.numerics.green_operator.ConvergenceFailure` is the
  *raising* form of this module's proposition, scoped to ``GreenOperator``'s
  own promise.  It keeps its Green-operator provenance; a caller who wants
  hard failure at a solver entry escalates the warning instead.
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Iterator
from dataclasses import dataclass, replace

__all__ = [
    "ESCALATION_FLAG",
    "ConvergenceWarning",
    "IterationBudget",
    "IterationRecord",
    "StoppingCriterion",
    "default_iteration_budget",
    "resolve_iteration_budget",
    "warn_if_unconverged",
]


class ConvergenceWarning(RuntimeWarning):
    r"""An iterative solve exhausted its budget; the answer is best-effort.

    Emitted by a public solver entry when the returned result carries
    ``converged = False`` — i.e. the loop hit its iteration cap with the
    stopping criterion unmet.  The result is still returned (it is the best
    iterate available), which is why this is a warning and not an exception;
    see the module docstring for the ERR-053 precedent behind that choice.

    The message states the budget that ran out, the tolerance that was not
    reached, and the last residual, so the reader can tell *how far* from
    converged the answer is — the distance between "one more sweep" and
    "diverging" is the whole diagnostic content.

    Escalate to a hard failure with :data:`ESCALATION_FLAG` (the category
    must be DOTTED — see the module docstring for why the short spelling
    silently was not a gate at all).
    """


#: The published CI escalation recipe, as a VALUE rather than a spelling.
#:
#: ``-W`` resolves an undotted category against ``builtins``, so the short
#: form ``error::ConvergenceWarning`` raises ``AttributeError`` and pytest
#: collects nothing.  Deriving the path from the class instead of retyping
#: it means a module move or a class rename cannot leave the documented
#: recipe pointing at something that does not exist — the failure mode this
#: constant exists to make unspellable (#340, 2026-08-09: the wrong spelling
#: had propagated by copy to four sites, one of them the emitted message).
#:
#: Gated by ``test_the_published_escalation_flag_actually_parses``, which
#: consumes this STRING through pytest's own parser rather than installing
#: the filter through the in-process API — the two are different claims, and
#: only the first one can fail in CI.
ESCALATION_FLAG = (
    f"-W error::{ConvergenceWarning.__module__}."
    f"{ConvergenceWarning.__qualname__}"
)


#: Fraction of a trajectory, measured from its END, used to estimate the rate.
#:
#: An iteration's EARLY residuals are dominated by the transient — whatever
#: content the initial guess happened to carry in the sub-dominant modes —
#: while the number that answers *"how many more sweeps?"* is the ASYMPTOTIC
#: rate, i.e. the dominant eigenvalue of the iteration operator, which only
#: the tail sees.  A half is the usual compromise: enough points to average
#: out per-step noise, few enough that the transient does not tilt the slope.
#: The fit always keeps at least two points, since a rate needs two.
_RATE_FIT_TAIL_FRACTION = 0.5


#: The slowest per-iteration contraction the DEFAULT budget promises to serve.
#:
#: `[M]` 2026-08-09 the representative worst case in the suite is the d=3
#: all-reflective box at ``rho = 0.9854`` (the undamped DD face sawtooth —
#: see ``derivations/sn_dd_face_transmission.py`` for why diamond alone
#: carries it).  ``0.986`` covers that with a hair of margin.
#:
#: ⚠ It deliberately does **NOT** cover ``Sigma_t/4`` (``rho = 0.99575``,
#: needing 5171 sweeps at ``1e-12``).  That is a stated limit, not an
#: oversight: serving ``0.996`` would make a genuinely *stuck* 1-D slab churn
#: ~6900 sweeps before admitting defeat, and the cost of a too-small budget
#: is now one warning that names the number to set — while the cost of a
#: too-large one is paid by every solve that was never going to converge.
_SERVED_RATE = 0.986


def _budget_from_law(
    *, log_initial: float, log_rate: float, tolerance: float
) -> int:
    r"""The geometric budget law, in ONE place.

    Smallest ``N`` such that the ``N``-th entry of ``r_k = e^{a} e^{bk}``
    lies strictly below ``tolerance``:

    .. math::

       N = \left\lfloor \frac{\ln \mathrm{tol} - a}{b} \right\rfloor + 2

    The ``+2`` rather than a ``ceil`` is not a fudge: it is what makes the
    expression correct with no branch when the crossing lands exactly on an
    integer (``ceil`` would return the index whose value EQUALS the
    tolerance, and a stopping test is strict).

    Both users of the law route through here —
    :meth:`StoppingCriterion.projected_iterations` fits ``(a, b)`` from an
    observed trajectory, and
    :func:`default_iteration_budget` evaluates the same law a priori at the
    served rate.  They are the posterior and the prior of one statement, and
    a second spelling of it is a twin waiting to drift.
    """
    return max(1, math.floor((math.log(tolerance) - log_initial) / log_rate) + 2)


def default_iteration_budget(
    tolerance: float, *, served_rate: float = _SERVED_RATE
) -> int:
    r"""How many iterations a tolerance needs, at the rate we promise to serve.

    ⭐ The point is not that this number is *right* — no single number can
    be, since the required count depends on a spectral radius the caller
    does not know.  The point is that it is **derived from a stated
    promise** instead of chosen, so it moves coherently with the tolerance
    and can be argued with.

    `[M]` 2026-08-09 the tree shipped five hardcoded budgets in SN plus a
    sixth in :class:`~orpheus.numerics.iteration.KEigenvalue`, in three
    ``(budget, tolerance)`` combinations.  The two SN families differed by
    **5x** where the law says the factor between ``1e-8`` and ``1e-12``
    should be ``ln(1e-12)/ln(1e-8) = 1.5``, and BOTH were **short** at d=3
    zero-leakage: 830 needed at ``1e-8`` against 200 shipped, 1441 at
    ``1e-12`` against 1000.  A constant cannot track a tolerance it does not
    know about.

    `[M]` 2026-08-09, measured end-to-end on the very configuration that
    opened this module — the d=3 all-reflective pure absorber, whose
    truncated exit is the founding defect described at the top of the file.
    Each row is a real converged solve, so ``NEEDED`` is observed, not fitted:

    =========  ======  =========  =========
    inner_tol  NEEDED  shipped    derived
    =========  ======  =========  =========
    1e-9         1007  1000  ✗       1471  ✓
    1e-12        1473  1000  ✗       1961  ✓
    1e-13        1631  1000  ✗       2125  ✓
    1e-15        2031  1000  ✗       2451  ✓
    =========  ======  =========  =========

    **The derived budget covers every row; the shipped constant covers
    none.**  Note how thin the first miss is — 1007 needed against 1000
    shipped, a shortfall of *seven sweeps*.  That is the margin on which the
    original gate read green for months, and it is the clearest statement of
    why a constant is the wrong shape here: nothing about ``1000`` knows
    that the answer moves by ~500 sweeps per four decades of tolerance.

    Parameters
    residual, and conservative whenever it is smaller (`[M]` the d=3 probe
    starts at ``3.7e-2``, so this over-budgets it by ~234 sweeps).  Erring
    long is the safe direction for a default precisely because the shortfall
    is what corrupts an answer, while the excess only costs time on a solve
    that was already failing — and :meth:`IterationRecord.projected_iterations`
    then tells that caller the number to set.

    Parameters
    ----------
    tolerance:
        The convergence tolerance the budget must serve.  Strictly positive
        and below one — a tolerance at or above the assumed unit initial
        residual is already met before the first iteration.
    served_rate:
        Contraction factor to size against; defaults to :data:`_SERVED_RATE`.
        Pass a measured ``rho`` to size a budget for a KNOWN problem.
    """
    if not 0.0 < tolerance < 1.0:
        raise ValueError(
            f"tolerance must lie in (0, 1) — a budget is sized against a "
            f"unit initial residual; got {tolerance!r}"
        )
    if not 0.0 < served_rate < 1.0:
        raise ValueError(
            f"served_rate must lie in (0, 1) to contract; got {served_rate!r}"
        )
    return _budget_from_law(
        log_initial=0.0, log_rate=math.log(served_rate), tolerance=tolerance
    )


def resolve_iteration_budget(max_iter: int | None, tolerance: float) -> int:
    """``None`` means *derive from the tolerance*; an int is a deliberate cap.

    Every entry point that takes an iteration budget performs exactly this
    resolution, so it is spelled once here rather than six times across two
    packages.  That matters beyond tidiness: a behavioural ``None`` default
    encodes an assumption about its surroundings, and an assumption restated
    at N call sites is an assumption that will differ at one of them
    (``[[lessons-L19]]``).  Here the assumption has a single home and a
    single docstring.

    Passing an explicit int is not deprecated and never will be — a caller
    who KNOWS their spectral radius, or who is deliberately starving a solve
    to measure its truncated exit, is exercising the API correctly.
    """
    return default_iteration_budget(tolerance) if max_iter is None else int(max_iter)


@dataclass(frozen=True)
class StoppingCriterion:
    r"""One quantity an iteration drives below one tolerance.

    A stopping test is a conjunction of these — SN's k-eigenvalue outer stops
    on ``dk < keff_tol`` **and** ``dphi < flux_tol``, and "which of the two is
    lagging?" is a question users ask and a flat ``converged`` cannot answer.
    Naming each one separately is what makes it answerable
    (:attr:`IterationRecord.binding_criterion`).

    The trajectory is the per-iteration history of the quantity, co-indexed
    with its siblings on the same :class:`IterationRecord`.  It holds
    MAGNITUDES (``|dk|``, ``||dphi||``, ``||r||``) — a negative entry is a
    producer bug and is refused at construction.  ``nan`` and ``inf`` are
    ACCEPTED: a diverging solve is a real state that must be recordable, and
    since ``nan < tol`` is ``False`` such a criterion correctly never clears.

    Attributes
    ----------
    name:
        The quantity's name in the caller's own vocabulary — ``"dk"``,
        ``"dphi"``, ``"residual"``.  It appears verbatim in
        :meth:`IterationRecord.report`, so it should read like the domain.
    trajectory:
        One entry per iteration, oldest first.  Empty means the level never
        iterated (see :attr:`cleared` for what that implies).
    tolerance:
        The threshold this quantity was judged against.  Non-negative.

        ``0.0`` is legal and means *unsatisfiable by a strict test* — a real
        production input, not a degenerate one: ``GreenOperator(tol=0)`` is
        how the unsatisfiability guard is exercised, and GMRES's
        exact-breakdown path reaches a literal ``0.0`` residual.  `[M]`
        2026-08-09 this class first shipped demanding a strictly positive
        tolerance and two production paths refused to build a record — the
        vv anti-pattern #16 shape, a type asserting more than its callers
        actually promise.  Nothing about the comparison changed: ``cleared``
        stays strict, so a zero tolerance never clears, exactly as the
        retired ``_claims_convergence`` behaved.
    """

    name: str
    trajectory: tuple[float, ...]
    tolerance: float

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("a stopping criterion must be named")
        if self.tolerance < 0.0 or math.isnan(self.tolerance):
            raise ValueError(
                f"{self.name}: tolerance must be non-negative, got "
                f"{self.tolerance!r}"
            )
        # `v < 0.0` is False for nan, so this admits nan/inf by construction.
        negative = [v for v in self.trajectory if v < 0.0]
        if negative:
            raise ValueError(
                f"{self.name}: a stopping criterion holds MAGNITUDES, so a "
                f"negative entry is a producer bug; got {negative[0]!r}"
            )

    # ── building one up, one iteration at a time ─────────────────────────

    @classmethod
    def reading(
        cls, name: str, value: float, tolerance: float,
    ) -> "StoppingCriterion":
        """ONE iteration's reading of this quantity — a trajectory of length 1.

        A stopping test is evaluated per iteration, but a criterion is a
        *trajectory*; this is the named bridge between the two, so a producer
        that measures ``|Δk|`` at one iterate says exactly that and the loop
        accumulates with :meth:`extended_with`.  Spelling the same thing as
        ``StoppingCriterion(name, (value,), tol)`` works and means the same —
        the named constructor exists so the per-iteration convention is
        greppable and so the realizers of
        :meth:`~orpheus.numerics.eigenvalue.EigenvalueSolver.measure_stopping_criteria`
        read like the domain::

            return (
                StoppingCriterion.reading("dk", abs(keff - keff_old), self.keff_tol),
                StoppingCriterion.reading("dphi", dphi, self.flux_tol),
            )
        """
        return cls(name=name, trajectory=(value,), tolerance=tolerance)

    def extended_with(self, reading: "StoppingCriterion") -> "StoppingCriterion":
        """This criterion's trajectory extended by a later reading of itself.

        Routes through :func:`~dataclasses.replace`, so ``__post_init__``
        re-fires on the result and the magnitude law is re-established for
        free rather than restated here (``coding-elegance`` Pattern 4 ∩
        Pattern 2 — an invariant with two spellings is an invariant that will
        drift).

        Refuses a reading of a DIFFERENT quantity, and refuses a tolerance
        that moved mid-solve: a trajectory is only comparable against one
        threshold, so a shifting tolerance silently invalidates
        :attr:`cleared`, :attr:`distance` and every projection derived from
        them.  Neither is reachable from the shipped loop; both are cheap to
        refuse and would otherwise be silent.
        """
        if reading.name != self.name:
            raise ValueError(
                f"cannot extend criterion {self.name!r} with a reading of "
                f"{reading.name!r} — trajectories are per-quantity"
            )
        if reading.tolerance != self.tolerance:
            raise ValueError(
                f"{self.name}: tolerance moved mid-solve "
                f"({self.tolerance!r} → {reading.tolerance!r}) — the whole "
                f"trajectory would no longer be judged against one threshold"
            )
        return replace(self, trajectory=self.trajectory + reading.trajectory)

    # ── what happened ────────────────────────────────────────────────────

    @property
    def n_iterations(self) -> int:
        """How many iterations this quantity was measured over."""
        return len(self.trajectory)

    @property
    def last(self) -> float | None:
        """The final value, or ``None`` if the level never iterated."""
        return self.trajectory[-1] if self.trajectory else None

    @property
    def cleared(self) -> bool:
        """Did the final value fall below the tolerance?

        Vacuously ``True`` for an empty trajectory, following ``all(())``:
        a criterion that was never measured is not one that FAILED, and
        conflating the two is precisely the misreading that turned "GMRES
        returned on the initial guess" into 44 phantom truncations.  A
        producer that considers zero iterations suspicious says so with
        :attr:`IterationRecord.min_iterations`, which is a statement about
        the LOOP, not about this quantity.
        """
        last = self.last
        return last is None or last < self.tolerance

    @property
    def distance(self) -> float:
        """How far the final value sits from clearing, as a ratio to tol.

        Below one means cleared; the LARGEST ratio across a record's criteria
        is the one that bound (whether it cleared last or failed worst), which
        is what makes one definition serve both cases with no branch.  Zero
        for an empty trajectory — nothing was measured, so nothing is blocking.

        Infinite against a zero tolerance, because a strict test can never
        clear one — which keeps ``cleared`` and ``distance < 1`` two
        spellings of one statement in every case, including that one.
        """
        last = self.last
        if last is None:
            return 0.0
        if self.tolerance == 0.0:
            return math.inf
        return last / self.tolerance

    # ── what it would take ───────────────────────────────────────────────

    def _log_fit(self) -> tuple[float, float] | None:
        r"""Least-squares fit of :math:`\ln r_k = a + b k` over the tail.

        Returns ``(a, b)`` with ``k`` indexed against the FULL trajectory, so
        the intercept is the fitted value at iteration zero and
        :meth:`projected_iterations` speaks in absolute iteration counts.

        ``None`` when fewer than two usable points remain — usable meaning
        finite and strictly positive, since a rate is a ratio of successive
        residuals and ``ln 0`` is not a number.
        """
        n = len(self.trajectory)
        if n < 2:
            return None
        first = min(int(n * (1.0 - _RATE_FIT_TAIL_FRACTION)), n - 2)
        points = [
            (float(k), math.log(v))
            for k, v in enumerate(self.trajectory)
            if k >= first and v > 0.0 and math.isfinite(v)
        ]
        if len(points) < 2:
            return None

        mean_index = sum(k for k, _ in points) / len(points)
        mean_log_residual = sum(y for _, y in points) / len(points)
        covariance = sum(
            (k - mean_index) * (y - mean_log_residual) for k, y in points
        )
        variance = sum((k - mean_index) ** 2 for k, _ in points)
        if variance == 0.0:
            return None

        log_rate = covariance / variance
        return mean_log_residual - log_rate * mean_index, log_rate

    @property
    def rate(self) -> float | None:
        r"""The observed geometric decay rate :math:`\rho` over the tail.

        ``rho < 1`` contracts, ``rho >= 1`` does not.  ``None`` when the
        trajectory is too short or too degenerate to support an estimate —
        which is the honest answer, not a default of 1.
        """
        fit = self._log_fit()
        return None if fit is None else math.exp(fit[1])

    def projected_iterations(self, tolerance: float | None = None) -> int | None:
        r"""Iterations the observed rate says are needed to clear ``tolerance``.

        This is the number that turns *"raise the budget"* into *"set
        ``max_inner=849``"*, and it is counted from iteration zero — so it is
        directly comparable to :attr:`IterationBudget.in_iterations`.

        ⛔ **NOT to :attr:`IterationBudget.limit`.**  This sentence read
        *"directly comparable to ``IterationRecord.budget``"* until
        2026-08-13, and it was the sentence the actionable advice rested on:
        the count here is in TRAJECTORY units, while the limit is in the
        caller's KNOB units, and the two coincide only when
        :attr:`IterationBudget.iterations_per_unit` is 1.  On the GMRES arm
        it is ``restart`` = ``n_dof``, so the advice printed an Arnoldi-step
        count against a restart-cycle knob (#349).  Route the conversion
        through :meth:`IterationBudget.covering`, which is the identity
        wherever the old sentence was true.

        Derived by extrapolating the fitted line to the crossing:

        .. math::

           \ln r_k = a + b k, \qquad
           N = \left\lfloor \frac{\ln \mathrm{tol} - a}{b} \right\rfloor + 2

        ⭐ The **intercept is fitted, not assumed**.  Reading the law as
        :math:`n = |\ln(\mathrm{tol}/r_0)| / |\ln\rho|` — with :math:`r_0` the
        first residual — treats the initial transient as if it lay on the
        asymptotic line.  `[M]` 2026-08-09 on the d=3 all-reflective control
        that offset is **403 of 1473 iterations, i.e. 27 %**, and a rate-only
        comparison of two splittings flipped SIGN on 4 of 5 configurations
        when re-measured as sweep counts.  A spectral radius predicts a RATE;
        only the fit predicts a COST.

        Returns ``None`` when no rate is estimable, and — deliberately — also
        when the trajectory is NOT decaying.  At :math:`\rho \ge 1` no budget
        suffices, and answering with a large number would be a lie: `[M]` the
        campaign found a mutated heterogeneous slab whose NEGATIVE dominant
        eigenvalue makes the increment sign-alternate, so its stop test is
        unsatisfiable *forever* and raising the budget is futile.  That is a
        different failure from a shortfall and must not be reported as one.

        Parameters
        ----------
        tolerance:
            Target to project against; defaults to this criterion's own.
            Passing a looser one answers *"what would I get for my budget?"*.
        """
        fit = self._log_fit()
        if fit is None:
            return None
        intercept, log_rate = fit
        if log_rate >= 0.0:
            return None
        if (self.tolerance if tolerance is None else tolerance) <= 0.0:
            return None  # no finite count reaches an exactly-zero residual
        return _budget_from_law(
            log_initial=intercept,
            log_rate=log_rate,
            tolerance=self.tolerance if tolerance is None else tolerance,
        )


@dataclass(frozen=True)
class IterationBudget:
    r"""A ceiling, the knob that set it, and **what one unit of it buys**.

    A budget and an iteration count are comparable only when they count the
    same thing.  For five of the tree's six producers they trivially do — a
    :class:`~orpheus.numerics.iteration.SourceIteration` pass, a
    power-iteration outer, a CP within-group pass, an MoC transport sweep are
    each one unit of their own ``max_iter``.  For
    :class:`~orpheus.numerics.iteration.KrylovAcceleration` they do **not**,
    and the gap is not a detail of bookkeeping — it is a factor of ``n_dof``.

    The GMRES asymmetry, measured
    =============================

    scipy's ``gmres(maxiter=m, restart=r)`` runs an OUTER loop of ``m``
    **restart cycles**, each an INNER Arnoldi loop of up to ``r`` steps, and
    under ``callback_type='pr_norm'`` the callback fires once per **inner**
    step.  So the trajectory this record carries is in Arnoldi steps while
    the caller's knob is in restart cycles, and one cycle buys ``r`` of them.

    `[M]` 2026-08-13, scipy 1.17.1, a dense ``60x60`` system, 12 of 12
    non-converging rows over ``m ∈ {1,2,3,5} x r ∈ {4,7,20}``:

    ===========  ===========  =========  =======
    ``maxiter``  ``restart``  callbacks  ``m*r``
    ===========  ===========  =========  =======
    1            4            4          4
    2            7            14         14
    3            20           60         60
    5            20           100        100
    ===========  ===========  =========  =======

    Exact in every row, with ``info == maxiter``; the two CONVERGING rows
    fall short (8 callbacks of 100, 8 of 300, ``info == 0``).  So the
    exchange rate is ``r``, exactly, and "the loop ran out" is
    ``callbacks >= m*r``.

    ⛔ **The defect this type exists to make unspellable** (#349, `[M]`
    reproduced at ``7b1618ab``): ``exhausted_budget`` was
    ``n_iterations >= budget`` over the raw pair, so a 1-group ``nx=16``
    GL(8) vacuum slab solved with ``inner_solver="krylov", max_inner=5``
    returned ``budget=5, n_iterations=91, status=CONVERGED`` — and
    ``exhausted_budget=True``.  The solve was healthy (its ``pr_residual``
    criterion cleared at ``9.4465e-11 < 1e-10``, zero warnings); the
    *property meaning "it stopped because it ran out"* was simply reading
    91 Arnoldi steps against 5 restart cycles.  Raising the knob 10x to
    ``max_inner=50`` did not fix it, because the two numbers were never
    commensurable.

    Why a TYPE and not a third int on the record
    ============================================

    ``.claude/rules/coding-standards.md`` makes this decidable: mint a type
    iff there are ≥2 non-isomorphic realizations of the concept AND a
    non-identity morphism is applied to it.  Both hold — restart cycles and
    fixed-point passes are the two realizations, and cycles → Arnoldi steps
    is the non-identity morphism, applied in :meth:`exhausted_by` and
    **inverted** in :meth:`covering`.  That inverse is the load-bearing half:
    the actionable advice (#340 N3) projects a count off the trajectory, in
    Arnoldi steps, and must hand the reader a number to type into a knob
    that takes cycles.  Before this type that conversion did not exist
    anywhere — an un-named welded operation, which is the failure mode
    Cardinal Rule 2 names.

    ⚠ What this does NOT do: it does not make a WRONG
    :attr:`iterations_per_unit` unspellable.  A producer can still pass the
    default ``1`` when its loop buys more, and no invariant here can catch
    that — only the producer knows.  What the type does is force the
    producer to have an opinion, and give the conversion a name that a
    reader can check.  The honest claim is *"unspellable without saying so"*,
    not *"impossible"*.

    ⭐ Note the deliberate NON-relationship to :func:`default_iteration_budget`
    despite the shared word.  That function answers *"how many ITERATIONS
    does this tolerance need?"* from the geometric rate law, so its output is
    in iteration units and is a correct :attr:`limit` only when
    :attr:`iterations_per_unit` is 1.  `[M]` on the Krylov path it is not:
    ``sn/solver.py:1174`` resolves ``max_inner=None`` through it and hands the
    result to scipy's ``maxiter``, so a derived budget of ~1471 becomes 1471
    *cycles* of ``n_dof`` steps each.  That over-provisions, which is the safe
    direction and is left alone here — but it is the same confusion, and it
    is why ``max_inner`` effectively never binds a Krylov inner.

    Attributes
    ----------
    limit:
        The ceiling, **in the caller's knob units** — the number the caller
        actually passed, so a diagnostic can quote it back verbatim.  Zero
        means *no budget*, which is how a direct solve (an LU resolvent) says
        it cannot be truncated; :attr:`exhausted_by` is then always ``False``.
    name:
        What the CALLER's API calls the knob that set :attr:`limit` —
        ``"max_inner"``, ``"max_outer"``, or the default ``"max_iter"``.
        Absorbed here from the record's retired ``budget_name`` field,
        because a limit and the knob that spells it are one fact: advice
        that says *"set X=N"* needs both halves or neither.

        `[M]` the shipped spellings are exactly
        ``{max_inner, max_outer, max_iter, max_inner_sweeps, params.max_inner}``.
        The field is a free ``str`` and the ONLY illegal value is the empty
        one — a whitelist would refuse
        :class:`~orpheus.numerics.green_operator.GreenOperator`, whose own
        knob IS ``max_iter`` (``green_operator.py:277``), the same shape as
        the ``GreenOperator(tol=0)`` refusal this campaign already had to
        back out (vv anti-pattern #16, in the direction that breaks
        production).  A level always has *some* knob, and advice with no
        token to type is not advice.

        ⭐ It is stated by the producer, not inferred by the reader, because
        advice needs a token the caller can actually type.  A diagnostic that
        says *"set max_outer=838"* about a starved INNER sends the reader to
        the wrong knob — and worse than wrong, because with a starved inner
        raising ``max_outer`` cannot help at all (the outer's stop test is
        entirely increments, which the starved inner suppresses; #340 F2).

        ⚠ It is NOT derivable from :attr:`IterationRecord.label`.  That
        string is chosen for humans, and reading a control decision off it is
        the stringly-typed dispatch
        :attr:`orpheus.sn.solution.IterationHistory._is_outer` already
        refuses.  Nor is it derivable from depth: the same
        :class:`~orpheus.numerics.iteration.SourceIteration` is the *inner*
        of an eigenvalue solve and the *whole* of a fixed-source one.  Only
        the site that CONSTRUCTED the level knows which of its own parameters
        supplied the cap, so that site states it.
    iterations_per_unit:
        How many RECORDED iterations one unit of :attr:`limit` buys.  ``1``
        (the default) is the identity conversion and the honest answer for
        every producer whose loop counter *is* its knob.  GMRES states its
        ``restart`` here.

        ⚠ It is not derivable and never will be.  Only the site that drove
        the loop knows how its library spends a budget unit — the same
        argument that puts :attr:`IterationRecord.iterations_run` on the
        producer rather than on the record.
    """

    limit: int = 0
    name: str = "max_iter"
    iterations_per_unit: int = 1

    def __post_init__(self) -> None:
        if self.limit < 0:
            raise ValueError(f"budget limit must be >= 0, got {self.limit}")
        if not self.name:
            raise ValueError(
                "budget name must name the caller's knob — it is what the "
                "advice tells the reader to set, and an empty one leaves "
                "nothing to type"
            )
        if self.iterations_per_unit < 1:
            raise ValueError(
                f"iterations_per_unit must be >= 1 — one unit of a budget "
                f"buys at least one iteration, or the budget cannot bind; "
                f"got {self.iterations_per_unit}"
            )

    @property
    def is_budgeted(self) -> bool:
        """Is there a cap at all?

        ``limit == 0`` is how a DIRECT level says it cannot be truncated — an
        LU back-substitution, a closed form.  The question is asked by
        :meth:`exhausted_by` and by every consumer that skips unbudgeted
        levels before advising a knob, so it is named once here rather than
        spelled ``limit > 0`` at each of them.
        """
        return self.limit > 0

    @property
    def in_iterations(self) -> int:
        """The ceiling in TRAJECTORY units — the only form a count may meet.

        This is the single site where the exchange rate is spent.  Every
        comparison against an iteration count goes through here, so there is
        no second spelling of the conversion to drift (Pattern 2).
        """
        return self.limit * self.iterations_per_unit

    def exhausted_by(self, n_iterations: int) -> bool:
        """Did a loop that ran ``n_iterations`` stop because it ran out?

        ``limit == 0`` means *unbudgeted*, and an unbudgeted loop can never
        be truncated however long it ran — that is what a ``DIRECT`` level
        (:class:`~orpheus.diffusion.solver` 's LU resolvent) relies on.
        """
        return self.is_budgeted and n_iterations >= self.in_iterations

    def covering(self, n_iterations: int) -> int:
        """The smallest setting of :attr:`name` that buys ``n_iterations``.

        The inverse of :attr:`in_iterations`, and the number that belongs in
        *"set {name}={this}"*.  Ceiling division, because a partial budget
        unit cannot be bought: 800 Arnoldi steps at 128 per cycle needs
        **7** cycles, not 6.25 and not 800.

        Rounds up to at least 1 whenever any iterations are wanted, so the
        advice never tells a reader to set a knob to zero.
        """
        if n_iterations <= 0:
            return 0
        return max(1, -(-n_iterations // self.iterations_per_unit))

    def __str__(self) -> str:
        """``max_inner=5`` — or ``max_inner=5 (x128 = 640 iterations)``.

        The conversion is printed ONLY when it is not the identity, so the
        five unit-consistent producers read exactly as they did before this
        type existed and the one asymmetric producer declares itself.
        """
        head = f"{self.name}={self.limit}"
        if self.iterations_per_unit == 1:
            return head
        return (
            f"{head} (x{self.iterations_per_unit} = "
            f"{self.in_iterations} iterations)"
        )


@dataclass(frozen=True)
class IterationRecord:
    r"""One LEVEL of an iterative solve: what it wanted, what it got, why it
    stopped — and, recursively, the same for every level beneath it.

    An SN k-eigenvalue solve is not one loop but a tree of them, and every
    diagnostic question worth asking is about a NODE of that tree: *which*
    level stalled, on *which* criterion, needing *what* budget.  A flat record
    answers none of them, and its ``converged`` silently comes to mean "the
    outermost one" — which is how a solve whose inner starved reports success
    (see the module docstring's measured 5.3x ``keff_tol``).

    ⭐ **Nothing here is stored; every verdict is derived.**  There is no
    ``converged`` field to set, so there is no ``converged=True`` to hardcode
    (#342) and no sixth spelling to drift (``... and iteration > 5``).  This
    is ``coding-elegance`` Pattern 4 applied to a boolean: the illegal state
    is not rejected at runtime, it is unrepresentable.

    Attributes
    ----------
    label:
        This level's name in the domain's vocabulary — ``"outer(power)"``,
        ``"inner(within-group g=1)"``, ``"krylov(gmres)"``.  It is what the
        reader sees first in :meth:`report`, so it should say which loop this
        is, not which function implements it.
    criteria:
        The conjunction this level stops on.  All co-indexed: one entry per
        criterion per iteration, enforced at construction, because they are
        measurements of the SAME iterations and a length mismatch means a
        producer dropped one.
    budget:
        The cap this level was given — the limit, the NAME of the knob that
        set it, and what one unit of that knob buys — so
        :attr:`exhausted_budget` is answerable without the caller having to
        remember what it passed.  See :class:`IterationBudget`, which carries
        the whole fact including the exchange rate that makes the comparison
        dimensionally sound (#349).

        ⛔ Until 2026-08-13 this was a bare ``int`` beside a separate
        ``budget_name: str``, and ``exhausted_budget`` compared it directly
        against :attr:`n_iterations`.  That is sound for five of the six
        producers and wrong for GMRES by a factor of ``restart``; the
        measurement and the reasoning are on :class:`IterationBudget`.  Both
        fields are absorbed there — a limit, the knob that spells it, and the
        rate that makes it comparable are ONE fact, and splitting them across
        loose ints is what let the mismatch ship.
    iterations_run:
        How many times the loop actually iterated, when that differs from
        the number of criterion measurements.  ``None`` (the default) means
        *the same*, which is the common case.

        ⭐ This is the ONE thing the record cannot derive, because the offset
        is a **per-producer fact**.  `[M]` 2026-08-09, both conventions ship
        in this tree: :class:`~orpheus.numerics.iteration.SourceIteration`
        measures the *difference* between successive iterates, so ``P``
        passes yield ``P - 1`` residuals and a run that exhausts
        ``max_iter=50`` records 49;
        :class:`~orpheus.numerics.iteration.KrylovAcceleration` gets one
        callback per inner Arnoldi step, so its count matches its own
        TRAJECTORY.  Inferring one rule would silently mis-read the other — which
        is exactly what happened before the record existed, where the same
        ``n_inner`` field was written as ``len(residuals)`` by one driver and
        ``len(residuals) + 1`` by the other, undocumented.

        ⛔ That last clause read *"so its counts match"* until 2026-08-13, and
        the missing word cost a merge cycle.  It is true of *count vs
        trajectory* — the only pair this attribute is about — and it READS as
        true of *count vs budget*, which is the pair it says nothing about
        and which for GMRES is off by ``restart``.  It is the sentence most
        likely to talk a reader out of suspecting #349.  The budget half of
        the question lives on :class:`IterationBudget`; this attribute is
        silent on it by construction, and now says so.

        It is an OBSERVATION, not a verdict.  The distinction this class
        rests on is that no *convergence judgement* may be stored; measured
        data plainly must be.  The construction invariant keeps the two
        consistent: you cannot measure a criterion more often than you
        iterated.
    min_iterations:
        Iterations below which this level refuses to claim convergence.  The
        home for a guard like SN's ``iteration <= 2``, and the honest way for
        a producer to say "zero iterations is not success here" without
        distorting what :attr:`StoppingCriterion.cleared` means.
    children:
        The levels nested inside one iteration of this one.
    """

    label: str
    criteria: tuple[StoppingCriterion, ...] = ()
    budget: IterationBudget = IterationBudget()
    iterations_run: int | None = None
    min_iterations: int = 0
    children: tuple[IterationRecord, ...] = ()

    def __post_init__(self) -> None:
        if not self.label:
            raise ValueError("an iteration record must be labelled")
        # ⭐ The budget's own invariants (non-negative limit, non-empty knob
        # name, a conversion that buys at least one iteration) live on
        # :class:`IterationBudget.__post_init__`, not here.  They moved with
        # the concept: a record is not the only thing that can hold a budget,
        # and restating them would be the Pattern-2 duplicate that
        # `coding-elegance` warns produces two spellings of one law.
        if self.min_iterations < 0:
            raise ValueError(
                f"{self.label}: min_iterations must be >= 0, got "
                f"{self.min_iterations}"
            )
        names = [criterion.name for criterion in self.criteria]
        if len(set(names)) != len(names):
            raise ValueError(
                f"{self.label}: criterion names must be unique, got {names}"
            )
        lengths = {criterion.n_iterations for criterion in self.criteria}
        if len(lengths) > 1:
            raise ValueError(
                f"{self.label}: criteria are co-indexed by iteration, so their "
                f"trajectories must be the same length; got "
                + ", ".join(
                    f"{c.name}={c.n_iterations}" for c in self.criteria
                )
            )
        if self.iterations_run is not None:
            measured = max(lengths, default=0)
            if self.iterations_run < measured:
                raise ValueError(
                    f"{self.label}: iterations_run={self.iterations_run} is "
                    f"fewer than the {measured} criterion measurements — a "
                    f"loop cannot measure more often than it iterates"
                )

    # ── the four states a loop can stop in ───────────────────────────────

    @property
    def n_iterations(self) -> int:
        """Iterations this level ran (zero if it never entered its loop).

        The producer's own count when it stated one, else the number of
        criterion measurements — see :attr:`iterations_run` for why only the
        producer can know the difference.
        """
        if self.iterations_run is not None:
            return self.iterations_run
        return self.criteria[0].n_iterations if self.criteria else 0

    @property
    def iterated(self) -> bool:
        """Did this level iterate at all?

        A direct solve — a factorisation, a closed form, a Krylov call that
        was satisfied by its initial guess — records ``False`` here and
        ``True`` for :attr:`converged`.  The two are different questions and
        collapsing them is what manufactures phantom truncations.
        """
        return self.n_iterations > 0

    @property
    def converged(self) -> bool:
        """Did THIS level meet all of its own criteria?

        Says nothing about the levels beneath it — that is
        :attr:`fully_converged`, and keeping them separate is what lets a
        caller tell "the outer stalled because an inner starved" from "the
        outer is genuinely non-critical".

        ⭐ **A level that RAN and measured nothing cannot claim convergence.**
        :attr:`StoppingCriterion.cleared` is vacuously ``True`` on an empty
        trajectory, which is the right reading when the loop never entered —
        GMRES returning on its initial guess *did* converge, and calling
        that a truncation is the misreading that produced 44 phantom rows in
        #340's audit.  It is the WRONG reading when the loop iterated and
        simply never got to measure: a `SourceIteration` given
        ``max_iter=1`` makes one pass and records no residual, because its
        stop compares SUCCESSIVE iterates.  There is no evidence there, and
        absence of evidence must not read as convergence.

        :attr:`iterated` is exactly the discriminator between those two, so
        the rule needs no special case beyond naming it.

        ⛔ That rule first shipped as ``any(criterion.n_iterations == 0 …)``,
        which is silent on the case of **no criteria at all** — ``any(())`` is
        ``False``, so a level that iterated while declaring nothing to measure
        fell through to ``all(())`` and claimed convergence.  That is the same
        vacuous-conjunction lie the power-iteration loop
        (:func:`~orpheus.numerics.eigenvalue.power_iteration`) refuses at the
        producer, left open one layer down;
        `[M]` 2026-08-09 it was found by writing the neighbouring test, not by
        review.  Since criteria are co-indexed the two spellings agree
        wherever there IS a criterion, so widening to "no criterion was
        measured" loses nothing and closes the hole.
        """
        if self.n_iterations < self.min_iterations:
            return False
        if self.iterated and not any(
            criterion.n_iterations > 0 for criterion in self.criteria
        ):
            return False
        return all(criterion.cleared for criterion in self.criteria)

    @property
    def exhausted_budget(self) -> bool:
        """Did this level stop because it ran out of iterations?

        The comparison is the budget's to make, not the record's: only
        :class:`IterationBudget` knows the exchange rate between the caller's
        knob and the trajectory this record counts (#349).
        """
        return self.budget.exhausted_by(self.n_iterations)

    @property
    def truncated(self) -> bool:
        """Budget gone, criterion unmet — the #340 defect, named.

        The returned iterate is mid-descent and arbitrary; any value asserted
        against it is green or red by luck.
        """
        return self.exhausted_budget and not self.converged

    @property
    def status(self) -> str:
        """The state, in one word, for a human: the four cases spelled aloud."""
        if not self.iterated:
            return "DIRECT"
        if self.converged:
            return "CONVERGED"
        if self.truncated:
            return "TRUNCATED"
        return "STOPPED"

    # ── the tree ─────────────────────────────────────────────────────────

    @property
    def fully_converged(self) -> bool:
        """Did this level AND every level beneath it converge?

        **The honest answer to "can I trust this number?"** — and the one a
        value gate must assert before asserting physics.
        """
        return self.converged and all(
            child.fully_converged for child in self.children
        )

    @property
    def first_failure(self) -> IterationRecord | None:
        """The deepest, earliest level that did not converge — the CAUSE.

        Children are searched before ``self``, so a starved inner is reported
        rather than the outer it starved.  This runs even when ``self``
        converged, because that is exactly the measured failure: an outer
        whose increment-only stop test is suppressed BY the starved inner
        reports success while carrying its error.
        """
        for child in self.children:
            found = child.first_failure
            if found is not None:
                return found
        return None if self.converged else self

    def walk(self) -> Iterator[IterationRecord]:
        """Every record in the tree, this level first, then depth-first."""
        yield self
        for child in self.children:
            yield from child.walk()

    # ── what bound, and what it would have taken ─────────────────────────

    @property
    def binding_criterion(self) -> StoppingCriterion | None:
        """The criterion furthest from clearing — the one that bound.

        If this level converged it is the one that cleared LAST (the
        rate-limiting one); if it did not, it is the one that failed worst.
        One ratio, :attr:`StoppingCriterion.distance`, answers both.
        """
        return max(self.criteria, key=lambda c: c.distance, default=None)

    @property
    def rate(self) -> float | None:
        """The binding criterion's observed decay rate."""
        binding = self.binding_criterion
        return None if binding is None else binding.rate

    def projected_iterations(self) -> int | None:
        """The budget that would have sufficed for the binding criterion.

        ``None`` when unknowable or when no budget suffices — see
        :meth:`StoppingCriterion.projected_iterations`.
        """
        binding = self.binding_criterion
        return None if binding is None else binding.projected_iterations()

    # ── the thing you paste into an issue ────────────────────────────────

    def report(self, *, _depth: int = 0) -> str:
        """The tree, human-readable, one level per indented block.

        Written to be pasted verbatim into a bug report: every number a
        reader needs to tell "one more sweep" from "diverging" is on the
        page, with no re-run.
        """
        pad = "  " * _depth
        # ⭐ Against ``in_iterations``, never the raw ``limit``.  The two differ
        # exactly when the knob's unit is not the trajectory's, and printing
        # the raw pair is what rendered a healthy GMRES leaf as the
        # self-contradicting ``CONVERGED (91/5 iterations)`` (#349).  The knob
        # is appended ONLY in that case, so every unit-consistent producer's
        # line is byte-identical to what it printed before this type existed —
        # including the ``(200/200 iterations)`` a committed gate pins.
        scale = (
            ""
            if self.budget.iterations_per_unit == 1
            else (
                f" [{self.budget.name}={self.budget.limit}"
                f" x {self.budget.iterations_per_unit}]"
            )
        )
        lines = [
            f"{pad}{self.label}: {self.status} "
            f"({self.n_iterations}/{self.budget.in_iterations} "
            f"iterations{scale})"
        ]
        binding = self.binding_criterion
        for criterion in self.criteria:
            tag = " <- binding" if criterion is binding else ""
            last = criterion.last
            # ⚠ An unmeasured criterion must NOT read "met" here.  `cleared`
            # is vacuously True on an empty trajectory and that is correct for
            # what it decides, but printing "met" beside a level whose own
            # status line says TRUNCATED gives a reader two statements that
            # cannot both be acted on.  The report's job is to be pasted into
            # a bug report unedited, so it says which of the two it is.
            if last is None:
                value, mark = "n/a", "not measured"
            else:
                value = f"{last:.3e}"
                mark = "met" if criterion.cleared else "MISSED"
            lines.append(
                f"{pad}  {criterion.name}: {value} vs tol "
                f"{criterion.tolerance:.3e}  {mark}{tag}"
            )
            if criterion is binding and not criterion.cleared:
                rate = criterion.rate
                needed = criterion.projected_iterations()
                if rate is None:
                    lines.append(f"{pad}    rate: not estimable")
                elif needed is None:
                    lines.append(
                        f"{pad}    rate: {rate:.6f} — NOT contracting; "
                        f"no budget suffices at this rate"
                    )
                else:
                    lines.append(
                        f"{pad}    rate: {rate:.6f} — needs ~{needed} "
                        f"iterations for tol {criterion.tolerance:.3e}"
                    )
        for child in self.children:
            lines.append(child.report(_depth=_depth + 1))
        return "\n".join(lines)


def warn_if_unconverged(
    record: IterationRecord,
    *,
    where: str,
    balance_defect: float | None = None,
) -> None:
    r"""Make a best-effort exit AUDIBLE at a public entry, in ANY family.

    The single ENTRY-LEVEL emission point for :class:`ConvergenceWarning`.
    Every public solver entry — :func:`~orpheus.sn.solver.solve_sn` and its
    four siblings, :func:`~orpheus.cp.solver.solve_cp`,
    :func:`~orpheus.moc.solver.solve_moc`,
    :func:`~orpheus.diffusion.solver.solve_diffusion_1d` — calls this after
    building its record, so a truncated solve announces itself once, in one
    voice, wherever it came from.

    ⛔ **There are TWO ConvergenceWarning emission points in this package and
    they must NOT be consolidated** (#340 N4.7, 2026-08-11).  This one is
    ENTRY-level: it speaks about an iteration TREE at the boundary where a
    caller receives a result.  The other,
    :meth:`~orpheus.numerics.iteration.KrylovAcceleration.solve`
    (``iteration.py``), is LEAF-level and mechanism-specific: it reports that
    SciPy's ``gmres`` returned ``info != 0`` for one inner solve.  A pass that
    "unifies the ConvergenceWarning emission points" would destroy the leaf's
    mechanism detail, or drown the entry's tree verdict in per-inner noise.
    Two levels, two warnings, both correct.

    ⭐ **Scope: "a truncated solve" means a solve ANY of whose levels did not
    converge** (#340 N6b, 2026-08-10).  The guard below is
    :attr:`IterationRecord.fully_converged`, so an outer that met its own
    criteria while standing on a TRUNCATED inner is AUDIBLE — which is the
    whole of the #340 headline defect (F1): the record always knew
    (``fully_converged`` was ``False``, ``first_failure`` named the inner) and
    nothing said so.

    ⛔ Until 2026-08-10 the guard was the TOP level only (``converged``) and
    `[M]` 20 tests in the shipped suite returned silently while standing on a
    starved inner.  Those 20 were adjudicated before the flip
    (`scratch/n6b_r2_adjudication.md`): 10 declare their truncation as the
    fixture and suppress this ONE category in-test, and 10 are now audible on
    purpose, tracked with measured budgets in
    `#352 <https://github.com/deOliveira-R/ORPHEUS/issues/352>`_.

    ⚠ The flip was originally scheduled to ride an N5 residual certificate
    that would separate a *corrupting* truncation from a *benign* one.  `[M]`
    that certificate was REFUTED by measurement: the raw defect's benign and
    corrupting populations overlap **634×** and it misses 15 of 16 corrupting
    cases, so it cannot gate.  The user's ruling (2026-08-10) was to widen the
    guard unconditionally and report the balance projection as a DIAGNOSTIC
    number rather than a threshold — a truncation the caller has not declared
    is worth saying out loud whether or not we can yet say how much it cost.

    Why a warning and not a raise: the ERR-053 / D-H.1e ruling (see this
    module's docstring) — legitimate callers harvest the residual history of a
    deliberately-truncated solve, and `[M]` an audit found zero production and
    zero ``examples/`` callers that depend on the answer of one.  Escalate in
    CI with :data:`ESCALATION_FLAG` — the DOTTED category; the short spelling
    does not parse, so it was never a gate (#340, 2026-08-09).

    ⭐ **It speaks for the level that FAILED, not for the entry that was
    called** (#340 N6, 2026-08-10).  ``where`` names the entry — the only
    thing the caller still supplies — and everything else is read off
    :attr:`IterationRecord.first_failure`, which searches CHILDREN FIRST.  So
    a starved inner is named, with ITS budget, ITS knob and ITS tolerance,
    even when the outer standing on it met its own criteria.  The three
    arguments this function used to take for those facts (``budget_name``,
    ``budget``, ``tol``) are retired: each entry was describing its own TOP
    level, so on a tree they described the wrong one, and the advice pointed
    at a knob that cannot help (F2).

    ⭐ **That retirement is also why there is no per-family "read this
    attribute" argument** (#340 N4.7, 2026-08-11).  Until this function left
    ``sn/solver.py`` its closing advice was the literal string
    ``solution.history.fully_converged`` — which presumes the caller named
    their local variable ``solution``, a fact the library cannot know, and
    which is simply wrong for the three families whose entries return a
    ``*Result``.  Threading a per-entry spelling in as a parameter would
    re-commit exactly the retired defect above: a fact asserted by the call
    site, free to drift from the object it describes.  The message names the
    ATTRIBUTE and its TYPE instead, which is true in every family and which
    each result type's own ``record`` docstring already points at.

    Parameters
    ----------
    record : IterationRecord
        The solve's iteration tree.  Read-only; nothing here mutates it.
    where : str
        The public entry's name, for the message's subject.  The only
        caller-supplied fact left, because it is the only one the record
        cannot know.
    balance_defect : float or None, optional
        ``‖R_g‖/‖Q_g‖`` for the RETURNED iterate, when the family computes
        one.  SN supplies it (see
        :func:`~orpheus.sn.solver._exit_balance_defect`); CP, MoC and
        diffusion do not yet, and pass ``None``.  ``None`` renders as an
        ABSENT clause rather than the words "unavailable", because an empty
        clause cannot be misread as a measurement.

    Notes
    -----
    The message carries the budget that ran out, the tolerance that was not
    reached, and how far away the last iterate was, because "one more sweep"
    and "diverging" want opposite responses from the reader.

    ⭐ And it carries the OBSERVED rate with the budget that rate projects,
    so the advice is a number rather than a direction.  *"Raise max_inner"*
    leaves the reader to guess, and a guess against a ``rho = 0.985`` mode is
    wrong by a factor of six; *"set max_inner=1343"* is actionable in one
    step.  The projection is fitted, not assumed — see
    :meth:`StoppingCriterion.projected_iterations` for why the intercept
    matters (`[M]` 27 % of the count on the d=3 control).

    ⚠ The projection is only as sharp as the tail it was measured on, and it
    reads LOW when the budget was cut inside the transient — the message says
    "so far" for exactly that reason.  It converges fast, though, and it
    converges from below, so a reader who follows the advice gets a bigger
    number next time rather than a wrong answer.  `[M]` 2026-08-09 on the d=3
    all-reflective absorber (``inner_tol=1e-13``, true count **1631**):

    ======  =========  =========
    budget  observed   projected
    ======  =========  =========
    50      0.956890         586
    200     0.985181        1618
    800     0.985388        1634
    4000    0.985399        1633
    ======  =========  =========

    At the *old* default of 200 the projection is already within **0.8 %** of
    the truth; only the deep-transient row under-reads.

    The **non-contracting** case is reported differently on purpose.  When
    the observed rate is ``>= 1`` no budget suffices, and telling that reader
    to raise the budget would send them down a road with no end: `[M]` #340
    found a configuration whose NEGATIVE dominant eigenvalue makes the
    increment sign-alternate, so its stop test is unsatisfiable forever.

    ``stacklevel=3`` attributes the warning to the entry's CALLER — frame 1 is
    this function, frame 2 the public entry, frame 3 user code.

    ⛔ **That is a PRECONDITION on the call site, and it was violated at two of
    the eight sites when N4.7 first landed** (`[M]` 2026-08-11; corrected the
    same day).  ``solve_sn_fixed_source`` dispatches to
    ``_solve_fixed_source_si`` / ``_solve_fixed_source_krylov``, and the call
    lived inside those private arms — so frame 3 was the entry's own
    ``return _solve_fixed_source_si(`` dispatch line, and the warning pointed
    the reader at ``sn/solver.py``, a file they did not write.  The fix was to
    hoist the call into the public entry (reading the record off the
    ``Solution`` being returned), **not** to pass a per-call ``stacklevel``: a
    frame COUNT is a fact about the call chain asserted at the call site, and
    it rots silently the moment a helper is interposed — the same defect class
    as the retired ``budget_name``/``budget``/``tol`` arguments above.

    ⟹ **Call this DIRECTLY from a public entry, one frame deep.** Anything
    else re-blames ``orpheus``.

    ⚠ And the obvious gate for it is **Mode-12 blind to half the class**.  A
    check of the form *"the warning is not attributed inside ``orpheus/``"*
    catches ``stacklevel`` being too SMALL (which lands in the package) and is
    structurally incapable of catching it being too LARGE (which lands in the
    caller's caller — still outside the package, still wrong).  A gate here
    must assert the exact ``filename`` **and** ``lineno``, not a containment
    predicate.
    """
    if record.fully_converged:
        return

    # ⭐ #340 N6: every fact below is read off the level that FAILED, which is
    # not necessarily the level the caller asked about.  ``first_failure``
    # searches CHILDREN FIRST, so a starved inner is named rather than the
    # outer it starved.
    #
    # ⛔ Until 2026-08-10 the budget, its NAME and the tolerance were passed in
    # by each entry point and described that entry's TOP level, while the rate
    # and projection below were already read off the record.  On any solve
    # whose inner starved, that welded one level's knob to another level's
    # trajectory: every number real, every pairing wrong, and — worst —
    # the result LOOKS level-correct.  It also advised the wrong knob, because
    # with a starved inner raising ``max_outer`` cannot help at all (the
    # outer's stop test is entirely increments, and the starved inner is what
    # suppresses them — F2).
    #
    # ⭐ The ``or`` arm is now PROVABLY dead, and that is the point of the
    # widened guard: ``fully_converged`` is ``self.converged and all(child
    # .fully_converged)`` and ``first_failure`` returns ``None`` iff
    # ``self.converged`` and every child's does — complementary by induction,
    # so past the guard above ``first_failure`` is never ``None``.  Under the
    # OLD ``converged`` guard the arm was also dead, but only incidentally (a
    # failing top level always names itself); the two predicates now coincide
    # exactly.  Kept so the level is a non-optional value for every read below
    # rather than a defended-against one.
    failing = record.first_failure or record
    level = "" if failing is record else f"{failing.label} "
    budget = failing.budget
    budget_name = budget.name

    # The criterion that ACTUALLY bound — the one furthest from clearing,
    # read straight off the record rather than reconstructed from a flat
    # list (#340 N2b-ii).
    #
    # ⛔ Until 2026-08-09 this rebuilt a criterion by hand: `flux_residuals`
    # when non-empty, else re-differencing `keff_history` into `|dk|`.  That
    # could only ever recover ONE of the outer's two criteria, and it picked
    # the wrong one exactly when it mattered.  `[M]` caught by the wide run,
    # not by review: on the mutated heterogeneous slab, whose NEGATIVE
    # dominant eigenvalue makes `dphi` alternate in sign forever, `|dk|` sat
    # at 3.3e-16 against a 1e-9 tolerance — so the projection dutifully
    # answered "you need 1 iteration" while the solve could never converge.
    # The stop-gap was a branch that detected "the recorded quantity cleared"
    # and refused to project; it is retired here, because the record carries
    # BOTH criteria and `binding_criterion` picks the one that failed.
    criterion = failing.binding_criterion
    distance = (
        "no criterion recorded"
        if criterion is None or criterion.last is None
        else f"last {criterion.name} {criterion.last:.3e}"
    )
    # A level with no criteria at all has no tolerance to quote, so the clause
    # is dropped rather than faked.
    against = (
        "" if criterion is None
        else f" without reaching tol={criterion.tolerance:.3e}"
    )

    rate = None if criterion is None else criterion.rate
    projected = None if criterion is None else criterion.projected_iterations()
    if criterion is not None and criterion.cleared:
        # Every criterion cleared, yet the level did not converge — so the
        # refusal came from the LOOP, not from a quantity: too few iterations
        # to claim (`min_iterations`), or a level that ran and measured
        # nothing.  Both are real states and neither is a budget problem, so
        # naming the state beats projecting a number.
        #
        # ⚠ This is NOT the retired stop-gap it replaces.  That branch fired
        # when the failing criterion was ABSENT from the history; this one
        # fires only when there is no failing criterion at all, which is a
        # statement the record can actually make.
        shortfall = failing.min_iterations
        advice = (
            f"every recorded criterion cleared, so the refusal is the loop's: "
            f"it ran {failing.n_iterations} of the {shortfall} "
            f"iterations required before convergence may be claimed. Read"
            if failing.n_iterations < shortfall
            else
            f"every recorded criterion cleared but the level measured "
            f"nothing — it ran {failing.n_iterations} iterations "
            f"without recording a value, so there is nothing to project "
            f"from. Read"
        )
    elif rate is None:
        advice = f"Raise {budget_name}, or read"
    elif projected is None:
        advice = (
            f"⛔ the observed rate is rho={rate:.6f} — this iteration is NOT "
            f"contracting, so NO {budget_name} suffices and raising it will "
            f"not help. Check the problem, not the budget. Or read"
        )
    else:
        # ⛔ ``projected`` is fitted over the TRAJECTORY, so it is a count of
        # recorded iterations — and ``budget_name`` names a knob that may not
        # take those.  Until 2026-08-13 both halves of "needs about N: set
        # X=N" printed the same N, which on the GMRES arm handed the reader
        # an Arnoldi-step count to type into a restart-cycle knob, over by a
        # factor of ``restart`` = ``n_dof`` (#349).  ``covering`` is the
        # conversion, and it is the identity for every unit-consistent
        # producer — so the two numbers still coincide wherever they should.
        setting = budget.covering(projected)
        advice = (
            f"At the rate observed so far (rho={rate:.6f}) this needs about "
            f"{projected} iterations: set {budget_name}={setting}. Or read"
        )

    # ⭐ Its OWN sentence, with its own subject, and deliberately not folded
    # into the clause above.  Every other number in this message belongs to
    # ``failing`` — the level ``first_failure`` named, which on an eigenvalue
    # solve is outer-iteration ONE's inner.  The balance defect belongs to
    # the RETURNED ITERATE, i.e. the last one.  Both are real; appending it
    # to the failing level's clause would read as one level's facts and be
    # the exact "every number real, every pairing wrong" defect N6a removed.
    # Absent (LD schemes, zero source, and every family but SN) it says
    # nothing rather than "unavailable" — an empty clause cannot be misread
    # as a measurement.
    balance = (
        "" if balance_defect is None else
        f"The returned iterate leaves a per-group balance defect of "
        f"‖R_g‖/‖Q_g‖ = {balance_defect:.3e} — a DIAGNOSTIC "
        f"magnitude, NOT a verdict: it tracks the error in keff better than "
        f"the raw residual does, but benign and corrupting solves overlap, "
        f"so weigh it and do not threshold it. "
    )
    warnings.warn(
        # ``{budget}`` renders as ``max_inner=200`` for every unit-consistent
        # producer — byte-identical to the retired ``{budget_name}={budget}``
        # over the raw int pair — and as ``max_inner=5 (x128 = 640
        # iterations)`` for GMRES, which is the one arm where the knob alone
        # does not say what the loop was allowed to do (#349).
        f"{where}: {level}hit {budget}{against} "
        f"({distance}). Returning a BEST-EFFORT iterate — it is mid-descent, "
        f"not the converged answer. {balance}{advice} "
        # ``fully_converged``, not ``converged``: this warning fires for ANY
        # level, and on a starved-inner solve the flat ``converged`` reads
        # True — so naming that one would send the reader to the single
        # predicate that cannot see what they were just warned about (#340
        # N6b).  Named as ATTRIBUTE-on-TYPE rather than as a variable path,
        # for the reason in this function's docstring.
        f"`fully_converged` on the IterationRecord this solve returned, "
        f"and handle it. "
        f"Silence this per-call with warnings.catch_warnings(); make it fatal "
        f"everywhere with {ESCALATION_FLAG}.",
        ConvergenceWarning,
        stacklevel=3,
    )
