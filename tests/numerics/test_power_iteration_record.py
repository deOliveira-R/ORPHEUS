"""The laws of the OUTER iteration record (#340, plan step N2b-i).

:func:`~orpheus.numerics.eigenvalue.power_iteration` no longer asks its solver
"are you converged?".  It asks what the solver MEASURED, accumulates the
readings into trajectories, and derives the verdict.  This module gates that
change at the level it lives — the loop and the protocol — rather than through
any one physics solver:

* the min-outer rule has exactly ONE home, and no realizer can spell a second;
* the derived ``converged`` agrees with the break-on-test semantics it
  replaced, INCLUDING at the budgets where the two could disagree;
* a solver that measures nothing is refused rather than certified;
* the subtree carries the children THIS loop caused, and a reused solver
  instance cannot contribute a stale one (#340 F12, fixed structurally);
* the headline: an outer whose own criteria clear while an inner starved
  reads ``converged=True`` and ``fully_converged=False``, and names the level.

⭐ The instrument is a SCRIPTED solver, not a converging physics solve.  The
claims here are about accumulation, verdict and tree assembly; a real solve
couples all three to physics that can mask an off-by-one in any of them.  The
end-to-end SN readings live in ``tests/sn/solve/test_convergence_contract.py``
(one contract, one home) — this file is the loop in isolation.

Runtime note (vv Mode 8): this module is COLLECTED, so pytest's assertion
rewriter keeps bare ``assert`` live under the canonical ``python -O``.
"""
from __future__ import annotations

import dataclasses
import inspect

import numpy as np
import pytest

from orpheus.numerics.convergence import IterationRecord, StoppingCriterion
from orpheus.numerics.eigenvalue import (
    MINIMUM_OUTER_ITERATIONS,
    PowerIterationOutcome,
    RecordingSolver,
    power_iteration,
)

pytestmark = pytest.mark.foundation


# ── the instrument ───────────────────────────────────────────────────────


class _ScriptedSolver:
    """An :class:`EigenvalueSolver` whose readings are a fixed script.

    Every physics method is a no-op; the ONLY live behaviour is the reading
    sequence, so a red here is a statement about the loop.  The script is
    read one entry per outer and the last entry repeats, which lets a test
    say "clears from iteration k onward" without sizing the list to the
    budget.
    """

    def __init__(
        self,
        script: list[float],
        *,
        tolerance: float = 1.0,
        names: tuple[str, ...] = ("r",),
    ) -> None:
        self._script = list(script)
        self._tolerance = tolerance
        self._names = names
        self.calls = 0

    # ⚠ Parameter NAMES match the protocol exactly.  A surrogate that stands
    # in for a contract has to honour it; pyright reads a renamed keyword as
    # a non-conformance, and it is right to — a caller passing by keyword
    # would break on the real thing.
    def initial_flux_distribution(self) -> np.ndarray:
        return np.ones(2)

    def compute_fission_source(
        self, flux_distribution: np.ndarray, keff: float,
    ) -> np.ndarray:
        return flux_distribution

    def solve_fixed_source(
        self, fission_source: np.ndarray, flux_distribution: np.ndarray,
    ) -> np.ndarray:
        return flux_distribution

    def compute_keff(self, flux_distribution: np.ndarray) -> float:
        return 1.0

    def measure_stopping_criteria(
        self, keff: float, keff_old: float,
        flux_distribution: np.ndarray, flux_old: np.ndarray,
    ) -> tuple[StoppingCriterion, ...]:
        value = self._script[min(self.calls, len(self._script) - 1)]
        self.calls += 1
        return tuple(
            StoppingCriterion.reading(name, value, self._tolerance)
            for name in self._names
        )


class _RecordingScriptedSolver(_ScriptedSolver):
    """A :class:`RecordingSolver` that appends one inner record per outer.

    ``inner_records`` deliberately ACCUMULATES across solves and is never
    reset — that is the un-hygienic realizer the loop must be robust to.
    """

    def __init__(self, script: list[float], *, inner_converges: bool = True,
                 **kw) -> None:
        super().__init__(script, **kw)
        self.inner_records: list[IterationRecord] = []
        self._inner_converges = inner_converges

    def solve_fixed_source(
        self, fission_source: np.ndarray, flux_distribution: np.ndarray,
    ) -> np.ndarray:
        self.inner_records.append(
            IterationRecord(
                label="inner(scripted)",
                criteria=(
                    StoppingCriterion(
                        name="residual",
                        trajectory=(1e-12,) if self._inner_converges else (1.0,),
                        tolerance=1e-8,
                    ),
                ),
                budget=1 if not self._inner_converges else 100,
                iterations_run=1,
            )
        )
        return flux_distribution


class _MuteSolver(_ScriptedSolver):
    """Reports NO criteria — the state the loop must refuse."""

    def measure_stopping_criteria(
        self, keff: float, keff_old: float,
        flux_distribution: np.ndarray, flux_old: np.ndarray,
    ) -> tuple[StoppingCriterion, ...]:
        return ()


# ── 1. the min-outer rule has one home ───────────────────────────────────


class TestTheMinOuterRuleIsStatedOnce:
    """`[M]` ``if iteration <= 2: return False`` was transcribed in all FIVE
    realizers of the protocol.  These gates make the sixth unspellable."""

    def test_the_constant_is_the_rule(self) -> None:
        assert MINIMUM_OUTER_ITERATIONS == 3

    @pytest.mark.parametrize(
        "module_path,class_name",
        [
            ("orpheus.numerics.iteration", "KEigenvalue"),
            ("orpheus.sn.solver", "SNSolver"),
            ("orpheus.cp.solver", "CPSolver"),
            ("orpheus.moc.core", "MOCSolver"),
            ("orpheus.diffusion.solver", "DiffusionSolver"),
        ],
    )
    def test_no_realizer_can_see_the_iteration_index(
        self, module_path: str, class_name: str,
    ) -> None:
        """No realizer's measurement method takes ``iteration``.

        The teeth are STRUCTURAL, not stylistic.  A realizer that cannot see
        the iteration index cannot re-introduce a per-solver min-outer guard,
        so the rule stays where the algorithm owns it.  Checking the five by
        signature also pins the retirement itself: a survivor still spelling
        ``converged(...)`` fails the first assertion.
        """
        import importlib

        cls = getattr(importlib.import_module(module_path), class_name)
        assert hasattr(cls, "measure_stopping_criteria"), (
            f"{class_name} has not migrated off the retired predicate"
        )
        assert not hasattr(cls, "converged"), (
            f"{class_name}.converged survived the N2b retirement — two "
            f"stopping surfaces is exactly the twin this step removed"
        )
        params = inspect.signature(cls.measure_stopping_criteria).parameters
        assert "iteration" not in params, (
            f"{class_name} can see the iteration index and could therefore "
            f"re-spell the min-outer guard locally"
        )

    def test_readings_that_clear_immediately_still_wait_for_the_rule(
        self,
    ) -> None:
        """A criterion below tolerance from the FIRST outer does not convert
        into a claim until the third.

        This is the rule as BEHAVIOUR rather than as a constant, and it is the
        gate that would red if the loop's ``n >= MINIMUM_OUTER_ITERATIONS``
        were dropped while the constant stayed.
        """
        for budget, expected in ((1, False), (2, False), (3, True)):
            solver = _ScriptedSolver([0.0], tolerance=1.0)
            outcome = power_iteration(solver, max_iter=budget)
            assert outcome.converged is expected, (
                f"max_iter={budget}: expected converged={expected}"
            )

    def test_the_loop_stops_the_instant_the_rule_allows(self) -> None:
        """...and not one outer later — a min-outer rule that also DELAYED a
        genuinely converged solve would be a cost, not a safeguard."""
        solver = _ScriptedSolver([0.0], tolerance=1.0)
        outcome = power_iteration(solver, max_iter=50)
        assert outcome.record.n_iterations == MINIMUM_OUTER_ITERATIONS


# ── 2. converged is derived, and agrees with what it replaced ────────────


class TestConvergedIsDerivedNotStored:
    """The campaign's founding rule, at the outer level (#342)."""

    def test_the_outcome_has_no_converged_FIELD(self) -> None:
        names = {f.name for f in dataclasses.fields(PowerIterationOutcome)}
        assert "converged" not in names, (
            "a stored converged is a transcription, and a transcription of a "
            "convergence verdict is #342"
        )
        assert "record" in names

    def test_it_cannot_be_handed_a_verdict(self) -> None:
        with pytest.raises(TypeError):
            PowerIterationOutcome(            # type: ignore[call-arg]
                keff=1.0, keff_history=[1.0], flux_distribution=np.ones(1),
                record=IterationRecord(label="outer(probe)"),
                converged=True,
            )

    def test_it_reads_through_to_the_record(self) -> None:
        """Two surfaces, one fact — cross-checked so neither can drift."""
        for script, budget in (([0.0], 10), ([5.0], 4), ([0.0], 2)):
            solver = _ScriptedSolver(script, tolerance=1.0)
            outcome = power_iteration(solver, max_iter=budget)
            assert outcome.converged is outcome.record.converged


class TestTheDerivedVerdictMatchesBreakOnTest:
    """`[M]` N2b-i replaced a flag the loop SET on breaking with a verdict
    DERIVED from the trajectories.  The two agree by argument; these are the
    budgets where the argument could have been wrong.

    (``plan-authoring`` §8 — an "equivalent by construction" claim owes a
    measurement at the boundary, not only in the middle.)
    """

    def test_a_generous_budget_converges_below_it(self) -> None:
        solver = _ScriptedSolver([9.0, 9.0, 9.0, 0.0], tolerance=1.0)
        outcome = power_iteration(solver, max_iter=100)
        assert outcome.converged is True
        assert outcome.record.status == "CONVERGED"
        assert outcome.record.n_iterations == 4

    def test_converging_on_the_LAST_allowed_outer_is_still_converged(
        self,
    ) -> None:
        """The exact case the retired ``len(keff_history) < max_outer``
        inference got WRONG, and the reason the flag was introduced.  The
        derived reading must not regress to it: here the budget is spent AND
        the criterion cleared, so ``exhausted_budget`` is True while
        ``truncated`` is False."""
        solver = _ScriptedSolver([9.0, 9.0, 9.0, 0.0], tolerance=1.0)
        outcome = power_iteration(solver, max_iter=4)
        assert outcome.converged is True
        assert outcome.record.exhausted_budget is True
        assert outcome.record.truncated is False
        assert outcome.record.status == "CONVERGED"

    def test_one_outer_short_is_truncated(self) -> None:
        solver = _ScriptedSolver([9.0, 9.0, 9.0, 0.0], tolerance=1.0)
        outcome = power_iteration(solver, max_iter=3)
        assert outcome.converged is False
        assert outcome.record.truncated is True

    def test_a_zero_budget_never_entered_its_loop(self) -> None:
        solver = _ScriptedSolver([0.0], tolerance=1.0)
        outcome = power_iteration(solver, max_iter=0)
        assert outcome.converged is False
        assert outcome.record.iterated is False
        assert outcome.record.status == "DIRECT"
        assert outcome.record.criteria == ()
        assert solver.calls == 0


# ── 3. what the loop refuses ─────────────────────────────────────────────


class TestTheLoopRefusesWhatItCannotJudge:

    def test_a_solver_that_measures_nothing_is_refused(self) -> None:
        """An empty conjunction is vacuously true, so silence would otherwise
        be certified as convergence at outer 3 — a #342-class claim assembled
        out of no evidence at all."""
        with pytest.raises(ValueError, match="no criteria"):
            power_iteration(_MuteSolver([0.0]), max_iter=10)

    def test_a_tolerance_that_moves_mid_solve_is_refused(self) -> None:
        """A trajectory judged against two thresholds makes ``cleared``,
        ``distance`` and every projection off it meaningless."""
        first = StoppingCriterion.reading("r", 1.0, 1e-6)
        with pytest.raises(ValueError, match="tolerance moved"):
            first.extended_with(StoppingCriterion.reading("r", 0.5, 1e-9))

    def test_a_reading_of_a_DIFFERENT_criterion_is_refused(self) -> None:
        first = StoppingCriterion.reading("dk", 1.0, 1e-6)
        with pytest.raises(ValueError, match="per-quantity"):
            first.extended_with(StoppingCriterion.reading("dphi", 0.5, 1e-6))


class TestTheAccumulationPreservesTheTrajectory:

    def test_reading_is_the_explicit_spelling(self) -> None:
        assert StoppingCriterion.reading("r", 0.25, 1e-3) == StoppingCriterion(
            name="r", trajectory=(0.25,), tolerance=1e-3,
        )

    def test_extension_appends_in_order(self) -> None:
        c = StoppingCriterion.reading("r", 3.0, 1.0)
        for v in (2.0, 1.0):
            c = c.extended_with(StoppingCriterion.reading("r", v, 1.0))
        assert c.trajectory == (3.0, 2.0, 1.0)
        assert c.last == 1.0

    def test_extension_re_fires_the_magnitude_INVARIANT(self) -> None:
        """``extended_with`` routes through ``replace`` precisely so the
        construction law re-runs rather than being restated (``coding-elegance``
        Pattern 4 ∩ 2).  A negative reading must be refused on the way IN, not
        discovered later by a consumer of ``rate``."""
        c = StoppingCriterion.reading("r", 1.0, 1.0)
        with pytest.raises(ValueError, match="MAGNITUDES"):
            c.extended_with(
                StoppingCriterion(name="r", trajectory=(-1.0,), tolerance=1.0)
            )

    def test_the_loop_accumulates_one_reading_per_outer_per_criterion(
        self,
    ) -> None:
        solver = _ScriptedSolver([9.0, 8.0, 7.0, 6.0], names=("a", "b"))
        outcome = power_iteration(solver, max_iter=4)
        assert [c.name for c in outcome.record.criteria] == ["a", "b"]
        for criterion in outcome.record.criteria:
            assert criterion.trajectory == (9.0, 8.0, 7.0, 6.0)
        assert outcome.record.n_iterations == 4


# ── 4. the subtree ───────────────────────────────────────────────────────


class TestTheSubtreeIsWhatTHISLoopCaused:

    def test_a_non_recording_solver_yields_a_leaf(self) -> None:
        solver = _ScriptedSolver([0.0])
        assert not isinstance(solver, RecordingSolver)
        assert power_iteration(solver, max_iter=5).record.children == ()

    def test_one_child_per_outer(self) -> None:
        solver = _RecordingScriptedSolver([9.0, 9.0, 9.0, 0.0])
        outcome = power_iteration(solver, max_iter=10)
        assert outcome.record.n_iterations == 4
        assert len(outcome.record.children) == 4

    def test_a_REUSED_solver_contributes_no_stale_child(self) -> None:
        """`[M]` #340 F12 — an instance accumulator with no reset double-counts
        across two solves.  The fix is structural: the loop slices from the
        length it saw at entry, so the tree is correct even for a realizer
        with no reset hook at all (which is what ``_RecordingScriptedSolver``
        deliberately is).

        Stated as the LAW — *each outcome carries exactly the children of its
        own outers* — rather than as two hardcoded counts.  The counts are not
        equal to each other: this solver's script pointer also persists, so
        the second solve clears sooner and legitimately runs fewer outers.
        `[M]` a first draft of this gate asserted ``4`` for both and reddened
        on the correct code, which is the more useful failure to have had:
        a per-outcome invariant cannot be satisfied by a lucky number.

        Without the slice the second outcome would carry EVERY child the
        instance had ever produced — and it would look plausible.
        """
        solver = _RecordingScriptedSolver([9.0, 9.0, 9.0, 0.0])
        first = power_iteration(solver, max_iter=10)
        second = power_iteration(solver, max_iter=10)

        for label, outcome in (("first", first), ("second", second)):
            assert len(outcome.record.children) == outcome.record.n_iterations, (
                f"{label}: {len(outcome.record.children)} children for "
                f"{outcome.record.n_iterations} outers"
            )
        assert len(solver.inner_records) == (
            first.record.n_iterations + second.record.n_iterations
        ), "fixture drift: this solver is supposed to be the un-hygienic one"
        assert second.record.children[0] is not first.record.children[0], (
            "the second solve re-served the first solve's opening child"
        )


class TestAnOuterCanConvergeOnAStarvedInner:
    """⭐ THE #340 headline, in isolation.

    An increment-only outer stop CANNOT see an upstream throttle: a truncated
    inner suppresses the very increments the outer reads, so the outer stalls
    and calls the stall convergence.  Before this step the returned object had
    no way to say so.  Now the outer's own verdict and the tree's verdict are
    separate questions with separate answers.
    """

    def test_the_outer_claims_convergence_and_the_tree_refuses_it(
        self,
    ) -> None:
        solver = _RecordingScriptedSolver([9.0, 9.0, 0.0],
                                          inner_converges=False)
        record = power_iteration(solver, max_iter=20).record

        assert record.converged is True, "the outer's OWN criteria did clear"
        assert record.fully_converged is False, (
            "the AND-fold over the tree must refuse what the level alone "
            "would claim — this is the whole of F6"
        )

    def test_the_OUTCOME_reports_the_level_and_not_the_fold(self) -> None:
        """The surface a caller actually touches must be the LEVEL's verdict.

        ⭐ This gate exists because a mutation found its absence.  `[M]` making
        :attr:`PowerIterationOutcome.converged` read ``record.fully_converged``
        instead of ``record.converged`` reddened **zero** of 249 tests: every
        other gate reads ``outcome.record.…`` directly, and the one that does
        cross-check the two surfaces
        (``test_it_reads_through_to_the_record``) uses fixtures where the
        level and the fold agree — so the distinction was described in three
        docstrings and asserted nowhere.

        A starved inner is the ONLY configuration that separates them, which
        is the same lesson as the campaign's pure-vs-transient fixture finding:
        before crediting a fixture with covering a distinction, ask what the
        two candidate implementations would BOTH produce on it.
        """
        solver = _RecordingScriptedSolver([9.0, 9.0, 0.0],
                                          inner_converges=False)
        outcome = power_iteration(solver, max_iter=20)

        assert outcome.record.fully_converged is False, "fixture drift"
        assert outcome.converged is True, (
            "PowerIterationOutcome.converged must answer for the OUTER level. "
            "Folding the children in here would silently re-answer a "
            "different question under the same name — and every caller "
            "reading it would inherit the swap with no diff to see."
        )

    def test_it_names_the_level_that_failed(self) -> None:
        solver = _RecordingScriptedSolver([9.0, 9.0, 0.0],
                                          inner_converges=False)
        failure = power_iteration(solver, max_iter=20).record.first_failure
        assert failure is not None
        assert failure.label == "inner(scripted)"
        assert failure.status == "TRUNCATED"

    def test_a_healthy_inner_leaves_both_verdicts_agreeing(self) -> None:
        """The positive control.  Without it the two gates above are satisfied
        by a ``fully_converged`` pinned False, and the pair would certify a
        tree that can never say yes."""
        solver = _RecordingScriptedSolver([9.0, 9.0, 0.0],
                                          inner_converges=True)
        record = power_iteration(solver, max_iter=20).record
        assert record.converged is True
        assert record.fully_converged is True
        assert record.first_failure is None


# ── 5. the report reads honestly ─────────────────────────────────────────


class TestTheReportDoesNotContradictItsOwnStatus:

    def test_an_unmeasured_criterion_does_not_read_as_MET(self) -> None:
        """``cleared`` is vacuously True on an empty trajectory and that is
        correct for the verdict it feeds.  Printing ``met`` beside a level
        whose status line says TRUNCATED hands the reader two statements that
        cannot both be acted on — and the report exists to be pasted into a
        bug report unedited."""
        record = IterationRecord(
            label="inner(probe)",
            criteria=(
                StoppingCriterion(name="residual", trajectory=(),
                                  tolerance=1e-8),
            ),
            budget=1,
            iterations_run=1,
        )
        text = record.report()
        assert record.status == "TRUNCATED"
        assert "not measured" in text
        assert "met" not in text.replace("not measured", "")
