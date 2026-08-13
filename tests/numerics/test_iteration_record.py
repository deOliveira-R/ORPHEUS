"""The DEFINING laws of the nested-iteration record (#340, plan step N1).

:class:`~orpheus.numerics.convergence.IterationRecord` is a rooted tree of
iteration levels carrying a **monotone AND-fold** (``fully_converged``), and
:class:`~orpheus.numerics.convergence.StoppingCriterion` is a magnitude
trajectory judged against one tolerance.  Both are math-bearing types, so
this module gates their *intrinsic laws* — the properties that define what
they ARE — rather than a sample of usage:

* the fold agrees with an independent flat traversal (recursion vs iteration);
* ``first_failure is None`` iff ``fully_converged`` (two derived surfaces,
  cross-checked against each other);
* ``cleared`` iff ``distance < 1`` (one statement, two spellings);
* the rate estimator recovers ``rho`` from a synthetic geometric decay, and
  the projection ROUND-TRIPS against a direct search of the same decay — a
  structurally-independent check (closed-form fit vs stepping the sequence).

⭐ The single most load-bearing gate is
:class:`TestConvergedIsDerivedNotStored`.  #342 was five sites each spelling
their own ``converged``, one of them a hardcoded ``True``; the remedy is not
"fix the five" but "make the sixth unspellable".  If a ``converged`` field
is ever added to either dataclass, that class of bug is back and this file
reds.

Runtime note (vv Mode 8): this module is COLLECTED, so pytest's assertion
rewriter keeps bare ``assert`` live under the canonical ``python -O``.
"""

from __future__ import annotations

import dataclasses
import math

import pytest

from orpheus.numerics.convergence import (
    IterationBudget,
    IterationRecord,
    StoppingCriterion,
)


def _geometric(r0: float, rho: float, n: int) -> tuple[float, ...]:
    """``n`` entries of the exact decay ``r_k = r0 * rho**k``."""
    return tuple(r0 * rho**k for k in range(n))


def _iterations_until_below(r0: float, rho: float, tol: float) -> int:
    """The SAME question the projection answers, by stepping the sequence.

    Structurally independent of :meth:`StoppingCriterion.projected_iterations`
    — that one fits a line in log space and solves for the crossing; this one
    walks the recurrence.  Neither can be wrong alone without the round-trip
    gate below reddening.
    """
    index, residual = 0, r0
    while residual >= tol:
        residual *= rho
        index += 1
    return index + 1  # entries 0..index inclusive


def _criterion(*, name: str = "residual", **kw) -> StoppingCriterion:
    return StoppingCriterion(name=name, **kw)


# ─── 1. The criterion is a magnitude trajectory ──────────────────────────


@pytest.mark.foundation
class TestStoppingCriterionRefusesWhatItCannotMean:
    """Construction guards, each with the positive leg that proves the guard
    is not simply always-raising (vv anti-pattern #11)."""

    def test_a_well_formed_criterion_is_accepted(self) -> None:
        """The positive leg.  Without it, a guard that rejects EVERYTHING
        would satisfy every negative test below."""
        criterion = _criterion(trajectory=(1e-2, 1e-4, 1e-6), tolerance=1e-5)
        assert criterion.n_iterations == 3
        assert criterion.last == 1e-6

    def test_a_negative_magnitude_is_a_producer_bug(self) -> None:
        with pytest.raises(ValueError, match="MAGNITUDES"):
            _criterion(trajectory=(1e-2, -1e-4), tolerance=1e-5)

    @pytest.mark.parametrize("tolerance", [-1e-8, math.nan])
    def test_a_negative_or_nan_tolerance_is_refused(
        self, tolerance: float
    ) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            _criterion(trajectory=(1e-2,), tolerance=tolerance)

    def test_a_ZERO_tolerance_is_legal_and_never_clears(self) -> None:
        """⛔ The guard first shipped demanding a strictly positive tolerance,
        and two PRODUCTION paths could then not build a record at all:
        ``GreenOperator(tol=0)`` is how the unsatisfiability guard is
        exercised, and GMRES's exact-breakdown path reaches a literal ``0.0``
        residual.  A type must not assert more than its callers promise (vv
        anti-pattern #16, in the direction that breaks production rather than
        the direction that yields a false red).

        Nothing about the comparison changed: the test stays strict, so a
        zero tolerance never clears — exactly as the retired
        ``_claims_convergence`` behaved.
        """
        criterion = _criterion(trajectory=(1e-2, 0.0), tolerance=0.0)
        assert criterion.cleared is False
        assert criterion.distance == math.inf  # keeps the two spellings tied
        assert criterion.projected_iterations() is None  # no finite count

    def test_an_unnamed_criterion_is_refused(self) -> None:
        with pytest.raises(ValueError, match="named"):
            StoppingCriterion(name="", trajectory=(1.0,), tolerance=1e-8)

    def test_a_diverged_nan_is_RECORDABLE_and_never_clears(self) -> None:
        """A diverging solve is a real state that must be representable.

        Rejecting ``nan`` at construction would make the one solve most in
        need of a diagnostic the one that cannot produce one.  Since
        ``nan < tol`` is ``False``, such a criterion correctly never clears
        — the honesty is free, no branch required.
        """
        criterion = _criterion(trajectory=(1e-2, math.nan), tolerance=1e-5)
        assert criterion.cleared is False
        assert criterion.rate is None  # not estimable from a non-finite tail

    def test_an_infinite_residual_is_recordable_and_never_clears(self) -> None:
        assert _criterion(trajectory=(math.inf,), tolerance=1e-5).cleared is False


@pytest.mark.foundation
class TestClearedAndDistanceAreOneStatement:
    """``cleared`` iff ``distance < 1`` — two spellings of one fact.

    They are derived independently (a comparison against ``tolerance`` vs a
    ratio to it), so agreement across the boundary cases is a real law and
    not a tautology.  The exact-equality row is the one that matters: a
    stopping test is strict, so hitting the tolerance exactly is NOT clearing
    it, and both spellings must say so.
    """

    @pytest.mark.parametrize(
        "last, tolerance",
        [(1e-9, 1e-8), (1e-8, 1e-8), (1e-7, 1e-8), (0.0, 1e-8),
         (0.0, 0.0), (1e-9, 0.0)],
    )
    def test_the_two_spellings_agree(self, last: float, tolerance: float) -> None:
        criterion = _criterion(trajectory=(1e-2, last), tolerance=tolerance)
        assert criterion.cleared == (criterion.distance < 1.0)

    def test_hitting_the_tolerance_exactly_does_not_clear_it(self) -> None:
        criterion = _criterion(trajectory=(1e-8,), tolerance=1e-8)
        assert criterion.cleared is False
        assert criterion.distance == 1.0

    def test_a_level_that_never_iterated_clears_vacuously(self) -> None:
        """``all(())`` is ``True``, and that is the correct reading.

        `[M]` #340's audit: reusing a production predicate that returns False
        for an EMPTY history turned "GMRES converged on the initial guess"
        into 44 phantom truncations of 90 rows.  A criterion that was never
        measured is not one that FAILED.
        """
        criterion = _criterion(trajectory=(), tolerance=1e-8)
        assert criterion.cleared is True
        assert criterion.last is None
        assert criterion.distance == 0.0  # not blocking, so it cannot bind


# ─── 2. The rate and the projection ──────────────────────────────────────


@pytest.mark.foundation
class TestTheRateIsRecoveredFromASyntheticDecay:
    """On an exact geometric decay the estimator must return ``rho`` itself.

    This is the pure-math anchor: the answer is known in closed form, so no
    solver, tolerance, or reference is involved.  The spread of rates covers
    the fast (``0.1``) through the campaign's measured d=3 worst (``0.9854``)
    to the near-stuck (``0.99575``, the ``Sigma_t/4`` configuration the
    served-rate default deliberately does NOT cover).
    """

    @pytest.mark.parametrize("rho", [0.1, 0.5, 0.9, 0.9854, 0.99575])
    def test_the_fitted_rate_is_the_generating_rate(self, rho: float) -> None:
        criterion = _criterion(
            trajectory=_geometric(3.7e-2, rho, 200), tolerance=1e-10
        )
        assert criterion.rate == pytest.approx(rho, rel=1e-9)

    def test_two_points_are_enough_for_a_rate(self) -> None:
        """The tail window must never starve itself below the two points a
        rate needs — ``min(int(n/2), n-2)`` is what guarantees it."""
        criterion = _criterion(trajectory=(1.0, 0.25), tolerance=1e-8)
        assert criterion.rate == pytest.approx(0.25)

    def test_one_point_has_no_rate(self) -> None:
        assert _criterion(trajectory=(1.0,), tolerance=1e-8).rate is None

    def test_no_iterations_has_no_rate(self) -> None:
        assert _criterion(trajectory=(), tolerance=1e-8).rate is None


@pytest.mark.foundation
class TestTheProjectionRoundTrips:
    """The projected count is the smallest N whose last entry clears — proved
    against a direct walk of the same recurrence.

    ⭐ This is the gate that pins the campaign's ``[M]`` lesson that a
    spectral radius predicts a RATE, not a COST.  The projection uses the
    FITTED intercept; the naive reading ``n = |ln(tol/r0)|/|ln rho|`` assumes
    the first residual lies on the asymptotic line.  On a pure geometric decay
    the two coincide by construction — which is exactly why the fixture is
    pure: it isolates the crossing arithmetic (the ``floor + 2``) from the
    intercept question, and any off-by-one in the crossing reds here.
    """

    @pytest.mark.parametrize("rho", [0.1, 0.5, 0.9, 0.9854])
    @pytest.mark.parametrize("tol", [1e-6, 1e-10])
    def test_the_projection_equals_a_direct_search(
        self, rho: float, tol: float
    ) -> None:
        r0 = 3.7e-2
        criterion = _criterion(trajectory=_geometric(r0, rho, 60), tolerance=tol)
        assert criterion.projected_iterations() == _iterations_until_below(
            r0, rho, tol
        )

    def test_the_projected_iterate_is_the_FIRST_one_below_tolerance(self) -> None:
        """Both halves of "smallest N": entry ``N-1`` clears and ``N-2`` does
        not.  Asserting only the first half would be satisfied by any
        over-estimate."""
        r0, rho, tol = 3.7e-2, 0.9854, 1e-10
        needed = _criterion(
            trajectory=_geometric(r0, rho, 60), tolerance=tol
        ).projected_iterations()
        assert needed is not None
        assert r0 * rho ** (needed - 1) < tol
        assert r0 * rho ** (needed - 2) >= tol

    def test_a_TRANSIENT_is_not_charged_to_the_asymptotic_rate(self) -> None:
        r"""⭐⭐ The intercept must be FITTED, not read off ``r0``.

        The pure-geometric rows above cannot see this: there ``r0`` lies on
        the asymptotic line by construction, so the naive law
        :math:`n = |\ln(\mathrm{tol}/r_0)|/|\ln\rho|` and the fitted crossing
        coincide exactly.  `[M]` 2026-08-09 mutation MA6 (substitute
        ``trajectory[0]`` for the fitted intercept) reddened **nothing**
        across the whole file — a measured blind spot, closed here.

        A real iteration starts off the asymptotic line: the initial guess
        carries content in the sub-dominant modes, which decays FAST and
        leaves the slow mode behind.  This fixture is that shape exactly —
        :math:`r_k = A\rho^k + B\sigma^k` with a large, fast transient
        (:math:`B/A = 50`, :math:`\sigma = 0.5 \ll \rho = 0.98`).

        `[M]` the two readings differ by **195 of 685 iterations = 28.5 %**,
        which is the same magnitude as the campaign's measured 27 % intercept
        share on the d=3 all-reflective control.  The reference is a direct
        walk of the true two-mode sequence — structurally independent of both
        the fit and the naive law.
        """
        A, rho, B, sigma, tol = 1.0, 0.98, 50.0, 0.5, 1e-6
        trajectory = tuple(A * rho**k + B * sigma**k for k in range(60))
        criterion = _criterion(trajectory=trajectory, tolerance=tol)

        index, residual = 0, trajectory[0]
        while residual >= tol:
            index += 1
            residual = A * rho**index + B * sigma**index
        truth = index + 1

        assert criterion.projected_iterations() == truth

        # ...and the naive reading is decisively wrong on the same data, so
        # the row above is a real discriminator and not a coincidence.
        naive = (
            math.floor(
                (math.log(tol) - math.log(trajectory[0]))
                / math.log(criterion.rate)  # type: ignore[arg-type]
            )
            + 2
        )
        assert naive > truth + 100

    def test_projecting_against_a_looser_tolerance_costs_less(self) -> None:
        """The "what would I get for my budget?" question."""
        criterion = _criterion(
            trajectory=_geometric(3.7e-2, 0.9854, 60), tolerance=1e-10
        )
        loose = criterion.projected_iterations(tolerance=1e-6)
        tight = criterion.projected_iterations()
        assert loose is not None and tight is not None
        assert loose < tight

    @pytest.mark.parametrize("rho", [1.0, 1.05])
    def test_a_NON_CONTRACTING_trajectory_has_NO_projection(
        self, rho: float
    ) -> None:
        """At ``rho >= 1`` no budget suffices, and a number would be a lie.

        `[M]` #340 found a mutated heterogeneous slab whose NEGATIVE dominant
        eigenvalue makes the flux increment sign-alternate, so ``dphi <
        flux_tol`` is unsatisfiable *forever*.  Raising the budget is futile
        there, and that is a structurally different failure from a shortfall
        — reporting "needs ~10^6" would send the reader down the wrong path.
        """
        criterion = _criterion(
            trajectory=_geometric(1e-3, rho, 40), tolerance=1e-10
        )
        assert criterion.rate is not None and criterion.rate >= 1.0
        assert criterion.projected_iterations() is None


# ─── 2b. The budget is a ceiling IN A UNIT (#349 / ERR-079) ──────────────


@pytest.mark.foundation
@pytest.mark.catches("ERR-079")
class TestTheBudgetOwnsItsExchangeRate:
    r"""The DEFINING laws of :class:`IterationBudget`, as a math object.

    It is a ceiling ``L`` in the caller's knob units together with a rate
    ``p`` converting one knob unit into recorded iterations, so it names the
    map :math:`L \mapsto L\,p` and its inverse
    :math:`n \mapsto \lceil n/p \rceil`.  These rows gate that pair — not a
    sample of usage — because the pair IS the type's reason to exist: before
    it, the forward map was welded into a comparison and the inverse did not
    exist at all, which is ERR-079.
    """

    @pytest.mark.parametrize("per_unit", [1, 2, 8, 144])
    @pytest.mark.parametrize("wanted", [1, 7, 63, 64, 65, 800])
    def test_covering_buys_what_was_asked_and_not_a_unit_more(
        self, per_unit: int, wanted: int,
    ) -> None:
        """``covering`` is the least ``L`` with ``L*p >= n`` — both halves.

        The upper half alone is satisfied by "recommend a huge number"; the
        lower half alone by "recommend zero".  Only the pair pins the ceiling
        division, and it is the pair a reader relies on when they type the
        advised setting: it must WORK, and it must not over-provision by a
        whole cycle (which on the SN arm is ``n_dof`` iterations).
        """
        budget = IterationBudget(0, "knob", iterations_per_unit=per_unit)
        setting = budget.covering(wanted)
        assert setting * per_unit >= wanted, "the advice must actually buy it"
        assert (setting - 1) * per_unit < wanted, "and not a unit more"

    def test_covering_round_trips_an_exact_multiple(self) -> None:
        """No off-by-one at the boundary, which is where ceiling division
        goes wrong: ``covering(k*p) == k`` exactly, never ``k+1``."""
        budget = IterationBudget(0, "knob", iterations_per_unit=8)
        for k in (1, 2, 5, 17):
            assert budget.covering(k * 8) == k

    def test_wanting_no_iterations_advises_no_budget(self) -> None:
        """The degenerate end, stated rather than left to the arithmetic:
        ``-(-0 // p)`` is 0 and the ``max(1, ...)`` floor must not turn a
        request for nothing into a recommendation of one."""
        assert IterationBudget(0, "knob", iterations_per_unit=8).covering(0) == 0

    def test_the_ceiling_binds_at_in_iterations_not_at_the_limit(self) -> None:
        """⛔ ERR-079 in one row: the boundary is ``L*p``, never ``L``.

        `[M]` the shipped defect read ``exhausted`` at 91 steps against a
        limit of 5; here the same shape is checked at every neighbouring
        count, so an off-by-one in either direction reds.
        """
        budget = IterationBudget(5, "max_inner", iterations_per_unit=30)
        assert budget.in_iterations == 150
        assert budget.exhausted_by(5) is False    # the retired reading
        assert budget.exhausted_by(149) is False
        assert budget.exhausted_by(150) is True
        assert budget.exhausted_by(151) is True

    @pytest.mark.parametrize("limit", [1, 2, 50, 1000])
    @pytest.mark.parametrize("ran", [0, 1, 49, 50, 51, 1001])
    def test_the_identity_rate_REPRODUCES_the_retired_int_comparison(
        self, limit: int, ran: int,
    ) -> None:
        """The compatibility law, and the reason five producers changed by
        nothing.

        ``iterations_per_unit=1`` must agree with the retired
        ``budget > 0 and n_iterations >= budget`` on EVERY input, so the
        carve is behaviour-preserving for `SourceIteration`, power iteration,
        CP, MoC and diffusion by construction rather than by inspection of
        their messages.  Only GMRES, which states a rate ≠ 1, moves.
        """
        retired = limit > 0 and ran >= limit
        assert IterationBudget(limit, "knob").exhausted_by(ran) is retired

    @pytest.mark.parametrize("ran", [0, 1, 10**6])
    def test_an_UNBUDGETED_level_is_never_exhausted_however_long_it_ran(
        self, ran: int,
    ) -> None:
        """``limit == 0`` is how a DIRECT solve says it cannot be truncated.

        Diffusion's inner is one LU back-substitution and relies on exactly
        this; a rate of any size must not resurrect a ceiling from a limit of
        zero (``0 * p == 0``, which a naive ``>=`` would read as *always*
        exhausted — the inverse of ERR-079 and the reason
        :attr:`is_budgeted` guards the comparison).
        """
        for per_unit in (1, 144):
            budget = IterationBudget(0, "knob", iterations_per_unit=per_unit)
            assert budget.is_budgeted is False
            assert budget.exhausted_by(ran) is False

    def test_it_renders_the_conversion_ONLY_when_it_is_not_the_identity(
        self,
    ) -> None:
        """``__str__`` is consumed by the ConvergenceWarning's subject line.

        The identity case must stay byte-identical to the retired
        ``f"{budget_name}={budget}"`` over the raw int pair — that is what
        keeps five families' messages unchanged — while the one arm where the
        knob does not say what the loop was allowed to do declares itself.
        """
        assert str(IterationBudget(200, "max_inner")) == "max_inner=200"
        assert str(IterationBudget(5, "max_inner", iterations_per_unit=144)) == (
            "max_inner=5 (x144 = 720 iterations)"
        )


# ─── 3. The record is a tree with a monotone fold ────────────────────────


def _leaf(label: str, *, last: float, tol: float, budget: int) -> IterationRecord:
    return IterationRecord(
        label=label,
        criteria=(_criterion(trajectory=(1e-1, last), tolerance=tol),),
        budget=IterationBudget(budget),
    )


def _three_level_tree(*, deep_last: float) -> IterationRecord:
    """outer -> two inners, the second of which owns a krylov leaf."""
    krylov = _leaf("krylov(gmres)", last=deep_last, tol=1e-12, budget=2)
    inner_b = IterationRecord(
        label="inner(within-group g=1)",
        criteria=(_criterion(trajectory=(1e-1, 1e-11), tolerance=1e-10),),
        budget=IterationBudget(200),
        children=(krylov,),
    )
    return IterationRecord(
        label="outer(power)",
        criteria=(_criterion(name="dk", trajectory=(1e-1, 1e-9), tolerance=1e-8),),
        budget=IterationBudget(500),
        children=(_leaf("inner(within-group g=0)", last=1e-11, tol=1e-10,
                        budget=200), inner_b),
    )


@pytest.mark.foundation
class TestTheFoldAgreesWithAnIndependentTraversal:
    """``fully_converged`` is the AND-fold over the tree.

    The property is defined by RECURSION; ``walk()`` flattens the tree by
    ITERATION.  Computing the same conjunction both ways is a genuine
    cross-check — a recursion that skipped a subtree, or a ``walk`` that
    missed a node, breaks the agreement without either being obviously wrong
    on its own.
    """

    @pytest.mark.parametrize("deep_last", [1e-13, 1e-3])
    def test_recursion_and_traversal_agree(self, deep_last: float) -> None:
        tree = _three_level_tree(deep_last=deep_last)
        assert tree.fully_converged == all(
            record.converged for record in tree.walk()
        )

    def test_walk_yields_every_node_exactly_once(self) -> None:
        tree = _three_level_tree(deep_last=1e-13)
        labels = [record.label for record in tree.walk()]
        assert len(labels) == 4
        assert len(set(labels)) == 4
        assert labels[0] == "outer(power)"  # self first, then depth-first

    def test_fully_converged_IMPLIES_converged(self) -> None:
        """Monotonicity of the fold: the root is one of its own conjuncts."""
        for deep_last in (1e-13, 1e-3):
            tree = _three_level_tree(deep_last=deep_last)
            assert not tree.fully_converged or tree.converged


@pytest.mark.foundation
class TestFirstFailureLocatesTheDEEPESTCause:
    """``first_failure`` and ``fully_converged`` are cross-checked opposites."""

    @pytest.mark.parametrize("deep_last", [1e-13, 1e-3])
    def test_no_failure_iff_fully_converged(self, deep_last: float) -> None:
        tree = _three_level_tree(deep_last=deep_last)
        assert (tree.first_failure is None) == tree.fully_converged

    def test_it_names_the_deepest_level_not_the_outermost(self) -> None:
        """The whole point: a starved krylov is reported, not the outer that
        merely inherited its failure."""
        tree = _three_level_tree(deep_last=1e-3)
        failure = tree.first_failure
        assert failure is not None
        assert failure.label == "krylov(gmres)"

    def test_a_parent_that_ALSO_failed_does_not_mask_its_child(self) -> None:
        """⭐ The row that gives the traversal ORDER its teeth.

        `[M]` 2026-08-09 mutation MA3 (check ``self`` before ``children``)
        left every other gate in this file GREEN, because in each fixture the
        only failing level was already the deepest — so both orders returned
        the same record and the ordering claim was unfalsifiable (vv Mode 12
        at the fixture, not the functional).  Discriminating requires a
        parent and a child failing TOGETHER: self-first names the parent,
        children-first names the cause.  A truncated inner starving its outer
        is the common shape in practice, so this is the realistic case, not a
        contrived one.
        """
        child = _leaf("inner(within-group)", last=1e-3, tol=1e-10, budget=2)
        parent = IterationRecord(
            label="outer(power)",
            criteria=(_criterion(name="dk", trajectory=(1e-1, 1e-2),
                                 tolerance=1e-8),),
            budget=IterationBudget(2),
            children=(child,),
        )
        assert parent.converged is False  # BOTH levels failed
        assert child.converged is False
        failure = parent.first_failure
        assert failure is not None
        assert failure.label == "inner(within-group)"

    def test_it_searches_children_even_when_SELF_converged(self) -> None:
        """`[M]` #340's measured defect, as a record.

        A 2-D all-reflective 2-group eigenvalue solve reported
        ``converged=True`` while a truncated inner had moved ``keff`` by
        **5.3x keff_tol**, because the outer stops on INCREMENTS and a
        starved inner suppresses exactly those.  The record must make that
        state legible rather than reproduce it: the outer's own verdict stays
        honestly ``True``, and the tree-wide verdict says ``False``.
        """
        tree = _three_level_tree(deep_last=1e-3)
        assert tree.converged is True  # its own increments did clear
        assert tree.fully_converged is False  # ...but the tree did not
        assert tree.first_failure is not None


# ─── 4. The four states, and the binding criterion ───────────────────────


@pytest.mark.foundation
class TestTheFourStatesAreExhaustiveAndDistinct:
    """A predicate makes one decision; a type classifies STATES.

    `[M]` #340's audit reused a production predicate as a detector and got 44
    of 90 rows wrong, because that predicate was also false for a state it
    was never asked about.  The remedy is one named property per state, so a
    consumer names the one it means — and one control per state, here.
    """

    @pytest.mark.parametrize(
        "criteria, budget, expected",
        [
            ((), 0, "DIRECT"),  # a factorisation: never iterated
            (((1e-1, 1e-11), 1e-10), 200, "CONVERGED"),
            (((1e-1, 1e-3), 1e-10), 2, "TRUNCATED"),  # budget gone, unmet
            (((1e-1, 1e-3), 1e-10), 200, "STOPPED"),  # unmet, budget left
        ],
    )
    def test_each_state_is_reached_and_named(
        self, criteria, budget: int, expected: str
    ) -> None:
        built = (
            ()
            if criteria == ()
            else (_criterion(trajectory=criteria[0], tolerance=criteria[1]),)
        )
        record = IterationRecord(
            label="level", criteria=built, budget=IterationBudget(budget),
        )
        assert record.status == expected

    def test_a_direct_solve_did_not_iterate_but_DID_converge(self) -> None:
        """The two questions are separate, and the Peierls/direct-solve
        families need exactly this distinction (``k_eff`` is ``float | None``
        because one type serves an iterating and a direct path).  Diffusion's
        inner is one LU back-substitution, so its record is a 1-level tree
        whose only honest reading is this one.

        The ``exhausted_budget`` row is load-bearing and was added after
        measurement: `[M]` 2026-08-09 mutation MA9 (drop the ``budget > 0``
        guard, so ``0 >= 0`` reads as exhausted) reddened **nothing** without
        it, because ``status`` short-circuits on ``iterated`` and
        ``truncated`` is masked by the vacuous ``converged``.  A level that
        never ran did not run OUT.
        """
        record = IterationRecord(label="direct(lu)")
        assert record.iterated is False
        assert record.converged is True
        assert record.exhausted_budget is False
        assert record.truncated is False

    def test_converging_ON_the_last_allowed_iteration_is_convergence(self) -> None:
        """Derived from the STATE, not from which way control flow left the
        loop.  A budget-exhausted exit whose criterion cleared is honest
        convergence; reporting it as truncation is the mirror-image lie.
        """
        record = IterationRecord(
            label="level",
            criteria=(_criterion(trajectory=(1e-1, 1e-11), tolerance=1e-10),),
            budget=IterationBudget(2),
        )
        assert record.exhausted_budget is True
        assert record.converged is True
        assert record.truncated is False

    def test_a_level_that_RAN_and_measured_NOTHING_cannot_claim(self) -> None:
        """⭐ Absence of evidence is not evidence of convergence.

        An empty trajectory is vacuously ``cleared``, which is right when
        the loop never entered (GMRES on its initial guess) and WRONG when
        it entered and never got to measure — a ``SourceIteration`` given
        ``max_iter=1`` makes one pass and records no residual, because its
        stop compares successive iterates.

        `[M]` this is the exact regression the retirement of
        ``_claims_convergence`` would otherwise have shipped: that predicate
        returned ``False`` on an empty history, so a one-pass SI solve read
        as unconverged; the naive replacement reads it as converged.  Both
        readings are wrong in the OTHER case, and ``iterated`` is the
        discriminator that makes one rule serve both.
        """
        unmeasured = _criterion(trajectory=(), tolerance=1e-10)
        ran = IterationRecord(
            label="inner", criteria=(unmeasured,),
            budget=IterationBudget(1), iterations_run=1,
        )
        assert ran.iterated is True
        assert ran.converged is False
        assert ran.truncated is True

        never_entered = IterationRecord(
            label="inner", criteria=(unmeasured,), budget=IterationBudget(200)
        )
        assert never_entered.iterated is False
        assert never_entered.converged is True  # nothing ran; nothing failed

    def test_min_iterations_withholds_a_premature_claim(self) -> None:
        """The home for SN's ``iteration <= 2`` guard — a statement about the
        LOOP, kept out of what ``cleared`` means for a quantity."""
        cleared = _criterion(trajectory=(1e-11,), tolerance=1e-10)
        assert IterationRecord(label="l", criteria=(cleared,)).converged is True
        assert (
            IterationRecord(
                label="l", criteria=(cleared,), min_iterations=3
            ).converged
            is False
        )


@pytest.mark.foundation
class TestTheBindingCriterionIsTheFurthestFromClearing:
    """One ratio serves both cases — no branch on converged-vs-not."""

    def test_when_converged_it_is_the_one_that_cleared_LAST(self) -> None:
        record = IterationRecord(
            label="outer(power)",
            criteria=(
                _criterion(name="dk", trajectory=(1e-2, 1e-12), tolerance=1e-8),
                _criterion(name="dphi", trajectory=(1e-2, 5e-9), tolerance=1e-8),
            ),
            budget=IterationBudget(500),
        )
        assert record.converged is True
        assert record.binding_criterion is not None
        assert record.binding_criterion.name == "dphi"

    def test_when_starved_it_is_the_one_that_failed_WORST(self) -> None:
        record = IterationRecord(
            label="outer(power)",
            criteria=(
                _criterion(name="dk", trajectory=(1e-2, 1e-7), tolerance=1e-8),
                _criterion(name="dphi", trajectory=(1e-2, 1e-3), tolerance=1e-8),
            ),
            budget=IterationBudget(2),
        )
        assert record.truncated is True
        assert record.binding_criterion is not None
        assert record.binding_criterion.name == "dphi"

    def test_a_level_with_no_criteria_has_no_binding_one(self) -> None:
        record = IterationRecord(label="direct(lu)")
        assert record.binding_criterion is None
        assert record.rate is None
        assert record.projected_iterations() is None


# ─── 5. The structural guarantee ─────────────────────────────────────────


@pytest.mark.foundation
class TestConvergedIsDerivedNotStored:
    """⭐ The load-bearing gate: there is no field to set wrongly.

    #342 was five sites each transcribing their own ``converged``, one of
    them a hardcoded ``True``, and a sixth (``... and iteration > 5``) turned
    up in ``derivations/``.  Fixing six sites leaves the seventh spellable.
    Removing the field does not.  (``coding-elegance`` Pattern 4 applied to a
    boolean; the same discipline as
    ``tests/sn/solve/test_convergence_contract.py``'s non-unpackable outcome.)
    """

    @pytest.mark.parametrize("cls", [IterationRecord, StoppingCriterion])
    def test_no_convergence_verdict_is_a_constructor_argument(self, cls) -> None:
        stored = {field.name for field in dataclasses.fields(cls)}
        assert not stored & {"converged", "fully_converged", "cleared",
                             "truncated", "status"}

    @pytest.mark.parametrize("cls", [IterationRecord, StoppingCriterion])
    def test_the_verdicts_are_read_only_properties(self, cls) -> None:
        assert isinstance(getattr(cls, "converged", None) or
                          getattr(cls, "cleared"), property)

    def test_a_verdict_cannot_be_overwritten_on_an_instance(self) -> None:
        record = IterationRecord(label="level")
        with pytest.raises((AttributeError, dataclasses.FrozenInstanceError)):
            record.converged = False  # type: ignore[misc]


@pytest.mark.foundation
class TestRecordConstructionEnforcesItsInvariants:
    """Guards, each paired with the positive leg (vv anti-pattern #11)."""

    def test_a_well_formed_record_is_accepted(self) -> None:
        record = IterationRecord(
            label="outer(power)",
            criteria=(
                _criterion(name="dk", trajectory=(1e-2, 1e-9), tolerance=1e-8),
                _criterion(name="dphi", trajectory=(1e-2, 1e-9), tolerance=1e-8),
            ),
            budget=IterationBudget(500),
        )
        assert record.n_iterations == 2

    def test_criteria_must_be_CO_INDEXED_by_iteration(self) -> None:
        """They measure the same iterations; a length mismatch means a
        producer dropped one, which would silently corrupt ``n_iterations``
        (read off the first criterion) and every state derived from it."""
        with pytest.raises(ValueError, match="co-indexed"):
            IterationRecord(
                label="outer",
                criteria=(
                    _criterion(name="dk", trajectory=(1e-2, 1e-9), tolerance=1e-8),
                    _criterion(name="dphi", trajectory=(1e-2,), tolerance=1e-8),
                ),
            )

    def test_criterion_names_must_be_unique(self) -> None:
        with pytest.raises(ValueError, match="unique"):
            IterationRecord(
                label="outer",
                criteria=(
                    _criterion(name="dk", trajectory=(1e-9,), tolerance=1e-8),
                    _criterion(name="dk", trajectory=(1e-9,), tolerance=1e-8),
                ),
            )

    def test_the_producer_may_state_a_count_that_EXCEEDS_its_measurements(
        self,
    ) -> None:
        """The offset is a per-producer fact and cannot be inferred.

        `[M]` both conventions ship: ``SourceIteration`` measures the
        DIFFERENCE between successive iterates, so an exhausted
        ``max_iter=50`` run records 49 residuals; ``KrylovAcceleration``
        gets one callback per iteration and records 50.  Before the record,
        the same ``n_inner`` field carried both — SI writing ``len`` and
        Krylov writing ``len + 1``, undocumented and exactly backwards.
        """
        criterion = _criterion(trajectory=tuple(1e-1 for _ in range(49)),
                               tolerance=1e-10)
        si_like = IterationRecord(
            label="inner(source-iteration)", criteria=(criterion,),
            budget=IterationBudget(50), iterations_run=50,
        )
        assert si_like.n_iterations == 50
        assert si_like.exhausted_budget is True  # would read False from len()
        assert si_like.truncated is True

        krylov_like = IterationRecord(
            label="inner(gmres)", criteria=(criterion,),
            budget=IterationBudget(50),
        )
        assert krylov_like.n_iterations == 49  # defaults to the measurements
        assert krylov_like.exhausted_budget is False

    def test_a_count_BELOW_the_measurements_is_refused(self) -> None:
        """You cannot measure a criterion more often than you iterated."""
        with pytest.raises(ValueError, match="fewer than"):
            IterationRecord(
                label="inner",
                criteria=(_criterion(trajectory=(1e-1, 1e-2), tolerance=1e-8),),
                iterations_run=1,
            )

    def test_an_unlabelled_record_is_refused(self) -> None:
        with pytest.raises(ValueError, match="labelled"):
            IterationRecord(label="")

    def test_an_unnamed_budget_knob_is_refused(self) -> None:
        """The advice has to give the reader a token to type (#340 N6).

        The knob name is what a ``ConvergenceWarning`` puts after "set", so
        an empty one degrades the message to ``hit =50`` — the same class of
        defect as an unlabelled record, and refused the same way.

        ⛔ The invariant MOVED on 2026-08-13 (#349): it used to live on
        ``IterationRecord.__post_init__`` and guard a ``budget_name: str``
        field, and it now lives on :class:`IterationBudget`, which owns the
        whole (limit, knob, exchange-rate) fact.  So the refusal fires at the
        BUDGET's construction, one frame earlier than it used to — which is
        the point of the move, since a budget can now be built and passed
        around without a record to hold it.
        """
        with pytest.raises(ValueError, match="budget name"):
            IterationBudget(50, "")

    def test_a_budget_unit_that_buys_no_iterations_is_refused(self) -> None:
        """The exchange rate must be able to bind (#349).

        A unit that buys zero iterations makes ``in_iterations`` zero for
        ANY limit, so ``exhausted_by`` would read ``True`` on a loop that
        ran once — the inverse of the defect this type was minted to fix,
        and worth refusing at the same door.  Negative is likewise nonsense.
        """
        for bad in (0, -1):
            with pytest.raises(ValueError, match="iterations_per_unit"):
                IterationBudget(50, "max_inner", iterations_per_unit=bad)

    def test_a_negative_limit_is_refused(self) -> None:
        """The positive leg of the pair lives in the state table above; this
        is the negative one, and it moved off the record with its concept."""
        with pytest.raises(ValueError, match="limit must be >= 0"):
            IterationBudget(-1)

    def test_the_knob_defaults_to_the_PRIMITIVE_s_own_parameter_name(
        self,
    ) -> None:
        """The positive leg — and the reason there is no whitelist.

        ⛔ `[M]` 2026-08-10 the three spellings actually shipped are
        ``{max_inner, max_outer, max_iter}``: a ``{max_inner, max_outer}``
        whitelist would refuse
        :class:`~orpheus.numerics.green_operator.GreenOperator`, whose own
        knob IS ``max_iter``.  That is ``vv`` anti-pattern #16 pointing the
        expensive way — a type asserting more than its callers promise —
        and this campaign already had to back one of those out
        (``GreenOperator(tol=0)``).  So the field is a free ``str`` whose
        default is the honest answer for a primitive nobody re-named.
        """
        assert IterationRecord(label="level").budget.name == "max_iter"
        assert IterationBudget(50).name == "max_iter"
        assert (
            IterationRecord(
                label="level", budget=IterationBudget(50, "max_inner"),
            ).budget.name
            == "max_inner"
        )

    def test_a_negative_min_iterations_is_refused(self) -> None:
        """⛔ This was parametrized over ``{"budget": -1}`` too, until #349
        made ``budget`` an :class:`IterationBudget` rather than an ``int``.

        The ``budget=-1`` row did not simply move — it stopped being a
        question this class can be ASKED.  A dataclass does not type-check
        its fields at runtime, so ``IterationRecord(budget=-1)`` now builds
        happily and fails later at the first ``.budget.exhausted_by(...)``;
        the guard against it is pyright, not a raise, and adding an
        ``isinstance`` here would be the harmful-stub anti-pattern
        (``coding-standards``: no runtime guard for a case the type system
        covers).  The VALUE invariant it used to assert lives on the type,
        in :meth:`test_a_negative_limit_is_refused` above — where it also
        covers budgets built without a record to hold them.
        """
        with pytest.raises(ValueError, match=">= 0"):
            IterationRecord(label="level", min_iterations=-1)

    def test_replace_RE_RUNS_the_invariant(self) -> None:
        """``dataclasses.replace`` routes back through ``__post_init__``, so
        a same-type-producing edit cannot bypass the law (``coding-elegance``
        Pattern 4 ∩ 2 — the invariant has exactly one home)."""
        record = IterationRecord(
            label="outer",
            criteria=(
                _criterion(name="dk", trajectory=(1e-2, 1e-9), tolerance=1e-8),
            ),
        )
        with pytest.raises(ValueError, match="co-indexed"):
            dataclasses.replace(
                record,
                criteria=record.criteria
                + (_criterion(name="dphi", trajectory=(1e-9,), tolerance=1e-8),),
            )


# ─── 6. The report is the artefact a user pastes into an issue ───────────


@pytest.mark.foundation
class TestTheReportCarriesWhatTheReaderMustAct_ON:
    """A bare "did not converge" cannot be acted on.

    The reader needs to tell "one more sweep" from "diverging", so the tree
    must name the level, the criterion, the distance, the rate, and the
    budget that would have sufficed — the same contract
    ``test_the_message_carries_budget_tolerance_and_distance`` pins for the
    single-level warning, lifted to the tree.
    """

    def test_it_names_level_criterion_distance_rate_and_the_budget_to_set(
        self,
    ) -> None:
        rho, r0, tol = 0.9854, 3.7e-2, 1e-10
        inner = IterationRecord(
            label="inner(within-group)",
            criteria=(
                _criterion(trajectory=_geometric(r0, rho, 200), tolerance=tol),
            ),
            budget=IterationBudget(200),
        )
        outer = IterationRecord(
            label="outer(power)",
            criteria=(
                _criterion(name="dk", trajectory=(1e-2, 1e-9), tolerance=1e-8),
            ),
            budget=IterationBudget(500),
            children=(inner,),
        )
        report = outer.report()

        assert "inner(within-group): TRUNCATED (200/200 iterations)" in report
        assert "binding" in report
        assert f"{rho:.6f}" in report  # the OBSERVED rate, not a guess
        needed = inner.projected_iterations()
        assert needed is not None and f"needs ~{needed}" in report
        assert "outer(power): CONVERGED" in report  # the honest per-level fact

    def test_a_non_contracting_level_says_NO_BUDGET_SUFFICES(self) -> None:
        record = IterationRecord(
            label="inner(stuck)",
            criteria=(
                _criterion(trajectory=_geometric(1e-3, 1.0, 40), tolerance=1e-10),
            ),
            budget=IterationBudget(40),
        )
        assert "no budget suffices" in record.report()

    def test_a_converged_tree_reports_no_rate_advice(self) -> None:
        """The other half of the pair — advice that always fires is noise,
        and noise gets filtered."""
        report = _three_level_tree(deep_last=1e-13).report()
        assert "needs ~" not in report
        assert "no budget suffices" not in report
