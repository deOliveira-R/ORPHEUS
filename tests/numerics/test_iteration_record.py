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

from orpheus.numerics.convergence import IterationRecord, StoppingCriterion


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

    @pytest.mark.parametrize("tolerance", [0.0, -1e-8])
    def test_a_nonpositive_tolerance_is_refused(self, tolerance: float) -> None:
        with pytest.raises(ValueError, match="strictly positive"):
            _criterion(trajectory=(1e-2,), tolerance=tolerance)

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
        [(1e-9, 1e-8), (1e-8, 1e-8), (1e-7, 1e-8), (0.0, 1e-8)],
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


# ─── 3. The record is a tree with a monotone fold ────────────────────────


def _leaf(label: str, *, last: float, tol: float, budget: int) -> IterationRecord:
    return IterationRecord(
        label=label,
        criteria=(_criterion(trajectory=(1e-1, last), tolerance=tol),),
        budget=budget,
    )


def _three_level_tree(*, deep_last: float) -> IterationRecord:
    """outer -> two inners, the second of which owns a krylov leaf."""
    krylov = _leaf("krylov(gmres)", last=deep_last, tol=1e-12, budget=2)
    inner_b = IterationRecord(
        label="inner(within-group g=1)",
        criteria=(_criterion(trajectory=(1e-1, 1e-11), tolerance=1e-10),),
        budget=200,
        children=(krylov,),
    )
    return IterationRecord(
        label="outer(power)",
        criteria=(_criterion(name="dk", trajectory=(1e-1, 1e-9), tolerance=1e-8),),
        budget=500,
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
            budget=2,
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
        record = IterationRecord(label="level", criteria=built, budget=budget)
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
            budget=2,
        )
        assert record.exhausted_budget is True
        assert record.converged is True
        assert record.truncated is False

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
            budget=500,
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
            budget=2,
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
            budget=500,
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

    def test_an_unlabelled_record_is_refused(self) -> None:
        with pytest.raises(ValueError, match="labelled"):
            IterationRecord(label="")

    @pytest.mark.parametrize("kw", [{"budget": -1}, {"min_iterations": -1}])
    def test_negative_counts_are_refused(self, kw) -> None:
        with pytest.raises(ValueError, match=">= 0"):
            IterationRecord(label="level", **kw)

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
            budget=200,
        )
        outer = IterationRecord(
            label="outer(power)",
            criteria=(
                _criterion(name="dk", trajectory=(1e-2, 1e-9), tolerance=1e-8),
            ),
            budget=500,
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
            budget=40,
        )
        assert "no budget suffices" in record.report()

    def test_a_converged_tree_reports_no_rate_advice(self) -> None:
        """The other half of the pair — advice that always fires is noise,
        and noise gets filtered."""
        report = _three_level_tree(deep_last=1e-13).report()
        assert "needs ~" not in report
        assert "no budget suffices" not in report
