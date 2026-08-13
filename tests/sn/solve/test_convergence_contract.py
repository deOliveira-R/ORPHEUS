"""The convergence CONTRACT: a best-effort exit is distinguishable and audible.

Issues #340 (a truncated exit was indistinguishable from a certified one at
the public entries) and #342 (``solve_sn`` hardcoded ``converged=True``).

The defect these gates exist to prevent is not a wrong number — it is a
**silently arbitrary** one.  A solve that runs out of iterations returns a
mid-descent iterate; if nothing distinguishes it from the converged answer,
a caller asserts physics against whatever the loop happened to be holding
when the budget expired.  That is how
``test_d3_pure_absorber_per_ordinate_psi_exact`` came to assert a
closed-form identity against the 999th of 1631 needed sweeps and pass, for
months, by luck (``scratch/d3_absorber_diagnosis.md``).

The contract has three parts, and each gets a gate below:

1. **The loop records why it stopped** — ``power_iteration`` returns a
   :class:`~orpheus.numerics.eigenvalue.PowerIterationOutcome` carrying
   ``converged``, instead of a triple that discarded it.  With the fact at
   the source, ``converged=True`` is not something a consumer can assert.
2. **The flag is honest at every entry** — forward and adjoint, eigenvalue
   and fixed-source, all read that fact (or the one shared residual
   predicate) rather than transcribing their own.
3. **A truncated exit is audible** — it emits
   :class:`~orpheus.numerics.convergence.ConvergenceWarning`, escalatable
   to a hard failure with
   :data:`~orpheus.numerics.convergence.ESCALATION_FLAG`.

⭐ Every gate here is written to be REDDENABLE, which for a boolean-flag
contract needs care: asserting ``converged is True`` on a solve that
converges is satisfied by the very hardcoded ``True`` this file exists to
forbid.  So each honesty gate is a PAIR — a converging configuration and a
deliberately-starved one — and it is the *starved* leg that has teeth.  A
regression to ``converged=True`` reds the starved leg of every pair.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture, make_mixture
from orpheus.geometry import BC
from orpheus.numerics.convergence import (
    ESCALATION_FLAG,
    ConvergenceWarning,
    IterationBudget,
    IterationRecord,
    StoppingCriterion,
    default_iteration_budget,
    warn_if_unconverged,
)
from orpheus.numerics.eigenvalue import PowerIterationOutcome
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import SNSolver, solve_sn, solve_sn_fixed_source
from orpheus.transport.mesh.axis import AxisMesh

_REFL = BC("reflective")


def _absorber_2g():
    """Pure absorber, 2 groups — no scattering iteration, so the reflective
    boundary coupling is the ONLY loop and its budget is the only knob."""
    return make_mixture(
        sig_t=np.array([0.8, 1.6]), sig_c=np.array([0.8, 1.6]),
        sig_f=np.array([0.0, 0.0]), nu=np.array([0.0, 0.0]),
        chi=np.zeros(2), sig_s=np.zeros((2, 2)),
    )


def _d3_reflective_axes():
    """The SI-hard configuration: all-reflective 3-D, [M] 1631 sweeps needed.

    Zero leakage plus zero scattering leaves only the boundary coupling, and
    only absorption damps it — the slow mode is the DD face sawtooth.  This
    is what makes the box a reliable *starvation* fixture: a modest
    ``max_inner`` cannot converge it, by measured construction rather than
    by a magic number.
    """
    return tuple(
        AxisMesh(edges=np.linspace(0.0, e, n + 1), bc_low=_REFL, bc_high=_REFL)
        for e, n in zip((1.0, 2.0, 3.0), (3, 4, 5))
    )


def _uniform_source(quad, cells, Q_g=(1.0, 0.5)):
    W = float(np.sum(quad.weights))
    return np.broadcast_to(
        (np.asarray(Q_g) / W)[None, :, *([None] * len(cells))],
        (quad.N, 2, *cells),
    ).copy()


def _eigenvalue(*, max_outer: int, keff_tol: float):
    """A 2-D all-reflective k-eigenvalue solve on a FISSILE mixture.

    Mixture "A" deliberately: "B" and "C" have ``SigF = 0``, so a k-solve on
    them divides by a zero production rate and yields ``nan`` — which the
    starved leg would then "pass" for entirely the wrong reason.  A control
    that is green because the physics is degenerate is not a control.
    """
    return solve_sn(
        {0: get_mixture("A", "2g")},
        tuple(
            AxisMesh(edges=np.linspace(0.0, e, n + 1),
                     bc_low=_REFL, bc_high=_REFL)
            for e, n in zip((2.0, 1.5), (4, 3))
        ),
        Quadrature.level_symmetric(sn_order=4),
        max_outer=max_outer, keff_tol=keff_tol, inner_tol=1e-10,
    )


def _fixed_source(max_inner: int, **kw):
    quad = Quadrature.level_symmetric(sn_order=4)
    return solve_sn_fixed_source(
        {0: _absorber_2g()}, _d3_reflective_axes(), quad,
        external_source=_uniform_source(quad, (3, 4, 5)),
        boundary_condition="reflective",
        inner_tol=1e-13, max_inner=max_inner, **kw,
    )


def _two_group_reflective_2d():
    """The cheap 2-D all-reflective box the #340 N6 tree fixtures share."""
    return tuple(
        AxisMesh(edges=np.linspace(0.0, e, n + 1), bc_low=_REFL, bc_high=_REFL)
        for e, n in zip((2.0, 1.5), (4, 3))
    )


def _nested_failure_eigenvalue():
    """BOTH levels fail, and every fact differs between them (#340 N6).

    The discriminating end-to-end fixture: the outer runs out of outers AND
    each inner runs out of inners, so ``first_failure`` is the child and the
    message must re-point to it.  `[M]` 2026-08-10 the seven candidate facts
    are pairwise distinct — knob (max_outer|max_inner), budget (3|20),
    tolerance (1e-12|1e-10), criterion (dphi|residual), last
    (6.308e-02|2.971e-01), rho (0.153338|0.962217), projection (17|586) — so
    each assertion in the consuming gate is individually attributable.

    `[M]` 0.18 s.  The d=3 all-reflective box discriminates identically and
    costs ~200x more; this one is the routine gate.
    """
    return solve_sn(
        {0: get_mixture("A", "2g")},
        _two_group_reflective_2d(),
        Quadrature.level_symmetric(sn_order=4),
        max_outer=3, keff_tol=1e-12, flux_tol=1e-12,
        max_inner=20, inner_tol=1e-10,
    )


def _starved_inner_converged_outer():
    """The outer MEETS its criteria while standing on 24 starved inners.

    `[M]` 2026-08-10, 0.31 s: ``converged=True``, ``fully_converged=False``,
    24 outers, 24 of 24 children TRUNCATED at 8/8.  This is #340 F1 in
    miniature — the increment-only outer stop test cannot see an upstream
    throttle, so the stall reads as convergence.

    It is the ONLY shape that separates the two guards, which is why both
    the commit-1 scope row and its commit-2 xfail sibling ride it.
    """
    return solve_sn(
        {0: get_mixture("A", "2g")},
        _two_group_reflective_2d(),
        Quadrature.level_symmetric(sn_order=4),
        max_outer=200, keff_tol=1e-8, flux_tol=1e-6,
        max_inner=8, inner_tol=1e-10,
    )


# ─── 0. The budget itself is derived, at every entry ─────────────────────


@pytest.mark.foundation
class TestEveryEntryDerivesItsBudget:
    """#340 N3: ``max_inner=None`` means *derive from the tolerance*.

    Six hardcoded constants shipped in three ``(budget, tolerance)``
    combinations, and `[M]` all of them were SHORT on the d=3 all-reflective
    absorber — the configuration whose truncated exit opened this issue.  A
    constant cannot track a tolerance it does not know about.

    These rows gate the WIRING (does each entry actually resolve?), which is
    a different claim from the law itself — that lives in
    ``tests/numerics/test_default_iteration_budget.py``.
    """

    @pytest.mark.parametrize("inner_tol", [1e-8, 1e-12])
    def test_the_solver_resolves_none_to_the_derived_budget(
        self, inner_tol: float
    ) -> None:
        from orpheus.sn.solver import SNSolver

        mesh = _sn_mesh_for_budget_probe()
        assert (
            SNSolver(mesh, inner_tol=inner_tol).max_inner
            == default_iteration_budget(inner_tol)
        )

    def test_an_explicit_budget_survives_untouched(self) -> None:
        """The deliberate-starvation path stays exactly as it was — `[M]` all
        7 truncations #340's audit found pass an explicit budget, so a
        resolution that second-guessed them would silently retune the lot."""
        from orpheus.sn.solver import SNSolver

        assert SNSolver(_sn_mesh_for_budget_probe(), max_inner=7).max_inner == 7

    def test_the_keigenvalue_posing_layer_derives_too(self) -> None:
        """The SIXTH constant, and the only one outside SN — a third
        ``(1000, 1e-8)`` pair no SN entry spelled.  #340's plan counted five
        and missed it, so it gets its own row rather than riding the others.
        """
        import inspect

        from orpheus.numerics.iteration import KEigenvalue

        assert (
            inspect.signature(KEigenvalue.__init__).parameters["max_inner"].default
            is None
        )

    @pytest.mark.parametrize(
        "entry",
        ["solve_sn", "solve_sn_adjoint", "solve_sn_fixed_source",
         "solve_sn_adjoint_fixed_source"],
    )
    def test_no_public_entry_still_ships_a_hardcoded_budget(
        self, entry: str
    ) -> None:
        """The retirement half: a constant left at ONE entry is the twin that
        reopens the whole defect, and it would be invisible to the rows above
        (which each exercise a different entry)."""
        import inspect

        import orpheus.sn.solver as solver_module

        default = (
            inspect.signature(getattr(solver_module, entry))
            .parameters["max_inner"]
            .default
        )
        assert default is None, f"{entry} still hardcodes max_inner={default}"


def _sn_mesh_for_budget_probe():
    """Smallest well-formed mesh — the budget rows care about resolution
    arithmetic, not about physics."""
    from orpheus.sn.solver import _as_sn_mesh

    return _as_sn_mesh(
        (AxisMesh(edges=np.linspace(0.0, 1.0, 3), bc_low=_REFL, bc_high=_REFL),),
        Quadrature.level_symmetric(sn_order=2),
        {0: _absorber_2g()},
        "reflective",
    )


@pytest.mark.foundation
class TestSNReportsBothCriteriaItStopsOn:
    """``SNSolver`` stops on ``dk`` AND ``dphi``; it must REPORT both.

    ⭐ This class exists because a mutation found its absence.  `[M]` patching
    :meth:`SNSolver.measure_stopping_criteria` to return only its first
    reading — which is *exactly* the pre-#340 state, where ``dphi`` was
    computed inside the predicate and discarded — reddened **zero** of 249
    tests.  Dropping the harder criterion makes convergence strictly EASIER,
    so every outcome-level gate stays green and the loss is invisible.

    The consequence is not hypothetical: with ``dphi`` unrecorded, the
    truncation warning could only ever project off ``|dk|``, and on a solve
    whose ``|dk|`` had cleared while ``dphi`` alternated in sign forever it
    answered "you need 1 more iteration".
    """

    def test_it_names_dk_and_dphi_against_their_own_tolerances(self) -> None:
        sn = _sn_mesh_for_budget_probe()
        solver = SNSolver(sn, keff_tol=1e-7, flux_tol=1e-6)
        phi = solver.initial_flux_distribution()

        readings = solver.measure_stopping_criteria(1.0, 1.0, phi, phi)

        assert [r.name for r in readings] == ["dk", "dphi"]
        assert [r.tolerance for r in readings] == [1e-7, 1e-6], (
            "each criterion must carry ITS OWN tolerance — a shared one would "
            "make `cleared` and every projection off it meaningless"
        )

    def test_each_reading_is_one_iteration_worth(self) -> None:
        """One reading per criterion per outer — what lets the loop
        concatenate them into co-indexed trajectories."""
        sn = _sn_mesh_for_budget_probe()
        solver = SNSolver(sn)
        phi = solver.initial_flux_distribution()

        readings = solver.measure_stopping_criteria(1.0, 1.0, phi, phi)

        for reading in readings:
            assert reading.n_iterations == 1, (
                f"{reading.name} reported {reading.n_iterations} values for "
                f"one iterate"
            )
        # An unchanged iterate reads zero on both — magnitudes, never signed.
        assert [r.last for r in readings] == [0.0, 0.0]

    def test_the_readings_are_MAGNITUDES_of_a_real_step(self) -> None:
        """A genuine step reports strictly positive, finite readings.

        The positive control for the two gates above: without it they are
        satisfied by a method pinned at zero, which would report two
        beautifully-named criteria that never move.
        """
        sn = _sn_mesh_for_budget_probe()
        solver = SNSolver(sn)
        phi = solver.initial_flux_distribution()

        readings = solver.measure_stopping_criteria(1.05, 1.0, phi, phi * 0.5)

        for reading in readings:
            assert reading.last is not None
            assert 0.0 < reading.last < np.inf, (
                f"{reading.name} read {reading.last!r} on a real step"
            )


# ─── 1. The primitive records the fact ───────────────────────────────────


@pytest.mark.foundation
class TestPowerIterationCarriesItsOutcome:
    """``power_iteration`` returns WHY it stopped, not just what it found."""

    def test_outcome_is_not_tuple_unpackable(self) -> None:
        """The old ``keff, hist, flux = power_iteration(...)`` idiom is DEAD.

        Not pedantry: that idiom is exactly how the flag got discarded at
        five call sites.  If the outcome destructured like the triple it
        replaced, every one of them would still compile while ignoring
        ``converged``, and the retirement would be cosmetic.
        """
        outcome = PowerIterationOutcome(
            keff=1.0, keff_history=[1.0], flux_distribution=np.ones(3),
            record=IterationRecord(label="outer(probe)"),
        )
        with pytest.raises(TypeError):
            _a, _b, _c = outcome           # type: ignore[misc]

    def test_the_outcome_cannot_be_HANDED_a_convergence_claim(self) -> None:
        """``converged`` is derived from the record, so it is not a ctor arg.

        The sharper half of #342's fix, and the reason this gate exists
        separately from the value gates below: the defect was a hand-written
        ``converged=True``, and a keyword that no longer EXISTS cannot be
        hand-written.  A gate that only checked the value would stay green
        against a re-introduced field defaulting to the optimistic answer.
        """
        with pytest.raises(TypeError):
            PowerIterationOutcome(
                keff=1.0, keff_history=[1.0], flux_distribution=np.ones(3),
                record=IterationRecord(label="outer(probe)"),
                converged=True,  # type: ignore[call-arg]
            )

    def test_a_starved_budget_reports_not_converged(self) -> None:
        """max_outer=1 cannot satisfy a criterion needing >=3 outers.

        The teeth of the whole file: ``MINIMUM_OUTER_ITERATIONS`` is 3, so
        this outcome is structurally unconverged and any hardcoded ``True``
        anywhere on the path reds here.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            sol = _eigenvalue(max_outer=1, keff_tol=1e-10)
        assert sol.history is not None
        assert sol.history.converged is False

    def test_a_generous_budget_reports_converged(self) -> None:
        """The positive control — the flag is not simply pinned False."""
        sol = _eigenvalue(max_outer=200, keff_tol=1e-8)
        assert sol.history is not None
        assert sol.history.converged is True


# ─── 2. The flag is honest at the fixed-source entry ─────────────────────


@pytest.mark.foundation
class TestFixedSourceFlagIsHonest:
    """The pair that would have caught the d3 defect the day it appeared."""

    def test_starved_budget_reports_not_converged(self) -> None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            sol = _fixed_source(max_inner=50)
        assert sol.history is not None
        assert sol.history.converged is False

    def test_sufficient_budget_reports_converged(self) -> None:
        """[M] 1631 sweeps needed; 4000 is ~2.5x headroom."""
        sol = _fixed_source(max_inner=4000)
        assert sol.history is not None
        assert sol.history.converged is True


# ─── 3. A truncated exit is AUDIBLE ──────────────────────────────────────


@pytest.mark.foundation
class TestTruncationIsAudible:
    """The warning exists, fires exactly when the flag is False, and says
    enough to act on."""

    def test_starved_solve_warns(self) -> None:
        with pytest.warns(ConvergenceWarning) as record:
            _fixed_source(max_inner=50)
        assert len(record) >= 1

    def test_converged_solve_is_SILENT(self) -> None:
        """The other half of the pair — a warning that always fires is noise,
        and noise gets filtered, which is how the next truncation goes
        unnoticed."""
        with warnings.catch_warnings():
            warnings.simplefilter("error", ConvergenceWarning)
            _fixed_source(max_inner=4000)   # must not raise

    def test_the_message_carries_budget_tolerance_and_distance(self) -> None:
        """A bare "did not converge" cannot be acted on.

        The reader needs to tell "one more sweep" from "diverging", so the
        message must name the budget that ran out, the tolerance that was
        missed, and how far the last iterate actually was.
        """
        with pytest.warns(ConvergenceWarning) as record:
            _fixed_source(max_inner=50)
        msg = str(record[0].message)
        assert "max_inner=50" in msg
        assert "tol=" in msg
        assert "last residual" in msg
        assert "BEST-EFFORT" in msg

    def test_the_message_names_the_BUDGET_TO_SET_not_just_a_direction(
        self,
    ) -> None:
        """⭐ #340 N3: "raise max_inner" is not advice, it is a direction.

        Against a ``rho = 0.985`` mode a reader's guess is wrong by a factor
        of six, so the message must carry the OBSERVED rate and the budget
        that rate projects.  The number is checked against the projection the
        same trajectory yields, so this gate cannot be satisfied by printing
        an arbitrary integer.

        ⚠ It reads LOW from inside the transient, and says "so far" for that
        reason — `[M]` on this fixture at budget 50 it projects 586 against a
        true 1631, sharpening to 1618 by budget 200 (0.8 %).  It converges
        from BELOW, so following the advice yields a bigger number next time
        rather than a wrong answer.
        """
        with pytest.warns(ConvergenceWarning) as record:
            sol = _fixed_source(max_inner=50)
        msg = str(record[0].message)

        assert sol.history is not None
        projected = StoppingCriterion(
            name="residual",
            trajectory=tuple(sol.history.flux_residuals),
            tolerance=1e-13,
        ).projected_iterations()
        assert projected is not None
        assert f"set max_inner={projected}" in msg
        assert "rho=" in msg

    def test_a_NON_CONTRACTING_solve_is_told_no_budget_will_help(self) -> None:
        """The other arm, and it must not be the same advice.

        `[M]` #340 found a configuration whose NEGATIVE dominant eigenvalue
        makes the increment sign-alternate, so its stop test is unsatisfiable
        *forever*: telling that reader to raise the budget sends them down a
        road with no end.  Driven here through the shared warning helper with
        a synthetic non-decaying history, because the arm is a property of
        the MESSAGE, not of any one solver configuration.
        """
        stalled = IterationRecord(
            label="inner(probe)",
            criteria=(
                StoppingCriterion(
                    name="residual",
                    trajectory=tuple(1e-3 for _ in range(40)),  # rho == 1
                    tolerance=1e-10,
                ),
            ),
            budget=IterationBudget(40, "max_inner"), iterations_run=40,
        )
        with pytest.warns(ConvergenceWarning) as record:
            warn_if_unconverged(stalled, where="probe")
        msg = str(record[0].message)
        assert "NOT" in msg and "contracting" in msg
        assert "set max_inner=" not in msg  # the advice must NOT be given

    def test_it_speaks_for_the_criterion_that_FAILED_not_one_that_cleared(
        self,
    ) -> None:
        """⛔ The regression for a lie this gate's own author shipped — and
        the SUCCESSOR to the stop-gap that first contained it.

        `[M]` 2026-08-09, caught by the wide run, not by review. The outer
        stops on ``dk AND dphi``, but only ``keff_history`` survived into the
        flat ``IterationHistory`` (#340 F2 — ``dphi`` was computed inside
        ``SNSolver.converged`` and discarded). On the mutated heterogeneous
        slab, whose NEGATIVE dominant eigenvalue makes ``dphi`` sign-alternate
        forever, ``|dk|`` sits at ``3.3e-16`` against a ``1e-9`` tolerance —
        so the freshly-added projection answered **"set max_outer=1"**.

        ⛔ The first fix was a REFUSAL: detect that the recorded criterion had
        cleared and decline to project ("no budget can honestly be projected
        from what is here"). That branch is RETIRED with N2b-ii, and this test
        is re-posed onto the successor rather than deleted — the situation it
        guarded (*the failing criterion is absent from the history*) is now
        unrepresentable, because the record carries every criterion the level
        was judged on.

        So the claim strengthens from "say nothing rather than something
        wrong" to **"say the right thing"**: ``binding_criterion`` is the one
        furthest from clearing, which IS the one that failed, and the message
        must speak for it.
        """

        # dk decayed to ~1.8e-15 against 1e-9 — cleared long ago.
        # dphi sits at 1e-3 against 1e-6 forever — rho == 1, never clears.
        record_ = IterationRecord(
            label="outer(power-iteration)",
            criteria=(
                StoppingCriterion(
                    name="dk",
                    trajectory=tuple(1e-3 * 0.5**k for k in range(40)),
                    tolerance=1e-9,
                ),
                StoppingCriterion(
                    name="dphi",
                    trajectory=tuple(1e-3 for _ in range(40)),
                    tolerance=1e-6,
                ),
            ),
            budget=IterationBudget(40, "max_outer"), iterations_run=40,
            min_iterations=3,
        )
        assert record_.converged is False, "fixture drift"

        with pytest.warns(ConvergenceWarning) as caught:
            warn_if_unconverged(record_, where="probe")
        msg = str(caught[0].message)

        assert "last dphi" in msg, (
            "the message must name the criterion that FAILED; naming the "
            "cleared one is how the 'set max_outer=1' lie was produced"
        )
        # ⭐ #340 N6 — and it must quote THAT criterion's tolerance.  Until
        # 2026-08-10 the tolerance was a hand-passed argument describing the
        # ENTRY, so this very message printed `tol=1.000e-09` (dk's) beside
        # `last dphi` (whose tolerance is 1e-6) — a mismatched pair that no
        # assertion here caught, because nothing asserted `tol=` at all.  The
        # negative leg is the one with teeth: it fails the moment the
        # tolerance stops being read off the binding criterion.
        assert "tol=1.000e-06" in msg, (
            "the tolerance must be the BINDING criterion's (dphi @ 1e-6), "
            "not whatever the caller happened to pass"
        )
        assert "1.000e-09" not in msg, (
            "dk's tolerance belongs to the criterion that CLEARED; quoting "
            "it beside 'last dphi' is the mismatched pair N6 removed"
        )
        # dphi is not contracting, so there is still no budget to promise —
        # but now that is a statement ABOUT dphi rather than a refusal to
        # speak at all.
        assert "NOT" in msg and "contracting" in msg
        assert "set max_outer=" not in msg
        assert "ALREADY cleared" not in msg, (
            "the retired stop-gap's wording must not survive its branch"
        )

    def test_it_names_the_LOOP_when_every_criterion_cleared(self) -> None:
        """The branch that replaced the stop-gap, and it is a different claim.

        A level can fail to converge with every criterion cleared — the loop
        refused, not a quantity: too few iterations to claim
        (``min_iterations``). That is a real state, it is not a budget
        problem, and naming it beats projecting a number at it.
        """

        record_ = IterationRecord(
            label="outer(power-iteration)",
            criteria=(
                StoppingCriterion(
                    name="dk", trajectory=(1e-12, 1e-14), tolerance=1e-9,
                ),
            ),
            budget=IterationBudget(2, "max_outer"), iterations_run=2,
            min_iterations=3,
        )
        assert record_.converged is False and record_.criteria[0].cleared

        with pytest.warns(ConvergenceWarning) as caught:
            warn_if_unconverged(record_, where="probe")
        msg = str(caught[0].message)
        assert "the refusal is the loop's" in msg
        assert "2 of the 3 iterations" in msg
        assert "set max_outer=" not in msg

    # ── #340 N6: the message speaks for the level that FAILED ────────────

    def test_the_message_speaks_for_the_FAILING_LEVEL_not_the_entry(
        self,
    ) -> None:
        """⭐ THE KEYSTONE of #340 N6, and the only gate that can see it.

        `[M]` 2026-08-10: all four audibility fixtures above route through
        ``solve_sn_fixed_source``, which is a **1-level tree**, so
        ``first_failure is record`` holds by *object identity* and six of the
        seven facts below are one object's attributes read twice.  No input
        can make the two sourcings differ there — the leaf rows are
        annihilated for this claim, not merely under-tested.  A synthetic
        two-level record is the cheapest shape that separates them.

        The claim is that EVERY fact is read off ``record.first_failure`` —
        not just the three that used to be caller-passed.  A partial
        re-point is *worse* than none: it welds the inner's knob and budget
        to the OUTER's criterion, rate and projection, which is a message
        where every number is real, every pairing is wrong, and the result
        reads as level-correct.  Hence the paired legs — each fact asserted
        PRESENT with the failing level's value and ABSENT with the top's.

        Trajectories are exact geometrics so ``rho`` and the projection are
        analytically determined (0.9 / 0.5) rather than solve-dependent; the
        assertions still RECOMPUTE them from the record, because a literal
        ``189`` would be a false red the day the rate-fit tail changes.
        """

        inner = IterationRecord(
            label="inner(within-group)",
            criteria=(
                StoppingCriterion(
                    name="residual",
                    trajectory=tuple(3.7e-2 * 0.9**k for k in range(30)),
                    tolerance=1e-10,
                ),
            ),
            budget=IterationBudget(30, "max_inner"), iterations_run=30,
        )
        outer = IterationRecord(
            label="outer(power-iteration)",
            criteria=(
                StoppingCriterion(
                    name="dk",
                    trajectory=tuple(1e-2 * 0.5**k for k in range(6)),
                    tolerance=1e-9,
                ),
                StoppingCriterion(
                    name="dphi",
                    trajectory=tuple(1e-3 * 0.5**k for k in range(6)),
                    tolerance=1e-7,
                ),
            ),
            budget=IterationBudget(6, "max_outer"), iterations_run=6,
            min_iterations=3, children=(inner,),
        )

        # ── leg 0: ACTIVATION.  Every fact must DIFFER between the two
        # levels, or the leg asserting it silently stops discriminating and
        # goes green forever (Mode 12, at the fixture rather than the
        # functional).  A fixture edit that collapses any pair reds HERE.
        assert outer.first_failure is inner, "fixture drift: the inner must fail"
        top, failing = outer.binding_criterion, inner.binding_criterion
        assert top is not None and failing is not None
        for fact, pair in {
            "knob": (outer.budget.name, inner.budget.name),
            "budget": (outer.budget, inner.budget),
            "tolerance": (top.tolerance, failing.tolerance),
            "criterion": (top.name, failing.name),
            "last value": (top.last, failing.last),
            "rate": (outer.rate, inner.rate),
            "projection": (outer.projected_iterations(),
                           inner.projected_iterations()),
        }.items():
            assert pair[0] != pair[1], (
                f"fixture drift: outer and failing inner agree on {fact} "
                f"({pair[0]!r}); the leg asserting it cannot discriminate"
            )

        with pytest.warns(ConvergenceWarning) as caught:
            warn_if_unconverged(outer, where="solve_sn")
        msg = str(caught[0].message)

        # level + knob + budget, in one substring so a mixed sourcing cannot
        # satisfy it by accident.
        assert f"{inner.label} hit {inner.budget}" in msg
        assert "max_outer" not in msg, "the top level's knob must not appear"
        assert outer.label not in msg, "the top level is not the failing one"

        assert f"tol={failing.tolerance:.3e}" in msg
        assert f"tol={top.tolerance:.3e}" not in msg

        assert f"last {failing.name}" in msg
        assert f"last {top.name}" not in msg

        assert f"rho={inner.rate:.6f}" in msg
        assert f"rho={outer.rate:.6f}" not in msg

        # ⭐ ``covering`` and not the raw projection (#349): the projection is
        # in TRAJECTORY units and the knob may not take those.  Here the inner
        # is a SourceIteration, so the conversion is the identity and this
        # assertion reads exactly as it did before — which is the point.  The
        # unit-carrying arm is gated directly on the type, in
        # ``tests/numerics/test_iteration_record.py``.
        projected = inner.projected_iterations()
        assert projected is not None, (
            "the fixture must project — an unprojectable level takes the "
            "'no budget suffices' arm and this assertion would be vacuous"
        )
        setting = inner.budget.covering(projected)
        assert f"set {inner.budget.name}={setting}" in msg

        # ``where`` is the ONE thing the caller still supplies, and it names
        # the ENTRY, not a level — it must NOT re-point.
        assert msg.startswith("solve_sn:")

    def test_a_LEAF_message_names_no_level_at_all(self) -> None:
        """The other half of the pair: on a 1-level tree nothing is prefixed.

        The level prefix exists to say "this is not the level you asked
        about".  On a leaf the failing level IS the record, so prefixing it
        would be noise — and the discriminator is an ``is`` identity, which
        no value comparison can stand in for.  Without this row the prefix
        could be emitted unconditionally and the whole suite would stay
        green (the tree row above asserts only that it IS present).
        """
        with pytest.warns(ConvergenceWarning) as caught:
            sol = _fixed_source(max_inner=50)
        msg = str(caught[0].message)

        assert sol.history is not None
        record = sol.history.record
        assert record.first_failure is record, "fixture drift: not a leaf"
        assert msg.startswith("solve_sn_fixed_source: hit "), (
            "a leaf has no level to disambiguate, so the message goes "
            f"straight to the budget; got {msg[:60]!r}"
        )
        assert record.label not in msg

    def test_a_level_with_NO_criteria_quotes_no_tolerance(self) -> None:
        """A level that measures nothing has no tolerance to name.

        MoC's inner is a fixed sweep count — no tolerance, no residual — so
        ``binding_criterion`` is ``None`` and the ``without reaching tol=``
        clause is DROPPED rather than faked with a zero.  Faking it would
        put a number in the reader's hands that no code ever compared
        anything against.
        """

        bare = IterationRecord(
            label="inner(moc-sweeps)", budget=IterationBudget(4, "max_inner"),
            iterations_run=4,
        )
        assert bare.binding_criterion is None and bare.converged is False

        with pytest.warns(ConvergenceWarning) as caught:
            warn_if_unconverged(bare, where="probe")
        msg = str(caught[0].message)

        assert "hit max_inner=4 (no criterion recorded)" in msg
        assert "tol=" not in msg, "there is no tolerance to quote"
        assert "rho=" not in msg, "there is no trajectory to fit"

    def test_a_real_nested_solve_reports_its_STARVED_INNER(self) -> None:
        """The same claim, end to end — what the synthetic row cannot see.

        The keystone above proves the CONSUMER's routing on a record built
        by hand.  It is structurally blind to whether production actually
        builds a nested record and stamps the caller-facing knob on the
        CHILD (``vv`` Mode 11).  This row runs the real eigenvalue entry and
        navigates the reference by hand through ``record.children[0]`` — a
        different route into the same object than the SUT's
        ``first_failure``, which is exactly the claim.

        `[M]` 2026-08-10, 0.18 s: outer TRUNCATED 3/3 binding ``dphi`` @
        1e-12 (rho 0.153338, projects 17); all three children TRUNCATED
        20/20 binding ``residual`` @ 1e-10 (rho 0.962217, projects 586).
        Deterministic — ``max_outer=3`` equals ``MINIMUM_OUTER_ITERATIONS``,
        so the loop always runs exactly three outers.

        ⚠ Asserts NOTHING about ``keff``: three starved inners make it an
        arbitrary mid-descent number (`[M]` 1.874, against a true 1.875).
        """
        with pytest.warns(ConvergenceWarning) as caught:
            sol = _nested_failure_eigenvalue()
        msg = str(caught[0].message)

        assert sol.history is not None
        record = sol.history.record
        child = record.children[0]
        assert record.first_failure is child, "fixture drift: the inner must fail"
        assert record.binding_criterion is not None

        assert f"{child.label} hit {child.budget}" in msg
        assert "max_outer" not in msg
        criterion = child.binding_criterion
        assert criterion is not None
        assert f"tol={criterion.tolerance:.3e}" in msg
        assert f"last {criterion.name}" in msg
        assert f"set max_inner={child.projected_iterations()}" in msg
        assert msg.startswith("solve_sn:")

    def test_the_tree_wide_truncation_is_audible(self) -> None:
        """⭐ THE #340 HEADLINE: a converged outer on a starved inner SPEAKS.

        This is F1 in one row.  The outer meets its own criteria, so the
        flat ``converged`` reads ``True`` and every increment-based stop
        test is satisfied — while 24 of 24 inners sat at their cap.  Until
        2026-08-10 that solve returned in total silence; the guard is now
        ``fully_converged``, so it does not.

        ⛔ **Three assertions, and the last two are why this row is not
        tautological.** ``pytest.warns`` alone would also pass if the
        fixture drifted so the OUTER stopped converging — in which case the
        OLD guard would have warned too and this row would be testing
        nothing (`vv-principles` #19: a positive reading that cannot
        discriminate).  The fixture-integrity pair below is inherited from
        the retired ``…_is_STILL_SILENT`` sibling, which is where it was
        load-bearing first; deleting that row without carrying them here
        would have quietly demoted this one.

        ⭐ And the third claim — that the message names the INNER — is what
        separates *a* warning from the *right* warning.  Its sibling
        ``test_a_real_nested_solve_reports_its_STARVED_INNER`` makes the
        same naming claim, but on a fixture whose outer ALSO failed, so the
        old guard already warned there.  Here the warning exists only
        because of the widening, which makes this the one row where the
        naming and the flip are tested together.

        `[M]` 2026-08-10: ``converged=True``, ``fully_converged=False``,
        24 outers, 24 of 24 children truncated at 8/8.
        """
        with pytest.warns(ConvergenceWarning) as caught:
            probe = _starved_inner_converged_outer()

        assert probe.history is not None
        assert probe.history.converged is True, (
            "fixture drift: the outer must CONVERGE, or this row stops "
            "testing the widened guard — the old one would have warned too"
        )
        assert probe.history.fully_converged is False, (
            "fixture drift: an inner must be starved, or this row "
            "degenerates into test_converged_solve_is_SILENT"
        )

        failing = probe.history.record.first_failure
        assert failing is not None and failing is not probe.history.record, (
            "fixture drift: the failure must be a CHILD, not the top level"
        )
        msg = str(caught[0].message)
        assert f"{failing.label} hit {failing.budget}" in msg
        assert "max_outer" not in msg, (
            "the advice must not point at the outer's knob — with a starved "
            "inner, raising max_outer cannot help at all (#340 F2)"
        )
        # ⛔ This asserted the literal ``solution.history.fully_converged``
        # until 2026-08-11 (#340 N4.7).  That string guessed the CALLER's local
        # variable name — a fact no library can know — and it was outright
        # wrong for the three families whose entries return a ``*Result``.  The
        # CLAIM is unchanged and is what this row still tests: the reader is
        # sent to the predicate that can SEE a starved inner, never to the flat
        # ``converged`` (which reads True on exactly this solve).
        assert "`fully_converged`" in msg and "IterationRecord" in msg, (
            "the message must send the reader to the predicate that can SEE "
            "this — the flat `converged` reads True on exactly this solve"
        )
        assert "solution.history" not in msg, (
            "the advice must not guess the caller's variable name; three of "
            "the four families have no `solution` at all (#340 N4.7)"
        )

    def test_the_level_facts_are_GONE_from_the_warning_signature(self) -> None:
        """The retirement half — a defaulted survivor is the twin.

        ``budget_name`` / ``budget`` / ``tol`` described the CALLER's top
        level.  Left on the signature with defaults, every call site would
        still compile while one of them quietly kept passing its own knob,
        and the re-sourcing would be cosmetic at that entry.  Set equality,
        not membership: it also catches a fact re-added under a new name.

        ⭐ Re-posed 2026-08-11 (#340 N4.7) when the helper moved to
        :mod:`orpheus.numerics.convergence` and became family-agnostic.  The
        parameter set grew by one, and the CLAIM had to grow with it or the
        set-equality would have been a false red.  The claim is now two-part
        and the second half is what keeps the first honest:

        1. Nothing a LEVEL knows may be a parameter — those come off
           ``record.first_failure``.
        2. Every parameter names something the record structurally CANNOT
           know: ``where`` (the entry, chosen by the caller) and
           ``balance_defect`` (a property of the RETURNED ITERATE, not of the
           iteration — see the helper's docstring for why it is deliberately
           not a record field).

        ⚠ So ``balance_defect`` is admitted by the second criterion, not
        grandfathered.  A future parameter that fails it — ``budget=``,
        ``tol=``, ``label=`` under any spelling — is the retired defect
        returning, which is exactly what the frozen set catches.
        """
        import inspect

        parameters = set(inspect.signature(warn_if_unconverged).parameters)
        assert parameters == {"record", "where", "balance_defect"}, (
            "the failing level's facts are read off the record, so the only "
            "things a caller still supplies are WHERE it was called and a "
            f"diagnostic about the returned iterate; got {sorted(parameters)}"
        )
        # The level-fact names, spelled out so the pin fails LOUDLY rather
        # than just differing by a set element nobody reads.
        level_facts = {
            "budget", "budget_name", "tol", "tolerance", "label",
            "min_iterations", "n_iterations", "iterations_run", "criterion",
            "criteria", "rate", "projected", "history",
        }
        assert not (parameters & level_facts), (
            "a level's own fact came back onto the signature: "
            f"{sorted(parameters & level_facts)} — read it off "
            "`record.first_failure` instead (#340 F2/N6a)"
        )

    def test_it_is_escalatable_to_an_error(self) -> None:
        """Prove the CATEGORY escalates rather than merely being emitted.

        ⚠ This proves the mechanism, NOT the published recipe — it installs
        the filter through the in-process API, so it is green for any
        spelling whatsoever.  The string is a separate claim and has its own
        gate below; read the two together.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("error", ConvergenceWarning)
            with pytest.raises(ConvergenceWarning):
                _fixed_source(max_inner=50)

    def test_the_published_escalation_flag_actually_parses(self) -> None:
        """The CI contract is a STRING, and a string can be unspellable.

        ⛔ 2026-08-09 (#340).  The recipe shipped as
        ``-W error::ConvergenceWarning`` at four sites — including inside
        the emitted warning message — and **it does not parse**: ``-W``
        resolves an undotted category against ``builtins``, so pytest dies
        at startup with ``AttributeError: module 'builtins' has no
        attribute 'ConvergenceWarning'`` and collects ZERO tests.  The CI
        contract was imaginary for as long as it was published, and the
        sibling gate above stayed green throughout because it never touched
        the string.

        So this gate consumes the STRING, through pytest's own parser.
        """
        from _pytest.config import UsageError, parse_warning_filter

        assert ESCALATION_FLAG.startswith("-W ")
        spec = ESCALATION_FLAG.removeprefix("-W ")

        # It parses, and resolves to the category we actually emit.
        assert parse_warning_filter(spec, escape=False)[2] is ConvergenceWarning

        # And the teeth: the short spelling that shipped must still be
        # rejected, so a future "simplification" back to it reds here
        # instead of silently disarming CI.
        with pytest.raises(UsageError):
            parse_warning_filter("error::ConvergenceWarning", escape=False)


# ─── 4. The advice names a knob the reader can actually type ─────────────


def _tiny_slab_inputs():
    """The smallest well-formed 1-D problem, shared by the knob sweep.

    Four cells, GL-4, two groups — `[M]` 0.01-0.20 s per entry.  The knob is
    stamped on every record whether or not the level truncates, so these
    configurations deliberately CONVERGE: the claim is about the label the
    producer wrote, not about a warning.
    """
    axes = (
        AxisMesh(edges=np.linspace(0.0, 1.0, 5), bc_low=_REFL, bc_high=_REFL),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    source = np.broadcast_to(
        (np.array([1.0, 0.5]) / float(np.sum(quad.weights)))[None, :, None],
        (quad.N, 2, 4),
    ).copy()
    return axes, quad, source


def _fixed_source_slab(inner_solver: str):
    axes, quad, source = _tiny_slab_inputs()
    return solve_sn_fixed_source(
        {0: _absorber_2g()}, axes, quad, external_source=source,
        boundary_condition="reflective", inner_solver=inner_solver,
    )


def _eigenvalue_slab(inner_solver: str):
    axes, quad, _ = _tiny_slab_inputs()
    return solve_sn(
        {0: get_mixture("A", "2g")}, axes, quad, inner_solver=inner_solver,
    )


def _adjoint_eigenvalue_slab():
    from orpheus.sn.solver import solve_sn_adjoint

    axes, quad, _ = _tiny_slab_inputs()
    return solve_sn_adjoint({0: get_mixture("A", "2g")}, axes, quad)


def _adjoint_fixed_source_slab():
    """⚠ This entry takes a DETECTOR RESPONSE shaped ``(ng, *spatial)`` —
    NOT the angular ``(N, ng, *spatial)`` its forward sibling takes.
    `[M]` 2026-08-10, passing the angular shape raises ``ValueError``."""
    from orpheus.sn.solver import solve_sn_adjoint_fixed_source

    axes, quad, _ = _tiny_slab_inputs()
    return solve_sn_adjoint_fixed_source(
        {0: _absorber_2g()}, axes, quad,
        detector_response=np.ones((2, 4)),
        boundary_condition="reflective",
    )


# The STARVED siblings of the two adjoint slabs.  The rows above deliberately
# CONVERGE (they gate the knob a producer stamped, not a warning); the
# attribution gates need the same entries with a level held open, because a
# silent solve says nothing about where a warning would have been blamed.
def _starved_adjoint_eigenvalue():
    from orpheus.sn.solver import solve_sn_adjoint

    axes, quad, _ = _tiny_slab_inputs()
    return solve_sn_adjoint(
        {0: get_mixture("A", "2g")}, axes, quad,
        max_outer=2, keff_tol=1e-14,
    )


def _starved_adjoint_fixed_source():
    from orpheus.sn.solver import solve_sn_adjoint_fixed_source

    axes, quad, _ = _tiny_slab_inputs()
    return solve_sn_adjoint_fixed_source(
        {0: _absorber_2g()}, axes, quad,
        detector_response=np.ones((2, 4)),
        boundary_condition="reflective",
        max_inner=2, inner_tol=1e-13,
    )


#: ``(id, run the entry, the entry itself, the knob its TOP level must name)``.
#:
#: The entry is carried as the FUNCTION, not a name to ``getattr`` and not a
#: tag to branch on: the row IS the call, so the knob-membership reference
#: below reads the signature of the very callable the row invoked, and every
#: argument is spelled explicitly rather than splatted from an untyped dict
#: (``coding-elegance`` anti-#4 — a tag plus a ``**kw`` splat is stringly-typed
#: dispatch, and pyright cannot check which parameter it lands on).
#:
#: `[M]` the fixed-source entries expose only ``max_inner``; the eigenvalue
#: entries expose BOTH ``max_outer`` and ``max_inner`` — which is exactly why
#: the membership leg alone cannot catch a SWAP, and why the second row below
#: exists.  Measured: under a swap, membership reds only the two fixed-source
#: rows; the four eigenvalue rows are caught by the exact-expectation leg
#: alone.
_KNOB_ROWS = [
    ("fixed_source-si",
     lambda: _fixed_source_slab("source_iteration"),
     solve_sn_fixed_source, "max_inner"),
    ("fixed_source-krylov",
     lambda: _fixed_source_slab("krylov"),
     solve_sn_fixed_source, "max_inner"),
    ("eigenvalue-si",
     lambda: _eigenvalue_slab("source_iteration"), solve_sn, "max_outer"),
    ("eigenvalue-krylov",
     lambda: _eigenvalue_slab("krylov"), solve_sn, "max_outer"),
    ("adjoint_eigenvalue", _adjoint_eigenvalue_slab, None, "max_outer"),
    ("adjoint_fixed_source", _adjoint_fixed_source_slab, None, "max_inner"),
]


def _entry_of(row_id: str, declared):
    """The entry function a row drives — resolved late for the adjoints.

    The two adjoint entries are imported inside their runners (the module
    keeps its import surface small), so their rows declare ``None`` and the
    function is fetched here rather than at table-definition time.
    """
    if declared is not None:
        return declared
    import orpheus.sn.solver as solver_module

    return getattr(
        solver_module,
        "solve_sn_adjoint" if row_id == "adjoint_eigenvalue"
        else "solve_sn_adjoint_fixed_source",
    )


@pytest.mark.foundation
class TestEveryEntryStampsItsCallerFacingKnob:
    """#340 N6: ``budget_name`` names a parameter of the entry that emitted it.

    The knob is the one fact in the warning that the record could not
    already answer, so it is the one a producer can silently FORGET — and
    the default (``"max_iter"``) is a plausible-looking string that is a
    parameter of **no** public SN entry.  A forgotten stamp therefore ships
    advice telling the reader to set something that does not exist, and no
    value gate anywhere can see it.

    The reference is ``inspect.signature`` of the public entry — the API
    itself, structurally independent of the record under test.
    """

    @pytest.mark.parametrize(
        "row_id,run,declared_entry,expected", _KNOB_ROWS,
        ids=[row[0] for row in _KNOB_ROWS],
    )
    def test_every_level_names_a_knob_this_entry_actually_has(
        self, row_id: str, run, declared_entry, expected: str
    ) -> None:
        import inspect

        solution = run()
        assert solution.history is not None
        entry = _entry_of(row_id, declared_entry)
        knobs = set(inspect.signature(entry).parameters)

        for level in solution.history.record.walk():
            assert level.budget.name in knobs, (
                f"{entry.__name__} returned a level advising `set "
                f"{level.budget.name}=...`, which is not one of its "
                f"parameters — the reader has nothing to type"
            )

    @pytest.mark.parametrize(
        "row_id,run,declared_entry,expected", _KNOB_ROWS,
        ids=[row[0] for row in _KNOB_ROWS],
    )
    def test_the_knob_is_the_RIGHT_one_at_each_depth(
        self, row_id: str, run, declared_entry, expected: str
    ) -> None:
        """The membership leg above is blind to a SWAP; this one is not.

        `[M]` 2026-08-10 ``solve_sn`` exposes BOTH ``max_outer`` and
        ``max_inner``, so a producer that stamped them the wrong way round
        satisfies "is a parameter of this entry" perfectly while sending
        every eigenvalue reader to the knob that cannot help them.  Measured
        under exactly that swap: membership reddened only the two
        fixed-source rows (whose entries have no ``max_outer`` to be wrong
        about); all four eigenvalue rows were caught here and nowhere else.
        """
        solution = run()
        assert solution.history is not None
        record = solution.history.record

        assert record.budget.name == expected, (
            f"{row_id}: the top level is governed by {expected}, not "
            f"{record.budget.name!r}"
        )
        for child in record.children:
            assert child.budget.name == "max_inner", (
                f"{row_id}: an inner level is governed by max_inner, not "
                f"{child.budget.name!r}"
            )


# ─── The exit balance projection (#340 N6b step 2) ───────────────────────


@pytest.mark.foundation
class TestExitBalanceDefect:
    r"""The number the warning reports, and what it is allowed to claim.

    :math:`R_g = \int_V \int_{4\pi} (A\psi - q)_g \, d\Omega \, dV`,
    reported as :math:`\lVert R_g \rVert / \lVert R_g(q) \rVert`.  It exists
    because the RAW residual cannot discriminate — `[M]` #340 N5 measured a
    634× overlap between benign and corrupting truncations, missing 15 of 16
    — while projecting onto the functional :math:`k` actually reads cuts
    that to 4.64×.

    ⚠ **4.64× is still an overlap, so nothing here asserts a THRESHOLD.**
    These rows pin the arithmetic, the plumbing and the message; the one
    thing they deliberately never do is assert that a defect above some
    value means the answer is wrong.  A gate that did would be the refuted
    N5 wearing a different name.
    """

    def test_the_projection_is_angle_then_volume(self) -> None:
        """Arithmetic, against a reference written independently here.

        Production composes two canonical reductions
        (``AngularField._integrate_angular_values`` then
        ``MaterialMesh.integrate_per_group``); this row re-expresses the
        same definition as one explicit ``einsum`` over the weights and
        volumes, so a transposed axis or a dropped measure has nowhere to
        hide.  Not a tolerance check — the two spellings are the same
        arithmetic and are expected to agree to a few ULP.
        """
        from orpheus.sn.solver import _balance_projection

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            sol = _starved_inner_converged_outer()
        sn_mesh = sol.mesh
        field = sol.angular_flux

        produced = _balance_projection(field, sn_mesh=sn_mesh)

        values = np.asarray(field.interior.values)      # (N, ng, *spatial)
        weights = np.asarray(sn_mesh.quad.weights)      # (N,)
        volumes = np.asarray(sn_mesh.volumes)           # spatial
        # An explicit double loop, not a vectorised restatement: the thing
        # most likely to be wrong is an AXIS, and a loop indexed by name
        # cannot inherit production's axis convention the way a clever
        # `einsum` can.
        reference = np.array([
            sum(
                float(weights[n]) * float(np.sum(values[n, g] * volumes))
                for n in range(values.shape[0])
            )
            for g in range(sn_mesh.ng)
        ])
        assert produced.shape == (sn_mesh.ng,)
        np.testing.assert_allclose(
            produced, reference, rtol=1e-13,
            err_msg="R_g must be Σ_n w_n Σ_i V_i x[n, g, i]",
        )

    def test_a_truncated_eigenvalue_solve_CARRIES_the_number(self) -> None:
        """End to end: the field is populated and the message quotes it.

        The two halves are one row on purpose — a populated field the
        message never prints, or a printed number the caller cannot read
        back, would each satisfy half the deliverable and neither is the
        deliverable.
        """
        with pytest.warns(ConvergenceWarning) as caught:
            sol = _starved_inner_converged_outer()

        assert sol.history is not None
        defect = sol.history.balance_defect
        assert defect is not None, (
            "a truncated eigenvalue exit must carry the balance defect"
        )
        assert np.isfinite(defect) and defect >= 0.0

        msg = str(caught[0].message)
        assert f"{defect:.3e}" in msg, (
            "the message must quote the number the history carries — not a "
            "separately-formatted one that could drift from it"
        )
        assert "DIAGNOSTIC" in msg and "NOT a verdict" in msg, (
            "the number must be labelled: it overlaps, so a reader who "
            "thresholds it will be wrong in both directions"
        )

    def test_it_is_its_OWN_sentence_not_the_failing_level_s(self) -> None:
        """⭐ The pairing claim — the defect belongs to the RETURNED iterate.

        ``first_failure`` is children-first, so on an eigenvalue solve every
        other number in the message belongs to outer-iteration ONE's inner,
        while the defect is measured on the LAST iterate.  Both are real;
        presenting them as one level's facts is the "every number real,
        every pairing wrong" defect N6a removed.  So the defect must appear
        in a sentence whose subject is the iterate, and must NOT sit inside
        the ``hit <knob>=<budget> …`` clause.
        """
        with pytest.warns(ConvergenceWarning) as caught:
            sol = _starved_inner_converged_outer()
        assert sol.history is not None
        msg = str(caught[0].message)
        defect_text = f"{sol.history.balance_defect:.3e}"

        head, _, tail = msg.partition("Returning a BEST-EFFORT iterate")
        assert defect_text not in head, (
            "the defect must not appear in the failing LEVEL's clause — "
            "that clause's numbers are the starved inner's, and this one "
            "is the returned iterate's"
        )
        assert "The returned iterate leaves" in tail and defect_text in tail

    def test_a_CARRYING_mesh_warns_WITHOUT_a_number(self) -> None:
        """⛔ The curvilinear exemption, pinned so it cannot drift silently.

        ``evaluate_residual`` REFUSES a bare System-A residual on a mesh
        with starting-direction levels, because it would omit ``r_B``.  The
        ``solve_sn`` exit assembles exactly that bare shape, so on a
        carrying mesh it reports no number — deliberately, and tracked as
        #354.

        ⭐ This row exists because the omission was found by a suite run,
        not by review: `[M]` 2026-08-10 the first wiring took the slice
        from 9 to **25** reds, all 16 new ones curvilinear solves raising
        out of ``evaluate_residual``.  Nothing in review caught it because
        ``_certify_within_group_exit`` calls the same function on the same
        meshes and never hits it — it is guarded on ``record.converged``
        and returns early on precisely the truncated solves this runs on.
        The complement of a guard reaches states its partner never visits.

        Both halves are asserted: the solve is still AUDIBLE (the omission
        must cost the number, never the warning), and the number is absent
        rather than wrong.
        """
        from orpheus.geometry import Mesh1D, Region, RegionMesh, StructuredGeometry

        geom = StructuredGeometry(
            geometry="SPH",
            regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
            bcs=(BC.reflective,),
        )
        mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=20),))
        with pytest.warns(ConvergenceWarning) as caught:
            sol = solve_sn(
                {0: get_mixture("A", "2g")}, mesh,
                Quadrature.gauss_legendre(n_ordinates=8),
                inner_solver="source_iteration",
                max_outer=3, max_inner=4, inner_tol=1e-12,
                keff_tol=1e-10, flux_tol=1e-9,
            )
        assert sol.mesh.radial_characteristic_field_space is not None, (
            "fixture drift: this row needs a CARRYING mesh, or it is "
            "silently re-testing the Cartesian path"
        )
        assert sol.history is not None
        assert sol.history.fully_converged is False
        assert sol.history.balance_defect is None, (
            "a bare System-A residual on a carrying mesh would omit r_B; "
            "the exit must report no number rather than a partial one (#354)"
        )
        assert "balance defect" not in str(caught[0].message), (
            "no number ⟹ no clause — an 'unavailable' string would read as "
            "a measurement to anyone skimming"
        )

    def test_a_CONVERGED_solve_neither_WARNS_nor_PAYS(self) -> None:
        """No warning, no number — and the second half is the cost claim.

        ``_exit_balance_defect`` returns ``None`` on a fully-converged tree
        BEFORE evaluating anything, so the happy path keeps exactly the
        cost it had before N6b.  That guard is the complement of
        ``_certify_within_group_exit``'s: the certificate fires when the
        solve claimed convergence and asserts, this fires when it did not
        and reports, and no solve pays for both forward applies.

        ⚠ ``balance_defect is None`` is therefore load-bearing, not
        incidental — a change that computed it unconditionally would leave
        every assertion here green while silently adding a forward apply
        to every converged solve in the suite.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("error", ConvergenceWarning)
            sol = _fixed_source(max_inner=2000)
        assert sol.history is not None
        assert sol.history.fully_converged is True
        assert sol.history.balance_defect is None, (
            "a converged exit is CERTIFIED, not diagnosed — reporting a "
            "defect here means the guard was dropped and every converged "
            "solve is now paying for a residual it does not need"
        )


@pytest.mark.foundation
class TestTheWarningBlamesTheCallerAtEverySNEntry:
    """⭐ ``stacklevel=3`` must point at the USER's line, at every SN entry.

    ⛔ **This was a LIVE defect until 2026-08-11** (#340 N4.7), and SN was the
    only family with it.  ``solve_sn_fixed_source`` dispatches to
    ``_solve_fixed_source_si`` / ``_solve_fixed_source_krylov``, and the
    emission sat inside those private arms — so frame 3 was the entry's own
    ``return _solve_fixed_source_si(`` line, and `[M]` the warning blamed
    ``orpheus/sn/solver.py:3541``.  A file the reader did not write.

    It survived because it is invisible to every other kind of gate: the
    warning still appeared, still named the right level, still carried the
    right budget and projection.  Only the *attribution* was wrong, and
    `[M]` ``grep -rn stacklevel tests/`` had exactly ONE hit before this
    class, unrelated.

    The fix was to hoist the call into the public entry — reading the record
    off the ``Solution`` about to be returned, which also collapsed two mirror
    emission points into one — and NOT to pass a per-call ``stacklevel``: a
    frame COUNT is a fact about the call chain asserted at the call site, and
    it rots the moment a helper is interposed.

    ⚠ The companion source gate
    (``tests/numerics/test_family_convergence_contract.py::
    TestEveryEntryCallsTheHelperCorrectly``) refuses a call from a private
    function outright.  `[M]` run against ``git show HEAD:`` before the fix, it
    fails on both arms — so the pair is: this class sees the SYMPTOM, that one
    forbids the CAUSE.
    """

    @pytest.mark.parametrize(
        "row_id,run",
        [
            ("solve_sn", lambda: _eigenvalue(max_outer=2, keff_tol=1e-14)),
            ("solve_sn_fixed_source[si]", lambda: _fixed_source(max_inner=2)),
            ("solve_sn_adjoint", _starved_adjoint_eigenvalue),
            ("solve_sn_adjoint_fixed_source", _starved_adjoint_fixed_source),
        ],
    )
    def test_the_warning_is_attributed_outside_orpheus(self, row_id, run):
        from pathlib import Path

        import orpheus

        root = Path(orpheus.__file__).parent.resolve()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            run()
        conv = [w for w in caught
                if issubclass(w.category, ConvergenceWarning)]
        assert conv, f"{row_id}: the fixture must starve a level"
        for w in conv:
            assert not Path(w.filename).resolve().is_relative_to(root), (
                f"{row_id}: blamed {w.filename}:{w.lineno}, inside orpheus/ — "
                "the emission is being made from a private helper again "
                "(#340 N4.7)"
            )

    def test_the_warning_names_the_EXACT_call_line(self) -> None:
        """⭐⭐ The sharp leg, and why containment alone is not enough.

        `[M]` a ``stacklevel`` that is too LARGE also escapes ``orpheus/`` — it
        lands on the caller's *caller*.  So the row above is structurally
        incapable of seeing an over-deep frame count; it is Mode-12 blind to
        half the failure class.  Measured by mutation: ``stacklevel=4`` leaves
        the containment rows GREEN and reds only this one.

        ⚠ The call below MUST be literal, one line under the ``f_lineno``
        read.  Routing it through ``_fixed_source(...)`` would pin the
        HELPER's line — a true statement about nothing.
        """
        import inspect

        quad = Quadrature.level_symmetric(sn_order=4)
        frame = inspect.currentframe()
        assert frame is not None, "CPython always provides one here"
        expected = frame.f_lineno + 3
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            solve_sn_fixed_source(
                {0: _absorber_2g()}, _d3_reflective_axes(), quad,
                external_source=_uniform_source(quad, (3, 4, 5)),
                boundary_condition="reflective",
                inner_tol=1e-13, max_inner=2,
            )
        conv = [w for w in caught
                if issubclass(w.category, ConvergenceWarning)]
        assert len(conv) == 1
        assert conv[0].filename == __file__
        assert conv[0].lineno == expected, (
            f"attributed to line {conv[0].lineno}, the call opens at "
            f"{expected} — stacklevel is off by "
            f"{conv[0].lineno - expected} lines' worth of frames"
        )
