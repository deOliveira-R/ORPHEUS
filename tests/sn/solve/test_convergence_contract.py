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
    IterationRecord,
    StoppingCriterion,
    default_iteration_budget,
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
            PowerIterationOutcome(              # type: ignore[call-arg]
                keff=1.0, keff_history=[1.0], flux_distribution=np.ones(3),
                record=IterationRecord(label="outer(probe)"),
                converged=True,
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
        from orpheus.sn.solution import IterationHistory
        from orpheus.sn.solver import _warn_if_unconverged

        stalled = IterationHistory(
            converged=False,
            flux_residuals=tuple(1e-3 for _ in range(40)),  # rho == 1
        )
        with pytest.warns(ConvergenceWarning) as record:
            _warn_if_unconverged(
                stalled, where="probe", budget_name="max_inner",
                budget=40, tol=1e-10,
            )
        msg = str(record[0].message)
        assert "NOT" in msg and "contracting" in msg
        assert "set max_inner=" not in msg  # the advice must NOT be given

    def test_it_REFUSES_to_project_off_a_criterion_that_already_cleared(
        self,
    ) -> None:
        """⛔ The regression for a lie this gate's own author shipped.

        `[M]` 2026-08-09, caught by the wide run, not by review. The outer
        stops on ``dk AND dphi``, but only ``keff_history`` survives into
        ``IterationHistory`` (#340 F2 — ``dphi`` is computed inside
        ``SNSolver.converged`` and discarded). On the mutated heterogeneous
        slab, whose NEGATIVE dominant eigenvalue makes ``dphi`` sign-alternate
        forever, ``|dk|`` sits at ``3.3e-16`` against a ``1e-9`` tolerance —
        so the freshly-added projection answered **"set max_outer=1"**.

        A confidently wrong number is worse than no number, and worse than
        the vague "raise the budget" it replaced. When the recorded criterion
        has cleared and the loop still did not converge, the binding one is
        not in this history, and the message must say so.
        """
        from orpheus.sn.solution import IterationHistory
        from orpheus.sn.solver import _warn_if_unconverged

        # |dk| decays to 3.3e-16, far below the 1e-9 it is judged against;
        # the loop is nonetheless unconverged, because dphi never cleared.
        keffs = tuple(1.0 + 1e-3 * 0.5**k for k in range(40))
        with pytest.warns(ConvergenceWarning) as record:
            _warn_if_unconverged(
                IterationHistory(converged=False, keff_history=keffs),
                where="probe", budget_name="max_outer", budget=40, tol=1e-9,
            )
        msg = str(record[0].message)
        assert "ALREADY cleared" in msg
        assert "set max_outer=" not in msg

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
