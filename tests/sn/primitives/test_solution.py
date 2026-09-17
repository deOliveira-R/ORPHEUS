"""Foundation tests for the Issue #197 PR-TYPED-5 :class:`Solution` type.

Pins the structural contract of :class:`Solution` +
:class:`IterationHistory` + :class:`SolutionDiff` — since step 3 of the
consumers campaign (2026-09-17) the carrier is (Problem, posing) + Strategy +
records, built here the way the mints build it (a real hub, a real posing):

* the state-on-domain law at construction (the outcome's state lives on the
  Problem's coupled space; a cross-hub pairing is refused);
* the flux members READ off the state (``angular_flux``, ``boundary_flux``,
  ``radial_characteristic`` by arity, ``scalar_flux`` derived and cached);
* the KIND is the outcome's TYPE (``EigenOutcome`` / ``SourceOutcome``; a
  source Solution has no ``keff`` attribute at all);
* :meth:`IterationHistory.dominance_ratio` for a recorded trajectory (the
  view type), and ``outcome.dominance_ratio()`` on the eigen kind;
* :meth:`SolutionBase.compare` field-by-field summary;
* :meth:`Solution.reaction_rate_density` :math:`\\sigma \\cdot \\phi` math;
* the ROLE axis (#276 A5): ``SolutionBase`` → {``Solution``,
  ``AdjointSolution``} — base non-instantiable, forward-physics trio
  structurally absent on the adjoint, ``importance`` alias, role-closed
  ``compare``.

These are foundation tests — software invariants of the typed return
contract, not physics claims about a solver, so they carry
``@pytest.mark.foundation`` per the V&V harness convention.
"""
from __future__ import annotations

from typing import Any, cast

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.numerics.convergence import (
    IterationBudget,
    IterationRecord,
    StoppingCriterion,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.timed_full_field import TimedFullField
from orpheus.numerics.coupled_system import CoupledField
from orpheus.numerics.gauge import ScaleGauge
from orpheus.numerics.outcome import (
    Certified,
    EigenOutcome,
    ExitCertificate,
    Measured,
    NotApplicable,
    SourceOutcome,
)
from orpheus.numerics.posing import SourcePosing
from orpheus.sn.splitting import Splitting, resolve_schedule
from orpheus.sn.solution import (
    AdjointSolution,
    IterationHistory,
    Solution,
    SolutionBase,
    SolutionDiff,
)

from tests.sn._test_helpers import placeholder_materials
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux

pytestmark = pytest.mark.foundation


def _history(
    *,
    converged: bool = True,
    keff_history: tuple[float, ...] = (),
    flux_residuals: tuple[float, ...] = (),
) -> IterationHistory:
    """Build a history that READS as requested, by stating a trajectory.

    ⛔ ``converged`` stopped being a field on 2026-08-09 (#340 N2b-ii):
    :class:`IterationHistory` is now a view over an
    :class:`~orpheus.numerics.convergence.IterationRecord` and every verdict
    is derived from what was measured.  So a test cannot assert the verdict
    into existence — it has to do what a producer does and supply a
    trajectory that clears (or does not) against a tolerance.

    That is deliberately a little more work here, and it is the whole point:
    the old spelling ``IterationHistory(converged=True)`` is exactly the
    hand-written claim #342 was, and making it unspellable in TESTS too is
    what stops a fixture from certifying a state the production types can no
    longer produce.

    Builds a LEAF record (a within-group solve), so ``n_inner`` is populated
    and ``n_outer`` is ``None``.
    """
    trajectory = tuple(flux_residuals) or ((0.0,) if converged else (1.0,))
    # A tolerance above the whole trajectory clears; ``0.0`` is legal and,
    # since ``cleared`` is a strict ``<``, never clears.
    tolerance = (max(trajectory) + 1.0) if converged else 0.0
    return IterationHistory(
        record=IterationRecord(
            label="inner(probe)",
            criteria=(
                StoppingCriterion(
                    name="residual", trajectory=trajectory,
                    tolerance=tolerance,
                ),
            ),
            iterations_run=len(trajectory),
        ),
        keff_history=keff_history,
    )


def _slab_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
    """Build a small slab :class:`SNMesh` for unit testing."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _quad8_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
    """GL(8) sibling — the ANGULAR-only discriminator: the (energy,
    spatial) marginal is identical, so only the ψ gate can refuse it."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    return SNMesh(
        mesh, Quadrature.gauss_legendre(n_ordinates=8),
        placeholder_materials(ng=ng),
    )


# ════════════════════════════════════════════════════════════════════
# IterationHistory
# ════════════════════════════════════════════════════════════════════


class TestIterationHistory:
    def test_a_record_that_measured_nothing_reads_empty(self) -> None:
        """An empty record reads empty everywhere, and claims convergence.

        ⛔ Was ``test_default_empty``, asserting the DEFAULTS of six fields.
        There are no defaults now (#340 N2b-ii): the record is required and
        every reading below is derived from it, so what this pins is the
        derivation's behaviour on the degenerate tree rather than a
        constructor's argument list.

        ``converged is True`` on a level that never entered its loop is
        deliberate and is the reading the retired ``_claims_convergence``
        got backwards: a direct solve — a factorisation, a Krylov call
        satisfied by its initial guess — did not FAIL to converge. The
        neighbouring state (a level that RAN and measured nothing) is the one
        that must not claim, and it is discriminated by ``iterated``.
        """
        h = IterationHistory(record=IterationRecord(label="inner(probe)"))
        assert h.keff_history == ()
        assert h.flux_residuals == ()
        assert h.n_inner == 0
        assert h.n_outer is None          # a leaf: no outer to report
        assert h.total_inner_iterations == 0
        assert h.converged is True
        assert h.fully_converged is True
        assert h.record.iterated is False
        assert h.dominance_ratio() is None
        assert h.latest_keff() is None
        assert h.latest_residual() is None

    def test_a_level_that_RAN_and_measured_nothing_does_not_claim(
        self,
    ) -> None:
        """The neighbour of the case above, and it must read the other way.

        Same empty criteria, but the level iterated — an SI capped at
        ``max_iter=1`` measures no difference at all. Vacuous truth is right
        for one of these two and wrong for the other, which is why the
        discriminator is named rather than a side being picked.
        """
        h = IterationHistory(
            record=IterationRecord(label="inner(probe)", iterations_run=1),
        )
        assert h.record.iterated is True
        assert h.converged is False

    def test_converged_answers_for_the_LEVEL_and_fully_converged_for_the_TREE(
        self,
    ) -> None:
        """Two questions, two answers — on the fixture that separates them.

        ⭐ Both are one-line delegations to the record, which is exactly why
        they need a gate here rather than being taken on trust: `[M]` the same
        conflation on ``PowerIterationOutcome`` reddened **0 of 249** tests,
        because every other gate reached through to the record directly and
        the fixtures in use had the level and the fold agreeing.  A delegating
        property is only gated by a fixture where the delegate's two readings
        DIFFER.

        The state below is the #340 headline: an outer whose own criteria
        cleared, standing on an inner that was truncated.
        """
        starved_inner = IterationRecord(
            label="inner(probe)",
            criteria=(
                StoppingCriterion(name="residual", trajectory=(1.0,),
                                  tolerance=1e-8),
            ),
            budget=IterationBudget(1), iterations_run=1,
        )
        h = IterationHistory(
            record=IterationRecord(
                label="outer(probe)",
                criteria=(
                    StoppingCriterion(name="dk", trajectory=(0.0,) * 3,
                                      tolerance=1.0),
                ),
                iterations_run=3, children=(starved_inner,) * 3,
            ),
        )

        assert h.converged is True, "the outer's OWN criteria did clear"
        assert h.fully_converged is False, (
            "the AND-fold over the tree must refuse what the level alone "
            "would claim — a value gate wants THIS one"
        )
        assert h.record.first_failure is not None
        assert h.record.first_failure.label == "inner(probe)"

    def test_the_inner_TOTAL_sums_the_children_not_the_outer_itself(
        self,
    ) -> None:
        """``total_inner_iterations`` is a sum over children, and it must be
        DISTINGUISHABLE from the outer's own count.

        ⭐ This gate exists because a mutation found its absence, and the
        near-miss is the instructive part. `[M]` replacing the sum with
        ``record.n_iterations`` reddened **0** — including
        ``test_si_convergence_rate``'s ``total_inner >= n_outer > 0``, which
        LOOKS like the guard for exactly this and cannot be: the mutation
        makes the two sides equal, and ``n >= n`` is true. **A ``>=`` between
        two quantities is blind to any mutation that collapses one onto the
        other.**

        So the fixture is built with children that ran MORE than once, which
        is the only shape where sum-over-children and count-of-outers differ.
        """
        inner = IterationRecord(
            label="inner(probe)",
            criteria=(
                StoppingCriterion(name="residual", trajectory=(0.0,) * 7,
                                  tolerance=1.0),
            ),
            iterations_run=7,
        )
        h = IterationHistory(
            record=IterationRecord(
                label="outer(probe)",
                criteria=(
                    StoppingCriterion(name="dk", trajectory=(0.0,) * 3,
                                      tolerance=1.0),
                ),
                iterations_run=3, children=(inner,) * 3,
            ),
        )

        assert h.n_outer == 3
        assert h.total_inner_iterations == 21, "3 outers x 7 inner passes"
        assert h.total_inner_iterations != h.n_outer, (
            "the fixture must SEPARATE the two, or this gate is a tautology"
        )
        assert h.n_inner is None, "no single inner solve to point at"

    def test_flux_residuals_stays_EMPTY_on_the_eigenvalue_path(self) -> None:
        """Deliberately not widened to ``dphi`` — and the choice is gated.

        ⭐ Also found blind by mutation. The outer DOES now record ``dphi``,
        which is the same quantity ``flux_residuals`` names, so returning it
        here is a tempting one-line "improvement". It would silently re-point
        every consumer that branches on ``if history.flux_residuals:`` — the
        DSA rate probes among them — from "this is a within-group solve" to
        "this is any solve at all".

        The eigenvalue path's criteria are not lost; they are one attribute
        away on :attr:`IterationHistory.record`. A deliberate design choice
        that nothing gates is an ungated coverage claim, so this is the gate.
        """
        h = IterationHistory(
            record=IterationRecord(
                label="outer(probe)",
                criteria=(
                    StoppingCriterion(name="dk", trajectory=(1e-3, 1e-6),
                                      tolerance=1.0),
                    StoppingCriterion(name="dphi", trajectory=(1e-2, 1e-5),
                                      tolerance=1.0),
                ),
                iterations_run=2,
                children=(IterationRecord(label="inner(probe)"),) * 2,
            ),
        )

        assert h.flux_residuals == (), (
            "the flat name means the WITHIN-GROUP residual; the outer's "
            "criteria reach the caller through .record"
        )
        assert h.latest_residual() is None
        # ...and the data is genuinely still available, one hop away.
        assert [c.name for c in h.record.criteria] == ["dk", "dphi"]

    def test_dominance_ratio_single_entry(self) -> None:
        h = _history(converged=True, keff_history=(1.234,))
        assert h.dominance_ratio() is None  # need ≥ 2

    def test_dominance_ratio_three_iters(self) -> None:
        h = _history(converged=True, keff_history=(1.0, 1.1, 1.10005))
        # |1.10005 - 1.1| / |1.1| ≈ 4.545e-5
        ratio = h.dominance_ratio()
        assert ratio is not None
        np.testing.assert_allclose(ratio, 5e-5 / 1.1, rtol=1e-10)

    def test_dominance_ratio_zero_prev(self) -> None:
        """A zero penultimate keff returns None to avoid divide-by-zero."""
        h = _history(converged=True, keff_history=(0.0, 1.0))
        assert h.dominance_ratio() is None

    def test_latest_keff(self) -> None:
        h = _history(converged=True, keff_history=(1.0, 1.1, 1.05))
        assert h.latest_keff() == 1.05

    def test_latest_residual(self) -> None:
        h = _history(converged=True, flux_residuals=(1e-3, 1e-6, 1e-9))
        assert h.latest_residual() == 1e-9

    def test_frozen(self) -> None:
        """IterationHistory is frozen — fields cannot be reassigned."""
        h = _history(converged=True, keff_history=(1.0,))
        with pytest.raises((AttributeError, Exception)):
            h.keff_history = (1.0, 1.1)  # type: ignore[misc]

    def test_keff_history_is_tuple(self) -> None:
        """The trajectory MUST be a tuple — frozen-with-mutable-list anti-pattern."""
        h = _history(converged=True, keff_history=(1.0, 1.1))
        assert isinstance(h.keff_history, tuple)
        # Pinning the contract: a list passed in survives but the field
        # is documented as a tuple; consumers MUST treat it as immutable.


# ════════════════════════════════════════════════════════════════════
# Solution — construction on the reshaped carrier (step 3, 2026-09-17)
# ════════════════════════════════════════════════════════════════════
#
# A Solution is (Problem, posing) + Strategy + records: ``mesh``, a kind-typed
# ``outcome`` (question + returned STATE + answer + gauge), ``strategy``,
# ``certificate``, ``record``.  The flux members are READ off the state.  The
# rows below build Solutions the way the mints do — on a real hub with a real
# posing — never a bag of hand-made fields.

_W4 = float(np.sum(Quadrature.gauss_legendre(n_ordinates=4).weights))  # ∫ 1 dΩ on the GL-4 rule


def _state(sn_mesh: SNMesh, fill: float = 1.0) -> CoupledField:
    """The returned state WHOLE: the one-system coupled field over the hub's full field."""
    member = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn_mesh.full_field_space)
    member.interior.values[:] = fill
    member.boundary.values[:] = fill
    return CoupledField(systems=(member,))


def _member(state: CoupledField) -> TimedFullField:
    """System A's composite (the member contract is structural; the fixture's is a TimedFullField)."""
    return cast(TimedFullField, state.systems[0])


def _state_with_scalar(sn_mesh: SNMesh, phi_values: np.ndarray) -> CoupledField:
    """A state whose derived scalar flux ∫ψ dΩ reproduces ``phi_values`` (ψ uniform in angle)."""
    state = _state(sn_mesh, 0.0)
    _member(state).interior.values[:] = (np.asarray(phi_values, dtype=float) / _W4)[None]
    return state


def _record(converged: bool = True) -> IterationRecord:
    return _history(converged=converged).record


def _certificate() -> ExitCertificate:
    return ExitCertificate(
        balance=NotApplicable("a fixture"), gauge=NotApplicable("a fixture"),
        rayleigh_gap=NotApplicable("a fixture"), admissibility=NotApplicable("a fixture"),
    )


def _strategy(sn_mesh: SNMesh) -> Splitting:
    return Splitting.from_schedule(sn_mesh.system, resolve_schedule(sn_mesh, "jacobi"))


def _eigen(sn_mesh: SNMesh, state: CoupledField | None = None, *, lam: float = 1.0, trajectory: tuple[float, ...] = ()) -> EigenOutcome:
    return EigenOutcome(
        posing=sn_mesh.eigen_posing, state=_state(sn_mesh) if state is None else state,
        lam=lam, trajectory=trajectory or (lam,), gauge=ScaleGauge(lambda s: 1.0, 1.0),
    )


def _source(sn_mesh: SNMesh, state: CoupledField | None = None) -> SourceOutcome:
    return SourceOutcome(
        posing=SourcePosing(sn_mesh.pencil.at(0.0), sn_mesh.system.space.zeros()),
        state=_state(sn_mesh) if state is None else state, gauge=sn_mesh.loss_kernel_gauge,
    )


def _solution(sn_mesh: SNMesh, outcome, *, cls: type[Any] = Solution, record: IterationRecord | None = None) -> Any:
    return cls(
        mesh=sn_mesh, outcome=outcome, strategy=_strategy(sn_mesh),
        certificate=_certificate(), record=_record() if record is None else record,
    )


class TestSolutionConstruction:
    def test_construct_fixed_source_reads_its_members_off_the_state(self) -> None:
        m = _slab_mesh()
        state = _state(m)
        sol = _solution(m, _source(m, state))
        assert sol.state is state
        assert sol.angular_flux is _member(state)
        assert sol.boundary_flux is _member(state).boundary
        assert sol.radial_characteristic is None
        assert sol.mesh is m
        assert not hasattr(sol, "keff"), "a source Solution has NO keff — the kind is the outcome's type"
        assert isinstance(sol.outcome, SourceOutcome)

    def test_construct_eigenvalue(self) -> None:
        m = _slab_mesh()
        sol = _solution(m, _eigen(m, lam=1.05, trajectory=(1.0, 1.05)))
        assert sol.outcome.keff == 1.05
        assert sol.outcome.trajectory == (1.0, 1.05)
        assert isinstance(sol.outcome, EigenOutcome)

    def test_scalar_flux_is_DERIVED_from_the_state(self) -> None:
        """φ = ∫ψ dΩ of the state's cell-average moment — one quantity, one representation."""
        m = _slab_mesh()
        sol = _solution(m, _source(m, _state(m, fill=1.5)))
        np.testing.assert_allclose(sol.scalar_flux.values, 1.5 * _W4, rtol=1e-14)
        assert sol.scalar_flux is sol.scalar_flux, "a cached reading — identity holds across reads"
        assert sol.scalar_flux.space == m.bulk_space

    def test_the_state_on_domain_law_refuses_a_cross_hub_pairing(self) -> None:
        """The outcome's state must live on THIS Problem's coupled space (content identity)."""
        m1, m8 = _slab_mesh(), _quad8_mesh()
        with pytest.raises(ValueError, match="state-on-domain"):
            _solution(m8, _source(m1))          # m1's state and question on m8
        with pytest.raises(ValueError, match="state-on-domain"):
            _solution(m8, _eigen(m1))

    def test_the_same_hub_cross_kind_replace_is_a_LEGAL_different_solve(self) -> None:
        """Both kinds' questions share the coupled space, so this is a different
        solve, not a refusal (the test-architect's declared null for arm A5)."""
        from dataclasses import replace
        m = _slab_mesh()
        sol = _solution(m, _eigen(m))
        other = replace(sol, outcome=_source(m))
        assert isinstance(other.outcome, SourceOutcome)

    def test_boundary_flux_delegates_to_the_states_member(self) -> None:
        m = _slab_mesh()
        state = _state(m)
        sol = _solution(m, _source(m, state))
        assert sol.boundary_flux is _member(state).boundary

    def test_frozen(self) -> None:
        m = _slab_mesh()
        sol = _solution(m, _eigen(m))
        with pytest.raises((AttributeError, Exception)):
            sol.outcome = _source(m)  # type: ignore[misc]


# ════════════════════════════════════════════════════════════════════
# Solution — the kind is the outcome's TYPE
# ════════════════════════════════════════════════════════════════════


class TestSolutionKind:
    def test_the_kind_is_the_outcomes_type(self) -> None:
        m = _slab_mesh()
        assert isinstance(_solution(m, _source(m)).outcome, SourceOutcome)
        assert isinstance(_solution(m, _eigen(m)).outcome, EigenOutcome)

    def test_no_field_of_the_carrier_is_optional_by_kind(self) -> None:
        import dataclasses, typing
        for f in dataclasses.fields(SolutionBase):
            assert "None" not in str(f.type) and typing.get_origin(f.type) is not typing.Union, f.name
        assert [f.name for f in dataclasses.fields(SolutionBase)] == ["mesh", "outcome", "strategy", "certificate", "record"]

    def test_an_eigen_outcome_refuses_a_source_question_at_the_sn_tier(self) -> None:
        m = _slab_mesh()
        with pytest.raises(TypeError, match="must be an EigenPosing"):
            EigenOutcome(posing=SourcePosing(m.pencil.lhs, m.system.space.zeros()), state=_state(m), lam=1.0, trajectory=(1.0,), gauge=ScaleGauge(lambda s: 1.0, 1.0))  # type: ignore[arg-type]


# ════════════════════════════════════════════════════════════════════
# Solution — diagnostics read the record and the outcome
# ════════════════════════════════════════════════════════════════════


class TestSolutionDiagnostics:
    def test_dominance_ratio_three_iterations(self) -> None:
        m = _slab_mesh()
        sol = _solution(m, _eigen(m, lam=1.10005, trajectory=(1.0, 1.1, 1.10005)))
        ratio = sol.outcome.dominance_ratio()
        assert ratio is not None
        np.testing.assert_allclose(ratio, 5e-5 / 1.1, rtol=1e-10)

    def test_dominance_ratio_of_a_one_point_trajectory_is_undefined(self) -> None:
        m = _slab_mesh()
        assert _solution(m, _eigen(m)).outcome.dominance_ratio() is None

    def test_converged_reads_the_record(self) -> None:
        m = _slab_mesh()
        assert _solution(m, _source(m), record=_record(converged=True)).converged() is True
        assert _solution(m, _source(m), record=_record(converged=False)).converged() is False

    def test_converged_asks_the_TREE_not_the_top_level(self) -> None:
        r"""An outer that met its own criteria while standing on a STARVED inner
        is not a trustworthy answer (#340 F1): ``Solution.converged()`` reads
        ``record.fully_converged`` — the tree-wide question."""
        outer = IterationRecord(
            label="outer(power-iteration)",
            budget=IterationBudget(50, "max_outer"),
            iterations_run=3,
            criteria=(StoppingCriterion(
                name="dk", tolerance=1e-6, trajectory=(1e-3, 1e-5, 1e-8)),),
            children=(IterationRecord(
                label="inner(source-iteration)",
                budget=IterationBudget(8, "max_inner"), iterations_run=8,
                criteria=(StoppingCriterion(
                    name="residual", tolerance=1e-10,
                    trajectory=(1e-2, 1e-3, 1e-4)),),
            ),),
        )
        assert outer.converged is True, "fixture drift: the outer must clear"
        assert outer.fully_converged is False, (
            "fixture drift: the inner must be starved, or this row "
            "degenerates into the leaf pair above"
        )
        m = _slab_mesh()
        sol = _solution(m, _eigen(m, trajectory=(1.0, 1.0, 1.0)), record=outer)
        assert sol.converged() is False, (
            "Solution.converged() must answer the TREE — an outer standing "
            "on a starved inner is not a trustworthy answer (#340 F1)"
        )

    def test_the_trajectory_replaces_keff_history_list(self) -> None:
        m = _slab_mesh()
        sol = _solution(m, _eigen(m, lam=1.04, trajectory=(1.0, 1.05, 1.04)))
        assert list(sol.outcome.trajectory) == [1.0, 1.05, 1.04]
        assert not hasattr(sol, "keff_history_list")

    def test_the_history_view_is_a_one_cycle_reading_of_the_new_members(self) -> None:
        """U2 → U6: the pre-step-3 ``IterationHistory`` reading, assembled from
        the record, the trajectory and the certificate's magnitudes — the
        ``float | None`` leak confined to the view that retires at U6."""
        m = _slab_mesh()
        sol = _solution(m, _eigen(m, lam=1.04, trajectory=(1.0, 1.04)))
        h = sol.history
        assert isinstance(h, IterationHistory)
        assert h.record is sol.record
        assert h.keff_history == (1.0, 1.04)
        assert h.balance_defect is None and h.gauge_correction is None, "NotApplicable reads None on the view"
        measured = ExitCertificate(balance=Measured(0.3), gauge=Measured(0.01), rayleigh_gap=NotApplicable("x"), admissibility=Certified(0.9, "x"))
        from dataclasses import replace
        h2 = replace(sol, certificate=measured).history
        assert h2.balance_defect == 0.3 and h2.gauge_correction == 0.01
        assert _solution(m, _source(m)).history.keff_history == ()


# ════════════════════════════════════════════════════════════════════
# Solution — reaction-rate accessor
# ════════════════════════════════════════════════════════════════════


class TestReactionRate:
    def test_reaction_rate_density_shape(self) -> None:
        m = _slab_mesh(nx=4, ng=2)
        phi_values = np.arange(m.ng * m.nx, dtype=float).reshape(m.ng, *m.spatial_shape) + 1.0
        sol = _solution(m, _source(m, _state_with_scalar(m, phi_values)))
        xs = np.full((m.ng, *m.spatial_shape), 0.5)
        rate = sol.reaction_rate_density(xs)
        assert rate.shape == (m.ng, *m.spatial_shape)

    def test_reaction_rate_density_math(self) -> None:
        r"""σ · φ at each cell — the named math reads as the formula."""
        m = _slab_mesh(nx=3, ng=2)
        phi_values = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])  # (ng=2, nx=3)
        sol = _solution(m, _source(m, _state_with_scalar(m, phi_values)))
        np.testing.assert_allclose(sol.scalar_flux.values, phi_values, rtol=1e-14)
        xs = np.full((m.ng, *m.spatial_shape), 0.0)
        xs[0, :] = 0.5
        xs[1, :] = 0.25
        rate = sol.reaction_rate_density(xs)
        np.testing.assert_allclose(rate[0, :], [0.5, 1.0, 1.5], rtol=1e-13)
        np.testing.assert_allclose(rate[1, :], [1.0, 1.25, 1.5], rtol=1e-13)

    def test_reaction_rate_density_zero_flux(self) -> None:
        m = _slab_mesh()
        sol = _solution(m, _source(m, _state(m, fill=0.0)))
        xs = np.full((m.ng, *m.spatial_shape), 0.5)
        rate = sol.reaction_rate_density(xs)
        np.testing.assert_array_equal(rate, np.zeros_like(rate))


# ════════════════════════════════════════════════════════════════════
# Solution — compare
# ════════════════════════════════════════════════════════════════════


class TestSolutionCompare:
    def test_compare_identical_within_tolerance(self) -> None:
        m = _slab_mesh()
        state = _state(m)
        sol_a = _solution(m, _eigen(m, state, lam=1.0))
        sol_b = _solution(m, _eigen(m, state, lam=1.0))
        diff = sol_a.compare(sol_b, rtol=1e-12)
        assert isinstance(diff, SolutionDiff)
        assert diff.keff_abs == Measured(0.0)
        assert diff.angular_flux_linf == 0.0
        assert diff.scalar_flux_linf == 0.0
        assert diff.within_tolerance is True

    def test_compare_different_keff(self) -> None:
        m = _slab_mesh()
        state = _state(m)
        sol_a = _solution(m, _eigen(m, state, lam=1.0))
        sol_b = _solution(m, _eigen(m, state, lam=1.001))
        diff = sol_a.compare(sol_b, rtol=1e-6)
        assert isinstance(diff.keff_abs, Measured)
        np.testing.assert_allclose(diff.keff_abs.value, 0.001, rtol=1e-12)
        assert diff.within_tolerance is False  # 0.001 > 1e-6 * 1.001

    def test_compare_different_flux(self) -> None:
        m = _slab_mesh()
        sol_a = _solution(m, _source(m, _state(m, fill=1.0)))
        sol_b = _solution(m, _source(m, _state(m, fill=1.1)))
        diff = sol_a.compare(sol_b, rtol=1e-12)
        np.testing.assert_allclose(diff.scalar_flux_linf, 0.1 * _W4, rtol=1e-12)
        np.testing.assert_allclose(diff.angular_flux_linf, 0.1, rtol=1e-12)
        assert diff.keff_abs == NotApplicable("the source kind carries no eigenvalue")
        assert diff.within_tolerance is False

    def test_compare_is_KIND_closed(self) -> None:
        """An eigen answer against a source answer is a different QUESTION —
        refused, where the pre-step-3 carrier silently skipped the eigenvalue
        channel (``keff_abs is None``)."""
        m = _slab_mesh()
        with pytest.raises(TypeError, match="kind mismatch"):
            _solution(m, _eigen(m)).compare(_solution(m, _source(m)))

    def test_compare_cross_mesh_rejected(self) -> None:
        """Two hubs on different quadratures realize different phase spaces —
        the layout gate (``same_phase_space``) rejects."""
        m1, m8 = _slab_mesh(), _quad8_mesh()
        sol_a = _solution(m1, _source(m1))
        sol_b = _solution(m8, _source(m8))
        with pytest.raises(ValueError, match="different discrete phase space"):
            sol_a.compare(sol_b)

    def test_compare_accepts_sibling_snmesh_over_same_constituents(self) -> None:
        """Two SNMesh wrappers over the same constituents realize ONE phase
        space (content identity) — ``compare`` must accept, and so must the
        state-on-domain law (the state's space equals the sibling's by content)."""
        geometry = Mesh1D(
            edges=np.linspace(0.0, 1.0, 5), mat_ids=np.zeros(4, dtype=int),
            coord=CoordSystem.CARTESIAN, bc_left=BC("vacuum"), bc_right=BC("vacuum"),
        )
        quad = Quadrature.gauss_legendre(n_ordinates=4)
        materials = placeholder_materials(ng=2)
        m1 = SNMesh(geometry, quad, materials)
        m2 = SNMesh(geometry, quad, materials)
        assert m1 is not m2 and m1.same_phase_space(m2)
        sol_a = _solution(m1, _source(m1))
        sol_b = _solution(m2, _source(m2))
        assert sol_a.compare(sol_b) is not None
        # and a state minted on the sibling is an element of THIS hub's space too
        _solution(m1, _source(m1, _state(m2)))


# ════════════════════════════════════════════════════════════════════
# The ROLE axis (#276 A5) — SolutionBase → {Solution, AdjointSolution}
# ════════════════════════════════════════════════════════════════════


class TestSolutionRoleAxis:
    """The A5 carrier ruling's structural contract: the role (forward vs
    adjoint) is a TYPE (the verb set differs), the kind a type PARAMETER (it
    does not) — two leaves, not four."""

    def test_base_not_instantiable(self) -> None:
        m = _slab_mesh()
        with pytest.raises(TypeError, match="not instantiable"):
            SolutionBase(mesh=m, outcome=_source(m), strategy=_strategy(m), certificate=_certificate(), record=_record())

    def test_adjoint_construction_shares_the_carrier(self) -> None:
        m = _slab_mesh()
        state = _state(m)
        adj = _solution(m, _source(m, state), cls=AdjointSolution)
        assert adj.angular_flux is _member(state)
        assert adj.boundary_flux is _member(state).boundary
        assert adj.mesh is m
        assert isinstance(adj.outcome, SourceOutcome) and not hasattr(adj, "keff")

    def test_adjoint_state_on_domain_enforced(self) -> None:
        m1, m8 = _slab_mesh(), _quad8_mesh()
        with pytest.raises(ValueError, match="state-on-domain"):
            _solution(m8, _source(m1), cls=AdjointSolution)

    def test_roles_are_siblings_not_subtypes(self) -> None:
        assert issubclass(Solution, SolutionBase)
        assert issubclass(AdjointSolution, SolutionBase)
        assert not issubclass(AdjointSolution, Solution)
        assert not issubclass(Solution, AdjointSolution)

    def test_forward_physics_structurally_absent_on_adjoint(self) -> None:
        forward_trio = ("homogenize", "condense", "reaction_rate_density")
        for name in forward_trio:
            assert hasattr(Solution, name), f"Solution must carry {name}"
            assert not hasattr(AdjointSolution, name), (
                f"AdjointSolution must NOT carry {name} — the forward-"
                "physics asymmetry is structural (#276 A5)"
            )
        assert hasattr(AdjointSolution, "importance")
        assert not hasattr(Solution, "importance")

    def test_importance_aliases_scalar_flux(self) -> None:
        m = _slab_mesh()
        adj = _solution(m, _source(m), cls=AdjointSolution)
        assert adj.importance is adj.scalar_flux

    def test_compare_role_mismatch_rejected(self) -> None:
        m = _slab_mesh()
        fwd = cast(SolutionBase, _solution(m, _source(m)))
        adj = cast(SolutionBase, _solution(m, _source(m), cls=AdjointSolution))
        with pytest.raises(TypeError, match="role mismatch"):
            fwd.compare(adj)
        with pytest.raises(TypeError, match="role mismatch"):
            adj.compare(fwd)

    def test_compare_same_role_adjoint_works(self) -> None:
        m = _slab_mesh()
        adj_a = _solution(m, _source(m), cls=AdjointSolution)
        adj_b = _solution(m, _source(m), cls=AdjointSolution)
        assert isinstance(adj_a.compare(adj_b), SolutionDiff)
