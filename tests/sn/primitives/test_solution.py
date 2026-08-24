"""Foundation tests for the Issue #197 PR-TYPED-5 :class:`Solution` type.

Pins the structural contract of :class:`Solution` +
:class:`IterationHistory` + :class:`SolutionDiff`:

* shape consistency at construction (delegates to the typed fields);
* mesh-binding identity across angular_flux / scalar_flux / boundary_flux;
* eigenvalue vs fixed-source discrimination via
  :meth:`Solution.is_eigenvalue` / :meth:`Solution.is_fixed_source`;
* :meth:`IterationHistory.dominance_ratio` for a recorded trajectory;
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

from typing import cast

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
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.timed_full_field import TimedFullField
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
    and ``n_outer`` is ``None``.  A test that needs the eigenvalue shape
    builds a record with children — see :func:`_outer_history`.
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


def _outer_history(
    *, keff_history: tuple[float, ...], n_outer: int,
) -> IterationHistory:
    """An EIGENVALUE-shaped history: an outer record over ``n_outer`` inners.

    The shape is what makes ``n_outer`` non-``None`` and ``n_inner`` ``None``
    — the view reads those off the tree's STRUCTURE (a level with children is
    an outer; a leaf is the solve), never off a field a producer set.
    """
    inner = IterationRecord(
        label="inner(probe)",
        criteria=(
            StoppingCriterion(name="residual", trajectory=(0.0,),
                              tolerance=1.0),
        ),
        iterations_run=1,
    )
    return IterationHistory(
        record=IterationRecord(
            label="outer(probe)",
            criteria=(
                StoppingCriterion(
                    name="dk", trajectory=(0.0,) * n_outer, tolerance=1.0,
                ),
            ),
            iterations_run=n_outer,
            children=(inner,) * n_outer,
        ),
        keff_history=keff_history,
    )


# ── Fixtures ─────────────────────────────────────────────────────────


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


def _make_fluxes(sn_mesh: SNMesh, fill: float = 1.0):
    """Build (state, scalar, boundary) for the given mesh.

    D-H.1b — :class:`Solution.angular_flux` is now a
    :class:`~orpheus.transport.timed_full_field.TimedFullField`
    composite. The triple return-tuple is preserved for fixture-shape
    stability; the first slot is now a TimedFullField (bulk angular
    flux + boundary trace + history), the third slot is the
    composite's boundary trace.
    """
    state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn_mesh.full_field_space)
    state.interior.values[:] = fill
    state.boundary.values[:] = fill
    phi = ScalarFlux(values=np.full((sn_mesh.ng, *sn_mesh.spatial_shape), fill), space=sn_mesh.bulk_space)
    return state, phi, state.boundary


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
# Solution — construction + mesh-identity contract
# ════════════════════════════════════════════════════════════════════


class TestSolutionConstruction:
    def test_construct_fixed_source(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol = Solution(
            angular_flux=psi,
            scalar_flux=phi,
            mesh=m,
        )
        assert sol.angular_flux is psi
        assert sol.scalar_flux is phi
        assert sol.boundary_flux is bf
        assert sol.mesh is m
        assert sol.keff is None
        assert sol.history is None

    def test_construct_eigenvalue(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        h = _outer_history(keff_history=(1.0, 1.05), n_outer=2)
        sol = Solution(
            angular_flux=psi, scalar_flux=phi,
            mesh=m, keff=1.05, history=h,
        )
        assert sol.keff == 1.05
        assert sol.history is h

    def test_space_content_angular_flux(self) -> None:
        """CS4b S3 (F2): a twin-carrier composite legitimately constructs
        (content agrees); a composite whose content differs (here: group
        structure) refuses on the space-content invariant."""
        m1 = _slab_mesh()
        state_twin = TimedFullField.zeros(
            interior=AngularFlux, boundary=AngularBoundaryFlux, space=_slab_mesh().full_field_space,
        )
        _, phi, bf = _make_fluxes(m1)
        Solution(angular_flux=state_twin, scalar_flux=phi, mesh=m1)
        # The discriminator is ANGULAR-ONLY (a GL8 carrier): the marginal
        # (energy, spatial) is IDENTICAL, so the φ gate cannot catch it —
        # this row is O6's own witness, not O7's ([M] the S3 battery
        # caught the first fixture (ng=3) being caught by O7 instead).
        state_foreign = TimedFullField.zeros(
            interior=AngularFlux, boundary=AngularBoundaryFlux,
            space=_quad8_mesh().full_field_space,
        )
        with pytest.raises(ValueError, match="angular_flux.interior.space"):
            Solution(
                angular_flux=state_foreign, scalar_flux=phi, mesh=m1,
            )

    def test_space_content_scalar_flux(self) -> None:
        """CS4b S3 (F2): φ must be ψ's non-angular marginal by CONTENT —
        a φ from a group-structure-differing carrier refuses."""
        m1 = _slab_mesh()
        m3 = _slab_mesh(ng=3)
        phi_foreign = ScalarFlux(values=np.zeros((m3.ng, *m3.spatial_shape)), space=m3.bulk_space)
        psi, _, bf = _make_fluxes(m1)
        with pytest.raises(ValueError, match="MARGINAL"):
            Solution(
                angular_flux=psi, scalar_flux=phi_foreign, mesh=m1,
            )

    def test_boundary_flux_delegates_to_angular_flux(self) -> None:
        """D-H.1b — ``sol.boundary_flux is sol.angular_flux.boundary``.

        :attr:`Solution.boundary_flux` is a property that returns
        :attr:`TimedFullField.boundary` — the composite's owned
        boundary trace.
        """
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol = Solution(angular_flux=psi, scalar_flux=phi, mesh=m)
        assert sol.boundary_flux is psi.boundary
        assert sol.boundary_flux is bf  # bf = state.boundary from _make_fluxes

    def test_frozen(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol = Solution(angular_flux=psi, scalar_flux=phi, mesh=m)
        with pytest.raises((AttributeError, Exception)):
            sol.keff = 1.0  # type: ignore[misc]


# ════════════════════════════════════════════════════════════════════
# Solution — discrimination
# ════════════════════════════════════════════════════════════════════


class TestSolutionDiscrimination:
    def test_fixed_source_no_keff(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol = Solution(angular_flux=psi, scalar_flux=phi, mesh=m)
        assert sol.is_fixed_source() is True
        assert sol.is_eigenvalue() is False

    def test_eigenvalue_has_keff(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol = Solution(
            angular_flux=psi, scalar_flux=phi,
            mesh=m, keff=1.234,
        )
        assert sol.is_eigenvalue() is True
        assert sol.is_fixed_source() is False

    def test_discrimination_methods_are_mutually_exclusive(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        for keff in (None, 1.0, 0.5, 2.0):
            sol = Solution(
                angular_flux=psi, scalar_flux=phi,
                mesh=m, keff=keff,
            )
            assert sol.is_eigenvalue() != sol.is_fixed_source()


# ════════════════════════════════════════════════════════════════════
# Solution — convergence diagnostics (delegating to history)
# ════════════════════════════════════════════════════════════════════


class TestSolutionDiagnostics:
    def test_dominance_ratio_three_iterations(self) -> None:
        """Solution.dominance_ratio delegates to IterationHistory."""
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        h = _history(converged=True, keff_history=(1.0, 1.1, 1.10005))
        sol = Solution(angular_flux=psi, scalar_flux=phi,
                       mesh=m, keff=1.10005, history=h)
        ratio = sol.dominance_ratio()
        assert ratio is not None
        np.testing.assert_allclose(ratio, 5e-5 / 1.1, rtol=1e-10)

    def test_dominance_ratio_no_history(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol = Solution(angular_flux=psi, scalar_flux=phi,
                       mesh=m)
        assert sol.dominance_ratio() is None

    def test_converged_no_history_is_not_a_convergence_claim(self) -> None:
        """A solution with NO history answers ``False``, not ``True``.

        ⛔ This gate asserted ``is True`` until 2026-08-08, under the
        docstring *"a solution without history is assumed converged"*.  That
        assumption is the #342 defect in its purest form: "nobody recorded
        whether this converged" is not evidence that it did, and the
        optimistic reading turned an ABSENCE of data into a positive claim
        — the one production branch in the tree that asserted convergence
        rather than reading it.

        The honest answer for an unrecorded solve is the same one a
        truncated solve gets, so a caller gating on this method treats both
        alike.  Re-posed, not deleted: the behaviour is still pinned, with
        the opposite expectation and the reason on the record.
        """
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol = Solution(angular_flux=psi, scalar_flux=phi,
                       mesh=m)
        assert sol.converged() is False

    def test_converged_with_history(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        h_yes = _history(converged=True)
        h_no = _history(converged=False)
        sol_yes = Solution(
            angular_flux=psi, scalar_flux=phi,
            mesh=m, history=h_yes,
        )
        sol_no = Solution(
            angular_flux=psi, scalar_flux=phi,
            mesh=m, history=h_no,
        )
        assert sol_yes.converged() is True
        assert sol_no.converged() is False

    def test_converged_asks_the_TREE_not_the_top_level(self) -> None:
        """A converged outer on a starved inner must answer ``False``.

        ⛔ The pair above cannot see this, and that is the point.  Both its
        histories are LEAVES, where ``converged`` and ``fully_converged``
        coincide by construction — so they stayed green through the
        2026-08-10 flip of :meth:`Solution.converged` from the top-level
        reading to the tree-wide one, and would stay green if it were
        flipped back.  This row is the discriminator they were missing.

        The shape is #340 F1: the outer meets its own criteria on
        increments that the starved inner suppressed, so the top level
        reports success while carrying the inner's error.
        """
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
        psi, phi, bf = _make_fluxes(m)
        sol = Solution(
            angular_flux=psi, scalar_flux=phi, mesh=m,
            history=IterationHistory(record=outer),
        )
        assert sol.converged() is False, (
            "Solution.converged() must answer the TREE — an outer standing "
            "on a starved inner is not a trustworthy answer (#340 F1)"
        )

    def test_keff_history_list(self) -> None:
        """The legacy keff_history accessor returns a plain list."""
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        h = _history(converged=True, keff_history=(1.0, 1.05, 1.04))
        sol = Solution(angular_flux=psi, scalar_flux=phi,
                       mesh=m, keff=1.04, history=h)
        lst = sol.keff_history_list()
        assert isinstance(lst, list)
        assert lst == [1.0, 1.05, 1.04]

    def test_keff_history_list_no_history(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol = Solution(angular_flux=psi, scalar_flux=phi,
                       mesh=m)
        assert sol.keff_history_list() == []


# ════════════════════════════════════════════════════════════════════
# Solution — reaction-rate accessor
# ════════════════════════════════════════════════════════════════════


class TestReactionRate:
    def test_reaction_rate_density_shape(self) -> None:
        """σ · φ has the same shape as φ (per-cell rate density)."""
        m = _slab_mesh(nx=4, ng=2)
        psi, _, bf = _make_fluxes(m)
        # Non-trivial flux + non-trivial XS to exercise the einsum
        phi_values = np.arange(
            m.ng * m.nx, dtype=float,
        ).reshape(m.ng, *m.spatial_shape) + 1.0
        phi = ScalarFlux(values=phi_values, space=m.bulk_space)
        sol = Solution(angular_flux=psi, scalar_flux=phi,
                       mesh=m)

        xs = np.full((m.ng, *m.spatial_shape), 0.5)
        rate = sol.reaction_rate_density(xs)
        assert rate.shape == (m.ng, *m.spatial_shape)

    def test_reaction_rate_density_math(self) -> None:
        r"""σ · φ at each cell — the named math reads as the formula."""
        m = _slab_mesh(nx=3, ng=2)
        psi, _, bf = _make_fluxes(m)
        phi_values = np.array([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
        ])  # (ng=2, nx=3)
        phi = ScalarFlux(values=phi_values, space=m.bulk_space)
        sol = Solution(angular_flux=psi, scalar_flux=phi,
                       mesh=m)

        # σ constant per group → rate[g,x] = σ[g] · φ[g,x]
        xs = np.full((m.ng, *m.spatial_shape), 0.0)
        xs[0, :] = 0.5
        xs[1, :] = 0.25
        rate = sol.reaction_rate_density(xs)

        np.testing.assert_allclose(
            rate[0, :], [0.5, 1.0, 1.5], rtol=1e-15,
        )
        np.testing.assert_allclose(
            rate[1, :], [1.0, 1.25, 1.5], rtol=1e-15,
        )

    def test_reaction_rate_density_zero_flux(self) -> None:
        """σ · 0 = 0 everywhere."""
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m, fill=0.0)
        sol = Solution(angular_flux=psi, scalar_flux=phi,
                       mesh=m)
        xs = np.full((m.ng, *m.spatial_shape), 0.5)
        rate = sol.reaction_rate_density(xs)
        np.testing.assert_array_equal(rate, np.zeros_like(rate))


# ════════════════════════════════════════════════════════════════════
# Solution — compare
# ════════════════════════════════════════════════════════════════════


class TestSolutionCompare:
    def test_compare_identical_within_tolerance(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol_a = Solution(angular_flux=psi, scalar_flux=phi,
                         mesh=m, keff=1.0)
        sol_b = Solution(angular_flux=psi, scalar_flux=phi,
                         mesh=m, keff=1.0)
        diff = sol_a.compare(sol_b, rtol=1e-12)
        assert isinstance(diff, SolutionDiff)
        assert diff.keff_abs == 0.0
        assert diff.angular_flux_linf == 0.0
        assert diff.scalar_flux_linf == 0.0
        assert diff.within_tolerance is True

    def test_compare_different_keff(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol_a = Solution(angular_flux=psi, scalar_flux=phi,
                         mesh=m, keff=1.0)
        sol_b = Solution(angular_flux=psi, scalar_flux=phi,
                         mesh=m, keff=1.001)
        diff = sol_a.compare(sol_b, rtol=1e-6)
        assert diff.keff_abs is not None
        np.testing.assert_allclose(diff.keff_abs, 0.001, rtol=1e-12)
        assert diff.within_tolerance is False  # 0.001 > 1e-6 * 1.001

    def test_compare_different_flux(self) -> None:
        m = _slab_mesh()
        psi, phi_a, bf = _make_fluxes(m, fill=1.0)
        _, phi_b, _ = _make_fluxes(m, fill=1.1)
        sol_a = Solution(angular_flux=psi, scalar_flux=phi_a,
                         mesh=m)
        sol_b = Solution(angular_flux=psi, scalar_flux=phi_b,
                         mesh=m)
        diff = sol_a.compare(sol_b, rtol=1e-12)
        # phi_a - phi_b = -0.1, L∞ = 0.1
        np.testing.assert_allclose(diff.scalar_flux_linf, 0.1, rtol=1e-12)
        assert diff.within_tolerance is False

    def test_compare_keff_only_one_eigenvalue(self) -> None:
        """When one is fixed-source, keff_abs is None."""
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol_a = Solution(angular_flux=psi, scalar_flux=phi,
                         mesh=m, keff=1.0)
        sol_b = Solution(angular_flux=psi, scalar_flux=phi,
                         mesh=m)  # fixed-source
        diff = sol_a.compare(sol_b)
        assert diff.keff_abs is None

    def test_compare_cross_mesh_rejected(self) -> None:
        """Two ``_slab_mesh()`` calls build DIFFERENT constituent objects
        (fresh Mesh1D + Quadrature + materials each), so the phase-space
        pairing guard (``SNMesh.is_same_phase_space``, P6 B2) rejects."""
        m1 = _slab_mesh()
        m2 = _slab_mesh()
        psi_a, phi_a, bf_a = _make_fluxes(m1)
        psi_b, phi_b, bf_b = _make_fluxes(m2)
        sol_a = Solution(angular_flux=psi_a, scalar_flux=phi_a, mesh=m1)
        sol_b = Solution(angular_flux=psi_b, scalar_flux=phi_b, mesh=m2)
        with pytest.raises(ValueError, match="different discrete phase space"):
            sol_a.compare(sol_b)

    def test_compare_accepts_sibling_snmesh_over_same_constituents(self) -> None:
        """The two-entry pattern (P6 B2): two ``SNMesh`` WRAPPERS built from
        the SAME constituent objects (geometry mesh, quadrature, materials)
        realize the same discrete phase space — ``compare`` must accept,
        even though ``mesh_a is not mesh_b``."""
        geometry = Mesh1D(
            edges=np.linspace(0.0, 1.0, 5),
            mat_ids=np.zeros(4, dtype=int),
            coord=CoordSystem.CARTESIAN,
            bc_left=BC("vacuum"),
            bc_right=BC("vacuum"),
        )
        quad = Quadrature.gauss_legendre(n_ordinates=4)
        materials = placeholder_materials(ng=2)
        m1 = SNMesh(geometry, quad, materials)
        m2 = SNMesh(geometry, quad, materials)
        assert m1 is not m2, "the wrappers must be distinct objects"
        assert m1.is_same_phase_space(m2), (
            "sibling SNMesh over the same constituents must share the phase space"
        )
        psi_a, phi_a, _ = _make_fluxes(m1)
        psi_b, phi_b, _ = _make_fluxes(m2)
        sol_a = Solution(angular_flux=psi_a, scalar_flux=phi_a, mesh=m1)
        sol_b = Solution(angular_flux=psi_b, scalar_flux=phi_b, mesh=m2)
        diff = sol_a.compare(sol_b)           # must NOT raise
        assert diff is not None, "compare must return a diff object"


# ════════════════════════════════════════════════════════════════════
# The ROLE axis (#276 A5) — SolutionBase → {Solution, AdjointSolution}
# ════════════════════════════════════════════════════════════════════


class TestSolutionRoleAxis:
    """The A5 carrier ruling's structural contract.

    The role (forward vs adjoint) is a TYPE, the problem kind
    (fixed-source vs eigenvalue) a property — and the forward-physics
    asymmetry is STRUCTURAL: the reaction-rate-preserving operations
    (``homogenize`` / ``condense`` / ``reaction_rate_density``) exist
    only on :class:`Solution`; :class:`AdjointSolution` carries the
    importance vocabulary instead.  These rows pin exactly that shape
    so a future "convenience" re-attachment of forward physics to the
    adjoint (or a role-less base instantiation) is a red test, not a
    silent drift.
    """

    def test_base_not_instantiable(self) -> None:
        """A role-less carrier is not a value that exists."""
        m = _slab_mesh()
        psi, phi, _ = _make_fluxes(m)
        with pytest.raises(TypeError, match="not instantiable"):
            SolutionBase(angular_flux=psi, scalar_flux=phi, mesh=m)

    def test_adjoint_construction_shares_the_carrier(self) -> None:
        """AdjointSolution constructs on the SAME carrier contract
        (fields, mesh-identity validation, boundary delegate,
        problem-kind discrimination) — the role changes semantics,
        never the carrier."""
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        adj = AdjointSolution(angular_flux=psi, scalar_flux=phi, mesh=m)
        assert adj.angular_flux is psi
        assert adj.scalar_flux is phi
        assert adj.boundary_flux is bf
        assert adj.mesh is m
        assert adj.is_fixed_source() and not adj.is_eigenvalue()

    def test_adjoint_space_content_enforced(self) -> None:
        """The base's space-content contract fires on the adjoint leaf too."""
        m1 = _slab_mesh()
        m3 = _slab_mesh(ng=3)
        phi_foreign = ScalarFlux(values=np.zeros((m3.ng, *m3.spatial_shape)), space=m3.bulk_space)
        psi, _, _ = _make_fluxes(m1)
        with pytest.raises(ValueError, match="MARGINAL"):
            AdjointSolution(angular_flux=psi, scalar_flux=phi_foreign, mesh=m1)

    def test_roles_are_siblings_not_subtypes(self) -> None:
        """AdjointSolution is NOT a Solution (and vice versa) — a
        subclass that removed capability would violate Liskov; the
        sibling-under-base shape is what makes the asymmetry legal."""
        assert issubclass(Solution, SolutionBase)
        assert issubclass(AdjointSolution, SolutionBase)
        assert not issubclass(AdjointSolution, Solution)
        assert not issubclass(Solution, AdjointSolution)

    def test_forward_physics_structurally_absent_on_adjoint(self) -> None:
        """The asymmetry is an ABSENT ATTRIBUTE, not a runtime refusal:
        Σ·φ* is not a reaction-rate density and an importance map has
        no reaction rate to preserve — the wrong physics is unspellable."""
        forward_trio = ("homogenize", "condense", "reaction_rate_density")
        for name in forward_trio:
            assert hasattr(Solution, name), f"Solution must carry {name}"
            assert not hasattr(AdjointSolution, name), (
                f"AdjointSolution must NOT carry {name} — the forward-"
                "physics asymmetry is structural (#276 A5)"
            )
        # ...and the vocabulary asymmetry points the other way:
        assert hasattr(AdjointSolution, "importance")
        assert not hasattr(Solution, "importance")

    def test_importance_aliases_scalar_flux(self) -> None:
        """One storage, two vocabularies — φ* IS the importance map."""
        m = _slab_mesh()
        psi, phi, _ = _make_fluxes(m)
        adj = AdjointSolution(angular_flux=psi, scalar_flux=phi, mesh=m)
        assert adj.importance is adj.scalar_flux

    def test_compare_role_mismatch_rejected(self) -> None:
        """Cross-role compare raises even on the SAME mesh — the check
        the mesh-identity guard cannot make (a forward flux and an
        importance map are different physical quantities).

        Leaf-typed calls are already a STATIC type error (``Self``), so
        the runtime guard is exercised the only way it can be reached:
        through base-typed values (``cast`` — pyright re-narrows a mere
        annotation back to the leaf) — exactly the untyped-caller
        surface it exists to defend.
        """
        m = _slab_mesh()
        psi, phi, _ = _make_fluxes(m)
        fwd = cast(SolutionBase, Solution(angular_flux=psi, scalar_flux=phi, mesh=m))
        adj = cast(
            SolutionBase, AdjointSolution(angular_flux=psi, scalar_flux=phi, mesh=m),
        )
        with pytest.raises(TypeError, match="role mismatch"):
            fwd.compare(adj)
        with pytest.raises(TypeError, match="role mismatch"):
            adj.compare(fwd)

    def test_compare_same_role_adjoint_works(self) -> None:
        """Adjoint-vs-adjoint compare flows through the shared base body."""
        m = _slab_mesh()
        psi, phi, _ = _make_fluxes(m)
        adj_a = AdjointSolution(angular_flux=psi, scalar_flux=phi, mesh=m)
        adj_b = AdjointSolution(angular_flux=psi, scalar_flux=phi, mesh=m)
        diff = adj_a.compare(adj_b)
        assert isinstance(diff, SolutionDiff)
        assert diff.within_tolerance
