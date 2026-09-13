r"""Consumers-campaign step 2 — PRE-carve anchors for the terminal object.

Campaign: ``.claude/plans/cs4c_binding_design.md`` §27.2 (R-cc6).
Verification plan: ``scratch/_consumers/test_architect_step2.md``.

What step 2 does, in the domain's terms
=======================================

A Problem builds its **terminal object** — for the eigenvalue kind the pencil
:math:`(A, F)` on ONE space — as its LAST step; a Strategy picks a point or a
path in :math:`\Lambda`, or a way to invert :math:`\mathcal A(\lambda)` there;
a Solution reproduces itself from the pair.  R-cc6 lands that in three
bit-identical-gated sub-steps:

(i)   split :class:`~orpheus.sn.coupled_system.WithinGroupSystem` into the
      Problem's record (``loss`` + ``space`` = the terminal object) and the
      Strategy's record (``implicit_operator`` + ``explicit_gains`` = ONE
      chosen splitting);
(ii)  adopt ONE posed ``F`` (gated in
      ``tests/sn/operators/test_step2_posed_fission_anchors.py``);
(iii) move the pencil onto the hub, so
      :func:`~orpheus.sn.coupled_system.build_within_group_system` runs once
      per **Problem** instead of once per **outer step**.

The three row kinds in this module, and why each exists
=======================================================

**RECORD** (``TestRecord*``) — green at the time of writing, describing the
pre-carve tree exactly.  Each is *designed to RED at the carve* and is then
**DELETED, never repaired**: its job is to make the API change LOUD, because
a ``strict`` xfail only flips on XPASS and is therefore SILENT when the carve
lands the API with the wrong semantics (``vv-principles`` Mode-8, fourth
class; the coda's pairing shape, ``lessons`` L82e).  ✅ Sub-step (i) landed
2026-09-13: ``TestRecordTheRecordsShape`` (three rows) is deleted, its
subject — the record's field set — now asserted by the ex-xfail row
:class:`TestRuledTheRecordSplits` as a permanent negative gate.

**THEOREM** (``TestLaw*``) — green before AND after.  These are the laws the
split makes load-bearing.  ``A = M − N`` was asserted once per RECORD
(``test_stage_separation.py`` :func:`test_reconstruction_identity_A_equals_M_minus_N`,
schedule-free); since the split the Strategy value is produced *per schedule*
(:meth:`~orpheus.sn.splitting.Splitting.from_schedule`), so the law is
asserted **per Strategy value** — this module's
:class:`TestLawTheSplittingLawHoldsPerStrategy`, the successor to the strict
xfail the split flipped (see the next paragraph).

**XFAIL(strict)** (``TestRuled*``) — the ruled post-carve behaviour, RED
until its sub-step lands; each is paired with the RECORD row that states the
pre-carve answer.  (i)'s marker is gone (XPASSed); (ii) and (iii) stand.

⭐ The row the split FLIPPED already shipped, and the split DEMOTED it
====================================================================

``test_stage_separation.py::test_driver_consumes_the_records_own_splitting[cart2d-gauss_seidel]``
was ``xfail(strict=True)`` with reason R7 and its module said *"the strict-xfail
set IS the campaign's todo list"*.  Sub-step (i) made it XPASS — and made it
a **tautology**, because after the split the driver consumes exactly the value
it is handed (``coding-standards``' single-sourcing demotion).  Its honest
successor is :class:`TestLawTheSplittingLawHoldsPerStrategy` below; the gate
was kept and its docstring re-scoped when the marker went.

⚠ Mode-12, MEASURED, and it decides where every (i) row reads
=============================================================

``[M]`` ``scratch/_consumers/probes2/p7_split_refutation.py`` on
:func:`cart2d_seedless`: the two splittings' implicit operators satisfy
``|M_jacobi·x − M_gs·x| = 0.000000e+00`` on the **bulk** and
``9.970929e-01`` on the **trace**.  They differ ONLY in the boundary block.
⟹ *any gate on sub-step (i) that reads a bulk-only functional is a provable
non-catcher for the splitting choice.*  Every row here that means to see the
choice asserts on the trace, and says so.

Marks
=====

``foundation`` — architecture invariants of the operator algebra; no theory
``:label:``, hence no ``verifies(...)`` (the verifies ⊥ level doctrine).
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

import orpheus.sn.coupled_system as _cs_mod
import orpheus.sn.solver as _solver_mod
from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry import BC, Mesh1D, Mesh2D
from orpheus.numerics.coupled_system import CoupledField
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.coupled_system import (
    WithinGroupSystem,
    build_within_group_system,
)
from orpheus.sn.solver import solve_sn, solve_sn_adjoint
from orpheus.sn.splitting import Splitting, resolve_schedule
from tests.sn.architecture._config import (
    cart2d_seedless,
    record_for,
    random_state,
    sphere_carrying,
    system_a,
)

pytestmark = pytest.mark.foundation

_SEED = 20260912

_SCHEDULES = ("jacobi", "gauss_seidel")

#: The two fields R-cc6 (i) moved OFF the Problem's record onto the
#: Strategy's value (:class:`~orpheus.sn.splitting.Splitting`, 2026-09-13).
_STRATEGY_FIELDS = frozenset({"implicit_operator", "explicit_gains"})


# ═══════════════════════════════════════════════════════════════════════
# Fissile fixtures — why NOT the neighbouring _config meshes
# ═══════════════════════════════════════════════════════════════════════
#
# ``lessons`` L26 / L7: a brief's "reuse the existing fixture" is a
# hypothesis.  ``[M]`` 2026-09-12 — every ``tests/sn/architecture/_config``
# mesh is built from ``_two_region_materials`` (``anisotropic_mixture``), which
# carries NO fission: ``solve_sn`` on ``slab_seedless()`` raises
# ``RuntimeError: leakage scale bridge is degenerate: the last inner solve's
# flux carries non-positive fission production``.  Those meshes are STRUCTURE
# fixtures (operator shapes, splitting laws) and are used as such above; every
# row here that SOLVES builds its own fissile mesh from the shipped
# cross-section library instead.


def _fissile_slab() -> "tuple[dict[int, Mixture], Mesh1D, Quadrature, int]":
    """A solvable 1-D slab: fuel | moderator, 2G, GL-8, reflective|vacuum.

    ``[M]`` ``probes2/p1_counts.py`` — converges in 4 outer steps at the
    default tolerances, so the ``n_outer > 1`` non-vacuity guard holds.
    """
    from orpheus.derivations.common.xs_library import get_mixture

    materials = {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, 9),
        mat_ids=np.array([0, 0, 0, 0, 1, 1, 1, 1]),
        bc_left=BC("reflective"), bc_right=BC("vacuum"),
    )
    return materials, mesh, Quadrature.gauss_legendre(n_ordinates=8), 0


def _fissile_cart2d() -> "tuple[dict[int, Mixture], Mesh2D, Quadrature, int]":
    """A solvable 2-D Cartesian box, heterogeneous, P1, ``nx != ny``, LS-4.

    The geometry choices mirror :func:`cart2d_seedless` for the SAME reasons
    (``nx != ny`` is the ``x↔y``-swap catcher; ``level_symmetric`` avoids the
    ERR-056 axis-aligned degeneracy) — only the materials differ, because the
    Strategy predicate needs an eigenvalue and ``_two_region_materials`` has
    no fission.  ``[M]`` converges in 17 outer steps per schedule, 6.50 s for
    the pair.
    """
    from orpheus.derivations.common.xs_library import get_mixture

    nx, ny = 4, 5
    mat_map = np.zeros((nx, ny), dtype=int)
    mat_map[nx // 2:, :] = 1
    materials = {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 3.0, ny + 1),
        mat_map=mat_map,
        bc_xmin=BC("reflective"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("reflective"), bc_ymax=BC("vacuum"),
    )
    return materials, mesh, Quadrature.level_symmetric(4), 1


# ═══════════════════════════════════════════════════════════════════════
# The route spy — one construction, two module bindings, self-asserting
# ═══════════════════════════════════════════════════════════════════════


class _BuildSpy:
    r"""Count calls to ``build_within_group_system`` through EVERY module
    binding, and assert its own installation.

    ``lessons`` L46e: a census plugin that rebinds only the defining module
    reports a confident ZERO, because consumers import the symbol into their
    own namespace.  ``[M]`` at ``b0fd3e7e`` the symbol is bound in **three**
    modules (``orpheus.sn.coupled_system``, ``orpheus.sn.solver``,
    ``orpheus.sn``); the two PRODUCTION call-site owners are the first two
    (``solver.py`` :2093/:2254/:2766/:3784/:4041 and ``coupled_system.py``
    :268).  Binding zero modules raises rather than returning 0.
    """

    def __init__(self, monkeypatch: pytest.MonkeyPatch) -> None:
        import sys

        self.calls = 0
        self.sites: list[str] = []
        original = _cs_mod.build_within_group_system

        def spy(*args: object, **kwargs: object) -> WithinGroupSystem:
            self.calls += 1
            frame = sys._getframe(1)
            self.sites.append(
                f"{frame.f_code.co_filename.rsplit('/', 1)[-1]}:"
                f"{frame.f_lineno} {frame.f_code.co_name}",
            )
            return original(*args, **kwargs)  # type: ignore[arg-type]

        bound = 0
        for module in list(sys.modules.values()):
            if module is None:
                continue
            try:
                if getattr(module, "build_within_group_system", None) is original:
                    monkeypatch.setattr(
                        module, "build_within_group_system", spy, raising=False,
                    )
                    bound += 1
            except Exception:  # pragma: no cover - defensive over exotic modules
                continue
        if bound == 0:
            raise RuntimeError(
                "the build spy bound 0 modules — the instrument is dead and "
                "every count it reports is a false zero (lessons L46e).",
            )
        self.bound_modules = bound

    def reset(self) -> None:
        self.calls = 0
        self.sites = []


# ═══════════════════════════════════════════════════════════════════════
# RECORD — today's tree, stated exactly.  DELETE these at the carve.
# ═══════════════════════════════════════════════════════════════════════


class TestRecordTheBuildRoute:
    r"""How many times per solve is the terminal object built, and by whom.

    This is the pre-carve anchor for sub-step (iii).  ``[M]``
    ``probes2/p1_counts.py`` / ``p2_counts2.py`` at ``b0fd3e7e``:

    ==================================  =========  ==========================
    solve                                n_outer    ``build_…`` calls
    ==================================  =========  ==========================
    1-D slab eigenvalue, SI                    4    4  (``solver.py:2093``)
    1-D slab eigenvalue, Krylov                4    4  (``solver.py:2254``)
    1-D sphere eigenvalue (carrying), SI       5    5  (``solver.py:2093``)
    2-D Cartesian eigenvalue, SI               3    3  (``solver.py:2093``)
    1-D slab FIXED SOURCE                    n/a    1  (``solver.py:3784``)
    1-D slab ADJOINT eigenvalue                4    1  (``solver.py:2766``)
    ==================================  =========  ==========================

    ⭐ Two readings the plan turns into design constraints: the eigen count is
    exactly ``n_outer`` (so the row asserts the MECHANISM, not a fixture
    number — a naive "move the call into a method" still reads ``n_outer``
    and this row still reds), and the ADJOINT and FIXED-SOURCE paths are
    ALREADY once-per-Problem, so sub-step (iii) changes the forward
    eigenvalue path only.
    """

    @pytest.mark.parametrize("inner_solver", ["source_iteration", "krylov"])
    def test_record_eigen_builds_once_per_OUTER_step(
        self, monkeypatch: pytest.MonkeyPatch, inner_solver: str,
    ) -> None:
        """RECORD — the forward eigenvalue path builds once per outer step.

        ⛔ DELETE at sub-step (iii).  The ruled successor is
        :meth:`TestRuledTheBuildIsOncePerProblem.test_ruled_eigen_builds_once`.
        """
        materials, mesh, quadrature, order = _fissile_slab()
        spy = _BuildSpy(monkeypatch)
        assert spy.bound_modules >= 2, (
            f"the spy bound only {spy.bound_modules} module(s); production "
            f"calls live in BOTH orpheus.sn.solver and "
            f"orpheus.sn.coupled_system."
        )
        spy.reset()
        solution = solve_sn(
            materials, mesh, quadrature, scattering_order=order,
            max_outer=200, inner_solver=inner_solver,
        )
        history = solution.history
        assert history is not None
        n_outer = len(history.keff_history)
        assert n_outer > 1, (
            f"non-vacuity: the fixture converged in {n_outer} outer step(s), "
            f"so 'once per outer' and 'once per Problem' are the same number "
            f"and this row cannot discriminate."
        )
        assert spy.calls == n_outer, (
            f"pre-carve the within-group system is built once per OUTER step: "
            f"expected {n_outer}, saw {spy.calls} at {spy.sites}."
        )

    def test_record_fixed_source_and_adjoint_already_build_once(
        self, monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """RECORD — ⭐ two of the three entries are ALREADY once-per-Problem.

        So sub-step (iii)'s blast radius is the forward eigenvalue path
        alone, and these two counts are a MUST-STAY-GREEN pin across it.
        """
        from orpheus.sn import solve_sn_fixed_source
        from orpheus.sn.solver import _as_sn_mesh

        materials, mesh, quadrature, order = _fissile_slab()
        spy = _BuildSpy(monkeypatch)

        spy.reset()
        adjoint = solve_sn_adjoint(
            materials, mesh, quadrature, scattering_order=order, max_outer=200,
        )
        adjoint_history = adjoint.history
        assert adjoint_history is not None
        n_outer_adj = len(adjoint_history.keff_history)
        assert n_outer_adj > 1, (
            f"non-vacuity: the adjoint converged in {n_outer_adj} outer "
            f"step(s), so once-per-outer and once-per-Problem coincide."
        )
        assert spy.calls == 1, (
            f"the adjoint poses ONCE ({spy.sites}); it converged in "
            f"{n_outer_adj} outer steps."
        )

        spy.reset()
        sn_mesh = _as_sn_mesh(mesh, quadrature, materials, scattering_order=order)
        source = np.ones(sn_mesh.angular_trial_space.shape)
        solve_sn_fixed_source(
            materials, mesh, quadrature, source, scattering_order=order,
        )
        assert spy.calls == 1, (
            f"the fixed-source entry poses ONCE; saw {spy.calls} at "
            f"{spy.sites}."
        )


class TestRecordTheStaleSigmaExposure:
    r"""HAZARD H1, stated as a measurement rather than a warning.

    ``[M]`` ``probes2/p5_hazard_H1.py``: a within-group system built BEFORE
    ``SNSolver.rebind_cross_sections`` and applied AFTER it is **bit-identical
    to its pre-rebind self** (``array_equal = True``), while a freshly built
    one moves ``max|Δ| = 3.7797926845799177``, ``max rel = 0.1619``.  Today
    the defect is unspellable only because ``build_within_group_system`` runs
    per outer step; sub-step (iii) removes exactly that protection.

    ⟹ the cached pencil and a non-mutating ``.at(σ)`` are ONE merge unit
    (plan F-5 / open ruling O-5).  This row is the CHARACTERIZATION that says
    what must become unspellable; it carries no ``verifies`` and asserts a
    one-sided fact about today's tree.
    """

    def test_record_a_reused_system_goes_stale_under_a_sigma_rebind(self) -> None:
        """RECORD — a reused system is stale after a σ rebind, silently."""
        from orpheus.sn.solver import SNSolver
        from tests.sn.architecture._config import slab_seedless

        sn_mesh = slab_seedless()
        solver = SNSolver(sn_mesh)
        reused = build_within_group_system(
            sn_mesh, solver.mat_xs,
            scattering_op=solver.scattering_op, n2n_op=solver.n2n_op,
        )
        state = random_state(reused, seed=_SEED)
        before = reused.loss.apply(state).to_flat().copy()

        solver.rebind_cross_sections(
            np.asarray(solver.mat_xs.total_cross_section) * 3.0,
        )
        after_reused = reused.loss.apply(state).to_flat()
        rebuilt = build_within_group_system(
            sn_mesh, solver.mat_xs,
            scattering_op=solver.scattering_op, n2n_op=solver.n2n_op,
        )
        after_fresh = rebuilt.loss.apply(state).to_flat()

        np.testing.assert_array_equal(
            before, after_reused,
            err_msg=(
                "the reused system MOVED under a σ rebind — if this reds, "
                "the rebind stopped being a silent staleness source and the "
                "H1 merge-unit ruling can be revisited."
            ),
        )
        defect = float(np.max(np.abs(after_reused - after_fresh)))
        scale = float(np.max(np.abs(after_fresh)))
        assert defect / scale > 1e-2, (
            f"the stale/fresh gap is only {defect / scale:.3e} relative — "
            f"the fixture no longer activates the σ_t path, so this "
            f"characterization has lost its subject."
        )


# ═══════════════════════════════════════════════════════════════════════
# THEOREM — green before AND after.  These are the step's real contracts.
# ═══════════════════════════════════════════════════════════════════════


def _splitting_image(
    record: WithinGroupSystem, schedule: str, sn_mesh: object, state: CoupledField,
) -> np.ndarray:
    """``(M − ΣN_i)·x`` for the Strategy value ``schedule`` labels.

    Stated on the ANGULAR bindings (the value's ``implicit``/``explicit``
    directly, not ``_within_group_si``): the law is about the OPERATORS, and
    the 2-D driver's moment re-binding is a separate claim that
    ``test_stage_separation.py`` already owns.
    """
    splitting = Splitting.from_schedule(
        record, resolve_schedule(sn_mesh, schedule),  # type: ignore[arg-type]
    )
    interior = system_a(state)
    accumulated = None
    for gain in splitting.explicit:
        image = gain.apply(interior)
        accumulated = image if accumulated is None else accumulated + image
    assert accumulated is not None
    return (
        CoupledField(systems=(splitting.implicit.apply(interior),)).to_flat()
        - CoupledField(systems=(accumulated,)).to_flat()
    )


class TestLawTheSplittingLawHoldsPerStrategy:
    r"""``A = M − N`` **for every Strategy value** — the successor to the R7
    strict xfail that sub-step (i) flips.

    Before the split ``test_stage_separation.py``'s law row was parametrized
    over the RECORD, which is schedule-free; the choice lived downstream in
    ``_select_si_splitting``.  Since the split the Strategy value is produced
    per schedule, so the law's denominator is the Strategy's value set — and
    the gate that flipped (*"the driver runs the objects the record
    advertises"*) is true BY CONSTRUCTION.  This is the row that keeps
    teeth; the value's own :meth:`~orpheus.sn.splitting.Splitting.law_residual`
    is the same statement as a method, asserted at ``0.0`` alongside (exact
    on the seedless arm — one flat operator sum; the carrying arm's block
    grid re-associates and is gated at nulp in ``test_stage_separation.py``).

    ``[M]`` ``probes2/p9_gs_law.py`` on :func:`cart2d_seedless`: the law is
    ``array_equal`` (``max|A − (M − N)| = 0.000000e+00``) for **both**
    schedules; the in-class ERR-056-family mutation (the lower boundary half
    kept implicit in ``M`` *and* lagged whole in ``N``) reads
    **``3.595068e+00``**.  Bit-exactness is EARNED here, not assumed: both
    sides are one flat operator sum over the same leaves.
    """

    @pytest.mark.parametrize("schedule", _SCHEDULES)
    def test_law_A_equals_M_minus_N_per_schedule(self, schedule: str) -> None:
        sn_mesh = cart2d_seedless()
        record = record_for(sn_mesh)
        state = random_state(record, seed=_SEED)
        loss_image = record.loss.apply(state).to_flat()
        splitting = Splitting.from_schedule(record, resolve_schedule(sn_mesh, schedule))
        np.testing.assert_array_equal(
            splitting.law_residual(state), np.zeros_like(loss_image),
            err_msg=f"Splitting.law_residual is not exactly zero for {schedule!r}",
        )
        np.testing.assert_array_equal(
            loss_image, _splitting_image(record, schedule, sn_mesh, state),
            err_msg=(
                f"A != M - N for inner_schedule={schedule!r}. Both sides are "
                f"one flat operator sum over the same leaves, so a non-zero "
                f"defect is an algebraic error, not FP drift (measured 0.0 "
                f"on both schedules at b0fd3e7e)."
            ),
        )

    def test_law_the_two_schedules_really_choose_different_splittings(
        self,
    ) -> None:
        r"""ACTIVATION leg — without it the row above is two copies of one
        measurement.

        ⚠ Mode-12: ``[M]`` the two implicit operators agree **exactly** on
        the BULK (``0.000000e+00``) and differ on the TRACE
        (``9.970929e-01``).  This row therefore reads the trace; a bulk
        assertion here would be a provable non-catcher.
        """
        sn_mesh = cart2d_seedless()
        record = record_for(sn_mesh)
        state = system_a(random_state(record, seed=_SEED))
        images = {}
        for schedule in _SCHEDULES:
            implicit = Splitting.from_schedule(
                record, resolve_schedule(sn_mesh, schedule),
            ).implicit
            out = implicit.apply(state)
            images[schedule] = (
                np.asarray(out.interior.values).copy(),
                np.asarray(out.boundary.values).copy(),
            )
        bulk_gap = float(
            np.max(np.abs(images["jacobi"][0] - images["gauss_seidel"][0])),
        )
        trace_gap = float(
            np.max(np.abs(images["jacobi"][1] - images["gauss_seidel"][1])),
        )
        np.testing.assert_array_equal(
            images["jacobi"][0], images["gauss_seidel"][0],
            err_msg=(
                "the two splittings moved the BULK — the declared Mode-12 "
                "blindness (the choice is a boundary-block choice) is false "
                "and every 'bulk-only is blind' note in this module must be "
                "re-measured."
            ),
        )
        assert trace_gap > 1e-2, (
            f"the two schedules produced trace-indistinguishable implicit "
            f"operators (gap {trace_gap:.3e}; bulk {bulk_gap:.3e}) — the "
            f"fixture stopped reaching the G-S arm, so the law row above is "
            f"two copies of one measurement."
        )


class TestLawTheStrategyPredicate:
    r"""The three-type gate's **Strategy** clause, run: *nothing on it moves
    the limit past tolerance*.

    ``[M]`` ``probes2/p7_split_refutation.py`` at
    ``keff_tol=1e-11 / flux_tol=1e-10 / inner_tol=1e-12`` on
    :func:`cart2d_seedless`: ``k_jacobi = 0.19266749853436`` vs
    ``k_gs = 0.19266749853430`` (``|Δk| = 6.140e-14``), ``max rel |Δφ| =
    4.061e-13`` (17 outer steps each, 6.50 s for the pair) — while the splitting objects differ by ``9.97e-01`` on the
    trace.  That pair IS the ``vv`` Mode-9 discriminator, and it is what makes
    ``inner_schedule`` a Strategy coordinate rather than Problem data.

    ⚠ Mode-9's own warning applies: a splitting must be verified on a config
    that BREAKS the degenerate coincidence.  This fixture is heterogeneous,
    2-group, mixed-BC, ``nx != ny``, level-symmetric — chosen by
    ``_config.cart2d_seedless`` precisely for the ERR-056 non-degeneracy.
    """

    def test_law_two_schedules_one_limit(self) -> None:
        materials, mesh, quadrature, order = _fissile_cart2d()
        results = {}
        for schedule in _SCHEDULES:
            solution = solve_sn(
                materials, mesh, quadrature, scattering_order=order,
                max_outer=400, inner_schedule=schedule,
                keff_tol=1e-11, flux_tol=1e-10, inner_tol=1e-12,
            )
            results[schedule] = (
                solution.keff, np.asarray(solution.scalar_flux.values).copy(),
            )
        delta_k = abs(results["jacobi"][0] - results["gauss_seidel"][0])
        scale = float(np.max(np.abs(results["jacobi"][1])))
        delta_phi = float(
            np.max(np.abs(results["jacobi"][1] - results["gauss_seidel"][1])),
        ) / scale
        assert delta_k < 1e-11, (
            f"|Δk| = {delta_k:.3e} between the two schedules exceeds the "
            f"keff tolerance they were both driven to — the splitting is "
            f"NOT limit-invariant and `inner_schedule` is not a Strategy "
            f"coordinate (measured 6.140e-14 on this fixture at b0fd3e7e)."
        )
        assert delta_phi < 1e-10, (
            f"max rel |Δφ| = {delta_phi:.3e} exceeds the flux tolerance "
            f"(measured 4.061e-13 on this fixture at b0fd3e7e)."
        )


class TestLawTheGaugeIsSigmaFree:
    r"""HAZARD H3 — ``loss_kernel_gauge`` is σ-FREE while ``loss`` is not.

    ``[M]`` ``probes2/p4_spaces_and_seeds.py``: after a ×3 σ_t rebind the hub
    returns the SAME ``LossKernelGauge`` object (it is a ``cached_property``
    over σ-free data, ``sn/mesh/augmented_mesh.py`` :1062), while a freshly
    built ``loss`` moves by ``max rel = 0.1619``.

    ⟹ deriving the gauge from the pencil would make it rebuild on every
    ``.at(σ)`` — invisible to every value test, because the rebuilt gauge has
    the same values.  This row pins the σ-freedom so that a later "derive it
    from the pencil" simplification reds on a COUNT, not on a value.
    """

    def test_law_the_gauge_survives_a_sigma_rebind_by_identity(self) -> None:
        from orpheus.sn.solver import SNSolver
        from tests.sn.architecture._config import slab_seedless

        sn_mesh = slab_seedless()
        solver = SNSolver(sn_mesh)
        gauge_before = sn_mesh.loss_kernel_gauge
        solver.rebind_cross_sections(
            np.asarray(solver.mat_xs.total_cross_section) * 3.0,
        )
        assert sn_mesh.loss_kernel_gauge is gauge_before, (
            "the loss-kernel gauge was rebuilt by a σ rebind — it is "
            "documented σ-free; if this reds, either the gauge gained a σ "
            "read or it was derived from the pencil (HAZARD H3)."
        )


# ═══════════════════════════════════════════════════════════════════════
# XFAIL(strict) — the ruled post-carve behaviour.  RED today, by design.
# ═══════════════════════════════════════════════════════════════════════


class TestRuledTheRecordSplits:
    """R-cc6 (i) — the Problem's record carries no splitting.

    ✅ LANDED 2026-09-13 (the consumers campaign's step 2, C1): the strict
    xfail this row carried XPASSed and was deleted; the row stays as the
    PERMANENT negative gate — a Strategy field re-appearing on the posed
    record is the weld this step removed.
    """

    def test_ruled_the_posed_record_carries_no_strategy_field(self) -> None:
        """The record ``build_within_group_system`` returns holds the loss,
        its space and its factors — never a splitting.

        Expressed over the FIELD SET rather than over a guessed class name:
        a row that named the value's class would be a guess wearing an
        assertion; the value is :class:`~orpheus.sn.splitting.Splitting`,
        and it is minted FROM this record, not stored on it.
        """
        record = build_within_group_system(
            cart2d_seedless(), cart2d_seedless().material_xs_field(),
        )
        names = frozenset(f.name for f in dataclasses.fields(record))
        assert names.isdisjoint(_STRATEGY_FIELDS), (
            f"the posed record carries {sorted(names & _STRATEGY_FIELDS)} "
            f"— the splitting is Strategy-side (a value minted from the "
            f"record's factors, never a field on the Problem's record)."
        )
        assert "factors" in names, (
            f"the posed record exposes no factors ({sorted(names)}) — the "
            f"Strategy value is minted from them (the split is not a field "
            f"partition: the factors are not recoverable from `loss`)."
        )


class TestRuledTheBuildIsOncePerProblem:
    """R-cc6 (iii) — the terminal object is built once per Problem."""

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "R-cc6 (iii) NOT LANDED — the forward eigenvalue path builds the "
            "within-group system once per OUTER step (measured: count == "
            "n_outer on slab/sphere/2-D). WHEN THIS XPASSES: the pencil is "
            "on the hub — delete this marker AND "
            "TestRecordTheBuildRoute.test_record_eigen_builds_once_per_OUTER_step."
        ),
    )
    def test_ruled_eigen_builds_once(
        self, monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """EXACTLY once, never ``> 0``.

        The brief's own measurable consequence: *a naive implementation that
        merely moves the call into a method (not a cached_property) still
        reads n_outer and FAILS*.  The non-vacuity guard below is what makes
        that true — on a fixture converging in one outer step the two claims
        coincide.
        """
        materials, mesh, quadrature, order = _fissile_slab()
        spy = _BuildSpy(monkeypatch)
        spy.reset()
        solution = solve_sn(
            materials, mesh, quadrature, scattering_order=order, max_outer=200,
        )
        history = solution.history
        assert history is not None
        n_outer = len(history.keff_history)
        if n_outer <= 1:
            pytest.fail(
                f"non-vacuity: the fixture converged in {n_outer} outer "
                f"step(s), so 'once per Problem' and 'once per outer' are "
                f"the same number and this row proves nothing.",
            )
        assert spy.calls == 1, (
            f"the terminal object was built {spy.calls} times over {n_outer} "
            f"outer steps, at {spy.sites}."
        )
