"""The #340 convergence contract, end-to-end, for EVERY production family.

Whether a solve can be trusted is a property of the **iteration**, not of the
**method** — so the contract is gated once, over all families, rather than
re-derived per package.  Three modules divide the #340 surface:

* ``tests/numerics/test_power_iteration_record.py`` — the loop and the
  protocol in isolation, against a SCRIPTED solver;
* ``tests/sn/solve/test_convergence_contract.py`` — SN end-to-end, including
  the warning text and the exit balance projection;
* **this module** — CP / MoC / diffusion end-to-end (plan step N4): each
  result carries its record, every level that can exhaust advises a knob its
  caller can actually reach, the tree has the SHAPE that family's design
  implies, and a deliberately starved run is AUDIBLE rather than silently
  certified.

⭐ The teeth are the starvation rows.  A gate that only checks a converging
solve cannot tell "``fully_converged`` answers the tree" from
"``fully_converged`` is stuck on True" — the positive reading does not
discriminate (``vv-principles`` #19).  Every family therefore appears twice:
once converging, once starved at a named level.

Runtime note (``vv-principles`` Mode 8): this module is COLLECTED, so pytest's
assertion rewriter keeps bare ``assert`` live under the canonical
``python -O``.  The shared knob reference lives in
``tests/_harness/predicates.py``, which is NOT collected and therefore uses
explicit ``raise``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.cp.solver import CPParams, solve_cp
from orpheus.derivations import get
from orpheus.diffusion.solver import solve_diffusion_1d
from orpheus.geometry import CoordSystem, Mesh1D
from orpheus.moc.solver import solve_moc
from orpheus.numerics.eigenvalue import RecordingSolver

from tests._harness.predicates import (
    assert_every_budgeted_level_names_a_reachable_knob,
    reachable_knobs,
)

pytestmark = pytest.mark.foundation


# ── fixtures: the smallest shipped case per family ──────────────────────

_SLAB_CASE = "cp_slab_2eg_2rg"
_MOC_CASE = "moc_cyl1D_1eg_1rg"   # 1-group: the cheap certified MoC solve


def _slab():
    """The shipped 2-group / 2-region slab, as CP's own tests build it."""
    case = get(_SLAB_CASE)
    gp = case.geom_params
    edges = np.concatenate([[0.0], np.cumsum(np.array(gp["thicknesses"]))])
    return case, Mesh1D(
        edges=edges, mat_ids=np.array(gp["mat_ids"]),
        coord=CoordSystem.CARTESIAN,
    )


def _ws_pin():
    """Single-region Wigner-Seitz cell, as ``tests/moc`` builds it."""
    r_cell = 3.6 / np.sqrt(np.pi)
    return Mesh1D(
        edges=np.array([0.0, r_cell]), mat_ids=np.array([0]),
        coord=CoordSystem.CYLINDRICAL,
    )


def _moc(**kw):
    mix = next(iter(get(_MOC_CASE).materials.values()))
    return solve_moc(
        {0: mix}, _ws_pin(), n_azi=8, n_polar=3, ray_spacing=0.05,
        max_outer=kw.pop("max_outer", 200), **kw,
    )


def _cp(**params):
    case, mesh = _slab()
    return solve_cp(case.materials, mesh, CPParams(**params))


def _diffusion(**kw):
    case, mesh = _slab()
    return solve_diffusion_1d(case.materials, mesh, **kw)


# ── the four rows: family, a converging run, and its entry point ────────

_CONVERGING = [
    ("cp-jacobi", lambda: _cp(solver_mode="jacobi"), solve_cp),
    ("cp-gauss-seidel", lambda: _cp(solver_mode="gauss_seidel"), solve_cp),
    ("moc", _moc, solve_moc),
    ("diffusion", _diffusion, solve_diffusion_1d),
]


class TestEveryFamilyHandsBackItsIterationTree:
    """The record reaches the caller — the F6 hole, closed."""

    @pytest.mark.parametrize(
        "row_id,run,entry", _CONVERGING, ids=[r[0] for r in _CONVERGING],
    )
    def test_the_result_carries_the_record(self, row_id, run, entry) -> None:
        record = run().record
        assert record.label == "outer(power-iteration)"
        assert record.n_iterations >= 1

    @pytest.mark.parametrize(
        "row_id,run,entry", _CONVERGING, ids=[r[0] for r in _CONVERGING],
    )
    def test_every_budgeted_level_advises_a_reachable_knob(
        self, row_id, run, entry,
    ) -> None:
        """The forgettable fact, gated against the API itself.

        `[M]` 2026-08-11 all three families shipped ``budget_name ==
        "max_iter"`` — the default, and a parameter of NO entry — because none
        passed ``budget_name`` to ``power_iteration``.  Latent only because
        the record was dropped before anyone could read it.
        """
        assert_every_budgeted_level_names_a_reachable_knob(run().record, entry)

    def test_cps_knob_is_reachable_THROUGH_its_params_dataclass(self) -> None:
        """The reference has to see into a knob BUNDLE, or it punishes CP.

        ``solve_cp`` exposes no ``max_outer``/``max_inner`` of its own — both
        are fields of :class:`CPParams`.  A reference built from
        ``signature().parameters`` alone would reject CP's correct stamp, so
        this row pins the dotted arm directly rather than leaving it implied
        by the rows above.
        """
        knobs = reachable_knobs(solve_cp)
        assert "params.max_outer" in knobs
        assert "params.max_inner" in knobs
        assert "max_outer" not in knobs, (
            "solve_cp does not expose a bare max_outer; if it ever does, the "
            "dotted arm above stops being the thing under test here"
        )


class TestTheTreeHasTheShapeEachDesignImplies:
    """Per-family structure — the part that is NOT common to all four."""

    def test_cp_jacobi_has_no_inner_level_at_all(self) -> None:
        """Deliberate, not missing: Jacobi LAGS the scattering source.

        Recording a child here would be a lie — nothing at that level
        converges.  Contrast diffusion, whose inner is exact and DOES record.
        """
        record = _cp(solver_mode="jacobi").record
        assert record.children == ()
        assert record.fully_converged is record.converged

    def test_cp_gauss_seidel_fans_out_per_group_with_a_real_criterion(
        self,
    ) -> None:
        """``ng`` children per outer, each judging the residual CP already
        computed and used to throw away (#340 N4)."""
        case, _ = _slab()
        ng = next(iter(case.materials.values())).ng
        record = _cp(solver_mode="gauss_seidel").record
        assert len(record.children) == ng * record.n_iterations
        for child in record.children:
            assert child.label.startswith("inner(within-group scatter, g=")
            assert [c.name for c in child.criteria] == ["dphi_g"]
            assert child.criteria[0].trajectory, (
                "a count-only child cannot claim convergence — capturing the "
                "residual is the whole point of this arm"
            )

    def test_diffusions_inner_is_DIRECT_because_it_is_exact(self) -> None:
        """One LU back-substitution: did not iterate, DID converge.

        The child exists so that "diffusion has no inner level" is
        distinguishable from "diffusion's inner was never recorded".
        """
        record = _diffusion().record
        assert len(record.children) == record.n_iterations
        for child in record.children:
            assert child.status == "DIRECT"
            assert child.converged is True
            assert child.budget == 0, (
                "an exact solve has no iteration budget, and a nonzero one "
                "would make it advise a knob that cannot help"
            )

    def test_moc_declares_BOTH_readings_its_inner_is_driven_by(self) -> None:
        """⭐ The Mode-12 pin — the reason this row exists at all.

        The obvious criterion is ``dphi``, and it is nearly blind to the mode
        this loop exists to converge.  `[M]` 2026-08-11, ``moc_cyl1D_1eg_1rg``
        from a cold boundary flux: ``‖Δφ‖`` is MACHINE ZERO from sweep 4 and
        clears ``1e-5`` at sweep **2**, while ``‖Δψ_b‖`` clears only at sweep
        **18** — the scalar flux is a volume moment, so a 1-group problem
        collapses it in one pass while the cyclic-track closure is untouched.
        Dropping ``dpsi_boundary`` would return with MoC's own convergence
        claim four orders unfinished, and no eigenvalue gate would notice.
        """
        record = _moc().record
        assert len(record.children) == record.n_iterations
        for child in record.children:
            assert child.label == "inner(transport sweeps)"
            assert [c.name for c in child.criteria] == [
                "dphi", "dpsi_boundary",
            ]

    @pytest.mark.parametrize(
        "row_id,cls_path",
        [
            ("cp", ("orpheus.cp.solver", "CPSolver")),
            ("moc", ("orpheus.moc.core", "MOCSolver")),
            ("diffusion", ("orpheus.diffusion.solver", "DiffusionSolver")),
        ],
    )
    def test_the_solver_declares_the_member_and_no_second_verdict(
        self, row_id, cls_path,
    ):
        """Without ``inner_records`` ``power_iteration`` takes its ``else ()``
        arm and the family's tree is silently one deep.

        ⚠ This is a CLASS-level ``hasattr``, a weaker question than the one
        production asks — ``power_iteration`` narrows with
        ``isinstance(solver, RecordingSolver)`` on a live instance, and
        ``issubclass`` against that Protocol raises ``TypeError`` outright
        because ``inner_records`` is a non-method member (measured while
        writing this module).  What actually pins the narrowing is
        CONSEQUENCE: were it failing, ``children`` would be empty and the
        three shape rows above would red.  This row's own job is the cheap
        structural half — the member exists, and the retired second stopping
        surface has not grown back.
        """
        import importlib

        module, name = cls_path
        cls = getattr(importlib.import_module(module), name)
        assert hasattr(cls, "inner_records")
        assert not hasattr(cls, "converged"), (
            "the N2b retirement: a solver must not carry a second stopping "
            "surface (pinned for all five realizers in "
            "test_power_iteration_record.py)"
        )

    def test_the_narrowing_production_uses_admits_a_live_solver(self) -> None:
        """The instance-level question, asked the way ``power_iteration`` asks
        it: ``isinstance``, on a real solver.

        Built directly rather than through ``solve_diffusion_1d`` because the
        entry point does not hand its solver back, and diffusion's is the
        cheapest of the three to construct standalone.
        """
        from orpheus.diffusion.augmented_mesh import DiffusionMesh
        from orpheus.diffusion.solver import DiffusionSolver

        case, mesh = _slab()
        solver = DiffusionSolver(DiffusionMesh(mesh, case.materials))
        assert isinstance(solver, RecordingSolver)
        assert list(solver.inner_records) == [], (
            "a fresh solver has run no inner solves; entries here would mean "
            "the accumulator is shared across instances"
        )


class TestAStarvedLevelIsAudibleInEveryFamily:
    """⭐ The teeth. A converging solve cannot show that ``fully_converged``
    answers the TREE rather than being stuck on ``True``."""

    def test_cp_gauss_seidel_starved_inner(self) -> None:
        record = _cp(solver_mode="gauss_seidel", max_inner=2).record
        failure = record.first_failure
        assert record.fully_converged is False
        assert failure is not None
        assert failure.label.startswith("inner(within-group scatter, g=")
        assert failure.status == "TRUNCATED"
        assert failure.budget_name == "params.max_inner"

    def test_moc_starved_inner(self) -> None:
        record = _moc(max_inner_sweeps=2).record
        failure = record.first_failure
        assert record.fully_converged is False
        assert failure is not None
        assert failure.label == "inner(transport sweeps)"
        assert failure.status == "TRUNCATED"
        assert failure.budget_name == "max_inner_sweeps"

    @pytest.mark.parametrize(
        "row_id,run",
        [
            ("cp-jacobi", lambda: _cp(solver_mode="jacobi", max_outer=3)),
            # `[M]` diffusion CONVERGES in exactly 3 outers here, and
            # MINIMUM_OUTER_ITERATIONS is 3, so NO budget can starve it —
            # the tolerance has to be the unreachable thing.  `1e-15` is
            # not enough either: an exact resolvent drives dk to ~1e-16
            # immediately.  `0.0` is StoppingCriterion's documented
            # unsatisfiable-by-a-strict-test input (both measured while
            # writing this row, which is why the comment exists).
            ("diffusion", lambda: _diffusion(max_outer=3, keff_tol=0.0)),
        ],
    )
    def test_a_starved_OUTER_is_audible(self, row_id, run) -> None:
        """Diffusion's inner is exact and cannot starve, so its only
        failable level is the outer — which is the true statement about the
        method, and still has to be reported."""
        record = run().record
        assert record.fully_converged is False
        assert record.first_failure is record, (
            "children are searched first, so the outer being named means no "
            "inner failed — for diffusion that is structural"
        )
        assert record.status == "TRUNCATED"
        assert record.budget_name in reachable_knobs(
            solve_cp if row_id == "cp-jacobi" else solve_diffusion_1d
        )
