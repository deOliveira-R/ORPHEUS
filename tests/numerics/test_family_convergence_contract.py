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


# ── the four rows: family, a shipped-defaults run, and its entry point ──
#
# ⛔ This was called ``_CONVERGING`` until 2026-08-11 (#340 N4.7) and the name
# was FALSE for one of its four rows.  `[M]` ``cp-gauss-seidel`` at shipped
# defaults reports ``fully_converged = False``: its thermal-group inner hits
# ``params.max_inner = 100`` at ``dphi_g = 8.387e-05`` against its own
# ``inner_tol = 1e-8`` (rho = 0.950674, needs ~279).  That is the live starved
# inner #340 N4 found, and N4.7 is what made it audible.
#
# Nothing was green-by-luck — neither consumer of this list asserts convergence,
# which is exactly why the wrong name survived review.  The rows are a
# representative run PER FAMILY; convergence is not among their claims, so the
# name no longer implies it.  Rows that must genuinely converge live in
# ``_SILENT_RUNS`` below, where the budget is raised to make it TRUE.
_FAMILY_RUNS = [
    ("cp-jacobi", lambda: _cp(solver_mode="jacobi"), solve_cp),
    ("cp-gauss-seidel", lambda: _cp(solver_mode="gauss_seidel"), solve_cp),
    ("moc", _moc, solve_moc),
    ("diffusion", _diffusion, solve_diffusion_1d),
]

# Rows that DO converge, for the silence half of the audibility pair (H6).
# ⚠ ``cp-gauss-seidel`` needs its budget raised past the shipped default to
# belong here at all — `[M]` 279 projected, so 400 clears with margin. Quoting
# the default would put a KNOWN-starved row in a list named for silence.
_SILENT_RUNS = [
    ("cp-jacobi", lambda: _cp(solver_mode="jacobi")),
    ("cp-gauss-seidel", lambda: _cp(solver_mode="gauss_seidel", max_inner=400)),
    ("moc", _moc),
    ("diffusion", _diffusion),
]

# Rows that starve a level ON PURPOSE — the teeth of every N4.7 gate.
#
# ⭐⭐ The composition is load-bearing, not a convenience. It carries BOTH
# starvation classes:
#
# * **CHILD-starved** (an inner truncates while the outer meets its own
#   criteria) — ``cp-gauss-seidel``, ``moc``;
# * **OUTER-starved** (the top level truncates) — ``cp-jacobi``, ``diffusion``.
#
# `[M]` 2026-08-11, by mutation: restoring the PRE-N6b guard
# (``if record.converged: return``) partitions this list EXACTLY — the two
# child-starved rows fall silent, the two outer-starved rows still fire. ⟹ a
# gate set whose starved fixtures all starve the OUTER is a **provable
# non-catcher for the entire N6b widening** (#340 F1: a converged outer
# standing on a truncated inner), because a failing top level always names
# itself. Any new row added here must state which class it belongs to.
#
# ⚠ ``cp-jacobi`` and ``diffusion`` cannot supply a child-starved row at all:
# Jacobi records no children by design (the scattering source is lagged) and
# diffusion's inner is an exact LU resolvent recorded as DIRECT. That is a fact
# about those methods, not a coverage gap — which is why the other two must.
_STARVED_RUNS = [
    # child-starved
    ("cp-gauss-seidel", lambda: _cp(solver_mode="gauss_seidel", max_inner=2)),
    ("moc", lambda: _moc(max_inner_sweeps=2)),
    # outer-starved
    ("cp-jacobi", lambda: _cp(solver_mode="jacobi", max_outer=3)),
    ("diffusion", lambda: _diffusion(max_outer=3, keff_tol=0.0)),
]


class TestEveryFamilyHandsBackItsIterationTree:
    """The record reaches the caller — the F6 hole, closed."""

    @pytest.mark.parametrize(
        "row_id,run,entry", _FAMILY_RUNS, ids=[r[0] for r in _FAMILY_RUNS],
    )
    def test_the_result_carries_the_record(self, row_id, run, entry) -> None:
        record = run().record
        assert record.label == "outer(power-iteration)"
        assert record.n_iterations >= 1

    @pytest.mark.parametrize(
        "row_id,run,entry", _FAMILY_RUNS, ids=[r[0] for r in _FAMILY_RUNS],
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


# ═══════════════════════════════════════════════════════════════════════
# N4.7 — the tree does not just ANSWER, it TELLS.  Six hazards, gated.
# ═══════════════════════════════════════════════════════════════════════
#
# N4 made a starved solve *answerable* in all four families; N4.7 made it
# *audible*.  The two are different properties and only the second one reaches
# a caller who did not think to look.
#
# ⚠ Scope split, deliberately: the SN legs of the attribution gates live in
# ``tests/sn/solve/test_convergence_contract.py``, which owns the SN fixtures.
# Each module gates the families it can build. The source-level gates (G1, G8)
# are family-agnostic and live here.


class TestTheTwoEmissionPointsStayTWO:
    """⛔ ``ConvergenceWarning`` is raised from exactly TWO places, at two
    different LEVELS, and they must never be consolidated.

    * ENTRY — :func:`~orpheus.numerics.convergence.warn_if_unconverged`:
      *"the tree this call returns is not fully converged"*, once per solve.
    * LEAF — :meth:`~orpheus.numerics.iteration.KrylovAcceleration.solve`:
      *"SciPy's gmres returned info != 0 for THIS inner solve"*, per inner.

    "Unify the ConvergenceWarning emission points" is a natural tidying
    instinct — the category IS shared — and acting on it either drowns the
    tree verdict in per-inner noise or discards the only report of a
    SciPy-level failure.  This is a SOURCE gate on purpose: a behavioural
    "both warnings still fire" test cannot see a refactor that keeps both
    MESSAGES while routing them through one function, whereas the
    ``(module, stacklevel)`` pair can.
    """

    def test_exactly_two_sites_and_they_are_the_named_two(self) -> None:
        import ast
        from pathlib import Path

        import orpheus

        root = Path(orpheus.__file__).parent
        found: set[tuple[str, int | None]] = set()
        for path in sorted(root.rglob("*.py")):
            tree = ast.parse(path.read_text(), filename=str(path))
            for node in ast.walk(tree):
                if not isinstance(node, ast.Call):
                    continue
                fn = node.func
                is_warn = (
                    (isinstance(fn, ast.Attribute) and fn.attr == "warn")
                    or (isinstance(fn, ast.Name) and fn.id == "warn")
                )
                if not is_warn:
                    continue
                names = {
                    a.id for a in node.args if isinstance(a, ast.Name)
                } | {
                    kw.value.id for kw in node.keywords
                    if isinstance(kw.value, ast.Name)
                }
                if "ConvergenceWarning" not in names:
                    continue
                level = next(
                    (
                        kw.value.value
                        for kw in node.keywords
                        if kw.arg == "stacklevel"
                        and isinstance(kw.value, ast.Constant)
                        and isinstance(kw.value.value, int)
                    ),
                    None,
                )
                found.add((str(path.relative_to(root)), level))

        assert found == {
            ("numerics/convergence.py", 3),
            ("numerics/iteration.py", 2),
        }, (
            "the ConvergenceWarning emission set changed. A THIRD site means a "
            "level nobody adjudicated; a MISSING one means a consolidation "
            "that destroyed either the tree verdict or the gmres mechanism "
            f"report (#340 N4.7). Found: {sorted(found)}"
        )

    def test_the_two_messages_are_distinguishable(self) -> None:
        """Same category, disjoint prefixes — so a reader (and a log grep) can
        tell which level is speaking without reading the whole sentence."""
        import warnings

        from orpheus.numerics.convergence import ConvergenceWarning

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            _cp(solver_mode="gauss_seidel", max_inner=2)
        entry = [
            str(w.message) for w in caught
            if issubclass(w.category, ConvergenceWarning)
        ]
        assert len(entry) == 1
        assert entry[0].startswith("solve_cp:")
        assert "BEST-EFFORT" in entry[0]
        # The leaf's own wording, pinned by SUBSTRING so this stays sensitive
        # to the phrase a consumer would match on (coding-standards: grep the
        # SHORTEST distinctive fragment, never the full sentence).
        assert "ERR-053" not in entry[0], (
            "the entry-level warning must not carry the LEAF's mechanism tag; "
            "if it does, the two levels have been merged"
        )


class TestTheWarningBlamesTheCALLER:
    """⭐ ``stacklevel=3`` must point at the user's line, not at ``orpheus``.

    `[M]` 2026-08-11 this was a LIVE defect at 2 of the 8 emission sites:
    ``solve_sn_fixed_source`` dispatches to two private arms and the call sat
    inside them, so frame 3 was the entry's own ``return _solve_...(``
    line and the warning blamed ``orpheus/sn/solver.py``.  Completely ungated
    before this class — ``grep -rn stacklevel tests/`` had one unrelated hit.
    """

    @pytest.mark.parametrize(
        "row_id,run", _STARVED_RUNS, ids=[r[0] for r in _STARVED_RUNS],
    )
    def test_the_warning_is_attributed_outside_orpheus(self, row_id, run):
        """The PORTABLE claim: survives being reached through a helper."""
        import warnings
        from pathlib import Path

        import orpheus
        from orpheus.numerics.convergence import ConvergenceWarning

        root = Path(orpheus.__file__).parent.resolve()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            run()
        conv = [w for w in caught
                if issubclass(w.category, ConvergenceWarning)]
        assert conv, "the fixture must starve a level"
        for w in conv:
            assert not Path(w.filename).resolve().is_relative_to(root), (
                f"{row_id}: the warning blamed {w.filename}:{w.lineno}, which "
                "is inside orpheus/ — the reader is sent to a file they did "
                "not write (#340 N4.7)"
            )

    def test_the_warning_names_the_EXACT_call_line(self) -> None:
        """⭐⭐ The SHARP leg, and the reason the portable one is not enough.

        `[M]` a ``stacklevel`` that is too LARGE escapes ``orpheus/`` too — it
        lands on the caller's *caller*.  The test above is therefore
        structurally incapable of catching an over-deep frame count (it is
        Mode-12 blind to exactly half the failure class), and only an exact
        ``lineno`` can see it.  Measured: ``stacklevel=4`` leaves the portable
        row GREEN while blaming the wrong line.

        ⚠ Every call here must be LITERAL, one line below its ``f_lineno``
        read.  Routing through a helper would pin the HELPER's line, which is
        a true statement about nothing.
        """
        import inspect
        import warnings

        from orpheus.numerics.convergence import ConvergenceWarning

        case, mesh = _slab()
        frame = inspect.currentframe()
        assert frame is not None, "CPython always provides one here"
        expected = frame.f_lineno + 3
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            solve_diffusion_1d(
                case.materials, mesh, max_outer=3, keff_tol=0.0)
        conv = [w for w in caught
                if issubclass(w.category, ConvergenceWarning)]
        assert len(conv) == 1
        assert conv[0].filename == __file__
        assert conv[0].lineno == expected, (
            f"attributed to line {conv[0].lineno}, the call is at {expected} — "
            "stacklevel is off by "
            f"{conv[0].lineno - expected} frames' worth of lines"
        )


class TestOneWarningPerStarvedSolve:
    """H4: exactly ONE, not one per truncated level.

    A warning that always fires is filtered, and the next truncation goes
    unnoticed (the MA4 lesson).  ``== 1`` is the whole gate — ``>= 1`` is what
    the SN module already asserts and it cannot see a per-level emission.
    `[M]` CP's tree carries up to ``ng`` truncated children and still warns
    once.
    """

    @pytest.mark.parametrize(
        "row_id,run", _STARVED_RUNS, ids=[r[0] for r in _STARVED_RUNS],
    )
    def test_exactly_one(self, row_id, run) -> None:
        import warnings

        from orpheus.numerics.convergence import ConvergenceWarning

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = run()
        conv = [w for w in caught
                if issubclass(w.category, ConvergenceWarning)]
        n_truncated = sum(
            1 for lvl in result.record.walk() if lvl.status == "TRUNCATED"
        )
        assert len(conv) == 1, (
            f"{row_id}: {len(conv)} warnings for {n_truncated} truncated "
            "levels — the entry speaks ONCE, for the level first_failure names"
        )


class TestAConvergedSolveIsSILENT:
    """H6, the other half of the audibility pair — and the half that keeps the
    warning worth reading.  Escalate the category to an error and require the
    converging run not to raise."""

    @pytest.mark.parametrize(
        "row_id,run", _SILENT_RUNS, ids=[r[0] for r in _SILENT_RUNS],
    )
    def test_a_converged_solve_emits_nothing(self, row_id, run) -> None:
        import warnings

        from orpheus.numerics.convergence import ConvergenceWarning

        with warnings.catch_warnings():
            warnings.simplefilter("error", ConvergenceWarning)
            result = run()
        assert result.record.fully_converged is True, (
            f"{row_id}: fixture drift — this row must CONVERGE or it stops "
            "testing silence and starts testing nothing"
        )


class TestTheAdviceIsFamilyAgnostic:
    """The advice must not guess the caller's variable name.

    ⛔ Until 2026-08-11 it read ``solution.history.fully_converged`` — a guess
    at a local name no library can know, and simply wrong for the three
    families whose entries return a ``*Result``.  SN gates this at
    ``tests/sn/solve/test_convergence_contract.py``; the claim is family-WIDE,
    so it is gated per family here too.
    """

    @pytest.mark.parametrize(
        "row_id,run", _STARVED_RUNS, ids=[r[0] for r in _STARVED_RUNS],
    )
    def test_it_names_the_attribute_and_its_type(self, row_id, run) -> None:
        import warnings

        from orpheus.numerics.convergence import ConvergenceWarning

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            run()
        msg = next(
            str(w.message) for w in caught
            if issubclass(w.category, ConvergenceWarning)
        )
        assert "`fully_converged`" in msg
        assert "IterationRecord" in msg
        assert "solution" not in msg, (
            f"{row_id}: the advice named a caller-side variable; there is no "
            "`solution` in this family at all (#340 N4.7)"
        )


class TestEveryEntryCallsTheHelperCorrectly:
    """The two call-site invariants, as a source gate.

    1. **The call is at a PUBLIC entry.** ``stacklevel=3`` is only correct one
       frame below the entry; wrapping the call in a private helper silently
       re-blames ``orpheus``.  This is the static form of the live defect the
       attribution class above measures dynamically.
    2. **``balance_defect`` is always EXPLICIT.**  An omitted argument cannot
       be told apart from a forgotten one, and SN's five sites DO have a
       number to pass — so "no balance defect" must be a stated claim rather
       than a silence.
    """

    def test_every_call_site_is_public_and_states_its_balance(self) -> None:
        import ast
        from pathlib import Path

        import orpheus

        root = Path(orpheus.__file__).parent
        sites: list[tuple[str, str, bool]] = []
        for pkg in ("sn", "cp", "moc", "diffusion"):
            path = root / pkg / "solver.py"
            tree = ast.parse(path.read_text(), filename=str(path))
            for fn in ast.walk(tree):
                if not isinstance(fn, ast.FunctionDef):
                    continue
                for node in ast.walk(fn):
                    if (
                        isinstance(node, ast.Call)
                        and isinstance(node.func, ast.Name)
                        and node.func.id == "warn_if_unconverged"
                    ):
                        explicit = any(
                            kw.arg == "balance_defect" for kw in node.keywords
                        )
                        sites.append((pkg, fn.name, explicit))

        # ⚠ SEVEN sites serve EIGHT (entry x inner-arm) combinations, and the
        # gap is the point: ``solve_sn_fixed_source`` has two inner arms
        # (source_iteration / krylov) and ONE emission, because #340 N4.7
        # hoisted the call out of both private arms into the entry. Counting
        # arms instead of sites is how this assertion was first written wrong.
        assert len(sites) == 7, (
            "expected 7 emission sites — 4 SN entries (solve_sn, "
            "solve_sn_adjoint, solve_sn_adjoint_fixed_source, "
            "solve_sn_fixed_source) + cp + moc + diffusion. A count of 8 "
            "usually means the fixed-source arms each warn again; a count of 6 "
            f"means an entry went silent. Found {len(sites)}: {sites}"
        )
        private = [(p, f) for p, f, _ in sites if f.startswith("_")]
        assert not private, (
            f"warn_if_unconverged is called from a PRIVATE function: {private}."
            " stacklevel=3 then blames the public entry's own dispatch line "
            "instead of the caller (#340 N4.7)"
        )
        silent = [(p, f) for p, f, ok in sites if not ok]
        assert not silent, (
            f"these sites omit balance_defect= instead of stating it: {silent}"
        )
