r"""S4 — the Strategy mints the Green (S\ :sub:`N` rows of the S4 gate spec).

Spec: ``scratch/numerics_untangling/packages/S4_gate_spec.md`` (rows G1–G3;
the numerics rows G4–G6 are ``tests/gates/numerics/test_strategy_mints_green.py``).

* G1 — the ROUTE (design §6.2, ``retirement-audit`` COUNTING-SPY): every
  driver a production entry constructs is constructed by
  :meth:`~orpheus.numerics.green_operator.GreenOperator.driver`, and every
  finalize is one :meth:`~orpheus.numerics.green_operator.GreenOperator.step`
  of the Green the last inner solve drove (ruling 69, #448). The retired
  spellings are gone. Value gates cannot see a route: ``[M]`` the 26
  production configurations are bit-identical before and after the carve
  (``S4_gates/baseline``).
* G2 — the RELAXATION face (rulings 68, Q3): the driven gains are DERIVED from
  N's terms through the realisation's representation, never handed in.
* G3 — the SPLITTING MAP (ruling 69): ``step`` is one application of
  :math:`M^{-1}(q + N\psi)` through the UN-windowed :math:`M^{-1}`, seeded by
  the iterate.
"""
from __future__ import annotations

import ast
import inspect
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from orpheus.numerics import iteration as iteration_mod
from orpheus.numerics.green_operator import GreenOperator
from orpheus.numerics.iteration import KrylovAcceleration, SourceIteration

pytestmark = [pytest.mark.foundation]

_HERE = "tests/gates/sn/solve/test_strategy_mints_green_sn.py"
_ORPHEUS = Path(iteration_mod.__file__).resolve().parents[1]
_GREEN_FILE = (_ORPHEUS / "numerics" / "green_operator.py").resolve()


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


# ── fixtures ────────────────────────────────────────────────────────────────


def _cfg(builder: str) -> dict:
    from tests.gates.sn.regression import _generate_snapshots as G

    return {
        "slab": lambda: G._slab_3region("2g", 8),
        "sphere": lambda: G._sphere_3region("2g", 8),
        "cart2d": lambda: G._cartesian_2d_het_si("2g"),
        "cart2d_p1": lambda: G._cartesian_2d_p1_aniso_het_si("2g"),
    }[builder]()


def _problem(builder: str):
    from orpheus.sn.solver import _as_problem

    cfg = _cfg(builder)
    return _as_problem(cfg["mesh"], cfg["quadrature"], cfg["materials"],
                       scattering_order=cfg["scattering_order"])


def _splitting(problem, schedule: str = "jacobi"):
    from orpheus.sn.splitting import Splitting, resolve_schedule

    return Splitting.from_schedule(problem.system, resolve_schedule(problem, schedule))


def _source(problem, seed: int):
    """A random per-ordinate source lifted by the production rhs builder (the
    coupled pair on a carrying mesh)."""
    from orpheus.sn.solver import _build_fixed_source_rhs

    rng = np.random.default_rng([seed, 4])
    shape = (problem.quad.n_ordinates, problem.ng, *problem.spatial_shape)
    return _build_fixed_source_rhs(rng.uniform(0.5, 1.5, shape), problem)


def _flux(problem, q, seed: int, *, windowed: bool = False):
    """A random flux iterate of the class the production driver holds (the
    moment composite when windowed, paired with ψ_B on a carrying mesh)."""
    from orpheus.sn.solver import _coupled_flux_state, _unwindowed_cold_start, _windowed_cold_start

    depth = (q.systems[0] if hasattr(q, "systems") else q).history_depth
    if windowed:
        cold = _windowed_cold_start(problem.system.factors.scattering, problem, history_depth=depth)
    else:
        cold = _unwindowed_cold_start(problem, history_depth=depth)
        if problem.system.is_coupled:
            cold = _coupled_flux_state(cold, problem)
    rng = np.random.default_rng([seed, 9])
    return type(cold).from_flat(rng.standard_normal(cold.to_flat().size), cold)


def _flat(x: Any) -> np.ndarray:
    return np.asarray(x.to_flat(), dtype=float)


# ═══════════════════════════════════════════════════════════════════════
# G1 — the route
# ═══════════════════════════════════════════════════════════════════════


def _in_green_driver() -> bool:
    """True iff a frame of ``GreenOperator.driver`` DEFINED IN green_operator.py
    is on the stack, keyed by file and QUALIFIED name: a driver built by a
    rebinding defined elsewhere is not recognised, and neither is a frame of
    the realisations' own ``driver`` (same file, same bare name; `[M]` keyed
    on the bare name the hand-built-driver arm read 0 reds once ruling 70 put
    ``Richardson.driver`` in that file)."""
    f = sys._getframe(2)
    while f is not None:
        code = f.f_code
        if code.co_qualname == "GreenOperator.driver" and Path(code.co_filename).resolve() == _GREEN_FILE:
            return True
        f = f.f_back
    return False


def _entry_runs() -> dict:
    """The production entries at loose budgets (the claim is the ROUTE), and the
    number of finalize steps each owes."""
    from orpheus.sn import solver as S

    slab, sph, c2d, c2dp1 = _cfg("slab"), _cfg("sphere"), _cfg("cart2d"), _cfg("cart2d_p1")
    loose = dict(max_outer=3, keff_tol=1e-3, flux_tol=1e-3, max_inner=40, inner_tol=1e-6)

    def eig(cfg, solver, **kw):
        return lambda: S.solve_sn(cfg["materials"], cfg["mesh"], cfg["quadrature"],
                                  inner_solver=solver, **{**loose, **kw})

    def fs(cfg, solver, src, **kw):
        return lambda: S.solve_sn_fixed_source(
            materials=cfg["materials"], mesh=cfg["mesh"], quadrature=cfg["quadrature"],
            external_source=src, scattering_order=cfg["scattering_order"],
            inner_solver=solver, max_inner=40, inner_tol=1e-6, **kw)

    ones = np.ones((8, 2, 8))
    return {
        "eig_si_slab": (eig(slab, "source_iteration"), 1),
        "eig_kr_slab": (eig(slab, "krylov"), 1),
        "eig_si_sphere": (eig(sph, "source_iteration"), 1),
        "eig_kr_sphere": (eig(sph, "krylov"), 1),
        "eig_si_c2d_gs": (eig(c2d, "source_iteration", inner_schedule="gauss_seidel"), 1),
        "fs_si_slab_dsa": (fs(slab, "source_iteration", ones, acceleration="dsa"), 0),
        "fs_kr_slab_dsa": (fs(slab, "krylov", ones, acceleration="dsa"), 0),
        "fs_si_c2d_windowed": (fs(c2dp1, "source_iteration", c2dp1["external_source"]), 1),
        "fs_kr_sphere": (fs(sph, "krylov", ones), 0),
        "adj_eig_slab": (lambda: S.solve_sn_adjoint(
            slab["materials"], slab["mesh"], slab["quadrature"], **loose), 0),
        "adj_fs_slab": (lambda: S.solve_sn_adjoint_fixed_source(
            slab["materials"], slab["mesh"], slab["quadrature"], np.ones((2, 8)),
            max_inner=40, inner_tol=1e-6), 0),
    }


_ENTRIES = list(_entry_runs())


@pytest.mark.parametrize("entry", _ENTRIES)
def test_g1a_every_production_driver_is_built_by_a_green(entry: str, monkeypatch) -> None:
    """Every ``SourceIteration`` / ``KrylovAcceleration`` a production entry
    constructs is constructed under ``GreenOperator.driver``, and at least one
    is (the activation leg). Reds on a site re-inlining
    ``SourceIteration(seeded_inverse(M), *gains)``, on ``KEigenvalue`` building
    its inner by hand, and on a ``driver`` rebound from another file."""
    seen: list[tuple[str, bool]] = []
    for cls in (SourceIteration, KrylovAcceleration):
        real = cls.__init__

        def spy(self, *a, _real=real, _name=cls.__name__, **k):
            seen.append((_name, _in_green_driver()))
            _real(self, *a, **k)

        monkeypatch.setattr(cls, "__init__", spy)
    thunk, _ = _entry_runs()[entry]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        thunk()
    _require(seen, f"[{entry}] Mode-11: no driver was constructed; the spy never fired")
    stray = [name for name, under in seen if not under]
    _require(not stray, f"[{entry}] {len(stray)} of {len(seen)} drivers built outside G.driver: {stray}")


@pytest.mark.parametrize("entry", _ENTRIES)
def test_g1b_every_finalize_is_one_step_of_the_green_the_last_inner_drove(entry: str, monkeypatch) -> None:
    """``GreenOperator.step`` runs once per eigenvalue solve and once per
    WINDOWED fixed-source solve (0 otherwise), on the SAME Green object whose
    driver ran last, so the finalize cannot drift from the iteration (#448).
    Reds on a finalize applying ``M.inverse()`` by hand (count 0) and on a
    finalize re-minting a Green (identity)."""
    drivers: list[int] = []
    steps: list[int] = []
    real_d, real_s = GreenOperator.driver, GreenOperator.step

    def spy_d(self, **kw):
        drivers.append(id(self))
        return real_d(self, **kw)

    def spy_s(self, *a, **kw):
        steps.append(id(self))
        return real_s(self, *a, **kw)

    monkeypatch.setattr(GreenOperator, "driver", spy_d)
    monkeypatch.setattr(GreenOperator, "step", spy_s)
    thunk, expected = _entry_runs()[entry]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        thunk()
    _require(len(steps) == expected, f"[{entry}] step ran {len(steps)}×, expected {expected}")
    if expected:
        _require(drivers and steps[-1] == drivers[-1],
                 f"[{entry}] the finalize stepped a Green other than the last one driven")


def test_g1c_the_retired_spellings_are_gone() -> None:
    """Retirement (``retirement-audit`` A–G). ``fixed_point_step`` is gone from
    ``orpheus.numerics.iteration`` and its ``__all__``; an AST pass over every
    file under ``orpheus/`` finds no ``fixed_point_step``, ``_within_group_si``
    or ``_within_group_krylov`` (the filter's control: it finds
    ``lagged_source``) and no ``SourceIteration(``/``KrylovAcceleration(`` call
    outside ``numerics/`` (three in ``sn/solver.py`` before S4);
    ``KEigenvalue``'s parameters start ``(posing, green)`` with no
    ``implicit``/``explicit``."""
    from orpheus.numerics.iteration import KEigenvalue

    _require(not hasattr(iteration_mod, "fixed_point_step"), "fixed_point_step survives")
    _require("fixed_point_step" not in getattr(iteration_mod, "__all__", ()), "still exported")
    retired = {"fixed_point_step", "_within_group_si", "_within_group_krylov"}
    found, control, ctor_sites, n_files = [], 0, [], 0
    for f in sorted(_ORPHEUS.rglob("*.py")):
        n_files += 1
        rel = f.relative_to(_ORPHEUS)
        for n in ast.walk(ast.parse(f.read_text())):
            name = (n.id if isinstance(n, ast.Name) else n.attr if isinstance(n, ast.Attribute)
                    else n.value if isinstance(n, ast.Constant) and isinstance(n.value, str) else None)
            if name in retired:
                found.append(f"{rel}:{n.lineno}:{name}")
            if name == "lagged_source":
                control += 1
            if (isinstance(n, ast.Call)
                    and getattr(n.func, "id", getattr(n.func, "attr", None)) in ("SourceIteration", "KrylovAcceleration")
                    and "numerics" not in rel.parts):
                ctor_sites.append(f"{rel}:{n.lineno}")
    _require(n_files > 100 and control > 0, f"census filter broken: files={n_files}, control={control}")
    _require(not found, f"retired spellings survive: {found}")
    _require(not ctor_sites, f"drivers constructed outside numerics/: {ctor_sites}")
    params = list(inspect.signature(KEigenvalue).parameters)
    _require(params[:2] == ["posing", "green"] and not ({"implicit", "explicit"} & set(params)),
             f"KEigenvalue signature {params}")


# ═══════════════════════════════════════════════════════════════════════
# G2 — the relaxation face
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("builder,schedule", [("slab", "jacobi"), ("sphere", "jacobi"),
                                              ("cart2d", "jacobi"), ("cart2d", "gauss_seidel")])
def test_g2a_the_driver_face_drives_ns_terms_through_the_representation(builder: str, schedule: str) -> None:
    """The SI Green the Strategy mints: its driver's gains are N's terms
    (``is``, in order) when unwindowed; when windowed every ``AngularLift`` term
    is re-bound on the moment iterate (the stage-separation predicate, field
    by field) and every other term is N's own object. Reds when the re-binding
    is skipped, HERE on the driver face and not only through ``apply`` (the
    orchestrator's Q3 (iii))."""
    from orpheus.sn.solver import _within_group_green_si
    from orpheus.transport.operators.angular_lift import AngularLift
    from tests.gates.sn.architecture.test_stage_separation import (
        _is_the_records_gain_rebound_on_the_iterate as rebound,
    )

    problem = _problem(builder)
    sp = _splitting(problem, schedule)
    green, windowed = _within_group_green_si(sp, problem)
    driver = green.driver(max_iter=1, tol=1e-12, budget_name="max_inner")
    _require(windowed == (builder == "cart2d"), f"windowing premise moved: windowed={windowed}")
    _require(tuple(map(id, driver.gains)) == tuple(map(id, green.driven_gains)),
             "the driver does not drive the Green's own derived gains")
    _require(len(driver.gains) == len(sp.explicit), "driven gain count ≠ N's term count")
    for d, r in zip(driver.gains, sp.explicit):
        if windowed and isinstance(r, AngularLift):
            _require(d is not r and rebound(d, r), f"{type(r).__name__} not re-bound on the moment iterate")
        else:
            _require(d is r, f"{type(r).__name__}: the driven gain is not N's term")
    if windowed:
        _require(any(isinstance(r, AngularLift) for r in sp.explicit),
                 "activation: N holds no AngularLift term to re-bind")


def test_g2b_extra_gains_enter_n() -> None:
    """``Splitting.green(extra_gains=(F,))`` puts F in N: it is the driver's last
    gain (``is``) and ``green.inner`` (M − N) applied to a random flux equals
    M·x − Σ Nᵢ·x − F·x. Reds when an extra gain reaches the driver only (the
    pre-S4 path): the promise face's N then omits it."""
    problem = _problem("slab")
    sp = _splitting(problem)
    F = problem.system.factors.fission
    green = sp.green(extra_gains=(F,))
    driver = green.driver(max_iter=1, tol=1e-12, budget_name="max_inner")
    _require(driver.gains[-1] is F, "the extra gain is not the driver's last gain")
    x = _flux(problem, _source(problem, 3), seed=3)
    rhs = sp.implicit.apply(x)
    for g in (*sp.explicit, F):
        rhs = rhs - g.apply(x)
    np.testing.assert_allclose(_flat(green.inner.apply(x)), _flat(rhs), rtol=1e-13, atol=1e-13)


@pytest.mark.catches("ERR-053")
@pytest.mark.parametrize("builder", ["slab", "sphere"])
def test_g2c_the_krylov_green_hands_gmres_the_values_members(builder: str) -> None:
    """The Krylov Green's driver: ``A is splitting.implicit``, its gains are N's
    terms (``is``), ``restart`` is the full ravel ``n_dof`` (ERR-053), and its
    preconditioner equals ``M⁻¹`` bitwise on a random source (no corrector).
    Reds on a clamped restart and on a re-derived gain."""
    from orpheus.sn.solver import _within_group_green_krylov

    problem = _problem(builder)
    sp = _splitting(problem)
    q = _source(problem, 5)
    n_dof = int(_flux(problem, q, seed=5).to_flat().size)
    drv = _within_group_green_krylov(sp, n_dof=n_dof).driver(max_iter=1, tol=1e-12, budget_name="max_inner")
    _require(isinstance(drv, KrylovAcceleration), f"driver is {type(drv).__name__}")
    _require(drv.A is sp.implicit, "GMRES operator is not the value's M")
    _require(tuple(map(id, drv.gains)) == tuple(map(id, sp.explicit)), "GMRES gains re-derived")
    _require(drv.restart == n_dof, f"restart {drv.restart} ≠ n_dof {n_dof} (ERR-053)")
    np.testing.assert_array_equal(_flat(drv.preconditioner.apply(q)), _flat(sp.implicit.inverse().apply(q)))


# ═══════════════════════════════════════════════════════════════════════
# G3 — the splitting map
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("builder", ["slab", "sphere", "cart2d"])
def test_g3a_step_is_one_application_of_the_splitting_map(builder: str) -> None:
    """``green.step(q, ψ)`` equals, bit for bit, ``M⁻¹(q + Σ gᵢψ)`` evaluated by
    hand through the UN-windowed ``M⁻¹`` seeded by ψ, with the Green's driven
    gains, on random q and ψ. On the windowed 2-D Green ψ is a moment iterate
    and the driver's preconditioner is the windowed product, so this row reds
    when ``step`` reuses the preconditioner (a moment composite comes back)."""
    from orpheus.sn.solver import _within_group_green_si

    problem = _problem(builder)
    sp = _splitting(problem)
    green, windowed = _within_group_green_si(sp, problem)
    q = _source(problem, 1)
    psi = _flux(problem, q, 2, windowed=windowed)
    rhs = q
    for g in green.driven_gains:
        rhs = rhs + g.apply(psi)
    hand = sp.implicit.inverse().apply(rhs, initial_guess=psi)
    stepped = green.step(q, psi)
    _require(type(stepped) is type(hand), f"step returned {type(stepped).__name__}, M⁻¹ returns {type(hand).__name__}")
    np.testing.assert_array_equal(_flat(stepped), _flat(hand))


def test_g3b_step_threads_the_iterate_as_the_seed(monkeypatch) -> None:
    """``step`` hands ``M⁻¹`` the iterate as its ``initial_guess``. ``[M]`` the
    26-configuration baseline is value-blind to a dropped seed (0 of 26
    differ), so this spy is the only catcher of a ``step`` calling
    ``apply(rhs)`` without it."""
    from orpheus.sn.operators.sweep_operator import SweepOperator

    problem = _problem("slab")
    green = _splitting(problem).green()
    seen: list[Any] = []
    real = SweepOperator.apply

    def spy(self, x, /, *, initial_guess=None):
        seen.append(initial_guess)
        return real(self, x, initial_guess=initial_guess)

    monkeypatch.setattr(SweepOperator, "apply", spy)
    q = _source(problem, 6)
    psi = _flux(problem, q, 4)
    green.step(q, psi)
    _require(seen and seen[-1] is psi, "step did not thread the iterate as the seed")


@pytest.mark.l1
def test_g3c_step_reproduces_a_converged_iterate() -> None:
    """The fixed point: driven to ``tol = 1e-13``, ``step(q, ψ*)`` is within
    ``1e-10`` of ψ* (relative, max norm); and ``step`` MOVES a random non-fixed
    point, so a ``step`` written as ``return psi`` reds the second leg."""
    problem = _problem("slab")
    green = _splitting(problem).green()
    q = _source(problem, 8)
    cold = 0.0 * _flux(problem, q, 0)
    psi_star, rec = green.driver(max_iter=5000, tol=1e-13, budget_name="max_inner").solve(q, initial_guess=cold)
    _require(rec.converged, "premise: the slab Green did not converge")
    gap = np.max(np.abs(_flat(green.step(q, psi_star)) - _flat(psi_star))) / np.max(np.abs(_flat(psi_star)))
    _require(gap <= 1e-10, f"step moved a converged iterate by {gap:.2e}")
    off = _flux(problem, q, 9)
    _require(not np.allclose(_flat(green.step(q, off)), _flat(off)), "step returns its argument")


# ═══════════════════════════════════════════════════════════════════════
# G2d / G9d — ruling 70 on the S_N realisations
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("schedule", ["jacobi", "gauss_seidel"])
def test_g2d_a_scaled_angular_lift_is_rebound_then_scaled(schedule: str) -> None:
    """qa F5: N = 0.9·S + N₂ₙ + B under the angular-moment window drives
    ``0.9·S.on_moment_domain()`` — the WORD re-bound, then scaled — and its
    action is exactly 0.9 times the unscaled Green's re-bound S on a moment
    iterate. Reds when a scaled term skips the re-binding (the pre-review
    ``c·w`` read as one opaque term: the moment iterate then reaches the
    angular S)."""
    from orpheus.numerics.green_operator import Richardson
    from orpheus.sn.solver import _angular_representation
    from orpheus.transport.operators.angular_lift import AngularLift
    from tests.gates.sn.architecture.test_stage_separation import (
        _is_the_records_gain_rebound_on_the_iterate as rebound,
    )

    problem = _problem("cart2d_p1")
    sp = _splitting(problem, schedule)
    S, *rest = sp.explicit
    _require(isinstance(S, AngularLift), f"premise: N's first term is {type(S).__name__}, not the angular lift S")
    representation = _angular_representation(problem.system.factors.scattering, problem)
    scaled_n = 0.9 * S
    for g in rest:
        scaled_n = scaled_n + g
    green = GreenOperator(sp.implicit, scaled_n, realisation=Richardson(representation=representation))
    unit = GreenOperator(sp.implicit, S, realisation=Richardson(representation=representation))
    gain = green.driven_gains[0]
    terms = getattr(gain, "terms", ())
    _require(len(terms) == 1 and terms[0][0] == 0.9, f"the scaled gain is not one term 0.9·w: {gain!r}")
    _require(rebound(terms[0][1], S), "the scaled lift's word was not re-bound on the moment iterate")
    q = _source(problem, 11)
    x = _flux(problem, q, 12, windowed=True)
    np.testing.assert_array_equal(_flat(gain.apply(x)), _flat(0.9 * unit.driven_gains[0].apply(x)))


def test_g9d_the_windowed_richardson_derives_its_preconditioner_from_m() -> None:
    """Ruling 70 on the S_N window: the Strategy's windowed Green's
    preconditioner is the window applied to the Green's OWN ``implicit_inverse``
    (a ``WindowedSweep`` whose sweep IS that object), never a sweep built
    beside it. Reds when the window wraps a re-derived ``M.inverse()``."""
    from orpheus.sn.operators.windowing import WindowedSweep
    from orpheus.sn.solver import _within_group_green_si

    problem = _problem("cart2d")
    green, windowed = _within_group_green_si(_splitting(problem, "gauss_seidel"), problem)
    _require(windowed, "premise: the 2-D Cartesian Green is not windowed")
    _require(isinstance(green.preconditioner, WindowedSweep), f"{type(green.preconditioner).__name__}")
    _require(green.preconditioner.sweep is green.implicit_inverse,
             "the windowed preconditioner does not wrap the Green's own M⁻¹")


def test_g9e_the_dual_of_a_strategy_green_carries_no_forward_realisation() -> None:
    """The dual split's Green is the default realisation of the DUAL split: a
    forward corrector (DSA, built from the forward low-order system) or a
    forward window is not carried into it. Structural on purpose: a corrector
    vanishes at convergence, so the dual's converged VALUE cannot see it.
    Reds when ``dual`` hands the forward realisation over (qa arm
    dual_keeps_realisation)."""
    from orpheus.numerics.iteration import NATIVE
    from orpheus.sn.acceleration.dsa import DSACorrection
    from orpheus.sn.solver import _within_group_green_si

    problem = _problem("slab")
    green, _ = _within_group_green_si(_splitting(problem), problem, corrector=DSACorrection.from_problem(problem))
    _require(green.realisation.corrector is not None, "premise: the forward Green carries no corrector")
    dual = green.dual()
    _require(dual.realisation is not green.realisation, "the dual kept the forward realisation")
    _require(getattr(dual.realisation, "corrector", None) is None, "the dual carries the forward corrector")
    _require(dual.realisation.representation is NATIVE, "the dual carries a forward representation")


# ═══════════════════════════════════════════════════════════════════════
# G5e — the promise face over the windowed realisation (qa F1, F2)
# ═══════════════════════════════════════════════════════════════════════

#: The promise tolerance of these rows, and the agreement band between
#: realisations: each answer is within tol of the solution in the true
#: residual; ``[M]`` qa 2026-10-21 (``S4_review/q1d_sn_apply.log``) the three
#: realisations agree to ≤ 1.7e-10 at tol 1e-10 (Jacobi and G-S), so the band
#: is ten times that.
_PROMISE_TOL = 1e-10
_AGREE = 2e-9


def _true_residual(sp, x, q) -> float:
    """‖(M − ΣNᵢ)x − q‖/‖q‖ from the splitting's own members (not the Green's
    ``inner``)."""
    r = sp.implicit.apply(x)
    for g in sp.explicit:
        r = r - g.apply(x)
    r = r - q
    return float(np.linalg.norm(_flat(r)) / np.linalg.norm(_flat(q)))


def _three_greens(schedule: str):
    from orpheus.sn.solver import _build_fixed_source_rhs, _within_group_green_krylov, _within_group_green_si

    cfg = _cfg("cart2d_p1")
    problem = _problem("cart2d_p1")
    q = _build_fixed_source_rhs(cfg["external_source"], problem)
    sp = _splitting(problem, schedule)
    windowed_green, windowed = _within_group_green_si(sp, problem)
    _require(windowed, "premise: the 2-D P1 Green is not windowed")
    greens = {
        "windowed": windowed_green,
        "default": sp.green(),
        "krylov": _within_group_green_krylov(sp, n_dof=int(q.to_flat().size)),
    }
    for g in greens.values():
        g.tol, g.max_iter = _PROMISE_TOL, 5000
    return sp, q, greens


@pytest.mark.l1
@pytest.mark.parametrize("schedule", ["jacobi", "gauss_seidel"])
def test_g5e_the_windowed_green_keeps_its_promise_cold(schedule: str) -> None:
    """qa F1: a windowed ``G.apply(q)`` from a cold start returns the
    FULL-ANGULAR solution (the iterate's class is q's, not the moment
    composite), its true residual below tol, and it agrees with the default
    (un-windowed) and the Krylov realisations within :data:`_AGREE`. Reds when
    apply returns the moment iterate, or when the cold start skips the window
    (the native M⁻¹q handed to a moment driver)."""
    sp, q, greens = _three_greens(schedule)
    answers = {}
    for tag, g in greens.items():
        x = g.apply(q)
        _require(_flat(x).size == _flat(q).size, f"[{tag}] apply returned a {type(x).__name__} of the wrong size")
        res = _true_residual(sp, x, q)
        _require(res < _PROMISE_TOL, f"[{tag}] true residual {res:.2e} ≥ tol")
        answers[tag] = _flat(x)
    ref = answers["krylov"]
    for tag in ("windowed", "default"):
        gap = np.max(np.abs(answers[tag] - ref)) / np.max(np.abs(ref))
        _require(gap <= _AGREE, f"[{tag}] differs from Krylov by {gap:.2e}")


@pytest.mark.l1
@pytest.mark.parametrize("schedule", ["jacobi", "gauss_seidel"])
def test_g5f_the_windowed_green_accepts_its_own_answer_as_a_seed(schedule: str) -> None:
    """qa F2: seeded with its own full-angular answer (a state of the split's
    own space), the windowed ``apply`` analyses the seed into the moment
    representation, does not raise, keeps the promise and returns the same
    answer within :data:`_AGREE`. Reds when the seed is handed to the moment
    driver un-analysed (a per-ordinate state where moments belong)."""
    sp, q, greens = _three_greens(schedule)
    g = greens["windowed"]
    first = g.apply(q)
    again = g.apply(q, initial_guess=first)
    res = _true_residual(sp, again, q)
    _require(res < _PROMISE_TOL, f"seeded true residual {res:.2e} ≥ tol")
    gap = np.max(np.abs(_flat(again) - _flat(first))) / np.max(np.abs(_flat(first)))
    _require(gap <= _AGREE, f"the seeded answer moved by {gap:.2e}")
