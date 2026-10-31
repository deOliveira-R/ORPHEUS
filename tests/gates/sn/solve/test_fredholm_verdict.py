r"""U0 rung 5 — every SN entry under the Fredholm verdict.

Spec ``scratch/membership_covector/u0/spec_u0.md`` §1 rung 5 (S1-S7, S9-S15; S8 and S16a/b are in ``tests/gates/sn/operators/test_predicted_loss_kernel.py``). Entries: eigen on ``solve_sn``,
``solve_sn_adjoint`` and the multiplying source's k-solve; source on ``solve_sn_fixed_source``,
``solve_sn_adjoint_fixed_source`` and ``solve_sn_multiplying_source`` (q on at(1)) — the orchestrator's list, 2026-10-31.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from orpheus.data.macro_xs.synthetic import make_mixture
from orpheus.mesh import AxisMesh
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import (
    solve_sn, solve_sn_adjoint, solve_sn_adjoint_fixed_source, solve_sn_fixed_source, solve_sn_multiplying_source,
)
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates._exterior import conditions_for
from tests.gates.sn.operators import test_loss_kernel_gauge as T

from tests.gates.transport._singular_loss_fixtures import (
    ABSORBER, FIS_FEEDS, FIS_ISOLATED, FIS_ISOLATED_KINF, FIS_REGULAR, PURE_1G, bordered_member, dense_of, null, slab,
)

from orpheus.numerics.fredholm import Gauge, NoSolution, Solve
from orpheus.transport import absorbing_chain as _ac  # noqa: F401  (the predicate the rows rest on)

_Refusal = NoSolution


def _verdict_of(solution):
    return solution.exit_report.solvability


_FRAGMENT = "no steady solution"
GL8 = Quadrature.gauss_legendre(n_ordinates=8)
pytestmark = [pytest.mark.l1, pytest.mark.rests_on(
    "tests/gates/transport/test_absorbing_chain_predictor.py",
    "tests/gates/numerics/test_fredholm_verdict.py",
    "tests/gates/sn/operators/test_predicted_loss_kernel.py",
)]


def _quiet(f):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return f()


# ── fixed source ─────────────────────────────────────────────────────────────────────────────────────────

@pytest.mark.catches("ERR-123")
@pytest.mark.parametrize("inner", [None, "krylov"])
def test_S1_a_positive_source_on_a_closed_pure_scatterer_is_refused(inner, monkeypatch):
    """[M] today 727.8 (SI) / 1.71e16 after 67.5 s (Krylov), admissibility 'always admitted: A is invertible'. The refusal
    precedes every sweep: a spy on the two private arms counts 0 calls."""
    import orpheus.sn.solver as S
    calls = []
    for arm in ("_solve_fixed_source_si", "_solve_fixed_source_krylov"):
        honest = getattr(S, arm)
        monkeypatch.setattr(S, arm, lambda *a, _h=honest, **k: (calls.append(1), _h(*a, **k))[1])
    with pytest.raises(_Refusal, match=_FRAGMENT):
        solve_sn_fixed_source({0: PURE_1G}, slab(), GL8, np.ones((8, 1, 8)), max_inner=200, inner_solver=inner)
    assert calls == []


def test_S2_a_closed_wall_carrying_inflow_is_refused():
    """The left co-vector's TRACE leg: zero bulk source, white(1.0) + uniform(1.0) inflow; the rhs is trace-only and pairs
    at 0.092 [M]. Today: max phi 20.5, max_inner exhausted, 'always admitted'."""
    wall = WallCondition.white(1.0) + WallCondition.uniform(1.0)
    assert isinstance(wall, WallCondition)
    with pytest.raises(_Refusal, match=_FRAGMENT):
        solve_sn_fixed_source({0: PURE_1G}, slab(None, 0), GL8, np.zeros((8, 1, 8)), max_inner=200,
                              boundary_conditions={0: wall})


@pytest.mark.parametrize("label,kw", [
    ("absorber", dict(materials={0: ABSORBER}, mesh=slab(), source=np.ones((8, 1, 8)), bcs={})),
    ("white0.5_inflow", dict(materials={0: PURE_1G}, mesh=slab(None, 0), source=np.zeros((8, 1, 8)),
                             bcs={0: WallCondition.white(0.5) + WallCondition.uniform(1.0)})),
])
def test_S3_an_admitted_posing_states_its_verdict(label, kw):
    """[M] today 'the pure-transport question is always admitted: A is invertible' on both (20.0 and 3.90)."""
    sol = _quiet(lambda: solve_sn_fixed_source(kw["materials"], kw["mesh"], GL8, kw["source"], boundary_conditions=kw["bcs"]))
    assert isinstance(_verdict_of(sol), Solve) and "absorption" in _verdict_of(sol).because
    assert "always admitted" not in str(sol.exit_report.admissibility)


def _one_d(c=1.0, cells=8, quad=GL8, stretch=1.0):
    axes = (AxisMesh(edges=T._graded(1.0, cells, stretch), role_low=T._R, role_high=T._R),)
    return axes, T._mesh((cells,), [(T._R, T._R)], c=c, ng=1, quad=quad, stretch=stretch)


@pytest.mark.parametrize("inner", [None, "krylov"])
def test_S4_a_balanced_source_is_gauged_to_the_minimum_G_norm_member(inner):
    """REFERENCE: dense bordered-G member on the production rhs. [M] today both arms return a member off it by a CONSTANT
    kernel component 2.39e-3 (a valid solution, not the gauged one). Arm ``trace_only_projector`` (today's
    _exit_gauge_trace scope) reds: the physical kernel is 1 on the bulk AND the trace."""
    from orpheus.sn.solver import _build_fixed_source_rhs
    axes, problem = _one_d()
    system, tmpl = T._loss_system(problem)
    n, nb = tmpl.to_flat().size, tmpl.interior.values.size
    A = dense_of(lambda e: T._apply_loss(system, tmpl, e), n)
    _, R = null(A)
    sv = np.linalg.svd(A, compute_uv=False)
    cond_regular = sv[0] / sv[-2]            # the condition of A on the complement of its one-dimensional kernel
    bal = np.linspace(-1.0, 1.0, 8) + 0.3 * np.sin(np.arange(8))
    bal -= bal.mean()
    src = np.broadcast_to(bal[None, None, :], (8, 1, 8)).copy()
    q = np.asarray(_build_fixed_source_rhs(src, problem).to_flat(), float)
    G = np.asarray(problem.full_field_space.metric_diagonal(), dtype=float)
    x_ref = bordered_member(A, q, R, (G[:, None] * R).T)
    sol = _quiet(lambda: solve_sn_fixed_source({0: T._mixture(1, 1.0)}, axes, GL8, src, max_inner=400, inner_tol=1e-13,
                                               inner_solver=inner, boundary_conditions=conditions_for(axes, T._CONDITIONS)))
    x = np.asarray(sol.angular_flux.to_flat())
    np.testing.assert_allclose(x, x_ref, rtol=0, atol=10 * 1e-13 * cond_regular * np.abs(x_ref).max())


# ── the union kernel on the d = 2 DD box ─────────────────────────────────────────────────────────────────

def _box(c, quad=None, ng=1):
    axes = tuple(AxisMesh(edges=T._graded(e, n), role_low=T._R, role_high=T._R) for e, n in ((1.0, 3), (2.0, 4)))
    return axes, quad or Quadrature.level_symmetric(sn_order=2)


def test_S5_the_DD_kernel_with_an_ordinary_source_is_gauged_as_today():
    """[M] c = 0.9: nullity 2, left vectors trace-only (mass 1.000), uniform source pairs 2.4e-15 -> Gauge, dim 2."""
    axes, quad = _box(0.9)
    n_ord = int(np.asarray(quad.weights).size)
    sol = _quiet(lambda: solve_sn_fixed_source({0: T._mixture(1, 0.9)}, axes, quad, np.ones((n_ord, 1, 3, 4)),
                                               boundary_conditions=conditions_for(axes, T._CONDITIONS)))
    assert isinstance(_verdict_of(sol), Gauge) and sol.outcome.gauge.dimension == 2


def test_S6_the_union_kernel_refuses_a_uniform_source():
    """[M] c = 1: nullity 3 = 2 + 1; uniform source pairs 0.522. Today max phi 2.90e3, max_inner exhausted, trace gauged."""
    axes, quad = _box(1.0)
    n_ord = int(np.asarray(quad.weights).size)
    with pytest.raises(_Refusal, match=_FRAGMENT):
        solve_sn_fixed_source({0: T._mixture(1, 1.0)}, axes, quad, np.ones((n_ord, 1, 3, 4)),
                              boundary_conditions=conditions_for(axes, T._CONDITIONS))


def test_S7_the_union_kernel_gauges_a_source_balanced_against_its_bulk_covector():
    """[M] L[:nb] has singular values (0.522, 1.9e-14, 1.5e-15): ONE bulk direction (the physical co-vector); a random
    bulk source with that direction removed pairs 1.2e-15 -> Gauge on DD (+) physical (G-orthogonal on the trace)."""
    axes, quad = _box(1.0)
    problem = T._mesh((3, 4), [(T._R, T._R)] * 2, c=1.0, ng=1, quad=quad)
    system, tmpl = T._loss_system(problem)
    n, nb = tmpl.to_flat().size, tmpl.interior.values.size
    L, _ = null(dense_of(lambda e: T._apply_loss(system, tmpl, e), n))
    u = np.linalg.svd(L[:nb], full_matrices=False)[0][:, 0]
    qb = np.random.default_rng(1).random(nb)
    qb -= u * (u @ qb)
    n_ord = int(np.asarray(quad.weights).size)
    sol = _quiet(lambda: solve_sn_fixed_source({0: T._mixture(1, 1.0)}, axes, quad, qb.reshape(n_ord, 1, 3, 4),
                                               boundary_conditions=conditions_for(axes, T._CONDITIONS)))
    assert isinstance(_verdict_of(sol), Gauge) and sol.outcome.gauge.dimension == 3


# ── eigen entries (q = the fission emission) ─────────────────────────────────────────────────────────────

@pytest.mark.catches("ERR-123")
@pytest.mark.verifies("fredholm-eigen-whole-pencil")
def test_S9_solve_sn_refuses_emission_into_the_class():
    """[M] today k = 1.25, max_inner exhausted, 'an eigen question has no admissibility to check'."""
    with pytest.raises(_Refusal, match=_FRAGMENT):
        solve_sn({0: FIS_FEEDS}, slab(), GL8, max_outer=60, max_inner=200)


@pytest.mark.verifies("fredholm-eigen-whole-pencil")
def test_S10_solve_sn_gauges_an_isolated_class():
    """[M] today k = 1.25 (the gauged value), no verdict recorded."""
    sol = _quiet(lambda: solve_sn({0: FIS_ISOLATED}, slab(), GL8))
    assert sol.outcome.lam == pytest.approx(FIS_ISOLATED_KINF, rel=1e-6)
    assert isinstance(_verdict_of(sol), Gauge)


#: An absorption-free class that CARRIES fission: the (n,2n) gain offsets the fission and the scatter loss exactly
#: (dyadic: removal 0.5 - 0.375 - 2 x 0.0625 = 0.0 in floats). Both eigen posings are refused on it.
FISSILE_N2N_BALANCED = make_mixture(sig_t=np.array([0.5]), sig_c=np.array([0.0]), sig_f=np.array([0.0625]),
                                    nu=np.array([2.5]), chi=np.array([1.0]), sig_s=np.array([[0.375]]),
                                    sig_2=np.array([[0.0625]]))


@pytest.mark.verifies("fredholm-eigen-whole-pencil")
def test_S11_solve_sn_adjoint_refuses_fission_carried_by_the_class():
    """The daggered posing's annihilators are flat(r), r in ker A, and the range of F^H is the production nu Sf:
    refused iff an absorption-free class carries fission. [M] FISSILE_N2N_BALANCED: NoSolution on both entries."""
    with pytest.raises(_Refusal, match=_FRAGMENT):
        solve_sn_adjoint({0: FISSILE_N2N_BALANCED}, slab(), GL8, max_outer=60, max_inner=200)


@pytest.mark.verifies("fredholm-eigen-whole-pencil")
def test_S11b_the_adjoint_refuses_exactly_where_the_forward_does():
    """The eigen verdict is the PENCIL's (the user's ruling, 2026-10-31: the posing derives it): the kernel must be the
    whole pencil's, F r = 0 on the modes AND F^T w = 0 on the annihilators, a condition symmetric under ``.H``. So
    FIS_FEEDS (fission emits into the absorption-free class {2}) is refused by BOTH eigen entries.

    [REFUTED 2026-10-31] the earlier row "the adjoint GAUGES a class that only receives fission neutrons": the
    daggered annihilators flat(r) pair to zero with the range of F^H whenever F r = 0, so that reading could never
    refuse, and gauging was not residual-neutral for A^H psi = mu F^H psi. [M] qa (scratch/membership_covector/u0/
    review/qa.md, F1): the adjoint then returned lambda = 0.9999999999999987 at every max_outer with a Rayleigh gap of
    0.769, against 0.0 with the gauge removed. Arm ``adjoint_reads_its_own_range`` (the daggered verdict paired with
    F^H only) reddens this row."""
    for entry in (solve_sn, solve_sn_adjoint):
        with pytest.raises(_Refusal, match=_FRAGMENT):
            entry({0: FIS_FEEDS}, slab(), GL8, max_outer=60, max_inner=200)


def test_S12_the_multiplying_source_k_solve_refuses_before_the_subcriticality_check():
    """DISCRIMINATION. [M] today SupercriticalSourceProblem after 288.9 s (the k-solve runs to its budget). Post-carve
    the eigen verdict fires first, so the row is fast."""
    with pytest.raises(_Refusal, match=_FRAGMENT) as info:
        solve_sn_multiplying_source({0: FIS_FEEDS}, slab(), GL8, np.ones((8, 2, 8)), max_inner=200)
    assert "SUBCRITICAL" not in str(info.value)


def test_S13_the_multiplying_source_is_solved_when_subcritical():
    """[M] today Asserted(k = 0.625) and max phi 233.2; kept beside the Solve verdict."""
    sol = _quiet(lambda: solve_sn_multiplying_source({0: FIS_REGULAR}, slab(), GL8, np.ones((8, 2, 8)), max_inner=200))
    assert isinstance(_verdict_of(sol), Solve)
    assert type(sol.exit_report.admissibility).__name__ == "Asserted"   # the subcriticality predicate, unchanged


# ── the adjoint fixed source: the LEFT projector ─────────────────────────────────────────────────────────

def test_S14_the_adjoint_fixed_source_refuses_a_positive_detector():
    """[M] today 727.8 silent, 'the pure-transport adjoint question has no admissibility to check'."""
    with pytest.raises(_Refusal, match=_FRAGMENT):
        solve_sn_adjoint_fixed_source({0: PURE_1G}, slab(), GL8, np.ones((1, 8)), max_inner=200)


@pytest.mark.parametrize("cells,quad,stretch", [(8, GL8, 1.0), (3, Quadrature.gauss_legendre(n_ordinates=4), 2.0)],
                         ids=["GL8_uniform", "GL4_graded"])
def test_S15_the_adjoint_member_is_annihilated_by_the_left_covectors(cells, quad, stretch):
    """THEOREM on the returned member x of A^H x = q (A^H = G^-1 A^T G): (i) the residual is at the solver's tolerance;
    (ii) w(x) = 0 canonically for every left co-vector w of A, i.e. x is G-orthogonal to ker A^H = sharp(w): the
    minimum-G-norm member of the ADJOINT coset. [M] cos(sharp(w), r) = 0.80 / 0.60, so the forward projector (today's
    sn/solver.py:2920 wiring; arm ``forward_projector``) gives w(x) = 7.7e-2 / 5.5e-2."""
    axes, problem = _one_d(cells=cells, quad=quad, stretch=stretch)
    system, tmpl = T._loss_system(problem)
    n = tmpl.to_flat().size
    A = dense_of(lambda e: T._apply_loss(system, tmpl, e), n)
    L, R = null(A)
    G = problem.full_field_space.metric_diagonal()
    V = np.diff(np.asarray(axes[0].edges))
    bal = np.random.default_rng(3).random(cells)
    bal -= (V @ bal) / V.sum()              # flat(r)(q) = 0: r is the constant, so the detector integrates to zero
    sv = np.linalg.svd(A, compute_uv=False)
    cond_regular = sv[0] / sv[-2]
    sol = _quiet(lambda: solve_sn_adjoint_fixed_source({0: T._mixture(1, 1.0)}, axes, quad, bal[None, :], inner_tol=1e-13,
                                                       max_inner=2000, boundary_conditions=conditions_for(axes, T._CONDITIONS)))
    x = np.asarray(sol.angular_flux.to_flat())
    w = L[:, 0]
    assert abs(w @ x) <= 10 * 1e-13 * cond_regular * np.linalg.norm(w) * np.linalg.norm(x)
