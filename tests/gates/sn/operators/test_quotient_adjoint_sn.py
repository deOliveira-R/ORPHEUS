r"""S6b on the S\ :sub:`N` carriers — degenerate traces, ruling 80, the adjoint modulo the kernel.

Ruling 80: every production S\ :sub:`N` operator with a degenerate end maps
``ker G_domain`` into ``ker G_codomain`` (``P_W A (I − P_V) = 0``), so today's
``.H`` is already the adjoint on the quotient ``V / ker G`` and no number moves.
The rows here pin that property leaf by leaf on the census's degenerate cases
(``scratch/numerics_untangling/packages/S6bc_census/report.md`` §2): the
cylinder on ``folded_product(2, 6)`` (μ_x = 0 bit-exactly, reflective and vacuum),
the sphere on GL-7 and the slab on GL-7 (a μ = 0 node). They also pin the trace
ranks against hand counts read from the quadrature, and the sweep's adjoint laws
modulo the kernel. Controls: the cylinder on ``folded_product(4, 8)`` and the
sphere on GL-8 are full rank.

``-O``-safe: every assertion is ``_require`` (a raise) or ``np.testing``.
"""

from __future__ import annotations

import collections

import numpy as np
import pytest

from orpheus.data.macro_xs.synthetic import get_mixture
from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.coupled_system import CoupledField, CoupledOperator
from orpheus.numerics.operator import AdjointOperator, OperatorProduct, OperatorSum, adjointable
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace
from orpheus.geometry.deck import MirrorDeck
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates._exterior import conditions_for

#: Exterior conditions of the tags this file's meshes carry (wave 2: a face's role is a deck or an int tag; the condition a tag
#: stands for is declared at the solve, keyword ``boundary_conditions``).
_CONDITIONS = {0: WallCondition.specular(1.0), 1: WallCondition.vacuum()}


pytestmark = [
    pytest.mark.foundation,
    pytest.mark.rests_on(
        "tests/gates/numerics/test_quotient_metric.py::test_B3_project_live_is_the_g_orthogonal_projector_onto_range_g",
        "tests/gates/numerics/test_quotient_metric.py::test_B1_metric_rank_is_the_hand_count",
    ),
]

_EPS = float(np.finfo(float).eps)


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _mesh(geom):
    return Mesher(geom).partition(CellsByCount.uniform_width(4)).mesh


#: the census's cases: (geometry, quadrature, number of faces carrying a trace)
_CASES = {
    "cyl_refl_fp26": (lambda: StructuredGeometry.cylinder((0.0, 1.0), (0,), outer=0),
                      lambda: Quadrature.folded_product(n_mu=2, n_phi=6), 1),
    "cyl_vac_fp26": (lambda: StructuredGeometry.cylinder((0.0, 1.0), (0,), outer=1),
                     lambda: Quadrature.folded_product(n_mu=2, n_phi=6), 1),
    "sph_vac_gl7": (lambda: StructuredGeometry.sphere((0.0, 3.0), (0,), outer=1),
                    lambda: Quadrature.gauss_legendre(n_ordinates=7), 1),
    "slab_refl_gl7": (lambda: StructuredGeometry.slab((0.0, 2.0), (0,), left=MirrorDeck(), right=1),
                      lambda: Quadrature.gauss_legendre(n_ordinates=7), 2),
    # full-rank controls
    "cyl_refl_fp48": (lambda: StructuredGeometry.cylinder((0.0, 1.0), (0,), outer=0),
                      lambda: Quadrature.folded_product(n_mu=4, n_phi=8), 1),
    "sph_vac_gl8": (lambda: StructuredGeometry.sphere((0.0, 3.0), (0,), outer=1),
                    lambda: Quadrature.gauss_legendre(n_ordinates=8), 1),
}
_DEGENERATE = ["cyl_refl_fp26", "cyl_vac_fp26", "sph_vac_gl7", "slab_refl_gl7"]


def _problem(name):
    from orpheus.sn.solver import _adjoint_posing_parts, _as_problem

    geom, rule, n_faces = _CASES[name]
    quad = rule()
    mesh = _mesh(geom())
    problem = _as_problem(mesh, quad, {0: get_mixture("A", "2g")}, boundary_conditions=conditions_for(mesh, _CONDITIONS))
    M, N, F, template, _ = _adjoint_posing_parts(problem)
    return problem, quad, n_faces, M, N, F, template


def _random_like(tmpl, seed):
    n = np.asarray(tmpl.to_flat()).size
    return type(tmpl).from_flat(np.random.default_rng(seed).standard_normal(n), tmpl)


def _flat(v):
    return np.asarray(v.to_flat(), dtype=float)


def _ends(op):
    V, W = getattr(op, "domain", None), getattr(op, "codomain", None)
    if not isinstance(V, FunctionSpace) or not isinstance(W, FunctionSpace):
        pytest.fail(f"fixture: {type(op).__name__} has an unbound end")
    return V, W


def _gnorm(space, v):
    return float(np.sqrt(max(space.inner_product(v, v), 0.0)))


# ═══════════════════════════════════════════════════════════════════════
# T1 — the trace ranks, against the quadrature's tangential count
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("name", list(_CASES))
def test_T1_the_trace_rank_is_the_non_tangential_count(name):
    """``trace.rank = size − ng · faces · #{n : μ_x,n = 0}``, the hand count read
    from the quadrature's direction cosines (independent of the trace metric):
    cylinder fp(2,6) 8 of 12, sphere GL-7 12 of 14, slab GL-7 24 of 28 (two faces);
    the controls full rank. The full field's rank is interior + trace, the coupled
    carrier's the sum over its systems, and each dual's rank its primal's. First
    red: a rank counting the element size."""
    problem, quad, n_faces, M, *_ = _problem(name)
    ff = problem.full_field_space
    tr = ff.trace_space if ff.trace_space is not None else pytest.fail("fixture: no trace block")
    interior = ff.interior_space if ff.interior_space is not None else pytest.fail("fixture: no interior block")
    ng = problem.ng
    tangential = int(np.sum(np.asarray(quad.mu_x, dtype=float) == 0.0))
    size = int(np.prod(tr.shape))
    _require(size == ng * n_faces * quad.n_ordinates, f"[{name}] trace size {size} is not ng·faces·N (fixture)")
    _require(tr.rank == size - ng * n_faces * tangential,
             f"[{name}] trace rank {tr.rank} vs hand {size - ng * n_faces * tangential}")
    if name in _DEGENERATE:
        _require(tangential > 0, f"[{name}] activation: no tangential ordinate")
    _require(ff.rank == interior.rank + tr.rank, f"[{name}] full-field rank is not the block sum")
    V = M.domain
    if not isinstance(V, FunctionSpace):
        pytest.fail(f"[{name}] fixture: M has no domain space")
    _require(V.rank == sum(s.rank for s in getattr(V, "systems", ())), f"[{name}] the coupled rank is not the system sum")
    _require(V.dual().rank == V.rank and tr.dual().rank == tr.rank, f"[{name}] a dual's rank differs")


# ═══════════════════════════════════════════════════════════════════════
# K80 — ruling 80: every degenerate-ended production leaf maps ker into ker
# ═══════════════════════════════════════════════════════════════════════


def _walk(op, depth=0):
    yield op
    if depth > 6:
        return
    if isinstance(op, CoupledOperator):
        for row in op.blocks:
            for b in row:
                if b is not None:
                    yield from _walk(b, depth + 1)
    elif isinstance(op, OperatorSum):
        for _, w in op.terms:
            yield from _walk(w, depth + 1)
    elif isinstance(op, OperatorProduct):
        for f in op.factors:
            yield from _walk(f, depth + 1)


#: leaf classes the census measured (report §2); the walk must reach each.
_REQUIRED = {"StreamingCollisionOperator", "ScatteringOperator", "N2NOperator", "FissionOperator",
             "SNBoundaryOperator", "CoupledSubstitutionOperator"}
_REQUIRED_CURVILINEAR = {"RadialCharacteristicSeeding", "RadialCharacteristicEmission", "SystemRestrictionOperator"}


@pytest.mark.parametrize("name", _DEGENERATE)
def test_K80_every_degenerate_ended_leaf_maps_the_kernel_into_the_kernel(name):
    """For every operator reachable from the loss, the Strategy pair ``M``, ``N``,
    the fission and the sweep ``M⁻¹`` (blocks, sum terms and word factors walked)
    whose domain or codomain is degenerate: ``‖P_W A (x − P_V x)‖ ≤ 1e-13 ‖A x‖``
    on a random ``x`` of the class the operator admits. `[M]` 2026-10-21 the worst
    reading per case ≤ 1.3e-16. Coverage: the census's leaf classes (and on the
    curvilinear cases the radial-characteristic blocks and the system restriction)
    are each reached at least once; operators with no admissible sample class are
    counted and must stay ≤ 1 (one coupled gain grid today). First red: the
    census's mutated boundary (a tangential slot feeding a live one, 0.70 at the
    leaf; battery arm ``boundary_leak``)."""
    problem, quad, _, M, N, F, t = _problem(name)
    src = M.apply(t)
    cands = [t, src, *t.systems, *src.systems]
    tops = [problem.pencil.lhs, M, N, F, M.inverse()]
    seen, reached, no_input, rows = set(), collections.Counter(), 0, []
    for top in tops:
        for op in _walk(top):
            if id(op) in seen:
                continue
            seen.add(id(op))
            V, W = getattr(op, "domain", None), getattr(op, "codomain", None)
            if not isinstance(V, FunctionSpace) or not isinstance(W, FunctionSpace):
                continue
            if V.rank == int(np.prod(V.shape)) and W.rank == int(np.prod(W.shape)):
                continue
            pair = None
            for c in cands:
                try:
                    x = _random_like(c, 7)
                    pair = (x, op.apply(x))
                    break
                except (TypeError, ValueError, AttributeError):
                    continue
            if pair is None:
                no_input += 1
                continue
            x, y = pair
            leak = W.project_live(op.apply(x - V.project_live(x)))
            r = float(np.linalg.norm(_flat(leak)) / max(np.linalg.norm(_flat(y)), 1e-300))
            reached[type(op).__name__] += 1
            rows.append((type(op).__name__, r))
    bad = [(n, r) for n, r in rows if r > 1e-13]
    _require(bad == [], f"[{name}] ker → live leaks: {bad}")
    required = _REQUIRED | (_REQUIRED_CURVILINEAR if not name.startswith("slab") else set())
    missing = required - set(reached)
    _require(not missing, f"[{name}] the walk did not reach {sorted(missing)} (reached {dict(reached)})")
    _require(no_input <= 1, f"[{name}] {no_input} degenerate-ended operators had no admissible sample")


# ═══════════════════════════════════════════════════════════════════════
# A — the sweep's adjoint laws modulo the kernel
# ═══════════════════════════════════════════════════════════════════════


def _source_sample(M, t, seed):
    r"""A source-classed sample WITH a kernel component (the sweep's admitted
    class): random coordinates on the class of ``M t``. ``M x`` alone has none
    (``M`` annihilates the tangential slots), which would make the rows blind."""
    tmpl = M.apply(t)
    n = np.asarray(tmpl.to_flat()).size
    return CoupledField.from_flat(np.random.default_rng(seed).standard_normal(n), tmpl)


@pytest.mark.parametrize("name", _DEGENERATE)
def test_A1_the_double_adjoint_of_the_sweep_is_the_sweep_modulo_the_kernel(name):
    """``A.H.H is A`` is the representative choice on ``V / ker G``; computed
    out, the double adjoint ``AdjointOperator(A.H)`` agrees with ``A`` after
    ``P_W`` (seminorm ≤ 1e-14 relative) and differs in the full norm, because the
    sweep passes tangential slots through as the identity. `[M]` 2026-10-21 on a
    random sample: full 0.269 (cylinder), 0.086 (sphere GL-7), 0.724 (slab GL-7);
    the census's 0.519 on the cylinder is the matrix measure of the same ker→ker
    block. Both legs asserted: the full-norm gap is the activation. First red: a
    sweep transpose that writes kernel input into live slots."""
    problem, quad, _, M, N, F, t = _problem(name)
    S = M.inverse()
    _, W = _ends(S)
    q = _source_sample(M, t, 5)
    a = S.apply(q)
    SH = S.H
    if not adjointable(SH):
        pytest.fail(f"[{name}] S.H is not adjointable")
    hh = AdjointOperator(SH).apply(q)  # the computed double adjoint (S.H.H is S by choice)
    d = hh - a
    full = float(np.linalg.norm(_flat(d)) / np.linalg.norm(_flat(a)))
    semi = _gnorm(W, d) / _gnorm(W, a)
    _require(full > 1e-3, f"[{name}] activation: the double adjoint equals the sweep in the full norm ({full:.1e})")
    _require(semi <= 1e-14, f"[{name}] A.H.H and A differ on the quotient: seminorm {semi:.2e}")


@pytest.mark.parametrize("name", _DEGENERATE)
def test_A2_sweep_reciprocity_on_degenerate_ends(name):
    """``⟨S q, y⟩_W = ⟨q, S.H y⟩_V`` with a kernel component in both samples,
    1e-13 relative (the pairing is the G-seminorm's form: it never reads the
    kernel). First red: a transpose sweep dropping a live completion (vv #37's
    one-sided partner is A3)."""
    problem, quad, _, M, N, F, t = _problem(name)
    S = M.inverse()
    q = _source_sample(M, t, 11)
    y = _source_sample(M, t, 12)  # [M] S.H admits the source class (it re-enters M's transpose)
    V, W = _ends(S)
    lhs = W.inner_product(S.apply(q), y)
    rhs = V.inner_product(q, S.H.apply(y))
    _require(abs(lhs - rhs) <= 1e-13 * max(1.0, abs(lhs)), f"[{name}] reciprocity {lhs!r} vs {rhs!r}")


@pytest.mark.parametrize("name", _DEGENERATE)
def test_A3_the_implicit_part_inverts_the_sweep_on_the_quotient(name):
    """``M (M⁻¹ q) = P q``: the sweep is an inverse on the quotient only
    (census §2: ``M·M⁻¹ e_k = 0`` on every dead slot). ``‖M M⁻¹ q − P q‖ ≤
    1e-14 ‖q‖`` on a random q, and on each dead coordinate ``M M⁻¹ e_k = 0``
    (``≤ 1e-14``); the full-norm gap ``‖M M⁻¹ q − q‖`` is the activation. First
    red: a sweep whose live block is not ``M⁻¹`` (the one-sided identity gate)."""
    problem, quad, _, M, N, F, t = _problem(name)
    S = M.inverse()
    q = _source_sample(M, t, 21)
    r = M.apply(S.apply(q))
    P = _ends(M)[1].project_live
    nq = float(np.linalg.norm(_flat(q)))
    _require(np.linalg.norm(_flat(r) - _flat(P(q))) <= 1e-14 * nq, f"[{name}] M M⁻¹ q != P q")
    _require(np.linalg.norm(_flat(r) - _flat(q)) > 1e-3 * nq, f"[{name}] activation: M M⁻¹ = I in the full norm")
    flat_q = _flat(q)
    dead = np.flatnonzero(_flat(P(q)) == 0.0)
    _require(dead.size > 0, f"[{name}] activation: no dead coordinate")
    for k in dead:
        e = np.zeros_like(flat_q)
        e[k] = 1.0
        ek = CoupledField.from_flat(e, q)
        _require(np.abs(_flat(M.apply(S.apply(ek)))).max() <= 1e-14, f"[{name}] M M⁻¹ e_{k} != 0")
