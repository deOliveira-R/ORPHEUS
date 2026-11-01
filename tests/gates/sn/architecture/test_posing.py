"""AC-a — no Strategy token reaches the Problem chain (consumers campaign step 2, C3b).

The chain ``SNProblem(...) → .system → .pencil → .eigen_posing`` is the PROBLEM's:
its members are determined by the generating data alone, so no callable on it
may accept a Strategy token (``_STRATEGY_TOKENS`` — the solver-choice vocabulary
of ``campaign_verification_plan.md`` §AC-a).  Non-tautological only since C3b
(pre-C3b the chain did not exist); the mutation that reds it: add
``inner_schedule: str = "jacobi"`` to ``build_within_group_system``.
"""
from __future__ import annotations

import inspect

import pytest

from orpheus.sn.coupled_system import build_within_group_system
from orpheus.sn.problem import SNProblem
import numpy as np
from orpheus.data.macro_xs.synthetic import get_mixture
from orpheus.geometry import StructuredGeometry
from orpheus.geometry.deck import MirrorDeck
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import _as_problem, _build_fixed_source_rhs
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates._exterior import conditions_for

pytestmark = pytest.mark.foundation

_STRATEGY_TOKENS = frozenset({
    "inner_solver", "inner_schedule", "max_iter", "max_inner", "tol", "inner_tol",
    "restart", "corrector", "preconditioner", "n_dof", "initial_guess",
})


def _require(cond: object, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


def _chain() -> list[tuple[str, object]]:
    """The callables on the Problem chain — a LIST, so a renamed member cannot
    silently empty the loop."""
    members = [
        ("SNProblem.__init__", SNProblem.__init__),
        ("SNProblem.from_axes", SNProblem.from_axes),
        ("SNProblem.from_material_mesh", SNProblem.from_material_mesh),
        ("SNProblem.with_scattering_order", SNProblem.with_scattering_order),
        ("SNProblem.with_cross_sections", SNProblem.with_cross_sections),
        ("build_within_group_system", build_within_group_system),
        ("SNProblem.system", SNProblem.__dict__["system"].func),
        ("SNProblem.pencil", SNProblem.__dict__["pencil"].func),
        ("SNProblem.eigen_posing", SNProblem.__dict__["eigen_posing"].func),
    ]
    return members


def test_ac_a_no_strategy_token_reaches_the_problem_chain() -> None:
    chain = _chain()
    _require(len(chain) >= 3, "the chain must be non-empty (a renamed member would empty the loop)")
    offenders = []
    for name, fn in chain:
        params = set(inspect.signature(fn).parameters)
        hit = params & _STRATEGY_TOKENS
        if hit:
            offenders.append((name, sorted(hit)))
    _require(not offenders, f"Strategy tokens on the Problem chain: {offenders}")


def test_ac_a_the_chain_is_determined_by_the_generating_data() -> None:
    """Two hubs built from equal data pose EQUAL records — the chain reads no
    solver state (the identity half of AC-a)."""
    from tests.gates.sn.operators.test_step2_posed_fission_anchors import _slab_hub

    a, b = _slab_hub(), _slab_hub()
    _require(a == b, "content-equal hubs")
    _require(a.system is a.system and a.pencil is a.pencil, "the members are cached per hub")
    _require(a.pencil.lhs is a.system.loss and a.pencil.rhs is a.system.production, "the pencil is over the hub's own record")
    _require(a.eigen_posing.pencil is a.pencil, "the posing is over the hub's own pencil")


# ── moved here by the layer DAG Phase 2 refresh (specs_2r §5) ──


_CONDITIONS = {0: WallCondition.vacuum()}


def _slab(length: float = 4.0):
    materials = {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
    geometry = StructuredGeometry.slab(
        (0.0, length / 2, length), (0, 1), left=MirrorDeck(), right=0,
    )
    mesh = Mesher(geometry).partition(CellsByCount.uniform_width(4)).mesh
    return _as_problem(mesh, Quadrature.gauss_legendre(8), materials, boundary_conditions=conditions_for(mesh, _CONDITIONS))


def _carrying_sphere():
    geometry = StructuredGeometry.sphere((0.0, 3.0), (0,), outer=0)
    mesh = Mesher(geometry).partition(CellsByCount.uniform_width(6)).mesh
    return _as_problem(mesh, Quadrature.gauss_legendre(8), {0: get_mixture("A", "2g")}, boundary_conditions=conditions_for(mesh, _CONDITIONS))


@pytest.mark.rests_on("tests/gates/numerics/test_outcome.py::TestLawTheSourceLivesOnTheOperatorsEnds::test_negative_a_foreign_source_is_REFUSED_at_the_posing")
class TestLawTheCoupledStateKnowsItsSpace:
    @pytest.mark.parametrize("hub_of, nx", [(_slab, 8), (_carrying_sphere, 6)], ids=["slab (1 system)", "carrying sphere (2 systems)"])
    def test_law_the_source_lives_on_the_pencils_domain(self, hub_of, nx: int) -> None:
        hub = hub_of()
        q = _build_fixed_source_rhs(np.ones((hub.quad.n_ordinates, 2, nx)), hub)
        posing = hub.source_posing(q)
        _require(posing.source.space == hub.pencil.at(1.0).domain, "the lifted source's DERIVED space equals the pencil's domain (content identity)")
        _require(posing.source.space == hub.system.space, "…and the builder-minted coupled space")
        _require(hub.system.space.zeros().space == hub.system.space, "the space's own zero element lives on it")

    def test_negative_a_cross_hub_source_is_REFUSED_at_the_posing(self) -> None:
        """The ends check is LIVE now: a source built on another Problem's space
        cannot be posed on this one (until step 3 the coupled state had no
        ``.space`` and this raised nothing)."""
        here, there = _slab(4.0), _slab(6.0)
        foreign = _build_fixed_source_rhs(np.ones((there.quad.n_ordinates, 2, 8)), there)
        with pytest.raises(ValueError, match="pose q on the operator's ends"):
            here.source_posing(foreign)
