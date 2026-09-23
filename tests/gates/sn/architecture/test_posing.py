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
