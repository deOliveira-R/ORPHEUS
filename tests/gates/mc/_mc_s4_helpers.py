r"""The MC builders of the S4 gates on the MC boundary table (wave 2, WF2.3 slice S4).

Not a test module (no ``test_`` prefix). Split on 2026-10-26 from ``tests/gates/transport/_wave2_s4_helpers.py`` (layer DAG ruling 27:
as import statements, that transport helper's MC accessors were a transport gate importing ``orpheus.mc``); the method-agnostic
builders (sites, the carry fixture, ``quiet``, ``bc``) stay there. The production module ``orpheus.mc.solver`` is imported at the top, as a module
object, and the accessors read their attributes at call time, so a monkeypatch of a module attribute reaches every gate.
"""

from __future__ import annotations

from typing import Any

from orpheus.mc import solver as mc_solver
from tests.gates.transport import _carry_cases as carry
from tests.gates.transport._wave2_s4_helpers import chart_mesh, quiet

#: A fragment MC's TOPOLOGY refusal must carry (section 4.3 delta 9: the torus, and the issue it is tracked under).
MC_TORUS_FRAGMENT = "torus"


def mc_table() -> Any:
    """``MC_REALIZERS``, the sealed MC realizer table (``orpheus/mc/solver.py``)."""
    return mc_solver.MC_REALIZERS


def mc_slab_mesh1d(*, left: Any, right: Any) -> Any:
    """The carried MC slab: breakpoints ``(0, 1.5, 3)``, materials ``(2, 0)`` (fuel, moderator), 6 + 6 cells; the faces' roles as given."""
    return chart_mesh("slab", left=left, right=right, mat_ids=(2, 0))


def mc_materials() -> dict[int, Any]:
    """``{2: fuel, 0: moderator}``: the carried mixtures under MC's ids (the HEAD capture's ``mc_mats``)."""
    mats = carry.materials()
    return {2: mats[1], 0: mats[0]}


def mc_head(mesh: Any, pitch: float, conditions: dict[int, Any] | None = None) -> Any:
    """An ``MCMesh`` over ``mesh`` with ``boundary_conditions=`` iff given."""
    return mc_solver.MCMesh(mesh, pitch=pitch) if conditions is None else mc_solver.MCMesh(mesh, pitch=pitch, boundary_conditions=conditions)


def mc_solve(geometry: Any) -> Any:
    """``solve_monte_carlo`` with the carried parameters (fixed seed), the geometry given; stdout silenced."""
    with quiet():
        return mc_solver.solve_monte_carlo(mc_materials(), mc_solver.MCParams(**carry.MC_PARAMS, geometry=geometry))
