r"""The CP builders of the S4 gates on the CP boundary table (wave 2, WF2.3 slice S4).

Not a test module (no ``test_`` prefix). Split on 2026-10-26 from ``tests/gates/transport/_wave2_s4_helpers.py`` (layer DAG ruling 27:
as import statements, that transport helper's CP accessors were a transport gate importing ``orpheus.cp``); the method-agnostic
builders (sites, the carry fixture, ``quiet``, ``bc``) stay there. The production module ``orpheus.cp.solver`` is imported at the top, as a module
object, and the accessors read their attributes at call time, so a monkeypatch of a module attribute reaches every gate.
"""

from __future__ import annotations

from typing import Any

from orpheus.cp import solver as cp_solver
from tests.gates.transport import _carry_cases as carry
from tests.gates.transport._wave2_s4_helpers import chart_mesh, quiet


def cp_table() -> Any:
    """``CP_REALIZERS``, the sealed CP realizer table (``orpheus/cp/solver.py``)."""
    return cp_solver.CP_REALIZERS


def cp_mesh(chart: str, *, outer: Any = 0, left: Any = None) -> Any:
    """The carried CP body of ``chart``: a slab has a mirror at the left face (CP's kernel images one) and ``outer`` on the right; a
    solid radial body has ``outer`` on its one face. ``left`` overrides the slab's left role (a tag or a deck)."""
    from orpheus.geometry.deck import MirrorDeck

    if chart == "slab":
        return chart_mesh(chart, left=left if left is not None else MirrorDeck(), right=outer)
    return chart_mesh(chart, outer=outer)


def cp_solve(mesh: Any, conditions: dict[int, Any] | None = None, params: Any = None) -> Any:
    """``solve_cp`` at ``CPParams()`` (what the HEAD capture ran) with ``boundary_conditions=`` iff given; stdout silenced."""
    kwargs: dict[str, Any] = {} if conditions is None else {"boundary_conditions": conditions}
    with quiet():
        return cp_solver.solve_cp(carry.materials(), mesh, params if params is not None else cp_solver.CPParams(), **kwargs)


def cp_head(mesh: Any, conditions: dict[int, Any] | None = None) -> Any:
    """A ``CPMesh`` over ``mesh`` with ``boundary_conditions=`` iff given (construction only: the table is read here)."""
    return cp_solver.CPMesh(mesh) if conditions is None else cp_solver.CPMesh(mesh, boundary_conditions=conditions)
