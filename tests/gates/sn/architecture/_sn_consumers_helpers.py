r"""Shared builders of the S2b gates that pose a problem through ``SNProblem`` and the ``solve_sn*`` entries (wave 2).

Not a test module (no ``test_`` prefix). Spec: ``.claude/plans/wave2_boundary_carve.md`` sections 2.10 and 6.2 (S2b gates
2 to 13). Every problem is spelled in the NEW vocabulary: a face carries a role (an ``int`` tag or a deck), the
conditions are a ``{tag: WallCondition}`` mapping handed to the entry as ``boundary_conditions=``.

It imports ``orpheus.transport.boundary_conditions`` at the top, as a module (layer DAG ruling 27; until 2026-10-26 it was
imported by string, because at the base tip of slice S2b ``import orpheus.sn`` raised), and the solver entries inside the
builders that read them. The bodies are carried in :mod:`tests.gates.transport._carry_cases`
(the same data the HEAD capture script reads), so the carry gate and the other gates pose the same problems.
"""

from __future__ import annotations

import types
from typing import Any

from orpheus.transport import boundary_conditions
from tests.gates.transport import _carry_cases as carry


def bc() -> types.ModuleType:
    """``orpheus.transport.boundary_conditions``."""
    return boundary_conditions


def mesh_of(body: str, *, left: Any = None, right: Any = None, outer: Any = None, inner: Any = None,
            breakpoints: tuple[float, ...] | None = None, mat_ids: tuple[int, ...] | None = None, cells: tuple[int, ...] | None = None) -> Any:
    """The mesh of a carried ``body`` with the given ROLES (a slab: ``left`` and ``right``; a radial body: ``outer`` and,
    when hollow, ``inner``). ``breakpoints``/``mat_ids``/``cells`` override the body's."""
    from orpheus.geometry import StructuredGeometry
    from orpheus.mesh import CellsByCount, Mesher

    b = carry.BODIES[body]
    bps = breakpoints if breakpoints is not None else b.breakpoints
    mats = mat_ids if mat_ids is not None else b.mat_ids
    n = cells if cells is not None else b.cells
    if b.coord == "CARTESIAN":
        geometry = StructuredGeometry.slab(bps, mats, left=left, right=right)
    elif b.coord == "CYLINDRICAL":
        geometry = StructuredGeometry.cylinder(bps, mats, outer=outer, inner=inner)
    else:
        geometry = StructuredGeometry.sphere(bps, mats, outer=outer, inner=inner)
    return Mesher(geometry).partition(tuple(CellsByCount.uniform_volume(k) for k in n)).mesh


def quadrature_of(body: str, ordinates: int | None = None) -> Any:
    return carry.quadrature(carry.BODIES[body], ordinates)


def problem(mesh: Any, quadrature: Any, conditions: Any = None, materials: Any = None) -> Any:
    """An ``SNProblem`` over ``mesh`` with ``conditions`` (``None`` omits the keyword, the empty default)."""
    from orpheus.sn.problem import SNProblem

    mats = materials if materials is not None else carry.materials()
    if conditions is None:
        return SNProblem(mesh, quadrature, mats)
    return SNProblem(mesh, quadrature, mats, boundary_conditions=conditions)


def wall(kind: str, value: float | None = None) -> Any:
    """``vacuum``, ``specular(a)``, ``white(a)``, ``uniform(v)`` as the shared exterior type."""
    w = bc().WallCondition
    match kind:
        case "vacuum":
            return w.vacuum()
        case "specular":
            return w.specular(value)
        case "white":
            return w.white(value)
        case "uniform":
            return w.uniform(value)
    raise ValueError(kind)


def roles_and_conditions(body: str, declaration: tuple) -> tuple[dict[str, Any], dict[int, Any]]:
    """The roles of the faces and the ``{tag: condition}`` of a carried ``(body, declaration)`` in the NEW spelling.

    ``mirror`` is a mirror deck on a Cartesian face and ``specular(1.0)`` on a curved one (rule R2: a reflective
    declaration on a curved face is a tag bound to the unit specular chain); a hollow radial body carries the reflective
    cavity its inner face has always meant as tag 1 bound to ``specular(1.0)`` (rule R6)."""
    from orpheus.geometry.deck import MirrorDeck

    b = carry.BODIES[body]
    kind = declaration[0]
    conditions: dict[int, Any] = {}
    if kind == "mirror" and b.coord == "CARTESIAN":
        role: Any = MirrorDeck()
    else:
        role = 0
        conditions[0] = wall("specular", 1.0) if kind == "mirror" else wall(kind, declaration[1] if len(declaration) > 1 else None)
    roles: dict[str, Any] = {}
    if b.coord == "CARTESIAN":
        roles = {"left": role, "right": role}
    else:
        roles = {"outer": role}
        if b.hollow:
            roles["inner"] = 1
            conditions[1] = wall("specular", 1.0)
    return roles, conditions


def build(body: str, declaration: tuple) -> tuple[Any, Any, dict[int, Any]]:
    """``(mesh, quadrature, conditions)`` of a carried ``(body, declaration)`` in the new spelling."""
    roles, conditions = roles_and_conditions(body, declaration)
    return mesh_of(body, **roles), quadrature_of(body), conditions
