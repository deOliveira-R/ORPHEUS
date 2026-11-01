r"""Shared builders of the S3 gates on the diffusion boundary table (wave 2, WF2.3 slice S3).

Not a test module (no ``test_`` prefix). Spec: ``.claude/plans/wave2_boundary_carve.md`` sections 2.7 (``ZeroFluxCondition``
and ``DIFFUSION_REALIZERS``), 4.2 (the diffusion column), 6.0 (the carry fixture) and 6.3 (the eight done-when gates).

How this module reads the production modules
============================================

The diffusion modules are imported at the top, as MODULES, and the accessors below read their attributes at call time, so a
monkeypatch of a module attribute reaches every gate (layer DAG ruling 27: an import whose target is known when the test is written
is an import statement, which the test-layer gate sees). Until 2026-10-26 they were imported by string at first use, because at the
base tip of slice S3 every ``orpheus.diffusion`` import raised ``ImportError`` (``augmented_mesh`` imported a name S2b removed from
``orpheus.transport.method``); that reason ended with the slice.

What it builds
==============

* the accessors of the production names (:func:`table`, :func:`zero_flux_condition`, :func:`minimal_space`);
* the faces' :class:`~orpheus.transport.boundary_conditions.FaceSite` of the three geometric classes (:func:`site`), the same
  keys ``tests/gates/transport/_admission_helpers.site`` builds;
* the carried diffusion problems of ``tests/gates/transport/_carry_cases.py`` in the NEW spelling (:func:`chart_mesh`,
  :func:`build`, :func:`solve`): the SAME charts, breakpoints, materials ``(1, 0)``, cell counts and solver settings the HEAD
  capture script used (``_capture_carry_head._chart_mesh`` (deleted in WF2.4; recover it with ``git show 197549f6:tests/gates/transport/_capture_carry_head.py``)), because a bit-equal ``k`` is meaningful only for the same problem;
* :func:`realized_action`, the action of a realized operator on a fixed vector: the realized operators are ``OperatorSum`` values
  with no ``__eq__``, so a gate reads them BY ACTION on a fixed, non-symmetric vector, and by type only where the type is
  the claim.
"""

from __future__ import annotations

import types
from typing import Any

import numpy as np

from orpheus.diffusion import augmented_mesh, boundary_realizer, method_space, solver, zero_flux
from orpheus.transport import boundary_conditions
from tests.gates.transport import _carry_cases as carry
from tests.gates.transport._wave2_s4_helpers import chart_mesh, summary  # noqa: F401  (re-exported: re-homed in layer DAG Phase 4)

#: A fixed non-symmetric probe of a face operator's input (one outflow current per group; two groups here, but the realized
#: scalar-multiple-of-identity is shape-agnostic).
PROBE = np.array([0.7, 1.9, 3.1, 0.2, 5.3, 2.3])


def bc() -> types.ModuleType:
    """``orpheus.transport.boundary_conditions``."""
    return boundary_conditions


def table() -> Any:
    """``DIFFUSION_REALIZERS``, the sealed diffusion realizer table."""
    return boundary_realizer.DIFFUSION_REALIZERS


def zero_flux_condition() -> Any:
    """The class ``ZeroFluxCondition`` (``orpheus/diffusion/zero_flux.py``)."""
    return zero_flux.ZeroFluxCondition


def minimal_space() -> Any:
    """The metadata-free method space the albedo-family realizers read nothing from."""
    return method_space.DiffusionMethodSpace.minimal()


def site(kind: str, face: str = "xmax") -> Any:
    """The three geometric sites of section 4.2 as ``FaceSite`` values (a face name overrides the default for ``CART``)."""
    from orpheus.geometry.coord import CoordSystem

    B = bc()
    match kind:
        case "CART":
            return B.FaceSite(face, CoordSystem.CARTESIAN, 1)
        case "RAD_OUTER":
            return B.FaceSite("xmax", CoordSystem.CYLINDRICAL, 1)
        case "RAD_INNER":
            return B.FaceSite("xmin", CoordSystem.CYLINDRICAL, 1)
    raise ValueError(kind)


SITES = ("CART", "RAD_OUTER", "RAD_INNER")
#: The sites with a trace slot, where a condition is REALIZED: the hollow inner face has none (#511), so the table admits only
#: ``specular(1.0)`` there and refuses the rest (the admission matrix's diffusion column pins that site, cell by cell).
TRACE_SITES = ("CART", "RAD_OUTER")


def roles_and_conditions(chart: str, condition: Any, *, mirror: bool = False) -> tuple[dict[str, Any], dict[int, Any]]:
    """The roles of the faces and the ``{tag: condition}`` of one wall ``condition`` on every outer face of ``chart``.

    ``mirror=True`` is the perfect mirror: a ``MirrorDeck`` on each slab face, and (rule R2) the tag bound to
    ``WallCondition.specular(1.0)`` on a curved face. Otherwise one tag, bound to ``condition``, on every outer face.
    """
    from orpheus.geometry.deck import MirrorDeck

    B = bc()
    coord, _ = carry.CHARTS[chart]
    if mirror and coord == "CARTESIAN":
        return {"left": MirrorDeck(), "right": MirrorDeck()}, {}
    bound = B.WallCondition.specular(1.0) if mirror else condition
    if coord == "CARTESIAN":
        return {"left": 0, "right": 0}, {0: bound}
    return {"outer": 0}, {0: bound}


def build(chart: str, declaration: carry.Declaration) -> tuple[Any, dict[int, Any]]:
    """``(mesh, conditions)`` of a carried ``(chart, declaration)`` in the new spelling (rules R2 and R4 of section 6.0)."""
    B = bc()
    match declaration:
        case ("vacuum",):
            roles, conditions = roles_and_conditions(chart, B.WallCondition.vacuum())
        case ("mirror",):
            roles, conditions = roles_and_conditions(chart, None, mirror=True)
        case ("albedo", a):  # the bare albedo of the old diffusion-only tests: R4, the diffuse exterior
            roles, conditions = roles_and_conditions(chart, B.WallCondition.white(a))
        case ("zero_flux",):
            roles, conditions = roles_and_conditions(chart, zero_flux_condition()())
        case _:
            raise ValueError(declaration)
    return chart_mesh(chart, **roles), conditions


def solve(mesh: Any, conditions: dict[int, Any] | None = None, materials: Any = None) -> Any:
    """``solve_diffusion_1d`` at its default tolerances (what the HEAD capture ran), with ``boundary_conditions=`` iff given."""
    mats = materials if materials is not None else carry.materials()
    if conditions is None:
        return solver.solve_diffusion_1d(mats, mesh)
    return solver.solve_diffusion_1d(mats, mesh, boundary_conditions=conditions)


def diffusion_mesh(mesh: Any, conditions: dict[int, Any] | None = None, materials: Any = None) -> Any:
    """A ``DiffusionMesh`` over ``mesh`` with ``boundary_conditions=`` iff given."""
    mats = materials if materials is not None else carry.materials()
    if conditions is None:
        return augmented_mesh.DiffusionMesh(mesh, mats)
    return augmented_mesh.DiffusionMesh(mesh, mats, boundary_conditions=conditions)


def realized_action(operator: Any) -> np.ndarray:
    """What a realized face operator does to :data:`PROBE` (the operators are read by action, never by ``==``)."""
    return np.asarray(operator.apply(PROBE))
