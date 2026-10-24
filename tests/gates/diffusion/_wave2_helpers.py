r"""Shared builders of the S3 gates on the diffusion boundary table (wave 2, WF2.3 slice S3).

Not a test module (no ``test_`` prefix). Spec: ``.claude/plans/wave2_boundary_carve.md`` sections 2.7 (``ZeroFluxCondition``
and ``DIFFUSION_REALIZERS``), 4.2 (the diffusion column), 6.0 (the carry fixture) and 6.3 (the eight done-when gates).

Why this module imports no diffusion module at the top
=====================================================

At the base tip of slice S3 every ``orpheus.diffusion.*`` import raises ``ImportError``: the package ``__init__`` imports the
solver, which imports ``augmented_mesh``, which imports ``resolve_boundary_conditions`` from ``orpheus.transport.method`` (S2b
removed it). A collection-time import would turn every gate into the same collection error and hide which gate fails for which
reason (``vv-principles`` anti-pattern 17c). So every diffusion module is read at first use, by the accessors below, inside each
test, and each gate's first red names its own line: the missing module ``orpheus.diffusion.zero_flux``, the missing table
``DIFFUSION_REALIZERS``, the missing keyword ``boundary_conditions``.

What it builds
==============

* the lazy accessors of the production names (:func:`table`, :func:`zero_flux_condition`, :func:`diffusion_module`);
* the faces' :class:`~orpheus.transport.boundary_conditions.FaceSite` of the three geometric classes (:func:`site`), the same
  keys ``tests/gates/transport/_admission_helpers.site`` builds;
* the carried diffusion problems of ``tests/gates/transport/_carry_cases.py`` in the NEW spelling (:func:`chart_mesh`,
  :func:`build`, :func:`solve`): the SAME charts, breakpoints, materials ``(1, 0)``, cell counts and solver settings the HEAD
  capture script used (``_capture_carry_head._chart_mesh``), because a bit-equal ``k`` is meaningful only for the same problem;
* :func:`realized_action`, the action of a realized operator on a fixed vector: the realized operators are ``OperatorSum`` values
  with no ``__eq__``, so a gate reads them BY ACTION on a fixed, non-symmetric vector, and by type only where the type is
  the claim.
"""

from __future__ import annotations

import importlib
import types
from typing import Any

import numpy as np

from tests.gates.transport import _carry_cases as carry

#: A fixed non-symmetric probe of a face operator's input (one outflow current per group; two groups here, but the realized
#: scalar-multiple-of-identity is shape-agnostic).
PROBE = np.array([0.7, 1.9, 3.1, 0.2, 5.3, 2.3])


def bc() -> types.ModuleType:
    """``orpheus.transport.boundary_conditions``."""
    return importlib.import_module("orpheus.transport.boundary_conditions")


def diffusion_module(name: str) -> types.ModuleType:
    """``orpheus.diffusion.<name>`` (read at first use: at the base tip the package ``__init__`` raises)."""
    return importlib.import_module(f"orpheus.diffusion.{name}")


def table() -> Any:
    """``DIFFUSION_REALIZERS``, the sealed diffusion realizer table."""
    return diffusion_module("boundary_realizer").DIFFUSION_REALIZERS


def zero_flux_condition() -> Any:
    """The class ``ZeroFluxCondition`` (``orpheus/diffusion/zero_flux.py``)."""
    return diffusion_module("zero_flux").ZeroFluxCondition


def minimal_space() -> Any:
    """The metadata-free method space the albedo-family realizers read nothing from."""
    return diffusion_module("method_space").DiffusionMethodSpace.minimal()


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


def chart_mesh(chart: str, *, left: Any = None, right: Any = None, outer: Any = None,
               breakpoints: tuple[float, ...] | None = None, mat_ids: tuple[int, ...] = (1, 0),
               cells: tuple[int, ...] = carry.OTHER_CELLS) -> Any:
    """The ``Mesh1D`` of a carried diffusion chart (``slab`` | ``cylinder`` | ``sphere``) with the given face ROLES.

    A slab takes ``left`` and ``right``; a solid radial body takes ``outer``. A role is an ``int`` tag or a deck.
    """
    from orpheus.geometry import StructuredGeometry
    from orpheus.mesh import CellsByCount, Mesher

    coord, default_breakpoints = carry.CHARTS[chart]
    bps = breakpoints if breakpoints is not None else default_breakpoints
    if coord == "CARTESIAN":
        geometry = StructuredGeometry.slab(bps, mat_ids, left=left, right=right)
    elif coord == "CYLINDRICAL":
        geometry = StructuredGeometry.cylinder(bps, mat_ids, outer=outer)
    else:
        geometry = StructuredGeometry.sphere(bps, mat_ids, outer=outer)
    return Mesher(geometry).partition(tuple(CellsByCount.uniform_volume(n) for n in cells)).mesh


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
    solver = diffusion_module("solver")
    mats = materials if materials is not None else carry.materials()
    if conditions is None:
        return solver.solve_diffusion_1d(mats, mesh)
    return solver.solve_diffusion_1d(mats, mesh, boundary_conditions=conditions)


def diffusion_mesh(mesh: Any, conditions: dict[int, Any] | None = None, materials: Any = None) -> Any:
    """A ``DiffusionMesh`` over ``mesh`` with ``boundary_conditions=`` iff given."""
    augmented = diffusion_module("augmented_mesh")
    mats = materials if materials is not None else carry.materials()
    if conditions is None:
        return augmented.DiffusionMesh(mesh, mats)
    return augmented.DiffusionMesh(mesh, mats, boundary_conditions=conditions)


def summary(admission: Any) -> tuple:
    """The pinned shape of an admission, the one ``tests/gates/transport/test_admission_matrix.summary`` reads: ``("R", faithful)``,
    ``("A", names, issues)`` or ``("F", axis name, issue)`` (the reason text is prose and is not pinned)."""
    B = bc()
    if isinstance(admission, B.Realized):
        return ("R", admission.faithful)
    if isinstance(admission, B.Approximated):
        return ("A", tuple(sorted(admission.names)), tuple(sorted(admission.issues)))
    if not isinstance(admission, B.Refused):
        raise TypeError(f"{admission!r} is not an admission")
    return ("F", admission.axis.name, admission.issue)


def realized_action(operator: Any) -> np.ndarray:
    """What a realized face operator does to :data:`PROBE` (the operators are read by action, never by ``==``)."""
    return np.asarray(operator.apply(PROBE))
