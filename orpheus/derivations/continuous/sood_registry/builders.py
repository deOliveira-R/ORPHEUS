r"""Builders: case → ``(materials, mesh, params)`` for production solvers.

Production solvers (:func:`orpheus.cp.solver.solve_cp`,
:func:`orpheus.sn.solver.solve_sn`) take ``materials: dict[int, Mixture]``
+ ``mesh: Mesh1D`` + per-solver params. This module provides ergonomic
helpers so a test can do:

>>> from orpheus.derivations.continuous.sood_registry import (
...     SOOD2003_CASES, build_materials, build_mesh, build_cp_params
... )
>>> case = SOOD2003_CASES["Ua-1-0-SL"]
>>> materials = build_materials(case)
>>> mesh = build_mesh(case, n_cells=64)
>>> # then solve_cp(materials, mesh, build_cp_params(case))

Each helper just unpacks fields off the case object. They exist so
production-solver consumer tests don't need to know the schema.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.mesh import Mesh1D

if TYPE_CHECKING:
    from .case import La13511Case


def build_materials(case: "La13511Case") -> dict[int, Mixture]:
    """Return the case's materials dict, ready for production solvers."""
    return case.materials


def build_mesh(case: "La13511Case", n_cells: int = 64) -> Mesh1D:
    """Build a :class:`Mesh1D` for ``case`` at ``n_cells`` refinement.

    Constructs the :class:`StructuredGeometry` via
    :meth:`La13511Case.to_geometry` (raises for ``case.geometry_kind ==
    "infinite"`` — infinite-medium cases have no spatial mesh and
    should be consumed via :func:`build_materials` alone) and pairs
    meshes every region with ``n_cells`` equal-volume cells
    (:meth:`CellsByCount.uniform_volume <orpheus.mesh.partition.CellsByCount.uniform_volume>`).
    """
    from orpheus.mesh import CellsByCount, Mesher

    return Mesher(case.to_geometry()).partition(CellsByCount.uniform_volume(n_cells)).mesh


def build_cp_params(case: "La13511Case", **kwargs):  # type: ignore[no-untyped-def]
    """Build a :class:`CPParams` with sensible defaults for ``case``.

    Imported lazily so importing :mod:`sood_registry` doesn't pull in
    the CP solver's transitive dependencies.
    """
    from orpheus.cp.solver import CPParams
    return CPParams(**kwargs)
