r"""Builders: case → ``(materials, mesh, params)`` for production solvers.

Production solvers (:func:`orpheus.cp.solver.solve_cp`,
:func:`orpheus.sn.solver.solve_sn`) take ``materials: dict[int, Mixture]``
+ ``mesh: Mesh1D`` + per-solver params. This module provides ergonomic
helpers so a test can do:

>>> from orpheus.derivations.continuous.sood_registry import (
...     LA13511_CASES, build_materials, build_mesh, build_cp_params
... )
>>> case = LA13511_CASES["Ua-1-0-SL"]
>>> materials = build_materials(case)
>>> mesh = build_mesh(case, n_cells=64)
>>> # then solve_cp(materials, mesh, build_cp_params(case))

Each helper just unpacks fields off the case object. They exist so
production-solver consumer tests don't need to know the schema.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry.mesh import Mesh1D

if TYPE_CHECKING:
    from .la13511 import La13511Case


def build_materials(case: "La13511Case") -> dict[int, Mixture]:
    """Return the case's materials dict, ready for production solvers."""
    return case.materials


def build_mesh(case: "La13511Case", n_cells: int = 64) -> Mesh1D:
    """Build a :class:`Mesh1D` for ``case`` at ``n_cells`` refinement.

    For ``case.geometry == "infinite"`` raises :class:`ValueError` —
    infinite-medium cases have no spatial mesh and should be consumed
    via :func:`build_materials` alone.
    """
    return case.geometry_spec.build(n_cells=n_cells)


def build_cp_params(case: "La13511Case", **kwargs):  # type: ignore[no-untyped-def]
    """Build a :class:`CPParams` with sensible defaults for ``case``.

    Imported lazily so importing :mod:`sood_registry` doesn't pull in
    the CP solver's transitive dependencies.
    """
    from orpheus.cp.solver import CPParams
    return CPParams(**kwargs)
