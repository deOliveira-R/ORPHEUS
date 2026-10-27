r"""The MoC builders of the S4 gates on the MoC boundary table (wave 2, WF2.3 slice S4).

Not a test module (no ``test_`` prefix). Split on 2026-10-26 from ``tests/gates/transport/_wave2_s4_helpers.py`` (layer DAG ruling 27:
as import statements, that transport helper's MoC accessors were a transport gate importing ``orpheus.moc``); the method-agnostic
builders (sites, the carry fixture, ``quiet``, ``bc``) stay there. The production modules ``orpheus.moc.geometry``, ``.quadrature`` and ``.solver`` are imported at the top, as module
objects, and the accessors read their attributes at call time, so a monkeypatch of a module attribute reaches every gate.
"""

from __future__ import annotations

from typing import Any

from orpheus.moc import geometry as moc_geometry
from orpheus.moc import quadrature as moc_quadrature
from orpheus.moc import solver as moc_solver
from tests.gates.transport import _carry_cases as carry
from tests.gates.transport._wave2_s4_helpers import chart_mesh, quiet

#: The reason text the MoC table's one named approximation carries (section 4.2, the MoC column; pinned exactly).
MOC_DOMAIN_NAME = "an equal-area square cell with reflecting flat faces, not the declared cylinder"


def moc_table() -> Any:
    """``MOC_REALIZERS``, the sealed MoC realizer table (``orpheus/moc/geometry.py``)."""
    return moc_geometry.MOC_REALIZERS


def moc_mesh1d(*, outer: Any = 0) -> Any:
    """The carried MoC cylinder: breakpoints ``(0, 1, 2)``, materials ``(1, 0)``, ``MOC_CELLS`` cells, ``outer`` on its one face."""
    return chart_mesh("cylinder", outer=outer, cells=carry.MOC_CELLS)


def moc_solve(mesh: Any, conditions: dict[int, Any] | None = None) -> Any:
    """``solve_moc`` at the carried ray settings with ``boundary_conditions=`` iff given; stdout silenced."""
    kwargs: dict[str, Any] = {} if conditions is None else {"boundary_conditions": conditions}
    with quiet():
        return moc_solver.solve_moc(carry.materials(), mesh, **carry.MOC_ARGS, **kwargs)


def moc_head(mesh: Any, conditions: dict[int, Any] | None = None) -> Any:
    """A ``MOCMesh`` over ``mesh`` (4 azimuthal angles, one polar) with ``boundary_conditions=`` iff given."""
    quadrature = moc_quadrature.MOCQuadrature.create(n_azi=4, n_polar=1)
    return moc_geometry.MOCMesh(mesh, quadrature) if conditions is None else moc_geometry.MOCMesh(mesh, quadrature, boundary_conditions=conditions)
