"""MoC refuses a mesh it would read as something else (#514), as a scope guard.

MoC reads any mesh as the solid Wigner–Seitz cylinder of a square pin cell, so a slab or a hollow cylinder would be traced as that
cylinder. The refusal is a declared scope boundary (``SCOPE-BOUNDARY[guard]``, #514), read before the declaration, so it stays a
``NotImplementedError`` after wave 2 (the specification's S3.13). The last row pins that ``solve_moc``'s default mesh is one MoC
admits (ERR-093).

Moved from ``tests/gates/mesh/test_dropped_laws_are_refused.py`` by the layer-DAG re-filing (Phase 1b, 2026-10-26): the
refusal is MoC's, so the row lives with MoC's gates; the bodies are unchanged. The two mesh builders the three methods share
live in ``tests/gates/transport/_dropped_laws_helpers.py``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import CoordSystem, MirrorDeck, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.moc.geometry import MOCMesh
from orpheus.moc.quadrature import MOCQuadrature
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates.transport._dropped_laws_helpers import _hollow, _slab

pytestmark = pytest.mark.foundation

_SPECULAR = {0: WallCondition.specular(1.0)}


# ── S3.13: MoC reads a solid Wigner–Seitz cylinder ──────────────────────────────────────────────────────────────────


class TestMOC:
    _QUADRATURE = MOCQuadrature.create(n_azi=4, n_polar=1)

    def test_a_cartesian_mesh_is_refused(self):
        """The scope guard (#514, ``SCOPE-BOUNDARY[guard]``) is read before the declaration: a slab is not a cylinder whatever its faces
        carry (here the mirror decks MoC would not realize either)."""
        mesh, conditions = _slab(MirrorDeck(), MirrorDeck())
        with pytest.raises(NotImplementedError, match="#514"):
            MOCMesh(mesh, self._QUADRATURE, boundary_conditions=conditions)

    def test_a_hollow_cylinder_is_refused(self):
        mesh = _hollow(CoordSystem.CYLINDRICAL)
        with pytest.raises(NotImplementedError, match="#514"):
            MOCMesh(mesh, self._QUADRATURE, boundary_conditions={1: WallCondition.specular(1.0), 0: WallCondition.specular(1.0)})

    def test_the_solid_pin_cell_builds(self):
        geometry = StructuredGeometry.cylinder(
            (0.0, 0.9, 1.1, 3.6 / np.sqrt(np.pi)), (2, 1, 0), outer=0,
        )
        mesh = Mesher(geometry).partition((
            CellsByCount.uniform_volume(2), CellsByCount.uniform_volume(1),
            CellsByCount.uniform_volume(2),
        )).mesh
        MOCMesh(mesh, self._QUADRATURE, boundary_conditions=_SPECULAR)

    @pytest.mark.catches("ERR-093")
    def test_the_default_mesh_builds(self):
        """``solve_moc``'s default mesh is one MoC admits.

        Its first red: the default was the Wigner-Seitz cell, whose model law
        is white, so ``solve_moc(materials)`` raised "MOC solver does not
        support boundary condition 'white'" before a single ray was traced.
        The default mesh now carries the tag 0 on its outer face and ``solve_moc`` binds it to the unit specular wall MoC links.
        """
        from orpheus.moc.solver import default_pin_cell_mesh

        mesh = default_pin_cell_mesh()
        assert dict(mesh.face_roles) == {"xmax": 0}
        MOCMesh(mesh, self._QUADRATURE, boundary_conditions=_SPECULAR)
