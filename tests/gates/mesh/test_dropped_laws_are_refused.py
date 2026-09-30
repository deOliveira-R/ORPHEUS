"""CP, MoC and MC refuse a declared law they would drop (#513, #514).

Each of the three reads only part of a mesh's declared boundary: CP and MC
resolve only the outer law, and MoC reads any mesh as the solid Wigner–Seitz
cylinder of a square pin cell. Before P1 step 2 a law they did not read was
dropped silently (``[M]`` 2026-09-25, the P1 specification's probes
``cp_slab_left_law.py`` and ``partial_law_readers.py``). Step 2 makes a
hollow body and a declared inner law expressible on a
:class:`~orpheus.geometry.StructuredGeometry`, so each method's refusal lands
with it (the specification's S3.12 to S3.14, moved into step 2 by the user's
ruling of 2026-09-29). Each refusal is a declared scope boundary. Since
P1 step 3b ``None`` is not a law (the mesh refuses it), so the rows that
admitted an undeclared left law are gone; each method's admitted pair is
the declared one.
"""
from __future__ import annotations

import pytest

from orpheus.cp.solver import CPMesh
from orpheus.geometry import BC, CoordSystem, StructuredGeometry
from orpheus.mc.solver import MCMesh
from orpheus.mesh import CellsByCount, Mesh1D, Mesher
from orpheus.moc.geometry import MOCMesh
from orpheus.moc.quadrature import MOCQuadrature

pytestmark = pytest.mark.foundation


def _slab(left, right) -> Mesh1D:
    geometry = StructuredGeometry.slab((0.0, 2.0), (0,), left=left, right=right)
    return Mesher(geometry).partition(CellsByCount.uniform_width(4)).mesh


def _hollow(coord: CoordSystem, inner, outer) -> Mesh1D:
    geometry = StructuredGeometry(
        coord=coord, breakpoints=(0.5, 2.0), mat_ids=(0,), boundaries=(inner, outer),
    )
    return Mesher(geometry).partition(CellsByCount.uniform_volume(4)).mesh


# ── S3.12: CP reads only the outer law ─────────────────────────────────


class TestCP:
    @pytest.mark.parametrize(
        "left, right", [(BC.vacuum, BC.white), (BC.white, BC.vacuum)],
        ids=["vacuum-white", "white-vacuum"],
    )
    def test_a_slab_whose_laws_differ_is_refused(self, left, right):
        with pytest.raises(NotImplementedError, match="#513"):
            CPMesh(_slab(left, right))

    @pytest.mark.parametrize("left", [BC.white], ids=["equal"])
    def test_a_slab_whose_left_law_is_read_builds(self, left):
        CPMesh(_slab(left, BC.white))

    @pytest.mark.parametrize(
        "coord", [CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL], ids=lambda c: c.name.lower(),
    )
    def test_a_hollow_body_with_an_inner_law_is_refused(self, coord):
        with pytest.raises(NotImplementedError, match="#513"):
            CPMesh(_hollow(coord, BC.white, BC.white))


# ── S3.13: MoC reads a solid Wigner–Seitz cylinder ─────────────────────


class TestMOC:
    _QUADRATURE = MOCQuadrature.create(n_azi=4, n_polar=1)

    def test_a_cartesian_mesh_is_refused(self):
        with pytest.raises(NotImplementedError, match="#514"):
            MOCMesh(_slab(BC.reflective, BC.reflective), self._QUADRATURE)

    def test_a_hollow_cylinder_is_refused(self):
        mesh = _hollow(CoordSystem.CYLINDRICAL, BC.reflective, BC.reflective)
        with pytest.raises(NotImplementedError, match="#514"):
            MOCMesh(mesh, self._QUADRATURE)

    def test_the_solid_pin_cell_builds(self):
        geometry = StructuredGeometry.wigner_seitz_pin_cell(
            boundaries=(BC.reflective,),
        )
        mesh = Mesher(geometry).partition((
            CellsByCount.uniform_volume(2), CellsByCount.uniform_volume(1),
            CellsByCount.uniform_volume(2),
        )).mesh
        MOCMesh(mesh, self._QUADRATURE)


# ── S3.14: MC reads periodic laws on a slab or a solid cylinder ───────


class TestMC:
    @pytest.mark.parametrize("left", [BC.vacuum, BC.reflective], ids=["vacuum", "reflective"])
    def test_a_non_periodic_left_law_is_refused(self, left):
        with pytest.raises(NotImplementedError, match="#513"):
            MCMesh(_slab(left, BC("periodic")), pitch=2.0)

    @pytest.mark.parametrize("left", [BC("periodic")], ids=["periodic"])
    def test_a_periodic_or_undeclared_left_law_builds(self, left):
        MCMesh(_slab(left, BC("periodic")), pitch=2.0)

    def test_a_hollow_cylinder_is_refused(self):
        """Its lookup would fill the cavity with the innermost material."""
        mesh = _hollow(CoordSystem.CYLINDRICAL, BC("periodic"), BC("periodic"))
        with pytest.raises(NotImplementedError, match="cavity"):
            MCMesh(mesh, pitch=4.0)
