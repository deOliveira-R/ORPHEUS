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

import re

import numpy as np
import pytest

import orpheus.cp.solver as cp_solver
from orpheus.cp.solver import CPMesh, solve_cp
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem, StructuredGeometry
from orpheus.geometry.boundary import AlbedoBoundary, ReflectiveBoundary, SpecularReturn
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
        "left, right",
        [
            (BC.vacuum, BC.white), (BC.white, BC.vacuum), (BC.white, BC.white),
            (AlbedoBoundary(0.7, SpecularReturn("x")), BC.white),
            (ReflectiveBoundary("y"), BC.white),
        ],
        ids=["vacuum-white", "white-vacuum", "white-white", "partial-specular-white",
             "wrong-axis-mirror-white"],
    )
    def test_a_slab_whose_left_law_is_not_the_mirror_is_refused(self, left, right):
        """CP's slab kernel images a mirror at the left face, so any other left
        law, even one equal to the right law, would be replaced. The
        ``white-white`` row was admitted until 2026-09-30, on the false claim
        that equal laws are what CP computes (``[M]``: CP declared white|white
        1.212883, SN reflective|white 1.212884, SN white|white 1.212537).

        The ``partial-specular-white`` row is the specular wall
        ``AlbedoBoundary(0.7, SpecularReturn("x"))``: it permutes ordinates
        like the mirror and returns 0.7 of the outflow, and CP would compute a
        full mirror in its place. Until the reflective cleanup (2026-10-01)
        the row was ``ReflectiveBoundary("x", 0.7)``, "the right class, the
        wrong law", the discriminator against a guard asking the law's type;
        the mirror now has no amplitude, so that law cannot be built (its
        successor is ``tests/gates/geometry/test_reflective_is_a_mirror.py``),
        and the row now discriminates a guard asking "does it permute
        ordinates" from one asking for the mirror. The ``wrong-axis-mirror-white`` row is a perfect mirror about
        y on the slab's x-face: the right kind, the wrong motion. A guard
        asking ``kind == "reflective"`` admitted it and CP returned the
        x-mirror's k to the last bit (the elegance review, 2026-09-30; SN
        refuses the same declaration). The message names the declared law, so
        each row is refused for its own law (the guard's reason, not only its
        issue number)."""
        with pytest.raises(
            NotImplementedError,
            match=rf"mirror at the left face.*{re.escape(repr(left))}.*#513",
        ):
            CPMesh(_slab(left, right))

    @pytest.mark.rests_on(
        "tests/gates/mesh/test_dropped_laws_are_refused.py::TestCP::"
        "test_a_slab_with_the_mirror_on_its_left_builds",
    )
    @pytest.mark.parametrize("right", [BC.white, BC.vacuum], ids=["white", "vacuum"])
    def test_the_slab_kernel_drops_the_left_law(self, right, monkeypatch):
        """Why the mirror is the only honest left law, as a RECORD: with the
        guard bypassed, CP's k does not depend on the left law at all, so a
        declared vacuum, white or partially reflecting left face computes the
        mirror's answer. ``[M]`` 2026-09-30, two groups, fuel (mixture A) on
        [0, 1] cm and moderator (mixture B) on [1, 2.5] cm, 8 cells: the four
        left laws give k = 1.0402206390766764 (white right) and
        1.2301755431553738 (vacuum right) to the last bit.

        This row reddens the day CP's slab kernel reads the left law; the
        guard ``_refuse_a_law_cp_drops`` (``SCOPE-BOUNDARY[guard]``, #513) is
        then to be widened to the laws it realizes, and this row re-posed."""
        materials = {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}

        def heterogeneous_slab(left):
            return Mesher(StructuredGeometry.slab(
                (0.0, 1.0, 2.5), (0, 1), left=left, right=right,
            )).partition(
                (CellsByCount.uniform_width(4), CellsByCount.uniform_width(4)),
            ).mesh

        honest = solve_cp(materials, heterogeneous_slab(BC.reflective)).keff
        monkeypatch.setattr(cp_solver, "_refuse_a_law_cp_drops", lambda mesh: None)
        for left in (BC.vacuum, BC.white, AlbedoBoundary(0.3, SpecularReturn("x"))):
            k = solve_cp(materials, heterogeneous_slab(left)).keff
            assert k == honest, (
                f"CP's k with the left law {left!r} is {k!r}, the mirror's is "
                f"{honest!r}: the slab kernel now reads the left law, so the "
                f"guard that admits only the mirror refuses laws CP realizes "
                f"(#513)."
            )

    @pytest.mark.parametrize("right", [BC.white, BC.vacuum], ids=["white", "vacuum"])
    def test_a_slab_with_the_mirror_on_its_left_builds(self, right):
        CPMesh(_slab(BC.reflective, right))

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
        geometry = StructuredGeometry.cylinder(
            (0.0, 0.9, 1.1, 3.6 / np.sqrt(np.pi)), (2, 1, 0), outer=BC.reflective,
        )
        mesh = Mesher(geometry).partition((
            CellsByCount.uniform_volume(2), CellsByCount.uniform_volume(1),
            CellsByCount.uniform_volume(2),
        )).mesh
        MOCMesh(mesh, self._QUADRATURE)

    @pytest.mark.catches("ERR-093")
    def test_the_default_mesh_builds(self):
        """``solve_moc``'s default mesh is one MoC admits.

        Its first red: the default was the Wigner-Seitz cell, whose model law
        is white, so ``solve_moc(materials)`` raised "MOC solver does not
        support boundary condition 'white'" before a single ray was traced.
        """
        from orpheus.moc.solver import default_pin_cell_mesh

        mesh = default_pin_cell_mesh()
        assert dict(mesh.face_laws) == {"xmax": BC.reflective}
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
