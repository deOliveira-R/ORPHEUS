"""CP, MoC and MC refuse a declaration they would drop (#513, #514), as admission or as a scope guard.

Each of the three reads only part of a mesh's declared boundary: CP resolves the outer face and a mirror at a slab's left face, MC
resolves the periodic faces of a unit cell, and MoC reads any mesh as the solid Wigner–Seitz cylinder of a square pin cell. Before P1
step 2 a law they did not read was dropped silently (``[M]`` 2026-09-25, the P1 specification's probes ``cp_slab_left_law.py`` and
``partial_law_readers.py``). Step 2 made a hollow body and a declared inner law expressible on a
:class:`~orpheus.geometry.StructuredGeometry`, so each method's refusal landed with it (the specification's S3.12 to S3.14, moved into
step 2 by the user's ruling of 2026-09-29). Each refusal is a declared scope boundary.

**The migration of wave 2 (slice S4; plan sections 2.8, 4.2, 6.4, rule R8).** The law half of each guard dissolved into the method's
admission table, so what was a ``NotImplementedError`` of the guard is a ``BoundaryRefused`` of the table (CP: the ``xmin`` face of a
Cartesian site admits a ``MirrorDeck`` only, ``Refused(PLUMBING, #513)``; MC: a slab admits the translation deck only, everything else
``Refused(PLUMBING, #654)``); the hollow-body half of each guard stays a guard (``NotImplementedError``, #514 for MoC and the cavity of
MC), and these rows keep it. A declaration is a deck (a symmetry of the domain) or a tag bound to a
:class:`~orpheus.transport.boundary_conditions.WallCondition`; a face has no axis, so a mirror about the wrong axis is not spellable
(rule R5).

Dispositions of the rows of the HEAD file, each with its reason (``retirement-audit`` D.14 and D.16):

* ``wrong-axis-mirror-white`` (a perfect mirror about y on the slab's x-face): DEAD. A ``MirrorDeck`` has no axis, so the subject (a
  guard asking ``kind == "reflective"`` admitting the wrong motion) cannot be constructed. Its replacement witness is the face-local
  condition law (a condition is invariant under the carrier's colouring symmetry or ``bind`` refuses it:
  ``tests/gates/transport/test_binding.py``).
* ``test_the_slab_kernel_drops_the_left_law`` (the RECORD that CP's ``k`` does not depend on the declared left law, measured by bypassing
  ``_refuse_a_law_cp_drops``): DEAD. The guard is an admission cell of the table, and a bypass of an admission would test the
  monkeypatch and not CP. The measured numbers it recorded (``[M]`` 2026-09-30, two groups, fuel A on [0, 1] cm and moderator B on
  [1, 2.5] cm, 8 cells: ``k = 1.0402206390766764`` with a white right face and ``1.2301755431553738`` with a vacuum one, to the last bit,
  for every left law) move to the docstring of ``tests/gates/cp/test_cp_realizer_table.py::test_cp_left_face_is_a_forced_mirror``, and the
  claim it witnessed (a left face the kernel replaces by a mirror is refused) is that gate and the CP rows of the admission pin.
* the MC rows "a non-periodic left law is refused" on a slab with a periodic right face: a periodic face needs a periodic partner (the
  geometry refuses a lone ``TranslationDeck`` with ``DeckError``), so the subject splits into the geometry's pairing rule (row below) and
  the table's refusal of every non-translation face (the other rows).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.cp.solver import CPMesh
from orpheus.geometry import CoordSystem, MirrorDeck, StructuredGeometry, TranslationDeck
from orpheus.geometry.deck import DeckError
from orpheus.mc.solver import MCMesh
from orpheus.mesh import CellsByCount, Mesh1D, Mesher
from orpheus.moc.geometry import MOCMesh
from orpheus.moc.quadrature import MOCQuadrature
from orpheus.transport.boundary_conditions import BoundaryRefused, RefusalAxis, WallCondition
from tests.gates._exterior import declare_sides

pytestmark = pytest.mark.foundation

_SPECULAR = {0: WallCondition.specular(1.0)}


def _slab(left, right) -> tuple[Mesh1D, dict]:
    """A slab with its faces declared BY VALUE (rule R1: equal conditions share a tag): the mesh and its ``{tag: condition}``."""
    roles, conditions = declare_sides(left, right)
    geometry = StructuredGeometry.slab((0.0, 2.0), (0,), left=roles[0], right=roles[1])
    return Mesher(geometry).partition(CellsByCount.uniform_width(4)).mesh, conditions


def _hollow(coord: CoordSystem, inner: int = 1, outer: int = 0) -> Mesh1D:
    """A hollow body [0.5, 2.0]: the inner face is the tag ``inner`` and the outer face the tag ``outer`` (rule R6)."""
    geometry = StructuredGeometry(coord=coord, breakpoints=(0.5, 2.0), mat_ids=(0,), face_roles=(inner, outer))
    return Mesher(geometry).partition(CellsByCount.uniform_volume(4)).mesh


# ── S3.12: CP reads only the outer face, and a mirror at a slab's left ──────────────────────────────────────────────


class TestCP:
    @pytest.mark.parametrize(
        "left, right",
        [
            (WallCondition.vacuum(), WallCondition.white()), (WallCondition.white(), WallCondition.vacuum()),
            (WallCondition.white(), WallCondition.white()), (WallCondition.specular(0.7), WallCondition.white()),
        ],
        ids=["vacuum-white", "white-vacuum", "white-white", "partial-specular-white"],
    )
    def test_a_slab_whose_left_face_is_not_the_mirror_is_refused(self, left, right):
        """CP's slab kernel images a mirror at the left face, so any other left condition, even one equal to the right one, would be
        replaced. The ``white-white`` row was admitted until 2026-09-30, on the false claim that equal laws are what CP computes
        (``[M]``: CP declared white|white 1.212883, SN reflective|white 1.212884, SN white|white 1.212537); here it is ONE tag bound to
        ``white`` on both faces (rule R1), refused at ``xmin``, the first face.

        The ``partial-specular-white`` row is the specular wall ``specular(0.7)``: it permutes ordinates like the mirror and returns 0.7
        of the outflow, and CP would compute a full mirror in its place; the row discriminates a table asking "does it permute
        ordinates" from one asking for the mirror. The refusal is the table's: ``BoundaryRefused`` whose admission is
        ``Refused(PLUMBING, #513)`` at the face ``xmin`` (each row is refused for its own face and issue, not only for its type)."""
        mesh, conditions = _slab(left, right)
        with pytest.raises(BoundaryRefused) as info:
            CPMesh(mesh, boundary_conditions=conditions)
        assert info.value.face == "xmin"
        assert info.value.admission.axis is RefusalAxis.PLUMBING and info.value.admission.issue == 513

    @pytest.mark.parametrize("right", [WallCondition.white(), WallCondition.vacuum()], ids=["white", "vacuum"])
    def test_a_slab_with_the_mirror_on_its_left_builds(self, right):
        mesh, conditions = _slab(MirrorDeck(), right)
        CPMesh(mesh, boundary_conditions=conditions)

    @pytest.mark.parametrize(
        "coord", [CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL], ids=lambda c: c.name.lower(),
    )
    def test_a_hollow_body_with_an_inner_condition_is_refused(self, coord):
        """The inner face of a hollow body carries a tag whose condition CP would drop: ``BoundaryRefused``, ``Refused(PLUMBING, #513)``
        at ``xmin`` (today's refusal, re-keyed by face; the same cell is pinned by ``test_cp_realizer_table``)."""
        with pytest.raises(BoundaryRefused) as info:
            CPMesh(_hollow(coord), boundary_conditions={1: WallCondition.white(), 0: WallCondition.white()})
        assert info.value.face == "xmin" and info.value.admission.issue == 513


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


# ── S3.14: MC reads the translation of a unit cell ──────────────────────────────────────────────────────────────────


class TestMC:
    @pytest.mark.parametrize(
        "faces", [(WallCondition.vacuum(), WallCondition.vacuum()), (MirrorDeck(), MirrorDeck())], ids=["vacuum", "mirror"],
    )
    def test_a_slab_whose_faces_are_not_the_translation_is_refused(self, faces):
        """MC applies the periodic wrap to every face: any other face (a wall, a mirror) is dropped. After the carve the table admits
        the translation deck on a slab and states everything else as ``Refused(PLUMBING, #654)`` (the registry held nothing else)."""
        mesh, conditions = _slab(*faces)
        with pytest.raises(BoundaryRefused) as info:
            MCMesh(mesh, pitch=2.0, boundary_conditions=conditions)
        assert info.value.admission.axis is RefusalAxis.PLUMBING and info.value.admission.issue == 654

    def test_a_lone_translation_face_is_not_a_slab(self):
        """The old row "a periodic right face with a non-periodic left law" cannot be declared: a translation deck needs its partner (the
        geometry refuses ``TranslationDeck`` against a tag, ``DeckError``), so MC never sees the dropped left law the guard was written for."""
        with pytest.raises(DeckError):
            StructuredGeometry.slab((0.0, 2.0), (0,), left=0, right=TranslationDeck())

    def test_a_translation_pair_builds(self):
        mesh, conditions = _slab(TranslationDeck(), TranslationDeck())
        MCMesh(mesh, pitch=2.0, boundary_conditions=conditions)

    def test_a_hollow_cylinder_is_refused(self):
        """Its lookup would fill the cavity with the innermost material (the scope guard, #513, is read before the declaration)."""
        mesh = _hollow(CoordSystem.CYLINDRICAL)
        with pytest.raises(NotImplementedError, match="cavity"):
            MCMesh(mesh, pitch=4.0, boundary_conditions={1: WallCondition.specular(1.0), 0: WallCondition.specular(1.0)})
