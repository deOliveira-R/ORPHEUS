r"""The bare 1-D mesh: its construction laws, its equality, its retirements (P1 step 3b).

``Mesh1D(coord, edges, volumes, region_ids, region_materials, face_roles)`` is what a 1-D
discretisation IS (the user's ruling of 2026-09-29, option A: the mesh holds
no geometry; it knows its cells, each cell's material, and each boundary
face's ROLE: a deck or an ``int`` tag, wave 2 S2b re-pose of the face laws). The mesher (:mod:`orpheus.mesh.mesher`) is the one call site that
lifts a geometry onto it; this file gates the constructor itself, its derived
views and its immutability (``TestTheValue``), and the constructor's coercion
(``TestMesh1D``, which moved here from ``tests/gates/geometry/test_geometry.py``
in the layer-DAG re-filing, `.claude/plans/layer_dag.md`, ruling 21). The five
rows ``TestMesh1D`` also held (widths, centres, total width, cell count,
immutability, through a mesher-built mesh) were deleted at ruling 30
(2026-10-26): ``TestTheValue::test_the_derived_views`` and
``::test_frozen_and_read_only`` reddened with each under its one mutation,
and the mesher route they added is ``test_mesher.py``'s claim.

Gate ids (``.claude/plans/reference_p1_spec.md`` §1.3, §1.3a): S3.8 (the
constructor's laws, re-posed on the bare constructor), S3.10 (the retired
fields and ``None``), S3.11 (the refusal rows ``TestMesh1D`` held while it
lived in ``geometry/test_geometry.py``, re-posed here). Every test is
``foundation``; claim kind THEOREM unless marked RECORD.

The volume law. Each stored volume must be the coordinate system's measure of
its cell, ``c (T(r_{j+1}) - T(r_j))``, within ``_volume_ulps(coord)`` ulp of
``c T(r_{j+1})``: an equal share ``m/n`` and the shell between the realised
edges differ by the rounding of ``T``, ``T^-1`` and the share. The control is a
volume from the WRONG coordinate system, about 1e15 ulp off.
"""
from __future__ import annotations

import math
import re
from collections.abc import Callable

import numpy as np
import pytest

import orpheus.mesh
import orpheus.mesh.structured as structured_module
from orpheus.geometry import CoordSystem, compute_areas_1d
from orpheus.geometry.deck import MirrorDeck
from orpheus.mesh import Mesh1D
from tests.gates.geometry._measure import EXPONENT as _EXPONENT

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/mesh/test_mesh1d.py"
_ONE_MEASURE = (
    "tests/gates/geometry/test_coord.py::TestTheOneMeasure::"
    "test_the_measure_is_c_times_the_difference_of_T"
)

_SLAB = CoordSystem.CARTESIAN
_CYLINDER = CoordSystem.CYLINDRICAL
_SPHERE = CoordSystem.SPHERICAL


def _mesh(coord=_SLAB, edges=(0.0, 0.5, 2.0), volumes=None, region_ids=None, region_materials=None, face_roles=None) -> Mesh1D:
    """A valid mesh, with the named fields replaced; volumes default to the
    coordinate system's measure of the cells."""
    e = np.asarray(edges, dtype=float)
    if volumes is None:
        volumes = coord.measure(e) if isinstance(coord, CoordSystem) else np.diff(e)
    if region_ids is None:
        region_ids = np.zeros(len(e) - 1, dtype=int)
    if region_materials is None:
        region_materials = tuple(range(int(np.max(region_ids)) + 1)) if len(region_ids) else (0,)
    if face_roles is None:
        # the vacuum tag is 0; a reflecting face is a mirror deck on a slab and a tag (1) on a curved body (rule R2)
        solid = isinstance(coord, CoordSystem) and coord is not _SLAB and e[0] == 0.0
        reflective = MirrorDeck() if coord is _SLAB else 1
        face_roles = {"xmax": 0} if solid else {"xmin": reflective, "xmax": 0}
    return Mesh1D(coord=coord, edges=edges, volumes=volumes, region_ids=region_ids, region_materials=region_materials,
                  face_roles=face_roles)


_ONE_ULP_ABOVE_ONE = float(np.nextafter(1.0, 2.0))


def _off_by(coord: CoordSystem, edges, j: int, k: float) -> np.ndarray:
    """The measure of the cells with volume ``j`` moved by ``k`` ulp of ``c T(r_{j+1})``."""
    e = np.asarray(edges, dtype=float)
    v = coord.measure(e).copy()
    unit = float(np.spacing(coord.measure_constant * e[j + 1] ** _EXPONENT[coord]))
    v[j] += k * unit
    return v


# ─────────────────────────────────────────────────────────────────────
# S3.8 / S3.11 — the construction laws, keyed and disjoint
# ─────────────────────────────────────────────────────────────────────

#: (case id, a callable building the mesh, error type, fragment).
_REFUSALS: list[tuple[str, Callable[[], object], type[Exception], str]] = [
    ("coord-string", lambda: _mesh(coord="SLB"), TypeError,  # type: ignore[arg-type]  # a refusal input
     "Mesh1D.coord is a CoordSystem member"),
    ("edge-string", lambda: _mesh(edges=("0", 1.0), volumes=[1.0], region_ids=[0]),
     TypeError, "Mesh1D.edges[0] must be a real number"),
    ("edge-bool", lambda: _mesh(edges=(False, True), volumes=[1.0], region_ids=[0]),
     TypeError, "Mesh1D.edges[0] must be a real number"),
    ("edge-nan", lambda: _mesh(edges=(0.0, math.nan, 2.0), volumes=[1.0, 1.0]),
     ValueError, "Mesh1D.edges[1] is NaN, which is not a number"),  # NaN: parse_real's own refusal since #405 P1 step 5
    # S3.11, re-posed: the non-monotone, equal-edge and too-few-edge rows
    ("edges-decreasing", lambda: _mesh(edges=(0.0, 2.0, 1.0), volumes=[2.0, 1.0]),
     ValueError, "are at least two strictly increasing positions"),
    ("edges-equal", lambda: _mesh(edges=(0.0, 1.0, 1.0), volumes=[1.0, 1.0]),
     ValueError, "are at least two strictly increasing positions"),
    ("one-edge", lambda: _mesh(edges=(0.0,), volumes=[], region_ids=[]),
     ValueError, "are at least two strictly increasing positions"),
    ("negative-radius", lambda: _mesh(coord=_SPHERE, edges=(-0.5, 1.0), volumes=[1.0]),
     ValueError, "a radial coordinate starts at r_0 >= 0"),
    ("volume-count", lambda: _mesh(volumes=[0.5]),
     ValueError, "Mesh1D.volumes: 2 cell(s) need 2 volume(s)"),
    ("volume-not-the-measure", lambda: _mesh(volumes=[0.5, 1.4]),
     ValueError, "is not the cartesian measure of the cell"),
    ("volume-wrong-coordinate", lambda: _mesh(coord=_SPHERE, volumes=_CYLINDER.measure(np.array([0.0, 0.5, 2.0]))),
     ValueError, "is not the spherical measure of the cell"),
    ("volume-non-finite", lambda: _mesh(volumes=[0.5, math.inf]),
     ValueError, "Mesh1D.volumes[1] is infinite"),  # the shared finite parse since #559
    # qa F4: on a cell one ulp wide, a zero or negative volume is within the
    # absolute band of its measure, so positivity is its own law
    ("volume-zero-on-a-thin-cell", lambda: _mesh(edges=(1.0, _ONE_ULP_ABOVE_ONE), volumes=[0.0], region_ids=[0]),
     ValueError, "a cell volume is positive"),
    ("volume-negative-on-a-thin-cell", lambda: _mesh(edges=(1.0, _ONE_ULP_ABOVE_ONE), volumes=[-1e-300], region_ids=[0]),
     ValueError, "a cell volume is positive"),
    # S3.11's material-id count and type rows moved to R7b2.0.2 (#405 P2 step 7b.2.0): the
    # mesh's labels and its region -> material map carry them (test_mesh1d_regions.py)
    # the face roles are a mapping over the face inventory (step 3c's FaceLaws, renamed FaceRoles at wave 2)
    ("slab-one-role", lambda: _mesh(face_roles={"xmax": 0}),
     ValueError, "Mesh1D.face_roles: this cartesian mesh has the boundary faces ('xmin', 'xmax')"),
    ("role-at-the-centre", lambda: _mesh(coord=_SPHERE, face_roles={"xmin": 0, "xmax": 0}),
     ValueError, "the centre r = 0 of a solid spherical body is an interior point"),
    ("hollow-one-role", lambda: _mesh(coord=_CYLINDER, edges=(0.5, 1.0, 2.0), face_roles={"xmax": 0}),
     ValueError, "a hollow cylindrical body has an inner surface"),
    # step 3c: the positional tuple is retired
    ("tuple-roles", lambda: _mesh(face_roles=(MirrorDeck(), 0)),
     TypeError, "Mesh1D.face_roles is a mapping from face name to role"),
    # S3.10: None is not a role
    ("none-role", lambda: _mesh(face_roles={"xmin": None, "xmax": 0}), TypeError,
     "a face role is a MirrorDeck, a TranslationDeck or an int tag, got None"),
    # S3.11, re-posed: the old ``bc_left`` type row
    ("string-role", lambda: _mesh(face_roles={"xmin": "vacuum", "xmax": 0}),
     TypeError, "a face role is a MirrorDeck, a TranslationDeck or an int tag, got 'vacuum'"),
    # wave 2: a reflecting curved wall is a response, never a deck (rule R2)
    ("deck-on-a-curved-wall", lambda: _mesh(coord=_SPHERE, face_roles={"xmax": MirrorDeck()}),
     ValueError, "a reflecting curved wall is a response; bind a specular exterior condition to a tag"),
]


class TestConstructionLaws:
    r"""S3.8's refusals. Each row names the clause it trips; one row asserts
    that no fragment appears in another row's message."""

    @pytest.mark.parametrize(
        "build, error, fragment",
        [pytest.param(b, e, f, id=c) for c, b, e, f in _REFUSALS],
    )
    def test_refusal(self, build, error, fragment):
        with pytest.raises(error, match=re.escape(fragment)):
            build()

    @pytest.mark.rests_on(f"{_HERE}::TestConstructionLaws::test_refusal")
    def test_refusal_fragments_are_disjoint(self):
        messages = {}
        for case, build, error, _ in _REFUSALS:
            with pytest.raises(error) as caught:
                build()
            messages[case] = str(caught.value)
        for case, *_, fragment in _REFUSALS:
            for other, message in messages.items():
                other_fragment = next(f for c, *_, f in _REFUSALS if c == other)
                if other_fragment == fragment or fragment in other_fragment:
                    continue
                assert fragment not in message, (case, other, message)

    @pytest.mark.parametrize(
        "coord, edges, faces",
        [
            pytest.param(_SLAB, (-1.0, 0.3, 2.0), ("xmin", "xmax"), id="slab"),
            pytest.param(_CYLINDER, (0.0, 0.3, 2.0), ("xmax",), id="solid-cylinder"),
            pytest.param(_SPHERE, (0.5, 0.9, 2.0), ("xmin", "xmax"), id="hollow-sphere"),
        ],
    )
    def test_the_positive_legs(self, coord, edges, faces):
        """A mesh at every legal boundary shape builds, its roles one per face,
        the faces as written here by hand from the topology law, and the
        boundary points paired with them in order."""
        mesh = _mesh(coord=coord, edges=edges)
        assert tuple(mesh.face_roles) == faces
        assert mesh.boundary_points == coord.boundary_points(edges[0], edges[-1])
        assert len(mesh.boundary_points) == len(faces)
        assert mesh.outer_role is mesh.face_roles["xmax"]


class TestTheVolumeLaw:
    r"""S3.8 (c) re-posed: each stored volume is the coordinate system's
    measure of its cell within ``_volume_ulps(coord) = 2p + 5`` ulp of
    ``c T(r_{j+1})`` (``p`` the measure exponent: 7 on a slab, 9 on a
    cylinder, 11 on a sphere).

    Claim kind: THEOREM with a DERIVED tolerance (the rounding count in
    ``_volume_ulps``'s docstring: ``p + 1.5`` ulp per edge, two edges, about
    2 for the subtraction, the constant and the share). Checked two ways:
    the band is two-edged at ``band -+ 1`` ulp (the addition of ``k`` ulp to
    a volume rounds by at most half an ulp, so ``band - 1`` is inside and
    ``band + 1`` outside); and the worst LEGAL producer output found by
    search is admitted. That worst case is the RECORD (libm-dependent:
    measured on macOS; glibc's ``cbrt`` is specified to 1 ulp only).
    The wrong-coordinate control sits about 1e15 ulp off.
    """

    @pytest.mark.rests_on(_ONE_MEASURE)
    @pytest.mark.parametrize("coord", [_SLAB, _CYLINDER, _SPHERE], ids=lambda c: c.name.lower())
    def test_the_band(self, coord):
        band = structured_module._volume_ulps(coord)
        assert band == 2 * _EXPONENT[coord] + 5
        edges = (0.3, 0.9, 2.0)
        _mesh(coord=coord, edges=edges, volumes=_off_by(coord, edges, 1, band - 1.0))
        _mesh(coord=coord, edges=edges, volumes=_off_by(coord, edges, 1, -(band - 1.0)))
        with pytest.raises(ValueError, match="measure of the cell"):
            _mesh(coord=coord, edges=edges, volumes=_off_by(coord, edges, 1, band + 1.0))
        with pytest.raises(ValueError, match=re.escape("Mesh1D.volumes[0]")):
            _mesh(coord=coord, edges=edges, volumes=_off_by(coord, edges, 0, -(band + 1.0)))

    @pytest.mark.parametrize(
        "coord, a, b, n, measured",
        [
            # [M] 2026-09-29, macOS libm: the elegance review's sphere case,
            # and the worst of 4000 random equal-volume intervals per
            # coordinate (scales 1e-3 to 1e3, n < 3000; the test-architect's
            # scratchpad probe worst.py).
            pytest.param(CoordSystem.SPHERICAL, 0.06675753124937545, 0.9507586186918889, 2400, 7.50390625,
                         id="sphere-elegance-review"),
            pytest.param(CoordSystem.SPHERICAL, 43.0446014117623, 81.72750086937991, 2677, 7.4599609375,
                         id="sphere-search"),
            pytest.param(CoordSystem.CYLINDRICAL, 5.8405877264467225, 23.15910405796302, 997, 6.623046875,
                         id="cylinder-search"),
            pytest.param(CoordSystem.CARTESIAN, 7.132951756503764e-05, 0.0027791855266102957, 2827, 3.08203125,
                         id="slab-search"),
        ],
    )
    def test_the_worst_legal_output_is_admitted(self, coord, a, b, n, measured):
        """The band's lower edge against the producers: the equal-volume
        partition with the largest measured gap builds (the type must not
        refuse correct output). The gap is re-measured here and must stay
        within 1 ulp of the recorded one (a libm with another ``cbrt`` can
        move it), or the row has lost its witness. The slab's worst, 3.08
        ulp, sits well inside its band of 7: the derivation is conservative
        there (no root, one rounding of the share)."""
        from orpheus.geometry import StructuredGeometry
        from orpheus.mesh import CellsByCount, Mesher

        geometry = StructuredGeometry.uniform_boundary(coord, (a, b), (0,), 0)
        mesh = Mesher(geometry).partition(CellsByCount.uniform_volume(n)).mesh
        unit = np.spacing(coord.measure_constant * coord.measure_coordinate(mesh.edges[1:]))
        gap = float(np.max(np.abs(mesh.volumes - coord.measure(mesh.edges)) / unit))
        band = structured_module._volume_ulps(coord)
        assert gap <= band
        assert gap >= measured - 1.0, f"the recorded worst case ({measured} ulp) now reads {gap} ulp"

    @pytest.mark.parametrize(
        "coord, wrong",
        [(_SLAB, _SPHERE), (_CYLINDER, _SPHERE), (_SPHERE, _CYLINDER), (_CYLINDER, _SLAB)],
        ids=lambda c: c.name.lower(),
    )
    def test_the_wrong_coordinate_is_refused(self, coord, wrong):
        """The control: volumes measured in another coordinate system."""
        edges = np.array([0.0, 0.5, 2.0])
        with pytest.raises(ValueError, match=f"is not the {coord.name.lower()} measure"):
            _mesh(coord=coord, edges=edges, volumes=wrong.measure(edges))

    def test_the_volumes_are_stored_not_recomputed(self):
        """An in-band equal share is stored as given, never replaced by the
        shell of the edges (ERR-020's invariant at the constructor)."""
        edges = (0.0, 0.5, 2.0)
        given = _off_by(_SPHERE, edges, 1, 3.0)
        mesh = _mesh(coord=_SPHERE, edges=edges, volumes=given)
        np.testing.assert_array_equal(mesh.volumes, given)
        assert not np.array_equal(mesh.volumes, _SPHERE.measure(np.asarray(edges)))


class TestTheValue:
    """The mesh is a frozen value with bitwise equality and no hash yet."""

    def test_equality_is_bitwise(self):
        a = _mesh(coord=_SPHERE, edges=(0.0, 0.5, 2.0))
        b = _mesh(coord=_SPHERE, edges=[0.0, 0.5, 2.0])
        assert a == b
        assert a != _mesh(coord=_SPHERE, volumes=_off_by(_SPHERE, (0.0, 0.5, 2.0), 1, 1.0))
        assert a != _mesh(coord=_SPHERE, region_ids=[0, 1])
        assert a != _mesh(coord=_SPHERE, face_roles={"xmax": 1})
        assert a != _mesh(coord=_SLAB, edges=(0.0, 0.5, 2.0))  # another coordinate system
        assert (a == (0.0, 0.5, 2.0)) is False

    def test_the_eq_hash_contract(self):
        """Content identity (#405 P1 step 5, 2026-10-02): two meshes built
        independently with equal cells and laws are equal and hash alike.
        Until step 5 this row was ``test_unhashable_until_step_5``, the
        RECORD that ``hash(mesh)`` raised, written to red when the step
        landed."""
        a, b = _mesh(), _mesh()
        assert a is not b
        assert a == b and hash(a) == hash(b)
        assert len({a, b}) == 1
        # The separating leg: an ``__eq__`` returning True everywhere would
        # pass the three lines above.
        other = _mesh(region_ids=[0, 1])
        assert other != a and len({a, other}) == 2

    def test_frozen_and_read_only(self):
        mesh = _mesh()
        with pytest.raises(AttributeError):
            mesh.edges = np.array([0.0, 1.0])  # type: ignore[misc]  # the frozen field is the subject
        for array in (mesh.edges, mesh.volumes, mesh.region_ids, mesh.mat_ids):
            assert not array.flags.writeable

    def test_the_inputs_are_copied(self):
        edges, volumes = np.array([0.0, 0.5, 2.0]), np.array([0.5, 1.5])
        mesh = _mesh(edges=edges, volumes=volumes)
        edges[1], volumes[0] = 0.25, 9.0
        np.testing.assert_array_equal(mesh.edges, [0.0, 0.5, 2.0])
        np.testing.assert_array_equal(mesh.volumes, [0.5, 1.5])

    @pytest.mark.parametrize(
        "edges",
        # A hollow mesh, and one starting at the origin: the TestMesh1D view rows deleted as this row's duplicates
        # (layer-DAG ruling 30) built origin meshes, and a slip wrong only when edges[0] == 0 reddens just the second
        # (`[M]` 2026-10-26, qa's centers_origin and widths_origin arms).
        [(0.5, 0.9, 2.0), (0.0, 1.0, 3.0, 6.0)],
        ids=("hollow", "from-the-origin"),
    )
    @pytest.mark.parametrize("coord", [_SLAB, _CYLINDER, _SPHERE], ids=lambda c: c.name.lower())
    def test_the_derived_views(self, coord, edges):
        """S3.8 (e): widths, centres and face areas from the edges."""
        edges = np.array(edges)
        mesh = _mesh(coord=coord, edges=edges)
        np.testing.assert_array_equal(mesh.widths, np.diff(edges))
        np.testing.assert_array_equal(mesh.centers, 0.5 * (edges[:-1] + edges[1:]))
        np.testing.assert_array_equal(mesh.areas, compute_areas_1d(coord, edges))
        assert mesh.N == len(edges) - 1 and mesh.total_width == edges[-1] - edges[0]

    def test_with_distinct_cell_ids(self):
        role = 3
        mesh = _mesh(coord=_SPHERE, edges=(0.5, 0.9, 2.0), region_ids=[0, 0], region_materials=(4,), face_roles={"xmin": 1, "xmax": role})
        relabelled = mesh.with_distinct_cell_ids()
        np.testing.assert_array_equal(relabelled.mat_ids, [0, 1])
        np.testing.assert_array_equal(relabelled.edges, mesh.edges)
        np.testing.assert_array_equal(relabelled.volumes, mesh.volumes)
        assert relabelled.face_roles["xmax"] == role


# ─────────────────────────────────────────────────────────────────────
# S3.10 — the retirements
# ─────────────────────────────────────────────────────────────────────


class TestTheRetirements:
    """S3.10: the retired fields, constructor and descriptor are gone (the
    witness of the retirement, ``retirement-audit`` D.18). ``None`` as a face
    law is refused by ``TestConstructionLaws`` (row ``none-role``); the SN
    ``boundary_condition=`` leg and the resolver's default are step 3c's: the
    two SN rows in ``tests/gates/sn/solve/test_no_entry_takes_a_boundary_condition.py``
    and the realization row in
    ``tests/gates/transport/test_realize_boundary_fills_no_default.py``."""

    @pytest.mark.parametrize(
        "name", ["bc_left", "bc_right", "precomputed_volumes", "from_geometry", "boundary_faces", "face_laws", "outer_law"],
    )
    def test_a_retired_mesh_attribute_is_gone(self, name):
        assert not hasattr(Mesh1D, name)
        assert not hasattr(_mesh(), name)

    @pytest.mark.parametrize("name", ["RegionMesh"])
    def test_a_retired_mesh_name_is_gone(self, name):
        assert not hasattr(orpheus.mesh, name)
        assert not hasattr(structured_module, name)

    def test_the_predecessor_subdivision_is_gone(self):
        import orpheus.mesh.factories as factories

        assert not hasattr(factories, "_subdivide_zone")

    def test_the_constructor_takes_exactly_the_discretisation(self):
        import dataclasses

        init_fields = [f.name for f in dataclasses.fields(Mesh1D) if f.init]
        # #405 P2 step 7b.2.0: the region labels and the region -> material map replace ``mat_ids``,
        # which is derived (R7b2.0.1 in ``test_mesh1d_regions.py``).
        assert init_fields == ["coord", "edges", "volumes", "region_ids", "region_materials", "face_roles"]


# ─────────────────────────────────────────────────────────────────────
# S3.10, step 3c's leg — no boundary law is filled in for the caller
# ─────────────────────────────────────────────────────────────────────


# The route gate that stood here (``TestOneBoundaryLawCheck``: the geometry and
# the mesh parse through one ``parse_boundary_laws``) was re-posed at step 3c,
# when the mesh's laws became a ``FaceLaws`` over ``face_inventory``: its
# successors are ``tests/gates/mesh/test_mesh2d_face_roles.py::TestTheOneElementParser``
# (every declared law of the five owners passes ``parse_boundary_law``) and
# ``::TestTheOneInventoryRule`` (both meshes derive their faces from
# ``face_inventory``).


# ═══════════════════════════════════════════════════════════════════════
# Mesh1D properties and validation
# ═══════════════════════════════════════════════════════════════════════

class TestMesh1D:
    """The constructor's coercion of its inputs. The construction refusals
    (edges, volumes, material ids, face roles) are the keyed table of this file
    (S3.8 and S3.11, re-posed at P1 step 3b from the rows this class held while
    it lived in ``tests/gates/geometry/test_geometry.py``); the derived views and
    immutability are ``TestTheValue``'s."""

    def test_coerces_to_float_and_int(self):
        """The bare constructor accepts lists; edges become float, region labels and mat_ids int."""
        mesh = Mesh1D(
            coord=CoordSystem.CARTESIAN, edges=[0, 1, 2], volumes=[1, 1],  # type: ignore[arg-type]  # the coercion is the subject
            region_ids=[0, 1], region_materials=(0, 1),  # type: ignore[arg-type]  # the coercion is the subject
            face_roles={"xmin": MirrorDeck(), "xmax": MirrorDeck()},
        )
        assert mesh.edges.dtype == float
        assert mesh.region_ids.dtype == int
        assert mesh.mat_ids.dtype == int
