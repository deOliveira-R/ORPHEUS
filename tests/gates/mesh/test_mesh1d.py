r"""The bare 1-D mesh: its construction laws, its equality, its retirements (P1 step 3b).

``Mesh1D(coord, edges, volumes, mat_ids, face_laws)`` is what a 1-D
discretisation IS (the user's ruling of 2026-09-29, option A: the mesh holds
no geometry; it knows its cells, each cell's material, and each boundary
face's law). The mesher (:mod:`orpheus.mesh.mesher`) is the one call site that
lifts a geometry onto it; this file gates the constructor itself.

Gate ids (``.claude/plans/reference_p1_spec.md`` §1.3, §1.3a): S3.8 (the
constructor's laws, re-posed on the bare constructor), S3.10 (the retired
fields and ``None``), S3.11 (the re-posed refusals of the old
``tests/gates/geometry/test_geometry.py::TestMesh1D`` rows). Every test is
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
from orpheus.geometry import BC, CoordSystem, compute_areas_1d
from orpheus.mesh import Mesh1D

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/mesh/test_mesh1d.py"
_ONE_MEASURE = (
    "tests/gates/mesh/test_partition.py::TestTheOneMeasure::"
    "test_the_measure_is_c_times_the_difference_of_T"
)

_SLAB = CoordSystem.CARTESIAN
_CYLINDER = CoordSystem.CYLINDRICAL
_SPHERE = CoordSystem.SPHERICAL
_EXPONENT = {_SLAB: 1, _CYLINDER: 2, _SPHERE: 3}


def _mesh(coord=_SLAB, edges=(0.0, 0.5, 2.0), volumes=None, mat_ids=None, face_laws=None) -> Mesh1D:
    """A valid mesh, with the named fields replaced; volumes default to the
    coordinate system's measure of the cells."""
    e = np.asarray(edges, dtype=float)
    if volumes is None:
        volumes = coord.measure(e) if isinstance(coord, CoordSystem) else np.diff(e)
    if mat_ids is None:
        mat_ids = np.zeros(len(e) - 1, dtype=int)
    if face_laws is None:
        solid = isinstance(coord, CoordSystem) and coord is not _SLAB and e[0] == 0.0
        face_laws = (BC.vacuum,) if solid else (BC.reflective, BC.vacuum)
    return Mesh1D(coord=coord, edges=edges, volumes=volumes, mat_ids=mat_ids, face_laws=face_laws)


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
    ("edge-string", lambda: _mesh(edges=("0", 1.0), volumes=[1.0], mat_ids=[0]),
     TypeError, "Mesh1D.edges[0] must be a real number"),
    ("edge-bool", lambda: _mesh(edges=(False, True), volumes=[1.0], mat_ids=[0]),
     TypeError, "Mesh1D.edges[0] must be a real number"),
    ("edge-nan", lambda: _mesh(edges=(0.0, math.nan, 2.0), volumes=[1.0, 1.0]),
     ValueError, "Mesh1D.edges must be finite"),
    # S3.11, re-posed: the non-monotone, equal-edge and too-few-edge rows
    ("edges-decreasing", lambda: _mesh(edges=(0.0, 2.0, 1.0), volumes=[2.0, 1.0]),
     ValueError, "are at least two strictly increasing positions"),
    ("edges-equal", lambda: _mesh(edges=(0.0, 1.0, 1.0), volumes=[1.0, 1.0]),
     ValueError, "are at least two strictly increasing positions"),
    ("one-edge", lambda: _mesh(edges=(0.0,), volumes=[], mat_ids=[]),
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
     ValueError, "Mesh1D.volumes must be finite"),
    # qa F4: on a cell one ulp wide, a zero or negative volume is within the
    # absolute band of its measure, so positivity is its own law
    ("volume-zero-on-a-thin-cell", lambda: _mesh(edges=(1.0, _ONE_ULP_ABOVE_ONE), volumes=[0.0], mat_ids=[0]),
     ValueError, "a cell volume is positive"),
    ("volume-negative-on-a-thin-cell", lambda: _mesh(edges=(1.0, _ONE_ULP_ABOVE_ONE), volumes=[-1e-300], mat_ids=[0]),
     ValueError, "a cell volume is positive"),
    # S3.11, re-posed: the wrong material-id count is a law of the bare
    # constructor again (option A: the mesh knows its cells' materials)
    ("mat-id-count", lambda: _mesh(mat_ids=[0, 1, 2]),
     ValueError, "Mesh1D.mat_ids: 2 cell(s) need 2 material id(s)"),
    ("mat-id-float", lambda: _mesh(mat_ids=[0, 1.0]), TypeError, "a material id is an int"),
    ("slab-one-law", lambda: _mesh(face_laws=(BC.vacuum,)),
     ValueError, "Mesh1D.face_laws: a slab has two boundary points"),
    ("law-at-the-centre", lambda: _mesh(coord=_SPHERE, face_laws=(BC.vacuum, BC.vacuum)),
     ValueError, "Mesh1D.face_laws: the centre r = 0 of a solid spherical body"),
    ("hollow-one-law", lambda: _mesh(coord=_CYLINDER, edges=(0.5, 1.0, 2.0), face_laws=(BC.vacuum,)),
     ValueError, "Mesh1D.face_laws: a hollow cylindrical body"),
    # S3.10: None is not a law
    ("none-law", lambda: _mesh(face_laws=(None, BC.vacuum)), TypeError, "None is not a boundary law"),
    # S3.11, re-posed: the old ``bc_left`` type row
    ("string-law", lambda: _mesh(face_laws=("vacuum", BC.vacuum)),
     TypeError, "must be a BC tag or a BoundaryTraceLaw instance"),
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
        "coord, edges",
        [
            pytest.param(_SLAB, (-1.0, 0.3, 2.0), id="slab"),
            pytest.param(_CYLINDER, (0.0, 0.3, 2.0), id="solid-cylinder"),
            pytest.param(_SPHERE, (0.5, 0.9, 2.0), id="hollow-sphere"),
        ],
    )
    def test_the_positive_legs(self, coord, edges):
        """A mesh at every legal boundary shape builds, its laws one per face."""
        mesh = _mesh(coord=coord, edges=edges)
        assert len(mesh.face_laws) == len(mesh.boundary_faces)
        assert mesh.boundary_faces == coord.boundary_points(edges[0], edges[-1])
        assert mesh.outer_law is mesh.face_laws[-1]


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

        geometry = StructuredGeometry.uniform_boundary(coord, (a, b), (0,), BC.vacuum)
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
        assert a != _mesh(coord=_SPHERE, mat_ids=[0, 1])
        assert a != _mesh(coord=_SPHERE, face_laws=(BC.reflective,))
        assert a != _mesh(coord=_SLAB, edges=(0.0, 0.5, 2.0))  # another coordinate system
        assert (a == (0.0, 0.5, 2.0)) is False

    def test_unhashable_until_step_5(self):
        """RECORD: content identity is P1 step 5's. When it lands this row
        reds on purpose; re-pose it as the eq/hash contract."""
        with pytest.raises(TypeError, match="unhashable"):
            hash(_mesh())

    def test_frozen_and_read_only(self):
        mesh = _mesh()
        with pytest.raises(AttributeError):
            mesh.edges = np.array([0.0, 1.0])  # type: ignore[misc]  # the frozen field is the subject
        for array in (mesh.edges, mesh.volumes, mesh.mat_ids):
            assert not array.flags.writeable

    def test_the_inputs_are_copied(self):
        edges, volumes = np.array([0.0, 0.5, 2.0]), np.array([0.5, 1.5])
        mesh = _mesh(edges=edges, volumes=volumes)
        edges[1], volumes[0] = 0.25, 9.0
        np.testing.assert_array_equal(mesh.edges, [0.0, 0.5, 2.0])
        np.testing.assert_array_equal(mesh.volumes, [0.5, 1.5])

    @pytest.mark.parametrize("coord", [_SLAB, _CYLINDER, _SPHERE], ids=lambda c: c.name.lower())
    def test_the_derived_views(self, coord):
        """S3.8 (e): widths, centres and face areas from the edges."""
        edges = np.array([0.5, 0.9, 2.0])
        mesh = _mesh(coord=coord, edges=edges)
        np.testing.assert_array_equal(mesh.widths, np.diff(edges))
        np.testing.assert_array_equal(mesh.centers, 0.5 * (edges[:-1] + edges[1:]))
        np.testing.assert_array_equal(mesh.areas, compute_areas_1d(coord, edges))
        assert mesh.N == 2 and mesh.total_width == 1.5

    def test_with_distinct_cell_ids(self):
        law = BC("albedo", {"albedo": 0.3})
        mesh = _mesh(coord=_SPHERE, edges=(0.5, 0.9, 2.0), mat_ids=[4, 4], face_laws=(BC.reflective, law))
        relabelled = mesh.with_distinct_cell_ids()
        np.testing.assert_array_equal(relabelled.mat_ids, [0, 1])
        np.testing.assert_array_equal(relabelled.edges, mesh.edges)
        np.testing.assert_array_equal(relabelled.volumes, mesh.volumes)
        assert relabelled.face_laws[1] is law


# ─────────────────────────────────────────────────────────────────────
# S3.10 — the retirements
# ─────────────────────────────────────────────────────────────────────


class TestTheRetirements:
    """S3.10: the retired fields, constructor and descriptor are gone (the
    witness of the retirement, ``retirement-audit`` D.18). ``None`` as a face
    law is refused by ``TestConstructionLaws`` (row ``none-law``); the SN
    ``boundary_condition=`` leg is step 3c's."""

    @pytest.mark.parametrize("name", ["bc_left", "bc_right", "precomputed_volumes", "from_geometry"])
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
        assert init_fields == ["coord", "edges", "volumes", "mat_ids", "face_laws"]


# ─────────────────────────────────────────────────────────────────────
# One boundary-law check for the geometry and the mesh
# ─────────────────────────────────────────────────────────────────────


class TestOneBoundaryLawCheck:
    """The geometry's laws and a mesh's face laws are parsed by ONE function,
    ``parse_boundary_laws`` (the elegance review of step 3b). Claim kind:
    THEOREM (a ROUTE gate): both constructors resolve the same function
    object, and a decoy of it moves both refusals."""

    def test_both_constructors_bind_the_one_function(self):
        import orpheus.geometry.structured_geometry as geometry_module

        assert structured_module.parse_boundary_laws is geometry_module.parse_boundary_laws

    def test_a_decoy_moves_both_refusals(self, monkeypatch):
        import sys

        import orpheus.geometry.structured_geometry as geometry_module
        from orpheus.geometry import StructuredGeometry

        original = geometry_module.parse_boundary_laws

        def decoy(laws, coord, r_0, r_R, where):
            raise ValueError(f"{where}: the decoy boundary check")

        rebound = 0
        for module in list(sys.modules.values()):
            if getattr(module, "parse_boundary_laws", None) is original:
                monkeypatch.setattr(module, "parse_boundary_laws", decoy)
                rebound += 1
        assert rebound >= 2
        with pytest.raises(ValueError, match="StructuredGeometry.* the decoy boundary check"):
            StructuredGeometry.slab((0.0, 1.0), (0,), left=BC.vacuum, right=BC.vacuum)
        with pytest.raises(ValueError, match="Mesh1D.face_laws: the decoy boundary check"):
            _mesh()
