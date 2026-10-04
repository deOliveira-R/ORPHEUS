r"""The 1-D mesh carries region labels (#405 P2 step 7b.2.0).

A mesh does not hold the geometry it came from: an external mesh (Gmsh,
Ansys) has no originating geometry and often no coordinate system, but it
does carry cell labels (physical volume groups). So :class:`Mesh1D` holds
``region_ids``, a region label on each cell, and ``region_materials``, the
material of each region, and ``mat_ids`` is DERIVED,
``region_materials[region_ids]`` (the user's ruling of 2026-10-03,
``.claude/plans/reference_cache.md``, "Step 7b.1 landed; step 7b.2 ruled").
Until this step the mesher composed its cell -> interval map with the
interval -> material map and kept only the composite, so the A|B|A regions 0
and 2, both material 0, could not be told apart and a per-region weight could
not be read.

The representation (spec §1.7b.2, "7b.2.0"): the labels are DENSE, ``0 ..
R-1``, and ``region_materials`` is a tuple indexed by label. A label is the
region's position, the same convention as a
:class:`~orpheus.numerics.mesh_free_function.RegionwiseConstant`'s rows and
the geometry's intervals; an importer renumbers external tags densely. Every
label has a material and every material entry labels a cell (no spectator
entry: it would change the digest and no reading).

Gate ids ``R7b2.0.<n>`` (``.claude/plans/reference_p2_spec.md`` §1.7b.2).
First red on the tree this file was written against (``78bbe024``): every row
fails, because ``Mesh1D`` takes ``mat_ids`` and has no ``region_ids``; the
content-identity roster (``test_content_identity_mesh.py``) and the
field-set pin (``test_mesh1d.py::TestTheRetirements``) were re-posed in the
same change. Every test is ``foundation``, claim kind THEOREM.
"""
from __future__ import annotations

import dataclasses
import re
from collections.abc import Callable

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, StructuredGeometry
from orpheus.mesh import AxisCoord, AxisMesh, CellsByCount, Mesh1D, Mesher, RadialAxisMesh
from orpheus.numerics.content import content_digest

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/mesh/test_mesh1d_regions.py"
_MESHER = "tests/gates/mesh/test_mesher.py::TestTheLift::test_every_cell_lies_in_exactly_one_interval"
_ROSTER = "tests/gates/mesh/test_content_identity_mesh.py::test_s5_3_population"

_SLAB, _CYLINDER, _SPHERE = CoordSystem.CARTESIAN, CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL


def _mesh(region_ids=(0, 1, 2), region_materials=(0, 1, 0), edges=(0.0, 0.5, 1.5, 2.0), coord=_SLAB) -> Mesh1D:
    e = np.asarray(edges, dtype=float)
    laws = {"xmax": BC.vacuum} if coord is not _SLAB and e[0] == 0.0 else {"xmin": BC.reflective, "xmax": BC.vacuum}
    return Mesh1D(
        coord=coord, edges=e, volumes=coord.measure(e), region_ids=region_ids,
        region_materials=region_materials, face_laws=laws,
    )


def _aba(coord: CoordSystem, cells_per_region: tuple[int, int, int] = (2, 4, 2)) -> tuple[StructuredGeometry, Mesh1D]:
    """The A|B|A body of the trajectory-resolvent cross-check: 0.5, 1.5, 2.0 cm, materials (0, 1, 0)."""
    g = StructuredGeometry.from_thicknesses(
        coord=coord, thicknesses=(0.5, 1.0, 0.5), mat_ids=(0, 1, 0),
        boundaries=(BC.vacuum,) if coord is not _SLAB else (BC.reflective, BC.vacuum),
    )
    return g, Mesher(g).partition(tuple(CellsByCount.uniform_width(n) for n in cells_per_region)).mesh


# ─────────────────────────────────────────────────────────────────────
# R7b2.0.1 — the constructor takes the labels; mat_ids is not stored
# ─────────────────────────────────────────────────────────────────────


def test_r7b2_0_1_the_constructor_takes_the_labels() -> None:
    """The init fields are the discretisation with labels; ``mat_ids`` is derived like ``widths``: a field
    that is not an init argument and not content (``compare=False``), one definition (the elegance
    review of 7b.2.0, S1: the class's derived quantities are all ``init=False, compare=False``)."""
    init_fields = [f.name for f in dataclasses.fields(Mesh1D) if f.init]
    assert init_fields == ["coord", "edges", "volumes", "region_ids", "region_materials", "face_laws"]
    derived = {f.name: f for f in dataclasses.fields(Mesh1D)}["mat_ids"]
    assert not derived.init and not derived.compare, "mat_ids is derived, never given and never content"
    with pytest.raises(TypeError):
        Mesh1D(  # the retired spelling is the subject
            coord=_SLAB, edges=np.array([0.0, 1.0]), volumes=np.array([1.0]),
            mat_ids=np.array([0]), face_laws={"xmin": BC.vacuum, "xmax": BC.vacuum},  # pyright: ignore[reportCallIssue]
        )


# ─────────────────────────────────────────────────────────────────────
# R7b2.0.2 — the labels are admitted, keyed refusals
# ─────────────────────────────────────────────────────────────────────

#: (case id, builder, error type, fragment of the message).
_REFUSALS: list[tuple[str, Callable[[], object], type[Exception], str]] = [
    ("label-count", lambda: _mesh(region_ids=(0, 1)), ValueError, "Mesh1D.region_ids: 3 cell(s) need 3 region label(s)"),
    ("label-float", lambda: _mesh(region_ids=(0, 1.0, 2)), TypeError, "a region label is an int"),
    ("label-bool", lambda: _mesh(region_ids=(0, True, 2)), TypeError, "a region label is an int"),
    ("label-negative", lambda: _mesh(region_ids=(0, -1, 2), region_materials=(0, 1, 0)), ValueError,
     "is a non-negative index"),
    ("label-without-material", lambda: _mesh(region_ids=(0, 1, 2), region_materials=(0, 1)), ValueError,
     "has no material"),
    ("material-without-cell", lambda: _mesh(region_ids=(0, 1, 1), region_materials=(0, 1, 0)), ValueError,
     "labels no cell"),
    ("material-float", lambda: _mesh(region_materials=(0, 1.0, 0)), TypeError, "a material id is an int"),
    ("materials-a-mapping", lambda: _mesh(region_materials={0: 0, 1: 1, 2: 0}), TypeError,  # type: ignore[arg-type]
     "(a tuple of material ids), got dict"),
    # The elegance review of 7b.2.0, B1: a dict given as labels was read as its keys, (0, 1, 2), a silently
    # wrong partition; a set and a one-pass iterator have no positions either.
    ("labels-a-mapping", lambda: _mesh(region_ids={0: 2, 1: 1, 2: 0}), TypeError, "sequence of int, got dict"),
    ("labels-a-set", lambda: _mesh(region_ids={0, 1, 2}), TypeError, "sequence of int, got set"),
    ("labels-an-iterator", lambda: _mesh(region_ids=iter((0, 1, 2))), TypeError, "sequence of int, got tuple_iterator"),
    ("materials-a-0d-array", lambda: _mesh(region_materials=np.array(5)), TypeError, "got ndarray"),
]


@pytest.mark.parametrize(
    "build, error, fragment", [r[1:] for r in _REFUSALS], ids=[r[0] for r in _REFUSALS],
)
def test_r7b2_0_2_refusal(build, error, fragment) -> None:
    with pytest.raises(error, match=re.escape(fragment)):
        build()


def test_r7b2_0_2_refusal_fragments_are_disjoint() -> None:
    """Each refusal raises its own message: no fragment matches another row's error."""
    messages = {}
    for case, build, error, _ in _REFUSALS:
        with pytest.raises(error) as info:
            build()
        messages[case] = str(info.value)
    for case, _, _, fragment in _REFUSALS:
        same = {c for c, m in messages.items() if fragment in m}
        siblings = {c for c, _, _, f in _REFUSALS if f == fragment}
        assert same == siblings, f"{case}: fragment {fragment!r} also matches {same - siblings}"


def test_r7b2_0_2_the_positive_legs() -> None:
    """Admitted: one region (an unlabelled mesh), a region of several cells, labels out of cell order, numpy ints."""
    one = _mesh(region_ids=(0, 0, 0), region_materials=(7,))
    assert one.mat_ids.tolist() == [7, 7, 7]
    unordered = _mesh(region_ids=(2, 0, 1), region_materials=(5, 6, 5))
    assert unordered.mat_ids.tolist() == [5, 5, 6]
    numpy_ints = _mesh(region_ids=np.array([0, 1, 2], dtype=np.int64), region_materials=tuple(np.array([0, 1, 0])))
    assert numpy_ints == _mesh()
    as_range = _mesh(region_ids=range(3), region_materials=range(3))  # a range is a sequence
    assert as_range == _mesh(region_materials=(0, 1, 2))


# ─────────────────────────────────────────────────────────────────────
# R7b2.0.3 — mat_ids is derived, read-only; the labels are a read-only copy
# ─────────────────────────────────────────────────────────────────────


def test_r7b2_0_3_mat_ids_is_the_map_of_the_labels() -> None:
    """``mat_ids == region_materials[region_ids]``, elementwise, on a labelling where the two differ."""
    mesh = _mesh(region_ids=(2, 1, 0), region_materials=(4, 9, 3))
    np.testing.assert_array_equal(mesh.mat_ids, [3, 9, 4])
    np.testing.assert_array_equal(mesh.mat_ids, np.asarray(mesh.region_materials)[mesh.region_ids])
    assert mesh.mat_ids.dtype.kind == "i" and mesh.region_ids.dtype.kind == "i"
    assert isinstance(mesh.region_materials, tuple)
    assert all(type(m) is int for m in mesh.region_materials)


def test_r7b2_0_3_read_only_copies() -> None:
    labels = np.array([0, 1, 2])
    mesh = _mesh(region_ids=labels)
    labels[0] = 1
    np.testing.assert_array_equal(mesh.region_ids, [0, 1, 2])
    for array in (mesh.region_ids, mesh.mat_ids):
        assert not array.flags.writeable


# ─────────────────────────────────────────────────────────────────────
# R7b2.0.4 — the mesher writes the labels it knows
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.rests_on(_MESHER)
@pytest.mark.parametrize("coord", [_SLAB, _CYLINDER, _SPHERE], ids=lambda c: c.name.lower())
def test_r7b2_0_4_the_mesher_labels_are_the_interval_indices(coord) -> None:
    """Each cell's label is the index of the geometry interval containing it (decided by containment,
    independently of the mesher's grouping), and the map is the geometry's ``mat_ids``."""
    g, mesh = _aba(coord)
    r = np.asarray(g.breakpoints)
    left, right = mesh.edges[:-1], mesh.edges[1:]
    contained = (r[None, :-1] <= left[:, None]) & (right[:, None] <= r[None, 1:])
    owner = np.argmax(contained, axis=1)
    np.testing.assert_array_equal(mesh.region_ids, owner)
    assert mesh.region_materials == tuple(g.mat_ids)
    np.testing.assert_array_equal(mesh.mat_ids, np.asarray(g.mat_ids)[owner])


@pytest.mark.rests_on(f"{_HERE}::test_r7b2_0_4_the_mesher_labels_are_the_interval_indices")
@pytest.mark.parametrize("coord", [_SLAB, _CYLINDER, _SPHERE], ids=lambda c: c.name.lower())
def test_r7b2_0_4_aba_regions_0_and_2_are_distinct(coord) -> None:
    """The X1 witness of the carve: regions 0 and 2 carry one material and remain two regions, so a
    per-region table with different rows 0 and 2 reads differently on their cells."""
    _, mesh = _aba(coord)
    fuel = mesh.mat_ids == 0
    assert sorted(set(mesh.region_ids[fuel].tolist())) == [0, 2]
    table = np.array([1.0, 0.0, 3.0])  # a per-region value, rows 0 and 2 different
    on_fuel = table[mesh.region_ids][fuel]
    assert set(on_fuel.tolist()) == {1.0, 3.0}, "regions 0 and 2 read as one"


@pytest.mark.parametrize("coord", [_SLAB, _SPHERE], ids=lambda c: c.name.lower())
def test_r7b2_0_4_refinement_keeps_the_labels(coord) -> None:
    """``Mesher.refine`` re-partitions every interval, so the labels stay the interval indices."""
    g, _ = _aba(coord)
    mesher = Mesher(g).partition(tuple(CellsByCount.uniform_width(n) for n in (1, 2, 1)))
    fine = mesher.refine(2).mesh
    np.testing.assert_array_equal(fine.region_ids, [0, 0, 1, 1, 1, 1, 2, 2])
    assert fine.region_materials == (0, 1, 0)


# ─────────────────────────────────────────────────────────────────────
# R7b2.0.5 — content identity: the labels and the map are content
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.rests_on(_ROSTER)
def test_r7b2_0_5_labels_are_content_with_mat_ids_fixed() -> None:
    """Two meshes with equal cells, laws and ``mat_ids`` but different labels are two values: a label is
    the region's position (the row a per-region table is read at), so relabelling is not a renaming."""
    a = _mesh(region_ids=(0, 1, 2), region_materials=(0, 1, 0))
    b = _mesh(region_ids=(2, 1, 0), region_materials=(0, 1, 0))
    np.testing.assert_array_equal(a.mat_ids, b.mat_ids)  # activation: the composite cannot tell them apart
    assert content_digest(a) != content_digest(b)
    assert a != b and len({a, b}) == 2


@pytest.mark.rests_on(_ROSTER)
def test_r7b2_0_5_the_map_is_content_while_the_mesh_holds_it() -> None:
    """While the region -> material map lives on the mesh (until #522 moves it to the system), it is
    content: equal values must be substitutable, and two meshes reading different ``mat_ids`` are not."""
    a, b = _mesh(region_materials=(0, 1, 0)), _mesh(region_materials=(0, 1, 2))
    assert content_digest(a) != content_digest(b) and a != b


def test_r7b2_0_5_the_mesher_and_the_constructor_agree() -> None:
    g, meshed = _aba(_SPHERE, (1, 1, 1))
    built = _mesh(coord=_SPHERE)
    assert meshed == built and content_digest(meshed) == content_digest(built)


# ─────────────────────────────────────────────────────────────────────
# R7b2.0.6 — the legacy adapter: one region per material run
# ─────────────────────────────────────────────────────────────────────


def test_r7b2_0_6_legacy_adapter_labels_material_runs() -> None:
    """``legacy_mesh_from_axes`` has no geometry, only a material per cell. It labels each maximal run of
    one material as its own region: the finest partition the material map determines that never
    merges two disjoint pieces (labelling by material id would fuse A|B|A's regions 0 and 2)."""
    from orpheus.mesh.axis import legacy_mesh_from_axes

    axis = AxisMesh(edges=np.linspace(0.0, 3.0, 7), bc_low=BC.reflective, bc_high=BC.vacuum)
    mesh = legacy_mesh_from_axes((axis,), mat_map=np.array([0, 0, 1, 1, 0, 0]))
    assert isinstance(mesh, Mesh1D), "one axis is a 1-D mesh"
    np.testing.assert_array_equal(mesh.region_ids, [0, 0, 1, 1, 2, 2])
    assert mesh.region_materials == (0, 1, 0)
    np.testing.assert_array_equal(mesh.mat_ids, [0, 0, 1, 1, 0, 0])

    radial = RadialAxisMesh(edges=np.array([0.0, 0.5, 1.5, 2.0]), coord=AxisCoord.RADIAL_SPHERICAL, bc_outer=BC.vacuum)
    one = legacy_mesh_from_axes((radial,))  # no map: one material, so one region
    assert isinstance(one, Mesh1D), "one axis is a 1-D mesh"
    np.testing.assert_array_equal(one.region_ids, [0, 0, 0])
    assert one.region_materials == (0,)


# ─────────────────────────────────────────────────────────────────────
# R7b2.0.7 — the polymorphic mints carry the labels
# ─────────────────────────────────────────────────────────────────────


def test_r7b2_0_7_with_distinct_cell_ids() -> None:
    """Each cell its own region and its own material, ``0 .. N-1``; cells and laws carried through."""
    mesh = _mesh(region_ids=(0, 0, 1), region_materials=(4, 4))
    distinct = mesh.with_distinct_cell_ids()
    np.testing.assert_array_equal(distinct.region_ids, [0, 1, 2])
    assert distinct.region_materials == (0, 1, 2)
    np.testing.assert_array_equal(distinct.mat_ids, [0, 1, 2])
    np.testing.assert_array_equal(distinct.edges, mesh.edges)
    assert dict(distinct.face_laws) == dict(mesh.face_laws)


def test_r7b2_0_7_replace_keeps_the_labels() -> None:
    """``dataclasses.replace(mesh, face_laws=...)``, the spelling 9 test sites use, keeps labels and map."""
    mesh = _mesh(region_ids=(2, 1, 0), region_materials=(4, 9, 3))
    other = dataclasses.replace(mesh, face_laws={"xmin": BC.vacuum, "xmax": BC.vacuum})
    np.testing.assert_array_equal(other.region_ids, mesh.region_ids)
    assert other.region_materials == mesh.region_materials
