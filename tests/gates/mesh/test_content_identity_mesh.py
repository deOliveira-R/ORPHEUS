r"""Content identity of the mesh layer: ``FaceLaws``, ``CellEdges``, ``Mesh1D``,
``Mesh2D`` (#405 P1 step 5).

Spec: ``.claude/plans/reference_p1_spec.md`` §1.5, the mesh legs of S5.2 and
S5.3, and S5.10 (``Mesh2D``: read-only copied arrays with signed zero
canonicalised, ``==`` on twins, ``hash``). The encoder and the mixin are gated
in ``tests/gates/numerics/test_content_identity.py``.

First red, measured on ``main`` ``1dc31163``: the module fails to import
(``orpheus.numerics.content`` does not exist). With the encoder alone added,
the rows red on the legacy behaviour each rejects: ``hash(FaceLaws)`` raises
(``Mapping`` defines ``__eq__``); ``FaceLaws == dict`` is True; ``hash(Mesh1D)``
raises (its hand-written ``__eq__`` nulls ``__hash__``); ``Mesh2D == Mesh2D``
raises ``ValueError`` (the generated ``__eq__`` compares ndarrays); ``Mesh2D``
aliases the caller's writeable arrays and keeps ``-0.0``.
"""

from __future__ import annotations

from typing import Any, cast

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, StructuredGeometry
from orpheus.geometry.boundary import AlbedoBoundary, ReflectiveBoundary, VacuumInflow
from orpheus.mesh import CellEdges, CellsByCount, FaceLaws, Mesh1D, Mesh2D, Mesher
from orpheus.numerics.content import content_digest
from tests.gates._content_identity_helpers import (
    Entry,
    check_equal_pair,
    check_perturbation,
    check_pickle,
    check_population,
    leg,
    pair_ids,
    param_id,
    perturbation_ids,
    require,
)

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/mesh/test_content_identity_mesh.py"
_ENCODER = "tests/gates/numerics/test_content_identity.py"
_S51 = f"{_ENCODER}::test_s5_1_digests_and_hashes_are_seed_stable"
_S54 = f"{_ENCODER}::TestS54EncoderCanonicalForms"

_SLAB, _CYL = CoordSystem.CARTESIAN, CoordSystem.CYLINDRICAL
_LAWS_1D = {"xmin": BC.reflective, "xmax": BC.vacuum}
_LAWS_2D = {"xmin": BC.reflective, "xmax": BC.vacuum, "ymin": BC.reflective, "ymax": BC.vacuum}
_UP = float(np.nextafter(1.0, 2.0))


def _face_laws(laws=None, inventory=("xmin", "xmax")) -> FaceLaws:
    return FaceLaws.over(inventory, dict(_LAWS_1D if laws is None else laws), "test")


def _m1(coord=_SLAB, edges=(0.5, 1.0, 2.0), volumes=None, region_ids=(0, 1), region_materials=(0, 1), laws=None) -> Mesh1D:
    """A hollow two-cell mesh (both faces carry a law on every coordinate system), one region per cell."""
    e = np.asarray(edges, dtype=float)
    v = coord.measure(e) if volumes is None else np.asarray(volumes, dtype=float)
    return Mesh1D(
        coord=coord, edges=edges, volumes=v, region_ids=np.asarray(region_ids),
        region_materials=tuple(region_materials), face_laws=dict(_LAWS_1D if laws is None else laws),
    )


def _volume_one_ulp_up() -> np.ndarray:
    v = _SLAB.measure(np.array([0.5, 1.0, 2.0])).copy()
    v[0] = np.nextafter(v[0], 1.0)
    return v


def _meshed() -> Mesh1D:
    """The same cells as :func:`_m1`, built by the mesher from a geometry."""
    geometry = StructuredGeometry(
        coord=_SLAB, breakpoints=(0.5, 1.0, 2.0), mat_ids=(0, 1),
        boundaries=(BC.reflective, BC.vacuum),
    )
    return Mesher(geometry).partition(CellsByCount.uniform_width(1)).mesh


def _m2(edges_x=(0.5, 1.0, 2.0), edges_y=(0.0, 1.0), mat_map=((0,), (1,)), laws=None, coord=_SLAB) -> Mesh2D:
    return Mesh2D(
        np.array(edges_x, dtype=float), np.array(edges_y, dtype=float), np.array(mat_map),
        face_laws=dict(_LAWS_2D if laws is None else laws), coord=coord,
    )


# ── The roster ──────────────────────────────────────────────────────────────

_FACE_LAWS = Entry(
    cls=FaceLaws, base=_face_laws, parts=("items",), fields_are_parts=False,
    perturb={
        "items": (
            leg("a law changed", lambda: _face_laws({"xmin": BC.reflective, "xmax": AlbedoBoundary(0.3)})),
            leg("the laws swapped between faces", lambda: _face_laws({"xmin": BC.vacuum, "xmax": BC.reflective})),
            leg("a face dropped", lambda: _face_laws({"xmax": BC.vacuum}, inventory=("xmax",))),
        ),
    },
    pairs=(
        ("two builds", _face_laws, _face_laws),
        ("item order", _face_laws, lambda: FaceLaws((("xmax", BC.vacuum), ("xmin", BC.reflective)))),
        ("typed laws built twice", lambda: _face_laws({"xmin": ReflectiveBoundary("x"), "xmax": VacuumInflow()}),
         lambda: _face_laws({"xmin": ReflectiveBoundary(axis="x"), "xmax": VacuumInflow()})),
    ),
)

_CELL_EDGES = Entry(
    cls=CellEdges, base=lambda: CellEdges(np.array([0.0, 0.5, 1.0])), parts=("edges",),
    perturb={"edges": (leg("an interior edge", lambda: CellEdges(np.array([0.0, 0.25, 1.0]))),)},
    # A list is the array_like input ``parse_positions`` admits at runtime.
    pairs=(("a list vs an array", lambda: CellEdges(cast(Any, [0.0, 0.5, 1.0])), lambda: CellEdges(np.array([0, 0.5, 1]))),),
)

_MESH1D = Entry(
    # #405 P2 step 7b.2.0 (R7b2.0.6): the region labels and the region -> material map replace the stored
    # material ids; ``mat_ids`` is derived from them (``region_materials[region_ids]``) and is not content.
    cls=Mesh1D, base=_m1, parts=("coord", "edges", "volumes", "region_ids", "region_materials", "face_laws"),
    fields_are_parts=False,  # widths, centers, areas are derived: not content
    perturb={
        "coord": (leg("cylinder", lambda: _m1(coord=_CYL), "volumes"),),
        "edges": (leg("an interior edge one ulp", lambda: _m1(edges=(0.5, _UP, 2.0), volumes=(0.5, 1.0))),),
        "volumes": (leg("one volume one ulp", lambda: _m1(volumes=_volume_one_ulp_up())),),
        # Two labels swapped, the map swapped with them (a declared co-moving part), so ``mat_ids``
        # reads (0, 1) on both. The pure witness that labels are content with the map held fixed is
        # R7b2.0.5 in test_mesh1d_regions.py (two regions of one material relabelled).
        "region_ids": (leg("two labels swapped, the map with them", lambda: _m1(region_ids=(1, 0), region_materials=(1, 0)), "region_materials"),),
        "region_materials": (leg("a region's material", lambda: _m1(region_materials=(0, 2))),),
        "face_laws": (leg("the outer law", lambda: _m1(laws={"xmin": BC.reflective, "xmax": AlbedoBoundary(0.3)})),),
    },
    pairs=(
        ("two builds", _m1, _m1),
        ("edges as a tuple vs an array", _m1, lambda: _m1(edges=np.array([0.5, 1.0, 2.0]))),
        ("the mesher vs the constructor", _m1, _meshed),
        ("a BC built, not the constant", _m1, lambda: _m1(laws={"xmin": BC("reflective"), "xmax": BC("vacuum")})),
    ),
)

_MESH2D = Entry(
    cls=Mesh2D, base=_m2, parts=("edges_x", "edges_y", "mat_map", "face_laws", "coord"),
    perturb={
        "edges_x": (leg("an interior edge", lambda: _m2(edges_x=(0.5, 1.5, 2.0))),),
        "edges_y": (leg("the top edge", lambda: _m2(edges_y=(0.0, 2.0))),),
        "mat_map": (leg("a material id", lambda: _m2(mat_map=((0,), (2,)))),),
        "face_laws": (leg("one law", lambda: _m2(laws={**_LAWS_2D, "ymax": AlbedoBoundary(0.3)})),),
        "coord": (leg("cylinder", lambda: _m2(coord=_CYL)),),
    },
    pairs=(
        ("two builds", _m2, _m2),
        ("a numpy-integer vs an int material map", _m2,
         lambda: _m2(mat_map=np.array([[0], [1]], dtype=np.int64))),
        ("signed zero on an edge", _m2, lambda: _m2(edges_y=(-0.0, 1.0))),
    ),
)

ROSTER: tuple[Entry, ...] = (_FACE_LAWS, _CELL_EDGES, _MESH1D, _MESH2D)


# ── S5.3 ─────────────────────────────────────────────────────────────────────


@pytest.mark.rests_on(_S51)
@pytest.mark.parametrize("entry", ROSTER, ids=lambda e: e.id)
def test_s5_3_population(entry: Entry) -> None:
    """First red with the encoder present: ``FaceLaws`` has no content parts
    (neither a dataclass nor the mixin); ``Mesh1D``'s compare fields include
    the derived ``widths``, ``centers``, ``areas``."""
    check_population(entry)


@pytest.mark.rests_on(f"{_HERE}::test_s5_3_population")
@pytest.mark.parametrize(
    "entry,part,the_leg",
    [pytest.param(e, p, l, id=param_id(e.id, p, l[0])) for e, p, l in perturbation_ids(ROSTER)],
)
def test_s5_3_a_moved_part_moves_the_digest(entry: Entry, part: str, the_leg) -> None:
    check_perturbation(entry, part, the_leg)


@pytest.mark.rests_on(f"{_HERE}::test_s5_3_population")
def test_s5_3_mesh1d_coord_alone_moves_the_digest() -> None:
    """The ``coord`` part, isolated. The roster's ``coord`` leg must also move
    ``volumes`` (a mesh's volumes are its coordinate system's measure), so a
    content that dropped ``coord`` stayed green there (``[M]`` post-carve
    battery, arm "Mesh1D.coord dropped": 0 reds). One cell ``[a, b]`` with
    ``a + b = 1/pi`` has slab measure ``b - a`` equal to the cylinder's
    ``pi (b^2 - a^2)`` in exact arithmetic (``[M]`` 1 ulp apart in floats,
    inside the volume band), so the two meshes share edges, volumes, ids and
    laws and differ in ``coord`` alone."""
    a = 0.05
    edges = np.array([a, 1.0 / np.pi - a])
    volumes = _SLAB.measure(edges)
    slab = _m1(coord=_SLAB, edges=edges, volumes=volumes, region_ids=(0,), region_materials=(0,))
    cylinder = _m1(coord=_CYL, edges=edges, volumes=volumes, region_ids=(0,), region_materials=(0,))
    for name in ("edges", "volumes", "region_ids", "mat_ids"):
        require(np.array_equal(getattr(slab, name), getattr(cylinder, name)), f"activation: {name} differs")
    require(dict(slab.face_laws) == dict(cylinder.face_laws), "activation: the laws differ")
    require(content_digest(slab) != content_digest(cylinder), "the digest does not read coord")
    require(slab != cylinder and len({slab, cylinder}) == 2, "== or a set does not read coord")


# ── S5.2 ─────────────────────────────────────────────────────────────────────


@pytest.mark.rests_on(_S51)
@pytest.mark.parametrize(
    "entry,pair", [pytest.param(e, p, id=param_id(e.id, p[0])) for e, p in pair_ids(ROSTER)]
)
def test_s5_2_equal_content_is_one_value(entry: Entry, pair) -> None:
    check_equal_pair(entry, pair)


@pytest.mark.rests_on(_S51)
@pytest.mark.parametrize("entry", ROSTER, ids=lambda e: e.id)
def test_s5_2_pickle_round_trip(entry: Entry) -> None:
    check_pickle(entry)


@pytest.mark.rests_on(_S51)
def test_s5_2_face_laws_are_not_a_dict() -> None:
    """A ``FaceLaws`` is a value of its own type: not equal to a plain ``dict``
    of the same items (``[M]`` today True, through ``Mapping.__eq__``), while
    its items read back as that dict. Re-poses every ``mesh.face_laws ==
    {...}`` assertion as ``dict(mesh.face_laws) == {...}``."""
    laws = _face_laws()
    require(not (laws == dict(_LAWS_1D)), "a FaceLaws is not equal to a dict")
    require(laws != dict(_LAWS_1D), "!= agrees")
    require(dict(laws) == dict(_LAWS_1D), "its items read back as the dict")


# ── S5.10: Mesh2D ────────────────────────────────────────────────────────────


@pytest.mark.rests_on(f"{_HERE}::test_s5_2_equal_content_is_one_value")
class TestS510Mesh2D:
    """``Mesh2D`` becomes a value: its arrays are read-only COPIES (``[M]``
    today ``np.asarray`` aliases the caller's array and leaves it writeable),
    signed zero is canonicalised in what it stores (``[M]`` today kept), and
    twins compare equal and hash (``[M]`` today ``==`` raises ``ValueError``,
    ``hash`` raises ``TypeError``)."""

    def test_the_arrays_are_read_only_copies(self) -> None:
        ex, ey, mm = np.array([0.5, 1.0, 2.0]), np.array([0.0, 1.0]), np.array([[0], [1]])
        mesh = Mesh2D(ex, ey, mm, face_laws=dict(_LAWS_2D))
        for name, stored, given in (("edges_x", mesh.edges_x, ex), ("edges_y", mesh.edges_y, ey), ("mat_map", mesh.mat_map, mm)):
            require(not stored.flags.writeable, f"{name} is writeable")
            require(not np.shares_memory(stored, given), f"{name} aliases the caller's array")
        before = content_digest(mesh)
        ex[1], mm[0, 0] = 1.25, 7
        require(mesh.edges_x[1] == 1.0 and mesh.mat_map[0, 0] == 0, "a caller's write moved the mesh")
        require(content_digest(mesh) == before, "a caller's write moved the digest")
        with pytest.raises(ValueError, match="read-only"):
            mesh.edges_x[1] = 1.25

    def test_signed_zero_is_stored_canonical(self) -> None:
        mesh = _m2(edges_y=(-0.0, 1.0))
        require(not bool(np.signbit(mesh.edges_y).any()), "-0.0 is stored")

    def test_twins_compare_equal_and_hash(self) -> None:
        a, b = _m2(), _m2()
        require(a == b, "twins compare equal")
        require(hash(a) == hash(b), "twins hash equal")
        require(not (a == _m2(mat_map=((1,), (1,)))), "a different material map is another mesh")
