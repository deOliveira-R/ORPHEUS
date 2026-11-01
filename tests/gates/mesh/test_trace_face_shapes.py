r"""``trace_face_shapes(axes)``: the traced-face list both trace spaces are laid out over (layer DAG Phase 4).

The REFERENCE is a hand-written literal per geometry: the faces carrying a trace, in ``face_labels`` order (axis
ascending, endpoint in the axis's own order), each with the per-axis cell counts of the OTHER axes. A solid radial
axis has no face at its pole and a hollow one has no trace at its inner surface (#511), so both read one face,
``xmax`` with endpoint ``"outer"``. The rows are the geometries a trace is built on: slab, solid and hollow sphere
and cylinder, a non-square 2-D Cartesian mesh (``nx != ny``, so a transposed shape is visible) and a synthetic
3-axis tuple (six faces).

Named first reds: reversing the label order (``slab`` and ``cart2d`` rows), a shape that keeps the face's own axis
(``cart2d``, ``3d``), and a label list that admits the hollow inner endpoint (``hollow_*`` rows).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import StructuredGeometry
from orpheus.geometry.deck import MirrorDeck
from orpheus.mesh import AxisMesh, CellsByCount, Mesh2D, Mesher, axes_from_legacy_mesh, trace_face_shapes

pytestmark = pytest.mark.foundation


def _radial(kind: str, inner: float) -> tuple:
    build = StructuredGeometry.sphere if kind == "sphere" else StructuredGeometry.cylinder
    roles = {"outer": 0} if inner == 0.0 else {"outer": 0, "inner": 0}
    mesh = Mesher(build((inner, 1.0), (0,), **roles)).partition(CellsByCount.uniform_width(4)).mesh
    return axes_from_legacy_mesh(mesh)


def _slab() -> tuple:
    geometry = StructuredGeometry.slab((0.0, 1.0), (0,), left=MirrorDeck(), right=MirrorDeck())
    return axes_from_legacy_mesh(Mesher(geometry).partition(CellsByCount.uniform_width(4)).mesh)


def _cart2d() -> tuple:
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 1.0, 4), edges_y=np.linspace(0.0, 1.0, 3), mat_map=np.zeros((3, 2), dtype=int),
        face_roles={face: MirrorDeck() for face in ("xmin", "xmax", "ymin", "ymax")},
    )
    return axes_from_legacy_mesh(mesh)


def _three_axes() -> tuple:
    return tuple(AxisMesh(np.linspace(0.0, 1.0, n + 1), MirrorDeck(), MirrorDeck()) for n in (2, 3, 5))


#: (axes builder, the hand literal: (face name, endpoint, shape) in the canonical order)
_CASES = {
    "slab": (_slab, (("xmin", "min", ()), ("xmax", "max", ()))),
    "sphere": (lambda: _radial("sphere", 0.0), (("xmax", "outer", ()),)),
    "hollow_sphere": (lambda: _radial("sphere", 0.5), (("xmax", "outer", ()),)),
    "cylinder": (lambda: _radial("cylinder", 0.0), (("xmax", "outer", ()),)),
    "hollow_cylinder": (lambda: _radial("cylinder", 0.5), (("xmax", "outer", ()),)),
    "cart2d": (_cart2d, (("xmin", "min", (2,)), ("xmax", "max", (2,)), ("ymin", "min", (3,)), ("ymax", "max", (3,)))),
    "3d": (_three_axes, (
        ("xmin", "min", (3, 5)), ("xmax", "max", (3, 5)), ("ymin", "min", (2, 5)), ("ymax", "max", (2, 5)),
        ("zmin", "min", (2, 3)), ("zmax", "max", (2, 3)),
    )),
}


@pytest.mark.parametrize("case", list(_CASES))
def test_the_traced_faces_are_the_hand_inventory(case: str) -> None:
    build, expected = _CASES[case]
    got = tuple((label.face_name, label.endpoint, shape) for label, shape in trace_face_shapes(build()))
    if got != expected:
        pytest.fail(f"{case}: trace_face_shapes reads {got}, the hand inventory is {expected}")


def test_every_shape_is_a_tuple_of_ints() -> None:
    """A shape is the codimension-1 cell count, an ``int`` tuple (a numpy scalar would hash and compare like an int
    and still break a ``FaceLayout`` key built from it)."""
    for case, (build, _) in _CASES.items():
        for label, shape in trace_face_shapes(build()):
            if type(shape) is not tuple or not all(type(n) is int for n in shape):
                pytest.fail(f"{case}/{label.face_name}: shape {shape!r} is not a tuple of int")
