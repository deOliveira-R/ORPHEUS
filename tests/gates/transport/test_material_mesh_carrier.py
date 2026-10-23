r"""``MaterialMesh`` is a ``BoundaryCarrier`` (wave 2, slice S2b, gate 21; plan section 2.10, decision 51).

``bind(carrier, conditions)`` takes the object that owns the faces and their colouring symmetry. SN and diffusion pass their
``MaterialMesh``; CP, MoC and MC pass their ``Mesh1D`` (S2a gate 14 and S4). So the material mesh grows the two members of the
Protocol: ``face_roles`` (a property, ``roles_of_axes(self.axes)``, a hollow radial axis's inner role included; ONE reader,
the one S1 made public) and ``symmetry`` (a property: one Cartesian ``AxisMesh`` passes its cell edges and ``mat_map`` values,
raw, to ``colouring_symmetry``, which merges the runs; a ``RadialAxisMesh`` the radial identity; anything else
``UncomputedSymmetries``).

The ladder::

    test_the_roles_are_the_axes_roles_inner_included     foundation: face_roles == roles_of_axes(axes)
    test_there_is_one_reader_of_the_roles                foundation: no second reader (AST), the docstring's "reads" is true
        ^
    test_the_symmetry_is_the_colouring_symmetry          edge: agrees with Mesh1D.symmetry on a slab (one fact, two readers)
    test_a_material_mesh_is_a_boundary_carrier           composition: the Protocol, and a real bind through it

First red at the base tip: ``AttributeError`` on ``MaterialMesh.face_roles`` (and the ``orpheus.transport`` parent import
chain raises ``ImportError`` first, since ``method.py`` still imports the retired ``BC``). Mutations: a ``face_roles`` that
omits the hollow inner role; a ``symmetry`` that ignores the materials (a symmetric and an asymmetric slab then agree);
a second ``roles_of_axes``-like reader in ``material_mesh.py``.
"""

from __future__ import annotations

import ast
import inspect
import os
from pathlib import Path

import pytest

import orpheus

pytestmark = pytest.mark.foundation
_THIS = "tests/gates/transport/test_material_mesh_carrier.py"
ROOT = Path(os.path.abspath(orpheus.__file__)).parent


def _materials(ids: tuple[int, ...]) -> dict:
    from orpheus.derivations.common.xs_library import get_mixture

    return {i: get_mixture(("A", "B", "C")[i % 3], "2g") for i in ids}


def _mesh(coord: str, breakpoints: tuple[float, ...], mats: tuple[int, ...], **roles):
    from orpheus.geometry import StructuredGeometry
    from orpheus.mesh import CellsByCount, Mesher

    build = {"slab": StructuredGeometry.slab, "cylinder": StructuredGeometry.cylinder}[coord]
    return Mesher(build(breakpoints, mats, **roles)).partition(tuple(CellsByCount.uniform_volume(2) for _ in mats)).mesh


def _material_mesh(coord, breakpoints, mats, **roles):
    from orpheus.transport.mesh.material_mesh import MaterialMesh

    mesh = _mesh(coord, breakpoints, mats, **roles)
    return MaterialMesh(mesh, _materials(tuple(sorted(set(mats))))), mesh


def test_the_roles_are_the_axes_roles_inner_included() -> None:
    r"""``MaterialMesh.face_roles`` equals ``roles_of_axes(axes)``: a slab's two faces, a solid cylinder's outer face only
    (the pole is not a face), a HOLLOW cylinder's inner face as well (the tag the axis holds and the trace inventory does
    not list, #511). The returned mapping is the carrier's ``Mapping[str, FaceRole]`` (decks are values).

    Red input at the base tip: no such member. Mutation: a ``face_roles`` listing the trace faces only.
    """
    from orpheus.geometry.deck import MirrorDeck
    from orpheus.mesh.axis import roles_of_axes

    slab, _ = _material_mesh("slab", (0.0, 1.0, 2.0), (0, 1), left=MirrorDeck(), right=3)
    assert dict(slab.face_roles) == dict(roles_of_axes(slab.axes)) == {"xmin": MirrorDeck(), "xmax": 3}
    solid, _ = _material_mesh("cylinder", (0.0, 1.0, 2.0), (0, 1), outer=2)
    assert dict(solid.face_roles) == {"xmax": 2}
    hollow, _ = _material_mesh("cylinder", (0.5, 1.0, 2.0), (0, 1), outer=0, inner=1)
    assert dict(hollow.face_roles) == dict(roles_of_axes(hollow.axes)) == {"xmin": 1, "xmax": 0}


@pytest.mark.rests_on(f"{_THIS}::test_the_roles_are_the_axes_roles_inner_included")
def test_there_is_one_reader_of_the_roles() -> None:
    r"""The property reads ``roles_of_axes`` (its docstring's "reads" is true), and ``roles_of_axes`` is DEFINED once in
    ``orpheus/`` (S1 made it public: do not add a second reader). The scanner is validated by a planted second definition.

    Red input at the base tip: ``MaterialMesh`` has no ``face_roles``. Mutation: a private ``_roles`` copy in
    ``material_mesh.py`` reading ``axis.role_low`` and friends.
    """
    from orpheus.transport.mesh.material_mesh import MaterialMesh

    def defined_in(source: str) -> int:
        return sum(1 for n in ast.walk(ast.parse(source)) if isinstance(n, ast.FunctionDef) and n.name == "roles_of_axes")

    assert defined_in("def roles_of_axes(axes):\n    pass\n") == 1, "the scanner cannot see a definition"
    files = sorted(ROOT.rglob("*.py"))
    assert len(files) > 100
    assert sum(defined_in(p.read_text()) for p in files) == 1, "roles_of_axes is defined more than once"
    prop = MaterialMesh.face_roles
    getter = getattr(prop, "fget", None) or getattr(prop, "func")
    assert "roles_of_axes" in inspect.getsource(getter), "MaterialMesh.face_roles does not read roles_of_axes"


@pytest.mark.rests_on(f"{_THIS}::test_the_roles_are_the_axes_roles_inner_included")
def test_the_symmetry_is_the_colouring_symmetry() -> None:
    r"""``MaterialMesh.symmetry`` is a ``Symmetries`` and, on a one-axis Cartesian mesh, EQUALS the ``Mesh1D.symmetry`` of the
    same slab (one fact, two readers; the material mesh passes cell edges and ``mat_map`` values, the mesh its region edges,
    and ``colouring_symmetry`` merges equal adjacent runs): a symmetric slab with the same tag on both faces has the mirror
    in its group (two elements), an asymmetric one (two materials) only the identity, and a symmetric slab with a mirror deck on each face
    (the same role at both ends) the mirror too (``[M]`` 2, 1, 2 elements on the Mesh1D of the same slabs). A radial mesh has the radial identity (one element); a 2-D Cartesian mesh is ``UncomputedSymmetries``.

    Red input at the base tip: no such member. Mutation: a ``symmetry`` that ignores ``mat_map`` (the asymmetric slab then
    has the mirror).
    """
    from orpheus.geometry import ColouringSymmetry, UncomputedSymmetries
    from orpheus.geometry.deck import MirrorDeck
    from orpheus.mesh import Mesh2D
    from orpheus.transport.mesh.material_mesh import MaterialMesh

    import numpy as np

    sym, sym_mesh = _material_mesh("slab", (0.0, 1.0, 2.0), (0, 0), left=0, right=0)
    asym, asym_mesh = _material_mesh("slab", (0.0, 1.0, 2.0), (0, 1), left=0, right=0)
    decks, decks_mesh = _material_mesh("slab", (0.0, 1.0, 2.0), (0, 0), left=MirrorDeck(), right=MirrorDeck())
    for carrier, mesh, elements in ((sym, sym_mesh, 2), (asym, asym_mesh, 1), (decks, decks_mesh, 2)):
        assert isinstance(carrier.symmetry, ColouringSymmetry) and len(carrier.symmetry.elements) == elements, (len(carrier.symmetry.elements), elements)
        assert carrier.symmetry == mesh.symmetry, "the material mesh and the mesh disagree on the colouring symmetry of one slab"
    radial, _ = _material_mesh("cylinder", (0.0, 1.0, 2.0), (0, 1), outer=0)
    assert isinstance(radial.symmetry, ColouringSymmetry) and len(radial.symmetry.elements) == 1

    plane = Mesh2D(np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0, 2.0]), np.zeros((2, 2), dtype=int), face_roles={f: 0 for f in ("xmin", "xmax", "ymin", "ymax")})
    assert isinstance(MaterialMesh(plane, _materials((0,))).symmetry, UncomputedSymmetries)


@pytest.mark.rests_on(f"{_THIS}::test_the_symmetry_is_the_colouring_symmetry")
def test_a_material_mesh_is_a_boundary_carrier() -> None:
    r"""``isinstance(material_mesh, BoundaryCarrier)`` (the Protocol is ``runtime_checkable``) and ``bind`` accepts it: the
    binding of a hollow body reads the inner tag from the carrier, and a lab-frame source on the symmetric slab's two-face
    tag is refused with the MATERIAL MESH's group (the invariance check reads the carrier's symmetry, not a second
    argument: ``bind`` has no ``symmetry`` parameter, S2a gate 14).

    Red input at the base tip: no such member. Mutation: a ``symmetry`` returning the identity group for the symmetric slab
    (the lab-frame condition then binds).
    """
    from tests.gates.transport._boundary_conditions_helpers import ProfileSampler

    from orpheus.transport import boundary_conditions as B

    hollow, _ = _material_mesh("cylinder", (0.5, 1.0, 2.0), (0, 1), outer=0, inner=1)
    assert isinstance(hollow, B.BoundaryCarrier)
    binding = B.bind(hollow, {0: B.WallCondition.vacuum(), 1: B.WallCondition.specular(1.0)})
    assert binding.tags["xmin"] == 1 and binding.tags["xmax"] == 0

    sym, _ = _material_mesh("slab", (0.0, 1.0, 2.0), (0, 0), left=0, right=0)
    lab = B.WallCondition.framed(B.LabFrameInflow(ProfileSampler(1.0)))
    with pytest.raises(B.NonInvariantConditionError):
        B.bind(sym, {0: lab})
    asym, _ = _material_mesh("slab", (0.0, 1.0, 2.0), (0, 1), left=0, right=0)
    B.bind(asym, {0: lab})
