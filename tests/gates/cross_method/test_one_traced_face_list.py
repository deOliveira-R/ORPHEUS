r"""``trace_face_shapes`` is the one spelling of the traced-face list both trace spaces read (layer DAG Phase 4).

Before Phase 4 the list ``[(face_name, face_shape) for label in face_labels(axes)]`` was spelled twice, in
``SNProblem.boundary_face_layout`` and ``DiffusionMesh._init_core``: a twin. Its retirement makes the mesh function
load-bearing (``retirement-audit`` D.18), and this file is its witness, two ways:

* the ROUTE: a decoy list (the honest one reversed) installed on ``orpheus.transport.phase_spaces`` moves BOTH a freshly
  posed SN slab's angular trace and a diffusion slab's scalar trace, and an activation counter proves each pose called
  it; the control leg (no decoy) keeps the honest order;
* the CENSUS: over ``orpheus/`` (AST), the face-layout and trace-space builders are called only where the list has its
  one reader: ``FaceLayout.from_named_shapes`` at 3 sites (the angular layout, ``ScalarTraceSpace.for_faces``'s own body
  and the radial-characteristic corner layout), ``ScalarTraceSpace.for_faces`` and
  ``AngularTraceSpace.from_quadrature_and_layout`` at 1 each, all in ``orpheus/transport/phase_spaces.py`` but the two
  numerics-internal sites. Positive control: the census finds ``phase_spaces.py``'s ``from_named_shapes`` call.

Cross-method because it poses both an SN problem and a diffusion mesh.
"""
from __future__ import annotations

import ast
import pathlib

import pytest

import orpheus
from orpheus.diffusion import DiffusionMesh
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.problem import SNProblem
from orpheus.transport import phase_spaces
from tests.gates.transport import _phase_space_cases as cases

pytestmark = pytest.mark.foundation

ORPHEUS_ROOT = pathlib.Path(orpheus.__file__).resolve().parent


def _poses():
    mm = cases.slab()
    mesh = cases.legacy_mesh(mm)
    sn = SNProblem(mesh, Quadrature.gauss_legendre(4), mm.materials)
    diffusion = DiffusionMesh(mesh, mm.materials)
    scalar_layout = diffusion.scalar_trace.layout
    if scalar_layout is None:
        pytest.fail("the diffusion scalar trace carries no face layout")
    return list(sn.angular_trace.layout.faces), list(scalar_layout.faces)


def test_the_control_poses_read_the_honest_order() -> None:
    angular, scalar = _poses()
    if angular != ["xmin", "xmax"] or scalar != ["xmin", "xmax"]:
        pytest.fail(f"control: the slab traces read {angular} and {scalar}")


@pytest.mark.rests_on("tests/gates/mesh/test_trace_face_shapes.py::test_the_traced_faces_are_the_hand_inventory", "tests/gates/cross_method/test_one_traced_face_list.py::test_the_control_poses_read_the_honest_order")
def test_both_traces_move_with_the_one_list(monkeypatch: pytest.MonkeyPatch) -> None:
    honest = phase_spaces.trace_face_shapes
    calls = []

    def decoy(axes):
        calls.append(len(axes))
        return tuple(reversed(honest(axes)))

    monkeypatch.setattr(phase_spaces, "trace_face_shapes", decoy)
    angular, scalar = _poses()
    if not calls:
        pytest.fail("activation: neither pose called the list through phase_spaces")
    if angular != ["xmax", "xmin"]:
        pytest.fail(f"the SN angular trace did not move with the list: {angular}")
    if scalar != ["xmax", "xmin"]:
        pytest.fail(f"the diffusion scalar trace did not move with the list: {scalar}")


_BUILDERS = ("from_named_shapes", "for_faces", "from_quadrature_and_layout")
_EXPECTED = {
    "from_named_shapes": {"transport/phase_spaces.py", "numerics/spaces/scalar_trace_space.py",
                          "numerics/spaces/radial_characteristic_space.py"},
    "for_faces": {"transport/phase_spaces.py"},
    "from_quadrature_and_layout": {"transport/phase_spaces.py"},
}


def _call_sites() -> dict[str, list[str]]:
    sites: dict[str, list[str]] = {name: [] for name in _BUILDERS}
    files = sorted(ORPHEUS_ROOT.rglob("*.py"))
    if len(files) < 100:
        pytest.fail(f"the census read {len(files)} files under {ORPHEUS_ROOT}")
    for path in files:
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr in _BUILDERS:
                sites[node.func.attr].append(path.relative_to(ORPHEUS_ROOT).as_posix())
    return sites


def test_the_trace_builders_are_called_only_beside_the_one_list() -> None:
    sites = _call_sites()
    if "transport/phase_spaces.py" not in sites["from_named_shapes"]:
        pytest.fail("positive control: the census does not find phase_spaces' own from_named_shapes call")
    for name, expected in _EXPECTED.items():
        if set(sites[name]) != expected or len(sites[name]) != len(expected):
            pytest.fail(f"{name} is called at {sorted(sites[name])}; the one-list sites are {sorted(expected)}")
