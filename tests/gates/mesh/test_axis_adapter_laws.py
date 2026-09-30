r"""The face laws of the adapter mesh ``SNProblem.from_axes`` synthesizes (P1 step 3b).

``SNProblem.from_axes`` stores the caller's axes and builds a legacy
:class:`~orpheus.mesh.Mesh1D` ADAPTER for the consumers still reading
``problem.mesh``. A mesh refuses ``None`` on a face, so the adapter carries
the laws SN computes (``orpheus/mesh/axis.py::_face_laws_of_axis``), two
guards with two retirements:

* an undeclared axis law becomes SN's reflective default
  (``ELEGANCE-DEBT[guard]`` #405, retires at P1 step 3c);
* the inner face of a hollow radial axis, which has no law slot, is the
  reflective cavity SN computes (``SCOPE-BOUNDARY[guard]`` #511).

Claim kind: THEOREM for each guard's declared behaviour; the reflective
inner face is the void cavity by ``tests/gates/mesh/test_hollow_inner_law.py
::test_a_reflective_inner_law_is_a_void_cavity``. First red: before the split
the adapter's hollow inner face had no witness (qa, review of step 3b).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem
from orpheus.mesh import AxisCoord, AxisMesh, Mesh1D, RadialAxisMesh
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.problem import SNProblem

pytestmark = pytest.mark.foundation

_VOID = "tests/gates/mesh/test_hollow_inner_law.py::test_a_reflective_inner_law_is_a_void_cavity"


def _adapter(axis) -> Mesh1D:
    """The 1-D adapter mesh ``SNProblem.from_axes`` synthesizes for ``axis``."""
    mesh = SNProblem.from_axes((axis,), Quadrature.gauss_legendre(8), {0: get_mixture("B", "2g")}).mesh
    if not isinstance(mesh, Mesh1D):
        raise AssertionError(f"a 1-D axis tuple gave a {type(mesh).__name__} adapter")
    return mesh


@pytest.mark.rests_on(_VOID)
@pytest.mark.parametrize("outer", [BC.vacuum, BC.reflective], ids=["vacuum", "reflective"])
def test_a_hollow_radial_axis_gets_a_reflective_inner_face(outer):
    axis = RadialAxisMesh(edges=np.array([0.5, 1.0, 2.0]), coord=AxisCoord.RADIAL_SPHERICAL, bc_outer=outer)
    mesh = _adapter(axis)
    assert mesh.coord is CoordSystem.SPHERICAL
    assert mesh.boundary_faces == (0.5, 2.0)
    assert mesh.face_laws[0] == BC.reflective
    assert mesh.face_laws[1] is outer


def test_a_solid_radial_axis_gets_its_outer_law_only():
    """The control: no inner face, so no reflective law is invented."""
    axis = RadialAxisMesh(edges=np.array([0.0, 1.0, 2.0]), coord=AxisCoord.RADIAL_SPHERICAL, bc_outer=BC.vacuum)
    assert _adapter(axis).face_laws == (BC.vacuum,)


@pytest.mark.parametrize(
    "low, high, expected",
    [
        pytest.param(None, BC.vacuum, (BC.reflective, BC.vacuum), id="undeclared-low"),
        pytest.param(BC.vacuum, None, (BC.vacuum, BC.reflective), id="undeclared-high"),
        pytest.param(None, None, (BC.reflective, BC.reflective), id="both-undeclared"),
        pytest.param(BC.vacuum, BC.vacuum, (BC.vacuum, BC.vacuum), id="declared-control"),
    ],
)
def test_an_undeclared_axis_law_becomes_reflective(low, high, expected):
    """The ELEGANCE-DEBT arm: an undeclared ``AxisMesh`` law is SN's
    reflective default on the adapter; a declared one is carried as given."""
    axis = AxisMesh(edges=np.linspace(0.0, 1.0, 3), bc_low=low, bc_high=high)
    assert _adapter(axis).face_laws == expected


def test_an_undeclared_radial_outer_law_becomes_reflective():
    axis = RadialAxisMesh(edges=np.array([0.5, 1.0, 2.0]), coord=AxisCoord.RADIAL_SPHERICAL)
    assert _adapter(axis).face_laws == (BC.reflective, BC.reflective)
