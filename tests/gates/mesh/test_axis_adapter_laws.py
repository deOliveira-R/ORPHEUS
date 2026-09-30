r"""The face laws of the adapter mesh ``SNProblem.from_axes`` synthesizes (P1 step 3b).

``SNProblem.from_axes`` stores the caller's axes and builds a legacy
:class:`~orpheus.mesh.Mesh1D` ADAPTER for the consumers still reading
``problem.mesh``. A mesh refuses ``None`` on a face, so the adapter carries
the laws SN computes (``orpheus/mesh/axis.py::_face_laws_of_axis``), which are:

* the declared laws of the axis, carried as given (the undeclared-law arm,
  ``ELEGANCE-DEBT[guard]`` #405, retired at P1 step 3c: an axis cannot be
  built without its laws, so there is nothing to default);
* the inner face of a hollow radial axis, which has no law slot, is the
  reflective cavity SN computes (``SCOPE-BOUNDARY[guard]`` #511).

Step 3c adds the axis classes' own laws: a law is required and parsed
(``orpheus/geometry/structured_geometry.py::parse_boundary_law``), so an
omitted law, ``None`` and a non-law object are refused at construction;
``with_uniform_bc`` (its one consumer was ``_apply_default_bcs``) is gone; and
the adapter round trip ``axes_from_legacy_mesh(legacy_mesh_from_axes(axes))``
keeps every declared law, in 1-D and 2-D.

Claim kind: THEOREM for each row; the reflective inner face is the void cavity
by ``tests/gates/mesh/test_hollow_inner_law.py
::test_a_reflective_inner_law_is_a_void_cavity``. First reds: before the split
the adapter's hollow inner face had no witness (qa, review of step 3b); the
pre-3c axis classes (laws defaulting to ``None``, no parse) red every
``TestTheAxisLaws`` refusal; a swapped pair in either adapter (``bc_low`` read
from ``xmax``, or the 2-D ``face_laws`` built with ``ymin``/``ymax`` crossed)
reds ``TestTheRoundTrip``.
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
    assert mesh.boundary_points == (0.5, 2.0)
    assert tuple(mesh.face_laws) == ("xmin", "xmax")
    assert mesh.face_laws["xmin"] == BC.reflective
    assert mesh.face_laws["xmax"] is outer


def test_a_solid_radial_axis_gets_its_outer_law_only():
    """The control: no inner face, so no reflective law is invented."""
    axis = RadialAxisMesh(edges=np.array([0.0, 1.0, 2.0]), coord=AxisCoord.RADIAL_SPHERICAL, bc_outer=BC.vacuum)
    laws = _adapter(axis).face_laws
    assert tuple(laws) == ("xmax",) and laws["xmax"] == BC.vacuum




# ─────────────────────────────────────────────────────────────────────
# Step 3c: an axis declares its laws
# ─────────────────────────────────────────────────────────────────────

_EDGES = np.array([0.0, 0.5, 2.0])
_ONE_PARSER = (
    "tests/gates/mesh/test_mesh2d_face_laws.py::TestTheOneElementParser::"
    "test_every_declared_law_passes_the_one_parser"
)


class TestTheAxisLaws:
    """The axis classes refuse an undeclared law (the pre-3c classes defaulted
    each law to ``None`` and parsed nothing, so every refusal row is red on
    them). The ``None`` and non-law rows of the five owners are also gated,
    frame by frame, by ``TestTheOneElementParser``."""

    def test_an_omitted_cartesian_law_is_refused(self):
        with pytest.raises(TypeError, match="bc_low"):
            AxisMesh(edges=_EDGES)  # type: ignore[call-arg]  # the refusal input
        with pytest.raises(TypeError, match="bc_high"):
            AxisMesh(edges=_EDGES, bc_low=BC.vacuum)  # type: ignore[call-arg]  # the refusal input

    def test_an_omitted_radial_law_is_refused(self):
        with pytest.raises(TypeError, match="bc_outer"):
            RadialAxisMesh(edges=_EDGES, coord=AxisCoord.RADIAL_SPHERICAL)  # type: ignore[call-arg]  # the refusal input

    @pytest.mark.rests_on(_ONE_PARSER)
    @pytest.mark.parametrize("slot", ["bc_low", "bc_high"])
    @pytest.mark.parametrize(
        "bad, fragment",
        [(None, "None is not a boundary law"), ("vacuum", "must be a BC tag or a BoundaryTraceLaw")],
        ids=["none", "non-law"],
    )
    def test_a_cartesian_law_that_is_not_a_law(self, slot, bad, fragment):
        laws = {"bc_low": BC.vacuum, "bc_high": BC.vacuum} | {slot: bad}
        with pytest.raises(TypeError, match=fragment) as info:
            AxisMesh(edges=_EDGES, bc_low=laws["bc_low"], bc_high=laws["bc_high"])
        assert f"AxisMesh.{slot}" in str(info.value)

    @pytest.mark.rests_on(_ONE_PARSER)
    @pytest.mark.parametrize(
        "bad, fragment",
        [(None, "None is not a boundary law"), (0, "must be a BC tag or a BoundaryTraceLaw")],
        ids=["none", "non-law"],
    )
    def test_a_radial_law_that_is_not_a_law(self, bad, fragment):
        with pytest.raises(TypeError, match=fragment) as info:
            RadialAxisMesh(edges=_EDGES, coord=AxisCoord.RADIAL_CYLINDRICAL, bc_outer=bad)
        assert "RadialAxisMesh.bc_outer" in str(info.value)

    def test_the_default_filling_verb_is_gone(self):
        """``with_uniform_bc`` retired with its one consumer ``_apply_default_bcs``
        (``retirement-audit`` G.24). First red: the verb re-added on either class
        or on the protocol."""
        from orpheus.mesh.axis import Axis1D

        for owner in (AxisMesh, RadialAxisMesh, Axis1D):
            assert not hasattr(owner, "with_uniform_bc"), owner.__name__

    def test_the_bc_table_carries_no_none(self):
        """The declared table is the laws given, by identity."""
        low, high = BC("albedo", {"albedo": 0.3}), BC.vacuum
        assert AxisMesh(edges=_EDGES, bc_low=low, bc_high=high).bc == {"min": low, "max": high}
        radial = RadialAxisMesh(edges=_EDGES, coord=AxisCoord.RADIAL_SPHERICAL, bc_outer=high)
        assert radial.bc["outer"] is high


def _round_trip(axes):
    from orpheus.mesh.axis import axes_from_legacy_mesh, legacy_mesh_from_axes

    return axes_from_legacy_mesh(legacy_mesh_from_axes(axes))


def _typed_law():
    from orpheus.geometry.boundary import VacuumInflow

    return VacuumInflow()


class TestTheRoundTrip:
    """``axes_from_legacy_mesh(legacy_mesh_from_axes(axes))`` keeps every
    declared law, the same object on the same endpoint. The laws on one tuple
    are pairwise distinct, so a swap of two faces cannot read equal. First red:
    a crossed pair in either adapter (e.g. the 2-D ``face_laws`` built with
    ``ymin`` and ``ymax`` exchanged, or ``bc_low`` read from ``xmax``)."""

    @staticmethod
    def _assert_same_laws(before, after):
        assert len(before) == len(after)
        for k, (a, b) in enumerate(zip(before, after)):
            assert type(a) is type(b), (k, type(a), type(b))
            np.testing.assert_array_equal(a.edges, b.edges)
            assert a.coord == b.coord
            assert a.endpoints == b.endpoints
            for endpoint in a.endpoints:
                assert b.bc[endpoint] is a.bc[endpoint], (k, endpoint, a.bc[endpoint], b.bc[endpoint])

    @pytest.mark.parametrize(
        "axis",
        [
            pytest.param(lambda: AxisMesh(edges=_EDGES, bc_low=BC.vacuum, bc_high=BC("albedo", {"albedo": 0.3})),
                         id="cartesian"),
            pytest.param(lambda: RadialAxisMesh(edges=_EDGES, coord=AxisCoord.RADIAL_SPHERICAL, bc_outer=_typed_law()),
                         id="solid-sphere"),
            pytest.param(lambda: RadialAxisMesh(edges=np.array([0.5, 1.0, 2.0]), coord=AxisCoord.RADIAL_CYLINDRICAL,
                                                bc_outer=BC.vacuum),
                         id="hollow-cylinder"),
        ],
    )
    def test_one_dimension(self, axis):
        axes = (axis(),)
        self._assert_same_laws(axes, _round_trip(axes))

    def test_two_dimensions(self):
        axes = (
            AxisMesh(edges=_EDGES, bc_low=BC.vacuum, bc_high=BC.reflective),
            AxisMesh(edges=np.array([0.0, 1.0, 2.0, 3.0]), bc_low=BC("albedo", {"albedo": 0.3}), bc_high=_typed_law()),
        )
        self._assert_same_laws(axes, _round_trip(axes))
