r"""The face roles of the axes, and of the adapter mesh ``SNProblem.from_axes`` synthesizes (P1 step 3b/3c, re-posed at wave 2 S2b).

``SNProblem.from_axes`` stores the caller's axes and builds a legacy :class:`~orpheus.mesh.Mesh1D` ADAPTER for the consumers still
reading ``problem.mesh``. Before wave 2 the adapter carried the LAWS SN computed (``_face_laws_of_axis``) and synthesized a
reflective cavity for the inner face of a hollow radial axis (``SCOPE-BOUNDARY[guard]`` #511). Both are gone: an axis declares a
face ROLE per endpoint (a deck or an ``int`` tag), the adapter carries it verbatim, and the hollow inner face is the TAG the axis
declares (``tests/gates/mesh/test_radial_inner_role.py``, which owns the inner-role and the synthesis claims; the two rows that
pinned the synthesized reflective cavity are retired with it and not repeated here).

What this module keeps, re-spelled on roles:

* ``TestTheAxisRoles``: an axis cannot be built without its roles, and a role that is not a deck or a tag is refused at
  construction (the pre-3c classes defaulted each law to ``None`` and parsed nothing, so every refusal row is red on them); the
  default-filling verb ``with_uniform_bc`` is gone; the declared table is the roles given, by identity;
* ``TestTheRoundTrip``: ``axes_from_legacy_mesh(legacy_mesh_from_axes(axes))`` keeps every declared role, the same object on
  the same endpoint, in 1-D and 2-D.

Claim kind: THEOREM for each row. First reds: a swapped pair in either adapter (``role_low`` read from ``xmax``, or the 2-D
``face_roles`` built with ``ymin``/``ymax`` crossed) reds ``TestTheRoundTrip``; the pre-3c axis classes red every refusal row of
``TestTheAxisRoles``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.deck import MirrorDeck
from orpheus.mesh import AxisCoord, AxisMesh, RadialAxisMesh

pytestmark = pytest.mark.foundation

_EDGES = np.array([0.0, 0.5, 2.0])
_ONE_PARSER = (
    "tests/gates/mesh/test_mesh2d_face_roles.py::TestTheOneElementParser::"
    "test_every_declared_role_passes_the_one_parser"
)
_ROLE_FRAGMENT = "a face role is a MirrorDeck, a TranslationDeck or an int tag"


class TestTheAxisRoles:
    """The axis classes refuse an undeclared or malformed role. The ``None`` and non-role rows of the five owners are also gated,
    frame by frame, by ``TestTheOneElementParser``."""

    def test_an_omitted_cartesian_role_is_refused(self):
        with pytest.raises(TypeError, match="role_low"):
            AxisMesh(edges=_EDGES)  # type: ignore[call-arg]  # the refusal input
        with pytest.raises(TypeError, match="role_high"):
            AxisMesh(edges=_EDGES, role_low=0)  # type: ignore[call-arg]  # the refusal input
        # control: the complete declaration builds
        assert AxisMesh(edges=_EDGES, role_low=0, role_high=0).roles == {"min": 0, "max": 0}

    def test_an_omitted_radial_role_is_refused(self):
        with pytest.raises(TypeError, match="role_outer"):
            RadialAxisMesh(edges=_EDGES, coord=AxisCoord.RADIAL_SPHERICAL)  # type: ignore[call-arg]  # the refusal input
        assert RadialAxisMesh(edges=_EDGES, coord=AxisCoord.RADIAL_SPHERICAL, role_outer=0).roles == {"outer": 0}

    @pytest.mark.rests_on(_ONE_PARSER)
    @pytest.mark.parametrize("slot", ["role_low", "role_high"])
    @pytest.mark.parametrize("bad", [None, "vacuum", 0.5], ids=["none", "string", "float"])
    def test_a_cartesian_role_that_is_not_a_role(self, slot, bad):
        roles = {"role_low": 0, "role_high": 0} | {slot: bad}
        with pytest.raises(TypeError, match=_ROLE_FRAGMENT) as info:
            AxisMesh(edges=_EDGES, role_low=roles["role_low"], role_high=roles["role_high"])
        assert f"AxisMesh.{slot}" in str(info.value)

    @pytest.mark.rests_on(_ONE_PARSER)
    @pytest.mark.parametrize("bad", [None, "vacuum", 0.5], ids=["none", "string", "float"])
    def test_a_radial_role_that_is_not_a_role(self, bad):
        with pytest.raises(TypeError, match=_ROLE_FRAGMENT) as info:
            RadialAxisMesh(edges=_EDGES, coord=AxisCoord.RADIAL_CYLINDRICAL, role_outer=bad)
        assert "RadialAxisMesh.role_outer" in str(info.value)

    def test_the_default_filling_verb_is_gone(self):
        """``with_uniform_bc`` retired with its one consumer ``_apply_default_bcs``
        (``retirement-audit`` G.24). First red: the verb re-added on either class
        or on the protocol."""
        from orpheus.mesh.axis import Axis1D

        for owner in (AxisMesh, RadialAxisMesh, Axis1D):
            assert not hasattr(owner, "with_uniform_bc"), owner.__name__

    def test_the_role_table_is_the_roles_given(self):
        """The declared table is the roles given, by identity (a deck object stays the same object; no ``None`` is ever stored)."""
        low, high = MirrorDeck(), 7
        axis = AxisMesh(edges=_EDGES, role_low=low, role_high=high)
        assert axis.roles["min"] is low and axis.roles["max"] == 7
        radial = RadialAxisMesh(edges=_EDGES, coord=AxisCoord.RADIAL_SPHERICAL, role_outer=high)
        assert radial.roles["outer"] == 7
        assert all(role is not None for role in (*axis.roles.values(), *radial.roles.values()))


def _round_trip(axes):
    from orpheus.mesh.legacy_adapter import axes_from_legacy_mesh, legacy_mesh_from_axes

    return axes_from_legacy_mesh(legacy_mesh_from_axes(axes))


class TestTheRoundTrip:
    """``axes_from_legacy_mesh(legacy_mesh_from_axes(axes))`` keeps every declared role, the same object on the same endpoint. The
    roles on one tuple are pairwise distinct, so a swap of two faces cannot read equal. First red: a crossed pair in either adapter
    (e.g. the 2-D ``face_roles`` built with ``ymin`` and ``ymax`` exchanged, or ``role_low`` read from ``xmax``)."""

    @staticmethod
    def _assert_same_roles(before, after):
        assert len(before) == len(after)
        for k, (a, b) in enumerate(zip(before, after)):
            assert type(a) is type(b), (k, type(a), type(b))
            np.testing.assert_array_equal(a.edges, b.edges)
            assert a.coord == b.coord
            assert a.endpoints == b.endpoints
            for endpoint in a.endpoints:
                assert b.roles[endpoint] == a.roles[endpoint], (k, endpoint, a.roles[endpoint], b.roles[endpoint])

    @pytest.mark.parametrize(
        "axis",
        [
            pytest.param(lambda: AxisMesh(edges=_EDGES, role_low=3, role_high=5), id="cartesian"),
            pytest.param(lambda: AxisMesh(edges=_EDGES, role_low=MirrorDeck(), role_high=5), id="cartesian-deck"),
            pytest.param(lambda: RadialAxisMesh(edges=_EDGES, coord=AxisCoord.RADIAL_SPHERICAL, role_outer=7), id="solid-sphere"),
            pytest.param(lambda: RadialAxisMesh(edges=np.array([0.5, 1.0, 2.0]), coord=AxisCoord.RADIAL_CYLINDRICAL,
                                                role_outer=2, role_inner=9),
                         id="hollow-cylinder"),
        ],
    )
    def test_one_dimension(self, axis):
        axes = (axis(),)
        self._assert_same_roles(axes, _round_trip(axes))

    def test_two_dimensions(self):
        axes = (
            AxisMesh(edges=_EDGES, role_low=3, role_high=MirrorDeck()),
            AxisMesh(edges=np.array([0.0, 1.0, 2.0, 3.0]), role_low=5, role_high=6),
        )
        self._assert_same_roles(axes, _round_trip(axes))
