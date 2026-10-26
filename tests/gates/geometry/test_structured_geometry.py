"""Foundation tests for :class:`StructuredGeometry`.

These tests pin the geometry value: the defining refusals of the value (its
breakpoints, its material ids, its derived boundary), the thickness
constructor, the retirement of the string kind tag, and the named
geometries. They are foundation-tier (software invariants): they verify the
code-shape contract, not a physics claim. The gate ids S2.1 to S2.5 are those
of the P1 verification specification (``.claude/plans/reference_p1_spec.md``
§1.2). The mesh built from a geometry (``TestMeshingAGeometry``) is
``tests/gates/mesh/test_mesher.py``, where it moved in the layer-DAG
re-filing (`.claude/plans/layer_dag.md`, ruling 21).
"""
from __future__ import annotations

import itertools
import math
from typing import Any

import numpy as np
import pytest

import orpheus.geometry
import orpheus.geometry.structured_geometry as structured_geometry_module
from orpheus.geometry import (
    CoordSystem,
    MirrorDeck,
    StructuredGeometry,
)


pytestmark = pytest.mark.foundation

_SLAB = CoordSystem.CARTESIAN
_CYLINDER = CoordSystem.CYLINDRICAL
_SPHERE = CoordSystem.SPHERICAL
_CURVILINEAR = (_CYLINDER, _SPHERE)


def _geometry(**overrides: Any) -> StructuredGeometry:
    """A valid two-interval slab, with the named fields replaced."""
    fields: dict[str, Any] = dict(
        coord=_SLAB,
        breakpoints=(0.0, 1.0, 2.0),
        mat_ids=(0, 1),
        face_roles=(0, 0),
    )
    fields.update(overrides)
    return StructuredGeometry(**fields)


# ─────────────────────────────────────────────────────────────────────
# S2.1 — the breakpoint laws: the value's defining refusals
# ─────────────────────────────────────────────────────────────────────

#: (case id, overrides, error type, the fragment that keys the refusal).
#: One row per clause of the constructor; S2.1's disjointness row reads
#: every message against every other row's fragment.
_BREAKPOINT_REFUSALS = [
    ("fewer-than-two", dict(breakpoints=(1.0,), mat_ids=()),
     ValueError, "at least 2 breakpoints"),
    ("equal-pair", dict(breakpoints=(0.0, 1.0, 1.0)),
     ValueError, "strictly increasing"),
    ("decreasing-pair", dict(breakpoints=(0.0, 2.0, 1.0)),
     ValueError, "strictly increasing"),
    ("nan", dict(breakpoints=(0.0, math.nan, 2.0)),
     ValueError, r"breakpoints\[1\] is NaN, which is not a number"),  # NaN: parse_real's own refusal since #405 P1 step 5
    ("inf", dict(breakpoints=(0.0, 1.0, math.inf)),
     ValueError, r"breakpoints\[2\] is infinite"),  # the shared finite parse since #559
    ("negative-radius-cylinder",
     dict(coord=_CYLINDER, breakpoints=(-0.5, 1.0), mat_ids=(0,)),
     ValueError, "starts at r_0 >= 0"),
    ("negative-radius-sphere",
     dict(coord=_SPHERE, breakpoints=(-0.5, 1.0), mat_ids=(0,)),
     ValueError, "starts at r_0 >= 0"),
    ("material-id-count", dict(mat_ids=(0,)),
     ValueError, "one material id per interval"),
    ("material-id-not-int", dict(mat_ids=(0, 1.0)),
     TypeError, "a material id is an int"),
]


class TestBreakpointLaws:
    """S2.1: the breakpoints and material ids a geometry admits."""

    @pytest.mark.parametrize(
        "overrides, error, fragment",
        [pytest.param(o, e, f, id=i) for i, o, e, f in _BREAKPOINT_REFUSALS],
    )
    def test_refusal(self, overrides, error, fragment):
        with pytest.raises(error, match=fragment):
            _geometry(**overrides)

    def test_refusal_fragments_are_disjoint(self):
        """Each refusal's message carries its own fragment and no other's."""
        messages = {}
        for case, overrides, error, _ in _BREAKPOINT_REFUSALS:
            with pytest.raises(error) as caught:
                _geometry(**overrides)
            messages[case] = str(caught.value)
        for case, _, _, fragment in _BREAKPOINT_REFUSALS:
            for other, message in messages.items():
                own = [f for c, _, _, f in _BREAKPOINT_REFUSALS if c == other]
                if fragment in own:
                    continue
                assert fragment not in message, (
                    f"the fragment of {case!r} ({fragment!r}) also appears in "
                    f"the refusal of {other!r}: {message!r}"
                )

    def test_a_slab_admits_any_origin(self):
        g = _geometry(breakpoints=(-1.0, 0.0, 2.0))
        assert g.breakpoints == (-1.0, 0.0, 2.0)

    def test_a_sphere_may_be_hollow(self):
        g = StructuredGeometry(
            coord=_SPHERE,
            breakpoints=(0.1, 1.0),
            mat_ids=(0,),
            face_roles=(1, 0),
        )
        assert g.breakpoints[0] == 0.1

    def test_adjacent_intervals_may_share_a_material(self):
        g = _geometry(mat_ids=(0, 0))
        assert g.mat_ids == (0, 0)

    def test_entries_are_canonicalised_to_float_and_int(self):
        """Integers and numpy scalars are parsed at the boundary, bit for bit."""
        g = _geometry(breakpoints=(0, np.float64(1.0), 2), mat_ids=(np.int64(3), 1))
        assert all(type(b) is float for b in g.breakpoints)
        assert all(type(m) is int for m in g.mat_ids)
        assert g.breakpoints == (0.0, 1.0, 2.0)
        assert g.mat_ids == (3, 1)

    def test_frozen(self):
        g = _geometry()
        with pytest.raises(AttributeError):
            g.coord = _SPHERE  # type: ignore[misc]

    def test_fields_are_keyword_only(self):
        """Two adjacent tuples of numbers cannot be swapped positionally."""
        with pytest.raises(TypeError):
            StructuredGeometry(_SLAB, (0.0, 1.0), (0,), (0, 0))  # type: ignore[misc]


# ─────────────────────────────────────────────────────────────────────
# S2.2 — the boundary is derived from the coordinate system and r_0
# ─────────────────────────────────────────────────────────────────────

#: (coord, r_0) → the boundary points of the geometry [r_0, 2.0].
_BOUNDARY_TABLE = [
    (_SLAB, 0.0, (0.0, 2.0)),
    (_SLAB, 0.5, (0.5, 2.0)),
    (_CYLINDER, 0.0, (2.0,)),
    (_CYLINDER, 0.5, (0.5, 2.0)),
    (_SPHERE, 0.0, (2.0,)),
    (_SPHERE, 0.5, (0.5, 2.0)),
]


class TestTheBoundaryIsDerived:
    """S2.2: one role per point of the region's topological boundary."""

    @pytest.mark.parametrize(
        "coord, r_0, points", _BOUNDARY_TABLE,
        ids=[f"{c.name.lower()}-r0={r}" for c, r, _ in _BOUNDARY_TABLE],
    )
    def test_boundary_points(self, coord, r_0, points):
        roles = tuple(0 for _ in points)
        g = StructuredGeometry(
            coord=coord, breakpoints=(r_0, 2.0), mat_ids=(0,), face_roles=roles,
        )
        assert g.boundary_points == points
        assert len(g.face_roles) == len(points)

    def test_a_role_at_the_centre_is_refused(self):
        """A DISCRIMINATION row: the old table refused this input too, for
        another reason ('requires 1 BC'); the new message names the centre."""
        with pytest.raises(ValueError) as caught:
            StructuredGeometry(
                coord=_SPHERE, breakpoints=(0.0, 2.0), mat_ids=(0,),
                face_roles=(0, 0),
            )
        assert "the centre r = 0" in str(caught.value)
        assert "requires 1 BC" not in str(caught.value)

    @pytest.mark.parametrize("coord", _CURVILINEAR, ids=lambda c: c.name.lower())
    def test_a_hollow_body_without_an_inner_role_is_refused(self, coord):
        with pytest.raises(ValueError, match="inner surface, which needs its own role"):
            StructuredGeometry(
                coord=coord, breakpoints=(0.5, 2.0), mat_ids=(0,),
                face_roles=(0,),
            )

    def test_a_slab_with_one_role_is_refused(self):
        with pytest.raises(ValueError, match="a slab has two boundary points"):
            _geometry(face_roles=(0,))

    def test_none_is_not_a_face_role(self):
        with pytest.raises(TypeError, match="a face role is a MirrorDeck, a TranslationDeck or an int tag, got None"):
            _geometry(face_roles=(None, 0))

    def test_face_roles_are_parsed_like_the_other_fields(self):
        """A list is canonicalised to a tuple; a non-sequence is refused."""
        assert _geometry(face_roles=[0, 1]).face_roles == (0, 1)
        with pytest.raises(TypeError, match="face_roles must be a sequence"):
            _geometry(face_roles=0)

    def test_a_role_is_a_deck_or_a_tag(self):
        """RE-POSED from ``test_a_law_is_a_tag_or_a_typed_law`` (wave 2, slice S1).

        A declaration is EITHER a deck (the domain continues across the face) or an ``int``
        tag (the domain ends there, and a problem binds the tag to a condition); a bare string
        is neither, and a response is not a role: the geometry names no response, and a
        function-valued inflow has no role spelling (it is bound to a tag where the method is
        known).
        """
        with pytest.raises(TypeError, match="a face role is a MirrorDeck, a TranslationDeck or an int tag"):
            _geometry(face_roles=("vacuum", 0))

    def test_a_deck_is_accepted(self):
        """The positive leg; without it the guard could reject every role."""
        g = StructuredGeometry(
            coord=_SLAB, breakpoints=(0.0, 2.0), mat_ids=(0,), face_roles=(MirrorDeck(), 3),
        )
        assert g.face_roles[0] == MirrorDeck() and g.face_roles[1] == 3


# ─────────────────────────────────────────────────────────────────────
# S2.3, S2.4 — breakpoints are stored; thicknesses are folded once
# ─────────────────────────────────────────────────────────────────────

#: A thickness list on which the sequential fold and a correctly rounded
#: cumulative sum differ (at the fourth breakpoint), so a constructor that
#: re-associated the sum would be seen.
_ASSOCIATION_SENSITIVE = (0.1, 0.2, 0.3, 0.4, 0.5)


class TestFromThicknesses:
    """S2.3: the thickness constructor is the left fold, bitwise."""

    def test_the_input_discriminates_the_association(self):
        """Positive control: on this input the two sums differ."""
        sequential = tuple(itertools.accumulate(_ASSOCIATION_SENSITIVE, initial=0.0))
        rounded = tuple(
            math.fsum(_ASSOCIATION_SENSITIVE[:k])
            for k in range(len(_ASSOCIATION_SENSITIVE) + 1)
        )
        assert sequential != rounded

    @pytest.mark.parametrize("r_0", [0.0, 0.25, -1.5])
    def test_breakpoints_are_the_left_fold(self, r_0):
        g = StructuredGeometry.from_thicknesses(
            coord=_SLAB,
            thicknesses=_ASSOCIATION_SENSITIVE,
            mat_ids=range(5),
            face_roles=(0, 0),
            r_0=r_0,
        )
        folded = [r_0]
        for t in _ASSOCIATION_SENSITIVE:
            folded.append(folded[-1] + t)
        assert g.breakpoints == tuple(folded)
        assert g.mat_ids == (0, 1, 2, 3, 4)

    @pytest.mark.parametrize("thickness", [0.0, -1.0], ids=["zero", "negative"])
    def test_a_non_positive_thickness_is_refused(self, thickness):
        """RE-POSED from ``TestRegion::test_zero_thickness_rejected`` and
        ``test_negative_thickness_rejected``: a thickness <= 0 is a
        non-increasing breakpoint pair."""
        with pytest.raises(ValueError, match="strictly increasing"):
            StructuredGeometry.from_thicknesses(
                coord=_SLAB, thicknesses=(1.0, thickness), mat_ids=(0, 1),
                face_roles=(0, 0),
            )


class TestFromThicknessesParsesLikeTheConstructor:
    """S1 of the elegance review: the thickness constructor admits exactly
    the scalars the constructor admits (no coercion of a string or a bool)."""

    @pytest.mark.parametrize(
        "thicknesses, r_0",
        [(("0.5", "2.0"), 0.0), ((True, True), 0.0), ((0.5, 2.0), "0.0")],
        ids=["string-thickness", "bool-thickness", "string-origin"],
    )
    def test_a_non_real_scalar_is_refused(self, thicknesses, r_0):
        with pytest.raises(TypeError, match="must be a real number"):
            StructuredGeometry.from_thicknesses(
                coord=_SLAB, thicknesses=thicknesses, mat_ids=(0, 1),
                face_roles=(0, 0), r_0=r_0,
            )

    def test_negative_zero_is_canonicalised(self):
        """-0.0 and 0.0 are one breakpoint; the stored bits are +0.0."""
        g = _geometry(breakpoints=(-0.0, 1.0, 2.0))
        assert math.copysign(1.0, g.breakpoints[0]) == 1.0


class TestBreakpointsAreStored:
    """S2.4: a breakpoint is never re-derived from widths."""

    _E = (0.0, 0.17, 0.45, 0.62, 1.0)

    def test_the_input_discriminates_a_re_derivation(self):
        """Positive control: re-adding the widths misses two breakpoints by
        one ULP (the census's site ``test_g_adjoint_reciprocity.py:244``)."""
        widths = np.diff(self._E).tolist()
        assert tuple(itertools.accumulate(widths, initial=0.0)) != self._E

    def test_breakpoints_round_trip_bitwise(self):
        g = StructuredGeometry(
            coord=_SLAB, breakpoints=self._E, mat_ids=(0, 1, 2, 3),
            face_roles=(0, 0),
        )
        assert g.breakpoints == self._E


# ─────────────────────────────────────────────────────────────────────
# S2.5 — the string kind tag retires
# ─────────────────────────────────────────────────────────────────────


class TestTheKindTagRetires:
    @pytest.mark.parametrize("tag", ["SLB", "CYL", "SPH", "HEX", "sph"])
    def test_a_string_coordinate_is_refused(self, tag):
        with pytest.raises(TypeError, match="must be a CoordSystem member"):
            _geometry(coord=tag)

    @pytest.mark.parametrize(
        "name", ["_GEOMETRY_TO_COORD", "_GEOMETRY_TO_N_ENDPOINTS", "Region"],
    )
    def test_the_tag_machinery_is_gone(self, name):
        assert not hasattr(structured_geometry_module, name)
        assert not hasattr(orpheus.geometry, name)

    @pytest.mark.parametrize("name", ["geometry", "regions", "bcs", "n_endpoints"])
    def test_the_old_fields_are_gone(self, name):
        assert not hasattr(_geometry(), name)


# ─────────────────────────────────────────────────────────────────────
# The factories
# ─────────────────────────────────────────────────────────────────────


class TestWignerSeitzPinCell:
    def test_default_construction(self):
        """The outer face is tag 0, which the Wigner-Seitz model binds to a white wall: part of the model, not a default."""
        g = StructuredGeometry.wigner_seitz_pin_cell(
            r_fuel=0.9, r_clad=1.1, pitch=3.6,
        )
        assert g.coord is _CYLINDER
        assert g.boundary_points == (g.breakpoints[-1],)
        assert g.face_roles == (0,)
        assert g.mat_ids == (2, 1, 0)

    def test_breakpoints_are_the_radii(self):
        """The factory states radii, so the radii are the breakpoints."""
        g = StructuredGeometry.wigner_seitz_pin_cell(
            r_fuel=0.9, r_clad=1.1, pitch=3.6,
        )
        assert g.breakpoints == (0.0, 0.9, 1.1, float(3.6 / np.sqrt(np.pi)))

    def test_extent_is_the_cell_radius(self):
        g = StructuredGeometry.wigner_seitz_pin_cell(
            r_fuel=0.9, r_clad=1.1, pitch=3.6,
        )
        assert g.domain_extent_cm == float(3.6 / np.sqrt(np.pi))


class TestPwrSlabHalfCell:
    def test_is_the_thickness_fold(self):
        g = StructuredGeometry.pwr_slab_half_cell(
            fuel_half=0.9, clad_thick=0.2, cool_thick=0.7,
        )
        assert g == StructuredGeometry.from_thicknesses(
            coord=_SLAB, thicknesses=(0.9, 0.2, 0.7), mat_ids=(2, 1, 0),
            face_roles=(MirrorDeck(), MirrorDeck()),
        )
