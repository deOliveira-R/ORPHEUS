r"""The geometry constructors that name each boundary point's law (P1 step 3a).

``StructuredGeometry.slab(breakpoints, mat_ids, *, left, right)``,
``.cylinder(..., *, outer, inner=None)`` and ``.sphere(..., *, outer,
inner=None)`` name the face each law belongs to, where the bare constructor
takes a position-coded tuple; ``.uniform_boundary(coord, breakpoints,
mat_ids, law)`` puts one law on the whole boundary whatever its point count;
``.from_homogeneous(width, boundary)`` is the finite slab a test uses for the
infinite medium (spec S4.1). Ruling 3 of 2026-09-29
(``.claude/plans/reference_cache.md``, "P1 step 3 opened").

Every constructor here is a SPELLING of the bare constructor, never a second
definition of the geometry: each one equals the bare constructor's value field
by field (the population is ``dataclasses.fields(StructuredGeometry)``, so a
field added later is compared without editing this file), and each refusal is
the bare constructor's own, reached through its one boundary check (the
message is compared as a string, which is the route witness: a second door
would word its own). Every test is ``foundation``; claim kind THEOREM unless
stated.
"""
from __future__ import annotations

import dataclasses
import math
import re
from collections.abc import Callable
from fractions import Fraction

import pytest

from orpheus.geometry import BC, CoordSystem, StructuredGeometry

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/geometry/test_named_face_constructors.py"
_GEOMETRY = "tests/gates/geometry/test_structured_geometry.py"
_S2_1 = f"{_GEOMETRY}::TestBreakpointLaws::test_refusal"
_S2_2 = f"{_GEOMETRY}::TestTheBoundaryIsDerived::test_boundary_points"

_SLAB = CoordSystem.CARTESIAN
_CYLINDER = CoordSystem.CYLINDRICAL
_SPHERE = CoordSystem.SPHERICAL
_RADIAL = {"cylinder": _CYLINDER, "sphere": _SPHERE}

#: Two different laws, so a constructor that swaps the faces is seen.
_LAW_A = BC.vacuum
_LAW_B = BC.reflective


def _typed_law():
    """A typed law (not a tag): a prescribed inflow, carried by identity."""
    from orpheus.geometry.boundary import ConstantInflowSource, PrescribedInflow

    return PrescribedInflow(source=ConstantInflowSource(value=2.5))


def _bare(coord, breakpoints, mat_ids, boundaries) -> StructuredGeometry:
    return StructuredGeometry(
        coord=coord, breakpoints=breakpoints, mat_ids=mat_ids, boundaries=boundaries,
    )


def _assert_same_value(named: StructuredGeometry, bare: StructuredGeometry) -> None:
    """Field by field over the TYPE's fields, laws by identity."""
    if type(named) is not StructuredGeometry:
        raise AssertionError(f"the named constructor built a {type(named).__name__}")
    for f in dataclasses.fields(StructuredGeometry):
        mine, theirs = getattr(named, f.name), getattr(bare, f.name)
        if f.name == "boundaries":
            if len(mine) != len(theirs) or any(
                m is not t for m, t in zip(mine, theirs, strict=True)
            ):
                raise AssertionError(f"boundaries {mine!r} are not {theirs!r} (by identity)")
        elif mine != theirs:
            raise AssertionError(f"field {f.name}: {mine!r} != {theirs!r}")


def _message(build: Callable[[], object], error: type[Exception]) -> str:
    with pytest.raises(error) as caught:
        build()
    return str(caught.value)


def _assert_disjoint(table) -> None:
    """Each row's fragment is in its own message and in no other row's."""
    messages = {case: _message(build, error) for case, build, error, _ in table}
    for case, _, _, fragment in table:
        for other, message in messages.items():
            other_fragment = next(f for c, _, _, f in table if c == other)
            if other_fragment == fragment:
                continue
            if fragment in message:
                raise AssertionError(
                    f"the fragment of {case!r} ({fragment!r}) is in {other!r}: {message!r}"
                )


# ─────────────────────────────────────────────────────────────────────
# slab
# ─────────────────────────────────────────────────────────────────────


class TestSlab:
    """``slab(breakpoints, mat_ids, *, left, right)``."""

    @pytest.mark.rests_on(_S2_1, _S2_2)
    @pytest.mark.parametrize("breakpoints", [(0.0, 1.0, 3.0), (-1.0, 2.0)], ids=["origin", "negative-origin"])
    def test_equals_the_bare_constructor(self, breakpoints):
        mat_ids = tuple(range(len(breakpoints) - 1))
        law = _typed_law()
        _assert_same_value(
            StructuredGeometry.slab(breakpoints, mat_ids, left=_LAW_A, right=law),
            _bare(_SLAB, breakpoints, mat_ids, (_LAW_A, law)),
        )
        _assert_same_value(
            StructuredGeometry.slab(breakpoints, mat_ids, left=_LAW_B, right=_LAW_A),
            _bare(_SLAB, breakpoints, mat_ids, (_LAW_B, _LAW_A)),
        )


_SLAB_REFUSALS = [
    ("no-left", lambda: StructuredGeometry.slab((0.0, 1.0), (0,), right=_LAW_A),  # type: ignore[arg-type, call-arg]  # a refusal input
     TypeError, "'left'"),
    ("no-right", lambda: StructuredGeometry.slab((0.0, 1.0), (0,), left=_LAW_A),  # type: ignore[arg-type, call-arg]  # a refusal input
     TypeError, "'right'"),
    ("positional-laws", lambda: StructuredGeometry.slab((0.0, 1.0), (0,), _LAW_A, _LAW_A),  # type: ignore[arg-type, call-arg]  # a refusal input
     TypeError, "positional argument"),
    ("none-law", lambda: StructuredGeometry.slab((0.0, 1.0), (0,), left=None, right=_LAW_A),  # type: ignore[arg-type, call-arg]  # a refusal input
     TypeError, "None is not a boundary law"),
    ("breakpoints", lambda: StructuredGeometry.slab((0.0, 0.0), (0,), left=_LAW_A, right=_LAW_A),
     ValueError, "strictly increasing"),
]


def test_negative_zero_is_one_position():
    """qa F4: ``-0.0`` is canonicalised to ``+0.0`` by the shared parser, so
    a slab from ``-0.0`` is the slab from ``0.0``, bit for bit."""
    g = StructuredGeometry.slab((-0.0, 1.0), (0,), left=_LAW_A, right=_LAW_A)
    assert math.copysign(1.0, g.breakpoints[0]) == 1.0
    assert g == StructuredGeometry.slab((0.0, 1.0), (0,), left=_LAW_A, right=_LAW_A)


class TestSlabRefusals:
    @pytest.mark.parametrize(
        "build, error, fragment", [pytest.param(b, e, f, id=c) for c, b, e, f in _SLAB_REFUSALS],
    )
    def test_refusal(self, build, error, fragment):
        with pytest.raises(error, match=re.escape(fragment)):
            build()

    @pytest.mark.rests_on(f"{_HERE}::TestSlabRefusals::test_refusal")
    def test_refusal_fragments_are_disjoint(self):
        _assert_disjoint(_SLAB_REFUSALS)

    def test_the_refusals_are_the_bare_constructors(self):
        """The route witness: the same defect, spelled through the bare
        constructor, gives the same message character for character."""
        assert _message(
            lambda: StructuredGeometry.slab((0.0, 1.0), (0,), left=None, right=_LAW_A),  # type: ignore[arg-type]  # a refusal input
            TypeError,
        ) == _message(lambda: _bare(_SLAB, (0.0, 1.0), (0,), (None, _LAW_A)), TypeError)
        assert _message(
            lambda: StructuredGeometry.slab((0.0, 0.0), (0,), left=_LAW_A, right=_LAW_A), ValueError,
        ) == _message(lambda: _bare(_SLAB, (0.0, 0.0), (0,), (_LAW_A, _LAW_A)), ValueError)


# ─────────────────────────────────────────────────────────────────────
# cylinder, sphere
# ─────────────────────────────────────────────────────────────────────


class TestRadial:
    """``cylinder`` and ``sphere``: ``outer`` always, ``inner`` exactly when hollow."""

    @pytest.mark.rests_on(_S2_1, _S2_2)
    @pytest.mark.parametrize("name", _RADIAL)
    def test_a_solid_body_equals_the_bare_constructor(self, name):
        law = _typed_law()
        build = getattr(StructuredGeometry, name)
        _assert_same_value(
            build((0.0, 0.5, 2.0), (0, 1), outer=law),
            _bare(_RADIAL[name], (0.0, 0.5, 2.0), (0, 1), (law,)),
        )

    @pytest.mark.rests_on(_S2_1, _S2_2)
    @pytest.mark.parametrize("name", _RADIAL)
    def test_a_hollow_body_equals_the_bare_constructor(self, name):
        """The inner law goes first (the bare order), whichever keyword is
        written first."""
        build = getattr(StructuredGeometry, name)
        _assert_same_value(
            build((0.5, 2.0), (0,), outer=_LAW_A, inner=_LAW_B),
            _bare(_RADIAL[name], (0.5, 2.0), (0,), (_LAW_B, _LAW_A)),
        )
        _assert_same_value(
            build((0.5, 2.0), (0,), inner=_LAW_A, outer=_LAW_B),
            _bare(_RADIAL[name], (0.5, 2.0), (0,), (_LAW_A, _LAW_B)),
        )

    @pytest.mark.parametrize("name", _RADIAL)
    def test_an_absent_inner_is_never_a_stored_none(self, name):
        g = getattr(StructuredGeometry, name)((0.0, 2.0), (0,), outer=_LAW_A)
        assert g.boundaries == (_LAW_A,)
        assert None not in g.boundaries


def _radial_refusals(name: str):
    build = getattr(StructuredGeometry, name)
    return [
        ("no-outer", lambda: build((0.0, 2.0), (0,)), TypeError, "'outer'"),
        ("positional-law", lambda: build((0.0, 2.0), (0,), _LAW_A), TypeError, "positional argument"),
        ("inner-on-a-solid-body", lambda: build((0.0, 2.0), (0,), outer=_LAW_A, inner=_LAW_B),
         ValueError, "the centre r = 0"),
        ("hollow-without-inner", lambda: build((0.5, 2.0), (0,), outer=_LAW_A),
         ValueError, "inner surface, which needs its own law"),
        ("none-outer", lambda: build((0.0, 2.0), (0,), outer=None), TypeError, "None is not a boundary law"),
        ("negative-radius", lambda: build((-0.5, 2.0), (0,), outer=_LAW_A),
         ValueError, "starts at r_0 >= 0"),
    ]


class TestRadialRefusals:
    @pytest.mark.parametrize("name", _RADIAL)
    @pytest.mark.parametrize(
        "case", [c for c, *_ in _radial_refusals("sphere")],
    )
    def test_refusal(self, name, case):
        _, build, error, fragment = next(r for r in _radial_refusals(name) if r[0] == case)
        with pytest.raises(error, match=re.escape(fragment)):
            build()

    @pytest.mark.rests_on(f"{_HERE}::TestRadialRefusals::test_refusal")
    @pytest.mark.parametrize("name", _RADIAL)
    def test_refusal_fragments_are_disjoint(self, name):
        _assert_disjoint(_radial_refusals(name))

    @pytest.mark.parametrize("name", _RADIAL)
    def test_an_explicit_none_inner_is_the_absent_keyword(self, name):
        """``inner=None`` is the absent keyword: on a hollow body it is
        refused exactly as an omitted ``inner`` is, and on a solid body it
        builds the one-law geometry."""
        build = getattr(StructuredGeometry, name)
        assert _message(lambda: build((0.5, 2.0), (0,), outer=_LAW_A, inner=None), ValueError) == (
            _message(lambda: build((0.5, 2.0), (0,), outer=_LAW_A), ValueError)
        )
        assert build((0.0, 2.0), (0,), outer=_LAW_A, inner=None).boundaries == (_LAW_A,)

    @pytest.mark.parametrize("name", _RADIAL)
    def test_the_refusals_are_the_bare_constructors(self, name):
        """The route witness: one boundary check, one message."""
        build, coord = getattr(StructuredGeometry, name), _RADIAL[name]
        assert _message(lambda: build((0.0, 2.0), (0,), outer=_LAW_A, inner=_LAW_B), ValueError) == (
            _message(lambda: _bare(coord, (0.0, 2.0), (0,), (_LAW_B, _LAW_A)), ValueError)
        )
        assert _message(lambda: build((0.5, 2.0), (0,), outer=_LAW_A), ValueError) == (
            _message(lambda: _bare(coord, (0.5, 2.0), (0,), (_LAW_A,)), ValueError)
        )


# ─────────────────────────────────────────────────────────────────────
# uniform_boundary
# ─────────────────────────────────────────────────────────────────────

#: (coord, breakpoints, the number of boundary points, S2.2's table).
_UNIFORM_TABLE = [
    (_SLAB, (0.0, 1.0, 2.0), 2),
    (_SLAB, (-1.0, 2.0), 2),
    (_CYLINDER, (0.0, 1.0, 2.0), 1),
    (_CYLINDER, (0.5, 2.0), 2),
    (_SPHERE, (0.0, 1.0, 2.0), 1),
    (_SPHERE, (0.5, 2.0), 2),
]


class TestUniformBoundary:
    """``uniform_boundary(coord, breakpoints, mat_ids, law)``: ``law`` at
    every boundary point, two on a slab or a hollow body, one on a solid
    cylinder or sphere."""

    @pytest.mark.rests_on(_S2_1, _S2_2)
    @pytest.mark.parametrize(
        "coord, breakpoints, n_points", _UNIFORM_TABLE,
        ids=[f"{c.name.lower()}-r0={b[0]}" for c, b, _ in _UNIFORM_TABLE],
    )
    def test_equals_the_bare_constructor(self, coord, breakpoints, n_points):
        law = _typed_law()
        mat_ids = tuple(range(len(breakpoints) - 1))
        _assert_same_value(
            StructuredGeometry.uniform_boundary(coord, breakpoints, mat_ids, law),
            _bare(coord, breakpoints, mat_ids, (law,) * n_points),
        )

    @pytest.mark.rests_on(f"{_HERE}::TestUniformBoundary::test_equals_the_bare_constructor")
    @pytest.mark.parametrize("coord", (_CYLINDER, _SPHERE), ids=lambda c: c.name.lower())
    @pytest.mark.parametrize(
        "r_0", [Fraction(1, 10**400), -0.0], ids=["fraction-below-the-smallest-float", "negative-zero"],
    )
    def test_hollowness_is_the_geometry_s_decision(self, coord, r_0):
        """The route witness for ONE definition of hollowness: an ``r_0``
        whose raw value is not ``0.0`` but whose parsed float is (a
        ``Fraction`` below the smallest subnormal), or the reverse spelling
        (``-0.0``), gives the solid geometry the bare constructor parses. A
        second hollowness test on the RAW input refused the Fraction row
        ([M] 2026-09-29, before the fix)."""
        g = StructuredGeometry.uniform_boundary(coord, (r_0, 1.0), (0,), _LAW_A)
        assert not g.is_hollow
        assert g.boundaries == (_LAW_A,)
        assert math.copysign(1.0, g.breakpoints[0]) == 1.0

    @pytest.mark.parametrize("breakpoints", [(-1.0, 1.0), (0.0, 1.0)], ids=["negative-origin", "origin"])
    def test_a_string_coordinate_gets_the_retirement_refusal(self, breakpoints):
        """elegance F8: the coordinate is parsed FIRST, so a retired string
        tag gets the keyed retirement refusal, the bare constructor's own,
        whatever the breakpoints (with ``(-1, 1)`` it used to raise an
        ``AttributeError`` from the breakpoint parser)."""
        with pytest.raises(TypeError, match=re.escape("kind tags ('SLB', 'CYL', 'SPH') are retired")):
            StructuredGeometry.uniform_boundary("SLB", breakpoints, (0,), _LAW_A)  # type: ignore[arg-type]  # a refusal input
        assert _message(
            lambda: StructuredGeometry.uniform_boundary("SLB", breakpoints, (0,), _LAW_A), TypeError,  # type: ignore[arg-type]  # a refusal input
        ) == _message(lambda: _bare("SLB", breakpoints, (0,), (_LAW_A, _LAW_A)), TypeError)

    def test_a_none_law_is_refused_by_the_bare_check(self):
        with pytest.raises(TypeError, match="None is not a boundary law"):
            StructuredGeometry.uniform_boundary(_SPHERE, (0.0, 1.0), (0,), None)  # type: ignore[arg-type]  # a refusal input


# ─────────────────────────────────────────────────────────────────────
# S4.1 — from_homogeneous
# ─────────────────────────────────────────────────────────────────────


class TestFromHomogeneous:
    r"""S4.1: ``from_homogeneous(width, boundary)`` is the slab ``[0, width]``
    of material 0 with ``boundary`` on both faces.

    Mutation witness: pass ``boundary`` to one face only (the other a fixed
    law) → the identity leg reds on the second face.
    """

    @pytest.mark.rests_on(f"{_HERE}::TestSlab::test_equals_the_bare_constructor")
    @pytest.mark.parametrize("width", [2.5, 3, 1e-3], ids=["float", "int", "thin"])
    def test_equals_the_named_slab(self, width):
        law = _typed_law()
        g = StructuredGeometry.from_homogeneous(width, law)
        _assert_same_value(g, StructuredGeometry.slab((0.0, width), (0,), left=law, right=law))
        assert g.breakpoints == (0.0, float(width))
        assert g.mat_ids == (0,)
        assert g.coord is _SLAB
        assert g.boundaries[0] is law and g.boundaries[1] is law

    @pytest.mark.parametrize("law", [BC.reflective, BC.vacuum, BC("albedo", {"albedo": 0.3})])
    def test_both_faces_carry_the_law(self, law):
        g = StructuredGeometry.from_homogeneous(1.0, law)
        assert g.boundaries == (law, law)


_HOMOGENEOUS_REFUSALS = [
    ("negative-zero", lambda: StructuredGeometry.from_homogeneous(-0.0, _LAW_B),
     ValueError, "the width is positive and finite"),
    ("zero", lambda: StructuredGeometry.from_homogeneous(0.0, _LAW_B),
     ValueError, "the width is positive and finite"),
    ("negative", lambda: StructuredGeometry.from_homogeneous(-1.0, _LAW_B),
     ValueError, "the width is positive and finite"),
    ("inf", lambda: StructuredGeometry.from_homogeneous(math.inf, _LAW_B),
     ValueError, "the width is positive and finite"),
    ("nan", lambda: StructuredGeometry.from_homogeneous(math.nan, _LAW_B),
     ValueError, "the width is positive and finite"),
    ("string", lambda: StructuredGeometry.from_homogeneous("1.0", _LAW_B),  # type: ignore[arg-type, call-arg]  # a refusal input
     TypeError, "must be a real number"),
    ("bool", lambda: StructuredGeometry.from_homogeneous(True, _LAW_B),
     TypeError, "must be a real number"),
    ("none-law", lambda: StructuredGeometry.from_homogeneous(1.0, None),  # type: ignore[arg-type, call-arg]  # a refusal input
     TypeError, "None is not a boundary law"),
]


class TestFromHomogeneousRefusals:
    @pytest.mark.parametrize(
        "build, error, fragment",
        [pytest.param(b, e, f, id=c) for c, b, e, f in _HOMOGENEOUS_REFUSALS],
    )
    def test_refusal(self, build, error, fragment):
        with pytest.raises(error, match=re.escape(fragment)):
            build()

    @pytest.mark.rests_on(f"{_HERE}::TestFromHomogeneousRefusals::test_refusal")
    def test_refusal_fragments_are_disjoint(self):
        _assert_disjoint(_HOMOGENEOUS_REFUSALS)
