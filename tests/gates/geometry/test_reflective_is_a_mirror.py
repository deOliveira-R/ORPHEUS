r"""``ReflectiveBoundary(axis)`` is the parameter-free deck mirror.

A reflective face is a symmetry plane of the problem: the domain is a
fundamental domain of the group the mirror generates, and the solution on the
far side is the mirror image of the solution on this side. Its geometry factor
is the Koopman operator of the mirror :math:`Q = I - 2\hat n\hat n^\top`, its
response factor is the identity, its source is zero, and its only datum is
which plane. A symmetry adds no physics, so it carries no amplitude: an albedo
on it was a category error (the fossil of the 2026-05 tensor-decomposition
framing, whose amplitude slot it was), and ERR-094 is what the fossil cost. A
partially specular wall is a response, spelled
``AlbedoBoundary(alpha, SpecularReturn(axis))``.

Plan of record: ``.claude/plans/boundary_law_ontology.md``, the fourth
exchange (why the albedo existed) and the tenth (this cleanup's scope).

The rows, each with its first red (``[M]`` 2026-10-01 at ``main`` ``0a5a23fa``, re-measured identical at ``88d30487``,
``.venv/bin/python -O -m pytest tests/gates/geometry/test_reflective_is_a_mirror.py``):

* the albedo is not a parameter: a second positional argument and an
  ``albedo=`` keyword are both a ``TypeError``, at any value including 1
  (red at HEAD: both construct);
* the law's only field is ``axis`` (red at HEAD: ``("axis", "albedo")``);
* the ``kind`` override retires: ``kind`` is the base's registry-key
  derivation, so it reads ``"reflective"`` on every axis because nothing else
  is spellable (red at HEAD: the override is in the class dict);
* the two amplitude checks ERR-094's review added retire with the field, and
  the certification still runs (red at HEAD: both in the class dict);
* the factors: ``G`` the mirror about ``axis``, ``R`` the unit scalar response
  that ``PeriodicBoundary``, the other deck law, also declares;
* value semantics: equality and hash on the axis; a mirror and a perfect
  specular wall, equal as matrices, are different values; pickle round-trips.

The factor, value-semantics and pickle rows are green at HEAD; their first red
is a mutation (named per row), since the state they reject is not spellable
there or here.
"""

from __future__ import annotations

import copy
import dataclasses
import pickle

import pytest

from orpheus.geometry.boundary import (
    AlbedoBoundary,
    BoundaryTraceLaw,
    PeriodicBoundary,
    ReflectiveBoundary,
    ScalarResponse,
    SelfPairedDeck,
    SpecularReturn,
    law_permutes_ordinates,
)
from orpheus.numerics.quadrature import Quadrature

pytestmark = pytest.mark.foundation

_AXES = ("x", "y", "z")

_MIRROR_MOTION = (
    "tests/gates/geometry/test_self_paired_deck.py::"
    "test_a_coordinate_mirror_constructs_and_permutes"
)
_INVOLUTION = (
    "tests/gates/geometry/test_bc_universal_invariants.py::"
    "TestReflectiveInvolutionInvariant::test_passes_for_gauss_legendre_x_axis"
)
_WALL_IS_THE_MIRROR_AT_ONE = (
    "tests/gates/geometry/test_reemission_closure.py::TestEquivalenceTheorems::"
    "test_specular_closure_equals_reflective"
)


# ─────────────────────────────────────────────────────────────────────
# 1. The albedo is not a parameter
# ─────────────────────────────────────────────────────────────────────


class TestTheAlbedoIsNotAParameter:
    """A second positional argument or an ``albedo=`` keyword is refused at
    construction, at every value. Including 1: an exemption for "the mirror's
    own value" would keep the slot and the category error with it."""

    @pytest.mark.parametrize("value", [0.5, 1.0, 0.0])
    def test_a_positional_albedo_is_a_type_error(self, value):
        with pytest.raises(TypeError, match="positional argument"):
            ReflectiveBoundary("x", value)  # type: ignore[call-arg]

    @pytest.mark.parametrize("value", [0.5, 1.0, 0.0])
    def test_an_albedo_keyword_is_a_type_error(self, value):
        with pytest.raises(TypeError, match="albedo"):
            ReflectiveBoundary(axis="x", albedo=value)  # type: ignore[call-arg]

    def test_the_registry_mint_refuses_it_too(self):
        """The polymorphic mint (``BoundaryTraceLaw.create``) reaches the same
        constructor; the parameter-free call returns the mirror."""
        assert BoundaryTraceLaw.create("reflective", axis="y") == ReflectiveBoundary("y")
        with pytest.raises(TypeError, match="albedo"):
            BoundaryTraceLaw.create("reflective", axis="x", albedo=1.0)

    def test_the_only_field_is_the_axis(self):
        """Equality, hash, ``dataclasses.replace`` and pickle all read the
        field tuple, so this one line is what makes "on the axis alone" true
        of each."""
        assert [f.name for f in dataclasses.fields(ReflectiveBoundary)] == ["axis"]

    def test_replace_cannot_reintroduce_it(self):
        with pytest.raises(TypeError, match="albedo"):
            dataclasses.replace(ReflectiveBoundary("x"), albedo=0.5)  # type: ignore[call-arg]


# ─────────────────────────────────────────────────────────────────────
# 2. kind, and the retired amplitude checks
# ─────────────────────────────────────────────────────────────────────


class TestKindIsTheRegistryKey:

    def test_the_kind_override_is_retired(self):
        """The override existed only to report ``"partial"`` for an albedo
        below 1. With the albedo gone the base's derivation (the registry key)
        is the whole answer, and an override would be a second spelling of
        it."""
        assert "kind" not in vars(ReflectiveBoundary)

    @pytest.mark.parametrize("axis", _AXES)
    def test_kind_is_reflective_on_every_axis(self, axis):
        """Mutation that reds it: an override returning anything but the key."""
        law = ReflectiveBoundary(axis)
        assert law.kind == "reflective"
        assert law.kind == type(law).key


class TestTheAmplitudeChecksRetire:
    """``assert_response_positive_if_declared`` and ``assert_submarkov`` were
    added at ERR-094's review to bound the albedo (ERR-043, ERR-046). Without
    the field they check nothing; ERR-043 and ERR-046 keep their catchers on
    ``WhiteBoundary`` and ``AlbedoBoundary``."""

    @pytest.mark.parametrize(
        "name", ["assert_submarkov", "assert_response_positive_if_declared"],
    )
    def test_the_override_is_gone(self, name):
        assert name not in vars(ReflectiveBoundary)

    @pytest.mark.parametrize(
        "axis, quadrature",
        [
            ("x", lambda: Quadrature.gauss_legendre(8)),
            ("x", lambda: Quadrature.lebedev(17)),
            ("y", lambda: Quadrature.lebedev(17)),
            ("z", lambda: Quadrature.lebedev(17)),
        ],
        ids=["x-GL8", "x-LEB17", "y-LEB17", "z-LEB17"],
    )
    @pytest.mark.rests_on(_INVOLUTION)
    def test_the_certification_still_runs_and_admits_the_mirror(self, axis, quadrature):
        """Positive leg: the mirror is realizable on every quadrature closed
        under it. Mutation that reds it: ``assert_realizable`` still calling a
        retired ``self.assert_submarkov()`` (an ``AttributeError``)."""
        ReflectiveBoundary(axis).assert_realizable(quadrature())

    def test_the_pairing_checks_survive(self):
        """The three pairing invariants (ERR-042, ERR-044, ERR-045) are the
        mirror's own and stay: guard against an over-deletion."""
        for name in (
            "assert_is_involutive",
            "assert_geometry_map_measure_preserving",
            "assert_reflection_maps_inflow_to_outflow",
            "assert_realizable",
        ):
            assert name in vars(ReflectiveBoundary), name


# ─────────────────────────────────────────────────────────────────────
# 3. The affine form's factors
# ─────────────────────────────────────────────────────────────────────


class TestTheFactors:

    @pytest.mark.rests_on(_MIRROR_MOTION)
    @pytest.mark.parametrize("axis", _AXES)
    def test_the_geometry_factor_is_the_mirror(self, axis):
        """Mutation that reds it: a mirror about a fixed axis."""
        law = ReflectiveBoundary(axis)
        assert law.geometry_map == SelfPairedDeck.mirror(axis=axis)
        assert law.geometry_map.permutes_ordinates
        assert law_permutes_ordinates(law)

    @pytest.mark.parametrize("axis", _AXES)
    def test_the_response_factor_is_the_identity(self, axis):
        """``R = I``, the same unit scalar response the other deck law
        (``PeriodicBoundary``) declares. Mutation that reds it: any amplitude
        other than 1 (the deleted field's default read through a constant)."""
        response = ReflectiveBoundary(axis).response_kernel
        assert response == ScalarResponse(1.0)
        assert response == PeriodicBoundary(axis).response_kernel
        assert response.amplitude == 1.0
        assert not response.is_zero


# ─────────────────────────────────────────────────────────────────────
# 4. Value semantics
# ─────────────────────────────────────────────────────────────────────


class TestValueSemantics:

    @pytest.mark.parametrize("axis", _AXES)
    def test_equal_on_the_axis(self, axis):
        a, b = ReflectiveBoundary(axis), ReflectiveBoundary(axis=axis)
        assert a == b
        assert hash(a) == hash(b)

    def test_different_axes_are_different_values(self):
        """Through the container, not by comparing hashes (a hash collision is
        legal): three axes, three members."""
        laws = [ReflectiveBoundary(a) for a in _AXES]
        assert len(set(laws)) == len(laws) == 3
        assert ReflectiveBoundary("x") != ReflectiveBoundary("y")

    @pytest.mark.rests_on(_WALL_IS_THE_MIRROR_AT_ONE)
    @pytest.mark.parametrize("axis", _AXES)
    def test_a_mirror_is_not_a_perfect_specular_wall(self, axis):
        """Equal as realized matrices (the SN realizer's shared deck kernel;
        ``test_reemission_closure.py``'s alpha = 1 equivalence pins that) and
        different values: one is a symmetry of the domain, the other a wall's
        constitutive response. Mutation that reds it: an ``__eq__`` comparing
        realized factors (``law_permutes_ordinates`` and the amplitude)."""
        mirror = ReflectiveBoundary(axis)
        wall = AlbedoBoundary(1.0, SpecularReturn(axis))
        assert mirror != wall
        assert wall != mirror
        assert len({mirror, wall}) == 2

    @pytest.mark.parametrize("axis", _AXES)
    def test_pickle_and_copy_round_trip(self, axis):
        law = ReflectiveBoundary(axis)
        for twin in (
            pickle.loads(pickle.dumps(law)), copy.copy(law), copy.deepcopy(law),
        ):
            assert type(twin) is ReflectiveBoundary
            assert twin == law
            assert hash(twin) == hash(law)
            assert twin.axis == axis
