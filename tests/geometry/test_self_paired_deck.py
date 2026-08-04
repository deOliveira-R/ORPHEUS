r"""The self-paired deck element's construction invariant (campaign step G5).

:class:`~orpheus.geometry.boundary.SelfPairedDeck` is :math:`G` for a face
paired with **itself** (Poincaré). Its whole content is a construction guard,
so these gates are that guard's positive and negative legs plus the properties
the guard makes derivable.

Why the guard is a CONJUNCTION, which is the thing most worth pinning
-------------------------------------------------------------------

The tempting spelling is "a self-paired face's pairing is an involution, so
guard on ``element_order() in (1, 2)``". That is the **converse** of the
insight and it is wrong. :math:`E(3)` has FOUR involution families:

===============  =====  ====  ===  ==========================
element          order   det   fix  pairing
===============  =====  ====  ===  ==========================
identity             1    +1    3  self-paired
reflection           2    −1    2  self-paired
**half-turn**        2    +1    1  face → **opposite**
**inversion**        2    −1    0  face → **opposite**
===============  =====  ====  ===  ==========================

The last two are exactly ``SpatialWrap``'s (deferred) job, so an order-only
guard admits the elements the type exists to exclude. The shipped guard is
``is_linear ∧ dim Fix ≥ d − 1``, which is strictly stronger, carries no
tolerance, and makes involution a **theorem** rather than a trusted premise:
a linear orthogonal :math:`Q` fixing a hyperplane pointwise is :math:`I` or a
reflection, hence :math:`Q^2 = I`.

Both clauses are load-bearing, and the witnesses below prove it in BOTH
directions — a glide passes the fixed-set clause and fails linearity; an
inversion passes linearity and fails the fixed-set clause. Neither clause is
implied by the other, so deleting either one admits a face-swapping element.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.boundary import (
    BoundaryGeometryMap,
    SelfPairedDeck,
    SpatialWrap,
)
from orpheus.geometry.transformation import RigidMotion

pytestmark = pytest.mark.foundation


def _glide() -> RigidMotion:
    """Reflection in x composed with a translation IN the mirror plane.

    ``fix = 2`` (it passes the fixed-set clause) but it is not linear, so only
    the linearity clause refuses it. The partner witness to ``inversion(3)``.
    """
    return RigidMotion.translation_by([0.0, 1.0, 0.0]) @ RigidMotion.reflection(
        normal=[1.0, 0.0, 0.0]
    )


def _half_turn() -> RigidMotion:
    """``diag(-1, -1, 1)`` — order 2, det +1, fixes only the z axis."""
    return RigidMotion.rotation_about_axis(axis=[0.0, 0.0, 1.0], angle=np.pi)


# ============================================================================
# G1a — POSITIVE: everything a self-paired face can legitimately carry
# ============================================================================


@pytest.mark.parametrize("dimension", [1, 2, 3])
def test_the_trivial_pairing_constructs_in_every_dimension(dimension: int) -> None:
    """``G = id`` is a self-pairing: it fixes everything, so it fixes the face."""
    deck = SelfPairedDeck.identity(dimension)
    assert deck.motion.dimension == dimension
    assert not deck.permutes_ordinates, (
        "the identity deck element relabels nothing, so it cannot permute the "
        "angular index"
    )


@pytest.mark.parametrize("axis", ["x", "y", "z"])
def test_a_coordinate_mirror_constructs_and_permutes(axis: str) -> None:
    """Every coordinate mirror is admissible and moves angle."""
    deck = SelfPairedDeck.mirror(axis=axis)
    assert deck.permutes_ordinates, (
        f"a mirror about {axis!r} exchanges the hemispheres, so realizing it "
        f"MUST permute the angular index — that is what makes G, not the "
        f"response, the carrier of the Gamma+ -> Gamma- crossing"
    )


def test_a_NON_coordinate_mirror_constructs() -> None:
    """`[M]` The plane :math:`x + y = 0` — which the retired ``axis: str``
    spelling could not name at all.

    This is the "tree constructs the object its type system cannot express"
    gap closed: a mirror is named by its NORMAL, and a normal is a vector, not
    a letter. Nothing in production installs a diagonal mirror today; the gate
    exists because the TYPE's admissible set is the claim, and a guard that
    silently refused a legitimate mirror would be a false constraint no
    consumer could see.
    """
    deck = SelfPairedDeck(RigidMotion.reflection(normal=[1.0, 1.0, 0.0]))
    assert deck.permutes_ordinates
    assert deck.domain_face("xmin") == "xmin"


# ============================================================================
# G1b/G1c — NEGATIVE, and the two clauses are separately load-bearing
# ============================================================================


def test_a_face_SWAPPING_involution_is_refused() -> None:
    r"""`[M]` The half-turn and the inversion are involutions and MUST still
    be refused — the whole reason the guard is not ``element_order``.

    Both satisfy :math:`g^2 = e`. Both map a face to its OPPOSITE, which is a
    genuine face pairing (``SpatialWrap``'s deferred job), not a self-pairing.
    An order-only guard admits both.
    """
    for name, motion in (("inversion", RigidMotion.inversion(3)),
                         ("half-turn", _half_turn())):
        assert motion.element_order() == 2, (
            f"{name} must be an involution for this gate to be the "
            f"discriminating one it claims to be"
        )
        with pytest.raises(ValueError, match="fix"):
            SelfPairedDeck(motion)


def test_an_affine_element_is_refused_even_when_it_fixes_a_plane() -> None:
    r"""`[M]` A glide fixes a 2-plane setwise yet is refused — the LINEARITY
    clause is the only thing that catches it.

    Partner witness to the previous test: that one passes linearity and fails
    the fixed-set clause; this one is the reverse. Neither clause implies the
    other, so both are load-bearing.
    """
    glide = _glide()
    assert glide.fixed_subspace_dimension == 2, (
        "the glide must pass the FIXED-SET clause, or it is not the witness "
        "this gate needs"
    )
    assert not glide.is_linear
    with pytest.raises(ValueError, match="LINEAR"):
        SelfPairedDeck(glide)


def test_a_mirror_with_an_OFFSET_is_refused() -> None:
    r"""`[M]` The carve's one genuinely invisible error class, closed by the
    TYPE because no gate could ever close it.

    ``on_directions`` drops the translation, so a mirror plane at the wrong
    POSITION produces a bit-identical ordinate permutation, a bit-identical
    realized image, and bit-identical snapshots — `vv-principles` Mode 12,
    designed-green at every tolerance in every regime. There is no measurement
    that distinguishes it, so the only available closure is to make it
    unspellable.
    """
    offset_mirror = RigidMotion.reflection(normal=[1.0, 0.0, 0.0], offset=2.5)
    plain = RigidMotion.reflection(normal=[1.0, 0.0, 0.0])
    probe = np.array([[0.3, -0.7, 0.6], [-1.0, 0.0, 0.0]])
    assert np.array_equal(
        offset_mirror.on_directions(probe), plain.on_directions(probe)
    ), (
        "if these ever differ, the offset is observable after all and this "
        "gate's rationale needs rewriting rather than its assertion"
    )
    with pytest.raises(ValueError, match="LINEAR"):
        SelfPairedDeck(offset_mirror)


@pytest.mark.parametrize(
    "name,motion",
    [
        ("quarter turn", RigidMotion.rotation_about_axis(
            axis=[0.0, 0.0, 1.0], angle=np.pi / 2)),
        ("pure translation", RigidMotion.translation_by([3.7, 0.0, 0.0])),
    ],
)
def test_non_involutions_are_refused(name: str, motion: RigidMotion) -> None:
    """A rotation and a translation are neither self-pairings nor involutions."""
    with pytest.raises(ValueError):
        SelfPairedDeck(motion)


def test_a_misnamed_axis_is_refused_at_CONSTRUCTION() -> None:
    """Not at realization, holding a quadrature it should never have reached.

    The retired ``SpecularMirror(axis="q")`` constructed happily and failed
    late with a bare ``ValueError`` raised from the quadrature's axis table —
    a diagnosis pointing at the wrong layer entirely.
    """
    with pytest.raises(ValueError, match="axis must be one of"):
        SelfPairedDeck.mirror(axis="q")


def test_an_axis_outside_the_dimension_is_refused() -> None:
    """A z-mirror has no meaning on a 1-D slab."""
    with pytest.raises(ValueError, match="out of range"):
        SelfPairedDeck.mirror(axis="z", dimension=1)


def test_a_non_motion_is_refused_with_a_type_error() -> None:
    """The deck element is the transformation, never a name for one."""
    with pytest.raises(TypeError, match="RigidMotion"):
        SelfPairedDeck("x")  # type: ignore[arg-type]


# ============================================================================
# G1e — the properties the guard makes DERIVABLE
# ============================================================================


@pytest.mark.parametrize("axis", ["x", "y", "z"])
def test_involution_is_a_derived_theorem_not_a_guard(axis: str) -> None:
    r"""Every admitted element satisfies :math:`g^2 = e` — asserted as a
    CONSEQUENCE of the guard, never assumed by it.

    This is the gate that would catch the guard being weakened back to
    something that admits a non-involution: the property is checked
    independently of the clauses that imply it.
    """
    for deck in (SelfPairedDeck.identity(), SelfPairedDeck.mirror(axis=axis)):
        motion = deck.motion
        assert motion.element_order() in (1, 2), (
            f"{motion!r} was admitted but is not an involution — the guard's "
            f"theorem (linear + fixes a hyperplane => Q^2 = I) is broken"
        )
        assert (motion @ motion) == RigidMotion.identity(motion.dimension)


#: Built INSIDE each test, never at import. A parametrize list that constructs
#: the SUT runs at COLLECTION time, so a guard regression that refuses a
#: legitimate element becomes a collection ERROR that takes the whole file
#: down — every other gate in it then reports nothing, and the diagnosis is
#: "1 error" instead of the named invariant that broke. `[M]` measured: an
#: over-strict-guard mutation produced exactly that, and no gate could isolate
#: it until these two were made lazy.
_ADMITTED = {
    "identity": lambda: SelfPairedDeck.identity(),
    "mirror-x": lambda: SelfPairedDeck.mirror(axis="x"),
}


@pytest.mark.parametrize("inhabitant", list(_ADMITTED))
def test_domain_face_is_self_paired_BY_CONSTRUCTION(inhabitant: str) -> None:
    """``domain_face(f) == f`` for every face and every admitted element.

    Self-pairing is the type's definition, so this cannot be a per-inhabitant
    choice that drifts from the guard. Swept over every canonical face name so
    a hard-coded single answer is not mistaken for the invariant.
    """
    deck = _ADMITTED[inhabitant]()
    for face in ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax"):
        assert deck.domain_face(face) == face


@pytest.mark.parametrize("inhabitant", list(_ADMITTED))
def test_it_conforms_to_the_geometry_map_protocol(inhabitant: str) -> None:
    """Runtime structural conformance, the same check the factor suite makes."""
    deck = _ADMITTED[inhabitant]()
    assert isinstance(deck, BoundaryGeometryMap)
    assert deck.is_adjointable


# ============================================================================
# G1f — CONTROL: the deferred half is untouched
# ============================================================================


def test_the_genuinely_paired_half_is_untouched() -> None:
    """``SpatialWrap`` still constructs and still answers the OPPOSITE face.

    A scope statement that cannot rot: if a later pass folds the periodic wrap
    into this type without building the surface-pair concept, this gate reds
    and says why.
    """
    wrap = SpatialWrap(axis="x")
    assert wrap.domain_face("xmin") == "xmax", (
        "SpatialWrap pairs two DISTINCT faces — that is exactly the case "
        "SelfPairedDeck refuses, and it must keep its own type until a "
        "surface-pair type exists"
    )
