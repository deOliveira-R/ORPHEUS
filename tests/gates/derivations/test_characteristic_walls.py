"""Gates for the characteristic reference's walls (:mod:`orpheus.derivations.continuous.characteristic.walls`).

P1 step (b), first rung, of the characteristic-reference campaign
(``.claude/plans/characteristic_reference_architecture.md``, "P1 step (b),
first rung: API sketch" and its "Ruled 2026-10-06" line; the ledger entry
"2026-10-06, the user, on P1 step (b)'s first rung"). Verification spec:
``scratch/characteristic_architecture/p1_verification_spec.md``; this file is
the rung's own rows (the factor table, the tag registry, the refusals).

The claims gated, each a THEOREM about the reading of a boundary law:

* **The factor table.** Every shipped law is read through its two factors
  (the deck, ``geometry_map``, and the response, ``response_kernel``) into a
  wall ``(breakpoint, specular, diffuse, partner)``. The expected walls are
  written by hand from the law's physics, not from its factors: a mirror
  returns along the reflected line, a matte wall returns isotropically, a wrap
  re-enters at the opposite wall, vacuum returns nothing.
* **The tag registry.** Each admitted tag kind reads the SAME walls as the
  typed law it names (the drift gate between the two spellings); the gate
  ranges over the registry itself, so a kind added without its typed twin
  here reds.
* **The refusals.** Each a ``NotImplementedError`` keyed to its own fragment;
  the fragments are pairwise disjoint over the refused inputs (one refusal
  cannot pass for another).

Declared blindness: a wall reads only amplitudes and the partner, so the
registry's axis and outward sign (the white law's hemisphere) are invisible
here; they matter to no reading of this reference.

Rows to move from ``foundation`` to ``l0`` when the archivist mints the label
``characteristic-closure`` on the rewritten page: the factor table and the
tag-registry drift gate (they define what the closure's amplitudes ARE).
"""
from __future__ import annotations

from typing import cast

import numpy as np
import pytest

from orpheus.derivations.continuous.characteristic import Wall, Walls
from orpheus.derivations.continuous.characteristic import walls as walls_module
from orpheus.geometry.boundary import (
    BC,
    AlbedoBoundary,
    BoundaryTraceLaw,
    ConstantInflowSource,
    IsotropicReturn,
    PeriodicBoundary,
    PrescribedInflow,
    ReflectiveBoundary,
    SpecularReturn,
    VacuumInflow,
    WhiteBoundary,
    ZeroFluxBoundary,
)
from orpheus.geometry.structured_geometry import StructuredGeometry

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/derivations/test_characteristic_walls.py::"
_FACTORS = "tests/gates/geometry/test_boundary_factors.py::"

_SLAB_BP = (-0.7, 0.3, 1.1, 2.0)        # n = 3: walls at breakpoints 0 and 3
_HOLLOW_BP = (0.4, 1.1, 2.0)            # n = 2: walls at 0 and 2
_SOLID_BP = (0.0, 0.3, 1.1, 2.0)        # n = 3: one wall, at 3


def _slab(left, right) -> StructuredGeometry:
    return StructuredGeometry.slab(_SLAB_BP, (0, 1, 2), left=left, right=right)


# ── the factor table ─────────────────────────────────────────────────────
#
# (id, left law, right law, expected (specular, diffuse, partner) at the left
# wall 0 and at the right wall 3). Written from the physics of each law.
_VAC = (0.0, 0.0)
_FACTOR_ROWS = [
    ("vacuum", VacuumInflow(), VacuumInflow(), (0.0, 0.0, 0), (0.0, 0.0, 3)),
    ("vacuum_tag", BC.vacuum, BC.vacuum, (0.0, 0.0, 0), (0.0, 0.0, 3)),
    ("reflective", ReflectiveBoundary(axis="x"), ReflectiveBoundary(axis="x"), (1.0, 0.0, 0), (1.0, 0.0, 3)),
    ("albedo_zero_is_vacuum", AlbedoBoundary(0.0), AlbedoBoundary(0.0), (0.0, 0.0, 0), (0.0, 0.0, 3)),
    ("albedo_specular", AlbedoBoundary(0.37, SpecularReturn(axis="x")),
     AlbedoBoundary(0.81, SpecularReturn(axis="x")), (0.37, 0.0, 0), (0.81, 0.0, 3)),
    ("albedo_isotropic", AlbedoBoundary(0.37, IsotropicReturn(axis="x", outward_sign=-1)),
     AlbedoBoundary(0.81, IsotropicReturn(axis="x", outward_sign=+1)), (0.0, 0.37, 0), (0.0, 0.81, 3)),
    ("white", WhiteBoundary(axis="x", outward_sign=-1, albedo=0.6),
     WhiteBoundary(axis="x", outward_sign=+1), (0.0, 0.6, 0), (0.0, 1.0, 3)),
    ("periodic", PeriodicBoundary(axis="x"), PeriodicBoundary(axis="x"), (1.0, 0.0, 3), (1.0, 0.0, 0)),
    ("prescribed_inflow_no_source_is_vacuum", PrescribedInflow(), PrescribedInflow(), (0.0, 0.0, 0), (0.0, 0.0, 3)),
    ("mixed_mirror_and_matte", ReflectiveBoundary(axis="x"), WhiteBoundary(axis="x", outward_sign=+1, albedo=0.25),
     (1.0, 0.0, 0), (0.0, 0.25, 3)),
]


@pytest.mark.parametrize(("left", "right", "at_left", "at_right"),
                         [r[1:] for r in _FACTOR_ROWS], ids=[r[0] for r in _FACTOR_ROWS])
@pytest.mark.rests_on(_FACTORS + "test_every_production_law_states_both_factors")
def test_every_law_reads_as_the_wall_its_physics_names(left, right, at_left, at_right) -> None:
    """A slab's two walls, each law read through its factors, against the hand-written wall.

    First reds: (a) the wrap's partner read as the wall itself (a periodic face
    as a mirror): the ``periodic`` row; (b) a Lambertian response read as
    specular: ``albedo_isotropic``, ``white``; (c) the vacuum arm returning a
    default amplitude 1: every vacuum row; (d) the white tag's albedo ignored:
    ``white`` (0.6 at the left wall).
    """
    walls = Walls.of(_slab(left, right))
    assert walls.n_regions == 3
    expected = (Wall(0, *at_left), Wall(3, *at_right))
    assert walls.walls == expected, f"{walls.walls} != {expected}"
    bp = np.array([0, 3])
    np.testing.assert_array_equal(walls.specular_at(bp), [at_left[0], at_right[0]])
    np.testing.assert_array_equal(walls.diffuse_at(bp), [at_left[1], at_right[1]])
    np.testing.assert_array_equal(walls.partner_at(bp), [at_left[2], at_right[2]])


@pytest.mark.rests_on(_HERE + "test_every_law_reads_as_the_wall_its_physics_names")
def test_a_body_has_a_wall_at_each_boundary_point_and_none_elsewhere() -> None:
    """Walls sit at breakpoints 0 and n of a hollow body and at n alone on a solid one, inner first.

    An interior breakpoint, a solid body's centre (breakpoint 0) and the
    absent-transit code n + 1 are no wall, and every lookup refuses them with
    ``IndexError``. First red: the per-wall lookup indexed by position instead
    of by breakpoint (the solid body's only wall would answer at 0).
    """
    hollow = Walls.of(StructuredGeometry.sphere(
        _HOLLOW_BP, (0, 1), inner=AlbedoBoundary(0.3, SpecularReturn(axis="x")),
        outer=AlbedoBoundary(0.6, SpecularReturn(axis="x"))))
    assert [w.breakpoint for w in hollow.walls] == [0, 2]
    np.testing.assert_array_equal(hollow.specular_at(np.array([[2, 0], [0, 2]])), [[0.6, 0.3], [0.3, 0.6]])
    solid = Walls.of(StructuredGeometry.sphere(_SOLID_BP, (0, 1, 2), outer=BC.reflective))
    assert solid.walls == (Wall(3, 1.0, 0.0, 3),)
    for walls, non_walls in ((hollow, (1, 3)), (solid, (0, 1, 2, 4))):
        for k in non_walls:
            for lookup in (walls.specular_at, walls.diffuse_at, walls.partner_at):
                with pytest.raises(IndexError):
                    lookup(np.array([k]))


# ── the tag registry ─────────────────────────────────────────────────────
#
# Each admitted kind with a parametrised tag and the typed law it must name at
# a wall of outward sign -1 (left) and +1 (right). Non-default parameters, so
# a registry that drops one cannot agree by default.
_TWINS = {
    "vacuum": (BC.vacuum, lambda sign: VacuumInflow()),
    "reflective": (BC.reflective, lambda sign: ReflectiveBoundary(axis="x")),
    "partial": (BC("partial", {"albedo": 0.37}), lambda sign: AlbedoBoundary(0.37, SpecularReturn(axis="x"))),
    "white": (BC.white, lambda sign: WhiteBoundary(axis="x", outward_sign=sign)),
    "periodic": (BC("periodic"), lambda sign: PeriodicBoundary(axis="x")),
}


def test_the_registry_admits_exactly_the_kinds_this_file_twins() -> None:
    """The drift gate ranges over the registry, not over a hand list: a new kind reds here until it is twinned."""
    assert set(walls_module.TAG_REGISTRY) == set(_TWINS)


@pytest.mark.parametrize("kind", sorted(_TWINS))
@pytest.mark.rests_on(_HERE + "test_every_law_reads_as_the_wall_its_physics_names",
                      _HERE + "test_the_registry_admits_exactly_the_kinds_this_file_twins")
def test_a_tag_reads_the_same_walls_as_the_law_it_names(kind: str) -> None:
    """``Walls.of`` of a tagged slab equals ``Walls.of`` of the slab carrying the typed twin, per kind.

    First reds: ``partial`` mapped to an albedo with no re-emission shape (it
    refuses); ``partial`` mapped to a mirror (1.0 against 0.37); ``white``
    mapped to a mirror (specular 1 against diffuse 1). The white tag takes no
    parameter (ruled 2026-10-06: a partial white wall is spelled
    ``WhiteBoundary(albedo=a)``; the parameterised tag is a refusal row).
    """
    tag, law = _TWINS[kind]
    by_tag = Walls.of(_slab(tag, tag))
    by_law = Walls.of(_slab(law(-1), law(+1)))
    assert by_tag == by_law, f"{kind}: {by_tag.walls} != {by_law.walls}"


def test_the_white_tag_without_an_albedo_is_a_full_return() -> None:
    """``BC("white")`` returns all of its outflow, as ``WhiteBoundary``'s default."""
    walls = Walls.of(_slab(BC.white, BC.white))
    assert walls.walls == (Wall(0, 0.0, 1.0, 0), Wall(3, 0.0, 1.0, 3))


@pytest.mark.parametrize("sign", [-1, +1])
def test_the_registry_hands_the_white_law_the_walls_outward_sign(sign: int) -> None:
    """The white tag's law carries the wall's outward sign, by content equality (a ``Walls`` cannot see it).

    The one registry datum no wall reading observes (the hemisphere the
    Lambertian averages over); first red: the sign dropped or negated.
    """
    assert walls_module.TAG_REGISTRY["white"](BC.white, sign) == WhiteBoundary(axis="x", outward_sign=sign)
    assert walls_module.TAG_REGISTRY["white"](BC.white, sign) != WhiteBoundary(axis="x", outward_sign=-sign)


# ── the refusals ─────────────────────────────────────────────────────────

_PREFIX = "the characteristic reference does not serve a"
_REFUSALS = [  # (id, geometry factory, fragment)
    ("albedo_with_unstated_shape", lambda: _slab(AlbedoBoundary(0.5), BC.vacuum), "unstated re-emission shape"),
    ("zero_flux", lambda: _slab(ZeroFluxBoundary(), BC.vacuum), "unstated re-emission shape"),
    ("inflow_source", lambda: _slab(PrescribedInflow(ConstantInflowSource(1.0)), BC.vacuum), "with an inflow source"),
    ("wrap_not_paired", lambda: _slab(PeriodicBoundary(axis="x"), BC.vacuum), "whose partner does not wrap back"),
    ("wrap_not_paired_right", lambda: _slab(BC.reflective, BC("periodic")), "whose partner does not wrap back"),
    ("wrap_on_hollow_sphere", lambda: StructuredGeometry.sphere(
        _HOLLOW_BP, (0, 1), inner=PeriodicBoundary(axis="x"), outer=PeriodicBoundary(axis="x")),
     "periodic wrap on a cylinder or sphere"),
    ("wrap_on_hollow_cylinder", lambda: StructuredGeometry.cylinder(
        _HOLLOW_BP, (0, 1), inner=BC("periodic"), outer=BC("periodic")), "periodic wrap on a cylinder or sphere"),
    ("wrap_on_solid_sphere", lambda: StructuredGeometry.sphere(_SOLID_BP, (0, 1, 2), outer=PeriodicBoundary(axis="x")),
     "periodic wrap on a cylinder or sphere"),
    ("unadmitted_tag_albedo", lambda: _slab(BC("albedo", {"albedo": 0.5}), BC.vacuum), "the reference admits the tag kinds"),
    ("unadmitted_tag_marshak", lambda: _slab(BC.vacuum, BC("marshak")), "the reference admits the tag kinds"),
    ("partial_without_albedo", lambda: _slab(BC("partial"), BC.vacuum), "takes exactly the parameters"),
    ("partial_with_an_extra_parameter", lambda: _slab(BC("partial", {"albedo": 0.5, "angle": 0.1}), BC.vacuum),
     "takes exactly the parameters"),
    ("vacuum_with_a_parameter", lambda: _slab(BC("vacuum", {"albedo": 0.5}), BC.vacuum), "takes exactly the parameters"),
    ("white_with_an_albedo", lambda: _slab(BC("white", {"albedo": 0.6}), BC.vacuum), "takes exactly the parameters"),
    ("reflective_with_a_parameter", lambda: _slab(BC.vacuum, BC("reflective", {"albedo": 1.0})),
     "takes exactly the parameters"),
    ("periodic_with_a_parameter", lambda: _slab(BC("periodic", {"shift": 0.0}), BC("periodic")),
     "takes exactly the parameters"),
    ("amplitude_above_one", lambda: _slab(AlbedoBoundary(1.5, SpecularReturn(axis="x")), BC.vacuum),
     "not a physical wall"),
    ("amplitude_negative", lambda: _slab(BC.vacuum, WhiteBoundary(axis="x", outward_sign=+1, albedo=-0.2)),
     "not a physical wall"),
]
_FRAGMENTS = sorted({r[2] for r in _REFUSALS})


@pytest.mark.parametrize(("make", "fragment"), [r[1:] for r in _REFUSALS], ids=[r[0] for r in _REFUSALS])
def test_an_unserved_wall_is_refused_by_its_own_fragment(make, fragment: str) -> None:
    """Each refused declaration raises ``NotImplementedError`` with its fragment and no other's.

    The geometry itself admits every one of these (the refusal is the
    reference's, not the geometry's). First reds, one per guard: the source
    check removed (``inflow_source`` reads as vacuum); the ``is_zero`` guard of
    the vacuum arm dropped (``albedo_with_unstated_shape`` reads specular 0.5);
    the wrap-back check removed; the radial-wrap check removed; the amplitude range
    check removed; an unknown tag falling through to vacuum; a tag's
    undeclared parameter dropped instead of refused (the white albedo, which
    production drops silently, is the case that motivated the rule).
    """
    geometry = make()
    with pytest.raises(NotImplementedError) as caught:
        Walls.of(geometry)
    message = str(caught.value)
    assert message.startswith(_PREFIX), message
    assert fragment in message, message
    others = [f for f in _FRAGMENTS if f != fragment and f in message]
    assert not others, f"the refusal also carries {others}: {message}"


def test_an_unadmitted_tag_is_refused_naming_every_admitted_kind() -> None:
    """The refusal lists the admitted kinds, read from the registry."""
    with pytest.raises(NotImplementedError) as caught:
        Walls.of(_slab(BC("marshak"), BC.vacuum))
    for kind in walls_module.TAG_REGISTRY:
        assert repr(kind) in str(caught.value), (kind, str(caught.value))


# ── the table's cells no shipped law reaches ─────────────────────────────


class _FactorLaw:
    """A stand-in law: just the two factors and the source the reader consults (no registration side effect)."""

    def __init__(self, deck, response) -> None:
        from orpheus.geometry.boundary import NoSource

        self.geometry_map, self.response_kernel, self.source = deck, response, NoSource()


def _cells():
    from orpheus.geometry.boundary import (
        LambertianReemission, PairedDeck, ScalarResponse, SelfPairedDeck, SpecularReemission,
    )
    mirror, wrap, identity = SelfPairedDeck.mirror(axis="x"), PairedDeck.wrap(axis="x"), SelfPairedDeck.identity()
    served = [
        ("identity_specular_reemission", identity, SpecularReemission(alpha=0.4), Wall(3, 0.4, 0.0, 3)),
        ("identity_lambertian", identity, LambertianReemission(alpha=0.7), Wall(3, 0.0, 0.7, 3)),
        ("mirror_full_return", mirror, ScalarResponse(1.0), Wall(3, 1.0, 0.0, 3)),
        ("wrap_full_return", wrap, ScalarResponse(1.0), Wall(3, 1.0, 0.0, 0)),
    ]
    refused = [
        ("mirror_partial_scalar", mirror, ScalarResponse(0.5)),
        ("mirror_specular_reemission", mirror, SpecularReemission(alpha=0.5)),
        ("mirror_lambertian", mirror, LambertianReemission(alpha=0.5)),
        ("wrap_partial_scalar", wrap, ScalarResponse(0.5)),
        ("wrap_lambertian", wrap, LambertianReemission(alpha=1.0)),
        ("identity_scalar_unstated", identity, ScalarResponse(0.5)),
    ]
    return served, refused


_SERVED, _REFUSED_CELLS = _cells()


@pytest.mark.parametrize(("deck", "response", "expected"), [c[1:] for c in _SERVED], ids=[c[0] for c in _SERVED])
@pytest.mark.rests_on(_HERE + "test_every_law_reads_as_the_wall_its_physics_names")
def test_a_served_factor_pair_reads_its_wall(deck, response, expected) -> None:
    """The reader's served cells on factor pairs built directly (the wall at breakpoint 3 of a 3-region body)."""
    assert walls_module._wall_of(cast(BoundaryTraceLaw, _FactorLaw(deck, response)), 3, opposite=0) == expected


@pytest.mark.parametrize(("deck", "response"), [c[1:] for c in _REFUSED_CELLS], ids=[c[0] for c in _REFUSED_CELLS])
def test_a_quotient_deck_returns_everything_and_a_shape_sits_on_the_identity(deck, response) -> None:
    """A mirror or a wrap admits only ScalarResponse(1); a re-emission shape only the identity deck.

    No shipped law reaches these cells; the first reds are a mirror read with
    its scalar amplitude (``mirror_partial_scalar`` returns 0.5 along the
    reflected line) and a shape read across a quotient deck.
    """
    with pytest.raises(NotImplementedError, match="unstated re-emission shape"):
        walls_module._wall_of(cast(BoundaryTraceLaw, _FactorLaw(deck, response)), 3, opposite=0)


def test_a_directly_built_walls_refuses_what_walls_of_refuses() -> None:
    """The invariants live on ``Walls`` itself: a wall off the ends, a repeated wall, an unpaired or radial wrap."""
    from orpheus.geometry.chart import Chart
    from orpheus.geometry.coord import CoordSystem

    slab, sphere = Chart(CoordSystem.CARTESIAN), Chart(CoordSystem.SPHERICAL)
    with pytest.raises(ValueError, match="walls sit at distinct breakpoints"):
        Walls((Wall(1, 0.0, 0.0, 1),), 3, slab)
    with pytest.raises(ValueError, match="walls sit at distinct breakpoints"):
        Walls((Wall(0, 0.0, 0.0, 0), Wall(0, 1.0, 0.0, 0)), 3, slab)
    with pytest.raises(NotImplementedError, match="whose partner does not wrap back"):
        Walls((Wall(0, 1.0, 0.0, 3), Wall(3, 1.0, 0.0, 3)), 3, slab)
    with pytest.raises(NotImplementedError, match="periodic wrap on a cylinder or sphere"):
        Walls((Wall(0, 1.0, 0.0, 2), Wall(2, 1.0, 0.0, 0)), 2, sphere)
    assert Walls((Wall(0, 1.0, 0.0, 3), Wall(3, 1.0, 0.0, 0)), 3, slab).partner_at(np.array([0, 3])).tolist() == [3, 0]
    with pytest.raises(NotImplementedError, match="not a physical wall"):
        Wall(3, 0.6, 1.0 + 1e-9, 3)
