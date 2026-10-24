"""Gates for the characteristic reference's walls (:mod:`orpheus.derivations.continuous.characteristic.walls`).

P1 step (b), first rung, of the characteristic-reference campaign
(``.claude/plans/characteristic_reference_architecture.md``, "P1 step (b),
first rung: API sketch" and its "Ruled 2026-10-06" line; the ledger entry
"2026-10-06, the user, on P1 step (b)'s first rung"), re-posed at wave 2, slice S5
(``.claude/plans/wave2_boundary_carve.md`` 2.9 and 6.5): the wall table no longer reads a typed boundary
law through its two factors; it reads a face ROLE (a deck, or an integer tag the specification's ``exterior``
states). Verification spec: ``scratch/characteristic_architecture/p1_verification_spec.md``; this file is the
rung's own rows (the role table, the refusals); the equivalence of the reading with the HEAD reading and with a hand
rule over 19 bodies is ``test_walls_of_reads_decks_and_amplitudes.py``.

The claims gated, each a THEOREM about the reading of a face role:

* **The role table.** Every role a reference serves is read into a wall ``(breakpoint, specular, diffuse,
  partner)``. The expected walls are written by hand from the physics, not from the reader: a mirror deck returns
  along the reflected line, a matte wall returns isotropically, a wrap deck re-enters at the opposite wall, vacuum
  returns nothing.
* **The refusals.** What the reader still refuses is the reader's own (an unstated tag; a directly built ``Walls``
  whose wrap does not close or sits on a radial chart). The refusals that moved are not here: a curved face's wrap and
  an unpaired wrap are geometry's (:class:`~orpheus.geometry.deck.DeckError`, ``test_curved_wraps_are_deck_errors.py``),
  a prescribed inflow, a zero flux and a combined wall have no spelling in ``WallAmplitudes`` (the specification
  refuses them where it is built, ``tests/gates/specification/test_exterior.py``), and an amplitude outside ``[0, 1]``
  is ``WallAmplitudes``'s parse refusal, a ``ValueError`` since ``Wall`` composes it (R-8b; it was a
  ``NotImplementedError`` while ``Wall`` held its own bounds).

Retired with the typed-law reader at S5, each with its reason (the classification owed in the commit body): the
string-tag registry rows (``test_the_registry_admits_exactly_the_kinds_this_file_twins``,
``test_a_tag_reads_the_same_walls_as_the_law_it_names``, ``test_the_white_tag_without_an_albedo_is_a_full_return``,
``test_the_registry_hands_the_white_law_the_walls_outward_sign`` and the ``unadmitted_tag_*``, ``*_with_a_parameter``
and ``partial_without_albedo`` refusal rows) are DEAD, their subject, the registry that parsed a tag into a typed
law, being deleted; ``albedo_with_unstated_shape``, ``zero_flux`` and ``inflow_source`` are DEAD BY TYPE (no field of
the specification can state them); the factor-pair cells (``test_a_served_factor_pair_reads_its_wall`` and
``test_a_quotient_deck_returns_everything_and_a_shape_sits_on_the_identity``) are DEAD with ``_wall_of``'s two-factor
match, their subject (a deck factor and a response factor that need not agree) not existing in a reader of face
roles; the ``wrap_not_paired*`` and ``wrap_on_*`` rows MOVED to geometry.

Declared blindness: a wall reads only amplitudes and the partner, so an exterior tag's number is invisible here.

Rows at ``l0`` with ``verifies("characteristic-closure")`` (the label on
``docs/theory/references/characteristic.rst``): the role table, which defines what the closure's amplitudes ARE.
Every other row is ``foundation``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.characteristic import Wall, Walls
from orpheus.geometry.deck import MirrorDeck, TranslationDeck
from orpheus.geometry.structured_geometry import StructuredGeometry
from orpheus.specification.exterior import WallAmplitudes

_HERE = "tests/gates/derivations/test_characteristic_walls.py::"

_SLAB_BP = (-0.7, 0.3, 1.1, 2.0)        # n = 3: walls at breakpoints 0 and 3
_HOLLOW_BP = (0.4, 1.1, 2.0)            # n = 2: walls at 0 and 2
_SOLID_BP = (0.0, 0.3, 1.1, 2.0)        # n = 3: one wall, at 3


def _slab(left, right) -> StructuredGeometry:
    return StructuredGeometry.slab(_SLAB_BP, (0, 1, 2), left=left, right=right)


# ── the role table ───────────────────────────────────────────────────────
#
# (id, left role, right role, exterior, expected (specular, diffuse, partner) at the left
# wall 0 and at the right wall 3). Written from the physics of each wall.
_MIRROR, _WRAP = MirrorDeck(), TranslationDeck()
_ROLE_ROWS = [
    ("vacuum", 0, 0, {0: WallAmplitudes.VACUUM}, (0.0, 0.0, 0), (0.0, 0.0, 3)),
    ("mirror_deck", _MIRROR, _MIRROR, {}, (1.0, 0.0, 0), (1.0, 0.0, 3)),
    ("specular_tags", 0, 1, {0: WallAmplitudes(0.37, 0.0), 1: WallAmplitudes(0.81, 0.0)}, (0.37, 0.0, 0), (0.81, 0.0, 3)),
    ("diffuse_tags", 0, 1, {0: WallAmplitudes(0.0, 0.37), 1: WallAmplitudes(0.0, 0.81)}, (0.0, 0.37, 0), (0.0, 0.81, 3)),
    ("white_full_return", 0, 1, {0: WallAmplitudes(0.0, 0.6), 1: WallAmplitudes.WHITE}, (0.0, 0.6, 0), (0.0, 1.0, 3)),
    ("wrap_deck", _WRAP, _WRAP, {}, (1.0, 0.0, 3), (1.0, 0.0, 0)),
    ("mixed_mirror_and_matte", _MIRROR, 0, {0: WallAmplitudes(0.0, 0.25)}, (1.0, 0.0, 0), (0.0, 0.25, 3)),
    ("shared_tag", 0, 0, {0: WallAmplitudes(0.0, 0.5)}, (0.0, 0.5, 0), (0.0, 0.5, 3)),
]


@pytest.mark.l0
@pytest.mark.verifies("characteristic-closure")
@pytest.mark.parametrize(("left", "right", "exterior", "at_left", "at_right"),
                         [r[1:] for r in _ROLE_ROWS], ids=[r[0] for r in _ROLE_ROWS])
def test_every_law_reads_as_the_wall_its_physics_names(left, right, exterior, at_left, at_right) -> None:
    """A slab's two walls, each face role read, against the hand-written wall.

    First reds: (a) the wrap's partner read as the wall itself (a wrap deck as
    a mirror): the ``wrap_deck`` row; (b) a diffuse amplitude read as
    specular: ``diffuse_tags``, ``white_full_return``; (c) the vacuum tag returning a
    default amplitude 1: the ``vacuum`` row; (d) the amplitudes of the wrong tag read
    (tag 0 for both faces): ``specular_tags``, ``diffuse_tags``, ``white_full_return``.
    """
    walls = Walls.of(_slab(left, right), exterior)
    assert walls.n_regions == 3
    expected = (Wall(0, WallAmplitudes(*at_left[:2]), at_left[2]), Wall(3, WallAmplitudes(*at_right[:2]), at_right[2]))
    assert walls.walls == expected, f"{walls.walls} != {expected}"
    bp = np.array([0, 3])
    np.testing.assert_array_equal(walls.specular_at(bp), [at_left[0], at_right[0]])
    np.testing.assert_array_equal(walls.diffuse_at(bp), [at_left[1], at_right[1]])
    np.testing.assert_array_equal(walls.partner_at(bp), [at_left[2], at_right[2]])


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_every_law_reads_as_the_wall_its_physics_names")
def test_a_body_has_a_wall_at_each_boundary_point_and_none_elsewhere() -> None:
    """Walls sit at breakpoints 0 and n of a hollow body and at n alone on a solid one, inner first.

    An interior breakpoint, a solid body's centre (breakpoint 0) and the
    absent-transit code n + 1 are no wall, and every lookup refuses them with
    ``IndexError``. First red: the per-wall lookup indexed by position instead
    of by breakpoint (the solid body's only wall would answer at 0).
    """
    hollow = Walls.of(
        StructuredGeometry.sphere(_HOLLOW_BP, (0, 1), inner=0, outer=1),
        {0: WallAmplitudes(0.3, 0.0), 1: WallAmplitudes(0.6, 0.0)},
    )
    assert [w.breakpoint for w in hollow.walls] == [0, 2]
    np.testing.assert_array_equal(hollow.specular_at(np.array([[2, 0], [0, 2]])), [[0.6, 0.3], [0.3, 0.6]])
    solid = Walls.of(StructuredGeometry.sphere(_SOLID_BP, (0, 1, 2), outer=0), {0: WallAmplitudes.MIRROR})
    assert solid.walls == (Wall(3, WallAmplitudes.MIRROR, 3),)
    for walls, non_walls in ((hollow, (1, 3)), (solid, (0, 1, 2, 4))):
        for k in non_walls:
            for lookup in (walls.specular_at, walls.diffuse_at, walls.partner_at):
                with pytest.raises(IndexError):
                    lookup(np.array([k]))


# ── the refusals ─────────────────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.parametrize(("make", "exterior", "unstated"), [
    (lambda: _slab(0, 0), {}, 0),
    (lambda: _slab(0, 1), {0: WallAmplitudes.VACUUM}, 1),
    (lambda: StructuredGeometry.sphere(_HOLLOW_BP, (0, 1), inner=2, outer=0), {0: WallAmplitudes.VACUUM}, 2),
], ids=["slab_no_entry", "slab_second_tag_missing", "hollow_inner_tag_missing"])
def test_a_tag_the_exterior_does_not_state_is_refused_by_number(make, exterior, unstated: int) -> None:
    """``Walls.of`` on a tag with no amplitudes is a ``ValueError`` naming the tag.

    The reader takes a bare mapping and admits it as the exact cover through the specification's own door
    (:func:`~orpheus.specification.exterior.admit_exterior`), so a missing tag is refused by number here as it is at
    construction of a specification (it would otherwise be a ``KeyError`` naming nothing). First red: the admission removed."""
    with pytest.raises(ValueError, match=rf"tag\(s\) \[{unstated}\]"):
        Walls.of(make(), exterior)


@pytest.mark.foundation
@pytest.mark.parametrize(("specular", "diffuse"), [(1.5, 0.0), (0.0, -0.2)], ids=["amplitude_above_one", "amplitude_negative"])
def test_an_unphysical_amplitude_is_refused_before_a_wall_exists(specular: float, diffuse: float) -> None:
    """An amplitude outside ``[0, 1]`` is ``WallAmplitudes``'s parse refusal (R-8b: a ``ValueError``, no longer a ``NotImplementedError``).

    Re-posed from the rows ``amplitude_above_one`` and ``amplitude_negative`` of the typed-law table: the refusal
    moved from ``Walls.of`` (which read a law holding the amplitude) to the value that holds it. The accepting
    control is the same amplitude within range."""
    WallAmplitudes(min(max(specular, 0.0), 1.0), min(max(diffuse, 0.0), 1.0))
    with pytest.raises(ValueError, match="not a physical wall") as caught:
        WallAmplitudes(specular, diffuse)
    assert type(caught.value) is ValueError


@pytest.mark.foundation
def test_a_directly_built_walls_refuses_what_walls_of_refuses() -> None:
    """The invariants live on ``Walls`` itself: a wall off the ends, a repeated wall, an unpaired or radial wrap."""
    from orpheus.geometry.chart import Chart
    from orpheus.geometry.coord import CoordSystem

    slab, sphere = Chart(CoordSystem.CARTESIAN), Chart(CoordSystem.SPHERICAL)
    vacuum, mirror = WallAmplitudes.VACUUM, WallAmplitudes.MIRROR
    with pytest.raises(ValueError, match="walls sit at distinct breakpoints"):
        Walls((Wall(1, vacuum, 1),), 3, slab)
    with pytest.raises(ValueError, match="walls sit at distinct breakpoints"):
        Walls((Wall(0, vacuum, 0), Wall(0, mirror, 0)), 3, slab)
    with pytest.raises(NotImplementedError, match="whose partner does not wrap back"):
        Walls((Wall(0, mirror, 3), Wall(3, mirror, 3)), 3, slab)
    with pytest.raises(NotImplementedError, match="periodic wrap on a cylinder or sphere"):
        Walls((Wall(0, mirror, 2), Wall(2, mirror, 0)), 2, sphere)
    assert Walls((Wall(0, mirror, 3), Wall(3, mirror, 0)), 3, slab).partner_at(np.array([0, 3])).tolist() == [3, 0]
    with pytest.raises(ValueError, match="not a physical wall"):
        Wall(3, WallAmplitudes(0.6, 1.0 + 1e-9), 3)
