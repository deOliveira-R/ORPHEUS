r"""A deck law cannot be composed; a response can.

A deck law (``ReflectiveBoundary``, ``PeriodicBoundary``) states that the domain
is a fundamental domain of a symmetry group: the solution across the face is
the solution on this side, moved by the deck element. It has no amplitude
(``R = I``), so a scaled deck is not a law: ``0.3 * ReflectiveBoundary("x")``
was a second spelling of the partial specular wall ``AlbedoBoundary(0.3,
SpecularReturn("x"))``, the fossil ERR-094 was paid for. A SUM containing a
deck is either a scaled deck in disguise (``0.5 R + 0.5 R``) or not sub-Markov
(``R + W`` returns twice the outflow), so ``LawSum`` refuses one too. Both
constructors refuse a leaf whose geometry factor is not the identity, in
``__post_init__``, so every construction route reaches the refusal: the eight
dunders on ``BoundaryTraceLaw``, the composers' own dunders, the direct
constructors and ``dataclasses.replace``. A route that skips construction
(``object.__new__``, unpickling) skips it too.

**A deck nested in a sum that is then scaled** (``0.5 * (R + W)``) refuses at
the INNER node's construction, ``R + W``, before the scaling is evaluated, and
needs no subtree walk: every constructor checks its direct children, every
child is a tree that was itself admitted, so by induction no admitted tree
contains a deck at any depth. The induction holds only while the refusal is in
the constructor and not in the dunders; ``test_replace_reaches_the_refusal``
is the row that fails if it moves.

Responses still compose (``AlbedoBoundary`` with and without a closure,
``WhiteBoundary``, ``VacuumInflow``, ``ZeroFluxBoundary``, ``PrescribedInflow``)
and realize as before. The discriminating response is ``AlbedoBoundary(1.0,
SpecularReturn("x"))``: it permutes ordinates exactly as the mirror does and
must compose, so a refusal keyed on "permutes ordinates" rather than on the
geometry factor reds it.

Refusal contract (the message is the gate): ``TypeError``; the message carries
``repr(deck)``, the word ``symmetry`` and the phrase ``cannot be scaled or
mixed``.

First reds (``[M]`` 2026-10-01 at ``main`` ``0a5a23fa``, re-measured identical at ``88d30487``, ``.venv/bin/python -O
-m pytest tests/gates/geometry/test_deck_laws_do_not_compose.py``): every
refusal row (every route composes at HEAD); the response controls are green.
"""

from __future__ import annotations

import dataclasses
import re
from dataclasses import dataclass
from typing import Callable

import pytest

from orpheus.geometry.boundary import (
    AlbedoBoundary,
    BoundaryTraceLaw,
    LawScaled,
    LawSum,
    PeriodicBoundary,
    PrescribedInflow,
    ReflectiveBoundary,
    ScalarResponse,
    SelfPairedDeck,
    SpecularReturn,
    VacuumInflow,
    WhiteBoundary,
    ZeroFluxBoundary,
)

pytestmark = pytest.mark.foundation


@dataclass(frozen=True)
class _SyntheticMirror(BoundaryTraceLaw):
    """A law no shipped class spells, whose geometry factor is a mirror.

    Unregistered (no ``key=``), so it never enters the production registry.
    It is the row that tells a refusal keyed on the FACTOR (the ruling) from
    one keyed on ``isinstance(law, (ReflectiveBoundary, PeriodicBoundary))``,
    which lets a third deck law (a rotation, a glide) compose."""

    axis: str = "z"

    @property
    def geometry_map(self) -> SelfPairedDeck:
        return SelfPairedDeck.mirror(axis=self.axis)

    @property
    def response_kernel(self) -> ScalarResponse:
        return ScalarResponse(1.0)

    @property
    def kind(self) -> str:
        return "synthetic_mirror"


_DECKS: dict[str, BoundaryTraceLaw] = {
    "reflective": ReflectiveBoundary("x"),
    "periodic": PeriodicBoundary("x"),
    "synthetic-mirror": _SyntheticMirror(),
}

#: Responses: geometry factor the identity. The specular wall is first because
#: it is the one a wrong predicate refuses.
_RESPONSES: dict[str, BoundaryTraceLaw] = {
    "albedo-specular-1": AlbedoBoundary(1.0, SpecularReturn("x")),
    "albedo-specular-0.4": AlbedoBoundary(0.4, SpecularReturn("x")),
    "albedo-bare": AlbedoBoundary(0.5),
    "white": WhiteBoundary("x", +1, 0.5),
    "vacuum": VacuumInflow(),
    "zero-flux": ZeroFluxBoundary(),
    "prescribed-inflow": PrescribedInflow(),
}

_OTHER = WhiteBoundary("x", +1, 0.5)

#: Every spelling that puts ``law`` into a composition node, as a direct child
#: or through one of the composers' own dunders.
_ROUTES: dict[str, Callable[[BoundaryTraceLaw], object]] = {
    "k*law": lambda d: 0.5 * d,
    "law*k": lambda d: d * 0.5,
    "1*law": lambda d: 1.0 * d,
    "law/k": lambda d: d / 2.0,
    "-law": lambda d: -d,
    "law+r": lambda d: d + _OTHER,
    "r+law": lambda d: _OTHER + d,
    "law-r": lambda d: d - _OTHER,
    "r-law": lambda d: _OTHER - d,
    "law+law": lambda d: d + d,
    "LawScaled(k,law)": lambda d: LawScaled(0.5, d),
    "LawSum(law,r)": lambda d: LawSum(d, _OTHER),
    "LawSum(r,law)": lambda d: LawSum(_OTHER, d),
    # the composers' own dunders, with the law as the right operand
    "LawSum+law": lambda d: LawSum(_OTHER, _OTHER) + d,
    "law+LawSum": lambda d: d + LawSum(_OTHER, _OTHER),
    "LawSum-law": lambda d: LawSum(_OTHER, _OTHER) - d,
    "LawScaled+law": lambda d: (0.5 * _OTHER) + d,
    "LawScaled-law": lambda d: (0.5 * _OTHER) - d,
    "law-LawScaled": lambda d: d - (0.5 * _OTHER),
    # nested: the deck under a sum that is then scaled
    "k*(law+r)": lambda d: 0.5 * (d + _OTHER),
    "k*(r+k*law)": lambda d: 0.5 * (_OTHER + 0.3 * d),
}

_MIRROR_FACTOR = (
    "tests/gates/geometry/test_reflective_is_a_mirror.py::TestTheFactors::"
    "test_the_geometry_factor_is_the_mirror"
)
_WRAP_FACTOR = (
    "tests/gates/geometry/test_paired_deck.py::TestTheGuardsPartitionE3::"
    "test_exactly_one_type_admits_each_motion"
)

_REFUSAL = re.compile(r"symmetry.*cannot be scaled or mixed|cannot be scaled or mixed.*symmetry", re.S)


def _assert_refused(route: str, deck: BoundaryTraceLaw) -> None:
    with pytest.raises(TypeError) as info:
        _ROUTES[route](deck)
    message = str(info.value)
    assert repr(deck) in message, (
        f"{route}: the refusal does not name the law {deck!r}: {message!r}"
    )
    assert _REFUSAL.search(message), (
        f"{route}: the refusal does not say a symmetry cannot be scaled or "
        f"mixed: {message!r}"
    )


# ─────────────────────────────────────────────────────────────────────
# 1. Every route refuses a deck law
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.rests_on(_MIRROR_FACTOR, _WRAP_FACTOR)
@pytest.mark.parametrize("route", list(_ROUTES))
@pytest.mark.parametrize("deck", list(_DECKS))
def test_a_deck_law_cannot_be_composed(deck, route):
    _assert_refused(route, _DECKS[deck])


@pytest.mark.rests_on(_MIRROR_FACTOR, _WRAP_FACTOR)
@pytest.mark.parametrize("deck", list(_DECKS))
def test_replace_reaches_the_refusal(deck):
    """``dataclasses.replace`` re-runs ``__post_init__`` and calls no dunder:
    green only when the refusal is in the constructors, which is what makes
    the nesting induction (module docstring) hold."""
    law = _DECKS[deck]
    for build in (
        lambda: dataclasses.replace(LawSum(_OTHER, _OTHER), a=law),
        lambda: dataclasses.replace(LawSum(_OTHER, _OTHER), b=law),
        lambda: dataclasses.replace(LawScaled(0.5, _OTHER), inner=law),
    ):
        with pytest.raises(TypeError, match=re.escape(repr(law))):
            build()


@pytest.mark.rests_on(
    "tests/gates/geometry/test_deck_laws_do_not_compose.py::"
    "test_a_deck_law_cannot_be_composed",
)
def test_the_nested_refusal_is_the_inner_node_s():
    """``0.5 * (R + W)``: the refusal is raised while building ``R + W``, so
    the composition nodes on the raising stack are ``LawSum`` instances only;
    the outer ``LawScaled`` is never built. Mutation that reds it: ``LawSum``
    admitting the deck and ``LawScaled`` walking its subtree instead (a
    ``LawScaled`` frame on the stack)."""
    deck = _DECKS["reflective"]
    with pytest.raises(TypeError) as info:
        _ = 0.5 * (deck + _OTHER)
    under_construction = [
        type(f.locals["self"]) for f in info.traceback
        if isinstance(f.locals.get("self"), (LawSum, LawScaled))
    ]
    assert under_construction, (
        "the refusal was not raised while a composition node was being built"
    )
    assert set(under_construction) == {LawSum}, under_construction


# ─────────────────────────────────────────────────────────────────────
# 2. Responses still compose — the same routes, a response operand
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("route", list(_ROUTES))
@pytest.mark.parametrize("response", list(_RESPONSES))
def test_a_response_composes_on_every_route(response, route):
    """The discrimination row of every refusal above: the same route with a
    response operand builds a tree. A refusal keyed on "permutes ordinates"
    reds the two ``albedo-specular`` columns."""
    tree = _ROUTES[route](_RESPONSES[response])
    assert isinstance(tree, (LawSum, LawScaled))


# RETIRED LEGS (wave 2, slice S3; the demotion check of each). The last three tests realized a tree through the descriptor-tree
# walker and the two method realizers, all three retired with the string-tagged laws:
#
# * ``test_a_deck_leaf_alone_is_still_a_law`` (a bare deck leaf realizes, the matrix of the perfect specular wall): the deck is
#   realized through the table's ONE ``DeckRealizer``: in SN "the mirror deck realizes the unit specular image on each face" is the
#   deck row of ``tests/gates/sn/architecture/test_sn_realizer_table.py`` (the test whose docstring opens "The same law for the
#   deck") and, in diffusion, ``test_diffusion_translation_is_topological`` (the mirror deck realizes the identity, ``A = 1``).
# * ``test_a_response_mix_realizes_through_sn_as_before`` (a response mix is the pointwise weighted sum of its leaves): the
#   functor on the response network, gated per method (SN: the ``RealizedSum`` witness of ``tests/gates/sn/architecture/test_wave2_carry.py``; diffusion: law L-L).
# * ``test_a_response_mix_realizes_through_diffusion_as_before`` (``0.3 * 1 + 0.7 * 0.5``): ``test_a_mixed_wall_realizes_as_the_sum_of_scalars``
#   of ``tests/gates/diffusion/test_diffusion_realizer_table.py``, on the rows ``(0.3, 0.5)``, ``(0.7, 0.1)`` and ``(0.3, 0.7)``.
