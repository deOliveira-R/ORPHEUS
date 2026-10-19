r"""Content identity of the ``Axis`` family (#405 P1 step 5, the 5a move).

``Axis`` carried a hand-written ``_identity_key`` tuple and a private
``_structural_bytes`` encoder (scalars by ``repr``, so ``1 != 1.0``) until
step 5 moved it onto the one content encoder,
:class:`~orpheus.numerics.content.ContentIdentity`. Its parts are its
dataclass fields: ``generator`` is provenance, not identity (the module
docstring's slot table, ruling of 2026-08-20), so it is declared
``field(compare=False)``, the encoder's one spelling of "not content". The
unit (the coordinates of the constant function, numerics-untangling rulings
39 and 41) is content on a MODAL axis only, and by the type: the kind is
``BasisKind.NODAL`` or ``Modal(unit)``, so a modal axis's ``kind`` part carries
its unit and a nodal one's carries nothing new (every nodal digest is the one
it was before the unit existed). The modal entries' ``kind`` part therefore
has a leg that moves only the unit.

The roster also holds :class:`~orpheus.numerics.content.FrozenMapping`, the
one frozen mapping value the content types hold (``Materials``, ``BC``
params, ``FaceLaws``).

The S5.2, S5.3 and pickle rows run the shared checkers of
``tests/gates/_content_identity_helpers.py`` over this roster; S5.9's
population row (``tests/gates/numerics/test_content_identity.py``) reads it,
so an ``Axis`` subclass added without an entry reds there.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from orpheus.numerics.axis import Axis, BasisKind, EnergyAxis, HarmonicAxis, LegendreAxis, Modal
from orpheus.numerics.content import FrozenMapping
from tests.gates._content_identity_helpers import (
    Entry,
    check_equal_pair,
    check_perturbation,
    check_pickle,
    check_population,
    leg,
    pair_ids,
    param_id,
    perturbation_ids,
)

pytestmark = pytest.mark.foundation

_W = np.array([0.5, 1.5])


def _axis(**kw) -> Axis:
    base: dict[str, Any] = dict(label="spatial", shape=(2,), weights=_W, kind=BasisKind.NODAL)
    return Axis(**(base | kw))


def _energy(**kw) -> EnergyAxis:
    base: dict[str, Any] = dict(label="energy", shape=(2,), kind=BasisKind.NODAL, edges=np.array([1.0, 2.0, 3.0]))
    return EnergyAxis(**(base | kw))


_E0 = np.array([1.0, 0.0])


def _harmonic(**kw) -> HarmonicAxis:
    base: dict[str, Any] = dict(label="angular", shape=(2,), weights=_W, kind=Modal(_E0))
    return HarmonicAxis(**(base | kw))


def _legendre(**kw) -> LegendreAxis:
    base: dict[str, Any] = dict(
        label="angular", shape=(2,), weights=_W, kind=Modal(_E0), spent_axis="z",
    )
    return LegendreAxis(**(base | kw))


def _read_only(array: np.ndarray) -> np.ndarray:
    array.flags.writeable = False
    return array


def _forged(value, name: str, part):
    """``value`` with one field overwritten past its constructor's laws."""
    object.__setattr__(value, name, part)
    return value


def _base_legs(build) -> dict:
    """A leg per base part, for any member's builder (the shape leg moves the
    weights with it, since the weights live on the shape; on a modal member it
    moves the kind too, since the unit lives on the shape). The kind leg mints
    the other kind lawfully; a modal member has a second kind leg that moves
    only the unit."""
    if build().kind is BasisKind.NODAL:
        return {
            "label": (leg("another role", lambda: build(label="other")),),
            "shape": (leg("another shape", lambda: build(shape=(3,), weights=np.array([0.5, 1.5, 2.0])), "weights"),),
            "kind": (leg("the other kind", lambda: build(kind=Modal(_E0))),),
            "weights": (leg("one weight", lambda: build(weights=np.array([0.5, 2.5]))),),
        }
    return {
        "label": (leg("another role", lambda: build(label="other")),),
        "shape": (leg(
            "another shape",
            lambda: build(shape=(3,), weights=np.array([0.5, 1.5, 2.0]), kind=Modal(np.array([1.0, 0.0, 0.0]))),
            "weights", "kind",
        ),),
        "kind": (
            leg("the other kind", lambda: build(kind=BasisKind.NODAL)),
            leg("the moments of one, not e0", lambda: build(kind=Modal(np.array([2.0, 0.0])))),
        ),
        "weights": (leg("one weight", lambda: build(weights=np.array([0.5, 2.5]))),),
    }


_AXIS = Entry(
    cls=Axis, base=_axis, parts=("label", "shape", "weights", "kind"),
    perturb=_base_legs(_axis),
    pairs=(
        ("two builds", _axis, _axis),
        ("weights as ints vs floats", lambda: _axis(weights=np.array([2, 3])),
         lambda: _axis(weights=np.array([2.0, 3.0]))),
        ("a generator is not content", _axis, lambda: _axis(generator=object())),
    ),
)

_ENERGY = Entry(
    cls=EnergyAxis, base=_energy, parts=("label", "shape", "weights", "kind", "edges"),
    perturb=_base_legs(_energy) | {
        "shape": (leg("another shape", lambda: _energy(shape=(3,), edges=np.array([1.0, 2.0, 3.0, 4.0])), "edges"),),
        # An EnergyAxis is NODAL with the counting measure BY CONSTRUCTION, so
        # no constructible instance differs in ``kind`` or ``weights``; the
        # legs forge the field on a built instance to witness that the
        # encoder still reads the part (it is content, inherited from Axis).
        "kind": (leg("kind forged", lambda: _forged(_energy(), "kind", Modal(np.array([1.0, 0.0])))),),
        "weights": (leg("weights forged", lambda: _forged(_energy(), "weights", _read_only(np.array([0.5, 1.5])))),),
        "edges": (leg("a group edge", lambda: _energy(edges=np.array([1.0, 2.5, 3.0]))),),
    },
    pairs=(("two builds", _energy, _energy),),
)

_HARMONIC = Entry(
    cls=HarmonicAxis, base=_harmonic, parts=("label", "shape", "weights", "kind"),
    perturb=_base_legs(_harmonic),
    pairs=(("two builds", _harmonic, _harmonic),),
)

_LEGENDRE = Entry(
    cls=LegendreAxis, base=_legendre,
    parts=("label", "shape", "weights", "kind", "spent_axis"),
    perturb=_base_legs(_legendre) | {
        "spent_axis": (leg("the other pole", lambda: _legendre(spent_axis="x")),),
    },
    pairs=(("two builds", _legendre, _legendre),),
)

_MODAL = Entry(
    cls=Modal, base=lambda: Modal(_E0), parts=("unit",),
    perturb={"unit": (leg("the moments of one, not e0", lambda: Modal(np.array([2.0, 0.0]))),)},
    pairs=(
        ("two builds", lambda: Modal(_E0), lambda: Modal(_E0)),
        ("a negative zero is zero", lambda: Modal(np.array([1.0, 0.0])), lambda: Modal(np.array([1.0, -0.0]))),
    ),
)

_FROZEN_MAPPING = Entry(
    cls=FrozenMapping, base=lambda: FrozenMapping({"a": 1.0, "b": 2.0}), parts=("items",),
    fields_are_parts=False,  # not a dataclass: its content is its items
    perturb={
        "items": (
            leg("a value", lambda: FrozenMapping({"a": 1.0, "b": 2.5})),
            leg("a key", lambda: FrozenMapping({"a": 1.0, "c": 2.0})),
            leg("an item dropped", lambda: FrozenMapping({"a": 1.0})),
        ),
    },
    pairs=(
        ("two builds", lambda: FrozenMapping({"a": 1.0, "b": 2.0}), lambda: FrozenMapping({"a": 1.0, "b": 2.0})),
        ("item order", lambda: FrozenMapping({"a": 1.0, "b": 2.0}), lambda: FrozenMapping({"b": 2.0, "a": 1.0})),
        ("an int vs a float value", lambda: FrozenMapping({"a": 1, "b": 2.0}), lambda: FrozenMapping({"a": 1.0, "b": 2.0})),
    ),
)

ROSTER: tuple[Entry, ...] = (_AXIS, _ENERGY, _HARMONIC, _LEGENDRE, _MODAL, _FROZEN_MAPPING)


@pytest.mark.parametrize("entry", ROSTER, ids=lambda e: e.id)
def test_s5_3_population(entry: Entry) -> None:
    check_population(entry)


@pytest.mark.parametrize(
    "entry,part,the_leg",
    [pytest.param(e, p, l, id=param_id(e.id, p, l[0])) for e, p, l in perturbation_ids(ROSTER)],
)
def test_s5_3_a_moved_part_moves_the_digest(entry: Entry, part: str, the_leg) -> None:
    check_perturbation(entry, part, the_leg)


@pytest.mark.parametrize(
    "entry,pair", [pytest.param(e, p, id=param_id(e.id, p[0])) for e, p in pair_ids(ROSTER)]
)
def test_s5_2_equal_content_is_one_value(entry: Entry, pair) -> None:
    """First red (the old ``repr`` scalar encoding): integer weights and their
    float twin gave different ``_structural_bytes``, so different space names."""
    check_equal_pair(entry, pair)


@pytest.mark.parametrize("entry", ROSTER, ids=lambda e: e.id)
def test_s5_2_pickle_round_trip(entry: Entry) -> None:
    check_pickle(entry)
