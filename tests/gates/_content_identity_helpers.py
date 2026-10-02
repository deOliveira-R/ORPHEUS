r"""Shared machinery of the content-identity gates (#405 P1 step 5, spec §1.5).

Not a test module (no ``test_`` prefix). Each owning tree (``tests/gates/data``,
``tests/gates/geometry``, ``tests/gates/mesh``) declares a ROSTER of
:class:`Entry` values, one per type whose identity is its content, and runs
the per-type legs of S5.2 (equal content: equal ``==``, ``hash`` and digest)
and S5.3 (a change to any part changes the digest) through the checkers here.
``tests/gates/numerics/test_content_identity.py`` imports the three rosters for
the population gate S5.9 and the route gate S5.6.

The population is the TYPE, never a hand list (spec S5.3, ``instrument-doctrine``
X2 STRUCTURAL-DENOMINATOR): an entry declares its expected part names, the gate
reads the part names off the production instance (``content_parts()`` for a
:class:`~orpheus.numerics.content.ContentIdentity`, the ``compare=True``
dataclass fields otherwise) and requires the perturbation table to cover
exactly them. A part added later reds S5.3 until it is perturbed.
"""

from __future__ import annotations

import dataclasses
import pickle
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from typing import Any

from orpheus.numerics.content import ContentIdentity, content_digest, encode

#: One perturbation: (leg id, builder of a value differing from the base in
#: this part, the OTHER parts the leg is allowed to move because a
#: construction law couples them, e.g. a mesh's coordinate system and its
#: volumes).
Leg = tuple[str, Callable[[], Any], tuple[str, ...]]

#: One equal-content pair: (leg id, builder of a, builder of b).
Pair = tuple[str, Callable[[], Any], Callable[[], Any]]


def require(cond: object, msg: str) -> None:
    """``-O``-safe assertion (the canonical runner strips a bare ``assert``)."""
    if not cond:
        raise AssertionError(msg)


def leg(leg_id: str, build: Callable[[], Any], *also: str) -> Leg:
    return (leg_id, build, tuple(also))


@dataclass(frozen=True)
class Entry:
    """One type whose identity is its content.

    ``fields_are_parts``: the part names must equal the ``compare=True``
    dataclass fields (S5.3's quantification over ``dataclasses.fields``);
    False for a type whose content is not its field list: a
    ``content_parts()`` override (``FrozenMapping`` and so ``FaceLaws``) or
    derived fields declared ``compare=False`` (``Mesh1D``'s ``widths``,
    ``centers``, ``areas``).
    ``content_identity``: the type carries the mixin (it has ``==`` and
    ``hash`` from the digest); False for a frozen-dataclass PART of a law
    (a reemission closure, an inflow source), which is content through the
    encoder's frozen-dataclass rule and keeps its dataclass equality.
    """

    cls: type
    base: Callable[[], Any]
    parts: tuple[str, ...]
    perturb: Mapping[str, tuple[Leg, ...]]
    pairs: tuple[Pair, ...] = ()
    fields_are_parts: bool = True
    content_identity: bool = True
    pickles: bool = True
    extra: Mapping[str, Any] = field(default_factory=dict)

    @property
    def id(self) -> str:
        return self.cls.__qualname__


def parts_of(value: Any) -> tuple[tuple[str, Any], ...]:
    """The named content parts the encoder reads for ``value``."""
    if isinstance(value, ContentIdentity):
        return tuple(value.content_parts())
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return tuple(
            (f.name, getattr(value, f.name)) for f in dataclasses.fields(value) if f.compare
        )
    raise AssertionError(
        f"{type(value).__qualname__} has no content parts: it is neither a "
        f"ContentIdentity nor a dataclass"
    )


def check_population(entry: Entry) -> None:
    """S5.3's denominator: the parts the production instance has, against the
    entry's expectation and its perturbation table."""
    base = entry.base()
    require(type(base) is entry.cls, f"{entry.id}: the base builder returned {type(base)}")
    require(
        isinstance(base, ContentIdentity) == entry.content_identity,
        f"{entry.id}: carries the ContentIdentity mixin = {isinstance(base, ContentIdentity)}, "
        f"the roster says {entry.content_identity}",
    )
    names = tuple(name for name, _ in parts_of(base))
    require(
        names == entry.parts,
        f"{entry.id}: the content parts are {names}, the roster expects {entry.parts} "
        f"(a part was added, removed, renamed or reordered: perturb it here)",
    )
    if entry.fields_are_parts:
        compare_fields = tuple(f.name for f in dataclasses.fields(entry.cls) if f.compare)
        require(
            names == compare_fields,
            f"{entry.id}: the content parts {names} are not the compare=True dataclass "
            f"fields {compare_fields}",
        )
    require(
        set(entry.perturb) == set(names),
        f"{entry.id}: the perturbation table covers {sorted(entry.perturb)}, the parts are "
        f"{sorted(names)}",
    )
    for name in names:
        require(len(entry.perturb[name]) >= 1, f"{entry.id}.{name}: no perturbation leg")


def check_perturbation(entry: Entry, part: str, the_leg: Leg) -> None:
    """S5.3, one leg: the moved value differs from the base in its digest, its
    ``==`` and a container, and only in ``part`` (plus declared co-movers), so
    the red names the part an encoder dropped."""
    leg_id, build, also = the_leg
    base, moved = entry.base(), build()
    require(type(moved) is type(base), f"{entry.id}.{part}[{leg_id}]: type changed")
    base_parts, moved_parts = dict(parts_of(base)), dict(parts_of(moved))
    require(
        encode(base_parts[part]) != encode(moved_parts[part]),
        f"{entry.id}.{part}[{leg_id}]: activation: the leg does not move the part's content",
    )
    for other in base_parts:
        if other == part or other in also:
            continue
        require(
            encode(base_parts[other]) == encode(moved_parts[other]),
            f"{entry.id}.{part}[{leg_id}]: the leg also moved the part {other!r} "
            f"(undeclared), so the row cannot attribute a red",
        )
    require(
        content_digest(moved) != content_digest(base),
        f"{entry.id}.{part}[{leg_id}]: the digest did not move: the encoder drops {part!r}",
    )
    require(not (moved == base), f"{entry.id}.{part}[{leg_id}]: == still holds")
    require(moved != base, f"{entry.id}.{part}[{leg_id}]: != disagrees with ==")
    if entry.content_identity:
        # Separation through the CONTAINER, never ``hash(a) != hash(b)``
        # (test-architect lessons §1, the identity family (b)).
        require(len({moved, base}) == 2, f"{entry.id}.{part}[{leg_id}]: a set merged them")


def check_equal_pair(entry: Entry, pair: Pair) -> None:
    """S5.2, one pair: two independently built values of equal content are one
    value to ``==``, ``hash``, a set and the digest."""
    leg_id, build_a, build_b = pair
    a, b = build_a(), build_b()
    tag = f"{entry.id}[{leg_id}]"
    require(type(a) is type(b) is entry.cls, f"{tag}: the pair is not two {entry.id}")
    require(a is not b, f"{tag}: activation: the builders returned one object")
    require(content_digest(a) == content_digest(b), f"{tag}: the digests differ")
    require(a == b and b == a, f"{tag}: equal content must compare equal, both ways")
    require(not (a != b), f"{tag}: != disagrees with ==")
    require(hash(a) == hash(b), f"{tag}: equal values must hash equal")
    require(len({a, b}) == 1, f"{tag}: a set of two equal values holds one member")


def check_pickle(entry: Entry) -> None:
    """The persisted round trip (ruling of 2026-10-01: the cache pickles its
    keys' values): the reloaded value is equal, hashes equal and digests equal."""
    base = entry.base()
    twin = pickle.loads(pickle.dumps(base))
    tag = f"{entry.id}[pickle]"
    require(type(twin) is type(base), f"{tag}: the type changed")
    require(twin == base, f"{tag}: the reloaded value is not equal")
    require(hash(twin) == hash(base), f"{tag}: the reloaded value hashes differently")
    require(content_digest(twin) == content_digest(base), f"{tag}: the digest moved")


def perturbation_ids(roster: tuple[Entry, ...]) -> list[tuple[Entry, str, Leg]]:
    return [
        (entry, part, the_leg)
        for entry in roster
        for part, legs in entry.perturb.items()
        for the_leg in legs
    ]


def pair_ids(roster: tuple[Entry, ...]) -> list[tuple[Entry, Pair]]:
    return [(entry, pair) for entry in roster for pair in entry.pairs]


def param_id(*words: str) -> str:
    """A pytest id with no whitespace, so ``-rf`` reports it whole."""
    return "-".join(w.replace(" ", "_") for w in words)
