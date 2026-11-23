"""The pinned debt and the fingerprint seed: lists of exceptions compared with what the gates measure, by equality.

Each debt list is a tuple of one entry type, and a violation that a check returns is a value of the type of the entry
that would pin it (an ``EdgeViolation`` is an ``EdgeDebt``), so ``pin`` compares two sets with nothing to parse. An
entry carries no wave: the wave of an entry is the wave of its concept, derived by ``wave_of``, the one function that
single-sources it. That a list never grows without a ruling is a review rule: every change of a list is a visible
edit of a tracked file.
"""

from __future__ import annotations

import dataclasses
import re
import tomllib
from collections.abc import Hashable, Mapping
from pathlib import Path
from typing import Any, TypeVar

from .population import ClassKey, Population
from .refs import EDGE_KINDS
from .table import Concept

ROLE_TESTS = ("arm", "twin_is_arm", "interface", "part_component_candidate", "part_no_core_reader", "component_no_core_reader")
#: The laws a kind-law entry may pin; the guards of the other laws are never pinned.
PINNABLE_LAWS = ("head_arm_kind", "apparatus_row")
#: Every law G4 states.
KIND_LAWS = ("head_arm_kind", "annex_record", "operator_subkind", "element_subkind", "covector_subkind", "apparatus_row")
_FINGERPRINT = re.compile(r"blake2b:[0-9a-f]{64}\Z")

T = TypeVar("T", bound=Hashable)


class DebtError(ValueError):
    """A debt or seed file that is not well formed."""


@dataclasses.dataclass(frozen=True, order=True)
class ClassDebt:
    """A class on the twin, weld or misplaced list."""

    key: ClassKey

    def __str__(self) -> str:
        return str(self.key)


@dataclasses.dataclass(frozen=True, order=True)
class RoleTestDebt:
    """A class that fails the structural test ``test`` of its role."""

    test: str
    key: ClassKey

    def __post_init__(self) -> None:
        if self.test not in ROLE_TESTS:
            raise DebtError(f"role test {self.test!r} is not one of {ROLE_TESTS}")

    def __str__(self) -> str:
        return f"{self.key} [{self.test}]"


@dataclasses.dataclass(frozen=True, order=True)
class EdgeDebt:
    """A reference between two classes that the phase order (or the layer order) refuses."""

    src: ClassKey
    dst: ClassKey
    kind: str

    def __post_init__(self) -> None:
        if self.kind not in EDGE_KINDS:
            raise DebtError(f"edge kind {self.kind!r} is not one of {EDGE_KINDS}")

    def __str__(self) -> str:
        return f"{self.src} -> {self.dst} [{self.kind}]"


@dataclasses.dataclass(frozen=True)
class CycleDebt:
    """A set of concepts that refer to one another in a cycle."""

    concepts: frozenset[str]

    def __post_init__(self) -> None:
        if len(self.concepts) < 2:
            raise DebtError(f"a cycle holds at least two concepts, got {sorted(self.concepts)}")

    def __str__(self) -> str:
        return f"concepts [{', '.join(sorted(self.concepts))}]"

    def __lt__(self, other: CycleDebt) -> bool:
        """The larger cycle first, then by the sorted concept ids (a set has no order of its own)."""
        return (-len(self.concepts), sorted(self.concepts)) < (-len(other.concepts), sorted(other.concepts))


@dataclasses.dataclass(frozen=True, order=True)
class KindLawDebt:
    """A class that violates the kind law ``law``."""

    law: str
    key: ClassKey

    def __post_init__(self) -> None:
        if self.law not in KIND_LAWS:
            raise DebtError(f"kind law {self.law!r} is not one of {KIND_LAWS}")

    def __str__(self) -> str:
        return f"{self.key} [{self.law}]"


@dataclasses.dataclass(frozen=True, order=True)
class OffAlgebraDebt:
    """A class with sub-kind ``operator`` that is not on the operator algebra."""

    key: ClassKey

    def __str__(self) -> str:
        return f"{self.key} [operator_off_algebra]"


@dataclasses.dataclass(frozen=True, order=True)
class CovectorOffLineageDebt:
    """A class with sub-kind ``covector`` that is not on the lineage of ``Covector``."""

    key: ClassKey

    def __str__(self) -> str:
        return f"{self.key} [covector_off_lineage]"


@dataclasses.dataclass(frozen=True)
class SeedEntry:
    """The fingerprint a class had when it was put on the twin or weld list."""

    key: ClassKey
    fingerprint: str

    def __post_init__(self) -> None:
        if not _FINGERPRINT.match(self.fingerprint):
            raise DebtError(f"fingerprint {self.fingerprint!r} is not 'blake2b:' and 64 hex digits")


Entry = ClassDebt | RoleTestDebt | EdgeDebt | CycleDebt | KindLawDebt | OffAlgebraDebt | CovectorOffLineageDebt


@dataclasses.dataclass(frozen=True)
class Debt:
    """The nine pinned lists."""

    twin: tuple[ClassDebt, ...] = ()
    weld: tuple[ClassDebt, ...] = ()
    role_test: tuple[RoleTestDebt, ...] = ()
    misplaced: tuple[ClassDebt, ...] = ()
    edge: tuple[EdgeDebt, ...] = ()
    cycle: tuple[CycleDebt, ...] = ()
    kind_law: tuple[KindLawDebt, ...] = ()
    operator_off_algebra: tuple[OffAlgebraDebt, ...] = ()
    covector_off_lineage: tuple[CovectorOffLineageDebt, ...] = ()

    def lists(self) -> dict[str, tuple[Entry, ...]]:
        """The lists by name, in file order."""
        return {f.name: getattr(self, f.name) for f in dataclasses.fields(self)}


def pin(measured: frozenset[T], pinned: frozenset[T]) -> tuple[frozenset[T], frozenset[T]]:
    """``(unpinned, stale)``: what is measured and not pinned, and what is pinned and no longer measured."""
    return measured - pinned, pinned - measured


def wave_of(entry: Entry, concepts: tuple[Concept, ...], population: Population) -> frozenset[str]:
    """The waves that own a debt entry: a set, empty when the entry's concept has no wave.

    A class entry has the wave of its declaration's concept; an edge, the wave of its source class's concept; a cycle,
    the waves of those of its concepts that have one.
    """
    wave = {c.id: c.wave for c in concepts}
    concept_of = {r.key: r.declaration.concept for r in population.rows if r.declaration is not None}

    def of_class(key: ClassKey) -> frozenset[str]:
        found = wave.get(concept_of.get(key, ""))
        return frozenset() if found is None else frozenset({found})

    if isinstance(entry, CycleDebt):
        return frozenset(w for c in entry.concepts if (w := wave.get(c)) is not None)
    if isinstance(entry, EdgeDebt):
        return of_class(entry.src)
    return of_class(entry.key)


# --- the loaders -----------------------------------------------------------------------------------------------------

def _entries(data: Mapping[str, Any], name: str, fields: tuple[str, ...]) -> list[dict[str, Any]]:
    rows = data.get(name, [])
    if not isinstance(rows, list):
        raise DebtError(f"'{name}' is an array of tables")
    for row in rows:
        if set(row) != set(fields):
            raise DebtError(f"[[{name}]] entry {dict(row)} has keys {sorted(row)}, the schema is {sorted(fields)}")
        for field in fields:
            if field != "concepts" and not isinstance(row[field], str):
                raise DebtError(f"[[{name}]] {field} {row[field]!r} is not a string")
    return rows


def _key(text: str) -> ClassKey:
    try:
        return ClassKey.parse(text)
    except ValueError as error:
        raise DebtError(str(error)) from error


def loads_debt(text: str) -> Debt:
    """The debt of a file spelled as TOML text; refuses any key or list outside the schema and any entry's extra field."""
    try:
        data = tomllib.loads(text)
    except tomllib.TOMLDecodeError as error:
        raise DebtError(f"the debt file is not TOML: {error}") from error
    names = tuple(f.name for f in dataclasses.fields(Debt))
    if set(data) - set(names):
        raise DebtError(f"the debt file has lists {sorted(set(data) - set(names))} outside {names}")
    cls = lambda name: tuple(ClassDebt(_key(r["class"])) for r in _entries(data, name, ("class",)))  # noqa: E731
    cycles = []
    for row in _entries(data, "cycle", ("concepts",)):
        concepts = row["concepts"]
        if not isinstance(concepts, list) or not all(isinstance(c, str) for c in concepts) or len(set(concepts)) != len(concepts):
            raise DebtError(f"[[cycle]] concepts {concepts!r} is not a list of distinct ids")
        cycles.append(CycleDebt(frozenset(concepts)))
    laws = []
    for row in _entries(data, "kind_law", ("law", "class")):
        if row["law"] not in PINNABLE_LAWS:
            raise DebtError(f"kind law {row['law']!r} is not pinnable; the pinnable laws are {PINNABLE_LAWS}")
        laws.append(KindLawDebt(row["law"], _key(row["class"])))
    return Debt(
        twin=cls("twin"), weld=cls("weld"),
        role_test=tuple(RoleTestDebt(r["test"], _key(r["class"])) for r in _entries(data, "role_test", ("test", "class"))),
        misplaced=cls("misplaced"),
        edge=tuple(EdgeDebt(_key(r["src"]), _key(r["dst"]), r["kind"]) for r in _entries(data, "edge", ("src", "dst", "kind"))),
        cycle=tuple(cycles), kind_law=tuple(laws),
        operator_off_algebra=tuple(OffAlgebraDebt(_key(r["class"])) for r in _entries(data, "operator_off_algebra", ("class",))),
        covector_off_lineage=tuple(CovectorOffLineageDebt(_key(r["class"])) for r in _entries(data, "covector_off_lineage", ("class",))))


def load_debt(path: Path) -> Debt:
    """The debt of the file at ``path``."""
    return loads_debt(Path(path).read_text())


def loads_seed(text: str) -> tuple[SeedEntry, ...]:
    """The fingerprint seed spelled as TOML text; one entry per class, and no list but ``[[seed]]``."""
    try:
        data = tomllib.loads(text)
    except tomllib.TOMLDecodeError as error:
        raise DebtError(f"the seed file is not TOML: {error}") from error
    if set(data) - {"seed"}:
        raise DebtError(f"the seed file has lists {sorted(set(data) - {'seed'})} besides [[seed]]")
    entries = tuple(SeedEntry(_key(r["class"]), r["fingerprint"]) for r in _entries(data, "seed", ("class", "fingerprint")))
    keys = [e.key for e in entries]
    if len(set(keys)) != len(keys):
        raise DebtError(f"the seed file has two entries for {sorted({str(k) for k in keys if keys.count(k) > 1})}")
    return entries


def load_seed(path: Path) -> tuple[SeedEntry, ...]:
    """The fingerprint seed of the file at ``path``."""
    return loads_seed(Path(path).read_text())
