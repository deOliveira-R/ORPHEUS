"""The concept table: one row per concept, loaded from ``docs/architecture/concept_table.toml``.

A concept is a mathematical object (``measure``, ``inverse``, ``solution``). Its row says what it is (``kind``), where
it sits on the order structure (``phase``), which package it belongs to (``layer``), the wave that removes the debt
against it (``wave``), the ruling that admitted it and the form of its head, the one site that defines it:

* ``ClassHead``: one class declares ``role: head`` for the concept, found on the tree;
* ``SumHead``: a named closed sum, a module-level union alias, whose leaves are the concept's arms;
* ``FunctionsHead``: the concept is realised by the entries of its ``functions`` field and has no class;
* ``MissingHead``: no head exists yet; the row carries the wave that creates it and may list its arms.

The status of a concept and its head class are derived, never stored. Every refusal that concerns one row is raised
here, at load; what needs two rows or the tree is a check of ``checks.py``.
"""

from __future__ import annotations

import dataclasses
import re
import tomllib
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from .declaration import KINDS
from .population import ClassKey, Population

#: The phases of the order structure, earliest first. The absence of a phase is not a member.
PHASE_ORDER = ("structure", "posing", "strategy", "solution")
WAVES = ("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8a", "W8b", "W9")
FUNCTIONS_HEAD_RULING = "O36"
_ID = re.compile(r"(EXT:)?[a-z][a-z0-9-]*\Z")
_LAYER = re.compile(r"[a-z][a-z0-9_]*\Z")
_RULING = re.compile(r"O[0-9]+\Z")
_FUNCTION = re.compile(r"orpheus(\.[A-Za-z_][A-Za-z0-9_]*)+\Z")
_THEORY = re.compile(r"(ref|eq):\S+\Z")
_KEYS = ("id", "name", "definition", "kind", "phase", "layer", "wave", "ruling", "sum", "functions_head", "missing", "arms",
         "functions", "theory")


class TableError(ValueError):
    """A row of the concept table that is not well formed, or a table that is not TOML."""


def rank(phase: str) -> int:
    """The position of a phase on ``PHASE_ORDER``."""
    return PHASE_ORDER.index(phase)


@dataclasses.dataclass(frozen=True)
class ClassHead:
    """The head is the one class that declares ``role: head`` for the concept."""


@dataclasses.dataclass(frozen=True)
class SumHead:
    """The head is a named closed sum: the module-level union alias ``alias`` of ``file``."""

    file: str
    alias: str

    def __post_init__(self) -> None:
        if not self.file.endswith(".py") or not self.alias.isidentifier():
            raise TableError(f"sum head '{self.file}:{self.alias}' is not '<file>.py:<Alias>'")


@dataclasses.dataclass(frozen=True)
class FunctionsHead:
    """The concept is realised by the entries of its ``functions`` field; it has no class, no alias and no wave."""


@dataclasses.dataclass(frozen=True)
class MissingHead:
    """The head does not exist yet; ``arms`` names the classes that already realise the concept."""

    arms: tuple[ClassKey, ...] = ()


Head = ClassHead | SumHead | FunctionsHead | MissingHead


@dataclasses.dataclass(frozen=True)
class Concept:
    """One row of the table. Invariants are refused at construction with ``TableError``."""

    id: str
    name: str
    definition: str
    kind: str
    phase: str | None
    layer: str
    wave: str | None
    ruling: str
    head: Head
    functions: tuple[str, ...] = ()
    theory: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        def refuse(field: str, value: object, allowed: str) -> None:
            raise TableError(f"concept {self.id!r}: {field} {value!r} is not {allowed}")

        if not _ID.match(self.id):
            refuse("id", self.id, "an id ('[a-z][a-z0-9-]*', optionally prefixed 'EXT:')")
        if not self.name.strip() or not self.definition.strip():
            refuse("name or definition", "", "empty")
        if self.kind not in KINDS:
            refuse("kind", self.kind, f"one of {KINDS}")
        if (self.phase is None) != (self.kind == "apparatus"):
            refuse("phase", self.phase, "absent exactly for kind 'apparatus'")
        if self.phase is not None and self.phase not in PHASE_ORDER:
            refuse("phase", self.phase, f"one of {PHASE_ORDER}")
        if not _LAYER.match(self.layer):
            refuse("layer", self.layer, "a package identifier")
        if self.wave is not None and self.wave not in WAVES:
            refuse("wave", self.wave, f"None or one of {WAVES}")
        if not _RULING.match(self.ruling):
            refuse("ruling", self.ruling, "a ruling 'O<number>'")
        if isinstance(self.head, MissingHead) and self.wave is None:
            refuse("wave", None, "absent on a row whose head is missing")
        if isinstance(self.head, FunctionsHead):
            if not self.functions:
                refuse("functions", (), "empty on a functions head")
            if self.wave is not None:
                refuse("wave", self.wave, "present on a functions head")
            if self.ruling != FUNCTIONS_HEAD_RULING:
                refuse("ruling", self.ruling, f"{FUNCTIONS_HEAD_RULING}, the ruling that admits a functions head")
        for entry in self.functions:
            if not _FUNCTION.match(entry):
                refuse("functions entry", entry, "a dotted name 'orpheus.<module>.<name>'")
        for label in self.theory:
            if not _THEORY.match(label):
                refuse("theory label", label, "'ref:<label>' or 'eq:<label>'")

    @property
    def slug(self) -> str:
        """The anchor of the concept on the generated page: ``concept-`` plus the id lowercased, ``:`` as ``-``."""
        return "concept-" + self.id.lower().replace(":", "-")


def _head_of(row: Mapping[str, Any], cid: str) -> Head:
    forms = [name for name in ("sum", "functions_head", "missing") if name in row]
    if len(forms) > 1:
        raise TableError(f"concept {cid!r}: two head forms {forms}")
    if "arms" in row and "missing" not in row:
        raise TableError(f"concept {cid!r}: 'arms' without 'missing'")
    if "sum" in row:
        file, sep, alias = str(row["sum"]).partition(":")
        if not sep:
            raise TableError(f"concept {cid!r}: sum {row['sum']!r} is not '<file>:<Alias>'")
        return SumHead(file, alias)
    if "functions_head" in row:
        if row["functions_head"] is not True:
            raise TableError(f"concept {cid!r}: functions_head is true or absent")
        return FunctionsHead()
    if "missing" in row:
        if row["missing"] is not True:
            raise TableError(f"concept {cid!r}: missing is true or absent")
        arms = row.get("arms", [])
        if not isinstance(arms, list):
            raise TableError(f"concept {cid!r}: arms is a list of '<file>::<qualname>'")
        return MissingHead(tuple(_parse_key(a, cid) for a in arms))
    return ClassHead()


def _parse_key(text: object, cid: str) -> ClassKey:
    try:
        return ClassKey.parse(str(text))
    except ValueError as error:
        raise TableError(f"concept {cid!r}: {error}") from error


def _text(row: Mapping[str, Any], field: str, cid: str) -> str:
    value = row.get(field)
    if not isinstance(value, str):
        raise TableError(f"concept {cid!r}: {field} {value!r} is not a string")
    return value


def _strings(row: Mapping[str, Any], field: str, cid: str) -> tuple[str, ...]:
    value = row.get(field, [])
    if not isinstance(value, list) or not all(isinstance(x, str) for x in value):
        raise TableError(f"concept {cid!r}: {field} {value!r} is not a list of strings")
    return tuple(value)


def loads_table(text: str) -> tuple[Concept, ...]:
    """The concepts of a table spelled as TOML text, in file order."""
    try:
        data = tomllib.loads(text)
    except tomllib.TOMLDecodeError as error:
        raise TableError(f"the concept table is not TOML: {error}") from error
    if set(data) - {"concept"}:
        raise TableError(f"the concept table has keys {sorted(set(data) - {'concept'})} besides [[concept]]")
    rows = data.get("concept", [])
    if not isinstance(rows, list):
        raise TableError("'concept' is an array of tables")
    out = []
    for row in rows:
        cid = row.get("id", "<no id>")
        unknown = sorted(set(row) - set(_KEYS))
        if unknown:
            raise TableError(f"concept {cid!r}: unknown keys {unknown}")
        for required in ("id", "name", "definition", "kind", "layer", "ruling"):
            if required not in row:
                raise TableError(f"concept {cid!r}: required key '{required}' is missing")
        phase, wave = row.get("phase"), row.get("wave")
        for field, value in (("phase", phase), ("wave", wave)):
            if value is not None and not isinstance(value, str):
                raise TableError(f"concept {cid!r}: {field} {value!r} is not a string")
        out.append(Concept(
            id=_text(row, "id", cid), name=_text(row, "name", cid), definition=_text(row, "definition", cid),
            kind=_text(row, "kind", cid), phase=phase, layer=_text(row, "layer", cid), wave=wave,
            ruling=_text(row, "ruling", cid), head=_head_of(row, cid), functions=_strings(row, "functions", cid),
            theory=_strings(row, "theory", cid)))
    return tuple(out)


def load_table(path: Path) -> tuple[Concept, ...]:
    """The concepts of the table file at ``path``."""
    return loads_table(Path(path).read_text())


def heads(population: Population, concepts: tuple[Concept, ...]) -> dict[str, tuple[ClassKey, ...]]:
    """For every concept id, the classes that declare ``role: head`` for it, in key order.

    The one derivation of a head; whether the count is right for the concept's head form is ``check_heads``'s.
    """
    out: dict[str, list[ClassKey]] = {c.id: [] for c in concepts}
    for row in population.rows:
        if row.declaration is not None and row.declaration.role == "head" and row.declaration.concept in out:
            out[row.declaration.concept].append(row.key)
    return {cid: tuple(keys) for cid, keys in out.items()}
