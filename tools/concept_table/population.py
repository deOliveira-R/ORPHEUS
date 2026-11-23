"""The population of the concept table: every class of ``orpheus/`` outside ``derivations/``, keyed and sorted.

A class is keyed by ``ClassKey(file, qualname)``. Its declaration is the trailing comment on the line that ends its
header (``declaration.py``). The population reads source text only: nothing is imported and nothing runs.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Mapping
from pathlib import Path
from typing import NamedTuple

from .declaration import Declaration, DeclarationError, parse_declaration, read_declarations


class DuplicateClass(ValueError):
    """Two classes with one key, or with one dotted name (the spelling the reference pass resolves by)."""


class ClassKey(NamedTuple):
    """A class of the tree: the repository-relative file and the qualified name inside it.

    ``qualname`` is ``Name`` for a module-level class, ``Outer.Inner`` for a class nested in a class and
    ``outer.<locals>.Name`` for a class local to a function. ``str(key)`` is ``file::qualname``.
    """

    file: str
    qualname: str

    def __str__(self) -> str:
        return f"{self.file}::{self.qualname}"

    @property
    def dotted(self) -> str:
        """The one definition of a class's dotted name.

        The file path without ``.py`` with ``/`` replaced by ``.``, a package ``__init__.py`` collapsing to the
        package path, then ``.`` and the qualname; for a function-local class the ``outer.<locals>.`` prefix is
        dropped, so that ``a/b.py::f.<locals>.C`` is ``a.b.C``.
        """
        parts = self.file.removesuffix(".py").split("/")
        if parts[-1] == "__init__":
            parts = parts[:-1]
        return ".".join([*parts, self.qualname.rsplit("<locals>.", 1)[-1]])

    @classmethod
    def parse(cls, text: str) -> ClassKey:
        """The key spelled ``file::qualname`` (the spelling of the debt file); refuses any other shape."""
        file, sep, qualname = text.partition("::")
        if not sep or not file.endswith(".py") or file.startswith("/") or not qualname:
            raise ValueError(f"class key {text!r} is not '<file>.py::<qualname>'")
        return cls(file, qualname)


@dataclasses.dataclass(frozen=True)
class ClassRow:
    """One class: its key, the line of its ``class`` keyword, the line that ends its header, its declaration.

    ``markers`` holds every ``# concept:`` comment on the header-end line; the declaration is the parse of the one
    marker when there is exactly one and it parses, otherwise ``None``.
    """

    key: ClassKey
    lineno: int
    header_end: int
    declaration: Declaration | None
    markers: tuple[str, ...]

    def __post_init__(self) -> None:
        if self.lineno < 1 or self.header_end < self.lineno:
            raise ValueError(f"{self.key}: header lines {self.lineno}..{self.header_end} are not a span of lines")
        if self.declaration is not None and len(self.markers) != 1:
            raise ValueError(f"{self.key}: a declaration needs exactly one marker, has {len(self.markers)}")


@dataclasses.dataclass(frozen=True)
class Population:
    """The classes of a tree, sorted by ``ClassKey``, with the markers that sit on no class header-end line.

    Invariants, refused at construction: keys are unique, dotted names are unique, ``rows`` is sorted by key.
    """

    rows: tuple[ClassRow, ...]
    orphan_markers: tuple[tuple[str, int], ...]
    files_read: int
    files_skipped_derivations: int

    def __post_init__(self) -> None:
        keys = [row.key for row in self.rows]
        if len(set(keys)) != len(keys):
            raise DuplicateClass(f"two classes with one key: {sorted({str(k) for k in keys if keys.count(k) > 1})}")
        dotted = [k.dotted for k in keys]
        clash = sorted({d for d in dotted if dotted.count(d) > 1})
        if clash:
            raise DuplicateClass(f"two classes with one dotted name: {clash}")
        if keys != sorted(keys):
            raise ValueError("population rows are not sorted by ClassKey")
        if self.files_read < 0 or self.files_skipped_derivations < 0:
            raise ValueError("file counts are not negative")

    @classmethod
    def from_sources(cls, sources: Mapping[str, str]) -> Population:
        """The population of a mapping ``{file: source text}``; every file given is read, ``derivations/`` included."""
        rows: list[ClassRow] = []
        orphans: list[tuple[str, int]] = []
        for file in sorted(sources):
            if not file.endswith(".py") or file.startswith("/"):
                raise ValueError(f"source key {file!r} is not a repository-relative '.py' path")
            try:
                found, orphan_lines = read_declarations(sources[file])
            except SyntaxError as error:
                raise ValueError(f"{file} is not Python source: {error}") from error
            orphans.extend((file, line) for line in orphan_lines)
            for qualname, (lineno, end, markers) in found.items():
                declaration = None
                if len(markers) == 1:
                    try:
                        declaration = parse_declaration(markers[0])
                    except DeclarationError:
                        declaration = None
                rows.append(ClassRow(ClassKey(file, qualname), lineno, end, declaration, tuple(markers)))
        return cls(tuple(sorted(rows, key=lambda r: r.key)), tuple(orphans), len(sources), 0)


def read_tree(root: Path) -> Population:
    """The population of ``root/orpheus/**/*.py``: sorted, without ``__pycache__`` and without ``derivations/``."""
    root = Path(root)
    if not (root / "orpheus").is_dir():
        raise ValueError(f"{root} holds no orpheus/ directory")
    sources: dict[str, str] = {}
    skipped = 0
    for path in sorted((root / "orpheus").rglob("*.py")):
        if "__pycache__" in path.parts:
            continue
        if "derivations" in path.relative_to(root).parts:
            skipped += 1
            continue
        sources[path.relative_to(root).as_posix()] = path.read_text()
    return dataclasses.replace(Population.from_sources(sources), files_skipped_derivations=skipped)
