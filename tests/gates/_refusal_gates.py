r"""The machinery the W0 refusal gates share (unit ``defects``, #709 and #711; the gate of #710 retired with N03).

A refusal gate has two instruments that an in-process ``pytest.raises`` cannot be:

* a **cold child** (``python -O -P -c``), which asserts that the refusal is a ``raise``
  and survives the ``-O`` the suite runs under (a ``monkeypatch`` in the pytest process
  never reaches the child);
* a **refusal-removal mutant**: a copy of the ``*.py`` files of ``orpheus/`` in which the
  first ``raise`` after the refusal's tag line is replaced by ``assert False, (``, so that
  under ``-O`` the refused input constructs again.

The bullets "Conventions" and "A leg that imports a patched copy of ``orpheus``" of the
specification (``.claude/plans/ontology_w0_spec.md`` 1.3) are what this file implements:
the child runs ``python -O -P -c`` with ``PYTHONPATH=<root>:<repository root>`` and prints
``orpheus.__file__`` first, then the sentinel ``IMPORTED`` after its last import; a mutant
whose child did not print the sentinel is never red evidence.

``import_root`` is the directory the PYTEST process imports ``orpheus`` from, not the
directory of this file: a run on a patched copy (``-P -o pythonpath=``) makes the
unmutated leg and every copy mutation read that copy.

A precondition that is missing (a tag line that does not exist yet, an attribute a
mutation needs) raises :class:`Uninstallable`, a distinct failure, never an
``AttributeError`` that reads as a refusal gate's own red.
"""

from __future__ import annotations

import os
import pathlib
import subprocess
import sys
from typing import NamedTuple

import orpheus

#: The repository root whose ``tests`` package a child imports (this file is ``tests/gates/_refusal_gates.py``).
REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]

SENTINEL = "IMPORTED"


class Uninstallable(Exception):
    """A mutation's precondition does not hold in this tree (no tag line, an absent attribute): not a gate's red."""


class Child(NamedTuple):
    returncode: int
    stdout: str
    stderr: str


def import_root() -> pathlib.Path:
    """The directory that holds the ``orpheus`` package this pytest process imported."""
    return pathlib.Path(orpheus.__file__).resolve().parents[1]


def copy_orpheus_py(source_root: pathlib.Path, destination: pathlib.Path) -> pathlib.Path:
    """Copy every ``*.py`` of ``source_root/orpheus`` (no ``__pycache__``) to ``destination/orpheus``; returns ``destination``.

    The tree is 846 MB on disk, almost all of it data; its ``*.py`` files are 9 to 12 MB.
    """
    package = source_root / "orpheus"
    for path in package.rglob("*.py"):
        if "__pycache__" in path.parts:
            continue
        target = destination / path.relative_to(source_root)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(path.read_bytes())
    return destination


def remove_refusal(root: pathlib.Path, relative_file: str, issue: int) -> None:
    """In the copy at ``root``, replace the first ``raise ValueError(`` / ``raise TypeError(`` after the tag line of ``issue`` by ``assert False, (``.

    The anchor is the TAG line ``ELEGANCE-DEBT[guard] #<issue>``, which must occur exactly once in the
    file (``raise ValueError(`` occurs 6 to 12 times in each of the three files). The call's arguments and
    its closing parenthesis stay, so the copy parses; under ``-O`` the assert is compiled out and the
    refused input constructs again.
    """
    path = root / relative_file
    text = path.read_text(encoding="utf-8")
    tag = f"ELEGANCE-DEBT[guard] #{issue}"
    if text.count(tag) != 1:
        raise Uninstallable(f"{relative_file}: the tag line {tag!r} occurs {text.count(tag)} times, not once")
    start = text.index(tag)
    candidates = [i for i in (text.find(k, start) for k in ("raise ValueError(", "raise TypeError(")) if i != -1]
    if not candidates:
        raise Uninstallable(f"{relative_file}: no raise after the tag line {tag!r}")
    at = min(candidates)
    end = at + len("raise ValueError(" if text.startswith("raise ValueError(", at) else "raise TypeError(")
    path.write_text(text[:at] + "assert False, (" + text[end:], encoding="utf-8")


def replace_once(root: pathlib.Path, relative_file: str, old: str, new: str) -> None:
    """In the copy at ``root``, replace the text ``old`` (which must occur exactly once) by ``new``."""
    path = root / relative_file
    text = path.read_text(encoding="utf-8")
    if text.count(old) != 1:
        raise Uninstallable(f"{relative_file}: the search text {old!r} occurs {text.count(old)} times, not once")
    path.write_text(text.replace(old, new), encoding="utf-8")


def run_child(imports: str, body: str, *, root: pathlib.Path | None = None) -> Child:
    """Run ``imports`` then ``body`` in ``python -O -P -c`` with ``orpheus`` read from ``root`` (default: where pytest reads it).

    The child prints ``orpheus.__file__``, runs ``imports``, prints :data:`SENTINEL`, runs ``body``.
    """
    root = import_root() if root is None else root
    code = f"import orpheus\nprint(orpheus.__file__)\n{imports}\nprint({SENTINEL!r})\n{body}\n"
    path = os.pathsep.join(dict.fromkeys([str(root), str(REPO_ROOT)]))
    done = subprocess.run(
        [sys.executable, "-O", "-P", "-c", code],
        env={**os.environ, "PYTHONPATH": path},
        capture_output=True, text=True, cwd=REPO_ROOT, check=False,
    )
    return Child(done.returncode, done.stdout, done.stderr)


def imported_from(child: Child, root: pathlib.Path) -> bool:
    """The child imported ``orpheus`` from under ``root`` and reached its sentinel (the mutant imported)."""
    lines = child.stdout.splitlines()
    if len(lines) < 2 or lines[1] != SENTINEL:
        return False
    return pathlib.Path(lines[0]).resolve().is_relative_to(root.resolve())


def refused(child: Child, fragment: str) -> bool:
    """The body was refused: non-zero exit and the message fragment on stderr."""
    return child.returncode != 0 and fragment in child.stderr


def constructed(child: Child, fragment: str) -> bool:
    """The body ran to its end: exit 0 and the refusal's fragment absent from stderr."""
    return child.returncode == 0 and fragment not in child.stderr
