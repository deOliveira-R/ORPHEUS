r"""G05t: ``face_streaming_normal`` lives in ``numerics/face_measure.py`` and nowhere else (N05, spec 7.6 row G05t, ON05-9).

A layout module holds layouts only: the codim-1 face measure ``face_streaming_normal`` moves out of ``numerics/face_layout.py`` (renamed ``coproduct_layout.py`` by the next commit) before the rename.
``ast`` over ``orpheus/`` and ``tests/``: ``numerics/face_measure.py`` defines it, ``numerics/coproduct_layout.py`` and ``numerics/face_layout.py`` define no function of that name, and every import of it
(``from <module> import face_streaming_normal`` and an attribute read ``<layout module>.face_streaming_normal``) names ``orpheus.numerics.face_measure``. The count of references is read and reported
in the failure message; the spec's figure for the base tip is 7 (2 in ``spaces/angular_trace_space.py``, 5 lines of ``tests/gates/numerics/test_face_streaming_normal.py``).

First red at the base tip: the function is defined in ``numerics/face_layout.py`` and ``face_measure.py`` does not exist. Mutation: re-define the function in the layout module (the planted-source row
reads it). Exact on any host: a parse.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

pytestmark = pytest.mark.foundation

ROOT = Path(__file__).resolve().parents[4]
NAME = "face_streaming_normal"
HOME = "orpheus.numerics.face_measure"
LAYOUT_MODULES = ("face_layout", "coproduct_layout")


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def defines(source: str) -> bool:
    return any(isinstance(n, ast.FunctionDef) and n.name == NAME for n in ast.walk(ast.parse(source)))


def references(source: str, package: str = "orpheus.numerics") -> list[tuple[int, str]]:
    """``(line, module named)`` of every import of the function and every read of it as an attribute of a layout module."""
    out: list[tuple[int, str]] = []
    for n in ast.walk(ast.parse(source)):
        if isinstance(n, ast.ImportFrom) and any(a.name == NAME for a in n.names):
            module = n.module or ""
            if n.level:
                module = f"{package}.{module}" if module else package
            out.append((n.lineno, module))
        if isinstance(n, ast.Attribute) and n.attr == NAME and isinstance(n.value, ast.Name) and n.value.id in LAYOUT_MODULES:
            out.append((n.lineno, f"<attribute of {n.value.id}>"))
    return out


def test_face_streaming_normal_lives_in_face_measure() -> None:
    planted = "def face_streaming_normal(a, b):\n    return a\n"
    _require(defines(planted), "the definition predicate does not read a planted definition")
    _require(references("from orpheus.numerics.face_layout import face_streaming_normal\nx = face_layout.face_streaming_normal\n") == [(1, "orpheus.numerics.face_layout"), (2, "<attribute of face_layout>")], "the reference predicate does not read the planted imports")
    measure = ROOT / "orpheus/numerics/face_measure.py"
    _require(measure.exists(), "orpheus/numerics/face_measure.py does not exist")
    _require(defines(measure.read_text()), "orpheus/numerics/face_measure.py does not define face_streaming_normal")
    for name in LAYOUT_MODULES:
        path = ROOT / f"orpheus/numerics/{name}.py"
        if path.exists():
            _require(not defines(path.read_text()), f"orpheus/numerics/{name}.py still defines face_streaming_normal")
    found: list[tuple[str, int, str]] = []
    for top in ("orpheus", "tests"):
        for path in sorted((ROOT / top).rglob("*.py")):
            if path == Path(__file__).resolve():
                continue
            rel = path.relative_to(ROOT).as_posix()
            found += [(rel, line, module) for line, module in references(path.read_text())]
    stray = [f for f in found if f[2] != HOME]
    _require(found, "no reference to face_streaming_normal was read (the census is empty)")
    _require(not stray, f"{len(stray)} of {len(found)} references to face_streaming_normal do not name {HOME}: {stray[:6]}")
