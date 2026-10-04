r"""Synthetic generator packages for the traced-memo gates (#405 P3): every mutation witness edits a COPY.

A gate never edits a tracked file (``process-discipline``, "Mutation-testing an uncommitted file"): each
test writes a fresh package under its ``tmp_path``, with a name unique to the test (so ``sys.modules`` never
holds a stale twin), puts its root at the FRONT of ``sys.path`` for the test only, and mutates the files it
wrote. Each generator appends one line to a counter file per run, so a gate can count the generations that
actually happened (X1: the activation is measured on the side of the generator, not reported by the memo).
"""
from __future__ import annotations

import importlib
import sys
import textwrap
import uuid
from pathlib import Path
from types import ModuleType

from . import _traced_memo_api as api

ALPHA = '''
"""The alpha module: a generator, its helpers, and code it never runs."""
import dataclasses

import numpy as np

from {memo_module} import traced_memo

SCALE = 3.0
UNREAD = 7.0
COUNTER = {counter!r}


def _count(tag):
    with open(COUNTER, "a") as handle:
        handle.write(tag + "\\n")


def helper(x):
    """Doc of a traced helper."""
    return x * SCALE


def unused(x):
    """Never called by any generator."""
    return x - 1.0


def nested_user(x):
    def inner(y):
        return y + 1.0
    squares = sum(v * v for v in (1.0, 2.0))
    return inner(x) + (lambda z: z * 0.5)(squares)


class Box:
    def method(self, x):
        return x + 2.0

    def other(self, x):
        return x - 2.0


@dataclasses.dataclass(frozen=True)
class Result:
    k: float
    field: np.ndarray
    n: int
    ok: bool


@dataclasses.dataclass(frozen=True)
class Spec:
    n: int

    def __post_init__(self):
        if self.n < 0:
            raise ValueError("Spec: n is a count")


@traced_memo
def generate(x: float) -> float:
    _count("generate")
    return helper(x) + Box().method(x) + nested_user(x) + float(np.sum(np.ones(2)))


@traced_memo
def solve(n: int, scale: float = 2.0, *, offset: float = 0.0) -> Result:
    _count("solve")
    field = np.linspace(0.0, scale, n) + offset
    return Result(float(field.sum()) / 3.0, field, n, True)


@traced_memo
def read_spec(spec: Spec) -> float:
    _count("read_spec")
    return float(spec.n) / 7.0


@traced_memo
def boom(x: float) -> float:
    _count("boom")
    raise ValueError(f"boom at {{x}}")


def as_array(value):
    return np.asarray(value, dtype=float)


@traced_memo(canonical={{"x": as_array}})
def total(x) -> float:
    _count("total")
    return float(np.sum(x)) + (0.5 if isinstance(x, np.ndarray) else 0.0)


@traced_memo
def environment_probe(name: str) -> str:
    import os
    return os.environ.get(name, "<unset>")


@traced_memo
def divide(x) -> float:
    _count("divide")
    return float((np.asarray(x) / np.asarray(3, dtype=np.asarray(x).dtype)).sum())
'''

BETA = '''
"""The beta module: a parent generator whose answer reads a child entry."""
from {memo_module} import traced_memo

from .alpha import _count

FACTOR = 5.0


@traced_memo
def inner(x: float) -> float:
    _count("inner")
    return x * FACTOR


@traced_memo
def outer(x: float) -> float:
    _count("outer")
    return inner(x) + 1.0


@traced_memo
def outer_twice(x: float) -> float:
    _count("outer_twice")
    return inner(x) + inner(x)
'''

GAMMA = '''
"""The gamma module: imported by nobody."""
NEVER = 1.0


def never(x):
    return x
'''

MARKER = '''
"""A module whose import writes a marker file: the validation-imports-nothing witness."""
from pathlib import Path

Path({marker!r}).write_text("imported")
from {memo_module} import traced_memo


@traced_memo
def plain(x: float) -> float:
    return x + 1.0
'''


DELTA = '''
"""The delta module: Python 3.12+ code objects generated from signatures (PEP 649 annotations, PEP 695 generics)."""
import typing


def annotated(x: float, y: "list[int]") -> float:
    return x + len(y)


def generic[T](value: T) -> T:
    return value


class Holder:
    width: float = 1.0

    def scaled[S](self, factor: S) -> S:
        return factor


def reads_annotations(x):
    hints = typing.get_type_hints(annotated)
    generic.__type_params__
    typing.get_type_hints(Holder)
    return Holder().scaled(generic(x)) + len(hints)
'''


class Package:
    """A fresh generator package under ``root``; ``edit`` mutates the copy, never a tracked file."""

    def __init__(self, tmp_path: Path) -> None:
        self.root = tmp_path / "src"
        self.name = f"memo_fixture_{uuid.uuid4().hex[:12]}"
        self.dir = self.root / self.name
        self.dir.mkdir(parents=True)
        self.counter = tmp_path / "counter.txt"
        self.counter.write_text("")
        self.marker = tmp_path / "marker.txt"
        fill = {"memo_module": api.MODULE, "counter": str(self.counter), "marker": str(self.marker)}
        (self.dir / "__init__.py").write_text("")
        for module, template in (("alpha", ALPHA), ("beta", BETA), ("gamma", GAMMA), ("marker", MARKER), ("delta", DELTA)):
            (self.dir / f"{module}.py").write_text(textwrap.dedent(template).format(**fill))
        sys.path.insert(0, str(self.root))
        importlib.invalidate_caches()

    def close(self) -> None:
        while str(self.root) in sys.path:
            sys.path.remove(str(self.root))
        for key in [k for k in sys.modules if k == self.name or k.startswith(self.name + ".")]:
            del sys.modules[key]

    def path(self, module: str) -> Path:
        return self.dir / f"{module}.py"

    def module(self, module: str) -> ModuleType:
        return importlib.import_module(f"{self.name}.{module}")

    def edit(self, module: str, old: str, new: str) -> None:
        """Replace exactly one occurrence of ``old`` (a gate whose edit misses its target must fail loudly)."""
        text = self.path(module).read_text()
        if text.count(old) != 1:
            raise AssertionError(f"edit target {old!r} occurs {text.count(old)} times in {module}.py")
        self.path(module).write_text(text.replace(old, new))

    def generations(self, tag: str | None = None) -> int:
        lines = self.counter.read_text().split()
        return len(lines) if tag is None else lines.count(tag)

    def relpath(self, module: str) -> str:
        return f"{self.name}/{module}.py"
