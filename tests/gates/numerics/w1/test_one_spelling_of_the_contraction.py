r"""G04l: the metric verbs of ``TensorSpace`` are calls of ``flat_on`` and ``sharp_on``; no other body spells the walk (N04, spec 7.5; law L10).

``ast`` over class ``TensorSpace`` in ``numerics/space.py``: ``flat``, ``apply_metric``, ``apply_inverse_metric`` and ``integral`` are each defined in its body and each
reaches ``flat_on`` or ``sharp_on`` through ``self`` calls (transitively through other methods of the class); ``flat_on`` and ``sharp_on`` are the only methods that call
``.apply``, ``.apply_inverse``, ``._walk`` or ``.apply_block``. The two leaf bodies of ``FunctionSpace`` (``_diagonal_apply_metric``, ``_diagonal_apply_inverse_metric``)
are the tagged legacy-leaf arm and outside the scope (class ``FunctionSpace`` is not read).

First red at the base tip: no class ``TensorSpace`` (the test names the missing class). Mutations, on source text through the gate's predicate (positive controls in the
file): inline a walk in ``apply_metric`` reddens; ``integral`` that stops delegating reddens; a second walk helper calling ``._walk`` reddens. Exact on any host.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

pytestmark = pytest.mark.foundation

ROOT = Path(__file__).resolve().parents[4]
VERBS = ("flat", "apply_metric", "apply_inverse_metric", "integral")
WALKERS = {"apply", "apply_inverse", "_walk", "apply_block"}
OWNERS = ("flat_on", "sharp_on")


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def findings(source: str) -> list[str]:
    cls = [n for n in ast.walk(ast.parse(source)) if isinstance(n, ast.ClassDef) and n.name == "TensorSpace"]
    if not cls:
        return ["class TensorSpace is not defined"]
    methods = {m.name: m for m in cls[0].body if isinstance(m, ast.FunctionDef)}
    self_calls = {
        name: {c.func.attr for c in ast.walk(m) if isinstance(c, ast.Call) and isinstance(c.func, ast.Attribute) and isinstance(c.func.value, ast.Name) and c.func.value.id == "self"}
        for name, m in methods.items()
    }

    def reaches(name: str, seen: frozenset[str] = frozenset()) -> bool:
        if name in OWNERS:
            return True
        if name in seen or name not in self_calls:
            return False
        return any(reaches(t, seen | {name}) for t in self_calls[name])

    out: list[str] = []
    for v in VERBS:
        if v not in methods:
            out.append(f"{v} is not defined in the body of TensorSpace")
        elif not reaches(v):
            out.append(f"{v} does not reach flat_on or sharp_on")
    for owner in OWNERS:
        if owner not in methods:
            out.append(f"{owner} is not defined in the body of TensorSpace")
    for name, m in methods.items():
        if name in OWNERS:
            continue
        for c in ast.walk(m):
            if isinstance(c, ast.Call) and isinstance(c.func, ast.Attribute) and c.func.attr in WALKERS:
                out.append(f"{name} spells a walk: .{c.func.attr}( at line {c.lineno}")
    return out


_GOOD = (
    "class TensorSpace:\n"
    "    def flat_on(self, x, labels): return self._m._walk(x)\n"
    "    def sharp_on(self, x, labels): return self._m._walk(x)\n"
    "    def flat(self, x): return self.flat_on(x, ())\n"
    "    def apply_metric(self, x): return self.flat_on(x, ())\n"
    "    def apply_inverse_metric(self, x): return self.sharp_on(x, ())\n"
    "    def integral(self, x): return self.flat(x)\n"
)


def test_the_instrument_reads_a_good_source_clean_and_each_mutation_red() -> None:
    _require(findings(_GOOD) == [], f"a clean synthetic source reads {findings(_GOOD)}")
    inline = _GOOD.replace("def apply_metric(self, x): return self.flat_on(x, ())", "def apply_metric(self, x): return self._m.apply(x)")
    _require(any("apply_metric" in f for f in findings(inline)), "an inlined walk in apply_metric is not found")
    stops = _GOOD.replace("def integral(self, x): return self.flat(x)", "def integral(self, x): return x.sum()")
    _require(any("integral does not reach" in f for f in findings(stops)), "an integral that stops delegating is not found")
    helper = _GOOD + "    def _again(self, x): return self._m._walk(x)\n"
    _require(any("_again" in f for f in findings(helper)), "a second walk helper is not found")


def test_the_metric_verbs_of_tensor_space_delegate() -> None:
    found = findings((ROOT / "orpheus/numerics/space.py").read_text())
    _require(not found, f"{len(found)} findings in class TensorSpace: {found}")
