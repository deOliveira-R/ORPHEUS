r"""G04g: the axis-built arm of ``FunctionSpace`` and the overlay are gone (N04, spec 7.5, 6.5 "Outcome").

``ast`` over the three source files the unit edits. In ``numerics/space.py``: class ``FunctionSpace`` has no ``axes`` field, no ``of_axes``, no ``__mul__``
(the class ``TensorSpace`` defines its own ``of_axes`` and ``__mul__``), and none of ``_axis_entries``, ``_entry_of_axis``, ``_require_positioned_over``,
``_overlay_forms``, ``_with_overlay``, ``_axes_metric``, ``_positioned_forms``, ``_tensor_product_factored_metric``, ``from_factors``, ``TensorProductSpace`` is defined anywhere
in the file; ``PositionedMetric`` is defined nowhere in ``numerics/metric.py``; ``_dressed`` nowhere in ``numerics/frame_faces.py``.

First red at the base tip ``cf924aaa``: the tip defines 10 of the 14 names and the field, constructor and dunder (the opening's count was 10 for the first
set; this gate adds ``_axes_metric``, ``_positioned_forms``, ``_tensor_product_factored_metric``, ``from_factors``, ``_dressed``). Mutations, run on a copy of the source text
through the gate's own predicate (positive controls inside the gate): append ``def _axis_entries(): pass`` reddens; re-adding the ``axes`` field to ``FunctionSpace`` reddens;
appending ``class PositionedMetric: pass`` to the ``metric.py`` text reddens. Exact on any host.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

pytestmark = pytest.mark.foundation

ROOT = Path(__file__).resolve().parents[4]
RETIRED_IN_SPACE = (
    "_axis_entries",
    "_entry_of_axis",
    "_require_positioned_over",
    "_overlay_forms",
    "_with_overlay",
    "_axes_metric",
    "_positioned_forms",
    "_tensor_product_factored_metric",
    "from_factors",
    "TensorProductSpace",
)


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _defined(tree: ast.AST) -> set[str]:
    out: set[str] = set()
    for n in ast.walk(tree):
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            out.add(n.name)
        elif isinstance(n, ast.AnnAssign) and isinstance(n.target, ast.Name):
            out.add(n.target.id)
        elif isinstance(n, ast.Assign):
            out |= {t.id for t in n.targets if isinstance(t, ast.Name)}
    return out


def _class_members(tree: ast.AST, name: str) -> set[str]:
    for n in ast.walk(tree):
        if isinstance(n, ast.ClassDef) and n.name == name:
            return _defined(ast.Module(body=n.body, type_ignores=[])) - {n.name}
    return set()


def _space_findings(source: str) -> list[str]:
    tree = ast.parse(source)
    found = [f"{n} is defined" for n in RETIRED_IN_SPACE if n in _defined(tree)]
    base = _class_members(tree, "FunctionSpace")
    found += [f"FunctionSpace.{n} is defined" for n in ("axes", "of_axes", "__mul__") if n in base]
    return found


def _requires_tensor_space_own(source: str) -> list[str]:
    own = _class_members(ast.parse(source), "TensorSpace")
    return [f"TensorSpace.{n} is not defined in its own body" for n in ("of_axes", "__mul__") if n not in own]


def test_the_instrument_finds_each_planted_name() -> None:
    clean = "class FunctionSpace:\n    x: int\nclass TensorSpace:\n    def of_axes(self): ...\n    def __mul__(self): ...\n"
    _require(_space_findings(clean) == [], "the predicate reads a clean synthetic source as dirty")
    _require(_space_findings(clean + "\ndef _axis_entries(): pass") == ["_axis_entries is defined"], "a planted _axis_entries is not found")
    _require("FunctionSpace.axes is defined" in _space_findings("class FunctionSpace:\n    axes: tuple = ()\n"), "a planted axes field is not found")


def test_space_py_carries_no_axes_arm_and_no_overlay() -> None:
    src = (ROOT / "orpheus/numerics/space.py").read_text()
    found = _space_findings(src) + _requires_tensor_space_own(src)
    _require(not found, f"{len(found)} findings in numerics/space.py: {found}")


def test_positioned_metric_and_dressed_are_gone() -> None:
    metric = ast.parse((ROOT / "orpheus/numerics/metric.py").read_text())
    _require("PositionedMetric" not in _defined(metric), "PositionedMetric is still defined in numerics/metric.py")
    _require("PositionedMetric" in _defined(ast.parse("class PositionedMetric: pass")), "the predicate does not find a planted PositionedMetric")
    faces = ast.parse((ROOT / "orpheus/numerics/frame_faces.py").read_text())
    _require("_dressed" not in _defined(faces), "_dressed is still defined in numerics/frame_faces.py")
