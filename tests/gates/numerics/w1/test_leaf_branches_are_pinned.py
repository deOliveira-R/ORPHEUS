r"""G04o: the legacy-leaf arm of ``FunctionSpace`` is a pinned set of members and no member reads ``self.axes`` (N04, spec 7.5).

``ast`` over class ``FunctionSpace`` in ``numerics/space.py``. (1) No member reads ``self.axes`` (13 members on the base tip: the 12 that compare it with ``None``
and ``axis``, which reads it into a local). (2) The members carrying the legacy-leaf tag ``ELEGANCE-DEBT[legacy-leaf]`` in their source segment are exactly ``PINNED``:
the 12 members the base tip tests ``self.axes is None`` in, measured here by the same parse on ``cf924aaa`` (``[M]``, 12 of 41 members); the unit keeps each as its
leaf arm with no branch and the tag. A 13th tagged member, a pinned member that loses its tag, or any read of ``self.axes`` reddens. G07c (N07) asserts the set empty. ``_overlay_forms`` and ``_with_overlay`` were among the base's 12 and are deleted by the unit (G04g), so ``PINNED`` holds the 10 that survive. The tag's bracket word and issue are N07's opening ruling (C5: ``ELEGANCE-DEBT[guard] #<issue>``), not this gate's: it accepts either spelling of a tag that names ``legacy-leaf``.

Not a measurement of the final set: the spec reads the pinned set "at the unit's last commit"; the builder may amend ``PINNED`` only by a ruling stated in the
commit (a member that the carve turns into a pure abstract method is removed from the set). First red at the base tip: 13 reads of ``self.axes`` and 0 tagged members.
Mutations on source text through the gate's predicate (positive controls in the file): add ``if self.axes is None: pass`` to a member; tag a 13th; drop a tag. Exact on any host.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.foundation

ROOT = Path(__file__).resolve().parents[4]
TAG_RE = re.compile(r"ELEGANCE-DEBT\[legacy-leaf\]|ELEGANCE-DEBT\[[a-z-]+\][^\n]*legacy-leaf")
PINNED = frozenset(
    {
        "__post_init__",
        "identity_key",
        "has_coordinate_cone",
        "_axis_index",
        "_axis_block",
        "without_axis",
        "with_axis",
        "positions_spanned_by",
        "unit",
        "_resolved_metric",
    }
)


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def read(source: str) -> tuple[set[str], set[str]]:
    """``(members reading self.axes, members whose source segment carries the tag)`` of class FunctionSpace."""
    lines = source.splitlines()
    cls = [n for n in ast.walk(ast.parse(source)) if isinstance(n, ast.ClassDef) and n.name == "FunctionSpace"][0]
    reads: set[str] = set()
    tagged: set[str] = set()
    for m in cls.body:
        if not isinstance(m, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        if any(isinstance(a, ast.Attribute) and a.attr == "axes" and isinstance(a.value, ast.Name) and a.value.id == "self" for a in ast.walk(m)):
            reads.add(m.name)
        segment = "\n".join(lines[m.lineno - 1 : m.end_lineno])
        if TAG_RE.search(segment):
            tagged.add(m.name)
    return reads, tagged


def test_the_instrument_reads_each_mutation() -> None:
    base = f"class FunctionSpace:\n    def unit(self):\n        # ELEGANCE-DEBT[legacy-leaf] N07\n        return 1\n    def other(self):\n        return 2\n"
    _require(read(base) == (set(), {"unit"}), f"a clean synthetic source reads {read(base)}")
    _require(read(base.replace("return 2", "return self.axes"))[0] == {"other"}, "a read of self.axes is not found")
    _require(read(base.replace("return 2", "return 2  # ELEGANCE-DEBT[legacy-leaf] N07"))[1] == {"unit", "other"}, "a 13th tagged member is not found")


def test_leaf_branches_are_pinned() -> None:
    reads, tagged = read((ROOT / "orpheus/numerics/space.py").read_text())
    _require(not reads, f"{len(reads)} members of FunctionSpace read self.axes: {sorted(reads)}")
    _require(tagged == PINNED, f"the legacy-leaf tagged members differ from the pinned set: missing {sorted(PINNED - tagged)}, extra {sorted(tagged - PINNED)}")
