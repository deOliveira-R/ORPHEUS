r"""G04o: the legacy-leaf arm of ``FunctionSpace`` is a pinned set of members and no member reads ``self.axes`` (N04, spec 7.5).

``ast`` over class ``FunctionSpace`` in ``numerics/space.py``. (1) No member reads ``self.axes`` (13 members on the base tip: the 12 that compare it with ``None``
and ``axis``, which reads it into a local). (2) The members carrying the legacy-leaf tag ``ELEGANCE-DEBT[guard] #759 legacy-leaf`` in their source segment are exactly ``PINNED``:
the 12 members the base tip tests ``self.axes is None`` in, measured here by the same parse on ``cf924aaa`` (``[M]``, 12 of 41 members); the unit keeps each as its
leaf arm with no branch and the tag. A 13th tagged member, a pinned member that loses its tag, or any read of ``self.axes`` reddens. G07c (N07) asserts the set empty. ``_overlay_forms`` and ``_with_overlay`` were among the base's 12 and are deleted by the unit (G04g), so ``PINNED`` holds the 10 that survive. The tag has one spelling, ``ELEGANCE-DEBT[guard] #759 legacy-leaf`` (issue #759, "FunctionSpace keeps a leaf branch for the arms that are not tensor spaces yet", retiring at N07); the earlier spelling ``ELEGANCE-DEBT[legacy-leaf]`` and a tag citing another issue are not tags of this set and a member carrying only one of them reddens as a lost tag.

Not a measurement of the final set: the spec reads the pinned set "at the unit's last commit"; the builder may amend ``PINNED`` only by a ruling stated in the
commit (a member that the carve turns into a pure abstract method is removed from the set). First red at the base tip: 13 reads of ``self.axes`` and 0 tagged members.
Mutations on source text through the gate's predicate (positive controls in the file): add ``if self.axes is None: pass`` to a member; tag a 13th; drop a tag. Exact on any host.

Second leg (ruling W3 of the N04 post-review pass): the doors outside ``numerics/space.py``. A door is an ``isinstance`` whose class argument is ``TensorSpace`` (or a tuple holding it) or a
``type(x) is TensorSpace`` test. By ``ast`` over ``orpheus/`` (outside ``space.py``, whose own two doors are ``TensorSpace.__mul__``, the binary-operator protocol, untagged, and
``TensorSpace.positions_spanned_by``, tagged in place) every door sits in a function pinned here by ``(file, qualified function)`` in one of two sets. ``TAGGED_DOORS`` tell a tensor space from a legacy leaf
and nothing else and carry ``ELEGANCE-DEBT[guard] #759 legacy-leaf`` on the door's own line or the two lines above it; N07's successor asserts it empty. ``UNTAGGED_DOORS`` are the doors that are not
leaf debt, each with its reason (a direct sum, the arm that stays at N07, reaches the False side, or the function refuses what is not a tensor space with a ``TypeError``). The census by ``ast`` at the unit's tip read 14 doors: 2 in ``space.py`` and 12 outside it, one of which (``BulkField.spatial_moments_per_axis_of``, class (iii)) was replaced by the ``axis`` verb every space has; the 11 that remain, each in its own function, are 1 tagged and 10 untagged with reasons.
A door in a function of neither set, a pinned function with no door, or a tag on an untagged door reddens. Mutation: a new untagged door in any file, and a tag dropped from the tagged door (planted rows in the file, and the mutations run on a copy).
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.foundation

ROOT = Path(__file__).resolve().parents[4]
TAG_RE = re.compile(r"ELEGANCE-DEBT\[guard\] #759 legacy-leaf")
DOOR_RE = TAG_RE
#: The doors outside ``space.py`` that tell a tensor space from a legacy leaf and nothing else: tagged, and asserted empty by N07's successor.
TAGGED_DOORS = frozenset({("orpheus/transport/fields/_bases.py", "BulkField._spatial_moment_tail_of")})
#: The doors that are not leaf debt, each with the reason it stays untagged. A tensor space is told from a direct sum (the second arm that stays at N07) or the
#: function refuses what is not a tensor space with a ``TypeError``, so the door survives the leaves. A new door anywhere is untagged and not here, and reddens.
UNTAGGED_DOORS = {
    ("orpheus/numerics/spaces/moment_head.py", "moment_head_axis"): "refuses a space that is not a one-axis moment head with a TypeError naming it; a direct sum is refused here after N07",
    ("orpheus/transport/fields/_bases.py", "BulkField.compose_spatial_moments"): "takes a tensor space and refuses anything else with a TypeError: a direct sum has no labelled factor to append a tail to",
    ("orpheus/transport/fields/_bases.py", "BulkField.ng"): "reads the energy axis off a labelled factor and raises a TypeError where there is none, which a direct sum has not",
    ("orpheus/transport/fields/_bases.py", "AngularField.N"): "reads the ordinate count off the leading labelled factor and raises a TypeError where there is none, which a direct sum has not",
    ("orpheus/transport/fields/_bases.py", "MomentField.head"): "takes the head from the labelled factors of a tensor space and refuses any other field space with a TypeError",
    ("orpheus/transport/fields/harmonic_moment_flux.py", "HarmonicMomentFlux.scalar_flux"): "an assert narrowing the type for the checker, not a guard: a moment field's space is a tensor space; it goes when the field's space is typed",
    ("orpheus/transport/frames/harmonic_frame.py", "HarmonicFrame.for_space"): "takes an angular tensor space and refuses anything else with a TypeError: a direct sum names no quadrature",
    ("orpheus/transport/frames/harmonic_frame.py", "HarmonicFrame.moment_space_on"): "takes an angular tensor space and refuses anything else with a TypeError: a direct sum names no cell axes",
    ("orpheus/transport/operators/_energy_conformity.py", "assert_energy_extent_conforms"): "the composite spaces (FullFieldSpace, CoupledSpace; the direct-sum arm from N05) land on the False side and pass unchecked, a measured inertness; it retires when the direct sum has its own energy check, not with the leaves",
    ("orpheus/transport/operators/fission.py", "FissionOperator.from_solver_data"): "takes a composite whose interior is a tensor space and refuses a non-tensor interior with a TypeError: a direct sum names no quadrature",
}
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
    base = f"class FunctionSpace:\n    def unit(self):\n        # ELEGANCE-DEBT[guard] #759 legacy-leaf: retires at N07\n        return 1\n    def other(self):\n        return 2\n"
    _require(read(base) == (set(), {"unit"}), f"a clean synthetic source reads {read(base)}")
    _require(read(base.replace("return 2", "return self.axes"))[0] == {"other"}, "a read of self.axes is not found")
    _require(read(base.replace("return 2", "return 2  # ELEGANCE-DEBT[guard] #759 legacy-leaf: retires at N07"))[1] == {"unit", "other"}, "a 13th tagged member is not found")
    for other in ("ELEGANCE-DEBT[legacy-leaf] N07", "ELEGANCE-DEBT[guard] #707 legacy-leaf: retires at N07"):
        _require(read(base.replace("ELEGANCE-DEBT[guard] #759 legacy-leaf: retires at N07", other))[1] == set(), f"the spelling {other!r} is accepted as the tag")


def test_leaf_branches_are_pinned() -> None:
    reads, tagged = read((ROOT / "orpheus/numerics/space.py").read_text())
    _require(not reads, f"{len(reads)} members of FunctionSpace read self.axes: {sorted(reads)}")
    _require(tagged == PINNED, f"the legacy-leaf tagged members differ from the pinned set: missing {sorted(PINNED - tagged)}, extra {sorted(tagged - PINNED)}")


def _names(node: ast.AST) -> list[str]:
    if isinstance(node, ast.Name):
        return [node.id]
    if isinstance(node, ast.Attribute):
        return [node.attr]
    if isinstance(node, ast.Tuple):
        return [n for e in node.elts for n in _names(e)]
    return []


def read_doors(source: str) -> tuple[set[str], set[str]]:
    """``(functions with an untagged door, functions with a tagged door)`` of a module: a door is ``isinstance(x, TensorSpace)`` (or a tuple holding it) or ``type(x) is TensorSpace``."""
    lines = source.splitlines()
    untagged: set[str] = set()
    tagged: set[str] = set()

    def is_door(n: ast.AST) -> bool:
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == "isinstance" and len(n.args) == 2:
            return "TensorSpace" in _names(n.args[1])
        if isinstance(n, ast.Compare) and isinstance(n.left, ast.Call) and getattr(n.left.func, "id", "") == "type":
            return any(isinstance(o, (ast.Is, ast.IsNot)) for o in n.ops) and any("TensorSpace" in _names(c) for c in n.comparators)
        return False

    def walk(node: ast.AST, stack: list[str], function: str | None) -> None:
        for child in ast.iter_child_nodes(node):
            if isinstance(child, ast.ClassDef):
                walk(child, stack + [child.name], function)
                continue
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
                walk(child, stack + [child.name], ".".join(stack + [child.name]))
                continue
            if is_door(child):
                function = function or "<module>"
                window = "\n".join(lines[max(child.lineno - 3, 0) : child.lineno])
                (tagged if DOOR_RE.search(window) else untagged).add(function)
            walk(child, stack, function)

    walk(ast.parse(source), [], None)
    return untagged, tagged


def test_the_door_instrument_reads_each_mutation() -> None:
    tag = "        # ELEGANCE-DEBT[guard] #759 legacy-leaf: a leaf has no axes; retires at N07\n"
    clean = f"class C:\n    def f(self, x):\n{tag}        if isinstance(x, TensorSpace):\n            return 1\n"
    _require(read_doors(clean) == (set(), {"C.f"}), f"a tagged door reads {read_doors(clean)}")
    _require(read_doors(clean.replace(tag, ""))[0] == {"C.f"}, "an untagged door is not found")
    _require(read_doors("def g(x):\n    return isinstance(x, (A, TensorSpace))\n")[0] == {"g"}, "a tuple door is not found")
    _require(read_doors("def g(x):\n    return type(x) is TensorSpace\n")[0] == {"g"}, "a type-is door is not found")
    _require(read_doors("def g(x):\n    return isinstance(x, Member)\n") == (set(), set()), "a non-door isinstance is found")
    _require(read_doors("ok = isinstance(x, TensorSpace)\n")[0] == {"<module>"}, "a door outside every function is not found")


def test_the_doors_outside_space_py_are_tagged_and_pinned() -> None:
    tagged_found: set[tuple[str, str]] = set()
    untagged_found: set[tuple[str, str]] = set()
    for path in sorted((ROOT / "orpheus").rglob("*.py")):
        rel = str(path.relative_to(ROOT))
        if rel == "orpheus/numerics/space.py":
            continue
        bad, good = read_doors(path.read_text())
        untagged_found |= {(rel, f) for f in bad}
        tagged_found |= {(rel, f) for f in good}
    _require(
        tagged_found == TAGGED_DOORS,
        f"the tagged doors differ from the pinned set: missing {sorted(TAGGED_DOORS - tagged_found)}, extra {sorted(tagged_found - TAGGED_DOORS)}",
    )
    _require(
        untagged_found == set(UNTAGGED_DOORS),
        f"the untagged doors differ from the pinned set with its reasons: missing {sorted(set(UNTAGGED_DOORS) - untagged_found)}, new {sorted(untagged_found - set(UNTAGGED_DOORS))}",
    )
