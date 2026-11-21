"""Gate G7 of unit ``memo-digest`` (specification part B 3.4, 4.2): ``_MEMO_SOURCES`` is the static import closure of the memo over its
two packages, the same set in every process on one tree, complete for a function-body import, and the memo pin is sensitive to a
byte moved between two adjacent files.

Tree ``T`` (``_memo_tree.build_tree``) holds copies of ``orpheus/__init__.py``, ``orpheus/memo/*.py`` and ``orpheus/values/*.py`` plus
``orpheus/values/stray.py``, a values module the memo does not import. Children: ordinary ones run ``-c`` with ``T`` first on
``sys.path`` and assert ``orpheus.__file__`` under ``T``; generating ones run the real ``_boot`` (``-P``, cwd outside ``T``).

Legs: (a) the order leg (two cold children, one importing ``stray`` first), (b) the generating leg (derived set equals the pinned files
under the two roots, both derived by the gate, the known member ``values/digest.py`` required in both), (c) the late-import leg, (d) the
boundary leg (a line moved across two adjacent files changes the pin). The arms (1) to (4) are copies of ``traced_memo.py`` edited
through ``_memo_tree.mutate`` (``Uninstallable`` unless the anchor occurs once): [R] at base, owed after the build.

First red at base: ``orpheus/memo`` does not exist (``build_tree`` raises ``Uninstallable`` on an empty tree); a scratch probe of the same
legs on the numerics layout reads the hand tuple of 4 files (commit body).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from tests.gates.memo import _memo_tree
from tests.gates.memo._memo_tree import build_tree, cold_child, generating_child, mutate

pytestmark = pytest.mark.foundation

_EIGHT = {"memo/__init__.py", "memo/_import_time.py", "memo/_traced_memo_boot.py", "memo/traced_memo.py",
          "values/__init__.py", "values/content.py", "values/digest.py", "values/scalars.py"}
_STRAY = {"stray.py": "STRAY = 1\n"}
_LIST = "import json; import orpheus.memo.traced_memo as t; print(json.dumps([str(p) for p in t._MEMO_SOURCES]))\n"


def _rel(tree: Path, paths: list[str]) -> set[str]:
    base = (tree / "orpheus").resolve()
    return {Path(p).resolve().relative_to(base).as_posix() for p in paths}


def _order_leg(tree: Path) -> tuple[list[str], list[str]]:
    first = json.loads(cold_child(tree, "import orpheus.values.stray\n" + _LIST).splitlines()[-1])
    second = json.loads(cold_child(tree, _LIST).splitlines()[-1])
    return first, second


_GENERATING = '''
derived = [os.path.realpath(str(p)) for p in traced_memo._MEMO_SOURCES]
roots = [os.path.realpath(os.path.join(R, "orpheus", d)) + os.sep for d in ("memo", "values")]
pinned = sorted(k for k in module.PINNED_SOURCES if k.startswith(tuple(roots)))
facts.update(derived=derived, pinned=pinned, sorted_ok=derived == sorted(derived),
             digest_in_both=all(os.path.realpath(os.path.join(R, "orpheus", "values", "digest.py")) in s for s in (derived, pinned)),
             running=traced_memo._RUNNING_MEMO_DIGEST is not None)
'''


def test_the_set_is_the_same_in_two_import_orders_and_is_the_eight_files(tmp_path: Path) -> None:
    """Leg (a). First red at the second pass's derivation (the loaded modules): 9 files against 8, the stray module among them."""
    tree = build_tree(tmp_path / "t", extra_values=_STRAY)
    first, second = _order_leg(tree)
    assert first == second and first == sorted(first), (first, second)
    assert _rel(tree, first) == _EIGHT and len(first) == 8
    assert "values/stray.py" not in _rel(tree, first)


@pytest.mark.rests_on('tests/gates/memo/test_memo_boot_loads_the_convention.py::test_a_generating_process_holds_one_module_object_of_the_convention_read_once', 'tests/gates/memo/test_memo_sources.py::test_the_set_is_the_same_in_two_import_orders_and_is_the_eight_files')
def test_the_set_equals_the_pinned_files_under_the_two_roots_in_a_generating_process(tmp_path: Path) -> None:
    """Leg (b). Both sides are derived by the gate (no literal of the eight paths); the known member is required in both."""
    facts = generating_child(build_tree(tmp_path / "t", extra_values=_STRAY), _GENERATING, tmp_path / "w")
    assert facts["derived"] == facts["pinned"] and facts["sorted_ok"] and facts["digest_in_both"] and facts["running"], facts


_LATE = '''
import importlib
units = os.path.realpath(os.path.join(R, "orpheus", "values", "units.py"))
facts["in_set"] = units in [os.path.realpath(str(p)) for p in traced_memo._MEMO_SOURCES]
facts["in_pinned"] = units in module.PINNED_SOURCES
with open(units, "w") as f:
    f.write("U = 2\\n")
import orpheus.values.units as u
facts["compiled_from_pin"] = u.U == 1 and module.PINNED_SOURCES[units] == b"U = 1\\n"
'''


@pytest.mark.rests_on('tests/gates/memo/test_memo_boot_loads_the_convention.py::test_a_generating_process_holds_one_module_object_of_the_convention_read_once', 'tests/gates/memo/test_memo_sources.py::test_the_set_is_the_same_in_two_import_orders_and_is_the_eight_files')
def test_a_function_body_import_is_in_the_set_and_pinned_and_compiled_from_the_pinned_bytes(tmp_path: Path) -> None:
    """Leg (c). First red: ``_pinned_memo_digest`` without the completion raises ``KeyError`` at the memo's import (the child exits 1)."""
    tree = build_tree(tmp_path / "t", extra_values={"units.py": "U = 1\n"})
    memo = tree / "orpheus" / "memo" / "traced_memo.py"
    memo.write_text(memo.read_text(encoding="utf-8") + "\n\ndef _planted_late_import() -> None:\n    from orpheus.values import units  # noqa: F401\n")
    facts = generating_child(tree, _LATE, tmp_path / "w")
    assert facts["in_set"] and facts["in_pinned"] and facts["compiled_from_pin"], facts


def _memo_digest_in(tree: Path) -> str:
    return cold_child(tree, "import orpheus.memo.traced_memo as t; print(t._memo_digest())\n").splitlines()[-1]


@pytest.mark.rests_on('tests/gates/memo/test_memo_sources.py::test_the_set_is_the_same_in_two_import_orders_and_is_the_eight_files')
def test_a_line_moved_across_two_adjacent_files_changes_the_memo_pin(tmp_path: Path) -> None:
    """Leg (d). The concatenated bytes of ``content.py`` + ``digest.py`` are equal in the two trees; the framed pin differs. First red:
    the bare concatenation (equal digests, measured 3.4)."""
    one = build_tree(tmp_path / "one")
    two = build_tree(tmp_path / "two")
    probe = "_PROBE_SWITCH = True\n"
    content = one / "orpheus" / "values" / "content.py"
    assert content.read_bytes().endswith(b"\n")
    content.write_bytes(content.read_bytes() + probe.encode())
    digest = two / "orpheus" / "values" / "digest.py"
    digest.write_bytes(probe.encode() + digest.read_bytes())
    assert (one / "orpheus/values/content.py").read_bytes() + (one / "orpheus/values/digest.py").read_bytes() == (
        (two / "orpheus/values/content.py").read_bytes() + (two / "orpheus/values/digest.py").read_bytes())
    assert _memo_digest_in(one) != _memo_digest_in(two)


# ---------------------------------------------------------------- the arms (each: honest first, then the mutant; [R] at base)
def test_arm_1_a_directory_listing_derivation_is_red_on_the_order_leg(tmp_path: Path) -> None:
    tree = build_tree(tmp_path / "t", extra_values=_STRAY)
    assert _rel(tree, _order_leg(tree)[0]) == _EIGHT
    listing = ("_MEMO_SOURCES = tuple(sorted(p for d in ('memo', 'values') for p in "
               "(Path(__file__).resolve().parent.parent / d).glob('*.py')))")
    mutate(tree / "orpheus" / "memo" / "traced_memo.py", "_MEMO_SOURCES = _memo_sources()", listing)
    assert "values/stray.py" in _rel(tree, _order_leg(tree)[1])


def test_arm_2_a_derived_set_missing_a_pinned_file_is_red_on_the_generating_leg(tmp_path: Path) -> None:
    ok = generating_child(build_tree(tmp_path / "ok", extra_values=_STRAY), _GENERATING, tmp_path / "w1")
    assert ok["derived"] == ok["pinned"]
    bad_tree = build_tree(tmp_path / "bad", extra_values=_STRAY, memo_edit=(
        "traced_memo.py", "_MEMO_SOURCES = _memo_sources()",
        '_MEMO_SOURCES = tuple(p for p in _memo_sources() if p.name != "scalars.py")'))
    bad = generating_child(bad_tree, _GENERATING, tmp_path / "w2")
    assert bad["derived"] != bad["pinned"]


def test_arm_3_a_pin_without_the_length_prefix_is_equal_across_the_boundary(tmp_path: Path) -> None:
    one = build_tree(tmp_path / "one", memo_edit=("traced_memo.py", 'len(blob).to_bytes(8, "big") + ', ""))
    two = build_tree(tmp_path / "two", memo_edit=("traced_memo.py", 'len(blob).to_bytes(8, "big") + ', ""))
    probe = b"_PROBE_SWITCH = True\n"
    content, digest = one / "orpheus/values/content.py", two / "orpheus/values/digest.py"
    content.write_bytes(content.read_bytes() + probe)
    digest.write_bytes(probe + digest.read_bytes())
    assert _memo_digest_in(one) == _memo_digest_in(two), "unframed: the arm must read EQUAL (the defect it re-introduces)"


def test_arm_4_a_pin_without_the_completion_raises_a_key_error_at_the_memo_import(tmp_path: Path) -> None:
    tree = build_tree(tmp_path / "t", extra_values={"units.py": "U = 1\n"}, memo_edit=("traced_memo.py", "pinned[real] = path.read_bytes()", "pass"))
    memo = tree / "orpheus" / "memo" / "traced_memo.py"
    memo.write_text(memo.read_text(encoding="utf-8") + "\n\ndef _planted_late_import() -> None:\n    from orpheus.values import units  # noqa: F401\n")
    with pytest.raises(RuntimeError, match="KeyError"):
        generating_child(tree, _LATE, tmp_path / "w")


def test_the_gate_cannot_run_on_an_empty_tree(tmp_path: Path) -> None:
    """The instrument refuses a vacuous population (runs at base through ``Uninstallable`` only when ``orpheus/memo`` is absent)."""
    if any((_memo_tree.ORPHEUS / "memo").glob("*.py")):
        pytest.skip("orpheus/memo exists: the refusal is not reachable")
    with pytest.raises(_memo_tree.Uninstallable, match="empty"):
        build_tree(tmp_path / "t")
