"""The cone store's recorder: the rows one capture writes, against sets read off the fixture source by hand.

Spec ``scratch/verification_tree/cone_store/spec_b1.md`` §3.1, rows R1-R9. The expected sets are the
spec's independent ground for the recorder (§5): they are written from the fixture's source
(``_cone_store_fixture.py``), never produced by any recorder. R7-R9 are AST censuses of
``tools/cone_store/``: one definition per quantity (X4), one recorder (ruling 15), one fork guard (ruling 38).
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

import tools.cone_store as cone  # noqa: F401  (red at collection until B1 lands, never vacuous)
from orpheus.numerics import traced_memo

from ._cone_store_fixture import REPO, Captured, capture_base, code_keys, pin_key

pytestmark = pytest.mark.foundation

_TOOL = REPO / "tools" / "cone_store"


@pytest.fixture(scope="module")
def base(tmp_path_factory) -> Captured:
    return capture_base(tmp_path_factory.mktemp("recorder") / "tree")


def _defs(*pairs: tuple[str, str]) -> set[tuple]:
    return {("DefPin", path, name, 0) for path, name in pairs}


def _modules(*paths: str) -> set[tuple]:
    return {("ModulePin", path, "", -1) for path in paths}


def test_r1_the_rows_of_test_b_are_the_hand_written_set(base: Captured) -> None:
    """R1. test_b's def and module rows, exactly: its own def, the three defs of ``pkg.user`` it enters (the
    cache is empty in its module's parent: repair (1)), the conftest helper the grandparent ran, and the
    module code of every first-party module its process imported. ``[M]`` first red: the single shared
    parent of ruling 38 as first written loses five of these (spec §0.1)."""
    expected = _defs(("tests/test_b.py", "test_b"), ("pkg/user.py", "_build"), ("pkg/user.py", "table"),
                     ("pkg/user.py", "above"), ("pkg/conf_dep.py", "conftest_helper"))
    expected |= _modules("pkg/__init__.py", "pkg/consts.py", "pkg/user.py", "pkg/conf_dep.py", "tests/__init__.py",
                         "tests/conftest.py", "tests/test_b.py")
    got = code_keys(base.store.pins("tests/test_b.py::test_b"))
    assert got == expected, f"missing {sorted(expected - got)}; extra {sorted(got - expected)}"


def test_r2_a_tests_helper_is_attributed_to_its_caller_only(base: Captured) -> None:
    """R2. ``tests/_helpers.py:helper`` is in test_d's rows (the ``tests/`` root) and not in test_plain's."""
    helper = ("DefPin", "tests/_helpers.py", "helper", 0)
    assert helper in code_keys(base.store.pins("tests/test_d.py::test_d"))
    assert helper not in code_keys(base.store.pins("tests/test_d.py::test_plain"))


def test_r3_a_data_file_is_a_row_and_a_tmp_file_is_not(base: Captured) -> None:
    """R3. test_reads_data holds a DataPin of ``tests/data/ref.txt`` by its bytes; test_tmp is captured as
    passed and pins nothing under its tmp path (the memo's refusal of a path observed then changed is a
    reason the store drops for test-owned paths, not a lost manifest)."""
    data = [p for p in base.store.pins("tests/test_h.py::test_reads_data") if type(p).__name__ == "DataPin"]
    ref = base.root / "tests/data/ref.txt"
    hits = [p for p in data if (base.root / p.path).resolve() == ref.resolve() or Path(p.path).resolve() == ref.resolve()]
    assert len(hits) == 1, f"DataPins: {data}"
    assert hits[0].digest == traced_memo._file_digest(str(ref))
    record = base.store.tests()["tests/test_h.py::test_tmp"]
    assert record.outcome == "passed"
    tmp_pins = [p for p in base.store.pins("tests/test_h.py::test_tmp")
                if type(p).__name__ in {"DataPin", "PresencePin", "ListingPin"} and "pytest-of" in str(p.path)]
    assert tmp_pins == []


def test_r4_a_warm_reader_holds_the_child_and_not_the_oracle(base: Captured) -> None:
    """R4. On the warm capture the reader holds ChildPin("pkg.oracle:oracle", ...) and no def of
    ``pkg/oracle.py``; the validation test (``bypass()``, in process) holds ``_oracle_body``."""
    reader = base.store.pins("tests/test_e.py::test_reader")
    children = [p for p in reader if type(p).__name__ == "ChildPin"]
    assert [c.function_id for c in children] == ["pkg.oracle:oracle"]
    assert ("DefPin", "pkg/oracle.py", "_oracle_body", 0) not in code_keys(reader)
    validation = code_keys(base.store.pins("tests/test_e.py::test_oracle_validation"))
    assert ("DefPin", "pkg/oracle.py", "_oracle_body", 0) in validation


def test_r5_twin_defs_are_told_apart_by_ordinal_and_a_lambda_by_its_def(base: Captured) -> None:
    """R5. ``pick`` of ordinal 1 (the ``else`` branch) ran, ordinal 0 did not; the lambda inside ``outer`` is
    pinned as ``outer``. The live def is the SECOND homonym so a key that drops the ordinal (every def read as
    ordinal 0) is seen (battery 2026-11-01: with the live def first, ORD was designed-green)."""
    keys = code_keys(base.store.pins("tests/test_f.py::test_twins"))
    assert ("DefPin", "pkg/twins.py", "pick", 1) in keys
    assert ("DefPin", "pkg/twins.py", "pick", 0) not in keys
    assert ("DefPin", "pkg/twins.py", "outer", 0) in keys


def test_r6_every_digest_is_the_memos(base: Captured) -> None:
    """R6. Every DefPin digest is ``function_digest(source, qualname, ordinal)`` and every ModulePin digest is
    ``skeleton_digest(source)`` of the file at capture; over every row of every test."""
    checked = 0
    for test in base.store.tests():
        for pin in base.store.pins(test):
            kind = type(pin).__name__
            if kind not in {"DefPin", "ModulePin"}:
                continue
            source = (base.root / pin.path).read_bytes()
            want = (traced_memo.function_digest(source, pin.qualname, pin.ordinal) if kind == "DefPin"
                    else traced_memo.skeleton_digest(source))
            assert pin.digest == want, f"{test}: {pin_key(pin)}"
            checked += 1
    assert checked > 50, f"only {checked} code rows checked"


def _tool_modules() -> list[tuple[Path, ast.Module]]:
    modules = [(path, ast.parse(path.read_text(), str(path))) for path in sorted(_TOOL.rglob("*.py"))]
    assert modules, f"no module under {_TOOL}"
    return modules


def _calls(tree: ast.AST) -> set[str]:
    """Dotted names of every call and attribute reference in ``tree``."""
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute):
            parts, cursor = [node.attr], node.value
            while isinstance(cursor, ast.Attribute):
                parts.append(cursor.attr)
                cursor = cursor.value
            if isinstance(cursor, ast.Name):
                parts.append(cursor.id)
            names.add(".".join(reversed(parts)))
        elif isinstance(node, ast.Name):
            names.add(node.id)
    return names


def _leaf_names(tree: ast.AST) -> set[str]:
    """The last segment of every dotted reference and import in ``tree`` (``self.boot.fork_blockers`` gives
    ``fork_blockers``), beside the full spellings."""
    full = _calls(tree) | _imported(tree)
    return full | {n.rsplit(".", 1)[-1] for n in full}


def _imported(tree: ast.AST) -> set[str]:
    found: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            found |= {f"{node.module}.{alias.name}" for alias in node.names}
        elif isinstance(node, ast.Import):
            found |= {alias.name for alias in node.names}
    return found


def test_r7_the_store_defines_no_digest_of_its_own() -> None:
    """R7. One definition per quantity (X4): ``tools/cone_store/`` computes no digest of a source. No
    ``hashlib``, ``ast.parse`` or ``ast.dump``, no def whose name holds ``digest``, and every ``digest`` name it
    uses is imported from the memo (``orpheus.numerics._traced_memo_boot`` or ``traced_memo``). Its def and
    skeleton digests come with the memo's pins (``collect_pins``), which some module of the package names.
    Census: every module of the package, the count asserted non-zero, and the known member naming
    ``collect_pins`` found. (Re-pointed in review round 1: the store no longer imports ``function_digest``; it
    reads the digests off the pins.)"""
    modules = _tool_modules()
    banned = {"hashlib", "ast.parse", "ast.dump", "blake2b", "sha256"}
    memo_modules = ("orpheus.numerics._traced_memo_boot.", "orpheus.numerics.traced_memo.")
    for path, tree in modules:
        where = path.relative_to(REPO)
        hit = (_calls(tree) | _imported(tree)) & banned
        assert not hit, f"{where} spells a digest of its own: {sorted(hit)}"
        own = [n.name for n in ast.walk(tree) if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and "digest" in n.name]
        assert not own, f"{where} defines {own}"
        foreign = [n for n in _imported(tree) if n.rsplit(".", 1)[-1].endswith("digest") and not n.startswith(memo_modules)]
        assert not foreign, f"{where} imports a digest from outside the memo: {foreign}"
    readers = [path for path, tree in modules if "collect_pins" in _leaf_names(tree)]
    assert readers, f"no module of {len(modules)} reads the memo's pins through collect_pins"


def test_r8_the_store_starts_no_second_recorder() -> None:
    """R8. No ``sys.monitoring`` tool of its own and no audit hook in ``tools/cone_store/``; it constructs the
    memo's ``Recording`` (ruling 15: one change-detection layer)."""
    modules = _tool_modules()
    banned = {"sys.monitoring.use_tool_id", "sys.monitoring.register_callback", "monitoring.use_tool_id",
              "monitoring.register_callback", "sys.addaudithook"}
    for path, tree in modules:
        hit = _calls(tree) & banned
        assert not hit, f"{path.relative_to(REPO)} starts a recorder of its own: {sorted(hit)}"
    assert any("Recording" in _calls(tree) for _, tree in modules), f"no module of {len(modules)} names Recording"


def test_r9_the_fork_guard_is_the_memos() -> None:
    """R9. ``tools/cone_store/`` asks the memo's one fork guard, ``fork_blockers``, and assembles none of its own:
    it names neither ``FORK_MEASURED`` nor ``_fork_unsafe`` (the guard's parts) and counts no threads (no
    ``ctypes``, ``proc_pidinfo``, ``ps``, ``threading.active_count``). (Re-pointed in review round 1: the guard
    became ``_traced_memo_boot.fork_blockers``.)"""
    modules = _tool_modules()
    names = set().union(*(_leaf_names(tree) for _, tree in modules))
    assert "fork_blockers" in names, "the memo's guard is not asked"
    banned = {"ctypes", "proc_pidinfo", "threading.active_count", "active_count", "FORK_MEASURED", "_fork_unsafe"}
    for path, tree in modules:
        flat = _calls(tree) | _imported(tree) | {n.rsplit(".", 1)[-1] for n in _imported(tree)} | \
            {n.rsplit(".", 1)[-1] for n in _calls(tree)}
        hit = flat & banned
        assert not hit, f"{path.relative_to(REPO)} assembles a guard of its own: {sorted(hit)}"
        literals = {n.value for n in ast.walk(tree) if isinstance(n, ast.Constant) and isinstance(n.value, str)}
        assert "ps" not in literals, f"{path.relative_to(REPO)} runs ps"
