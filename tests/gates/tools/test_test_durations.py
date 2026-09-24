"""The instrument behind the #405 duration tables must read pytest's JUnit output correctly.

``tools/test_durations.py`` turns the JUnit XML of a sharded CI run into a
per-test duration table.  A table that mislabels a test, drops one, or reads a
timeout as a pass would send the campaign after the wrong tests, so each
reading is checked against an independent instrument: the node ids against
pytest's own ``--collect-only``, the outcomes against a synthetic file whose
every outcome is known, and the timeout against a test that sleeps past it.
"""
from __future__ import annotations

import pathlib
import subprocess
import sys

import pytest

from tools.test_durations import deal, gate_files, read_junit, summary_markdown, tree_of

pytestmark = pytest.mark.foundation

SYNTHETIC = '''
import time, pytest
def test_plain(): pass
class TestNested:
    @pytest.mark.parametrize("x", [1, 2])
    def test_param(self, x): pass
def test_fail(): assert False
def test_skip(): pytest.skip("not here")
def test_sleeps_past_the_timeout(): time.sleep(3)
'''

EXPECTED_OUTCOMES = {
    "test_plain": "passed",
    "TestNested::test_param[1]": "passed",
    "TestNested::test_param[2]": "passed",
    "test_fail": "failed",
    "test_skip": "skipped",
    "test_sleeps_past_the_timeout": "timeout",
}


@pytest.fixture(scope="module")
def synthetic_run(tmp_path_factory) -> tuple[pathlib.Path, pathlib.Path, set[str]]:
    """Run the synthetic file with a 1 s per-test timeout; return (root, JUnit XML, collected ids)."""
    root = tmp_path_factory.mktemp("durations")
    (root / "tests" / "gates" / "demo").mkdir(parents=True)
    (root / "tests" / "gates" / "demo" / "test_demo.py").write_text(SYNTHETIC)
    junit = root / "out.xml"
    base = [sys.executable, "-O", "-m", "pytest", "-q", "-p", "no:cacheprovider"]
    subprocess.run([*base, "--timeout=1", f"--junitxml={junit}", "tests/gates/demo"], cwd=root, capture_output=True)
    collected = subprocess.run([*base, "--collect-only", "tests/gates/demo"], cwd=root, capture_output=True, text=True)
    ids = {line for line in collected.stdout.splitlines() if "::" in line}
    return root, junit, ids


def test_node_ids_are_the_ones_pytest_collects(synthetic_run) -> None:
    root, junit, collected = synthetic_run
    read = {case.node_id for case in read_junit(junit, root)}
    if read != collected or len(collected) != len(EXPECTED_OUTCOMES):
        pytest.fail(f"read {sorted(read)}\ncollected {sorted(collected)}")


def test_every_outcome_is_read_including_a_timeout(synthetic_run) -> None:
    root, junit, _ = synthetic_run
    prefix = "tests/gates/demo/test_demo.py::"
    got = {case.node_id.removeprefix(prefix): case.outcome for case in read_junit(junit, root)}
    if got != EXPECTED_OUTCOMES:
        pytest.fail(f"outcomes {got}, expected {EXPECTED_OUTCOMES}")


def test_the_timed_out_test_carries_its_time(synthetic_run) -> None:
    root, junit, _ = synthetic_run
    timed_out = [c for c in read_junit(junit, root) if c.outcome == "timeout"]
    if len(timed_out) != 1 or timed_out[0].seconds < 1.0:
        pytest.fail(f"expected one timeout of at least 1 s, got {timed_out}")


def test_a_classname_with_no_module_file_is_refused(tmp_path) -> None:
    junit = tmp_path / "orphan.xml"
    junit.write_text('<testsuites><testsuite><testcase classname="tests.gates.gone.test_x" name="t" time="1"/>'
                     "</testsuite></testsuites>")
    with pytest.raises(ValueError, match="no test module"):
        read_junit(junit, tmp_path)


def test_dealing_partitions_the_gate_files() -> None:
    files = gate_files()
    if not files or "tests/gates/test_layer_imports.py" not in files:
        pytest.fail(f"the gate-file listing lost a known member ({len(files)} files)")
    groups = deal(files, 40)
    dealt = [f for g in groups for f in g]
    if sorted(dealt) != files or len(groups) != 40 or min(map(len, groups)) == 0:
        pytest.fail(f"{len(groups)} groups over {len(files)} files are not a partition into non-empty shards")


def test_the_summary_lists_the_slowest_first(synthetic_run) -> None:
    root, junit, _ = synthetic_run
    text = summary_markdown(read_junit(junit, root), top=2)
    rows = [line for line in text.splitlines() if line.startswith("| ") and "`" in line]
    if not rows or "test_sleeps_past_the_timeout" not in rows[0]:
        pytest.fail(f"the slowest case does not head the list:\n{text}")


@pytest.mark.parametrize(
    ("node_id", "tree"),
    [
        ("tests/gates/sn/sweep/test_x.py::test_a", "sn"),
        ("tests/gates/test_layer_imports.py::test_no_forbidden_imports[orpheus/sn/solver.py]", "(root)"),
        ("tests/gates/test_docstring_xrefs.py::TestK::test_path[tests/gates/sn/x.py]", "(root)"),
        ("tests/gates/derivations/test_y.py::TestK::test_b[1.0-a/b]", "derivations"),
    ],
)
def test_the_tree_is_read_from_the_path_not_the_parametrize_id(node_id: str, tree: str) -> None:
    """A parametrize id may contain ``/``; the tree comes from the file path alone.

    `[M]` 2026-09-23, the first CI run: splitting the whole node id on ``/``
    put three root-level tests in trees named after their parametrize ids.
    """
    if tree_of(node_id) != tree:
        pytest.fail(f"tree_of({node_id!r}) = {tree_of(node_id)!r}, expected {tree!r}")
