"""The cone store's ``select``: which tests a working-tree change selects, and positive controls 1-3.

Spec ``scratch/verification_tree/cone_store/spec_b1.md`` §3.3, rows S1-S10 and P1-P3. Every expected set
is written from the fixture's source. ``select`` runs in a fresh process (``_cone_store_fixture.select``).
"""
from __future__ import annotations

import os
import time
from pathlib import Path

import pytest

import tools.cone_store as cone

from ._cone_store_fixture import (BASE_IDS, BASELINE, MEMO_IDS, Captured, build, capture, capture_base, copy_captured,
                                  edit, ids_of, paths, select, store_of)

pytestmark = pytest.mark.foundation

ALL = BASE_IDS | MEMO_IDS


@pytest.fixture(scope="module")
def base(tmp_path_factory) -> Captured:
    return capture_base(tmp_path_factory.mktemp("selection") / "tree")


@pytest.fixture
def tree(base: Captured, tmp_path: Path) -> Path:
    return copy_captured(base.root, tmp_path / "tree")


def test_s1_no_edit_selects_only_the_unknown_and_the_not_passed(tree: Path) -> None:
    """S1. With no edit, ``select`` lists exactly the tests that did not pass at capture (``test_fail``); no
    base test spawns. The negative leg of every over-selection mutant."""
    assert ids_of(tree) == ALL, "the fixture's collected ids are not the hand-written set"
    selected, _ = select(tree)
    assert selected == BASELINE, f"extra {sorted(selected - BASELINE)}"


def test_s2_a_new_test_is_selected_alone(tree: Path) -> None:
    """S2. A test added after capture is selected (it is new), and no other test is: the new def binds a name no
    captured run read. Re-posed at B2 (#689), INVERTED: under the module skeleton the module's other tests were
    selected too (the skeleton moved); a statement pin sees only the statements a run observes."""
    with (tree / "tests/test_d.py").open("a") as handle:
        handle.write("\n\ndef test_new():\n    assert True\n")
    selected, _ = select(tree)
    assert selected == BASELINE | {"tests/test_d.py::test_new"}, f"got {sorted(selected)}"


def test_s3_a_deleted_test_is_not_listed_and_its_rows_are_pruned(tree: Path) -> None:
    """S3. ``tests/test_f.py`` deleted: ``select`` does not list ``test_twins``, and the next capture prunes it."""
    (tree / "tests/test_f.py").unlink()
    selected, _ = select(tree)
    assert "tests/test_f.py::test_twins" not in selected
    capture(tree)
    assert "tests/test_f.py::test_twins" not in store_of(tree).tests()


def test_s4_a_test_that_failed_at_capture_is_always_selected(tree: Path) -> None:
    """S4. ``test_fail`` failed at capture (its cone stops at the failure) and is selected with no edit."""
    assert store_of(tree).tests()["tests/test_g.py::test_fail"].outcome == "failed"
    selected, _ = select(tree)
    assert "tests/test_g.py::test_fail" in selected


def test_s5_another_environment_selects_everything(tree: Path) -> None:
    """S5 (ruling 34). The capture ran with ``OMP_NUM_THREADS=1``; ``select`` with ``OMP_NUM_THREADS=2`` finds no
    row under its key and lists every test."""
    selected, _ = select(tree, OMP_NUM_THREADS="2")
    assert selected == ALL, f"missing {sorted(ALL - selected)}"


@pytest.mark.parametrize("component", ["interpreter", "blas", "distributions", "threads"])
def test_s6_each_component_moves_the_environment_key(component: str, monkeypatch: pytest.MonkeyPatch) -> None:
    """S6. ``environment_key()`` changes when one of its four components changes, each varied at its source
    (spec §1 surface): ``traced_memo._python_identity``, ``traced_memo._blas``, ``importlib.metadata.distributions``,
    ``OMP_NUM_THREADS``. One id per component: a four-part key is four claims."""
    import importlib.metadata

    from orpheus.numerics import traced_memo

    monkeypatch.setenv("OMP_NUM_THREADS", "1")
    before = cone.environment_key()
    if component == "interpreter":
        monkeypatch.setattr(traced_memo, "_python_identity", lambda: "3.99.0|cpython-399")
    elif component == "blas":
        monkeypatch.setattr(traced_memo, "_blas", lambda: "openblas")
    elif component == "distributions":
        real = list(importlib.metadata.distributions())
        monkeypatch.setattr(importlib.metadata, "distributions", lambda **kw: real[:-1])
    else:
        monkeypatch.setenv("OMP_NUM_THREADS", "2")
    assert cone.environment_key() != before


def test_s7_an_edit_to_a_def_nobody_ran_selects_nothing(tree: Path) -> None:
    """S7 (def grain). ``pkg/user.py:untouched`` ran in no test; editing its body selects only the baseline."""
    edit(tree, "pkg/user.py", "def untouched():\n    return 0", "def untouched():\n    return 1")
    selected, _ = select(tree)
    assert selected == BASELINE, f"extra {sorted(selected - BASELINE)}"


@pytest.mark.parametrize("what", ["def", "file"])
def test_s8_a_deleted_def_or_file_selects_its_readers(tree: Path, what: str) -> None:
    """S8. Deleting ``above`` selects ``test_b`` and ``test_fail``; deleting ``tests/_helpers.py`` selects ``test_d``."""
    if what == "def":
        edit(tree, "pkg/user.py", "def above(x):\n    return x > TOL\n", "")
        expected = {"tests/test_b.py::test_b", "tests/test_g.py::test_fail"}
    else:
        (tree / "tests/_helpers.py").unlink()
        expected = {"tests/test_d.py::test_d"}
    selected, _ = select(tree)
    assert expected <= selected, f"missing {sorted(expected - selected)}"
    if what == "file":
        assert "tests/test_d.py::test_plain" not in selected


def test_s9_a_data_edit_selects_its_reader_only(tree: Path) -> None:
    """S9. An edit to ``tests/data/ref.txt`` selects ``test_reads_data`` and nothing beyond the baseline."""
    (tree / "tests/data/ref.txt").write_text("0.25\n")
    selected, _ = select(tree)
    assert selected == BASELINE | {"tests/test_h.py::test_reads_data"}


def test_s10_an_edit_to_a_tests_own_body_selects_it_alone(tree: Path) -> None:
    """S10. An edit to ``test_b``'s own body selects ``test_b`` and nothing beyond the baseline."""
    edit(tree, "tests/test_b.py", "assert user.above(1.0)", "assert user.above(2.0)")
    selected, _ = select(tree)
    assert selected == BASELINE | {"tests/test_b.py::test_b"}


def test_p1_control_1_module_grain_selects_the_readers_of_the_constant(tree: Path) -> None:
    """P1 (positive control 1). Mutating ``TOL`` selects every test whose executed code reads it (on this tree every
    importer of ``pkg.user`` reads it, so module and statement grain agree here; the exact control is
    ``test_cone_store_statement_grain.py`` SG1, SG2, SG4): ``test_a[t0]`` (``_build`` at collection), ``test_b``, ``test_fail`` and ``test_nested`` (``above``),
    and on this tree no other."""
    edit(tree, "pkg/consts.py", "TOL = 1e-15", "TOL = 1e-14")
    readers = {"tests/test_a.py::test_a[t0]", "tests/test_b.py::test_b", "tests/test_g.py::test_fail",
               "tests/test_m.py::test_nested"}
    selected, _ = select(tree)
    assert selected == BASELINE | readers


def test_p2_control_2_an_oracle_edit_selects_its_validation_and_its_reader(tree: Path) -> None:
    """P2 (positive control 2). Mutating ``_oracle_body`` selects ``test_oracle_validation`` (its DefPin) and
    ``test_reader`` (its ChildPin: the entry's manifest no longer validates), and no other test. Route leg:
    ``select`` regenerated nothing (the memo store's files and their contents are unchanged)."""
    memo = paths(tree)["memo"]
    before = {p: p.read_bytes() for p in memo.rglob("*") if p.is_file()}
    assert before, "the memo store is empty: the reader never generated"
    edit(tree, "pkg/oracle.py", "return float(n) * 0.5", "return float(n) * 0.5 + 0.0")
    selected, _ = select(tree)
    assert selected == BASELINE | MEMO_IDS, f"got {sorted(selected)}"
    after = {p: p.read_bytes() for p in memo.rglob("*") if p.is_file()}
    assert after == before, "select regenerated the memo entry"


@pytest.mark.parametrize("kind", ["def-docstring", "module-docstring", "comment", "blank-line"])
def test_p3_control_3_a_prose_or_layout_edit_selects_nothing(tree: Path, kind: str) -> None:
    """P3 (positive control 3). A docstring on ``_build``, a module docstring on ``pkg/user.py``, a comment in
    ``above`` and a blank line between two defs each select nothing beyond the baseline."""
    if kind == "def-docstring":
        edit(tree, "pkg/user.py", "def _build(n):\n", 'def _build(n):\n    """The table."""\n')
    elif kind == "module-docstring":
        text = (tree / "pkg/user.py").read_text()
        (tree / "pkg/user.py").write_text('"""The user module."""\n' + text)
    elif kind == "comment":
        edit(tree, "pkg/user.py", "    return x > TOL", "    # strictly above\n    return x > TOL")
    else:
        edit(tree, "pkg/user.py", "\n\ndef above(x):", "\n\n\ndef above(x):")
    selected, _ = select(tree)
    assert selected == BASELINE, f"extra {sorted(selected - BASELINE)}"
    assert os.path.exists(tree / "pkg/user.py")


def test_s11_a_test_that_reads_a_source_as_text_is_selected_by_a_docstring_edit(tmp_path: Path) -> None:
    """S11 (added by the battery, 2026-11-01). ``test_reads_source`` imports ``pkg.user`` (so its module code ran
    and is pinned) and reads ``pkg/user.py`` as TEXT, asserting a docstring is absent. A docstring added to
    ``_build`` changes its verdict, so it must be selected. The store drops the DataPin of a ``.py`` file that is
    also a pinned module (to keep the def grain against pytest's assertion-rewriting read of each test file):
    that drop under-selects every census test reading the source it imported (the elegance-debt ledger, the
    layer linter, the docstring cross-reference gate). First red: today's build."""
    root = build(tmp_path / "tree", memo=False, reads_source=True)
    capture(root)
    edit(root, "pkg/user.py", "def _build(n):\n", 'def _build(n):\n    """The table."""\n')
    selected, _ = select(root)
    assert "tests/test_k.py::test_reads_source" in selected


# ── Review round 1: qa's under-selections (qa_b1.md §1, U3-U5, U7-U9), each repaired by the builder ──

PKG_USER_IMPORTERS = frozenset({"tests/test_a.py::test_a[t0]", "tests/test_b.py::test_b", "tests/test_g.py::test_fail",
                                "tests/test_m.py::test_nested"})


def test_s12_a_new_package_that_shadows_an_imported_module_selects_its_importers(tree: Path) -> None:
    """S12 (qa U5). ``pkg/consts/__init__.py`` created beside ``pkg/consts.py``: the package is found first, so
    every test whose process imported ``pkg.consts`` now reads another ``TOL``. Seen only through the presence
    pins of the names that would shadow each imported module (the import system's own probes are not recorded).
    First red: SHADOW-DROP."""
    (tree / "pkg/consts").mkdir()
    (tree / "pkg/consts/__init__.py").write_text("TOL = 1.0\n")
    selected, _ = select(tree)
    assert PKG_USER_IMPORTERS <= selected, f"missing {sorted(PKG_USER_IMPORTERS - selected)}"


def test_s13_a_new_conftest_in_a_test_directory_selects_every_test(tmp_path: Path) -> None:
    """S13 (qa U3; re-posed after qa's second pass). A ``conftest.py`` created in ``tests/sub/`` (an autouse
    fixture that raises) selects ``tests/sub/test_x.py::test_x``, and, because the store promises full-run
    verdicts and a full run loads every conftest session-wide, every other test too: the set of conftests is in
    the environment key. First red: PRESENCE-DROP and KEY-CONFTESTS-DROP together (either alone leaves the other
    route)."""
    root = build(tmp_path / "tree", memo=False, subdir=True)
    capture(root)
    (root / "tests/sub/conftest.py").write_text(
        "import pytest\n\n\n@pytest.fixture(autouse=True)\ndef boom():\n    raise RuntimeError('new conftest')\n")
    selected, _ = select(root)
    assert "tests/sub/test_x.py::test_x" in selected
    assert BASE_IDS <= selected, f"missing {sorted(BASE_IDS - selected)}"


def test_s14_a_change_to_the_pytest_settings_selects_everything(tree: Path) -> None:
    """S14 (qa U4). A ``filterwarnings`` entry added to ``pyproject.toml``'s pytest settings: pytest reads the ini
    before any recorder starts, so its bytes are part of the environment key and every test is selected.
    First red: ENV-INI."""
    edit(tree, "pyproject.toml", 'testpaths = ["tests"]', 'testpaths = ["tests"]\nfilterwarnings = ["error::DeprecationWarning"]')
    selected, _ = select(tree)
    assert ALL <= selected, f"missing {sorted(ALL - selected)}"


def test_s15_a_module_whose_collection_errors_is_listed_and_recorded(tree: Path) -> None:
    """S15 (qa U8). A syntax error in ``tests/test_d.py``: ``select`` lists the module's id, and a capture records
    it as a row with outcome ``collection-error`` instead of losing it. First red: COLLERR-DROP."""
    (tree / "tests/test_d.py").write_text("def test_d(:\n    pass\n")
    selected, _ = select(tree)
    assert "tests/test_d.py" in selected, sorted(selected)
    capture(tree)
    record = store_of(tree).tests().get("tests/test_d.py")
    assert record is not None and record.outcome == "collection-error", record


def test_s16_a_test_that_reads_and_rewrites_a_repository_file_is_cone_unknown(tmp_path: Path) -> None:
    """S16 (qa U7). ``test_rw`` reads ``tests/data/counter.txt`` and writes it back plus one: no state of that file
    pins its verdict, so the test is cone-unknown and selected with no edit (a ``tmp_path`` file, outside the
    repository, stays dropped: R3). First red: RW-DROP."""
    root = build(tmp_path / "tree", memo=False, read_rewrite=True)
    capture(root)
    assert store_of(root).tests()["tests/test_rw.py::test_rw"].unknown
    assert not store_of(root).tests()["tests/test_h.py::test_tmp"].unknown
    selected, _ = select(root)
    assert "tests/test_rw.py::test_rw" in selected


def test_s17_an_edit_made_while_the_capture_runs_makes_the_test_cone_unknown(tmp_path: Path) -> None:
    """S17 (qa U9). ``test_race`` runs ``above``, then waits; the gate flips ``>`` to ``>=`` in ``above`` (still
    true for 1.0) while it waits. The pins read the edited bytes after the test ran, so they would agree with
    the edited tree; the late write must instead make the cone unknown, and ``select`` lists the test.
    First red: LATE-DROP."""
    import threading

    root = build(tmp_path / "tree", memo=False, race=True)
    mark = tmp_path / "race_mark"
    failure: list[BaseException] = []

    def run() -> None:
        try:
            capture(root, "tests/test_race.py", CONE_RACE_MARK=str(mark))
        except BaseException as error:  # noqa: BLE001 - re-raised below in the gate's thread
            failure.append(error)

    worker = threading.Thread(target=run)
    worker.start()
    deadline = time.monotonic() + 120
    while not mark.exists() and time.monotonic() < deadline and worker.is_alive():
        time.sleep(0.05)
    assert mark.exists(), "precondition: test_race never ran"
    edit(root, "pkg/user.py", "return x > TOL", "return x >= TOL")
    Path(str(mark) + ".edited").write_text("edited")
    worker.join(180)
    assert not failure, failure
    record = store_of(root).tests()["tests/test_race.py::test_race"]
    assert record.outcome == "passed"
    assert any("changed during the capture" in reason for reason in record.unknown), record.unknown
    selected, _ = select(root)
    assert "tests/test_race.py::test_race" in selected


# ── qa's second pass (qa_b1_delta.md, D1-D5): the store promises FULL-RUN verdicts ──

def test_s18_an_off_path_conftests_session_hook_selects_a_test_elsewhere(tmp_path: Path) -> None:
    """S18 (D2). ``tests/zother/conftest.py``'s ``pytest_configure`` is edited to rebind ``tests._helpers.helper``.
    In a full run that conftest is loaded before any test, so ``tests/test_d.py::test_d`` (another directory)
    turns red; the grandparent loads every conftest a full run loads, so the edit selects it. First red:
    GP-CONFTESTS-DROP (the grandparent loads only pytest's initial conftests)."""
    root = build(tmp_path / "tree", memo=False, offpath=True)
    capture(root)
    edit(root, "tests/zother/conftest.py", "    pass\n",
         "    import tests._helpers as helpers\n\n    helpers.helper = lambda: 8\n")
    selected, _ = select(root)
    assert "tests/test_d.py::test_d" in selected, sorted(selected)


def test_s19_a_new_conftest_in_a_test_star_directory_selects_every_test(tmp_path: Path) -> None:
    """S19 (D1). ``tests/test_extra/conftest.py`` created (pytest loads ``<testpath>/test*/conftest.py`` as an
    initial conftest, session-wide): every test is selected through the conftest set in the environment key.
    First red: KEY-CONFTESTS-DROP."""
    root = build(tmp_path / "tree", memo=False)
    capture(root)
    (root / "tests/test_extra").mkdir()
    (root / "tests/test_extra/conftest.py").write_text("import pkg.user\n\npkg.user.TOL = 10.0\n")
    selected, _ = select(root)
    assert BASE_IDS <= selected, f"missing {sorted(BASE_IDS - selected)}"


def test_s20_a_user_ignore_collect_hook_that_fills_a_cache_stays_recorded(tmp_path: Path) -> None:
    """S20 (D3). ``tests/conftest.py::pytest_ignore_collect`` calls a ``functools.cache`` helper that reads
    ``tests/data/flag.txt``; ``test_flag`` calls it again (a hit in the forked child). Only pytest's own
    ``pytest_ignore_collect`` implementations run unrecorded, so the hook's read is a dependency and an edit to
    the file selects ``test_flag``. First red: IGNORE-ALL (every implementation silenced)."""
    root = build(tmp_path / "tree", memo=False, user_ignore=True)
    capture(root)
    (root / "tests/data/flag.txt").write_text("changed\n")
    selected, _ = select(root)
    assert "tests/test_flag.py::test_flag" in selected, sorted(selected)


def test_s21_a_pytest11_entry_point_under_an_unchanged_version_selects_everything(tmp_path: Path) -> None:
    """S21 (D4). The local distribution ``localplug 1.0`` gains ``[pytest11] localplug = localplug_mod`` with its
    name and version unchanged: pytest loads it before any recorder starts, so the entry points are in the
    environment key and every test is selected. First red: KEY-PYTEST11-DROP."""
    root = build(tmp_path / "tree", memo=False, localplug=True)
    capture(root)
    (root / "localplug-1.0.dist-info/entry_points.txt").write_text("[pytest11]\nlocalplug = localplug_mod\n")
    selected, _ = select(root)
    assert BASE_IDS <= selected, f"missing {sorted(BASE_IDS - selected)}"


def test_s22_pytest_plugins_in_the_environment_selects_everything(tmp_path: Path) -> None:
    """S22 (D5). ``PYTEST_PLUGINS=localplug_mod`` at select (absent at capture): every test is selected through the
    environment key. First red: KEY-PYTEST-PLUGINS-DROP."""
    root = build(tmp_path / "tree", memo=False, localplug=True)
    capture(root)
    selected, _ = select(root, PYTEST_PLUGINS="localplug_mod")
    assert BASE_IDS <= selected, f"missing {sorted(BASE_IDS - selected)}"


def test_s23_a_reader_of_its_own_stores_entry_is_selected_by_the_generated_code_only(tmp_path: Path) -> None:
    """S23. ``test_own_store`` reads ``own_oracle(4)`` from a store under its ``tmp_path``; the generation imports
    ``pkg.deep``, which the test's process never imports. With no edit it is not selected: its child is pinned by
    the entry's own pins, not by a ``ChildPin`` judged against the default store, where the entry is absent
    (``[M]`` 2026-11-03, B3 acceptance X1.3: 70 such readers selected with no edit). An edit to ``_deep_body``
    selects it (the generated code is held, not dropped with the child), resolved through ``sys.path`` since no
    code of ``pkg/deep.py`` ran in the test's process; that file's identity is recorded (qa F3: it was written as
    missing, disabling stage one for it). First reds: the ``ChildPin`` (no-edit leg), a cone that drops the pins of
    files its own process never ran (edit leg), and the missing state (qa F3)."""
    root = capture_base(tmp_path / "tree", own_store=True).root
    own = "tests/test_own.py::test_own_store"
    pins = store_of(root).pins(own)
    if children := [pin for pin in pins if type(pin).__name__ == "ChildPin"]:
        pytest.fail(f"the reader of its own store holds a ChildPin: {children}")
    deep = store_of(root).identities().get("pkg/deep.py")
    if deep is None or deep.ident is None or deep.digest is None:
        pytest.fail(f"pkg/deep.py, read only by the generation, is recorded without its identity: {deep}")
    selected, _ = select(root)
    assert selected == BASELINE, f"extra {sorted(selected - BASELINE)}"
    edit(root, "pkg/deep.py", "return float(n) * 0.25", "return float(n) * 0.25 + 0.0")
    selected, _ = select(root)
    assert selected == BASELINE | {own}, f"got {sorted(selected)}"
