"""The cone store at statement grain (increment B2, issue #689): positive control 1 exact, and recall at the store.

Spec: ``scratch/verification_tree/cone_store/b2/gate_spec.md`` §1 (the store rows) and §5 (control 1 on the production
tree). The memo-level rows, one per shape of the design's census, are
``tests/gates/numerics/test_traced_memo_statement_pins.py``; this file asks the same questions of ``select``.

Before B2 (B1, `e9df7ed9`) a test's cone pinned every first-party module whose code ran (at collection or in the test) by its skeleton, so an
edit to one constant selects every test whose process imported the module. After B2 it selects the tests whose
executed code reads the constant. The precision rows here are Red at B1 (`e9df7ed9`); the recall rows are green and stay green.

Expected sets are written from the fixture's sources by hand (control 1 on the fixture), or measured by an instrument
that is not the store (control 1 on the production tree: each test run alone in a fresh interpreter under a monitor
that records whether the reader's code object started).
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import textwrap
import time
from pathlib import Path

import pytest

from ._cone_store_fixture import REPO, _trusted_age_s, capture, edit, environment, ids_of, paths, select, store_of

pytestmark = pytest.mark.foundation

#: The statement-grain tree. ``pkg/consts.py`` binds two constants; ``test_alpha`` reads one, ``test_beta`` the other,
#: ``test_neutral`` imports the same module chain and reads neither. The other modules carry one recall shape each.
_TREE: dict[str, str] = {
    "pyproject.toml": """
        [tool.pytest.ini_options]
        addopts = "-p no:cacheprovider"
        testpaths = ["tests"]
    """,
    "pkg/__init__.py": "",
    "pkg/consts.py": """
        ALPHA_READ_TOL = 1e-15
        BETA_READ_TOL = 2e-15
    """,
    "pkg/user.py": """
        from pkg import consts


        def above_alpha(x):
            return x > consts.ALPHA_READ_TOL


        def below_beta(x):
            return x < consts.BETA_READ_TOL


        def neutral():
            return 0
    """,
    "pkg/registry.py": """
        GAMMA_TABLE = {}
    """,
    "pkg/fill.py": """
        from pkg.registry import GAMMA_TABLE

        GAMMA_TABLE.update(gamma_entry=3)
    """,
    "pkg/k.py": """
        DELTA_SOURCE_LIMIT = 4
    """,
    "pkg/facade.py": """
        from pkg.k import DELTA_SOURCE_LIMIT as DELTA_PUBLIC_LIMIT
    """,
    "pkg/sigs.py": """
        EPSILON_SCALE = 2


        def scaled(x, scale=EPSILON_SCALE):
            return x * scale
    """,
    "tests/__init__.py": "",
    "tests/test_s.py": """
        from pkg import user


        def test_alpha():
            assert user.above_alpha(1.0)


        def test_beta():
            assert user.below_beta(0.0)


        def test_neutral():
            assert user.neutral() == 0
    """,
    "tests/test_t.py": """
        from pkg import fill, registry


        def test_table():
            assert registry.GAMMA_TABLE["gamma_entry"] == 3
    """,
    "tests/test_r.py": """
        from pkg import facade


        def test_reexport():
            assert facade.DELTA_PUBLIC_LIMIT == 4
    """,
    "pkg/confmod.py": """
        CONF_LIMIT = 6
    """,
    "pkg/late.py": """
        LATE_CONST = 9
    """,
    "tests/conftest.py": """
        import pytest

        from pkg.confmod import CONF_LIMIT


        @pytest.fixture
        def conf_limit():
            return CONF_LIMIT
    """,
    "tests/test_q.py": """
        import pytest


        @pytest.fixture
        def late_value():
            from pkg import late

            return late.LATE_CONST


        def test_conf(conf_limit):
            assert conf_limit == 6


        def test_late(late_value):
            assert late_value == 9
    """,
    # Code under the repository outside the package and ``tests/`` (here ``devtools/``, a stand-in for ``tools/``,
    # which a fixture cannot name: the cone store itself is imported as ``tools``).
    "devtools/__init__.py": "",
    "devtools/helper_tool.py": """
        def tool_value():
            return 4
    """,
    # qa's final F1: an unread table naming a def; deleting the def makes the module raise at import.
    "pkg/reg.py": """
        def reg_g():
            return 2


        _REG_TABLE = [reg_g]
        REG_VALUE = 1
    """,
    "tests/test_reg.py": """
        def test_reg():
            from pkg import reg  # inside the test: a module that raises at import fails the test, not collection

            assert reg.REG_VALUE == 1
    """,
    "tests/test_tool.py": """
        from devtools.helper_tool import tool_value


        def test_tool():
            assert tool_value() == 4
    """,
    "tests/test_p.py": """
        import pytest

        QA_CASES_LIST = [1, 2]
        QA_MESHES = [4, 8]


        @pytest.mark.parametrize("x", QA_CASES_LIST, ids=["first", "second"])
        def test_param(x):
            assert x < 3


        @pytest.fixture(params=QA_MESHES, ids=["m0", "m1"])
        def qa_cells_zz(request):
            return request.param


        def test_cells(qa_cells_zz):
            if qa_cells_zz % 4:
                raise AssertionError
    """,
    "tests/test_d.py": """
        from pkg import sigs


        def test_default():
            assert sigs.scaled(3) == 6
    """,
}

ALPHA = "tests/test_s.py::test_alpha"
BETA = "tests/test_s.py::test_beta"
NEUTRAL = "tests/test_s.py::test_neutral"
TABLE = "tests/test_t.py::test_table"
REEXPORT = "tests/test_r.py::test_reexport"
DEFAULT = "tests/test_d.py::test_default"
PARAM = ("tests/test_p.py::test_param[first]", "tests/test_p.py::test_param[second]")
CELLS = ("tests/test_p.py::test_cells[m0]", "tests/test_p.py::test_cells[m1]")
CONF = "tests/test_q.py::test_conf"
LATE = "tests/test_q.py::test_late"
TOOL = "tests/test_tool.py::test_tool"
REG = "tests/test_reg.py::test_reg"
ALL = frozenset({ALPHA, BETA, NEUTRAL, TABLE, REEXPORT, DEFAULT, *PARAM, *CELLS, CONF, LATE, TOOL, REG})

_THIS = "tests/gates/tools/test_cone_store_statement_grain.py"
_SG0 = f"{_THIS}::test_sg0_no_edit_selects_nothing"
_GT = "tests/gates/tools/test_cone_store_ground_truth.py"


def _build(root: Path) -> Path:
    for relative, source in _TREE.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(textwrap.dedent(source).lstrip("\n"))
    time.sleep(_trusted_age_s())
    return root


@pytest.fixture(scope="module")
def captured(tmp_path_factory) -> Path:
    root = _build(tmp_path_factory.mktemp("statement_grain") / "tree")
    capture(root)
    return root


@pytest.fixture
def tree(captured: Path, tmp_path: Path) -> Path:
    from ._cone_store_fixture import copy_captured

    return copy_captured(captured, tmp_path / "tree")


@pytest.fixture(params=["rows", "rows-and-columns"])
def old_form(captured: Path, tmp_path: Path, request: pytest.FixtureRequest) -> Path:
    """A copy of the captured tree whose store is rewritten in the form B1 wrote: no ``format`` meta, a skeleton
    ``ModulePin`` row in place of every ``StatementsPin`` and, in the second form, a ``files`` table without the
    ``digest`` column (B1's schema)."""
    import sqlite3

    from ._cone_store_fixture import copy_captured

    root = copy_captured(captured, tmp_path / "old")
    with sqlite3.connect(paths(root)["store"]) as db:
        db.execute("DELETE FROM meta WHERE key = 'format'")
        rows = db.execute("SELECT pin_id, row FROM pins").fetchall()
        rewritten = 0
        for pin_id, row in rows:
            fields = json.loads(row)
            if fields[0] == "StatementsPin":
                db.execute("UPDATE pins SET row = ? WHERE pin_id = ?",  # rows are unique: keep each one's digest
                           (json.dumps(["ModulePin", fields[1], hashlib.sha256(row.encode()).hexdigest()]), pin_id))
                rewritten += 1
        if request.param == "rows-and-columns":
            db.execute("CREATE TABLE files_b1(path TEXT PRIMARY KEY, ident TEXT, settled INTEGER NOT NULL)")
            db.execute("INSERT INTO files_b1 SELECT path, ident, settled FROM files")
            db.execute("DROP VIEW IF EXISTS rows")
            db.execute("DROP TABLE files")
            db.execute("ALTER TABLE files_b1 RENAME TO files")
            columns = [row[1] for row in db.execute("PRAGMA table_info(files)")]
            if "digest" in columns:
                pytest.fail("the B1 files table still has a digest column: the fixture is vacuous")
    if not rewritten:
        pytest.fail("the captured store holds no StatementsPin row to rewrite: the fixture is vacuous")
    return root


@pytest.mark.rests_on(f"{_THIS}::test_sg0_no_edit_selects_nothing")
def test_sg5_a_store_in_the_old_form_is_emptied_and_every_test_is_new(old_form: Path) -> None:
    """SG5. A store written before B2 (no ``format`` meta, ``ModulePin`` rows the memo no longer decodes) opens
    without raising, and ``select`` lists every collected test as new. Reddening mutation: drop the format check in
    ``ConeStore.__init__`` (``KeyError: 'ModulePin'`` in the pin decoder)."""
    selected, _ = select(old_form)
    assert selected == ALL, f"missing {sorted(ALL - selected)}"
    # The store is usable again: a capture writes it (each file with its digest), and a select then lists nothing.
    # Reddening mutation for the column form: empty the old tables instead of dropping them (the B1 ``files`` table
    # has no ``digest`` column, so the capture's write fails).
    capture(old_form)
    selected, _ = select(old_form)
    assert selected == frozenset(), f"after a capture into the reset store: {sorted(selected)}"


@pytest.mark.rests_on(f"{_THIS}::test_sg0_no_edit_selects_nothing")
def test_sg6_an_import_that_starts_raising_stays_selected_after_a_select(tree: Path) -> None:
    """SG6 (qa's final F1). Deleting ``reg_g``, which only the unread ``_REG_TABLE`` names, makes ``pkg/reg.py`` raise
    at import: only the import outcome sees it. ``test_reg`` is selected, and selected AGAIN by a second ``select`` with
    nothing changed between: the first select refreshes the file's stored identity (its statement pin was clean), and
    the refresh must keep the digest of the bytes it judged, never the capture's. Reddening mutation: the refresh
    stores the capture's digest (the second select vouches for the old bytes and reads the import set fresh)."""
    edit(tree, "pkg/reg.py", "def reg_g():\n    return 2\n", "")
    # Settle the edit before the selects: a file written within ``TRUSTED_AGE_NS`` is refreshed as unsettled, and an
    # unsettled identity never passes stage one, so the second select would read the file whatever the refresh kept
    # (``[M]`` 2026-11-02: without the wait the reddening mutation stayed green).
    time.sleep(_trusted_age_s())
    first, _ = select(tree)
    second, _ = select(tree)
    assert REG in first, f"the first select missed the import that now raises: {sorted(first)}"
    assert REG in second, f"the second select read the broken import set as fresh: {sorted(second)}"


def test_sg0_no_edit_selects_nothing(tree: Path) -> None:
    """SG0 (the negative leg of every precision row). With no edit, ``select`` lists nothing: every test passed at
    capture and none spawns. Green at B1 and after B2. Population leg: pytest collects exactly the hand-written ids
    (X2)."""
    assert ids_of(tree) == ALL, f"collected ids differ from the hand-written set: {sorted(ids_of(tree) ^ ALL)}"
    selected, _ = select(tree)
    assert selected == frozenset(), f"extra {sorted(selected)}"


@pytest.mark.rests_on(_SG0)
def test_sg1_control_1_exact_an_edit_to_a_read_constant_selects_its_readers_only(tree: Path) -> None:
    """SG1 (positive control 1 exact, on the fixture). ``ALPHA_READ_TOL`` is read by ``above_alpha`` only: editing it
    selects ``test_alpha`` and no other test, though ``test_beta`` and ``test_neutral`` import the same module.
    Red at B1 (`e9df7ed9`): the skeleton of ``pkg/consts.py`` moves, and all three tests of ``test_s.py`` are selected."""
    edit(tree, "pkg/consts.py", "ALPHA_READ_TOL = 1e-15", "ALPHA_READ_TOL = 1e-14")
    selected, _ = select(tree)
    assert selected == {ALPHA}, f"extra {sorted(selected - {ALPHA})}; missing {sorted({ALPHA} - selected)}"


@pytest.mark.rests_on(_SG0)
def test_sg2_control_1_exact_the_sibling_constant(tree: Path) -> None:
    """SG2. Editing ``BETA_READ_TOL`` selects ``test_beta`` alone. Red at B1 (`e9df7ed9`), as SG1."""
    edit(tree, "pkg/consts.py", "BETA_READ_TOL = 2e-15", "BETA_READ_TOL = 3e-15")
    selected, _ = select(tree)
    assert selected == {BETA}, f"extra {sorted(selected - {BETA})}; missing {sorted({BETA} - selected)}"


@pytest.mark.rests_on(_SG0)
@pytest.mark.parametrize(("file", "old", "new", "readers"), [
    pytest.param("pkg/fill.py", "GAMMA_TABLE.update(gamma_entry=3)", "GAMMA_TABLE.update(gamma_entry=5)", {TABLE},
                 id="effect-statement"),
    pytest.param("pkg/k.py", "DELTA_SOURCE_LIMIT = 4", "DELTA_SOURCE_LIMIT = 5", {REEXPORT}, id="reexport"),
    pytest.param("pkg/sigs.py", "EPSILON_SCALE = 2", "EPSILON_SCALE = 3", {DEFAULT}, id="default-argument"),
    # qa's S2a and S2b'' (ruling R1): a test module's constant read only by a decorator; ids held fixed, so the
    # test is not "new" and only an import-time pin of the test module can select it.
    pytest.param("tests/test_p.py", "QA_CASES_LIST = [1, 2]", "QA_CASES_LIST = [1, 5]", {PARAM[1]},
                 id="qa-S2a-parametrize-constant"),
    pytest.param("tests/test_p.py", "QA_MESHES = [4, 8]", "QA_MESHES = [4, 9]", {CELLS[1]},
                 id="qa-S2b-fixture-params-constant"),
    # The quiet sections of the capture (``unrecorded``: pytest's own collection plumbing) must not swallow
    # first-party code: a conftest's module code (its import, run while the conftest is loaded) and an import made
    # inside a fixture at setup, each read by one test.
    pytest.param("pkg/confmod.py", "CONF_LIMIT = 6", "CONF_LIMIT = 7", {CONF}, id="conftest-module-code"),
    pytest.param("pkg/late.py", "LATE_CONST = 9", "LATE_CONST = 8", {LATE}, id="fixture-time-import"),
    # First-party is every source under the repository (outside the virtual environment and the bytecode cache), not
    # only ``<package>/`` and ``tests/``: a body edit to a def under ``devtools/`` that a test ran selects it.
    # Reddening mutation: ``Reading.first_party`` back to the two roots (the def is dropped from every cone).
    pytest.param("devtools/helper_tool.py", "    return 4", "    return 2 + 2", {TOOL}, id="code-outside-package-and-tests"),
])
def test_sg3_recall_at_the_store(tree: Path, file: str, old: str, new: str, readers: set[str]) -> None:
    """SG3 (recall at the store; green at B1 and after B2). An import-time edit the test observes selects it: a bare
    expression statement filling a registry, a constant read through a renamed re-export, a default argument. The
    edit is one only an import-time pin sees (no def that ran changes). Today the selection is exactly the readers,
    since each module is imported by one test file only; the row asserts containment, not equality, so the effect
    statement's module grain (a declared over-selection) cannot redden it."""
    edit(tree, file, old, new)
    selected, _ = select(tree)
    assert readers <= selected, f"missing {sorted(readers - selected)}"


# ── control 1 on the production tree (slow) ────────────────────────────────────────────────────────────────────

_AXIS = "orpheus/mesh/axis.py"
_READER = "face_outflow_ordinates"
_MODULE = "tests/gates/sn/primitives/test_axis_primitive.py"
#: The sample: three tests that call the reader directly, six that import ``axis.py`` and (by the monitor's reading,
#: asserted below) may or may not reach it. Every one imports ``axis.py``, so at B1 an edit to its constant selects all.
_SAMPLE = tuple(f"{_MODULE}::{name}" for name in (
    "test_f0_1_axismesh_endpoints_min_max",
    "test_f0_1_axismesh_monotonicity_validated",
    "test_f0_5_face_shape_2d_each_face",
    "test_f0_6_face_outflow_ordinates_matches_inline_expression_1d",
    "test_f0_6_face_outflow_ordinates_2d_LS4",
    "test_f0_6_face_outflow_ordinates_radial",
    "test_f0_7_synthetic_3d_face_outflow_axis_beyond_quad_dim_is_empty",
    "test_f0_8_n_unknowns_flat_slab",
    "test_f0_8_n_unknowns_flat_2d_cartesian",
))

_MONITOR = r'''
"""Records whether orpheus/mesh/axis.py:face_outflow_ordinates started, from this plugin's import (before the test
module is collected) to the end of the session. Not the store's recorder: a tool id of its own, one code object."""
import os, sys
from pathlib import Path

_TOOL = 5
_HIT = []


def _start(code, offset):
    if code.co_qualname == "face_outflow_ordinates" and code.co_filename.endswith(os.path.join("orpheus", "mesh", "axis.py")):
        _HIT.append(code.co_filename)
    return sys.monitoring.DISABLE


sys.monitoring.use_tool_id(_TOOL, "gt-axis")
sys.monitoring.register_callback(_TOOL, sys.monitoring.events.PY_START, _start)
sys.monitoring.set_events(_TOOL, sys.monitoring.events.PY_START)


def pytest_unconfigure(config):
    Path(os.environ["GT_AXIS_OUT"]).write_text("1" if _HIT else "0")
'''


def _clone(target: Path) -> Path:
    """A copy-on-write clone (APFS ``cp -c``) of what a run of the sample reads: the package, the tests, the tools and
    the configuration. The production files are never edited; the clone is."""
    target.mkdir(parents=True)
    for name in ("orpheus", "tests", "tools"):
        subprocess.run(["cp", "-c", "-R", str(REPO / name), str(target / name)], check=True)
    for name in ("pyproject.toml",):
        shutil.copy2(REPO / name, target / name)
    for cache in target.rglob("__pycache__"):
        shutil.rmtree(cache, ignore_errors=True)
    time.sleep(_trusted_age_s())
    return target


def _readers_by_monitor(root: Path, plugin_dir: Path) -> frozenset[str]:
    """The tests of the sample whose run (collection included) started the reader, each run alone in a fresh
    interpreter."""
    found: set[str] = set()
    for test_id in _SAMPLE:
        out = plugin_dir / "hit.txt"
        out.unlink(missing_ok=True)
        env = environment(root, GT_AXIS_OUT=str(out))
        env["PYTHONPATH"] = os.pathsep.join([str(plugin_dir), str(root)])
        run = subprocess.run([sys.executable, "-O", "-m", "pytest", "-q", "-p", "gt_axis_monitor", "-p", "no:cacheprovider",
                              test_id], cwd=root, env=env, capture_output=True, text=True, timeout=600)
        if run.returncode != 0 or not out.exists():
            pytest.fail(f"the monitor's run of {test_id} exited {run.returncode}:\n{run.stdout[-2000:]}\n{run.stderr[-2000:]}")
        if out.read_text() == "1":
            found.add(test_id)
    return frozenset(found)


def _cone(root: Path, command: str, *arguments: str) -> subprocess.CompletedProcess[str]:
    """``python -O -m tools.cone_store <command>`` on the clone, with ``orpheus`` as the first-party package (the
    fixture's ``cli`` names ``pkg``)."""
    where = paths(root)
    argv = [sys.executable, "-O", "-m", "tools.cone_store", command, "--repo", str(root), "--store", str(where["store"]),
            "--package", "orpheus", "--memo-root", str(where["memo"]), *arguments]
    run = subprocess.run(argv, cwd=root, env=environment(root), capture_output=True, text=True, timeout=3600)
    if run.returncode != 0:
        pytest.fail(f"{' '.join(argv[:6])} exited {run.returncode}:\n{run.stdout[-3000:]}\n{run.stderr[-3000:]}")
    return run


def _select_ids(root: Path) -> frozenset[str]:
    run = _cone(root, "select", *_SAMPLE)
    return frozenset(line.strip() for line in run.stdout.splitlines() if "::" in line)


@pytest.mark.slow
@pytest.mark.rests_on(f"{_GT}::test_t2_a_mechanism_sample_of_the_real_suite_contains_its_alone_cones", _SG0)
def test_sg4_control_1_exact_on_the_production_tree(tmp_path: Path) -> None:
    """SG4 (positive control 1 exact, production tree). In a clone of the repository, the sample is captured, then
    ``_OUTFLOW_COSINE_TOL`` (``orpheus/mesh/axis.py``, read only by ``face_outflow_ordinates``) is edited: ``select``
    lists exactly the sample's tests whose run started that def, by an independent monitor. The sample holds at least
    one reader and one test that imports ``axis.py`` and never reaches the reader, or exactness is vacuous. RED
    at B1: every test of the sample was selected (each pinned ``axis.py``'s skeleton)."""
    root = _clone(tmp_path / "clone")
    plugin_dir = tmp_path / "plugin"
    plugin_dir.mkdir()
    (plugin_dir / "gt_axis_monitor.py").write_text(_MONITOR)
    readers = _readers_by_monitor(root, plugin_dir)
    non_readers = frozenset(_SAMPLE) - readers
    assert readers and non_readers, f"the sample does not discriminate: readers {sorted(readers)}"
    _cone(root, "capture", *_SAMPLE)
    assert _select_ids(root) == frozenset(), "the sample is selected with no edit (the negative leg)"
    edit(root, _AXIS, "_OUTFLOW_COSINE_TOL = 1e-15", "_OUTFLOW_COSINE_TOL = 2e-15")
    selected = _select_ids(root)
    assert selected == readers, (f"{len(selected)} of {len(_SAMPLE)} selected; extra (non-readers) "
                                 f"{sorted(selected - readers)}; missing {sorted(readers - selected)}")


# ── qa's final F4 / elegance F1-F2: an import set's rows in the store's spelling ──────────────────────────────

_SPELLING_TREE = {
    "pyproject.toml": """
        [tool.pytest.ini_options]
        addopts = "-p no:cacheprovider"
        testpaths = ["tests"]
        pythonpath = ["src"]
    """,
    "src/libmod.py": """
        def lib_g():
            return 2


        _LIB_TABLE = [lib_g]
        LIB_VALUE = 1
    """,
    "tests/__init__.py": "",
    "tests/test_lib.py": """
        def test_lib():
            import libmod  # inside the test: an import that raises fails the test, not collection

            assert libmod.LIB_VALUE == 1
    """,
}
LIB = "tests/test_lib.py::test_lib"


@pytest.fixture(scope="module", params=["prepend", "importlib"])
def spelling_captured(tmp_path_factory, request: pytest.FixtureRequest) -> Path:
    """The tree captured under one of pytest's import modes. Under ``prepend`` pytest puts the rootdir ahead of ``src``
    on ``sys.path``, so the memo already spells ``libmod`` as ``src/libmod.py``; under ``importlib`` it does not, so the
    memo spells it ``libmod.py`` (relative to ``src``, the first entry holding it) and only the store's relocation
    makes the row a repository file (``[M]`` ``ta/probe_sg7_kinds.py``, 2026-11-02)."""
    root = tmp_path_factory.mktemp(f"spelling_{request.param}") / "tree"
    for relative, source in _SPELLING_TREE.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        text = textwrap.dedent(source).lstrip("\n")
        if relative == "pyproject.toml" and request.param == "importlib":
            text = text.replace('addopts = "-p no:cacheprovider"', 'addopts = "-p no:cacheprovider --import-mode=importlib"')
        path.write_text(text)
    time.sleep(_trusted_age_s())
    capture(root)
    return root


@pytest.mark.parametrize(("edit_kind", "expected"), [
    pytest.param("none", frozenset(), id="no-edit"),
    pytest.param("comment", frozenset(), id="comment"),
    pytest.param("import-raises", frozenset({LIB}), id="import-raises"),
])
def test_sg7_a_module_under_another_sys_path_entry_is_judged_in_the_stores_spelling(
        spelling_captured: Path, tmp_path: Path, edit_kind: str, expected: frozenset[str]) -> None:
    """SG7 (qa's final F4). ``libmod`` lives in ``src/``, a ``sys.path`` entry inside the repository other than its
    root (``pythonpath``), so the memo pins it as ``libmod.py`` and the store as ``src/libmod.py``. The import set's
    rows are re-kept in the store's spelling at capture, so a row is vouched by the file's stage one: no edit and a
    comment edit select nothing; deleting ``lib_g`` (named by the unread ``_LIB_TABLE``) selects ``test_lib``.
    Reddening mutation: skip ``pins.relocated_rows`` (under ``importlib`` the row ``libmod.py`` names no file of the
    store); measure under the selector's ``sys.path`` (the comment leg); derive a root from the pin path (the cone is
    unknown)."""
    from ._cone_store_fixture import copy_captured

    tree = copy_captured(spelling_captured, tmp_path / "tree")
    if edit_kind == "comment":
        edit(tree, "src/libmod.py", "LIB_VALUE = 1", "LIB_VALUE = 1  # the value")
    elif edit_kind == "import-raises":
        edit(tree, "src/libmod.py", "def lib_g():\n    return 2\n", "")
    selected, _ = select(tree)
    assert selected == expected, f"got {sorted(selected)}"


def test_sg8_the_capture_baseline_covers_every_repository_source(tmp_path: Path) -> None:
    """SG8 (elegance's final F0). ``recent_identities`` takes the baseline identity of every file a cone can pin that is
    not settled at the capture's start: a file under ``tools/`` (first-party since B2's ``Reading``) is in it, as are
    the package's and the tests'; ``.git``, ``__pycache__`` and the bytecode cache are not. Reddening mutation: walk only
    ``<package>/`` and ``tests/`` (a ``tools/`` file written just before a capture then makes every cone that ran it
    unknown)."""
    from tools.cone_store.driver import recent_identities

    for relative in ("pkg/a.py", "tests/test_a.py", "tools/helper.py", ".git/HEAD", "pkg/__pycache__/a.cpython.pyc"):
        (tmp_path / relative).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / relative).write_text("x = 1\n")
    found = recent_identities(tmp_path, time.time_ns())
    assert {"pkg/a.py", "tests/test_a.py", "tools/helper.py"} <= set(found), sorted(found)
    assert not [p for p in found if p.startswith(".git") or "__pycache__" in p], sorted(found)


_F5_NODE = "tests/gates/test_layer_imports.py::test_the_points_and_the_group_action_are_geometry"


@pytest.mark.slow
def test_sg9_the_capture_pins_none_of_its_own_code_on_the_production_tree(tmp_path: Path) -> None:
    """SG9 (the main agent's F5, ``[M]`` 533 of 533 cones of ``test_layer_imports.py`` pinned six ``worker.py`` defs and
    the statements of ``cone``, ``pins``, ``store`` and ``worker``). One real test captured into a store of its own: no
    pin of its cone names ``tools/cone_store/worker.py``, and no statement pin names ``tools/cone_store/{cone,pins,
    store}.py`` (the capture's code runs in every capture and in no full run). Positive control: the cone holds the
    test's own def. Reddening mutations: keep the worker's code objects in the recording, or import ``cone`` (and so
    ``pins`` and ``store``) before the recording stops. Since B3 both cones are asked (qa's F6c: the strip of the own
    cone had no witness, the fixture trees holding no worker); ``[M]`` 2026-11-03 on a clone of the tree whose
    ``Recorded.cone_of`` keeps the worker's code: red on the selection cone (six worker ``DefPin``); with the strip
    skipped for the own cone alone (``run`` reads it through ``cone.cone_of``): red on the own cone (one worker
    ``DefPin``, its ``PresencePin`` and ``StatementsPin``), green before this leg."""
    import tools.cone_store as cone

    store = tmp_path / "f5.sqlite"
    argv = [sys.executable, "-O", "-m", "tools.cone_store", "capture", "--repo", str(REPO), "--store", str(store),
            "--memo-root", str(tmp_path / "memo"), _F5_NODE]
    run = subprocess.run(argv, cwd=REPO, env=dict(os.environ, OMP_NUM_THREADS="1"), capture_output=True, text=True,
                         timeout=1800)
    if run.returncode != 0:
        pytest.fail(f"capture exited {run.returncode}:\n{run.stdout[-3000:]}\n{run.stderr[-3000:]}")
    opened = cone.open_store(store)
    # Both cones of the test (B3, qa's F6c): the selection cone and the own cone, which a second recording started
    # inside the worker's ``run`` reads, so the worker's code that runs the test's protocol is inside it too.
    for which, pins in (("selection", opened.pins(_F5_NODE)), ("own", opened.own_pins(_F5_NODE))):
        paths_of = [(type(p).__name__, str(getattr(p, "path", ""))) for p in pins]
        assert ("DefPin", "tests/gates/test_layer_imports.py") in paths_of, f"the {which} cone lacks the test's own def"
        worker = [p for p in paths_of if p[1] == "tools/cone_store/worker.py"]
        readers = [p for p in paths_of if p[0] == "StatementsPin"
                   and p[1] in {"tools/cone_store/cone.py", "tools/cone_store/pins.py", "tools/cone_store/store.py"}]
        assert not worker and not readers, f"the capture's own code is in the {which} cone: {worker + readers}"



# ── the import outcome's roots in the importing process's precedence (the main agent's cause (b), 2026-11-02) ──────

_PRECEDENCE_TREE = {
    "pyproject.toml": """
        [tool.pytest.ini_options]
        addopts = "-p no:cacheprovider"
        testpaths = ["tests"]
    """,
    "tests/__init__.py": "",
    "tests/test_o.py": """
        import orpheus.data.citation  # an ORPHEUS module the capture has not imported: its root joins the set


        def test_o():
            assert orpheus.data.citation is not None
    """,
}
O = "tests/test_o.py::test_o"


def test_sg10_a_tree_package_shadowed_by_the_repositorys_own_is_measured_in_the_runs_precedence(tmp_path: Path) -> None:
    """SG10. The tree's tests are a package ``tests`` under its rootdir, and the run also imports an ORPHEUS module, whose
    root (this repository) holds a ``tests`` package of its own. The import outcome re-imports ``tests.test_o`` with the
    roots in the IMPORTING process's ``sys.path`` precedence, so it finds the tree's package: the cone is known and holds
    an import-outcome pin, and no edit selects nothing. Reddening mutation: roots in the import set's order (ORPHEUS's
    root first: ``No module named 'tests.test_o'``, the cone unknown, the test always selected)."""
    root = tmp_path / "tree"
    for relative, source in _PRECEDENCE_TREE.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(textwrap.dedent(source).lstrip("\n"))
    time.sleep(_trusted_age_s())
    capture(root)
    record = store_of(root).tests()[O]
    assert not record.unknown, f"the cone is unknown: {record.unknown}"
    assert any(type(p).__name__ == "ImportOutcomePin" for p in store_of(root).pins(O)), "no import-outcome pin"
    selected, _ = select(root)
    assert selected == frozenset(), f"got {sorted(selected)}"


# ── the guarded probe of a name new to a module (B2 gate spec row U2, a recall row since ruling 43, R4) ──────────────

_PROBE_TREE = {
    "pyproject.toml": """
        [tool.pytest.ini_options]
        addopts = "-p no:cacheprovider"
        testpaths = ["tests"]
    """,
    "pkg/__init__.py": "",
    "pkg/target.py": """
        PRESENT_VALUE = 1
    """,
    "tests/__init__.py": "",
    "tests/test_u2.py": """
        from pkg import target


        def test_probe():
            assert getattr(target, "U2_NEW_NAME", 0) == 0
    """,
}
PROBE = "tests/test_u2.py::test_probe"


def test_u2_a_guarded_probe_of_a_name_new_to_its_module_selects_the_prober(tmp_path: Path) -> None:
    """U2 (``scratch/verification_tree/cone_store/b2/gate_spec.md``, "now recall"; the theory page's boundary list).
    ``test_probe`` reads ``U2_NEW_NAME`` through ``getattr(module, name, default)``, and ``pkg/target.py`` does not bind
    it at capture. Adding ``U2_NEW_NAME = 1`` to the module changes the test's verdict, and ``select`` lists it: the
    name is an entry name of the run (a string token of its code), and the entry names are stored whole, so the new
    statement binding it is observed. Precision leg: adding a statement binding a name the run never names selects
    nothing. Negative leg: no edit selects nothing.

    ``[M]`` 2026-11-03: green. First red ``[M]``: on a clone of the tree whose ``collect_pins`` stores each module's
    entry names restricted to the identifiers of its source at capture (ruling 43's R4 fallback), the edit selects
    nothing."""
    root = tmp_path / "tree"
    for relative, source in _PROBE_TREE.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(textwrap.dedent(source).lstrip("\n"))
    time.sleep(_trusted_age_s())
    capture(root)
    assert ids_of(root) == {PROBE}
    selected, _ = select(root)
    assert selected == frozenset(), f"no edit: {sorted(selected)}"
    edit(root, "pkg/target.py", "PRESENT_VALUE = 1\n", "PRESENT_VALUE = 1\nU2_OTHER_VALUE = 2\n")
    selected, _ = select(root)
    assert selected == frozenset(), f"a name the run never names: {sorted(selected)}"
    edit(root, "pkg/target.py", "U2_OTHER_VALUE = 2\n", "U2_OTHER_VALUE = 2\nU2_NEW_NAME = 1\n")
    selected, _ = select(root)
    assert selected == {PROBE}, f"the probed name, new to the module: {sorted(selected)}"
