"""The cone store's attribution of collection-time code: repair (1) of ruling 38, and the item-6 separation.

Spec ``scratch/verification_tree/cone_store/spec_b1.md`` §0.1 and §3.2, rows C1-C6 and H1. A grandparent
imports pytest and the conftests; per test module it forks a collecting parent that collects that module
only; that parent forks one child per test. ``[M]`` 2026-11-01, the single shared parent of ruling 38 as
first written loses, for ``test_b``, the defs ``_build`` and ``table`` and three modules' code (§0.1), which
is the first red of C2, C3, C5 and C6.
"""
from __future__ import annotations

from pathlib import Path

import pytest

import tools.cone_store as cone  # noqa: F401  (red at collection until B1 lands)

from ._cone_store_fixture import (BASE_IDS, BASELINE, MEMO_IDS, Captured, capture, capture_base, copy_captured, edit,
                                  in_suite_cones, select, store_of)

pytestmark = pytest.mark.foundation

ALL = BASE_IDS | MEMO_IDS
PKG_USER_READERS = frozenset({"tests/test_a.py::test_a[t0]", "tests/test_b.py::test_b", "tests/test_g.py::test_fail",
                              "tests/test_m.py::test_nested"})


@pytest.fixture(scope="module")
def base(tmp_path_factory) -> Captured:
    return capture_base(tmp_path_factory.mktemp("collection") / "tree")


@pytest.fixture
def tree(base: Captured, tmp_path: Path) -> Path:
    return copy_captured(base.root, tmp_path / "tree")


def test_c1_a_parametrize_fill_is_attributed_to_its_module(tree: Path) -> None:
    """C1. An edit to the body of ``_build`` (run only while ``test_a.py`` was collected, filling ``table(3)``)
    selects ``test_a[t0]``."""
    edit(tree, "pkg/user.py", "return [i * TOL for i in range(n)]", "return [i * TOL + 0 for i in range(n)]")
    selected, _ = select(tree)
    assert "tests/test_a.py::test_a[t0]" in selected


def test_c2_a_cache_filled_by_another_modules_collection_does_not_hide_the_fill(tree: Path) -> None:
    """C2. The same edit selects ``test_b``, of another module collected after ``test_a.py``. ``[M]`` first red:
    the single shared parent (spec §0.1)."""
    edit(tree, "pkg/user.py", "return [i * TOL for i in range(n)]", "return [i * TOL + 0 for i in range(n)]")
    selected, _ = select(tree)
    assert "tests/test_b.py::test_b" in selected


def test_c3_a_constant_edit_selects_every_importer_of_its_reader(tree: Path) -> None:
    """C3 (control 1, also P1). An edit to ``TOL`` selects exactly the tests whose code reads it, beside the baseline:
    ``test_a[t0]``, ``test_b``, ``test_fail``, ``test_nested``. On this tree every importer of ``pkg.user`` reads
    ``TOL`` (``_build`` and ``above``), so the row cannot tell statement grain from module grain; the exact control
    at statement grain is ``test_cone_store_statement_grain.py`` SG1, SG2 and SG4 (B2, #689)."""
    edit(tree, "pkg/consts.py", "TOL = 1e-15", "TOL = 2e-15")
    selected, _ = select(tree)
    assert selected == PKG_USER_READERS | BASELINE, f"extra {sorted(selected - PKG_USER_READERS - BASELINE)}; missing {sorted(PKG_USER_READERS - selected)}"


@pytest.mark.parametrize("change", ["body", "statement"])
def test_c4_what_the_grandparent_ran_is_attributed_to_every_test(tree: Path, change: str) -> None:
    """C4. ``tests/conftest.py`` runs ``conftest_helper`` in the grandparent; an edit to its body, or to ``FLAG``
    (an import-time statement of ``pkg/conf_dep.py`` the helper reads), selects every test of the tree. Re-posed at
    B2 (#689): ``FLAG`` was unread, and the row passed only through the module skeleton; under statement pins that
    edit selects nothing, so the helper now reads ``FLAG`` and the row asserts attribution of an OBSERVED statement."""
    if change == "body":
        edit(tree, "pkg/conf_dep.py", "return 3", "return 1 + 2")
    else:
        edit(tree, "pkg/conf_dep.py", "FLAG = True", "FLAG = False")
    selected, _ = select(tree)
    assert selected == ALL, f"missing {sorted(ALL - selected)}"


def test_c5_a_package_init_is_attributed_to_every_test_under_it(tree: Path) -> None:
    """C5. An import-time effect added to ``tests/__init__.py`` selects every test. ``[M]`` first red: the single
    parent attributed package ``__init__`` code to the first module only, 31 of 31 real samples (spec §0.1).
    Re-posed at B2 (#689): the edit was ``MARK = 1``, an unread binding that statement pins rightly ignore; an
    effect (a warnings filter every test runs under) is what every test observes."""
    (tree / "tests/__init__.py").write_text("import warnings\n\nwarnings.simplefilter('default')\n")
    selected, _ = select(tree)
    assert selected == ALL, f"missing {sorted(ALL - selected)}"


def test_c6_the_rows_do_not_depend_on_the_order_modules_are_captured(base: Captured, tmp_path: Path) -> None:
    """C6. A capture with the test modules given in reverse order writes the same (test, pin) rows as the
    forward capture: repair (1) is order-independent; the single parent is not."""
    reverse = copy_captured(base.root, tmp_path / "reverse")
    files = sorted({test.split("::")[0] for test in ALL}, reverse=True)
    capture(reverse, *files)
    forward, backward = base.store, store_of(reverse)
    assert set(backward.tests()) == set(forward.tests())
    for test in forward.tests():
        assert backward.pins(test) == forward.pins(test), test


def test_h1_the_fork_capture_sees_a_fill_the_in_suite_monitor_cannot(base: Captured) -> None:
    """H1 (brief item 6). ``test_c1`` fills ``hidden.cached_f(4)`` and ``test_c2`` hits it in one process.
    The store's ``test_c2`` rows hold ``_fill_body`` (its child forks from a parent whose cache is empty);
    the in-suite monitor's cone of ``test_c2`` does not (a declared blindness of that instrument, asserted
    so no in-suite mode can claim equivalence silently). First red: NOFORK."""
    fill = ("DefPin", "pkg/hidden.py", "_fill_body", 0)
    from ._cone_store_fixture import code_keys

    assert fill in code_keys(base.store.pins("tests/test_c.py::test_c2"))
    monitor = in_suite_cones(base.root, "tests/test_c.py")
    assert ("pkg/hidden.py", "_fill_body") in monitor["tests/test_c.py::test_c1"], "the monitor saw no fill at all"
    assert ("pkg/hidden.py", "_fill_body") not in monitor["tests/test_c.py::test_c2"]
