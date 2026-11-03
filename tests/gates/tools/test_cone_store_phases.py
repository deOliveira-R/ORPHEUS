"""The cone store's two cones per test (B3): the SELECTION cone ``select`` judges, and the OWN cone placement reads.

Spec ``scratch/verification_tree/cone_store/b3/spec_b3.md`` §2.2, rows PH1-PH5, re-posed for the representation
that landed: a test's selection cone is the ``cones`` relation (everything its process ran, the grandparent's
conftest imports and its module's collection included, ruling 38), its own cone the ``own_cones`` relation (a second
recording started at pytest's protocol for the test, on its own ``sys.monitoring`` tool, so a def the conftests already
ran is reported again when the test runs it), and the grandparent's cone is stored once under ``store.GRANDPARENT``.

The fixture's root conftest runs ``pkg.conf_dep.conftest_helper`` at import. ``test_p`` calls it again in its own run;
``test_b`` never does. The expected sets are read off the fixture's sources by hand.
"""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from tools.cone_store import store as store_module

from ._cone_store_fixture import Captured, capture_base, code_keys, copy_captured, edit, pin_key, select

pytestmark = pytest.mark.foundation

HELPER = ("DefPin", "pkg/conf_dep.py", "conftest_helper", 0)
TEST_P = "tests/test_p.py::test_p"
TEST_B = "tests/test_b.py::test_b"


def nodes(pins) -> set[tuple]:
    """The nodes a cone names (:func:`pin_key`: kind, path, qualname, ordinal), the run's import outcome excepted.
    A ``StatementsPin``'s ``names`` and ``digest`` are a projection on the entry names of the run that recorded it,
    and an ``ImportOutcomePin`` is that run's import set, so two cones of one process name the same module with
    different rows: containment is asked of the nodes, never of the rows."""
    return {pin_key(p) for p in pins if type(p).__name__ != "ImportOutcomePin"}


@pytest.fixture(scope="module")
def base(tmp_path_factory) -> Captured:
    return capture_base(tmp_path_factory.mktemp("phases") / "tree", memo=False, phases=True)


def test_ph1_a_def_the_conftest_ran_is_in_the_own_cone_of_a_test_that_runs_it_again(base: Captured) -> None:
    """PH1. ``test_p``'s own cone holds ``conftest_helper`` although the grandparent ran it first, and so does its
    selection cone. First red: one recorder for both cones (a code object reports once per tool, so a def the
    conftests entered is never reported again, and the own cone lacks it)."""
    own = code_keys(base.store.own_pins(TEST_P))
    assert ("DefPin", "tests/test_p.py", "test_p", 0) in own, "the own cone lost the test's own function (control)"
    assert HELPER in own
    assert HELPER in code_keys(base.store.pins(TEST_P))


def test_ph2_a_def_only_the_conftest_ran_is_in_the_selection_cone_and_not_the_own(base: Captured) -> None:
    """PH2. ``test_b`` never calls ``conftest_helper``: its own cone lacks it, its selection cone holds it. First red:
    the own cone written as the selection cone (the phases not told apart)."""
    own = code_keys(base.store.own_pins(TEST_B))
    assert ("DefPin", "pkg/user.py", "above", 0) in own, "the own cone lost a def the test ran (control)"
    assert HELPER not in own
    assert HELPER in code_keys(base.store.pins(TEST_B))


def test_ph3_selection_reads_the_whole_selection_cone(base: Captured, tmp_path: Path) -> None:
    """PH3, selection soundness is unchanged. (a) Every test's selection cone contains its own cone and the
    grandparent's share, pin for pin (the union property: the selection recording saw everything either did).
    (b) An edit to ``conftest_helper``'s body selects ``test_b``, which reaches the def only through the conftest,
    as well as ``test_p``. First red of (b): ``select`` judging the own cones (``test_b`` drops out).

    ``[M]`` 2026-11-03, the one-off independent reference (not a gate: it needs B2's tools from git): the B2 tools
    of ``83db308f`` and the B3 tools captured the same 7 real tests (``test_concept_order.py``, two of
    ``test_iteration.py``) with selection rows equal for 7 of 7 tests (1860-1943 rows each)."""
    stored = base.store
    share = nodes(stored.share(store_module.GRANDPARENT) or ())
    assert HELPER in share, "the share is empty or lost the conftest's def (control)"
    for test_id in stored.tests():
        selection = nodes(stored.pins(test_id))
        assert nodes(stored.own_pins(test_id)) <= selection, (test_id, sorted(nodes(stored.own_pins(test_id)) - selection))
        assert share <= selection, (test_id, sorted(share - selection))
    tree = copy_captured(base.root, tmp_path / "tree")
    edit(tree, "pkg/conf_dep.py", "return 3 if FLAG else 0", "return 1 + 2 if FLAG else 0")
    selected, _ = select(tree)
    assert {TEST_P, TEST_B} <= selected, f"missing {sorted({TEST_P, TEST_B} - selected)}"


def test_ph4_the_grandparents_cone_is_stored_as_a_set(base: Captured) -> None:
    """PH4. The grandparent's cone is stored (not only its count): it holds ``conftest_helper``, and it is contained
    in every test's selection cone (by node). First red: the share not recorded (``share`` returns ``None``)."""
    share = base.store.share(store_module.GRANDPARENT)
    assert share is not None, "no grandparent share stored"
    assert HELPER in code_keys(share)
    for test_id in base.store.tests():
        assert nodes(share) <= nodes(base.store.pins(test_id)), test_id


def test_ph5_a_store_of_the_previous_format_is_emptied_when_opened(tmp_path: Path) -> None:
    """PH5. A store whose ``format`` is 2 (B2: no own cones, no shares) is emptied when opened, so a B2 row is never
    read as a test with an empty own cone. First red: ``FORMAT`` left at 2."""
    assert store_module.FORMAT != 2
    path = tmp_path / "old.sqlite"
    with store_module.ConeStore(path) as database:
        database.record(store_module.StoredTest("tests/test_x.py::test_x", "key", "fork", "passed", ()), (), {})
        database.set_meta("format", 2)
    with store_module.ConeStore(path) as database:
        assert database.tests() == {}
        assert database.meta("format") == store_module.FORMAT
    tables = {name for (name,) in sqlite3.connect(path).execute("SELECT name FROM sqlite_master WHERE type='table'")}
    assert {"own_cones", "shares"} <= tables
