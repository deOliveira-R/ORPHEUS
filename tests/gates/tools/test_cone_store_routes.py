"""The cone store's routes: the fork guard at each level, the fresh-interpreter fallback, and "cone unknown".

Spec ``scratch/verification_tree/cone_store/spec_b1.md`` §3.5, rows G1-G4 and U1-U9. Each test's route
(``fork`` or ``fresh``) is a column of its record, which these rows read: a fallback that silently forks
anyway, or a fork that silently degrades to fresh, gives the same rows and can only be seen there.
"""
from __future__ import annotations

from pathlib import Path

import pytest

import tools.cone_store as cone

from ._cone_store_fixture import (BASE_IDS, BASELINE, SPAWN_AT_COLLECTION_IDS, SPAWNER_IDS, Captured, build, capture,
                                  capture_base, code_keys, paths, select, store_of)

pytestmark = pytest.mark.foundation


@pytest.fixture(scope="module")
def plain(tmp_path_factory) -> Path:
    root = build(tmp_path_factory.mktemp("routes-plain") / "tree", memo=False)
    capture(root)
    return root


def _names_only(pins) -> frozenset[tuple]:
    """Code rows by node, without digests: a variant's conftest has another skeleton digest."""
    return code_keys(pins)


def test_g1_a_thread_in_the_conftest_makes_every_capture_fresh(plain: Path, tmp_path: Path) -> None:
    """G1. A daemon thread started by ``tests/conftest.py`` leaves the grandparent with 2 OS threads: every test
    is captured by route ``fresh``, and its rows name the same nodes as the plain tree's fork capture."""
    root = build(tmp_path / "tree", memo=False, thread_at_conftest=True)
    capture(root)
    records, reference = store_of(root).tests(), store_of(plain)
    assert set(records) == BASE_IDS
    assert {t: r.route for t, r in records.items()} == dict.fromkeys(BASE_IDS, "fresh")
    for test in BASE_IDS:
        assert _names_only(store_of(root).pins(test)) == _names_only(reference.pins(test)), test


def test_g2_a_thread_started_by_one_module_makes_only_that_module_fresh(tmp_path: Path) -> None:
    """G2. ``tests/test_thread.py`` starts a thread when it is imported: its collecting parent holds 2 OS threads,
    so its tests are ``fresh``, and every other module's tests are ``fork``."""
    root = build(tmp_path / "tree", memo=False, thread_at_module=True)
    capture(root)
    routes = {t: r.route for t, r in store_of(root).tests().items()}
    assert routes.pop("tests/test_thread.py::test_thread") == "fresh"
    assert routes == dict.fromkeys(BASE_IDS, "fork")


def test_g3_an_unmeasured_platform_makes_every_capture_fresh(tmp_path: Path) -> None:
    """G3. With ``traced_memo.fork_platform`` reporting linux/x86_64/openblas (not in ``FORK_MEASURED``), set by
    the fixture's conftest, every test is captured ``fresh``."""
    root = build(tmp_path / "tree", memo=False, unmeasured_platform=True)
    capture(root)
    assert {t: r.route for t, r in store_of(root).tests().items()} == dict.fromkeys(BASE_IDS, "fresh")


def test_g4_the_fork_route_and_the_fresh_route_write_the_same_rows(plain: Path, tmp_path: Path) -> None:
    """G4. A capture forced fresh for every test (``--route fresh``) writes the same rows as the fork capture
    (repair (1)'s promise; the two routes share the recorder and nothing else)."""
    fresh_store = tmp_path / "fresh.sqlite"
    capture(plain, "--route", "fresh", "--store", str(fresh_store))
    forked, fresh = store_of(plain), store_of(plain, fresh_store)
    assert {t: r.route for t, r in fresh.tests().items()} == dict.fromkeys(BASE_IDS, "fresh")
    assert {t: r.route for t, r in forked.tests().items()} == dict.fromkeys(BASE_IDS, "fork")
    for test in BASE_IDS:
        assert forked.pins(test) == fresh.pins(test), test


@pytest.fixture(scope="module")
def spawners(tmp_path_factory) -> Path:
    root = build(tmp_path_factory.mktemp("routes-spawn") / "tree", memo=False, spawners=True)
    capture(root)
    return root


@pytest.mark.parametrize("api", sorted(SPAWNER_IDS))
def test_u1_to_u5_a_test_that_starts_a_process_is_cone_unknown(spawners: Path, api: str) -> None:
    """U1-U5. A test starting a process by ``subprocess.run``, ``os.system``, ``multiprocessing`` (spawn),
    ``os.fork`` or ``os.posix_spawn`` is marked unknown, and ``select`` lists it with no edit."""
    test = SPAWNER_IDS[api]
    assert store_of(spawners).tests()[test].unknown, f"{test} is not cone unknown"
    selected, _ = select(spawners)
    assert test in selected


def test_u6_a_process_started_at_collection_marks_its_module_unknown(spawners: Path) -> None:
    """U6. ``tests/test_spawn_c.py`` starts a process at module level: both its tests are unknown, and no test of
    a module that starts nothing."""
    records = store_of(spawners).tests()
    assert all(records[t].unknown for t in SPAWN_AT_COLLECTION_IDS)
    assert not any(records[t].unknown for t in BASE_IDS), sorted(t for t in BASE_IDS if records[t].unknown)


def test_u7_the_stores_own_forks_are_not_the_tests(plain: Path) -> None:
    """U7. No test of the plain tree is unknown: the store's own forks are not a test's spawn."""
    records = store_of(plain).tests()
    assert set(records) == BASE_IDS
    assert not any(r.unknown for r in records.values()), {t: r.unknown for t, r in records.items() if r.unknown}
    selected, _ = select(plain)
    assert selected == BASELINE


@pytest.fixture(scope="module")
def memo_base(tmp_path_factory) -> Captured:
    return capture_base(tmp_path_factory.mktemp("routes-memo") / "tree")


def test_u8_a_memo_generation_is_not_the_tests_spawn(memo_base: Captured) -> None:
    """U8. On the cold capture ``test_reader`` missed, the memo generated ``oracle(2)`` in its own process, and
    ``test_reader`` was not unknown (its ChildPin covers what the generation read)."""
    memo = paths(memo_base.root)["memo"]
    assert any(p.is_file() for p in memo.rglob("*")), "precondition: the cold capture generated nothing"
    assert memo_base.cold_unknown["tests/test_e.py::test_reader"] == ()


def test_u9_an_unknown_test_carries_its_reason_and_status_counts_it(spawners: Path) -> None:
    """U9. Each unknown test records why (the program or the event), and ``status`` counts the unknown tests:
    5 spawners plus 2 tests of the module that spawns at collection."""
    records = store_of(spawners).tests()
    reasons = records[SPAWNER_IDS["os.system"]].unknown
    assert reasons and all(isinstance(r, str) and r for r in reasons)
    assert any("true" in r or "os.system" in r for r in reasons), reasons
    status = cone.status(spawners, paths(spawners)["store"])
    assert status.unknown_tests == len(SPAWNER_IDS) + len(SPAWN_AT_COLLECTION_IDS)
    assert status.tests == len(records)
