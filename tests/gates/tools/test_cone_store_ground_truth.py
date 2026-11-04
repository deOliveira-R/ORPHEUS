"""Positive control 4: the store's cone contains each test's cone captured alone in a fresh interpreter.

Spec ``scratch/verification_tree/cone_store/spec_b1.md`` §3.6, rows T1-T4 and CU-OVER. Containment, 0 misses,
over every pin kind (defs, import-time statements, data, memo children). The alone capture is the store's own recorder
run by ``--route fresh`` on one node id, so route is the only axis that differs (X4): the recorder's own
ground is R1-R5. ``[M]`` first red of T1: the single shared parent loses five nodes of ``test_b`` (§0.1).

What containment cannot see (spec §7 finding 1): over-selection. CU-OVER records it; nothing bounds it.
"""
from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

import tools.cone_store as cone

from ._cone_store_fixture import BASE_IDS, MEMO_IDS, REPO, Captured, capture, capture_base, paths, store_of

pytestmark = pytest.mark.foundation

#: T2's sample, by mechanism (spec §3.6). The spawner is asserted unknown in both captures and excluded
#: from containment.
T2_SAMPLE = {
    "parametrize-fill": (
        "tests/gates/numerics/test_generating_measure.py::test_legendre_matches_numpy_leggauss[5]",
        "tests/gates/numerics/test_generating_measure.py::test_jacobi_matches_scipy[0.0-0.0-2]",
        "tests/gates/numerics/test_generating_measure.py::test_weights_sum_to_the_measures_mass[1-legendre-measure0]",
        "tests/gates/numerics/test_generating_measure.py::test_is_symmetric_is_derived_correctly[legendre-measure0-True]",
    ),
    "in-module-cache": (
        "tests/gates/numerics/test_invariance.py::TestR2ImportGraph::test_a2_symmetry_does_not_import_manifold_or_measure_at_runtime",
        "tests/gates/numerics/test_invariance.py::TestR2ImportGraph::test_a4_manifold_no_longer_duck_types_the_group",
        "tests/gates/numerics/test_invariance.py::TestR2OneClosure::test_b3_predicate_equals_certificate_existence_for_finite_groups[lebedev(5)]",
    ),
    "plain": (
        "tests/gates/mesh/test_mesh1d_regions.py::test_r7b2_0_1_the_constructor_takes_the_labels",
        "tests/gates/mesh/test_mesh1d_regions.py::test_r7b2_0_2_refusal[label-count]",
        "tests/gates/mesh/test_mesh1d_regions.py::test_r7b2_0_2_refusal[label-float]",
    ),
    "tests-helper": (
        "tests/gates/mesh/test_face_roles.py::TestFaceRolesOver::test_coord_is_required_and_keyword_only",
        "tests/gates/mesh/test_face_roles.py::TestFaceRolesOver::test_tags_build_in_inventory_order[slab]",
    ),
    "memo-children": (
        "tests/gates/derivations/test_traced_memo_clients.py::test_m4_14_a_point_reading_pins_the_blocks_as_children_and_assembles_none",
        "tests/gates/derivations/test_traced_memo_clients.py::test_m4_9_the_exact_medium_is_read_in_the_asking_process",
    ),
}
T2_SPAWNER = "tests/gates/tools/test_test_durations.py::test_node_ids_are_the_ones_pytest_collects"
T2_IDS = tuple(t for group in T2_SAMPLE.values() for t in group)

#: Round 2's sample N: 200 tests collected (``[M]`` 2026-11-01).
T3_FILES = ("tests/gates/numerics/test_invariance.py", "tests/gates/numerics/test_measure.py",
            "tests/gates/numerics/test_space_of_axes.py")


def _names(pin) -> frozenset[str]:
    """The entry names a ``StatementsPin`` keeps (``.cache/sets``); a missing set reads as none, so it cannot hide a
    miss."""
    from orpheus.numerics import traced_memo

    return frozenset(traced_memo._read_kept_set(pin.names) or ())


def _observed(repo: Path, path: str, names: frozenset[str]) -> frozenset[int]:
    """The import-time statements of ``path`` (relative to ``repo``, as it reads now) that a run with entry
    names ``names`` observes: the selection a ``StatementsPin`` projects (``_import_time.observed_digests``)."""
    from orpheus.numerics import _import_time, traced_memo

    return frozenset(_import_time.selection(traced_memo._source_index((repo / path).read_bytes()).statements, names))


def _misses(store, alone, test: str, repo: Path = REPO) -> list:
    """The alone rows the store's cone does not contain, for a store captured over ``repo``. A ``StatementsPin`` (B2, #689) is contained when the store
    holds one for the same file whose observed statements include every statement the alone run observes there. The
    entry names are one set per run, over every file, so the alone run's names may hold a name the store's lack that
    binds nothing in this file (``[M]`` 2026-11-03: two T3 runs read such a name set as a miss, the pre-fix tree's
    included, while the alone run's projection equalled the store's in 10 of 10 printed pins); comparing the names
    read that as a miss, comparing the observed statements does not. Every other pin is contained by identity."""
    held = set(store.pins(test))
    names: dict[str, frozenset[str]] = {}
    for pin in held:
        if type(pin).__name__ == "StatementsPin":
            names[pin.path] = names.get(pin.path, frozenset()) | _names(pin)
    missed = []
    for pin in alone.pins(test):
        if type(pin).__name__ == "StatementsPin":
            if pin.path not in names or not _observed(repo, pin.path, _names(pin)) <= _observed(repo, pin.path, names[pin.path]):
                missed.append(pin)
        elif pin not in held:
            missed.append(pin)
    return sorted(missed, key=repr)


def _excess(store, alone, test: str) -> int:
    return len(set(store.pins(test)) - set(alone.pins(test)))


def _real_cli(command: str, store: Path, *arguments: str, timeout: float) -> subprocess.CompletedProcess[str]:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    argv = [sys.executable, "-O", "-m", "tools.cone_store", command, "--repo", str(REPO), "--store", str(store), *arguments]
    run = subprocess.run(argv, cwd=REPO, env=env, capture_output=True, text=True, timeout=timeout)
    if run.returncode != 0:
        raise AssertionError(f"{' '.join(argv)} exited {run.returncode}:\n{run.stdout[-3000:]}\n{run.stderr[-3000:]}")
    return run


@pytest.fixture(scope="module")
def fixture_pair(tmp_path_factory) -> tuple[Captured, object]:
    """The base tree captured by fork, and each of its tests captured alone, fresh, into a second store."""
    base = capture_base(tmp_path_factory.mktemp("truth") / "tree")
    alone_store = base.root.parent / "alone.sqlite"
    for test in sorted(BASE_IDS | MEMO_IDS):
        capture(base.root, "--route", "fresh", "--store", str(alone_store), test)
    return base, store_of(base.root, alone_store)


def test_t1_every_fixture_cone_contains_its_alone_cone(fixture_pair) -> None:
    """T1. For every test of the fixture, the store's rows contain the alone rows, 0 misses."""
    base, alone = fixture_pair
    store = base.store
    assert set(alone.tests()) == BASE_IDS | MEMO_IDS
    misses = {t: _misses(store, alone, t, base.root) for t in alone.tests()}
    assert not any(misses.values()), {t: m for t, m in misses.items() if m}
    print({t: _excess(store, alone, t) for t in sorted(alone.tests())})


def test_t4_the_fork_and_the_alone_run_give_the_same_outcomes_on_the_fixture(fixture_pair) -> None:
    """T4 (fixture half). Outcomes agree for every test; ``test_c3`` (green only when no earlier test of its
    module ran in its process) is the input that reddens a capture that does not fork per test."""
    base, alone = fixture_pair
    forked, fresh = base.store.tests(), alone.tests()
    assert forked["tests/test_c.py::test_c3"].outcome == "passed"
    assert {t: r.outcome for t, r in forked.items()} == {t: r.outcome for t, r in fresh.items()}


@pytest.fixture(scope="module")
def real_pair(tmp_path_factory) -> tuple[object, object, float]:
    where = tmp_path_factory.mktemp("truth-real")
    started = time.monotonic()
    _real_cli("capture", where / "store.sqlite", *T2_IDS, T2_SPAWNER, timeout=120)
    for test in (*T2_IDS, T2_SPAWNER):
        _real_cli("capture", where / "alone.sqlite", "--route", "fresh", test, timeout=120)
    elapsed = time.monotonic() - started
    return cone.open_store(where / "store.sqlite"), cone.open_store(where / "alone.sqlite"), elapsed


def test_t2_a_mechanism_sample_of_the_real_suite_contains_its_alone_cones(real_pair) -> None:
    """T2. 15 real tests by mechanism (parametrize fill, in-module cache, plain, a ``tests/`` helper, memo
    children, a spawner): every non-spawner's store rows contain its alone rows, 0 misses. Activation leg:
    the 14 non-spawners were captured by ``fork`` (a real conftest that started a thread pool would make every
    route ``fresh``, and fresh against fresh is green for every design). Budget: 120 s per command."""
    store, alone, elapsed = real_pair
    records = store.tests()
    assert {t: records[t].route for t in T2_IDS} == dict.fromkeys(T2_IDS, "fork")
    assert records[T2_SPAWNER].unknown and alone.tests()[T2_SPAWNER].unknown
    misses = {t: _misses(store, alone, t) for t in T2_IDS}
    assert not any(misses.values()), {t: m[:10] for t, m in misses.items() if m}
    print(f"T2 capture time {elapsed:.1f} s; excess per test", {t: _excess(store, alone, t) for t in T2_IDS})


def test_t2_a_real_memo_reader_holds_what_its_own_stores_entries_depended_on(real_pair) -> None:
    """T2, reads leg: ``test_m4_14`` reads a point value from a memo store under its own tmp path, an entry whose
    children are the solve and one block per group. Its cone holds no ``ChildPin``, which would be judged against the
    default store where the entry is absent (spec §7 finding 4: the test was always selected, until the B3 acceptance
    fix of 2026-11-03), and holds the block assembly's ``LineRule.transport``, which ran only in a block's generating
    process: the child is pinned by what its entry depended on, two children deep."""
    store, _, _ = real_pair
    pins = store.pins(T2_SAMPLE["memo-children"][0])
    if children := [p for p in pins if type(p).__name__ == "ChildPin"]:
        pytest.fail(f"test_m4_14 holds ChildPins of its own store: {children[:3]}")
    transport = ("DefPin", "orpheus/derivations/continuous/characteristic/assembly.py", "LineRule.transport")
    if transport not in {(type(p).__name__, getattr(p, "path", ""), getattr(p, "qualname", "")) for p in pins}:
        pytest.fail("test_m4_14 does not hold LineRule.transport, which its entries' generations ran")


def test_t4_the_fork_and_the_alone_run_give_the_same_outcomes_on_the_sample(real_pair) -> None:
    """T4 (real half). Outcomes agree for every T2 test."""
    store, alone, _ = real_pair
    assert {t: store.tests()[t].outcome for t in T2_IDS} == {t: alone.tests()[t].outcome for t in T2_IDS}


@pytest.mark.slow
def test_t3_round_2s_sample_n_contains_its_alone_cones(tmp_path: Path) -> None:
    """T3. Round 2's sample N, 200 tests, at the grain round 2 omitted (module code, ``tests/``, reads): 0 misses
    (:func:`_misses`). The activation leg: some alone pin observes a statement by name, beyond those every run
    observes by effect, so a store that lost names would read as a miss. Budget about 45 min."""
    store_path, alone_path = tmp_path / "store.sqlite", tmp_path / "alone.sqlite"
    _real_cli("capture", store_path, *T3_FILES, timeout=600)
    store = cone.open_store(store_path)
    tests = sorted(store.tests())
    assert len(tests) == 200, len(tests)
    for test in tests:
        _real_cli("capture", alone_path, "--route", "fresh", test, timeout=120)
    alone = cone.open_store(alone_path)
    by_name = sum(1 for t in tests for pin in alone.pins(t) if type(pin).__name__ == "StatementsPin"
                  and _observed(REPO, pin.path, _names(pin)) - _observed(REPO, pin.path, frozenset()))
    assert by_name, "no alone pin observes a statement by name: the statement leg of the predicate cannot fire"
    misses = {t: _misses(store, alone, t) for t in tests}
    assert not any(misses.values()), {t: m[:10] for t, m in misses.items() if m}


def test_cu_over_the_grandparents_share_is_recorded(fixture_pair) -> None:
    """CU-OVER (RECORD). ``status`` reports how many pins the grandparent attributes to every test. On the
    fixture it holds at least ``conftest_helper`` and the conftest's statement pin (the positive control of the
    count); on the real tree it is the 91 modules the root conftest imports (spec §0.2), printed here only."""
    base, _ = fixture_pair
    status = cone.status(base.root, paths(base.root)["store"])
    assert status.grandparent_pins >= 2, status
    print(f"fixture grandparent pins: {status.grandparent_pins}")
