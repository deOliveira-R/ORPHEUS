"""B5: ``place`` reading a test's reads by channel, with the tiered home (rulings 47, 48 and 49).

Spec ``scratch/verification_tree/cone_store/b5/spec_b5.md`` §3.4-§3.6, rows V1-V10 (placement by channel), D1-D3
(the data home) and A1-A9 (the apparatus). The fixture is ``_READS`` (``build(root, reads=True)``): a layer
``pkg/rd`` ranked ``low < mid, rex < high`` (:data:`READS_ORDER`), ``pkg/other`` another package. Every expectation is
read off those sources by hand.

The home is tiered (ruling 49): def nodes decide it when the test has any; else statement nodes (the channels
``direct``, ``decorator``, ``string``, ``via_apparatus``, and ``module``, a read of the module object, C8); else data nodes (C13: ``layer`` when the reads cover the
layer alone, ``package`` when they walk the package); else the apparatus (C14, the test-tree file). Each test entry
carries ``via`` (layer node to its sorted channels among ``def``, ``statement``, ``data``), ``items`` (node, a layer
node name or an apparatus file path, to its sorted read names), ``home_union``, ``statement_straddle`` (node to kind,
as ``foreign``), ``placed_by``, ``unresolved`` (the sorted dotted spellings) and ``apparatus_items`` (the apparatus items it
rests on, sorted ``[path, name]`` pairs). The result carries ``apparatus_supports`` (proposed premise edges
``{upper, lower, node, c}`` over test functions), each declared ``rests_on`` row gains ``c_items`` (containment over
defs and items), and the summary gains ``homes_union``, ``union_moves``, ``statement_straddlers``,
``above_layer_statement``, ``pure_fraction_statement`` and ``null_pure_fraction_statement``.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from ._cone_store_fixture import READS_ORDER, build, capture, paths

pytestmark = pytest.mark.foundation

PREFIX = "tests/rd/"
T = "tests/rd/test_rd.py::"
A = "tests/rd/test_audit.py::"
TEST_RD = "tests/rd/test_rd.py"
HELPER = "tests/rd/_helper.py"


@pytest.fixture(scope="module")
def captured(tmp_path_factory) -> Path:
    root = build(tmp_path_factory.mktemp("placement_reads") / "tree", memo=False, reads=True)
    capture(root, "tests/rd")
    return root


def _place(root: Path, declared: tuple[tuple[str, str], ...] = (), store: Path | None = None) -> dict:
    from tools.cone_store import placement

    from tests.gates.numerics.test_concept_order import import_kinds

    layer = placement.Layer("pkg/rd", "pkg", ())
    return placement.place(store if store is not None else paths(root)["store"], PREFIX, layer, READS_ORDER, import_kinds(root / "pkg/rd"),
                           declared=declared)


@pytest.fixture(scope="module")
def placed(captured: Path) -> dict:
    return _place(captured)


def row(placed: dict, name: str) -> dict:
    test_id = name if "::" in name else T + name
    assert test_id in placed["tests"], f"{test_id} not placed or listed"
    return placed["tests"][test_id]


# ── group V, placement by channel ──


def test_v1_a_value_read_from_a_higher_node_does_not_move_a_def_home(placed: dict) -> None:
    """V1 (ruling 49). ``test_input_value`` runs ``low.f`` and reads ``high.HIGH``: home ``low`` by def;
    ``high`` is a statement node, the union home, and a statement straddle (``high`` imports nothing); purity is
    read from defs alone, so the test is pure. First red: the union home."""
    entry = row(placed, "test_input_value")
    assert entry["home"] == "low" and entry["placed_by"] == "def"
    assert entry["via"] == {"low": ["def", "statement"], "high": ["statement"]}
    assert entry["home_union"] == "high"
    assert "high" in entry["statement_straddle"]
    assert entry["pure_rt"] is True


def test_v2_a_reexporter_is_not_a_node(placed: dict) -> None:
    """V2 (C7, ruled: pinned, not a node). ``test_reexport`` reads ``TOL`` through ``rex``: home ``low`` by
    statement; ``rex`` is absent from ``via``. First red: home at the re-exporter (``rex`` ranks above ``low``)."""
    entry = row(placed, "test_reexport")
    assert entry["home"] == "low" and entry["placed_by"] == "statement"
    assert "rex" not in entry["via"]


def test_v3_a_single_file_audit_homes_by_data(placed: dict) -> None:
    """V3 (C13). ``test_audit_one`` opens ``pkg/rd/high.py`` as text, runs no layer code and reads no binding:
    home ``high`` by data. First red: data in no channel (the test stays unplaced)."""
    entry = row(placed, A + "test_audit_one")
    assert entry["home"] == "high" and entry["placed_by"] == "data"
    assert entry["via"] == {"high": ["data"]}


#: The def homes of the fixture, by hand: the tests that run a layer def, and the node of the highest-ranked one.
_DEF_HOMES = {T + "test_def_and_assign": "low", T + "test_input_value": "low", T + "test_uses_premise": "low",
              T + "test_superinstruction": "high", T + "test_closure_unbound": "high", T + "test_two_cells": "high",
              T + "test_fast_pair": "low", A + "test_audit_run": "mid"}


def test_v4_every_def_home_is_kept(placed: dict) -> None:
    """V4 (AC1 in the fixture). Every test that runs a layer def keeps its def home, whatever it reads beside it.
    First red: any channel admitted into the def tier."""
    by_def = {t: e["home"] for t, e in placed["tests"].items() if e.get("placed_by") == "def"}
    print(f"V4: {len(by_def)} def-placed tests")
    assert by_def == _DEF_HOMES


def test_v5_a_statement_node_carries_its_items(placed: dict) -> None:
    """V5. ``test_global``: ``items == {low: [TOL]}``, placed by statement. First red: the node without its items."""
    entry = row(placed, "test_global")
    assert entry["placed_by"] == "statement" and entry["home"] == "low"
    assert entry["items"] == {"low": ["TOL"]}


def test_v6_the_null_is_shuffled_per_channel(placed: dict) -> None:
    """V6 (C5; attack a21). The null is read per channel, three seeded draws each, on one denominator: the
    statement-placed tests. The mutation it guards, shuffling ``(node, def)`` labels only, leaves every
    statement-placed test unshuffled, so its null equals the measured fraction on every seed; the row asserts that
    the statement null differs from the measured fraction on at least one seed. A RECORD of the fixture with its
    seeds, as PL16 is: ``[M]`` 2026-11-07 (the orchestrator's first green run of the build) the measured statement
    fraction is 23/24 (``test_two_reads_apart``, ``mid`` and ``high`` unrelated, is the one straddler) and the null
    reads 24/24, 22/24 and 24/24 (seeds 1, 2, 3; 22/23, 23/23, 21/23 and 23/23 before
    ``test_rebound`` joined the fixture). A null may read above the measured fraction: on one straddler a
    shuffle can map its two nodes onto one. First red: ``(node, def)`` labels only shuffled (the null is 23/24 on
    every seed)."""
    summary = placed["summary"]
    assert len(summary["null_pure_fraction"]) == 3, summary["null_pure_fraction"]
    measured, null = summary["pure_fraction_statement"], summary["null_pure_fraction_statement"]
    n = 24
    assert len(null) == 3, null
    for fraction in (measured, *null.values()):
        assert abs(fraction * n - round(fraction * n)) < 1e-4, f"{fraction} is not a fraction of {n}: the fixture changed"
    assert round(measured * n) == 23, f"measured {measured}"
    assert any(value != measured for value in null.values()), f"the statement null equals the measured {measured}: {null}"
    assert {seed: round(value * n) for seed, value in null.items()} == {"1": 24, "2": 22, "3": 24}, null


_PLACE_DRIVER = r'''
import json, sys
from tools.cone_store import placement
order, kinds = json.loads(sys.argv[2]), json.loads(sys.argv[3])
found = placement.place(sys.argv[1], "tests/rd/", placement.Layer("pkg/rd", "pkg", ()), [tuple(t) for t in order],
                        {(a, b): frozenset(k) for a, b, k in kinds})
print("B5-PLACE " + json.dumps(found, sort_keys=True))
'''


def _place_in_process(root: Path, hash_seed: str) -> str:
    """The placement JSON computed in a fresh interpreter with ``PYTHONHASHSEED`` set: a set iterated into the report
    reads in another order under another hash seed, which a second call in this process cannot show. The order and
    the import kinds are handed over as JSON, so the child imports nothing of the test package (the fixture root,
    first on its path, holds a ``tests`` package of its own)."""
    import subprocess
    import sys

    from tests.gates.numerics.test_concept_order import import_kinds

    from ._cone_store_fixture import REPO, environment

    kinds = [[a, b, sorted(k)] for (a, b), k in sorted(import_kinds(root / "pkg/rd").items())]
    driver = root.parent / f"{root.name}.place_driver.py"
    driver.write_text(_PLACE_DRIVER)
    run = subprocess.run([sys.executable, "-O", str(driver), str(paths(root)["store"]), json.dumps(READS_ORDER),
                          json.dumps(kinds)], cwd=REPO, env=environment(root, PYTHONHASHSEED=hash_seed),
                         capture_output=True, text=True, timeout=600)
    line = next((l for l in run.stdout.splitlines() if l.startswith("B5-PLACE ")), None)
    if run.returncode != 0 or line is None:
        raise AssertionError(f"place driver exited {run.returncode}:\n{run.stdout[-2000:]}\n{run.stderr[-2000:]}")
    return line[len("B5-PLACE "):]


def test_v7_the_result_is_byte_identical_across_runs(captured: Path, placed: dict) -> None:
    """V7. A second placement gives the same JSON, in this process and in two fresh interpreters under two hash
    seeds. First reds: an unseeded shuffle of the null (the in-process leg); a set iterated into the report (the
    hash-seed legs; a two-element set coincides under one seed in two, so two seeds are read)."""
    here = json.dumps(placed, sort_keys=True)
    assert json.dumps(_place(captured), sort_keys=True) == here
    for seed in ("1", "2"):
        assert _place_in_process(captured, seed) == here, f"PYTHONHASHSEED={seed}"


def test_v8_the_union_reading_is_reported_beside_the_tiered_one(placed: dict) -> None:
    """V8 (C5, ruling 49). The summary reports ``homes_union``, ``union_moves`` (the placed tests whose union home is
    not their home; among the def- and statement-placed ones, by hand: ``test_input_value``, def ``low`` and
    union ``high``, and ``test_audit_run``, def ``mid`` and union ``high`` through its data reads) and ``statement_straddlers``;
    ``pure`` stays computed from defs alone. First red: union purity substituted for def purity."""
    summary, tests = placed["summary"], placed["tests"]
    moved = sorted(t for t, e in tests.items() if e.get("placed_by") in ("def", "statement")
                   and e["home_union"] != (e["home"] or "+".join(e["tie"] or ())))
    assert moved == [A + "test_audit_run", T + "test_input_value"], moved
    assert summary["union_moves"] >= len(moved)
    assert isinstance(summary["homes_union"], dict) and summary["homes_union"]
    assert T + "test_two_reads_apart" in summary["statement_straddlers"]
    assert tests[T + "test_input_value"]["pure_rt"] is True
    assert T + "test_input_value" not in summary["straddlers"]


def test_v9_an_above_layer_read_is_counted_apart(placed: dict) -> None:
    """V9. ``test_other`` reads ``pkg/other/o.py:OTHER`` and runs no ``pkg/other`` code: ``above_layer_statement``
    counts it, ``above_layer`` (defs) does not. First red: the two counts merged."""
    assert placed["summary"]["above_layer_statement"] == 1
    assert row(placed, "test_other")["above_layer"] == {}


def test_v10_an_unresolved_read_is_counted(placed: dict) -> None:
    """V10 (C11). The fixture ``late`` binds ``low.LATE`` with ``setattr``; the test reads it, and no statement binds
    it: ``unresolved == ["pkg.rd.low.LATE"]``. The test is still placed by its read of ``low``. First red: the counter
    dropped (a silent zero)."""
    entry = row(placed, "test_unresolved")
    assert entry["unresolved"] == ["pkg.rd.low.LATE"]
    assert entry["home"] == "low"
    assert all(not e.get("unresolved") for t, e in placed["tests"].items() if t != T + "test_unresolved")


def test_v11_a_store_of_b4s_form_is_refused_naming_the_own_reads(captured: Path, tmp_path: Path) -> None:
    """V11 (the archivist's reading, 2026-11-07). A store of form 3 (B4's: own cones, no ``own_reads``) is refused
    with a message naming the own reads, BEFORE it is opened: the copy still holds every test row afterwards (opening
    a store of another form drops its tables). First reds: the old fixed message (it names the own cones, which a
    form-3 store has); the refusal removed (the store is opened, emptied, and placed as zero tests)."""
    import sqlite3
    import shutil

    store = tmp_path / "b4_form.sqlite"
    shutil.copy2(paths(captured)["store"], store)
    with sqlite3.connect(store) as db:
        db.execute("UPDATE meta SET value = '3' WHERE key = 'format'")
        db.execute("DROP TABLE own_reads")
        before = db.execute("SELECT COUNT(*) FROM tests").fetchone()[0]
    assert before >= 30, before
    with pytest.raises(ValueError, match="own reads"):
        _place(captured, store=store)
    with sqlite3.connect(store) as db:
        assert db.execute("SELECT COUNT(*) FROM tests").fetchone()[0] == before


# ── group D, the data home ──


def test_d1_an_audit_of_the_layer_alone_homes_at_the_layer(placed: dict) -> None:
    """D1 (C13). ``test_audit_layer`` opens every file of ``pkg/rd`` and nothing else, runs no layer code, reads no
    binding: ``placed_by == "layer"``, no node home. First red: home at the highest-ranked file read (``high``)."""
    entry = row(placed, A + "test_audit_layer")
    assert entry["placed_by"] == "layer" and entry["home"] is None


def test_d2_an_audit_of_the_package_is_a_census_finding(placed: dict) -> None:
    """D2 (C13, ruled: listed apart). ``test_audit_package`` opens every file under ``pkg/``: ``placed_by ==
    "package"``. First red: ``layer``."""
    entry = row(placed, A + "test_audit_package")
    assert entry["placed_by"] == "package" and entry["home"] is None


def test_d3_data_never_outranks_a_def(placed: dict) -> None:
    """D3 (ruling 49). ``test_audit_run`` runs ``mid.g`` and opens every file of the layer: home ``mid`` by def; the
    data reads are listed under ``via``. First red: the data tier above the def tier (``high``, or ``layer``)."""
    entry = row(placed, A + "test_audit_run")
    assert entry["home"] == "mid" and entry["placed_by"] == "def"
    assert entry["via"]["high"] == ["data"] and "data" in entry["via"]["mid"]


# ── group A, the apparatus ──


def test_a1_a_premise_test_homes_at_its_apparatus(placed: dict) -> None:
    """A1 (ruling 48). ``test_premise`` reads only ``_G``, a constant of its own file built from no layer binding:
    placed by the apparatus, home the file. First red: test-tree bindings resolved to nothing (unplaced)."""
    entry = row(placed, "test_premise")
    assert entry["placed_by"] == "apparatus" and entry["home"] == TEST_RD
    assert entry["items"] == {TEST_RD: ["_G"]}


def test_a2_a_case_table_delivers_its_layer_node(placed: dict) -> None:
    """A2 (C6 as refined). ``test_case_table`` reads ``CASES = (low.Rec(1.0),)``: resolution continues through the
    test-tree binding to ``low:Rec`` (``via_apparatus``), so the home is ``low`` by statement, and the test rests on
    the apparatus item. First red: closure stopped at the apparatus (the test homes at its file)."""
    entry = row(placed, "test_case_table")
    assert entry["home"] == "low" and entry["placed_by"] == "statement"
    assert [TEST_RD, "CASES"] in entry["apparatus_items"]


def test_a3_resolution_stops_at_the_first_layer_binding(placed: dict) -> None:
    """A3 (C6). ``test_helper_table`` reads ``_helper.TABLE = (BUILT,)``, ``BUILT`` a ``mid`` constant built from
    ``low``: home ``mid``; ``low`` is not a node (its reads are ``closure``, pins only); the test rests on
    ``(_helper.py, TABLE)``. First red: closure followed through the layer binding (``low`` joins ``via``)."""
    entry = row(placed, "test_helper_table")
    assert entry["home"] == "mid" and "low" not in entry["via"]
    assert [HELPER, "TABLE"] in entry["apparatus_items"]


def test_a4_a_layer_constants_builder_is_provenance(placed: dict) -> None:
    """A4 (ruled side of attack a6). ``test_legendre_like`` reads ``mid.BUILT`` directly: ``via == {mid:
    [statement]}``, ``low`` (which built ``BUILT``) absent. First red: closure everywhere."""
    entry = row(placed, "test_legendre_like")
    assert entry["via"] == {"mid": ["statement"]} and entry["home"] == "mid"


def test_a5_a_helpers_reads_are_the_tests(placed: dict) -> None:
    """A5. ``test_helper_def`` calls ``_helper.make()``, which reads ``high.HIGH``; the test runs no layer def: home
    ``high`` by statement, and it rests on ``(_helper.py, make)``. First red: a helper's reads attributed to the
    helper (the test would home at the helper file)."""
    entry = row(placed, "test_helper_def")
    assert entry["home"] == "high" and entry["placed_by"] == "statement"
    assert [HELPER, "make"] in entry["apparatus_items"]


def test_a6_a_user_of_a_premise_is_proposed_to_rest_on_it(placed: dict) -> None:
    """A6 (C14). ``test_uses_premise`` reads ``_G`` and runs ``low.f``: the edge ``test_uses_premise ->
    test_premise`` is proposed at the node ``tests/rd/test_rd.py`` with ``c == 1.0``; ``test_global`` and
    ``test_case_table`` (which never read ``_G``) are not proposed. First red: containment over defs only (a premise
    test has no defs)."""
    edges = {(e["upper"], e["lower"]): e for e in placed["apparatus_supports"]}
    edge = edges.get((T + "test_uses_premise", T + "test_premise"))
    assert edge is not None and edge["c"] == 1.0 and edge["node"] == TEST_RD, placed["apparatus_supports"]
    assert (T + "test_global", T + "test_premise") not in edges
    assert (T + "test_case_table", T + "test_premise") not in edges


def test_a7_a_hook_is_no_tests_apparatus(placed: dict) -> None:
    """A7 (C4). The conftest's hook runs in every test: no test's ``items`` names it. The conftest's fixture ``box``
    is an apparatus item of ``test_instance_attr`` only, the one test requesting it. First red: C4 dropped."""
    holders = {t: e["items"].get("tests/rd/conftest.py") for t, e in placed["tests"].items()
               if e.get("items", {}).get("tests/rd/conftest.py")}
    assert holders == {T + "test_instance_attr": ["box"]}, holders


def test_a8_the_test_itself_is_not_its_item(placed: dict) -> None:
    """A8. ``test_global``'s items name no binding of its own file (its own def is the test, not apparatus). First
    red: the test's own statement read as an item."""
    assert TEST_RD not in row(placed, "test_global")["items"]


def test_a9_a_declared_edge_onto_a_premise_is_measured_over_items(captured: Path) -> None:
    """A9 (AC10 in the fixture). ``rests_on`` declared from ``test_uses_premise`` onto ``test_premise``: the row is
    measured, ``c_items == 1.0`` (``_G`` is contained), where the defs-only measure has nothing to read (the lower
    test ran no def). First red: declared edges measured over defs only (the row is missing)."""
    declared = ((T + "test_uses_premise", T + "test_premise"),)
    found = [r for r in _place(captured, declared)["declared_rests_on"]
             if (r["upper"], r["lower"]) == (T + "test_uses_premise", T + "test_premise")]
    assert len(found) == 1 and found[0]["c_items"] == 1.0, found
