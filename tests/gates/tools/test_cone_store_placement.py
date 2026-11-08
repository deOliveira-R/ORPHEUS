"""``place``: a layer's tests homed against the layer's concept order, read from each test's OWN cone (B3).

Spec ``scratch/verification_tree/cone_store/b3/spec_b3.md`` §2.3-§2.5 and §4, rows PL1-PL11, and the
orchestrator's PL13-PL15 (PL12, the capture instrument counted apart, retired with the instrument's exception
2026-11-03: the recorder no longer records its own apparatus, ``_APPARATUS`` in ``_traced_memo_boot``). The fixture layer is
``pkg/lay`` with the order :data:`PLACEMENT_ORDER` (``_cone_store_fixture.py``): ``beta`` imports ``zeta`` and ``s``
at top level, ``alpha`` imports ``beta``, ``gamma`` imports ``zeta`` lazily and ``beta`` under ``TYPE_CHECKING``,
``beta2`` imports nothing; ``pkg/up`` is the upstream layer, ``pkg/top`` a higher one; ``tests/lay/conftest.py`` runs
``alpha.fd`` (hence ``fb``, ``fa``, ``_fa_inner``) at import, so those defs are in every test's inherited floor. Every
expectation below is read off those sources by hand.

The interface is ``tools.cone_store.placement.place(store_path, layer_prefix, layer, order, kinds, *, threshold, min_defs,
seeds)``, ``layer`` a ``placement.Layer(root, package, upstream_layers)`` as the orchestrator wrote it (2026-11-03). ``result["tests"][test_id]`` holds
``upstream_layers`` (upstream module to def count) and ``above_layer`` (package to def count), ``helpers``, ``elsewhere`` and
``unplaced`` (``None`` or a reason starting ``no layer read``, ``not passed`` or ``cone unknown``), and, when placed,
``home`` (a node, ``None`` for a tie), ``tie`` (the tied nodes or ``None``), ``pure_rt``, ``pure_lazy``, ``pure_typ``
and ``foreign`` (node to ``lazy``, ``typing``, ``none`` or ``tie``); ``result["supports"]`` is a list of ``{"upper",
"lower", "c"}`` over test functions (a parametrized function's cases are one node, its id without the brackets);
``result["duplicates"]`` a list of ``{"a", "b", ...}``; ``result["summary"]["null_pure_fraction"]`` the shuffled
null per seed.
"""
from __future__ import annotations

import json
import shutil
import sqlite3
from pathlib import Path

import pytest

from ._cone_store_fixture import PLACEMENT_ORDER, SUBJECT_ORDER, build, capture, paths

pytestmark = pytest.mark.foundation

PREFIX = "tests/lay/"
CODE_ROOT = "pkg/lay"
UPSTREAM = ("pkg/up",)
T = "tests/lay/test_pl.py::"


@pytest.fixture(scope="module")
def captured(tmp_path_factory) -> Path:
    root = build(tmp_path_factory.mktemp("placement") / "tree", memo=False, placement=True)
    capture(root, "tests/lay")
    return root


def _place(root: Path, store: Path | None = None, declared: tuple[tuple[str, str], ...] = ()) -> dict:
    from tools.cone_store import placement

    from tests.gates.numerics.test_concept_order import import_kinds

    layer = placement.Layer(CODE_ROOT, "pkg", UPSTREAM)
    return placement.place(store if store is not None else paths(root)["store"], PREFIX, layer, PLACEMENT_ORDER,
                           import_kinds(root / CODE_ROOT), declared=declared)


@pytest.fixture(scope="module")
def placed(captured: Path) -> dict:
    return _place(captured)


def row(placed: dict, name: str) -> dict:
    return placed["tests"][T + name]


def test_pl1_the_home_is_the_maximal_rank_node(placed: dict) -> None:
    """PL1. ``test_d_only`` runs ``alpha``, ``beta``, ``zeta`` (two defs in ``zeta``): its home is ``alpha``, the top
    of the order, neither the alphabetical extreme nor the module with most defs. ``test_low`` is homed at ``zeta``.
    First red: a home by name or by def count."""
    assert row(placed, "test_d_only")["home"] == "alpha"
    assert row(placed, "test_low")["home"] == "zeta"


def test_pl2_purity_is_closure_under_the_runtime_import_graph(placed: dict) -> None:
    """PL2. ``test_d_only`` is pure (``beta`` and ``zeta`` are ancestors of ``alpha`` by top-level imports).
    ``test_none`` adds ``beta2``, below ``alpha`` in the order but imported by nothing ``alpha`` imports: a straddler
    of kind ``none``. First red: closure by rank order (which calls ``test_none`` pure)."""
    assert row(placed, "test_d_only")["pure_rt"] is True and row(placed, "test_d_only")["foreign"] == {}
    none = row(placed, "test_none")
    assert none["home"] == "alpha" and none["pure_rt"] is False
    assert none["foreign"] == {"beta2": "none"}


def test_pl3_each_omitted_edge_kind_explains_its_foreign_node(placed: dict) -> None:
    """PL3. ``test_lazy_typ`` runs ``gamma``, ``zeta`` (reached by ``gamma``'s lazy import) and ``beta`` (named by
    ``gamma`` under ``TYPE_CHECKING``): kinds ``lazy`` and ``typing``, and the three closures read impure, impure,
    pure. First red: a kind walk that does not enter defs or does not read ``TYPE_CHECKING``."""
    lazy_typ = row(placed, "test_lazy_typ")
    assert lazy_typ["home"] == "gamma"
    assert lazy_typ["foreign"] == {"zeta": "lazy", "beta": "typing"}
    assert (lazy_typ["pure_rt"], lazy_typ["pure_lazy"], lazy_typ["pure_typ"]) == (False, False, True)


def test_pl4_two_maximal_nodes_of_one_tier_are_a_tie(placed: dict) -> None:
    """PL4. ``test_tie`` runs ``beta`` and ``beta2``, one tier: both are its home, flagged, never broken silently.
    First red: a tie broken by name or by order of appearance."""
    tie = row(placed, "test_tie")
    assert tie["home"] is None and tie["tie"] == ["beta", "beta2"]
    assert tie["pure_rt"] is False and tie["foreign"].get("beta") == "tie"


def test_pl5_the_inherited_floor_does_not_move_a_home(placed: dict) -> None:
    """PL5. The conftest ran ``alpha.fd`` in every test's process; ``test_low`` runs only ``zeta``. Its home is
    ``zeta``. First red: placement reading the selection cone (every home becomes ``alpha``)."""
    assert row(placed, "test_low")["home"] == "zeta"
    assert row(placed, "test_low")["pure_rt"] is True


def test_pl6_a_def_both_the_floor_and_the_test_ran_is_the_tests(placed: dict) -> None:
    """PL6. ``test_d_only`` runs ``fd``, ``fb``, ``fa``, all also in the floor: it is placed, at ``alpha``. First red:
    placement by subtracting the floor (the test would be left with no layer def)."""
    d_only = row(placed, "test_d_only")
    assert d_only["unplaced"] is None and d_only["home"] == "alpha"


def test_pl7_other_layers_are_counted_never_foreign(placed: dict) -> None:
    """PL7. ``test_layers`` runs ``zeta``, ``pkg/up/u.py`` (an upstream layer) and ``pkg/top/t.py`` (a higher one):
    home ``zeta``, pure, the first counted under ``upstream_layers`` by its module (spec §2.5), the second under
    ``above_layer``. First reds: either read as a foreign node of the layer; the upstream count kept per layer root."""
    layers = row(placed, "test_layers")
    assert layers["home"] == "zeta" and layers["pure_rt"] is True and layers["foreign"] == {}
    assert set(layers["upstream_layers"]) == {"pkg/up/u"}, layers["upstream_layers"]
    assert any("pkg/top" in key for key in layers["above_layer"]), layers["above_layer"]
    assert not any("pkg/top" in key for key in layers["upstream_layers"])


@pytest.mark.parametrize(("name", "reason"), [("test_nodef", "no layer read"), ("test_fails", "not passed"),
                                              ("test_spawn", "cone unknown")])
def test_pl8_an_unplaced_test_is_listed_with_its_reason(placed: dict, name: str, reason: str) -> None:
    """PL8. A test that read nothing of the layer or of the apparatus (B5: a test running no layer def but reading a
    binding is placed by that read), a test that failed and a test that started a process are kept in the result,
    unplaced, with their reason. First red: dropping them (a ``KeyError``), or placing a failed test's partial cone."""
    unplaced = row(placed, name)
    assert unplaced.get("home") is None
    assert unplaced["unplaced"] is not None and unplaced["unplaced"].startswith(reason), unplaced["unplaced"]


def _supports(placed: dict) -> set[tuple[str, str]]:
    return {(s["upper"], s["lower"]) for s in placed["supports"]}


def test_pl9_supports_by_fractional_containment(placed: dict) -> None:
    """PL9. Own ``pkg`` defs: ``test_s_lower`` 20 (``f0``-``f19``; two parametrized cases, one node);
    ``test_s_upper95`` holds 19 of them (and ``fb``, ``fa``, ``_fa_inner``), so ``c`` = 0.95, proposed; ``test_s_upper90``
    holds 18 and nothing else, ``c`` = 0.9, not (it is no superset of ``test_s_upper95`` either: a fixture where the
    two contained each other at 0.95 makes a 2-cycle, and the reduction then drops ``test_s_super``'s edge too);
    ``test_s_super`` holds all of both, but its edge to ``test_s_lower`` is implied through ``test_s_upper95`` and
    reduced away; ``test_tiny``'s cone has one def and is nobody's lower; ``test_dup1`` and ``test_dup2`` have equal
    cones, a duplicate candidate. First reds: ``>`` for ``>=``; no collapse of cases; no transitive reduction; no
    minimum cone size; equal cones not grouped."""
    lower, up95, up90, sup = (T + n for n in ("test_s_lower", "test_s_upper95", "test_s_upper90", "test_s_super"))
    found = _supports(placed)
    assert (up95, lower) in found
    assert next(s["c"] for s in placed["supports"] if (s["upper"], s["lower"]) == (up95, lower)) == pytest.approx(0.95)
    assert (up90, lower) not in found
    assert (sup, up95) in found and (sup, lower) not in found
    assert not any(l == T + "test_tiny" for _, l in found)
    assert not any("[" in u or "[" in l for u, l in found)
    assert any({d["a"], d["b"]} == {T + "test_dup1", T + "test_dup2"} for d in placed["duplicates"])


def test_pl10_the_result_is_byte_identical_across_runs(captured: Path, placed: dict) -> None:
    """PL10. A second run gives the same JSON, the shuffled null included (three seeded draws). First red: an unseeded
    shuffle, or an unordered set serialised."""
    again = _place(captured)
    assert json.dumps(again, sort_keys=True) == json.dumps(placed, sort_keys=True)
    assert len(placed["summary"]["null_pure_fraction"]) == 3


def test_pl16_the_shuffled_null_falls_below_the_measured_purity(placed: dict) -> None:
    """PL16 (qa's F3 on B3, 2026-11-03). On the fixture the placement carries information: each seeded draw of the
    null reads a smaller pure fraction than the measured one, the draws differ from one another, and both read on one
    denominator, the placed tests. A RECORD of the fixture with its seeds: ``[M]`` 2026-11-03 the measured fraction is
    9/14 and the null 5/14, 4/14 and 7/14 (seeds 1, 2, 3). First red ``[M]``: the shuffle replaced by the identity
    (null = measured, 9/14 three times); one seed for every draw (the draws equal)."""
    summary = placed["summary"]
    measured, null, n = summary["pure_fraction"], summary["null_pure_fraction"], summary["placed"]
    assert n == 14 and len(null) == 3, f"placed {n}, draws {sorted(null)}: the fixture changed"
    for fraction in (measured, *null.values()):
        assert abs(fraction * n - round(fraction * n)) < 1e-4, f"{fraction} is not a fraction of the {n} placed tests"
    assert round(measured * n) == 9, f"measured {measured}"
    above = {seed: value for seed, value in null.items() if value >= measured}
    assert not above, f"the null reads at or above the measured {measured}: {above}"
    assert len(set(null.values())) == len(null), f"the seeded draws coincide: {null}"


def test_pl11_a_store_without_own_cones_is_refused(captured: Path, tmp_path: Path) -> None:
    """PL11. A store of the previous form (B2: no own cones) is refused with a message naming the own cones, before
    it is opened (opening empties it) and never read through the selection cone. First red: a silent fallback (the
    store opened, emptied, and placed as zero tests, or read raw)."""
    store = tmp_path / "b2_form.sqlite"
    shutil.copy2(paths(captured)["store"], store)
    with sqlite3.connect(store) as db:
        db.execute("UPDATE meta SET value = '2' WHERE key = 'format'")
    with pytest.raises(ValueError, match="own cones"):
        _place(captured, store)


def test_pl13_an_import_cycle_closes_exactly_whatever_the_visiting_order() -> None:
    """PL13 (the elegance review of B3, 2026-11-03). ``ancestors`` on a cycle: ``a`` and ``b`` import each other
    and ``c`` imports ``a``, so ``c``'s ancestors are ``{a, b}`` and each of ``a``, ``b`` has the other. The
    memoised walk it replaces answered ``{b, c}`` or ``{}`` by visiting order. The live witness: the numerics
    graph's ``traced_memo`` <-> ``_traced_memo_boot`` cycle, where the top-level import must be in the widest
    closure of ``traced_memo``."""
    from tests.gates.numerics.test_concept_order import import_kinds
    from tools.cone_store import placement

    kinds = {("a", "b"): frozenset({"runtime"}), ("b", "a"): frozenset({"lazy"}), ("a", "c"): frozenset({"runtime"})}
    for nodes in (["a", "b", "c"], ["c", "b", "a"], ["b", "c", "a"]):
        widest = placement.ancestors(kinds, nodes, ("runtime", "lazy"))
        assert widest == {"a": frozenset({"b"}), "b": frozenset({"a"}), "c": frozenset({"a", "b"})}, nodes
    live = placement.ancestors(import_kinds(), ["traced_memo", "_traced_memo_boot"], ("runtime", "lazy", "typing"))
    assert "_traced_memo_boot" in live["traced_memo"] and "traced_memo" in live["_traced_memo_boot"]


def test_pl14_an_import_kind_with_no_closure_is_refused(tmp_path: Path) -> None:
    """PL14. A kind the concept order's walk gains and this placement has no closure for is refused by name, never
    read as ``none`` (the kinds are spelled in the gate and in ``placement.CLOSURES``)."""
    from tools.cone_store import placement

    with pytest.raises(ValueError, match="import kinds this placement has no closure for"):
        placement.place(tmp_path / "absent.sqlite", PREFIX, placement.Layer(CODE_ROOT), PLACEMENT_ORDER,
                        {("alpha", "beta"): frozenset({"conditional"})})


def test_pl15_each_declared_rests_on_edge_is_measured_on_function_cones(captured: Path) -> None:
    """PL15 (spec §2.4, the threshold re-judged on the declared population). Declared edges, as pytest resolves them:
    ``test_s_upper90`` over ``test_s_lower`` holds 18 of its 20 defs, ``c`` = 0.9; ``test_s_super`` over the CASE
    ``test_s_lower[1]`` collapses to the function, ``c`` = 1.0; the case ``test_s_lower[2]`` over ``test_tiny``
    (one def, ``fa``, none of the ``s`` defs) is ``c`` = 0, and the same edge declared from the other case is one
    edge; a target no stored test has is counted, not measured. First reds: no case collapse (the ``[1]`` target
    unmatched, the two ``test_s_lower`` cases two edges); ``>`` for ``>=`` in the counts; an unmatched target given a
    ``c``."""
    declared = ((T + "test_s_upper90", T + "test_s_lower"), (T + "test_s_super", T + "test_s_lower[1]"),
                (T + "test_s_lower[1]", T + "test_tiny"), (T + "test_s_lower[2]", T + "test_tiny"),
                (T + "test_low", T + "test_nowhere"))
    found = _place(captured, declared=declared)
    measured = {(row["upper"], row["lower"]): row["c"] for row in found["declared_rests_on"]}
    assert measured == {(T + "test_s_upper90", T + "test_s_lower"): pytest.approx(0.9),
                        (T + "test_s_super", T + "test_s_lower"): 1.0,
                        (T + "test_s_lower", T + "test_tiny"): 0.0}
    assert found["declared_counts"] == {"1.0": 1, "0.95": 1, "0.9": 2, "0.8": 2, "below": 1}
    assert found["declared_unmeasured"] == {"target not judged": 1}


# ── ruling 50: a support's lower SUBJECT is contained too ──
#
# The fixture is ``_SUBJECT`` (``build(root, subject=True)``): layer ``pkg/sj``, order ``a, b < mid, mid2 < top``; every
# containment below is counted off its sources by hand (``_cone_store_fixture.py``, the comment above ``SUBJECT_ORDER``).

SJ = "tests/sj/test_sj.py::"


@pytest.fixture(scope="module")
def subject_captured(tmp_path_factory) -> Path:
    root = build(tmp_path_factory.mktemp("subject") / "tree", memo=False, subject=True)
    capture(root, "tests/sj")
    return root


def _place_subject(root: Path, declared: tuple[tuple[str, str], ...] = (), threshold: float = 0.95) -> dict:
    from tools.cone_store import placement

    from tests.gates.numerics.test_concept_order import import_kinds

    return placement.place(paths(root)["store"], "tests/sj/", placement.Layer("pkg/sj", "pkg", ()), SUBJECT_ORDER,
                           import_kinds(root / "pkg/sj"), declared=declared, threshold=threshold)


@pytest.fixture(scope="module")
def subject_placed(subject_captured: Path) -> dict:
    return _place_subject(subject_captured)


def _proposed(placed: dict) -> dict[tuple[str, str], dict]:
    """Every proposed support, along the order or against it, by (upper, lower)."""
    return {(r["upper"], r["lower"]): r for r in [*placed["supports"], *placed["supports_against_order"]]}


def test_pl17_an_upper_that_never_runs_the_lowers_subject_is_no_support(subject_placed: dict) -> None:
    """PL17 (ruling 50). ``test_upper_misses`` contains ``test_lower``'s whole cone at 40/41 (0.976, above 0.95) and
    none of its subject (``mid.m0``, the lower's only def at its home ``mid``): no support is proposed, and the pair
    is absent from the counts at 0.95 (3: scenarios A, B and C's full containments). First red: the cone's
    containment alone (the predicate before ruling 50)."""
    assert subject_placed["tests"][SJ + "test_lower"]["home"] == "mid"
    assert (SJ + "test_upper_misses", SJ + "test_lower") not in _proposed(subject_placed)
    assert subject_placed["support_counts"]["0.95"] == 3, subject_placed["support_counts"]


def test_pl18_an_upper_that_runs_the_subject_is_a_support_with_its_c_home(subject_placed: dict) -> None:
    """PL18 (ruling 50). ``test_upper_runs`` contains ``test_lower``'s cone at 41/41 and its subject at 1/1: proposed,
    with ``c == 1.0`` and ``c_home == 1.0`` on the record. First red: ``c_home`` missing from the record; a subject
    read from the upper's home instead of the lower's (the upper's ``top`` defs are not in the lower: 0, no support)."""
    record = _proposed(subject_placed).get((SJ + "test_upper_runs", SJ + "test_lower"))
    assert record is not None, sorted(_proposed(subject_placed))
    assert record["c"] == 1.0 and record["c_home"] == 1.0


def test_pl19_a_declared_edge_carries_c_home_and_none_for_a_lower_not_placed_by_defs(subject_captured: Path) -> None:
    """PL19 (ruling 50: declared rows carry ``c_home`` and are not filtered). Declared ``test_upper_misses`` over
    ``test_lower``: kept, ``c`` = 40/41 and ``c_home`` = 0.0; ``test_upper_runs`` over ``test_lower``: ``c_home`` = 1.0;
    ``test_upper_misses`` over ``test_reads_a_constant`` (placed by statement, no def): ``c`` and ``c_home`` are
    ``None`` (measured by ``c_items`` alone). First reds: the subject's containment replaced (``c_home`` 1.0 on the
    first row); a ``c_home`` given to a lower with no def home (a value, or a division by an empty subject)."""
    declared = ((SJ + "test_upper_misses", SJ + "test_lower"), (SJ + "test_upper_runs", SJ + "test_lower"),
                (SJ + "test_upper_misses", SJ + "test_reads_a_constant"))
    rows = {(r["upper"], r["lower"]): r for r in _place_subject(subject_captured, declared)["declared_rests_on"]}
    misses = rows[(SJ + "test_upper_misses", SJ + "test_lower")]
    assert misses["c"] == pytest.approx(40 / 41, abs=1e-6), misses
    assert misses["c_home"] == 0.0, misses
    assert rows[(SJ + "test_upper_runs", SJ + "test_lower")]["c_home"] == 1.0
    constant = rows[(SJ + "test_upper_misses", SJ + "test_reads_a_constant")]
    assert constant["c"] is None and constant["c_home"] is None, constant


def test_pl20_a_tied_lowers_subject_spans_every_tied_node(subject_placed: dict) -> None:
    """PL20 (ruling 50, a tie). ``test_tie_lower`` is tied at ``mid`` and ``mid2``; its subject is ``{p0, n0}``.
    ``test_tie_partial`` runs ``p0`` and not ``n0``: cone 41/42, subject 1/2, no support. ``test_tie_full`` runs both:
    proposed, ``c_home == 1.0``. First red: the subject read at one tied node only (``mid``: ``test_tie_partial``
    then holds it whole and is proposed)."""
    lower = subject_placed["tests"][SJ + "test_tie_lower"]
    assert lower["home"] is None and lower["tie"] == ["mid", "mid2"]
    found = _proposed(subject_placed)
    assert (SJ + "test_tie_partial", SJ + "test_tie_lower") not in found
    assert found[(SJ + "test_tie_full", SJ + "test_tie_lower")]["c_home"] == 1.0


def test_pl21_a_support_against_the_order_carries_its_c_home(subject_placed: dict) -> None:
    """PL21 (qa F6). ``test_against_upper`` homes at ``side``, which imports nothing, and contains
    ``test_against_lower`` (home ``mid``) whole: the support is listed against the order, with ``c_home == 1.0`` on
    its record. First red: ``c_home`` written on the along-order records only."""
    against = {(r["upper"], r["lower"]): r for r in subject_placed["supports_against_order"]}
    record = against.get((SJ + "test_against_upper", SJ + "test_against_lower"))
    assert record is not None, sorted(against)
    assert record["c_home"] == 1.0, record


def test_pl22_a_declared_lower_with_a_cone_and_no_def_home_has_c_and_no_c_home(subject_captured: Path) -> None:
    """PL22 (qa F6). ``test_reads_and_runs_above`` is placed by statement and ran one def ABOVE the layer
    (``pkg/sjx/o.o0``), so declared under ``test_against_upper`` it has a cone, ``c == 0.0``, and no subject,
    ``c_home is None``. First red: a ``c_home`` defaulted to ``c`` for a lower with no def home."""
    declared = ((SJ + "test_against_upper", SJ + "test_reads_and_runs_above"),)
    found = _place_subject(subject_captured, declared)
    assert found["tests"][SJ + "test_reads_and_runs_above"]["placed_by"] == "statement"
    (row,) = found["declared_rests_on"]
    assert row["c"] == 0.0, row
    assert row["c_home"] is None, row


def test_pl23_the_subject_is_applied_at_every_counted_threshold(subject_placed: dict) -> None:
    """PL23 (qa F6). The counts at :data:`THRESHOLDS`, by hand: at 1.0 and 0.95 the full containments of scenarios A,
    B and C (3); at 0.9 and 0.8 also ``test_edge_upper`` over ``test_edge_lower`` (cone 0.98, subject 0.9 exactly);
    every other pair holding at 0.8 has a subject below 0.8 (scenario E's 2/3, B's 1/2, A's 0). First red: the
    subject compared with the proposal threshold (0.95) instead of the count's own (0.9 and 0.8 read 3)."""
    assert subject_placed["support_counts"] == {"1.0": 3, "0.95": 3, "0.9": 4, "0.8": 4}


def test_pl24_a_subject_contained_at_exactly_the_threshold_holds(subject_captured: Path) -> None:
    """PL24 (qa F6, the inclusive ``>=``). Placed at threshold 0.9, ``test_edge_upper`` over ``test_edge_lower``
    (subject 9/10, exactly 0.9 in binary as in the literal) is proposed with ``c_home == 0.9``. First red: ``>`` for
    ``>=`` on the subject."""
    found = _place_subject(subject_captured, threshold=0.9)
    record = _proposed(found).get((SJ + "test_edge_upper", SJ + "test_edge_lower"))
    assert record is not None, sorted(_proposed(found))
    assert record["c_home"] == 0.9, record


def test_pl25_a_functions_subject_spans_the_homes_of_all_its_cases(subject_captured: Path) -> None:
    """PL25 (qa F5, ruling 4 with ruling 50). ``test_two_homes``'s cases home at ``mid`` and ``top``; the function's
    subject is its defs at both, ``{v0, v1, h0}``, where ``v1`` (run by the case homed at ``top``) is at the other
    case's home. ``test_two_homes_upper`` runs ``v0`` and ``h0`` and not ``v1``: cone 42/43, subject 2/3, no support,
    and a declared row reads ``c_home`` 2/3. First red: each case's subject read at its own home (``{v0, h0}``: 1.0,
    a support)."""
    declared = ((SJ + "test_two_homes_upper", SJ + "test_two_homes"),)
    found = _place_subject(subject_captured, declared)
    assert {found["tests"][SJ + f"test_two_homes[{k}]"]["home"] for k in (1, 2)} == {"mid", "top"}
    assert (SJ + "test_two_homes_upper", SJ + "test_two_homes") not in _proposed(found)
    (row,) = found["declared_rests_on"]
    assert row["c_home"] == pytest.approx(2 / 3, abs=1e-6), row
