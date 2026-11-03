"""``place``: each test of a layer placed in the layer's concept order, read from its OWN cone (B3, ruling 22).

The layer's code is cut into NODES at the concept order's grain (``tests/gates/numerics/test_concept_order.py``: a
module directly under the layer's root, a package directly under it being one node), and the order is a sequence of
tiers, upstream first. For a test ``t``, ``own(t)`` is the set of nodes holding a def its own cone ran
(:meth:`ConeStore.own_pins`, what pytest's protocol for that test entered; the conftests' floor is not in it unless
the test runs it too).

* The HOME is the node of ``own(t)`` of highest rank; two nodes of one tier at that rank are a TIE, reported as the
  tie set and never broken.
* The test is PURE iff every other node of ``own(t)`` is an ancestor of the home in the RUNTIME import graph (a
  top-level import, at any depth). The closure is the import graph, never the rank order: under the rank order every
  lower node is upstream, so purity would only mean "no tie", which a shuffled null reads the same.
* A foreign node (in ``own(t)``, not the home, not upstream of it at runtime) carries the KIND of the import edge
  that first explains it: ``lazy`` (upstream once imports inside a def count), ``typing`` (once ``TYPE_CHECKING``
  imports count too), ``none`` (no import path: the test itself composes the two), or ``tie``. Purity is reported
  under the three cumulative closures side by side (:data:`CLOSURES`), so the difference between them counts the
  straddles each omitted edge kind explains.
* Code of the layers below (``upstream_layers``) is counted per test by its module, never a straddle; code of any other
  first-party package is ABOVE-LAYER, a finding of its own; test helpers are counted, and so is code in no region.
* The memo's recorder (``orpheus/numerics/_traced_memo_boot.py``) is an ordinary node: the code its hooks run on a
  run's events is its declared apparatus, which it never records (``_APPARATUS``), so a test's own cone holds a
  recorder def only when the test ran one.
* A test is UNPLACED, with its reason, when its cone is unknown, it did not pass, or its own cone holds no layer def.

Proposed ``executes`` supports (S17: containment is fractional): over the test FUNCTIONS (a parametrised test's
cases are one claim, ruling 4; their own defs are united), ``c(U ⊇ L) = |D(L) ∩ D(U)| / |D(L)|`` over own defs
under the first-party package. ``U -> L`` is proposed when ``c >= threshold``, ``|D(L)| >= min_defs`` and the two
cones are not equal (equal cones of at least ``min_defs`` defs are DUPLICATE CANDIDATES, with no verdict: ruling 21
decides on kill sets; two smaller equal cones are below the size a containment is read at).
Fractional containment can cycle, so functions that contain each other at the threshold are a NEAR-DUPLICATE group,
also listed with no verdict, and the relation is transitively reduced over its condensation (each group one node).
A pair whose lower home is neither the upper home nor upstream of it under the widest closure is listed apart,
against the order, for the orchestrator's reading (S18).

The declared ``rests_on`` edges (``declared``, each (upper test id, target id) as pytest resolved the marker:
``tests/_harness/rests_on_manifest.py``, the one resolver) are measured on the same function cones: each target
collapses to its functions (a case to its function, a class or a file to every function under it), and each edge
reports ``c(U ⊇ L)``, with the cumulative counts at :data:`THRESHOLDS` and the count below the lowest, so the
threshold is judged on the layer's whole declared population. An edge whose upper or target holds no stored function
that passed with a known cone, or whose lower ran no def, is counted by that reason, never given a ``c``.

The null shuffles the node label of each distinct (node, def) across the layer with a fixed seed and re-reads purity;
a placement that the null matches carries no information (S13). The output is a plain dict, byte-identical across runs
when dumped with sorted keys.
"""
from __future__ import annotations

import random
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Any

import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import breadth_first_order, connected_components

from orpheus.numerics import traced_memo as memo

from tools.cone_store.store import FORMAT, ConeStore, stored_format

#: The import edge kinds, each with the output key of the cumulative closure that first admits it. A closure admits
#: its kind and every kind before it: ``runtime`` reads top-level imports; ``lazy`` adds imports inside a def;
#: ``typing`` adds imports under ``TYPE_CHECKING``. The graph is the concept order's (``import_kinds`` beside
#: ``TIERS`` in ``tests/gates/numerics/test_concept_order.py``), passed in: this package parses no source
#: (``test_cone_store_recorder.py`` R7), and one walk serves the gate and the placement. A kind outside this table
#: is refused (:func:`place`), so a kind the gate gains cannot silently read as ``none``.
CLOSURES: dict[str, str] = {"runtime": "pure_rt", "lazy": "pure_lazy", "typing": "pure_typ"}
#: The containment thresholds the report counts proposals at, beside the one it proposes at.
THRESHOLDS = (1.0, 0.95, 0.9, 0.8)


# ── where a def is ──


@dataclass(frozen=True)
class InLayer:
    node: str


@dataclass(frozen=True)
class Upstream:
    """A def of a layer below, by that layer's node (its root itself for the layer's own ``__init__``)."""

    module: str


@dataclass(frozen=True)
class Above:
    package: str


@dataclass(frozen=True)
class Helper:
    pass


@dataclass(frozen=True)
class Elsewhere:
    """First-party code in no region above: the layer's own ``__init__`` (no node), a file outside the package."""


Region = InLayer | Upstream | Above | Helper | Elsewhere


@dataclass(frozen=True)
class Layer:
    """Where a layer's code is, and what surrounds it: ``root`` (``orpheus/numerics``, relative to the repository),
    the package it belongs to, and the roots of the layers below it."""

    root: str
    package: str = "orpheus"
    upstream_layers: tuple[str, ...] = ("orpheus/geometry", "orpheus/values")

    def node_of(self, path: str) -> str | None:
        """The node of a layer file (a module or a package directly under the root); ``None`` for a file outside
        the layer and for the layer's own ``__init__``."""
        return _node_under(self.root, path)

    def region(self, path: str) -> Region:
        if (node := self.node_of(path)) is not None:
            return InLayer(node)
        for upstream in self.upstream_layers:
            if path.startswith(upstream.rstrip("/") + "/"):
                node = _node_under(upstream, path)
                return Upstream(upstream if node is None else f"{upstream.rstrip('/')}/{node}")
        parts = PurePosixPath(path).parts
        if parts[:1] == ("tests",):
            return Helper()
        if parts[:1] == (self.package,) and not path.startswith(self.root.rstrip("/") + "/"):
            return Above("/".join(parts[:2]))
        return Elsewhere()


def _node_under(root: str, path: str) -> str | None:
    """The node of ``path`` under ``root`` at the concept order's grain (a module or a package directly under the
    root); ``None`` for a file outside the root and for the root's own ``__init__``."""
    parts, head = PurePosixPath(path).parts, PurePosixPath(root).parts
    if parts[: len(head)] != head or len(parts) <= len(head):
        return None
    name = parts[len(head)].removesuffix(".py")
    return None if name == "__init__" else name


#: The numerics layer, the census's first (ruling 5).
NUMERICS = Layer("orpheus/numerics")


@dataclass(frozen=True)
class OwnDefs:
    """A test's own defs, sorted by region: the layer's per node, and the counts elsewhere."""

    layer: dict[str, frozenset[str]]
    upstream: Counter[str] = field(default_factory=Counter)
    above: Counter[str] = field(default_factory=Counter)
    helpers: int = 0
    elsewhere: int = 0
    first_party: frozenset[tuple[str, str]] = frozenset()


def own_defs(pins: Iterable[Any], layer: Layer) -> OwnDefs:
    by_node: dict[str, set[str]] = defaultdict(set)
    upstream: Counter[str] = Counter()
    above: Counter[str] = Counter()
    helpers = elsewhere = 0
    first_party: set[tuple[str, str]] = set()
    for pin in pins:
        if not isinstance(pin, memo.DefPin):
            continue
        match layer.region(pin.path):
            case InLayer(node):
                by_node[node].add(pin.qualname)
                first_party.add((pin.path, pin.qualname))
            case Upstream(module):
                upstream[module] += 1
                first_party.add((pin.path, pin.qualname))
            case Above(package):
                above[package] += 1
                first_party.add((pin.path, pin.qualname))
            case Helper():
                helpers += 1
            case Elsewhere():
                elsewhere += 1
    return OwnDefs({node: frozenset(q) for node, q in by_node.items()}, upstream, above, helpers, elsewhere,
                   frozenset(first_party))


# ── graphs ──


def _reach(successors: Mapping[str, Iterable[str]], nodes: Sequence[str]) -> dict[str, frozenset[str]]:
    """Each node's reachable set in a directed graph, the node itself excluded: one breadth-first order per node
    over a sparse adjacency, so a cycle is walked once, whatever the visiting order."""
    index = {node: i for i, node in enumerate(nodes)}
    rows, cols = [], []
    for node, targets in successors.items():
        for target in targets:
            rows.append(index[node])
            cols.append(index[target])
    graph = csr_array((np.ones(len(rows)), (rows, cols)), shape=(len(nodes), len(nodes)))
    reached = {}
    for node, i in index.items():
        order = breadth_first_order(graph, i, directed=True, return_predecessors=False)
        reached[node] = frozenset(nodes[j] for j in order if j != i)
    return reached


def ancestors(kinds: Mapping[tuple[str, str], frozenset[str]], nodes: Iterable[str],
              admitted: Iterable[str]) -> dict[str, frozenset[str]]:
    """Each node's ancestors in the import graph whose edges carry a kind in ``admitted``: reachability over the
    reversed graph (importer to imported), so an import cycle is closed exactly."""
    admitted = frozenset(admitted)
    up: dict[str, set[str]] = defaultdict(set)
    for (imported, importer), found in kinds.items():
        if found & admitted:
            up[importer].add(imported)
    ordered = sorted(set(nodes) | {n for edge in kinds for n in edge})
    reached = _reach(up, ordered)
    return {node: reached[node] for node in nodes}


def _groups(relation: set[tuple[str, str]]) -> dict[str, str]:
    """Each function's strongly connected group in ``relation``, named by its least member."""
    nodes = sorted({n for pair in relation for n in pair})
    if not nodes:
        return {}
    index = {node: i for i, node in enumerate(nodes)}
    rows = [index[u] for u, _ in relation]
    cols = [index[l] for _, l in relation]
    graph = csr_array((np.ones(len(rows)), (rows, cols)), shape=(len(nodes), len(nodes)))
    _, labels = connected_components(graph, directed=True, connection="strong")
    least: dict[int, str] = {}
    for node, label in zip(nodes, labels):
        least.setdefault(int(label), node)  # nodes are sorted, so the first seen is the least
    return {node: least[int(label)] for node, label in zip(nodes, labels)}


# ── placement ──


def _reading(nodes: set[str], rank: Mapping[str, int], closures: Mapping[str, Mapping[str, frozenset[str]]]) -> dict[str, Any]:
    """Home, tie and purity of one test's layer nodes."""
    top = max(rank[node] for node in nodes)
    homes = sorted(node for node in nodes if rank[node] == top)
    foreign: dict[str, str] = {}
    for node in sorted(nodes - set(homes)):
        foreign[node] = next((kind for kind in CLOSURES if all(node in closures[kind][home] for home in homes)),
                             "none")
    foreign = {node: kind for node, kind in foreign.items() if kind != "runtime"}
    if len(homes) > 1:
        foreign.update({home: "tie" for home in homes})
    reading: dict[str, Any] = {"home": homes[0] if len(homes) == 1 else None, "tie": homes if len(homes) > 1 else None,
                               "foreign": foreign}
    for kind, key in CLOSURES.items():
        reading[key] = len(homes) == 1 and all(node in closures[kind][homes[0]] for node in nodes - set(homes))
    return reading


def place(store_path: str | Path, layer_prefix: str, layer: Layer, order: Sequence[Sequence[str]],
          kinds: Mapping[tuple[str, str], Iterable[str]], *, threshold: float = 0.95, min_defs: int = 3,
          seeds: Sequence[int] = (1, 2, 3), declared: Iterable[tuple[str, str]] = ()) -> dict[str, Any]:
    """The placement of every stored test whose id starts with ``layer_prefix``, against ``order`` (tiers of node
    names, upstream first) over ``layer``'s code and its import graph ``kinds``, (imported node, importer node) to
    the kinds of import that make the edge, and the containment of each ``declared`` ``rests_on`` edge (module
    docstring)."""
    graph = {edge: frozenset(found) for edge, found in kinds.items()}
    unknown_kinds = set().union(*graph.values()) - set(CLOSURES) if graph else set()
    if unknown_kinds:
        raise ValueError(f"cone_store place: import kinds this placement has no closure for: {sorted(unknown_kinds)}")
    rank = {node: i for i, tier in enumerate(order) for node in tier}
    listed = list(CLOSURES)
    closures = {kind: ancestors(graph, rank, listed[: listed.index(kind) + 1]) for kind in CLOSURES}

    if not Path(store_path).is_file():
        raise ValueError(f"cone_store place: no store at {store_path}")
    found = stored_format(store_path)
    if found != FORMAT:  # read before opening: ConeStore empties a store of another form
        raise ValueError(f"cone_store place: {store_path} is a store of form {found}, not {FORMAT}: it holds no own "
                         "cones (written before B3), so there is nothing to place; capture again")
    with ConeStore(store_path) as store:
        stored = {test_id: test for test_id, test in store.tests().items() if test_id.startswith(layer_prefix)}
        owned = {test_id: own_defs(store.own_pins(test_id), layer) for test_id in sorted(stored)}
    if not stored:
        raise ValueError(f"cone_store place: no stored test starts with {layer_prefix!r}")

    tests: dict[str, Any] = {}
    placed: dict[str, set[str]] = {}
    for test_id in sorted(stored):
        test, defs = stored[test_id], owned[test_id]
        nodes = set(defs.layer)
        entry: dict[str, Any] = {
            "upstream_layers": dict(sorted(defs.upstream.items())), "above_layer": dict(sorted(defs.above.items())),
            "helpers": defs.helpers, "elsewhere": defs.elsewhere, "unplaced": None}
        if test.unknown:
            entry["unplaced"] = "cone unknown: " + "; ".join(test.unknown)
        elif test.outcome != "passed":
            entry["unplaced"] = f"not passed: {test.outcome}"
        elif not nodes:
            entry["unplaced"] = "no layer def"
        if entry["unplaced"] is None:
            entry.update(_reading(nodes, rank, closures))
            placed[test_id] = nodes
        tests[test_id] = entry

    def pure_fraction(per_test: Mapping[str, set[str]]) -> float:
        readings = [_reading(nodes, rank, closures)["pure_rt"] for nodes in per_test.values() if nodes]
        return round(sum(readings) / len(readings), 6) if readings else 0.0

    null = {}
    labelled = sorted({(node, qualname) for test_id in placed for node, names in owned[test_id].layer.items()
                       for qualname in names})
    for seed in seeds:
        labels = [node for node, _ in labelled]
        random.Random(seed).shuffle(labels)
        moved = dict(zip(labelled, labels))
        shuffled = {test_id: {moved[(node, q)] for node, names in owned[test_id].layer.items() for q in names}
                    for test_id in placed}
        null[str(seed)] = pure_fraction(shuffled)

    supports = _supports({test_id: owned[test_id] for test_id in placed}, tests, closures["typing"], threshold,
                         min_defs)
    readings = [tests[test_id] for test_id in placed]
    summary = {
        "tests": len(tests), "placed": len(placed),
        "unplaced": dict(sorted(Counter(entry["unplaced"].split(":")[0] for entry in tests.values()
                                        if entry["unplaced"]).items())),
        "pure": {kind: sum(r[key] for r in readings) for kind, key in CLOSURES.items()},
        "pure_fraction": pure_fraction(placed), "null_pure_fraction": null, "ties": sum(1 for r in readings if r["tie"]),
        "foreign_kinds": dict(sorted(Counter(k for r in readings for k in r["foreign"].values()).items())),
        "homes": dict(sorted(Counter(r["home"] or "+".join(r["tie"]) for r in readings).items())),
        "above_layer_tests": sum(1 for entry in tests.values() if entry["above_layer"]),
        "straddlers": sorted((test_id for test_id in placed if not tests[test_id]["pure_rt"]),
                             key=lambda t: (sorted(set(tests[t]["foreign"].values())), t)),
    }
    judged = [test_id for test_id, entry in tests.items()
              if entry["unplaced"] is None or entry["unplaced"].startswith("no layer def")]
    rests_on = _declared(_functions({test_id: owned[test_id] for test_id in judged}), declared)
    return {"layer": layer_prefix, "code_root": layer.root, "threshold": threshold, "min_defs": min_defs,
            "summary": summary, "tests": tests, **supports, **rests_on}


def _function(test_id: str) -> str:
    """The test function of a test id: a parametrised case is one claim with its function (ruling 4)."""
    return test_id.split("[")[0]


def _functions(owned: Mapping[str, OwnDefs]) -> dict[str, set[tuple[str, str]]]:
    """Each test function's own first-party defs, its cases united."""
    functions: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for test_id, defs in owned.items():
        functions[_function(test_id)] |= defs.first_party
    return functions


def _declared(functions: Mapping[str, set[tuple[str, str]]], declared: Iterable[tuple[str, str]]) -> dict[str, Any]:
    """The containment of each declared ``rests_on`` edge, its cumulative counts, and the edges not measured by
    reason (module docstring)."""
    rows: list[dict[str, Any]] = []
    unmeasured: Counter[str] = Counter()
    for upper, target in sorted({(_function(u), t) for u, t in declared}):
        lowers = sorted(f for f in functions if f == _function(target) or f.startswith(target + "::"))
        if upper not in functions:
            unmeasured["upper not judged"] += 1
        elif not lowers:
            unmeasured["target not judged"] += 1
        for lower in lowers if upper in functions else ():
            if lower == upper:  # a class or file target holding the upper itself
                continue
            if not functions[lower]:
                unmeasured["lower ran no def"] += 1
                continue
            c = len(functions[lower] & functions[upper]) / len(functions[lower])
            rows.append({"upper": upper, "lower": lower, "target": target, "c": round(c, 6)})
    counts = {str(t): sum(1 for row in rows if row["c"] >= t) for t in THRESHOLDS}
    counts["below"] = sum(1 for row in rows if row["c"] < THRESHOLDS[-1])
    return {"declared_rests_on": rows, "declared_counts": counts, "declared_unmeasured": dict(sorted(unmeasured.items()))}


def _supports(owned: Mapping[str, OwnDefs], tests: Mapping[str, Any], widest: Mapping[str, frozenset[str]],
              threshold: float, min_defs: int) -> dict[str, Any]:
    """The proposed ``executes`` supports between test functions, the duplicate and near-duplicate candidates, and
    the counts at each of :data:`THRESHOLDS` (module docstring). ``widest`` is each node's ancestors under the
    widest closure, the order a support is read along."""
    functions = _functions(owned)
    homes: dict[str, set[str]] = defaultdict(set)
    for test_id in owned:
        reading = tests[test_id]
        homes[_function(test_id)] |= {reading["home"]} if reading["home"] else set(reading["tie"] or ())
    by_def: dict[tuple[str, str], list[str]] = defaultdict(list)
    for function, cone in functions.items():
        for item in cone:
            by_def[item].append(function)
    containment: dict[tuple[str, str], float] = {}
    for lower, cone in functions.items():
        if len(cone) < min_defs:
            continue
        shared: Counter[str] = Counter(upper for item in cone for upper in by_def[item] if upper != lower)
        for upper, count in shared.items():
            containment[(upper, lower)] = count / len(cone)
    duplicates = sorted({tuple(sorted((u, l))) for (u, l), c in containment.items()
                         if c == 1.0 and functions[u] == functions[l]})
    duplicate_pairs = {pair for a, b in duplicates for pair in ((a, b), (b, a))}
    counts = {str(t): sum(1 for pair, c in containment.items() if c >= t and pair not in duplicate_pairs)
              for t in THRESHOLDS}
    relation = {pair for pair, c in containment.items() if c >= threshold and pair not in duplicate_pairs}
    group = _groups(relation)
    members: dict[str, set[str]] = defaultdict(set)
    for function, name in group.items():
        members[name].add(function)
    near = sorted(tuple(sorted(found)) for found in members.values() if len(found) > 1)
    condensed: dict[str, set[str]] = defaultdict(set)
    for upper, lower in relation:
        if group[upper] != group[lower]:
            condensed[group[upper]].add(group[lower])
    reach = _reach(condensed, sorted(set(group.values())))
    reduced = sorted((upper, lower) for upper, lower in relation if group[upper] != group[lower]
                     and not any(group[lower] in reach[middle] for middle in condensed[group[upper]]))
    proposals, against = [], []
    for upper, lower in reduced:
        along = any(lh == uh or lh in widest.get(uh, frozenset()) for lh in homes[lower] for uh in homes[upper])
        record = {"upper": upper, "lower": lower, "c": round(containment[(upper, lower)], 6)}
        (proposals if along else against).append(record)
    return {"supports": proposals, "supports_against_order": against,
            "duplicates": [{"a": a, "b": b, "defs": len(functions[a])} for a, b in duplicates],
            "near_duplicates": [list(found) for found in near], "support_counts": counts}
