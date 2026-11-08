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
* A test is UNPLACED, with its reason, when its cone is unknown, it did not pass, or it read nothing of the layer or
  of the apparatus.

What a test READ (B5, rulings 47 to 49; :mod:`tools.cone_store.reads`, the ``own_reads`` rows) places a test that ran no
layer def: a test of a module's declared structure (a dataclass's fields, an enum's members, a constant) reads a
statement that ran at collection, and a source audit opens a layer file as data. The home is TIERED (ruling 49): the
defs decide it when the test ran any layer def; else the layer nodes whose bindings it read (``statement``, by the
channels whose :attr:`~tools.cone_store.reads.Channel.is_node` holds); else the files it opened as data (``data``;
``layer`` when it opened every node, ``package`` when it opened files of the package outside the layer too, a census
finding); else the APPARATUS (ruling 48): a test-tree file whose bindings it read or whose defs it ran, its own file
first, a tie otherwise. :func:`home_of` returns the home as one value of the :data:`Home` sum, matched once. Running a
module's code is a dependence on its behaviour; reading a value with none of its module's code run is a dependence on an
input, so a statement read never moves a def home. Every channel is reported under ``via`` whatever homes the test, with
the union home (every channel's nodes, highest rank wins: the refuted reading, kept) and the statement straddle (a
def-placed test's statement nodes its home does not explain at runtime). Purity, the straddlers, the null and the
proposed supports are read over the def-placed tests, as in B3; the statement-placed ones have their own pure fraction
and null, over (node, item) labels; a fraction over no test is ``None``. A test homed at an apparatus file is the
proposed lower of every test function whose items at that file contain its own (``apparatus_supports``), and each
declared ``rests_on`` edge is measured over defs and read items too (``c_items``).

Proposed ``executes`` supports (S17: containment is fractional): over the test FUNCTIONS (a parametrised test's
cases are one claim, ruling 4; their own defs are united), ``c(U ⊇ L) = |D(L) ∩ D(U)| / |D(L)|`` over own defs
under the first-party package. ``U -> L`` is proposed when ``c >= threshold``, ``|D(L)| >= min_defs`` and the two
cones are not equal (equal cones of at least ``min_defs`` defs are DUPLICATE CANDIDATES, with no verdict: ruling 21
decides on kill sets; two smaller equal cones are below the size a containment is read at). The lower's SUBJECT, its
own defs at its home nodes (each tied node for a tie, its cases united), is contained at the same threshold, ``c_home(U ⊇ L) =
|S(L) ∩ D(U)| / |S(L)| >= threshold`` (ruling 50): the 5 % the whole cone may miss cannot all fall on the code the
lower is about, since an upper that never runs the lower's subject cannot be reached by a defect in it.
Fractional containment can cycle, so functions that contain each other at the threshold are a NEAR-DUPLICATE group,
also listed with no verdict, and the relation is transitively reduced over its condensation (each group one node).
A pair whose lower home is neither the upper home nor upstream of it under the widest closure is listed apart,
against the order, for the orchestrator's reading (S18).

The declared ``rests_on`` edges (``declared``, each (upper test id, target id) as pytest resolved the marker:
``tests/_harness/rests_on_manifest.py``, the one resolver) are measured on the same function cones: each target
collapses to its functions (a case to its function, a class or a file to every function under it), and each edge
reports ``c(U ⊇ L)`` and, for a lower placed by its defs, ``c_home`` (``None`` otherwise), with the cumulative
counts over ``c`` at :data:`THRESHOLDS` and the count below the lowest, so the threshold is judged on the layer's
whole declared population, every test that passed with a known cone included, whatever homes it. An edge whose
upper or target holds no such function, or whose lower ran no def and read nothing, is counted by that reason. A
lower that ran no def but read items has no ``c``: it is measured by ``c_items`` alone (over defs and read items,
B5), and counted apart (``declared_items_only``).

The null shuffles the node label of each distinct (node, def) across the layer with a fixed seed and re-reads purity;
a placement that the null matches carries no information (S13). The output is a plain dict, byte-identical across runs
when dumped with sorted keys.
"""
from __future__ import annotations

import random
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence, Set
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Any, ClassVar

import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import breadth_first_order, connected_components

from orpheus.numerics import traced_memo as memo

from tools.cone_store.reads import MODULE_ITEM, Channel, Read, in_test_tree, is_hook
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
#: What a store of an earlier form lacks that placement reads (``store.FORMAT``'s history).
_LACKS: dict[int | None, str] = {2: "own cones (written before B3)", 3: "own reads (written before B5)"}


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
        if in_test_tree(path):
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
    """A test's own defs, sorted by region: the layer's per node, the test-tree defs per file (its own def, the code
    nested in it and the hooks excepted: the apparatus it ran, ruling 48), and the counts elsewhere."""

    layer: dict[str, frozenset[tuple[str, str]]]
    upstream: Counter[str] = field(default_factory=Counter)
    above: Counter[str] = field(default_factory=Counter)
    helpers: dict[str, frozenset[str]] = field(default_factory=dict)
    elsewhere: int = 0
    first_party: frozenset[tuple[str, str]] = frozenset()

    def at(self, nodes: Iterable[str]) -> frozenset[tuple[str, str]]:
        """The layer defs at ``nodes``, each as (path, qualname)."""
        return frozenset().union(*(self.layer.get(node, frozenset()) for node in nodes))


def own_defs(pins: Iterable[Any], layer: Layer, test_id: str) -> OwnDefs:
    by_node: dict[str, set[tuple[str, str]]] = defaultdict(set)
    upstream: Counter[str] = Counter()
    above: Counter[str] = Counter()
    helpers: dict[str, set[str]] = defaultdict(set)
    elsewhere = 0
    first_party: set[tuple[str, str]] = set()
    own = _test_qualname(test_id)
    for pin in pins:
        if not isinstance(pin, memo.DefPin):
            continue
        match layer.region(pin.path):
            case InLayer(node):
                by_node[node].add((pin.path, pin.qualname))
                first_party.add((pin.path, pin.qualname))
            case Upstream(module):
                upstream[module] += 1
                first_party.add((pin.path, pin.qualname))
            case Above(package):
                above[package] += 1
                first_party.add((pin.path, pin.qualname))
            case Helper() if not _is_own(pin.qualname, own, pin.path, test_id) and not is_hook(pin.qualname):
                helpers[pin.path].add(pin.qualname)
            case Helper():
                pass
            case Elsewhere():
                elsewhere += 1
    return OwnDefs({node: frozenset(q) for node, q in by_node.items()}, upstream, above,
                   {path: frozenset(q) for path, q in helpers.items()}, elsewhere, frozenset(first_party))


@dataclass(frozen=True)
class OwnReads:
    """What a test's own run read (B5, :mod:`tools.cone_store.reads`), sorted by channel: the layer's nodes it read a
    binding of (``statement``, node to the names read, a module read adding the node with no name), the layer's nodes
    it opened as data (``data``) and whether it opened a file of the package outside the layer too, the apparatus it
    read (test-tree file to items: the bindings it read there and the test-tree defs it ran, its own def and the
    hooks excepted), its reads of the package above the layer, and the names no statement binds."""

    statement: dict[str, frozenset[str]]
    data: frozenset[str]
    data_outside: bool
    apparatus: dict[str, frozenset[str]]
    above: int
    unresolved: tuple[str, ...]


def own_reads(reads: Iterable[Read], pins: Iterable[Any], layer: Layer, defs: OwnDefs) -> OwnReads:
    """A test's reads sorted by channel, with the files it opened as data from its pins and the apparatus defs it ran
    from ``defs``."""
    statement: dict[str, set[str]] = defaultdict(set)
    apparatus: dict[str, set[str]] = defaultdict(set, {path: set(q) for path, q in defs.helpers.items()})
    data: set[str] = set()
    above = 0
    data_outside = False
    unresolved = sorted(read.name for read in reads if read.channel is Channel.UNRESOLVED and read.name)
    for read in reads:
        if not read.channel.is_node:
            continue
        match layer.region(read.path):
            case InLayer(node):
                statement[node] |= set() if read.name is None else {read.name}
            case Helper() if read.name is not None:
                apparatus[read.path].add(read.name)
            case Above():
                above += 1
            case _:
                pass
    for pin in pins:
        match pin:
            case memo.DataPin(path=path) if (node := layer.node_of(path)) is not None:
                data.add(node)
            case memo.DataPin(path=path) if path.startswith(layer.package + "/") \
                    and not path.startswith(layer.root.rstrip("/") + "/"):
                data_outside = True  # the layer's own ``__init__`` is in the layer, and no node
            case _:
                pass
    return OwnReads({node: frozenset(names) for node, names in statement.items()}, frozenset(data), data_outside,
                    {path: frozenset(items) for path, items in apparatus.items()}, above, tuple(unresolved))


def _test_qualname(test_id: str) -> str:
    """The qualified name of a test's function: ``f.py::C::m[case]`` is ``C.m``."""
    return ".".join(test_id.split("[")[0].split("::")[1:])


def _is_own(qualname: str, own: str, path: str, test_id: str) -> bool:
    """Whether a def is the test's own (its function, or code nested in it): not apparatus (A8)."""
    return path == test_id.split("::")[0] and (qualname == own or qualname.startswith(own + ".<locals>."))


# ── how a test is homed ──


@dataclass(frozen=True)
class DefHome:
    """Homed by the layer nodes whose defs it ran."""

    nodes: frozenset[str]
    by: ClassVar[str] = "def"


@dataclass(frozen=True)
class StatementHome:
    """Homed by the layer nodes whose bindings it read, having run no layer def."""

    nodes: frozenset[str]
    by: ClassVar[str] = "statement"


@dataclass(frozen=True)
class DataHome:
    """Homed by the layer nodes whose files it opened as data, having run and read none."""

    nodes: frozenset[str]
    by: ClassVar[str] = "data"


@dataclass(frozen=True)
class ScopeHome:
    """It opened every node's file (``layer``), and files of the package outside the layer too (``package``, a census
    finding, C13): no node is the home."""

    by: str


@dataclass(frozen=True)
class ApparatusHome:
    """Homed at the test-tree files whose bindings it read or whose defs it ran (ruling 48), its own file first; two
    other files are a tie (C14)."""

    files: frozenset[str]
    by: ClassVar[str] = "apparatus"


Home = DefHome | StatementHome | DataHome | ScopeHome | ApparatusHome


def home_of(test_id: str, defs: OwnDefs, reads: OwnReads, layer_nodes: Iterable[str]) -> Home | None:
    """How a test is homed (ruling 49): by its defs when it ran any layer def, else by the bindings it read, else by
    the files it opened as data (the whole layer when it opened every node, the package when it opened more), else
    by the apparatus; ``None`` when it read nothing of either."""
    if defs.layer:
        return DefHome(frozenset(defs.layer))
    if reads.statement:
        return StatementHome(frozenset(reads.statement))
    if reads.data:
        if reads.data == set(layer_nodes):
            return ScopeHome("package" if reads.data_outside else "layer")
        return DataHome(reads.data)
    if reads.apparatus:
        own = test_id.split("::")[0]
        return ApparatusHome(frozenset({own}) if own in reads.apparatus else frozenset(reads.apparatus))
    return None


def _home_nodes(reading: Mapping[str, Any]) -> frozenset[str]:
    """A reading's home as a node set: the home, or every tied node."""
    return frozenset([reading["home"]] if reading["home"] else reading["tie"] or ())


def _label(reading: Mapping[str, Any]) -> str:
    """A reading's home, or its tie joined by ``+``."""
    return reading["home"] or "+".join(reading["tie"])


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
        lacks = _LACKS.get(found, "rows of the form this placement reads")
        raise ValueError(f"cone_store place: {store_path} is a store of form {found}, not {FORMAT}: it holds no "
                         f"{lacks}, so there is nothing to place; capture again")
    with ConeStore(store_path) as store:
        stored = {test_id: test for test_id, test in store.tests().items() if test_id.startswith(layer_prefix)}
        pins = {test_id: store.own_pins(test_id) for test_id in sorted(stored)}
        owned = {test_id: own_defs(pins[test_id], layer, test_id) for test_id in sorted(stored)}
        read = {test_id: own_reads(store.own_reads(test_id), pins[test_id], layer, owned[test_id])
                for test_id in sorted(stored)}
    if not stored:
        raise ValueError(f"cone_store place: no stored test starts with {layer_prefix!r}")

    tests: dict[str, Any] = {}
    placed: dict[str, set[str]] = {}  # placed by defs: the population every def-only reading is taken over
    by_statement: dict[str, set[str]] = {}
    homes: dict[str, Home] = {}
    for test_id in sorted(stored):
        test, defs, reads = stored[test_id], owned[test_id], read[test_id]
        nodes, said = set(defs.layer), set(reads.statement)
        entry: dict[str, Any] = {
            "upstream_layers": dict(sorted(defs.upstream.items())), "above_layer": dict(sorted(defs.above.items())),
            "helpers": sum(map(len, defs.helpers.values())), "elsewhere": defs.elsewhere, "unplaced": None}
        if test.unknown:
            entry["unplaced"] = "cone unknown: " + "; ".join(test.unknown)
        elif test.outcome != "passed":
            entry["unplaced"] = f"not passed: {test.outcome}"
        else:
            entry.update(_channels(defs, reads, nodes | said | reads.data, rank, closures))
            home = home_of(test_id, defs, reads, rank)
            if home is None:
                entry["unplaced"] = "no layer read"
            else:
                homes[test_id] = home
                entry["placed_by"] = home.by
                match home:
                    case DefHome(found):
                        entry.update(_reading(set(found), rank, closures))
                        placed[test_id] = nodes
                        entry["statement_straddle"] = _straddle(said - nodes, entry, closures)
                    case StatementHome(found):
                        entry.update(_reading(set(found), rank, closures))
                        by_statement[test_id] = said
                        entry["statement_straddle"] = dict(entry["foreign"])
                    case DataHome(found):
                        entry.update(_reading(set(found), rank, closures))
                    case ScopeHome():
                        entry.update(_outside())
                    case ApparatusHome(files):
                        entry.update(_outside(files))
        tests[test_id] = entry

    def pure_fraction(per_test: Mapping[str, set[str]]) -> float | None:
        """The pure fraction; ``None`` over no test (a population with no member has no fraction)."""
        readings = [_reading(nodes, rank, closures)["pure_rt"] for nodes in per_test.values() if nodes]
        return round(sum(readings) / len(readings), 6) if readings else None

    def null(labelled_of: Mapping[str, Mapping[str, Iterable[str]]]) -> dict[str, float | None]:
        """The pure fraction with the node label of each distinct (node, item) shuffled across the layer (S13)."""
        labelled = sorted({(node, item) for per in labelled_of.values() for node, items in per.items()
                           for item in items})
        found = {}
        for seed in seeds:
            labels = [node for node, _ in labelled]
            random.Random(seed).shuffle(labels)
            moved = dict(zip(labelled, labels))
            found[str(seed)] = pure_fraction({test_id: {moved[(node, item)] for node, items in per.items()
                                                        for item in items} for test_id, per in labelled_of.items()})
        return found

    statement_items = {test_id: {node: sorted(read[test_id].statement[node]) or [MODULE_ITEM] for node in said}
                       for test_id, said in by_statement.items()}
    function_homes = _function_homes({test_id: tests[test_id] for test_id in placed})
    subjects = _subjects({test_id: owned[test_id] for test_id in placed}, function_homes)
    supports = _supports({test_id: owned[test_id] for test_id in placed}, function_homes, subjects,
                         closures["typing"], threshold, min_defs)
    readings = [tests[test_id] for test_id in placed]
    union_homes = {test_id: e["home_union"] for test_id, e in tests.items() if e.get("home_union")}
    by_node = {test_id for test_id, home in homes.items() if isinstance(home, DefHome | StatementHome)}
    summary = {
        "tests": len(tests), "placed": len(placed),
        "placed_by": dict(sorted(Counter(e["placed_by"] for e in tests.values() if e.get("placed_by")).items())),
        "unplaced": dict(sorted(Counter(entry["unplaced"].split(":")[0] for entry in tests.values()
                                        if entry["unplaced"]).items())),
        "pure": {kind: sum(r[key] for r in readings) for kind, key in CLOSURES.items()},
        "pure_fraction": pure_fraction(placed),
        "null_pure_fraction": null({t: {node: {name for _, name in defs} for node, defs in owned[t].layer.items()}
                                    for t in placed}),
        "pure_fraction_statement": pure_fraction(by_statement),
        "null_pure_fraction_statement": null(statement_items),
        "ties": sum(1 for r in readings if r["tie"]),
        "foreign_kinds": dict(sorted(Counter(k for r in readings for k in r["foreign"].values()).items())),
        "homes": dict(sorted(Counter(_label(r) for r in readings).items())),
        "homes_union": dict(sorted(Counter(union_homes.values()).items())),
        "union_moves": sum(1 for test_id, home in union_homes.items()
                           if test_id in by_node and home != _label(tests[test_id])),
        "above_layer_tests": sum(1 for entry in tests.values() if entry["above_layer"]),
        "above_layer_statement": sum(1 for test_id in tests if read[test_id].above),
        "straddlers": sorted((test_id for test_id in placed if not tests[test_id]["pure_rt"]),
                             key=lambda t: (sorted(set(tests[t]["foreign"].values())), t)),
        "statement_straddlers": sorted(t for t, e in tests.items() if e.get("statement_straddle")),
    }
    # every test that passed with a known cone, whatever homes it: a lower that read nothing is counted by reason
    judged = [test_id for test_id, test in stored.items() if not test.unknown and test.outcome == "passed"]
    claims = _functions({test_id: owned[test_id] for test_id in judged},
                        {test_id: _items(read[test_id]) for test_id in judged})
    rests_on = _declared(_functions({test_id: owned[test_id] for test_id in judged}), declared, claims, subjects)
    return {"layer": layer_prefix, "code_root": layer.root, "threshold": threshold, "min_defs": min_defs,
            "summary": summary, "tests": tests, **supports, **rests_on,
            "apparatus_supports": _apparatus_supports(homes, read, judged, threshold)}


def _outside(files: Iterable[str] = ()) -> dict[str, Any]:
    """The reading of a test no layer node homes: the test-tree files it is homed at (none for the layer or the
    package), never pure."""
    homes = sorted(files)
    return {"home": homes[0] if len(homes) == 1 else None, "tie": homes if len(homes) > 1 else None, "foreign": {},
            **{key: False for key in CLOSURES.values()}}


def _channels(defs: OwnDefs, reads: OwnReads, every: set[str], rank: Mapping[str, int],
              closures: Mapping[str, Mapping[str, frozenset[str]]]) -> dict[str, Any]:
    """Every channel of a test, reported whatever homes it: each layer node's channels, the items read per node
    and per apparatus file, the union home (every channel's nodes, highest rank wins: ruling 49's refuted reading,
    kept), the unresolved names and the apparatus items it rests on."""
    via: dict[str, list[str]] = {}
    for node in sorted(every):
        via[node] = [channel for channel, held in (("def", node in defs.layer), ("statement", node in reads.statement),
                                                   ("data", node in reads.data)) if held]
    items = {node: sorted(names) for node, names in sorted(reads.statement.items()) if names}
    items |= {path: sorted(names) for path, names in sorted(reads.apparatus.items())}
    union = _reading(every, rank, closures) if every else None
    return {"via": via, "items": items,
            "home_union": None if union is None else _label(union),
            "unresolved": list(reads.unresolved),
            "apparatus_items": sorted([path, name] for path, names in reads.apparatus.items() for name in names)}


def _straddle(nodes: set[str], entry: Mapping[str, Any], closures: Mapping[str, Mapping[str, frozenset[str]]]) -> dict[str, str]:
    """The statement nodes of a def-placed test that its def home does not explain at runtime, each with the kind of
    the first import edge that does (``none``: no path), as ``foreign`` reads them."""
    homes = _home_nodes(entry)
    found = {node: next((kind for kind in CLOSURES if all(node in closures[kind][home] for home in homes)), "none")
             for node in sorted(nodes - homes)}
    return {node: kind for node, kind in found.items() if kind != "runtime"}


def _items(reads: OwnReads) -> set[tuple[str, str]]:
    """A test's read items as claim members beside its defs: (layer node or apparatus file, name)."""
    return {(f"@{node}", name) for node, names in reads.statement.items() for name in names} | \
        {(path, name) for path, names in reads.apparatus.items() for name in names}


def _apparatus_supports(homes: Mapping[str, Home], read: Mapping[str, OwnReads], judged: Iterable[str],
                        threshold: float) -> list[dict[str, Any]]:
    """The proposed premise edges (C14, ruling 48): a test function homed at an apparatus file is the lower of every
    other function whose items at that file contain its own, ``c = |I(P) ∩ I(U)| / |I(P)|``, over functions."""
    items: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    premises: dict[str, set[str]] = defaultdict(set)
    for test_id in judged:
        function = _function(test_id)
        for path, names in read[test_id].apparatus.items():
            items[function][path] |= names
        if isinstance(home := homes.get(test_id), ApparatusHome):
            premises[function] |= home.files
    edges = []
    for lower, nodes in sorted(premises.items()):
        for node in sorted(nodes):
            premise = items[lower][node]
            if not premise:
                continue
            for upper in sorted(items):
                if upper == lower or not items[upper].get(node):
                    continue
                c = _contained(premise, items[upper][node])
                if c is not None and c >= threshold:
                    edges.append({"upper": upper, "lower": lower, "node": node, "c": _rounded(c)})
    return edges


def _function(test_id: str) -> str:
    """The test function of a test id: a parametrised case is one claim with its function (ruling 4)."""
    return test_id.split("[")[0]


def _functions(owned: Mapping[str, OwnDefs],
               items: Mapping[str, set[tuple[str, str]]] | None = None) -> dict[str, set[tuple[str, str]]]:
    """Each test function's own first-party defs, its cases united; with ``items``, its read items beside them."""
    functions: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for test_id, defs in owned.items():
        functions[_function(test_id)] |= defs.first_party | (items or {}).get(test_id, set())
    return functions


def _function_homes(readings: Mapping[str, Mapping[str, Any]]) -> dict[str, frozenset[str]]:
    """Each test function's home nodes, its cases united (ruling 4: a parametrised test is one claim)."""
    homes: dict[str, set[str]] = defaultdict(set)
    for test_id, reading in readings.items():
        homes[_function(test_id)] |= _home_nodes(reading)
    return {function: frozenset(nodes) for function, nodes in homes.items()}


def _subjects(owned: Mapping[str, OwnDefs],
              homes: Mapping[str, frozenset[str]]) -> dict[str, frozenset[tuple[str, str]]]:
    """Each def-placed test function's subject (ruling 50): the function's own defs at its home nodes, its cases
    united."""
    subjects: dict[str, frozenset[tuple[str, str]]] = {}
    for test_id, defs in owned.items():
        function = _function(test_id)
        subjects[function] = subjects.get(function, frozenset()) | defs.at(homes[function])
    return subjects


def _contained[T](part: Set[T], cone: Set[T]) -> float | None:
    """``|part ∩ cone| / |part|``, the fraction of ``part`` inside ``cone``; ``None`` for an empty part (a
    containment over no member has no value)."""
    return len(part & cone) / len(part) if part else None


def _rounded(value: float | None) -> float | None:
    """A containment as the report writes it."""
    return None if value is None else round(value, 6)


def _declared(functions: Mapping[str, set[tuple[str, str]]], declared: Iterable[tuple[str, str]],
              claims: Mapping[str, set[tuple[str, str]]],
              subjects: Mapping[str, frozenset[tuple[str, str]]]) -> dict[str, Any]:
    """The containment of each declared ``rests_on`` edge over defs (``c``), of the lower's subject (``c_home``, for
    a lower placed by its defs, ruling 50) and over defs and read items (``c_items``, B5), its cumulative counts over
    ``c``, and the edges not measured by reason (module docstring). A lower test that ran no def but read items is
    measured by ``c_items`` alone (``c`` is ``None``)."""
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
            low, high = claims.get(lower, set()), claims.get(upper, set())
            if not functions[lower] and not low:
                unmeasured["lower ran no def and read nothing"] += 1
                continue
            rows.append({"upper": upper, "lower": lower, "target": target,
                         "c": _rounded(_contained(functions[lower], functions[upper])),
                         "c_home": _rounded(_contained(subjects.get(lower, frozenset()), functions[upper])),
                         "c_items": _rounded(_contained(low, high))})
    counts = {str(t): sum(1 for row in rows if row["c"] is not None and row["c"] >= t) for t in THRESHOLDS}
    counts["below"] = sum(1 for row in rows if row["c"] is not None and row["c"] < THRESHOLDS[-1])
    return {"declared_rests_on": rows, "declared_counts": counts,
            "declared_items_only": sum(1 for row in rows if row["c"] is None),  # measured by ``c_items`` alone
            "declared_unmeasured": dict(sorted(unmeasured.items()))}


def _supports(owned: Mapping[str, OwnDefs], homes: Mapping[str, frozenset[str]],
              subjects: Mapping[str, frozenset[tuple[str, str]]], widest: Mapping[str, frozenset[str]],
              threshold: float, min_defs: int) -> dict[str, Any]:
    """The proposed ``executes`` supports between test functions, the duplicate and near-duplicate candidates, and
    the counts at each of :data:`THRESHOLDS` (module docstring). ``widest`` is each node's ancestors under the
    widest closure, the order a support is read along."""
    functions = _functions(owned)
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

    def c_home(pair: tuple[str, str]) -> float | None:
        upper, lower = pair
        return _contained(subjects[lower], functions[upper])

    def holds(pair: tuple[str, str], at: float) -> bool:
        """``U -> L`` at ``at``: the lower's cone and its subject both contained (ruling 50)."""
        home = c_home(pair)
        return containment[pair] >= at and home is not None and home >= at

    duplicates = sorted({tuple(sorted((u, l))) for (u, l), c in containment.items()
                         if c == 1.0 and functions[u] == functions[l]})
    duplicate_pairs = {pair for a, b in duplicates for pair in ((a, b), (b, a))}
    counts = {str(t): sum(1 for pair in containment if holds(pair, t) and pair not in duplicate_pairs)
              for t in THRESHOLDS}
    relation = {pair for pair in containment if holds(pair, threshold) and pair not in duplicate_pairs}
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
        record = {"upper": upper, "lower": lower, "c": round(containment[(upper, lower)], 6),
                  "c_home": _rounded(c_home((upper, lower)))}
        (proposals if along else against).append(record)
    return {"supports": proposals, "supports_against_order": against,
            "duplicates": [{"a": a, "b": b, "defs": len(functions[a])} for a, b in duplicates],
            "near_duplicates": [list(found) for found in near], "support_counts": counts}
