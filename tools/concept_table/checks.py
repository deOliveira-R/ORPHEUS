"""The checks of the concept table: pure functions over values, each returning the frozenset of its violations.

A violation is a value of the type of the debt entry that would pin it (``EdgeDebt``, ``ClassDebt``, ...), except for the
checks of the table, the declarations and the heads, whose violations are messages. No check takes the debt as an
argument: the gate compares what a check measures with the pinned list by ``debt.pin``. Every message names the class
or the concept at fault and the file that fixes it.
"""

from __future__ import annotations

import ast
import collections
import functools
import hashlib
import io
import re
import tokenize
from collections.abc import Callable, Collection, Iterable, Mapping
from pathlib import Path

from .declaration import parse_declaration
from .debt import (ClassDebt, CovectorOffLineageDebt, CycleDebt, Debt, EdgeDebt, KindLawDebt, OffAlgebraDebt, RoleTestDebt,
                   SeedEntry, wave_of)
from .population import ClassKey, ClassRow, Population
from .refs import Edge, alias_leaves, read_function_edges
from .table import ClassHead, Concept, FunctionsHead, MissingHead, SumHead, heads, rank

TABLE_FILE = "docs/architecture/concept_table.toml"
CORE_ROLES = ("head", "arm", "weld", "twin")
OPERATOR_PROTOCOL = ClassKey("orpheus/numerics/operator.py", "LinearOperator")
MEMBER_CLASS = ClassKey("orpheus/numerics/space.py", "Member")
COVECTOR_CLASS = ClassKey("orpheus/numerics/space.py", "Covector")
_LABEL = re.compile(r"^\.\. _([^:\n][^\n]*?):\s*$")
_MATH_LABEL = re.compile(r"^\s+:label:\s+(\S+)\s*$")


# --- shared helpers --------------------------------------------------------------------------------------------------

def package(key: ClassKey) -> str:
    """The top-level package of a class's file: ``orpheus/sn/x.py`` is ``sn`` and ``orpheus/plotting.py`` is ``plotting``."""
    return key.file.split("/")[1].removesuffix(".py")


def may_reference(a: str, b: str, forbidden: Mapping[str, Collection[str]]) -> bool:
    """Whether package ``a`` may reference package ``b``: ``a == b``, or ``b`` is not forbidden to ``a``.

    A package with no key in ``forbidden`` may reference every package. The one predicate of the layer order, used by
    the layer check and by the placement check.
    """
    return a == b or b not in forbidden.get(a, ())


@functools.lru_cache(maxsize=2048)
def _tree(source: str) -> ast.Module:
    return ast.parse(source)


def _class_node(source: str, qualname: str) -> ast.ClassDef:
    """The ``ClassDef`` named by a qualified name (``Outer.Inner`` or ``outer.<locals>.Name``) in ``source``."""
    node: ast.AST = _tree(source)
    for part in qualname.split("."):
        if part == "<locals>":
            continue
        found = None
        for child in ast.walk(node):
            if isinstance(child, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == part and child is not node:
                found = child
                break
        if found is None:
            raise ValueError(f"no class {qualname!r} in the source")
        node = found
    if not isinstance(node, ast.ClassDef):
        raise ValueError(f"{qualname!r} is not a class")
    return node


def _members(source: str, qualname: str) -> frozenset[str]:
    """The non-dunder functions and the annotated fields of a class body."""
    out = set()
    for n in _class_node(source, qualname).body:
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and not n.name.startswith("__"):
            out.add(n.name)
        if isinstance(n, ast.AnnAssign) and isinstance(n.target, ast.Name):
            out.add(n.target.id)
    return frozenset(out)


def _imported_protocol_names(source: str) -> tuple[frozenset[str], frozenset[str]]:
    """The local names bound to ``typing.Protocol`` (or ``typing_extensions.Protocol``) and to the two modules."""
    protocols, modules = set(), set()
    for n in ast.walk(_tree(source)):
        if isinstance(n, ast.ImportFrom) and n.module in ("typing", "typing_extensions") and n.level == 0:
            protocols.update(a.asname or a.name for a in n.names if a.name == "Protocol")
        elif isinstance(n, ast.Import):
            modules.update(a.asname or a.name for a in n.names if a.name in ("typing", "typing_extensions"))
    return frozenset(protocols), frozenset(modules)


def is_protocol(source: str, qualname: str) -> bool:
    """Whether a class has ``typing.Protocol`` as a base, resolved through the module's imports and never by substring.

    The base is a bare name imported from ``typing``, a subscript ``Protocol[...]`` of one, or ``typing.Protocol``.
    """
    protocols, modules = _imported_protocol_names(source)
    for base in _class_node(source, qualname).bases:
        base = base.value if isinstance(base, ast.Subscript) else base
        if isinstance(base, ast.Name) and base.id in protocols:
            return True
        if isinstance(base, ast.Attribute) and base.attr == "Protocol" and isinstance(base.value, ast.Name) and base.value.id in modules:
            return True
    return False


def _declared(population: Population) -> dict[ClassKey, ClassRow]:
    return {r.key: r for r in population.rows if r.declaration is not None}


def _inheritors(edges: Iterable[Edge]) -> dict[ClassKey, frozenset[ClassKey]]:
    """The transitive ancestors of every class over the ``inherits`` edges."""
    parents: dict[ClassKey, set[ClassKey]] = collections.defaultdict(set)
    for e in edges:
        if e.kind == "inherits":
            parents[e.src].add(e.dst)
    out: dict[ClassKey, frozenset[ClassKey]] = {}
    for start in list(parents):
        seen: set[ClassKey] = set()
        stack = list(parents[start])
        while stack:
            k = stack.pop()
            if k not in seen:
                seen.add(k)
                stack.extend(parents.get(k, ()))
        out[start] = frozenset(seen)
    return out


def _tarjan(graph: Mapping[str, Collection[str]]) -> list[frozenset[str]]:
    """The strongly connected components of at least two nodes (iterative, so no recursion limit applies)."""
    index: dict[str, int] = {}
    low: dict[str, int] = {}
    stack: list[str] = []
    on: set[str] = set()
    out: list[frozenset[str]] = []
    counter = 0
    nodes = sorted(set(graph) | {w for ws in graph.values() for w in ws})
    for root in nodes:
        if root in index:
            continue
        work = [(root, iter(sorted(graph.get(root, ()))))]
        index[root] = low[root] = counter
        counter += 1
        stack.append(root)
        on.add(root)
        while work:
            v, successors = work[-1]
            advanced = False
            for w in successors:
                if w not in index:
                    index[w] = low[w] = counter
                    counter += 1
                    stack.append(w)
                    on.add(w)
                    work.append((w, iter(sorted(graph.get(w, ())))))
                    advanced = True
                    break
                if w in on:
                    low[v] = min(low[v], index[w])
            if advanced:
                continue
            work.pop()
            if work:
                parent = work[-1][0]
                low[parent] = min(low[parent], low[v])
            if low[v] == index[v]:
                component = []
                while True:
                    w = stack.pop()
                    on.discard(w)
                    component.append(w)
                    if w == v:
                        break
                if len(component) > 1:
                    out.append(frozenset(component))
    return out


# --- G0: the table ---------------------------------------------------------------------------------------------------

def doc_labels(docs_root: Path) -> dict[str, str]:
    """The labels defined in ``docs/**/*.rst``, lowercased, each mapped to its kind.

    ``"ref"`` is an explicit target (``.. _label:``), ``"eq"`` the label of a math directive (``:label:``), and ``"both"``
    a name that the docs define as both, which ``:ref:`` and ``:eq:`` each resolve.
    """
    out: dict[str, str] = {}
    for path in sorted(Path(docs_root).rglob("*.rst")):
        if "_build" in path.parts:
            continue
        for line in path.read_text().split("\n"):
            if (m := _LABEL.match(line)) is not None:
                name, kind = m.group(1).strip().lower(), "ref"
            elif (m := _MATH_LABEL.match(line)) is not None:
                name, kind = m.group(1).lower(), "eq"
            else:
                continue
            out[name] = kind if out.get(name, kind) == kind else "both"
    return out


def function_resolver(root: Path) -> Callable[[str], bool]:
    """A resolver ``name -> bool`` for the dotted names of a ``functions`` field, reading ``root/orpheus`` by ``ast``.

    The longest dotted prefix that is a file ``orpheus/<path>.py`` or ``orpheus/<path>/__init__.py`` is the module; the
    rest is ``<name>`` (a module-level ``def``, ``class`` or assignment) or ``<Class>.<member>`` (a ``def``, an annotated
    field or an assignment of the class body). Nothing is imported and no re-export or base class is followed.
    """
    root = Path(root)

    @functools.lru_cache(maxsize=None)
    def tree(path: Path) -> ast.Module:
        return ast.parse(path.read_text())

    def defines(body: list[ast.stmt], name: str, *, members: bool) -> ast.stmt | None:
        kinds = (ast.FunctionDef, ast.AsyncFunctionDef) if members else (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
        for n in body:
            if isinstance(n, kinds) and n.name == name:
                return n
            if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == name for t in n.targets):
                return n
            if isinstance(n, ast.AnnAssign) and isinstance(n.target, ast.Name) and n.target.id == name:
                return n
        return None

    def resolve(name: str) -> bool:
        parts = name.split(".")
        for i in range(len(parts), 0, -1):
            base = root.joinpath(*parts[:i])
            path = base.with_suffix(".py") if base.with_suffix(".py").is_file() else base / "__init__.py"
            if not path.is_file():
                continue
            rest = parts[i:]
            if not rest or len(rest) > 2:
                return False
            node = defines(tree(path).body, rest[0], members=False)
            if node is None:
                return False
            if len(rest) == 1:
                return True
            return isinstance(node, ast.ClassDef) and defines(node.body, rest[1], members=True) is not None
        return False

    return resolve


def check_table(concepts: tuple[Concept, ...], known_layers: Collection[str], expected: int, labels: Mapping[str, str],
                resolve: Callable[[str], bool]) -> frozenset[str]:
    """The facts of the table that need two rows or the tree; every single-row fact was refused when the table was loaded.

    Ids are unique; the row count is ``expected``; each layer is a package of the tree; the slugs are unique; each
    ``functions`` entry resolves; each ``theory`` label is defined in ``docs/`` with the right kind and sits under one
    concept only.
    """
    out: set[str] = set()
    ids = collections.Counter(c.id for c in concepts)
    out.update(f"concept {i!r}: id {i!r} is not unique; fix it in {TABLE_FILE}" for i, n in ids.items() if n > 1)
    if len(concepts) != expected:
        out.add(f"the table has {len(concepts)} rows, not {expected}; fix it in {TABLE_FILE}")
    slugs: dict[str, str] = {}
    owner: dict[str, str] = {}
    for c in concepts:
        if c.layer not in known_layers:
            out.add(f"concept {c.id!r}: layer {c.layer!r} is not a package of orpheus/; fix it in {TABLE_FILE}")
        if slugs.setdefault(c.slug, c.id) != c.id:
            out.add(f"concept {c.id!r}: slug {c.slug!r} is also the slug of {slugs[c.slug]!r}; fix it in {TABLE_FILE}")
        for entry in c.functions:
            if not resolve(entry):
                out.add(f"concept {c.id!r}: functions entry {entry!r} is not defined in the tree; fix it in {TABLE_FILE}")
        for label in c.theory:
            kind, _, name = label.partition(":")
            if labels.get(name.lower()) not in (kind, "both"):
                out.add(f"concept {c.id!r}: theory label {label!r} is not a defined {'equation ' if kind == 'eq' else ''}label of docs/; fix it in {TABLE_FILE}")
            if owner.setdefault(label, c.id) != c.id:
                out.add(f"concept {c.id!r}: theory label {label!r} is also under {owner[label]!r}; fix it in {TABLE_FILE}")
    return frozenset(out)


# --- G1: every class declares one concept ----------------------------------------------------------------------------

def check_declarations(population: Population, concept_ids: Collection[str]) -> frozenset[str]:
    """Every class has exactly one marker on its header-end line; it parses; it is canonical; its concept is a row.

    Every id of its ``welded`` list is a row; no marker sits on a line that ends no class header; no class lies under
    ``derivations/``.
    """
    ids = set(concept_ids)
    out: set[str] = set()
    for row in population.rows:
        where = f"{row.key}"
        if "derivations" in row.key.file.split("/"):
            out.add(f"{where}: a class under derivations/ is outside the population")
        if not row.markers:
            out.add(f"{where}: no declaration; add the trailing comment on line {row.header_end} of {row.key.file}")
            continue
        if len(row.markers) > 1:
            out.add(f"{where}: {len(row.markers)} declarations")
            continue
        try:
            declaration = parse_declaration(row.markers[0])
        except ValueError as error:
            out.add(f"{where}: unparsable ({error})")
            continue
        if declaration.concept not in ids:
            out.add(f"{where}: concept {declaration.concept!r} is not in the table")
        out.update(f"{where}: welded concept {w!r} is not in the table" for w in declaration.welded[1:] if w not in ids)
        if row.markers[0].split("#")[1].strip() != declaration.render()[2:].strip():
            out.add(f"{where}: the declaration is not canonical; write '{declaration.render()}'")
    out.update(f"{file}:{line}: a concept comment on no class header-end line" for file, line in population.orphan_markers)
    return frozenset(out)


# --- G2: every concept has one head ----------------------------------------------------------------------------------

def check_heads(population: Population, concepts: tuple[Concept, ...], sources: Mapping[str, str]) -> frozenset[str]:
    """Every concept has the head its form says: exactly one class (``ClassHead``), none beside a sum, none for a functions head,
    none for a missing head (a missing concept whose head now exists is stale), and the arms and leaves it names are population
    classes that declare the concept."""
    out: set[str] = set()
    declared = _declared(population)
    found = heads(population, concepts)
    for c in concepts:
        have = found[c.id]
        names = ", ".join(k.qualname for k in have)
        head = c.head
        if isinstance(head, ClassHead):
            if not have:
                out.add(f"concept {c.id!r}: 0 heads; declare one class 'role: head' for it")
            elif len(have) > 1:
                out.add(f"concept {c.id!r}: {len(have)} heads ({names})")
        elif isinstance(head, SumHead):
            if have:
                out.add(f"concept {c.id!r}: head beside a sum ({names})")
            leaves = alias_leaves(sources, head.file, head.alias)
            if len(leaves) < 2:
                out.add(f"concept {c.id!r}: sum alias {head.alias!r} not found in {head.file} as a union of at least two classes")
            for leaf in leaves:
                row = declared.get(leaf)
                if row is None or row.declaration is None or row.declaration.concept != c.id or row.declaration.role != "arm":
                    out.add(f"concept {c.id!r}: sum leaf {leaf} does not declare an arm of the concept")
        elif isinstance(head, FunctionsHead):
            if have:
                out.add(f"concept {c.id!r}: a functions head with a head class ({names})")
        elif isinstance(head, MissingHead):
            if have:
                out.add(f"concept {c.id!r}: missing but head exists ({names}); drop 'missing' in {TABLE_FILE}")
            for arm in head.arms:
                row = declared.get(arm)
                if row is None or row.declaration is None or row.declaration.concept != c.id:
                    out.add(f"concept {c.id!r}: arm {arm} is not a class that declares the concept")
    return frozenset(out)


# --- G3: the structural role tests -----------------------------------------------------------------------------------

def check_roles(population: Population, edges: Iterable[Edge], concepts: tuple[Concept, ...],
                sources: Mapping[str, str]) -> frozenset[RoleTestDebt]:
    """The classes that fail the structural test of their role.

    ``arm``: inherits a head, arm, weld or twin of its concept, is a leaf of its concept's sum, is listed in its concept's
    arms, or conforms to its concept's Protocol head. ``twin``: none of those (``twin_is_arm`` otherwise). ``interface``:
    has ``typing.Protocol`` as a base. ``part``: every reference to it comes from its concept (class, a class whose welded
    list holds the concept, or a function of a module defining one); an outside reader that includes a core class of the
    concept makes it ``part_component_candidate``, otherwise ``part_no_core_reader``. ``component``: a head, arm, weld or
    twin of its concept refers to it.
    """
    edges = tuple(edges)
    declared = _declared(population)
    ancestors = _inheritors(edges)
    listed: dict[str, set[ClassKey]] = collections.defaultdict(set)
    for c in concepts:
        if isinstance(c.head, MissingHead):
            listed[c.id].update(c.head.arms)
        elif isinstance(c.head, SumHead):
            listed[c.id].update(alias_leaves(sources, c.head.file, c.head.alias))
    head_of = {c: ks[0] for c, ks in heads(population, concepts).items() if ks}
    readers: dict[ClassKey, set[ClassKey]] = collections.defaultdict(set)
    for e in edges:
        readers[e.dst].add(e.src)
    function_readers: dict[ClassKey, set[str]] = collections.defaultdict(set)
    for fe in read_function_edges(sources, population):
        function_readers[fe.dst].add(fe.src_file)
    concept_files: dict[str, set[str]] = collections.defaultdict(set)
    for r in declared.values():
        d = r.declaration
        if d is not None:
            for cid in (d.concept, *d.welded[1:]):
                concept_files[cid].add(r.key.file)

    def inside(src: ClassKey, cid: str) -> bool:
        row = declared.get(src)
        return row is not None and row.declaration is not None and cid in (row.declaration.concept, *row.declaration.welded)

    out: set[RoleTestDebt] = set()
    for k, row in declared.items():
        d = row.declaration
        assert d is not None
        cid = d.concept
        lineage = any(a in declared and declared[a].declaration.concept == cid and declared[a].declaration.role in CORE_ROLES  # type: ignore[union-attr]
                      for a in ancestors.get(k, ()))
        is_listed = k in listed.get(cid, ())
        head = head_of.get(cid)
        conforms = False
        if head is not None and head != k and is_protocol(sources[head.file], head.qualname):
            wanted = _members(sources[head.file], head.qualname)
            conforms = bool(wanted) and wanted <= _members(sources[k.file], k.qualname)
        if d.role == "arm" and not (lineage or is_listed or conforms):
            out.add(RoleTestDebt("arm", k))
        if d.role == "twin" and (lineage or is_listed or conforms):
            out.add(RoleTestDebt("twin_is_arm", k))
        if d.role == "interface" and not is_protocol(sources[k.file], k.qualname):
            out.add(RoleTestDebt("interface", k))
        core = [s for s in readers.get(k, ()) if s in declared and declared[s].declaration.concept == cid  # type: ignore[union-attr]
                and declared[s].declaration.role in CORE_ROLES]  # type: ignore[union-attr]
        if d.role == "component" and not core:
            out.add(RoleTestDebt("component_no_core_reader", k))
        if d.role == "part":
            outside = [s for s in readers.get(k, ()) if not inside(s, cid)]
            outside += [f for f in function_readers.get(k, ()) if f not in concept_files.get(cid, ())]
            if outside:
                out.add(RoleTestDebt("part_component_candidate" if core else "part_no_core_reader", k))
    return frozenset(out)


# --- G4: the kind laws -----------------------------------------------------------------------------------------------

def _operator_hits(population: Population, edges: Iterable[Edge], sources: Mapping[str, str]) -> frozenset[ClassKey]:
    """The classes of the operator algebra: ``LinearOperator`` itself, its inheritors, and the classes whose body defines
    every non-dunder member of the Protocol. Empty when the population does not hold ``LinearOperator``."""
    keys = {r.key for r in population.rows}
    if OPERATOR_PROTOCOL not in keys:
        return frozenset()
    ancestors = _inheritors(edges)
    wanted = _members(sources[OPERATOR_PROTOCOL.file], OPERATOR_PROTOCOL.qualname)
    return frozenset(
        k for k in keys if k == OPERATOR_PROTOCOL or OPERATOR_PROTOCOL in ancestors.get(k, ())
        or (bool(wanted) and wanted <= _members(sources[k.file], k.qualname)))


def _declared_sub(row: ClassRow) -> tuple[str, str | None, str]:
    d = row.declaration
    assert d is not None
    return d.kind, d.subkind, d.role


def check_kinds(population: Population, concepts: tuple[Concept, ...], edges: Iterable[Edge],
                sources: Mapping[str, str]) -> frozenset[KindLawDebt]:
    """The kind laws, each violation a ``KindLawDebt`` naming its law.

    ``head_arm_kind``: a head or arm row has its concept's kind. ``annex_record``: an annex row is a record.
    ``apparatus_row``: kind apparatus only on rows of a concept of kind apparatus. ``operator_subkind``: over the rows
    that are not interfaces, every class of the operator algebra has sub-kind operator. ``element_subkind``: every class
    that inherits ``Member`` has sub-kind member or covector. ``covector_subkind``: ``Covector`` and its inheritors have
    sub-kind covector.
    """
    edges = tuple(edges)
    kind_of = {c.id: c.kind for c in concepts}
    declared = _declared(population)
    ancestors = _inheritors(edges)
    hits = _operator_hits(population, edges, sources)
    out: set[KindLawDebt] = set()
    for k, row in declared.items():
        kind, sub, role = _declared_sub(row)
        d = row.declaration
        assert d is not None
        concept_kind = kind_of.get(d.concept)
        if role in ("head", "arm") and concept_kind is not None and kind != concept_kind:
            out.add(KindLawDebt("head_arm_kind", k))
        if role == "annex" and kind != "record":
            out.add(KindLawDebt("annex_record", k))
        if kind == "apparatus" and concept_kind is not None and concept_kind != "apparatus":
            out.add(KindLawDebt("apparatus_row", k))
        if role != "interface" and k in hits and sub != "operator":
            out.add(KindLawDebt("operator_subkind", k))
        if MEMBER_CLASS in ancestors.get(k, ()) and sub not in ("member", "covector"):
            out.add(KindLawDebt("element_subkind", k))
        if role != "interface" and (k == COVECTOR_CLASS or COVECTOR_CLASS in ancestors.get(k, ())) and sub != "covector":
            out.add(KindLawDebt("covector_subkind", k))
    return frozenset(out)


def check_off_algebra(population: Population, edges: Iterable[Edge], sources: Mapping[str, str]) -> frozenset[OffAlgebraDebt]:
    """The converse of Law 3: a class that is not an interface, has sub-kind operator and is not on the operator algebra."""
    hits = _operator_hits(population, tuple(edges), sources)
    return frozenset(OffAlgebraDebt(k) for k, row in _declared(population).items()
                     if _declared_sub(row)[1] == "operator" and _declared_sub(row)[2] != "interface" and k not in hits)


def check_covector_lineage(population: Population, edges: Iterable[Edge]) -> frozenset[CovectorOffLineageDebt]:
    """The converse of Law 6: a class that is not an interface, has sub-kind covector and is neither ``Covector`` nor its inheritor."""
    ancestors = _inheritors(tuple(edges))
    return frozenset(
        CovectorOffLineageDebt(k) for k, row in _declared(population).items()
        if _declared_sub(row)[1] == "covector" and _declared_sub(row)[2] != "interface"
        and k != COVECTOR_CLASS and COVECTOR_CLASS not in ancestors.get(k, ()))


# --- G5, G5L, G6, G7 -------------------------------------------------------------------------------------------------

def _judged(edges: Iterable[Edge], population: Population) -> list[tuple[Edge, str, str]]:
    """The class-to-class edges whose two ends are declared classes of two different concepts, with the two concepts."""
    declared = _declared(population)
    out = []
    for e in edges:
        a, b = declared.get(e.src), declared.get(e.dst)
        if a is not None and b is not None and a.declaration and b.declaration and a.declaration.concept != b.declaration.concept:
            out.append((e, a.declaration.concept, b.declaration.concept))
    return out


def _phase_violations(judged: list[tuple[Edge, str, str]], concepts: tuple[Concept, ...]) -> set[Edge]:
    phase = {c.id: c.phase for c in concepts}
    out = set()
    for e, a, b in judged:
        pa, pb = phase.get(a), phase.get(b)
        if pa is not None and pb is not None and rank(pb) > rank(pa):
            out.add(e)
    return out


def check_phase(edges: Iterable[Edge], population: Population, concepts: tuple[Concept, ...]) -> frozenset[EdgeDebt]:
    """A reference from a class of one concept to a class of another whose phase is later on ``PHASE_ORDER``.

    An end whose concept has no phase (an apparatus) exempts the edge; the converse, an earlier phase, is legal.
    """
    return frozenset(EdgeDebt(e.src, e.dst, e.kind) for e in _phase_violations(_judged(edges, population), concepts))


def check_layer(edges: Iterable[Edge], population: Population, forbidden: Mapping[str, Collection[str]]) -> frozenset[EdgeDebt]:
    """A class-to-class reference, of any kind, from a package that may not reference the package of its target."""
    keys = {r.key for r in population.rows}
    return frozenset(EdgeDebt(e.src, e.dst, e.kind) for e in edges
                     if e.src in keys and e.dst in keys and not may_reference(package(e.src), package(e.dst), forbidden))


def check_placement(population: Population, concepts: tuple[Concept, ...], forbidden: Mapping[str, Collection[str]]) -> frozenset[ClassDebt]:
    """The class declaring ``head`` for a concept sits in the package of the concept's layer; every other row sits in that
    package or in one that may reference it."""
    layer = {c.id: c.layer for c in concepts}
    out = set()
    for k, row in _declared(population).items():
        d = row.declaration
        assert d is not None
        want = layer.get(d.concept)
        if want is None:
            continue
        here = package(k)
        if (d.role == "head" and here != want) or (d.role != "head" and not may_reference(here, want, forbidden)):
            out.add(ClassDebt(k))
    return frozenset(out)


def check_cycles(edges: Iterable[Edge], population: Population, concepts: tuple[Concept, ...]) -> frozenset[CycleDebt]:
    """The cycles of the concept graph (T11 as amended): the strict inheritance half and the quotient half.

    The inheritance half: the concept graph of the judged ``inherits`` edges is acyclic. The quotient half: the graph of
    the judged edges that are not phase violations has its strongly connected components of at least two concepts as the
    cycles; each is a ``CycleDebt``.
    """
    judged = _judged(edges, population)
    violating = _phase_violations(judged, concepts)

    def graph(chosen: Iterable[tuple[Edge, str, str]]) -> dict[str, set[str]]:
        g: dict[str, set[str]] = collections.defaultdict(set)
        for _, a, b in chosen:
            g[a].add(b)
        return g

    inheritance = _tarjan(graph(j for j in judged if j[0].kind == "inherits"))
    quotient = _tarjan(graph(j for j in judged if j[0] not in violating))
    return frozenset(CycleDebt(c) for c in (*inheritance, *quotient))


# --- the debt lists, the wave rule, the unweld rule ------------------------------------------------------------------

def measure_twin(population: Population) -> frozenset[ClassDebt]:
    """The classes whose declared role is ``twin``."""
    return frozenset(ClassDebt(k) for k, r in _declared(population).items() if r.declaration and r.declaration.role == "twin")


def measure_weld(population: Population) -> frozenset[ClassDebt]:
    """The classes whose declared role is ``weld``, and the heads that carry a welded list."""
    return frozenset(ClassDebt(k) for k, r in _declared(population).items()
                     if r.declaration and (r.declaration.role == "weld" or (r.declaration.role == "head" and r.declaration.welded)))


def check_wave_rule(debt: Debt, concepts: tuple[Concept, ...], population: Population) -> frozenset:
    """The entries of the nine lists whose wave set is empty: an entry that no wave owns is debt nobody removes."""
    return frozenset(e for entries in debt.lists().values() for e in entries if not wave_of(e, concepts, population))


def _position(lines: list[str], row: int, col: int) -> int:
    """The character column of a UTF-8 byte offset of ``ast``."""
    return len(lines[row - 1].encode()[:col].decode())


_SKIPPED = {tokenize.COMMENT, tokenize.NL, tokenize.NEWLINE, tokenize.INDENT, tokenize.DEDENT}


def fingerprint(source: str, qualname: str) -> str:
    """The BLAKE2b digest of the token stream of a class statement, as ``"blake2b:" + hexdigest``.

    The lines from the first decorator (or the ``class`` line) to the last line of the class are tokenized; comment, newline,
    indent and dedent tokens are dropped and so is every string token inside a docstring (a string-constant expression
    statement first in the body of the class or of any function in it). A trailing declaration comment, a docstring edit and
    a comment edit leave it unchanged; one added statement changes it.
    """
    node = _class_node(source, qualname)
    lines = source.splitlines(keepends=True)
    start = min([node.lineno, *(d.lineno for d in node.decorator_list)])
    text = "".join(lines[start - 1:node.end_lineno])
    spans = []
    for n in ast.walk(node):
        if isinstance(n, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)) and n.body:
            first = n.body[0]
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) and isinstance(first.value.value, str):
                assert first.end_lineno is not None and first.end_col_offset is not None
                spans.append(((first.lineno, _position(lines, first.lineno, first.col_offset)),
                              (first.end_lineno, _position(lines, first.end_lineno, first.end_col_offset))))
    digest = hashlib.blake2b(digest_size=32)
    for tok in tokenize.generate_tokens(io.StringIO(text).readline):
        if tok.type in _SKIPPED:
            continue
        a = (tok.start[0] + start - 1, tok.start[1])
        z = (tok.end[0] + start - 1, tok.end[1])
        if tok.type == tokenize.STRING and any(s <= a and z <= e for s, e in spans):
            continue
        digest.update(f"{tokenize.tok_name[tok.type]}\x1f{tok.string}\x1e".encode())
    return "blake2b:" + digest.hexdigest()


def check_unweld(population: Population, seed: Iterable[SeedEntry],
                 sources: Mapping[str, str]) -> tuple[frozenset[ClassDebt], frozenset[ClassDebt]]:
    """``(uncovered, relabelled)``: a class leaves the twin or weld list only with code motion.

    Arm 1: every class declared ``twin``, ``weld`` or a welded head has a seed entry. Arm 2: for every seed entry whose
    class is declared none of those today, either no class with that key exists or its fingerprint differs from the seed's;
    otherwise the class was relabelled without code motion. A class that stays on a list may change freely.
    """
    seed = tuple(seed)
    listed = {d.key for d in measure_twin(population) | measure_weld(population)}
    seeded = {s.key for s in seed}
    uncovered = frozenset(ClassDebt(k) for k in listed if k not in seeded)
    present = {r.key for r in population.rows}
    relabelled = frozenset(
        ClassDebt(s.key) for s in seed
        if s.key not in listed and s.key in present and fingerprint(sources[s.key.file], s.key.qualname) == s.fingerprint)
    return uncovered, relabelled
