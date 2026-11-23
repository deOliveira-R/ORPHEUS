"""The reference pass: which class of ``orpheus/`` refers to which, read from the syntax tree alone.

An edge is a reference from the class that owns a piece of code to a class named in it. The kinds, strongest first:
``inherits`` (a base of the class header), ``top`` (a name bound by a module-level import, also inside a module-level
``try`` or a module-level ``if`` that is not a type-checking guard, or defined in the same module), ``late`` (a name
imported inside a function body of the owner) and ``typing`` (a name imported under ``if TYPE_CHECKING:``, or a string
that is one identifier). Of the edges between one pair of classes the strongest kind is kept.

A module-level class owns every reference in its body, nested definitions included; a class nested in a class (``Outer.Inner``) is a
class of the population too, owns its own body besides, and is the end of a reference spelled ``Outer.Inner``; a module-level function, assignment
or annotated assignment is owned by the pseudo-owner of its file (``FunctionEdge``), and the module's ``__all__`` list is
no reference, being strings. A generic subscript is not a base: ``class A(Base[Arg])`` inherits ``Base`` and only refers
to ``Arg``. A module-level union alias of plain names (``A | B``, ``Union[A, B]``, ``A | None``) is expanded: a name that
is such an alias yields one reference to each leaf class. Nothing is imported and nothing runs.
"""

from __future__ import annotations

import ast
import collections
import dataclasses
import functools
import re
from collections.abc import Mapping

from .declaration import walk_classes
from .population import ClassKey, Population

EDGE_KINDS = ("inherits", "top", "late", "typing")
_RANK = {kind: rank for rank, kind in enumerate(EDGE_KINDS)}
_IDENT = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
_ALIAS_DEPTH = 4
_CHASE_DEPTH = 6


@dataclasses.dataclass(frozen=True)
class Edge:
    """A reference from class ``src`` to class ``dst`` of one kind in ``EDGE_KINDS``."""

    src: ClassKey
    dst: ClassKey
    kind: str

    def __post_init__(self) -> None:
        if self.kind not in EDGE_KINDS:
            raise ValueError(f"edge kind {self.kind!r} is not one of {EDGE_KINDS}")

    def __str__(self) -> str:
        return f"{self.src} -> {self.dst} [{self.kind}]"


@dataclasses.dataclass(frozen=True)
class FunctionEdge:
    """A reference from a module-level function or assignment of ``src_file`` to class ``dst``; the checks never judge it."""

    src_file: str
    dst: ClassKey
    kind: str

    def __post_init__(self) -> None:
        if self.kind not in EDGE_KINDS:
            raise ValueError(f"edge kind {self.kind!r} is not one of {EDGE_KINDS}")


def is_type_checking_guard(test: ast.expr) -> bool:
    """``if TYPE_CHECKING:`` or ``if typing.TYPE_CHECKING:``, exactly; ``if not TYPE_CHECKING:`` is not a guard."""
    return (isinstance(test, ast.Name) and test.id == "TYPE_CHECKING") or (
        isinstance(test, ast.Attribute) and test.attr == "TYPE_CHECKING")


def module_name(file: str) -> str:
    """The dotted module of a file: ``orpheus/a/b.py`` is ``orpheus.a.b`` and ``orpheus/a/__init__.py`` is ``orpheus.a``."""
    parts = file.removesuffix(".py").split("/")
    return ".".join(parts[:-1] if parts[-1] == "__init__" else parts)


def _leaves(x: ast.expr) -> list[str] | None:
    """The names of a union spelled ``A | B``, ``Union[A, B]`` or ``A | None``; ``None`` if a leaf is not a plain name."""
    if isinstance(x, ast.BinOp) and isinstance(x.op, ast.BitOr):
        left, right = _leaves(x.left), _leaves(x.right)
        return None if left is None or right is None else left + right
    if isinstance(x, ast.Subscript) and ast.unparse(x.value).split(".")[-1] == "Union":
        elements = x.slice.elts if isinstance(x.slice, ast.Tuple) else [x.slice]
        out: list[str] = []
        for element in elements:
            leaves = _leaves(element)
            if leaves is None:
                return None
            out += leaves
        return out
    if isinstance(x, ast.Name):
        return [x.id]
    if isinstance(x, ast.Constant) and x.value is None:
        return []
    return None


# A resolved name: (module, class name). A binding: (source module, imported symbol or None, kind).
_Cls = tuple[str, str]
_Binding = tuple[str, "str | None", str]


class _Pass:
    """The reference pass over one mapping of sources."""

    def __init__(self, sources: Mapping[str, str]) -> None:
        self.files = {module_name(f): f for f in sorted(sources)}
        try:
            self.trees = {m: ast.parse(sources[f]) for m, f in self.files.items()}
        except SyntaxError as error:
            raise ValueError(f"the reference pass cannot parse a source: {error}") from error
        self.defs = {m: {n.name for n in ast.walk(t) if isinstance(n, ast.ClassDef)} for m, t in self.trees.items()}
        # Classes nested in a class, by qualified name: ``Outer.Inner`` is a class of its own, an edge end and an owner.
        self.nested = {m: {q: c for q, c in walk_classes(t) if "." in q and "<locals>" not in q} for m, t in self.trees.items()}
        self.alias: dict[str, dict[str, list[str]]] = {}
        for m, tree in self.trees.items():
            for node in tree.body:
                target = value = None
                if isinstance(node, ast.Assign) and len(node.targets) == 1:
                    target, value = node.targets[0], node.value
                elif isinstance(node, ast.AnnAssign):
                    target, value = node.target, node.value
                elif isinstance(node, ast.TypeAlias):
                    target, value = node.name, node.value
                if isinstance(target, ast.Name) and value is not None:
                    leaves = _leaves(value)
                    if leaves is not None and len(leaves) >= 2:
                        self.alias.setdefault(m, {})[target.id] = leaves
        self.table = {m: self._module_table(m) for m in self.files}

    # --- name resolution -------------------------------------------------------------------------------------

    def _is_package(self, m: str) -> bool:
        return self.files[m].endswith("__init__.py")

    def _resolve_from(self, m: str, node: ast.ImportFrom) -> str:
        if node.level:
            base = m.split(".") if self._is_package(m) else m.split(".")[:-1]
            base = base[: len(base) - (node.level - 1)] if node.level > 1 else base
            return ".".join(base + ([node.module] if node.module else []))
        return node.module or ""

    def _module_table(self, m: str) -> dict[str, _Binding]:
        """The names bound by the imports of the module level: name -> (module, symbol or None, kind)."""
        table: dict[str, _Binding] = {}

        def add(node: ast.stmt, kind: str) -> None:
            if isinstance(node, ast.ImportFrom):
                source = self._resolve_from(m, node)
                for a in node.names:
                    table[a.asname or a.name] = (source, a.name, kind)
            elif isinstance(node, ast.Import):
                for a in node.names:
                    if a.asname:
                        table[a.asname] = (a.name, None, kind)
                    else:
                        table[a.name.split(".")[0]] = (a.name.split(".")[0], None, kind)

        def imports_below(node: ast.stmt, kind: str) -> None:
            for k in ast.walk(node):
                if isinstance(k, (ast.Import, ast.ImportFrom)):
                    add(k, kind)

        for node in self.trees[m].body:
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                add(node, "top")
            elif isinstance(node, ast.If):
                imports_below(node, "typing" if is_type_checking_guard(node.test) else "top")
            elif isinstance(node, ast.Try):
                imports_below(node, "top")
        return table

    def _chase(self, module: str, symbol: str, depth: int = 0) -> _Cls | None:
        """The class a ``(module, symbol)`` import names, followed through re-exports, or ``None``."""
        if depth > _CHASE_DEPTH or module not in self.files:
            return None
        if symbol in self.defs[module]:
            return (module, symbol)
        if symbol in self.table[module]:
            source, imported, _ = self.table[module][symbol]
            return None if imported is None else self._chase(source, imported, depth + 1)
        return None

    def _target(self, m: str, local: dict[str, _Binding], name: str, chain: tuple[str, ...] = ()) -> tuple[_Cls, str] | None:
        """The class a name (with an attribute chain) used in module ``m`` refers to, and the kind of its binding."""
        if name in local:
            source, symbol, kind = local[name]
        elif name in self.table[m]:
            source, symbol, kind = self.table[m][name]
        elif name in self.defs[m]:
            return ((m, name), "top")
        else:
            return None
        if symbol is not None:
            found = self._chase(source, symbol)
            if found:
                return (found, kind)
            submodule = f"{source}.{symbol}"
            if submodule in self.files and chain:
                found = self._chase(submodule, chain[0])
                if found:
                    return (found, kind)
            return None
        if source in self.files and chain:
            found = self._chase(source, chain[0])
            if found:
                return (found, kind)
        return None

    def _alias_targets(self, m: str, local: dict[str, _Binding], name: str, depth: int = 0) -> list[tuple[_Cls, str]]:
        """A name that is a module-level union alias, here or imported: the classes of its leaves, each with its kind."""
        if depth > _ALIAS_DEPTH:
            return []
        if name in self.alias.get(m, {}) and name not in self.defs[m]:
            out: list[tuple[_Cls, str]] = []
            for leaf in self.alias[m][name]:
                found = self._target(m, {}, leaf)
                if found:
                    out.append((found[0], "top"))
                else:
                    out += self._alias_targets(m, {}, leaf, depth + 1)
            return out
        binding = local.get(name) or self.table[m].get(name)
        if binding and binding[1] is not None and binding[0] in self.files:
            return [(cls, binding[2]) for cls, _ in self._alias_targets(binding[0], {}, binding[1], depth + 1)]
        return []

    # --- the scan --------------------------------------------------------------------------------------------

    def _scan(self, owner: _Cls, m: str, node: ast.AST, local: dict[str, _Binding], edges: collections.Counter) -> None:
        for n in ast.walk(node):
            if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)):
                for k in ast.walk(n):
                    if isinstance(k, ast.ImportFrom):
                        source = self._resolve_from(m, k)
                        for a in k.names:
                            local[a.asname or a.name] = (source, a.name, "late")
                            hit = self._chase(source, a.name)  # the import is itself a reference to the class it binds
                            if hit and hit != owner:
                                edges[(owner, hit, "late")] += 1
        for n in ast.walk(node):
            found = None
            if isinstance(n, ast.Name):
                found = self._target(m, local, n.id)
            elif isinstance(n, ast.Attribute) and isinstance(n.value, ast.Name):
                found = self._target(m, local, n.value.id, (n.attr,))
            elif isinstance(n, ast.Constant) and isinstance(n.value, str) and len(n.value) < 200:
                for word in set(_IDENT.findall(n.value)):
                    if word in self.defs[m] or word in self.table[m] or word in local:
                        hit = self._target(m, local, word)
                        stripped = n.value.strip()
                        for char in "'\" |[],.":
                            stripped = stripped.replace(char, "")
                        if (hit and word[0].isupper() and stripped.isidentifier()) or (hit and n.value.strip() == word):
                            if hit[0] != owner:
                                edges[(owner, hit[0], "typing")] += 1
                continue
            if found and found[0] != owner:
                edges[(owner, found[0], found[1])] += 1
            if found and isinstance(n, ast.Attribute):  # ``Outer.Inner`` refers to the nested class too
                inner = (found[0][0], f"{found[0][1]}.{n.attr}")
                if inner[1] in self.nested.get(inner[0], ()) and inner != owner:
                    edges[(owner, inner, found[1])] += 1
            if found is None and isinstance(n, ast.Name):
                for cls, kind in self._alias_targets(m, local, n.id):
                    if cls != owner:
                        edges[(owner, cls, kind)] += 1

    def _inherits(self, owner: _Cls, m: str, cls: ast.ClassDef, counted: collections.Counter) -> None:
        for base in cls.bases:
            base = base.value if isinstance(base, ast.Subscript) else base
            for k in ast.walk(base):
                found = None
                if isinstance(k, ast.Name):
                    found = self._target(m, {}, k.id)
                elif isinstance(k, ast.Attribute) and isinstance(k.value, ast.Name):
                    found = self._target(m, {}, k.value.id, (k.attr,))
                if found and found[0] != owner:
                    counted[(owner, found[0], "inherits")] += 1

    def edges(self) -> dict[tuple[_Cls, _Cls], str]:
        """The strongest kind of every (owner, class) pair; a function owner is ``(module, '<functions>')``."""
        counted: collections.Counter = collections.Counter()
        for m, tree in self.trees.items():
            for cls in (n for n in tree.body if isinstance(n, ast.ClassDef)):
                owner = (m, cls.name)
                self._inherits(owner, m, cls, counted)
                self._scan(owner, m, cls, {}, counted)
                for qualname, inner in self.nested[m].items():
                    if qualname.startswith(cls.name + "."):  # a nested class owns its own body, beside the class that holds it
                        self._inherits((m, qualname), m, inner, counted)
                        self._scan((m, qualname), m, inner, {}, counted)
            functions = (m, "<functions>")
            for node in tree.body:
                if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "__all__" for t in node.targets):
                    continue
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Assign, ast.AnnAssign)):
                    self._scan(functions, m, node, {}, counted)
        best: dict[tuple[_Cls, _Cls], str] = {}
        for (src, dst, kind), _ in counted.items():
            if (src, dst) not in best or _RANK[kind] < _RANK[best[(src, dst)]]:
                best[(src, dst)] = kind
        return best


@functools.lru_cache(maxsize=4)
def _pass_of(items: tuple[tuple[str, str], ...]) -> _Pass:
    return _Pass(dict(items))


def _pass(sources: Mapping[str, str]) -> _Pass:
    return _pass_of(tuple(sorted(sources.items())))


def _keys(population: Population) -> dict[str, ClassKey]:
    return {row.key.dotted: row.key for row in population.rows}


def read_edges(sources: Mapping[str, str], population: Population) -> tuple[Edge, ...]:
    """The class-to-class edges of ``sources``, sorted; an edge with an end outside the population is dropped."""
    keys = _keys(population)
    out = []
    for (src, dst), kind in _pass(sources).edges().items():
        if src[1] == "<functions>":
            continue
        a, b = keys.get(f"{src[0]}.{src[1]}"), keys.get(f"{dst[0]}.{dst[1]}")
        if a is not None and b is not None and a != b:
            out.append(Edge(a, b, kind))
    return tuple(sorted(out, key=lambda e: (e.src, e.dst, e.kind)))


def read_function_edges(sources: Mapping[str, str], population: Population) -> tuple[FunctionEdge, ...]:
    """The references owned by module-level functions and assignments, sorted; the checks count them and do not judge them."""
    keys = _keys(population)
    pass_ = _pass(sources)
    out = []
    for (src, dst), kind in pass_.edges().items():
        if src[1] != "<functions>":
            continue
        b = keys.get(f"{dst[0]}.{dst[1]}")
        if b is not None:
            out.append(FunctionEdge(pass_.files[src[0]], b, kind))
    return tuple(sorted(out, key=lambda e: (e.src_file, e.dst, e.kind)))


def alias_leaves(sources: Mapping[str, str], file: str, alias: str) -> tuple[ClassKey, ...]:
    """The classes that the module-level union alias ``alias`` of ``file`` names, as the pass resolved them through imports.

    Empty when ``file`` is not in ``sources`` or ``alias`` is not a union alias of at least two names there.
    """
    pass_ = _pass(sources)
    module = module_name(file)
    if module not in pass_.files or alias not in pass_.alias.get(module, {}) or alias in pass_.defs[module]:
        return ()
    leaves = ((pass_.files[m], name) for (m, name), _ in pass_._alias_targets(module, {}, alias))
    return tuple(ClassKey(*leaf) for leaf in dict.fromkeys(leaves))
