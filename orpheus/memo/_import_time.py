r"""What a run can observe of a module's import-time code: its statements, and which of them a run selects.

A module's import-time code is the statements its body runs when it is imported: the top-level statements and,
recursively, the statements of every class body, with every def body set aside (a def body runs only when called,
and the traced memo pins it by its own digest when it does). A run observes such a statement in one of two ways:

* **by name.** The statement binds a name (an assignment, a def, a class, an import), and code the run executed
  reads that name. The names a run read are its *entry names*: the names of every code object it started
  (:func:`code_names`), except the module-level code of first-party modules, which is read here statically, per
  statement. The selection is closed (:func:`closure`): a selected statement's reads are entry names too, so a
  constant derived from another (``B = A * 2``) selects both, and an ``ImportFrom`` follows a re-export.
* **by effect.** The statement does something no binding names: a call outside a short list of pure callees, a
  store into a subscript or an attribute, a ``del``, a ``raise`` or an ``assert``, a class with a base or a keyword
  that is not a builtin or a trusted import (its creation can run a first-party ``__init_subclass__``,
  ``__set_name__`` or metaclass), an import (the imported module's own effects). An effect is selected for every run that imported its module.

Within a class, the header (decorators, bases, keywords) is selected with any statement of its body, and a selected
header selects every non-def statement of its body and every statement binding a dunder: a field iterated by
``dataclasses.fields`` is read by no name, and the interpreter and C libraries read ``__hash__``,
``__match_args__`` or ``__array_priority__``, never Python code by name.

Names are not qualified by module, so every module binding a name the run read is selected: a collision
over-selects. The rule and the 42 fixtures it was attacked with are ruling 43 of
``.claude/plans/verification_tree_refinement.md``; what it cannot see is listed on
:doc:`/theory/verification/cone_store`.

This module imports nothing from ``orpheus``: the traced memo reads it, and it is one of the memo's own sources.
"""
from __future__ import annotations

import ast
import builtins
import copy
import re
import types
from collections.abc import Callable, Iterable, Iterator, Mapping
from typing import NamedTuple

DEF = (ast.FunctionDef, ast.AsyncFunctionDef)
#: The anonymous expressions, each with the name of the code object it compiles to. Each opens a scope of its own:
#: a name stored inside one binds nothing at module or class level, and a lambda's body does not run at import.
ANONYMOUS_NAME: Mapping[type[ast.AST], str] = {
    ast.Lambda: "<lambda>", ast.GeneratorExp: "<genexpr>", ast.ListComp: "<listcomp>",
    ast.SetComp: "<setcomp>", ast.DictComp: "<dictcomp>",
}
_IDENTIFIER = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")

#: Callees whose call at import does nothing but return a value, by their last name; pure only when the callee's
#: root is a builtin the module does not rebind, or a binding imported from one of :data:`_PURE_MODULES` (so a
#: first-party ``cache`` or ``field`` is an effect).
_PURE_CALLEES = frozenset({
    "frozenset", "tuple", "set", "dict", "list", "int", "float", "complex", "str", "bytes", "bool", "range", "len",
    "sorted", "reversed", "enumerate", "zip", "map", "filter", "min", "max", "sum", "abs", "round", "repr",
    "isinstance", "issubclass", "getattr", "hasattr", "object", "slice", "field", "dataclass", "TypeVar",
    "ParamSpec", "TypeVarTuple", "NewType", "namedtuple", "NamedTuple", "cast", "overload", "final", "abstractmethod",
    "total_ordering", "wraps", "partial", "cache", "lru_cache", "cached_property", "property", "staticmethod",
    "classmethod", "Path", "compile", "array", "asarray", "zeros", "ones", "full", "empty", "eye", "arange",
    "linspace", "logspace", "sqrt", "exp", "log", "finfo", "iinfo", "dtype", "float64", "int64", "complex128",
    "mark", "parametrize", "fixture", "param", "skipif", "xfail", "usefixtures", "Fraction", "Decimal", "auto",
    "itemgetter", "attrgetter", "MappingProxyType", "Symbol", "symbols", "IndexedBase", "Interval", "Rational",
    "Integer", "Function", "mpf", "mpc",
})
_PURE_MODULES = frozenset({
    "builtins", "dataclasses", "functools", "typing", "typing_extensions", "collections", "enum", "abc", "pathlib",
    "re", "math", "operator", "types", "fractions", "decimal", "numpy", "pytest", "sympy", "mpmath",
})
_BUILTINS = frozenset(dir(builtins))


class Statement(NamedTuple):
    """One import-time statement of a module, as the selection reads it.

    ``digest``: the digest of its normalised source, every def body read as ``pass`` (a class statement's own body is
    its member statements, not part of it). ``binds``: the names it binds. ``reads``: the names it reads at import,
    with the names an ``ImportFrom`` takes from its module and the identifier tokens of its strings. ``effect``:
    whether it does more than bind (selected for every run that imported the module). ``imports``: it is an import,
    whose reads join the closure only when it is selected by name. ``defines``: it is a def. ``owner``: the index of
    the class statement whose body holds it (``None``: the module's body). ``first`` and ``last``: its lines, its
    decorators included."""

    digest: str
    binds: tuple[str, ...]
    reads: tuple[str, ...]
    effect: bool
    imports: bool
    defines: bool
    owner: int | None
    first: int
    last: int

    @property
    def member(self) -> bool:
        """A non-def statement of a class body: selected with its class's header."""
        return self.owner is not None and not self.defines

    @property
    def dunder(self) -> bool:
        """A class-body statement binding a dunder: selected with its class's header."""
        return self.owner is not None and any(n.startswith("__") and n.endswith("__") for n in self.binds)


def string_names(constants: Iterable[object]) -> Iterator[str]:
    """The identifier tokens of every string among ``constants``, recursing through tuples and frozensets."""
    for value in constants:
        if isinstance(value, str):
            yield from _IDENTIFIER.findall(value)
        elif isinstance(value, (tuple, frozenset)):
            yield from string_names(value)


def code_names(code: types.CodeType) -> Iterator[str]:
    """The names a code object can read by name: its ``co_names`` (globals and attributes) and the identifier tokens
    of its string constants (``getattr(x, "name")``, ``attrgetter("a.b")``, a ``match`` class pattern's keywords)."""
    yield from code.co_names
    yield from string_names(code.co_consts)


# ── one statement ─────────────────────────────────────────────────────────────────────────────────────


def span(node: ast.AST) -> range:
    """The lines of a statement, its decorators included."""
    first = min([node.lineno, *(d.lineno for d in getattr(node, "decorator_list", []))])  # type: ignore[attr-defined]
    return range(first, node.end_lineno + 1)  # type: ignore[attr-defined]


def _evaluated(node: ast.FunctionDef | ast.AsyncFunctionDef | ast.ClassDef | ast.Lambda) -> list[ast.expr]:
    """What a def, a class or a lambda statement evaluates when it runs, beside its body: a def's decorators and
    defaults, a lambda's defaults, a class's decorators, bases and keywords (annotations are lazy, PEP 649)."""
    if isinstance(node, ast.ClassDef):
        return [*node.decorator_list, *node.bases, *(k.value for k in node.keywords)]
    defaults = [*node.args.defaults, *(d for d in node.args.kw_defaults if d is not None)]
    return defaults if isinstance(node, ast.Lambda) else [*node.decorator_list, *defaults]


def _annotations(node: ast.FunctionDef | ast.AsyncFunctionDef) -> list[ast.expr]:
    """A def's annotations, read (by name) when they are evaluated, not when the def runs."""
    arguments = node.args
    every = (*arguments.posonlyargs, *arguments.args, *arguments.kwonlyargs, arguments.vararg, arguments.kwarg)
    return [a.annotation for a in every if a is not None and a.annotation is not None] + \
        ([node.returns] if node.returns is not None else [])


def _reads_binds(statement: ast.stmt) -> tuple[set[str], set[str]]:
    """The names ``statement`` reads at import and the names it binds in the namespace whose body holds it. A def's
    body and a split class's body are not walked (the def's decorators, defaults and annotations are); a class nested
    in a compound statement is walked whole, its body's bindings being its own attributes, not the module's."""
    reads: set[str] = set()
    binds: set[str] = set()

    def walk(node: ast.AST, binding: bool) -> None:
        match node:
            case ast.FunctionDef() | ast.AsyncFunctionDef():
                if binding:
                    binds.add(node.name)
                for part in (*_evaluated(node), *_annotations(node)):
                    walk(part, False)
                return
            case ast.ClassDef():
                if binding:
                    binds.add(node.name)
                for part in _evaluated(node):
                    walk(part, False)
                if node is not statement:
                    for member in node.body:
                        walk(member, False)
                return
            case ast.Name():
                if isinstance(node.ctx, (ast.Store, ast.Del)):
                    if binding:
                        binds.add(node.id)
                else:
                    reads.add(node.id)
            case ast.Attribute():
                reads.add(node.attr)
            case ast.Constant() if isinstance(node.value, str):
                reads.update(string_names((node.value,)))
            case ast.ImportFrom():
                for alias in node.names:
                    binds.add(alias.asname or alias.name)
                    reads.add(alias.name)
            case ast.Import():
                for alias in node.names:
                    binds.add((alias.asname or alias.name).split(".")[0])
            case ast.ExceptHandler() if node.name and binding:
                binds.add(node.name)
            case ast.MatchAs() | ast.MatchStar() if node.name and binding:
                binds.add(node.name)
            case ast.MatchMapping() if node.rest and binding:
                binds.add(node.rest)
            case ast.NamedExpr():
                walk(node.target, True)  # a walrus binds in the enclosing namespace, even inside a comprehension
                walk(node.value, False)
                return
        for child in ast.iter_child_nodes(node):
            walk(child, binding and not isinstance(node, (*ANONYMOUS_NAME, ast.comprehension)))

    walk(statement, True)
    return reads, binds


def _trusted_roots(tree: ast.Module) -> frozenset[str]:
    """The names whose call can be pure in this module: a builtin the module never rebinds, or a name bound only by
    imports from one of :data:`_PURE_MODULES`."""
    imported: set[str] = set()
    rebound: set[str] = set()
    for node in ast.walk(tree):
        match node:
            case ast.ImportFrom():
                package = "" if node.level else (node.module or "").split(".")[0]
                for alias in node.names:
                    (imported if package in _PURE_MODULES else rebound).add(alias.asname or alias.name)
            case ast.Import():
                for alias in node.names:
                    package = alias.name.split(".")[0]
                    (imported if package in _PURE_MODULES else rebound).add((alias.asname or alias.name).split(".")[0])
            case ast.FunctionDef() | ast.AsyncFunctionDef() | ast.ClassDef():
                rebound.add(node.name)
            case ast.Name() if isinstance(node.ctx, (ast.Store, ast.Del)):
                rebound.add(node.id)
    return frozenset((imported | _BUILTINS) - rebound)


def _root(node: ast.expr) -> str | None:
    """The name an expression is reached from through attributes and subscripts alone (``np`` of ``np.zeros``,
    ``Generic`` of ``Generic[T]``); ``None`` when a call or anything else stands in the chain."""
    while isinstance(node, (ast.Attribute, ast.Subscript)):
        node = node.value
    return node.id if isinstance(node, ast.Name) else None


def _pure_callee(callee: ast.expr, trusted: frozenset[str]) -> bool:
    """Whether calling ``callee`` at import only returns a value: its root is trusted and either its last name is a
    pure callee (``np.zeros``; not ``np.random.default_rng(0).normal``) or it is a pytest mark (``pytest.mark.slow``,
    which attaches a mark and nothing else)."""
    if _root(callee) not in trusted:
        return False
    last = callee.attr if isinstance(callee, ast.Attribute) else callee.id if isinstance(callee, ast.Name) else None
    chain: list[str] = []
    while isinstance(callee, ast.Attribute):
        chain.append(callee.attr)
        callee = callee.value
    return last in _PURE_CALLEES or "mark" in chain[1:]


#: Pure callees that CALL a function they are given: ``map``'s and ``filter``'s first argument, every ``key=``.
_CALLS_ITS_ARGUMENT = frozenset({"map", "filter", "sorted", "min", "max"})


def _calls_untrusted(call: ast.Call, trusted: frozenset[str]) -> bool:
    """Whether a call to a pure higher-order callee calls a function whose root is not trusted (``map(register, X)``
    calls ``register``; a lambda's body is not walked, so a lambda counts as untrusted)."""
    callee = call.func
    last = callee.attr if isinstance(callee, ast.Attribute) else callee.id if isinstance(callee, ast.Name) else None
    if last not in _CALLS_ITS_ARGUMENT:
        return False
    called = [k.value for k in call.keywords if k.arg == "key"]
    if last in ("map", "filter") and call.args:
        called.append(call.args[0])
    return any(_root(f) not in trusted for f in called)


def _import_time_nodes(statement: ast.stmt, split: bool) -> Iterator[ast.AST]:
    """The nodes evaluated when ``statement`` runs at import: a def's decorators and defaults (annotations are lazy),
    a split class's decorators, bases and keywords, never a def's or a lambda's body."""
    def walk(node: ast.AST) -> Iterator[ast.AST]:
        yield node
        match node:
            case ast.FunctionDef() | ast.AsyncFunctionDef() | ast.Lambda():
                children: list[ast.AST] = list(_evaluated(node))
            case ast.ClassDef() if split and node is statement:
                children = list(_evaluated(node))
            case _:
                children = list(ast.iter_child_nodes(node))
        for child in children:
            yield from walk(child)
    return walk(statement)


def _effect(statement: ast.stmt, trusted: frozenset[str], split: bool) -> bool:
    """Whether ``statement`` does more than bind names (the module docstring's list)."""
    if isinstance(statement, (ast.Import, ast.ImportFrom)):
        return True
    for node in _import_time_nodes(statement, split):
        match node:
            case ast.Expr() if not isinstance(node.value, ast.Constant):
                return True
            case ast.Delete() | ast.Raise() | ast.Assert() | ast.Import() | ast.ImportFrom() | ast.Global():
                return True
            case ast.Subscript() | ast.Attribute() if isinstance(node.ctx, (ast.Store, ast.Del)):
                return True
            case ast.ClassDef() if any(_root(e) not in trusted for e in (*node.bases, *(k.value for k in node.keywords))):
                return True
            case ast.Call() if not _pure_callee(node.func, trusted) or _calls_untrusted(node, trusted):
                return True
            case ast.FunctionDef() | ast.AsyncFunctionDef() | ast.ClassDef():
                # a decorator is called with the def or class, whether or not it is written as a call
                if any(not _pure_callee(d.func if isinstance(d, ast.Call) else d, trusted) for d in node.decorator_list):
                    return True
    return False


class _Hollowed(ast.NodeTransformer):
    """Every def body read as ``pass``: the def's own pin digests the body."""

    def visit_FunctionDef(self, node: ast.FunctionDef) -> ast.AST:
        node.body = [ast.Pass()]
        return node

    visit_AsyncFunctionDef = visit_FunctionDef  # type: ignore[assignment]


def _hollowed(statement: ast.stmt) -> ast.stmt:
    """``statement`` as its digest reads it, the parsed tree untouched: every def body read as ``pass``, and a split
    class's body too (its member statements are statements of their own)."""
    if isinstance(statement, (*DEF, ast.ClassDef)):
        clone = copy.copy(statement)  # decorators, bases, defaults and annotations hold no def
        clone.body = [ast.Pass()]
        return clone
    if not any(isinstance(node, DEF) for node in ast.walk(statement)):
        return statement
    return _Hollowed().visit(copy.deepcopy(statement))


def statements(tree: ast.Module, digest: Callable[[ast.AST], str]) -> tuple[Statement, ...]:
    """The import-time statements of a parsed module (docstrings already stripped), in source order: the module's
    body, and after each class statement of it, recursively, its body's statements. ``digest`` normalises a node (the
    traced memo's, so a statement and a def are digested alike)."""
    trusted = _trusted_roots(tree)
    found: list[Statement] = []

    def add(body: list[ast.stmt], owner: int | None) -> None:
        for node in body:
            reads, binds = _reads_binds(node)
            index = len(found)
            lines = span(node)
            found.append(Statement(digest(_hollowed(node)), tuple(sorted(binds)), tuple(sorted(reads)),
                                   _effect(node, trusted, split=True), isinstance(node, (ast.Import, ast.ImportFrom)),
                                   isinstance(node, DEF), owner, lines.start, lines.stop - 1))
            if isinstance(node, ast.ClassDef):
                add(node.body, index)

    add(tree.body, None)
    return tuple(found)


def holding(statements: tuple[Statement, ...], line: int) -> Statement | None:
    """The innermost statement whose lines hold ``line`` (``None``: none). A class's member statements follow it, so
    the last holder is the innermost."""
    held = [s for s in statements if s.first <= line <= s.last]
    return held[-1] if held else None


# ── a run's selection ─────────────────────────────────────────────────────────────────────────────────


def selection(statements: tuple[Statement, ...], names: frozenset[str]) -> frozenset[int]:
    """The statements of one module a run with entry names ``names`` observes: every effect, every statement binding
    one of ``names``, and the class rules (a header with any statement of its body; a selected header's members and
    dunders), to a fixpoint."""
    chosen = {i for i, s in enumerate(statements) if s.effect or not names.isdisjoint(s.binds)}
    while True:
        grown = {owner for i in chosen if (owner := statements[i].owner) is not None} - chosen
        grown |= {i for i, s in enumerate(statements)
                  if i not in chosen and s.owner in chosen and (s.member or s.dunder)}
        if not grown:
            return frozenset(chosen)
        chosen |= grown


def closure(modules: Mapping[str, tuple[Statement, ...]], seed: Iterable[str]) -> frozenset[str]:
    """The entry names of a run over its first-party modules: ``seed`` (the names its executed code read), grown by
    the reads of every statement it selects until nothing grows. An import's reads grow it only when the import is
    selected by name (a re-export the run reads), not as an effect."""
    names = frozenset(seed)
    while True:
        grown: set[str] = set()
        for rows in modules.values():
            for i in selection(rows, names):
                s = rows[i]
                if not (s.imports and names.isdisjoint(s.binds)):
                    grown.update(s.reads)
        if grown <= names:
            return names
        names |= grown


def observed_digests(statements: tuple[Statement, ...], names: frozenset[str]) -> tuple[str, ...]:
    """The digests, in source order, of the statements of one module a run with entry names ``names`` observes."""
    return tuple(statements[i].digest for i in sorted(selection(statements, names)))
