"""The independence legs of an L4 corroboration harness: an old family compared with the new code that replaces it.

A corroboration row compares an old spelling with its replacement during a
migration (``vv-principles``: L4, no correctness content; its value is that a
move which changes an answer shows up). The comparison means something only
while the two sides are independent, so every row first asserts it, in two
legs, and turns RED with "delete it in the migration commit" when either
fails (``retirement-audit`` D.14: once the old side reaches the new code the
row compares the new code with itself through a facade).

**The new side** is a :class:`NewSide`: the module prefixes of the new code,
the entry points the runtime leg counts (``Class.method``, each class defined
in exactly one of the new modules), and the label the refusal names.

**The static leg** (:func:`assert_independent`) reads each old module's own
source by AST and collects every dotted name it can reach without running it:

- an import, absolute or relative (resolved against the module's package),
  late, inside a function, or under ``TYPE_CHECKING`` alike;
- an attribute chain on a name an import bound (``continuous.characteristic.X``
  after ``from orpheus.derivations import continuous``), and the same chain
  spelled ``getattr(<bound name>, "<literal>")``;
- a string constant in an expression position that spells a new module, alone
  or before a ``:`` (``importlib.import_module("<new module>")``,
  ``sys.modules["<new module>"]``, an entry-point string
  ``"<new module>:Name"``); docstrings are skipped, since naming a module in
  prose imports nothing.

A name matches when it is a new module, lies under one, or is a name a parent
package of a new module re-exports from it (read from the parent's namespace
at check time, so a re-export added later is caught without editing a list).

A transitive leg (``assert_closure_independent``, walking every first-party
module the old modules import) served the characteristic reference's
corroboration file and was deleted with it in step (e2) of the
characteristic-reference campaign (the user's ruling of 2026-10-10: an
uncalled helper is deleted); the per-module leg is the one P0's file uses.

What the static legs still miss, by construction (qa ``qa_c/static_probe.py``
measured the shapes): a module name COMPUTED at run time (an f-string, a
concatenation, ``getattr(package, "char" + "acteristic")``, a ``getattr`` whose
first argument is not a name an import bound); a module reached only through
an attribute chain on a module object obtained at run time
(``sys.modules[...]`` with a computed key, a function's return value);
``exec``/``eval``; an object of the new code handed IN to the old code (an
argument, a registry entry, a pickled value) and called there; and anything a
third-party or C module does. The runtime leg exists for those.

**The runtime leg** (:func:`spy`, :func:`without`) counts, while the old
spelling is evaluated, every code object of the new modules that starts
(``sys.monitoring``, ``PY_START`` set per code object with
``set_local_events``, so no code outside the new modules pays for it), plus
the named entry points, so that a refusal names what was reached at the door
as well as the function that ran. The old evaluation runs under
:func:`orpheus.numerics.traced_memo.bypass`: a ``@traced_memo`` call otherwise
runs in a fresh interpreter on a miss and runs nothing on a hit, so an
in-process spy would read 0 either way. What it cannot see: code of the new
modules that runs in another process, and C code the new modules call.
It reads only what runs inside its window: an old value computed earlier and
cached in this process (a ``functools.cache``, a ``cached_property`` on a
shared object) passes it vacuously, so a row evaluates its old side afresh
(``[M]`` 2026-10-09: with the old Garcia solve cached across two rows, the
second row stayed green behind the first row's refusal).
"""
from __future__ import annotations

import ast
import importlib
import importlib.util
import inspect
import pkgutil
import sys
import types
from collections.abc import Callable, Iterable, Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any, TypeVar

import pytest

from orpheus.numerics.traced_memo import bypass

T = TypeVar("T")

_DELETE = "delete it in the migration commit (retirement-audit D.14)"


@dataclass(frozen=True)
class NewSide:
    """The new code of a migration: its module prefixes, the entry points the runtime leg counts, and its name in a refusal."""

    modules: tuple[str, ...]
    entry_points: tuple[str, ...]
    label: str


def family(package: str) -> tuple[str, ...]:
    """A package and every module under it, by name (``pkgutil.walk_packages``)."""
    root = importlib.import_module(package)
    return (package, *(info.name for info in pkgutil.walk_packages(root.__path__, package + ".")))


def _under(name: str, prefixes: Iterable[str]) -> bool:
    return any(name == p or name.startswith(p + ".") for p in prefixes)


def _docstrings(tree: ast.Module) -> set[int]:
    """The ``id`` of every docstring constant in ``tree``."""
    found = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)) and node.body:
            first = node.body[0]
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) and isinstance(first.value.value, str):
                found.add(id(first.value))
    return found


def _dotted(node: ast.expr) -> list[str] | None:
    """``a.b.c`` as ``["a", "b", "c"]``, or None when the chain does not end in a name."""
    parts: list[str] = []
    while isinstance(node, ast.Attribute):
        parts.append(node.attr)
        node = node.value
    if not isinstance(node, ast.Name):
        return None
    parts.append(node.id)
    return parts[::-1]


def references(module_name: str) -> dict[str, set[str]]:
    """The dotted names a module's own source reaches without running, by shape: ``import``, ``attribute``, ``string``."""
    module = importlib.import_module(module_name)
    return references_in(inspect.getsource(module), module.__package__ or module_name)


def references_in(source: str, package: str) -> dict[str, set[str]]:
    """:func:`references` of a source text whose relative imports resolve against ``package``."""
    tree = ast.parse(source)
    found: dict[str, set[str]] = {"import": set(), "attribute": set(), "string": set()}
    bound: dict[str, str] = {}                                           # a local name -> the module path it binds
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            base = importlib.util.resolve_name("." * node.level + (node.module or ""), package) if node.level else node.module
            if base is None:
                continue
            found["import"].add(base)
            for alias in node.names:
                found["import"].add(f"{base}.{alias.name}")
                bound[alias.asname or alias.name] = f"{base}.{alias.name}"
        elif isinstance(node, ast.Import):
            for alias in node.names:
                found["import"].add(alias.name)
                if alias.asname:
                    bound[alias.asname] = alias.name
                else:
                    head = alias.name.split(".")[0]
                    bound[head] = head

    def resolved(expression: ast.expr) -> str | None:
        parts = _dotted(expression)
        return ".".join([bound[parts[0]], *parts[1:]]) if parts and parts[0] in bound else None

    docstrings = _docstrings(tree)
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute):
            if (name := resolved(node)) is not None:
                found["attribute"].add(name)
        elif (isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "getattr"
              and len(node.args) >= 2 and isinstance(node.args[1], ast.Constant) and isinstance(node.args[1].value, str)):
            if (name := resolved(node.args[0])) is not None:
                found["attribute"].add(f"{name}.{node.args[1].value}")
        elif isinstance(node, ast.Constant) and isinstance(node.value, str) and id(node) not in docstrings:
            found["string"].add(node.value.split(":")[0].strip())
    return found


def _reexported(new: NewSide) -> set[str]:
    """Every ``parent.name`` a parent package of a new module binds to an object of the new code, read at check time."""
    names = set()
    for module in new.modules:
        parts = module.split(".")
        for depth in range(1, len(parts)):
            parent_name = ".".join(parts[:depth])
            if _under(parent_name, new.modules):
                continue
            parent = importlib.import_module(parent_name)
            for name, value in vars(parent).items():
                origin = getattr(value, "__module__", None) or getattr(value, "__name__", None)
                if isinstance(origin, str) and _under(origin, new.modules):
                    names.add(f"{parent_name}.{name}")
    return names


def _hits(names: Iterable[str], new: NewSide, reexported: set[str]) -> list[str]:
    return sorted(n for n in names if _under(n, new.modules) or _under(n, reexported))


def _refusal(subject: str, found: Mapping[str, set[str]], new: NewSide, reexported: set[str]) -> str | None:
    hits = {shape: h for shape, names in found.items() if (h := _hits(names, new, reexported))}
    if not hits:
        return None
    return (f"{subject} now imports {new.label} ({hits}): this corroboration row compares {new.label} "
            f"with itself; {_DELETE}")


def assert_independent(old_modules: Iterable[str], new: NewSide) -> None:
    """The static leg: no old module's own source names the new code (else the row is a tautology; delete it)."""
    reexported = _reexported(new)
    for module_name in old_modules:
        if (message := _refusal(module_name, references(module_name), new, reexported)) is not None:
            raise AssertionError(message)


def assert_source_independent(subject: str, source: str, package: str, new: NewSide) -> None:
    """:func:`assert_independent` on a source text (its relative imports resolved against ``package``): the leg's control."""
    if (message := _refusal(subject, references_in(source, package), new, _reexported(new))) is not None:
        raise AssertionError(message)


def _entry_class(new: NewSide, class_name: str) -> type:
    """The one class named ``class_name`` defined in a new module."""
    owners = set()
    for prefix in new.modules:
        for module_name in family(prefix) if hasattr(importlib.import_module(prefix), "__path__") else (prefix,):
            value = getattr(importlib.import_module(module_name), class_name, None)
            if isinstance(value, type) and _under(value.__module__, new.modules):
                owners.add(value)
    if len(owners) != 1:
        raise LookupError(f"the entry point's class {class_name} names {len(owners)} classes of {new.label}, not one")
    return owners.pop()


def _counted(original: Any, counts: dict[str, int], key: str) -> Any:
    """``original`` (a function, a classmethod, a staticmethod or any descriptor) with every call counted under ``key``."""
    if isinstance(original, classmethod):
        function = original.__func__

        def counted_class(cls: type, *args: Any, **kwargs: Any) -> Any:
            counts[key] += 1
            return function(cls, *args, **kwargs)

        return classmethod(counted_class)
    if isinstance(original, staticmethod):
        static = original.__func__

        def counted_static(*args: Any, **kwargs: Any) -> Any:
            counts[key] += 1
            return static(*args, **kwargs)

        return staticmethod(counted_static)

    def counted(self: Any, *args: Any, **kwargs: Any) -> Any:
        counts[key] += 1
        return original.__get__(self, type(self))(*args, **kwargs)

    return counted


def _code_objects(new: NewSide) -> dict[types.CodeType, str]:
    """Every code object defined in the new modules (functions, methods, descriptors' functions, nested code), by name."""
    codes: dict[types.CodeType, str] = {}

    def add(code: types.CodeType, module_name: str) -> None:
        if code in codes:
            return
        codes[code] = f"{module_name.rsplit('.', 1)[-1]}.{code.co_qualname}"
        for constant in code.co_consts:
            if isinstance(constant, types.CodeType):
                add(constant, module_name)

    def visit(value: Any, module_name: str, seen: set[int]) -> None:
        if id(value) in seen:
            return
        seen.add(id(value))
        if isinstance(value, (classmethod, staticmethod)):
            value = value.__func__
        if isinstance(value, property):
            for accessor in (value.fget, value.fset, value.fdel):
                if accessor is not None:
                    visit(accessor, module_name, seen)
            return
        for attribute in ("func", "__wrapped__"):                          # cached_property, a traced memo
            inner = getattr(value, attribute, None) if not isinstance(value, type) else None
            if callable(inner):
                visit(inner, module_name, seen)
        if isinstance(value, types.FunctionType) and value.__module__ == module_name:
            add(value.__code__, module_name)
        elif isinstance(value, type) and value.__module__ == module_name:
            for member in vars(value).values():
                visit(member, module_name, seen)

    for prefix in new.modules:
        names = family(prefix) if hasattr(importlib.import_module(prefix), "__path__") else (prefix,)
        for module_name in names:
            seen: set[int] = set()
            for value in vars(importlib.import_module(module_name)).values():
                visit(value, module_name, seen)
    return codes


@contextmanager
def spy(monkeypatch: pytest.MonkeyPatch, new: NewSide) -> Iterator[dict[str, int]]:
    """A counting spy on the new side: its named entry points (``Class.method``) and every code object of its modules.

    The dict maps each entry point, and each code object of the new modules
    that has started (``module.qualname``), to its call count. Each entry point
    must be defined on its class itself (patching an inherited method would
    count calls the class never routes through it).
    """
    codes = _code_objects(new)                                           # before the entry points are wrapped
    counts: dict[str, int] = {}
    for entry in new.entry_points:
        class_name, name = entry.split(".")
        cls = _entry_class(new, class_name)
        if name not in vars(cls):
            raise LookupError(f"{entry} is not defined on {cls.__qualname__} itself")
        counts[entry] = 0
        monkeypatch.setattr(cls, name, _counted(vars(cls)[name], counts, entry))
    monitoring = sys.monitoring
    tool = next((t for t in (5, 4, 3) if monitoring.get_tool(t) is None), None)
    if tool is None:
        raise RuntimeError("no free sys.monitoring tool id for the runtime leg")
    monitoring.use_tool_id(tool, "corroboration runtime leg")

    def started(code: types.CodeType, offset: int) -> None:
        name = codes.get(code)
        if name is not None:
            counts[name] = counts.get(name, 0) + 1

    monitoring.register_callback(tool, monitoring.events.PY_START, started)
    try:
        for code in codes:
            monitoring.set_local_events(tool, code, monitoring.events.PY_START)
        yield counts
    finally:
        for code in codes:
            monitoring.set_local_events(tool, code, 0)
        monitoring.register_callback(tool, monitoring.events.PY_START, None)
        monitoring.free_tool_id(tool)


def without(counts: Mapping[str, int], new: NewSide, old: Callable[..., T], *args: Any) -> T:
    """Evaluate the old spelling in this process (memo bypassed) and require that it ran no code of the new side."""
    before = dict(counts)
    with bypass():
        value = old(*args)
    reached = {k: counts[k] - before.get(k, 0) for k in counts if counts[k] != before.get(k, 0)}
    if reached:
        raise AssertionError(
            f"the old spelling called {new.label} ({reached}): this corroboration row compares {new.label} "
            f"with itself; {_DELETE}"
        )
    return value
