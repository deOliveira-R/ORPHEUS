"""What a test's own run READ, module-qualified (B5, rulings 47 to 49): the bindings it loaded, resolved to the statement
that binds each one, beside the defs it ran (the own cone, B3).

A test of a module's declared structure (a dataclass's fields, an enum's members, a Protocol by ``isinstance``, a
constant) runs none of the module's code: its subject is a statement that ran at collection, which the own recording
does not see (``collect_pins`` mints statement pins for files with a code object in the recording only). So each
executed first-party code object is read from its bytecode: every LOAD chain (``m.a.b``: a load, then its attribute
loads) is resolved where the code's own namespace binds its root (a local import, a cell bound by an import, the
module's globals) and then attribute by attribute, until it ends on a binding that is not a module. ``co_names`` is not
the read set: it holds globals and attributes in one unordered tuple, and cannot pair an attribute with its receiver.

A name is resolved to the statement that BINDS it. A re-export (``from x import y`` in a module the chain passed
through) is followed while the next module binds the same object (``is``), so a ``try/except ImportError`` fallback or a
``TYPE_CHECKING`` import stops at the module whose statement really binds it; each module passed through is recorded
(``reexport``), never a node. Module-level and class-body code ran at import, so it is not kept: it is compiled again
from the module's bytes and read with the same instruction reader, and an instruction is mapped to its statement by
its line (:func:`~orpheus.numerics._import_time.holding`, the innermost statement holding it, so a class's header
holds its bases and decorators and its member statements hold the body). This package parses no source (R7).

The channels of a :class:`Read` (:class:`Channel`):

* ``direct``: a chain of code the test ran.
* ``decorator``: a chain in the decorators or defaults of a statement on the test's qualified name (its classes'
  headers, then its def): evaluated at collection, the values reach the test as arguments.
* ``string``: an identifier token of a string constant of code the test ran, bound by a module that code read
  (``getattr(m, "X")``: pytest's assertion rewriter stores the constant apart from its receiver).
* ``module``: a chain that ends on a module object (``hasattr(m, "Old")``, ``m.__file__``), with no bound name; its
  pin names the tokens the code tried on it, so a later binding of a retired name stales it.
* ``via_apparatus``: a binding reached through the statement of a binding under ``tests/`` (a fixture, a case table, a
  test-module constant): the apparatus delivers a value, whose subject is the first binding outside the test tree.
* ``closure``: a binding read by the statement of a binding outside the test tree (what built ``LEGENDRE``): the
  value's provenance, pinned, never a node. Followed one statement deep.
* ``reexport``: a module a chain passed through by identity.
* ``ambiguous``: a file in ``sys.modules`` under two names whose relative imports resolve apart; the read is kept and
  counted.
* ``unresolved``: a name bound at run time that no statement binds (``setattr`` into a module, ``globals()[k] = v``),
  its dotted spelling as the name (C11); counted, never dropped. Its bare name joins the module's pin, so a statement
  that binds it later stales the pin. Met while following a binding's provenance (``closure``), it is not the test's
  read, and only the pin keeps it.

A fused instruction (``LOAD_FAST_LOAD_FAST``, ``STORE_FAST_LOAD_FAST``: Python 3.14 joins two adjacent local
accesses) is read as its two halves, the table derived from the opcode names (:data:`_FUSED`), so an instruction the
interpreter gains this way cannot fall through unread.

Not read, by declaration: pytest hook implementations (``pytest_*``: the floor every test runs, F1), module-level
code that ran during the test (an import's own reads are the imported module's provenance), a class attribute read
through an instance whose class the test never names, and a name built at run time.
"""
from __future__ import annotations

import dis
import importlib.util
import inspect
import os
import sys
import types
from collections import defaultdict
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path, PurePosixPath
from typing import Any, NamedTuple


class Channel(StrEnum):
    """How a read was reached (module docstring)."""

    DIRECT = "direct"
    DECORATOR = "decorator"
    STRING = "string"
    MODULE = "module"
    VIA_APPARATUS = "via_apparatus"
    CLOSURE = "closure"
    REEXPORT = "reexport"
    AMBIGUOUS = "ambiguous"
    UNRESOLVED = "unresolved"

    @property
    def is_node(self) -> bool:
        """Whether the read makes a node of the test (``closure`` and ``reexport`` are pinned only)."""
        return self in _NODE

    @property
    def is_pinned(self) -> bool:
        """Whether the read names a statement its module's pin is projected on (an unresolved or ambiguous read is
        counted, and names none)."""
        return self not in (Channel.UNRESOLVED, Channel.AMBIGUOUS)


_NODE = frozenset({Channel.DIRECT, Channel.DECORATOR, Channel.STRING, Channel.MODULE, Channel.VIA_APPARATUS})

#: Names the import system binds on a module object, which no statement binds: reading one reads the module.
MODULE_DUNDERS = frozenset({"__file__", "__name__", "__path__", "__spec__", "__loader__", "__package__", "__dict__",
                            "__cached__", "__builtins__", "__doc__", "__all__", "__annotations__"})
#: The spelling of a read of the module object (``Read.name is None``) in a row and among a node's read items.
MODULE_ITEM = ""

_LOAD_GLOBAL = frozenset({"LOAD_GLOBAL", "LOAD_NAME", "LOAD_FROM_DICT_OR_GLOBALS"})
_LOAD_LOCAL = frozenset({"LOAD_FAST", "LOAD_FAST_BORROW", "LOAD_FAST_CHECK", "LOAD_FAST_AND_CLEAR", "LOAD_DEREF",
                         "LOAD_FROM_DICT_OR_DEREF"})
_STORE = frozenset({"STORE_FAST", "STORE_DEREF", "STORE_NAME", "STORE_GLOBAL"})
_SKIP = frozenset({"CACHE", "NOP", "EXTENDED_ARG", "NOT_TAKEN"})
#: Each fused instruction and the two it joins: an opcode whose name is two opcode names joined (``argval`` is the
#: pair of their arguments).
_FUSED: dict[str, tuple[str, str]] = {
    name: ("_".join(parts[:i]), "_".join(parts[i:])) for name in dis.opmap for parts in [name.split("_")]
    for i in range(1, len(parts)) if "_".join(parts[:i]) in dis.opmap and "_".join(parts[i:]) in dis.opmap}
_MISSING = object()
#: ``inspect.CO_HAS_DOCSTRING`` (Python 3.14): the code's first constant is its docstring.
_CO_HAS_DOCSTRING: int = getattr(inspect, "CO_HAS_DOCSTRING")


def in_test_tree(path: str) -> bool:
    """Whether a repository-relative path is under ``tests/`` (apparatus, ruling 48)."""
    return PurePosixPath(path).parts[:1] == ("tests",)


def is_hook(qualname: str) -> bool:
    """Whether a def is a pytest hook implementation (``pytest_*``): the floor every test runs, never a read."""
    return qualname.rpartition(".")[2].startswith("pytest_")


class Read(NamedTuple):
    """One read of a test's own run: the repository-relative ``path`` of the module whose statement binds ``name``
    (``None``: the module object itself), reached by ``channel``."""

    path: str
    name: str | None
    channel: Channel

    def to_row(self) -> tuple[str, str, str]:
        return self.path, MODULE_ITEM if self.name is None else self.name, str(self.channel)

    @classmethod
    def from_row(cls, row: Any) -> Read:
        path, name, channel = row
        return cls(path, None if name == MODULE_ITEM else name, Channel(channel))


class Chain(NamedTuple):
    """A load chain of one code object: ``root`` loaded by ``scope`` (``global`` or ``local``), then ``attrs``, at
    ``line``."""

    scope: str
    root: str
    attrs: tuple[str, ...]
    line: int


class Imported(NamedTuple):
    """What an import statement bound a name to: the module ``module`` at relative ``level``, then the names taken
    from it (``from m import a`` takes ``a``; ``import a.b.c as d`` takes ``b``, ``c``), with ``star`` for
    ``from m import *``; ``fromlist`` tells the two forms apart."""

    module: str
    level: int
    taken: tuple[str, ...]
    fromlist: bool


class CodeReads(NamedTuple):
    """What one code object's instructions read: its load chains, the names its imports bind, its star imports."""

    chains: tuple[Chain, ...]
    bound: dict[str, Imported]
    stars: tuple[Imported, ...]


def _operations(code: types.CodeType) -> Iterator[tuple[str, Any, int]]:
    """``code``'s instructions as (opname, argument, line), each fused instruction as its two halves."""
    for ins in dis.get_instructions(code):
        if ins.opname in _SKIP:
            continue
        line = ins.positions.lineno if ins.positions is not None and ins.positions.lineno is not None else 0
        if (halves := _FUSED.get(ins.opname)) is not None:
            for opname, argument in zip(halves, ins.argval):
                yield opname, argument, line
        else:
            yield ins.opname, ins.argval, line


def instructions(code: types.CodeType) -> CodeReads:
    """The load chains of ``code`` (never of the code objects it holds), the names its imports bind, and its star
    imports. ``getattr(m, "c")``/``hasattr(m, "c")`` with a constant extends the chain on ``m`` by ``c``: the callee
    loads, then its receiver, then the constant (``probe`` is the step reached)."""
    chains: list[Chain] = []
    bound: dict[str, Imported] = {}
    stars: list[Imported] = []
    current: Chain | None = None
    probe: str | None = None  # ``callee`` after ``getattr``/``hasattr`` loads, ``receiver`` after the next load
    pending: Imported | None = None
    recent: list[Any] = []  # the constants just loaded (an import's level and fromlist)

    def close() -> None:
        nonlocal current
        if current is not None:
            chains.append(current)
        current = None

    for opname, argument, line in _operations(code):
        if opname in ("LOAD_CONST", "LOAD_SMALL_INT"):
            if probe == "receiver" and current is not None and isinstance(argument, str):
                current = current._replace(attrs=current.attrs + (argument,))
                close()
            probe = None
            recent = [*recent, argument][-2:]
            continue
        if opname == "IMPORT_NAME":
            close()
            level, fromlist = (recent[-2], recent[-1]) if len(recent) >= 2 else (0, None)
            pending = Imported(argument, level if isinstance(level, int) else 0, (), fromlist is not None)
            if fromlist == ("*",):
                stars.append(pending)
        elif opname == "IMPORT_FROM" and pending is not None:
            pending = pending._replace(taken=pending.taken + (argument,))
        elif opname in _STORE and pending is not None:
            bound[argument] = pending
            pending = pending._replace(taken=()) if pending.fromlist else None
        elif opname == "POP_TOP" and pending is not None and pending.fromlist:
            pending = None
        elif opname in _LOAD_GLOBAL or opname in _LOAD_LOCAL:
            close()
            probe = "receiver" if probe == "callee" else None
            if isinstance(argument, str):
                current = Chain("global" if opname in _LOAD_GLOBAL else "local", argument, (), line)
                if argument in ("getattr", "hasattr") and opname in _LOAD_GLOBAL:
                    probe = "callee"
        elif opname == "LOAD_ATTR" and current is not None:
            current = current._replace(attrs=current.attrs + (argument,))
        elif opname not in ("PUSH_NULL", "SWAP", "COPY"):
            close()
            probe = None
        recent = []
    close()
    return CodeReads(tuple(chains), bound, tuple(stars))


def _body_codes(code: types.CodeType) -> Iterator[types.CodeType]:
    """``code`` and every class body it holds, recursively: the code that runs at import (a def's body, a lambda's
    and an annotation scope's do not run there)."""
    yield code
    for constant in code.co_consts:
        if isinstance(constant, types.CodeType) and not constant.co_flags & inspect.CO_OPTIMIZED:
            yield from _body_codes(constant)


def _enclosing(qualname: str) -> list[str]:
    """The qualified names of the defs enclosing a nested def, outermost first (``f.<locals>.g`` is in ``f``)."""
    pieces = qualname.split(".<locals>.")
    return [".<locals>.".join(pieces[:i]) for i in range(1, len(pieces))]


@dataclass(frozen=True, eq=False)  # one per file, compared by identity
class _Source:
    """One first-party module as the resolution reads it: its object(s), its statement index, the chains of the code
    that ran when it was imported (compiled again from its bytes), each with the index of the statement holding its
    line, and the import bindings of its module body."""

    real: str
    path: str
    module: types.ModuleType
    names: frozenset[str]
    index: Any  # the memo's statement index (``traced_memo._source_index``)
    chains: tuple[tuple[int | None, Chain], ...]
    bound: dict[str, Imported]
    stars: tuple[Imported, ...]

    def header_chains(self, statement: int) -> list[Chain]:
        """The chains a statement evaluates itself: a def's decorators and defaults, a class's header, any other
        statement whole (a statement of its body holds its own)."""
        return [chain for holder, chain in self.chains if holder == statement]


class _Unbound(NamedTuple):
    """A name a module holds that no statement binds (set at run time): its module, bare and dotted."""

    source: _Source
    name: str
    dotted: str


class Sources:
    """The interpreter's first-party modules as they are after a test, each read once: the reusable half of a
    resolution (:func:`resolve`)."""

    def __init__(self, memo: types.ModuleType, first_party: Any, relative: Any) -> None:
        self.memo = memo
        self.first_party = first_party
        self.relative = relative
        self.by_file: dict[str, list[types.ModuleType]] = {}
        for module in list(sys.modules.values()):
            file = getattr(module, "__file__", None)
            if isinstance(file, str):
                self.by_file.setdefault(os.path.realpath(file), []).append(module)
        self.cache: dict[str, _Source | None] = {}

    def source(self, module: types.ModuleType) -> _Source | None:
        file = getattr(module, "__file__", None)
        return self.of_file(os.path.realpath(file)) if isinstance(file, str) else None

    def of_file(self, real: str) -> _Source | None:
        if real not in self.cache:
            self.cache[real] = self._read(real) if self.first_party(real) else None
        return self.cache[real]

    def _read(self, real: str) -> _Source | None:
        modules = self.by_file.get(real, [])
        if not modules:
            return None
        try:
            data = Path(real).read_bytes()
        except OSError:
            return None
        index = self.memo._source_index(data)
        statements = index.statements
        position = {id(statement): i for i, statement in enumerate(statements)}
        codes = [instructions(code) for code in _body_codes(compile(data, real, "exec", dont_inherit=True))]

        def holder(line: int) -> int | None:
            held = self.memo._import_time.holding(statements, line)
            return None if held is None else position[id(held)]

        chains = tuple((holder(chain.line), chain) for code in codes for chain in code.chains)
        return _Source(real, self.relative(real), modules[0], frozenset(m.__name__ for m in modules), index, chains,
                       codes[0].bound, codes[0].stars)

    def index_of(self, path: str) -> Any:
        return next((source.index for source in self.cache.values() if source is not None and source.path == path),
                    None)

    @staticmethod
    def absolute(imported: Imported, within: types.ModuleType) -> str:
        if not imported.level:
            return imported.module
        package = within.__name__ if hasattr(within, "__path__") else within.__name__.rpartition(".")[0]
        return importlib.util.resolve_name("." * imported.level + imported.module, package)

    def imported(self, imported: Imported, within: types.ModuleType) -> tuple[Any, str | None]:
        """What an import binding names: ``(module, None)`` for a module, ``(module, name)`` for a name taken from
        a module (``None`` for the first: nothing first-party)."""
        name = self.absolute(imported, within)
        if not imported.fromlist:
            dotted = name if imported.taken else name.split(".")[0]
            return sys.modules.get(dotted), None
        module = sys.modules.get(name)
        taken = imported.taken[-1] if imported.taken else None
        if module is None or taken is None:
            return None, None
        submodule = sys.modules.get(f"{name}.{taken}")
        if submodule is not None and module.__dict__.get(taken, _MISSING) is submodule:
            return submodule, None
        return module, taken


class RunReads(NamedTuple):
    """The reads of one test's own run, and per module path read the names its pin is projected on (the names read
    in it, the tokens tried on it as an object, the names its provenance met unbound)."""

    reads: frozenset[Read]
    names: Mapping[str, frozenset[str]]


def resolve(sources: Sources, codes: list[types.CodeType], test: tuple[str, str] | None) -> RunReads:
    """The reads of the executed ``codes``; ``test`` is the executed test's (real file, qualified name), whose
    statements' decorators and defaults are read as ``decorator``."""
    return _Run(sources).read(codes, test)


class _Run:
    """One resolution: the reads it has recorded so far, and the tokens tried on each module read as an object."""

    def __init__(self, sources: Sources) -> None:
        self.sources = sources
        self.reads: set[Read] = set()
        self.tried: dict[str, set[str]] = defaultdict(set)

    # ── a binding ──

    def binding(self, module: types.ModuleType, name: str, passed: list[tuple[_Source, str]]) -> Any:
        """Where ``module.name`` is bound: ``(source, name)`` for a statement's binding, a module object, an
        :class:`_Unbound`, or ``None`` (not first-party, or not bound). ``passed`` gains each (module, name) a
        re-export was followed through."""
        source = self.sources.source(module)
        if source is None:
            return None
        if name in MODULE_DUNDERS:
            return module
        value = module.__dict__.get(name, _MISSING)
        if isinstance(value, types.ModuleType):
            return value
        statements = [s for s in source.index.statements if s.owner is None and name in s.binds]
        # the import that may have bound it, then (when no statement names it) each star import
        candidates: list[tuple[Any, str | None]] = []
        if name in source.bound:
            candidates.append(self.sources.imported(source.bound[name], module))
        if not statements:
            candidates += [(sys.modules.get(self.sources.absolute(star, module)), name) for star in source.stars]
        for target, taken in candidates:
            if value is _MISSING or target is None or taken is None:
                continue
            if target.__dict__.get(taken, _MISSING) is value:  # a re-export, followed while the object is the same
                if self.sources.source(target) is None:
                    return None  # the object is another package's: not a first-party read
                passed.append((source, name))
                return self.binding(target, taken, passed)
        if any(not s.imports for s in statements):
            return source, name
        if value is not _MISSING and not statements:
            return _Unbound(source, name, f"{module.__name__}.{name}")
        return None

    # ── a chain ──

    def chain(self, module: types.ModuleType, chain: Chain, local: Mapping[str, Imported]) -> Any:
        """The end of ``chain`` read in ``module``'s code, whose local imports are ``local``; returns what
        :meth:`binding` returns, with the passed re-exports recorded."""
        passed: list[tuple[_Source, str]] = []
        if chain.scope == "local":
            if chain.root not in local:
                return None
            target, taken = self.sources.imported(local[chain.root], module)
            if target is None:
                return None
            current = target if taken is None else self.binding(target, taken, passed)
        else:
            current = self.binding(module, chain.root, passed)
        for attribute in chain.attrs:
            if not isinstance(current, types.ModuleType):
                break
            found = self.binding(current, attribute, passed)
            if found is None:
                if (source := self.sources.source(current)) is not None:
                    self.tried[source.path].add(attribute)
                break
            current = found
        # the reading code's own module binding its root by an import is how the code reached the name, not a re-export
        self.reads.update(Read(source.path, name, Channel.REEXPORT) for source, name in passed
                          if not (source.module is module and name == chain.root))
        return current

    def record(self, found: Any, channel: Channel) -> tuple[_Source, str] | None:
        """Record the end of a chain under ``channel``; the binding it names, when it names one."""
        match found:
            case types.ModuleType():
                if (source := self.sources.source(found)) is not None:
                    self.reads.add(Read(source.path, None, Channel.MODULE if channel is Channel.DIRECT else channel))
            case _Unbound(source, name, dotted):
                self.tried[source.path].add(name)  # a statement binding it later stales the module's pin
                if channel.is_node:
                    self.reads.add(Read(source.path, dotted, Channel.UNRESOLVED))
            case (source, name):
                self.reads.add(Read(source.path, name, channel))
                if len(source.names) > 1:
                    self.reads.add(Read(source.path, None, Channel.AMBIGUOUS))
                return source, name
        return None

    def follow(self, start: tuple[_Source, str], apparatus: bool) -> None:
        """The statement of a resolved binding: through apparatus (a binding under ``tests/``) to the first binding
        outside the test tree (``via_apparatus``), and one statement deep from a binding outside it (``closure``)."""
        work, seen = [(start, apparatus)], {start}
        while work:
            (source, name), through = work.pop()
            binders = [i for i, s in enumerate(source.index.statements)
                       if s.owner is None and not s.imports and name in s.binds]
            for chain in (chain for i in binders for chain in source.header_chains(i)):
                channel = Channel.VIA_APPARATUS if through else Channel.CLOSURE
                found = self.record(self.chain(source.module, chain, source.bound), channel)
                if found is None or found in seen or not through:
                    continue
                seen.add(found)
                work.append((found, in_test_tree(found[0].path)))

    # ── the run ──

    def read(self, codes: list[types.CodeType], test: tuple[str, str] | None) -> RunReads:
        cells: dict[tuple[str, str], dict[str, Imported]] = {}  # (file, binding def) to the cells its imports bind
        own: list[tuple[types.ModuleType, types.CodeType, CodeReads]] = []
        for code in codes:
            real = os.path.realpath(code.co_filename)
            if not self.sources.first_party(real) or is_hook(code.co_qualname) or code.co_name == "<module>":
                continue
            if not code.co_flags & inspect.CO_OPTIMIZED and "<locals>" not in code.co_qualname:
                continue  # a class body that ran at import
            modules = self.sources.by_file.get(real)
            if not modules:
                continue
            found = instructions(code)
            cells[(real, code.co_qualname)] = {name: imported for name, imported in found.bound.items()
                                               if name in code.co_cellvars}
            own.append((modules[0], code, found))
        resolved: list[tuple[_Source, str]] = []
        for module, code, found in own:
            real = os.path.realpath(code.co_filename)
            local = {name: imported for outer in _enclosing(code.co_qualname)
                     for name, imported in cells.get((real, outer), {}).items()} | found.bound
            reached: list[types.ModuleType] = []
            for chain in found.chains:
                end = self.chain(module, chain, local)
                if isinstance(end, types.ModuleType):
                    reached.append(end)
                if (binding := self.record(end, Channel.DIRECT)) is not None:
                    resolved.append(binding)
            constants = code.co_consts[1:] if code.co_flags & _CO_HAS_DOCSTRING else code.co_consts
            tokens = frozenset(self.sources.memo._import_time.string_names(constants))  # the docstring excluded
            for target in reached:
                if (source := self.sources.source(target)) is None:
                    continue
                self.tried[source.path] |= tokens
                for token in sorted(tokens):
                    if token in MODULE_DUNDERS or token not in target.__dict__:
                        continue
                    if (binding := self.record(self.binding(target, token, []), Channel.STRING)) is not None:
                        resolved.append(binding)
        if test is not None:
            resolved += self.decorators(*test)
        for binding in resolved:
            self.follow(binding, in_test_tree(binding[0].path))
        names: dict[str, set[str]] = defaultdict(set)
        for read in self.reads:
            if read.channel.is_pinned:
                names[read.path] |= set() if read.name is None else {read.name}
        for path, tried in self.tried.items():
            names[path] |= tried
        return RunReads(frozenset(self.reads), {path: frozenset(found) for path, found in names.items()})

    def decorators(self, real: str, qualname: str) -> list[tuple[_Source, str]]:
        """The bindings read by the decorators and defaults of each statement on ``qualname`` in ``real`` (its
        classes' headers, then its def)."""
        source = self.sources.of_file(real)
        if source is None:
            return []
        found: list[tuple[_Source, str]] = []
        owner: int | None = None
        for part in (p for p in qualname.split(".") if p != "<locals>"):
            hits = [i for i, s in enumerate(source.index.statements)
                    if s.owner == owner and not s.imports and part in s.binds]
            if not hits:
                break
            owner = hits[-1]  # the last binding is the one Python keeps
            for chain in source.header_chains(owner):
                if (end := self.record(self.chain(source.module, chain, source.bound), Channel.DECORATOR)) is not None:
                    found.append(end)
        return found
