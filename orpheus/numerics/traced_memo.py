r"""The traced memo: a pure function of content, memoised on disk and keyed on the code that ran (#405 P3).

A reference reading costs seconds to minutes and changes only when the code that produces it changes. The
memo stores each call's answer under ``.cache/references/`` and serves it again for as long as nothing that
ran to produce it has changed. What ran is not declared, it is RECORDED: on a miss the call is generated
in a fresh interpreter traced by :mod:`sys.monitoring` from its first line, and the entry keeps the
:class:`Manifest` of that run beside the answer.

* **The key** is the digest of the function's identity, its signature-bound arguments (defaults applied,
  each declared canonical form put in), the type tree of those arguments, and the platform tag. The
  argument digest is :func:`~orpheus.numerics.content.encode`'s, which follows ``==`` (``8 == 8.0``); the
  type tree separates what that identifies and the function can tell apart.
* **The manifest** pins what ran: each first-party def by the digest of its normalised source (docstrings,
  comments and layout excluded), each first-party module that ran by its skeleton (the module with every
  function body removed: imports, constants, decorators, signatures, class attributes), each third-party
  distribution by its version, the interpreter, each data file the run opened by its bytes, and each memo
  the run called by its key and payload digest (validated recursively).
* **Validation** re-hashes the manifest against the files on the current ``sys.path``; it imports and runs
  nothing. A changed pin makes the entry :class:`Stale`, an unreadable or inconsistent entry
  :class:`Corrupt`; neither is ever served, and the next call regenerates it.
* **The payload** is JSON (every float as ``float.hex``, so it crosses bit for bit) and one ``.npz`` for
  the arrays, never a pickle; each load hands out fresh, read-only arrays.

Why a fresh process (the user's ruling of 2026-10-04, ``.claude/plans/reference_cache.md``, "P3 rulings"):
a trace records only code that RUNS, so anything served from memory (a ``functools.cache`` hit, a
``cached_property`` computed earlier, a monkeypatched function) would hide a dependency from it. A process
born for one call has nothing in memory, so the trace is complete without any discipline on how the rest
of the package memoises. The cost is one interpreter start per miss, about 0.65 s ``[M]`` 2026-10-04.

What no manifest can see: native state (a C library's own globals) and a file a C library opens itself
(HDF5 through ``h5py``), which raises no audit event; a generator reading such a file takes the file's
digest as an argument. A test that monkeypatches anything a generation runs reads under :func:`bypass`.
"""
from __future__ import annotations

import ast
import contextlib
import contextvars
import dataclasses
import functools
import hashlib
import importlib
import importlib.metadata
import inspect
import io
import json
import os
import pickle
import platform
import shutil
import site
import subprocess
import sys
import sysconfig
import types
import typing
import uuid
from collections.abc import Callable, Generator, Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any, NamedTuple

import numpy as np

from orpheus.numerics.content import ContentIdentity, _rebuild, encode

__all__ = [
    "Absent", "ChildPin", "Corrupt", "DataPin", "DefPin", "DistributionPin", "Hit", "Manifest", "ModulePin",
    "Stale", "TracedMemo", "Unencodable", "Unpinnable", "Verdict", "bypass", "cache_root", "decode_payload",
    "default_root", "encode_payload", "function_digest", "platform_tag", "skeleton_digest", "trace_call",
    "traced_memo", "validate",
]

#: The entry format; an entry of another schema is stale.
_SCHEMA = 1


def _digest(data: bytes) -> str:
    return hashlib.blake2b(data, digest_size=32).hexdigest()


class Unpinnable(RuntimeError):
    """A traced code object no manifest can pin: a site-packages file of no distribution, or code compiled
    under a file's name at a line where the file holds no such def (validation could never re-hash it)."""


class Unencodable(TypeError):
    """A return value the payload cannot write exactly; it is refused, never written lossily."""


# ── the normalised source ───────────────────────────────────────────────────────────────────────────

_DEF = (ast.FunctionDef, ast.AsyncFunctionDef)
#: Statements whose bodies hold module-level defs (a def under ``if``/``try``/``with`` is still top level).
_BLOCK = (ast.If, ast.Try, ast.With)


def _without_docstrings(tree: ast.Module) -> ast.Module:
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, *_DEF)):
            first = node.body[0] if node.body else None
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) and isinstance(first.value.value, str):
                node.body = node.body[1:] or [ast.Pass()]
    return tree


def _ast_digest(node: ast.AST) -> str:
    """The digest of ``node``'s structure: line and column attributes excluded, so layout never moves it."""
    return _digest(ast.dump(node, include_attributes=False).encode())


def _block_bodies(node: ast.stmt) -> list[ast.stmt]:
    return [*getattr(node, "body", []), *getattr(node, "orelse", []), *getattr(node, "finalbody", []),
            *(s for h in getattr(node, "handlers", []) for s in h.body)]


@dataclasses.dataclass(frozen=True)
class _SourceText:
    """One source file, parsed once: its tree, the digest of each outermost def, and its skeleton's digest."""

    tree: ast.Module
    defs: Mapping[str, str]
    skeleton: str


@functools.cache
def _source_text(data: bytes) -> _SourceText:
    """The parse of a source text, cached by its BYTES (an edited file is another key, never a stale read)."""
    tree = ast.parse(data)
    stripped = _without_docstrings(ast.parse(data))
    defs: dict[str, str] = {}

    def collect(body: list[ast.stmt], prefix: str) -> None:
        for node in body:
            if isinstance(node, _DEF):
                defs[prefix + node.name] = _ast_digest(node)  # a later def of one name wins, as Python binds it
            elif isinstance(node, ast.ClassDef):
                collect(node.body, f"{prefix}{node.name}.")
            elif isinstance(node, _BLOCK):
                collect(_block_bodies(node), prefix)

    collect(stripped.body, "")
    for node in ast.walk(stripped):
        if isinstance(node, _DEF):
            node.body = [ast.Pass()]
    return _SourceText(tree, types.MappingProxyType(defs), _ast_digest(stripped))


def function_digest(source: bytes | str, qualname: str) -> str | None:
    """The digest of the outermost def ``qualname`` (``function`` or ``Class.method``) with its docstrings
    stripped; ``None`` when the source holds no such def."""
    return _source_text(source.encode() if isinstance(source, str) else source).defs.get(qualname)


def skeleton_digest(source: bytes | str) -> str:
    """The digest of the module with every function body replaced by ``pass`` and every docstring stripped:
    what a traced function reads by name without its own body changing (constants, imports, decorators,
    signatures with their defaults, class attributes, dataclass fields)."""
    return _source_text(source.encode() if isinstance(source, str) else source).skeleton


#: Code objects Python generates from a SIGNATURE or a class body rather than from a def body: lazy
#: annotations (PEP 649, ``__annotate__``) and the scopes of type parameters (PEP 695). Their source is the
#: signature, which the skeleton pins; they never pin the body of the def they belong to (reading a def's
#: annotations does not run it).
_SIGNATURE_SCOPES = ("__annotate__", "<generic parameters of ", "<TypeVar bound of ", "<TypeVar constraint of ",
                     "<TypeAlias ", "<type parameters of ")
#: Anonymous code objects and the expression that compiles to each.
_ANONYMOUS: Mapping[str, type[ast.AST]] = {
    "<lambda>": ast.Lambda, "<genexpr>": ast.GeneratorExp, "<listcomp>": ast.ListComp,
    "<setcomp>": ast.SetComp, "<dictcomp>": ast.DictComp,
}


def _span(node: ast.AST) -> range:
    first = min([node.lineno, *(d.lineno for d in getattr(node, "decorator_list", []))])  # type: ignore[attr-defined]
    return range(first, node.end_lineno + 1)  # type: ignore[attr-defined]


def _holds(scope: ast.AST, line: int, name: str) -> bool:
    """Whether ``scope`` holds the source of a code object named ``name`` first reported at ``line``."""
    if name.startswith(_SIGNATURE_SCOPES):
        return True
    for node in ast.walk(scope):
        if isinstance(node, (*_DEF, ast.ClassDef)) and node.name == name and line in _span(node):
            return True
        if name in _ANONYMOUS and isinstance(node, _ANONYMOUS[name]) and line in _span(node):
            return True
    return False


class _ModuleLevel:
    """The code object belongs to no def: a class body, a module-level lambda, a signature; the skeleton pins it."""


_MODULE_LEVEL = _ModuleLevel()


def _outermost_def(tree: ast.Module, line: int, name: str) -> str | _ModuleLevel | None:
    """The qualified name of the outermost def whose source holds the code object ``name`` first reported at
    ``line``; :data:`_MODULE_LEVEL` for module-level code; ``None`` when the file holds no such code there."""

    def search(body: list[ast.stmt], prefix: str) -> str | _ModuleLevel | None:
        for node in body:
            if isinstance(node, _DEF) and line in _span(node):
                if name.startswith(_SIGNATURE_SCOPES) and line < node.body[0].lineno:
                    return _MODULE_LEVEL  # the def's own signature
                return prefix + node.name if _holds(node, line, name) else None
            if isinstance(node, ast.ClassDef) and line in _span(node):
                inner = search(node.body, f"{prefix}{node.name}.")
                if isinstance(inner, str):
                    return inner
                return _MODULE_LEVEL if _holds(node, line, name) else None
            if isinstance(node, _BLOCK):
                inner = search(_block_bodies(node), prefix)
                if inner is not _MODULE_LEVEL:
                    return inner
        return _MODULE_LEVEL

    found = search(tree.body, "")
    if found is _MODULE_LEVEL and not _holds(tree, line, name):
        return None
    return found


# ── where a traced file belongs ─────────────────────────────────────────────────────────────────────


def _real(path: str) -> str:
    return os.path.realpath(path)


@functools.cache
def _site_dirs() -> tuple[str, ...]:
    paths = sysconfig.get_paths()
    dirs = {*site.getsitepackages(), paths["purelib"], paths["platlib"]}
    return tuple(sorted({_real(d) + os.sep for d in dirs}, key=len, reverse=True))


@functools.cache
def _stdlib_dirs() -> tuple[str, ...]:
    paths = sysconfig.get_paths()
    return tuple({_real(paths["stdlib"]) + os.sep, _real(paths["platstdlib"]) + os.sep})


@functools.cache
def _distributions_by_top_level() -> Mapping[str, tuple[str, ...]]:
    return {top: tuple(sorted(names)) for top, names in importlib.metadata.packages_distributions().items()}


def _recording_distribution(site_dir: str, relative: str) -> tuple[str, bool] | None:
    """The distribution whose ``RECORD`` lists ``relative`` (a path under ``site_dir``) and whether it is an
    editable install. Searching the raw ``RECORD`` texts costs 0.23 s; building the whole file index through
    :mod:`importlib.metadata` costs 9.4 to 10.9 s per process ``[M]`` 2026-10-04 (24 484 files)."""
    for info in Path(site_dir).glob("*.dist-info"):
        record = info / "RECORD"
        if record.is_file() and f"\n{relative}," in "\n" + record.read_text():
            distribution = importlib.metadata.Distribution.at(info)
            direct = json.loads(distribution.read_text("direct_url.json") or "{}")
            return distribution.metadata["Name"], bool(direct.get("dir_info", {}).get("editable"))
    return None


class Dropped(NamedTuple):
    """Generated code (``<string>``, ``<frozen …>``): pinned by the interpreter version and the skeleton."""

    filename: str


class Interpreter(NamedTuple):
    """A standard-library file: pinned by the interpreter version."""

    path: str


class Distribution(NamedTuple):
    """A file of an installed distribution: pinned by the distribution's version."""

    name: str


class Source(NamedTuple):
    """A file pinned by its content: first-party source, or an editable install's own file."""

    path: str


Origin = Dropped | Interpreter | Distribution | Source


def origin(filename: str) -> Origin:
    """Where a traced file belongs, and so how the manifest pins it."""
    if filename.startswith("<"):
        return Dropped(filename)
    real = _real(filename)
    for site_dir in _site_dirs():
        if real.startswith(site_dir):
            relative = real[len(site_dir):]
            top = relative.split(os.sep)[0].removesuffix(".py")
            if names := _distributions_by_top_level().get(top):
                return Distribution(names[0])
            recorded = _recording_distribution(site_dir, Path(relative).as_posix())
            if recorded is None:
                raise Unpinnable(f"{filename}: in site-packages and in no installed distribution, so no version pins it")
            name, editable = recorded
            return Source(real) if editable else Distribution(name)  # an editable version never moves
    if any(real.startswith(d) for d in _stdlib_dirs()):
        return Interpreter(real)
    return Source(real)


def _relative(real: str) -> str:
    """``real`` relative to the ``sys.path`` entry that holds it: how validation finds it again."""
    for entry in sys.path:
        root = _real(entry or os.getcwd()) + os.sep
        if real.startswith(root):
            return real[len(root):]
    raise Unpinnable(f"{real}: a traced source file on no sys.path entry, so validation could not find it")


def _locate(relative: str) -> Path | None:
    for entry in sys.path:
        candidate = Path(entry or os.getcwd()) / relative
        if candidate.is_file():
            return candidate
    return None


def python_identity() -> str:
    return f"{sys.version}|{sys.implementation.cache_tag}"


@functools.cache
def platform_tag() -> str:
    """The operating system, the machine, the ``-O`` level (a bare ``assert`` in a generator is stripped under
    it) and the BLAS (Accelerate and OpenBLAS differ in the last bits, #504)."""
    config = np.show_config(mode="dicts")
    blas = config.get("Build Dependencies", {}).get("blas", {}).get("name", "unknown") if isinstance(config, dict) else "unknown"
    return f"{sys.implementation.cache_tag}-{sys.platform}-{platform.machine()}-O{sys.flags.optimize}-{blas}"


# ── the manifest ────────────────────────────────────────────────────────────────────────────────────


class _Checkout:
    """The files on the current ``sys.path`` and the store, as one validation reads them (no read outlives it)."""

    def __init__(self, store: Store) -> None:
        self.store = store
        self._texts: dict[str, bytes | None] = {}

    def text(self, relative: str) -> bytes | None:
        if relative not in self._texts:
            path = _locate(relative)
            self._texts[relative] = path.read_bytes() if path is not None else None
        return self._texts[relative]


class DefPin(NamedTuple):
    """A first-party def that ran, by its normalised source."""

    path: str
    qualname: str
    digest: str

    def stale(self, checkout: _Checkout) -> str | None:
        text = checkout.text(self.path)
        if text is None:
            return None  # the module's own pin reports the file gone
        found = function_digest(text, self.qualname)
        if found is None:
            return f"function {self.path}:{self.qualname} gone"
        return f"function {self.path}:{self.qualname} changed" if found != self.digest else None


class ModulePin(NamedTuple):
    """A first-party module whose body ran or which holds a def that ran, by its skeleton."""

    path: str
    digest: str

    def stale(self, checkout: _Checkout) -> str | None:
        text = checkout.text(self.path)
        if text is None:
            return f"module {self.path} gone"
        return f"module {self.path} skeleton changed" if skeleton_digest(text) != self.digest else None


class DistributionPin(NamedTuple):
    """A third-party distribution whose code ran, by its version."""

    name: str
    version: str

    def stale(self, checkout: _Checkout) -> str | None:
        try:
            now: str | None = importlib.metadata.version(self.name)
        except importlib.metadata.PackageNotFoundError:
            now = None
        return f"distribution {self.name}: {self.version} != {now}" if now != self.version else None


class DataPin(NamedTuple):
    """A file the run opened for reading that no other pin covers, by its bytes."""

    path: str
    digest: str

    def stale(self, checkout: _Checkout) -> str | None:
        path = Path(self.path) if os.path.isabs(self.path) else _locate(self.path)
        if path is None or not path.is_file():
            return f"data {self.path} gone"
        return f"data {self.path} changed" if _digest(path.read_bytes()) != self.digest else None


class ChildPin(NamedTuple):
    """A memo the run called, hit or generated, by its key and the digest of the payload it served."""

    function_id: str
    key: str
    payload_digest: str

    def stale(self, checkout: _Checkout) -> str | None:
        verdict = checkout.store.lookup(self.function_id, self.key)
        if not isinstance(verdict, Hit):
            return f"child {self.function_id}/{self.key[:12]}: {verdict}"
        if verdict.payload_digest != self.payload_digest:
            return f"child {self.function_id}/{self.key[:12]}: its payload {self.payload_digest[:12]} became {verdict.payload_digest[:12]}"
        return None


_PINS: Mapping[str, Callable[..., Any]] = {
    "functions": DefPin, "modules": ModulePin, "distributions": DistributionPin, "children": ChildPin, "data": DataPin,
}


@dataclasses.dataclass(frozen=True)
class Manifest:
    """What one generation ran, pinned so that a later checkout can be compared with it without running it."""

    functions: tuple[DefPin, ...]
    modules: tuple[ModulePin, ...]
    distributions: tuple[DistributionPin, ...]
    python: str
    children: tuple[ChildPin, ...]
    data: tuple[DataPin, ...] = ()

    def __post_init__(self) -> None:
        for field, pin in _PINS.items():  # each row becomes its pin, so every pin can validate itself
            object.__setattr__(self, field, tuple(pin(*row) for row in getattr(self, field)))

    def pins(self) -> Iterator[DefPin | ModulePin | DistributionPin | ChildPin | DataPin]:
        yield from (*self.distributions, *self.modules, *self.functions, *self.data, *self.children)

    def to_json(self) -> dict[str, Any]:
        return {"python": self.python} | {field: [list(pin) for pin in getattr(self, field)] for field in _PINS}

    @classmethod
    def from_json(cls, data: Mapping[str, Any]) -> Manifest:
        return cls(python=str(data["python"]), **{field: data[field] for field in _PINS})


def _data_pins(opened: Sequence[str], store: Store) -> tuple[DataPin, ...]:
    """Every file the run opened for reading that no other pin covers: not code (``.py``, ``.pyc``), not the
    standard library or a distribution (versioned), not an entry of the store (a child pins it)."""
    excluded = (_real(str(store.root)) + os.sep, *_site_dirs(), *_stdlib_dirs())
    pins = set()
    for name in opened:
        real = _real(name)
        if real.endswith((".py", ".pyc")) or real.startswith(excluded) or not os.path.isfile(real):
            continue
        try:
            label = _relative(real)
        except Unpinnable:
            label = real
        pins.add(DataPin(label, _digest(Path(real).read_bytes())))
    return tuple(sorted(pins))


def _manifest(codes: set[types.CodeType], children: Sequence[ChildPin], opened: Sequence[str], store: Store) -> Manifest:
    """The manifest of a run that executed ``codes``, called ``children`` and opened ``opened``."""
    by_file: dict[str, list[types.CodeType]] = {}
    distributions: set[str] = set()
    for code in codes:
        match origin(code.co_filename):
            case Distribution(name):
                distributions.add(name)
            case Source(path):
                by_file.setdefault(path, []).append(code)
            case Dropped() | Interpreter():
                pass
    functions: set[DefPin] = set()
    modules: set[ModulePin] = set()
    for real, file_codes in by_file.items():
        relative = _relative(real)
        text = _source_text(Path(real).read_bytes())
        modules.add(ModulePin(relative, text.skeleton))
        for code in file_codes:
            if code.co_name == "<module>":
                continue
            match _outermost_def(text.tree, code.co_firstlineno, code.co_name):
                case str(qualname):
                    functions.add(DefPin(relative, qualname, text.defs[qualname]))
                case None:
                    raise Unpinnable(f"{real}:{code.co_firstlineno} {code.co_qualname}: no def of that name at its line")
                case _:
                    pass  # module-level code: the skeleton pins it
    return Manifest(
        tuple(sorted(functions)), tuple(sorted(modules)),
        tuple(sorted(DistributionPin(d, importlib.metadata.version(d)) for d in distributions)),
        python_identity(), tuple(sorted(set(children))), _data_pins(opened, store),
    )


def validate(manifest: Manifest, store: Store | None = None) -> tuple[str, ...]:
    """Every reason ``manifest`` no longer describes the checkout (empty: it does). Imports and runs nothing."""
    checkout = _Checkout(store if store is not None else Store.current())
    reasons = [f"python: {manifest.python!r} != {python_identity()!r}"] if manifest.python != python_identity() else []
    return (*reasons, *(reason for pin in manifest.pins() if (reason := pin.stale(checkout)) is not None))


# ── tracing ─────────────────────────────────────────────────────────────────────────────────────────

#: The memos a generation in progress has called, so that its manifest pins them as children.
_CHILDREN: contextvars.ContextVar[list[ChildPin] | None] = contextvars.ContextVar("traced_memo_children", default=None)


@contextlib.contextmanager
def _tracing() -> Generator[set[types.CodeType]]:
    """Every code object that starts while the block runs (each reported once: the callback disables itself)."""
    monitoring = sys.monitoring
    tool = next((t for t in (monitoring.PROFILER_ID, 3, 4, monitoring.OPTIMIZER_ID) if monitoring.get_tool(t) is None), None)
    if tool is None:
        raise RuntimeError("traced_memo: no free sys.monitoring tool id")
    monitoring.use_tool_id(tool, "traced_memo")
    codes: set[types.CodeType] = set()

    def started(code: types.CodeType, _offset: int) -> object:
        codes.add(code)
        return monitoring.DISABLE

    monitoring.register_callback(tool, monitoring.events.PY_START, started)
    monitoring.set_events(tool, monitoring.events.PY_START)
    try:
        yield codes
    finally:
        monitoring.set_events(tool, 0)
        monitoring.register_callback(tool, monitoring.events.PY_START, None)
        monitoring.free_tool_id(tool)
        monitoring.restart_events()


def trace_call(function: Callable[..., Any], *args: Any, **kwargs: Any) -> tuple[Any, Manifest]:
    """Call ``function`` in THIS process under the tracer: its value and the manifest of what ran. Code that
    ran earlier in the process and is served from memory now is not in it, which is why a generation runs in
    a fresh process; this is the tracer's own instrument."""
    children: list[ChildPin] = []
    token = _CHILDREN.set(children)
    try:
        with _tracing() as codes:
            value = function(*args, **kwargs)
        return value, _manifest(codes, children, (), Store.current())
    finally:
        _CHILDREN.reset(token)


# ── the payload ─────────────────────────────────────────────────────────────────────────────────────


def encode_payload(value: Any) -> tuple[Any, dict[str, np.ndarray]]:
    """``(JSON tree, arrays)``: floats as ``float.hex``, numpy scalars with their dtype, arrays by name, tuples,
    and frozen dataclasses by their constructor fields. Anything else is :class:`Unencodable`."""
    arrays: dict[str, np.ndarray] = {}

    def node(v: Any) -> Any:
        if v is None or isinstance(v, (bool, str)):
            return {"v": v}
        if isinstance(v, np.generic):  # before ``float``: ``np.float64`` is a ``float``
            if v.dtype.kind not in "biuf":
                raise Unencodable(f"a {v.dtype} scalar has no exact payload")
            return {"np": v.dtype.str, "a": node(np.asarray(v))["a"]}
        if isinstance(v, int):
            return {"i": str(v)}
        if isinstance(v, float):
            return {"f": float.hex(v)}
        if isinstance(v, np.ndarray):
            if v.dtype.kind not in "biufc":
                raise Unencodable(f"a {v.dtype} array has no exact payload")
            name = f"a{len(arrays)}"
            arrays[name] = v
            return {"a": name}
        if isinstance(v, tuple):
            return {"t": [node(x) for x in v]}
        if dataclasses.is_dataclass(v) and not isinstance(v, type):
            cls = type(v)
            return {"dc": f"{cls.__module__}:{cls.__qualname__}",
                    "fields": {f.name: node(getattr(v, f.name)) for f in dataclasses.fields(v) if f.init}}
        raise Unencodable(f"a {type(v).__module__}.{type(v).__qualname__} has no exact payload")

    return node(value), arrays


def _returned_types(function: Callable[..., Any]) -> set[str]:
    """The dataclasses ``function``'s return annotation names, transitively through their fields: the only
    types a payload of it may name (a payload never imports a type its function does not declare)."""
    found: set[str] = set()

    def visit(annotation: Any) -> None:
        for member in typing.get_args(annotation) or (annotation,):
            if isinstance(member, type) and dataclasses.is_dataclass(member):
                name = f"{member.__module__}:{member.__qualname__}"
                if name not in found:
                    found.add(name)
                    for hint in typing.get_type_hints(member).values():
                        visit(hint)

    visit(inspect.get_annotations(function, eval_str=True).get("return"))
    return found


def _resolve(dotted: str) -> Any:
    module, _, qualname = dotted.partition(":")
    value: Any = importlib.import_module(module)
    for part in qualname.split("."):
        value = getattr(value, part)
    return value


def decode_payload(tree: Any, arrays: Mapping[str, np.ndarray], allowed: set[str]) -> Any:
    """The value :func:`encode_payload` wrote: every array a fresh read-only copy, every dataclass rebuilt
    through its constructor (so its laws re-run), and only a dataclass named in ``allowed``."""

    def value(v: Mapping[str, Any]) -> Any:
        match v:
            case {"v": plain}:
                return plain
            case {"i": text}:
                return int(text)
            case {"f": text}:
                return float.fromhex(text)
            case {"np": dtype, "a": name}:
                return np.dtype(dtype).type(arrays[name][()])
            case {"a": name}:
                array = np.array(arrays[name])
                array.setflags(write=False)
                return array
            case {"t": items}:
                return tuple(value(x) for x in items)
            case {"dc": dotted, "fields": fields}:
                if dotted not in allowed:
                    raise ValueError(f"the payload names {dotted}, which the function does not return")
                return _resolve(dotted)(**{k: value(x) for k, x in fields.items()})
            case _:
                raise ValueError(f"an unknown payload node {sorted(v)}")

    return value(tree)


def _payload_digest(tree: Any, npz: bytes) -> str:
    return _digest(json.dumps(tree, sort_keys=True).encode() + b"\0" + npz)


# ── the store ───────────────────────────────────────────────────────────────────────────────────────


@dataclasses.dataclass(frozen=True)
class Hit:
    """A valid entry: its directory, its payload's digest and its payload tree."""

    directory: Path
    payload_digest: str
    tree: Any


@dataclasses.dataclass(frozen=True)
class Absent:
    """No entry."""


@dataclasses.dataclass(frozen=True)
class Stale:
    """An entry whose manifest no longer describes the checkout, with every reason."""

    reasons: tuple[str, ...]


@dataclasses.dataclass(frozen=True)
class Corrupt:
    """An entry that cannot be read, or whose payload does not match its digest."""

    reason: str


Verdict = Hit | Absent | Stale | Corrupt


def default_root() -> Path:
    """``.cache/references`` at the repository root (gitignored)."""
    return Path(__file__).resolve().parents[2] / ".cache" / "references"


_ROOT: contextvars.ContextVar[Path | None] = contextvars.ContextVar("traced_memo_root", default=None)
_BYPASS: contextvars.ContextVar[bool] = contextvars.ContextVar("traced_memo_bypass", default=False)


@contextlib.contextmanager
def cache_root(path: Path) -> Generator[Path]:
    """Read and write the store under ``path`` while the block runs."""
    token = _ROOT.set(Path(path))
    try:
        yield Path(path)
    finally:
        _ROOT.reset(token)


@contextlib.contextmanager
def bypass() -> Generator[None]:
    """Every memoised call runs in THIS process and reads and writes nothing: the spelling for a test that
    monkeypatches anything a generation runs (a patch never reaches a fresh process, so without the bypass the
    test would read the honest entry instead of its patched answer)."""
    token = _BYPASS.set(True)
    try:
        yield
    finally:
        _BYPASS.reset(token)


@dataclasses.dataclass(frozen=True)
class Store:
    """The entries under one root: ``<root>/<function id>/<key>/{entry.json, payload.npz}``."""

    root: Path

    @classmethod
    def current(cls) -> Store:
        return cls(_ROOT.get() or default_root())

    def directory(self, function_id: str, key: str) -> Path:
        return self.root / function_id / key

    def lookup(self, function_id: str, key: str) -> Verdict:
        directory = self.directory(function_id, key)
        if not (directory / "entry.json").exists():
            return Absent()
        try:
            entry = json.loads((directory / "entry.json").read_text())
            npz = (directory / "payload.npz").read_bytes()
            if entry["schema"] != _SCHEMA:
                return Stale((f"schema {entry['schema']} != {_SCHEMA}",))
            if _payload_digest(entry["payload"], npz) != entry["payload_digest"]:
                return Corrupt("the payload does not match its digest")
            manifest = Manifest.from_json(entry["manifest"])
        except (OSError, ValueError, KeyError, TypeError) as error:
            return Corrupt(f"{type(error).__name__}: {error}")
        if reasons := validate(manifest, self):
            return Stale(reasons)
        return Hit(directory, entry["payload_digest"], entry["payload"])

    def load(self, hit: Hit, allowed: set[str]) -> Any:
        with np.load(hit.directory / "payload.npz", allow_pickle=False) as npz:
            arrays = {name: npz[name] for name in npz.files}
        return decode_payload(hit.tree, arrays, allowed)

    def write(self, function_id: str, key: str, manifest: Manifest, value: Any) -> None:
        """Stage the entry in a sibling directory, then move it into place in one ``os.replace``: a reader sees
        the whole entry or none of it."""
        tree, arrays = encode_payload(value)
        buffer = io.BytesIO()
        np.savez(buffer, allow_pickle=False, **arrays)
        npz = buffer.getvalue()
        directory = self.directory(function_id, key)
        staging = directory.parent / f".staging-{uuid.uuid4().hex}"
        staging.mkdir(parents=True)
        (staging / "payload.npz").write_bytes(npz)
        (staging / "entry.json").write_text(json.dumps({
            "schema": _SCHEMA, "function": function_id, "key": key, "manifest": manifest.to_json(),
            "payload": tree, "payload_digest": _payload_digest(tree, npz),
        }))
        if directory.exists():
            shutil.rmtree(directory, ignore_errors=True)
        os.replace(staging, directory)


# ── the key and the process boundary ────────────────────────────────────────────────────────────────


def _type_tree(value: Any) -> Any:
    """The types inside a top-level argument that its content digest identifies (an array's dtype, a
    container's element types); a :class:`ContentIdentity` value's constructor canonicalises its own."""
    if isinstance(value, np.ndarray):
        return ("ndarray", value.dtype.str)
    if isinstance(value, (tuple, list)):
        return (type(value).__name__, tuple(_type_tree(x) for x in value))
    if isinstance(value, Mapping) and not isinstance(value, ContentIdentity):
        return ("mapping", tuple(sorted((repr(k), _type_tree(v)) for k, v in value.items())))
    return f"{type(value).__module__}.{type(value).__qualname__}"


class _ConstructorPickler(pickle.Pickler):
    """Every dataclass crosses into the generating process through its constructor, so its construction
    runs, and is traced, there (census finding F1: a default unpickle restores the built object, and in 2 of
    2 traced trajectory runs its construction was absent from the trace)."""

    def reducer_override(self, obj: Any) -> Any:
        if dataclasses.is_dataclass(obj) and not isinstance(obj, (type, ContentIdentity)):
            return (_rebuild, (type(obj), {f.name: getattr(obj, f.name) for f in dataclasses.fields(obj) if f.init}))
        return NotImplemented


def _constructor_pickle(value: Any) -> bytes:
    buffer = io.BytesIO()
    _ConstructorPickler(buffer, protocol=pickle.HIGHEST_PROTOCOL).dump(value)
    return buffer.getvalue()


#: The generating process's first code: the tracer, the audit hook on ``open``, then this module. It is run
#: by path with ``-P``, so its own directory is not put on ``sys.path`` (``orpheus/numerics/operator.py``
#: would otherwise shadow the standard library's ``operator``).
_BOOT = Path(__file__).with_name("_traced_memo_boot.py")


@dataclasses.dataclass(frozen=True)
class _Job:
    """One generation, as the generating process receives it."""

    memo: TracedMemo
    args: tuple[Any, ...]
    kwargs: dict[str, Any]
    key: str
    root: str


def _generate_here(job: _Job, codes: set[types.CodeType], stop_recording: Callable[[], None], opened: Sequence[str]) -> dict[str, Any]:
    """The generating process's body (called by the boot script): run the call, stop recording, write the
    entry; the result crosses back to the parent as ``{"ok": True}`` or ``{"ok": False, "error": ...}``."""
    store = Store(Path(job.root))
    _ROOT.set(store.root)
    children: list[ChildPin] = []
    _CHILDREN.set(children)
    try:
        value = job.memo.__wrapped__(*job.args, **job.kwargs)
        stop_recording()
        store.write(job.memo.function_id, job.key, _manifest(codes, children, opened, store), value)
        return {"ok": True}
    except BaseException as error:  # noqa: BLE001 - every failure crosses back to the caller
        stop_recording()
        try:
            pickle.dumps(error)
        except Exception:  # noqa: BLE001 - an exception that cannot cross is described instead
            error = RuntimeError(f"{type(error).__name__}: {error}")
        return {"ok": False, "error": error}


class TracedMemo:
    """A pure function of content-identified arguments, memoised on disk and generated in a fresh process."""

    def __init__(self, function: Callable[..., Any], canonical: Mapping[str, Callable[[Any], Any]] | None = None) -> None:
        functools.update_wrapper(self, function)
        self.__wrapped__ = function
        self.function_id = f"{function.__module__}:{function.__qualname__}"
        self._signature = inspect.signature(function)
        self._canonical = dict(canonical or {})

    def __reduce__(self) -> Any:
        return (_resolve, (self.function_id,))

    def __get__(self, instance: object, owner: type | None = None) -> Any:
        """A memoised method binds like a function: ``instance.method(x)`` is the memo called with
        ``(instance, x)``, so the instance's content is part of the key (it must have content identity)."""
        return self if instance is None else types.MethodType(self, instance)

    def _arguments(self, args: tuple[Any, ...], kwargs: dict[str, Any]) -> inspect.BoundArguments:
        """The call bound to the signature, defaults applied and every declared canonical form put in: the key
        and the generating process see the same arguments."""
        bound = self._signature.bind(*args, **kwargs)
        bound.apply_defaults()
        for name, canonical in self._canonical.items():
            if bound.arguments.get(name) is not None:
                bound.arguments[name] = canonical(bound.arguments[name])
        return bound

    def _key(self, bound: inspect.BoundArguments) -> str:
        arguments = dict(bound.arguments)
        return _digest(encode((self.function_id, arguments, _type_tree(arguments), platform_tag())))

    def key(self, *args: Any, **kwargs: Any) -> str:
        return self._key(self._arguments(args, kwargs))

    def lookup(self, *args: Any, **kwargs: Any) -> Verdict:
        """The entry's verdict for this call, without generating anything."""
        return Store.current().lookup(self.function_id, self.key(*args, **kwargs))

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        if _BYPASS.get():
            return self.__wrapped__(*args, **kwargs)
        bound = self._arguments(args, kwargs)
        key = self._key(bound)
        store = Store.current()
        verdict = store.lookup(self.function_id, key)
        if not isinstance(verdict, Hit):
            self._generate(_Job(self, bound.args, bound.kwargs, key, str(store.root)))
            verdict = store.lookup(self.function_id, key)
            if not isinstance(verdict, Hit):
                raise RuntimeError(f"traced_memo: the entry just generated for {self.function_id} does not validate: {verdict}")
        if (children := _CHILDREN.get()) is not None:
            children.append(ChildPin(self.function_id, key, verdict.payload_digest))
        return store.load(verdict, _returned_types(self.__wrapped__))

    def _generate(self, job: _Job) -> None:
        """Run ``job`` in a fresh interpreter with this process's ``-O`` level, ``sys.path`` and working
        directory, and its environment without any ``ORPHEUS_*`` variable (a withdrawn generator therefore
        always refuses there, and no variable selects an answer no key holds)."""
        optimize = ["-" + "O" * sys.flags.optimize] if sys.flags.optimize else []
        run = subprocess.run(
            [sys.executable, *optimize, "-P", str(_BOOT)],
            input=pickle.dumps({"sys_path": list(sys.path), "job": _constructor_pickle(job)}),
            capture_output=True, cwd=os.getcwd(),
            env={k: v for k, v in os.environ.items() if not k.startswith("ORPHEUS_")},
        )
        try:
            result = pickle.loads(run.stdout)
        except Exception as error:  # noqa: BLE001 - the process died before it could answer
            raise RuntimeError(f"traced_memo: the process generating {self.function_id} failed:\n{run.stderr.decode()[-4000:]}") from error
        if not result["ok"]:
            error = result["error"]
            error.add_note(f"raised in the process generating {self.function_id}")
            raise error


def traced_memo(function: Callable[..., Any] | None = None, *, canonical: Mapping[str, Callable[[Any], Any]] | None = None) -> Any:
    """Memoise ``function`` on disk (:class:`TracedMemo`); ``canonical`` maps a parameter to the function that
    puts its argument in canonical form (``{"radii": as_float_array}``), for the key and the generation alike."""
    if function is None:
        return lambda f: TracedMemo(f, canonical)
    return TracedMemo(function, canonical)
