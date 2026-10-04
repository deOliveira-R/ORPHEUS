r"""The traced memo: a pure function of content, memoised on disk and keyed on what ran (#405 P3).

A reference reading costs seconds to minutes and changes only when something that produces it changes. The
memo stores each call's answer under ``.cache/references/`` and serves it again for as long as nothing the
call depended on has changed. The dependencies are not declared, they are RECORDED: on a miss the call is
generated in a fresh interpreter recording from its first line (:mod:`~orpheus.numerics._traced_memo_boot`),
and the entry keeps the :class:`Manifest` of that run beside the answer.

* **The key** is the digest of the function's identity, its signature-bound arguments (defaults applied,
  each declared canonical form put in) and the platform tag, all through
  :func:`~orpheus.numerics.content.encode_exact`: the arguments' content walked exactly, every leaf with its
  type and bits. Content equality follows ``==`` (``8 == 8.0``, ``-0.0 == 0.0``), and a function can tell
  those apart, so a key that followed ``==`` could serve one call another's answer.
* **The manifest** is a set of :data:`Pin`\ s, each able to say whether the checkout still matches it:
  each first-party def that ran by the digest of its normalised source (docstrings, comments and layout
  excluded), each first-party module that ran by its skeleton (every function body removed: imports,
  constants, decorators, signatures, class attributes), each third-party distribution by its version, the
  interpreter, each file the run read and each directory it listed (a missing one included: the run
  depended on its absence), the working directory when the run read a relative path, each memo entry the
  run read (validated recursively), and the memo's own source (it wrote the entry).
* **Validation** re-hashes the manifest against the current files; it imports and runs nothing. A changed
  pin makes the entry :class:`Stale`, an unreadable or inconsistent entry :class:`Corrupt`; neither is
  served, and the next call regenerates it.
* **The payload** is exact: JSON with every float as ``float.hex`` and every array in the entry's ``.npz``,
  admitting only the types it writes exactly (a subclass of one is refused, never written as its base);
  never a pickle. Every load hands out fresh, read-only arrays. An entry is ONE file, replaced atomically.

Why a fresh process (the user's ruling of 2026-10-04, ``.claude/plans/reference_cache.md``, "P3 rulings"):
a recording sees only code that RUNS, so anything served from memory (a ``functools.cache`` hit, a
``cached_property`` computed earlier, a monkeypatched function) would hide a dependency from it. A process
born for one call has nothing in memory, so the recording is complete without any discipline on how the
rest of the package memoises. The cost is one interpreter start per miss, about 0.65 s ``[M]`` 2026-10-04.

**The generator contract** (the user's ruling of 2026-10-04, after qa's third review: a recorder written in
Python cannot be airtight against arbitrary Python, and each review found another way a deliberately
adversarial generator could hide a dependency). A memoised function:

* reads data only by opening real files (``open``, ``Path.read_*``, ``np.load``), never through a loader
  (``pkgutil.get_data``) or a module's absence (``try: import accel``, ``find_spec(...) is None``);
* depends on a file's presence and its bytes, never on its other metadata (its size or time from ``stat``,
  whether it is a link);
* reads no environment variable (the M5.1 census holds this for ``orpheus/``), starts no process, and
  branches on no other memo's ``lookup``;
* takes arguments with a constructor form: no ``InitVar``, no mapping beside ``dict`` and ``FrozenMapping``,
  no sparse matrix without stored arrays (each refused when keyed).

SCOPE-BOUNDARY[guard] machinery: a recorder below Python (a system-call tracer) for what the contract
excludes. ruling: the user, 2026-10-04 (``.claude/plans/reference_cache.md``, "P3 closed by contract").
revisit: a generator that needs what the contract excludes, or a cold rebuild (P5: every entry regenerated
without the cache and compared) that finds a stale entry. That cold rebuild is the contract's witness.

What no manifest can see, and how each is held: a probe the recorder does not wrap (``os.access``,
``os.scandir``'s entries); native state (a C library's own globals) and a file a
C library opens itself (HDF5 through ``h5py``) are invisible, so a generator reading such a file takes the
file's digest as an argument. A run may start only a declared machine query (``uname``), pinned by its executable's bytes; any other
started process is refused, since its own reads and children are unrecorded. A run sees only a declared
environment, which is pinned. A ``bypass()`` entered in one thread holds for every thread of the process.
A test that monkeypatches anything a generation runs reads under :func:`bypass`.
"""
from __future__ import annotations

import ast
import contextlib
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
import site
import subprocess
import sys
import sysconfig
import traceback
import types
import typing
import uuid
import zipfile
from collections.abc import Callable, Generator, Iterable, Mapping
from pathlib import Path
from typing import Any, Generic, NamedTuple, ParamSpec, TypeVar, overload

import numpy as np

from orpheus.numerics._traced_memo_boot import Recording, note_child, spawning, unrecorded
from orpheus.numerics.content import ContentIdentity, _rebuild, constructor_arguments, encode_exact

__all__ = [
    "Absent", "ChildPin", "Corrupt", "DataPin", "DefPin", "DistributionPin", "EnvironmentPin", "Hit", "InterpreterPin", "PresencePin",
    "ListingPin", "Manifest", "MemoPin", "ModulePin", "Pin", "Stale", "TracedMemo", "Unencodable",
    "Unpinnable", "Verdict", "WorkingDirectoryPin", "bypass", "cache_root", "decode_payload", "default_root",
    "encode_payload", "function_digest", "origin", "platform_tag", "skeleton_digest", "trace_call",
    "traced_memo", "validate",
]


def _digest(data: bytes) -> str:
    return hashlib.blake2b(data, digest_size=32).hexdigest()


class Unpinnable(RuntimeError):
    """A run depends on something no manifest can pin: a site-packages file of no distribution, code compiled
    under a file's name at a line where the file holds no such code, or a process the run started."""


class Unencodable(TypeError):
    """A return value the payload cannot write exactly; it is refused, never written lossily."""


# ── the normalised source ───────────────────────────────────────────────────────────────────────────

_DEF = (ast.FunctionDef, ast.AsyncFunctionDef)


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


def _span(node: ast.AST) -> range:
    first = min([node.lineno, *(d.lineno for d in getattr(node, "decorator_list", []))])  # type: ignore[attr-defined]
    return range(first, node.end_lineno + 1)  # type: ignore[attr-defined]


class _Def(NamedTuple):
    """An outermost def: its qualified name, its rank among the defs of that name (0 first), its span with
    decorators, its node and the digest of its normalised source."""

    qualname: str
    ordinal: int
    span: range
    node: ast.FunctionDef | ast.AsyncFunctionDef
    digest: str


@dataclasses.dataclass(frozen=True)
class _SourceText:
    """One source file, parsed once: its outermost defs and its skeleton's digest."""

    tree: ast.Module
    defs: tuple[_Def, ...]
    skeleton: str

    def digest(self, qualname: str, ordinal: int = -1) -> str | None:
        """The digest of the ``ordinal``-th def named ``qualname`` (``-1``: the last, the one Python binds)."""
        named = [d for d in self.defs if d.qualname == qualname]
        return named[ordinal].digest if -len(named) <= ordinal < len(named) else None

    def enclosing(self, line: int) -> _Def | None:
        """The outermost def whose span holds ``line`` (outermost defs never overlap)."""
        return next((d for d in self.defs if line in d.span), None)


@functools.cache
def _source_text(data: bytes) -> _SourceText:
    """The parse of a source text, cached by its BYTES (an edited file is another key, never a stale read).

    ONE walk finds every outermost def: it descends every statement except a def's body, so a def under any
    block (``if``, ``for``, ``while``, ``match``, ``try``, ``with``) is found; a class adds its name to the
    qualified names of the defs in its body."""
    tree = _without_docstrings(ast.parse(data))
    found: list[_Def] = []
    seen: dict[str, int] = {}

    def walk(parent: ast.AST, prefix: str) -> None:
        for node in ast.iter_child_nodes(parent):
            if isinstance(node, _DEF):
                qualname = prefix + node.name
                seen[qualname] = seen.get(qualname, -1) + 1
                found.append(_Def(qualname, seen[qualname], _span(node), node, _ast_digest(node)))
            elif isinstance(node, ast.ClassDef):
                walk(node, f"{prefix}{node.name}.")
            elif isinstance(node, (ast.stmt, ast.excepthandler, ast.match_case)):
                walk(node, prefix)

    walk(tree, "")
    skeleton = _without_docstrings(ast.parse(data))
    for node in ast.walk(skeleton):
        if isinstance(node, _DEF):
            node.body = [ast.Pass()]
    return _SourceText(tree, tuple(found), _ast_digest(skeleton))


def function_digest(source: bytes | str, qualname: str, ordinal: int = -1) -> str | None:
    """The digest of the ``ordinal``-th outermost def ``qualname`` (``function`` or ``Class.method``; ``-1``,
    the default, is the last, which Python binds), docstrings stripped; ``None`` when there is none."""
    return _source_text(source.encode() if isinstance(source, str) else source).digest(qualname, ordinal)


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


# ── where a file a run touched belongs ──────────────────────────────────────────────────────────────


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
    """Generated code (``<string>``, ``<frozen …>``): pinned by the interpreter and the skeleton."""

    filename: str


class Interpreter(NamedTuple):
    """A standard-library file: pinned by the interpreter."""

    path: str


class Distribution(NamedTuple):
    """A file of an installed distribution: pinned by the distribution's version."""

    name: str


class Source(NamedTuple):
    """A file pinned by its content: first-party source, an editable install's own file, a data file."""

    path: str


Origin = Dropped | Interpreter | Distribution | Source


def origin(filename: str) -> Origin:
    """Where a file a run touched belongs, and so how the manifest pins it."""
    if filename.startswith("<"):
        return Dropped(filename)
    real = _real(filename)
    for site_dir in _site_dirs():
        if real.startswith(site_dir):
            relative = real[len(site_dir):]
            top = relative.split(os.sep)[0].removesuffix(".py")
            if top.endswith(".dist-info"):  # a distribution's own metadata, read by importlib.metadata
                return Distribution(importlib.metadata.Distribution.at(Path(site_dir) / top).metadata["Name"])
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


def _pin_path(real: str) -> str:
    """How a source file is pinned: relative to the ``sys.path`` entry that holds it, so validation finds the
    file the same import would load now; absolute when no entry holds it (an import through an editable
    install's finder, qa finding 12 of 2026-10-04)."""
    for entry in sys.path:
        root = _real(entry or os.getcwd()) + os.sep
        if real.startswith(root):
            return real[len(root):]
    return real


def _locate(path: str) -> Path | None:
    if os.path.isabs(path):
        return Path(path) if os.path.isfile(path) else None
    for entry in sys.path:
        candidate = Path(entry or os.getcwd()) / path
        if candidate.is_file():
            return candidate
    return None


def _python_identity() -> str:
    return f"{sys.version}|{sys.implementation.cache_tag}"


@functools.cache
def platform_tag() -> str:
    """The operating system, the machine, the ``-O`` level (a bare ``assert`` in a generator is stripped under
    it) and the BLAS (Accelerate and OpenBLAS differ in the last bits, #504)."""
    config = np.show_config(mode="dicts")
    blas = config.get("Build Dependencies", {}).get("blas", {}).get("name", "unknown") if isinstance(config, dict) else "unknown"
    return f"{sys.implementation.cache_tag}-{sys.platform}-{platform.machine()}-O{sys.flags.optimize}-{blas}"


# ── the manifest ────────────────────────────────────────────────────────────────────────────────────

#: The memo's own source and the code its write path runs after the recording stops (the payload's
#: constructor form is :func:`~orpheus.numerics.content.constructor_arguments`): an edit to any of it makes
#: every entry stale (qa finding 8, and its second review, 2026-10-04).
_MEMO_SOURCES = tuple(Path(__file__).resolve().with_name(name) for name in ("traced_memo.py", "_traced_memo_boot.py", "content.py"))

#: The programs a run may start, each pinned by its executable's bytes: queries whose answer the machine fixes.
#: ``uname`` is run by ``platform.processor()`` during the imports of the P3 clients (``[M]`` 2026-10-04); it is
#: admitted by its system path, never by its name (a script named ``uname`` earlier on ``PATH`` is refused).
#: SCOPE-BOUNDARY[guard] machinery: a recorder for a started program's own reads and children (a traced
#: sub-process). ruling: the orchestrator, #405 P3, from qa's second review (a ``cat`` of an unpinned file, a
#: ``shell=True`` and a ``/usr/bin/env python3`` each served a stale value when every program was admitted
#: by its bytes). revisit: when a generator must start a program outside this table.
_ADMITTED_PROGRAMS = frozenset(os.path.realpath(p) for p in ("/usr/bin/uname", "/bin/uname") if os.path.isfile(p))

#: The environment a generating process receives, and nothing else: its values are pinned
#: (:class:`EnvironmentPin`), and ``PYTHONHASHSEED`` is fixed so that iteration over a set of strings is one
#: order in every generation. A variable outside this list (an ``ORPHEUS_*`` switch, a shell's own) never
#: reaches a generation, so no answer can follow a value no key holds (qa's second review, 2026-10-04).
_CHILD_ENVIRONMENT = ("PATH", "HOME", "TMPDIR", "LANG", "LC_ALL", "LC_CTYPE", "OMP_NUM_THREADS",
                      "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
_FIXED_ENVIRONMENT = {"PYTHONHASHSEED": "0"}
#: The declared variables a generation's answer can follow, and so the ones its pin holds. ``PATH`` is passed
#: (``uname`` is found through it) but not pinned: the one program a generation may start is admitted by its
#: system path and pinned by its bytes, and pinning ``PATH`` made every entry stale when a virtual environment
#: was activated (qa's third review, 2026-10-04).
_PINNED_ENVIRONMENT = tuple(k for k in _CHILD_ENVIRONMENT if k != "PATH") + tuple(_FIXED_ENVIRONMENT)

#: The deepest chain of generations one call may start: a cycle through distinct keys (``f(n)`` calling
#: ``f(n + 1)``) has no repeated key for the ancestry to catch.
_MAX_GENERATION_DEPTH = 16


def _child_environment(environment: Mapping[str, str] | None = None) -> dict[str, str]:
    """What a generation receives from ``environment`` (default: this process's): the declared variables and
    the fixed ones. A generating process projects to itself; ``trace_call`` projects its caller's."""
    source = os.environ if environment is None else environment
    return {k: source[k] for k in _CHILD_ENVIRONMENT if k in source} | _FIXED_ENVIRONMENT


class _Checkout:
    """The source files as one validation reads them (no read outlives it), and the store children live in."""

    def __init__(self, store: Store) -> None:
        self.store = store
        self._texts: dict[str, bytes | None] = {}

    def text(self, path: str) -> bytes | None:
        if path not in self._texts:
            found = _locate(path)
            self._texts[path] = found.read_bytes() if found is not None else None
        return self._texts[path]


def _file_digest(path: str) -> str:
    """The digest of a file's bytes, or ``"absent"``: a run that found no file there depended on that."""
    return _digest(Path(path).read_bytes()) if os.path.isfile(path) else "absent"


def _presence(path: str) -> str:
    """What a probe of ``path`` finds: a file, a directory, something else, or nothing."""
    if os.path.isfile(path):
        return "file"
    if os.path.isdir(path):
        return "directory"
    return "other" if os.path.lexists(path) else "absent"


def _listing_digest(path: str) -> str:
    return _digest("\0".join(sorted(os.listdir(path))).encode()) if os.path.isdir(path) else "absent"


def _memo_digest() -> str:
    return _digest(b"".join(p.read_bytes() for p in _MEMO_SOURCES))


class DefPin(NamedTuple):
    """A first-party def that ran: its file, its qualified name, its rank among the defs of that name, and
    the digest of its normalised source."""

    path: str
    qualname: str
    ordinal: int
    digest: str

    def stale(self, checkout: _Checkout) -> str | None:
        text = checkout.text(self.path)
        found = None if text is None else function_digest(text, self.qualname, self.ordinal)
        if found is None:
            return f"function {self.path}:{self.qualname} gone"
        return f"function {self.path}:{self.qualname} changed" if found != self.digest else None


class ModulePin(NamedTuple):
    """A first-party module that ran, by its skeleton."""

    path: str
    digest: str

    def stale(self, checkout: _Checkout) -> str | None:
        text = checkout.text(self.path)
        if text is None:
            return f"module {self.path} gone"
        return f"module {self.path} skeleton changed" if skeleton_digest(text) != self.digest else None


class DistributionPin(NamedTuple):
    """A third-party distribution whose code ran or whose file was read, by its version."""

    name: str
    version: str

    def stale(self, checkout: _Checkout) -> str | None:
        try:
            now: str | None = importlib.metadata.version(self.name)
        except importlib.metadata.PackageNotFoundError:
            now = None
        return f"distribution {self.name}: {self.version} != {now}" if now != self.version else None


class InterpreterPin(NamedTuple):
    """The interpreter that ran (its version and bytecode tag): the standard library and generated code."""

    identity: str

    def stale(self, checkout: _Checkout) -> str | None:
        return f"python: {self.identity!r} != {_python_identity()!r}" if self.identity != _python_identity() else None


class DataPin(NamedTuple):
    """A file the run opened for reading, by its absolute path and its bytes (``"absent"``: it was missing)."""

    path: str
    digest: str

    def stale(self, checkout: _Checkout) -> str | None:
        return f"data {self.path} changed" if _file_digest(self.path) != self.digest else None


class ListingPin(NamedTuple):
    """A directory the run listed, by its absolute path and the names in it (``"absent"``: it was missing)."""

    path: str
    digest: str

    def stale(self, checkout: _Checkout) -> str | None:
        return f"listing {self.path} changed" if _listing_digest(self.path) != self.digest else None


class PresencePin(NamedTuple):
    """A path the run probed (``exists``, ``isfile``, ``stat``) without opening it, by what was there."""

    path: str
    kind: str

    def stale(self, checkout: _Checkout) -> str | None:
        return f"presence {self.path}: {self.kind} became {_presence(self.path)}" if _presence(self.path) != self.kind else None


class EnvironmentPin(NamedTuple):
    """The environment the generating process received (the declared variables only), by its digest."""

    digest: str

    def stale(self, checkout: _Checkout) -> str | None:
        return "the generation's environment changed" if _environment_digest() != self.digest else None


class WorkingDirectoryPin(NamedTuple):
    """The working directory a run started in, when it read a relative path (from another, the same path is
    another file); the starting one, since a run that changes directory does so deterministically from it."""

    path: str

    def stale(self, checkout: _Checkout) -> str | None:
        return f"working directory {self.path} != {os.getcwd()}" if os.getcwd() != self.path else None


class ChildPin(NamedTuple):
    """A memo entry the run read, hit or generated, by its key and the digest of the payload it served."""

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


class MemoPin(NamedTuple):
    """The memo's own source, which wrote the entry after the recording stopped (qa finding 8, 2026-10-04)."""

    digest: str

    def stale(self, checkout: _Checkout) -> str | None:
        return "the traced memo's own source changed" if _memo_digest() != self.digest else None


Pin = (DefPin | ModulePin | DistributionPin | InterpreterPin | DataPin | PresencePin | ListingPin | EnvironmentPin
       | WorkingDirectoryPin | ChildPin | MemoPin)


def _environment_digest(environment: Mapping[str, str] | None = None) -> str:
    """The digest of the pinned variables of ``environment`` (default: what this process would hand a
    generation now)."""
    received = _child_environment(environment)
    return _digest(json.dumps(sorted((k, received[k]) for k in _PINNED_ENVIRONMENT if k in received)).encode())
_PIN_KINDS: Mapping[str, Callable[..., Pin]] = {kind.__name__: kind for kind in typing.get_args(Pin)}


@dataclasses.dataclass(frozen=True)
class Manifest:
    """What one generation depended on, pinned so that a later checkout can be compared with it without
    running it. JSON rows are tagged by their pin's kind and parsed once, in :meth:`from_json`."""

    pins: frozenset[Pin]

    def __post_init__(self) -> None:
        strays = [p for p in self.pins if not isinstance(p, typing.get_args(Pin))]
        if strays:
            raise TypeError(f"Manifest: every member is a pin, got {strays[:3]!r}")

    def of(self, kind: type[Any]) -> tuple[Any, ...]:
        """The pins of one kind, sorted."""
        return tuple(sorted(p for p in self.pins if isinstance(p, kind)))

    @property
    def functions(self) -> tuple[DefPin, ...]:
        return self.of(DefPin)

    @property
    def modules(self) -> tuple[ModulePin, ...]:
        return self.of(ModulePin)

    @property
    def distributions(self) -> tuple[DistributionPin, ...]:
        return self.of(DistributionPin)

    @property
    def children(self) -> tuple[ChildPin, ...]:
        return self.of(ChildPin)

    @property
    def data(self) -> tuple[DataPin, ...]:
        return self.of(DataPin)

    def replacing(self, kind: type[Any], pins: Iterable[Pin]) -> Manifest:
        """This manifest with every pin of ``kind`` replaced by ``pins``."""
        return Manifest(frozenset({p for p in self.pins if not isinstance(p, kind)} | set(pins)))

    def to_json(self) -> list[list[Any]]:
        return sorted([type(pin).__name__, *pin] for pin in self.pins)

    @classmethod
    def from_json(cls, rows: list[list[Any]]) -> Manifest:
        return cls(frozenset(_PIN_KINDS[kind](*fields) for kind, *fields in rows))


def _def_pin(text: _SourceText, path: str, code: types.CodeType) -> DefPin | None:
    """The pin of the outermost def whose source holds ``code``; ``None`` for module-level code (a module or
    class body, a module-level lambda, a signature), which the skeleton pins."""
    line, name = code.co_firstlineno, code.co_name
    found = text.enclosing(line)
    if name == "<module>" or (found is not None and name.startswith(_SIGNATURE_SCOPES) and line < found.node.body[0].lineno):
        return None
    if found is None:
        if _holds(text.tree, line, name):
            return None
    elif _holds(found.node, line, name):
        return DefPin(path, found.qualname, found.ordinal, found.digest)
    raise Unpinnable(f"{path}:{line} {code.co_qualname}: no def of that name at its line")


def _manifest(recording: Recording, store: Store) -> Manifest:
    """The manifest of a recorded run: every code object it started, file it read, directory it listed and
    memo entry it read, each pinned by its origin."""
    pins: set[Pin] = {InterpreterPin(_python_identity()), MemoPin(_memo_digest()), EnvironmentPin(_environment_digest(recording.environment))}
    pins |= {_program_pin(program) for program in recording.spawned}
    pins |= {ChildPin(*child) for child in recording.children}
    by_file: dict[str, list[types.CodeType]] = {}
    for code in recording.codes:
        match origin(code.co_filename):
            case Distribution(name):
                pins.add(DistributionPin(name, importlib.metadata.version(name)))
            case Source(path):
                by_file.setdefault(path, []).append(code)
            case Dropped() | Interpreter():
                pass
    for real, file_codes in by_file.items():
        path = _pin_path(real)
        text = _source_text(Path(real).read_bytes())
        pins.add(ModulePin(path, text.skeleton))
        pins |= {pin for code in file_codes if (pin := _def_pin(text, path, code)) is not None}
    store_root = _real(str(store.root)) + os.sep
    for opened in recording.opened:
        if _real(opened).startswith(store_root):
            continue  # a store entry is a child
        match origin(opened):
            case Distribution(name):
                pins.add(DistributionPin(name, importlib.metadata.version(name)))
            case Source():
                pins.add(DataPin(opened, _file_digest(opened)))
            case Dropped() | Interpreter():
                pass
    pins |= {PresencePin(probed, _presence(probed)) for probed in recording.probed if not _real(probed).startswith(store_root)}
    pins |= {ListingPin(listed, _listing_digest(listed)) for listed in recording.listed}
    if recording.relative:
        pins.add(WorkingDirectoryPin(recording.directory))
    return Manifest(frozenset(pins))


def _program_pin(program: str) -> DataPin:
    """A program the run started, pinned by its executable's bytes, when it is one of
    :data:`_ADMITTED_PROGRAMS` (a query whose answer the machine fixes); any other start is refused: a
    Python process (a process pool, ``env python3``) ran code no recording saw, a shell or another program
    read files and started processes no recording saw, and a ``fork`` or an ``os.system`` line names no
    program at all."""
    if program in _ADMITTED_PROGRAMS:
        return DataPin(program, _file_digest(program))
    raise Unpinnable(
        f"the run started a process ({program}) whose own reads and children no recording sees; a generation may "
        f"start only {sorted(_ADMITTED_PROGRAMS)} (a scope boundary of the traced memo)"
    )


def validate(manifest: Manifest, store: Store | None = None) -> tuple[str, ...]:
    """Every reason ``manifest`` no longer describes the checkout (empty: it does). Imports and runs nothing."""
    checkout = _Checkout(store if store is not None else Store.current())
    ordered = sorted(manifest.pins, key=lambda pin: (type(pin).__name__, tuple(pin)))
    return tuple(reason for pin in ordered if (reason := pin.stale(checkout)) is not None)


def trace_call(function: Callable[..., Any], *args: Any, **kwargs: Any) -> tuple[Any, Manifest]:
    """Call ``function`` in THIS process under the recorder a generation uses: its value and the manifest of
    what it depended on. Code that ran earlier in the process and is served from memory now is not in it,
    which is why a generation runs in a fresh process; this is the recorder's own instrument."""
    recording = Recording().start()
    try:
        value = function(*args, **kwargs)
    finally:
        recording.stop()
    return value, _manifest(recording, Store.current())


# ── the payload ─────────────────────────────────────────────────────────────────────────────────────

#: The array and scalar kinds the payload writes exactly: boolean, signed and unsigned integer, real, complex.
_EXACT_KINDS = "biufc"


def _dotted(cls: type) -> str:
    return f"{cls.__module__}:{cls.__qualname__}"


def encode_payload(value: Any, allowed: set[str]) -> tuple[Any, dict[str, np.ndarray]]:
    """``(JSON tree, arrays)``: floats as ``float.hex``, numpy scalars with their dtype, arrays by name, tuples,
    and the dataclasses named in ``allowed`` by their constructor form. A type is admitted EXACTLY: a subclass
    of one (a ``NamedTuple``, an ``IntEnum``, a masked array) would lose what makes it the subclass (qa finding
    3 of 2026-10-04)."""
    arrays: dict[str, np.ndarray] = {}

    def node(v: Any) -> Any:
        kind = type(v)
        if v is None or kind in (bool, str):
            return {"v": v}
        if kind is int:
            return {"i": str(v)}
        if kind is float:
            return {"f": float.hex(v)}
        if isinstance(v, np.generic) and kind is v.dtype.type and v.dtype.kind in _EXACT_KINDS:
            return {"np": v.dtype.str, "a": node(np.asarray(v))["a"]}
        if isinstance(v, np.ndarray) and kind is np.ndarray and v.dtype.kind in _EXACT_KINDS:
            name = f"a{len(arrays)}"
            arrays[name] = v
            return {"a": name}
        if isinstance(v, tuple) and kind is tuple:
            return {"t": [node(x) for x in v]}
        if dataclasses.is_dataclass(v) and _dotted(kind) in allowed:
            return {"dc": _dotted(kind), "fields": {k: node(x) for k, x in constructor_arguments(v).items()}}
        raise Unencodable(
            f"a {kind.__module__}.{kind.__qualname__} has no exact payload (admitted exactly: None, bool, int, float, "
            f"str, numpy scalars and arrays of kind {_EXACT_KINDS!r}, tuples, and the dataclasses the function's "
            f"return annotation names: {sorted(allowed)})"
        )

    return node(value), arrays


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
                return _rebuild(_resolve(dotted), {k: value(x) for k, x in fields.items()})
            case _:
                raise ValueError(f"an unknown payload node {sorted(v)}")

    return value(tree)


def _payload_digest(tree: Any, arrays: Mapping[str, np.ndarray]) -> str:
    """The digest of a payload's CONTENT (its tree, each array's name, dtype, shape and bytes), not of the
    container that stores it."""
    parts = [json.dumps(tree, sort_keys=True).encode()]
    for name in sorted(arrays):
        array = np.ascontiguousarray(arrays[name])
        parts.append(f"{name}|{array.dtype.str}|{array.shape}|".encode() + array.tobytes())
    return _digest(b"\0".join(parts))


# ── the store ───────────────────────────────────────────────────────────────────────────────────────


@dataclasses.dataclass(frozen=True)
class Hit:
    """A valid entry: the digest of its payload, its tree and its arrays, as verified (qa finding 7: the value
    served is the one the digest was checked on, never a second read)."""

    payload_digest: str
    tree: Any
    arrays: Mapping[str, np.ndarray]


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

#: The member of an entry's ``.npz`` that holds its JSON (the key, the manifest, the payload tree, its digest).
_ENTRY = "__entry__"


def default_root() -> Path:
    """``.cache/references`` at the repository root (gitignored)."""
    return Path(__file__).resolve().parents[2] / ".cache" / "references"


@dataclasses.dataclass(frozen=True)
class Store:
    """The entries under one root, one file each: ``<root>/<function id>/<key>.npz``."""

    root: Path

    @classmethod
    def current(cls) -> Store:
        target = _target()
        return target if isinstance(target, Store) else cls(default_root())

    def path(self, function_id: str, key: str) -> Path:
        return self.root / function_id / f"{key}.npz"

    def lookup(self, function_id: str, key: str) -> Verdict:
        with unrecorded():
            return self._lookup(function_id, key)

    def _lookup(self, function_id: str, key: str) -> Verdict:
        try:
            data = self.path(function_id, key).read_bytes()
        except FileNotFoundError:
            return Absent()
        try:
            with np.load(io.BytesIO(data), allow_pickle=False) as npz:
                arrays = {name: npz[name] for name in npz.files}
            entry = json.loads(arrays.pop(_ENTRY).tobytes())
            if _payload_digest(entry["payload"], arrays) != entry["payload_digest"]:
                return Corrupt("the payload does not match its digest")
            manifest = Manifest.from_json(entry["manifest"])
        except (OSError, ValueError, KeyError, TypeError, zipfile.BadZipFile) as error:
            return Corrupt(f"{type(error).__name__}: {error}")
        if reasons := validate(manifest, self):
            return Stale(reasons)
        return Hit(entry["payload_digest"], entry["payload"], arrays)

    def write(self, function_id: str, key: str, manifest: Manifest, value: Any, allowed: set[str]) -> None:
        """Write the entry to a temporary file beside it, then move it into place with one ``os.replace``: a
        reader sees the old entry or the new one, never a part of either (qa finding 10 of 2026-10-04)."""
        tree, arrays = encode_payload(value, allowed)
        entry = {"function": function_id, "key": key, "manifest": manifest.to_json(),
                 "payload": tree, "payload_digest": _payload_digest(tree, arrays)}
        path = self.path(function_id, key)
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
        with temporary.open("wb") as handle:
            np.savez(handle, allow_pickle=False, **arrays, **{_ENTRY: np.frombuffer(json.dumps(entry).encode(), np.uint8)})
        os.replace(temporary, path)


class InProcess:
    """The target of a call under :func:`bypass`: it runs in this process and reads and writes nothing."""


#: Where memoised calls go, innermost last: a :class:`Store`, or :class:`InProcess` under :func:`bypass`.
#: Process state, not context: a memo called from a worker thread reads the same target (qa finding 2).
_TARGETS: list[Store | InProcess] = []


def _target() -> Store | InProcess:
    return _TARGETS[-1] if _TARGETS else Store(default_root())


@contextlib.contextmanager
def _targeting(target: Store | InProcess) -> Generator[None]:
    _TARGETS.append(target)
    try:
        yield
    finally:
        _TARGETS.pop()


@contextlib.contextmanager
def cache_root(path: Path) -> Generator[Path]:
    """Read and write the store under ``path`` while the block runs."""
    with _targeting(Store(Path(path))):
        yield Path(path)


@contextlib.contextmanager
def bypass() -> Generator[None]:
    """Every memoised call runs in THIS process and reads and writes nothing: the spelling for a test that
    monkeypatches anything a generation runs (a patch never reaches a fresh process, so without the bypass the
    test would read the honest entry instead of its patched answer)."""
    with _targeting(InProcess()):
        yield


# ── the key and the process boundary ────────────────────────────────────────────────────────────────


class _ConstructorPickler(pickle.Pickler):
    """Every dataclass crosses into the generating process through its constructor, so its construction
    runs, and is recorded, there (census finding F1: a default unpickle restores the built object, and in 2
    of 2 traced trajectory runs its construction was absent from the trace)."""

    def reducer_override(self, obj: Any) -> Any:
        if dataclasses.is_dataclass(obj) and not isinstance(obj, (type, ContentIdentity)):
            return (_rebuild, (type(obj), constructor_arguments(obj)))
        return NotImplemented


def _constructor_pickle(value: Any) -> bytes:
    buffer = io.BytesIO()
    _ConstructorPickler(buffer, protocol=pickle.HIGHEST_PROTOCOL).dump(value)
    return buffer.getvalue()


#: The generating process's first code, run by path with ``-P`` (``orpheus/numerics/operator.py`` would
#: otherwise shadow the standard library's ``operator`` from the script's own directory).
_BOOT = Path(__file__).with_name("_traced_memo_boot.py")

#: The calls a generating process is answering, outermost first: a memo calling itself with one key would
#: otherwise start interpreters without end (qa finding 13).
_ANCESTRY: tuple[tuple[str, str], ...] = ()


@dataclasses.dataclass(frozen=True)
class _Job:
    """One generation, as the generating process receives it."""

    memo: TracedMemo[..., Any]
    args: tuple[Any, ...]
    kwargs: dict[str, Any]
    key: str
    root: str
    ancestry: tuple[tuple[str, str], ...]


def _failure(error: BaseException) -> dict[str, Any]:
    """An exception as it crosses back: pickled, and described in case the parent cannot rebuild it."""
    described = "".join(traceback.format_exception(error))
    try:
        pickled: bytes | None = pickle.dumps(error)
    except Exception:  # noqa: BLE001 - an exception that cannot cross is described instead
        pickled = None
    return {"ok": False, "error": pickled, "described": described}


def _generate_here(job_bytes: bytes, recording: Recording) -> bytes:
    """The generating process's body (called by the boot script, recording since its first line): run the
    call, stop recording, write the entry; the answer is a pickled ``{"ok": ...}``."""
    global _ANCESTRY
    try:
        job: _Job = pickle.loads(job_bytes)
        _TARGETS.append(Store(Path(job.root)))
        _ANCESTRY = job.ancestry
        value = job.memo.__wrapped__(*job.args, **job.kwargs)
    except BaseException as error:  # noqa: BLE001 - every failure crosses back to the caller
        recording.stop()
        return pickle.dumps(_failure(error))
    recording.stop()
    store = Store(Path(job.root))
    try:
        store.write(job.memo.function_id, job.key, _manifest(recording, store), value, job.memo.returned_types)
    except BaseException as error:  # noqa: BLE001 - an unpinnable run or an unencodable value crosses back
        return pickle.dumps(_failure(error))
    return pickle.dumps({"ok": True})


P = ParamSpec("P")
R = TypeVar("R")


class TracedMemo(Generic[P, R]):
    """A pure function of content-identified arguments, memoised on disk and generated in a fresh process."""

    def __init__(self, function: Callable[P, R], canonical: Mapping[str, Callable[[Any], Any]] | None = None) -> None:
        functools.update_wrapper(self, function)
        self.__wrapped__ = function
        self.function_id = f"{function.__module__}:{function.__qualname__}"
        self._signature = inspect.signature(function)
        self._canonical = dict(canonical or {})
        if unknown := sorted(set(self._canonical) - set(self._signature.parameters)):
            raise TypeError(f"traced_memo({self.function_id}): canonical names no parameter of it: {unknown}")

    def __reduce__(self) -> Any:
        return (_resolve, (self.function_id,))

    def __get__(self, instance: object, owner: type | None = None) -> Any:
        """A memoised method binds like a function: ``instance.method(x)`` is the memo called with
        ``(instance, x)``, so the instance's content is part of the key (it must have content identity)."""
        return self if instance is None else types.MethodType(self, instance)

    @functools.cached_property
    def returned_types(self) -> set[str]:
        """The dataclasses the return annotation names, transitively through their fields: the only types a
        payload of this function may name (a payload never imports a type its function does not declare)."""
        found: set[str] = set()

        def visit(annotation: Any) -> None:
            for member in typing.get_args(annotation) or (annotation,):
                if isinstance(member, type) and dataclasses.is_dataclass(member) and _dotted(member) not in found:
                    found.add(_dotted(member))
                    for hint in typing.get_type_hints(member).values():
                        visit(hint)

        visit(inspect.get_annotations(self.__wrapped__, eval_str=True).get("return"))
        return found

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
        return _digest(encode_exact((self.function_id, dict(bound.arguments), platform_tag())))

    def key(self, *args: P.args, **kwargs: P.kwargs) -> str:
        return self._key(self._arguments(args, kwargs))

    def lookup(self, *args: P.args, **kwargs: P.kwargs) -> Verdict:
        """The entry's verdict for this call, without generating anything."""
        return Store.current().lookup(self.function_id, self.key(*args, **kwargs))

    def __call__(self, *args: P.args, **kwargs: P.kwargs) -> R:
        store = _target()
        if isinstance(store, InProcess):
            return self.__wrapped__(*args, **kwargs)
        bound = self._arguments(args, kwargs)
        key = self._key(bound)
        if (self.function_id, key) in _ANCESTRY:
            raise RecursionError(f"traced_memo: {self.function_id} calls itself with the same arguments")
        if len(_ANCESTRY) >= _MAX_GENERATION_DEPTH:
            raise RecursionError(f"traced_memo: {self.function_id} would start a generation {_MAX_GENERATION_DEPTH + 1} deep")
        verdict = store.lookup(self.function_id, key)
        if not isinstance(verdict, Hit):
            self._generate(_Job(self, bound.args, bound.kwargs, key, str(store.root), (*_ANCESTRY, (self.function_id, key))))
            verdict = store.lookup(self.function_id, key)
            if not isinstance(verdict, Hit):
                raise RuntimeError(f"traced_memo: the entry just generated for {self.function_id} does not validate: {verdict}")
        note_child((self.function_id, key, verdict.payload_digest))
        return decode_payload(verdict.tree, verdict.arrays, self.returned_types)

    def _generate(self, job: _Job) -> None:
        """Run ``job`` in a fresh interpreter with this process's ``-O`` level, ``sys.path`` and working
        directory, and the declared environment only (:data:`_CHILD_ENVIRONMENT`: a withdrawn generator's
        ``ORPHEUS_*`` switch never reaches it, and no variable selects an answer no pin holds). The answer comes back on the
        process's standard output, which the boot script reserves for it (qa finding 9)."""
        optimize = ["-" + "O" * sys.flags.optimize] if sys.flags.optimize else []
        with spawning():
            run = subprocess.run(
                [sys.executable, *optimize, "-P", str(_BOOT)],
                input=pickle.dumps({"sys_path": list(sys.path), "job": _constructor_pickle(job)}),
                capture_output=True, cwd=os.getcwd(),
                env=_child_environment(),
            )
        try:
            result = pickle.loads(run.stdout)
        except Exception as error:  # noqa: BLE001 - the process died before it could answer
            raise RuntimeError(f"traced_memo: the process generating {self.function_id} failed:\n{run.stderr.decode()[-4000:]}") from error
        if result["ok"]:
            return
        try:
            error = pickle.loads(result["error"]) if result["error"] is not None else None
        except Exception:  # noqa: BLE001 - pickled there, not rebuildable here (qa finding 11)
            error = None
        if not isinstance(error, BaseException):
            raise RuntimeError(f"traced_memo: the process generating {self.function_id} raised:\n{result['described']}")
        error.add_note(f"raised in the process generating {self.function_id}")
        raise error


@overload
def traced_memo(function: Callable[P, R], /) -> TracedMemo[P, R]: ...
@overload
def traced_memo(*, canonical: Mapping[str, Callable[[Any], Any]]) -> Callable[[Callable[P, R]], TracedMemo[P, R]]: ...
def traced_memo(function: Callable[P, R] | None = None, /, *, canonical: Mapping[str, Callable[[Any], Any]] | None = None) -> Any:
    """Memoise ``function`` on disk (:class:`TracedMemo`); ``canonical`` maps a parameter to the function that
    puts its argument in canonical form (``{"radii": as_float_array}``), for the key and the generation alike."""
    if function is None:
        return lambda f: TracedMemo(f, canonical)
    return TracedMemo(function, canonical)
