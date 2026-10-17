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
* does not depend on a set's iteration order, which follows the hash seed: the key orders a set by its
  elements, and the generation iterates under ``PYTHONHASHSEED=0``, not the caller's seed.

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
import atexit
import contextlib
import copy
import csv
import dataclasses
import functools
import importlib
import importlib.metadata
import inspect
import io
import json
import os
import pickle
import platform
import site
import stat
import subprocess
import sys
import tempfile
import threading
import sysconfig
import traceback
import types
import typing
import uuid
import warnings
import zipfile
from collections.abc import Callable, Generator, Iterable, Mapping
from pathlib import Path
from typing import Any, Generic, NamedTuple, ParamSpec, TypeVar, overload

import numpy as np

from orpheus.numerics import _traced_memo_boot
from orpheus.numerics._traced_memo_boot import (
    CACHE, CONTENT, FORK_MEASURED, LISTING, Kept, Witness, PRESENCE, Observation, Recording, _read_frame, _write_frame, digest as _digest, file_digest as _file_digest, listing_digest as _listing_digest, presence as _presence, frame, note_child,
    spawning, unrecorded, write_atomically,
)
from orpheus.numerics.content import ContentIdentity, _rebuild, constructor_arguments, encode_exact

__all__ = [
    "Absent", "ChildPin", "Corrupt", "DataPin", "DefPin", "DistributionPin", "EnvironmentPin", "Hit", "InterpreterPin", "PresencePin",
    "ListingPin", "Manifest", "MemoPin", "ModulePin", "Pin", "Stale", "TracedMemo", "Unencodable",
    "TemplateRefused", "Unpinnable", "Verdict", "WorkingDirectoryPin", "bypass", "cache_root", "decode_payload", "default_root",
    "encode_payload", "function_digest", "origin", "platform_tag", "skeleton_digest", "trace_call",
    "traced_memo", "validate",
]


class Unpinnable(RuntimeError):
    """A run depends on something no manifest can pin: a site-packages file of no distribution, code compiled
    under a file's name at a line where the file holds no such code, or a process the run started."""


class TemplateRefused(RuntimeWarning):
    """A fresh fork template refused its first job, so this process generates every later miss in a fresh
    interpreter: still correct, slower, and the forked route goes unexercised. Warned once, when it happens,
    with the template's reasons. While it was silent, a full run could lose the forked-route rows to skips and
    print no cause (``[M]`` 2026-10-16: three extra skips in one run, not reproduced, cause unrecovered)."""


class Unencodable(TypeError):
    """A return value the payload cannot write exactly; it is refused, never written lossily."""


# ── the normalised source ───────────────────────────────────────────────────────────────────────────

_DEF = (ast.FunctionDef, ast.AsyncFunctionDef)
#: The nodes whose fields hold statements: a statement, an ``except`` handler and a ``match`` case. A docstring
#: and a def are statements, so every walk that looks for either descends these and never an expression.
_HOLDS_STATEMENTS = (ast.stmt, ast.excepthandler, ast.match_case)


def _without_docstrings(tree: ast.Module) -> ast.Module:
    """``tree`` with the docstring of the module, of every class and of every def removed.

    A docstring is the first statement of a module, class or def body, and every class and def is a
    statement, so the walk descends statements only (and the handlers and cases that hold statements),
    never an expression: walking every expression node cost a cheap miss about 0.2 s (``[M]`` #602)."""
    def strip(node: ast.AST) -> None:
        if isinstance(node, (ast.Module, ast.ClassDef, *_DEF)):
            first = node.body[0] if node.body else None
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) and isinstance(first.value.value, str):
                # A body that was only a docstring keeps a ``pass`` AT THE DOCSTRING'S LINES, so the first line of
                # a def's body is still a line of the file (``_def_pin`` reads it; digests exclude locations).
                node.body = node.body[1:] or [ast.copy_location(ast.Pass(), first)]
        for child in ast.iter_child_nodes(node):
            if isinstance(child, _HOLDS_STATEMENTS):
                strip(child)

    strip(tree)
    return tree


def _ast_digest(node: ast.AST) -> str:
    """The digest of ``node``'s structure: line and column attributes excluded, so layout never moves it."""
    return _digest(ast.dump(node, include_attributes=False).encode())


def _span(node: ast.AST) -> range:
    first = min([node.lineno, *(d.lineno for d in getattr(node, "decorator_list", []))])  # type: ignore[attr-defined]
    return range(first, node.end_lineno + 1)  # type: ignore[attr-defined]


#: Code objects Python generates from a SIGNATURE or a class body rather than from a def body: lazy
#: annotations (PEP 649, ``__annotate__``) and the scopes of type parameters (PEP 695). Their source is the
#: signature, which the skeleton pins; they never pin the body of the def they belong to (reading a def's
#: annotations does not run it).
_SIGNATURE_SCOPES = ("__annotate__", "<generic parameters of ", "<TypeVar bound of ", "<TypeVar constraint of ",
                     "<TypeAlias ", "<type parameters of ")
#: The anonymous expressions, each with the name of the code object it compiles to.
_ANONYMOUS_NAME: Mapping[type[ast.AST], str] = {
    ast.Lambda: "<lambda>", ast.GeneratorExp: "<genexpr>", ast.ListComp: "<listcomp>",
    ast.SetComp: "<setcomp>", ast.DictComp: "<dictcomp>",
}


class _IndexedDef(NamedTuple):
    """An outermost def: its qualified name, its rank among the defs of that name (0 first), its span with
    decorators (first and last line, inclusive), the line of its stripped body's first statement, and the
    digest of its normalised source."""

    qualname: str
    ordinal: int
    first: int
    last: int
    body_line: int
    digest: str


class _Scope(NamedTuple):
    """A node a code object can be compiled from (a def, a class, a lambda, a comprehension): the name of that
    code object and the node's span (inclusive)."""

    name: str
    first: int
    last: int


@dataclasses.dataclass(frozen=True)
class _SourceIndex:
    """What a manifest and a validation read from one source text: the skeleton's digest, every outermost def,
    and every node a code object can be compiled from. It holds no tree, so it is written to disk and read
    back (:func:`_source_index`); :func:`_parsed_index` builds it from the text."""

    skeleton: str
    defs: tuple[_IndexedDef, ...]
    scopes: tuple[_Scope, ...]

    def digest(self, qualname: str, ordinal: int = -1) -> str | None:
        """The digest of the ``ordinal``-th def named ``qualname`` (``-1``: the last, the one Python binds)."""
        named = [d for d in self.defs if d.qualname == qualname]
        return named[ordinal].digest if -len(named) <= ordinal < len(named) else None

    def enclosing(self, line: int) -> _IndexedDef | None:
        """The outermost def whose span holds ``line`` (outermost defs never overlap)."""
        return next((d for d in self.defs if d.first <= line <= d.last), None)

    def holds(self, within: _IndexedDef | None, line: int, name: str) -> bool:
        """Whether the source of a code object named ``name`` first reported at ``line`` lies in ``within``'s
        span (``None``: anywhere in the module). A node lies in an outermost def's span exactly when the def's
        tree holds it, since outermost defs never overlap (``[M]`` 2026-10-15, the elegance review: 24 262 of
        24 262 nodes over the tracked files)."""
        if name.startswith(_SIGNATURE_SCOPES):
            return True
        return any(s.name == name and s.first <= line <= s.last
                   and (within is None or within.first <= s.first and s.last <= within.last) for s in self.scopes)

    def to_json(self) -> bytes:
        return json.dumps({"skeleton": self.skeleton, "defs": self.defs, "scopes": self.scopes}).encode()

    @classmethod
    def from_json(cls, body: bytes) -> _SourceIndex:
        data = json.loads(body)
        return cls(data["skeleton"], tuple(_IndexedDef(*d) for d in data["defs"]), tuple(_Scope(*r) for r in data["scopes"]))


@functools.cache
def _parsed_index(data: bytes) -> _SourceIndex:
    """The index of a source text, built from its parse and cached in memory by its BYTES (an edited file is
    another key, never a stale read).

    ONE walk finds every outermost def: it descends every statement except a def's body, so a def under any
    block (``if``, ``for``, ``while``, ``match``, ``try``, ``with``) is found; a class adds its name to the
    qualified names of the defs in its body. The same walk returns the skeleton: a shallow copy of the
    statement spine with every outermost def's body read as ``pass`` (a def nested in a body vanishes with that
    body), so the parsed tree is never changed. Every def is digested: the index is built once per content
    and then read back from disk."""
    tree = _without_docstrings(ast.parse(data))
    defs: list[_IndexedDef] = []
    seen: dict[str, int] = {}

    def hollowed(node: ast.AST, prefix: str) -> ast.AST:
        if isinstance(node, _DEF):
            qualname = prefix + node.name
            seen[qualname] = seen.get(qualname, -1) + 1
            span = _span(node)
            defs.append(_IndexedDef(qualname, seen[qualname], span.start, span.stop - 1, node.body[0].lineno, _ast_digest(node)))
            hollow = copy.copy(node)
            hollow.body = [ast.Pass()]
            return hollow
        inner = f"{prefix}{node.name}." if isinstance(node, ast.ClassDef) else prefix
        clone = copy.copy(node)
        for field, value in ast.iter_fields(node):
            if isinstance(value, list) and any(isinstance(v, _HOLDS_STATEMENTS) for v in value):
                setattr(clone, field, [hollowed(child, inner) for child in value])
        return clone

    skeleton = _ast_digest(hollowed(tree, ""))
    scopes = tuple(
        _Scope(node.name if isinstance(node, (*_DEF, ast.ClassDef)) else _ANONYMOUS_NAME[type(node)],
               _span(node).start, _span(node).stop - 1)
        for node in ast.walk(tree) if isinstance(node, (*_DEF, ast.ClassDef, *_ANONYMOUS_NAME))
    )
    return _SourceIndex(skeleton, tuple(defs), scopes)


def function_digest(source: bytes | str, qualname: str, ordinal: int = -1) -> str | None:
    """The digest of the ``ordinal``-th outermost def ``qualname`` (``function`` or ``Class.method``; ``-1``,
    the default, is the last, which Python binds), docstrings stripped; ``None`` when there is none."""
    return _parsed_index(source.encode() if isinstance(source, str) else source).digest(qualname, ordinal)


def skeleton_digest(source: bytes | str) -> str:
    """The digest of the module with every function body replaced by ``pass`` and every docstring stripped:
    what a traced function reads by name without its own body changing (constants, imports, decorators,
    signatures with their defaults, class attributes, dataclass fields)."""
    return _parsed_index(source.encode() if isinstance(source, str) else source).skeleton


#: Where the source indexes live, beside the store and the bytecode.
_SOURCE_INDEXES = Path(CACHE) / "sources"


class _CorruptIndex(NamedTuple):
    """An index file that does not hold the index of the key it is named by (a parse error, a body that does
    not match its digest, a header for other bytes) or cannot be read. It is rebuilt, never read."""

    reason: str


def _read_index(path: Path, header: Mapping[str, str]) -> _SourceIndex | _CorruptIndex | None:
    """The index stored at ``path`` if it is the index ``header`` names; ``None`` when there is none."""
    try:
        raw = path.read_bytes()
    except FileNotFoundError:
        return None
    except OSError as error:
        return _CorruptIndex(f"{type(error).__name__}: {error}")
    try:
        head, body = raw.split(b"\n", 1)
        if json.loads(head) != {**header, "body": _digest(body)}:
            return _CorruptIndex("the header names other bytes, or the body does not match its digest")
        return _SourceIndex.from_json(body)
    except (ValueError, TypeError, KeyError) as error:
        return _CorruptIndex(f"{type(error).__name__}: {error}")


@functools.cache
def _source_index(data: bytes) -> _SourceIndex:
    """The index of a source text: from memory, else from disk, else parsed.

    An index on disk is named by the digest of (the text's bytes, the interpreter, the memo's source), so a
    change to the normalisation is a change of name. Only a GENERATING process writes one: its memo was
    compiled from bytes it pinned (:data:`_RUNNING_MEMO_DIGEST`), so the name it writes under is the code that
    built the index. Any other process may run its memo from a timestamp ``.pyc`` or from a text edited after
    it was read, so it names indexes by the memo's source on disk, reads the indexes generations wrote under
    that name (the code that wrote the entries it validates), and parses in memory when there is none, never
    writing (qa's round-two review of #602: a process running old code wrote an index under new code's name,
    and a later process read the old normalisation's digests as its own). ``[M]`` #602: 128 indexes read back
    in 9.6 ms against 225 ms of parsing; building them costs 0.56 s once per content."""
    parts = (_digest(data), _python_identity(), _INDEX_MEMO_DIGEST)
    header: dict[str, str] = dict(zip(("bytes", "python", "memo"), parts, strict=True))
    path = _SOURCE_INDEXES / f"{_digest(chr(0).join(parts).encode())}.json"
    with unrecorded():
        found = _read_index(path, header)
        if isinstance(found, _SourceIndex):
            return found
        index = _parsed_index(data)
        if _RUNNING_MEMO_DIGEST is not None:
            body = index.to_json()
            write_atomically(path, json.dumps({**header, "body": _digest(body)}).encode() + b"\n" + body)
        return index


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


class _Installed:
    """The installed environment as one checkout reads it: which distribution owns each file under the
    site-packages directories, and each distribution's version. Built once per checkout (one manifest or one
    validation), each part on first use.

    Ownership is a property of the FILE, never of its top-level name: the distribution whose ``RECORD`` lists
    it. A namespace top such as ``sphinxcontrib`` is claimed by several distributions, and the top-level map
    that preceded this table pinned ``sphinxcontrib/tikz.py`` to ``sphinxcontrib-applehelp``, so an upgrade of
    ``sphinxcontrib-tikz`` left valid every entry that ran its code (``[M]`` 2026-10-15, the test-architect,
    #602). Reading the raw ``RECORD`` bytes of every distribution builds the table for 30 446 paths in 18 to
    20 ms; the top-level map cost 0.22 s. A ``RECORD`` row is parsed as CSV, since a path may be quoted and
    hold a comma. What is installed as a distribution is kept for the process (:data:`_DISTRIBUTIONS`): a find
    parses metadata across ``sys.path`` (1245 lookups for 6 distributions cost a cheap miss 0.64 s, ``[M]`` #602)."""

    def __init__(self) -> None:
        self._listed: dict[str, dict[str, Path]] | None = None
        self._shared: set[tuple[str, str]] = set()
        self._owners: dict[Path, tuple[str, bool]] = {}

    def _table(self) -> dict[str, dict[str, Path]]:
        if self._listed is None:
            self._listed = {}
            for site_dir in _site_dirs():
                listed = self._listed.setdefault(site_dir, {})
                for info in Path(site_dir).glob("*.dist-info"):
                    record = info / "RECORD"
                    if record.is_file():
                        for row in csv.reader(io.StringIO(record.read_text(encoding="utf-8"))):
                            if row and listed.setdefault(row[0], info) != info:
                                self._shared.add((site_dir, row[0]))
        return self._listed

    def owner(self, site_dir: str, relative: str) -> tuple[str, bool] | None:
        """The name of the distribution whose ``RECORD`` lists ``relative`` (a path under ``site_dir``, in
        ``/`` form) and whether it is an editable install; ``None`` when no ``RECORD`` lists it. A file that
        several ``RECORD`` files list is refused: pinning one owner would leave an upgrade of another unseen
        (qa's round-two review of #602; 0 of 30 446 paths in the environment)."""
        info = self._table().get(site_dir, {}).get(relative)
        if info is None:
            return None
        if (site_dir, relative) in self._shared:
            raise Unpinnable(f"{site_dir}{relative}: listed by the RECORD of several installed distributions, so no one version pins it")
        if info not in self._owners:
            distribution = importlib.metadata.Distribution.at(info)
            direct = json.loads(distribution.read_text("direct_url.json") or "{}")
            self._owners[info] = distribution.metadata["Name"], bool(direct.get("dir_info", {}).get("editable"))
        return self._owners[info]

    def installed(self, name: str) -> tuple[str, str] | None:
        """What is installed as the distribution ``name``: its version and the digest of its ``RECORD`` (the
        hash of every file it installed, so a rebuild reinstalled under the same version is another install,
        qa's review of #602 cut E); ``None`` when none is installed. From the process's :data:`_DISTRIBUTIONS`."""
        return _DISTRIBUTIONS((name, tuple(sys.path), os.getcwd()))


def _find_installed(key: tuple[str, tuple[str, ...], str], witness: Witness) -> tuple[str, str] | None:
    """What is installed as the distribution ``key[0]`` under ``sys.path`` ``key[1]`` and working directory
    ``key[2]``: its version and the digest of its ``RECORD``, or ``None``. The witness stats every ``sys.path``
    directory (an install or an uninstall adds or removes a ``.dist-info`` there) before the find, then the
    found directory's ``METADATA``, ``PKG-INFO`` and ``RECORD`` before they are read, named from the directory
    rather than from what ``RECORD`` lists (the elegance review of #602's hit cost, A6). A find that names no
    directory is never kept."""
    name, entries, directory = key
    for entry in entries:
        witness.stat(entry or directory)
    try:
        distribution = importlib.metadata.distribution(name)
    except importlib.metadata.PackageNotFoundError:
        return None
    info = getattr(distribution, "_path", None)  # a PathDistribution's directory; the stdlib names it nowhere else
    if info is None:
        witness.unkeepable()
    else:
        for file in ("METADATA", "PKG-INFO", "RECORD"):
            witness.stat(os.path.join(info, file))
    record = distribution.read_text("RECORD")
    return distribution.version, _digest(record.encode()) if record is not None else "absent"


#: What validations have found installed, kept for the process (``[M]`` 2026-10-17: 38 142 finds held 47.7 of
#: the 292 profiled seconds of the warm three-file run, one per distribution per validation).
_DISTRIBUTIONS = Kept(_find_installed)


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


def origin(filename: str, installed: _Installed) -> Origin:
    """Where a file a run touched belongs, and so how the manifest pins it. ``installed`` is the ownership
    table of the checkout asking: building one reads every ``RECORD``, so it is built once per checkout and
    passed, never rebuilt per file."""
    if filename.startswith("<"):
        return Dropped(filename)
    real = _real(filename)
    for site_dir in _site_dirs():
        if real.startswith(site_dir):
            relative = real[len(site_dir):]
            top = relative.split(os.sep)[0]
            if top.endswith(".dist-info"):  # a distribution's own metadata, read by importlib.metadata
                return Distribution(importlib.metadata.Distribution.at(Path(site_dir) / top).metadata["Name"])
            owned = installed.owner(site_dir, Path(relative).as_posix())
            if owned is None:
                raise Unpinnable(f"{filename}: in site-packages and listed by no installed distribution's RECORD, so no version pins it")
            name, editable = owned
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


def _python_identity() -> str:
    return f"{sys.version}|{sys.implementation.cache_tag}"


@functools.cache
def _blas() -> str:
    """The BLAS numpy was built against (Accelerate and OpenBLAS differ in the last bits, #504)."""
    config = np.show_config(mode="dicts")
    return config.get("Build Dependencies", {}).get("blas", {}).get("name", "unknown") if isinstance(config, dict) else "unknown"


def platform_tag() -> str:
    """The operating system, the machine, the ``-O`` level (a bare ``assert`` in a generator is stripped under
    it) and the BLAS."""
    return f"{sys.implementation.cache_tag}-{sys.platform}-{platform.machine()}-O{sys.flags.optimize}-{_blas()}"


def fork_platform() -> tuple[str, str, str]:
    """The platform as the fork protocol's table names it (:data:`~orpheus.numerics._traced_memo_boot.FORK_MEASURED`):
    the operating system, the machine and the BLAS."""
    return sys.platform, platform.machine(), _blas()


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
    the fixed ones. The parent projects its own and hands it to the generation on the job (``_Job.environment``),
    so the pin never reads the child's ``os.environ``, which PEP 538 rewrites under ``LANG=C``; ``trace_call``
    projects its caller's."""
    source = os.environ if environment is None else environment
    return {k: source[k] for k in _CHILD_ENVIRONMENT if k in source} | _FIXED_ENVIRONMENT


def _read_source(key: tuple[str, tuple[str, ...], str], witness: Witness) -> bytes | None:
    """The bytes of the source a pin's path ``key[0]`` names under ``sys.path`` ``key[1]`` and working directory
    ``key[2]``: the absolute path itself, or the first entry holding a regular file there (the one the same import
    would load); ``None`` when none does. Every candidate up to the one read is witnessed, a missing one
    included, so a file appearing on an earlier entry is a change and the shadowing file is read."""
    path, entries, directory = key
    candidates = (path,) if os.path.isabs(path) else tuple(os.path.join(entry or directory, path) for entry in entries)
    for candidate in candidates:
        found = witness.stat(candidate)
        if found is not None and stat.S_ISREG(found.st_mode):
            return Path(candidate).read_bytes()
    return None


#: The pinned sources validations have read, kept for the process. Without it every validation re-read every
#: pinned file: ``[M]`` 2026-10-17, the warm three-file run (system, reading, assembly; not slow) made 5346
#: lookups, 832 287 opens and 1.4 M stats over a few hundred files, and the validations held 111 of its 292
#: profiled seconds, so the warm run cost what no memo did (222 s against 221 s; 155 s after).
_SOURCES = Kept(_read_source)
#: A data pin's file and a listing pin's directory, kept for the process.
_FILE_DIGESTS = Kept(_file_digest)
_LISTINGS = Kept(_listing_digest)


class _Checkout:
    """One validation's view of the checkout: one state per source file (every pin of a file is judged against
    one read of it; this does not reach across files or across children, each of which validates with its own
    checkout), the installed environment, and the store children live in."""

    def __init__(self, store: Store) -> None:
        self.store = store
        self._texts: dict[str, bytes | None] = {}
        self.installed = _Installed()

    def text(self, path: str) -> bytes | None:
        if path not in self._texts:
            self._texts[path] = _SOURCES((path, tuple(sys.path), os.getcwd()))
        return self._texts[path]



def _read_memo_digest(paths: tuple[Path, ...], witness: Witness) -> str:
    for path in paths:
        witness.stat(str(path))
    return _digest(b"".join(path.read_bytes() for path in paths))


_MEMO_DIGESTS = Kept(_read_memo_digest)


def _memo_digest() -> str:
    """The digest of the memo's source on disk now (kept for the process while its three files are unchanged)."""
    return _MEMO_DIGESTS(_MEMO_SOURCES)


#: The memo's source as this process imported it: it built every index this process writes, so it keys them
#: (an edit to the memo's source after the import changes the next process's indexes, not this one's).
def _pinned_memo_digest() -> str | None:
    """The digest of the bytes this process compiled the memo from, when it compiled them through the
    pinned-source loader (a generating process); ``None`` in any other process, which may run the memo from a
    timestamp ``.pyc`` or from a text edited after it was read, and so cannot name its own code."""
    pinned = _traced_memo_boot.PINNED_SOURCES
    if pinned is None:
        return None
    return _digest(b"".join(pinned[os.path.realpath(p)] for p in _MEMO_SOURCES))


#: The memo's source this process RAN, named by its pinned bytes in a generating process (``None``
#: elsewhere): the memo pin of every entry a generation writes, and the name of every index it writes.
_RUNNING_MEMO_DIGEST = _pinned_memo_digest()
#: The memo's source a process names indexes by: the code it ran in a generation; elsewhere the source on disk at
#: import, which names the indexes generations wrote from it (no recording is active when a process that is not
#: generating imports the memo, so this read is nobody's dependency).
_INDEX_MEMO_DIGEST = _RUNNING_MEMO_DIGEST if _RUNNING_MEMO_DIGEST is not None else _memo_digest()


class DefPin(NamedTuple):
    """A first-party def that ran: its file, its qualified name, its rank among the defs of that name, and
    the digest of its normalised source."""

    path: str
    qualname: str
    ordinal: int
    digest: str

    def stale(self, checkout: _Checkout) -> str | None:
        text = checkout.text(self.path)
        found = None if text is None else _source_index(text).digest(self.qualname, self.ordinal)
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
        return f"module {self.path} skeleton changed" if _source_index(text).skeleton != self.digest else None


class DistributionPin(NamedTuple):
    """A third-party distribution whose code ran or whose file was read, by its version and its ``RECORD``."""

    name: str
    version: str
    record: str

    def stale(self, checkout: _Checkout) -> str | None:
        now = checkout.installed.installed(self.name)
        if now is None:
            return f"distribution {self.name}: {self.version} is no longer installed"
        if now != (self.version, self.record):
            return f"distribution {self.name}: {self.version} (RECORD {self.record[:12]}) != {now[0]} (RECORD {now[1][:12]})"
        return None


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
        return f"data {self.path} changed" if _FILE_DIGESTS(self.path) != self.digest else None


class ListingPin(NamedTuple):
    """A directory the run listed, by its absolute path and the names in it (``"absent"``: it was missing)."""

    path: str
    digest: str

    def stale(self, checkout: _Checkout) -> str | None:
        return f"listing {self.path} changed" if _LISTINGS(self.path) != self.digest else None


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

    def with_compiled_sources(self, sources: Mapping[str, bytes]) -> Manifest:
        """This manifest, holding also the whole bytes each first-party source was compiled from: a process that
        keeps running that code (a template) is stale when any of them differs on disk, an edit to a def body
        that never ran included, which the skeleton and def pins do not see."""
        return Manifest(self.pins | {DataPin(path, _digest(text)) for path, text in sources.items()})

    def replacing(self, kind: type[Any], pins: Iterable[Pin]) -> Manifest:
        """This manifest with every pin of ``kind`` replaced by ``pins``."""
        return Manifest(frozenset({p for p in self.pins if not isinstance(p, kind)} | set(pins)))

    def to_json(self) -> list[list[Any]]:
        return sorted([type(pin).__name__, *pin] for pin in self.pins)

    @classmethod
    def from_json(cls, rows: list[list[Any]]) -> Manifest:
        return cls(frozenset(_PIN_KINDS[kind](*fields) for kind, *fields in rows))


def _def_pin(index: _SourceIndex, path: str, code: types.CodeType) -> DefPin | None:
    """The pin of the outermost def whose source holds ``code``; ``None`` for module-level code (a module or
    class body, a module-level lambda, a signature), which the skeleton pins."""
    line, name = code.co_firstlineno, code.co_name
    found = index.enclosing(line)
    if name == "<module>" or (found is not None and name.startswith(_SIGNATURE_SCOPES) and line < found.body_line):
        return None
    if index.holds(found, line, name):
        return None if found is None else DefPin(path, found.qualname, found.ordinal, found.digest)
    raise Unpinnable(f"{path}:{line} {code.co_qualname}: no def of that name at its line")


def _unchanged_by_the_run(recording: Recording, kind: Observation, path: str, seen: str) -> str:
    """What the run saw at ``path``, refused when the run itself changed it afterwards: the state it saw no longer
    exists, so no state a later lookup can find pins the answer (qa's third review of #602 cut E, ``[M]``
    2026-10-16: a generator that read ``n`` and wrote ``n + 1`` was served as a Hit, while a cold run returned
    ``n + 1``)."""
    if kind.changed_by(path, recording.written) and kind.observe(path) != seen:
        raise Unpinnable(f"{path}: the run observed it and then changed it, so no state of it pins the answer")
    return seen


def _manifest(recording: Recording, store: Store) -> Manifest:
    """The manifest of a recorded run: every code object it started, file it read, directory it listed and
    memo entry it read, each pinned by its origin."""
    # The memo pin names the code that ran in a generation; trace_call, the in-process instrument, re-reads it.
    running = _RUNNING_MEMO_DIGEST if _RUNNING_MEMO_DIGEST is not None else _memo_digest()
    pins: set[Pin] = {InterpreterPin(_python_identity()), MemoPin(running), EnvironmentPin(_environment_digest(recording.environment))}
    pins |= {_program_pin(program) for program in recording.spawned}
    pins |= {ChildPin(*child) for child in recording.children}
    checkout = _Checkout(store)
    by_file: dict[str, list[types.CodeType]] = {}
    for code in recording.codes:
        match origin(code.co_filename, checkout.installed):
            case Distribution(name):
                pins.add(_distribution_pin(checkout, name))
            case Source(path):
                by_file.setdefault(path, []).append(code)
            case Dropped() | Interpreter():
                pass
    for real, file_codes in by_file.items():
        path = _pin_path(real)
        index = _source_index(_compiled_bytes(recording, real))
        pins.add(ModulePin(path, index.skeleton))
        pins |= {pin for code in file_codes if (pin := _def_pin(index, path, code)) is not None}
    store_root = _real(str(store.root)) + os.sep
    for opened, read in recording.opened.items():
        if _real(opened).startswith(store_root):
            continue  # a store entry is a child
        match origin(opened, checkout.installed):
            case Distribution(name):
                pins.add(_distribution_pin(checkout, name))
            case Source():
                pins.add(DataPin(opened, _unchanged_by_the_run(recording, CONTENT, opened, read)))
            case Dropped() | Interpreter():
                pass
    pins |= {PresencePin(probed, _unchanged_by_the_run(recording, PRESENCE, probed, kind))
             for probed, kind in recording.probed.items() if not _real(probed).startswith(store_root)}
    pins |= {ListingPin(listed, _unchanged_by_the_run(recording, LISTING, listed, seen))
             for listed, seen in recording.listed.items()}
    if recording.relative:
        pins.add(WorkingDirectoryPin(recording.directory))
    return Manifest(frozenset(pins))


def _distribution_pin(checkout: _Checkout, name: str) -> DistributionPin:
    """The pin of a distribution whose code ran or whose file was read, by what is installed."""
    installed = checkout.installed.installed(name)
    if installed is None:
        raise Unpinnable(f"the distribution {name} owns a file the run touched, but no version of it is installed")
    return DistributionPin(name, *installed)


def _compiled_bytes(recording: Recording, real: str) -> bytes:
    """The bytes a source file's code ran from. A generating process compiles every first-party module from
    the bytes its loader read and records them (``_traced_memo_boot``), so they are pinned even when the file
    changed during the run, or within the second a timestamp ``.pyc`` cannot tell apart; code of such a file
    that the loader did not compile (an ``exec`` of a file read by hand, ``runpy``) has no bytes to pin and is
    refused. :func:`trace_call`, the in-process instrument, has no loader and reads the file after the run."""
    if recording.sources is None:
        return Path(real).read_bytes()
    if real not in recording.sources:
        raise Unpinnable(
            f"{real}: code of this file ran without being imported (an exec or runpy of a source file), so no "
            "bytes were recorded for it and a re-read after the run could pin bytes that never ran"
        )
    return recording.sources[real]


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



def _resolve_memo(function_id: str) -> TracedMemo[..., Any]:
    """The memo named ``function_id``: its module attribute, or the memo a decorator stacked above it wraps
    (``functools.wraps`` keeps the chain in ``__wrapped__``; :func:`inspect.unwrap` follows it and raises on a cycle).
    The attribute alone would be the outer wrapper, whose ``__wrapped__`` is the memo again, so the child would call
    the memo instead of its body."""
    memo = inspect.unwrap(_resolve(function_id), stop=lambda f: isinstance(f, TracedMemo))
    if not isinstance(memo, TracedMemo) or memo.function_id != function_id:
        raise TypeError(f"traced_memo: {function_id} names no memo, directly or through __wrapped__")
    return memo

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
    return Path(CACHE) / "references"


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
        written = io.BytesIO()
        np.savez(written, allow_pickle=False, **arrays, **{_ENTRY: np.frombuffer(json.dumps(entry).encode(), np.uint8)})
        write_atomically(self.path(function_id, key), written.getvalue())


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
    #: The environment the asking process hands the generation (:func:`_child_environment`), which its entry
    #: pins. Not the generating process's own ``os.environ``: under ``LANG=C`` the interpreter coerces the
    #: locale at its start (PEP 538) and adds ``LC_CTYPE=C.UTF-8``, so that environment differs from the one
    #: the asking process validates against, and every entry would be born stale.
    environment: dict[str, str]


class _Generated(NamedTuple):
    """A generating process's answer: the entry is written."""


class _Failed(NamedTuple):
    """A generating process's answer: the call, or the writing of its entry, raised. The exception is pickled,
    and described in case the parent cannot rebuild it."""

    error: bytes | None
    described: str


class _Refused(NamedTuple):
    """A template's answer: it will not fork (what its imports depended on changed, or it holds more than one
    OS thread); it exits after answering."""

    reasons: tuple[str, ...]


class _Died(NamedTuple):
    """A template's answer: the child it forked ended without answering, with this wait status."""

    status: int


def _failure(error: BaseException) -> _Failed:
    described = "".join(traceback.format_exception(error))
    try:
        pickled: bytes | None = pickle.dumps(error)
    except Exception:  # noqa: BLE001 - an exception that cannot cross is described instead
        pickled = None
    return _Failed(pickled, described)


def _generate_here(job_bytes: bytes, recording: Recording) -> bytes:
    """The generating process's body (called by the boot script, recording since its first line): run the
    call, stop recording, write the entry; the answer is a pickled :class:`_Generated` or :class:`_Failed`."""
    global _ANCESTRY
    try:
        job: _Job = pickle.loads(job_bytes)
        _TARGETS.append(Store(Path(job.root)))
        _ANCESTRY = job.ancestry
        recording.environment = job.environment
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
    return pickle.dumps(_Generated())


#: How many children one template forks before it is replaced: its resident size grew about 16 kB per fork in
#: the protocol (``[M]`` 2026-10-16, 99 to 156 MB over 3600 forks).
_FORKS_PER_TEMPLATE = 256
def _forks_here() -> bool:
    """Whether this process's misses are forked from a template: on a measured platform, never inside a
    generating process (a nested miss spawns its own fresh interpreter, as it always has), and not after a
    fresh template refused to fork."""
    return _traced_memo_boot.PINNED_SOURCES is None and _SLOT.refused is None and fork_platform() in FORK_MEASURED


class _Launch(NamedTuple):
    """How a generating process is started, by either route: the interpreter, the ``-O`` level, ``sys.path``,
    the working directory and the environment the job is handed. A template serves a job iff it was started
    with an equal launch, since the generation's manifest pins that state as it was at the template's start."""

    executable: str
    optimize: int
    sys_path: tuple[str, ...]
    directory: str
    environment: tuple[tuple[str, str], ...]

    @classmethod
    def of(cls, job: _Job) -> _Launch:
        return cls(sys.executable, sys.flags.optimize, tuple(sys.path), os.getcwd(), tuple(sorted(job.environment.items())))

    def argv(self, *mode: str) -> list[str]:
        return [self.executable, *(["-" + "O" * self.optimize] if self.optimize else []), "-P", str(_BOOT), *mode]

    def envelope(self) -> bytes:
        return frame(pickle.dumps({"sys_path": list(self.sys_path)}))

    def options(self) -> dict[str, Any]:
        return {"cwd": self.directory, "env": dict(self.environment)}


def _answer(data: bytes | None, errors: str) -> _Generated | _Failed | _Refused | _Died:
    """A frame from a generating process or a template, read back; a missing or unreadable one raises with the
    tail of the process's standard error."""
    try:
        if data is not None:
            return pickle.loads(data)
    except Exception as error:  # noqa: BLE001 - the process ended inside its answer
        raise RuntimeError(f"traced_memo: a generating process answered unreadably:\n{errors}") from error
    raise RuntimeError(f"traced_memo: a generating process ended without answering:\n{errors}")


class _Template:
    """A template process (``_traced_memo_boot.template``) and the launch it was started with."""

    def __init__(self, launch: _Launch) -> None:
        self.launch = launch
        self.forks = 0
        self.errors = tempfile.TemporaryFile()
        with spawning():
            self.process = subprocess.Popen(launch.argv("--template"), stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                            stderr=self.errors, **launch.options())
        assert self.process.stdin is not None and self.process.stdout is not None
        os.write(self.process.stdin.fileno(), launch.envelope())

    def generate(self, job: bytes) -> _Generated | _Failed | _Refused | _Died | None:
        """The template's answer to ``job``; ``None`` when it has exited (it ended at a previous refusal, or was
        killed)."""
        assert self.process.stdin is not None and self.process.stdout is not None
        try:
            _write_frame(self.process.stdin.fileno(), job)
        except BrokenPipeError:
            return None
        data = _read_frame(self.process.stdout)
        if data is None:
            return None
        self.forks += 1
        return _answer(data, self.stderr_tail())

    def stderr_tail(self) -> str:
        self.errors.seek(0)
        return self.errors.read().decode(errors="replace")[-4000:]

    def release(self) -> None:
        """Close this process's ends of the pipes and the error file, without waiting for the template."""
        for pipe in (self.process.stdin, self.process.stdout):
            if pipe is not None:
                pipe.close()
        self.errors.close()

    def close(self) -> None:
        self.release()
        try:
            self.process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            self.process.kill()
            self.process.wait()


class _TemplateSlot:
    """This process's template, if it has started one; the lock that gives one miss at a time its pipes (two
    threads missing at once interleaved their frames and read each other's answers: qa's review of #602 cut E,
    ``[M]`` 14 of 40 calls raised); and why a fresh template refused its first job, if one has. The refusal's
    cause holds for this process (an OS thread its imports started, say) and for a child forked from it, so
    their misses spawn from then on instead of starting a template per miss."""

    def __init__(self) -> None:
        self.template: _Template | None = None
        self.lock = threading.Lock()
        self.refused: tuple[str, ...] | None = None

    def forget(self) -> None:
        """In a child of this process (``os.fork``): the template is the parent's, so the child closes its copies
        of the pipes (the template reads its end of input only when every holder of its stdin has closed it) and
        starts its own; the lock may have been held by a thread the fork did not copy."""
        if self.template is not None:
            self.template.release()
        self.template, self.lock = None, threading.Lock()

    def close(self) -> None:
        if self.template is not None:
            self.template.close()
            self.template = None

    def answer(self, launch: _Launch, job: bytes) -> _Generated | _Failed | None:
        """The answer of a child forked from the template, starting or replacing the template as the launch, its
        fork count or its refusal requires; ``None`` when no template forked the job (the miss then spawns)."""
        with self.lock:
            for _ in range(2):
                if self.template is not None and (self.template.launch != launch or self.template.forks >= _FORKS_PER_TEMPLATE):
                    self.close()
                if self.template is None:
                    self.template = _Template(launch)
                fresh = self.template.forks == 0
                match self.template.generate(job):
                    case _Generated() | _Failed() as answer:
                        return answer
                    case _Died(status):
                        raise RuntimeError(f"traced_memo: a forked generating process died (wait status {status}):\n{self.template.stderr_tail()}")
                    case _Refused(reasons) if fresh:
                        self.refused = reasons
                        self.close()
                        warnings.warn(f"traced_memo: a fresh fork template refused ({'; '.join(reasons)}); this process "
                                      "now generates every miss in a fresh interpreter", TemplateRefused, stacklevel=2)
                        return None
                    case _Refused() | None:
                        self.close()
            return None


_SLOT = _TemplateSlot()
os.register_at_fork(after_in_child=_SLOT.forget)
atexit.register(_SLOT.close)

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
        return (_resolve_memo, (self.function_id,))

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
            self._generate(_Job(self, bound.args, bound.kwargs, key, str(store.root), (*_ANCESTRY, (self.function_id, key)),
                                _child_environment()))
            verdict = store.lookup(self.function_id, key)
            if not isinstance(verdict, Hit):
                raise RuntimeError(f"traced_memo: the entry just generated for {self.function_id} does not validate: {verdict}")
        note_child((self.function_id, key, verdict.payload_digest))
        return decode_payload(verdict.tree, verdict.arrays, self.returned_types)

    def _generate(self, job: _Job) -> None:
        """Run ``job`` in a generating process: a child forked from this process's template where forking is
        measured (:func:`_forks_here`), else a fresh interpreter. Either is started by the job's
        :class:`_Launch`: this process's ``-O`` level, ``sys.path`` and working directory, and the declared
        environment only (:data:`_CHILD_ENVIRONMENT`: a withdrawn generator's ``ORPHEUS_*`` switch never reaches
        it, and no variable selects an answer no pin holds). The answer comes back as one frame on the
        process's standard output, which the boot script reserves for it (qa finding 9)."""
        job_bytes, launch = _constructor_pickle(job), _Launch.of(job)
        answer = _SLOT.answer(launch, job_bytes) if _forks_here() else None
        if answer is None:
            with spawning():
                run = subprocess.run(launch.argv(), input=launch.envelope() + frame(job_bytes), capture_output=True,
                                     **launch.options())
            answer = _answer(_read_frame(io.BytesIO(run.stdout)), run.stderr.decode(errors="replace")[-4000:])
        match answer:
            case _Generated():
                return
            case _Failed(pickled, described):
                try:
                    error = pickle.loads(pickled) if pickled is not None else None
                except Exception:  # noqa: BLE001 - pickled there, not rebuildable here (qa finding 11)
                    error = None
                if not isinstance(error, BaseException):
                    raise RuntimeError(f"traced_memo: the process generating {self.function_id} raised:\n{described}")
                error.add_note(f"raised in the process generating {self.function_id}")
                raise error
            case _:
                raise RuntimeError(f"traced_memo: the process generating {self.function_id} answered {answer!r}")


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
