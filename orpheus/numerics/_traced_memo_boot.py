"""The traced memo's recorder, and the first code its generating process runs (#405 P3).

This module imports only the standard library, so the generating process can start recording before any
``orpheus`` code runs: it is run by path with ``-P`` (``python -P _traced_memo_boot.py``), and its
``main`` starts a :class:`Recording` on its first line, then loads the job and
:mod:`orpheus.numerics.traced_memo`. Importing it runs nothing (a roster walk over the package imports it),
and :func:`orpheus.numerics.traced_memo.trace_call` records through the same :class:`Recording`: one
definition of what a generation depends on.

A recording holds, for the time it runs, across every thread of the process:

* every code object that started (``sys.monitoring`` ``PY_START``; each reports once, then disables itself);
* every path opened for reading and every path whose presence was probed (``os.stat``, ``os.lstat``, and so
  ``os.path.exists``, ``isfile``, ``isdir`` and ``pathlib``'s ``exists``), each made absolute at the moment
  it was touched, a bytes path decoded, a missing path included (the run then depended on its absence);
* every directory listed (``os.listdir``, ``os.scandir``);
* every process the run started, by its program or, when none can be named, by its event;
* the working directory the run started in, and whether it read a relative path;
* the memo entries the run read, which :mod:`~orpheus.numerics.traced_memo` adds as it serves them.

The import system's own reads, probes and listings are not recorded: they are how it finds a module, and the
module it finds is pinned where its code runs. Recording them would make every new file in any ``sys.path``
directory invalidate every entry. The import system is the frozen bootstrap and the ``importlib`` package
(``importlib.metadata`` lists every ``sys.path`` directory to find the installed distributions, which are
pinned by version). The exemption is decided by the IMMEDIATE caller (the import-system frame that made the
call), never by an import being somewhere on the stack: a module body that lists a
directory while it is imported is the run's own dependency (qa's second review, 2026-10-04).

Audit hooks cannot be removed, so one hook is installed per process and serves the recordings active at the
moment of each event; the ``os.stat`` probe is installed while a recording is active and removed after.
"""
import hashlib
import importlib.machinery
import marshal
import os
import sys
import threading

#: The recordings active in this process, innermost last; every event is offered to each of them.
_ACTIVE: list["Recording"] = []
_LOCK = threading.Lock()
_HOOKED = False
#: Set while the memo itself starts a generating process, so that spawn is not the run's own.
_SPAWNING = threading.local()
#: Set while the memo itself validates or loads an entry: its reads of the files a CHILD entry pinned are not
#: the run's dependencies (the child pin covers them), and recording them would pin every source file of the
#: child by its bytes, so a comment edit would make the parent stale (the archivist, 2026-10-04: 105 of 105
#: of a solve's modules pinned as data on the reading that read it).
_QUIET = threading.local()
#: ``os.stat`` and ``os.lstat`` as they were before the first recording replaced them.
_STAT = {"stat": os.stat, "lstat": os.lstat}

_READ_EVENTS = {"open"}
_LIST_EVENTS = {"os.listdir", "os.scandir"}
_SPAWN_EVENTS = {"subprocess.Popen", "_posixsubprocess.fork_exec", "os.fork", "os.forkpty", "os.posix_spawn",
                 "os.spawn", "os.exec", "os.system"}


def _path(raw: object) -> tuple[str, bool] | None:
    """A path made absolute now, and whether it was relative; ``None`` for a file descriptor."""
    if isinstance(raw, int):
        return None
    try:
        text = os.fsdecode(raw if raw is not None else ".")  # type: ignore[arg-type]
    except TypeError:
        return None
    return os.path.abspath(text), not os.path.isabs(text)


#: The import system's own code: its frozen bootstrap and the ``importlib`` package (``importlib.metadata``
#: lists every ``sys.path`` directory to find the installed distributions, which are pinned by version).
#: Both spellings of its directory, computed once at import: the check runs inside the ``os.stat`` probe, so it
#: must not itself call ``os.stat`` (``realpath`` does).
_IMPORTLIB = tuple({os.path.dirname(f) + os.sep for f in (__import__("importlib").__file__,
                                                            os.path.realpath(__import__("importlib").__file__))})


def _by_importer(depth: int) -> bool:
    """Whether the Python frame ``depth`` levels up, the one that made the call, belongs to the import system."""
    try:
        frame = sys._getframe(depth)
    except ValueError:
        return False
    filename = frame.f_code.co_filename
    return filename.startswith(("<frozen importlib", *_IMPORTLIB))


def _offer(field: str, raw: object) -> None:
    found = _path(raw)
    if found is None:
        return
    path, relative = found
    for recording in list(_ACTIVE):
        getattr(recording, field).add(path)
        recording.relative |= relative


def _program(event: str, args: tuple[object, ...]) -> str:
    """The absolute path of the program a spawn event runs, or the event itself when none can be named."""
    import shutil

    name: object = None
    if event == "subprocess.Popen":
        executable, argv = (args + (None, None))[:2]
        name = executable if executable is not None else (argv[0] if isinstance(argv, (list, tuple)) and argv else argv)
    elif event == "_posixsubprocess.fork_exec" and args and isinstance(args[0], (list, tuple)):
        # (the executable candidates, the argument list, ...): argv[0] is only a label; the program is the
        # first candidate that exists.
        found = next((os.fsdecode(c) for c in args[0] if os.path.isfile(os.fsdecode(c))), None)
        return os.path.realpath(found) if found is not None else event
    if isinstance(name, (str, bytes, os.PathLike)):
        found = shutil.which(os.fsdecode(name))
        if found is not None:
            return os.path.realpath(found)
    return event


def _audit(event: str, args: tuple[object, ...]) -> None:
    if not _ACTIVE or getattr(_QUIET, "on", False):
        return
    if event in _READ_EVENTS:
        mode = args[1] if len(args) > 1 and isinstance(args[1], str) else "r"
        if not any(c in mode for c in "wax+") and not _by_importer(2):
            _offer("opened", args[0] if args else None)
    elif event in _LIST_EVENTS and not _by_importer(2):
        _offer("listed", args[0] if args else None)
    elif event in _SPAWN_EVENTS and not getattr(_SPAWNING, "on", False):
        program = _program(event, args)
        for recording in list(_ACTIVE):
            recording.spawned.add(program)


def _probing(name: str):
    original = _STAT[name]

    def probe(path, *args, **kwargs):
        if _ACTIVE and not getattr(_QUIET, "on", False) and not _by_importer(2):
            _offer("probed", path)
        return original(path, *args, **kwargs)

    return probe


class Recording:
    """What one run depends on, recorded while it runs (module docstring)."""

    def __init__(self) -> None:
        self.codes: set = set()
        self.opened: set[str] = set()
        self.probed: set[str] = set()
        self.listed: set[str] = set()
        self.spawned: set[str] = set()
        self.relative = False
        self.directory = os.getcwd()
        self.environment: dict[str, str] = dict(os.environ)
        self.children: set[tuple[str, str, str]] = set()
        #: The bytes each first-party module was compiled from, by real path, when the pinned-source loader
        #: is installed (a generating process); ``None`` in a process without it (``trace_call``).
        self.sources: dict[str, bytes] | None = None
        self._tool: int | None = None

    def start(self) -> "Recording":
        global _HOOKED
        monitoring = sys.monitoring
        with _LOCK:
            if not _HOOKED:
                sys.addaudithook(_audit)
                _HOOKED = True
            free = [t for t in (monitoring.PROFILER_ID, 3, 4, monitoring.OPTIMIZER_ID) if monitoring.get_tool(t) is None]
            if not free:
                raise RuntimeError("traced_memo: no free sys.monitoring tool id")
            self._tool = free[0]
            monitoring.use_tool_id(self._tool, "traced_memo")
            codes = self.codes

            def started(code, _offset):
                codes.add(code)
                return monitoring.DISABLE

            monitoring.register_callback(self._tool, monitoring.events.PY_START, started)
            monitoring.set_events(self._tool, monitoring.events.PY_START)
            if not _ACTIVE:
                os.stat, os.lstat = _probing("stat"), _probing("lstat")
            self.directory = os.getcwd()
            self.environment = dict(os.environ)  # what the run RECEIVED; a run may set a variable later
            _ACTIVE.append(self)
        return self

    def stop(self) -> None:
        monitoring = sys.monitoring
        with _LOCK:
            if self._tool is None:
                return
            monitoring.set_events(self._tool, 0)
            monitoring.register_callback(self._tool, monitoring.events.PY_START, None)
            monitoring.free_tool_id(self._tool)
            monitoring.restart_events()
            self._tool = None
            _ACTIVE.remove(self)
            if not _ACTIVE:
                os.stat, os.lstat = _STAT["stat"], _STAT["lstat"]


def note_child(child: tuple[str, str, str]) -> None:
    """A memo entry a run read: every active recording pins it."""
    for recording in list(_ACTIVE):
        recording.children.add(child)


class unrecorded:
    """The memo's own reading of the store (a lookup and its validation), which no run depends on."""

    def __enter__(self) -> None:
        self._outer = getattr(_QUIET, "on", False)
        _QUIET.on = True

    def __exit__(self, *_exc: object) -> None:
        _QUIET.on = self._outer


class spawning:
    """The memo's own start of a generating process, which is not a dependency of the run that asked."""

    def __enter__(self) -> None:
        _SPAWNING.on = True

    def __exit__(self, *_exc: object) -> None:
        _SPAWNING.on = False


#: The repository's cache directory (gitignored), the one root of the memo's three stores: the entries
#: (``references``), the source indexes (``sources``) and the pinned bytecode (``bytecode``).
CACHE = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))), ".cache")
#: Where the pinned-source loader keeps bytecode, one directory per interpreter and ``-O`` level. An entry is
#: named by the digest of (path, interpreter, ``-O`` level, source bytes), so sharing it across stores and
#: processes is sound.
_BYTECODE = os.path.join(CACHE, "bytecode", f"{sys.implementation.cache_tag}-O{sys.flags.optimize}")
#: The bytes each first-party module of this process was compiled from, by real path, once the pinned-source
#: loader is installed (a generating process); ``None`` in any other process.
PINNED_SOURCES: dict[str, bytes] | None = None


def write_atomically(path, data: bytes) -> None:
    """Write ``data`` to a temporary file beside ``path``, then move it into place with one ``os.replace``: a
    reader sees the old file or the new one, never a part of either. The memo's three stores write this way."""
    directory = os.path.dirname(os.fspath(path))
    os.makedirs(directory, exist_ok=True)
    partial = f"{os.fspath(path)}.{os.getpid()}.{threading.get_ident()}.part"
    with open(partial, "wb") as out:
        out.write(data)
    os.replace(partial, path)


def _bytecode_key(path: str, source: bytes) -> bytes:
    head = f"{path}\0{sys.implementation.cache_tag}\0{sys.flags.optimize}\0".encode()
    return hashlib.blake2b(head + source, digest_size=32).digest()


class _PinnedSourceLoader(importlib.machinery.SourceFileLoader):
    """A first-party module compiled from the bytes this process read, never from a timestamp ``.pyc``.

    CPython accepts a timestamp ``.pyc`` when the source's whole-second mtime and its size match, so a
    same-size edit within the second of the previous write ran the OLD code while the manifest pinned the NEW
    bytes: a regenerated entry served the old answer and validated (the test-architect's probe, 2026-10-15,
    #602). This loader compiles from the bytes it reads, offers those bytes to every active recording (the
    manifest pins exactly them), and keeps its bytecode apart from ``__pycache__``: an entry carries its key
    and the digest of its code, and one that does not match both is recompiled, never run. Its own reads are
    the memo's, not the run's (``unrecorded``)."""

    def get_code(self, fullname: str):  # type: ignore[override]
        with unrecorded():
            path = self.get_filename(fullname)
            source = self.get_data(path)
            if PINNED_SOURCES is not None:
                PINNED_SOURCES[os.path.realpath(path)] = source
            key = _bytecode_key(path, source)
            entry = os.path.join(_BYTECODE, key.hex())
            try:
                with open(entry, "rb") as stored:
                    data = stored.read()
                if data[:32] == key and data[32:64] == hashlib.blake2b(data[64:], digest_size=32).digest():
                    return marshal.loads(data[64:])
            except (OSError, EOFError, ValueError, TypeError):
                pass
            code = self.source_to_code(source, path)
            payload = marshal.dumps(code)
            try:
                write_atomically(entry, key + hashlib.blake2b(payload, digest_size=32).digest() + payload)
            except OSError:
                pass
            return code


def _vendor_dirs() -> tuple[str, ...]:
    """The standard library and the site-packages directories: code there is pinned by the interpreter or by
    a distribution's version, never by its bytes."""
    import site
    import sysconfig

    paths = sysconfig.get_paths()
    dirs = {paths["stdlib"], paths["platstdlib"], paths["purelib"], paths["platlib"], *site.getsitepackages()}
    return tuple(sorted({os.path.realpath(d) + os.sep for d in dirs}, key=len, reverse=True))


class _PinnedSourceFinder:
    """The first entry of ``sys.meta_path`` in a generating process: it finds a module through every finder
    after it, in order (so a module only an editable install's finder knows is covered), and gives a
    first-party source module the pinned-source loader."""

    vendor: tuple[str, ...] = ()

    @classmethod
    def find_spec(cls, fullname, path=None, target=None):
        for finder in sys.meta_path:
            if finder is cls or not hasattr(finder, "find_spec"):
                continue
            spec = finder.find_spec(fullname, path, target)
            if spec is not None:
                break
        else:
            return None
        if type(spec.loader) is importlib.machinery.SourceFileLoader and spec.origin is not None:
            with unrecorded():
                vendored = os.path.realpath(spec.origin).startswith(cls.vendor)
            if not vendored:
                spec.loader = _PinnedSourceLoader(fullname, spec.origin)
        return spec


def _install_pinned_sources(recording: "Recording") -> None:
    """Compile every first-party module this process imports from now on through the pinned-source loader,
    and record its bytes on ``recording``.

    The modules already imported were compiled before the loader existed: this script, and what the
    interpreter's start imported (``site`` runs an editable install's finder module, which the manifest pins
    by its bytes). Their bytes are recorded now, the closest this process can come to the moment they were
    compiled."""
    global PINNED_SOURCES
    PINNED_SOURCES = recording.sources = {}
    with unrecorded():
        _PinnedSourceFinder.vendor = _vendor_dirs()
        started = {os.path.realpath(__file__)} | {
            os.path.realpath(f) for module in list(sys.modules.values())
            if isinstance(f := getattr(module, "__file__", None), str) and f.endswith(".py")
        }
        for path in started:
            with open(path, "rb") as source:
                recording.sources[path] = source.read()
    sys.meta_path.insert(0, _PinnedSourceFinder)


def main() -> None:
    """Record from the first line, then run the job; answer through the descriptor the parent reads.

    Standard output is redirected to standard error before the job runs, so a ``print`` in a generator
    cannot corrupt the answer (qa finding 9 of 2026-10-04).
    """
    # This file runs as ``__main__``; the memo imports it by name. Without the alias that import would load a
    # second copy whose list of active recordings is empty, so every memo entry the run read, and the memo's
    # own start of a nested generating process, would be offered to nobody.
    sys.modules["orpheus.numerics._traced_memo_boot"] = sys.modules[__name__]
    recording = Recording().start()
    _install_pinned_sources(recording)
    import pickle  # after the recorder starts, so its own import is recorded

    result_fd = os.dup(1)
    os.dup2(2, 1)
    envelope = pickle.loads(sys.stdin.buffer.read())
    sys.path[:] = envelope["sys_path"]
    from orpheus.numerics import traced_memo

    answer = traced_memo._generate_here(envelope["job"], recording)
    with os.fdopen(result_fd, "wb") as channel:
        channel.write(answer)
    # The entry is written and the answer sent; the process ends without tearing down its objects (the parsed
    # sources and the imported modules took 0.33 s to free per miss, ``[M]`` #602). Nothing after this line
    # could change what the parent reads.
    sys.stdout.flush()  # a generator's print, which went to standard error (qa finding 9)
    sys.stderr.flush()
    os._exit(0)


if __name__ == "__main__":  # an import (a roster walk over the package) runs nothing
    main()
