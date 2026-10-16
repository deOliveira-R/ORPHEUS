"""The traced memo's recorder, and the first code its generating process runs (#405 P3).

This module imports only the standard library, so the generating process can start recording before any
``orpheus`` code runs: it is run by path with ``-P`` (``python -P _traced_memo_boot.py``), and its
``main`` starts a :class:`Recording` on its first line, then loads the job and
:mod:`orpheus.numerics.traced_memo`. Importing it runs nothing (a roster walk over the package imports it),
and :func:`orpheus.numerics.traced_memo.trace_call` records through the same :class:`Recording`: one
definition of what a generation depends on.

A recording holds, for the time it runs, across every thread of the process:

* every code object that started (``sys.monitoring`` ``PY_START``; each reports once, then disables itself);
* every path opened for reading, with the digest of its bytes, and every path whose presence was probed
  (``os.stat``, ``os.lstat``, and so ``os.path.exists``, ``isfile``, ``isdir`` and ``pathlib``'s ``exists``),
  with what was there, each made absolute at the moment it was touched, a bytes path decoded, a missing path
  included (the run then depended on its absence);
* every directory listed (``os.listdir``, ``os.scandir``), with the digest of its names;
* every path the run wrote, truncated, renamed or removed: its own output.

Each input is observed when the run first touches it, never after the run: an entry pins what the run saw, so a
file, a flag or a directory changed while the run used it leaves the entry born stale instead of serving the
old answer under the new state (qa's reviews of #602 cut E, ``[M]`` 2026-10-16). A path the run touches only
after it has written it is its own output, not an input. A path the run observed and THEN changed is refused
by the manifest: the state the run saw no longer exists, so no later state pins the answer.
* every process the run started, by its program or, when none can be named, by its event;
* the working directory the run started in, and whether it read a relative path;
* the memo entries the run read, which :mod:`~orpheus.numerics.traced_memo` adds as it serves them.

The import system's own reads, probes and listings are not recorded: they are how it finds a module, and the
module it finds is pinned where its code runs. Recording them would make every new file in any ``sys.path``
directory invalidate every entry. The import system is the frozen bootstrap and the ``importlib`` package
(``importlib.metadata`` lists every ``sys.path`` directory to find the installed distributions, which are
pinned by version and RECORD). The exemption is decided by the IMMEDIATE caller (the import-system frame that made the
call), never by an import being somewhere on the stack: a module body that lists a
directory while it is imported is the run's own dependency (qa's second review, 2026-10-04).

Audit hooks cannot be removed, so one hook is installed per process and serves the recordings active at the
moment of each event; the ``os.stat`` probe is installed while a recording is active and removed after.
"""
import copy
import hashlib
import importlib.machinery
import marshal
import os
import struct
import sys
import threading
from collections.abc import Callable
from typing import NamedTuple

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

#: An open reads a path's prior bytes, changes them, or both (:func:`_access`).
_OPEN_EVENTS = {"open"}
#: The events by which the run changes a path without opening it.
_CHANGE_EVENTS = {"os.remove", "os.rename", "os.truncate"}
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
#: lists every ``sys.path`` directory to find the installed distributions, which are pinned by version and RECORD).
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


def digest(data: bytes) -> str:
    """The memo's one digest: of a source, a data file, a manifest's payload, an index's key."""
    return hashlib.blake2b(data, digest_size=32).hexdigest()


def file_digest(path: str) -> str:
    """The digest of a regular file's bytes, or ``"absent"``: a run that found no regular file there depended on
    that. A directory, a device or a pipe is never read: ``/dev/urandom`` has no end (qa's second review of #602
    cut E, ``[M]`` 2026-10-16: a run that read it hung)."""
    if not os.path.isfile(path):
        return "absent"
    try:
        with open(path, "rb") as file:
            return digest(file.read())
    except FileNotFoundError:
        return "absent"


def presence(path: str) -> str:
    """What a probe of ``path`` finds: a file, a directory, something else, or nothing."""
    if os.path.isfile(path):
        return "file"
    if os.path.isdir(path):
        return "directory"
    return "other" if os.path.lexists(path) else "absent"


def listing_digest(path: str) -> str:
    """The digest of a directory's names, or ``"absent"``."""
    return digest("\0".join(sorted(os.listdir(path))).encode()) if os.path.isdir(path) else "absent"


class Observation(NamedTuple):
    """One kind of input: the :class:`Recording` field that maps each path to what the run saw there, what is
    seen at a path, and whether the run's own writes change what is seen there."""

    field: str
    observe: Callable[[str], str]
    changed_by: Callable[[str, set[str]], bool]


#: A file's bytes, a path's presence, a directory's names. A write into a directory changes its names.
CONTENT = Observation("opened", file_digest, lambda path, written: path in written)
PRESENCE = Observation("probed", presence, lambda path, written: path in written)
LISTING = Observation("listed", listing_digest,
                      lambda path, written: path in written or any(os.path.dirname(w) == path for w in written))


def _observe(kind: Observation, raw: object) -> None:
    """An input the run touched, observed at its first touch (the module docstring). A path the run has already
    written is its own output and is not observed."""
    found = _path(raw)
    if found is None:
        return
    path, relative = found
    seen: str | None = None
    for recording in list(_ACTIVE):
        observed = getattr(recording, kind.field)
        if path not in observed and path not in recording.written:
            if seen is None:
                with unrecorded():
                    seen = kind.observe(path)
            observed[path] = seen
        recording.relative |= relative


def _access(flags: object) -> tuple[bool, bool]:
    """Whether an ``open`` event reads the path's prior bytes, and whether it changes them, from its ``O_*``
    flags, which both of its emitters carry (``io.open``'s mode string is redundant with them: ``"r+"`` is
    ``O_RDWR``, ``"w+"`` adds ``O_TRUNC``, ``"a+"`` ``O_APPEND``). A truncating or exclusive open reads nothing
    that was there. The elegance review of #602 cut E (``[M]`` 2026-10-16): an ``"r+"`` open was taken for a
    write alone, so what it read was never pinned."""
    flags = flags if isinstance(flags, int) else 0
    access = flags & os.O_ACCMODE
    reads = access != os.O_WRONLY and not flags & (os.O_TRUNC | os.O_EXCL)
    writes = access != os.O_RDONLY or bool(flags & (os.O_CREAT | os.O_TRUNC | os.O_APPEND))
    return reads, writes


def _note_written(raw: object) -> None:
    """A path the run wrote, truncated, renamed or removed: its own output. A generator that writes a scratch
    file, reads it back and deletes it reads only its own output (qa's second review of #602 cut E, ``[M]``
    2026-10-16)."""
    found = _path(raw)
    if found is None:
        return
    for recording in list(_ACTIVE):
        recording.written.add(found[0])


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
    if event in _OPEN_EVENTS and not _by_importer(2):
        path, _, flags = (*args, None, None, 0)[:3]
        reads, writes = _access(flags)
        if reads:
            _observe(CONTENT, path)  # first: what it read is an input, refused if the run then changed it
        if writes:
            _note_written(path)
    elif event in _CHANGE_EVENTS and not _by_importer(2):
        for path in args[: 2 if event == "os.rename" else 1]:
            _note_written(path)
    elif event in _LIST_EVENTS and not _by_importer(2):
        _observe(LISTING, args[0] if args else None)
    elif event in _SPAWN_EVENTS and not getattr(_SPAWNING, "on", False):
        program = _program(event, args)
        for recording in list(_ACTIVE):
            recording.spawned.add(program)


def _probing(name: str):
    original = _STAT[name]

    def probe(path, *args, **kwargs):
        if _ACTIVE and not getattr(_QUIET, "on", False) and not _by_importer(2):
            _observe(PRESENCE, path)
        return original(path, *args, **kwargs)

    return probe


#: What a recording holds about its run, the fields :meth:`Recording.snapshot` copies and
#: :meth:`Recording.restore` puts back: one list, so a field is carried by both or by neither.
_RUN_STATE = ("codes", "opened", "written", "probed", "listed", "spawned", "children", "relative", "directory", "environment")


class Recording:
    """What one run depends on, recorded while it runs (module docstring)."""

    def __init__(self) -> None:
        self.codes: set = set()
        #: Each input path, with what the run found at its first touch (the module docstring).
        self.opened: dict[str, str] = {}
        self.probed: dict[str, str] = {}
        self.listed: dict[str, str] = {}
        #: Each path the run wrote, truncated, renamed or removed.
        self.written: set[str] = set()
        self.spawned: set[str] = set()
        self.relative = False
        self.directory = os.getcwd()
        self.environment: dict[str, str] = dict(os.environ)
        self.children: set[tuple[str, str, str]] = set()
        #: The bytes each first-party module was compiled from, by real path, when the pinned-source loader
        #: is installed (a generating process); ``None`` in a process without it (``trace_call``). The
        #: process's record, not the run's: every recording of the process shares it.
        self.sources: dict[str, bytes] | None = PINNED_SOURCES
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
            recording = self

            def started(code, _offset):
                recording.codes.add(code)
                return monitoring.DISABLE

            monitoring.register_callback(self._tool, monitoring.events.PY_START, started)
            monitoring.set_events(self._tool, monitoring.events.PY_START)
            if not _ACTIVE:
                os.stat, os.lstat = _probing("stat"), _probing("lstat")
            self.directory = os.getcwd()
            self.environment = dict(os.environ)  # what the run RECEIVED; a run may set a variable later
            _ACTIVE.append(self)
        return self

    def snapshot(self) -> "Recording":
        """A copy of what this recording holds now, not started: the state a forked generating child is put back
        to (:meth:`restore`), and the run whose manifest a template re-validates before each fork
        (:func:`template`)."""
        held = Recording()
        held._take(self)
        return held

    def restore(self, held: "Recording") -> None:
        """Put back what :meth:`snapshot` held, and re-enable the events of every code object, so a code object
        that already reported (and disabled itself) reports again when it next starts."""
        self._take(held)
        sys.monitoring.restart_events()

    def _take(self, other: "Recording") -> None:
        for name in _RUN_STATE:
            setattr(self, name, copy.copy(getattr(other, name)))

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


def frame(data: bytes) -> bytes:
    """``data`` as one length-prefixed message, the unit of every pipe between the memo and its processes."""
    return struct.pack(">Q", len(data)) + data


def _read_frame(stream) -> bytes | None:
    """One length-prefixed message from ``stream``; ``None`` at its end, or when it ends inside the message."""
    head = stream.read(8)
    if len(head) < 8:
        return None
    (size,) = struct.unpack(">Q", head)
    data = stream.read(size)
    return data if len(data) == size else None


def _write_frame(fd: int, data: bytes) -> None:
    view = memoryview(frame(data))
    while view:
        view = view[os.write(fd, view):]


def _libproc_thread_counter() -> Callable[[], int]:
    """A function counting this process's OS threads through libproc, native ones included (a BLAS pool, which
    ``threading`` cannot see)."""
    import ctypes

    class TaskInfo(ctypes.Structure):  # libproc's struct proc_taskinfo
        _fields_ = [(name, ctypes.c_uint64) for name in ("virtual", "resident", "user", "system", "threads_user",
                                                          "threads_system")] + \
                   [(name, ctypes.c_int32) for name in ("policy", "faults", "pageins", "cow_faults", "sent", "received",
                                                         "mach", "unix", "csw", "threadnum", "running", "priority")]

    libproc = ctypes.CDLL("/usr/lib/libproc.dylib")
    pid = os.getpid()

    def count() -> int:
        info = TaskInfo()
        if libproc.proc_pidinfo(pid, 4, ctypes.c_uint64(0), ctypes.byref(info), ctypes.sizeof(info)) != ctypes.sizeof(info):
            raise OSError("proc_pidinfo(PROC_PIDTASKINFO) failed")
        return info.threadnum

    return count


#: The platforms on which a template may fork its generating children, as (operating system, machine, BLAS),
#: each with the counter of its OS threads: one table, the one door.
#: SCOPE-BOUNDARY[guard] machinery: the fork-safety protocol of #602 run on another platform (thousands of forks
#: per thread setting, each child's answer compared with a fresh interpreter's), and an OS-thread counter for it.
#: ruling: the user, 2026-10-15, "land it only if that gate holds"; measured 2026-10-16 on darwin/arm64/Accelerate
#: only (32 400 forks, 0 failures). revisit: when the protocol runs on Linux/OpenBLAS (CI); until then every
#: other platform spawns a fresh interpreter per miss.
FORK_MEASURED: dict[tuple[str, str, str], Callable[[], Callable[[], int]]] = {
    ("darwin", "arm64", "accelerate"): _libproc_thread_counter,
}


def _fork_unsafe(count_threads: Callable[[], int]) -> tuple[str, ...]:
    """Why this process must not fork now: it holds more than one OS thread (a fork copies one thread and every
    lock the others held; a fork after a BLAS call crashed its child, ``[M]`` 2026-10-16, 20 of 20). Empty when
    it may fork."""
    threads = count_threads()
    return () if threads == 1 else (f"{threads} OS threads",)


def _boot():
    """The first code of a generating process and of a template: record from the first line, compile every
    first-party module from the bytes it pins, reserve the descriptor the parent reads, read the envelope
    (``sys.path``), and import the memo. Both start with the argument list a fresh process has, so no generator
    can read which of the two it runs in (qa's review of #602 cut E: ``len(sys.argv)`` was 2 forked, 1 fresh)."""
    # This file runs as ``__main__``; the memo imports it by name. Without the alias that import would load a
    # second copy whose list of active recordings is empty, so every memo entry the run read, and the memo's
    # own start of a nested generating process, would be offered to nobody.
    sys.modules["orpheus.numerics._traced_memo_boot"] = sys.modules[__name__]
    del sys.argv[1:]
    recording = Recording().start()
    _install_pinned_sources(recording)
    import pickle  # after the recorder starts, so its own import is recorded

    # Standard output goes to standard error, so a ``print`` in a generator cannot corrupt the answer (qa
    # finding 9 of 2026-10-04).
    answers = os.dup(1)
    os.dup2(2, 1)
    jobs = sys.stdin.buffer
    envelope = _read_frame(jobs)
    if envelope is None:
        os._exit(0)
    sys.path[:] = pickle.loads(envelope)["sys_path"]
    from orpheus.numerics import traced_memo

    return recording, answers, jobs, traced_memo


def _end() -> None:
    """End without tearing down the process's objects (the parsed sources and the imported modules took 0.33 s
    to free per miss, ``[M]`` #602): the answer is sent, and nothing after it could change what the parent reads."""
    sys.stdout.flush()  # a generator's print, which went to standard error (qa finding 9)
    sys.stderr.flush()
    os._exit(0)


def template() -> None:
    """A template process (#602, cut E): it starts as a generating process does (:func:`_boot`), then forks one
    generating child per job read from the parent.

    The child is the generating process a fresh interpreter would be, without paying the interpreter's start
    and the memo's imports again (about 0.5 s per miss). It inherits the recording as it stood when the imports
    finished, with the events of every code object re-enabled (:meth:`Recording.restore`), so its manifest
    holds what a fresh process's would, and its standard input is ``/dev/null``, as a fresh process's is at its
    end. Before each fork the template checks that the manifest of its imports still validates, the bytes of
    every first-party source it compiled included (an edit to a def body the imports never ran, a data file, a
    listing or a distribution they read), and that it runs one OS thread; if not, it answers
    :class:`~orpheus.numerics.traced_memo._Refused` and exits, and the parent starts another. A parent that
    ends closes the job pipe, and the template exits at its end.

    Protocol, on file descriptors: the parent writes frames on standard input (the envelope, then one job per
    frame); the template answers one frame per job on the descriptor that was standard output: the child's
    answer, :class:`~orpheus.numerics.traced_memo._Refused`, or :class:`~orpheus.numerics.traced_memo._Died`."""
    recording, answers, jobs, traced_memo = _boot()
    import pickle

    ready = recording.snapshot()
    with unrecorded():
        count_threads = FORK_MEASURED[traced_memo.fork_platform()]()
        # What the imports depended on, pinned as a generation pins it, with the whole bytes of every source
        # they compiled: module state built from a state that has changed since would run old values under a
        # manifest naming the new ones (the test-architect's design note, #602 cut E).
        store = traced_memo.Store.current()
        imported = traced_memo._manifest(ready, store).with_compiled_sources(PINNED_SOURCES or {})
    while (job := _read_frame(jobs)) is not None:
        with unrecorded():
            reasons = traced_memo.validate(imported, store) or _fork_unsafe(count_threads)
        if reasons:
            _write_frame(answers, pickle.dumps(traced_memo._Refused(reasons)))
            break
        read_end, write_end = os.pipe()
        with spawning():
            pid = os.fork()
        if pid == 0:
            os.close(read_end)
            os.dup2(os.open(os.devnull, os.O_RDONLY), 0)
            recording.restore(ready)
            answer = traced_memo._generate_here(job, recording)
            with os.fdopen(write_end, "wb") as channel:
                channel.write(answer)
            _end()
        os.close(write_end)
        with os.fdopen(read_end, "rb") as channel:
            answer = channel.read()
        _, status = os.waitpid(pid, 0)
        _write_frame(answers, answer if answer else pickle.dumps(traced_memo._Died(status)))
    _end()


def main() -> None:
    """A fresh generating process: start (:func:`_boot`), run the one job, answer it, end."""
    recording, answers, jobs, traced_memo = _boot()
    job = _read_frame(jobs)
    if job is not None:
        _write_frame(answers, traced_memo._generate_here(job, recording))
    _end()


if __name__ == "__main__":  # an import (a roster walk over the package) runs nothing
    template() if sys.argv[1:] == ["--template"] else main()
