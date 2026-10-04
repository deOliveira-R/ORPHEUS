"""The traced memo's recorder, and the first code its generating process runs (#405 P3).

This module imports only the standard library, so the generating process can start recording before any
``orpheus`` code runs: it is run by path with ``-P`` (``python -P _traced_memo_boot.py``), and its
``main`` starts a :class:`Recording` on its first line, then loads the job and
:mod:`orpheus.numerics.traced_memo`. Importing it runs nothing (a roster walk over the package imports it),
and :func:`orpheus.numerics.traced_memo.trace_call` records through the same :class:`Recording`: one
definition of what a generation depends on.

A recording holds, for the time it runs, across every thread of the process:

* every code object that started (``sys.monitoring`` ``PY_START``; each reports once, then disables itself);
* every path opened for reading, made absolute at the moment it was opened (a later change of working
  directory cannot move it), a bytes path decoded, and a path that does not exist included (the run then
  depended on its absence);
* every directory listed (``os.listdir``, ``os.scandir``);
* every process the run started: the program, or the event when no program can be named (a ``fork``, an
  ``os.system`` command line);
* the memo entries the run read, which :mod:`~orpheus.numerics.traced_memo` adds as it serves them.

Audit hooks cannot be removed, so one hook is installed per process and serves the recordings active at the
moment of each event.
"""
import os
import sys
import threading

#: The recordings active in this process, innermost last; every event is offered to each of them.
_ACTIVE: list["Recording"] = []
_LOCK = threading.Lock()
_HOOKED = False
#: Set while the memo itself starts a generating process, so that spawn is not the run's own.
_SPAWNING = threading.local()

_READ_EVENTS = {"open"}
_LIST_EVENTS = {"os.listdir", "os.scandir"}
_SPAWN_EVENTS = {"subprocess.Popen", "os.fork", "os.forkpty", "os.posix_spawn", "os.spawn", "os.exec", "os.system"}


def _path(raw: object) -> tuple[str, bool] | None:
    """A path made absolute now, and whether it was relative (the run then depended on the working
    directory); ``None`` for a file descriptor."""
    if isinstance(raw, int):
        return None
    try:
        text = os.fsdecode(raw if raw is not None else ".")  # type: ignore[arg-type]
    except TypeError:
        return None
    return os.path.abspath(text), not os.path.isabs(text)


def _offer(field: str, raw: object) -> None:
    found = _path(raw)
    if found is None:
        return
    path, relative = found
    for recording in list(_ACTIVE):
        getattr(recording, field).add(path)
        recording.relative |= relative


def _importing() -> bool:
    """Whether the import system is the caller: its directory listings (``FileFinder`` filling its cache) are
    how it finds a module, and the module it found is pinned where its code runs. Pinning them would make
    every new file in any ``sys.path`` directory invalidate every entry."""
    frame = sys._getframe(2)
    while frame is not None:
        if frame.f_code.co_filename.startswith("<frozen importlib"):
            return True
        frame = frame.f_back
    return False


def _audit(event: str, args: tuple[object, ...]) -> None:
    if not _ACTIVE:
        return
    if event in _READ_EVENTS:
        mode = args[1] if len(args) > 1 and isinstance(args[1], str) else "r"
        if not any(c in mode for c in "wax+"):
            _offer("opened", args[0] if args else None)
    elif event in _LIST_EVENTS and not _importing():
        _offer("listed", args[0] if args else None)
    elif event in _SPAWN_EVENTS and not getattr(_SPAWNING, "on", False):
        program = _program(event, args)
        for recording in list(_ACTIVE):
            recording.spawned.add(program)


def _program(event: str, args: tuple[object, ...]) -> str:
    """The absolute path of the program a spawn event runs, or the event itself when none can be named."""
    import shutil

    if event == "subprocess.Popen":
        executable, argv = (args + (None, None))[:2]
        name = executable if executable is not None else (argv[0] if isinstance(argv, (list, tuple)) and argv else argv)
        if isinstance(name, (str, bytes, os.PathLike)):
            found = shutil.which(os.fsdecode(name))
            if found is not None:
                return os.path.realpath(found)
    return event


class Recording:
    """What one run depends on, recorded while it runs (module docstring)."""

    def __init__(self) -> None:
        self.codes: set = set()
        self.opened: set[str] = set()
        self.listed: set[str] = set()
        self.spawned: set[str] = set()
        self.relative = False
        self.children: set[tuple[str, str, str]] = set()
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


def note_child(child: tuple[str, str, str]) -> None:
    """A memo entry a run read: every active recording pins it."""
    for recording in list(_ACTIVE):
        recording.children.add(child)


class spawning:
    """The memo's own start of a generating process, which is not a dependency of the run that asked."""

    def __enter__(self) -> None:
        _SPAWNING.on = True

    def __exit__(self, *_exc: object) -> None:
        _SPAWNING.on = False


def main() -> None:
    """Record from the first line, then run the job; answer through the result descriptor the parent passed.

    Standard output is redirected to standard error before the job runs, so a ``print`` in a generator
    cannot corrupt the answer (qa finding 9 of 2026-10-04).
    """
    # This file runs as ``__main__``; the memo imports it by name. Without the alias that import would load a
    # second copy whose list of active recordings is empty, so every memo entry the run read, and the memo's
    # own start of a nested generating process, would be offered to nobody.
    sys.modules["orpheus.numerics._traced_memo_boot"] = sys.modules[__name__]
    recording = Recording().start()
    import pickle  # after the recorder starts, so its own import is recorded

    result_fd = os.dup(1)
    os.dup2(2, 1)
    envelope = pickle.loads(sys.stdin.buffer.read())
    sys.path[:] = envelope["sys_path"]
    from orpheus.numerics import traced_memo

    answer = traced_memo._generate_here(envelope["job"], recording)
    with os.fdopen(result_fd, "wb") as channel:
        channel.write(answer)


if __name__ == "__main__":  # an import (a roster walk over the package) runs nothing
    main()
