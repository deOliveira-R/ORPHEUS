"""The processes of a capture: ``python -m tools.cone_store.worker <role> <config.json> [arguments]``.

A capture is three levels of process (ruling 38 as corrected), one plugin class each, chosen once in :func:`main`:

* the GRANDPARENT (:class:`Grandparent`, role ``gp``) starts the memo's ``Recording`` in a
  ``pytest_load_initial_conftests`` hook wrapper, ahead of the root conftest, lets pytest import the conftests, and
  forks one COLLECTING PARENT per test module;
* the collecting parent (:class:`Collector`) is that fork continuing as an ordinary pytest session on the module's
  arguments: it collects the module, and runs each test in a CHILD forked from it;
* the child runs pytest's own protocol for one test, stops the recording, reads it as a cone
  (:func:`tools.cone_store.cone.cone_of`) and writes it to the results directory.

A child also records the test's OWN cone: a second ``Recording``, started when pytest's protocol for the test starts
and stopped when it returns. A recording reports each code object the first time it starts, per recording (each has
its own ``sys.monitoring`` tool), so a def the conftests ran is in the own cone again when the test runs it. The
inherited recording is the selection cone (``select`` judges it); the own one is what placement reads (B3).

The recording is inherited across both forks, so a child's cone is the grandparent's (conftest imports) plus its
module's collection plus its own run: no cache or ``sys.modules`` entry filled by an EARLIER module can hide code
from a later one, because the later module's parent was forked from the grandparent, not from the earlier module's.

The memo's fork guard (``fork_blockers``) decides before every fork whether the process may fork. A process that
may not fork runs the test in a FRESH interpreter (:class:`FreshChild`, role ``fresh``: the same pytest session on
one node id, recorded from its first line). A grandparent that may not fork spawns a :class:`FreshCollector` per
module (role ``collector_fresh``, no recorder, every test fresh); so does ``--route fresh`` (:func:`fresh_route`).

Every spawn and fork here is wrapped in the memo's ``spawning()``, so it is not recorded as the test's. Role
``list`` (:class:`Lister`) collects and writes the node ids pytest would run and the modules whose collection failed.
"""
from __future__ import annotations

import contextlib
import functools
import json
import os
import subprocess
import sys
import traceback
from collections.abc import Callable
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from tools.cone_store import _recorder, conftests
from tools.cone_store.plan import Config, pycache_dir, worker_command, worker_environment

#: No terminal reporter: its progress and width code ran inside each test's recorded protocol, differently on the
#: fork and fresh routes (each test's entry names took ``wcwidth`` and the progress line's names on one route only).
#: A capture's output is its log, never read for a verdict, which comes from the reports.
_FLAGS = ["-p", "no:cacheprovider", "-p", "no:terminal"]


class Lister:
    """Writes the node ids a collection yields to ``<results>/ids.json`` and the modules whose collection failed to
    ``<results>/errors.json``."""

    def __init__(self, config: Config) -> None:
        self.config = config
        self.errors: list[str] = []

    def pytest_collectreport(self, report: Any) -> None:
        if report.failed:
            self.errors.append(report.nodeid)

    def pytest_sessionfinish(self, session: Any) -> None:
        Path(self.config.results, "ids.json").write_text(json.dumps([item.nodeid for item in session.items]))
        Path(self.config.results, "errors.json").write_text(json.dumps(sorted(set(self.errors))))
        Path(self.config.results, "conftests.json").write_text(
            json.dumps(conftests.relative(self.config.repo, conftests.of_config(session.config))))


class Worker:
    """What every capture role shares: the configuration, the results directory and the way to start a process."""

    def __init__(self, config: Config, config_path: str) -> None:
        self.config = config
        self.config_path = config_path
        self.boot: ModuleType = _recorder.boot()

    def result_path(self, nodeid: str) -> Path:
        return Path(self.config.results, f"{self.boot.digest(nodeid.encode())[:40]}.json")

    def write(self, nodeid: str, route: str, outcome: str, cone: dict[str, Any] | None, unknown: tuple[str, ...] = (),
              own: dict[str, Any] | None = None, reads: list[list[Any]] | None = None) -> None:
        """The result of one test: its selection cone, its own cone under ``own`` (``None`` when no run was
        recorded, a crash) and what its own run read under ``reads`` (:class:`~tools.cone_store.reads.Read` rows)."""
        path = self.result_path(nodeid)
        body = cone or {"pins": [], "files": {}, "unknown": []}
        body["unknown"] = [*body["unknown"], *unknown]
        temporary = path.with_suffix(f".{os.getpid()}.tmp")
        temporary.write_text(json.dumps({"id": nodeid, "route": route, "outcome": outcome, **body, "own": own,
                                         "reads": reads or []}))
        os.replace(temporary, path)

    def spawn(self, role: str, *arguments: str) -> None:
        """A worker process of ``role``, waited for. Its output is the capture's log, not the user's."""
        with self.boot.spawning():
            subprocess.run(worker_command(role, Path(self.config_path), *arguments), cwd=self.config.repo,
                           env=worker_environment(Path(self.config.repo)),
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def fresh(self, nodeid: str) -> None:
        """Run one test in its own interpreter (a crash is recorded: the cone is unknown)."""
        self.spawn("fresh", nodeid)
        if not self.result_path(nodeid).exists():
            self.write(nodeid, "fresh", "crashed", None, ("fresh interpreter produced no result",))

    def finish(self) -> None:
        """Called when this process's pytest session returns."""

    @staticmethod
    def flush() -> None:
        sys.stdout.flush()
        sys.stderr.flush()


class FreshCollector(Worker):
    """A session over one module whose every test runs in its own interpreter."""

    def pytest_runtestloop(self, session: Any) -> bool | None:
        if session.config.option.collectonly:
            return None
        for item in session.items:
            self.fresh(item.nodeid)
        return True


def fresh_route(config: Config, config_path: str, _arguments: list[str]) -> int:
    """``--route fresh``: no recorder here; a fresh collector per module, a fresh interpreter per test."""
    worker = Worker(config, config_path)
    for _module, module_args in config.groups:
        worker.spawn("collector_fresh", *module_args)
    return 0


def _within(stretch: Callable[[], Any], function: Callable[..., Any]) -> Callable[..., Any]:
    """``function`` run inside ``stretch()``: the memo's ``unrecorded()`` or its ``recorded()`` (the one switch of
    what a recording takes, the boot module's)."""

    @functools.wraps(function)
    def within(*args: Any, **kwargs: Any) -> Any:
        with stretch():
            return function(*args, **kwargs)

    return within


def _quiet(method: Callable[..., Any]) -> Callable[..., Any]:
    """A role's own bookkeeping, run under the memo's ``unrecorded()``: the fork guard, the forks and the route each
    role takes are the capture's machinery, not the test's, and code started while it runs is recorded by no run.
    Without it a test's entry names depended on the route that captured it (the fork guard's ``ctypes``, the
    after-fork hooks), so the fork and fresh routes wrote different rows (G4). This module's own code objects are
    removed from the recording before it is read (:meth:`Recorded.read`), so the wrapper adds none."""

    @functools.wraps(method)
    def quietly(self: Worker, *args: Any, **kwargs: Any) -> Any:
        return _within(self.boot.unrecorded, method)(self, *args, **kwargs)

    return quietly


class Recorded(Worker):
    """A process recorded from before the conftests: the grandparent, a collecting parent, a fresh child."""

    def __init__(self, config: Config, config_path: str) -> None:
        super().__init__(config, config_path)
        self.recording: Any = None
        self.memo: Any = None
        self.store: Any = None
        self.reports: list[Any] = []
        self.exit = contextlib.ExitStack()

    def adopt(self, other: Recorded) -> Recorded:
        """Continue another recorded process's recording as this role (the fork turns a grandparent into a collector)."""
        self.recording, self.memo, self.store, self.exit = other.recording, other.memo, other.store, other.exit
        return self

    def start(self) -> None:
        from _pytest import config as pytest_config
        from _pytest import fixtures as pytest_fixtures
        from _pytest import main as pytest_main
        from _pytest import python as pytest_python
        from _pytest.assertion import rewrite
        from orpheus.numerics._traced_memo_boot import Recording  # the aliased module: imports no package

        # Pytest's own reads are not the test's dependencies: the assertion rewriter reads each test module's source
        # to compile it, and collection lists each directory to find the test files in it (a conftest that appears is
        # seen through pytest's probe of its path, a new test file changes no existing test). They run under the memo's
        # own `unrecorded()`, which the audit hook honours, so a read of a source by a TEST (a census, a docstring gate)
        # still is a dependency.
        # Pytest's sort of a module's items by fixture scope (``reorder_items``) is its scheduling, and it compares
        # items only when there are several: a collecting parent sorts the module's whole list, a fresh child its one
        # test, so the sort's code entered a test's cone on the fork route only (G4). A conftest's own
        # ``pytest_collection_modifyitems`` stays recorded.
        for module, name in ((rewrite, "_rewrite_test"), (pytest_main, "scandir"), (pytest_python, "scandir"),
                             (pytest_main, "_in_venv"), (pytest_fixtures, "reorder_items")):
            setattr(module, name, self.unrecorded_call(getattr(module, name)))
        # Distribution discovery lists every ``sys.path`` entry for metadata directories. Which distributions are
        # installed is the environment key's fact (``environment.distributions``), not a test's: a listing of the
        # repository root or of ``tests/`` taken for it made every file added or deleted there select every test.
        # Both metadata modules a process may load find distributions through ``FastPath.children``.
        for metadata in _metadata_modules():
            metadata.FastPath.children = self.unrecorded_call(metadata.FastPath.children)
        # Pytest's initial conftest discovery, in the grandparent, globs ``<testpath>/test*`` and asks ``is_dir()`` of
        # each entry, so every test FILE's presence was recorded there and attributed to every test: deleting any
        # test file selected the whole suite (S8[file]). The glob is silenced; the loading of each directory's
        # conftests (``_loadconftestmodules``: its probe of ``conftest.py`` and the import) stays recorded.
        manager = pytest_config.PytestPluginManager
        manager._try_load_conftest = self.unrecorded_call(manager._try_load_conftest)
        manager._loadconftestmodules = self.recorded_call(manager._loadconftestmodules)
        self.recording = Recording().start()
        from orpheus.numerics import traced_memo  # after the recorder starts, so its import is recorded

        self.memo = traced_memo
        if self.config.memo_root is not None:
            self.exit.enter_context(traced_memo.cache_root(Path(self.config.memo_root)))
        self.store = traced_memo.Store(Path(self.config.memo_root)) if self.config.memo_root else traced_memo.Store.current()

    def unrecorded_call(self, function: Callable[..., Any]) -> Callable[..., Any]:
        """``function`` run under the memo's ``unrecorded()``."""
        return _within(self.boot.unrecorded, function)

    def recorded_call(self, function: Callable[..., Any]) -> Callable[..., Any]:
        """``function`` run with recording on, inside an ``unrecorded_call`` (the memo's ``recorded()``)."""
        return _within(self.boot.recorded, function)

    def blockers(self) -> tuple[str, ...]:
        """Why this process must not fork now: the memo's guard, asked of the platform as it is now."""
        return self.boot.fork_blockers(self.memo.fork_platform())

    @pytest.hookimpl(wrapper=True, tryfirst=True)
    def pytest_load_initial_conftests(self, early_config: Any, parser: Any, args: list[str]) -> Any:
        self.start()
        result = yield
        self.load_every_conftest(early_config)
        self.silence_pytest_walk(early_config)
        self.conftests_loaded(early_config, args)
        return result

    def load_every_conftest(self, early_config: Any) -> None:
        """Import, recorded and before any test, every conftest a full run loads (:mod:`tools.cone_store.conftests`),
        through pytest's own conftest loading: their module code and session hooks act on every test of a full run,
        so they are in every test's cone (qa's second pass, D2: an off-path ``pytest_configure``)."""
        manager = early_config.pluginmanager
        with self.boot.unrecorded():  # the store's own walk: the environment key carries the set it finds
            found = conftests.of_config(early_config)
        for path in found:
            manager._loadconftestmodules(
                path.parent, early_config.known_args_namespace.importmode, early_config.rootpath,
                consider_namespace_packages=early_config.getini("consider_namespace_packages"))

    def silence_pytest_walk(self, early_config: Any) -> None:
        """Pytest's collection walk asks every entry of every directory it descends whether to ignore it, even for a
        session over one module, so each sibling file's presence would be every test's dependency (S8[file]). Only
        PYTEST'S own ``pytest_ignore_collect`` implementations run unrecorded; a conftest's or a plugin's stays
        recorded (qa's second pass, D3: a user hook that fills a cache)."""
        for impl in early_config.pluginmanager.hook.pytest_ignore_collect.get_hookimpls():
            if getattr(impl.plugin, "__name__", "").startswith("_pytest."):
                impl.function = self.unrecorded_call(impl.function)

    def conftests_loaded(self, early_config: Any, args: list[str]) -> None:
        """The conftests are imported and recorded."""

    def pytest_runtest_logreport(self, report: Any) -> None:
        self.reports.append(report)

    def fork(self, target: Callable[[], None]) -> int:
        """Run ``target`` in a forked child, which exits when it returns; the wait status."""
        self.flush()
        with self.boot.spawning():
            pid = os.fork()
        if pid == 0:
            code = 1
            try:
                target()
                code = 0
            except BaseException:
                traceback.print_exc()
            finally:
                self.flush()
                os._exit(code)
        return os.waitpid(pid, 0)[1]

    def run(self, item: Any, route: str) -> None:
        """Pytest's own protocol for ``item`` under a second recording (the test's own cone), then both recordings
        read as cones: the own one after the inherited one stopped, so the reading's imports are in neither. The own
        cone is what the own run did (:class:`~tools.cone_store.cone.ReadCones`): a memo entry the test read from
        another store is in the selection cone by the pins it was validated on, which another process ran."""
        self.reports.clear()
        protocol = self.recorded_call(item.ihook.pytest_runtest_protocol)
        own = self.boot.Recording().start()
        try:
            protocol(item=item, nextitem=None)  # the test: recorded by both
        finally:
            own.stop()
        outcome = "failed" if any(r.failed for r in self.reports) else \
            "skipped" if any(r.skipped for r in self.reports) else "passed"
        selection = self.read()
        function = getattr(item, "function", None)
        code = getattr(function, "__code__", None)
        test = None if code is None else (os.path.realpath(code.co_filename), function.__qualname__)
        cones = self.cone_of(own, test)
        self.write(item.nodeid, route, outcome, selection.to_json(self.memo), own=cones.ran.to_json(self.memo),
                   reads=[list(read.to_row()) for read in cones.reads])

    def read(self) -> Any:
        """The recording, stopped, read as its selection cone. The modules that read the cone are imported only now, after the
        recording stopped (imported during it, their import-time statements were pinned in every cone). A test of
        the cone store imports them itself, recorded."""
        self.recording.stop()
        return self.cone_of(self.recording).selection

    def cone_of(self, recording: Any, test: tuple[str, str] | None = None) -> Any:
        """A stopped recording read as its two cones. The capture's own code is not the test's: this module runs in every
        capture and in no full run, so its code objects are removed before the reading (`[M]` 2026-11-02, a capture
        of ``tests/gates/test_layer_imports.py``: six of its defs were pinned in each of 533 cones)."""
        mine = Recorded.cone_of.__code__.co_filename
        for key in [key for key in recording.codes if key[0] == mine]:
            del recording.codes[key]
        from tools.cone_store import cone

        return cone.cone_of(recording, self.store, self.config, self.memo, self.boot, test)


class FreshChild(Recorded):
    """A fresh interpreter over one node id, recorded from its first line. It starts as the grandparent does (no
    path argument, so the same conftests are discovered the same way) and takes its node id as the argument once
    they are imported, as a collecting parent takes its module: the two routes record the same."""

    def __init__(self, config: Config, config_path: str, nodeid: str) -> None:
        super().__init__(config, config_path)
        self.nodeid = nodeid

    @_quiet
    def conftests_loaded(self, early_config: Any, args: list[str]) -> None:
        args.append(self.nodeid)

    @_quiet
    def pytest_runtestloop(self, session: Any) -> bool | None:
        if session.config.option.collectonly:
            return None
        for item in session.items:
            if item.nodeid == self.nodeid:
                self.run(item, "fresh")
        return True


class Collector(Recorded):
    """A session over one module, forked from the grandparent; a child is forked per test."""

    def finish(self) -> None:
        self.flush()
        os._exit(0)  # a fork of the grandparent ends here, not in the grandparent's caller

    @_quiet
    def pytest_runtestloop(self, session: Any) -> bool | None:
        if session.config.option.collectonly:
            return None
        for item in session.items:
            if self.blockers():
                self.fresh(item.nodeid)
            else:
                status = self.fork(lambda item=item: self.run(item, "fork"))
                if not self.result_path(item.nodeid).exists():
                    self.write(item.nodeid, "fork", "crashed", None, (f"child died: wait status {status}",))
        return True


class Grandparent(Recorded):
    """The first process: records the conftests' imports, then forks a collecting parent per module."""

    collector: Collector | None = None

    def finish(self) -> None:
        if self.collector is not None:
            self.collector.finish()

    @_quiet
    def conftests_loaded(self, early_config: Any, args: list[str]) -> None:
        base = list(args)
        if not self.blockers():
            self.fork(self.share)
        for _module, module_args in self.config.groups:
            self.flush()
            if self.blockers():
                self.spawn("collector_fresh", *module_args)
                continue
            with self.boot.spawning():
                pid = os.fork()
            if pid == 0:
                self.collector = Collector(self.config, self.config_path)
                self.collector.adopt(self)
                early_config.pluginmanager.register(self.collector)
                args[:] = [*base, *module_args]
                return  # this process goes on as pytest's session over this module
            os.waitpid(pid, 0)
        self.flush()
        os._exit(0)

    def share(self) -> None:
        """A child of the grandparent that writes the cone the grandparent attributes to every test."""
        Path(self.config.results, "gp.json").write_text(json.dumps(self.read().to_json(self.memo)))


def _session(plugin: Callable[[Config, str], Worker]) -> Callable[[Config, str, list[str]], int]:
    def run(config: Config, config_path: str, rest: list[str]) -> int:
        worker = plugin(config, config_path)
        pytest.main([*_FLAGS, *rest], plugins=[worker])
        worker.finish()
        return 0

    return run


def _grandparent(config: Config, config_path: str, rest: list[str]) -> int:
    return (fresh_route if config.fresh_only else _session(Grandparent))(config, config_path, rest)


def _fresh_child(config: Config, config_path: str, rest: list[str]) -> int:
    pytest.main(_FLAGS, plugins=[FreshChild(config, config_path, rest[-1])])
    return 0


def _list(config: Config, _config_path: str, rest: list[str]) -> int:
    pytest.main(["--collect-only", "--continue-on-collection-errors", *_FLAGS, *rest], plugins=[Lister(config)])
    return 0


def _metadata_modules() -> list[Any]:
    """The distribution-metadata modules this process can load: the standard library's and, when installed, the
    ``importlib_metadata`` backport (a plugin may import either)."""
    import importlib.metadata

    found: list[Any] = [importlib.metadata]
    try:
        import importlib_metadata
    except ImportError:
        return found
    return [*found, importlib_metadata]


#: The role named on the command line, to the function that runs it: the one place a role is chosen.
ROLES: dict[str, Callable[[Config, str, list[str]], int]] = {
    "gp": _grandparent,
    "collector_fresh": _session(FreshCollector),
    "fresh": _fresh_child,
    "list": _list,
}


def main(argv: list[str]) -> int:
    role, config_path, *rest = argv
    config = Config.load(config_path)
    sys.pycache_prefix = pycache_dir(config.repo)  # this process's bytecode only; no interpreter a test starts inherits it
    # No bytecode is ever written, so the prefix stays empty and every capture compiles each module the same way: the
    # import system's cache-write code (``_write_atomic`` and its strings) ran only in a capture that met a cold
    # cache, and put its names in that test's entry names alone (the routes' rows differed, G4 and C6).
    sys.dont_write_bytecode = True
    return ROLES[role](config, config_path, rest)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
