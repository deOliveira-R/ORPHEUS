"""The processes of a capture: ``python -m tools.cone_store.worker <role> <config.json> [arguments]``.

A capture is three levels of process (ruling 38 as corrected), one plugin class each, chosen once in :func:`main`:

* the GRANDPARENT (:class:`Grandparent`, role ``gp``) starts the memo's ``Recording`` in a
  ``pytest_load_initial_conftests`` hook wrapper, ahead of the root conftest, lets pytest import the conftests, and
  forks one COLLECTING PARENT per test module;
* the collecting parent (:class:`Collector`) is that fork continuing as an ordinary pytest session on the module's
  arguments: it collects the module, and runs each test in a CHILD forked from it;
* the child runs pytest's own protocol for one test, stops the recording, reads it as a cone
  (:func:`tools.cone_store.cone.cone_of`) and writes it to the results directory.

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

_FLAGS = ["-p", "no:cacheprovider", "-q", "--color=no"]


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

    def write(self, nodeid: str, route: str, outcome: str, cone: dict[str, Any] | None, unknown: tuple[str, ...] = ()) -> None:
        path = self.result_path(nodeid)
        body = cone or {"pins": [], "files": {}, "unknown": []}
        body["unknown"] = [*body["unknown"], *unknown]
        temporary = path.with_suffix(f".{os.getpid()}.tmp")
        temporary.write_text(json.dumps({"id": nodeid, "route": route, "outcome": outcome, **body}))
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


class Recorded(Worker):
    """A process recorded from before the conftests: the grandparent, a collecting parent, a fresh child."""

    def __init__(self, config: Config, config_path: str) -> None:
        super().__init__(config, config_path)
        self.recording: Any = None
        self.memo: Any = None
        self.cone: Any = None
        self.store: Any = None
        self.reports: list[Any] = []
        self.exit = contextlib.ExitStack()

    def adopt(self, other: Recorded) -> Recorded:
        """Continue another recorded process's recording as this role (the fork turns a grandparent into a collector)."""
        self.recording, self.memo, self.cone, self.store, self.exit = other.recording, other.memo, other.cone, other.store, other.exit
        return self

    def start(self) -> None:
        from _pytest import config as pytest_config
        from _pytest import main as pytest_main
        from _pytest import python as pytest_python
        from _pytest.assertion import rewrite
        from orpheus.numerics._traced_memo_boot import Recording  # the aliased module: imports no package

        # Pytest's own reads are not the test's dependencies: the assertion rewriter reads each test module's source
        # to compile it, and collection lists each directory to find the test files in it (a conftest that appears is
        # seen through pytest's probe of its path, a new test file changes no existing test). They run under the memo's
        # own `unrecorded()`, which the audit hook honours, so a read of a source by a TEST (a census, a docstring gate)
        # still is a dependency.
        for module, name in ((rewrite, "_rewrite_test"), (pytest_main, "scandir"), (pytest_python, "scandir"),
                             (pytest_main, "_in_venv")):
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

        from tools.cone_store import cone

        self.memo, self.cone = traced_memo, cone
        if self.config.memo_root is not None:
            self.exit.enter_context(traced_memo.cache_root(Path(self.config.memo_root)))
        self.store = traced_memo.Store(Path(self.config.memo_root)) if self.config.memo_root else traced_memo.Store.current()

    def unrecorded_call(self, function: Callable[..., Any]) -> Callable[..., Any]:
        """``function`` run under the memo's ``unrecorded()``."""
        boot = self.boot

        def quiet(*args: Any, **kwargs: Any) -> Any:
            with boot.unrecorded():
                return function(*args, **kwargs)

        return quiet

    def recorded_call(self, function: Callable[..., Any]) -> Callable[..., Any]:
        """``function`` run with recording on, inside an ``unrecorded_call`` (the memo's ``unrecorded()`` keeps the
        outer state on its instance; the flag it sets is the boot module's ``_QUIET``)."""
        quiet_flag = self.boot._QUIET

        def loud(*args: Any, **kwargs: Any) -> Any:
            outer = getattr(quiet_flag, "on", False)
            quiet_flag.on = False
            try:
                return function(*args, **kwargs)
            finally:
                quiet_flag.on = outer

        return loud

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
        """Pytest's own protocol for ``item``, then the recording read as a cone."""
        self.reports.clear()
        item.ihook.pytest_runtest_protocol(item=item, nextitem=None)
        outcome = "failed" if any(r.failed for r in self.reports) else \
            "skipped" if any(r.skipped for r in self.reports) else "passed"
        self.write(item.nodeid, route, outcome, self.read().to_json(self.memo))

    def read(self) -> Any:
        self.recording.stop()
        return self.cone.cone_of(self.recording, self.store, self.config, self.memo, self.boot)


class FreshChild(Recorded):
    """A fresh interpreter over one node id, recorded from its first line. It starts as the grandparent does (no
    path argument, so the same conftests are discovered the same way) and takes its node id as the argument once
    they are imported, as a collecting parent takes its module: the two routes record the same."""

    def __init__(self, config: Config, config_path: str, nodeid: str) -> None:
        super().__init__(config, config_path)
        self.nodeid = nodeid

    def conftests_loaded(self, early_config: Any, args: list[str]) -> None:
        args.append(self.nodeid)

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
        """A child of the grandparent that writes how many pins the grandparent attributes to every test."""
        cone = self.read()
        Path(self.config.results, "gp.json").write_text(json.dumps({"pins": len(cone.pins), "unknown": list(cone.unknown)}))


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
    return ROLES[role](config, config_path, rest)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
