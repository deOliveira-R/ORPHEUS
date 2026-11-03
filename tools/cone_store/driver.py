"""``capture``: run tests under the grandparent/parent/child scheme and store their cones."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
import time
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from orpheus.numerics._traced_memo_boot import Witness

from tools.cone_store.collect import collect_ids, normalise
from tools.cone_store.cone import Cone, Reading
from tools.cone_store.environment import environment_key
from tools.cone_store.plan import Config, worker_command, worker_environment
from tools.cone_store.store import GRANDPARENT, ConeStore, StoredTest


@dataclass(frozen=True)
class CaptureReport:
    captured: int
    routes: dict[str, int]
    unknown: int
    not_passed: int
    pruned: int
    seconds: float

    def __str__(self) -> str:
        routes = ", ".join(f"{route} {count}" for route, count in sorted(self.routes.items()))
        return (f"captured {self.captured} tests ({routes}); {self.unknown} cone-unknown, {self.not_passed} not passed, "
                f"{self.pruned} pruned, {self.seconds:.1f} s")


def module_groups(repo: Path, expected: Sequence[str], arguments: Sequence[str]) -> list[tuple[str, tuple[str, ...]]]:
    """For each module holding an expected test, in order of appearance, the pytest arguments that collect it:
    the module's path, or the node ids the caller named in it."""
    named = [normalise(repo, a) for a in arguments]
    modules = list(dict.fromkeys(test_id.split("::")[0] for test_id in expected))
    groups = []
    for module in modules:
        mine = [a for a in named if a.split("::")[0] == module]
        whole = not named or any(a == module for a in mine) or any(
            module.startswith(a.split("::")[0].rstrip("/") + "/") for a in named if "::" not in a and a != module)
        groups.append((module, (module,) if whole or not mine else tuple(dict.fromkeys(mine))))
    return groups


def recent_identities(repo: Path, since_ns: int) -> dict[str, list[int]]:
    """The identity, at ``since_ns``, of every file a cone can pin (:meth:`Reading.kept`) that is not settled against
    it (written within the memo's ``TRUSTED_AGE_NS`` before, or after): the baseline a worker compares a file with,
    to tell a file written during the capture from one written just before it. The walk is the repository's, pruned
    where ``Reading`` keeps nothing (the bytecode cache) and at ``__pycache__`` and ``.git`` (git's own store, read by
    git, a process no test's cone admits); a file with no baseline that a capture finds unsettled makes its cone
    unknown, never a wrong one. Elegance's final pass on B2 (F0): this walked only the package and ``tests/`` while
    ``Reading`` counted ``tools/`` first-party, so a ``tools/`` file written just before a capture made every test
    that ran it unknown."""
    reading = Reading(repo)
    found: dict[str, list[int]] = {}
    for directory, names, files in os.walk(reading.repo):
        names[:] = [name for name in names
                    if name not in ("__pycache__", ".git") and reading.kept(os.path.join(directory, name) + os.sep)]
        for file in files:
            witness = Witness(reference_ns=since_ns)
            path = os.path.join(directory, file)
            stat = witness.stat(path)
            if stat is not None and not witness.settled:
                found[reading.relative(path)] = list(witness.seen[0][1] or ())
    return found


def capture(repo: str | Path, store: str | Path, package: str = "orpheus", memo_root: str | Path | None = None,
            ids: Sequence[str] = (), route: str = "auto", since_ns: int = 0) -> CaptureReport:
    """Run the tests ``ids`` names (default: all) and record each test's cone in ``store``."""
    started = time.monotonic()
    since_ns = since_ns or time.time_ns()
    repo = Path(repo).resolve()
    recent = recent_identities(repo, since_ns)
    collected = collect_ids(repo, package, ids)
    expected = collected.ids
    groups = module_groups(repo, expected, ids)
    key = environment_key(repo, collected.conftests)
    with tempfile.TemporaryDirectory(prefix="cone_capture_") as results:
        config = Path(results, "config.json")
        Config(str(repo), package, None if memo_root is None else str(Path(memo_root).resolve()), results,
               tuple(groups), route == "fresh", since_ns, recent).dump(config)
        run = subprocess.run(worker_command("gp", config), cwd=repo, env=worker_environment(repo),
                             capture_output=True, text=True)
        if run.returncode != 0:
            raise RuntimeError(f"cone_store: the grandparent exited {run.returncode}:\n{run.stdout[-3000:]}{run.stderr[-3000:]}")
        found = {}
        for path in Path(results).glob("*.json"):
            if path.name not in ("config.json", "gp.json", "ids.json", "errors.json"):
                data = json.loads(path.read_text())
                found[data["id"]] = data
        share = json.loads(Path(results, "gp.json").read_text()) if Path(results, "gp.json").exists() else None
    from orpheus.numerics import traced_memo

    with ConeStore(store) as database:
        for module in collected.errors:  # a module that does not import is recorded, not lost: select always lists it
            database.record(StoredTest(module, key, "none", "collection-error", ("collection error",)), (), {})
        for test_id in expected:
            data = found.get(test_id)
            if data is None:
                database.record(StoredTest(test_id, key, "none", "uncaptured", ("no worker produced a result",)), (), {})
                continue
            cone = Cone.from_json(traced_memo, data)
            own = Cone.from_json(traced_memo, data["own"]) if data["own"] is not None else Cone((), {}, ())
            database.record(StoredTest(test_id, key, data["route"], data["outcome"], cone.unknown), cone.pins, cone.files,
                            own.pins, own.files)
        pruned = _prune(database, [*expected, *collected.errors],
                        [module for module, args in groups if args == (module,)] + collected.errors, whole=not ids)
        # No grandparent cone (the fresh route, a grandparent that died): the stored one is another capture's.
        grandparent = Cone.from_json(traced_memo, share) if share is not None else Cone((), {}, ())
        database.record_share(GRANDPARENT, grandparent.pins, grandparent.files)
        database.set_meta("pruned", pruned)
        database.collect_garbage()
        stored = [database.tests()[t] for t in expected]
    routes = Counter(t.route for t in stored)
    return CaptureReport(len(stored), dict(routes), sum(1 for t in stored if t.unknown),
                         sum(1 for t in stored if t.outcome != "passed"), pruned, time.monotonic() - started)


def _prune(database: ConeStore, expected: Sequence[str], modules: Sequence[str], *, whole: bool) -> int:
    """Forget stored tests that no longer exist: those of a re-captured module that it no longer holds, and,
    after a capture of everything, those of a module that is gone. Returns how many."""
    keep = set(expected)
    touched = set(modules)
    gone = [test_id for test_id in database.tests()
            if test_id not in keep and (whole or test_id.split("::")[0] in touched)]
    return database.delete(gone)
