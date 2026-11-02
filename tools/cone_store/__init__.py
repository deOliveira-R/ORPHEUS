"""The cone store: per test, the first-party code and data the test's verdict depends on, and from it the
tests a working-tree change requires.

    python -m tools.cone_store {capture,select,status} --repo R --store S [--package P] [--memo-root M] [ids]

``capture`` runs the tests under a grandparent process that records (the memo's ``Recording``) from before the
root conftest is imported, forks one collecting parent per test module, and forks one child per test; each child's
cone is the recording's pins (``orpheus.numerics.traced_memo.collect_pins``): defs and module skeletons of code
under ``<repo>/<package>/`` and ``<repo>/tests/``, data, presence and listing pins of paths under the repository,
memo children read and distributions. ``select`` lists the collected tests whose cone is stale, new, unknown,
captured under another environment key, or not passing at capture, and a module that failed to collect. The store
is SQLite (``.cache/cone_store/cones.sqlite`` by default); see :mod:`tools.cone_store.store` for its layout and
:mod:`tools.cone_store.worker` for the processes. A process that may not fork (the memo's ``fork_blockers``: an
unmeasured platform, more than one OS thread) runs the test in a fresh interpreter instead.

SCOPE BOUNDARIES (what a green ``select`` does not say):

* SCOPE-BOUNDARY[guard] machinery: a recorder of environment reads (no audit event fires on ``os.environ``), #687.
  ruling: the verification-tree plan, "The cone store -- design", rulings 34 and 42 (the key names what it names).
  revisit: a test whose verdict follows a variable the key does not name. The key holds the interpreter, the
  platform and BLAS, the installed distributions and their ``pytest11`` entry points, ``PYTHONHASHSEED`` and the
  thread-count variables the memo fixes, ``PYTEST_ADDOPTS``, ``PYTEST_PLUGINS``, ``PYTEST_DISABLE_PLUGIN_AUTOLOAD``,
  ``PYTHONWARNINGS``, ``MKL_CBWR``, the set of conftests a full run loads, and the bytes of ``pyproject.toml``,
  ``pytest.ini``, ``setup.cfg`` and ``tox.ini``; any other variable is invisible.
* SCOPE-BOUNDARY[guard] machinery: the memo's pinned-source loader, extended to pytest's assertion-rewriting
  loader, so the bytes that ran are recorded (#688). ruling: the orchestrator, 2026-11-01, review round 1 (a B1
  scope boundary, a B2 candidate).
  revisit: when a source is edited between a test's run and the read of its bytes. A pin reads its bytes after the
  test ran, in the capture child. A file written within ``TRUSTED_AGE_NS`` before the capture began, or after it,
  is stored unsettled; a first-party file whose identity at the end of a test differs from the one it had when the
  capture began makes that test cone-unknown. A file edited and restored, or edited after its test and before the
  test is judged, within one timestamp tick of the stat, is not seen.
* SCOPE-BOUNDARY[guard] machinery: per-statement pins (B2, ruling 33, #689), so an edit selects the tests that ran the
  edited statement. ruling: the plan's ruling 33 (def grain in B1). revisit: B2. Until then over-selection is not
  bounded: whatever the grandparent imports (the root conftest's import closure) is attributed to every test, so an
  edit to the skeleton of any of those modules selects the whole suite, a new file selects every test that probed
  or listed its directory, and agreement with an alone run (containment) passes for every superset.
* SCOPE-BOUNDARY[guard] machinery: a recorder of the reads and children of a process a test starts (#690). ruling: the
  memo's own contract (``traced_memo``, "What no manifest can see"). revisit: a suite that starts processes. A test
  that starts a process, or whose recorded code cannot be pinned, is cone-unknown and always selected.

Importing this package imports nothing of ``orpheus``: the capture workers must start the recorder first.
"""
from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

__all__ = ["capture", "environment_key", "open_store", "select", "status"]


def capture(repo: str | Path, store: str | Path, package: str = "orpheus", memo_root: str | Path | None = None,
            ids: Sequence[str] = (), route: str = "auto", since_ns: int = 0) -> Any:
    from tools.cone_store import driver

    return driver.capture(repo, store, package, memo_root, ids, route, since_ns)


def select(repo: str | Path, store: str | Path, package: str = "orpheus", memo_root: str | Path | None = None,
           ids: Sequence[str] = ()) -> Any:
    from tools.cone_store import selection

    return selection.select(repo, store, package, memo_root, ids)


def status(repo: str | Path, store: str | Path) -> Any:
    from tools.cone_store import selection

    return selection.status(repo, store)


def open_store(path: str | Path) -> Any:
    """The store at ``path`` for reading: ``tests()`` (test id to record with ``route``, ``outcome``, ``unknown``),
    ``pins(test_id)`` (the memo's pin objects, paths relative to the repository) and ``identities()``."""
    from tools.cone_store import store

    return store.ConeStore(path)


def environment_key() -> str:
    from tools.cone_store import environment

    return environment.environment_key()
