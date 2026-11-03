"""The cone store: per test, the first-party code and data the test's verdict depends on, and from it the
tests a working-tree change requires.

    python -m tools.cone_store {capture,select,status,place} --repo R --store S [--package P] [--memo-root M] [ids]

``capture`` runs the tests under a grandparent process that records (the memo's ``Recording``) from before the
root conftest is imported, forks one collecting parent per test module, and forks one child per test; each child's
cone is the recording's pins (``orpheus.numerics.traced_memo.collect_pins``): the defs that ran and the import-time
statements the test could observe (``orpheus.numerics._import_time``: by a name its executed code read, or by
effect) of first-party code (any source under the repository that is not the interpreter's or a distribution's,
``tools/`` included), the outcome of importing the modules it
imported, data, presence and listing pins of paths under the repository, memo children read and distributions. ``select`` lists the collected tests whose cone is stale, new, unknown,
captured under another environment key, or not passing at capture, and a module that failed to collect. Each test
also keeps its OWN cone (what pytest's protocol for that test ran, recorded by a second recording), and ``place``
reads it to place a layer's tests in its concept order (:mod:`tools.cone_store.placement`). The store
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
* SCOPE-BOUNDARY[guard] machinery: a recorder of name lookups in module and class namespaces (the import-time
  statements a test observes are found by the names its executed code could read, ``co_names`` and the identifier
  tokens of its strings, closed over the statements; #689). ruling: the verification-tree plan, ruling 43, after qa's
  42 attack fixtures. revisit: an edit that changes a verdict while every pin of the test is valid. Not seen: a name
  built from data and absent from every string of the code that ran; iteration over a namespace (``vars(module)``,
  a module's ``__dict__``); C code reading a non-dunder class attribute. (A guarded probe of a name new to a module,
  ``getattr(module, name, default)``, is seen: the entry names are stored whole, gate U2.) Over-selection is bounded but not small: the names the grandparent's code read (pytest's and the
  conftests') are entry names of every test, and every effect statement of an imported module is observed.
* SCOPE-BOUNDARY[guard] machinery: a recorder of the reads and children of a process a test starts (#690). ruling: the
  memo's own contract (``traced_memo``, "What no manifest can see"). revisit: a suite that starts processes. A test
  that starts a process, or whose recorded code cannot be pinned, is cone-unknown and always selected.

Importing this package imports nothing of ``orpheus``: the capture workers must start the recorder first.
"""
from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

__all__ = ["capture", "environment_key", "open_store", "place", "select", "status"]


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


def place(store: str | Path, layer_prefix: str, code_root: str, order: Sequence[Sequence[str]], kinds: Any,
          threshold: float = 0.95, declared: Sequence[tuple[str, str]] = ()) -> Any:
    """Place the stored tests under ``layer_prefix`` in ``order`` over the layer at ``code_root`` (the numerics
    layer's surroundings, :data:`tools.cone_store.placement.NUMERICS`, with that root), and measure the
    ``declared`` ``rests_on`` edges (upper test id, target id)."""
    from dataclasses import replace

    from tools.cone_store import placement

    return placement.place(store, layer_prefix, replace(placement.NUMERICS, root=code_root), order, kinds,
                           threshold=threshold, declared=declared)


def open_store(path: str | Path) -> Any:
    """The store at ``path`` for reading: ``tests()`` (test id to record with ``route``, ``outcome``, ``unknown``),
    ``pins(test_id)`` (the memo's pin objects, paths relative to the repository) and ``identities()``."""
    from tools.cone_store import store

    return store.ConeStore(path)


def environment_key() -> str:
    from tools.cone_store import environment

    return environment.environment_key()
