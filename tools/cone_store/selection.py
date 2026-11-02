"""``select`` and ``status``: which tests a working-tree change requires, from the store alone.

A test is selected when it is new (not in the store), when the store holds it under another environment key, when
its cone is unknown, when it did not pass at capture (its cone stops at the failure), or when any pin of its cone is
stale; a module that failed to collect is selected by its path. A pin is judged in two stages. Stage one stats the
pinned FILE: a file whose identity (device, inode, size, modification and change times) is unchanged since
capture, and was settled then (older than the memo's ``TRUSTED_AGE_NS`` before the capture began), is fresh without
being read. Stage two judges every pin on the file with the memo's own ``Pin.stale``, against the memo's own
``_Checkout``. After a stage two that found every pin on the file clean the file's stored identity is refreshed, so
the next ``select`` stops at stage one; a file with a stale pin is judged again each time.
"""
from __future__ import annotations

import contextlib
import sys
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from orpheus.numerics import traced_memo
from orpheus.numerics._traced_memo_boot import FileIdentity, Witness

from tools.cone_store import pins as pin_rows
from tools.cone_store.collect import collect_ids
from tools.cone_store.environment import environment_key
from tools.cone_store.store import ConeStore


@dataclass(frozen=True)
class Selection:
    """The node ids to run, and for each the reasons."""

    ids: tuple[str, ...]
    reasons: dict[str, tuple[str, ...]] = field(default_factory=dict)

    def __iter__(self) -> Iterator[str]:
        return iter(self.ids)

    def __contains__(self, test_id: object) -> bool:
        return test_id in self.ids

    def __len__(self) -> int:
        return len(self.ids)


@contextlib.contextmanager
def _importable(repo: Path) -> Iterator[None]:
    """The memo resolves a pin's path through ``sys.path`` (a ``DefPin`` is relative to an entry of the process that
    ran it), and a memo child's manifest names paths relative to the process that generated it."""
    sys.path.insert(0, str(repo))
    try:
        yield
    finally:
        sys.path.remove(str(repo))


def _absolute(repo: Path, pin: Any) -> Any:
    """The pin as the memo judges it: a pin that names a path other than by a ``sys.path`` entry (data, presence,
    listing) names it absolutely."""
    return pin if isinstance(pin, (traced_memo.DefPin, traced_memo.ModulePin)) else \
        pin_rows.relocated(pin, lambda p: str(repo / p.path))


def _stale_pins(database: ConeStore, repo: Path, memo_store: Any, wanted: set[int]) -> dict[int, str]:
    """The reason each stale pin among ``wanted`` is stale (fresh pins are absent)."""
    checkout = traced_memo._Checkout(memo_store)
    table = database.pin_table()
    identities = database.identities()
    by_file: dict[str, list[int]] = {}
    stale: dict[int, str] = {}

    def judge(pin_id: int) -> None:
        reason = _absolute(repo, table[pin_id][0]).stale(checkout)
        if reason is not None:
            stale[pin_id] = reason

    for pin_id in wanted:
        file = table[pin_id][1]
        if file is None:
            judge(pin_id)
        else:
            by_file.setdefault(file, [])
    for pin_id, (_, file) in table.items():
        if file in by_file:
            by_file[file].append(pin_id)
    for file, pin_ids in by_file.items():
        ident, settled = identities.get(file, (None, False))
        witness = Witness()
        found = witness.stat(str(repo / file))  # before the read it guards
        if ident is not None and settled and found is not None and tuple(FileIdentity.of(found)) == ident:
            continue  # stage one: the file is the one captured
        for pin_id in pin_ids:
            judge(pin_id)
        if found is not None and not any(pin_id in stale for pin_id in pin_ids):
            database.refresh(file, tuple(FileIdentity.of(found)), witness.settled)  # the stage-one cache
    return {pin_id: reason for pin_id, reason in stale.items() if pin_id in wanted}


def select(repo: str | Path, store: str | Path, package: str = "orpheus", memo_root: str | Path | None = None,
           ids: Sequence[str] = ()) -> Selection:
    """The tests among ``ids`` (default: all collected) that the working tree requires."""
    repo = Path(repo).resolve()
    collected = collect_ids(repo, package, ids)
    key = environment_key(repo, collected.conftests)
    reasons: dict[str, list[str]] = {module: [f"collection error in {module}"] for module in collected.errors}
    with ConeStore(store) as database, _importable(repo):
        stored = database.tests()
        memo_store = traced_memo.Store(Path(memo_root).resolve()) if memo_root is not None else traced_memo.Store.current()
        candidates: set[str] = set()
        for test_id in collected.ids:
            test = stored.get(test_id)
            if test is None:
                reasons[test_id] = ["new: not in the store"]
            elif test.env_key != key:
                reasons[test_id] = ["environment key differs from the capture's"]
            elif test.unknown:
                reasons[test_id] = [f"cone unknown: {reason}" for reason in test.unknown]
            elif test.outcome != "passed":
                reasons[test_id] = [f"did not pass at capture ({test.outcome})"]
            else:
                candidates.add(test_id)
        if candidates:
            pairs = [(t, p) for t, p in database.cones() if t in candidates]
            stale = _stale_pins(database, repo, memo_store, {p for _, p in pairs})
            for test_id, pin_id in pairs:
                if pin_id in stale:
                    reasons.setdefault(test_id, []).append(stale[pin_id])
    chosen = tuple([*collected.errors, *(test_id for test_id in collected.ids if test_id in reasons)])
    return Selection(chosen, {test_id: tuple(reasons[test_id]) for test_id in chosen})


@dataclass(frozen=True)
class Status:
    path: str
    tests: int
    routes: dict[str, int]
    unknown: dict[str, tuple[str, ...]]
    unknown_tests: int
    collection_errors: int
    not_passed: int
    other_environment: int
    pins: int
    rows: int
    grandparent_pins: int | None
    pruned: int | None

    def __str__(self) -> str:
        routes = ", ".join(f"{route} {count}" for route, count in sorted(self.routes.items())) or "none"
        lines = [f"store {self.path}",
                 f"tests {self.tests} (routes: {routes}); {self.unknown_tests} cone-unknown, {self.collection_errors} "
                 f"collection errors, {self.not_passed} not passed, {self.other_environment} under another environment key",
                 f"pins {self.pins} distinct, {self.rows} rows; grandparent pins {self.grandparent_pins}; "
                 f"pruned by the last capture {self.pruned}"]
        lines += [f"  unknown {test_id}: {'; '.join(reasons)}" for test_id, reasons in sorted(self.unknown.items())]
        return "\n".join(lines)


def status(repo: str | Path, store: str | Path) -> Status:
    """What the store holds (one collection, for the conftest set the environment key carries)."""
    key = environment_key(repo, collect_ids(Path(repo).resolve(), "orpheus").conftests)
    with ConeStore(store) as database:
        tests = database.tests()
        routes: dict[str, int] = {}
        for test in tests.values():
            routes[test.route] = routes.get(test.route, 0) + 1
        unknown = {t.test_id: t.unknown for t in tests.values() if t.unknown}
        errors = sum(1 for t in tests.values() if t.outcome == "collection-error")
        return Status(
            str(database.path), len(tests), routes, unknown, len(unknown), errors,
            sum(1 for t in tests.values() if t.outcome != "passed"), sum(1 for t in tests.values() if t.env_key != key),
            database.count("pins"), database.count("cones"), database.meta("grandparent_pins"), database.meta("pruned"))
