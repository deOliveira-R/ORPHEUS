"""The SQLite store: one row per (test, pin), the pins shared between tests, the file identities shared between pins.

Layout. ``tests`` holds a test id with the environment key it was captured under, its capture route (``fork`` or
``fresh``), its outcome and its cone-unknown reasons (a JSON list of strings, never a bare flag). ``pins`` holds
each distinct pin once, as the memo's own row (``[kind, *fields]``, :mod:`tools.cone_store.pins`), with the file
it judges when it judges one. ``files`` holds, per file, the identity it had at capture, whether that identity
was settled (older than the memo's ``TRUSTED_AGE_NS`` before the capture began) and, for a source, the digest of its
bytes then, which a pin naming the file's bytes elsewhere (an ``ImportOutcomePin``'s row) is compared with while
stage one vouches for the file. ``cones`` joins tests to pins.
The view ``rows`` is the flat reading, one row per test and pin. Paths are relative to the repository.

A test has two cones, two relations over the same pins (B3). ``cones`` is its SELECTION cone: everything its process
ran, the grandparent's conftest imports and its module's collection included (ruling 38), which is what ``select``
judges. ``own_cones`` is its OWN cone: what pytest's protocol for that one test ran (setup, call, teardown), read
from a second recording started at the protocol, so a def the conftests also ran is in it when the test runs it
too; placement reads only this one. ``own_reads`` holds what that run READ (B5, :mod:`tools.cone_store.reads`): one
row per binding it reached, by the repository path of the module whose statement binds it, the bound name (empty for a
read of the module object) and the channel it was reached by; placement reads these rows, never the name sets a pin
keeps in the capturing checkout's cache. ``shares`` holds a cone common to many tests under a scope name: the
grandparent's (:data:`GRANDPARENT`), the floor every test of a capture inherits.

A file's identity is a fact about the FILE: it is written when any test's capture reads the file, and it is
trusted (stage one of ``select``) only while every pin on the file was judged against those bytes. When a capture
finds a file with another identity than the stored one, or the stored one was not settled, the file is stored
unsettled, so the next ``select`` judges every pin on it (stage two) and refreshes the identity when all are clean.
The refresh takes the digest of the bytes it judged with the identity: an identity vouches for the bytes it was
taken on, never for the bytes of the capture, so a pin that remembers the capture's bytes (an import set's row)
still sees an edit that every pin on the file was clean against (qa's final pass on B2, F1).
"""
from __future__ import annotations

import json
import sqlite3
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, NamedTuple

from orpheus.numerics import traced_memo as memo

from tools.cone_store import pins as pin_rows
from tools.cone_store.reads import Read

_SCHEMA = """
CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY, value TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS tests(
    test_id TEXT PRIMARY KEY, env_key TEXT NOT NULL, route TEXT NOT NULL, outcome TEXT NOT NULL,
    unknown TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS files(path TEXT PRIMARY KEY, ident TEXT, settled INTEGER NOT NULL, digest TEXT);
CREATE TABLE IF NOT EXISTS pins(
    pin_id INTEGER PRIMARY KEY, row TEXT NOT NULL UNIQUE, file TEXT REFERENCES files(path));
CREATE INDEX IF NOT EXISTS pins_by_file ON pins(file);
CREATE TABLE IF NOT EXISTS cones(
    test_id TEXT NOT NULL REFERENCES tests(test_id) ON DELETE CASCADE, pin_id INTEGER NOT NULL REFERENCES pins(pin_id),
    PRIMARY KEY(test_id, pin_id)) WITHOUT ROWID;
CREATE INDEX IF NOT EXISTS cones_by_pin ON cones(pin_id);
CREATE TABLE IF NOT EXISTS own_cones(
    test_id TEXT NOT NULL REFERENCES tests(test_id) ON DELETE CASCADE, pin_id INTEGER NOT NULL REFERENCES pins(pin_id),
    PRIMARY KEY(test_id, pin_id)) WITHOUT ROWID;
CREATE TABLE IF NOT EXISTS own_reads(
    test_id TEXT NOT NULL REFERENCES tests(test_id) ON DELETE CASCADE, path TEXT NOT NULL, name TEXT NOT NULL,
    channel TEXT NOT NULL, PRIMARY KEY(test_id, path, name, channel)) WITHOUT ROWID;
CREATE TABLE IF NOT EXISTS shares(
    scope TEXT NOT NULL, pin_id INTEGER NOT NULL REFERENCES pins(pin_id), PRIMARY KEY(scope, pin_id)) WITHOUT ROWID;
CREATE VIEW IF NOT EXISTS rows AS
    SELECT c.test_id AS test_id, t.env_key AS env_key, p.row AS pin, p.file AS file, f.ident AS ident, f.settled AS settled
    FROM cones c JOIN tests t USING(test_id) JOIN pins p USING(pin_id) LEFT JOIN files f ON f.path = p.file;
"""

#: The form of the rows this code writes and reads: 4 since B5, whose ``own_reads`` hold what each test's own run read
#: (3 was B3, whose ``own_cones`` and ``shares`` hold each test's own cone and the grandparent's; 2 was B2, #689: ``StatementsPin`` and ``ImportOutcomePin`` replaced the module
#: skeleton's ``ModulePin``, and ``files`` kept each source's digest). A store of another form holds rows this code
#: cannot read as it means them (a B2 store has no own cones, and placement would read none), so it is emptied when
#: opened, and every test is then new (selected, and captured again).
FORMAT = 4

#: The scope of the grandparent's cone in ``shares``: what the conftests' import ran, inherited by every test.
GRANDPARENT = "<grandparent>"

class FileState(NamedTuple):
    """A file as a capture, or a clean ``select``, last saw it: its identity (a ``FileIdentity`` as a tuple, ``None``
    when the file was missing), whether that identity was settled, and the digest of the bytes it was taken on
    (``None`` for a file only data, presence or listing pins judge, whose bytes no other pin names)."""

    ident: tuple[int, ...] | None
    settled: bool
    digest: str | None = None


_MISSING = FileState(None, False)


def stored_format(path: str | Path) -> int | None:
    """The form of the store at ``path`` without opening it as a :class:`ConeStore` (which empties a store of
    another form); ``None`` when there is no store there."""
    if not Path(path).exists():
        return None
    db = sqlite3.connect(f"file:{Path(path)}?mode=ro", uri=True)
    try:
        found = db.execute("SELECT value FROM meta WHERE key = 'format'").fetchone()
    except sqlite3.OperationalError:
        return None
    finally:
        db.close()
    return None if found is None else json.loads(found[0])


@dataclass(frozen=True)
class StoredTest:
    test_id: str
    env_key: str
    route: str
    outcome: str
    unknown: tuple[str, ...]


class ConeStore:
    """The store at ``path`` (created, with its directory, when absent)."""

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.db = sqlite3.connect(self.path)
        self.db.executescript(_SCHEMA)
        self.db.execute("PRAGMA foreign_keys = ON")
        if self.meta("format") != FORMAT:  # another form's tables are dropped, not emptied: their columns may differ
            with self.db:
                self.db.execute("DROP VIEW IF EXISTS rows")
                for table in ("own_reads", "own_cones", "shares", "cones", "pins", "files", "tests", "meta"):
                    self.db.execute(f"DROP TABLE IF EXISTS {table}")
            self.db.executescript(_SCHEMA)
            self.set_meta("format", FORMAT)

    def close(self) -> None:
        self.db.close()

    def __enter__(self) -> ConeStore:
        return self

    def __exit__(self, *_exc: object) -> None:
        self.close()

    # ── meta ──

    def meta(self, key: str, default: Any = None) -> Any:
        found = self.db.execute("SELECT value FROM meta WHERE key = ?", (key,)).fetchone()
        return default if found is None else json.loads(found[0])

    def set_meta(self, key: str, value: Any) -> None:
        with self.db:
            self.db.execute("INSERT INTO meta VALUES (?, ?) ON CONFLICT(key) DO UPDATE SET value = excluded.value",
                            (key, json.dumps(value)))

    # ── writing ──

    def record(self, test: StoredTest, pins: Iterable[memo.Pin], files: Mapping[str, FileState],
               own: Iterable[memo.Pin] = (), own_files: Mapping[str, FileState] | None = None,
               reads: Iterable[Read] = ()) -> None:
        """Replace everything the store holds about ``test.test_id``: its selection cone ``pins``, its own cone
        ``own`` and what its own run read, ``reads``. ``files`` is the state of each file the selection cone judges,
        ``own_files`` of each file the own cone judges; where both name a file the selection's state is kept (one
        capture, one identity per file)."""
        with self.db:
            self.db.execute("DELETE FROM cones WHERE test_id = ?", (test.test_id,))
            self.db.execute("DELETE FROM own_cones WHERE test_id = ?", (test.test_id,))
            self.db.execute("DELETE FROM own_reads WHERE test_id = ?", (test.test_id,))
            self.db.execute("INSERT OR REPLACE INTO tests VALUES (?, ?, ?, ?, ?)",
                            (test.test_id, test.env_key, test.route, test.outcome, json.dumps(list(test.unknown))))
            states = {**(own_files or {}), **files}
            for table, cone in (("cones", pins), ("own_cones", own)):
                for pin_id in self._pin_ids(cone, states):
                    self.db.execute(f"INSERT OR IGNORE INTO {table} VALUES (?, ?)", (test.test_id, pin_id))
            self.db.executemany("INSERT OR IGNORE INTO own_reads VALUES (?, ?, ?, ?)",
                                [(test.test_id, *read.to_row()) for read in reads])

    def record_share(self, scope: str, pins: Iterable[memo.Pin], files: Mapping[str, FileState]) -> None:
        """Replace the cone held under ``scope`` (:data:`GRANDPARENT`)."""
        with self.db:
            self.db.execute("DELETE FROM shares WHERE scope = ?", (scope,))
            for pin_id in self._pin_ids(pins, files):
                self.db.execute("INSERT OR IGNORE INTO shares VALUES (?, ?)", (scope, pin_id))

    def _pin_ids(self, pins: Iterable[memo.Pin], files: Mapping[str, FileState]) -> list[int]:
        """Each pin's id, the pin and the state of the file it judges written first (inside the caller's transaction)."""
        pins = list(pins)
        for path in {f for pin in pins if (f := pin_rows.file_of(pin)) is not None}:
            self._note_file(path, files.get(path, _MISSING))
        ids = []
        for pin in pins:
            row = pin_rows.encode(pin)
            self.db.execute("INSERT OR IGNORE INTO pins(row, file) VALUES (?, ?)", (row, pin_rows.file_of(pin)))
            ids.append(self.db.execute("SELECT pin_id FROM pins WHERE row = ?", (row,)).fetchone()[0])
        return ids

    def _note_file(self, path: str, state: FileState) -> None:
        old = self.db.execute("SELECT ident, settled FROM files WHERE path = ?", (path,)).fetchone()
        encoded = None if state.ident is None else json.dumps(list(state.ident))
        trusted = state.settled and (old is None or (bool(old[1]) and old[0] == encoded))
        self.db.execute("INSERT INTO files VALUES (?, ?, ?, ?) ON CONFLICT(path) DO UPDATE SET ident = excluded.ident, "
                        "settled = excluded.settled, digest = excluded.digest", (path, encoded, int(trusted), state.digest))

    def delete(self, test_ids: Iterable[str]) -> int:
        """Forget these tests (their pins stay until :meth:`collect_garbage`); the number deleted."""
        deleted = 0
        with self.db:
            for test_id in test_ids:
                self.db.execute("DELETE FROM cones WHERE test_id = ?", (test_id,))
                self.db.execute("DELETE FROM own_cones WHERE test_id = ?", (test_id,))
                self.db.execute("DELETE FROM own_reads WHERE test_id = ?", (test_id,))
                deleted += self.db.execute("DELETE FROM tests WHERE test_id = ?", (test_id,)).rowcount
        return deleted

    def collect_garbage(self) -> None:
        with self.db:
            self.db.execute("DELETE FROM pins WHERE pin_id NOT IN "
                            "(SELECT pin_id FROM cones UNION SELECT pin_id FROM own_cones UNION SELECT pin_id FROM shares)")
            self.db.execute("DELETE FROM files WHERE path NOT IN (SELECT file FROM pins WHERE file IS NOT NULL)")

    def refresh(self, path: str, state: FileState) -> None:
        """A file's state after every pin on it was judged clean against its current bytes (stage one's cache)."""
        assert state.ident is not None  # narrows: a refresh follows a stat that found the file
        with self.db:
            self.db.execute("UPDATE files SET ident = ?, settled = ?, digest = ? WHERE path = ?",
                            (json.dumps(list(state.ident)), int(state.settled), state.digest, path))

    # ── reading ──

    def tests(self) -> dict[str, StoredTest]:
        return {row[0]: StoredTest(row[0], row[1], row[2], row[3], tuple(json.loads(row[4])))
                for row in self.db.execute("SELECT test_id, env_key, route, outcome, unknown FROM tests")}

    def pins(self, test_id: str) -> list[memo.Pin]:
        """The pins of one test's selection cone, as the memo's own objects with paths relative to the repository,
        in row order."""
        return [pin_rows.decode(row) for (row,) in self.db.execute(
            "SELECT p.row FROM cones c JOIN pins p USING(pin_id) WHERE c.test_id = ? ORDER BY p.row", (test_id,))]

    def own_pins(self, test_id: str) -> list[memo.Pin]:
        """The pins of one test's own cone (what its protocol ran), in row order."""
        return [pin_rows.decode(row) for (row,) in self.db.execute(
            "SELECT p.row FROM own_cones c JOIN pins p USING(pin_id) WHERE c.test_id = ? ORDER BY p.row", (test_id,))]

    def own_reads(self, test_id: str) -> list[Read]:
        """What one test's own run read (:class:`~tools.cone_store.reads.Read`), sorted."""
        return [Read.from_row(row) for row in self.db.execute(
            "SELECT path, name, channel FROM own_reads WHERE test_id = ? ORDER BY path, name, channel", (test_id,))]

    def share(self, scope: str) -> list[memo.Pin] | None:
        """The cone held under ``scope``, in row order; ``None`` when none was stored."""
        if self.db.execute("SELECT 1 FROM shares WHERE scope = ? LIMIT 1", (scope,)).fetchone() is None:
            return None
        return [pin_rows.decode(row) for (row,) in self.db.execute(
            "SELECT p.row FROM shares s JOIN pins p USING(pin_id) WHERE s.scope = ? ORDER BY p.row", (scope,))]

    def identities(self) -> dict[str, FileState]:
        """Each file's state: its identity at capture (or at the last clean ``select``), whether it was settled, and the digest of the bytes it was taken on."""
        return {path: FileState(None if ident is None else tuple(json.loads(ident)), bool(settled), digest)
                for path, ident, settled, digest in self.db.execute("SELECT path, ident, settled, digest FROM files")}

    def cones(self) -> Iterator[tuple[str, int]]:
        """Every (test id, pin id)."""
        yield from self.db.execute("SELECT test_id, pin_id FROM cones")

    def pin_table(self, pin_ids: Iterable[int] | None = None) -> dict[int, tuple[memo.Pin, str | None]]:
        """Pin id to (pin, the file it judges): every pin, or the given ones."""
        wanted = None if pin_ids is None else set(pin_ids)
        return {pin_id: (pin_rows.decode(row), file) for pin_id, row, file in self.db.execute(
            "SELECT pin_id, row, file FROM pins") if wanted is None or pin_id in wanted}

    def count(self, table: str) -> int:
        return self.db.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
