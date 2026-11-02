"""The SQLite store: one row per (test, pin), the pins shared between tests, the file identities shared between pins.

Layout. ``tests`` holds a test id with the environment key it was captured under, its capture route (``fork`` or
``fresh``), its outcome and its cone-unknown reasons (a JSON list of strings, never a bare flag). ``pins`` holds
each distinct pin once, as the memo's own row (``[kind, *fields]``, :mod:`tools.cone_store.pins`), with the file
it judges when it judges one. ``files`` holds, per file, the identity it had at capture and whether that identity
was settled (older than the memo's ``TRUSTED_AGE_NS`` before the capture began). ``cones`` joins tests to pins.
The view ``rows`` is the flat reading, one row per test and pin. Paths are relative to the repository.

A file's identity is a fact about the FILE: it is written when any test's capture reads the file, and it is
trusted (stage one of ``select``) only while every pin on the file was judged against those bytes. When a capture
finds a file with another identity than the stored one, or the stored one was not settled, the file is stored
unsettled, so the next ``select`` judges every pin on it (stage two) and refreshes the identity when all are clean.
"""
from __future__ import annotations

import json
import sqlite3
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from orpheus.numerics import traced_memo as memo

from tools.cone_store import pins as pin_rows

_SCHEMA = """
CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY, value TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS tests(
    test_id TEXT PRIMARY KEY, env_key TEXT NOT NULL, route TEXT NOT NULL, outcome TEXT NOT NULL,
    unknown TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS files(path TEXT PRIMARY KEY, ident TEXT, settled INTEGER NOT NULL);
CREATE TABLE IF NOT EXISTS pins(
    pin_id INTEGER PRIMARY KEY, row TEXT NOT NULL UNIQUE, file TEXT REFERENCES files(path));
CREATE INDEX IF NOT EXISTS pins_by_file ON pins(file);
CREATE TABLE IF NOT EXISTS cones(
    test_id TEXT NOT NULL REFERENCES tests(test_id) ON DELETE CASCADE, pin_id INTEGER NOT NULL REFERENCES pins(pin_id),
    PRIMARY KEY(test_id, pin_id)) WITHOUT ROWID;
CREATE INDEX IF NOT EXISTS cones_by_pin ON cones(pin_id);
CREATE VIEW IF NOT EXISTS rows AS
    SELECT c.test_id AS test_id, t.env_key AS env_key, p.row AS pin, p.file AS file, f.ident AS ident, f.settled AS settled
    FROM cones c JOIN tests t USING(test_id) JOIN pins p USING(pin_id) LEFT JOIN files f ON f.path = p.file;
"""

#: A file's identity at capture (a ``FileIdentity`` as a tuple, ``None`` when the file was missing) and whether it was settled.
FileState = tuple[tuple[int, ...] | None, bool]


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

    def record(self, test: StoredTest, pins: Iterable[memo.Pin], files: Mapping[str, FileState]) -> None:
        """Replace everything the store holds about ``test.test_id``; ``files`` is the state of each file its pins judge."""
        with self.db:
            self.db.execute("DELETE FROM cones WHERE test_id = ?", (test.test_id,))
            self.db.execute("INSERT OR REPLACE INTO tests VALUES (?, ?, ?, ?, ?)",
                            (test.test_id, test.env_key, test.route, test.outcome, json.dumps(list(test.unknown))))
            pins = list(pins)
            for path in {f for pin in pins if (f := pin_rows.file_of(pin)) is not None}:
                self._note_file(path, *files.get(path, (None, False)))
            for pin in pins:
                row = pin_rows.encode(pin)
                self.db.execute("INSERT OR IGNORE INTO pins(row, file) VALUES (?, ?)", (row, pin_rows.file_of(pin)))
                (pin_id,) = self.db.execute("SELECT pin_id FROM pins WHERE row = ?", (row,)).fetchone()
                self.db.execute("INSERT OR IGNORE INTO cones VALUES (?, ?)", (test.test_id, pin_id))

    def _note_file(self, path: str, ident: tuple[int, ...] | None, settled: bool) -> None:
        old = self.db.execute("SELECT ident, settled FROM files WHERE path = ?", (path,)).fetchone()
        encoded = None if ident is None else json.dumps(list(ident))
        trusted = settled and (old is None or (bool(old[1]) and old[0] == encoded))
        self.db.execute("INSERT INTO files VALUES (?, ?, ?) ON CONFLICT(path) DO UPDATE SET ident = excluded.ident, "
                        "settled = excluded.settled", (path, encoded, int(trusted)))

    def delete(self, test_ids: Iterable[str]) -> int:
        """Forget these tests (their pins stay until :meth:`collect_garbage`); the number deleted."""
        deleted = 0
        with self.db:
            for test_id in test_ids:
                self.db.execute("DELETE FROM cones WHERE test_id = ?", (test_id,))
                deleted += self.db.execute("DELETE FROM tests WHERE test_id = ?", (test_id,)).rowcount
        return deleted

    def collect_garbage(self) -> None:
        with self.db:
            self.db.execute("DELETE FROM pins WHERE pin_id NOT IN (SELECT pin_id FROM cones)")
            self.db.execute("DELETE FROM files WHERE path NOT IN (SELECT file FROM pins WHERE file IS NOT NULL)")

    def refresh(self, path: str, ident: tuple[int, ...], settled: bool) -> None:
        """A file's identity after every pin on it was judged clean against its current bytes (stage one's cache)."""
        with self.db:
            self.db.execute("UPDATE files SET ident = ?, settled = ? WHERE path = ?",
                            (json.dumps(list(ident)), int(settled), path))

    # ── reading ──

    def tests(self) -> dict[str, StoredTest]:
        return {row[0]: StoredTest(row[0], row[1], row[2], row[3], tuple(json.loads(row[4])))
                for row in self.db.execute("SELECT test_id, env_key, route, outcome, unknown FROM tests")}

    def pins(self, test_id: str) -> list[memo.Pin]:
        """The pins of one test, as the memo's own objects with paths relative to the repository, in row order."""
        return [pin_rows.decode(row) for (row,) in self.db.execute(
            "SELECT p.row FROM cones c JOIN pins p USING(pin_id) WHERE c.test_id = ? ORDER BY p.row", (test_id,))]

    def identities(self) -> dict[str, FileState]:
        """Each file's identity at capture and whether it was settled."""
        return {path: (None if ident is None else tuple(json.loads(ident)), bool(settled))
                for path, ident, settled in self.db.execute("SELECT path, ident, settled FROM files")}

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
