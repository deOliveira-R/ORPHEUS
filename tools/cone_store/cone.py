"""A recorded run read as a cone: the memo's pins, restricted and relocated for the store.

The pins come from :func:`orpheus.numerics.traced_memo.collect_pins`, which returns every refusal beside the pins
(ruling 40). A refusal is a value of the memo's closed :class:`~orpheus.numerics.traced_memo.Refusal` set, read in
one ``match`` (:func:`_unknown`): a dependency no pin can name makes the cone UNKNOWN (an unadmitted process, a file
of the repository the run read and then changed, code that ran without being imported, a distribution file no
``RECORD`` lists), and is dropped when it concerns something that is not the repository's (a temporary file the run
made and read back).

Kept: defs and module skeletons of first-party code (a file under ``<repo>/<package>/`` or ``<repo>/tests/``),
data, presence and listing pins of paths under the repository, memo children read, distributions. The listing of
the directory of every first-party module that ran, and of every ``sys.path`` entry under the repository, is added
beside the recorded ones: the import system's own stats and listings are not recorded (they are how a module is
found), so a new module or package that would shadow an imported one is seen only through those directories.
Dropped by decision: the interpreter, environment, working-directory and memo-source pins (the environment key and
the first-party digests cover them), and everything in the bytecode cache directory (:func:`plan.pycache_dir`).
"""
from __future__ import annotations

import importlib.machinery
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any

from tools.cone_store import pins as pin_rows
from tools.cone_store.plan import Config, pycache_dir
from tools.cone_store.store import FileState


@dataclass(frozen=True)
class Cone:
    """The pins of a run (paths relative to the repository), the state of each file they judge, and why the cone
    is unknown (empty: it is known)."""

    pins: tuple[Any, ...]
    files: dict[str, FileState]
    unknown: tuple[str, ...]

    def to_json(self, memo: ModuleType) -> dict[str, Any]:
        return {"pins": memo.Manifest(frozenset(self.pins)).to_json(),
                "files": {path: [None if ident is None else list(ident), settled] for path, (ident, settled) in self.files.items()},
                "unknown": list(self.unknown)}

    @classmethod
    def from_json(cls, memo: ModuleType, data: dict[str, Any]) -> Cone:
        return cls(tuple(memo.Manifest.from_json(data["pins"]).pins),
                   {path: (None if ident is None else tuple(ident), settled) for path, (ident, settled) in data["files"].items()},
                   tuple(data["unknown"]))


class Reading:
    """Where paths belong, for one repository: the roots of first-party code, and what the store keeps."""

    def __init__(self, repo: Path | str, package: str) -> None:
        self.repo = os.path.realpath(repo)
        self.roots = (f"{self.repo}{os.sep}{package}{os.sep}", f"{self.repo}{os.sep}tests{os.sep}")
        self.cache = pycache_dir(repo) + os.sep

    def first_party(self, real: str) -> bool:
        return real.startswith(self.roots)

    def kept(self, real: str) -> bool:
        """Under the repository and not in the bytecode cache."""
        return real.startswith(self.repo + os.sep) and not real.startswith(self.cache)

    def tool_owned(self, real: str) -> bool:
        """The bytecode cache directory itself or a directory above it (the capture creates them): not the repository's."""
        return (self.cache).startswith(real + os.sep)

    def relative(self, real: str) -> str:
        return os.path.relpath(real, self.repo)


def cone_of(recording: Any, store: Any, config: Config, memo: ModuleType, boot: ModuleType) -> Cone:
    """The cone of a recording that has been stopped. ``store`` is the memo's ``Store`` (children live there).
    A file's identity is settled when it was last written more than ``TRUSTED_AGE_NS`` before ``config.since_ns``,
    the capture's start, through the memo's own :class:`Witness`."""
    reading = Reading(config.repo, config.package)
    reals = {memo._real(name) for name, _ in recording.codes if not name.startswith("<")}
    code_files = {real for real in reals if reading.first_party(real)}
    data_files = {real for opened in recording.opened if reading.kept(real := memo._real(opened))}
    states: dict[str, FileState] = {}
    for real in code_files | data_files:  # each file's identity is read BEFORE its pins read its bytes
        witness = boot.Witness(reference_ns=config.since_ns)
        found = witness.stat(real)
        states[real] = (None if found is None else tuple(boot.FileIdentity.of(found)), witness.settled)

    collected = memo.collect_pins(recording, store)
    by_pin_path = {memo._pin_path(real): real for real in code_files}

    def real_of(path: str) -> str:
        return path if os.path.isabs(path) else by_pin_path.get(path, "")

    def to_store(pin: pin_rows.PathPin) -> str | None:
        match pin:
            case memo.DefPin() | memo.ModulePin():
                real = real_of(pin.path)
                return reading.relative(real) if reading.first_party(real) else None
            case _:
                real = memo._real(pin.path)
                return reading.relative(real) if reading.kept(real) and not reading.tool_owned(real) else None

    kept: list[Any] = []
    for pin in collected.pins:
        match pin:
            case memo.InterpreterPin() | memo.EnvironmentPin() | memo.WorkingDirectoryPin() | memo.MemoPin():
                continue
            case _:
                if (moved := pin_rows.relocated(pin, to_store)) is not None:
                    kept.append(moved)
    kept += [pin_rows.relocated(pin, to_store) for pin in shadow_pins(code_files, reading, memo, boot)]
    pins = tuple(sorted({pin for pin in kept if pin is not None}, key=pin_rows.encode))
    judged = {pin_rows.file_of(pin) for pin in pins}
    files = {reading.relative(real): state for real, state in states.items() if reading.relative(real) in judged}
    unknown = [reason for refusal in collected.refusals
               if (reason := _unknown(refusal, memo, reading, real_of)) is not None]
    unknown += _changed_during_capture(files, config)
    return Cone(pins, files, tuple(unknown))


def shadow_pins(code_files: set[str], reading: Reading, memo: ModuleType, boot: ModuleType) -> list[Any]:
    """What a new file would have to be, to change which file an import resolved to: for each first-party module that
    ran, the presence of every name that would be found before it. In its directory, a package directory or an
    extension module of the same name; and, for a top-level name, the same names in the other ``sys.path`` entries
    under the repository. The import system's own stats are not recorded (they are how a module is found), so a
    new module or package that shadows an imported one would otherwise be seen by no pin. Pins of presence, not of
    the directory's listing: a listing changes when ``__pycache__`` appears, which shadows nothing."""
    entries = [real for entry in sys.path if os.path.isdir(real := os.path.realpath(entry or os.getcwd())) and reading.kept(real + os.sep)]
    suffixes = ("", *importlib.machinery.EXTENSION_SUFFIXES, ".py")
    candidates: set[str] = set()
    for real in code_files:
        root = next((entry for entry in entries if real.startswith(entry + os.sep)), None)
        if root is None:
            continue
        parts = os.path.relpath(real, root).split(os.sep)
        parts[-1] = parts[-1].removesuffix(".py")
        if parts[-1] == "__init__":
            parts.pop()
        for depth, name in enumerate(parts):
            directories = [os.path.join(root, *parts[:depth])] + ([e for e in entries if e != root] if depth == 0 else [])
            candidates |= {os.path.join(directory, name + suffix) for directory in directories for suffix in suffixes}
    return [memo.PresencePin(path, boot.presence(path)) for path in sorted(candidates)]


def _unknown(refusal: Any, memo: ModuleType, reading: Reading, real_of: Any) -> str | None:
    """Why a cone is unknown, from one refusal of ``collect_pins``; ``None`` for a refusal that is not about the
    repository (a temporary file the run wrote, a path outside the roots)."""
    match refusal:
        case memo.UnadmittedProcess(program):
            return f"process started: {program}"
        case memo.ChangedByRun(path):
            return f"{reading.relative(memo._real(path))} was read and then changed by the run" \
                if reading.kept(memo._real(path)) else None
        case memo.NoRecordedBytes(path) | memo.SourceGone(path):
            return str(refusal) if reading.first_party(memo._real(path)) else None
        case memo.NoDefAtLine(path):
            return str(refusal) if reading.first_party(real_of(path)) else None
        case memo.AmbiguousOwner() | memo.UnownedSitePackage() | memo.UninstalledDistribution():
            return str(refusal)
    raise TypeError(f"cone_store: a refusal of a kind it does not know: {refusal!r}")


def _changed_during_capture(files: dict[str, FileState], config: Config) -> list[str]:
    """A file that is not settled against the capture's start (written within ``TRUSTED_AGE_NS`` before it, or
    after it) and whose identity is not the one it had when the capture began was written during the capture: the
    pins read its bytes after the test ran, so they may pin bytes that did not run."""
    return [f"{path} changed during the capture" for path, (ident, settled) in sorted(files.items())
            if not settled and ident is not None and config.recent.get(path) != list(ident)]
