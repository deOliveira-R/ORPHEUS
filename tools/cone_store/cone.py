"""A recorded run read as a cone: the memo's pins, restricted and relocated for the store.

The pins come from :func:`orpheus.numerics.traced_memo.collect_pins`, which returns every refusal beside the pins
(ruling 40). A refusal is a value of the memo's closed :class:`~orpheus.numerics.traced_memo.PinRefusal` set, read in
one ``match`` (:func:`_unknown`): a dependency no pin can name makes the cone UNKNOWN (an unadmitted process, a file
of the repository the run read and then changed, code that ran without being imported, a distribution file no
``RECORD`` lists, an imported module whose import cannot be measured again), and is dropped when it concerns something that is not the repository's (a temporary file the run
made and read back).

Kept: defs and the observed import-time statements of first-party code (a source under the repository that is not
the interpreter's or a distribution's; one projection per module, on the run's entry names), the run's import outcome, data, presence and listing pins of paths under the repository, memo children read from the memo's store, distributions. A memo entry read from another store (a test's own ``cache_root``) arrives as the pins its entry was validated on,
children flattened in and placed where the read found them (``collect_pins``), kept and dropped like the run's own;
those pins are in the selection cone and kept out of the run's own cone (:class:`ReadCones`) unless the run ran them too,
since another process ran them. The listing of
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
from typing import Any, NamedTuple

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
                "files": {path: [None if state.ident is None else list(state.ident), state.settled, state.digest]
                          for path, state in self.files.items()},
                "unknown": list(self.unknown)}

    @classmethod
    def from_json(cls, memo: ModuleType, data: dict[str, Any]) -> Cone:
        return cls(tuple(memo.Manifest.from_json(data["pins"]).pins),
                   {path: FileState(None if ident is None else tuple(ident), settled, digest)
                    for path, (ident, settled, digest) in data["files"].items()},
                   tuple(data["unknown"]))


class Reading:
    """Where paths belong, for one repository: which code is first-party, and what the store keeps."""

    def __init__(self, repo: Path | str) -> None:
        self.repo = os.path.realpath(repo)
        self.cache = pycache_dir(repo) + os.sep
        self.environment = tuple({os.path.realpath(p) + os.sep for p in (sys.prefix, sys.base_prefix)})

    def first_party(self, real: str) -> bool:
        """Code under the repository that is not the interpreter's or an installed distribution's (a virtual
        environment inside the checkout) and not the bytecode cache: the package, the tests, and every other source a
        test can run (``tools/``, whose tests are under ``tests/gates/tools``; B1 kept only the package and the tests,
        so an edit to ``tools/`` selected none of the tests that ran it)."""
        return self.kept(real) and not real.startswith(self.environment)

    def kept(self, real: str) -> bool:
        """Under the repository and not in the bytecode cache."""
        return real.startswith(self.repo + os.sep) and not real.startswith(self.cache)

    def tool_owned(self, real: str) -> bool:
        """The bytecode cache directory itself or a directory above it (the capture creates them): not the repository's."""
        return (self.cache).startswith(real + os.sep)

    def relative(self, real: str) -> str:
        return os.path.relpath(real, self.repo)


class ReadCones(NamedTuple):
    """One recording read as two cones. ``ran`` holds what the run did: the code it ran, what it read, probed and
    listed. ``selection`` holds also what the run depends on only through a memo entry it read from another store
    (``CollectedPins.through_reads``), which the generating process ran: ``select`` judges a test by its selection
    cone, ``place`` homes it by its own, which is the ``ran`` cone of the own recording. Both carry the same unknown
    reasons."""

    ran: Cone
    selection: Cone


def cone_of(recording: Any, store: Any, config: Config, memo: ModuleType, boot: ModuleType) -> ReadCones:
    """The cones of a recording that has been stopped. ``store`` is the memo's ``Store`` (children live there).
    A file's identity is settled when it was last written more than ``TRUSTED_AGE_NS`` before ``config.since_ns``,
    the capture's start, through the memo's own :class:`Witness`."""
    reading = Reading(config.repo)
    reals = {memo._real(name) for name, _ in recording.codes if not name.startswith("<")}
    code_files = {real for real in reals if reading.first_party(real)}
    data_files = {real for opened in recording.opened if reading.kept(real := memo._real(opened))}

    def state_of(real: str, code: bool) -> FileState:
        witness = boot.Witness(reference_ns=config.since_ns)
        found = witness.stat(real)
        return FileState(None if found is None else tuple(boot.FileIdentity.of(found)), witness.settled,
                         boot.file_digest(real) if found is not None and code else None)

    # each file's identity is read BEFORE its pins read its bytes
    states = {real: state_of(real, real in code_files) for real in code_files | data_files}

    collected = memo.collect_pins(recording, store)
    by_pin_path = {memo._pin_path(real): real for real in reals}

    def real_of(path: str) -> str:
        return path if os.path.isabs(path) else by_pin_path.get(path, "")

    def row_to_store(path: str) -> str:
        """An import-set row's path as the store spells the module's other pins: relative to the repository for
        first-party code (its ``StatementsPin``'s file, which stage one vouches for), the real path otherwise (#691:
        not the repository's, so ``select`` never mistakes it for a repository file of the same relative name)."""
        real = real_of(path)
        return reading.relative(real) if reading.first_party(real) else (real or path)

    def to_store(pin: pin_rows.PathPin) -> str | None:
        match pin:
            case memo.DefPin() | memo.StatementsPin():
                real = real_of(pin.path)
                return reading.relative(real) if reading.first_party(real) else None
            case _:
                real = memo._real(pin.path)
                return reading.relative(real) if reading.kept(real) and not reading.tool_owned(real) else None

    def for_store(pin: Any) -> Any:
        match pin:
            case memo.InterpreterPin() | memo.EnvironmentPin() | memo.WorkingDirectoryPin() | memo.MemoPin():
                return None
            case memo.ImportOutcomePin():
                return pin_rows.relocated_rows(pin, row_to_store)
            case _:
                return pin_rows.relocated(pin, to_store)

    ran = {for_store(pin) for pin in [*collected.ran, *shadow_pins(code_files, reading, memo, boot)]} - {None}
    # a pin the run ran and a read's entry also holds is the run's: relocation spells the two alike
    through_reads = {for_store(pin) for pin in collected.through_reads} - {None} - ran
    pins = tuple(sorted(ran | through_reads, key=pin_rows.encode))
    judged = {pin_rows.file_of(pin) for pin in pins}
    # A file only a child read from another store names (its generation read it, the run did not): its state is taken
    # now, after the generation read it; the state carries the digest it was taken on, so a change in between is a
    # stale pin, never a vouched one. Without it the store would record the file as missing (qa F3, 2026-11-03).
    held = {reading.relative(real) for real in states}
    code = {pin.path for pin in pins if isinstance(pin, (memo.DefPin, memo.StatementsPin))}
    states |= {os.path.join(reading.repo, relative): state_of(os.path.join(reading.repo, relative), relative in code)
               for relative in judged - held if relative is not None}
    files = {reading.relative(real): state for real, state in states.items() if reading.relative(real) in judged}
    unknown = [reason for refusal in collected.refusals
               if (reason := _unknown(refusal, memo, reading, real_of)) is not None]
    unknown += _changed_during_capture(files, config)
    own_pins = tuple(pin for pin in pins if pin in ran)
    own_judged = {pin_rows.file_of(pin) for pin in own_pins}
    return ReadCones(Cone(own_pins, {path: state for path, state in files.items() if path in own_judged}, tuple(unknown)),
                     Cone(pins, files, tuple(unknown)))


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
        case memo.SourceGone(path) if not os.path.isabs(path):
            # a foreign child's source, relative to a ``sys.path`` entry of the read: which entry held it, and so
            # whether it was the repository's, cannot be told once it is gone (qa F5, 2026-11-03)
            return str(refusal)
        case memo.NoRecordedBytes(path) | memo.SourceGone(path) | memo.UnnamedModule(path):
            return str(refusal) if reading.first_party(memo._real(path)) else None
        case memo.NoDefAtLine(path):
            return str(refusal) if reading.first_party(real_of(path)) else None
        case memo.AmbiguousOwner() | memo.UnownedSitePackage() | memo.UninstalledDistribution() | memo.UnmeasurableImport():
            return str(refusal)
    raise TypeError(f"cone_store: a refusal of a kind it does not know: {refusal!r}")


def _changed_during_capture(files: dict[str, FileState], config: Config) -> list[str]:
    """A file that is not settled against the capture's start (written within ``TRUSTED_AGE_NS`` before it, or
    after it) and whose identity is not the one it had when the capture began was written during the capture: the
    pins read its bytes after the test ran, so they may pin bytes that did not run."""
    return [f"{path} changed during the capture" for path, state in sorted(files.items())
            if not state.settled and state.ident is not None and config.recent.get(path) != list(state.ident)]

