r"""The recorder never records its own apparatus, and records the rest of its file like any other code.

The audit hook, the ``os.stat`` probe, what they call (``_observe``, ``_path``, ``_by_importer``) and the
``unrecorded`` stretch the probe enters run inside every recording, before any quiet stretch can begin. Recorded,
they put the same six recorder defs in every cone the cone store read, and homed every test with no def of its own
at the recorder (the archivist, ``[M]`` 2026-11-03: 10 of 10 own cones of the B3 sample). They are declared at
their definition (``_apparatus``) and ``started`` skips them. First red: ``_APPARATUS`` emptied (the leak back).
"""
from __future__ import annotations

import os
import subprocess
import sys
import types
from pathlib import Path

import pytest

from orpheus.numerics import _traced_memo_boot

pytestmark = pytest.mark.foundation

_REPO = Path(__file__).resolve().parents[3]


def _recorded_boot_defs(run) -> tuple[set[str], _traced_memo_boot.Recording]:
    recording = _traced_memo_boot.Recording().start()
    try:
        run()
    finally:
        recording.stop()
    boot = os.path.realpath(_traced_memo_boot.__file__)
    return {code.co_qualname for (path, code) in recording.codes if os.path.realpath(path) == boot}, recording


def test_a_probed_path_records_no_apparatus(tmp_path):
    """A run that probes a path and reads a file: the probe fired (the path is observed, the activation leg), and
    no def of the recorder is among the codes the run started."""
    target = tmp_path / "probed.txt"
    target.write_text("x")

    def run():
        assert os.path.exists(target)
        with open(target) as handle:
            handle.read()

    names, recording = _recorded_boot_defs(run)
    if str(target) not in recording.probed:
        pytest.fail(f"the stat probe did not observe {target}: the gate does not reach the hook path")
    if names:
        pytest.fail(f"the recorder recorded its own apparatus: {sorted(names)}")


def test_the_rest_of_the_recorder_is_recorded_like_any_code():
    """A run that calls a recorder def outside the apparatus records it: the file is not excluded wholesale, so a
    test of the recorder still depends on what it ran."""
    names, _ = _recorded_boot_defs(lambda: _traced_memo_boot.digest(b"x"))
    if "digest" not in names:
        pytest.fail(f"the recorder's own `digest` was not recorded: {sorted(names)}")


def test_every_hook_the_recording_installs_is_declared_apparatus():
    """The callbacks a recording installs (the audit hook and the ``os.stat``/``os.lstat`` probes) are declared,
    read off the installed objects rather than a list of names."""
    recording = _traced_memo_boot.Recording().start()
    try:
        installed = (os.stat, os.lstat)
    finally:
        recording.stop()
    for hook in (*installed, _traced_memo_boot._audit):
        code = hook.__code__
        if (code.co_filename, code) not in _traced_memo_boot._APPARATUS:
            pytest.fail(f"{code.co_qualname} is installed by a recording but not declared apparatus")


# ── a nested recording, every hook path (qa's F2 on B3, 2026-11-03) ─────────────────────────────────────────────────

def _nest(code: types.CodeType) -> set[tuple[str, types.CodeType]]:
    """``code`` and every body nested in it, by (file, code): what one ``@_apparatus`` declares."""
    found, bodies = set(), [code]
    while bodies:
        body = bodies.pop()
        found.add((body.co_filename, body))
        bodies.extend(const for const in body.co_consts if isinstance(const, types.CodeType))
    return found


#: The declarations by (file, code): the defs ``@_apparatus`` was applied to, each standing for its nest (the bodies
#: nested in it are declared with it), the population the arm row ranges over. ``_probing``'s probe is declared when
#: a recording first replaces ``os.stat``, so one probe is made here (its code object is the one every probe shares).
_traced_memo_boot._probing("stat")
_NESTED = {key for _, code in _traced_memo_boot._APPARATUS for key in _nest(code) - {(code.co_filename, code)}}
_DECLARED = sorted(_traced_memo_boot._APPARATUS - _NESTED, key=lambda key: key[1].co_qualname)


#: The memo entry :func:`_nested_run` notes as read.
_CHILD = _traced_memo_boot.ReadChild("/apparatus-store", "apparatus", "child", "digest", frozenset(), (), "/")


def _nested_run(tmp_path) -> tuple[_traced_memo_boot.Recording, _traced_memo_boot.Recording]:
    """An outer recording around an inner one that starts and stops, as the cone store's worker nests a test's own
    recording inside the inherited one, with every hook path taken inside both: a ``stat`` probe, a read, a write,
    a listing, a spawn, a memo child noted, a ``spawning`` stretch, and an ``unrecorded`` stretch holding a
    ``recorded`` one (so ``unrecorded.__exit__`` runs with recording on). The inner recording is constructed before
    the outer starts, so no def of the recorder is the run's own code."""
    boot = _traced_memo_boot
    read, written = tmp_path / "read.txt", tmp_path / "written.txt"
    read.write_text("x")
    inner = boot.Recording()
    outer = boot.Recording().start()  # replaces os.stat: ``_probing`` and the ``@_apparatus`` on its probe run here
    try:
        inner.start()
        try:
            os.stat(read)
            with open(read) as handle:
                handle.read()
            written.write_text("y")
            os.listdir(tmp_path)
            subprocess.run([sys.executable, "-c", ""], check=True)
            boot.note_child(_CHILD)
            with boot.spawning():
                pass
            with boot.unrecorded():
                with boot.recorded():
                    pass
        finally:
            inner.stop()
    finally:
        outer.stop()
    return outer, inner


def _boot_codes(recording) -> set[str]:
    boot = os.path.realpath(_traced_memo_boot.__file__)
    return {code.co_qualname for (path, code) in recording.codes if os.path.realpath(path) == boot}


def test_a_nested_recording_records_no_apparatus_on_any_hook_path(tmp_path):
    """F2. The outer recording of :func:`_nested_run` starts no def of the recorder, and every hook path fired in
    both recordings (the activation legs: the probed, read, written, listed, spawned and child fields).

    First red ``[M]`` 2026-11-03 (a ``-p`` plugin rebinding ``_APPARATUS`` without one declaration, refusing its
    re-declaration): 13 of the 17 declarations redden this row, among them the four qa found unwitnessed
    (``Recording.start``, ``Recording.stop``, ``note_child``, ``spawning.__enter__``), and ``spawning.__exit__``,
    ``unrecorded.__exit__``. The other four (``_audit``, ``_access``, ``_note_written``, ``_program``) run only
    inside the audit hook, which CPython calls with monitoring paused (no ``__cantrace__``): un-declared, they are
    still never recorded, so no run can witness them today; the arm row below makes the hook traceable for them."""
    outer, inner = _nested_run(tmp_path)
    read, written = str(tmp_path / "read.txt"), str(tmp_path / "written.txt")
    for name, recording in (("outer", outer), ("inner", inner)):
        missing = [field for field, fired in (
            ("probed", read in recording.probed), ("opened", read in recording.opened),
            ("written", written in recording.written), ("listed", str(tmp_path) in recording.listed),
            ("spawned", bool(recording.spawned)), ("children", _CHILD in recording.children))
            if not fired]
        if missing:
            pytest.fail(f"the {name} recording saw no {missing}: the run does not reach those hook paths")
    names = _boot_codes(outer)
    if names:
        pytest.fail(f"the outer recording recorded the recorder's apparatus: {sorted(names)}")


def _in_a_fresh_interpreter(entry: str, *args: str) -> None:
    """Run this module's ``entry(*args)`` in a fresh interpreter and fail on a non-zero exit, with what it printed.

    Two declarations are reached only by the FIRST recording of a process: ``Recording.start`` installs the
    ``os.stat`` probes (``_probing``, whose probe ``_apparatus`` declares) only while no recording is active. In
    process, a recording already active around the test (the cone store's capture child records every test) hides
    both, and the arm row read red under capture while green in a plain run (``[M]`` 2026-11-03, the
    memo-coexistence gate M1: ``[_apparatus]`` and ``[_probing]``). A fresh interpreter also keeps the traceable
    audit hook these runs set from reaching a recording around the test."""
    run = subprocess.run(
        [sys.executable, "-O", "-c",
         f"import sys; from tests.gates.numerics.test_traced_memo_apparatus import {entry}; "
         f"sys.exit({entry}(*sys.argv[1:]))", *args],
        cwd=_REPO, capture_output=True, text=True, timeout=300)
    if run.returncode != 0:
        pytest.fail(run.stdout + run.stderr[-3000:])


def _traceable(directory: str) -> int:
    """:func:`_nested_run` with the audit hook made traceable (``__cantrace__``): exit 0 when the outer recording
    started no def of the recorder."""
    setattr(_traced_memo_boot._audit, "__cantrace__", True)
    outer, _ = _nested_run(Path(directory))
    if names := _boot_codes(outer):
        print(f"with the audit hook traceable, the outer recording recorded the recorder's apparatus: {sorted(names)}")
        return 1
    return 0


def test_with_the_audit_hook_traceable_a_nested_recording_records_no_apparatus(tmp_path):
    """The audit hook and the defs only it calls run with monitoring paused today, so the row above cannot see a
    body of theirs left undeclared. Made traceable, it can. First red ``[M]`` 2026-11-03: ``_program``'s generator
    expression, recorded while ``@_apparatus`` declared a def's own code object and not the bodies nested in it."""
    _in_a_fresh_interpreter("_traceable", str(tmp_path))


def _arm(qualname: str, directory: str) -> int:
    """The arm of :func:`test_each_declaration_is_reached_by_the_nested_run`: the declaration named ``qualname``
    un-declared with its nest (and refused if declared again, as the ``@_apparatus`` on ``_probing``'s probe
    re-declares it on each call), the audit hook made traceable, and :func:`_nested_run` run. Exit 0 when the outer
    recording recorded it."""
    declared = next(key for key in _DECLARED if key[1].co_qualname == qualname)
    refused = _nest(declared[1])

    class Without(set):
        def add(self, key):
            if key not in refused:
                super().add(key)

    _traced_memo_boot._APPARATUS = Without(_traced_memo_boot._APPARATUS - refused)
    # CPython runs an audit hook with monitoring paused unless the hook sets ``__cantrace__``, so ``_audit`` and the
    # defs only it calls (``_access``, ``_note_written``, ``_program``) are never started under a recording today.
    # The arm makes the hook traceable, so their declarations are witnessed for the day it is.
    setattr(_traced_memo_boot._audit, "__cantrace__", True)
    outer, _ = _nested_run(Path(directory))
    if declared in outer.codes:
        return 0
    print(f"{qualname} un-declared, and the nested run never recorded it: no hook path reaches it "
          f"(recorded: {sorted(_boot_codes(outer))})")
    return 1


@pytest.mark.rests_on("tests/gates/numerics/test_traced_memo_apparatus.py::"
                      "test_a_nested_recording_records_no_apparatus_on_any_hook_path")
@pytest.mark.parametrize("declared", _DECLARED, ids=[code.co_qualname for _, code in _DECLARED])
def test_each_declaration_is_reached_by_the_nested_run(declared, tmp_path):
    """Every declaration is load-bearing: with that one def un-declared, the outer recording of :func:`_nested_run`
    records it. A declaration no hook path reaches reddens here, as a dead declaration or a path the run lost. The
    arm runs in a fresh interpreter (:func:`_in_a_fresh_interpreter`).
    ``[M]`` 2026-11-03: 16 declarations, 16 green; a case reddens when ``_nested_run`` drops the step that reaches
    its def (``recorded`` replaced by an empty stretch: ``unrecorded.__exit__`` red)."""
    if not _DECLARED:
        pytest.fail("no declared apparatus: the population is empty")
    _in_a_fresh_interpreter("_arm", declared[1].co_qualname, str(tmp_path))
