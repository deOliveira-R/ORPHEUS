"""The traced memo beside the cone store's capture, and the witness of the memo's ``collect_pins`` split.

Spec ``scratch/verification_tree/cone_store/spec_b1.md`` §3.7 and §6, rows M1, M2, and M3 (the orchestrator's
addition, 2026-11-01). The capture activates the memo's ``Recording`` so ``note_child`` and the audit hook
reach it; M1 runs the memo's own suites under the capture, M2 nests a ``trace_call`` inside a captured test.

M3 needs no B1 and runs today; M1 and M2 import ``tools.cone_store`` inside the row, so they error (never
pass) until B1 lands.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

from orpheus.numerics import _traced_memo_boot as boot
from orpheus.numerics import traced_memo

from ._cone_store_fixture import REPO, code_keys

pytestmark = pytest.mark.foundation


def _two_refusals(scratch: Path) -> None:
    """Start a program the memo does not admit, then read a file and change it: two distinct refusals."""
    subprocess.run(["/bin/echo", "x"], capture_output=True, check=True)
    path = scratch / "seen.txt"
    path.read_text()
    path.write_text("changed")


def test_m3_collect_pins_returns_every_refusal_and_manifest_raises_the_first(tmp_path: Path) -> None:
    """M3. A recording with two distinct refusals (an unadmitted program, a file observed then changed) gives
    both from ``collect_pins`` as typed ``Refusal`` values, in the run's event order (programs before files, as
    ``_manifest`` read them): ``UnadmittedProcess`` naming ``/bin/echo``, then ``ChangedByRun`` naming the file.
    Each ``str`` is the message the memo raised before the split. The other pins are kept (the def that ran).
    ``_manifest`` raises the first, as an ``Unpinnable`` carrying that refusal."""
    seen = tmp_path / "seen.txt"
    seen.write_text("before")
    recording = boot.Recording().start()
    try:
        _two_refusals(tmp_path)
    finally:
        recording.stop()
    store = traced_memo.Store(tmp_path / "store")
    collected = traced_memo.collect_pins(recording, store)
    kinds = [type(r) for r in collected.refusals]
    assert kinds == [traced_memo.UnadmittedProcess, traced_memo.ChangedByRun], collected.refusals
    process, changed = collected.refusals
    assert isinstance(process, traced_memo.UnadmittedProcess) and isinstance(changed, traced_memo.ChangedByRun)
    assert os.path.realpath(process.program) == os.path.realpath("/bin/echo"), process
    assert os.path.realpath(changed.path) == os.path.realpath(seen), changed
    assert str(process).startswith(f"the run started a process ({process.program}) whose own reads and children no recording sees")
    assert str(changed) == f"{changed.path}: the run observed it and then changed it, so no state of it pins the answer"
    defs = {p.qualname for p in collected.pins if isinstance(p, traced_memo.DefPin)}
    assert "_two_refusals" in defs, sorted(defs)
    with pytest.raises(traced_memo.Unpinnable) as raised:
        traced_memo._manifest(recording, store)
    assert raised.value.refusal == process
    assert str(raised.value) == str(process)


def test_m3_an_unowned_site_package_is_a_refusal_never_raised() -> None:
    """M3, qa's third leg. A code object whose file lies in site-packages and is listed by no installed
    distribution's ``RECORD`` reaches ``origin()`` inside ``collect_pins`` and comes back as an
    ``UnownedSitePackage`` refusal; ``collect_pins`` does not raise. The code object is compiled under a file
    name in the real site directory that no file holds (nothing is written there)."""
    site_dir = traced_memo._site_dirs()[0]
    filename = os.path.join(site_dir, "cone_store_unowned_probe_module.py")
    assert not os.path.exists(filename), f"precondition: {filename} exists"
    code = compile("value = 1\n", filename, "exec")
    recording = boot.Recording()
    recording.codes[(filename, code)] = code
    collected = traced_memo.collect_pins(recording, traced_memo.Store(Path(site_dir) / "never-written"))
    unowned = [r for r in collected.refusals if isinstance(r, traced_memo.UnownedSitePackage)]
    assert [r.path for r in unowned] == [filename], collected.refusals


def test_m3_a_source_the_run_deleted_is_a_refusal_never_raised(tmp_path: Path) -> None:
    """M3, the deleted-source leg (M1's ``test_m2_2`` flip, 2026-11-01). A run imports a module, runs its def,
    and deletes the module's file. In a process without the pinned-source loader (a capture child,
    ``trace_call``) ``collect_pins`` reads that source after the run: it returns a ``SourceGone`` naming the
    file, keeps the other pins (the driving def that ran), and does not raise. ``_manifest`` raises it."""
    import importlib
    import sys

    package = tmp_path / "pkg_gone"
    package.mkdir()
    module = package / "gone_module.py"
    module.write_text("def f():\n    return 1\n")
    sys.path.insert(0, str(package))
    try:
        recording = boot.Recording().start()
        try:
            _run_and_delete(importlib, module)
        finally:
            recording.stop()
    finally:
        sys.path.remove(str(package))
        sys.modules.pop("gone_module", None)
    store = traced_memo.Store(tmp_path / "store")
    collected = traced_memo.collect_pins(recording, store)
    gone = [r for r in collected.refusals if isinstance(r, traced_memo.SourceGone)]
    assert [os.path.realpath(r.path) for r in gone] == [os.path.realpath(module)], collected.refusals
    defs = {p.qualname for p in collected.pins if isinstance(p, traced_memo.DefPin)}
    assert "_run_and_delete" in defs, sorted(defs)
    with pytest.raises(traced_memo.Unpinnable) as raised:
        traced_memo._manifest(recording, store)
    assert isinstance(raised.value.refusal, traced_memo.SourceGone), raised.value.refusal


def _run_and_delete(importlib, module: Path) -> None:
    """Import ``gone_module``, run its def, delete its file."""
    imported = importlib.import_module("gone_module")
    assert imported.f() == 1
    module.unlink()


def test_m2_a_nested_trace_call_is_seen_by_the_capture_and_unpolluted_by_it(tmp_path: Path) -> None:
    """M2. ``test_nested`` calls ``trace_call(user.above, 1.0)`` inside a captured child: the capture's rows hold
    ``above``, and the inner manifest holds ``above`` and no def of ``tools/cone_store`` (asserted inside the
    fixture test, whose outcome must be passed)."""
    import tools.cone_store as cone  # noqa: F401

    from ._cone_store_fixture import build, capture, store_of

    root = build(tmp_path / "tree", memo=False)
    capture(root, "tests/test_m.py")
    store = store_of(root)
    assert store.tests()["tests/test_m.py::test_nested"].outcome == "passed"
    assert ("DefPin", "pkg/user.py", "above", 0) in code_keys(store.pins("tests/test_m.py::test_nested"))


_MEMO_SUITES = (*sorted(str(p.relative_to(REPO)) for p in (REPO / "tests/gates/numerics").glob("test_traced_memo_*.py")),
                "tests/gates/derivations/test_traced_memo_clients.py")


def _plain_outcomes(files: tuple[str, ...]) -> dict[str, str]:
    run = subprocess.run([sys.executable, "-O", "-m", "pytest", "-q", "-rA", "--color=no", "-p", "no:cacheprovider", *files],
                         cwd=REPO, env=dict(os.environ, OMP_NUM_THREADS="1"), capture_output=True, text=True, timeout=3600)
    found = re.findall(r"(PASSED|FAILED|ERROR|SKIPPED|XFAIL|XPASS) (\S+::\S+)", run.stdout)
    return {test: outcome.lower() for outcome, test in found}


@pytest.mark.slow
def test_m1_the_memos_own_suites_give_the_same_outcomes_under_the_capture(tmp_path: Path) -> None:
    """M1 (RECORD). The memo's own test files give the same outcome per test captured by the store as run
    plainly. A red names a production change for the orchestrator (spec §6.1: ``os.stat`` stays probed while any
    recording is active), never a workaround in ``tools/``."""
    import tools.cone_store as cone

    assert len(_MEMO_SUITES) >= 2, _MEMO_SUITES
    plain = _plain_outcomes(_MEMO_SUITES)
    assert plain, "the plain run reported no outcome"
    store_path = tmp_path / "store.sqlite"
    run = subprocess.run([sys.executable, "-O", "-m", "tools.cone_store", "capture", "--repo", str(REPO), "--store",
                          str(store_path), *_MEMO_SUITES], cwd=REPO, env=dict(os.environ, OMP_NUM_THREADS="1"),
                         capture_output=True, text=True, timeout=3600)
    assert run.returncode == 0, run.stderr[-3000:]
    captured = {t: r.outcome for t, r in cone.open_store(store_path).tests().items()}
    differ = {t: (plain.get(t), captured.get(t)) for t in set(plain) | set(captured) if plain.get(t) != captured.get(t)}
    if differ:
        pytest.fail("outcomes differ (plain, captured):\n" + "\n".join(f"{t}: {v}" for t, v in sorted(differ.items())))
