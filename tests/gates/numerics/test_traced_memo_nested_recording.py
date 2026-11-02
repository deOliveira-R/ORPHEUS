r"""An enclosing recording does not take the memo's own reads as its run's dependencies.

``trace_call`` builds its manifest after its own recording stops, and building it re-reads the source of every
file that ran (in a process without the pinned-source loader). A recording still active around the call (the
cone store's capture of a test that calls ``trace_call``) was offered those reads, and recorded the file's whole
bytes as the caller's data dependency, so a docstring edit to it re-selected the caller (found by the cone
store's builder, 2026-11-01). The re-read is the memo's own, made under ``unrecorded``.
"""
from __future__ import annotations

import importlib
import sys

import pytest

from orpheus.numerics import _traced_memo_boot, traced_memo

pytestmark = pytest.mark.foundation


@pytest.fixture
def traced_module(tmp_path, monkeypatch):
    (tmp_path / "nested_subject.py").write_text("def value():\n    return 3\n")
    monkeypatch.syspath_prepend(str(tmp_path))
    yield tmp_path / "nested_subject.py"
    sys.modules.pop("nested_subject", None)


def _call() -> int:
    return importlib.import_module("nested_subject").value()


def test_the_inner_manifest_still_pins_the_file(traced_module):
    """The inner call's own manifest pins the subject (the activation leg: the re-read happens)."""
    _, manifest = traced_memo.trace_call(_call)
    if not any(getattr(pin, "path", "") == "nested_subject.py" for pin in manifest.pins):
        pytest.fail(f"the inner manifest does not pin nested_subject.py: {sorted(map(str, manifest.pins))}")


@pytest.mark.rests_on("tests/gates/numerics/test_traced_memo_nested_recording.py::test_the_inner_manifest_still_pins_the_file")
def test_an_enclosing_recording_is_not_offered_the_manifests_re_read(traced_module):
    """Around a ``trace_call``, an outer recording holds no read of the subject's source: the import's read is
    the loader's, and the manifest's re-read is the memo's."""
    outer = _traced_memo_boot.Recording().start()
    try:
        traced_memo.trace_call(_call)
    finally:
        outer.stop()
    subject = str(traced_module.resolve())
    read = [path for path in outer.opened if traced_memo._real(path) == subject]
    if read:
        pytest.fail(f"the enclosing recording took the memo's re-read of the source as a dependency: {read}")


def test_the_memos_digest_of_its_own_source_is_offered_to_no_recording(monkeypatch):
    """The memo reads its own three source files to name itself (``_memo_digest``, run when it is imported). A
    recording active then (the cone store's capture records from before the memo's import) is not offered the read,
    which would make the memo's whole bytes a dependency of every test. A fresh keeper forces the read to happen."""
    monkeypatch.setattr(traced_memo, "_MEMO_DIGESTS", _traced_memo_boot.Kept(traced_memo._read_memo_digest))
    recording = _traced_memo_boot.Recording().start()
    try:
        traced_memo._memo_digest()
    finally:
        recording.stop()
    sources = {traced_memo._real(str(path)) for path in traced_memo._MEMO_SOURCES}
    read = sorted(path for path in recording.opened if traced_memo._real(path) in sources)
    if read:
        pytest.fail(f"a recording was offered the memo's read of its own source: {read}")
