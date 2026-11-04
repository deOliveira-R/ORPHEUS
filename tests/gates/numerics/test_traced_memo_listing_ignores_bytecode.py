r"""A directory listing leaves the interpreter's bytecode cache out.

``__pycache__`` appears in a package directory whenever something imports from it and writes bytecode, so a listing
that counted it moved the listing pin of every package directory after any ordinary run: the cone store then
re-selected the whole suite (2026-11-01). A real name added still moves it. Leaving it out is an accepted limit of the
generator contract (``_LISTING_IGNORED``): a run that lists subdirectories without filtering ``__pycache__`` is not
re-generated when one appears.
"""
from __future__ import annotations

import os

import pytest

from orpheus.numerics import _traced_memo_boot as boot

pytestmark = pytest.mark.foundation


def test_a_bytecode_cache_leaves_the_listing_unchanged_and_a_new_module_moves_it(tmp_path):
    package = tmp_path / "listed"
    package.mkdir()
    (package / "mod.py").write_text("X = 1\n")
    before = boot.listing_digest(str(package))
    (package / "__pycache__").mkdir()
    if boot.listing_digest(str(package)) != before:
        pytest.fail("creating __pycache__ moved the directory's listing digest")
    (package / "new.py").write_text("Y = 2\n")
    if boot.listing_digest(str(package)) == before:
        pytest.fail("adding new.py left the listing digest unchanged (the activation leg)")


@pytest.mark.catches("ERR-125")
def test_a_walk_records_no_listing_of_a_bytecode_cache(tmp_path):
    """A run that walks a tree lists each directory in it, ``__pycache__`` included. Its listing changes whenever
    bytecode is written, so it re-selected every source-walking test on every ``select`` (``[M]`` 2026-11-03, B3
    acceptance X1.3: 10 of 229 tests selected with no edit). A directory of an ignored name is not an input: its
    listing is not recorded, while its parent's still is (the activation leg)."""
    package = tmp_path / "walked"
    (package / "__pycache__").mkdir(parents=True)
    (package / "mod.py").write_text("X = 1\n")
    recording = boot.Recording().start()
    try:
        walked = [names for _, _, names in os.walk(package)]
    finally:
        recording.stop()
    assert walked
    if str(package) not in recording.listed:
        pytest.fail(f"the walked directory's own listing was not recorded: {sorted(recording.listed)}")
    if cached := [path for path in recording.listed if os.path.basename(path) == "__pycache__"]:
        pytest.fail(f"the listing of a bytecode cache was recorded: {cached}")
