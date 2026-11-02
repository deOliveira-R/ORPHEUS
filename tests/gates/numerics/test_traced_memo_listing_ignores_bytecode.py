r"""A directory listing leaves the interpreter's bytecode cache out.

``__pycache__`` appears in a package directory whenever something imports from it and writes bytecode, so a listing
that counted it moved the listing pin of every package directory after any ordinary run: the cone store then
re-selected the whole suite (2026-11-01). A real name added still moves it. Leaving it out is an accepted limit of the
generator contract (``_LISTING_IGNORED``): a run that lists subdirectories without filtering ``__pycache__`` is not
re-generated when one appears.
"""
from __future__ import annotations

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
