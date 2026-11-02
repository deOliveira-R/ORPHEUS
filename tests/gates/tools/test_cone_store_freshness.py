"""The cone store's freshness check in two stages, and its route (spec ``vv-principles`` #26).

Spec ``scratch/verification_tree/cone_store/spec_b1.md`` §3.4, rows F1-F6. Stage one compares each pinned
file's identity (device, inode, size, modification and change times) with the identity at capture; stage
two parses only the files stage one could not clear. The observable is the set of ``.py`` files under the
fixture root that the selecting process OPENED (an audit hook in a fresh process,
``_cone_store_fixture.select``): AST-ONLY gives every right answer and opens every file, so only these rows
see it. Each row builds its own aged tree without the memo part (the memo's validation of a child entry
reads ``pkg/oracle.py`` through its own reads, which is not the store's stage two).
"""
from __future__ import annotations

import os
import time
from pathlib import Path

import pytest

import tools.cone_store as cone  # noqa: F401  (red at collection until B1 lands)

from ._cone_store_fixture import BASELINE, _trusted_age_s, build, capture, edit, select, store_of

pytestmark = pytest.mark.foundation

PKG_USER_ABOVE = frozenset({"tests/test_b.py::test_b", "tests/test_g.py::test_fail", "tests/test_m.py::test_nested"})


@pytest.fixture
def captured(tmp_path: Path) -> Path:
    root = build(tmp_path / "tree", memo=False)
    capture(root)
    return root


def test_f1_an_unchanged_checkout_is_judged_without_reading_a_source(captured: Path) -> None:
    """F1. With no file changed since capture and every identity settled, ``select`` opens no source file.
    Activation leg: after touching ``pkg/user.py`` (its change time moves), exactly that file is opened."""
    unsettled = [path for path, (_, settled) in store_of(captured).identities().items() if not settled]
    assert unsettled == [], f"precondition: the aged tree has unsettled identities {unsettled}"
    selected, opened = select(captured)
    assert selected == BASELINE
    assert opened == frozenset(), f"stage two ran on an unchanged checkout: {sorted(opened)}"
    os.utime(captured / "pkg/user.py")
    selected, opened = select(captured)
    assert opened == {"pkg/user.py"}, sorted(opened)
    assert selected == BASELINE


def test_f2_an_identity_unsettled_at_capture_is_read_again(tmp_path: Path) -> None:
    """F2. A file written less than ``TRUSTED_AGE_NS`` before the capture recorded its identity is marked
    unsettled, and ``select`` opens it at stage two even though its identity is unchanged (git's racily
    clean rule; the memo's ``Witness.settled``)."""
    root = build(tmp_path / "tree", memo=False)
    user = root / "pkg/user.py"
    user.write_text(user.read_text())
    capture(root, "tests/test_b.py")
    identities = store_of(root).identities()
    assert "pkg/user.py" in identities, sorted(identities)
    identity, settled = identities["pkg/user.py"]
    assert not settled, "precondition: the capture took longer than TRUSTED_AGE_NS, so no identity is unsettled"
    _, opened = select(root)
    assert "pkg/user.py" in opened


def test_f3_an_edit_that_keeps_size_and_modification_time_is_seen(captured: Path) -> None:
    """F3. ``>`` swapped for ``<`` in ``above`` (same size), then the old access and modification times put back
    with ``os.utime``: the change time moved, so the edit selects every test that ran ``above``."""
    user = captured / "pkg/user.py"
    old = user.stat()
    edit(captured, "pkg/user.py", "return x > TOL", "return x < TOL")
    os.utime(user, ns=(old.st_atime_ns, old.st_mtime_ns))
    new = user.stat()
    assert (new.st_size, new.st_mtime_ns) == (old.st_size, old.st_mtime_ns), "precondition: size or mtime moved"
    selected, _ = select(captured)
    assert PKG_USER_ABOVE <= selected, f"missing {sorted(PKG_USER_ABOVE - selected)}"


def _comment_and_age(root: Path) -> None:
    edit(root, "pkg/user.py", "    return x > TOL", "    # strictly above\n    return x > TOL")
    time.sleep(_trusted_age_s())


def test_f4_a_comment_edit_passes_stage_one_and_stops_at_stage_two(captured: Path) -> None:
    """F4. A comment-only edit moves the file's identity (stage one fires) and not its digests (stage two
    clears): the file is opened once and no row goes stale."""
    _comment_and_age(captured)
    selected, opened = select(captured)
    assert opened == {"pkg/user.py"}, sorted(opened)
    assert selected == BASELINE


def test_f5_a_cleared_file_has_its_identity_refreshed(captured: Path) -> None:
    """F5 (ruled 2026-11-01: refresh, as git refreshes its index). After F4's ``select`` the store records the
    file's new identity, so a second ``select``, in a new process, opens nothing."""
    _comment_and_age(captured)
    select(captured)
    selected, opened = select(captured)
    assert opened == frozenset(), f"the identity was not refreshed: {sorted(opened)}"
    assert selected == BASELINE


def test_f6_a_def_moved_down_the_file_stays_fresh(captured: Path) -> None:
    """F6. Three comment lines inserted after the imports move every def of ``pkg/user.py`` to new lines with
    bodies unchanged: no row is stale (keys are qualname and ordinal, digests strip layout)."""
    edit(captured, "pkg/user.py", "from pkg.consts import TOL\n", "from pkg.consts import TOL\n# a\n# b\n# c\n")
    selected, _ = select(captured)
    assert selected == BASELINE, f"extra {sorted(selected - BASELINE)}"
