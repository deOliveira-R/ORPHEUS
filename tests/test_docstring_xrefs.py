"""Gate: no docstring names a Python object that does not exist.

An unresolved Python-domain cross-reference (``:func:`` / ``:class:`` /
``:meth:`` / ``:mod:`` / ``:attr:``) renders as **plain text** and draws no
``-W`` warning. The natural next move — a nitpicky (``-n``) build — does not
help either, and the reason is structural rather than a config gap:

**Sphinx can only diagnose what it RENDERS.** A docstring in a module that is
never ``automodule``'d is never read, so its roles never become nodes and
nothing can warn at any severity. This tree carries 45 live ``automodule``
directives against 263 modules in ``orpheus/``, and ``tests/`` is never
rendered at all while holding ~500 fully-qualified xrefs.

The project's knowledge graph inherits the same blind spot for the same
reason: Nexus marks a reference ``unresolved`` when it finds no target NODE,
and a target node exists only if Sphinx rendered the target — so live symbols
in un-rendered modules land in that bucket beside genuinely dead ones, and it
cannot be gated on.

So this gate resolves each target by **importing** it, via
``tools/check_docstring_xrefs.py``. Render coverage is irrelevant.

It must be a subprocess: the checker imports arbitrary modules while probing,
and doing that inside the pytest process would pollute ``sys.modules`` for
every test that runs after it.

When this gate reds, the fix is per-site adjudication, not a blanket delete —
see ``.claude/rules/coding-standards.md`` "Retire as you go": a past-tense
mention is correct history and becomes a ``literal``; a present-tense claim or
an imperative instruction is a MUST-FIX. Do NOT silence a hit by adding to the
checker's ``ALLOWLIST``, which ships empty deliberately.
"""

from __future__ import annotations

import pathlib
import subprocess
import sys

import pytest

REPO_ROOT = pathlib.Path(__file__).parent.parent
CHECKER = REPO_ROOT / "tools" / "check_docstring_xrefs.py"


@pytest.mark.foundation
@pytest.mark.parametrize("tree", ["orpheus", "tests"])
def test_no_docstring_names_a_missing_object(tree: str) -> None:
    """Every fully-qualified xref in ``tree`` resolves against the live namespace.

    Parametrised by tree so a production regression and a test-tree regression
    are distinguishable in the report rather than sharing one red row.
    """
    result = subprocess.run(
        [sys.executable, str(CHECKER), tree],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )
    assert result.returncode == 0, (
        f"Dead Sphinx cross-references in {tree}/ — each one names a Python "
        f"object that no longer exists. No Sphinx build can see these (see this "
        f"module's docstring for why), so this gate is the only thing standing "
        f"between a module move and a docstring that lies.\n\n"
        f"Adjudicate each site (repoint / past-tense literal / rewrite), do not "
        f"blanket-delete, and do not add to the checker's ALLOWLIST.\n\n"
        f"{result.stdout}{result.stderr}"
    )
