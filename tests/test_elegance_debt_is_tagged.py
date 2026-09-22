r"""Every ``ELEGANCE-DEBT[guard]`` token in ``orpheus/`` carries its issue and its
retiring change — the debt ledger the ``coding-standards`` rule promises.

**The rule** (``.claude/rules/coding-standards.md``, "A guard is elegance debt
— tag it, and name what retires it", ``[R]`` user 2026-09-07): a runtime guard
is a protection *today* and never the target state; every guard that lands
carries the token ``ELEGANCE-DEBT[guard]``, the issue number that prices its
retirement, and ONE sentence naming the structural change that makes the
guarded mistake unspellable. ``grep -rn "ELEGANCE-DEBT" orpheus/`` is then the
ledger — and a ledger is only as honest as its entries, which is what this
gate reads.

**Three legs, per ``plan-authoring`` §6c (a gate lands with the case it
catches):**

* **structure** — every token line carries ``#<issue>`` and a non-empty
  sentence within three lines that names a change (the predicate below);
* **non-vacuity** — at least ONE tagged guard exists in ``orpheus/`` (`[M]`
  the token count was **0** before CS4c step 6 item 6.3 landed the first,
  ``FullField.require_member`` / #457 — a gate over an empty population is
  vacuously green, so this file lands in the SAME commit as the first token);
* **positive control, inside the body** — a fixture string carrying the token
  WITHOUT an issue number is classified as a violation by the same predicate
  that scans the tree, so a broken predicate cannot read the tree as clean
  (``code-search``: validate a filter against a known member before trusting
  its negative).

Python ``re`` over ``pathlib.rglob``, never the shell (``grep`` here is
ugrep, and a completeness claim is re-run in Python). Foundation mark:
a corpus invariant, no physics claim.
"""

from __future__ import annotations

import pathlib
import re

import pytest

pytestmark = pytest.mark.foundation

_TOKEN = re.compile(r"ELEGANCE-DEBT\[guard\]")
_ISSUE = re.compile(r"#\d+")
#: A sentence that NAMES a change: "retires when …", "until …", "when …".
_RETIRING_CHANGE = re.compile(r"\b(retires?|retired|until|when)\b", re.IGNORECASE)


def _repo_root() -> pathlib.Path:
    return pathlib.Path(__file__).resolve().parents[1]


def _violations_in(text: str, *, label: str) -> list[str]:
    """Every token in ``text`` that lacks an issue number on its line or a
    retiring-change sentence within ±3 lines — the ONE predicate both the
    tree scan and the positive control run."""
    lines = text.splitlines()
    out: list[str] = []
    for i, line in enumerate(lines):
        if not _TOKEN.search(line):
            continue
        if not _ISSUE.search(line):
            out.append(f"{label}:{i + 1}: token without an issue number")
            continue
        window = " ".join(lines[max(0, i - 3): i + 4])
        if not _RETIRING_CHANGE.search(window):
            out.append(
                f"{label}:{i + 1}: token names no retiring change within ±3 lines"
            )
    return out


def _tagged_guards() -> tuple[list[str], list[str]]:
    """``(token sites, violations)`` over every ``orpheus/**/*.py``."""
    root = _repo_root() / "orpheus"
    sites: list[str] = []
    violations: list[str] = []
    for path in sorted(root.rglob("*.py")):
        text = path.read_text(encoding="utf-8")
        rel = path.relative_to(_repo_root()).as_posix()
        for i, line in enumerate(text.splitlines()):
            if _TOKEN.search(line):
                sites.append(f"{rel}:{i + 1}")
        violations.extend(_violations_in(text, label=rel))
    return sites, violations


def test_the_predicate_flags_a_tag_without_an_issue_number():
    """POSITIVE CONTROL — the same predicate the tree scan runs must classify
    an untagged-by-issue token as a violation, and a well-formed one as
    clean. Without this leg an empty ``violations`` list is compatible with a
    predicate that matches nothing."""
    bad = "    ELEGANCE-DEBT[guard] — retires when B is bound on its own end.\n"
    good = (
        "    ELEGANCE-DEBT[guard] #457 — a runtime guard is a protection;\n"
        "    it retires when every leaf is bound on the end it acts on.\n"
    )
    if not _violations_in(bad, label="control"):
        pytest.fail("the predicate did not flag a token without an issue number")
    if _violations_in(good, label="control"):
        pytest.fail(f"the predicate flagged a well-formed tag: {_violations_in(good, label='control')}")
    missing_change = "    ELEGANCE-DEBT[guard] #999 — a protection for now.\n"
    if not _violations_in(missing_change, label="control"):
        pytest.fail("the predicate did not flag a token that names no retiring change")


def test_every_tagged_guard_carries_its_issue_and_its_retiring_change():
    """The ledger: ≥ 1 tagged guard in ``orpheus/`` (non-vacuity), and every
    one carries ``#<issue>`` on its line and names the change that retires it
    within three lines."""
    sites, violations = _tagged_guards()
    if not sites:
        pytest.fail(
            "no ELEGANCE-DEBT[guard] token in orpheus/ — the ledger is empty, "
            "so this gate is vacuous; either a guard lost its tag or the last "
            "guard retired and this file should retire with it"
        )
    if violations:
        pytest.fail(
            "ELEGANCE-DEBT[guard] ledger violations:\n  " + "\n  ".join(violations)
        )


# ------------------------------------------------------------ scope boundaries

#: A guard at the declared edge of machinery a user ruling chose not to build
#: (``coding-standards``, "A guard is debt, or it is a declared scope boundary",
#: ``[R]`` user 2026-09-22). It carries three fields within three lines: the
#: machinery that would derive the refused property, the ruling that defers
#: it, and what would overturn that ruling.
_BOUNDARY = re.compile(r"SCOPE-BOUNDARY\[guard\]")
_BOUNDARY_FIELDS = ("machinery:", "ruling:", "revisit:")


def _boundary_violations_in(text: str, *, label: str) -> list[str]:
    """Every boundary token in ``text`` missing a field within ±3 lines — the
    ONE predicate both the tree scan and the positive control run."""
    lines = text.splitlines()
    out: list[str] = []
    for i, line in enumerate(lines):
        if not _BOUNDARY.search(line):
            continue
        window = " ".join(lines[max(0, i - 3): i + 4]).lower()
        missing = [f for f in _BOUNDARY_FIELDS if f not in window]
        if missing:
            out.append(f"{label}:{i + 1}: boundary token without {', '.join(missing)}")
    return out


def _boundary_guards() -> tuple[list[str], list[str]]:
    root = _repo_root() / "orpheus"
    sites: list[str] = []
    violations: list[str] = []
    for path in sorted(root.rglob("*.py")):
        text = path.read_text(encoding="utf-8")
        rel = path.relative_to(_repo_root()).as_posix()
        sites += [f"{rel}:{i + 1}" for i, line in enumerate(text.splitlines()) if _BOUNDARY.search(line)]
        violations.extend(_boundary_violations_in(text, label=rel))
    return sites, violations


def test_the_boundary_predicate_flags_a_tag_missing_a_field():
    """POSITIVE CONTROL — a boundary token with all three fields is clean, and
    one missing any field is a violation, by the predicate the scan runs."""
    good = (
        "    SCOPE-BOUNDARY[guard] — machinery: the Gröbner-basis engine.\n"
        "    ruling: D0.1, the user, 2026-08-31. revisit: a request outside.\n"
    )
    if _boundary_violations_in(good, label="control"):
        pytest.fail(f"the predicate flagged a well-formed boundary: {_boundary_violations_in(good, label='control')}")
    for field in _BOUNDARY_FIELDS:
        bad = good.replace(field, "")
        if not _boundary_violations_in(bad, label="control"):
            pytest.fail(f"the predicate did not flag a boundary missing {field!r}")


def test_every_scope_boundary_names_its_machinery_ruling_and_revisit_condition():
    """The boundary ledger: ≥ 1 declared boundary in ``orpheus/`` (the orbit
    catalogue's door, non-vacuity), and every one names the machinery that
    would derive the refused property, the ruling, and what would overturn it."""
    sites, violations = _boundary_guards()
    if not any("numerics/manifold.py" in s for s in sites):
        pytest.fail("the orbit catalogue's door carries no SCOPE-BOUNDARY[guard] token; the ledger's known member is missing")
    if violations:
        pytest.fail("SCOPE-BOUNDARY[guard] ledger violations:\n  " + "\n  ".join(violations))
