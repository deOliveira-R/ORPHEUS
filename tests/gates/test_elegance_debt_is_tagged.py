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

from tests.gates import _refusal_gates as rg

pytestmark = pytest.mark.foundation

_TOKEN = re.compile(r"ELEGANCE-DEBT\[guard\]")
_ISSUE = re.compile(r"#\d+")
#: A sentence that NAMES a change: "retires when …", "until …", "when …".
_RETIRING_CHANGE = re.compile(r"\b(retires?|retired|until|when)\b", re.IGNORECASE)


def _repo_root() -> pathlib.Path:
    return pathlib.Path(__file__).resolve().parents[2]


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


def _tagged_guards(
    root: pathlib.Path | None = None,
) -> tuple[list[str], list[str], list[tuple[str, str]]]:
    """``(token sites, violations, token lines)`` over every ``*.py`` under ``root``.

    ``root`` is an ``orpheus/`` directory (default: the repository's); a mutation is a copy of it
    (``plan`` 5.3: every walk takes a root path). A token line is ``(site, text)``; the wave leg
    reads those, so the token is walked once.
    """
    root = _repo_root() / "orpheus" if root is None else root
    sites: list[str] = []
    violations: list[str] = []
    token_lines: list[tuple[str, str]] = []
    for path in sorted(root.rglob("*.py")):
        text = path.read_text(encoding="utf-8")
        rel = path.relative_to(root.parent).as_posix()
        for i, line in enumerate(text.splitlines()):
            if _TOKEN.search(line):
                sites.append(f"{rel}:{i + 1}")
                token_lines.append((f"{rel}:{i + 1}", line))
        violations.extend(_violations_in(text, label=rel))
    return sites, violations, token_lines


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
    sites, violations, _ = _tagged_guards()
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
    if not any("geometry/manifold.py" in s for s in sites):
        pytest.fail("the orbit catalogue's door carries no SCOPE-BOUNDARY[guard] token; the ledger's known member is missing")
    if violations:
        pytest.fail("SCOPE-BOUNDARY[guard] ledger violations:\n  " + "\n  ".join(violations))


# ------------------------------------------------------------------ the wave leg

#: The issues of the W0 units that tag a refusal (``ontology_w0_spec.md`` 5.3, decision 9): each token line
#: naming one carries ``#<issue> (W<digit>[ab]?)``, the wave that retires the refusal, written ONCE, in the tag.
#: The wave that retires a refusal deletes the refusal and its tag in one diff and removes its issue from this set.
_WAVE_TAGGED_ISSUES = frozenset({709, 710, 711})


def _wave_violations(token_lines: list[tuple[str, str]]) -> list[str]:
    """Every token line naming an issue of the set without ``#<issue> (W<digit>[ab]?)`` right after the number, and every
    issue of the set with no token line at all. The wave token sits INSIDE parentheses after the number: a bare
    ``\\bW3\\b`` would be met by the tag's own sentence ``Retires when W3`` and a wrong wave would pass."""
    out: list[str] = []
    for issue in sorted(_WAVE_TAGGED_ISSUES):
        mine = [(site, line) for site, line in token_lines if re.search(rf"#{issue}(?!\d)", line)]
        if not mine:
            out.append(f"#{issue}: no ELEGANCE-DEBT[guard] token line names it")
        for site, line in mine:
            if not re.search(rf"#{issue}(?!\d)\s*\(W\d[ab]?\)", line):
                out.append(f"{site}: the tag of #{issue} carries no wave in the form '#{issue} (W<digit>)'")
    return out


def test_the_wave_predicate_reads_the_parenthesised_wave_and_nothing_else():
    """POSITIVE CONTROL, inside the body: the wave in the sentence only is a violation; a well-formed tag is clean;
    ``#7090`` is not ``#709``; an issue with no token line is a violation."""
    tags = lambda *lines: [(f"p:{i}", ln) for i, ln in enumerate(lines)]  # noqa: E731
    well = [
        "# ELEGANCE-DEBT[guard] #709 (W3): x. Retires when y.",
        "# ELEGANCE-DEBT[guard] #710 (W1): x. Retires when y.",
        "# ELEGANCE-DEBT[guard] #711 (W8a): x. Retires when y.",
    ]
    if _wave_violations(tags(*well)):
        pytest.fail(f"the predicate flagged well-formed tags: {_wave_violations(tags(*well))}")
    sentence_only = [well[0].replace("#709 (W3)", "#709"), *well[1:]]
    if len(_wave_violations(tags(*sentence_only))) != 1:
        pytest.fail("the predicate did not flag a wave written in the sentence only")
    wrong_form = [well[0].replace("(W3)", "(3)"), *well[1:]]
    if len(_wave_violations(tags(*wrong_form))) != 1:
        pytest.fail("the predicate did not flag a wave that is not W<digit>")
    prefix_only = ["# ELEGANCE-DEBT[guard] #7090 (W3): x. Retires when y.", *well[1:]]
    if not any(v.startswith("#709:") for v in _wave_violations(tags(*prefix_only))):
        pytest.fail("the predicate read #7090 as #709")
    if len(_wave_violations(tags(*well[1:]))) != 1:
        pytest.fail("the predicate did not flag an issue with no token line")


def _imported_orpheus() -> pathlib.Path:
    """The ``orpheus/`` directory pytest imported (the repository's, except in a run on a patched copy)."""
    return rg.import_root() / "orpheus"


def test_every_refusal_tag_of_w0_carries_the_wave_that_retires_it():
    """First red: 0 of the 3 tag lines exist today. The leg makes the tag carry its wave and nothing more; it does not
    force the retirement (the wave's own done-when does)."""
    _, _, token_lines = _tagged_guards(_imported_orpheus())
    violations = _wave_violations(token_lines)
    if violations:
        pytest.fail("wave leg violations:\n  " + "\n  ".join(violations))


def test_mutation_deleting_a_wave_from_the_tag_is_one_violation_naming_its_issue(tmp_path):
    """Mutation of the wave leg: copy every ``orpheus/**/*.py`` (no ``__pycache__``) to ``<tmp>/orpheus``, require the
    unmutated copy green (a reading red before the mutation proves nothing), delete ``(W3)`` from the #709 tag line
    of ``transport/fields/_bases.py`` in the copy, and require exactly one violation, naming 709."""
    rg.copy_orpheus_py(rg.import_root(), tmp_path)
    copy = tmp_path / "orpheus"
    _, _, lines = _tagged_guards(copy)
    before = _wave_violations(lines)
    if before:
        raise RuntimeError(f"UNINSTALLABLE: the unmutated copy is already red, so the mutation proves nothing: {before}")
    bases = copy / "transport" / "fields" / "_bases.py"
    text = bases.read_text(encoding="utf-8")
    tag = [ln for ln in text.splitlines() if "ELEGANCE-DEBT[guard] #709" in ln]
    if len(tag) != 1 or tag[0].count("(W3)") != 1:
        raise RuntimeError(f"UNINSTALLABLE: the #709 tag line is not unique with one '(W3)': {tag}")
    bases.write_text(text.replace(tag[0], tag[0].replace("(W3)", "")), encoding="utf-8")
    _, _, mutated = _tagged_guards(copy)
    after = _wave_violations(mutated)
    if len(after) != 1 or "#709" not in after[0]:
        pytest.fail(f"expected exactly one violation naming 709, got {after}")
