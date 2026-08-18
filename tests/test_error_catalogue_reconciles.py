"""The error catalogue, its markers, and its index describe one set.

Three surfaces name ``ERR-NNN``, and each is maintained by a different
action:

* ``docs/theory/verification/error_catalog.rst`` — one
  ``.. error-entry::`` per defect, written by hand when a bug is caught;
* ``tests/`` — ``@pytest.mark.catches("ERR-NNN")``, written when the
  regression gate is added;
* ``.claude/skills/vv-principles/error_index.md`` — GENERATED from the
  knowledge graph and injected into the ``vv-principles`` skill.

Before the 2026-08-17 move the first two were reconciled by hand against
a markdown file, and nothing kept it true. These are pure text checks —
no graph, no venv, no build — so they run in the ordinary suite and fail
in seconds rather than at the next Sphinx build.

⚠ Each test below is a SEPARATE arm on purpose. A single "the three sets
agree" assertion is green on arrival and reds off whichever arm the
mutation happens to reach first, certifying the other two
(``vv-principles`` #17, the granularity trap).
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
CATALOG = REPO_ROOT / "docs" / "theory" / "verification" / "error_catalog.rst"
INDEX = REPO_ROOT / ".claude" / "skills" / "vv-principles" / "error_index.md"
TESTS = REPO_ROOT / "tests"

ERR_RE = re.compile(r"\bERR-\d{3}\b")


def _catalogue_ids() -> set[str]:
    """Ids DECLARED as entries — the `.. error-entry::` argument only.

    Deliberately not every `ERR-NNN` in the page: entries cross-reference
    each other constantly, and counting a mention as a declaration would
    make the page trivially agree with itself.
    """
    return set(
        re.findall(r"^\.\. error-entry:: (ERR-\d{3})\s*$",
                   CATALOG.read_text(encoding="utf-8"), re.M)
    )


def _marker_ids() -> set[str]:
    found: set[str] = set()
    for path in TESTS.rglob("test_*.py"):
        for call in re.findall(r"catches\(([^)]*)\)", path.read_text(encoding="utf-8")):
            found.update(ERR_RE.findall(call))
    return found


def _index_ids() -> set[str]:
    return set(re.findall(r"^\| (ERR-\d{3}) \|", INDEX.read_text(encoding="utf-8"), re.M))


@pytest.mark.foundation
def test_every_marker_names_a_declared_entry():
    """Arm 1. A marker with no entry READS as coverage and is not one.

    It is also a hard build failure: nexus warns per unresolved marker
    once the project has declared anything, and the canonical gate is
    ``sphinx-build -E -W``.
    """
    orphans = _marker_ids() - _catalogue_ids()
    assert not orphans, (
        f"{len(orphans)} `catches` marker id(s) name no `.. error-entry::` in "
        f"{CATALOG.relative_to(REPO_ROOT)}: {sorted(orphans)}. Declare them, or "
        "if the tag names something other than a catalogued defect (a mutation, "
        "a failure-mode family) it does not belong in `catches`."
    )


@pytest.mark.foundation
def test_every_declared_entry_has_a_catching_test():
    """Arm 2. A catalogued defect nothing pins is an unguarded regression.

    `vv-principles` defines a `catches` marker as a coverage CLAIM, so
    the absence of one is the honest signal that the claim was never
    made.
    """
    uncaught = _catalogue_ids() - _marker_ids()
    assert not uncaught, (
        f"{len(uncaught)} catalogued defect(s) have no `@pytest.mark.catches`: "
        f"{sorted(uncaught)}. Either add the gate, or say in the entry why no "
        "test can exist."
    )


@pytest.mark.foundation
def test_the_generated_index_matches_the_corpus():
    """Arm 3. The index is DERIVED; a mismatch means it went stale.

    It is injected into `vv-principles` and read by five preloading
    agents, so a stale index misinforms them with no other symptom —
    the injection succeeds either way.
    """
    catalogue, index = _catalogue_ids(), _index_ids()
    assert index == catalogue, (
        f"index/corpus disagree — only in index: {sorted(index - catalogue)}; "
        f"only in corpus: {sorted(catalogue - index)}. Regenerate with "
        "`python -m tools.verification.generate_error_index`."
    )


@pytest.mark.foundation
def test_the_id_sequence_has_no_gaps_or_duplicates():
    """Arm 4. The next free id must be unambiguous.

    The skill tells an agent to append the next sequential id; a gap
    makes "next" ambiguous and a duplicate silently merges two defects
    into one node, since the id is the graph key.
    """
    ids = sorted(_catalogue_ids())
    numbers = [int(i.split("-")[1]) for i in ids]
    assert numbers == list(range(1, len(numbers) + 1)), (
        "ERR ids are not a contiguous 1..N run — "
        f"got {len(numbers)} ids spanning {numbers[0]}..{numbers[-1]}"
    )
    raw = re.findall(r"^\.\. error-entry:: (ERR-\d{3})\s*$",
                     CATALOG.read_text(encoding="utf-8"), re.M)
    assert len(raw) == len(set(raw)), "a `.. error-entry::` id is declared twice"
