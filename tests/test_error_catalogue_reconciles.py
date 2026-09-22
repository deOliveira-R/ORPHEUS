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

A fourth surface is the catalogue's own prose: each entry names the tests
that catch it by path (``tests/<dir>/test_x.py::Class::test_y``). Arms 5
and 6 hold that prose to the tree, because a re-layout moves test files
and nothing else notices (`[M]` 2026-09-22: 30 of 99 cited paths and 12
of 77 cited functions were dead after the ``sn/`` re-layout). The
convention the arm enforces: a path citation is a claim that the test
exists NOW; a retired test is named in history by module and directory
("the ``test_projection_operators`` module, then under
``tests/numerics/``"), never by a path.

⚠ Each test below is a SEPARATE arm on purpose. A single "the three sets
agree" assertion is green on arrival and reds off whichever arm the
mutation happens to reach first, certifying the other two
(``vv-principles`` #17, the granularity trap).
"""

from __future__ import annotations

import ast
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
    """Ids named by a ``catches(...)`` CALL — a decorator, a ``pytestmark``
    entry, a ``pytest.param(marks=...)`` — read by an AST pass.

    Not a text grep: a test that fixtures the harness's ID checker carries
    the literal text ``catches("ERR-999")`` inside a string, and a grep
    counted it as a marker, so arm 1 read red on ``main`` from 2026-09-20
    (``399b285f``) until 2026-09-22. A membership question is parsed.
    Arm 2 is this census's positive control: a census that missed real
    markers would leave declared entries uncaught and redden it.
    """
    found: set[str] = set()
    for path in TESTS.rglob("test_*.py"):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "catches"
            ):
                for arg in node.args:
                    if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                        found.update(ERR_RE.findall(arg.value))
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


#: A test cited by path: ``tests/<…>.py`` then up to two ``::`` segments,
#: each an identifier that may carry ONE brace alternation
#: (``test_si_krylov_eigenvalue_equivalence_{sphere,cylinder}``).
_CITED_TEST_RE = re.compile(
    r"(tests/[A-Za-z0-9_./-]+\.py)"
    r"((?:::[A-Za-z0-9_]*(?:\{[A-Za-z0-9_, ]+\})?[A-Za-z0-9_]*)*)"
)


def _expand_braces(segment: str) -> list[str]:
    """``a_{x, y}_b`` → ``[a_x_b, a_y_b]``; a plain name → ``[name]``."""
    m = re.fullmatch(r"([A-Za-z0-9_]*)\{([^}]*)\}([A-Za-z0-9_]*)", segment)
    if m is None:
        return [segment]
    head, alts, tail = m.groups()
    return [head + alt.strip() + tail for alt in alts.split(",")]


def _defined_names(path: Path) -> tuple[set[str], dict[str, set[str]]]:
    """Module-level names (defs, classes, assignments) and each class's methods,
    read by an AST pass (a membership question is parsed, never grepped)."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    top: set[str] = set()
    methods: dict[str, set[str]] = {}
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            top.add(node.name)
        elif isinstance(node, ast.ClassDef):
            top.add(node.name)
            methods[node.name] = {
                m.name for m in node.body
                if isinstance(m, (ast.FunctionDef, ast.AsyncFunctionDef))
            }
        elif isinstance(node, ast.Assign):
            top.update(t.id for t in node.targets if isinstance(t, ast.Name))
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            top.add(node.target.id)
    return top, methods


def _unresolved_test_citations(text: str) -> list[str]:
    """Every ``tests/…py[::A[::B]]`` in ``text`` that names no file, or names a
    class, function or module-level name the file does not define."""
    cache: dict[Path, tuple[set[str], dict[str, set[str]]]] = {}
    missing: list[str] = []
    for m in _CITED_TEST_RE.finditer(text):
        rel, tail = m.group(1), m.group(2)
        line = text.count("\n", 0, m.start()) + 1
        path = REPO_ROOT / rel
        if not path.is_file():
            missing.append(f"line {line}: {rel} (no such file)")
            continue
        segments = [s for s in tail.split("::") if s][:2]
        if not segments:
            continue
        if path not in cache:
            cache[path] = _defined_names(path)
        top, methods = cache[path]
        if len(segments) == 1:
            wanted = [(name,) for name in _expand_braces(segments[0])]
        else:
            wanted = [
                (cls, fn)
                for cls in _expand_braces(segments[0])
                for fn in _expand_braces(segments[1])
            ]
        for name in wanted:
            found = (
                name[0] in top if len(name) == 1
                else name[1] in methods.get(name[0], set())
            )
            if not found:
                missing.append(f"line {line}: {rel}::{'::'.join(name)}")
    return missing


@pytest.mark.foundation
def test_every_test_the_catalogue_cites_exists():
    """Arm 5. An entry's "caught by" prose names a test that exists NOW.

    `vv-principles` reads the catalogue as the record of which gate pins
    which defect; a citation to a moved or deleted test sends the reader
    to nothing while the entry reads as covered. The ``catches`` markers
    are arm 2's business; this arm is the PROSE.
    """
    missing = _unresolved_test_citations(CATALOG.read_text(encoding="utf-8"))
    assert not missing, (
        f"{len(missing)} test citation(s) in {CATALOG.relative_to(REPO_ROOT)} "
        "resolve to nothing:\n  " + "\n  ".join(missing) + "\nRe-point each to "
        "the test carrying the entry's `catches` marker, or, for a retired "
        "gate, name it in history by module and directory, not by path."
    )


@pytest.mark.foundation
def test_the_citation_checker_can_fail():
    """Arm 6. Arm 5's positive control (X1): a green arm 5 is evidence
    only if the checker reddens on a dead file and on a dead name, and
    resolves a live one, brace alternation included."""
    text = (
        "``tests/test_error_catalogue_reconciles.py::test_the_citation_checker_can_fail`` "
        "``tests/no_such_directory/test_absent.py`` "
        "``tests/test_error_catalogue_reconciles.py::test_the_{citation,missing}_checker_can_fail``"
    )
    got = _unresolved_test_citations(text)
    assert got == [
        "line 1: tests/no_such_directory/test_absent.py (no such file)",
        "line 1: tests/test_error_catalogue_reconciles.py::test_the_missing_checker_can_fail",
    ], got
