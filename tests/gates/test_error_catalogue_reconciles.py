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
``tests/gates/numerics/``"), never by a path.

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

REPO_ROOT = Path(__file__).resolve().parents[2]
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


def _open_ids() -> set[str]:
    """Ids of entries whose defect is still OPEN: the entry's ``**Status:**`` line reads OPEN and names its issue.

    An open defect has no catcher, by the meaning of one: a test that reddens
    under the defect would hold the suite red until it is fixed. Such an
    entry says so in its status, with the issue that tracks the fix, and is
    exempt from arm 2 until its status changes (ERR-110, #546).
    """
    text = CATALOG.read_text(encoding="utf-8")
    blocks = re.split(r"^\.\. error-entry:: ", text, flags=re.M)[1:]
    found = set()
    for block in blocks:
        entry = block.split("\n", 1)[0].strip()
        status = re.search(r"\*\*Status:\*\*(.*?)(?:\n\s*\n)", block, re.S)
        if status and "OPEN" in status.group(1) and re.search(r"#\d+", status.group(1)):
            found.add(entry)
    return found


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

    An entry whose catcher AND subject are both gone carries a ``**Retired:**`` paragraph and is exempt
    (``_retired_ids``, read from the tree; a catalogued defect cannot be re-pinned once its code is gone).
    """
    uncaught = _uncaught_ids(_catalogue_ids(), _marker_ids(), _open_ids(), _retired_ids())
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
        "``tests/gates/test_error_catalogue_reconciles.py::test_the_citation_checker_can_fail`` "
        "``tests/no_such_directory/test_absent.py`` "
        "``tests/gates/test_error_catalogue_reconciles.py::test_the_{citation,missing}_checker_can_fail``"
    )
    got = _unresolved_test_citations(text)
    assert got == [
        "line 1: tests/no_such_directory/test_absent.py (no such file)",
        "line 1: tests/gates/test_error_catalogue_reconciles.py::test_the_missing_checker_can_fail",
    ], got


#: The line an entry carries when every one of its catchers is withdrawn:
#: ``**Status:** dormant — every catcher is withdrawn under #506`` (several
#: issues are listed ``#506, #512``).
_DORMANT_RE = re.compile(
    r"^\s*\*\*Status:\*\* dormant — every catcher is withdrawn under "
    r"(#\d+(?:, #\d+)*)\s*$",
    re.M,
)


def _entry_bodies(text: str) -> dict[str, str]:
    """``ERR-NNN -> body``: the text from its ``.. error-entry::`` line up to
    the next entry (or the end of the page)."""
    heads = list(re.finditer(r"^\.\. error-entry:: (ERR-\d{3})\s*$", text, re.M))
    return {
        m.group(1): text[m.end(): (heads[i + 1].start() if i + 1 < len(heads) else len(text))]
        for i, m in enumerate(heads)
    }


def _dormancy_violations(
    dormant_errors: dict[str, list[int]], bodies: dict[str, str]
) -> list[str]:
    """Entries whose dormancy line disagrees with the audit's ``dormant_errors``.

    ``dormant_errors`` maps each entry whose catchers are ALL withdrawn to
    its withdrawal issues (derived once, in the audit). Such an entry must
    carry the dormancy line naming exactly those issues; an entry carrying
    the line while not dormant (a running catcher, or no withdrawn one) is a
    stale dormancy.
    """
    violations: list[str] = []
    for err, body in sorted(bodies.items()):
        line = _DORMANT_RE.search(body)
        if err in dormant_errors:
            wanted = ", ".join(f"#{n}" for n in dormant_errors[err])
            if line is None:
                violations.append(
                    f"{err}: every catcher is withdrawn ({wanted}) but the entry "
                    "lacks the line "
                    f"`**Status:** dormant — every catcher is withdrawn under {wanted}`"
                )
            elif line.group(1) != wanted:
                violations.append(
                    f"{err}: the dormancy line names {line.group(1)}, its catchers' "
                    f"withdrawals are {wanted}"
                )
        elif line is not None:
            violations.append(
                f"{err}: carries a dormancy line but is not dormant (a running "
                "catcher, or no withdrawn one) — a stale dormancy"
            )
    return violations


@pytest.mark.foundation
def test_a_dormant_entry_says_so_and_only_a_dormant_one():
    """Arm 7. An entry whose every catcher is withdrawn says so, in its body.

    Arm 2 counts a ``catches`` CLAIM whatever the test's other markers
    are, so a catcher marked ``@pytest.mark.withdrawn`` (skipped unless
    ``ORPHEUS_RUN_WITHDRAWN`` names its issue) keeps arm 2 green while
    nothing pins the defect. This arm asks the RUNNING question, through
    pytest's own marker resolution (the audit's ``dormant_errors``, derived
    from the registry, so no AST twin of module/class/function
    resolution), and holds the catalogue's prose to it both ways: a
    dormant entry carries ``**Status:** dormant — every catcher is
    withdrawn under #N``, and an entry with a running catcher does not.

    Its positive control is in its own body (X1): a synthetic payload with
    one dormant entry and no line, and one non-dormant entry carrying a
    stale line, must both classify as violations, and the corrected pair
    as clean.
    """
    synthetic = {"ERR-900": [506]}
    broken = {
        "ERR-900": "\n   :title: x\n\n   **Bug:** y\n",
        "ERR-901": "\n   :title: x\n\n   **Status:** dormant — every catcher is withdrawn under #506\n",
    }
    fixed = {
        "ERR-900": "\n   :title: x\n\n   **Status:** dormant — every catcher is withdrawn under #506\n",
        "ERR-901": "\n   :title: x\n\n   **Bug:** y\n",
    }
    got = _dormancy_violations(synthetic, broken)
    assert [v.split(":")[0] for v in got] == ["ERR-900", "ERR-901"], got
    assert _dormancy_violations(synthetic, fixed) == []

    from tests._harness.audit import audit_payload

    violations = _dormancy_violations(
        audit_payload()["dormant_errors"],
        _entry_bodies(CATALOG.read_text(encoding="utf-8")),
    )
    assert not violations, (
        f"{len(violations)} catalogue entr(y/ies) disagree with their catchers' "
        "withdrawal state:\n  " + "\n  ".join(violations)
    )



@pytest.mark.foundation
def test_an_open_entry_names_its_issue_and_has_no_catcher_yet():
    """Arm 2's exemption reads the status line, and only an OPEN one with an issue.

    Its positive control is ERR-110 (#546, CP's vacuum k), open since
    2026-10-10 and documented by strict xfails. An open entry that has a
    catcher is no longer open: its status is stale, and arm 2 needs no
    exemption for it.
    """
    opened = _open_ids()
    assert "ERR-110" in opened, f"the status-line parse misses ERR-110 (OPEN, #546); it found {sorted(opened)}"
    stale = opened & _marker_ids()
    assert not stale, f"entries marked OPEN carry a `catches` marker, so their status is stale: {sorted(stale)}"



# ---------------------------------------------------------------------------
# Gate D3 (specification ontology_w0_spec.md 3.1): an entry "retired with its code"
# ---------------------------------------------------------------------------

#: The paragraph an entry carries when its catcher and its subject are both gone:
#: ``**Retired:** with its code: <file>:<qualname>[, <file>:<qualname>...]`` then one sentence of date and reason. The paragraph runs to
#: the next blank line, so a wrapped paragraph is read whole; the pairs are the tokens ``path.py:qualname`` (a date sentence has none).
_RETIRED_PARA_RE = re.compile(r"^[ \t]*\*\*Retired:\*\* with its code:(.*?)(?:\n[ \t]*\n|\Z)", re.M | re.S)
_RETIRED_PAIR_RE = re.compile(r"[\w/]+\.py:\w+(?:\.\w+)*")


def _retired_pairs(text: str) -> dict[str, list[str]]:
    """``ERR-NNN -> [file:qualname, ...]`` for every entry that carries a ``**Retired:**`` paragraph."""
    found: dict[str, list[str]] = {}
    for err, body in _entry_bodies(text).items():
        paragraph = _RETIRED_PARA_RE.search(body)
        if paragraph is not None:
            found[err] = _RETIRED_PAIR_RE.findall(paragraph.group(1))
    return found


def _is_defined(root: Path, pair: str) -> bool:
    """Is ``file:qualname`` defined in ``root/file``? Read by the ONE definition finder of the retirement engine (a second finder would
    be a second spelling of one predicate, and a bare-name search would let a later homonym un-retire the entry)."""
    from tests.gates._retired._finder import definitions

    rel, qualname = pair.split(":", 1)
    path = root / rel
    return path.is_file() and qualname in definitions(path.read_text(encoding="utf-8"))


def _row_pairs() -> set[str]:
    """``file:qualname`` of every ``Retired`` and ``Moved`` row of the engine's rows modules: the one home of the retired set."""
    from tests.gates._retired import Moved, Retired, _walk

    return {f"{row.file}:{row.qualname}" for row in _walk.load_rows()[0] if isinstance(row, (Retired, Moved))}


def _retired_ids(
    root: Path = REPO_ROOT, text: str | None = None, markers: set[str] | None = None, rows: set[str] | None = None
) -> set[str]:
    """Ids of the entries that are RETIRED WITH THEIR CODE: the ``**Retired:**`` paragraph names at least one pair, none of the named
    pairs is defined under ``root``, and the entry carries no ``catches`` marker (an entry with both is stale, never retired).

    Every pair must be the ``file:qualname`` of a ``Retired`` or ``Moved`` row (``rows``, default the engine's): a pair that is no row
    is an ``AssertionError``, because a pair naming an absent file reads as a gone definition and a misspelt path would retire its
    entry for good. The rows were red at their base, which proves their paths existed; the engine's walk decides the code is gone."""
    text = CATALOG.read_text(encoding="utf-8") if text is None else text
    markers = _marker_ids() if markers is None else markers
    rows = _row_pairs() if rows is None else rows
    unknown = sorted(pair for pairs in _retired_pairs(text).values() for pair in pairs if pair not in rows)
    assert not unknown, f"`**Retired:**` pairs that are no Retired or Moved row of tests/gates/_retired: {unknown}"
    return {
        err
        for err, pairs in _retired_pairs(text).items()
        if pairs and err not in markers and not any(_is_defined(root, pair) for pair in pairs)
    }


def _uncaught_ids(catalogue: set[str], markers: set[str], opened: set[str], retired: set[str]) -> set[str]:
    return catalogue - markers - opened - retired


@pytest.mark.foundation
def test_the_retired_state_of_an_entry_is_read_from_the_tree_and_the_markers():
    """D3's positive control (X1), one synthetic entry per shape: a retired line naming a gone function retires; naming a defined
    one does not; a wrapped two-line paragraph is read whole (a pair on its second line that is still defined keeps the entry live);
    a retired line plus a marker is stale; an entry with no retired line, and one whose line names no pair, retire nothing; an absent
    file is a gone definition, but only for a pair that is a row; a pair that is no row is
    refused. Red input: a reader that takes the first line only, that searches names bare, or that ignores markers."""
    here = "tests/gates/test_error_catalogue_reconciles.py"
    gone, live = f"{here}:no_such_function_xyz", f"{here}:_retired_ids"

    def entry(err: str, *lines: str) -> str:
        return f".. error-entry:: {err}\n\n   :title: x\n\n" + "".join(f"   {line}\n" for line in lines) + "\n"

    text = "".join([
        entry("ERR-900", f"**Retired:** with its code: {live}. 2026-11-20, a date and a reason."),
        entry("ERR-901", f"**Retired:** with its code: {gone}. 2026-11-20, a date and a reason."),
        entry("ERR-902", f"**Retired:** with its code: {gone},", f"{here}:also_gone. 2026-11-20 and a reason."),
        entry("ERR-903", f"**Retired:** with its code: {gone},", f"{live}. 2026-11-20 and a reason."),
        entry("ERR-904", f"**Retired:** with its code: {gone}. 2026-11-20."),
        entry("ERR-905", "**Bug:** nothing retired here."),
        entry("ERR-906", "**Retired:** with its code: nothing to parse here."),
        entry("ERR-907", "**Retired:** with its code: tests/gates/no_such_file_at_all.py:thing."),
    ])
    rows = {live, gone, f"{here}:also_gone", "tests/gates/no_such_file_at_all.py:thing"}
    assert _retired_ids(REPO_ROOT, text, markers={"ERR-904"}, rows=rows) == {"ERR-901", "ERR-902", "ERR-907"}
    assert set(_retired_pairs(text)) == {"ERR-90%d" % k for k in (0, 1, 2, 3, 4, 6, 7)}
    # the negative control: a misspelt path (a pair that is no row) is refused, not read as a gone definition
    with pytest.raises(AssertionError, match="no Retired or Moved row"):
        _retired_ids(REPO_ROOT, text, markers={"ERR-904"}, rows=rows - {"tests/gates/no_such_file_at_all.py:thing"})
    assert _uncaught_ids({"ERR-1", "ERR-2"}, set(), set(), {"ERR-2"}) == {"ERR-1"}


@pytest.mark.foundation
@pytest.mark.rests_on("tests/gates/test_error_catalogue_reconciles.py::test_the_retired_state_of_an_entry_is_read_from_the_tree_and_the_markers")
def test_err_118_is_retired_with_its_code_and_no_retired_entry_is_stale():
    """Dead-A, gate D3. ERR-118's catcher called ``select_quadrature`` and its subject was ``QuadratureSpec.structural_flags``; both
    retire with the selector, so the entry carries the ``**Retired:**`` paragraph naming the two ``file:qualname`` pairs, and no entry
    that carries one is stale (a pair still defined, or a marker still present). Red input in today's tree: the entry carries no such
    paragraph (and both pairs are defined)."""
    text = CATALOG.read_text(encoding="utf-8")
    pairs = _retired_pairs(text).get("ERR-118", [])
    assert pairs == [
        "orpheus/numerics/quadrature/registry.py:select_quadrature",
        "orpheus/numerics/quadrature/registry.py:QuadratureSpec.structural_flags",
    ], pairs
    stale = set(_retired_pairs(text)) - _retired_ids()
    assert not stale, f"entries carry a `**Retired:**` paragraph but are not retired (a pair is defined, or a marker remains): {sorted(stale)}"
    assert "ERR-118" in _retired_ids()


def _retired_root(tmp_path: Path, append: str = "") -> Path:
    """A root holding only the CURRENT ``registry.py`` at its relative path, with ``append`` added at its end."""
    rel = Path("orpheus/numerics/quadrature/registry.py")
    target = tmp_path / rel
    target.parent.mkdir(parents=True)
    target.write_text((REPO_ROOT / rel).read_text(encoding="utf-8") + append, encoding="utf-8")
    return tmp_path


@pytest.mark.foundation
@pytest.mark.rests_on("tests/gates/test_error_catalogue_reconciles.py::test_err_118_is_retired_with_its_code_and_no_retired_entry_is_stale")
@pytest.mark.parametrize(
    "append, markers",
    [
        ("\n\ndef select_quadrature():\n    return None\n", set()),
        ("\n\nclass QuadratureSpec:\n    def structural_flags(self):\n        return None\n", set()),
        ("", {"ERR-118"}),
    ],
    ids=["the-selector-redefined", "the-flags-accessor-redefined", "a-marker-returns"],
)
def test_err_118_leaves_the_retired_set_when_a_pair_is_redefined_or_a_marker_returns(tmp_path: Path, append: str, markers: set[str]):
    """Gate D3's mutation, one arm per leg of the predicate (``vv-principles`` #17a). The root is a copy of ``registry.py`` alone, at its
    relative path, so the tracked file is never touched: the unmutated copy must retire ERR-118 (the precondition: a reading that is
    red before the mutation proves nothing), and re-adding ONE named definition, or a ``catches`` marker, must un-retire it. Red input
    in today's tree: the precondition (the entry carries no ``**Retired:**`` paragraph, the copy still defines both pairs)."""
    honest = _retired_ids(_retired_root(tmp_path / "honest"), markers=set())
    assert "ERR-118" in honest, "precondition: the unmutated copy retires ERR-118"
    mutant = _retired_ids(_retired_root(tmp_path / "mutant", append), markers=markers)
    assert "ERR-118" not in mutant
