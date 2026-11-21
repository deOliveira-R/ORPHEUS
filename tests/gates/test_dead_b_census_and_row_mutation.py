"""The gates of sub-unit dead-B (numerics rest and eigenvalue) that the shared retirement engine does not assert.

The engine (``tests/gates/test_retirements.py``) asserts every row and every ``Kept`` entry of ``tests/gates/_retired/rows_w0_dead_b.py``
green on the real tree, and reads ``hasattr``/``getattr``/``setattr`` string accesses, so this file holds only:

* D4 (ontology W0 specification 3.2), no silent skip: no line of ``orpheus``, ``tests``, ``tools``, ``examples`` or ``docs`` spells the
  retired Rayleigh-quotient iteration (the rows module and the error catalogue excepted), and its six-test class is gone.
* the mutation of the ``direct_eigenvalue`` row (specification 1.5): a synthetic root holding the row's ``file`` without the definition
  and its ``home`` with it reads green, then (i) the definition returning to the old file and (ii) a ``monkeypatch.setattr`` on
  ``from orpheus.numerics import eigenvalue as m`` each turn exactly that row red, through the engine's seam ``_walk.walk(root, rows=...)``.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

from tests.gates._retired import Moved, _walk
from tests.gates._retired import rows_w0_dead_b as dead_b

_REPO = Path(__file__).resolve().parents[2]
_ENGINE = "tests/gates/test_retirements.py"

# The name is spelled from two halves so that no line of this file carries it (D4's grep has no exclusion for this file).
_RQI = "rayleigh_quotient" + "_iteration"
_RQI_CLASS = "TestRayleighQuotient" + "Iteration"


def _key(row: object) -> tuple[type, object]:
    return (type(row), row)


def _scan_for(pattern: re.Pattern[str], trees: tuple[str, ...], excluded: tuple[str, ...]) -> tuple[int, list[str]]:
    hits: list[str] = []
    seen = 0
    for tree in trees:
        for path in sorted((_REPO / tree).rglob("*")):
            rel = path.relative_to(_REPO).as_posix()
            if not path.is_file() or path.suffix not in {".py", ".rst", ".md"} or "_build" in path.parts:
                continue
            if any(rel == e or rel.startswith(e + "/") for e in excluded):
                continue
            seen += 1
            for n, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
                if pattern.search(line):
                    hits.append(f"{rel}:{n}")
    return seen, hits


@pytest.mark.foundation
def test_no_line_spells_the_retired_rqi_and_its_test_class_is_gone() -> None:
    """Specification 3.2 D4. Population: the ``.py``, ``.rst`` and ``.md`` files of ``orpheus``, ``tests``, ``tools``, ``examples`` and
    ``docs`` (without ``docs/_build``), excluding ``tests/gates/_retired`` and ``docs/theory/verification/error_catalog.rst``. The
    pattern is a plain substring, as the specification's ``git grep`` predicate; the positive control is the pattern read over a line that
    holds the name and over one that embeds it in a longer identifier; the file count is asserted non-zero. Red input at the base: ``orpheus/numerics/eigenvalue.py``
    (the definition), ``tests/gates/numerics/test_eigenvalue.py`` (the six tests, bound by a ``getattr`` that skips silently) and the
    docs lines of the tense table."""
    pattern = re.compile(re.escape(_RQI))  # a plain substring, as the specification's git grep predicate (X2)
    assert pattern.search(f"x = getattr(_eig, '{_RQI}', None)") and pattern.search(f"def test_{_RQI}_x")
    seen, hits = _scan_for(
        pattern,
        ("orpheus", "tests", "tools", "examples", "docs"),
        ("tests/gates/_retired", "docs/theory/verification/error_catalog.rst"),
    )
    assert seen > 1000, seen
    assert not hits, f"{len(hits)} of {seen} files' lines spell the retired name: {hits[:12]}"
    tree = ast.parse((_REPO / "tests/gates/numerics/test_eigenvalue.py").read_text(encoding="utf-8"))
    classes = {n.name for n in ast.walk(tree) if isinstance(n, ast.ClassDef)}
    assert "TestStandardEigenproblem" in classes, sorted(classes)  # the positive control: the AST filter finds a class that exists
    assert _RQI_CLASS not in classes, "the six tests skip silently once the symbol is gone (vv-principles mode 8(6))"


@pytest.mark.foundation
@pytest.mark.rests_on(
    f"{_ENGINE}::test_a_moved_row_goes_red_when_its_old_definition_returns_and_when_its_home_loses_the_name",
    f"{_ENGINE}::test_a_use_planted_in_the_copy_turns_exactly_its_name_row_red",
)
def test_the_direct_eigenvalue_row_goes_red_when_its_definition_returns_and_when_a_setter_names_it(tmp_path: Path) -> None:
    """The mutation of dead-B's row, the root-path mechanism (specification 1.3 and 1.5): the unmutated root is green for the row
    (required first); (i) the definition planted back in the root's ``numerics/eigenvalue.py`` turns exactly the row red by the
    definition leg; (ii) on a fresh root a test file holding ``from orpheus.numerics import eigenvalue as m`` and
    ``monkeypatch.setattr(m, "direct_eigenvalue", f)`` turns exactly the row red by the setter leg, in that file. The row is read
    from the rows module; the synthetic root is built from the row's own fields, so this gate does not check the spelling of those
    fields (gate 1 reads them on the tree)."""
    row = next(r for r in dead_b.ROWS if type(r) is Moved and r.qualname == "direct_eigenvalue")

    def root(sub: str) -> Path:
        base = tmp_path / sub
        for rel, body in ((row.file, "def power_iteration(): pass\n"), (row.home, f"def {row.home_name}(): pass\n")):
            (base / rel).parent.mkdir(parents=True, exist_ok=True)
            (base / rel).write_text(body, encoding="utf-8")
        return base

    def read(base: Path) -> list[tuple[str, int]]:
        return _walk.walk(base, rows=[row])[_key(row)]

    first = root("definition")
    assert read(first) == [], "the unmutated root is not green: a reading red before the mutation proves nothing"
    (first / row.file).write_text("def power_iteration(): pass\ndef direct_eigenvalue(): pass\n", encoding="utf-8")
    assert [f for f, _ in read(first)] == [row.file]

    second = root("setter")
    assert read(second) == []
    planted = "tests/gates/numerics/test_planted_spy.py"
    (second / planted).write_text(
        "from orpheus.numerics import eigenvalue as m\n\n\ndef test_x(monkeypatch):\n    monkeypatch.setattr(m, 'direct_eigenvalue', None)\n",
        encoding="utf-8",
    )
    assert {f for f, _ in read(second)} == {planted}
