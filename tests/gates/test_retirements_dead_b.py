"""The done-when gates of sub-unit dead-B (numerics rest and eigenvalue) of the ontology W0 specification (part A, 3.2 D4, 3.7, Appendix D-B).

The shared engine (``tests/gates/test_retirements.py``) reads the rows of every ``rows_*.py`` but asserts a verdict on the real tree
only for wave 2, so the verdict of this sub-unit's rows (``tests/gates/_retired/rows_w0_dead_b.py``) is asserted here, over the rows
module of this sub-unit alone. Four gates:

* every ``Retired`` and ``Moved`` row of the module reads green on the tree and every ``Kept`` entry reads green, with the input counts
  (12 ``Retired``, 7 ``Moved``, 8 ``Kept``) and a control row the walk must read red. First red at the base: 19 of 19 rows
  (each definition is still present; each ``Moved`` home is absent), 0 of 8 ``Kept``.
* D4, no silent skip: no line of ``orpheus``, ``tests``, ``tools``, ``examples`` or ``docs`` spells the retired
  Rayleigh-quotient iteration (the rows module and the error catalogue excepted), and its six-test class is gone. First red at the
  base: the definition, its tests and its prose.
* the ``hasattr``-guarded string tuple of ``tests/gates/homogeneous/test_homogeneous.py``: the engine does not read ``hasattr``
  (decision E1 (a)), so after ``dominant_eigenpair`` and ``direct_eigenvalue`` leave ``orpheus.numerics.eigenvalue`` a surviving
  tuple entry would be a silent no-op. First red at the base: the two entries (lines 269 and 270).
* the mutation of the ``direct_eigenvalue`` row (specification 1.5, "Mutation of dead-B's row"): a synthetic root that holds the row's
  ``file`` without the definition and its ``home`` with it reads green, then (i) the definition returning to the old file and (ii) a
  ``monkeypatch.setattr(m, "direct_eigenvalue", f)`` with ``m`` bound by ``from orpheus.numerics import eigenvalue as m`` each turn
  exactly that row red. This gate is green at the base: it mutates the row, not the tree.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

from tests.gates._retired import Kept, Moved, Retired, RetiredName, _walk
from tests.gates._retired import rows_w0_dead_b as dead_b

_REPO = Path(__file__).resolve().parents[2]
_THIS = "tests/gates/test_retirements_dead_b.py"
_ENGINE = "tests/gates/test_retirements.py"
_CONTROL = RetiredName("LinearOperator")  # a name the five trees use today: the walk must read it red (X1)

# The name is spelled from two halves so that no line of this file carries it (D4's grep has no exclusion for this file).
_RQI = "rayleigh_quotient" + "_iteration"
_RQI_CLASS = "TestRayleighQuotient" + "Iteration"


def _key(row: object) -> tuple[type, object]:
    return (type(row), row)


@pytest.fixture(scope="module")
def tree_walk() -> dict:
    return _walk.walk(_REPO, rows=[*dead_b.ROWS, *dead_b.KEPT, _CONTROL])


@pytest.mark.foundation
@pytest.mark.rests_on(
    f"{_ENGINE}::test_a_kept_entry_is_green_while_defined_and_red_when_its_definition_or_file_is_gone",
    f"{_ENGINE}::test_the_home_leg_is_green_for_a_home_that_defines_the_name_and_red_for_one_that_lacks_it_or_is_absent",
    f"{_ENGINE}::test_the_positive_control_of_a_definition_row_finds_each_of_nine_planted_lines",
)
def test_every_dead_b_row_is_green_on_the_tree_and_every_kept_entry_stays_defined(tree_walk: dict) -> None:
    """Specification 3.7 done-when 1, for this sub-unit. Counts first (12 ``Retired``, 7 ``Moved``, 8 ``Kept``, the engine's
    ``Row`` types), then the control row read red, then every row green. Red input: any of the 19 definitions still present, any of
    the 7 homes missing its name, any reference by import resolution, ancestor package, ``getattr`` or ``setattr`` string, or a
    ``Kept`` definition gone."""
    by_type = {t.__name__: sum(1 for r in dead_b.ROWS if type(r) is t) for t in (Retired, Moved)}
    assert by_type == {"Retired": 12, "Moved": 7}, by_type
    assert len(dead_b.KEPT) == 8 and all(type(k) is Kept for k in dead_b.KEPT)
    assert tree_walk[_key(_CONTROL)], "the control reads green: the walk scanned no file that uses it"
    assert len(tree_walk) == len(dead_b.ROWS) + len(dead_b.KEPT) + 1, len(tree_walk)
    red = {k: tree_walk[k] for k in (_key(r) for r in dead_b.ROWS) if tree_walk[k]}
    print(f"dead-B rows green: {len(dead_b.ROWS) - len(red)} of {len(dead_b.ROWS)}")
    assert not red, f"{len(dead_b.ROWS) - len(red)} of {len(dead_b.ROWS)} rows green; red: " + "; ".join(
        f"{k[1][1]}: {v[:2]}" for k, v in sorted(red.items(), key=repr)
    )
    gone = {k: tree_walk[k] for k in (_key(e) for e in dead_b.KEPT) if tree_walk[k]}
    assert not gone, f"kept definitions gone: {sorted(map(repr, gone))}"


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
def test_no_line_spells_the_retired_rayleigh_quotient_iteration_and_its_test_class_is_gone() -> None:
    """Specification 3.2 D4. Population: the ``.py``, ``.rst`` and ``.md`` files of ``orpheus``, ``tests``, ``tools``, ``examples`` and
    ``docs`` (without ``docs/_build``), excluding ``tests/gates/_retired`` and ``docs/theory/verification/error_catalog.rst``. The
    positive control is the pattern read over the name with a following word character and over a line that holds the name (the
    filter finds a known member); the file count is asserted non-zero. Red input at the base: ``orpheus/numerics/eigenvalue.py``
    (the definition), ``tests/gates/numerics/test_eigenvalue.py`` (the six tests, bound by a ``getattr`` that skips silently) and the
    docs lines of the tense table."""
    pattern = re.compile(r"(?<![A-Za-z0-9_])" + re.escape(_RQI) + r"(?![A-Za-z0-9_])")
    assert pattern.search(f"x = getattr(_eig, '{_RQI}', None)") and not pattern.search(f"{_RQI}_x")
    seen, hits = _scan_for(
        pattern,
        ("orpheus", "tests", "tools", "examples", "docs"),
        ("tests/gates/_retired", "docs/theory/verification/error_catalog.rst"),
    )
    assert seen > 1000, seen
    assert not hits, f"{len(hits)} of {seen} files' lines spell the retired name: {hits[:12]}"
    tree = ast.parse((_REPO / "tests/gates/numerics/test_eigenvalue.py").read_text(encoding="utf-8"))
    classes = {n.name for n in ast.walk(tree) if isinstance(n, ast.ClassDef)}
    assert "TestRayleighQuotientInversion" not in classes  # a control: the filter distinguishes a neighbouring class name
    assert _RQI_CLASS not in classes, "the six tests skip silently once the symbol is gone (vv-principles mode 8(6))"


@pytest.mark.foundation
def test_no_hasattr_guarded_string_in_the_homogeneous_gate_names_a_moved_eigenvalue_function() -> None:
    """After the move neither ``orpheus.numerics.eigenvalue`` nor the homogeneous namespace holds ``dominant_eigenpair`` or
    ``direct_eigenvalue``; a surviving string entry of the ``hasattr``-guarded tuple of ``tests/gates/homogeneous/test_homogeneous.py``
    would make its spy a silent no-op, and the engine does not read ``hasattr`` (decision E1 (a)). Population: the string constants
    (non-docstring) of that file, by ``ast``. The control is the same walk over a planted source holding both strings. Red input at
    the base: the two tuple entries (lines 269 and 270)."""
    names = {"dominant_eigenpair", "direct_eigenvalue"}

    def found(source: str) -> list[int]:
        return sorted(n.lineno for n in ast.walk(ast.parse(source)) if isinstance(n, ast.Constant) and n.value in names)

    assert len(found('T = ("dominant_eigenpair", "direct_eigenvalue")\n')) == 2
    hits = found((_REPO / "tests/gates/homogeneous/test_homogeneous.py").read_text(encoding="utf-8"))
    assert not hits, f"string constants naming a moved function at lines {hits}"


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
    from the rows module, so a misspelt file, qualname or home in the rows reddens the unmutated precondition."""
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
