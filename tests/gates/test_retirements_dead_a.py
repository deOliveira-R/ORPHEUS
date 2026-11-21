"""The gates of sub-unit dead-A of unit ``dead`` (specification ``ontology_w0_spec.md`` section 3, Appendix D-A.7).

The rows are data, ``tests/gates/_retired/rows_w0_dead_a.py``; the walk that reads them is the shared engine's
(``tests/gates/test_retirements.py``), and this file adds only what is dead-A's: that its 25 ``Retired`` and ``Moved`` rows read green
on the real tree and its 30 ``Kept`` entries read green, and the mutations of the specification's 1.5 and 3.1 ("run by dead-A") that
prove each leg can fail on those rows.

Red input in today's (base) tree: every ``Retired`` and ``Moved`` row of the module is defined there (25 of 25), so
``test_every_dead_a_row_is_gone_and_every_kept_entry_stands`` and every mutation test, whose precondition is that the unmutated
copy reads green for the row it mutates, are red until the unit lands. The ``Kept`` legs are green on both sides of the landing.

Every mutation runs on a COPY of the ``.py`` files of the five trees under a temporary directory (the ``orpheus/`` directory is
hundreds of megabytes of data files; the copy keeps the ``.py`` files only), so the tracked tree is never touched.
"""

from __future__ import annotations

import ast
import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest

_REPO = Path(__file__).resolve().parents[2]
_ENGINE = "tests/gates/test_retirements.py"
_THIS = "tests/gates/test_retirements_dead_a.py"
_TREES = ("orpheus", "tests", "tools", "examples", "student_resources")
_CONTROL = "LinearOperator"  # a name the five trees use today: the row a walk that scanned no file cannot read red


@pytest.fixture(scope="module")
def dead_a() -> SimpleNamespace:
    from tests.gates._retired import Kept, Moved, Retired, RetiredName, _walk, rows_w0_dead_a

    return SimpleNamespace(
        Kept=Kept, Moved=Moved, Retired=Retired, RetiredName=RetiredName, walk=_walk.walk,
        rows=tuple(rows_w0_dead_a.ROWS), kept=tuple(rows_w0_dead_a.KEPT),
        selector=rows_w0_dead_a.SELECTOR, selection_error=rows_w0_dead_a.SELECTION_ERROR,
    )


def _red(verdicts: dict) -> set[tuple[type, object]]:
    return {key for key, findings in verdicts.items() if findings}


@pytest.fixture(scope="module")
def real_walk(dead_a: SimpleNamespace) -> dict:
    """One walk of the real tree over the dead-A rows and the control row."""
    return dead_a.walk(_REPO, rows=[*dead_a.rows, *dead_a.kept, dead_a.RetiredName(_CONTROL)])


@pytest.mark.foundation
@pytest.mark.rests_on(
    f"{_ENGINE}::test_the_positive_control_of_a_definition_row_finds_each_of_nine_planted_lines",
    f"{_ENGINE}::test_the_home_leg_is_green_for_a_home_that_defines_the_name_and_red_for_one_that_lacks_it_or_is_absent",
)
def test_every_dead_a_row_is_gone_and_every_kept_entry_stands(dead_a: SimpleNamespace, real_walk: dict) -> None:
    """The population is counted first (14 ``Retired``, 11 ``Moved``, 30 ``Kept``, as the specification's 1.5 counts them), the control
    row reads red (a walk that scanned nothing would read every other row green), then the red set must be empty: no retired
    definition is defined, no reference to one survives by import resolution, each moved name is defined at its home and gone from its
    file, and each kept definition is still defined. Red input at the base tip: the 25 ``Retired`` and ``Moved`` rows are all defined
    (the ``Kept`` entries are green)."""
    counts = {name: sum(1 for r in dead_a.rows if type(r).__name__ == name) for name in ("Retired", "Moved")}
    assert counts == {"Retired": 14, "Moved": 11} and len(dead_a.kept) == 30, (counts, len(dead_a.kept))
    assert real_walk[(dead_a.RetiredName, dead_a.RetiredName(_CONTROL))], "the control reads green: the walk scanned no file that uses it"
    assert len(real_walk) == 14 + 11 + 30 + 1, len(real_walk)
    red = _red(real_walk) - {(dead_a.RetiredName, dead_a.RetiredName(_CONTROL))}
    print(f"dead-A rows green: {55 - len(red)} of 55")
    assert not red, f"{55 - len(red)} of 55 dead-A entries green; red: " + "; ".join(f"{k[1]}: {real_walk[k][:2]}" for k in sorted(red, key=repr))


@pytest.mark.foundation
def test_the_kept_entries_are_defined_on_the_real_tree(dead_a: SimpleNamespace, real_walk: dict) -> None:
    """The keep leg alone (a test-only input, green today and after the unit): the 30 entries are still defined. Split from the test above
    so that the rows' red at the base tip does not hide this leg's reading. Input count beside the finding count: 30 and 0."""
    kept_keys = [(dead_a.Kept, entry) for entry in dead_a.kept]
    assert len(kept_keys) == 30
    red = [key for key in kept_keys if real_walk[key]]
    assert not red, red


# ---------------------------------------------------------------------------
# Mutations: one copy of the five trees' .py files, one edit per arm, restored from bytes after each
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def tree_copy(tmp_path_factory: pytest.TempPathFactory) -> Path:
    root = tmp_path_factory.mktemp("dead_a_root")

    def only_python(directory: str, names: list[str]) -> list[str]:
        """The names ``copytree`` must skip: every file that is not a ``.py``, and every ``__pycache__``."""
        base = Path(directory)
        return [n for n in names if n == "__pycache__" or ((base / n).is_file() and not n.endswith(".py"))]

    for tree in _TREES:
        if (_REPO / tree).is_dir():
            shutil.copytree(_REPO / tree, root / tree, ignore=only_python)
    return root


def _cut_definition(source: str, qualname: str) -> str:
    """``source`` without the definition ``qualname`` (decorators included): its line is the engine's finder's (the one qualname
    walker), the span is the ``ast`` node that starts on that line."""
    from tests.gates._retired._finder import definitions

    line = definitions(source)[qualname]
    node = next(
        n for n in ast.walk(ast.parse(source))
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef, ast.Assign, ast.AnnAssign)) and n.lineno == line
    )
    lines = source.splitlines(keepends=True)
    first = min([node.lineno, *(d.lineno for d in getattr(node, "decorator_list", []))])
    del lines[first - 1:node.end_lineno]
    return "".join(lines)


def _walk_copy(dead_a: SimpleNamespace, root: Path, edits: dict[str, str]) -> dict:
    """The verdicts on ``root`` with each ``{relative path: new text}`` of ``edits`` applied, the original bytes restored afterwards."""
    originals = {rel: (root / rel).read_bytes() if (root / rel).exists() else None for rel in edits}
    try:
        for rel, text in edits.items():
            (root / rel).parent.mkdir(parents=True, exist_ok=True)
            (root / rel).write_text(text, encoding="utf-8")
        return dead_a.walk(root, rows=[*dead_a.rows, *dead_a.kept])
    finally:
        for rel, data in originals.items():
            if data is None:
                (root / rel).unlink(missing_ok=True)
            else:
                (root / rel).write_bytes(data)


def _gains(dead_a: SimpleNamespace, root: Path, edits: dict[str, str], target: tuple[type, object]) -> set[tuple[type, object]]:
    """The rows that the edit turns red, after requiring ``target`` GREEN on the unmutated copy (a reading that is red before the
    mutation proves nothing, X1): a mutation test asserts the red set of its OWN row and no other."""
    before = dead_a.walk(root, rows=[*dead_a.rows, *dead_a.kept])
    assert not before[target], f"precondition: {target[1]} reads red on the unmutated copy: {before[target][:2]}"
    return _red(_walk_copy(dead_a, root, edits)) - _red(before)


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_kept_entries_are_defined_on_the_real_tree")
def test_deleting_a_kept_definition_reddens_exactly_its_entry(dead_a: SimpleNamespace, tree_copy: Path) -> None:
    """Spec 1.5, D1': in the copy, ``DiscreteMeasure.symmetry_groups`` (row 58, kept) is cut from ``numerics/measure.py``; the red
    ``Kept`` set must be exactly that entry. The mutation root holds the ``.py`` files of all five trees: on a root of one file 47 of the
    50 ``Kept`` entries of the whole unit are red whether or not the definition is deleted, so "some red" proves nothing there. Red
    input: the entry itself, once its definition is gone (a test-only input; the entry is green in today's tree)."""
    rel = "orpheus/numerics/measure.py"
    target = (dead_a.Kept, next(k for k in dead_a.kept if k.qualname == "DiscreteMeasure.symmetry_groups"))
    source = (tree_copy / rel).read_text(encoding="utf-8")
    assert _gains(dead_a, tree_copy, {rel: _cut_definition(source, "DiscreteMeasure.symmetry_groups")}, target) == {target}


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_every_dead_a_row_is_gone_and_every_kept_entry_stands")
def test_redefining_the_selector_reddens_exactly_its_row(dead_a: SimpleNamespace, tree_copy: Path) -> None:
    """The definition leg: a ``def`` of the selector's name appended to the copy of ``registry.py`` reddens the ``Retired`` row of
    that name and no other. Red input at the base tip: the precondition (the row reads red there, the definition exists)."""
    rel = "orpheus/numerics/quadrature/registry.py"
    row = dead_a.selector
    target = (dead_a.Retired, row)
    source = (tree_copy / rel).read_text(encoding="utf-8")
    assert _gains(dead_a, tree_copy, {rel: source + f"\n\ndef {row.qualname}():\n    return None\n"}, target) == {target}


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_every_dead_a_row_is_gone_and_every_kept_entry_stands")
def test_importing_a_retired_name_through_its_package_reddens_its_row(dead_a: SimpleNamespace, tree_copy: Path) -> None:
    """The reference leg through an ANCESTOR package (specification 1.5, decision E-D1): a planted test file importing the selector's
    error class from ``orpheus.numerics.quadrature`` reddens the row of that name and no other (the package
    ``__init__`` re-exports it today; after the unit it binds nothing of the name, and the leftover counts). Red input at the base tip:
    the precondition."""
    row = dead_a.selection_error
    target = (dead_a.Retired, row)
    plant = "tests/gates/numerics/test_planted_dead_a_reference.py"
    # built at run time from the row: this file spells no retired name, so the walk (which scans it) does not read its own plant
    text = f"from orpheus.numerics.quadrature import {row.qualname}\n\nvalue = {row.qualname}\n"
    assert _gains(dead_a, tree_copy, {plant: text}, target) == {target}


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_every_dead_a_row_is_gone_and_every_kept_entry_stands")
def test_a_moved_name_missing_at_its_home_reddens_its_row(dead_a: SimpleNamespace, tree_copy: Path) -> None:
    """The home leg: ``UNIT_INTERVAL`` (row 27, a ``Moved`` row) cut from its home ``tests/gates/geometry/_manifold_fixtures.py`` in the
    copy reddens that row and no other. Red input at the base tip: the precondition (the home does not exist yet, so the row reads red)."""
    home = "tests/gates/geometry/_manifold_fixtures.py"
    row = next(r for r in dead_a.rows if isinstance(r, dead_a.Moved) and r.home_name == "UNIT_INTERVAL")
    target = (dead_a.Moved, row)
    source = (tree_copy / home).read_text(encoding="utf-8") if (tree_copy / home).exists() else ""
    assert source, "precondition: the home exists in the copy"
    assert _gains(dead_a, tree_copy, {home: _cut_definition(source, "UNIT_INTERVAL")}, target) == {target}
