"""The mutations of unit ``roles`` (specification ``ontology_w0_spec.md`` 4.2 R1), built at run time from the rows of
``tests/gates/_retired/rows_w0_roles.py`` and read by the shared engine's walk.

The root is ``orpheus/numerics/operator.py`` alone, copied to its relative path in a temporary directory (the file the roles rows
bind). The unmutated walk of that root is required green first: a red reading before the mutation proves nothing (X1). Then:

* a class ``LinearOperator`` with a class attribute ``block_role`` is appended: the red set is exactly the ``Retired`` row
  ``LinearOperator.block_role`` and the ``RetiredName`` row ``block_role``;
* the line defining ``D2`` is deleted: the red set is exactly the ``Kept`` entry of ``D2``.

The two mutations are the entries' own fields, so a rows module whose rows are renamed reddens them; no planted source is stored.
This module is separate from ``test_retirements.py`` (the engine's own gate, not edited by a unit) and reads the same walk.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from tests.gates._retired import Kept, Retired, RetiredName
from tests.gates._retired import rows_w0_roles
from tests.gates._retired._walk import walk

pytestmark = pytest.mark.foundation

_REPO = Path(__file__).resolve().parents[2]
_OPERATOR = "orpheus/numerics/operator.py"


def _entries() -> list:
    return [*rows_w0_roles.ROWS, *rows_w0_roles.KEPT]


def _red(verdicts: dict) -> set:
    return {key for key, found in verdicts.items() if found}


def _copy_operator(root: Path) -> Path:
    target = root / _OPERATOR
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text((_REPO / _OPERATOR).read_text(encoding="utf-8"), encoding="utf-8")
    return target


def test_the_copy_of_operator_py_alone_reads_green_before_any_mutation(tmp_path: Path) -> None:
    """The control of both mutations below: 25 rows and the one ``Kept`` entry, none red on the unmutated copy."""
    _copy_operator(tmp_path)
    assert _red(walk(tmp_path, rows=_entries())) == set()


def test_a_block_role_planted_on_a_linear_operator_reddens_its_two_rows(tmp_path: Path) -> None:
    """Red set exactly ``{Retired(operator.py, LinearOperator.block_role), RetiredName(block_role)}``."""
    retired = next(r for r in rows_w0_roles.ROWS if type(r) is Retired and r.qualname == "LinearOperator.block_role")
    named = next(r for r in rows_w0_roles.ROWS if type(r) is RetiredName and r.name == "block_role")
    target = _copy_operator(tmp_path)
    assert _red(walk(tmp_path, rows=_entries())) == set()
    target.write_text(target.read_text(encoding="utf-8") + "\n\nclass LinearOperator:\n    block_role = None\n", encoding="utf-8")
    assert _red(walk(tmp_path, rows=_entries())) == {(Retired, retired), (RetiredName, named)}


def test_deleting_the_definition_of_d2_reddens_exactly_its_kept_entry(tmp_path: Path) -> None:
    """Red set exactly ``{Kept(operator.py, D2)}``."""
    kept = next(k for k in rows_w0_roles.KEPT if type(k) is Kept and k.qualname == "D2")
    target = _copy_operator(tmp_path)
    assert _red(walk(tmp_path, rows=_entries())) == set()
    lines = target.read_text(encoding="utf-8").split("\n")
    kept_lines = [i for i, line in enumerate(lines) if line.startswith(f"{kept.qualname} = ")]
    assert len(kept_lines) == 1
    del lines[kept_lines[0]]
    target.write_text("\n".join(lines), encoding="utf-8")
    assert _red(walk(tmp_path, rows=_entries())) == {(Kept, kept)}
