r"""**#708** -- a mutation of this unit's rows through the engine's seam (``tests.gates._retired._walk.walk(root, rows=...)``).

The engine (``tests/gates/test_retirements.py``) asserts every row of ``rows_w0_angular_mass.py`` green on the real tree; this file holds only
what the engine does not assert: the module's shape (``test_the_module_holds_five_retired_rows_and_no_kept``) and that the walk can
redden one of this unit's rows.

* ``test_a_planted_bulk_field_contains_reddens_its_row``: a ``.py``-only mini-root holding a copy of ``_bases.py`` at the same relative path;
  the unmutated root must read the five rows green (precondition), the copy with ``PLANTED_BULK_FIELD_CONTAINS`` of ``_controls.py`` inserted
  after the end of the ``BulkField`` class docstring must put the ``BulkField.CONTAINS`` row in the red set, exactly that row. The planted
  text spells a real retired name, so it lives in ``tests/gates/_retired/_controls.py`` (excluded from the engine's scan), never here.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from tests.gates._retired import _controls, _walk, rows_w0_angular_mass

pytestmark = pytest.mark.foundation

_REPO = Path(__file__).resolve().parents[2]
_BASES = "orpheus/transport/fields/_bases.py"
#: the end of the ``BulkField`` class docstring; occurs once in ``_bases.py`` (measured at the base tip)
_DOCSTRING_END = '    ``ScalarFlux``, ...).\n    """\n'


def _verdicts(root: Path) -> dict:
    return _walk.walk(root, rows=[*rows_w0_angular_mass.ROWS, *rows_w0_angular_mass.KEPT])


def _red_rows(verdicts: dict) -> set:
    return {key[1] for key, findings in verdicts.items() if findings}


def test_the_module_holds_five_retired_rows_and_no_kept():
    assert len(rows_w0_angular_mass.ROWS) == 5 and rows_w0_angular_mass.KEPT == ()


@pytest.mark.rests_on("tests/gates/test_retirements.py::test_every_row_and_kept_entry_of_the_package_is_green_on_the_tree")
def test_a_planted_bulk_field_contains_reddens_its_row(tmp_path):
    source = (_REPO / _BASES).read_text()
    assert source.count(_DOCSTRING_END) == 1, "the insertion anchor must occur exactly once"
    target = tmp_path / _BASES
    target.parent.mkdir(parents=True)
    target.write_text(source)
    unmutated = _red_rows(_verdicts(tmp_path))
    assert not unmutated, f"the unmutated root must read the five rows green first: {sorted(unmutated)}"
    target.write_text(source.replace(_DOCSTRING_END, _DOCSTRING_END + _controls.PLANTED_BULK_FIELD_CONTAINS))
    red = _red_rows(_verdicts(tmp_path))
    assert red == {r for r in rows_w0_angular_mass.ROWS if r.qualname == "BulkField.CONTAINS"}, sorted(red)
