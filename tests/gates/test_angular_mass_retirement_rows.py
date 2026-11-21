r"""**#708** -- the five rows of ``tests/gates/_retired/rows_w0_angular_mass.py`` read green on the tree, and the engine can redden them.

The shared engine (``tests/gates/test_retirements.py``) asserts, on the real tree, only the wave-2 rows; a rows module of a later unit
is read by ``load_rows`` but nothing asserts it green. This file is that assertion for the rows of sub-unit ``angular-mass``
(``ontology_w0_spec.md`` 5.7.3 E6): the cross-class source dunder ``AngularSourceSink.__add__`` / ``__radd__`` / ``CONTAINS`` and
``BulkField.CONTAINS`` / ``BulkField.__add__`` are no longer defined and nothing imports, reads or sets them by that definition site.

* ``test_the_five_rows_read_green_on_the_tree``: RECORD of the tree. Red input in today's tree (the base tip ``954c4a4f``): all 5 rows
  are defined (findings at ``angular_source_sink.py`` lines 114, 116, 149 and ``_bases.py`` lines 373, 375), 5 of 5 red. Its positive
  control is the engine's own control row ``RetiredName('LinearOperator')``, which must read red (a walk that scanned no file reads every
  row green).
* ``test_a_planted_bulk_field_contains_reddens_its_row``: THEOREM of the walk on the rows. Mechanism: a ``.py``-only mini-root holding a copy of
  ``_bases.py`` at the same relative path; the unmutated root must read the five rows green (precondition; it fails today because
  ``_bases.py`` still defines both ``BulkField`` rows), the copy with ``PLANTED_BULK_FIELD_CONTAINS`` of ``_controls.py`` inserted after the end of
  the ``BulkField`` class docstring must put the ``BulkField.CONTAINS`` row in the red set, exactly that row. The planted text spells a real
  retired name, so it lives in ``tests/gates/_retired/_controls.py`` (excluded from the engine's scan), never here.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from tests.gates._retired import RetiredName, _controls, _walk, rows_w0_angular_mass

pytestmark = pytest.mark.foundation

_REPO = Path(__file__).resolve().parents[2]
_BASES = "orpheus/transport/fields/_bases.py"
#: the end of the ``BulkField`` class docstring; occurs once in ``_bases.py`` (measured at the base tip)
_DOCSTRING_END = '    ``ScalarFlux``, ...).\n    """\n'
_CONTROL = RetiredName("LinearOperator")  # a name the five trees use today: it must read red


def _verdicts(root: Path) -> dict:
    return _walk.walk(root, rows=[*rows_w0_angular_mass.ROWS, *rows_w0_angular_mass.KEPT, _CONTROL])


def _red_rows(verdicts: dict) -> set:
    return {key[1] for key, findings in verdicts.items() if findings and key[1] != _CONTROL}


def test_the_module_holds_five_retired_rows_and_no_kept():
    assert len(rows_w0_angular_mass.ROWS) == 5 and rows_w0_angular_mass.KEPT == ()


@pytest.mark.rests_on("tests/gates/test_retirements.py::test_the_real_walk_reads_every_row_of_the_package")
def test_the_five_rows_read_green_on_the_tree():
    verdicts = _verdicts(_REPO)
    assert verdicts[(RetiredName, _CONTROL)], "the control reads green: the walk scanned no file that uses it"
    red = _red_rows(verdicts)
    assert not red, f"{5 - len(red)} of 5 rows green; red: " + "; ".join(
        f"{r.qualname}: {verdicts[(type(r), r)][:2]}" for r in sorted(red)
    )


@pytest.mark.rests_on("tests/gates/test_angular_mass_retirement_rows.py::test_the_five_rows_read_green_on_the_tree")
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
