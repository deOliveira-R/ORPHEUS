"""One table of the axis letters in production: :data:`orpheus.geometry.faces.AXIS_NAMES`.

Until layer DAG WP-B (2026-10-22) the point groups held a transcription of the
face naming's axis table (``AXIS_LETTER`` in the symmetry module, equal by
hand), and three of their constructors and one basis re-spelled the literal
inline. Two spellings of one convention agree only by transcription (X4 of the
instrument doctrine), so WP-B left one. This gate keeps it one: no module under
``orpheus/`` other than ``geometry/faces.py`` spells the letters as a sequence
literal. Tests are outside its population on purpose: a test listing the three
axes is an independent oracle, not a second production table.
"""

from __future__ import annotations

import pathlib
import re

import pytest

from orpheus.geometry.faces import AXIS_INDEX, AXIS_NAMES

ORPHEUS = pathlib.Path(__file__).resolve().parents[3] / "orpheus"

#: A tuple, list or set literal of the three letters, either quote style, with
#: or without a trailing comma.
_TABLE = re.compile(r"""[\(\[\{]\s*(['"])x\1\s*,\s*(['"])y\2\s*,\s*(['"])z\3\s*,?\s*[\)\]\}]""")


def _spellings(root: pathlib.Path) -> list[str]:
    return sorted(
        f"{path.relative_to(root)}:{i}"
        for path in root.rglob("*.py")
        for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1)
        if _TABLE.search(line)
    )


@pytest.mark.foundation
def test_the_pattern_finds_each_literal_shape() -> None:
    """Positive control (X1): every shape the census claims to see, it sees."""
    for literal in ('("x", "y", "z")', "['x', 'y', 'z']", '{"x","y","z",}', 'Literal["x", "y", "z"]'):
        assert _TABLE.search(literal), literal
    assert not _TABLE.search('("x", "y")') and not _TABLE.search('"xyz"')


@pytest.mark.foundation
def test_production_spells_the_axis_letters_once() -> None:
    spellings = _spellings(ORPHEUS)
    assert [s.split(":")[0] for s in spellings] == ["geometry/faces.py"], spellings


@pytest.mark.foundation
def test_the_inverse_is_the_table_read_backwards() -> None:
    assert AXIS_NAMES == ("x", "y", "z")
    assert all(AXIS_NAMES[AXIS_INDEX[a]] == a for a in AXIS_NAMES)
    assert sorted(AXIS_INDEX.values()) == list(range(len(AXIS_NAMES)))
