r"""The gates of the shared retirement engine (specification ``ontology_w0_spec.md`` section 1.5, unit ``roles``, first slice).

The engine is ``tests/gates/_retired/`` (the five row types in ``__init__``, the definition finder ``_finder``, the walk
``_walk``) and its data are the ``rows_*.py`` modules of that package. This file is its done-when: each test below names the input
that reddens it (X1) and each is RED at the base tip because the engine does not exist there (the fixture :func:`engine` fails with
``ImportError``, ``cannot import name '_retired'``, on every node; a module-level import would read as one collection error and ``FAILED = 0``).

**The seam this file pins (choices of the gate author, section 1.5 not fixing the spelling).** A builder conforms to these names:

* ``tests.gates._retired``: the NamedTuples ``Retired``, ``Moved``, ``RetiredName``, ``RetiredModule``, ``Kept`` (fields as 1.5).
* ``tests.gates._retired._finder.definitions(source: str) -> dict[str, int]``: qualname to line of every definition of ``source``.
* ``tests.gates._retired._walk.load_rows(rows_dir: Path | None = None) -> tuple[tuple, tuple]``: the pair ``(rows, kept)``, the ``ROWS``
  and the ``KEPT`` of every ``rows_*.py`` of ``rows_dir`` (default: the package directory, resolved from ``_walk.__file__``).
* ``tests.gates._retired._walk.walk(root: Path, rows: Iterable | None = None, rows_dir: Path | None = None) -> dict``: the per-row
  verdicts. ``rows`` (rows and ``Kept`` entries mixed) overrides the loader. The key of a verdict is
  ``(type(row), row)`` (decision C12: ``RetiredName('x') == RetiredModule('x')`` as tuples); the value is a list of
  ``(file, line)`` pairs, ``file`` the posix path relative to ``root``, one pair per finding (decision C9: one finding per
  ``(row, file, line)``), ``line`` ``0`` for a finding that has no line (an absent file, a home lacking its name, a module file
  that exists).
* ``tests.gates._retired._walk.GateError``: raised by ``walk`` for an invalid row or entry, with the message fragments ``empty``,
  ``.py`` and ``duplicate`` (case-insensitive), and ``unparsable`` for a file under the root that does not parse (decision C10).

**Not pinned, by choice.** The planted sources of the controls are built at run time here from synthetic names (``p``, ``m``,
``NAME``) and, for the mutations, from the fields of the rows of ``rows_wave2``, so this file spells no retired name
(``_controls.py`` of 1.5 is then not needed by any gate). A mutation is a mini-root built at run time, the copy to which the row's
fields are applied: the real tree is walked once (done-when 2).

Kinds (``vv-principles``): the controls and mutations are THEOREMs of the walk (each leg reads a shape or it does not); the real-tree
rows are RECORDs of the tree today (green once the wave-2 retirement landed); the structural census (done-when 7) is a THEOREM of
the engine's source.

Ladder (``rests_on``): the finder and the gate errors are the foundations; the controls rest on the finder; each mutation rests
on the control of its row type; the real-tree reading rests on every control.
"""

from __future__ import annotations

import ast
import importlib.util
import inspect
import shutil
import textwrap
from pathlib import Path
from types import SimpleNamespace

import pytest

_REPO = Path(__file__).resolve().parents[2]
_THIS = "tests/gates/test_retirements.py"
_WAVE2_TEST = _REPO / "tests" / "gates" / "test_wave2_retirement.py"
_ROW_TYPE_NAMES = ("Retired", "Moved", "RetiredName", "RetiredModule")


@pytest.fixture(scope="module")
def engine() -> SimpleNamespace:
    """The engine, imported inside the fixture so that a missing engine reddens each test on its own line."""
    from tests.gates import _retired
    from tests.gates._retired import _finder, _walk, rows_wave2

    return SimpleNamespace(
        R=_retired, finder=_finder, walk=_walk.walk, load_rows=_walk.load_rows, load_rows_by_module=_walk.load_rows_by_module, GateError=_walk.GateError,
        walk_module=_walk, rows_wave2=rows_wave2, package_dir=Path(_walk.__file__).resolve().parent,
    )


# ---------------------------------------------------------------------------
# Planting: run-time synthetic roots, the planted lines computed from the planted text
# ---------------------------------------------------------------------------


class Plant:
    """A synthetic root. ``add(rel, text)`` writes ``text`` (a line ending ``  #*`` is a planted line); ``planted`` is the
    set of ``(rel, line)`` pairs of the marked lines, computed from the text and never typed."""

    MARK = "#*"

    def __init__(self, root: Path) -> None:
        self.root = root
        self.planted: set[tuple[str, int]] = set()

    def add(self, rel: str, text: str) -> None:
        lines = textwrap.dedent(text).strip("\n").split("\n")
        out = []
        for number, line in enumerate(lines, start=1):
            if line.rstrip().endswith(self.MARK):
                self.planted.add((rel, number))
                line = line.rstrip()[: -len(self.MARK)].rstrip()
            out.append(line)
        path = self.root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(out) + "\n", encoding="utf-8")


def _findings(verdicts: dict, key: tuple) -> list[tuple[str, int]]:
    if key not in verdicts:
        raise AssertionError(f"no verdict for {key!r}; keys: {sorted(map(repr, verdicts))}")
    return [tuple(f) for f in verdicts[key]]


def _assert_exactly(found: list[tuple[str, int]], planted: set[tuple[str, int]]) -> None:
    """One finding per planted line: set equality AND no pair reported twice (never a count of names, C9)."""
    assert len(found) == len(set(found)), f"a (file, line) reported twice: {sorted(found)}"
    assert set(found) == planted, f"missing {sorted(planted - set(found))}; unexpected {sorted(set(found) - planted)}"


_M = "orpheus/p/m.py"


def _positive_retired(plant: Plant) -> None:
    """The nine planted statements of the positive control of a ``Retired`` row ``(orpheus/p/m.py, NAME)`` (1.5), under ``tests/``
    and ``orpheus/p/``; ``orpheus/p/m.py`` defines nothing, so only references can redden the row."""
    plant.add(_M, '"""Defines nothing."""')
    plant.add("tests/u1.py", """
        import orpheus.p.m as m
        m.NAME  #*
        getattr(m, "NAME")  #*
        """)
    plant.add("tests/u2.py", """
        from orpheus.p import m
        m.NAME  #*
        """)
    plant.add("tests/u3.py", """
        import orpheus.p.m as m
        def test(monkeypatch):
            monkeypatch.setattr(m, "NAME", f)  #*
        """)
    plant.add("tests/u4.py", "from orpheus.p.m import NAME  #*")
    plant.add("tests/u5.py", "from orpheus.p import NAME  #*")
    plant.add("tests/u6.py", 'FIXTURE = "from orpheus.p.m import NAME\\n"  #*')
    plant.add("orpheus/p/__init__.py", """
        from .m import NAME  #*
        __all__ = ["NAME"]  #*
        """)


def _negative_retired(plant: Plant) -> None:
    """Statements that must read nothing (1.5): a homonym ancestor, a different object, an unbound receiver, ``hasattr`` on a receiver
    that is not the row's module, a fixture string of another module."""
    plant.add(_M, '"""Defines nothing."""')
    plant.add("orpheus/p/other.py", "NAME = 1")
    plant.add("orpheus/p/__init__.py", "from .other import NAME")
    plant.add("tests/n1.py", """
        import numpy as np
        NAME = 1
        f(NAME=True)
        np.polynomial.laguerre
        """)
    plant.add("tests/n2.py", """
        import orpheus.p.m as m
        import numpy as np
        hasattr(np, "NAME")
        getattr(x, "NAME")
        from q import NAME
        S = "from q import NAME\\n"
        """)
    plant.add("tests/n3.py", "from orpheus.p import NAME")


@pytest.fixture
def retired_row(engine: SimpleNamespace):
    return engine.R.Retired(_M, "NAME")


# ---------------------------------------------------------------------------
# Foundation: the finder, the gate errors, the verdict key
# ---------------------------------------------------------------------------


@pytest.mark.foundation
def test_the_finder_reads_every_definition_shape_it_must(engine: SimpleNamespace) -> None:
    r"""The one finder reads a ``def``, a ``class``, a class attribute, a module ``Assign``, an ``AnnAssign``, a tuple target
    (C1), a nested qualname (C2), a definition under ``if``/``try``/``with`` (C3) and ``self.x = ...`` in a method (C5), and
    reads nothing for a name a function merely assigns. Red input: any engine whose finder drops one shape (1.5, "five shapes").
    The expected map is computed from the planted text: line numbers are found, not typed."""
    source = textwrap.dedent("""
        def fn(): pass
        class K:
            attr = 1
            def meth(self):
                self.slot = 2
                local = 3
            class Inner: pass
        top = 1
        ann: int = 2
        a, b = 1, 2
        if True:
            under_if = 1
        try:
            under_try = 1
        except Exception:
            pass
        with open:
            under_with = 1
        """).strip("\n")
    found = engine.finder.definitions(source)
    lines = source.split("\n")

    def at(fragment: str) -> int:
        return next(i for i, line in enumerate(lines, start=1) if fragment in line)

    expected = {
        "fn": at("def fn"), "K": at("class K"), "K.attr": at("attr = 1"), "K.meth": at("def meth"), "K.slot": at("self.slot"),
        "K.Inner": at("class Inner"), "top": at("top = 1"), "ann": at("ann: int"), "a": at("a, b"), "b": at("a, b"),
        "under_if": at("under_if"), "under_try": at("under_try"), "under_with": at("under_with"),
    }
    assert found == expected, (sorted(set(found) ^ set(expected)), found)
    assert "local" not in found and "K.meth.local" not in found


def _rows_dir_with(tmp_path: Path, body: str) -> Path:
    rows_dir = tmp_path / "rows_dir"
    rows_dir.mkdir()
    (rows_dir / "rows_planted.py").write_text(
        "from tests.gates._retired import Kept, Moved, Retired, RetiredModule, RetiredName\n" + textwrap.dedent(body), encoding="utf-8"
    )
    return rows_dir


_VALID_ROWS = """
    ROWS = (Retired('orpheus/p/m.py', 'f'), RetiredName('g'), RetiredModule('orpheus.p.q'),
            Moved('orpheus/p/m.py', 'h', 'orpheus/p/home.py', 'h'))
    KEPT = (Kept('orpheus/p/m.py', 'k', 'a wave'),)
    """

_INVALID_ROWS = {
    "empty name": ("RetiredName('')", "empty"),
    "empty qualname": ("Retired('orpheus/p/m.py', '')", "empty"),
    "empty module": ("RetiredModule('')", "empty"),
    "empty home_name": ("Moved('orpheus/p/m.py', 'f', 'orpheus/p/h.py', '')", "empty"),
    "file not py": ("Retired('orpheus/p/m.txt', 'f')", r"\.py"),
    "path leading space": ("Retired(' orpheus/p/m.py', 'f')", "repository path"),
    "path absolute": ("Retired('/abs/orpheus/p/m.py', 'f')", "repository path"),
    "path dot slash": ("Retired('./orpheus/p/m.py', 'f')", "repository path"),
    "path outside the five trees": ("Retired('scratch/m.py', 'f')", "repository path"),
    "home path leading space": ("Moved('orpheus/p/m.py', 'f', ' orpheus/p/h.py', 'f')", "repository path"),
    "kept inside ROWS": ("Kept('orpheus/p/m.py', 'f', 'a wave')", "ROWS"),
    "retired and moved of one definition": ("Retired('orpheus/p/m.py', 'f'), Moved('orpheus/p/m.py', 'f', 'orpheus/p/h.py', 'f')", "duplicate"),
    "home not py": ("Moved('orpheus/p/m.py', 'f', 'orpheus/p/h', 'f')", r"\.py"),
    "duplicate row": ("RetiredName('g'), RetiredName('g')", "duplicate"),
    "duplicate moved": ("Moved('orpheus/p/m.py', 'f', 'orpheus/p/h.py', 'f'), Moved('orpheus/p/m.py', 'f', 'orpheus/p/h.py', 'f')", "duplicate"),
}


@pytest.mark.foundation
@pytest.mark.parametrize("kind", sorted(_INVALID_ROWS))
def test_an_invalid_row_in_a_rows_module_is_a_gate_error(engine: SimpleNamespace, tmp_path: Path, kind: str) -> None:
    r"""An empty field, a ``file`` or ``home`` that does not end in ``.py`` and a duplicate cannot be made unspellable on a
    ``NamedTuple`` (``__post_init__`` is never called, measured), so the walk raises ``GateError``: shown once per kind (done-when 6),
    each through a rows module read from a rows directory. Positive control: the valid module of the same shape walks clean (a
    ``GateError`` that fires on everything is not a gate)."""
    row, fragment = _INVALID_ROWS[kind]
    root = tmp_path / "root"
    root.mkdir()
    (tmp_path / "good").mkdir()
    (tmp_path / "bad").mkdir()
    engine.walk(root, rows_dir=_rows_dir_with(tmp_path / "good", _VALID_ROWS))  # the control: valid rows, no error
    bad = _rows_dir_with(tmp_path / "bad", f"ROWS = ({row},)\nKEPT = ()\n")
    with pytest.raises(engine.GateError, match=f"(?i){fragment}"):
        engine.walk(root, rows_dir=bad)


@pytest.mark.foundation
def test_an_empty_keeper_a_kept_collision_and_an_unparsable_file_are_gate_errors(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""The remaining gate errors of 1.5: an empty ``keeper``; a ``RetiredName`` whose name is the last component of a ``Kept``
    qualname; a file under the root that does not parse (C10: the engine of the review skipped it, a silent zero). Positive
    control: the same walk on a parsable root with a ``Kept`` entry and a distinct name reads clean."""
    R = engine.R
    root = tmp_path / "root"
    Plant(root).add("tests/ok.py", "x = 1")
    Plant(root).add("orpheus/p/k.py", "def kept_fn(): pass")
    engine.walk(root, rows=[R.Kept("orpheus/p/k.py", "kept_fn", "a wave"), R.RetiredName("other")])
    with pytest.raises(engine.GateError, match=r"(?i)empty"):
        engine.walk(root, rows=[R.Kept("orpheus/p/k.py", "kept_fn", "")])
    with pytest.raises(engine.GateError, match=r"(?i)kept"):
        engine.walk(root, rows=[R.Kept("orpheus/p/k.py", "C.kept_fn", "a wave"), R.RetiredName("kept_fn")])
    broken = tmp_path / "broken"
    Plant(broken).add("tests/ok.py", "x = 1")
    (broken / "tests" / "bad.py").write_text("def (:\n", encoding="utf-8")
    with pytest.raises(engine.GateError, match=r"(?i)unparsable"):
        engine.walk(broken, rows=[R.RetiredName("anything")])


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_finder_reads_every_definition_shape_it_must")
def test_a_verdict_key_carries_the_type_of_its_row(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""``RetiredName('x') == RetiredModule('x')`` as tuples (C12), so a mapping keyed by the bare row would merge them. Input: both
    rows in one list; the walk returns two keys, ``(type, row)`` each. Red input: an engine keyed by the bare tuple returns
    one key. (A ``Retired`` and a ``Kept`` of equal fields are no second control: they are refused as one definition named twice.)"""
    R = engine.R
    verdicts = engine.walk(tmp_path, rows=[R.RetiredName("x"), R.RetiredModule("x")])
    assert set(verdicts) == {(engine.R.RetiredName, R.RetiredName("x")), (engine.R.RetiredModule, R.RetiredModule("x"))}


# ---------------------------------------------------------------------------
# Controls: one finding per planted line, positive and negative on separate roots
# ---------------------------------------------------------------------------


@pytest.mark.foundation
def test_the_kept_and_the_rows_are_typed_apart_and_one_definition_is_one_entry(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""``load_rows`` returns ``(rows, kept)``: a ``Retired`` inside a ``KEPT`` is a gate error, and a ``Kept`` and a ``Retired`` of one
    ``(file, qualname)`` (a gate permanently red with a misleading verdict) are refused as one definition named twice. Positive
    control: a valid module returns the two tuples typed apart."""
    R = engine.R
    (tmp_path / "good").mkdir()
    rows, kept = engine.load_rows(_rows_dir_with(tmp_path / "good", _VALID_ROWS))
    assert len(rows) == 4 and len(kept) == 1 and all(isinstance(k, R.Kept) for k in kept)
    bodies = {
        "retired inside KEPT": ("ROWS = ()\nKEPT = (Retired('orpheus/p/m.py', 'f'),)\n", "KEPT"),
    }
    for label, (body, fragment) in bodies.items():
        with pytest.raises(engine.GateError, match=fragment):
            (tmp_path / label.replace(" ", "_")).mkdir()
            engine.load_rows(_rows_dir_with(tmp_path / label.replace(" ", "_"), body))
    root = tmp_path / "root"
    root.mkdir()
    with pytest.raises(engine.GateError, match=r"(?i)duplicate"):
        engine.walk(root, rows=[R.Retired("orpheus/p/m.py", "f"), R.Kept("orpheus/p/m.py", "f", "a wave")])


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_finder_reads_every_definition_shape_it_must")
@pytest.mark.parametrize("moved", [False, True], ids=["retired", "moved"])
def test_the_positive_control_of_a_definition_row_finds_each_of_nine_planted_lines(
    engine: SimpleNamespace, tmp_path: Path, moved: bool,
) -> None:
    r"""The nine statements of 1.5 (an ``Attribute`` and a ``getattr`` through ``import p.m as m``, an ``Attribute`` through
    ``from p import m``, the ``monkeypatch.setattr`` setter, ``from p.m import NAME``, the ancestor import ``from p import NAME``, the
    fixture string, the re-export ``from .m import NAME`` and the ``__all__`` member) are each reported once, for a ``Retired`` row
    and for a ``Moved`` row (whose home defines ``NAME``, so the home leg adds nothing). Red input: an engine without the
    ``Attribute`` leg (9 first reds today that no other leg finds), without the ancestor leg, the setter leg, the string leg or
    the ``__all__`` leg. Compared as set equality over planted ``(file, line)`` with no pair twice."""
    plant = Plant(tmp_path)
    _positive_retired(plant)
    if moved:
        plant.add("orpheus/p/home.py", "NAME = 1")
        row = engine.R.Moved(_M, "NAME", "orpheus/p/home.py", "NAME")
    else:
        row = engine.R.Retired(_M, "NAME")
    assert len(plant.planted) == 9
    verdicts = engine.walk(tmp_path, rows=[row])
    _assert_exactly(_findings(verdicts, (type(row), row)), plant.planted)


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_finder_reads_every_definition_shape_it_must")
@pytest.mark.parametrize("moved", [False, True], ids=["retired", "moved"])
def test_the_negative_control_of_a_definition_row_reads_nothing(engine: SimpleNamespace, tmp_path: Path, moved: bool) -> None:
    r"""On its own root (the negative needs ``p/__init__.py`` to bind ``NAME`` to ``.other``, the positive to ``.m``): a module-level
    ``NAME = 1`` and ``f(NAME=True)`` in a file importing nothing from ``p.m``, ``np.polynomial.laguerre`` against a row
    ``(m.py, laguerre)``, ``getattr(x, "NAME")`` with ``x`` unbound, ``hasattr(np, "NAME")`` with ``np`` another module, ``from q import NAME`` and its fixture
    string with ``q`` not an ancestor, and an ancestor that re-exports a DIFFERENT object under the row's name (the homonym). Red
    input: an engine that matches an import by spelling alone, reads ``hasattr`` off any receiver, or counts every ancestor import."""
    plant = Plant(tmp_path)
    _negative_retired(plant)
    rows: list = [engine.R.Retired(_M, "NAME"), engine.R.Retired(_M, "laguerre")]
    if moved:
        plant.add("orpheus/p/home.py", "NAME = 1")
        rows = [engine.R.Moved(_M, "NAME", "orpheus/p/home.py", "NAME"), engine.R.Retired(_M, "laguerre")]
    verdicts = engine.walk(tmp_path, rows=rows)
    for row in rows:
        _assert_exactly(_findings(verdicts, (type(row), row)), set())


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_finder_reads_every_definition_shape_it_must")
def test_the_leftover_of_a_retirement_is_found_through_a_package_that_binds_nothing(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""The tenth control (decision E-D1, 1.5): ``from p import NAME`` where ``p/__init__.py`` binds nothing named ``NAME`` is the
    leftover of a partial retirement and is reported (1 of 1). Red input: the earlier rule ("an ancestor counts only when it binds
    the name to the row's module") reads 0 here, exactly the state the leg exists to catch (``[M]`` 1.5)."""
    plant = Plant(tmp_path)
    plant.add(_M, '"""Defines nothing."""')
    plant.add("orpheus/p/__init__.py", '"""Binds nothing."""')
    plant.add("tests/left.py", "from orpheus.p import NAME  #*")
    row = engine.R.Retired(_M, "NAME")
    _assert_exactly(_findings(engine.walk(tmp_path, rows=[row]), (engine.R.Retired, row)), plant.planted)


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_finder_reads_every_definition_shape_it_must")
def test_an_ancestor_chain_is_read_through_an_absolute_link_and_a_homonym_through_one(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""Decision C18 (links read in both forms): the module is ``orpheus/p/sub/m.py``; ``orpheus/p/__init__.py`` re-exports ``NAME``
    by the ABSOLUTE ``from orpheus.p.sub import NAME`` and ``orpheus/p/sub/__init__.py`` by the relative ``from .m import NAME``,
    so ``from orpheus.p.sub import NAME`` (the ancestor ``p.sub`` binds the row's object) and ``from orpheus.p import NAME``
    both reach the row's object (each reported). On a second root the absolute link is
    ``from orpheus.p.other import NAME`` (a different object): not reported. Red input: an engine that reads only relative links
    reads the second chain as unbound-or-other for both roots."""
    row = engine.R.Retired("orpheus/p/sub/m.py", "NAME")
    reach = Plant(tmp_path / "reach")
    reach.add("orpheus/p/sub/m.py", '"""Defines nothing."""')
    reach.add("orpheus/p/sub/__init__.py", "from .m import NAME  #*")
    reach.add("orpheus/p/__init__.py", "from orpheus.p.sub import NAME  #*")
    reach.add("tests/user.py", "from orpheus.p import NAME  #*")
    _assert_exactly(_findings(engine.walk(tmp_path / "reach", rows=[row]), (engine.R.Retired, row)), reach.planted)
    other = Plant(tmp_path / "other")
    other.add("orpheus/p/sub/m.py", '"""Defines nothing."""')
    other.add("orpheus/p/other.py", "NAME = 1")
    other.add("orpheus/p/__init__.py", "from orpheus.p.other import NAME")
    other.add("tests/user.py", "from orpheus.p import NAME")
    _assert_exactly(_findings(engine.walk(tmp_path / "other", rows=[row]), (engine.R.Retired, row)), set())


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_finder_reads_every_definition_shape_it_must")
def test_the_string_leg_reads_a_relative_fixture_import_by_its_suffix_for_a_definition_row_only(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""Decision C22 and the string leg of 1.5: the fixture string ``"from ..m import NAME"`` (a relative import whose module is a
    dotted suffix of the row's module at component boundaries) is found for a ``Retired`` row, and on a row of ``p/xm.py`` (the
    row's module ends in the characters ``m`` but not in the component ``.m``) it is not, and a ``RetiredName`` row never reads a fixture string (another gate's planted
    control must be able to spell a retired name). Red input: an engine with no suffix rule, a character suffix, or a
    ``RetiredName`` walk over strings."""
    plant = Plant(tmp_path)
    plant.add(_M, '"""Defines nothing."""')
    plant.add("tests/s1.py", 'S = "from ..m import NAME\\n"  #*')
    rows = [engine.R.Retired(_M, "NAME"), engine.R.RetiredName("NAME")]
    verdicts = engine.walk(tmp_path, rows=rows)
    _assert_exactly(_findings(verdicts, (engine.R.Retired, rows[0])), plant.planted)
    _assert_exactly(_findings(verdicts, (engine.R.RetiredName, rows[1])), set())
    # the boundary: on a row of ``orpheus/p/xm.py`` the string ``from ..m import NAME`` names ``m``, and ``orpheus.p.xm`` ends in the
    # characters ``m`` but not in the component ``.m`` (the character rule is red here, the boundary rule reads nothing)
    other = Plant(tmp_path / "boundary")
    other.add("orpheus/p/xm.py", '"""Defines nothing."""')
    other.add("tests/s2.py", 'S = "from ..m import NAME\\n"')
    row = engine.R.Retired("orpheus/p/xm.py", "NAME")
    _assert_exactly(_findings(engine.walk(tmp_path / "boundary", rows=[row]), (engine.R.Retired, row)), set())


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_finder_reads_every_definition_shape_it_must")
def test_a_definition_in_the_rows_own_file_is_the_first_leg_and_an_absent_file_defines_nothing(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""Leg (a): ``class C: NAME = 1`` in the row's file reddens ``(file, "C.NAME")`` at the line of the assignment, ``(file, "NAME")``
    reads nothing there, and a row whose file is absent from the root is green (1.5: a one-file mutation root is legal for a
    ``Retired`` row). Red input: an engine that reads names and not definition sites."""
    plant = Plant(tmp_path)
    plant.add(_M, """
        class C:
            NAME = 1  #*
        """)
    nested, flat, absent = engine.R.Retired(_M, "C.NAME"), engine.R.Retired(_M, "NAME"), engine.R.Retired("orpheus/p/missing.py", "C.NAME")
    verdicts = engine.walk(tmp_path, rows=[nested, flat, absent])
    _assert_exactly(_findings(verdicts, (engine.R.Retired, nested)), plant.planted)
    _assert_exactly(_findings(verdicts, (engine.R.Retired, flat)), set())
    _assert_exactly(_findings(verdicts, (engine.R.Retired, absent)), set())


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_finder_reads_every_definition_shape_it_must")
def test_the_home_leg_is_green_for_a_home_that_defines_the_name_and_red_for_one_that_lacks_it_or_is_absent(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""Decision C11 / E-D2: a ``Moved`` row is red on a root that lacks its ``home`` (an absent home reads as a home that lacks
    ``home_name``). Three roots, each its own ``tmp_path``: home defines ``home_name`` (green), home defines something else (red,
    ``(home, 0)``), home absent (red, ``(home, 0)``). Red input: an engine with no home leg, or one that reads an absent home as
    green."""
    row = engine.R.Moved(_M, "f", "orpheus/p/home.py", "f")
    for label, home_text, expected in (
        ("defines", "def f(): pass", set()), ("lacks", "def g(): pass", {("orpheus/p/home.py", 0)}), ("absent", None, {("orpheus/p/home.py", 0)}),
    ):
        root = tmp_path / label
        plant = Plant(root)
        plant.add(_M, '"""Defines nothing."""')
        if home_text is not None:
            plant.add("orpheus/p/home.py", home_text)
        _assert_exactly(_findings(engine.walk(root, rows=[row]), (engine.R.Moved, row)), expected)


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_finder_reads_every_definition_shape_it_must")
def test_the_positive_control_of_a_name_row_finds_each_of_seven_planted_lines_and_its_negative_none(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""``RetiredName('NAME')``: ``getattr(x, "NAME")``, ``x.NAME``, ``__all__ = ["NAME"]``, ``monkeypatch.setattr(m, "NAME", f)``,
    ``from m import NAME``, ``NAME = 1`` and ``f(NAME=True)`` are each found (7 of 7, a ``Store`` and a keyword included). Its
    negative on its own root (the fixture string ``"from m import NAME"``, a bare ``"NAME"`` in a tuple) reads nothing; the ``hasattr``
    string is read as an identifier like ``getattr``'s (``test_hasattr_reads_...`` below). Red input: an engine that reads the strings
    of a fixture, or misses a store, a keyword or an ``__all__`` member."""
    positive = Plant(tmp_path / "pos")
    positive.add("tests/rn.py", """
        getattr(x, "NAME")  #*
        x.NAME  #*
        __all__ = ["NAME"]  #*
        monkeypatch.setattr(m, "NAME", f)  #*
        from m import NAME  #*
        NAME = 1  #*
        f(NAME=True)  #*
        """)
    row = engine.R.RetiredName("NAME")
    assert len(positive.planted) == 7
    _assert_exactly(_findings(engine.walk(tmp_path / "pos", rows=[row]), (engine.R.RetiredName, row)), positive.planted)
    negative = Plant(tmp_path / "neg")
    negative.add("tests/rn.py", """
        S = "from m import NAME\\n"
        T = ("NAME", "other")
        """)
    _assert_exactly(_findings(engine.walk(tmp_path / "neg", rows=[row]), (engine.R.RetiredName, row)), set())


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_positive_control_of_a_name_row_finds_each_of_seven_planted_lines_and_its_negative_none")
def test_a_name_row_reads_the_binding_forms_that_are_not_name_nodes(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""``RetiredName('NAME')`` also reads the stores that are not ``Name`` nodes: ``except E as NAME``, ``global NAME``, the module
    of ``from NAME import x``, a ``match`` capture, ``__all__ += ["NAME"]`` and an annotated ``__all__: list[str] = ["NAME"]``
    (6 of 6, one root each line of text). Red input: an engine that reads ``Name`` nodes and keywords only."""
    plant = Plant(tmp_path)
    plant.add("tests/rn.py", """
        try:
            pass
        except Exception as NAME:  #*
            pass
        def f():
            global NAME  #*
        from NAME import x  #*
        match 1:
            case NAME:  #*
                pass
        __all__ = []
        __all__ += ["NAME"]  #*
        __all__: list[str] = ["NAME"]  #*
        """)
    row = engine.R.RetiredName("NAME")
    _assert_exactly(_findings(engine.walk(tmp_path, rows=[row]), (engine.R.RetiredName, row)), plant.planted)


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_positive_control_of_a_definition_row_finds_each_of_nine_planted_lines")
def test_a_local_name_bound_by_two_imports_reads_both_bindings(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""C6 read strictly: ``from orpheus.p import m`` at the top and ``from q import m`` inside a function bind one local name to two
    targets, and ``m.NAME`` is read against both, so the row is red (1 of 1). Red input: an engine whose alias map keeps the last
    binding (``ast.walk`` visits the function body last) reads 0."""
    plant = Plant(tmp_path)
    plant.add(_M, '"""Defines nothing."""')
    plant.add("tests/two.py", """
        from orpheus.p import m
        m.NAME  #*
        def f():
            from q import m
        """)
    row = engine.R.Retired(_M, "NAME")
    _assert_exactly(_findings(engine.walk(tmp_path, rows=[row]), (engine.R.Retired, row)), plant.planted)


def _module_positive(plant: Plant) -> None:
    plant.add("tests/a.py", """
        import orpheus.p.m  #*
        import orpheus.p.m.sub  #*
        """)
    plant.add("tests/b.py", """
        from typing import TYPE_CHECKING
        if TYPE_CHECKING:
            from orpheus.p.m import x  #*
            __import__("orpheus.p.m")  #*
        def f():
            from orpheus.p import m  #*
            import importlib
            importlib.import_module("orpheus.p.m")  #*
        """)
    plant.add("orpheus/p/user.py", "from . import m  #*")


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_finder_reads_every_definition_shape_it_must")
def test_the_positive_control_of_a_module_row_finds_each_of_seven_planted_lines_and_its_negative_none(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""``RetiredModule('orpheus.p.m')``: ``import p.m``, ``import p.m.sub``, ``from p.m import x``, ``from p import m``, a relative
    ``from . import m`` in a file of the package, ``importlib.import_module("p.m")`` and ``__import__("p.m")``, spread over the
    three positions (top level, function body, ``TYPE_CHECKING``), 7 of 7. Its negative on its own root (``import p.m2``,
    ``import p.mx``, the string ``"p.m"`` as a tuple member) reads nothing. Red input: an engine that matches by character prefix
    (``m2``, ``mx``), reads no relative import, misses a position, or has no string leg."""
    positive = Plant(tmp_path / "pos")
    _module_positive(positive)
    row = engine.R.RetiredModule("orpheus.p.m")
    assert len(positive.planted) == 7
    _assert_exactly(_findings(engine.walk(tmp_path / "pos", rows=[row]), (engine.R.RetiredModule, row)), positive.planted)
    negative = Plant(tmp_path / "neg")
    negative.add("tests/n.py", """
        import orpheus.p.m2
        import orpheus.p.mx
        T = ("orpheus.p.m", "other")
        """)
    _assert_exactly(_findings(engine.walk(tmp_path / "neg", rows=[row]), (engine.R.RetiredModule, row)), set())




@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_finder_reads_every_definition_shape_it_must")
def test_the_legs_with_no_other_control_find_each_of_their_planted_lines(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""The legs of 1.5 that no control above plants: the dotted-qualname setter (``getattr(C, "attr")`` and
    ``monkeypatch.setattr(C, "attr", v)`` with ``C`` imported from the row's module, a row ``(m.py, "C.attr")``), the dotted-target
    setter (``monkeypatch.setattr("orpheus.p.m.NAME", v)``), the C19 chain (``import orpheus.p.m`` then ``orpheus.p.m.NAME``) and
    the submodule string (``importlib.import_module("orpheus.p.m.sub")``). Negatives on the same root: ``getattr(D, "attr")`` for a
    class ``D`` of another module, and a homonym class ``C`` re-exported by the ancestor package (``from .other import C``), which
    must not redden the row. Red input: an engine that lacks one of the four legs reads 0 of that leg's planted lines."""
    plant = Plant(tmp_path)
    plant.add(_M, '"""Defines nothing."""')
    plant.add("tests/d1.py", """
        from orpheus.p.m import C
        getattr(C, "attr")  #*
        def test(monkeypatch):
            monkeypatch.setattr(C, "attr", 1)  #*
        """)
    plant.add("tests/d2.py", """
        from other import D
        getattr(D, "attr")
        """)
    plant.add("tests/d3.py", """
        def test(monkeypatch):
            monkeypatch.setattr("orpheus.p.m.NAME", 1)  #*
        """)
    plant.add("tests/d4.py", """
        import orpheus.p.m
        orpheus.p.m.NAME  #*
        """)
    plant.add("tests/d5.py", """
        import importlib
        importlib.import_module("orpheus.p.m.sub")  #*
        """)
    rows = [engine.R.Retired(_M, "C.attr"), engine.R.Retired(_M, "NAME")]
    verdicts = engine.walk(tmp_path, rows=rows)
    expected_dotted = {("tests/d1.py", 2), ("tests/d1.py", 4)}
    _assert_exactly(_findings(verdicts, (engine.R.Retired, rows[0])), expected_dotted)
    _assert_exactly(_findings(verdicts, (engine.R.Retired, rows[1])), {("tests/d3.py", 2), ("tests/d4.py", 2)})
    module_row = engine.R.RetiredModule("orpheus.p.m")
    _assert_exactly(_findings(engine.walk(tmp_path, rows=[module_row]), (engine.R.RetiredModule, module_row)),
                    {("tests/d1.py", 1), ("tests/d4.py", 1), ("tests/d5.py", 2), (_M, 0)})
    homonym = Plant(tmp_path / "homonym")
    homonym.add(_M, '"""Defines nothing."""')
    homonym.add("orpheus/p/other.py", "class C: pass")
    homonym.add("orpheus/p/__init__.py", "from .other import C")
    homonym.add("tests/h.py", 'from orpheus.p import C\ngetattr(C, "attr")\n')
    _assert_exactly(_findings(engine.walk(tmp_path / "homonym", rows=[rows[0]]), (engine.R.Retired, rows[0])), set())


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_finder_reads_every_definition_shape_it_must")
def test_the_ancestor_binding_follows_every_link_and_a_dead_link_counts(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""Decision E-D1 on the shapes a package can bind a name by. Each root plants the row ``(orpheus/p/m.py, NAME)``, a package
    ``orpheus/p/__init__.py``, a consumer ``tests/u.py: from orpheus.p import NAME`` and reports (1) a link through a sibling module
    that is gone (``from .api import NAME``, ``api.py`` empty: the leftover of a retirement, reported at the consumer), (2) a chain
    through a sibling that re-exports the row (``api.py: from .m import NAME``: reported at the consumer and at ``api.py``), (3) the
    explicit re-export idiom ``from .m import NAME as NAME`` (reported at the consumer and the package line), (4) a link to a sibling
    that defines its own ``NAME`` (a different object: nothing reported), (5) a function-body import in the package, which is no module-level binding (the consumer is
    reported), (6) an external link ``from numpy import NAME`` (another object: nothing). Red input: an engine that reads a sibling
    chain, an ``as`` that repeats the name or a nested import as a different object reads 0 for roots 1 to 3 and 5."""
    row = engine.R.Retired(_M, "NAME")

    def root(label: str, init: str, extra: dict[str, str]) -> tuple[Path, Plant]:
        path = tmp_path / label
        plant = Plant(path)
        plant.add(_M, '"""Defines nothing."""')
        plant.add("orpheus/p/__init__.py", init)
        for rel, text in extra.items():
            plant.add(rel, text)
        plant.add("tests/u.py", "from orpheus.p import NAME  #*")
        return path, plant

    cases = [
        ("dead_link", "from .api import NAME", {"orpheus/p/api.py": '"""Empty."""'}, False),
        ("sibling_chain", "from .api import NAME", {"orpheus/p/api.py": "from .m import NAME  #*"}, False),
        ("as_same_name", "from .m import NAME as NAME  #*", {}, False),
        ("own_definition", "from .api import NAME", {"orpheus/p/api.py": "NAME = 1"}, True),
        ("external", "from numpy import NAME", {}, True),
        ("nested_only", "def f():\n    from .other import NAME", {"orpheus/p/other.py": "NAME = 1"}, False),
    ]
    for label, init, extra, reads_nothing in cases:
        path, plant = root(label, init, extra)
        found = _findings(engine.walk(path, rows=[row]), (engine.R.Retired, row))
        _assert_exactly(found, set() if reads_nothing else plant.planted)


# The controls of one rule for when an ancestor package counts (decision E-D1), each a function of the walk it reads so that the
# mutation tests below run the very same control against a mutated copy of the engine.

_ANCESTOR_USES = {
    "attribute_import_as": "import orpheus.p as g\ng.NAME  #*",
    "attribute_from_import": "from orpheus import p\np.NAME  #*",
    "getattr_through_ancestor": 'import orpheus.p as g\ngetattr(g, "NAME")  #*',
    "setattr_through_ancestor": 'import orpheus.p as g\ndef test(monkeypatch):\n    monkeypatch.setattr(g, "NAME", v)  #*',
}


def _control_ancestor_use(walk, R, tmp_path: Path, use: str, homonym: bool) -> None:
    """The row ``(orpheus/p/m.py, NAME)``; the ancestor package ``orpheus.p`` re-exports ``NAME`` from the row's module (or, with
    ``homonym``, defines a different ``NAME`` of its own) and a consumer reads ``NAME`` off the PACKAGE. Reported (one finding per
    planted line, the package's own import included) or, for the homonym, nothing."""
    plant = Plant(tmp_path)
    plant.add(_M, '"""Defines nothing."""')
    if homonym:
        plant.add("orpheus/p/other.py", "NAME = 1")
        plant.add("orpheus/p/__init__.py", "from .other import NAME")
    else:
        plant.add("orpheus/p/__init__.py", "from .m import NAME  #*")
    plant.add("tests/u.py", _ANCESTOR_USES[use])
    row = R.Retired(_M, "NAME")
    _assert_exactly(_findings(walk(tmp_path, rows=[row]), (R.Retired, row)), set() if homonym else plant.planted)


def _control_absent_link_module(walk, R, tmp_path: Path) -> None:
    """``orpheus/p/__init__.py: from .api import NAME`` with NO ``api.py`` under the root: the link leads to a module of the five
    trees that is gone, the binding is "none", the ancestor counts and the consumer's leftover import is reported."""
    plant = Plant(tmp_path)
    plant.add(_M, '"""Defines nothing."""')
    plant.add("orpheus/p/__init__.py", "from .api import NAME")
    plant.add("tests/u.py", "from orpheus.p import NAME  #*")
    row = R.Retired(_M, "NAME")
    _assert_exactly(_findings(walk(tmp_path, rows=[row]), (R.Retired, row)), plant.planted)


def _control_try_star(walk, R, tmp_path: Path) -> None:
    """The ancestor ``__init__`` binds ``NAME`` to a DIFFERENT object inside an ``except*`` handler at module level: the binding is
    "other" (the homonym), so the consumer's ``from orpheus.p import NAME`` reads nothing. An engine that does not descend into
    ``ast.TryStar`` sees no binding, counts the ancestor, and reports the consumer."""
    plant = Plant(tmp_path)
    plant.add(_M, '"""Defines nothing."""')
    plant.add("orpheus/p/other.py", "NAME = 1")
    plant.add("orpheus/p/__init__.py", "try:\n    pass\nexcept* ValueError:\n    from .other import NAME")
    plant.add("tests/u.py", "from orpheus.p import NAME  #*")
    row = R.Retired(_M, "NAME")
    _assert_exactly(_findings(walk(tmp_path, rows=[row]), (R.Retired, row)), set())


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_finder_reads_every_definition_shape_it_must")
@pytest.mark.parametrize("use", sorted(_ANCESTOR_USES))
def test_an_undotted_name_read_off_an_ancestor_package_is_reported_and_a_homonym_ancestor_is_not(
    engine: SimpleNamespace, tmp_path: Path, use: str
) -> None:
    r"""One rule for when an ancestor counts, on every leg (E-D1): ``g.NAME`` (``import orpheus.p as g`` or ``from orpheus import p``),
    ``getattr(g, "NAME")`` and ``monkeypatch.setattr(g, "NAME", v)`` with ``g`` the package that re-exports the row's ``NAME`` are
    each reported (2 findings per root: the use and the package's import); the same uses where the package defines its own,
    different ``NAME`` read nothing. Red input: an attribute leg or an undotted setter leg that reads the row's own module alone
    reads 1 of the 2 planted lines (``test_a_walk_that_restores_the_own_module_alone_...`` runs that mutation)."""
    _control_ancestor_use(engine.walk, engine.R, tmp_path / "positive", use, homonym=False)
    _control_ancestor_use(engine.walk, engine.R, tmp_path / "homonym", use, homonym=True)


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_finder_reads_every_definition_shape_it_must")
def test_a_link_to_a_module_file_that_is_absent_binds_nothing_and_an_except_star_block_is_read(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""Two shapes of ``_Tree.binding``: a chain link ``from .api import NAME`` whose ``api.py`` does not exist is "none" (the
    existing dead-link case plants an EMPTY ``api.py``, which takes the other branch), and a binding inside ``except*`` is read.
    Red input: classing the absent module as "other" reads 0 for the first; not descending into ``ast.TryStar`` reports the consumer
    of the second."""
    _control_absent_link_module(engine.walk, engine.R, tmp_path / "absent")
    _control_try_star(engine.walk, engine.R, tmp_path / "star")


def _control_hasattr(walk, R, tmp_path: Path) -> None:
    """``hasattr(g, "NAME")`` reads a retired attribute by a string and answers False once it is gone. Positive root: ``g`` is the
    row's module by import (a ``Retired`` row reports the line), and a ``RetiredName`` row reads the string as an identifier, as for
    ``getattr``. Negative root: ``hasattr`` on another module and on an unbound receiver reads nothing for the ``Retired`` row."""
    plant = Plant(tmp_path / "pos")
    plant.add(_M, '"""Defines nothing."""')
    plant.add("tests/h.py", 'import orpheus.p.m as g\nhasattr(g, "NAME")  #*')
    retired, named = R.Retired(_M, "NAME"), R.RetiredName("NAME")
    verdicts = walk(tmp_path / "pos", rows=[retired, named])
    _assert_exactly(_findings(verdicts, (R.Retired, retired)), plant.planted)
    _assert_exactly(_findings(verdicts, (R.RetiredName, named)), plant.planted)
    negative = Plant(tmp_path / "neg")
    negative.add(_M, '"""Defines nothing."""')
    negative.add("tests/h.py", 'import numpy as np\nhasattr(np, "NAME")\nhasattr(x, "NAME")')
    _assert_exactly(_findings(walk(tmp_path / "neg", rows=[retired]), (R.Retired, retired)), set())


def _control_class_qualified(walk, R, tmp_path: Path) -> None:
    """A row ``(m.py, "Owner.member")`` is read by ``<expr>.member`` where ``<expr>`` resolves by the file's imports to ``p.m.Owner`` or
    to ``<ancestor>.Owner``: the class by a from-import, by a module alias, by the plain dotted chain, and through an ancestor package
    that re-exports it (4 positives). Negatives: a homonym ``Owner`` of another module, and an instance read ``x.member`` with ``x`` not
    resolvable by an import (the documented limit)."""
    plant = Plant(tmp_path / "pos")
    plant.add(_M, '"""Defines nothing."""')
    plant.add("orpheus/p/__init__.py", "from .m import Owner")
    plant.add("tests/c1.py", "from orpheus.p.m import Owner\nOwner.member()  #*")
    plant.add("tests/c2.py", "import orpheus.p.m as mm\nmm.Owner.member  #*")
    plant.add("tests/c3.py", "import orpheus.p.m\norpheus.p.m.Owner.member  #*")
    plant.add("tests/c4.py", "from orpheus.p import Owner\nOwner.member  #*")
    row = R.Retired(_M, "Owner.member")
    _assert_exactly(_findings(walk(tmp_path / "pos", rows=[row]), (R.Retired, row)), plant.planted)
    negative = Plant(tmp_path / "neg")
    negative.add(_M, '"""Defines nothing."""')
    negative.add("orpheus/p/other.py", "class Owner: pass")
    negative.add("orpheus/p/__init__.py", "from .other import Owner")
    negative.add("tests/n1.py", "from orpheus.p.other import Owner\nOwner.member")
    negative.add("tests/n2.py", "from orpheus.p import Owner\nOwner.member")
    negative.add("tests/n3.py", "from orpheus.p.m import Owner\nx = Owner()\nx.member")
    _assert_exactly(_findings(walk(tmp_path / "neg", rows=[row]), (R.Retired, row)), set())


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_positive_control_of_a_definition_row_finds_each_of_nine_planted_lines")
def test_hasattr_reads_a_retired_attribute_by_its_string(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""``hasattr(g, "NAME")`` is a reader by string beside ``getattr``, ``setattr`` and ``delattr`` (``hasattr`` returns False silently
    once the attribute is gone, so an unread one hides a leftover). Red input: an engine whose string accessors omit ``hasattr``."""
    _control_hasattr(engine.walk, engine.R, tmp_path)


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_positive_control_of_a_definition_row_finds_each_of_nine_planted_lines")
def test_a_dotted_row_is_read_by_a_class_qualified_attribute(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""``Owner.member`` read off the class resolved by import (module, alias or ancestor), reported; a homonym class and an instance
    read are not. Red input: an engine that reads a dotted row by its definition and its string setters alone."""
    _control_class_qualified(engine.walk, engine.R, tmp_path)


def _mutated_walk(engine: SimpleNamespace, tmp_path: Path, old: str, new: str):
    """``walk`` of a COPY of ``_walk.py`` under ``tmp_path`` with ``old`` replaced by ``new`` (the tracked file is never touched)."""
    source = Path(engine.walk_module.__file__).read_text(encoding="utf-8")
    assert source.count(old) == 1, f"the mutation site must occur once, found {source.count(old)}: {old!r}"
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "_walk_mutant.py"
    path.write_text(source.replace(old, new), encoding="utf-8")
    spec = importlib.util.spec_from_file_location("_walk_mutant", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.walk


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_an_undotted_name_read_off_an_ancestor_package_is_reported_and_a_homonym_ancestor_is_not")
@pytest.mark.parametrize("use", sorted(_ANCESTOR_USES))
def test_a_walk_that_restores_the_own_module_alone_misses_the_ancestor_reads(engine: SimpleNamespace, tmp_path: Path, use: str) -> None:
    """Mutation on a copy: the readers of the attribute and undotted setter legs revert to ``{module}``; each positive reddens."""
    walk = _mutated_walk(engine, tmp_path / "copy", "readers = {module, *ancestors}", "readers = {module}")
    with pytest.raises(AssertionError, match="missing"):
        _control_ancestor_use(walk, engine.R, tmp_path / "run", use, homonym=False)


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_hasattr_reads_a_retired_attribute_by_its_string")
def test_a_walk_without_hasattr_among_its_string_accessors_misses_the_hasattr_read(engine: SimpleNamespace, tmp_path: Path) -> None:
    """Mutation on a copy: ``hasattr`` removed from the string accessors; the positive of the hasattr control reddens."""
    walk = _mutated_walk(engine, tmp_path / "copy", '"delattr", "hasattr")', '"delattr")')
    with pytest.raises(AssertionError, match="missing"):
        _control_hasattr(walk, engine.R, tmp_path / "run")


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_a_dotted_row_is_read_by_a_class_qualified_attribute")
def test_a_walk_that_reads_no_class_qualified_attribute_misses_the_dotted_read(engine: SimpleNamespace, tmp_path: Path) -> None:
    """Mutation on a copy: the attribute leg is switched off for a dotted row; the class-qualified positives redden."""
    walk = _mutated_walk(engine, tmp_path / "copy", "if owners & f.resolve(value):", "if not dotted_qualname and owners & f.resolve(value):")
    with pytest.raises(AssertionError, match="missing"):
        _control_class_qualified(walk, engine.R, tmp_path / "run")


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_a_link_to_a_module_file_that_is_absent_binds_nothing_and_an_except_star_block_is_read")
def test_mutations_of_the_absent_link_and_the_except_star_descent_redden_their_controls(engine: SimpleNamespace, tmp_path: Path) -> None:
    """Mutations on copies: an absent link module classed "other" misses the leftover import; ``ast.TryStar`` dropped from the
    descent reports the consumer of the homonym."""
    absent = _mutated_walk(engine, tmp_path / "c1", 'return "none" if module.split(".")[0] in SCAN_TREES else "other"', 'return "other"')
    with pytest.raises(AssertionError, match="missing"):
        _control_absent_link_module(absent, engine.R, tmp_path / "r1")
    star = _mutated_walk(engine, tmp_path / "c2", "(ast.If, ast.Try, ast.With, ast.TryStar)", "(ast.If, ast.Try, ast.With)")
    with pytest.raises(AssertionError, match="unexpected"):
        _control_try_star(star, engine.R, tmp_path / "r2")


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_finder_reads_every_definition_shape_it_must")
@pytest.mark.parametrize("planted_rel", ["orpheus/p/m.py", "orpheus/p/m/__init__.py", "orpheus/p/m/stray.py"])
def test_the_absence_leg_reports_a_leftover_module_file_or_package_and_is_silent_on_a_root_with_neither(
    engine: SimpleNamespace, tmp_path: Path, planted_rel: str,
) -> None:
    r"""Decision E-D6: ``RetiredModule('orpheus.p.m')`` is red when ``orpheus/p/m.py`` or ``orpheus/p/m/__init__.py`` exists under
    the root, or any ``.py`` under the directory ``orpheus/p/m/`` (``stray.py`` with no ``__init__.py``: a namespace package, which the
    deleted homes test, ``path.exists()`` on the directory, read red), though nothing imports it (1 of 1 each, one root each); a root
    holding none is silent. Red input: an engine
    with the import legs only (the wave-2 shape of the check lived in a separate path test that the port deletes)."""
    row = engine.R.RetiredModule("orpheus.p.m")
    planted = Plant(tmp_path / "planted")
    planted.add(planted_rel, '"""A leftover."""  #*')
    found = _findings(engine.walk(tmp_path / "planted", rows=[row]), (engine.R.RetiredModule, row))
    assert len(found) == 1 and found[0][0] == planted_rel, found
    clean = Plant(tmp_path / "clean")
    clean.add("orpheus/p/other.py", "x = 1")
    _assert_exactly(_findings(engine.walk(tmp_path / "clean", rows=[row]), (engine.R.RetiredModule, row)), set())


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_finder_reads_every_definition_shape_it_must")
def test_a_kept_entry_is_green_while_defined_and_red_when_its_definition_or_file_is_gone(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""``Kept(file, qualname, keeper)`` is red when its definition is gone or its file is absent from the root (D1', 1.5). Three
    roots: defined (green), file present without the definition (red, ``(file, 0)``), file absent (red). Red input: an engine
    that does not read ``Kept`` entries, or reads a ``Kept`` whose file is absent as green."""
    entry = engine.R.Kept("orpheus/p/k.py", "g", "a wave")
    for label, text, expected in (("defined", "def g(): pass", set()), ("gone", "def h(): pass", {("orpheus/p/k.py", 0)}), ("absent", None, {("orpheus/p/k.py", 0)})):
        plant = Plant(tmp_path / label)
        plant.add("tests/benign.py", "x = 1")
        if text is not None:
            plant.add("orpheus/p/k.py", text)
        _assert_exactly(_findings(engine.walk(tmp_path / label, rows=[entry]), (engine.R.Kept, entry)), expected)


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_positive_control_of_a_definition_row_finds_each_of_nine_planted_lines")
def test_the_walk_reads_the_five_trees_and_not_its_own_files_or_the_rows_directory(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""The scan population (1.5, decision E-D3): ``orpheus``, ``tests``, ``tools``, ``examples`` and ``student_resources``, minus
    ``tests/gates/_retired/`` and this file. A user planted in any of the five trees is read (5 of 5 reported); the same user
    planted at the root outside the five trees, under ``tests/gates/_retired/`` or as ``tests/gates/test_retirements.py`` is not
    (0 of 3: a package planted at the root is the reason the control of 1.5 reported 7 of 9). Red input: an engine that scans
    every ``.py`` under the root, or its own planted sources."""
    plant = Plant(tmp_path)
    plant.add(_M, '"""Defines nothing."""')
    user = "from orpheus.p.m import NAME"
    for tree in ("orpheus", "tests", "tools", "examples", "student_resources"):
        plant.add(f"{tree}/leak.py", user + "  #*")
    for rel in ("leak_outside/leak.py", "tests/gates/_retired/leak.py", _THIS):
        plant.add(rel, user)
    row = engine.R.Retired(_M, "NAME")
    _assert_exactly(_findings(engine.walk(tmp_path, rows=[row]), (engine.R.Retired, row)), plant.planted)


# ---------------------------------------------------------------------------
# Done-when 4: one mutation per row type, exactly the mutated row red, after the same walk read it green
# ---------------------------------------------------------------------------


def _first(engine: SimpleNamespace, kind: type):
    return next(r for r in engine.rows_wave2.ROWS if type(r) is kind)


def _mutation_rows(engine: SimpleNamespace) -> list:
    """Every row of ``rows_wave2`` plus the two types it has none of, so that a mutation is read against 65 + 3 entries."""
    R = engine.R
    return [
        *engine.rows_wave2.ROWS,
        R.Retired("orpheus/q/rt.py", "Owner.attr"),
        R.Moved("orpheus/q/mv.py", "fn", "orpheus/q/home.py", "fn"),
        R.Kept("orpheus/q/kp.py", "alpha", "a consumer"),
        R.Kept("orpheus/q/kp.py", "beta", "a consumer"),
    ]


def _base_root(engine: SimpleNamespace, root: Path) -> Plant:
    """The copy: the ``Retired`` row's file (its owner class without the attribute), the
    ``Moved`` file and its home, the ``Kept`` file, and a benign test."""
    plant = Plant(root)
    plant.add("orpheus/q/rt.py", "class Owner:\n    pass")
    plant.add("orpheus/q/mv.py", "x = 1")
    plant.add("orpheus/q/home.py", "def fn(): pass")
    plant.add("orpheus/q/kp.py", "def alpha(): pass\ndef beta(): pass")
    plant.add("tests/benign.py", "x = 1")
    return plant


def _red_keys(verdicts: dict) -> set:
    return {key for key, found in verdicts.items() if found}


def _mutate_and_read(engine: SimpleNamespace, tmp_path: Path, mutate) -> tuple[set, dict]:
    """The unmutated walk (all entries green, required before the mutation, decision E-D2) then the mutated walk on the same root."""
    _base_root(engine, tmp_path)
    rows = _mutation_rows(engine)
    before = engine.walk(tmp_path, rows=rows)
    if _red_keys(before):
        raise AssertionError(f"the unmutated copy is not green (a reading red before the mutation proves nothing): {sorted(map(repr, _red_keys(before)))}")
    mutate(tmp_path)
    after = engine.walk(tmp_path, rows=rows)
    return _red_keys(after), after


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_positive_control_of_a_definition_row_finds_each_of_nine_planted_lines")
def test_a_definition_planted_in_the_copy_turns_exactly_its_retired_row_red(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""Mutation of a ``Retired`` row: the copy's file of the synthetic ``Retired`` row gains the definition the row names (a class
    attribute), mechanism: a mini-root copy, mutated text, walk on that root. Red set
    must be exactly that row. Red input of the unmutated walk: none (required green first)."""
    row = next(r for r in _mutation_rows(engine) if type(r) is engine.R.Retired and r.file == "orpheus/q/rt.py")

    def mutate(root: Path) -> None:
        (root / row.file).write_text("class Owner:\n    attr = None\n", encoding="utf-8")

    red, after = _mutate_and_read(engine, tmp_path, mutate)
    assert red == {(engine.R.Retired, row)}, sorted(map(repr, red))
    assert (row.file, 2) in {tuple(f) for f in after[(engine.R.Retired, row)]}


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_positive_control_of_a_name_row_finds_each_of_seven_planted_lines_and_its_negative_none")
def test_a_use_planted_in_the_copy_turns_exactly_its_name_row_red(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""Mutation of a ``RetiredName`` row: a test file of the copy assigns the row's name (built from the row's field). Red set
    exactly that row; the planted file is among its findings."""
    row = _first(engine, engine.R.RetiredName)

    def mutate(root: Path) -> None:
        (root / "tests" / "planted_use.py").write_text(f"{row.name} = 1\n", encoding="utf-8")

    red, after = _mutate_and_read(engine, tmp_path, mutate)
    assert red == {(engine.R.RetiredName, row)}, sorted(map(repr, red))
    assert ("tests/planted_use.py", 1) in {tuple(f) for f in after[(engine.R.RetiredName, row)]}


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_positive_control_of_a_module_row_finds_each_of_seven_planted_lines_and_its_negative_none")
def test_an_import_planted_in_the_copy_turns_exactly_its_module_row_red(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""Mutation of a ``RetiredModule`` row: a test file of the copy imports the row's module; and, separately on a fresh copy, a
    file at the module's path appears (the absence leg). Each turns exactly that row red."""
    row = _first(engine, engine.R.RetiredModule)

    def by_import(root: Path) -> None:
        (root / "tests" / "planted_use.py").write_text(f"import {row.module}\n", encoding="utf-8")

    red, after = _mutate_and_read(engine, tmp_path / "import", by_import)
    assert red == {(engine.R.RetiredModule, row)}, sorted(map(repr, red))
    assert ("tests/planted_use.py", 1) in {tuple(f) for f in after[(engine.R.RetiredModule, row)]}

    def by_file(root: Path) -> None:
        target = root / (row.module.replace(".", "/") + ".py")
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("x = 1\n", encoding="utf-8")

    red, _ = _mutate_and_read(engine, tmp_path / "file", by_file)
    assert red == {(engine.R.RetiredModule, row)}, sorted(map(repr, red))


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_home_leg_is_green_for_a_home_that_defines_the_name_and_red_for_one_that_lacks_it_or_is_absent")
def test_a_moved_row_goes_red_when_its_old_definition_returns_and_when_its_home_loses_the_name(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""Mutation of a ``Moved`` row (the root holds its ``file`` AND its ``home``, decision E-D2): (i) the definition returns to
    the old file; (ii) on a fresh copy, the home loses the name. Each turns exactly the ``Moved`` row red."""
    row = next(r for r in _mutation_rows(engine) if type(r) is engine.R.Moved)

    def returns(root: Path) -> None:
        (root / row.file).write_text(f"def {row.qualname}(): pass\n", encoding="utf-8")

    red, _ = _mutate_and_read(engine, tmp_path / "returns", returns)
    assert red == {(engine.R.Moved, row)}, sorted(map(repr, red))

    def loses(root: Path) -> None:
        (root / row.home).write_text("y = 1\n", encoding="utf-8")

    red, _ = _mutate_and_read(engine, tmp_path / "loses", loses)
    assert red == {(engine.R.Moved, row)}, sorted(map(repr, red))


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_a_kept_entry_is_green_while_defined_and_red_when_its_definition_or_file_is_gone")
def test_deleting_a_kept_definition_in_the_copy_turns_exactly_that_kept_entry_red(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""Mutation of a ``Kept`` entry: the copy's file loses the definition of one of two entries; the red ``Kept`` set is exactly that
    entry (the unmutated copy's is empty, required first). Red input: an engine whose ``Kept`` leg reads the file and not the
    qualname."""
    entries = [r for r in _mutation_rows(engine) if type(r) is engine.R.Kept]

    def mutate(root: Path) -> None:
        (root / entries[0].file).write_text(f"def {entries[1].qualname}(): pass\n", encoding="utf-8")

    red, _ = _mutate_and_read(engine, tmp_path, mutate)
    assert red == {(engine.R.Kept, entries[0])}, sorted(map(repr, red))


# ---------------------------------------------------------------------------
# Done-when 5: the rows are the package's, never the root's
# ---------------------------------------------------------------------------


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_a_verdict_key_carries_the_type_of_its_row")
def test_the_default_rows_are_read_from_the_package_directory_and_never_from_the_root(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""Decision E-D4: the walk on a root that holds no rows module (an empty root, and a root that plants a ``rows_planted.py`` under
    ``tests/gates/_retired/``) reads the full row count of the package directory, counted here by importing each ``rows_*.py`` of
    the package directory by file location (a count independent of ``load_rows``), and none of the planted row. Red input: an engine
    that globs under ``root`` reads 0 rows on both roots. The count is also asserted non-zero and at least wave 2's 65."""
    modules = sorted(engine.package_dir.glob("rows_*.py"))
    assert modules, engine.package_dir
    expected = 0
    for path in modules:
        spec = importlib.util.spec_from_file_location(f"_independent_{path.stem}", path)
        assert spec is not None and spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        import sys

        sys.modules[spec.name] = module
        try:
            spec.loader.exec_module(module)
        finally:
            sys.modules.pop(spec.name, None)
        expected += len(module.ROWS) + len(module.KEPT)
    assert expected >= 65, expected
    empty = tmp_path / "empty"
    empty.mkdir()
    assert len(engine.walk(empty)) == expected
    planted = tmp_path / "planted"
    Plant(planted).add("tests/gates/_retired/rows_planted.py", "from tests.gates._retired import RetiredName\nROWS = (RetiredName('only_in_the_root'),)\nKEPT = ()\n")
    verdicts = engine.walk(planted)
    assert len(verdicts) == expected
    assert (engine.R.RetiredName, engine.R.RetiredName("only_in_the_root")) not in verdicts


# ---------------------------------------------------------------------------
# Done-when 2: the real tree, wave 2
# ---------------------------------------------------------------------------


_REAL_CONTROL = "LinearOperator"  # a retired-name row that is RED on the real tree: a name the five trees use today


def _assert_every_row_green(engine: SimpleNamespace, root: Path, rows_dir: Path | None = None, verdicts: dict | None = None) -> int:
    """The one gate over a tree: every row of every ``rows_*.py`` of ``rows_dir`` (default: the package) and every ``Kept`` entry is
    green on ``root``. The failure message lists each red entry with the file name of its rows module and its first findings.
    ``verdicts`` is a walk already made of these entries (the real-tree fixture, which carries a control row besides). Returns the
    number of entries read."""
    by_module = engine.load_rows_by_module(rows_dir)
    entries = [(name, entry) for name, (rows, kept) in by_module.items() for entry in (*rows, *kept)]
    if verdicts is None:
        verdicts = engine.walk(root, rows_dir=rows_dir)
    red = [(name, entry, verdicts[(type(entry), entry)]) for name, entry in entries if verdicts[(type(entry), entry)]]
    assert not red, f"{len(entries) - len(red)} of {len(entries)} entries green; red: " + "; ".join(
        f"{name}: {entry}: {found[:3]}" for name, entry, found in red
    )
    return len(entries)


@pytest.fixture(scope="module")
def real_walk(engine: SimpleNamespace) -> dict:
    """One walk of the real tree with the default rows and one control row, ``RetiredName(_REAL_CONTROL)``, that must read red: a walk
    that scanned no file would read every row green, and the control is the row it cannot read green (X1)."""
    rows, kept = engine.load_rows()
    return engine.walk(_REPO, rows=[*rows, *kept, engine.R.RetiredName(_REAL_CONTROL)])


@pytest.mark.foundation
@pytest.mark.rests_on(
    f"{_THIS}::test_the_positive_control_of_a_definition_row_finds_each_of_nine_planted_lines",
    f"{_THIS}::test_the_positive_control_of_a_name_row_finds_each_of_seven_planted_lines_and_its_negative_none",
    f"{_THIS}::test_the_positive_control_of_a_module_row_finds_each_of_seven_planted_lines_and_its_negative_none",
    f"{_THIS}::test_the_absence_leg_reports_a_leftover_module_file_or_package_and_is_silent_on_a_root_with_neither",
    f"{_THIS}::test_hasattr_reads_a_retired_attribute_by_its_string",
    f"{_THIS}::test_a_dotted_row_is_read_by_a_class_qualified_attribute",
)
def test_every_row_and_kept_entry_of_the_package_is_green_on_the_tree(engine: SimpleNamespace, real_walk: dict) -> None:
    r"""Every row of every rows module of the package (wave 2 and each unit's ``rows_w0_<unit>.py``) and every ``Kept`` entry read green on
    the real tree: ``k of N`` is printed. RECORD of the tree. The wave-2 statement stays: 65 rows (2 ``Retired``, 60 ``RetiredName``, 3
    ``RetiredModule``; ``KEPT = ()``, decisions D5 and R2), the three module rows being the three paths the deleted path test named
    (decision E-D6), each module's file and package directory absent. The input count stands beside the finding count (X1): the walk
    holds one verdict per entry and the control, which reads red (a walk that scanned no file would read it green). Red input: a name,
    module or definition of any row still present anywhere in the five trees; the control below plants one in a rows module."""
    R = engine.R
    rows = list(engine.rows_wave2.ROWS)
    by_type = {name: sum(1 for r in rows if type(r).__name__ == name) for name in _ROW_TYPE_NAMES}
    assert by_type == {"Retired": 2, "Moved": 0, "RetiredName": 60, "RetiredModule": 3}, by_type
    assert tuple(engine.rows_wave2.KEPT) == ()
    assert real_walk[(R.RetiredName, R.RetiredName(_REAL_CONTROL))], "the control reads green: the walk scanned no file that uses it"
    total = _assert_every_row_green(engine, _REPO, verdicts=real_walk)
    assert total and len(real_walk) == total + 1, (len(real_walk), total)
    print(f"rows green: {total} of {total}; wave-2 rows green: {len(rows)} of {len(rows)}")
    for row in rows:
        if isinstance(row, R.RetiredModule):
            base = _REPO / row.module.replace(".", "/")
            assert not base.with_suffix(".py").exists() and not (base / "__init__.py").exists(), row


_PLANTED_HEAD = "from tests.gates._retired import Kept, Moved, Retired, RetiredModule, RetiredName\n"


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_every_row_and_kept_entry_of_the_package_is_green_on_the_tree")
@pytest.mark.parametrize(
    "body",
    ["ROWS = (RetiredName('LinearOperator'),)\nKEPT = ()\n", "ROWS = ()\nKEPT = (Kept('orpheus/p/not_a_module.py', 'not_defined', 'a wave'),)\n"],
    ids=["a red row", "a kept definition that does not exist"],
)
def test_the_real_tree_gate_fails_on_a_planted_rows_module_and_names_it(engine: SimpleNamespace, tmp_path: Path, body: str) -> None:
    r"""The gate body, called with a rows directory holding a copy of ``rows_wave2.py`` (green: 65 of 65, no raise) and a planted
    ``rows_w0_planted.py`` whose one entry is red on the real tree, fails and names ``rows_w0_planted.py`` and not ``rows_wave2.py``.
    Red input: a gate that asserts the rows of wave 2 alone."""
    rows_dir = tmp_path / "rows"
    rows_dir.mkdir()
    shutil.copy(engine.package_dir / "rows_wave2.py", rows_dir / "rows_wave2.py")
    assert _assert_every_row_green(engine, _REPO, rows_dir) == 65
    (rows_dir / "rows_w0_planted.py").write_text(_PLANTED_HEAD + body, encoding="utf-8")
    with pytest.raises(AssertionError, match=r"rows_w0_planted\.py") as raised:
        _assert_every_row_green(engine, _REPO, rows_dir)
    assert "rows_wave2.py" not in str(raised.value), str(raised.value)


# ---------------------------------------------------------------------------
# The port of the wave-2 file (the edits of 1.5 that this file's gates see)
# ---------------------------------------------------------------------------

_DELETED_NAMES = frozenset({
    "retired_uses", "_findings", "_package_of", "_python_files", "_THIS", "_ORPHEUS", "RETIRED_PATHS", "RETIRING_NAMES",
    "_OUTSIDE_TREES", "_TEST_HOMONYMS", "_HOMONYMS", "_PLANTED", "RETIRED_MODULES", "forbidden_imports",
    "test_nothing_under_orpheus_names_the_retired_vocabulary", "test_no_test_example_or_tool_names_the_retired_vocabulary",
    "test_every_excluded_homonym_survives_where_it_is_said_to", "test_the_census_sees_every_shape_it_must",
    "test_the_census_reads_nothing_in_a_clean_source", "test_the_walk_finds_a_planted_module_and_misses_a_clean_one",
    "test_the_three_retiring_homes_do_not_exist",
})
_KEPT_TESTS = frozenset({
    "test_the_capture_scripts_are_deleted_and_their_fixtures_stay", "test_the_cold_import_list_of_the_layer_gate_names_no_retired_module",
    "test_the_sn_method_space_carries_no_index_sets",
})


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_every_row_and_kept_entry_of_the_package_is_green_on_the_tree")
def test_the_wave2_file_keeps_only_its_checks_that_are_not_name_walks() -> None:
    r"""Section 1.5: the wave-2 file's own name and module walks, their planted controls, the homes test (decision E-D6) and the
    names only they read are deleted, three tests stay (the capture scripts, the cold-import list, the SN method space), its
    docstring says the census is the shared engine's, and no kept test names a deleted test in a ``rests_on`` line (the five lines
    that did would redden ``test_every_rests_on_target_resolves_but_the_listed_rows``). Population: the top-level definitions and
    assignments of the file, by ``ast``. Red input at the base tip: all 21 deleted names are defined there."""
    tree = ast.parse(_WAVE2_TEST.read_text(encoding="utf-8"))
    defined: set[str] = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            defined.add(node.name)
        elif isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            defined |= {t.id for t in targets if isinstance(t, ast.Name)}
        elif isinstance(node, ast.ImportFrom):
            defined |= {a.asname or a.name for a in node.names}
    assert not (defined & _DELETED_NAMES), f"deleted names still defined: {sorted(defined & _DELETED_NAMES)}"
    tests = {n.name for n in tree.body if isinstance(n, ast.FunctionDef) and n.name.startswith("test_")}
    assert tests == _KEPT_TESTS, sorted(tests ^ _KEPT_TESTS)
    docstring = ast.get_docstring(tree) or ""
    assert "shared engine" in docstring and "tests/gates/test_retirements.py" in docstring, docstring[:200]
    source = _WAVE2_TEST.read_text(encoding="utf-8")
    assert "test_the_census_reads_nothing_in_a_clean_source" not in source.split('"""', 2)[2], "a rests_on line still names a deleted test"
    assert "rows_wave2" in source


# ---------------------------------------------------------------------------
# Done-when 7: no second definition of import resolution
# ---------------------------------------------------------------------------


def _alpha_normal_form(fn: ast.AST) -> str:
    """The ``ast.dump`` of a function with its name, docstring and every local (argument or ``Name``) renamed by first occurrence."""
    node = ast.parse(ast.unparse(fn)).body[0]
    if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
        raise TypeError(node)
    node.name = "_"
    if node.body and isinstance(node.body[0], ast.Expr) and isinstance(node.body[0].value, ast.Constant) and isinstance(node.body[0].value.value, str):
        node.body = node.body[1:] or [ast.Pass()]
    node.returns = None
    seen: dict[str, str] = {}
    for sub in ast.walk(node):
        if isinstance(sub, ast.arg):
            sub.annotation = None
            sub.arg = seen.setdefault(sub.arg, f"v{len(seen)}")
        elif isinstance(sub, ast.Name):
            sub.id = seen.setdefault(sub.id, f"v{len(seen)}")
        elif isinstance(sub, ast.AnnAssign):
            pass
    return ast.dump(node)


def _import_resolution_census(source: str, scan_functions: dict[str, ast.AST]) -> dict[str, object]:
    """Of an engine source: the calls into ``imported_names`` / ``forbidden_imports`` and the functions α-equivalent to either of the two
    functions of ``tests.gates._import_scan`` (``scan_functions``, name to node)."""
    tree = ast.parse(source)
    calls = [
        n for n in ast.walk(tree)
        if isinstance(n, ast.Call) and (
            (isinstance(n.func, ast.Name) and n.func.id in scan_functions)
            or (isinstance(n.func, ast.Attribute) and n.func.attr in scan_functions)
        )
    ]
    normal = {name: _alpha_normal_form(fn) for name, fn in scan_functions.items()}
    twins = [
        n.name for n in ast.walk(tree)
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and n.name not in scan_functions and _alpha_normal_form(n) in normal.values()
    ]
    return {"calls": len(calls), "twins": twins}


def _scan_functions() -> dict[str, ast.AST]:
    from tests.gates import _import_scan

    tree = ast.parse(inspect.getsource(_import_scan))
    return {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in ("imported_names", "forbidden_imports")}


@pytest.mark.foundation
def test_the_census_of_import_resolution_finds_a_planted_copy_and_a_planted_call_and_reads_a_clean_source() -> None:
    r"""The positive control of the census below (X1): a source that copies ``imported_names`` under another name with renamed
    locals has one twin and no call; a source that calls it has one call and no twin; a source that does neither has neither.
    Red input: a census that compares names, or that never counts a call."""
    fns = _scan_functions()
    assert set(fns) == {"imported_names", "forbidden_imports"}
    renamed = ast.parse(ast.unparse(fns["imported_names"])).body[0]
    assert isinstance(renamed, ast.FunctionDef)
    renamed.name = "resolve_everything"
    for sub in ast.walk(renamed):
        if isinstance(sub, ast.Name):
            sub.id += "_x"
        elif isinstance(sub, ast.arg):
            sub.arg += "_x"
    copy = ast.unparse(renamed)
    assert _import_resolution_census(copy, fns) == {"calls": 0, "twins": ["resolve_everything"]}
    call = "from tests.gates._import_scan import imported_names\n\ndef g(src):\n    return imported_names(src, 'p')\n"
    assert _import_resolution_census(call, fns) == {"calls": 1, "twins": []}
    assert _import_resolution_census("def h(x):\n    return x + 1\n", fns) == {"calls": 0, "twins": []}


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_census_of_import_resolution_finds_a_planted_copy_and_a_planted_call_and_reads_a_clean_source")
def test_the_engine_calls_the_one_import_resolver_and_defines_no_second(engine: SimpleNamespace) -> None:
    r"""AST census of the engine's non-rows files (``__init__``, ``_finder``, ``_walk``): at least one call into
    ``tests.gates._import_scan`` (``imported_names`` or ``forbidden_imports``) and no function α-equivalent (locals renamed) to either
    of the two functions of that module. Predicate stated: this bounds a COPY or a rename of the resolver; it does not bound a
    different spelling of relative-import arithmetic, which the controls above read by behaviour. Red input: an engine that carries
    its own resolution (the review engines read ``ImportFrom.level`` themselves)."""
    fns = _scan_functions()
    files = [p for p in sorted(engine.package_dir.glob("*.py")) if not p.name.startswith("rows_")]
    assert {"__init__.py", "_finder.py", "_walk.py"} <= {p.name for p in files}, [p.name for p in files]
    calls = 0
    for path in files:
        result = _import_resolution_census(path.read_text(encoding="utf-8"), fns)
        assert result["twins"] == [], f"{path.name} defines a second import resolution: {result['twins']}"
        calls += int(result["calls"])  # type: ignore[call-overload]
    assert calls >= 1, "the engine never calls tests.gates._import_scan"
