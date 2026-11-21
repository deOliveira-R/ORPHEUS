r"""The gates of the shared retirement engine (specification ``ontology_w0_spec.md`` section 1.5, unit ``roles``, first slice).

The engine is ``tests/gates/_retired/`` (the five row types in ``__init__``, the definition finder ``_finder``, the walk
``_walk``) and its data are the ``rows_*.py`` modules of that package. This file is its done-when: each test below names the input
that reddens it (X1) and each is RED at the base tip because the engine does not exist there (the fixture :func:`engine` fails with
``ModuleNotFoundError`` on every node; a module-level import would read as one collection error and ``FAILED = 0``).

**The seam this file pins (choices of the gate author, section 1.5 not fixing the spelling).** A builder conforms to these names:

* ``tests.gates._retired``: the NamedTuples ``Retired``, ``Moved``, ``RetiredName``, ``RetiredModule``, ``Kept`` (fields as 1.5).
* ``tests.gates._retired._finder.definitions(source: str) -> dict[str, int]``: qualname to line of every definition of ``source``.
* ``tests.gates._retired._walk.load_rows(rows_dir: Path | None = None) -> list``: the ``ROWS`` then ``KEPT`` of every
  ``rows_*.py`` of ``rows_dir`` (default: the package directory, resolved from ``_walk.__file__``).
* ``tests.gates._retired._walk.walk(root: Path, rows: Iterable | None = None, rows_dir: Path | None = None) -> dict``: the per-row
  verdicts. ``rows`` (rows and ``Kept`` entries mixed) overrides the loader. The key of a verdict is
  ``(type(row).__name__, row)`` (decision C12: ``RetiredName('x') == RetiredModule('x')`` as tuples); the value is a list of
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
        R=_retired, finder=_finder, walk=_walk.walk, load_rows=_walk.load_rows, GateError=_walk.GateError,
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
    """Statements that must read nothing (1.5): a homonym ancestor, a different object, an unbound receiver, ``hasattr``, a fixture
    string of another module."""
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
        hasattr(m, "NAME")
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
    rows in one list; the walk returns two keys, ``(type name, row)`` each. Red input: an engine keyed by the bare tuple returns
    one key. Positive control: a ``Retired`` and a ``Kept`` of equal fields are also two keys."""
    R = engine.R
    verdicts = engine.walk(tmp_path, rows=[R.RetiredName("x"), R.RetiredModule("x")])
    assert set(verdicts) == {("RetiredName", R.RetiredName("x")), ("RetiredModule", R.RetiredModule("x"))}
    both = engine.walk(tmp_path, rows=[R.Retired("a.py", "b"), R.Kept("a.py", "b", "k")])
    assert {k[0] for k in both} == {"Retired", "Kept"} and len(both) == 2


# ---------------------------------------------------------------------------
# Controls: one finding per planted line, positive and negative on separate roots
# ---------------------------------------------------------------------------


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
    _assert_exactly(_findings(verdicts, (type(row).__name__, row)), plant.planted)


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_finder_reads_every_definition_shape_it_must")
@pytest.mark.parametrize("moved", [False, True], ids=["retired", "moved"])
def test_the_negative_control_of_a_definition_row_reads_nothing(engine: SimpleNamespace, tmp_path: Path, moved: bool) -> None:
    r"""On its own root (the negative needs ``p/__init__.py`` to bind ``NAME`` to ``.other``, the positive to ``.m``): a module-level
    ``NAME = 1`` and ``f(NAME=True)`` in a file importing nothing from ``p.m``, ``np.polynomial.laguerre`` against a row
    ``(m.py, laguerre)``, ``getattr(x, "NAME")`` with ``x`` unbound, ``hasattr(m, "NAME")``, ``from q import NAME`` and its fixture
    string with ``q`` not an ancestor, and an ancestor that re-exports a DIFFERENT object under the row's name (the homonym). Red
    input: an engine that matches an import by spelling alone, reads ``hasattr``, or counts every ancestor import."""
    plant = Plant(tmp_path)
    _negative_retired(plant)
    rows: list = [engine.R.Retired(_M, "NAME"), engine.R.Retired(_M, "laguerre")]
    if moved:
        plant.add("orpheus/p/home.py", "NAME = 1")
        rows = [engine.R.Moved(_M, "NAME", "orpheus/p/home.py", "NAME"), engine.R.Retired(_M, "laguerre")]
    verdicts = engine.walk(tmp_path, rows=rows)
    for row in rows:
        _assert_exactly(_findings(verdicts, (type(row).__name__, row)), set())


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
    _assert_exactly(_findings(engine.walk(tmp_path, rows=[row]), ("Retired", row)), plant.planted)


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
    _assert_exactly(_findings(engine.walk(tmp_path / "reach", rows=[row]), ("Retired", row)), reach.planted)
    other = Plant(tmp_path / "other")
    other.add("orpheus/p/sub/m.py", '"""Defines nothing."""')
    other.add("orpheus/p/other.py", "NAME = 1")
    other.add("orpheus/p/__init__.py", "from orpheus.p.other import NAME")
    other.add("tests/user.py", "from orpheus.p import NAME")
    _assert_exactly(_findings(engine.walk(tmp_path / "other", rows=[row]), ("Retired", row)), set())


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_finder_reads_every_definition_shape_it_must")
def test_the_string_leg_reads_a_relative_fixture_import_by_its_suffix_for_a_definition_row_only(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""Decision C22 and the string leg of 1.5: the fixture string ``"from ..m import NAME"`` (a relative import whose module is a
    dotted suffix of the row's module at component boundaries) is found for a ``Retired`` row, ``"from ..xm import NAME"`` (a
    suffix NOT at a component boundary) is not, and a ``RetiredName`` row never reads a fixture string (another gate's planted
    control must be able to spell a retired name). Red input: an engine with no suffix rule, a character suffix, or a
    ``RetiredName`` walk over strings."""
    plant = Plant(tmp_path)
    plant.add(_M, '"""Defines nothing."""')
    plant.add("tests/s1.py", 'S = "from ..m import NAME\\n"  #*')
    plant.add("tests/s2.py", 'S = "from ..xm import NAME\\n"')
    rows = [engine.R.Retired(_M, "NAME"), engine.R.RetiredName("NAME")]
    verdicts = engine.walk(tmp_path, rows=rows)
    _assert_exactly(_findings(verdicts, ("Retired", rows[0])), plant.planted)
    _assert_exactly(_findings(verdicts, ("RetiredName", rows[1])), set())


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
    _assert_exactly(_findings(verdicts, ("Retired", nested)), plant.planted)
    _assert_exactly(_findings(verdicts, ("Retired", flat)), set())
    _assert_exactly(_findings(verdicts, ("Retired", absent)), set())


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
        _assert_exactly(_findings(engine.walk(root, rows=[row]), ("Moved", row)), expected)


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_finder_reads_every_definition_shape_it_must")
def test_the_positive_control_of_a_name_row_finds_each_of_seven_planted_lines_and_its_negative_none(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""``RetiredName('NAME')``: ``getattr(x, "NAME")``, ``x.NAME``, ``__all__ = ["NAME"]``, ``monkeypatch.setattr(m, "NAME", f)``,
    ``from m import NAME``, ``NAME = 1`` and ``f(NAME=True)`` are each found (7 of 7, a ``Store`` and a keyword included). Its
    negative on its own root (``hasattr(m, "NAME")``, the fixture string ``"from m import NAME"``, a bare ``"NAME"`` in a tuple)
    reads nothing. Red input: an engine that reads strings, reads ``hasattr`` (every absence witness of a retirement would
    redden), or misses a store, a keyword or an ``__all__`` member."""
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
    _assert_exactly(_findings(engine.walk(tmp_path / "pos", rows=[row]), ("RetiredName", row)), positive.planted)
    negative = Plant(tmp_path / "neg")
    negative.add("tests/rn.py", """
        hasattr(m, "NAME")
        S = "from m import NAME\\n"
        T = ("NAME", "other")
        """)
    _assert_exactly(_findings(engine.walk(tmp_path / "neg", rows=[row]), ("RetiredName", row)), set())


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
    _assert_exactly(_findings(engine.walk(tmp_path / "pos", rows=[row]), ("RetiredModule", row)), positive.planted)
    negative = Plant(tmp_path / "neg")
    negative.add("tests/n.py", """
        import orpheus.p.m2
        import orpheus.p.mx
        T = ("orpheus.p.m", "other")
        """)
    _assert_exactly(_findings(engine.walk(tmp_path / "neg", rows=[row]), ("RetiredModule", row)), set())


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_finder_reads_every_definition_shape_it_must")
@pytest.mark.parametrize("planted_rel", ["orpheus/p/m.py", "orpheus/p/m/__init__.py"])
def test_the_absence_leg_reports_a_leftover_module_file_or_package_and_is_silent_on_a_root_with_neither(
    engine: SimpleNamespace, tmp_path: Path, planted_rel: str,
) -> None:
    r"""Decision E-D6: ``RetiredModule('orpheus.p.m')`` is red when ``orpheus/p/m.py`` or ``orpheus/p/m/__init__.py`` exists under
    the root, though nothing imports it (1 of 1 each, one root each); a root holding neither is silent. Red input: an engine
    with the import legs only (the wave-2 shape of the check lived in a separate path test that the port deletes)."""
    row = engine.R.RetiredModule("orpheus.p.m")
    planted = Plant(tmp_path / "planted")
    planted.add(planted_rel, '"""A leftover."""  #*')
    found = _findings(engine.walk(tmp_path / "planted", rows=[row]), ("RetiredModule", row))
    assert len(found) == 1 and found[0][0] == planted_rel, found
    clean = Plant(tmp_path / "clean")
    clean.add("orpheus/p/other.py", "x = 1")
    _assert_exactly(_findings(engine.walk(tmp_path / "clean", rows=[row]), ("RetiredModule", row)), set())


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
        _assert_exactly(_findings(engine.walk(tmp_path / label, rows=[entry]), ("Kept", entry)), expected)


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
    _assert_exactly(_findings(engine.walk(tmp_path, rows=[row]), ("Retired", row)), plant.planted)


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
    assert red == {("Retired", row)}, sorted(map(repr, red))
    assert (row.file, 2) in {tuple(f) for f in after[("Retired", row)]}


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_positive_control_of_a_name_row_finds_each_of_seven_planted_lines_and_its_negative_none")
def test_a_use_planted_in_the_copy_turns_exactly_its_name_row_red(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""Mutation of a ``RetiredName`` row: a test file of the copy assigns the row's name (built from the row's field). Red set
    exactly that row; the planted file is among its findings."""
    row = _first(engine, engine.R.RetiredName)

    def mutate(root: Path) -> None:
        (root / "tests" / "planted_use.py").write_text(f"{row.name} = 1\n", encoding="utf-8")

    red, after = _mutate_and_read(engine, tmp_path, mutate)
    assert red == {("RetiredName", row)}, sorted(map(repr, red))
    assert ("tests/planted_use.py", 1) in {tuple(f) for f in after[("RetiredName", row)]}


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_positive_control_of_a_module_row_finds_each_of_seven_planted_lines_and_its_negative_none")
def test_an_import_planted_in_the_copy_turns_exactly_its_module_row_red(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""Mutation of a ``RetiredModule`` row: a test file of the copy imports the row's module; and, separately on a fresh copy, a
    file at the module's path appears (the absence leg). Each turns exactly that row red."""
    row = _first(engine, engine.R.RetiredModule)

    def by_import(root: Path) -> None:
        (root / "tests" / "planted_use.py").write_text(f"import {row.module}\n", encoding="utf-8")

    red, after = _mutate_and_read(engine, tmp_path / "import", by_import)
    assert red == {("RetiredModule", row)}, sorted(map(repr, red))
    assert ("tests/planted_use.py", 1) in {tuple(f) for f in after[("RetiredModule", row)]}

    def by_file(root: Path) -> None:
        target = root / (row.module.replace(".", "/") + ".py")
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("x = 1\n", encoding="utf-8")

    red, _ = _mutate_and_read(engine, tmp_path / "file", by_file)
    assert red == {("RetiredModule", row)}, sorted(map(repr, red))


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_home_leg_is_green_for_a_home_that_defines_the_name_and_red_for_one_that_lacks_it_or_is_absent")
def test_a_moved_row_goes_red_when_its_old_definition_returns_and_when_its_home_loses_the_name(engine: SimpleNamespace, tmp_path: Path) -> None:
    r"""Mutation of a ``Moved`` row (the root holds its ``file`` AND its ``home``, decision E-D2): (i) the definition returns to
    the old file; (ii) on a fresh copy, the home loses the name. Each turns exactly the ``Moved`` row red."""
    row = next(r for r in _mutation_rows(engine) if type(r) is engine.R.Moved)

    def returns(root: Path) -> None:
        (root / row.file).write_text(f"def {row.qualname}(): pass\n", encoding="utf-8")

    red, _ = _mutate_and_read(engine, tmp_path / "returns", returns)
    assert red == {("Moved", row)}, sorted(map(repr, red))

    def loses(root: Path) -> None:
        (root / row.home).write_text("y = 1\n", encoding="utf-8")

    red, _ = _mutate_and_read(engine, tmp_path / "loses", loses)
    assert red == {("Moved", row)}, sorted(map(repr, red))


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
    assert red == {("Kept", entries[0])}, sorted(map(repr, red))


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
    assert ("RetiredName", engine.R.RetiredName("only_in_the_root")) not in verdicts


# ---------------------------------------------------------------------------
# Done-when 2: the real tree, wave 2
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def real_walk(engine: SimpleNamespace) -> dict:
    """One walk of the real tree (about 40 s) with the default rows; every row-reading test below shares it."""
    return engine.walk(_REPO)


@pytest.mark.foundation
@pytest.mark.rests_on(
    f"{_THIS}::test_the_positive_control_of_a_definition_row_finds_each_of_nine_planted_lines",
    f"{_THIS}::test_the_positive_control_of_a_name_row_finds_each_of_seven_planted_lines_and_its_negative_none",
    f"{_THIS}::test_the_positive_control_of_a_module_row_finds_each_of_seven_planted_lines_and_its_negative_none",
    f"{_THIS}::test_the_absence_leg_reports_a_leftover_module_file_or_package_and_is_silent_on_a_root_with_neither",
)
def test_every_wave2_row_is_green_on_the_tree(engine: SimpleNamespace, real_walk: dict) -> None:
    r"""The 65 wave-2 rows (2 ``Retired`` for the two retired properties, 60 ``RetiredName``, 3 ``RetiredModule``; ``KEPT = ()``, decisions
    D5 and R2) read green on the real tree: ``k of 65`` is printed, and the three counts are the specification's. RECORD of the tree.
    The three module rows are the three paths the deleted path test named (decision E-D6): each module's file and package directory
    are absent. Red input: a wave-2 name, module or definition still present anywhere in the five trees."""
    R = engine.R
    rows = list(engine.rows_wave2.ROWS)
    by_type = {name: sum(1 for r in rows if type(r).__name__ == name) for name in _ROW_TYPE_NAMES}
    assert by_type == {"Retired": 2, "Moved": 0, "RetiredName": 60, "RetiredModule": 3}, by_type
    assert tuple(engine.rows_wave2.KEPT) == ()
    red = {key: real_walk[key] for key in ((type(r).__name__, r) for r in rows) if real_walk[key]}
    print(f"wave-2 rows green: {len(rows) - len(red)} of {len(rows)}")
    assert not red, f"{len(rows) - len(red)} of {len(rows)} wave-2 rows green; red: " + "; ".join(f"{k[1]}: {v[:3]}" for k, v in sorted(red.items(), key=repr))
    for row in rows:
        if isinstance(row, R.RetiredModule):
            base = _REPO / row.module.replace(".", "/")
            assert not base.with_suffix(".py").exists() and not (base / "__init__.py").exists(), row


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_every_wave2_row_is_green_on_the_tree")
def test_the_real_walk_reads_every_row_of_the_package(engine: SimpleNamespace, real_walk: dict) -> None:
    """The input count beside the finding count (X1): the real walk holds one verdict per entry of ``load_rows()``, non-zero."""
    entries = engine.load_rows()
    assert entries and len(real_walk) == len(entries), (len(real_walk), len(entries))


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
@pytest.mark.rests_on(f"{_THIS}::test_every_wave2_row_is_green_on_the_tree")
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
