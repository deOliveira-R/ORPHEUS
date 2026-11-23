"""Machinery gates of the concept table: planted fixtures over synthetic inputs.

The concept table (``docs/architecture/concept_table.toml``) says which concept each class of ``orpheus/``
realises, and ``tools/concept_table/`` reads a declaration (a trailing comment on the line that ends a class
header), the table, the debt and the references between classes. This file tests that machinery on synthetic
sources and synthetic tables, so each check is shown to read both a clean input and a violating one. It reads
neither the table nor the declarations of the tree; the gates that do are in ``test_concept_table.py``.

Every check is pure and returns the frozenset of its violations, typed as the debt entry that would pin it, so a
planted fixture is a value and never a file edited on disk.

The names this file takes from ``tools.concept_table`` and is the contract it holds the builder to (the
specification fixes those marked S; the rest are fixed here):

* ``declaration``: ``ROLES``, ``KINDS``, ``SUBKINDS``, ``Declaration``, ``DeclarationError``,
  ``parse_declaration``, ``header_end``, ``read_declarations`` (S, the reviewed kit file).
* ``population``: ``ClassKey`` (S, with ``dotted``), ``Population.from_sources`` (S), ``read_tree`` (S),
  ``DuplicateClass`` (S); a population holds ``rows`` (each with ``key``, ``declaration``, ``markers``) and
  ``orphan_markers``, ``files_read``, ``files_skipped_derivations`` (S).
* ``refs``: ``read_edges(sources, population)`` giving ``Edge(src, dst, kind)`` values, ``alias_leaves`` (S),
  ``is_type_checking_guard``.
* ``table``: ``loads_table(text)`` and ``load_table(path)`` giving concepts, ``TableError`` (a ``ValueError``).
* ``debt``: ``loads_debt(text)``, ``load_debt(path)``, ``DebtError`` (a ``ValueError``), ``pin`` (S), ``wave_of`` (S),
  and the entry types ``ClassDebt``, ``RoleTestDebt``, ``EdgeDebt``, ``CycleDebt``, ``KindLawDebt``,
  ``OffAlgebraDebt``, ``CovectorOffLineageDebt`` (S).
* ``checks``: ``check_table``, ``check_declarations``, ``check_heads``, ``check_roles``, ``check_kinds``,
  ``check_off_algebra``, ``check_covector_lineage``, ``check_phase``, ``check_layer``, ``check_placement``,
  ``check_cycles`` (S), ``may_reference`` (S), ``fingerprint(source, qualname)``, ``check_unweld(population,
  seed, sources)`` returning the pair (uncovered, relabelled), ``check_wave_rule(debt, concepts, population)``.
"""

from __future__ import annotations

import ast
import dataclasses
import functools
import importlib
import types
from pathlib import Path

import pytest

pytestmark = pytest.mark.foundation

ROOT = Path(__file__).resolve().parents[2]
_MODULES = ("declaration", "population", "refs", "table", "debt", "checks")


@functools.lru_cache(maxsize=None)
def ct() -> types.SimpleNamespace:
    """The modules of ``tools.concept_table``, imported by the test that needs them, so each test fails alone."""
    return types.SimpleNamespace(**{m: importlib.import_module(f"tools.concept_table.{m}") for m in _MODULES})


# --- synthetic inputs ----------------------------------------------------------------------------------------

FORBIDDEN = {
    "values": frozenset(),
    "numerics": frozenset({"sn", "transport"}),
    "geometry": frozenset({"numerics", "sn", "transport"}),
    "transport": frozenset({"sn"}),
    "sn": frozenset(),
}


def toml_row(cid: str, kind: str = "element", phase: str | None = "structure", layer: str = "numerics", *,
             wave: str | None = "W1", ruling: str = "O23", extra: str = "") -> str:
    lines = ["[[concept]]", f'id = "{cid}"', f'name = "{cid}"', f'definition = "the {cid}"', f'kind = "{kind}"']
    if phase is not None:
        lines.append(f'phase = "{phase}"')
    lines.append(f'layer = "{layer}"')
    if wave is not None:
        lines.append(f'wave = "{wave}"')
    lines.append(f'ruling = "{ruling}"')
    return "\n".join(lines) + ("\n" + extra if extra else "") + "\n"


def decl(cid: str, role: str, kind: str, welded: str = "") -> str:
    return f"# concept: {cid} | role: {role} | kind: {kind}" + (f" | welded: {welded}" if welded else "")


def klass(name: str, concept: str, role: str, kind: str = "element", bases: str = "", welded: str = "",
          body: str = "    pass") -> str:
    return f"class {name}({bases}):  {decl(concept, role, kind, welded)}\n{body}\n\n"


def key(file: str, qualname: str):
    return ct().population.ClassKey(file, qualname)


# --- the declaration grammar and the reader (2.3) ------------------------------------------------------------

def test_reader_positive_control_and_single_line_header():
    """A planted one-class source with a marker reads exactly one marker; the same class without it reads zero."""
    d = ct().declaration
    src = f"class A:  {decl('x', 'arm', 'set')}\n    pass\n"
    found, orphans = d.read_declarations(src)
    assert len(found["A"][2]) == 1 and orphans == []
    assert d.read_declarations("class A:\n    pass\n")[0]["A"][2] == []


def test_reader_decorated_class_marker_on_decorator_line_is_an_orphan():
    d = ct().declaration
    src = f"@dataclass  {decl('x', 'arm', 'set')}\nclass A:\n    pass\n"
    found, orphans = d.read_declarations(src)
    assert found["A"][2] == [] and orphans == [1]


def test_reader_multiline_header_marks_only_the_last_line():
    d = ct().declaration
    good = f"class A(\n    B,\n    C,\n):  {decl('x', 'arm', 'set')}\n    pass\n"
    early = f"class A(  {decl('x', 'arm', 'set')}\n    B,\n):\n    pass\n"
    assert d.read_declarations(good)[0]["A"] == (1, 4, [f"{decl('x', 'arm', 'set')}"])
    found, orphans = d.read_declarations(early)
    assert found["A"][2] == [] and orphans == [1]


def test_reader_keeps_a_trailing_comment_and_the_parser_stops_at_the_second_hash():
    d = ct().declaration
    src = f"class A(B):  {decl('x', 'arm', 'set')}  # noqa: E501\n    pass\n"
    marker = d.read_declarations(src)[0]["A"][2][0]
    assert d.parse_declaration(marker) == d.Declaration("x", "arm", "set", None, ())


def test_reader_function_local_and_nested_classes():
    d = ct().declaration
    src = ("def f():\n    class L:  " + decl("x", "part", "record") + "\n        pass\n"
           "class Outer:  " + decl("x", "head", "record") + "\n    class Inner:  " + decl("x", "part", "record") + "\n        pass\n")
    found, _ = d.read_declarations(src)
    assert set(found) == {"f.<locals>.L", "Outer", "Outer.Inner"} and all(len(v[2]) == 1 for v in found.values())


def test_reader_marker_above_the_class_is_an_orphan():
    d = ct().declaration
    found, orphans = d.read_declarations(f"{decl('x', 'arm', 'set')}\nclass A:\n    pass\n")
    assert found["A"][2] == [] and orphans == [1]


def test_declaration_render_and_parse_are_inverse_and_noncanonical_text_is_not_canonical():
    d = ct().declaration
    for dec in (d.Declaration("x", "arm", "set", None, ()), d.Declaration("EXT:citation", "head", "element", "member", ()),
                d.Declaration("a", "weld", "morphism", "operator", ("a", "b")), d.Declaration("a", "head", "set", None, ("a", "b", "c"))):
        assert d.parse_declaration(dec.render()) == dec
    sloppy = "#concept: x|role: arm|kind: set"
    assert d.parse_declaration(sloppy).render() == "# concept: x | role: arm | kind: set" != sloppy


@pytest.mark.parametrize("args", [
    ("x", "boss", "set", None, ()), ("x", "arm", "thing", None, ()), ("x", "arm", "set", "operator", ()),
    ("x", "arm", "morphism", "member", ()), ("X1", "arm", "set", None, ()), ("x", "weld", "set", None, ()),
    ("x", "arm", "set", None, ("x", "y")), ("x", "weld", "set", None, ("y", "x")), ("x", "weld", "set", None, ("x",)),
    ("x", "weld", "set", None, ("x", "y", "y")),
])
def test_declaration_invariants_are_refused(args):
    d = ct().declaration
    with pytest.raises(d.DeclarationError):
        d.Declaration(*args)


@pytest.mark.parametrize("text", [
    "# concept: x | kind: set | role: arm", "# concept: x | role: arm", "# concept: x | role: arm | kind: set | note: y",
    "# concept: x | role: arm | role: arm | kind: set",
])
def test_parse_refuses_wrong_fields_order_and_repeats(text):
    d = ct().declaration
    with pytest.raises(d.DeclarationError):
        d.parse_declaration(text)


# --- the population (2.2, 2.7) ---------------------------------------------------------------------------------

def test_classkey_dotted_is_one_rule():
    k = ct().population.ClassKey
    assert k("orpheus/sn/loss_representation/__init__.py", "LossRepresentation").dotted == "orpheus.sn.loss_representation.LossRepresentation"
    assert k("orpheus/memo/_traced_memo_boot.py", "_libproc_thread_counter.<locals>.TaskInfo").dotted == "orpheus.memo._traced_memo_boot.TaskInfo"
    assert str(k("orpheus/a.py", "A")) == "orpheus/a.py::A"


def test_population_is_sorted_and_refuses_duplicates():
    p = ct().population
    pop = p.Population.from_sources({"orpheus/b.py": klass("Z", "x", "arm") + klass("A", "x", "arm"), "orpheus/a.py": klass("M", "x", "arm")})
    assert [str(r.key) for r in pop.rows] == ["orpheus/a.py::M", "orpheus/b.py::A", "orpheus/b.py::Z"]
    with pytest.raises(p.DuplicateClass):
        p.Population.from_sources({"orpheus/a/__init__.py": "class X:\n    pass\n", "orpheus/a.py": "class X:\n    pass\n"})
    with pytest.raises(Exception):
        dataclasses.replace(pop, rows=tuple(reversed(pop.rows)))


def test_population_records_orphans_and_markers():
    p = ct().population
    pop = p.Population.from_sources({"orpheus/a.py": f"@d  {decl('x', 'arm', 'set')}\nclass A:\n    pass\nclass B:  {decl('x', 'arm', 'set')}\n    pass\n"})
    assert pop.orphan_markers == (("orpheus/a.py", 1),)
    by = {r.key.qualname: r for r in pop.rows}
    assert by["A"].declaration is None and by["B"].declaration is not None and len(by["B"].markers) == 1


def test_read_tree_skips_derivations_and_sorts(tmp_path):
    """The positive control of the skip: the planted derivations file is counted as skipped and absent from the rows."""
    for rel, src in {"orpheus/b.py": "class B:\n    pass\n", "orpheus/a.py": "class A:\n    pass\n",
                     "orpheus/derivations/d.py": "class D:\n    pass\n", "orpheus/__pycache__/x.py": "class P:\n    pass\n"}.items():
        (tmp_path / rel).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / rel).write_text(src)
    pop = ct().population.read_tree(tmp_path)
    assert [r.key.qualname for r in pop.rows] == ["A", "B"]
    assert pop.files_skipped_derivations == 1 and pop.files_read == 2


def test_read_tree_on_the_tree_finds_the_known_classes_and_none_under_derivations():
    pop = ct().population.read_tree(ROOT)
    keys = {str(r.key) for r in pop.rows}
    assert pop.files_skipped_derivations > 0
    assert "orpheus/sn/loss_representation/__init__.py::LossRepresentation" in keys
    assert "orpheus/memo/_traced_memo_boot.py::_libproc_thread_counter.<locals>.TaskInfo" in keys
    assert not any("/derivations/" in k for k in keys)


# --- the reference pass (2.6) ----------------------------------------------------------------------------------

def _edges(sources):
    m = ct()
    pop = m.population.Population.from_sources(sources)
    return {(str(e.src), str(e.dst), e.kind) for e in m.refs.read_edges(sources, pop)}


def test_refs_control_1_inheritance_across_packages_on_the_tree():
    m = ct()
    pop = m.population.read_tree(ROOT)
    sources = {str(r.key.file): (ROOT / r.key.file).read_text() for r in pop.rows}
    got = {(str(e.src), str(e.dst), e.kind) for e in m.refs.read_edges(sources, pop)}
    assert any(s.endswith("::FullFieldSpace") and d.endswith("::FunctionSpace") and k == "inherits" for s, d, k in got)
    assert any(s.endswith("::SNProblem") and d.endswith("::FunctionSpace") for s, d, k in got)


def test_refs_control_2_generic_subscript_is_not_a_base():
    got = _edges({"orpheus/a.py": "class Base:\n    pass\nclass Arg:\n    pass\nclass A(Base[Arg]):\n    pass\n"})
    assert ("orpheus/a.py::A", "orpheus/a.py::Base", "inherits") in got
    assert ("orpheus/a.py::A", "orpheus/a.py::Arg", "top") in got
    assert ("orpheus/a.py::A", "orpheus/a.py::Arg", "inherits") not in got


def test_refs_control_3_all_list_is_not_a_reference():
    """Neither a class edge nor a function-owned edge: an ``__all__`` string belongs to the module, not to a class, so the
    class half alone cannot see the skip (removing it read green there and red only on the tree's role gate)."""
    m = ct()
    sources = {"orpheus/a.py": '__all__ = ["B"]\nclass B:\n    pass\n'}
    got = _edges(sources)
    assert not any(d.endswith("::B") for _, d, _ in got)
    pop = m.population.Population.from_sources(sources)
    assert not [e for e in m.refs.read_function_edges(sources, pop) if e.dst.qualname == "B"]
    positive = {"orpheus/a.py": "class B:\n    pass\nx = B\n"}
    assert [e.dst.qualname for e in m.refs.read_function_edges(positive, m.population.Population.from_sources(positive))] == ["B"], \
        "positive control: a function-owned reference to B is read"


def test_refs_nested_class_is_an_edge_end_and_an_owner():
    """``Outer.Inner`` is a class of the population: a reference to it from another file is an edge to it, and its own body's
    references are its own, beside those the class that holds it owns."""
    sources = {"orpheus/a.py": "class Base:\n    pass\nclass Other:\n    pass\nclass Outer:\n    class Inner(Base):\n        x: Other\n",
               "orpheus/b.py": "from orpheus.a import Outer\nclass User:\n    x: Outer.Inner\n"}
    got = _edges(sources)
    assert ("orpheus/b.py::User", "orpheus/a.py::Outer.Inner", "top") in got
    assert ("orpheus/a.py::Outer.Inner", "orpheus/a.py::Base", "inherits") in got
    assert ("orpheus/a.py::Outer.Inner", "orpheus/a.py::Other", "top") in got
    assert ("orpheus/a.py::Outer", "orpheus/a.py::Other", "top") in got


def test_class_node_is_looked_up_by_the_qualname_walk():
    c = ct().checks
    src = "class A:\n    def m(self):\n        class B:\n            pass\n"
    assert c._class_node(src, "A").lineno == 1 and c._class_node(src, "A.m.<locals>.B").lineno == 3
    with pytest.raises(ValueError):
        c._class_node(src, "A.B")  # a class B inside a method is not a class B of A


def test_refs_control_4_union_alias_expands_to_its_leaves():
    src = ("class Modal:\n    pass\nclass FromFrame:\n    pass\nAxisKind = Modal | FromFrame\n"
           "class Axis:\n    kind: AxisKind\n")
    got = _edges({"orpheus/a.py": src})
    assert {("orpheus/a.py::Axis", "orpheus/a.py::Modal", "top"), ("orpheus/a.py::Axis", "orpheus/a.py::FromFrame", "top")} <= got
    m = ct()
    leaves = m.refs.alias_leaves({"orpheus/a.py": src}, "orpheus/a.py", "AxisKind")
    assert {l.qualname for l in leaves} == {"Modal", "FromFrame"}


def test_refs_kinds_and_the_strongest_kind_of_a_pair():
    src = ("from typing import TYPE_CHECKING\nclass B:\n    pass\nclass C:\n    pass\n"
           "class A(B):\n    x: B\n    def f(self):\n        from orpheus.c import C\n")
    got = _edges({"orpheus/a.py": src, "orpheus/c.py": "class C:\n    pass\n"})
    assert ("orpheus/a.py::A", "orpheus/a.py::B", "inherits") in got and ("orpheus/a.py::A", "orpheus/a.py::B", "top") not in got
    assert ("orpheus/a.py::A", "orpheus/c.py::C", "late") in got
    typed = _edges({"orpheus/a.py": "from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n    from orpheus.c import C\nclass A:\n    x: C\n",
                    "orpheus/c.py": "class C:\n    pass\n"})
    assert ("orpheus/a.py::A", "orpheus/c.py::C", "typing") in typed


def test_refs_control_5_not_type_checking_is_a_top_import():
    got = _edges({"orpheus/a.py": "from typing import TYPE_CHECKING\nif not TYPE_CHECKING:\n    from orpheus.c import C\nclass A:\n    x: C\n",
                  "orpheus/c.py": "class C:\n    pass\n"})
    assert ("orpheus/a.py::A", "orpheus/c.py::C", "top") in got and ("orpheus/a.py::A", "orpheus/c.py::C", "typing") not in got


def test_refs_control_6_guard_witness_three_copies_agree():
    """The pass's predicate and the two other walkers' predicates agree on every ``ast.If`` of the tree and on the planted ``not``."""
    from tests.gates.numerics.test_concept_order import _is_type_checking_guard as b
    from tests.gates.test_layer_imports import _is_type_checking_guard as a

    mine = ct().refs.is_type_checking_guard
    planted = ast.parse("if not TYPE_CHECKING:\n    pass").body[0].test  # type: ignore[attr-defined]
    assert not (mine(planted) or a(planted) or b(planted))
    n = guards = 0
    for p in sorted((ROOT / "orpheus").rglob("*.py")):
        if "derivations" in p.parts or "__pycache__" in p.parts:
            continue
        for x in ast.walk(ast.parse(p.read_text())):
            if isinstance(x, ast.If):
                n += 1
                guards += bool(a(x.test))
                assert bool(mine(x.test)) == bool(a(x.test)) == bool(b(x.test)), (p, x.lineno)
    assert n > 0 and guards > 0, "positive control: the walk reads if-statements and finds guards"


# --- the loaders (G0's single-row half) ------------------------------------------------------------------------

def test_loader_reads_a_clean_table_and_derives_heads_forms():
    m = ct()
    text = toml_row("a") + toml_row("m", "apparatus", None, "memo", wave=None) + toml_row("d", "morphism", "structure", "values", wave=None, ruling="O36", extra='functions_head = true\nfunctions = ["orpheus.values.digest.digest"]') + toml_row("g", wave="W2", extra="missing = true")
    rows = m.table.loads_table(text)
    assert [c.id for c in rows] == ["a", "m", "d", "g"] and rows[1].phase is None
    assert isinstance(rows[0].head, m.table.ClassHead) and isinstance(rows[2].head, m.table.FunctionsHead) and isinstance(rows[3].head, m.table.MissingHead)


@pytest.mark.parametrize("text", [
    toml_row("a").replace('ruling = "O23"\n', ""),                                  # a row without a ruling
    toml_row("a", phase="none"),                                                    # the absence is the key's absence
    toml_row("m", "apparatus", "structure"),                                        # an apparatus with a phase
    toml_row("a", phase=None),                                                      # a non-apparatus without a phase
    toml_row("a", extra='status = "done"'),                                         # an unknown key
    toml_row("a", extra='sum = "orpheus/a.py:A"\nmissing = true'),                  # two head forms
    toml_row("a", wave=None, extra="missing = true"),                               # missing without a wave
    toml_row("a", extra='arms = ["orpheus/a.py::A"]'),                              # arms without missing
    toml_row("a", wave=None, ruling="O23", extra='functions_head = true\nfunctions = ["orpheus.a.f"]'),  # ruling other than O36
    toml_row("a", wave=None, ruling="O36", extra="functions_head = true"),          # empty functions
    toml_row("a", wave="W0"),                                                       # W0 is never a row's wave
    toml_row("a", kind="thing"),
    toml_row("A b"),
])
def test_loader_refuses_single_row_violations(text):
    with pytest.raises(ct().table.TableError):
        ct().table.loads_table(text)


def test_debt_loader_refuses_schema_violations():
    m = ct()
    assert m.debt.loads_debt('[[twin]]\nclass = "orpheus/a.py::A"\n').twin == (m.debt.ClassDebt(key("orpheus/a.py", "A")),)
    for bad in ('[[twin]]\nclass = "orpheus/a.py::A"\nwave = "W1"\n', '[[nosuch]]\nclass = "orpheus/a.py::A"\n',
                '[[role_test]]\ntest = "nosuch"\nclass = "orpheus/a.py::A"\n', '[[kind_law]]\nlaw = "annex_record"\nclass = "orpheus/a.py::A"\n',
                '[[edge]]\nsrc = "orpheus/a.py::A"\ndst = "orpheus/b.py::B"\nkind = "sideways"\n'):
        with pytest.raises(m.debt.DebtError):
            m.debt.loads_debt(bad)


def test_pin_returns_unpinned_and_stale():
    pin = ct().debt.pin
    assert pin(frozenset({1, 2}), frozenset({2, 3})) == (frozenset({1}), frozenset({3}))
    assert pin(frozenset({1}), frozenset({1})) == (frozenset(), frozenset())


# --- G0: check_table --------------------------------------------------------------------------------------------

def _table_findings(text: str, *, expected: int, labels=None, resolve=lambda n: True, layers=("numerics", "memo", "values")):
    m = ct()
    return m.checks.check_table(m.table.loads_table(text), frozenset(layers), expected, labels or {}, resolve)


def test_check_table_clean_and_planted_violations():
    clean = toml_row("a") + toml_row("b")
    assert not _table_findings(clean, expected=2)
    assert _table_findings(clean, expected=3), "the count is EXPECTED_CONCEPTS"
    assert _table_findings(toml_row("a") + toml_row("a"), expected=2), "duplicate id"
    assert _table_findings(toml_row("EXT:citation") + toml_row("ext-citation"), expected=2), "one slug for two ids"
    assert _table_findings(toml_row("a", layer="nowhere"), expected=1), "layer not in the tree"
    labelled = toml_row("a", extra='theory = ["ref:real"]') + toml_row("b", extra='theory = ["ref:real"]')
    assert _table_findings(labelled, expected=2, labels={"real": "ref"}), "a label under two concepts"
    assert not _table_findings(toml_row("a", extra='theory = ["ref:real"]') + toml_row("b"), expected=2, labels={"real": "ref"})
    assert _table_findings(toml_row("a", extra='theory = ["ref:nosuch"]'), expected=1, labels={"real": "ref"}), "undefined label"
    assert _table_findings(toml_row("a", extra='theory = ["eq:real"]'), expected=1, labels={"real": "ref"}), "eq: needs an equation label"
    assert _table_findings(toml_row("a", extra='functions = ["orpheus.numerics.space.FunctionSpace.nosuch_member"]'), expected=1, resolve=lambda n: False)


def test_label_scan_and_function_resolver_controls():
    c = ct().checks
    labels = c.doc_labels(ROOT / "docs")
    assert labels.get("architecture-conceptual-view") == "ref" and labels.get("spaces-axis-product") == "eq"
    assert "spaces-the-axiz" not in labels
    resolve = c.function_resolver(ROOT)
    assert resolve("orpheus.numerics.space.FunctionSpace")
    assert not resolve("orpheus.numerics.nosuch.Thing") and not resolve("orpheus.numerics.space.FunctionSpace.nosuch_member")


# --- G1 ----------------------------------------------------------------------------------------------------------

def _decls(sources: dict[str, str], ids=("x", "y")):
    m = ct()
    return m.checks.check_declarations(m.population.Population.from_sources(sources), frozenset(ids))


def test_check_declarations_planted():
    assert not _decls({"orpheus/a.py": klass("A", "x", "arm")})
    assert _decls({"orpheus/a.py": "class A:\n    pass\n"}), "undeclared"
    assert _decls({"orpheus/a.py": f"@d  {decl('x', 'arm', 'set')}\nclass A:\n    pass\n"}), "marker on the decorator line"
    assert _decls({"orpheus/a.py": f"class A(\n    B,  {decl('x', 'arm', 'set')}\n):\n    pass\n"}), "marker on the first header line"
    assert _decls({"orpheus/a.py": klass("A", "nosuch", "arm")}), "unknown concept"
    assert _decls({"orpheus/a.py": klass("A", "x", "weld", "element", welded="x, nosuch")}), "welded id not a row"
    assert not _decls({"orpheus/a.py": klass("A", "x", "weld", "element", welded="x, y")})
    assert _decls({"orpheus/a.py": f"{decl('x', 'arm', 'set')}\nclass A:\n    pass\n"}), "marker above the class"
    assert _decls({"orpheus/a.py": "class A:  #concept: x|role: arm|kind: set\n    pass\n"}), "non-canonical spelling"
    assert _decls({"orpheus/derivations/d.py": klass("D", "x", "arm")}), "a class under derivations"


def test_check_declarations_orphan_marker_arm_alone():
    """A declared class plus a stray marker on a body line: the class is not the finding, only the orphan arm can be."""
    found = _decls({"orpheus/a.py": klass("A", "x", "arm", body=f"    x = 1  {decl('x', 'arm', 'set')}")})
    assert any("on no class header-end line" in f for f in found), found
    assert not any("no declaration" in f for f in found)


# --- G2 ----------------------------------------------------------------------------------------------------------

def _heads(table_text: str, sources: dict[str, str]):
    m = ct()
    return m.checks.check_heads(m.population.Population.from_sources(sources), m.table.loads_table(table_text), sources)


def test_check_heads_planted():
    one = {"orpheus/numerics/a.py": klass("H", "x", "head") + klass("A", "x", "arm")}
    assert not _heads(toml_row("x"), one)
    assert _heads(toml_row("x"), {"orpheus/numerics/a.py": klass("H", "x", "head") + klass("H2", "x", "head")}), "two heads"
    assert _heads(toml_row("x"), {"orpheus/numerics/a.py": klass("A", "x", "arm")}), "a class head concept with no head"
    assert not _heads(toml_row("x", extra="missing = true"), {"orpheus/numerics/a.py": klass("A", "x", "arm")})
    assert _heads(toml_row("x", extra="missing = true"), one), "missing but a head exists"
    fn = toml_row("x", wave=None, ruling="O36", extra='functions_head = true\nfunctions = ["orpheus.numerics.a.f"]')
    assert _heads(fn, {"orpheus/numerics/a.py": klass("H", "x", "head")}), "a functions head concept given a head class"
    sumrow = toml_row("x", extra='sum = "orpheus/numerics/a.py:U"')
    ok = {"orpheus/numerics/a.py": klass("A", "x", "arm") + klass("B", "x", "arm") + "U = A | B\n"}
    assert not _heads(sumrow, ok)
    assert _heads(sumrow, {"orpheus/numerics/a.py": klass("A", "x", "arm") + klass("B", "x", "arm") + "class U:\n    pass\n"}), "alias is a class"
    assert _heads(sumrow, {"orpheus/numerics/a.py": klass("A", "x", "arm") + klass("B", "y", "arm") + "U = A | B\n"}), "a leaf of another concept"
    wrong_role = _heads(sumrow, {"orpheus/numerics/a.py": klass("A", "x", "arm") + klass("B", "x", "part") + "U = A | B\n"})
    assert any("sum leaf" in f for f in wrong_role), "a leaf of the concept that is not declared an arm"
    assert _heads(sumrow, {"orpheus/numerics/a.py": klass("A", "x", "arm") + klass("B", "x", "arm") + klass("H", "x", "head") + "U = A | B\n"}), "head beside a sum"
    miss = toml_row("x", extra='missing = true\narms = ["orpheus/numerics/a.py::A", "orpheus/numerics/a.py::Z"]')
    assert _heads(miss, {"orpheus/numerics/a.py": klass("A", "x", "arm")}), "an arm key that is no class"


# --- G3 ----------------------------------------------------------------------------------------------------------

def _roles(sources: dict[str, str], table_text: str = toml_row("x") + toml_row("y")):
    m = ct()
    pop = m.population.Population.from_sources(sources)
    return m.checks.check_roles(pop, frozenset(m.refs.read_edges(sources, pop)), m.table.loads_table(table_text), sources)


def _role_tests(found) -> set[tuple[str, str]]:
    return {(r.test, r.key.qualname) for r in found}


def test_check_roles_planted():
    head = klass("H", "x", "head")
    assert not _roles({"orpheus/numerics/a.py": head + klass("A", "x", "arm", bases="H")})
    assert ("arm", "A") in _role_tests(_roles({"orpheus/numerics/a.py": head + klass("A", "x", "arm")})), "an arm inheriting nothing"
    assert ("twin_is_arm", "T") in _role_tests(_roles({"orpheus/numerics/a.py": head + klass("T", "x", "twin", bases="H")}))
    assert ("interface", "I") in _role_tests(_roles({"orpheus/numerics/a.py": head + klass("I", "x", "interface")})), "no Protocol base"
    proto = "from typing import Protocol\nclass ProtocolError(Exception):\n    pass\n"
    assert ("interface", "I") in _role_tests(_roles({"orpheus/numerics/a.py": proto + head + klass("I", "x", "interface", bases="ProtocolError")})), "a substring is not the Protocol"
    assert not _role_tests(_roles({"orpheus/numerics/a.py": proto + head + klass("I", "x", "interface", bases="Protocol")}))
    part = klass("P", "x", "part")
    other = klass("O", "y", "head", body="    p: P")
    assert ("part_no_core_reader", "P") in _role_tests(_roles({"orpheus/numerics/a.py": head + part + other})), "a part read by another concept"
    assert ("part_component_candidate", "P") in _role_tests(_roles({"orpheus/numerics/a.py": head + part + other + klass("A", "x", "arm", bases="H", body="    p: P")}))
    assert not _role_tests(_roles({"orpheus/numerics/a.py": head + part + klass("A", "x", "arm", bases="H", body="    p: P")}))
    comp = klass("C", "x", "component")
    assert ("component_no_core_reader", "C") in _role_tests(_roles({"orpheus/numerics/a.py": head + comp}))
    assert not _role_tests(_roles({"orpheus/numerics/a.py": head + comp + klass("A", "x", "arm", bases="H", body="    c: C")}))


# --- G4 ----------------------------------------------------------------------------------------------------------

def _kinds(sources: dict[str, str], table_text: str):
    m = ct()
    pop = m.population.Population.from_sources(sources)
    edges = frozenset(m.refs.read_edges(sources, pop))
    concepts = m.table.loads_table(table_text)
    return m.checks.check_kinds(pop, concepts, edges, sources), m.checks.check_off_algebra(pop, edges, sources), m.checks.check_covector_lineage(pop, edges)


def _laws(found) -> set[tuple[str, str]]:
    return {(k.law, k.key.qualname) for k in found}


def test_check_kinds_planted_laws():
    t = toml_row("x", "element") + toml_row("r", "record") + toml_row("a", "apparatus", None, "memo", wave=None)
    base = klass("H", "x", "head", "element") + klass("N", "r", "annex", "record")
    assert not _kinds({"orpheus/numerics/a.py": base}, t)[0]
    assert ("head_arm_kind", "A") in _laws(_kinds({"orpheus/numerics/a.py": base + klass("A", "x", "arm", "set", bases="H")}, t)[0]), "law 1"
    assert ("annex_record", "N") in _laws(_kinds({"orpheus/numerics/a.py": klass("H", "x", "head", "element") + klass("N", "r", "annex", "element")}, t)[0]), "law 2"
    assert ("apparatus_row", "H") in _laws(_kinds({"orpheus/numerics/a.py": klass("H", "x", "head", "apparatus")}, t)[0]), "law 4"
    member = "from orpheus.numerics.space import Member\n"
    law5 = {"orpheus/numerics/space.py": klass("Member", "x", "head", "element", "member") + "", "orpheus/numerics/b.py": "from orpheus.numerics.space import Member\n" + klass("F", "x", "arm", "element", bases="Member")}
    assert ("element_subkind", "F") in _laws(_kinds(law5, t)[0]), "law 5, the guard"


def test_check_off_algebra_and_covector_lineage_planted():
    t = toml_row("x", "morphism") + toml_row("y", "element")
    op = "class LinearOperator:  " + decl("x", "weld", "morphism.operator", "x, y") + "\n    def apply(self):\n        pass\n"
    ok = {"orpheus/numerics/operator.py": op, "orpheus/numerics/b.py": "from orpheus.numerics.operator import LinearOperator\n" + klass("B", "x", "arm", "morphism.operator", bases="LinearOperator")}
    _, off, cov = _kinds(ok, t)
    assert not off and not cov
    bad = dict(ok)
    bad["orpheus/numerics/c.py"] = klass("Q", "x", "arm", "morphism.operator")
    assert {o.key.qualname for o in _kinds(bad, t)[1]} == {"Q"}, "converse of law 3: an operator off the algebra"
    missing = dict(ok)
    missing["orpheus/numerics/b.py"] = "from orpheus.numerics.operator import LinearOperator\n" + klass("B", "x", "arm", "morphism", bases="LinearOperator")
    assert ("operator_subkind", "B") in _laws(_kinds(missing, t)[0]), "forward half of law 3"
    cv = {"orpheus/numerics/space.py": klass("Covector", "y", "head", "element.covector"), "orpheus/numerics/b.py": klass("W", "y", "arm", "element.covector")}
    assert {c.key.qualname for c in _kinds(cv, t)[2]} == {"W"}, "converse of law 6"
    wrong = {"orpheus/numerics/space.py": klass("Covector", "y", "head", "element.member")}
    assert ("covector_subkind", "Covector") in _laws(_kinds(wrong, t)[0]), "forward half of law 6"


# --- G5, G5L, G6, G7 ---------------------------------------------------------------------------------------------

PHASED = (toml_row("s", phase="structure") + toml_row("p", phase="posing") + toml_row("o", phase="solution")
          + toml_row("m", "apparatus", None, "memo", wave=None))
SRC5 = {"orpheus/numerics/a.py": klass("S", "s", "head") + klass("P", "p", "head") + klass("O", "o", "head"),
        "orpheus/memo/m.py": klass("M", "m", "head", "apparatus")}


@functools.lru_cache(maxsize=None)
def _k() -> dict:
    return {n: key(f, n) for f, n in (("orpheus/numerics/a.py", "S"), ("orpheus/numerics/a.py", "P"), ("orpheus/numerics/a.py", "O"), ("orpheus/memo/m.py", "M"))}


def _phase(edges):
    m = ct()
    pop = m.population.Population.from_sources(SRC5)
    return m.checks.check_phase(tuple(m.refs.Edge(*e) for e in edges), pop, m.table.loads_table(PHASED))


def test_check_phase_planted_and_the_reversed_reading():
    m = ct()
    K = _k()
    assert _phase([(K["S"], K["O"], "top")]) == frozenset({m.debt.EdgeDebt(K["S"], K["O"], "top")}), "structure referencing solution"
    assert not _phase([(K["O"], K["S"], "top")]), "negative reading: solution to structure is legal"
    assert not _phase([(K["S"], K["M"], "top"), (K["M"], K["O"], "top")]), "an apparatus end is exempt"
    assert not _phase([(K["S"], K["S"], "top")])


def test_check_layer_planted_and_may_reference():
    m = ct()
    assert m.checks.may_reference("numerics", "numerics", FORBIDDEN) and m.checks.may_reference("sn", "numerics", FORBIDDEN)
    assert not m.checks.may_reference("numerics", "sn", FORBIDDEN) and m.checks.may_reference("unknownpkg", "sn", FORBIDDEN)
    pop = m.population.Population.from_sources({"orpheus/numerics/a.py": klass("A", "s", "head"), "orpheus/sn/b.py": klass("B", "s", "arm")})
    ka, kb = key("orpheus/numerics/a.py", "A"), key("orpheus/sn/b.py", "B")
    assert m.checks.check_layer((m.refs.Edge(ka, kb, "top"),), pop, FORBIDDEN) == frozenset({m.debt.EdgeDebt(ka, kb, "top")})
    assert not m.checks.check_layer((m.refs.Edge(kb, ka, "inherits"),), pop, FORBIDDEN)


def test_check_placement_planted():
    m = ct()
    t = m.table.loads_table(toml_row("x", layer="numerics"))
    good = m.population.Population.from_sources({"orpheus/numerics/a.py": klass("H", "x", "head") + klass("A", "x", "arm")})
    assert not m.checks.check_placement(good, t, FORBIDDEN)
    far = m.population.Population.from_sources({"orpheus/sn/a.py": klass("H", "x", "head")})
    assert {v.key.qualname for v in m.checks.check_placement(far, t, FORBIDDEN)} == {"H"}, "a head outside its layer"
    blocked = m.population.Population.from_sources({"orpheus/numerics/a.py": klass("H", "x", "head"), "orpheus/geometry/b.py": klass("A", "x", "arm")})
    assert {v.key.qualname for v in m.checks.check_placement(blocked, t, FORBIDDEN)} == {"A"}, "a row whose package may not import the layer"
    above = m.population.Population.from_sources({"orpheus/numerics/a.py": klass("H", "x", "head"), "orpheus/sn/b.py": klass("A", "x", "arm")})
    assert not m.checks.check_placement(above, t, FORBIDDEN), "a row above the layer may import it: not a violation"


def _cycles(edges, text=toml_row("a") + toml_row("b") + toml_row("c", phase="solution")):
    m = ct()
    src = {"orpheus/numerics/a.py": klass("A", "a", "head") + klass("B", "b", "head") + klass("C", "c", "head")}
    pop = m.population.Population.from_sources(src)
    ks = {n: key("orpheus/numerics/a.py", n) for n in "ABC"}
    return m.checks.check_cycles(tuple(m.refs.Edge(ks[s], ks[d], k) for s, d, k in edges), pop, m.table.loads_table(text))


def test_check_cycles_planted():
    m = ct()
    assert _cycles([("A", "B", "top"), ("B", "A", "top")]) == frozenset({m.debt.CycleDebt(frozenset({"a", "b"}))}), "green on G5, red on G7"
    assert not _cycles([("A", "B", "top")]), "negative reading: one direction only"
    assert _cycles([("A", "B", "inherits"), ("B", "A", "inherits")]), "the quotient half also sees an inheritance cycle inside one phase"
    assert not _cycles([("A", "C", "top"), ("C", "A", "top")]), "an edge rising in phase is G5's, and is not counted in the quotient"


def _inheritance_cycles(edges, text=toml_row("a") + toml_row("b") + toml_row("c", phase="solution")):
    m = ct()
    src = {"orpheus/numerics/a.py": klass("A", "a", "head") + klass("B", "b", "head") + klass("C", "c", "head")}
    pop = m.population.Population.from_sources(src)
    ks = {n: key("orpheus/numerics/a.py", n) for n in "ABC"}
    return m.checks.check_inheritance_cycles(tuple(m.refs.Edge(ks[s], ks[d], k) for s, d, k in edges), pop)


def test_check_inheritance_cycles_planted():
    """The strict half of T11 has a witness of its own, and no debt entry can pin it."""
    m = ct()
    assert _inheritance_cycles([("A", "B", "inherits"), ("B", "A", "inherits")]) == frozenset({frozenset({"a", "b"})})
    assert not _inheritance_cycles([("A", "B", "inherits")]), "negative reading: one direction only"
    assert not _inheritance_cycles([("A", "B", "top"), ("B", "A", "top")]), "a cycle of other kinds is the quotient's"
    # A cycle across two phases: the quotient half drops the phase-violating edge and is silent, so only this half sees it.
    cross = [("A", "C", "inherits"), ("C", "A", "inherits")]
    assert _inheritance_cycles(cross) == frozenset({frozenset({"a", "c"})})
    assert not _cycles(cross), "the quotient half does not see the cross-phase inheritance cycle"
    # A pinned quotient set equal to the inheritance component does not absorb it: the strict half returns a set that no
    # pin touches, so the gate's assertion that it is empty stays red.
    quotient = _cycles([("A", "B", "top"), ("B", "A", "top")])
    assert m.debt.pin(quotient, quotient) == (frozenset(), frozenset())
    assert _inheritance_cycles([("A", "B", "inherits"), ("B", "A", "inherits")]), "still red beside the pin"


# --- debt equality, the wave rule and the unweld rule ------------------------------------------------------------

def test_wave_of_is_a_frozenset_of_the_entry_concepts_waves():
    m = ct()
    t = m.table.loads_table(toml_row("a", wave="W2") + toml_row("b", wave=None) + toml_row("c", wave="W5"))
    src = {"orpheus/numerics/a.py": klass("A", "a", "head") + klass("B", "b", "head") + klass("C", "c", "head")}
    pop = m.population.Population.from_sources(src)
    ka, kb, kc = (key("orpheus/numerics/a.py", n) for n in "ABC")
    w = lambda e: m.debt.wave_of(e, t, pop)
    assert w(m.debt.ClassDebt(ka)) == frozenset({"W2"}) and w(m.debt.ClassDebt(kb)) == frozenset()
    assert w(m.debt.EdgeDebt(kc, ka, "top")) == frozenset({"W5"})
    assert w(m.debt.CycleDebt(frozenset({"a", "b", "c"}))) == frozenset({"W2", "W5"})
    d = m.debt.loads_debt('[[twin]]\nclass = "orpheus/numerics/a.py::B"\n')
    assert m.checks.check_wave_rule(d, t, pop) == frozenset({m.debt.ClassDebt(kb)}), "a concept with no wave reddens the rule"


def _unweld(old_body: str, new_body: str, role_new: str):
    m = ct()
    old = klass("T", "x", "twin", body=old_body)
    new = klass("T", "x", role_new, body=new_body)
    fp = m.checks.fingerprint(old, "T")
    seed = (m.debt.SeedEntry(key("orpheus/a.py", "T"), fp),)
    src = {"orpheus/a.py": new}
    return m.checks.check_unweld(m.population.Population.from_sources(src), seed, src)


def test_unweld_needs_code_motion():
    m = ct()
    k = key("orpheus/a.py", "T")
    assert _unweld("    pass", "    pass", "arm")[1] == frozenset({m.debt.ClassDebt(k)}), "relabel with an unchanged body"
    assert _unweld('    """one line."""\n    x = 1', '    """other line."""\n    x = 1', "arm")[1] == frozenset({m.debt.ClassDebt(k)}), "docstring-only edit is not motion"
    assert not _unweld("    x = 1", "    x = 1\n    y = 2", "arm")[1], "one added statement is motion"
    assert not _unweld('    """one line."""\n    x = 1', '    """other line."""\n    x = 1', "twin")[1], "a twin that stays may change freely"
    fp = m.checks.fingerprint(klass("T", "x", "twin"), "T")
    assert fp == m.checks.fingerprint(klass("T", "y", "twin", "set"), "T"), "the declaration comment is not in the fingerprint"
    assert fp.startswith("blake2b:")
    src = {"orpheus/a.py": klass("T", "x", "twin")}
    uncovered, _ = m.checks.check_unweld(m.population.Population.from_sources(src), (), src)
    assert uncovered == frozenset({m.debt.ClassDebt(k)}), "arm 1: a listed class with no seed entry"
    assert m.checks.check_unweld(m.population.Population.from_sources({}), (m.debt.SeedEntry(k, fp),), {})[1] == frozenset(), "a deleted class is green"
    seed = (m.debt.SeedEntry(k, m.checks.fingerprint(klass("T", "x", "twin", body="    def f(self) -> 'T':\n        return T()"), "T")),)
    for file, name in (("orpheus/b.py", "T"), ("orpheus/a.py", "U"), ("orpheus/b.py", "U")):
        src = {file: klass(name, "x", "arm", body=f"    def f(self) -> 'T':\n        return {name}()")}
        got = m.checks.check_unweld(m.population.Population.from_sources(src), seed, src)[1]
        assert got == frozenset({m.debt.ClassDebt(key(file, name))}), f"an unchanged body is found as {file}::{name}"
    src = {"orpheus/b.py": klass("U", "x", "arm", body="    def f(self) -> 'T':\n        return U() or None")}
    assert not m.checks.check_unweld(m.population.Population.from_sources(src), seed, src)[1], "a moved, renamed and edited class is motion"
