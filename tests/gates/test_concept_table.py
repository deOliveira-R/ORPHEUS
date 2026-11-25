"""Tree gates of the concept table: every class declares its concept, and the laws of the ontology hold of the table.

Reads the tracked table ``docs/architecture/concept_table.toml``, the pinned debt ``docs/architecture/concept_debt.toml``,
the fingerprints ``docs/architecture/concept_seed.toml`` and the declaration (a trailing comment on the line that ends
each class header) of every class of ``orpheus/`` outside ``derivations/``. The gates are G0 (the table is well formed),
G1 (every class declares one concept), G2 (every concept has one head), G3 (the structural role tests), G4 (the kind
laws), G5 (the phase order of references), G5L (the layer order read per class), G6 (placement), G7 (cycles), the debt
equality with its wave rule, and the unweld rule (a class leaves the twin or weld list only with code motion).

A debt list is compared with what the gates measure by equality (``pin``): a new violation is unpinned and a pinned entry
that no longer violates is stale, and both are red. That a list never grows without a ruling is a review rule: every
change of a list is a visible edit of a tracked file.

Mutations are in-memory edits of the sources, the edges or the debt, each opening with a precondition that raises
``Uninstallable`` when the edit did not bite, so a mutation that changed nothing cannot read as a red. The names taken
from ``tools.concept_table`` are listed in ``test_concept_table_machinery``.
"""

from __future__ import annotations

import dataclasses
import functools
import importlib
import re
import subprocess
import sys
import types
import pytest

from tests.gates.test_concept_table_machinery import ct, key
from tests.gates.test_layer_imports import FORBIDDEN_EDGES
from tools.concept_table.paths import DEBT_FILE, ROOT, SEED_FILE, TABLE_FILE

pytestmark = pytest.mark.foundation

EXPECTED_CONCEPTS = 88


class Uninstallable(Exception):
    """A mutation whose precondition does not hold in today's tree: it would not bite."""


def _sources() -> dict[str, str]:
    return ct().population.read_sources(ROOT)[0]


@functools.lru_cache(maxsize=None)
def world() -> types.SimpleNamespace:
    m = ct()
    sources = _sources()
    pop = m.population.Population.from_sources(sources)
    edges = frozenset(m.refs.read_edges(sources, pop))
    return types.SimpleNamespace(
        sources=sources, pop=pop, edges=edges, concepts=m.table.load_table(ROOT / TABLE_FILE),
        debt=m.debt.load_debt(ROOT / DEBT_FILE), seed=m.debt.load_seed(ROOT / SEED_FILE))


def _known_layers() -> frozenset[str]:
    return frozenset(p.name for p in (ROOT / "orpheus").iterdir() if p.is_dir() and p.name != "__pycache__")


def _declaration_of(sources: dict[str, str], k) -> tuple[int, str]:
    """(line number, marker text) of the declaration of class ``k``; ``Uninstallable`` when it has none."""
    found, _ = ct().declaration.read_declarations(sources[k.file])
    markers = found[k.qualname][2]
    if len(markers) != 1:
        raise Uninstallable(f"{k} has {len(markers)} markers, a mutation of its declaration cannot bite")
    return found[k.qualname][1], markers[0]


def edit(sources: dict[str, str], file: str, qualname: str, **changes: object) -> dict[str, str]:
    """A copy of ``sources`` with the declaration of one class replaced by itself with ``changes`` (``None`` removes it)."""
    k = key(file, qualname)
    line, marker = _declaration_of(sources, k)
    lines = sources[file].split("\n")
    assert marker in lines[line - 1]
    if changes.get("remove"):
        lines[line - 1] = lines[line - 1].replace("  " + marker, "")
    else:
        new = dataclasses.replace(ct().declaration.parse_declaration(marker), **changes)  # type: ignore[arg-type]
        lines[line - 1] = lines[line - 1].replace(marker, new.render())
    if "\n".join(lines) == sources[file]:
        raise Uninstallable(f"editing {k} left the source unchanged")
    return {**sources, file: "\n".join(lines)}


def _measure(sources: dict[str, str]):
    m = ct()
    pop = m.population.Population.from_sources(sources)
    return pop, frozenset(m.refs.read_edges(sources, pop))


def _red(findings, needle: str = "") -> bool:
    return any(needle in str(f) for f in findings) and bool(findings)


# --- G0 ---------------------------------------------------------------------------------------------------------

def test_g0_table_is_well_formed_and_has_89_rows():
    m, w = ct(), world()
    assert len(w.concepts) == EXPECTED_CONCEPTS
    findings = m.checks.check_table(w.concepts, _known_layers(), EXPECTED_CONCEPTS, m.checks.doc_labels(ROOT / "docs"), m.checks.function_resolver(ROOT))
    assert not findings, "\n".join(sorted(map(str, findings)))


def test_g0_mutations_redden_the_table():
    m, w = ct(), world()
    labels, resolve = m.checks.doc_labels(ROOT / "docs"), m.checks.function_resolver(ROOT)
    tensor = next(c for c in w.concepts if c.id == "tensor")
    eq = next((t for t in tensor.theory if t.kind == "eq"), None)
    if eq is None:
        raise Uninstallable("`tensor` carries no eq: theory label")
    broken = tuple(dataclasses.replace(c, theory=tuple(m.table.TheoryLink(t.kind, t.label + "x") if t == eq else t for t in c.theory)) if c is tensor else c for c in w.concepts)
    assert m.checks.check_table(broken, _known_layers(), EXPECTED_CONCEPTS, labels, resolve), "an undefined equation label"
    dropped = tuple(c for c in w.concepts if c.id != "tensor")
    assert m.checks.check_table(dropped, _known_layers(), EXPECTED_CONCEPTS, labels, resolve), "87 rows"
    text = (ROOT / TABLE_FILE).read_text()

    def row_edit(cid: str, old: str, new: str) -> str:
        head = text.find(f'id = "{cid}"')
        end = text.find("[[concept]]", head)
        end = len(text) if end < 0 else end
        if head < 0 or old not in text[head:end]:
            raise Uninstallable(f"the row of {cid!r} has no {old!r}")
        return text[:head] + text[head:end].replace(old, new, 1) + text[end:]

    ruling_line = next(l for l in text[text.index('id = "tensor"'):].split("\n") if l.startswith("ruling = "))
    with pytest.raises(m.table.TableError):
        m.table.loads_table(row_edit("tensor", ruling_line + "\n", ""))  # the ruling is required
    with pytest.raises(m.table.TableError):
        m.table.loads_table(row_edit("memo", 'kind = "apparatus"', 'kind = "apparatus"\nphase = "structure"'))  # an apparatus has no phase


def test_g0_concept_table_is_loadable_with_the_closed_head_forms():
    w = world()
    m = ct()
    forms = {type(c.head).__name__ for c in w.concepts}
    assert forms <= {"ClassHead", "SumHead", "FunctionsHead", "MissingHead"} and len(forms) == 4


# --- G1 ---------------------------------------------------------------------------------------------------------

def test_g1_every_class_declares_one_concept():
    m, w = ct(), world()
    findings = m.checks.check_declarations(w.pop, frozenset(c.id for c in w.concepts))
    assert not findings, f"{len(findings)} of {len(w.pop.rows)} classes:\n" + "\n".join(sorted(map(str, findings))[:40])
    assert len(w.pop.rows) > 0 and not w.pop.orphan_markers


def test_g1_population_by_read_tree_equals_the_sources_walk():
    w = world()
    assert [str(r.key) for r in ct().population.read_tree(ROOT).rows] == [str(r.key) for r in w.pop.rows]


def test_g1_mutations():
    m, w = ct(), world()
    ids = frozenset(c.id for c in w.concepts)
    k = key("orpheus/numerics/space.py", "TensorProductSpace")
    stripped = edit(w.sources, k.file, k.qualname, remove=True)
    pop, _ = _measure(stripped)
    found = m.checks.check_declarations(pop, ids)
    assert len(found) == 1 and "TensorProductSpace" in str(next(iter(found)))
    changed = edit(w.sources, "orpheus/cp/solver.py", "CPResult", concept="nosuch")
    assert _red(m.checks.check_declarations(_measure(changed)[0], ids), "nosuch")


# --- G2 ---------------------------------------------------------------------------------------------------------

def test_g2_every_concept_has_one_head():
    m, w = ct(), world()
    findings = m.checks.check_heads(w.pop, w.concepts, w.sources)
    assert not findings, "\n".join(sorted(map(str, findings)))


def test_g2_mutations():
    m, w = ct(), world()
    flipped = edit(w.sources, "orpheus/numerics/space.py", "TensorProductSpace", role="arm")
    pop, _ = _measure(flipped)
    assert _red(m.checks.check_heads(pop, w.concepts, flipped), "tensor"), "the head of `tensor` flipped to arm"
    promoted = edit(w.sources, "orpheus/numerics/green_operator.py", "Richardson", role="head")
    assert _red(m.checks.check_heads(_measure(promoted)[0], w.concepts, promoted), "Richardson"), "a head beside a missing concept"
    cat = {c.id: c for c in w.concepts}
    assert isinstance(cat["algorithm"].head, ct().table.MissingHead), "precondition of the control"
    heads = [c.id for c in w.concepts if isinstance(c.head, ct().table.ClassHead)]
    assert len(heads) == 66, "control: 66 class-head concepts, one head each"


# --- G3 ---------------------------------------------------------------------------------------------------------

def test_g3_role_tests_equal_the_pinned_debt():
    m, w = ct(), world()
    unpinned, stale = m.debt.pin(m.checks.check_roles(w.pop, w.edges, w.concepts, w.sources), frozenset(w.debt.role_test))
    assert not unpinned and not stale, f"unpinned {sorted(map(str, unpinned))}; stale {sorted(map(str, stale))}"


@pytest.mark.parametrize("file,qualname,role,test", [
    ("orpheus/diffusion/augmented_mesh.py", "DiffusionMesh", "twin", "twin_is_arm"),
    ("orpheus/cp/solver.py", "CPResult", "arm", "arm"),
    ("orpheus/data/energy_grid.py", "WithinGroupSpectrum", "part", "part_no_core_reader"),
    ("orpheus/cp/solver.py", "CPDecks", "component", "component_no_core_reader"),
    ("orpheus/data/materials.py", "Materials", "part", "part_component_candidate"),
])
def test_g3_mutation_flips_a_role(file, qualname, role, test):
    m, w = ct(), world()
    flipped = edit(w.sources, file, qualname, role=role)
    pop, edges = _measure(flipped)
    got = m.checks.check_roles(pop, edges, w.concepts, flipped)
    new = got - frozenset(w.debt.role_test)
    assert any(r.test == test and r.key.qualname == qualname for r in new), (test, sorted(map(str, new)))


# --- G4 ---------------------------------------------------------------------------------------------------------

def test_g4_kind_laws_and_both_converse_lists_equal_the_debt():
    m, w = ct(), world()
    u, s = m.debt.pin(m.checks.check_kinds(w.pop, w.concepts, w.edges, w.sources), frozenset(w.debt.kind_law))
    assert not u and not s, f"kind_law unpinned {sorted(map(str, u))}; stale {sorted(map(str, s))}"
    u, s = m.debt.pin(m.checks.check_off_algebra(w.pop, w.edges, w.sources), frozenset(w.debt.operator_off_algebra))
    assert not u and not s, f"operator_off_algebra unpinned {sorted(map(str, u))}; stale {sorted(map(str, s))}"
    u, s = m.debt.pin(m.checks.check_covector_lineage(w.pop, w.edges), frozenset(w.debt.covector_off_lineage))
    assert not u and not s, f"covector_off_lineage unpinned {sorted(map(str, u))}; stale {sorted(map(str, s))}"


@pytest.mark.parametrize("file,qualname,changes,law", [
    ("orpheus/cp/solver.py", "CPSpecularChain", {"kind": "element"}, "annex_record"),
    ("orpheus/numerics/assembled_operator.py", "SparseAssembledOperator", {"subkind": None}, "operator_subkind"),
    ("orpheus/values/content.py", "ContentIdentity", {"kind": "apparatus"}, "apparatus_row"),
    ("orpheus/numerics/green_operator.py", "Richardson", {"kind": "morphism"}, "head_arm_kind"),
    ("orpheus/numerics/field.py", "Field", {"subkind": None}, "element_subkind"),
    ("orpheus/numerics/space.py", "Covector", {"subkind": "member"}, "covector_subkind"),
])
def test_g4_mutation_reddens_its_law(file, qualname, changes, law):
    m, w = ct(), world()
    mutated = edit(w.sources, file, qualname, **changes)
    pop, edges = _measure(mutated)
    new = m.checks.check_kinds(pop, w.concepts, edges, mutated) - frozenset(w.debt.kind_law)
    assert any(k.law == law and k.key.qualname == qualname for k in new), (law, sorted(map(str, new)))


def test_g4_apparatus_mutation_is_also_a_head_arm_kind_violation():
    """Specification 2.8, G4: the ``ContentIdentity`` mutation (kind apparatus on a row of a concept that is not an apparatus)
    breaks the apparatus law and, the concept having a different kind, the head-or-arm law, so a test asserts both."""
    m, w = ct(), world()
    mutated = edit(w.sources, "orpheus/values/content.py", "ContentIdentity", kind="apparatus")
    pop, edges = _measure(mutated)
    new = m.checks.check_kinds(pop, w.concepts, edges, mutated) - frozenset(w.debt.kind_law)
    assert {k.law for k in new if k.key.qualname == "ContentIdentity"} >= {"apparatus_row", "head_arm_kind"}, sorted(map(str, new))


def test_g4_converse_mutations_name_the_class():
    m, w = ct(), world()
    op = edit(w.sources, "orpheus/numerics/metric.py", "HilbertMetric", subkind="operator")
    pop, edges = _measure(op)
    assert any(v.key.qualname == "HilbertMetric" for v in m.checks.check_off_algebra(pop, edges, op) - frozenset(w.debt.operator_off_algebra))
    cov = edit(w.sources, "orpheus/numerics/field.py", "Field", subkind="covector")
    pop, edges = _measure(cov)
    assert any(v.key.qualname == "Field" for v in m.checks.check_covector_lineage(pop, edges) - frozenset(w.debt.covector_off_lineage))


def test_g4_controls_law_1_passes_axis_and_law_3_hit_set_is_not_empty():
    m, w = ct(), world()
    assert not any(k.key.qualname == "Axis" for k in m.checks.check_kinds(w.pop, w.concepts, w.edges, w.sources))
    on_algebra = [r for r in w.pop.rows if r.declaration and r.declaration.subkind == "operator" and r.declaration.role != "interface"]
    assert len(on_algebra) > len(w.debt.operator_off_algebra) > 0


# --- G5, G5L ----------------------------------------------------------------------------------------------------

def test_g5_phase_order_equals_the_pinned_edges():
    m, w = ct(), world()
    u, s = m.debt.pin(m.checks.check_phase(tuple(w.edges), w.pop, w.concepts), frozenset(w.debt.edge))
    assert not u and not s, f"unpinned {sorted(map(str, u))}; stale {sorted(map(str, s))}"


def test_g5_mutation_and_negative_control():
    m, w = ct(), world()
    fs, ev, cid = (key("orpheus/numerics/space.py", "FunctionSpace"), key("orpheus/numerics/eigenvalue.py", "EigenvalueSolver"),
                   key("orpheus/values/content.py", "ContentIdentity"))
    base = frozenset(w.debt.edge)
    rising = m.refs.Edge(fs, ev, "top")
    if rising in w.edges:
        raise Uninstallable("the edge is already in the reference pass")
    got = m.checks.check_phase(tuple(w.edges | {rising}), w.pop, w.concepts)
    assert m.debt.EdgeDebt(fs, ev, "top") in got - base
    falling = m.refs.Edge(ev, cid, "top")
    assert not (m.checks.check_phase(tuple(w.edges | {falling}), w.pop, w.concepts) - base), "solution to structure is legal"
    assert m.debt.EdgeDebt(ev, cid, "top") not in m.checks.check_phase(tuple(w.edges | {falling}), w.pop, w.concepts)
    assert m.debt.EdgeDebt(cid, ev, "top") in m.checks.check_phase(tuple(w.edges | {m.refs.Edge(cid, ev, "top")}), w.pop, w.concepts), "the reverse reads red"


def test_g5l_no_class_edge_runs_against_forbidden_edges():
    m, w = ct(), world()
    assert len(w.edges) > 0
    assert not m.checks.check_layer(tuple(w.edges), w.pop, FORBIDDEN_EDGES)


def test_g5l_mutation_and_positive_control_of_the_zero():
    m, w = ct(), world()
    fs, sn = key("orpheus/numerics/space.py", "FunctionSpace"), key("orpheus/sn/problem.py", "SNProblem")
    if m.refs.Edge(fs, sn, "top") in w.edges:
        raise Uninstallable("the edge is already in the reference pass")
    got = m.checks.check_layer(tuple(w.edges | {m.refs.Edge(fs, sn, "top")}), w.pop, FORBIDDEN_EDGES)
    assert got == frozenset({m.debt.EdgeDebt(fs, sn, "top")})
    pkg = lambda k: k.file.split("/")[1]
    reversed_reads = sum(1 for e in w.edges if not m.checks.may_reference(pkg(e.dst), pkg(e.src), FORBIDDEN_EDGES))
    assert reversed_reads > 0, "the reversed predicate fires, so the predicate can fire"


# --- G6 ---------------------------------------------------------------------------------------------------------

def test_g6_placement_equals_the_pinned_misplaced_rows():
    m, w = ct(), world()
    u, s = m.debt.pin(m.checks.check_placement(w.pop, w.concepts, FORBIDDEN_EDGES), frozenset(w.debt.misplaced))
    assert not u and not s, f"unpinned {sorted(map(str, u))}; stale {sorted(map(str, s))}"


@pytest.mark.parametrize("file,qualname,moved", [
    ("orpheus/numerics/space.py", "TensorProductSpace", "orpheus/sn/space.py"),
    ("orpheus/numerics/factor.py", "Factor", "orpheus/sn/factor.py"),
    ("orpheus/numerics/space.py", "FunctionSpace", "orpheus/geometry/space.py"),
])
def test_g6_mutation_moves_a_row_out_of_its_layer(file, qualname, moved):
    m, w = ct(), world()
    if moved in w.sources:
        raise Uninstallable(f"{moved} exists")
    src = dict(w.sources)
    src[moved] = src.pop(file)
    pop, _ = _measure(src)
    new = m.checks.check_placement(pop, w.concepts, FORBIDDEN_EDGES)
    assert any(v.key.file == moved for v in new - frozenset(w.debt.misplaced)), qualname


# --- G7 ---------------------------------------------------------------------------------------------------------

def test_g7_cycles_equal_the_pinned_cycles():
    m, w = ct(), world()
    u, s = m.debt.pin(m.checks.check_cycles(tuple(w.edges), w.pop, w.concepts), frozenset(w.debt.cycle))
    assert not u and not s, f"unpinned {sorted(map(str, u))}; stale {sorted(map(str, s))}"


def test_g7_mutation_adds_a_fourth_component():
    m, w = ct(), world()
    a, b = key("orpheus/values/content.py", "ContentIdentity"), key("orpheus/values/content.py", "_ByEquality")
    extra = {m.refs.Edge(a, b, "top"), m.refs.Edge(b, a, "top")}
    if extra & w.edges:
        raise Uninstallable("one of the edges is already in the reference pass")
    got = m.checks.check_cycles(tuple(w.edges | extra), w.pop, w.concepts)
    assert m.debt.CycleDebt(frozenset({"digest", "value"})) in got - frozenset(w.debt.cycle)


def test_g7_inheritance_half_is_acyclic_on_the_tree():
    m, w = ct(), world()
    assert not m.checks.check_inheritance_cycles(tuple(w.edges), w.pop), "no pin can absorb an inheritance cycle"


def test_g7_inheritance_control_moves_rolepair():
    """Positive control of the strict half: moving ``RolePair`` to concept ``fluxsrc`` closes a cycle of inheritance edges,
    and the half returns exactly that component."""
    m, w = ct(), world()
    inherits = tuple(e for e in w.edges if e.kind == "inherits")
    assert inherits, "positive control: the tree has inheritance edges"
    k = key("orpheus/transport/fields/_bases.py", "RolePair")
    moved = edit(w.sources, k.file, k.qualname, concept="fluxsrc")
    pop, edges = _measure(moved)
    assert m.checks.check_inheritance_cycles(tuple(edges), pop) == frozenset({frozenset({"fluxsrc", "member"})})


# --- debt equality, the wave rule, the unweld rule ----------------------------------------------------------------

def test_debt_twin_and_weld_lists_equal_the_declared_roles():
    m, w = ct(), world()
    for measured, pinned, name in ((m.checks.measure_twin(w.pop), w.debt.twin, "twin"), (m.checks.measure_weld(w.pop), w.debt.weld, "weld")):
        u, s = m.debt.pin(measured, frozenset(pinned))
        assert not u and not s, f"{name}: unpinned {sorted(map(str, u))[:10]}; stale {sorted(map(str, s))[:10]}"


def test_debt_wave_rule_every_entry_has_a_wave():
    m, w = ct(), world()
    assert not m.checks.check_wave_rule(w.debt, w.concepts, w.pop)


def test_debt_stale_and_unpinned_entries_are_red():
    m, w = ct(), world()
    first = w.debt.edge[0]
    u, s = m.debt.pin(m.checks.check_phase(tuple(w.edges), w.pop, w.concepts), frozenset(w.debt.edge) - {first})
    assert first in u, "deleting a pinned edge leaves its violation unpinned"
    fake = m.debt.EdgeDebt(key("orpheus/numerics/eigenvalue.py", "EigenvalueSolver"), key("orpheus/values/content.py", "ContentIdentity"), "top")
    u, s = m.debt.pin(m.checks.check_phase(tuple(w.edges), w.pop, w.concepts), frozenset(w.debt.edge) | {fake})
    assert s == frozenset({fake}), "an entry that no longer violates is stale"


def test_unweld_rule_on_the_tree_and_its_mutations():
    m, w = ct(), world()
    uncovered, relabelled = m.checks.check_unweld(w.pop, w.seed, w.sources)
    assert not uncovered and not relabelled
    assert len(w.seed) > 0
    k = key("orpheus/cp/solver.py", "CPResult")
    mutated = edit(w.sources, k.file, k.qualname, role="arm")
    pop, _ = _measure(mutated)
    assert m.debt.ClassDebt(k) in m.checks.check_unweld(pop, w.seed, mutated)[1], "relabel without code motion"
    moved_to = "orpheus/cp/_moved.py"
    assert moved_to not in w.sources
    node = m.checks._class_node(mutated[k.file], k.qualname)
    lines = mutated[k.file].splitlines(keepends=True)
    body = "".join(lines[min([node.lineno, *(d.lineno for d in node.decorator_list)]) - 1:node.end_lineno])
    moved = {**w.sources, k.file: mutated[k.file].replace(body, ""), moved_to: body}
    pop, _ = _measure(moved)
    assert m.debt.ClassDebt(key(moved_to, k.qualname)) in m.checks.check_unweld(pop, w.seed, moved)[1], "relabel and move to another file"
    renamed = {**w.sources, k.file: re.sub(r"\bCPResult\b", "CPOutcome", mutated[k.file])}
    pop, _ = _measure(renamed)
    assert m.debt.ClassDebt(key(k.file, "CPOutcome")) in m.checks.check_unweld(pop, w.seed, renamed)[1], "relabel and rename"
    seed = tuple(s for s in w.seed if s.key != k)
    assert len(seed) == len(w.seed) - 1, "precondition: CPResult has a seed entry"
    assert m.debt.ClassDebt(k) in m.checks.check_unweld(w.pop, seed, w.sources)[0], "arm 1: a listed class with no seed entry"


def test_seed_fingerprints_reproduce_on_the_tree():
    m, w = ct(), world()
    listed = {d.key for d in m.checks.measure_twin(w.pop) | m.checks.measure_weld(w.pop)}
    on_list = [s for s in w.seed if s.key in listed]
    bad = [str(s.key) for s in on_list if m.checks.fingerprint(w.sources[s.key.file], s.key.qualname) != s.fingerprint]
    assert len(on_list) > 0 and not bad, f"{len(bad)} of {len(on_list)} fingerprints of listed classes differ: {bad[:5]}"
    assert (ROOT / SEED_FILE).read_text() == importlib.import_module("tools.concept_table.__main__").reseed(ROOT), "concept_seed.toml is not what `python -m tools.concept_table seed` writes"


def test_g2_g3_control_4_modal_axis_edges_and_components():
    """The concept half of reference-pass control 4: the alias expansion is what gives these readers a core reader."""
    m, w = ct(), world()
    names = {(e.src.qualname, e.dst.qualname, e.kind) for e in w.edges}
    assert {("Axis", "Modal", "top"), ("Axis", "FromFrame", "top"), ("InflowTerms", "FaceFrameInflow", "top")} <= names
    failing = {r.key.qualname for r in m.checks.check_roles(w.pop, w.edges, w.concepts, w.sources)}
    assert not failing & {"Modal", "FromFrame", "FaceFrameInflow"}


# --- the report -------------------------------------------------------------------------------------------------

def test_report_first_line():
    """The line is parsed, not matched: the counts are printed for the reader and asserted only where each equals the
    population it counts, so a correctly declared new class does not redden CI (specification 2.7, 2.11 item 1)."""
    m, w = ct(), world()
    out = subprocess.run([sys.executable, "-m", "tools.concept_table", "report"], cwd=ROOT, capture_output=True, text=True, check=True).stdout
    words = out.splitlines()[0].split()
    fields = {name: int(n) for name, n in zip(words[0::2], words[1::2], strict=True)}
    lists = w.debt.lists()
    assert fields["classes"] == fields["declared"] == len(w.pop.rows)
    assert fields["concepts"] == EXPECTED_CONCEPTS
    assert {n: fields[n] for n in lists} == {n: len(entries) for n, entries in lists.items()}
    assert len(lists) == 9
