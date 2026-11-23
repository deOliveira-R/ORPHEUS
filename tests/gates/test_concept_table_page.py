"""G8, the generated concept page is the table: ``docs/architecture/concept_table.rst`` is ``render`` of the tracked inputs.

The page is written by ``tools/concept_table/page.py`` from three tracked inputs (the table
``docs/architecture/concept_table.toml``, the pinned debt ``concept_debt.toml`` and the classes of ``orpheus/`` read
by ``read_tree``) and committed. This file asserts, one claim per test:

* the committed page equals ``render(*load_inputs(ROOT))`` byte for byte, and ``--check`` agrees under two hash seeds;
* ``render`` is invariant under any permutation of the table rows, the debt entries, the members of a cycle and the
  order the sources are read in;
* the member bullets of the page are, as a set of ``(file, qualname)``, the classes of ``orpheus/`` outside
  ``derivations/`` read by this file's own ``ast`` walk, which shares no code with ``read_tree`` or ``walk_classes``;
* every concept has its ``.. _<slug>:`` label on a heading, no class or function carries a Sphinx role, and every
  theory link of the table is on the page with its explicit link text (the conservation of the replaced hand-written
  table: 107 occurrences on 102 distinct targets, 36 ``functions`` entries; the multiset equality with the old table
  was checked once by the builder and is recorded in its commit message, and no tracked gate reads ``scratch/``);
* the enumerations the page shows are the code's own (``tuple(ROLE_INFO) == ROLES``);
* the generator is registered in ``docs/conf.py``, the toctree and the ``gates`` workflow name the page, and the hand
  written table of ``conceptual_view.rst`` is replaced by a pointer.

Each mutation is an in-memory edit of an input (or an edit of a temporary copy, for ``--check``), and opens with a
precondition that raises ``Uninstallable`` when the edit did not change the input, so a mutation that bit nothing
cannot read as a red. The tracked tree is never edited. The names taken from ``tools.concept_table.page`` are
``load_inputs(root)``, ``render(table, debt, population)``, ``main(argv)`` and ``ROLE_INFO``. ``page`` is imported by
the test that needs it, so each test fails alone with the module missing.

Reddening mutations, one per assertion (each is exercised by the test it is named in):

* (a) a field of one concept row changed in a copy of the table (``test_page_reddens_on_a_concept_field``);
* (b) a class dropped from a copy of the population (``test_member_lists_redden_on_a_dropped_class``);
* (c) a debt entry removed from a copy of the debt (``test_page_reddens_on_a_removed_debt_entry``).

The same three are run through ``python -m tools.concept_table.page --check`` on a temporary copy of the tree in
``test_check_mode_exits_nonzero_on_each_mutation``.
"""

from __future__ import annotations

import ast
import dataclasses
import functools
import importlib
import re
import shutil
import subprocess
import sys
import types
from pathlib import Path

import pytest

from tools.concept_table.paths import DEBT_FILE, ROOT, TABLE_FILE

pytestmark = pytest.mark.foundation

PAGE_FILE = "docs/architecture/concept_table.rst"
MODULE = "tools.concept_table.page"
# The conservation counts of the replaced hand-written table (specification C 3.3), over the tracked table.
THEORY_OCCURRENCES_CARRIED = 107
THEORY_TARGETS_CARRIED = 102
FUNCTIONS_ENTRIES = 36
BULLET = re.compile(r"^- ``(?P<qualname>[^`]+)`` \(``(?P<file>[^`]+\.py)``\) ``(?P<kind>[a-z.]+)``(?: \[welded: (?P<welded>[^\]]*)\])?$")
PYTHON_ROLES = re.compile(r":(?:py:)?(?:class|func|meth|mod|attr|obj|data|exc):`")


class Uninstallable(Exception):
    """A mutation whose precondition does not hold in today's tree: it would not bite."""


def page() -> types.ModuleType:
    """The generator, imported by the test that needs it."""
    return importlib.import_module(MODULE)


@functools.lru_cache(maxsize=None)
def inputs():
    """``(table, debt, population)`` read from the tracked tree by the generator's own loader."""
    return page().load_inputs(ROOT)


def committed() -> str:
    path = ROOT / PAGE_FILE
    if not path.is_file():
        raise AssertionError(f"{PAGE_FILE} is not committed: run python -m {MODULE}")
    return path.read_text(encoding="utf-8")


def first_difference(a: str, b: str) -> str:
    for n, (x, y) in enumerate(zip(a.split("\n"), b.split("\n")), 1):
        if x != y:
            return f"line {n}: committed {x!r}, rendered {y!r}"
    return f"lengths differ: {len(a)} against {len(b)}"


def rendered(table=None, debt=None, population=None) -> str:
    t, d, p = inputs()
    return page().render(t if table is None else table, d if debt is None else debt, p if population is None else population)


# --- the page is the table --------------------------------------------------------------------------------------


def test_page_equals_render_of_the_tracked_inputs() -> None:
    """The committed page is ``render`` of the tracked inputs, byte for byte, UTF-8, ``\\n`` line ends, one final newline.

    Reddens: the page is not committed, or any edit of an input that the page shows (the three mutations below are
    the named instances, ``a`` a concept field, ``b`` a class of the population, ``c`` a debt entry).
    """
    text = rendered()
    raw = (ROOT / PAGE_FILE).read_bytes() if (ROOT / PAGE_FILE).is_file() else b""
    assert raw.decode("utf-8") == committed() and b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n")
    assert text == committed(), (
        f"{PAGE_FILE} is stale: run python -m {MODULE}; first difference, {first_difference(committed(), text)}")


def test_page_opens_with_the_generator_marker() -> None:
    """Its first line is an ``.rst`` comment naming the module that generates it (a generated file says so).

    Reddens: a first line that does not name ``tools.concept_table.page``.
    """
    first = committed().split("\n", 1)[0]
    assert first.startswith("..") and MODULE in first, first


@pytest.mark.parametrize("field,edit", [
    ("definition", lambda c: dataclasses.replace(c, definition=c.definition + " (edited)")),
    ("name", lambda c: dataclasses.replace(c, name=c.name + " edited")),
    ("ruling", lambda c: dataclasses.replace(c, ruling="O1" if c.ruling != "O1" else "O2")),
])
def test_page_reddens_on_a_concept_field(field, edit) -> None:
    """Mutation (a): a field of the concept ``digest`` changed in a copy of the table makes ``render`` differ from the page.

    Control: the unmutated copy renders to the committed page (the reading is not trivially red). Precondition: the
    edited row differs from the original, so a mutation that bit nothing cannot pass.
    """
    table, _, _ = inputs()
    original = next(c for c in table if c.id == "digest")
    mutant = edit(original)
    if mutant == original:
        raise Uninstallable(f"editing {field} of digest changed nothing")
    changed = tuple(mutant if c.id == "digest" else c for c in table)
    assert rendered(table=tuple(table)) == committed()
    assert rendered(table=changed) != committed(), f"the page does not show the {field} of a concept"


@pytest.mark.parametrize("name", ["role_test", "misplaced", "edge", "cycle", "kind_law", "operator_off_algebra", "covector_off_lineage"])
def test_page_reddens_on_a_removed_debt_entry(name) -> None:
    """Mutation (c): a debt entry removed from a copy of the debt makes ``render`` differ (every list but ``twin`` and ``weld`` is shown by name).

    Precondition: the list is non-empty, so the removal bites. Control: the unmutated debt renders to the page.
    """
    _, debt, _ = inputs()
    entries = getattr(debt, name)
    if not entries:
        raise Uninstallable(f"debt list {name} is empty, nothing to remove")
    mutant = dataclasses.replace(debt, **{name: entries[1:]})
    assert rendered(debt=debt) == committed()
    assert rendered(debt=mutant) != committed(), f"the page does not show the entries of the {name} list"


@pytest.mark.parametrize("name", ["twin", "weld"])
def test_page_counts_the_twin_and_weld_lists(name) -> None:
    """The ``twin`` and ``weld`` lists are shown by count only (their membership is the member roles); the count is on the page.

    Reddens: a removal from either list that leaves the page unchanged.
    """
    _, debt, _ = inputs()
    entries = getattr(debt, name)
    if not entries:
        raise Uninstallable(f"debt list {name} is empty")
    assert rendered(debt=dataclasses.replace(debt, **{name: entries[1:]})) != committed()


# --- invariance (F6) --------------------------------------------------------------------------------------------


def test_render_is_invariant_under_permutation_of_its_inputs() -> None:
    """``render(sigma . inputs) == render(inputs)`` for table rows, debt entries, cycle members and reversed source order.

    Reddens: an implementation that iterates a list or a ``frozenset`` in input order. The reversed-sources arm
    tests the population's sorting and not ``render`` (the population sorts), the other arms test ``render``.
    """
    from tools.concept_table.debt import CycleDebt
    from tools.concept_table.population import Population, read_sources

    table, debt, population = inputs()
    sources, skipped = read_sources(ROOT)
    reversed_sources = Population.from_sources(dict(reversed(list(sources.items()))), skipped)
    assert reversed_sources == population
    shuffled = dataclasses.replace(debt, **{
        f.name: tuple(reversed(getattr(debt, f.name))) for f in dataclasses.fields(debt)})
    shuffled = dataclasses.replace(shuffled, cycle=tuple(
        CycleDebt(frozenset(sorted(c.concepts, reverse=True))) for c in shuffled.cycle))
    if shuffled == debt:
        raise Uninstallable("the permutation of the debt changed nothing")
    assert rendered(table=tuple(reversed(table)), debt=shuffled, population=reversed_sources) == rendered()


@pytest.mark.parametrize("seed", ["0", "1"])
def test_check_mode_is_hash_seed_independent(seed) -> None:
    """``python -m tools.concept_table.page --check`` exits 0 under ``PYTHONHASHSEED`` 0 and 1 (the build hook runs it in a fresh process).

    Reddens: a dependence on set iteration order (the 39-concept cycle iterates differently under the two seeds), or a stale page.
    """
    done = subprocess.run([sys.executable, "-m", MODULE, "--check"], cwd=ROOT, capture_output=True, text=True,
                          env={"PYTHONHASHSEED": seed, "PATH": "/usr/bin:/bin"})
    assert done.returncode == 0, done.stdout + done.stderr


# --- the member lists are the population, read independently ----------------------------------------------------


def ast_classes(root: Path) -> frozenset[tuple[str, str]]:
    """``(file, qualname)`` of every class of ``root/orpheus`` outside ``derivations/``, by this file's own walk.

    The qualname is Python's: ``Outer.Inner`` for a nested class, ``f.<locals>.Name`` for a class local to a function.
    Shares no code with ``read_tree`` or ``walk_classes`` (X4: the page's input and the gate's oracle are two walkers).
    """
    found: set[tuple[str, str]] = set()

    def visit(node: ast.AST, prefix: str, file: str) -> None:
        for child in ast.iter_child_nodes(node):
            if isinstance(child, ast.ClassDef):
                found.add((file, prefix + child.name))
                visit(child, prefix + child.name + ".", file)
            elif isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
                visit(child, prefix + child.name + ".<locals>.", file)
            else:
                visit(child, prefix, file)

    for path in (root / "orpheus").rglob("*.py"):
        rel = path.relative_to(root)
        if "__pycache__" in rel.parts or "derivations" in rel.parts:
            continue
        visit(ast.parse(path.read_text()), "", rel.as_posix())
    return frozenset(found)


def bullets(text: str) -> list[tuple[str, str]]:
    return [(m["file"], m["qualname"]) for line in text.split("\n") if (m := BULLET.match(line))]


def test_walk_control_finds_known_members() -> None:
    """Positive control of the independent walk: it finds a module-level class and the one function-local class (``TaskInfo``)."""
    walked = ast_classes(ROOT)
    assert any(q == "TaskInfo" or q.endswith(".<locals>.TaskInfo") for _, q in walked)
    assert sum("<locals>" in q for _, q in walked) == 1
    assert len(walked) == 676


def test_member_lists_equal_the_population() -> None:
    """The ``(file, qualname)`` of the page's member bullets equal the classes walked by ``ast``: 676 of 676, none extra, none twice.

    Reddens: a class missing from the page, a class listed that is not in the tree, a class listed under two concepts.
    """
    listed = bullets(committed())
    walked = ast_classes(ROOT)
    assert len(listed) == len(set(listed)), "a class is listed twice"
    assert set(listed) == walked, (sorted(set(listed) - walked)[:5], sorted(walked - set(listed))[:5])


def test_member_lists_redden_on_a_dropped_class() -> None:
    """Mutation (b): a class dropped from a copy of the population leaves the rendered member bullets short of the ``ast`` walk.

    The dropped class is a ``part``, ``component`` or ``arm`` row that no ``missing`` head names as an arm, so ``render``
    still succeeds. Precondition: the copy has one row fewer. Control: the unmutated population renders the walked set.
    """
    table, _, population = inputs()
    named = {k for c in table if hasattr(c.head, "arms") for k in c.head.arms}
    victim = next(r for r in population.rows if r.declaration is not None and r.declaration.role in ("part", "component", "arm")
                  and r.key not in named)
    mutant = dataclasses.replace(population, rows=tuple(r for r in population.rows if r.key != victim.key))
    if len(mutant.rows) != len(population.rows) - 1:
        raise Uninstallable("the population copy did not lose a row")
    walked = ast_classes(ROOT)
    assert set(bullets(rendered())) == walked
    assert set(bullets(rendered(population=mutant))) != walked


def test_every_concept_has_its_label_on_a_heading() -> None:
    """Each concept has ``.. _<Concept.slug>:`` once, followed by a title whose text is the concept id (slug already begins ``concept-``).

    Reddens: a missing label, a doubled prefix (``concept-concept-``), a label on a table cell and not a title.
    """
    lines = committed().split("\n")
    table, _, _ = inputs()
    assert not any("concept-concept-" in line for line in lines)
    for c in table:
        at = [n for n, line in enumerate(lines) if line == f".. _{c.slug}:"]
        assert len(at) == 1, (c.id, at)
        rest = [line for line in lines[at[0] + 1:at[0] + 5] if line.strip()]
        assert rest[0] == c.id and set(rest[1]) <= set("-=~^\"'`#*+") and len(rest[1]) >= len(c.id), (c.id, rest[:2])


def test_page_label_and_pointer_to_the_debt_heading() -> None:
    """The page carries ``architecture-concept-table`` and links ``architecture-hand-rolled-debt``, which ``conceptual_view.rst`` defines above its debt heading."""
    assert ".. _architecture-concept-table:" in committed()
    assert ":ref:`" in committed() and "architecture-hand-rolled-debt" in committed()
    view = (ROOT / "docs/architecture/conceptual_view.rst").read_text()
    assert re.search(r"^\.\. _architecture-hand-rolled-debt:\n\n+What is still hand-rolled", view, re.M)


def test_classes_and_functions_carry_no_sphinx_role() -> None:
    """Classes and ``functions`` entries are literal names with file paths; no Python-domain role appears (222 of 676 classes resolve, an unresolved role is a silent dead link).

    Reddens: any ``:class:``, ``:func:``, ``:meth:``, ``:mod:`` ... role on the page.
    """
    text = committed()
    assert not PYTHON_ROLES.search(text), PYTHON_ROLES.search(text)
    table, _, _ = inputs()
    for c in table:
        for entry in c.functions:
            assert f"``{entry}``" in text, entry


def test_theory_links_are_conserved_on_the_page() -> None:
    """Conservation (3.3): the table carries 107 theory occurrences on 102 distinct targets, each on one concept, each on the page with its explicit link text.

    Pinned counts are the measurement of the integration; a change of a count is a ruled edit of the table. The multiset
    equality with the replaced table was checked once by the builder (recorded in the commit message).
    Reddens: a target on two concepts, a target missing from the page, a ``ref`` without link text (a label on a paragraph fails ``-W``).
    """
    table, _, _ = inputs()
    owner: dict[str, set[str]] = {}
    for c in table:
        for t in c.theory:
            owner.setdefault(t, set()).add(c.id)
    assert sum(len(c.theory) for c in table) == THEORY_OCCURRENCES_CARRIED
    assert len(owner) == THEORY_TARGETS_CARRIED
    assert all(len(v) == 1 for v in owner.values()), sorted(t for t, v in owner.items() if len(v) > 1)
    assert sum(len(c.functions) for c in table) == FUNCTIONS_ENTRIES
    text = committed()
    for t in owner:
        kind, label = t.split(":", 1)
        shown = f":ref:`{label} <{label}>`" if kind == "ref" else f":eq:`{label}`"
        assert shown in text, t


def test_enumerations_are_the_codes_own() -> None:
    """``tuple(ROLE_INFO) == ROLES``, and the page shows every role, kind, phase, gloss and gate from the code's own constants (never re-typed).

    Reddens: a role added to ``ROLES`` with no gloss, a gloss that falls out of the page.
    """
    from tools.concept_table.declaration import KINDS, ROLES, SUBKINDS
    from tools.concept_table.table import PHASE_ORDER

    mod = page()
    assert tuple(mod.ROLE_INFO) == ROLES
    text = committed()
    for role, (gloss, gate) in mod.ROLE_INFO.items():
        assert gloss in text and gate in text, role
    for word in (*KINDS, *PHASE_ORDER, *(s for subs in SUBKINDS.values() for s in subs)):
        assert word in text, word


# --- the unit's wiring ------------------------------------------------------------------------------------------


def test_generator_is_registered_in_conf() -> None:
    """``docs/conf.py`` ``_GENERATORS`` holds ``("tools.concept_table.page", "concept table", "builder-inited")`` (the mechanism of the verification matrix)."""
    tree = ast.parse((ROOT / "docs/conf.py").read_text())
    value = next(n.value for n in tree.body if isinstance(n, ast.Assign) and any(
        isinstance(t, ast.Name) and t.id == "_GENERATORS" for t in n.targets))
    assert ("tools.concept_table.page", "concept table", "builder-inited") in ast.literal_eval(value)


def test_page_is_in_the_toctree_after_the_conceptual_view() -> None:
    """``docs/architecture/index.rst`` lists ``concept_table`` right after ``conceptual_view`` (a page in no toctree warns, and ``-W`` fails)."""
    lines = [line.strip() for line in (ROOT / "docs/architecture/index.rst").read_text().split("\n")]
    assert lines.index("concept_table") == lines.index("conceptual_view") + 1


def test_gate_runs_in_the_gates_workflow() -> None:
    """The ``gates`` workflow step that runs ``test_concept_table.py`` also runs this file."""
    lines = [x for x in (ROOT / ".github/workflows/gates.yml").read_text().split("\n") if "tests/gates/test_concept_table.py" in x]
    assert len(lines) == 1 and "tests/gates/test_concept_table_page.py" in lines[0]


def test_hand_written_table_is_replaced_by_a_pointer() -> None:
    """``conceptual_view.rst`` keeps the heading ``The concept table`` over a pointer, has no ``list-table`` under it, defines ``architecture-concept-table-outside``, and ``layering.rst`` points there.

    Reddens (positive control of the regex: 1 on the base text): the old ``list-table`` under the heading.
    """
    view = (ROOT / "docs/architecture/conceptual_view.rst").read_text()
    assert re.search(r"^The concept table\n-+\n", view, re.M)
    assert len(re.findall(r"^The concept table\n-+\n\n\.\. list-table::", view, re.M)) == 0
    assert ":doc:`concept_table`" in view
    assert re.search(r"^\.\. _architecture-concept-table-outside:", view, re.M)
    assert re.search(r"^\.\. _architecture-conceptual-view:", view, re.M)
    layering = (ROOT / "docs/architecture/layering.rst").read_text()
    assert layering.count("recorded by name in the concept") == 0
    assert layering.count(":ref:`architecture-concept-table-outside`") == 1


# --- --check on a temporary copy of the tree --------------------------------------------------------------------


@pytest.fixture(scope="module")
def tree_copy(tmp_path_factory) -> Path:
    """A copy of what the generator reads: the ``.py`` files of ``orpheus/``, ``tools/``, the table, the debt and the page."""
    root = tmp_path_factory.mktemp("page_tree")
    for top in ("orpheus", "tools"):
        for path in (ROOT / top).rglob("*.py"):
            rel = path.relative_to(ROOT)
            if "__pycache__" in rel.parts:
                continue
            (root / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, root / rel)
    for rel in (TABLE_FILE, DEBT_FILE, PAGE_FILE, "docs/architecture/concept_seed.toml"):
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / rel, root / rel)
    return root


def run_check(root: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run([sys.executable, "-m", MODULE, "--check"], cwd=root, capture_output=True, text=True,
                          env={"PYTHONHASHSEED": "0", "PATH": "/usr/bin:/bin"})


def rewrite(root: Path, rel: str, old: str, new: str):
    path = root / rel
    original = path.read_text()
    if original.count(old) != 1:
        raise Uninstallable(f"{old!r} occurs {original.count(old)} times in {rel}")
    path.write_text(original.replace(old, new))
    return lambda: path.write_text(original)


def test_check_mode_exits_nonzero_on_each_mutation(tree_copy) -> None:
    """``--check`` exits 0 on the copy, and non-zero after (a) a concept field edited, (b) a class dropped from a copied source, (c) a debt entry removed.

    Each edit is made on the temporary copy and restored from its own text; the tracked files are never touched.
    Precondition of each arm: ``rewrite`` found its anchor exactly once. Control: the unedited copy exits 0.
    """
    assert run_check(tree_copy).returncode == 0, run_check(tree_copy).stdout
    arms = {
        "a": (TABLE_FILE, 'definition = "the content digest: one hash of one encoding"',
              'definition = "the content digest: one hash of one encoding (edited)"'),
        "b": ("orpheus/data/cell_xs.py", "class CellXS", "class _GoneXS_"),
        "c": (DEBT_FILE, '[[misplaced]]\nclass = "orpheus/data/cell_xs.py::CellXS"\n\n', ""),
    }
    for name, (rel, old, new) in arms.items():
        restore = rewrite(tree_copy, rel, old, new)
        try:
            done = run_check(tree_copy)
        finally:
            restore()
        assert done.returncode != 0, f"arm ({name}) left --check green"
    assert run_check(tree_copy).returncode == 0
