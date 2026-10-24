r"""The retirement of wave 2 (WF2.4, unit R): the typed-law vocabulary and its three homes are gone, and nothing names them.

The carve (slices S1 to S5) replaced the law classes of ``orpheus/geometry/boundary/`` (16 modules), the realization seam
``orpheus/transport/boundary_realization.py`` and the SN shim ``orpheus/sn/boundary/bound_compat.py``. A retirement that only
deletes loses coverage and leaves orphans (``retirement-audit`` C and G), so this file is the done-when of the deletion commit:
every row below is RED at the base tip (``8c5578c5`` plus the S5 checkpoint) and green exactly when the retirement is complete.
Plan: ``.claude/plans/wave2_boundary_carve.md`` section 7 (7.1 census, 7.2 dependency table, 7.3 the edits).

The claims, each with its kind (``vv-principles``) and the input that reddens it (X1):

* the control (RECORD of the instrument): the census sees each shape it must (a module path in each import spelling, a retired
  name as a definition, an attribute read, an annotation, an import) and reads nothing in a clean source. Green today; reddens
  on a scanner that reads nothing.
* the three retiring paths do not exist (THEOREM of the retirement: a path is present or not). Red input: the tree today holds
  17 files under them.
* nothing under ``orpheus/`` imports a retired module or defines, reads or imports a retiring name (population: every ``.py``
  file under ``orpheus/``, AST; the count is printed). Red input: the retiring files themselves and ``orpheus/sn/boundary/``.
* nothing under ``tests/``, ``examples/``, ``tools/`` and ``student_resources/`` does either; there is no exemption (ruling of
  2026-10-24: the two frozen capture scripts are deleted, so they count as live users until then).
  The first red of this row IS the precondition list of WF2.4: the test files that still import the retired API and are not
  the 7.3 list (they are re-spelled or deleted in the deletion commit, each with its ``retirement-audit`` D.14 demotion check).
* the layer gate's cold-import list names no retired module (a retired module imported in a fresh interpreter is an
  ``ImportError``, so the gate would turn red by the deletion alone).
* ERR-041's index sets retire: ``SNMethodSpace`` carries no ``inflow_indices`` / ``outflow_indices`` (the typed operators read
  their half-spaces off the trace ladder, and ERR-041's invariant is carried there by
  ``test_angular_trace_space.py::TestTheInflowHalfSpaceIsTheIncomingOrdinates``).

**Two populations, said plainly.** A bare-name scan over a name that has a live homonym is red forever, so the names are split.
:data:`RETIRING_NAMES` have no surviving homonym under ``orpheus/`` (a ``def`` census at the base tip decides, ``[M]`` 2026-10-24).
The names that SURVIVE under the same spelling are not scanned: ``BoundaryError``, ``BoundaryResponseNotPositiveError`` and
``SubmarkovViolationError`` (new homes in ``orpheus/transport/boundary_conditions/errors.py``), ``stamp_boundary_role`` (new home
``realization.py``), ``_attenuated_kernel_operator`` (moved into ``realized_operators.py``), ``specular_albedos`` and
``require_vacuum`` (``orpheus/derivations/common/reference_body.py``), ``outer_law`` (a property of the MoC lattice reference),
``boundaries`` (a local of the Peierls slab), ``bc_outer`` (a local of the SN sweep), ``bc_kind`` (a parameter of test helpers),
``face_inventory`` (``orpheus/geometry/faces.py``) and ``face_laws`` (a mapping attribute of the diffusion operators: the
retiring ``face_laws`` was the geometry property). Their OLD homes are inside the retiring paths, which the path row and the
module-import row cover. The exclusion is asserted live (:data:`_HOMONYMS`): a name listed there that no longer occurs in its
survivor is a stale excuse and reddens.

**No exemption, and the two capture scripts.** ``tests/gates/transport/_capture_carry_head.py`` and
``tests/gates/derivations/_capture_reference_carry_head.py`` spelled the retired API by design (frozen records of how
``carry_head.npz`` and ``reference_carry_head.json`` were captured from ``5d04c296``). The ruling of 2026-10-24 deletes them: the
fixtures stay, each reading gate's docstring names the capture commit and the script's git path (recoverable with ``git show``),
and until the deletion the scripts are live users the census reads like any other file.
"""

from __future__ import annotations

import ast
import dataclasses
from pathlib import Path

import pytest

from tests.gates._import_scan import forbidden_imports

_THIS = "tests/gates/test_wave2_retirement.py"
_REPO = Path(__file__).resolve().parents[2]
_ORPHEUS = _REPO / "orpheus"

#: The three retiring homes, as paths (the package is a directory) and as dotted modules.
RETIRED_PATHS: tuple[Path, ...] = (
    _ORPHEUS / "geometry" / "boundary",
    _ORPHEUS / "transport" / "boundary_realization.py",
    _ORPHEUS / "sn" / "boundary" / "bound_compat.py",
)
RETIRED_MODULES: tuple[str, ...] = (
    "orpheus.geometry.boundary",
    "orpheus.transport.boundary_realization",
    "orpheus.sn.boundary.bound_compat",
)

#: The 7.2 symbols with no surviving homonym under ``orpheus/`` (section 7.2, rows R1 to R20; the docstring says which names are not here).
RETIRING_NAMES: frozenset[str] = frozenset(
    {
        # R1, R6: the certification verbs and the law-side realizability
        "assert_specular_pairing_valid", "assert_specular_pairing_involutive", "assert_specular_pairing_measure_preserving",
        "assert_specular_pairing_maps_inflow_to_outflow", "_specular_pairing", "_axis_cosines",
        "assert_realizable", "law_permutes_ordinates", "geometry_map", "response_kernel",
        # R3, R4, R5: the string registries and the tag bridge
        "TAG_REGISTRY", "_law_of", "_wall_of", "specular_albedo", "to_alpha", "BC_REGISTRY",
        "BOUNDARY_OPERATOR_REGISTRY", "_law_from_tag", "resolve_boundary_conditions", "realize_boundary_law",
        # R7, R8, R9: the shim, the recursion, the law algebra
        "_BoundBoundaryOperator", "realize_recursively", "BoundaryRealizer", "SNBoundaryRealizer", "DiffusionBoundaryRealizer",
        "LawSum", "LawScaled", "LawNode", "_admit_as_operand",
        # R10, R11, R12: the response factors, the source, the decks
        "ScalarResponse", "LambertianReemission", "SpecularReemission", "ReemissionClosure", "BoundaryResponseKernel",
        "NoSource", "ConstantInflowSource", "InflowSourceSpec",
        "SelfPairedDeck", "PairedDeck", "BoundaryGeometryMap", "_mirror_motion",
        # R13, R14, R16: the law classes, their parsers and the tag
        "WhiteBoundary", "ReflectiveBoundary", "PeriodicBoundary", "VacuumInflow", "PrescribedInflow", "AlbedoBoundary",
        "SpecularReturn", "IsotropicReturn", "ZeroFluxBoundary", "BoundaryTraceLaw",
        "parse_boundary_law", "parse_boundary_laws", "FaceLaws", "BC",
        # R20: the SN builders that the typed realized operators replaced
        "_narrowed_zero_operator", "_deck_kernel", "_specular_crossed_diagnosis", "_checked_angular_average",
        # the ERR-041 row of 7.3: the guard, its error and the index sets it cross-checked
        "VacuumAppliedToOutgoingTraceError", "inflow_indices", "outflow_indices",
    }
)

#: The names a test file may spell for a reason that is not the retiring symbol (a generic local), removed from the TESTS census only.
_TEST_HOMONYMS: frozenset[str] = frozenset({"geometry_map", "response_kernel"})

#: Names that survive under the same spelling, and the one survivor each (the exclusion is asserted live).
_HOMONYMS: dict[str, str] = {
    "BoundaryError": "orpheus/transport/boundary_conditions/errors.py",
    "BoundaryResponseNotPositiveError": "orpheus/transport/boundary_conditions/errors.py",
    "SubmarkovViolationError": "orpheus/transport/boundary_conditions/errors.py",
    "stamp_boundary_role": "orpheus/transport/boundary_conditions/realization.py",
    "_attenuated_kernel_operator": "orpheus/sn/boundary/realized_operators.py",
    "specular_albedos": "orpheus/derivations/common/reference_body.py",
    "require_vacuum": "orpheus/derivations/common/reference_body.py",
    "outer_law": "orpheus/derivations/continuous/mms/moc.py",
    "boundaries": "orpheus/derivations/continuous/peierls_nystrom/slab.py",
    "bc_outer": "orpheus/sn/loss_representation/__init__.py",
    "face_inventory": "orpheus/geometry/faces.py",
    "face_laws": "orpheus/diffusion/operators.py",
}

_CAPTURE_SCRIPTS: tuple[str, ...] = (
    "tests/gates/transport/_capture_carry_head.py",
    "tests/gates/derivations/_capture_reference_carry_head.py",
)
_CARRY_FIXTURES: tuple[str, ...] = ("tests/gates/transport/carry_head.npz", "tests/gates/derivations/reference_carry_head.json")

_OUTSIDE_TREES: tuple[str, ...] = ("tests", "examples", "tools", "student_resources")


def retired_uses(source: str, names: frozenset[str], package: str | None = None) -> list[str]:
    """What ``source`` spells of the retired vocabulary: an import of a retired module (any spelling, any depth, relative
    imports resolved against ``package``), and a retiring name as a ``def``, a ``class``, a ``Name``, an ``Attribute``, a keyword
    argument or an imported name. A docstring or comment naming one is not a use and is not found (AST, not text)."""
    found = [f"import {name}" for name in forbidden_imports(source, RETIRED_MODULES, package) if name != "<relative>"]
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Name) and node.id in names:
            found.append(f"{node.lineno}: name {node.id}")
        elif isinstance(node, ast.Attribute) and node.attr in names:
            found.append(f"{node.lineno}: attribute {node.attr}")
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and node.name in names:
            found.append(f"{node.lineno}: def {node.name}")
        elif isinstance(node, ast.keyword) and node.arg in names:
            found.append(f"{node.value.lineno}: keyword {node.arg}")
        elif isinstance(node, ast.ImportFrom):
            found += [f"{node.lineno}: import {a.name}" for a in node.names if a.name in names]
    return found


def _package_of(path: Path) -> str:
    """The dotted package a relative import in ``path`` is read from (``orpheus/a/b.py`` is ``orpheus.a``)."""
    return ".".join(path.relative_to(_REPO).parts[:-1])


def _python_files(*roots: Path) -> list[Path]:
    return sorted(p for root in roots for p in root.rglob("*.py") if "__pycache__" not in p.parts and "_build" not in p.parts)


def _findings(files: list[Path], names: frozenset[str]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for path in files:
        rel = path.relative_to(_REPO).as_posix()
        if rel == _THIS:
            continue
        hits = retired_uses(path.read_text(encoding="utf-8"), names, _package_of(path))
        if hits:
            out[rel] = hits
    return out


# ---------------------------------------------------------------------------
# The control: the census sees every shape it must
# ---------------------------------------------------------------------------

_PLANTED: tuple[tuple[str, str], ...] = (
    ("from orpheus.geometry.boundary import X\n", "import orpheus.geometry.boundary"),
    ("import orpheus.transport.boundary_realization\n", "import orpheus.transport.boundary_realization"),
    ("from orpheus.sn.boundary import bound_compat\n", "import orpheus.sn.boundary.bound_compat"),
    ("from orpheus.transport import boundary_realization\n", "import orpheus.transport.boundary_realization"),
    ("def f():\n    from orpheus.geometry.boundary._tag import Y\n", "import orpheus.geometry.boundary._tag"),
    ("from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n    from orpheus.geometry.boundary import BC\n", "import orpheus.geometry.boundary"),
    ("from .boundary import BC\n", "import orpheus.geometry.boundary"),  # relative, read from package orpheus.geometry
    ("from orpheus.geometry import BC\n", "import BC"),
    ("class LawSum: pass\n", "def LawSum"),
    ("x = obj.response_kernel\n", "attribute response_kernel"),
    ("y: BoundaryTraceLaw\n", "name BoundaryTraceLaw"),
    ("f(inflow_indices=1)\n", "keyword inflow_indices"),
)


@pytest.mark.foundation
@pytest.mark.parametrize(("source", "expected"), _PLANTED, ids=[e for _, e in _PLANTED])
def test_the_census_sees_every_shape_it_must(source: str, expected: str) -> None:
    r"""The positive control of every row below (X1): one planted source per shape reads the shape it plants, a docstring or comment
    naming a retired name reads nothing, and a clean source reads nothing. Passes today; reddens on a scanner that reads nothing
    or that misses one import spelling (the relative import is resolved against its package, which the control names)."""
    hits = retired_uses(source, RETIRING_NAMES, package="orpheus.geometry")
    assert any(expected in h for h in hits), (source, hits)


@pytest.mark.foundation
def test_the_census_reads_nothing_in_a_clean_source() -> None:
    clean = '"""BoundaryTraceLaw, LawSum and orpheus.geometry.boundary, named in prose."""\n# BC inflow_indices\nfrom orpheus.transport.boundary_conditions import WallCondition\nx = WallCondition\n'
    assert retired_uses(clean, RETIRING_NAMES, package="orpheus.mesh") == []


@pytest.mark.foundation
def test_the_walk_finds_a_planted_module_and_misses_a_clean_one(tmp_path: Path) -> None:
    r"""The positive control of the FILE walk (X1): a planted module in a throwaway ``orpheus/`` tree that defines a retired class
    and imports a retired module is found by ``_python_files`` plus ``retired_uses`` (both shapes read), a clean sibling is not, and
    a ``__pycache__`` copy is not walked. The census rows read the real tree through the same two functions."""
    package = tmp_path / "orpheus" / "geometry"
    (package / "__pycache__").mkdir(parents=True)
    (package / "planted.py").write_text("from orpheus.geometry.boundary import BC\nclass LawSum: pass\n")
    (package / "clean.py").write_text("from orpheus.transport.boundary_conditions import WallCondition\n")
    (package / "__pycache__" / "stale.py").write_text("class LawSum: pass\n")
    files = _python_files(tmp_path / "orpheus")
    assert [f.name for f in files] == ["clean.py", "planted.py"], files
    hits = {f.name: retired_uses(f.read_text(), RETIRING_NAMES, "orpheus.geometry") for f in files}
    assert hits["clean.py"] == []
    assert any("import orpheus.geometry.boundary" in h for h in hits["planted.py"]) and any("def LawSum" in h for h in hits["planted.py"])


# ---------------------------------------------------------------------------
# The retirement itself
# ---------------------------------------------------------------------------


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_census_reads_nothing_in_a_clean_source")
def test_the_three_retiring_homes_do_not_exist() -> None:
    r"""``orpheus/geometry/boundary/`` (the package, 16 modules), ``orpheus/transport/boundary_realization.py`` and
    ``orpheus/sn/boundary/bound_compat.py`` are absent: a path is present or it is not.

    Population: the three paths (printed in the message with the files still under them). Positive control: the successor home
    ``orpheus/transport/boundary_conditions/`` and the surviving SN module ``orpheus/sn/boundary/realizer.py`` exist (the
    realizer keeps its path), so the check is not reading a moved root. Red input at the base tip: 17 files under the three paths.
    """
    assert (_ORPHEUS / "transport" / "boundary_conditions" / "__init__.py").is_file()
    assert (_ORPHEUS / "sn" / "boundary" / "realizer.py").is_file()
    present = {
        p.relative_to(_REPO).as_posix(): sorted(f.name for f in _python_files(p)) if p.is_dir() else [p.name]
        for p in RETIRED_PATHS
        if p.exists()
    }
    assert not present, f"retiring homes still present ({sum(map(len, present.values()))} files): {present}"


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_census_reads_nothing_in_a_clean_source")
def test_nothing_under_orpheus_names_the_retired_vocabulary() -> None:
    r"""No module under ``orpheus/`` imports a retired module (any spelling) or defines, reads or imports a retiring name
    (:data:`RETIRING_NAMES`, the populations split in the module docstring).

    Population (X2): every ``.py`` file under ``orpheus/`` parsed with ``ast``, the count printed. Positive control: the planted
    sources of ``test_the_census_sees_every_shape_it_must``. Red input at the base tip: the 17 retiring files, and
    ``orpheus/sn/boundary/{errors,realizer}.py`` and ``orpheus/sn/mesh/method_space.py`` (the ERR-041 guard, its error and the two
    index sets). Mutation: re-adding ``from orpheus.geometry.boundary import BC`` under ``TYPE_CHECKING`` anywhere in ``orpheus/``.
    """
    files = _python_files(_ORPHEUS)
    assert len(files) >= 300, f"the census read only {len(files)} files under orpheus/"
    findings = _findings(files, RETIRING_NAMES)
    assert not findings, (
        f"{len(findings)} of {len(files)} files under orpheus/ still name the retired vocabulary: "
        + "; ".join(f"{k} ({len(v)}: {v[0]})" for k, v in sorted(findings.items()))
    )


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_census_reads_nothing_in_a_clean_source")
def test_no_test_example_or_tool_names_the_retired_vocabulary() -> None:
    r"""No file under ``tests/``, ``examples/``, ``tools/`` or ``student_resources/`` imports a retired module or names a retiring
    name; no file is exempt (the two frozen capture scripts are deleted by ruling, so they are live users until then).

    This is the precondition of WF2.4 as a gate: a retired test is deleted and a test of a surviving subject is re-spelled onto
    the successor, so the only way this row is green is that every retiring-symbol user is gone, the two capture scripts included. The red
    message lists them, and a file the plan's 7.3 list does not name is a slice defect, not a deletion.

    Population: every ``.py`` file under the four trees, minus this file; the count is printed. Positive control:
    the planted sources. ``geometry_map`` and ``response_kernel`` are not read in tests (a generic local), see :data:`_TEST_HOMONYMS`.
    """
    files = _python_files(*(_REPO / tree for tree in _OUTSIDE_TREES if (_REPO / tree).exists()))
    assert len(files) >= 500, f"the census read only {len(files)} files"
    findings = _findings(files, RETIRING_NAMES - _TEST_HOMONYMS)
    assert not findings, (
        f"{len(findings)} of {len(files)} files outside orpheus/ still name the retired vocabulary: "
        + "; ".join(f"{k} ({len(v)}: {v[0]})" for k, v in sorted(findings.items()))
    )


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_census_reads_nothing_in_a_clean_source")
def test_the_capture_scripts_are_deleted_and_their_fixtures_stay() -> None:
    r"""The two frozen capture scripts no longer exist (ruling of 2026-10-24: no exemption remains), and the fixtures they were run
    to produce stay (the carry gates read them). Red input at the base tip: both scripts are present. Green when both are deleted
    and both fixtures are files."""
    present = [rel for rel in _CAPTURE_SCRIPTS if (_REPO / rel).exists()]
    assert not present, f"the frozen capture scripts still exist: {present}"
    for rel in _CARRY_FIXTURES:
        assert (_REPO / rel).is_file(), rel


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_census_reads_nothing_in_a_clean_source")
def test_every_excluded_homonym_survives_where_it_is_said_to() -> None:
    r"""The names the census does not scan are each read in their one survivor (:data:`_HOMONYMS`), so the exclusion is a
    measured fact and not an excuse. Reddens on a stale entry (a survivor that no longer holds the name)."""
    for name, rel in sorted(_HOMONYMS.items()):
        tree = ast.parse((_REPO / rel).read_text(encoding="utf-8"))
        spelled = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)}
        spelled |= {n.attr for n in ast.walk(tree) if isinstance(n, ast.Attribute)}
        spelled |= {n.name for n in ast.walk(tree) if isinstance(n, (ast.FunctionDef, ast.ClassDef))}
        spelled |= {a.name for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) for a in n.names}
        assert name in spelled, f"{name} is listed as surviving in {rel} and is not spelled there"


@pytest.mark.foundation
def test_the_cold_import_list_of_the_layer_gate_names_no_retired_module() -> None:
    r"""``test_entry_point_imports_in_a_fresh_interpreter`` imports each listed entry in a cold interpreter; a retired module in
    the list is an ``ImportError`` the day it is deleted. The list is read as the string constants of that test's decorators (AST).

    Positive control: the same list holds the live entry ``orpheus.mesh`` and at least 30 entries. Red input at the base tip: the
    three entries ``orpheus.geometry.boundary``, ``orpheus.transport.boundary_realization`` and ``orpheus.sn.boundary.bound_compat``.
    """
    tree = ast.parse((_REPO / "tests" / "gates" / "test_layer_imports.py").read_text(encoding="utf-8"))
    fn = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "test_entry_point_imports_in_a_fresh_interpreter")
    entries = [c.value for d in fn.decorator_list for c in ast.walk(d) if isinstance(c, ast.Constant) and isinstance(c.value, str)]
    assert len(entries) >= 30 and "orpheus.mesh" in entries, entries
    retired = [e for e in entries if e in RETIRED_MODULES]
    assert not retired, f"the cold-import list names retired modules: {retired}"


@pytest.mark.foundation
def test_the_sn_method_space_carries_no_index_sets() -> None:
    r"""ERR-041's index sets are retired: ``SNMethodSpace`` has no ``inflow_indices`` and no ``outflow_indices`` field. The typed
    operators read their half-spaces off the trace ladder (``trace.inflow_space(face)``), so the pair had 0 production readers
    outside the zero map's own cross-check, and the invariant that guard checked is carried by the ladder's own witness
    (``tests/gates/numerics/test_angular_trace_space.py::TestTheInflowHalfSpaceIsTheIncomingOrdinates``).

    Read from the dataclass fields (a name grep would read a homonym local). Positive control: the same class still carries
    ``quadrature``, ``face`` and ``trace``. Red input at the base tip: both fields exist.
    """
    from orpheus.sn.mesh.method_space import SNMethodSpace

    names = {f.name for f in dataclasses.fields(SNMethodSpace)}
    assert {"quadrature", "face", "trace"} <= names, sorted(names)
    assert not {"inflow_indices", "outflow_indices"} & names, sorted(names)
