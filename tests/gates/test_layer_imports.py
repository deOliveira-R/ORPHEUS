"""Import-linter: enforce the L0 / L1 / (input) / L2 / L3 layer contract.

Layers (per plan §P3.0):
  (values) (values/)   — what a value is and when two values are the same
                         (layer_dag.md ruling 6); imports nothing from orpheus.
  L0  (derivations/)   — Branch-1 references (SymPy / mpmath);
                         beside L1 and the input layer, below L2: it
                         may import numerics, geometry, mesh and data.
  L1  (numerics/)      — math primitives; knows no neutrons.
  (input) (geometry/, data/) — geometry + nuclear data.
  (mesh)  (mesh/)      — the discretisation overlay on the geometry:
                         imports geometry and numerics; geometry,
                         data and numerics never import it.
  L2  (transport/)     — transport vocabulary; method-agnostic.
  L3  (sn/, pn/, moc/, cp/, mc/, diffusion/, kinetics/,
       fuel/, thermal_hydraulics/, homogeneous/) — one method's
       machinery; method-specific.
  L4  (plotting.py)    — orchestration; consumes everything.

Imports flow only L3 -> L2 -> L1 (and L3 -> L0 references, and any
layer above the input layer -> input). L1 imports no input package except
geometry, an edge pending layer_dag.md P1; numerics -> data is forbidden
(layer_dag.md P2). Forbidden edges raise the parametrised test.

Tolerances:
  - TYPE_CHECKING imports are edges like any other (retired at layer-DAG
    P3: 7 of 7 tolerated imports were transport -> sn/diffusion, 3 of them
    dead, 0 from numerics; a typing-only import of a method type from L2
    was a concept in the wrong layer).
  - WHITELIST entries cover transitional exemptions retired by
    later Phase 3 steps; each whitelist entry carries a comment
    naming its retirement trigger.
"""
from __future__ import annotations

import ast
import pathlib
import subprocess
import sys
from collections.abc import Iterable

import pytest

ORPHEUS_ROOT = pathlib.Path(__file__).resolve().parents[2] / "orpheus"

# Layer assignment.
# Layer DAG ruling 6: what a value is (the admitted reals) and when two values are
# the same (content identity); the lowest package, importing nothing from orpheus.
VALUES_PACKAGES: frozenset[str] = frozenset({"values"})
L0_PACKAGES: frozenset[str] = frozenset({"derivations"})
INPUT_PACKAGES: frozenset[str] = frozenset({"geometry", "data"})
MESH_PACKAGES: frozenset[str] = frozenset({"mesh"})
# #405 P1 step 8: the specification composes materials, geometry and a question;
# above data, geometry and numerics; a sibling of mesh (neither imports the other); below transport, L3 and derivations.
SPECIFICATION_PACKAGES: frozenset[str] = frozenset({"specification"})
# #405 P2 step 1: the reference package (readings, published and reference solutions,
# the certificates); above the specification, below transport, L3 and derivations,
# which writes reference solutions and so may import it.
REFERENCE_PACKAGES: frozenset[str] = frozenset({"reference"})
# P2 (layer_dag.md ruling 2): data sits strictly above numerics; numerics names nothing in data.
DATA_PACKAGES: frozenset[str] = frozenset({"data"})
L2_PACKAGES: frozenset[str] = frozenset({"transport"})
L3_PACKAGES: frozenset[str] = frozenset(
    {
        "sn",
        "pn",
        "moc",
        "cp",
        "mc",
        "diffusion",
        "kinetics",
        "fuel",
        "thermal_hydraulics",
        "homogeneous",
    }
)

#: Every top-level name under ``orpheus/``: its packages and its modules (``plotting``).
ORPHEUS_TOP_LEVEL: frozenset[str] = frozenset(
    path.stem if path.suffix == ".py" else path.name
    for path in ORPHEUS_ROOT.iterdir()
    if (path.is_dir() and (path / "__init__.py").is_file()) or (path.suffix == ".py" and path.name != "__init__.py")
)

FORBIDDEN_EDGES: dict[str, frozenset[str]] = {
    "values": ORPHEUS_TOP_LEVEL - VALUES_PACKAGES,
    "numerics": DATA_PACKAGES | MESH_PACKAGES | SPECIFICATION_PACKAGES | REFERENCE_PACKAGES | L2_PACKAGES | L3_PACKAGES,
    "geometry": MESH_PACKAGES | SPECIFICATION_PACKAGES | REFERENCE_PACKAGES | L2_PACKAGES | L3_PACKAGES,
    "data": MESH_PACKAGES | SPECIFICATION_PACKAGES | REFERENCE_PACKAGES | L2_PACKAGES | L3_PACKAGES,
    "mesh": SPECIFICATION_PACKAGES | REFERENCE_PACKAGES | L2_PACKAGES | L3_PACKAGES,
    "specification": MESH_PACKAGES | REFERENCE_PACKAGES | L2_PACKAGES | L3_PACKAGES | L0_PACKAGES,
    "reference": MESH_PACKAGES | L2_PACKAGES | L3_PACKAGES | L0_PACKAGES,
    "transport": L3_PACKAGES,
    "sn": L3_PACKAGES - {"sn"},
    "pn": L3_PACKAGES - {"pn"},
    "moc": L3_PACKAGES - {"moc"},
    "cp": L3_PACKAGES - {"cp"},
    "mc": L3_PACKAGES - {"mc"},
    "diffusion": L3_PACKAGES - {"diffusion"},
    "kinetics": L3_PACKAGES - {"kinetics"},
    "fuel": L3_PACKAGES - {"fuel"},
    "thermal_hydraulics": L3_PACKAGES - {"thermal_hydraulics"},
    "homogeneous": L3_PACKAGES - {"homogeneous"},
    "derivations": L2_PACKAGES | L3_PACKAGES,
}

WHITELIST: frozenset[tuple[str, str]] = frozenset(
    {
        # RETIRE_IN_P3_FOLLOWUP — MMS source uses MOCMesh / MOCQuadrature
        # at module level; move to test side or import only L2 primitives.
        ("derivations/continuous/mms/moc.py", "moc"),
        # RETIRE_IN_P3_FOLLOWUP — the non-vacuum MMS reference lazily builds
        # its prescribed-inflow source from transport vocabulary
        # (AngularBoundarySourceSink / AngularSourceSink / TimedFullField) inside
        # function bodies; move to the test side or import only L2 primitives.
        ("derivations/continuous/mms/sn.py", "transport"),
    }
)


def _iter_python_modules(root: pathlib.Path) -> Iterable[pathlib.Path]:
    for p in root.rglob("*.py"):
        if "__pycache__" in p.parts:
            continue
        yield p


def _top_level_package(rel_path: pathlib.Path) -> str:
    return rel_path.parts[0]


def _absolute_module(rel_path: pathlib.PurePosixPath, node: ast.ImportFrom) -> str | None:
    """The absolute dotted name an ``ImportFrom`` reads, relative imports resolved.

    ``from ..transport import x`` in ``orpheus/mesh/axis.py`` is
    ``orpheus.transport``: a relative import of level ``n`` climbs ``n - 1``
    packages above the one holding the module. Without this the linter reads
    the bare ``transport`` and drops it at the ``orpheus.`` prefix test.
    """
    if node.level == 0:
        return node.module
    package = ("orpheus", *rel_path.parts[:-1])
    if node.level - 1 > len(package) - 1:
        return None
    base = package[: len(package) - (node.level - 1)]
    return ".".join((*base, node.module) if node.module else base)


def _is_orpheus_module(dotted: str) -> bool:
    """Whether ``dotted`` names a module or package file under ``orpheus/``."""
    parts = dotted.split(".")
    if parts[0] != "orpheus" or len(parts) < 2:
        return False
    path = ORPHEUS_ROOT.joinpath(*parts[1:])
    return path.with_suffix(".py").is_file() or (path / "__init__.py").is_file()


def _imports_of_source(
    rel_path: pathlib.PurePosixPath, src: str,
) -> list[tuple[str, bool]]:
    """Parse imports, marking TYPE_CHECKING-guarded ones."""
    tree = ast.parse(src, filename=str(rel_path))
    results: list[tuple[str, bool]] = []

    def _visit(node: ast.AST, in_tc: bool) -> None:
        if (
            isinstance(node, ast.If)
            and isinstance(node.test, ast.Name)
            and node.test.id == "TYPE_CHECKING"
        ):
            for child in node.body:
                _visit(child, True)
            for child in node.orelse:
                _visit(child, in_tc)
            return
        if isinstance(node, ast.ImportFrom):
            module = _absolute_module(rel_path, node)
            if module:
                # A name that is itself a module (`from orpheus import sn`,
                # `from orpheus.numerics import operator`) is an edge to THAT
                # module; recording only `node.module` would let the package
                # root stand in for it and hide the edge.
                submodules = [f"{module}.{a.name}" for a in node.names if _is_orpheus_module(f"{module}.{a.name}")]
                other_names = len(submodules) < len(node.names)
                if other_names or not submodules:
                    results.append((module, in_tc))
                results.extend((name, in_tc) for name in submodules)
        elif isinstance(node, ast.Import):
            for alias in node.names:
                results.append((alias.name, in_tc))
        for child in ast.iter_child_nodes(node):
            _visit(child, in_tc)

    _visit(tree, False)
    return results


def _check_source(
    rel: pathlib.PurePosixPath, src: str, whitelist: frozenset[tuple[str, str]] | None = None,
) -> list[str]:
    """The layer violations of one module's source, ``rel`` relative to ``orpheus/``.

    ``whitelist`` is the exemption set read; the vacuous-row guard passes the
    set minus one row to ask whether that row still names a live edge. The
    default ``None`` reads the module's ``WHITELIST`` at CALL time, so a probe
    that rebinds ``WHITELIST`` reaches every caller.
    """
    if whitelist is None:
        whitelist = WHITELIST
    src_pkg = _top_level_package(rel)
    if src_pkg not in FORBIDDEN_EDGES:
        return []
    forbidden = FORBIDDEN_EDGES[src_pkg]
    rel_str = rel.as_posix()
    violations: list[str] = []
    for module_name, _is_tc in _imports_of_source(rel, src):
        if not module_name.startswith("orpheus."):
            continue
        tgt_pkg = module_name.split(".")[1]
        if tgt_pkg not in forbidden:
            continue
        if (rel_str, tgt_pkg) in whitelist:
            continue
        violations.append(
            f"{rel} imports {module_name} "
            f"(forbidden: {src_pkg} → {tgt_pkg})"
        )
    return violations


def _check_module(module_path: pathlib.Path) -> list[str]:
    rel = pathlib.PurePosixPath(module_path.relative_to(ORPHEUS_ROOT).as_posix())
    return _check_source(rel, module_path.read_text())


_ALL_MODULES = sorted(_iter_python_modules(ORPHEUS_ROOT))


@pytest.mark.foundation
@pytest.mark.parametrize("module_path", _ALL_MODULES, ids=str)
def test_no_forbidden_imports(module_path: pathlib.Path) -> None:
    violations = _check_module(module_path)
    assert not violations, "\n".join(violations)


@pytest.mark.foundation
@pytest.mark.parametrize("target", ["orpheus.numerics.space", "orpheus.geometry.chart", "orpheus.plotting"])
def test_values_imports_nothing_from_orpheus(target: str) -> None:
    """Ruling 6: ``values`` is the bottom of the package graph. A values module importing any other
    orpheus package or module is a violation; the row reads it on a synthetic source, so the gate's
    red exists in the tree (X1), and ``ORPHEUS_TOP_LEVEL`` is checked to hold the names it must."""
    assert {"numerics", "geometry", "plotting", "values"} <= ORPHEUS_TOP_LEVEL, sorted(ORPHEUS_TOP_LEVEL)
    violations = _check_source(pathlib.PurePosixPath("values/probe.py"), f"import {target}\n")
    assert violations and "forbidden: values" in violations[0], violations


# ---------------------------------------------------------------------------
# The points and the group action live in geometry (layer_dag.md rulings 1, 8)
# ---------------------------------------------------------------------------
#
# The point groups (`symmetry`), the exact realization of C_n they generate
# (`roots_of_unity`) and the point sets they act on (the `Manifold` closed sum)
# are geometry: numerics (measures, quadratures, bases) imports them, never the
# reverse. Until WP-B (2026-10-22) all three were numerics modules, which put
# geometry -> numerics edges under the rigid-motion core.

# `manifold` is in the set although the numerics half kept a module of its own
# (`index_set`): a numerics module of that stem could only be a re-export shim.
_GEOMETRY_OWNED_NUMERICS_STEMS: frozenset[str] = frozenset({"symmetry", "roots_of_unity", "manifold"})


def _modules_defining_class(root: pathlib.Path, name: str) -> list[str]:
    """Every module under ``root`` whose top level defines ``class <name>`` (AST)."""
    return sorted(
        str(path.relative_to(ORPHEUS_ROOT))
        for path in _iter_python_modules(root)
        if any(
            isinstance(node, ast.ClassDef) and node.name == name
            for node in ast.parse(path.read_text(encoding="utf-8")).body
        )
    )


@pytest.mark.foundation
def test_the_points_and_the_group_action_are_geometry() -> None:
    """No numerics module is the point-group, root-of-unity or manifold module (a
    re-export shim of the old path included), none defines the
    ``Manifold`` ABC, and geometry holds the point groups. The class finder is checked on
    its known member first (X1): geometry's own ``RigidMotion``."""
    assert _modules_defining_class(ORPHEUS_ROOT / "geometry", "RigidMotion") == [
        "geometry/transformation.py"
    ]
    misplaced = sorted(
        str(p.relative_to(ORPHEUS_ROOT))
        for p in _iter_python_modules(ORPHEUS_ROOT / "numerics")
        if p.stem in _GEOMETRY_OWNED_NUMERICS_STEMS
    )
    assert not misplaced, f"geometry concepts still in numerics: {misplaced}"
    assert _modules_defining_class(ORPHEUS_ROOT / "numerics", "Manifold") == []
    assert _modules_defining_class(ORPHEUS_ROOT / "geometry", "Manifold") == ["geometry/manifold.py"]
    assert (ORPHEUS_ROOT / "geometry" / "symmetry.py").is_file()
    assert (ORPHEUS_ROOT / "geometry" / "roots_of_unity.py").is_file()


_INDEX_SET = "orpheus.numerics.index_set"


def _geometry_imports_of_the_index_sets(rel: pathlib.PurePosixPath, src: str) -> list[str]:
    """The targets of ``src`` naming the finite index sets' module, every import kind."""
    targets = _imported_targets(rel, src) + _dynamic_imports(src)
    return [t for t in targets if t == _INDEX_SET or t.startswith(_INDEX_SET + ".")]


@pytest.mark.foundation
def test_geometry_never_imports_the_finite_index_sets() -> None:
    """Ruling 8: the ``Manifold`` closed sum's geometric members name no numerics member.

    ``IndexSet`` and ``EnergyGroups`` subclass geometry's ``Manifold`` from
    :mod:`orpheus.numerics.index_set`, the allowed direction; each member answers
    ``ambient_dim`` itself, so no geometry module needs the index sets at any import
    kind, a string literal handed to ``importlib.import_module`` included. The predicate
    is checked first on the five shapes it must refuse (X1)."""
    probe = pathlib.PurePosixPath("geometry/probe.py")
    for src in (
        "from orpheus.numerics.index_set import IndexSet\n",
        "def f():\n    from orpheus.numerics.index_set import ENERGY\n",
        "from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n    import orpheus.numerics.index_set\n",
        "from orpheus.numerics import index_set\n",
        "import importlib\ndef f():\n    return importlib.import_module('orpheus.numerics.index_set')\n",
    ):
        assert _geometry_imports_of_the_index_sets(probe, src), src
    offenders = [
        f"{path.relative_to(ORPHEUS_ROOT)}: {target}"
        for path in _iter_python_modules(ORPHEUS_ROOT / "geometry")
        for target in _geometry_imports_of_the_index_sets(
            pathlib.PurePosixPath(path.relative_to(ORPHEUS_ROOT).as_posix()),
            path.read_text(encoding="utf-8"),
        )
    ]
    assert not offenders, offenders


# ---------------------------------------------------------------------------
# The numerics -> geometry back-edge, and the discipline that makes it safe
# ---------------------------------------------------------------------------
#
# `numerics/__init__.py` imports `measure`, whose supports are geometry's
# (`geometry.manifold`, with the point groups of `geometry.symmetry`; until
# layer DAG WP-B of 2026-10-22 `numerics/__init__` imported `symmetry` itself,
# a numerics module importing `geometry.transformation`). That is a genuine
# package CYCLE: importing `orpheus.numerics` runs `orpheus.geometry.__init__`, which
# imports `orpheus.geometry.boundary`, which imports back into
# `orpheus.numerics` — while `orpheus.numerics.__init__` is still
# mid-execution.
#
# It resolves, and the reason is precise: a partially-initialised package can
# serve `from orpheus.numerics.measure import X` (a SUBMODULE import, resolved
# through the module system) but NOT `from orpheus.numerics import X` (an
# ATTRIBUTE lookup on a package whose body has not reached that line yet).
# Every one of geometry's numerics imports happens to be the first form.
#
# "Happens to be" is not a guarantee, and the failure is a hard ImportError at
# interpreter start-up, not a subtle wrong answer. These two gates make the
# discipline explicit: one structural, one end-to-end.

_INPUT_PACKAGE_ROOTS = sorted(INPUT_PACKAGES | MESH_PACKAGES | SPECIFICATION_PACKAGES)


@pytest.mark.foundation
@pytest.mark.parametrize("package", _INPUT_PACKAGE_ROOTS)
def test_input_layer_imports_numerics_only_by_submodule(package: str) -> None:
    """An input-layer package may not import the `numerics` PACKAGE itself.

    `from orpheus.numerics.measure import DiscreteMeasure` is fine;
    `from orpheus.numerics import DiscreteMeasure` is not, because the second
    form is an attribute lookup that fails while `numerics/__init__` is still
    executing — which is exactly the state it is in when `measure` reaches
    across to `geometry.manifold`.
    """
    offenders: list[str] = []
    for module_path in _iter_python_modules(ORPHEUS_ROOT / package):
        source = module_path.read_text(encoding="utf-8")
        tree = ast.parse(source, filename=str(module_path))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                # `level > 0` is a relative import — never cross-package.
                if node.level == 0 and node.module == "orpheus.numerics":
                    names = ", ".join(a.name for a in node.names)
                    offenders.append(
                        f"{module_path.relative_to(ORPHEUS_ROOT)}:{node.lineno} "
                        f"`from orpheus.numerics import {names}` — use the "
                        f"submodule (`from orpheus.numerics.<mod> import ...`)"
                    )
            elif isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name == "orpheus.numerics":
                        offenders.append(
                            f"{module_path.relative_to(ORPHEUS_ROOT)}:"
                            f"{node.lineno} `import orpheus.numerics`"
                        )
    assert not offenders, (
        "input-layer -> numerics imports must be SUBMODULE-level; these are "
        "package-level and will deadlock the numerics->geometry back-edge:\n"
        + "\n".join(offenders)
    )


@pytest.mark.foundation
@pytest.mark.parametrize(
    "entry",
    [
        "orpheus",
        "orpheus.numerics",
        "orpheus.geometry.symmetry",
        # layer DAG WP-B (2026-10-22): the points and the group action moved
        # down into geometry; the face naming left face_layout for geometry;
        # the finite index sets stayed in numerics, in their own module.
        "orpheus.geometry.faces",
        "orpheus.geometry.roots_of_unity",
        "orpheus.numerics.index_set",
        "orpheus.numerics.face_layout",
        # R2 of #434 (2026-09-03): the four entries the V1 cycle killed that
        # the list above could not see — `[M]` with `manifold -> symmetry` at
        # module scope AND symmetry still importing the axis table from
        # manifold, 6 of 9 entry points died while `import orpheus` stayed
        # green (order-dependent; plan-authoring §6d).
        "orpheus.geometry.manifold",
        "orpheus.numerics.measure",
        # (``orpheus.numerics.invariance`` was an entry until S5, 2026-10-19,
        # which deleted the module: invariance is the measure's composition.)
        "orpheus.numerics.quadrature.registry",
        "orpheus.geometry",
        "orpheus.geometry.transformation",
        # P1 step 1 of #405 (2026-09-25): the mesh became its own package,
        # importing geometry (whose boundary package now holds `BC`) and
        # numerics; each moved module, and the boundary package, cold.
        "orpheus.geometry.boundary",
        # Layer DAG WP-E1 (ruling 7): the realization seam left geometry for
        # transport, and the SN shim for sn; each new home, cold.
        "orpheus.transport.boundary_realization",
        "orpheus.sn.boundary.bound_compat",
        "orpheus.mesh",
        "orpheus.mesh.structured",
        "orpheus.mesh.factories",
        "orpheus.mesh.axis",
        # the layer-DAG campaign: the legacy-mesh adapter left the axis
        # module for its own, above axis and structured.
        "orpheus.mesh.legacy_adapter",
        "orpheus.sn.solver",
        # step 2 of the consumers campaign (2026-09-13): the Strategy value
        # module — solver -> splitting -> operators/loss_representation, no
        # runtime edge back to the hub (plan-authoring §6d).
        "orpheus.sn.splitting",
        # step 3 U1 (2026-09-17): the Solution-tier types — the gauges (the
        # section that picked the representative) and the kind-typed outcomes
        # + the exit report; numerics-tier, importing posing/pencil/
        # operator only (plan-authoring §6d; consumers_step3_design.md §6.3).
        "orpheus.numerics.gauge",
        "orpheus.numerics.outcome",
        # step 2 C3b (2026-09-13): the terminal-object types — numerics-tier,
        # importing only the operator algebra (plan §6.3: no new package edge).
        "orpheus.numerics.pencil",
        "orpheus.numerics.posing",
        # #405 P1 steps 6 and 7 (2026-10-02): the mesh-free functions (S6.21),
        # the question values (S7.13) and the one real-number parser (#559).
        "orpheus.numerics.mesh_free_function",
        "orpheus.numerics.question",
        "orpheus.values",
        "orpheus.values.content",
        "orpheus.values.scalars",
        # #405 P1 step 8 (S8.3): the specification and its keys.
        "orpheus.specification",
        "orpheus.data.cells",
        "orpheus.geometry.extent",
        # #405 P2 step 1: the enclosure and the reference package.
        "orpheus.numerics.enclosure",
        # #405 P2 step 3: the observables.
        "orpheus.numerics.observable",
        "orpheus.reference",
    ],
)
def test_entry_point_imports_in_a_fresh_interpreter(entry: str) -> None:
    """Each entry point imports cleanly from a COLD interpreter.

    The end-to-end companion to the structural gate above, and the one that
    cannot be fooled: a cycle that resolves only because some other module
    happened to be imported first will fail here. It must be a subprocess —
    in-process the modules are already in ``sys.modules`` and every import is
    a no-op.
    """
    result = subprocess.run(
        [sys.executable, "-c", f"import {entry}"],
        capture_output=True,
        text=True,
        cwd=ORPHEUS_ROOT.parent,
    )
    assert result.returncode == 0, (
        f"`import {entry}` failed in a fresh interpreter:\n{result.stderr}"
    )


@pytest.mark.foundation
def test_values_loads_nothing_else_of_orpheus_in_a_fresh_interpreter() -> None:
    """Ruling 6 at run time: a cold ``import orpheus.values.content`` loads only ``orpheus`` and
    ``orpheus.values``. The AST row above cannot see an ``importlib`` call or an import a module's
    code performs when it runs; this one reads ``sys.modules`` after the import, so it can. The
    loaded set must contain the module imported (the instrument's control), and nothing outside
    the package (``[M]`` 2026-10-22: 4 modules; 64 when the encoder lived in numerics)."""
    probe = "import sys, orpheus.values.content; print(*sorted(m for m in sys.modules if m.split('.')[0] == 'orpheus'))"
    result = subprocess.run([sys.executable, "-c", probe], capture_output=True, text=True, cwd=ORPHEUS_ROOT.parent)
    assert result.returncode == 0, result.stderr
    loaded = set(result.stdout.split())
    assert "orpheus.values.content" in loaded, loaded
    outside = sorted(m for m in loaded if m != "orpheus" and m.split(".")[1] != "values")
    assert not outside, f"importing orpheus.values.content loaded modules outside orpheus.values: {outside}"


# ---------------------------------------------------------------------------
# The linter's own law, on synthetic sources
# ---------------------------------------------------------------------------
#
# Every parametrised row above reads a real module, so a row the tree never
# exercises has no witness: the `geometry -> mesh` row, for instance, is only
# ever green, and a linter that lost it would stay green. These legs feed the
# pure `_check_source` one import each, and each forbidden leg must name the
# edge it forbids.

_FORBIDDEN_LEGS = [
    ("mesh/x.py", "from orpheus.transport.fields import X\n", "mesh → transport"),
    ("geometry/__init__.py", "from orpheus.mesh.structured import Mesh1D\n", "geometry → mesh"),
    ("data/x.py", "import orpheus.mesh\n", "data → mesh"),
    ("numerics/x.py", "from orpheus.mesh import Mesh1D\n", "numerics → mesh"),
    # A relative import climbing out of the package: `from ..transport`
    # in `orpheus/mesh/x.py` reads `orpheus.transport`.
    ("mesh/x.py", "from ..transport import fields\n", "mesh → transport"),
    # A package named as an imported name off the root: `from orpheus import sn`.
    ("transport/x.py", "from orpheus import sn\n", "transport → sn"),
    # #405 P1 step 8 (S8.3): the tempting spelling of the group-count rule.
    ("specification/x.py", "from orpheus.transport.mesh import MaterialMesh\n", "specification → transport"),
    ("specification/x.py", "from orpheus.mesh import Mesher\n", "specification → mesh"),
    ("specification/x.py", "from orpheus.derivations import get\n", "specification → derivations"),
    ("specification/x.py", "from orpheus.sn.problem import SNProblem\n", "specification → sn"),
    ("data/x.py", "from orpheus.specification import Specification\n", "data → specification"),
    ("geometry/x.py", "from orpheus.specification import Specification\n", "geometry → specification"),
    ("numerics/x.py", "from orpheus.specification import Specification\n", "numerics → specification"),
    ("mesh/x.py", "from orpheus.specification import Specification\n", "mesh → specification"),
    # A lazy L0 -> L3 import: a function-body import is an edge like any other.
    ("derivations/x.py", "def f():\n    from orpheus.cp.solver import CPParams\n", "derivations → cp"),
    # The retired WHITELIST row: re-adding it turns this red.
    (
        "derivations/continuous/sood_registry/builders.py",
        "def f():\n    from orpheus.cp.solver import CPParams\n",
        "derivations → cp",
    ),
    # layer_dag.md P2: numerics names nothing in data, guarded or not, and L1
    # has no TYPE_CHECKING tolerance at all.
    ("numerics/x.py", "from orpheus.data.energy_grid import EnergyGrid\n", "numerics → data"),
    (
        "numerics/x.py",
        "from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n    from orpheus.data.macro_xs.mixture import Mixture\n",
        "numerics → data",
    ),
    (
        "numerics/x.py",
        "from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n    from orpheus.sn.problem import SNProblem\n",
        "numerics → sn",
    ),
    # Layer DAG P3: a TYPE_CHECKING import is an edge like any other.
    ("transport/x.py", "from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n    from orpheus.sn.problem import SNProblem\n", "transport → sn"),
    ("numerics/x.py", "from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n    from orpheus.transport.mesh.material_mesh import MaterialMesh\n", "numerics → transport"),
    # Layer DAG WP-E1 (ruling 7): the realization seam is transport's, and a
    # geometry module that reached for it again would re-weld the declaration
    # to its realization.
    ("geometry/boundary/x.py", "from orpheus.transport.boundary_realization import BoundaryRealizer\n", "geometry → transport"),
    ("geometry/boundary/x.py", "from orpheus.sn.boundary.bound_compat import _BoundBoundaryOperator\n", "geometry → sn"),
]

_ADMITTED_LEGS = [
    ("mesh/x.py", "from orpheus.geometry.boundary import BC\n"),
    ("mesh/x.py", "from orpheus.numerics.measure import DiscreteMeasure\n"),
    ("mesh/x.py", "from .structured import Mesh1D\n"),
    ("transport/x.py", "from orpheus.mesh import Mesh1D\n"),
    ("derivations/x.py", "from orpheus.mesh import Mesh1D\n"),
    ("specification/x.py", "from orpheus.data.cells import CellCoefficient\n"),
    ("specification/x.py", "from orpheus.geometry.extent import GeometryExtent\n"),
    ("specification/x.py", "from orpheus.numerics.question import Eigen\n"),
    ("derivations/x.py", "from orpheus.specification import Specification\n"),
    ("transport/x.py", "from orpheus.specification import Specification\n"),
    ("sn/x.py", "from orpheus.specification import Specification\n"),
    # layer_dag.md P2: the energy arm's view and rule live in data, over the numerics axis.
    ("data/x.py", "from orpheus.numerics.axis import EnergyAxis\n"),
    ("sn/x.py", "from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n    from orpheus.transport.mesh.material_mesh import MaterialMesh\n"),
    # Layer DAG WP-E1: the realization reads the declarations, and both
    # method realizers read the realization.
    ("transport/x.py", "from orpheus.geometry.boundary import BoundaryTraceLaw, LawScaled, LawSum\n"),
    ("sn/x.py", "from orpheus.transport.boundary_realization import stamp_boundary_role\n"),
    ("diffusion/x.py", "from orpheus.transport.boundary_realization import stamp_boundary_role\n"),
]


@pytest.mark.foundation
@pytest.mark.parametrize(("rel", "source", "edge"), _FORBIDDEN_LEGS)
def test_the_linter_refuses_each_forbidden_edge(rel: str, source: str, edge: str) -> None:
    violations = _check_source(pathlib.PurePosixPath(rel), source)
    assert len(violations) == 1 and f"(forbidden: {edge})" in violations[0], violations


@pytest.mark.foundation
@pytest.mark.parametrize(("rel", "source"), _ADMITTED_LEGS)
def test_the_linter_admits_each_allowed_edge(rel: str, source: str) -> None:
    assert _check_source(pathlib.PurePosixPath(rel), source) == []


# ---------------------------------------------------------------------------
# Geometry declares; it never builds an operator
# ---------------------------------------------------------------------------
#
# Layer DAG ruling 1 (the user, 2026-10-22): geometry is declarative. A
# geometry module that imports the operator algebra is a module that BUILDS
# or names operators, which is realization, and realization is transport's
# (ruling 7). The package-level gate cannot see this: geometry -> numerics is
# an admitted edge, and `numerics.operator` is one numerics submodule among
# many. Until WP-E1 two geometry modules did it, `boundary/_realizer.py` (the
# realization seam) and `boundary/_bound_compat.py` (the SN shim, a
# `LinearOperator` subclass); both moved out, and this gate is their
# retirement's witness. A TYPE_CHECKING import is refused too: the last one,
# `BoundaryTraceLaw.realize -> LinearOperator`, was a hook with 0 overrides
# and 0 callers that only raised, and it retired with WP-E1's review (a
# declaration does not know what it is realized to; the realizer does).

_OPERATOR_ALGEBRA = "orpheus.numerics.operator"


def _operator_algebra_imports(
    rel: pathlib.PurePosixPath, src: str,
) -> list[str]:
    """The imports of the operator algebra in one module's source, at every
    import kind (top level, function body, ``TYPE_CHECKING``)."""
    return [
        module_name
        for module_name, _is_tc in _imports_of_source(rel, src)
        if module_name == _OPERATOR_ALGEBRA or module_name.startswith(_OPERATOR_ALGEBRA + ".")
    ]


@pytest.mark.foundation
def test_no_geometry_module_imports_the_operator_algebra() -> None:
    geometry_root = ORPHEUS_ROOT / "geometry"
    modules = sorted(_iter_python_modules(geometry_root))
    # The denominator is stated, so an empty walk cannot read as a clean tree.
    assert len(modules) > 0, f"no modules under {geometry_root}"
    offenders = [
        f"{module_path.relative_to(ORPHEUS_ROOT)} imports {target}"
        for module_path in modules
        for target in _operator_algebra_imports(
            pathlib.PurePosixPath(module_path.relative_to(ORPHEUS_ROOT).as_posix()),
            module_path.read_text(encoding="utf-8"),
        )
    ]
    assert not offenders, (
        "geometry is declarative (layer DAG rulings 1 and 7): no geometry "
        "module imports orpheus.numerics.operator, at any import kind; a realization "
        "belongs in orpheus.transport.boundary_realization or a method "
        "package:\n" + "\n".join(offenders)
    )


_OPERATOR_ALGEBRA_REFUSED_LEGS = [
    ("geometry/boundary/x.py", "from orpheus.numerics.operator import BlockRole\n"),
    ("geometry/boundary/x.py", "def f():\n    from orpheus.numerics.operator import LinearOperator\n"),
    ("geometry/boundary/x.py", "import orpheus.numerics.operator\n"),
    # The module named as an imported name, not as the `from` path.
    ("geometry/boundary/x.py", "from orpheus.numerics import operator\n"),
    (
        "geometry/boundary/x.py",
        "from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n    from orpheus.numerics.operator import LinearOperator\n",
    ),
]

_OPERATOR_ALGEBRA_ADMITTED_LEGS = [
    ("geometry/boundary/x.py", "from orpheus.numerics.operator_names import X\n"),
    ("geometry/boundary/x.py", "from orpheus.numerics.measure import DiscreteMeasure\n"),
]


@pytest.mark.foundation
@pytest.mark.parametrize(("rel", "source"), _OPERATOR_ALGEBRA_REFUSED_LEGS)
def test_the_operator_algebra_predicate_refuses_each_runtime_shape(rel: str, source: str) -> None:
    assert _operator_algebra_imports(pathlib.PurePosixPath(rel), source) == [_OPERATOR_ALGEBRA]


@pytest.mark.foundation
@pytest.mark.parametrize(("rel", "source"), _OPERATOR_ALGEBRA_ADMITTED_LEGS)
def test_the_operator_algebra_predicate_admits_each_other_shape(rel: str, source: str) -> None:
    assert _operator_algebra_imports(pathlib.PurePosixPath(rel), source) == []


# ---------------------------------------------------------------------------
# The closed references import only the interface vocabulary of numerics
# ---------------------------------------------------------------------------
#
# The user's ruling of 2026-10-06 (`.claude/plans/characteristic_reference_
# architecture.md`, rulings ledger): a continuous reference is closed, purpose
# built and limited in scope, while production numerics is versatile and
# changes as production grows, so churn must stay in production. A reference
# under `derivations/continuous/` or `derivations/common/` imports from
# `orpheus.numerics` only the vocabulary it shares with its consumers (the
# question, the observables, the answer function and the traced memo); content
# identity and the scalar parsers are `orpheus.values`, below every package, so
# a reference imports them freely; its mathematics is the reference kernel's
# (`orpheus/derivations/common/dense_pencil.py`, `quadrature.py`), on numpy,
# scipy and mpmath. A change to this set is a contract change.
#
# `content` and `scalars` were members until layer DAG ruling 6 (2026-10-22)
# moved them to `orpheus.values`; a reference still normalises a key argument
# with the rule the content encoder applies (#405 P3 step 2c, C1).

REFERENCE_INTERFACE_NUMERICS: frozenset[str] = frozenset(
    {"question", "observable", "mesh_free_function", "traced_memo"}
)

_CLOSED_REFERENCE_ROOTS: tuple[str, ...] = ("derivations/continuous/", "derivations/common/")

# The production machinery a closed reference never imports at all (the user's
# ruling of 2026-10-06, "Add mesh and methods"): the discretization overlay,
# the transport vocabulary and every method package. Its input and interface
# vocabulary (data, geometry, the specification, the reference package) stays
# importable.
REFERENCE_FORBIDDEN_PACKAGES: frozenset[str] = MESH_PACKAGES | L2_PACKAGES | L3_PACKAGES

# Exempt by name, with the reason (the user's rulings on the gate's scope,
# 2026-10-06): these two MMS harnesses pose a production method, so their
# subject is a production object and they follow production's changes by
# nature. (The Sood registry's builder module, exempt until layer DAG P4 of
# 2026-10-22, was deleted with its entry.)
# `derivations/discrete/` is outside the roots for the same reason (an algebra
# of record of a production discretization).
_REFERENCE_INSULATION_EXEMPT: frozenset[str] = frozenset(
    {
        "derivations/continuous/mms/sn.py",  # the MMS harness that poses SN problems
        "derivations/continuous/mms/moc.py",  # the MMS harness that poses MoC problems
    }
)


def _imported_targets(rel_path: pathlib.PurePosixPath, src: str) -> list[str]:
    """Every dotted name a module's imports bind, at any depth, ``TYPE_CHECKING`` included.

    ``from orpheus.numerics import question`` binds ``orpheus.numerics.question``
    and ``from orpheus import numerics`` binds ``orpheus.numerics``: the imported
    NAME is part of the target when the module is a package the gate rules on, so
    a package-level import cannot pass for an interface one, nor the reverse.
    """
    targets: list[str] = []
    for node in ast.walk(ast.parse(src, filename=str(rel_path))):
        if isinstance(node, ast.ImportFrom):
            module = _absolute_module(rel_path, node)
            if not module:
                continue
            if module in ("orpheus", "orpheus.numerics"):
                targets.extend(f"{module}.{alias.name}" for alias in node.names)
            else:
                targets.append(module)
        elif isinstance(node, ast.Import):
            targets.extend(alias.name for alias in node.names)
    return targets


def _numerics_submodule(module_name: str) -> str | None:
    """The `numerics` submodule a dotted name reads, ``""`` for the package itself."""
    if module_name == "orpheus.numerics":
        return ""
    if module_name.startswith("orpheus.numerics."):
        return module_name.split(".")[2]
    return None


def _dynamic_imports(src: str) -> list[str]:
    """Module names passed as string LITERALS to ``importlib.import_module`` or ``__import__``.

    A name built at run time (concatenated, formatted) is not decidable from the
    source and is not seen: the gate's declared blind spot.
    """
    names: list[str] = []
    for node in ast.walk(ast.parse(src)):
        if not (isinstance(node, ast.Call) and node.args):
            continue
        func = node.func
        called = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", "")
        argument = node.args[0]
        if called in {"import_module", "__import__"} and isinstance(argument, ast.Constant) and isinstance(argument.value, str):
            names.append(argument.value)
    return names


def _reference_insulation_violations(
    rel: pathlib.PurePosixPath, src: str, exempt: frozenset[str] | None = None,
) -> list[str]:
    """The non-interface `numerics` imports of one closed-reference module.

    ``exempt`` is the exemption set read; the vacuous-row guard passes the set
    minus one entry to ask whether that entry still exempts a live import. The
    default ``None`` reads ``_REFERENCE_INSULATION_EXEMPT`` at CALL time.
    """
    if exempt is None:
        exempt = _REFERENCE_INSULATION_EXEMPT
    rel_str = rel.as_posix()
    if not rel_str.startswith(_CLOSED_REFERENCE_ROOTS) or rel_str in exempt:
        return []
    imported = _imported_targets(rel, src) + _dynamic_imports(src)
    violations: list[str] = []
    for module_name in imported:
        if module_name.startswith("orpheus.") and module_name.split(".")[1] in REFERENCE_FORBIDDEN_PACKAGES:
            violations.append(
                f"{rel} imports {module_name} (a closed reference never imports production "
                f"machinery: {sorted(REFERENCE_FORBIDDEN_PACKAGES)})"
            )
            continue
        submodule = _numerics_submodule(module_name)
        if submodule is None or submodule in REFERENCE_INTERFACE_NUMERICS:
            continue
        violations.append(
            f"{rel} imports {module_name} (a closed reference may import from numerics only "
            f"{sorted(REFERENCE_INTERFACE_NUMERICS)})"
        )
    return violations


@pytest.mark.foundation
def test_closed_references_import_only_the_numerics_interface() -> None:
    violations = [
        violation
        for module_path in _iter_python_modules(ORPHEUS_ROOT / "derivations")
        for violation in _reference_insulation_violations(
            pathlib.PurePosixPath(module_path.relative_to(ORPHEUS_ROOT).as_posix()),
            module_path.read_text(encoding="utf-8"),
        )
    ]
    assert not violations, "\n".join(violations)


_INSULATION_REFUSED_LEGS = [
    ("derivations/common/x.py", "from orpheus.numerics.eigenvalue import dominant_eigenpair\n"),
    ("derivations/continuous/x/y.py", "import orpheus.numerics.quadrature.registry\n"),
    ("derivations/common/x.py", "from ...numerics.eigenvalue import dominant_eigenpair\n"),
    ("derivations/common/x.py", "def f():\n    from orpheus.numerics.pencil import OperatorPencil\n"),
    ("derivations/common/x.py", "from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n    from orpheus.numerics.operator import LinearOperator\n"),
    ("derivations/common/x.py", "from orpheus.numerics import power_iteration\n"),
    ("derivations/common/x.py", "import importlib\nm = importlib.import_module('orpheus.numerics.eigenvalue')\n"),
    ("derivations/common/x.py", "from orpheus import numerics\n"),
    ("derivations/common/x.py", "from ... import numerics\n"),
]

_INSULATION_ADMITTED_LEGS = [
    ("derivations/common/x.py", "from orpheus.numerics.question import Eigen\n"),
    ("derivations/common/x.py", "from orpheus.numerics import question\n"),
    ("derivations/continuous/x/y.py", "from orpheus.numerics.traced_memo import traced_memo\n"),
    ("derivations/continuous/mms/sn.py", "from orpheus.numerics.quadrature import Quadrature\n"),
    ("derivations/discrete/sn/x.py", "from orpheus.numerics.eigenvalue import dominant_eigenpair\n"),
    ("derivations/common/x.py", "import numpy as np\nimport scipy.linalg\n"),
    ("derivations/continuous/mms/moc.py", "from orpheus.moc.geometry import MOCMesh\n"),
    ("derivations/common/x.py", "from orpheus.geometry import StructuredGeometry\nfrom orpheus.data.materials import Materials\n"),
    ("derivations/common/x.py", "from orpheus.specification import Specification\nfrom orpheus.reference import ReferenceSolution\n"),
]


_INSULATION_REFUSED_PRODUCTION_LEGS = [
    ("derivations/continuous/x/y.py", "from orpheus.mesh import Mesh1D\n"),
    ("derivations/common/x.py", "def f():\n    from orpheus.transport.fields import X\n"),
    ("derivations/continuous/x/y.py", "from orpheus.cp import CPParams\n"),
    ("derivations/continuous/x/y.py", "from orpheus import mesh\n"),
    ("derivations/common/x.py", "from ...sn.problem import SNProblem\n"),
    # The retired Sood builder: re-adding its exempt entry turns this red.
    ("derivations/continuous/sood_registry/builders.py", "from orpheus.cp import CPParams\n"),
]


@pytest.mark.foundation
@pytest.mark.parametrize(("rel", "source"), _INSULATION_REFUSED_LEGS)
def test_reference_insulation_refuses_each_import_shape(rel: str, source: str) -> None:
    violations = _reference_insulation_violations(pathlib.PurePosixPath(rel), source)
    assert len(violations) == 1 and "may import from numerics only" in violations[0], violations


@pytest.mark.foundation
@pytest.mark.parametrize(("rel", "source"), _INSULATION_REFUSED_PRODUCTION_LEGS)
def test_reference_insulation_refuses_production_machinery(rel: str, source: str) -> None:
    violations = _reference_insulation_violations(pathlib.PurePosixPath(rel), source)
    assert len(violations) == 1 and "never imports production machinery" in violations[0], violations


@pytest.mark.foundation
@pytest.mark.parametrize(("rel", "source"), _INSULATION_ADMITTED_LEGS)
def test_reference_insulation_admits_the_interface_and_the_exempt(rel: str, source: str) -> None:
    assert _reference_insulation_violations(pathlib.PurePosixPath(rel), source) == []


# ---------------------------------------------------------------------------
# Every exemption row names a live edge
# ---------------------------------------------------------------------------
#
# An exemption that outlives its edge is a hole with no witness: a WHITELIST
# row for a deleted module, or an insulation exemption for a module that no
# longer imports production machinery, admits the next import spelled there
# while every other test stays green. Each row is checked by removing it alone
# and requiring its own module to redden the predicate it exempts. Layer DAG
# P4 (2026-10-22), whose deletion of the Sood builders left two such rows.


@pytest.mark.foundation
def test_every_exemption_row_names_a_live_edge() -> None:
    stale: list[str] = []
    for rel, pkg in sorted(WHITELIST):
        path = ORPHEUS_ROOT / rel
        if not path.is_file():
            stale.append(f"WHITELIST row {(rel, pkg)}: {rel} does not exist")
            continue
        violations = _check_source(
            pathlib.PurePosixPath(rel), path.read_text(encoding="utf-8"), WHITELIST - {(rel, pkg)},
        )
        if not any(v.endswith(f" → {pkg})") for v in violations):
            stale.append(f"WHITELIST row {(rel, pkg)}: {rel} imports nothing from orpheus.{pkg}")
    for rel in sorted(_REFERENCE_INSULATION_EXEMPT):
        path = ORPHEUS_ROOT / rel
        if not path.is_file():
            stale.append(f"insulation exemption {rel!r}: the file does not exist")
            continue
        if not _reference_insulation_violations(
            pathlib.PurePosixPath(rel), path.read_text(encoding="utf-8"), _REFERENCE_INSULATION_EXEMPT - {rel},
        ):
            stale.append(f"insulation exemption {rel!r}: the module imports nothing the exemption admits")
    assert not stale, "stale exemption rows:\n" + "\n".join(stale)


# ---------------------------------------------------------------------------
# Inside a package, the module graph is a DAG (the open cycles are listed)
# ---------------------------------------------------------------------------
#
# The package linter above sees edges BETWEEN packages; a cycle between two
# modules of one package is invisible to it. Inside a package the same law
# holds: a module imports only modules whose vocabulary is smaller, so the
# graph of the package's direct children (its modules and subpackages) is
# acyclic at every import kind (top level, function body, TYPE_CHECKING). A
# cycle is a concept in the wrong module. An import of the package facade
# itself (`from orpheus.mesh import X` inside `orpheus/mesh/`) is resolved
# through the facade's own imports, so it cannot hide an edge.
#
# The cycles still open are listed EXACTLY: a new cycle is red, and a listed
# cycle that no longer exists is red too (the row is retired with the cycle).
# [M] 2026-10-22 (layer DAG rounds 2-4): the mesh, data, transport and
# diffusion cycles were broken, and the sn cycle shrank from 8 members to 5;
# these remain.

_OPEN_INTRA_PACKAGE_CYCLES: dict[str, frozenset[frozenset[str]]] = {
    # The traced-memo boot pair (numerics_untangling.md: by design, a bootstrap).
    "numerics": frozenset({frozenset({"_traced_memo_boot", "traced_memo"})}),
    # The SN hub caches values built from above it (`system`, the loss-kernel
    # gauge): the discretization welded to the posing layer (#529).
    "sn": frozenset({frozenset({"coupled_system", "loss_representation", "operators", "problem", "sweep"})}),
}


def _package_children(package_dir: pathlib.Path) -> set[str]:
    modules = {p.stem for p in package_dir.glob("*.py") if p.stem != "__init__"}
    subpackages = {p.name for p in package_dir.iterdir() if (p / "__init__.py").is_file()}
    return modules | subpackages


def _intra_package_graph(package_dir: pathlib.Path) -> dict[str, set[str]]:
    """The import graph of one package's direct children, at every import kind."""
    package = package_dir.name
    children = _package_children(package_dir)

    def targets(path: pathlib.Path) -> tuple[set[str], bool]:
        rel = pathlib.PurePosixPath(path.relative_to(ORPHEUS_ROOT).as_posix())
        found: set[str] = set()
        facade = False
        for module_name, _is_tc in _imports_of_source(rel, path.read_text(encoding="utf-8")):
            parts = module_name.split(".")
            if parts[:2] != ["orpheus", package]:
                continue
            if len(parts) == 2:
                facade = True
            elif parts[2] in children:
                found.add(parts[2])
        return found, facade

    facade_targets, _ = targets(package_dir / "__init__.py")
    graph: dict[str, set[str]] = {child: set() for child in children}
    for path in package_dir.rglob("*.py"):
        rel_parts = path.relative_to(package_dir).with_suffix("").parts
        if rel_parts == ("__init__",):
            continue
        child = rel_parts[0]
        found, facade = targets(path)
        if facade:
            found |= facade_targets
        graph[child] |= found - {child}
    return graph


def _cycles(graph: dict[str, set[str]]) -> frozenset[frozenset[str]]:
    """The strongly connected components of ``graph`` with more than one node."""
    index: dict[str, int] = {}
    low: dict[str, int] = {}
    stack: list[str] = []
    on_stack: set[str] = set()
    found: set[frozenset[str]] = set()

    def visit(node: str) -> None:
        index[node] = low[node] = len(index)
        stack.append(node)
        on_stack.add(node)
        for nxt in sorted(graph[node]):
            if nxt not in index:
                visit(nxt)
                low[node] = min(low[node], low[nxt])
            elif nxt in on_stack:
                low[node] = min(low[node], index[nxt])
        if low[node] == index[node]:
            component = set()
            while True:
                member = stack.pop()
                on_stack.discard(member)
                component.add(member)
                if member == node:
                    break
            if len(component) > 1:
                found.add(frozenset(component))

    for node in sorted(graph):
        if node not in index:
            visit(node)
    return frozenset(found)


@pytest.mark.foundation
def test_the_cycle_finder_finds_a_planted_cycle() -> None:
    """Positive control for the predicate the tree scan runs."""
    graph = {"axis": {"structured"}, "structured": {"face_laws"}, "face_laws": {"axis"}, "leaf": set()}
    assert _cycles(graph) == frozenset({frozenset({"axis", "structured", "face_laws"})})
    assert _cycles({"a": {"b"}, "b": set()}) == frozenset()


@pytest.mark.foundation
def test_every_package_module_graph_is_acyclic_but_the_listed_cycles() -> None:
    packages = sorted(p for p in ORPHEUS_ROOT.iterdir() if (p / "__init__.py").is_file())
    # The denominator is stated, so an empty walk cannot read as a clean tree.
    assert len(packages) > 10, [p.name for p in packages]
    new: list[str] = []
    stale: list[str] = []
    for package_dir in packages:
        found = _cycles(_intra_package_graph(package_dir))
        listed = _OPEN_INTRA_PACKAGE_CYCLES.get(package_dir.name, frozenset())
        new += [f"{package_dir.name}: {sorted(c)}" for c in found - listed]
        stale += [f"{package_dir.name}: {sorted(c)}" for c in listed - found]
    assert not new, (
        "a module cycle inside a package (a concept sits in the wrong module): " + "; ".join(new)
    )
    assert not stale, (
        "a listed open cycle no longer exists; retire its _OPEN_INTRA_PACKAGE_CYCLES row: " + "; ".join(stale)
    )
