"""Import-linter: enforce the L0 / L1 / (input) / L2 / L3 layer contract.

Layers (per plan §P3.0):
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
layer -> input). Forbidden edges raise the parametrised test.

Tolerances:
  - TYPE_CHECKING imports of L3-only types inside L1/L2 modules are
    permitted (string annotations don't create runtime edges).
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
L0_PACKAGES: frozenset[str] = frozenset({"derivations"})
L1_PACKAGES: frozenset[str] = frozenset({"numerics"})
INPUT_PACKAGES: frozenset[str] = frozenset({"geometry", "data"})
MESH_PACKAGES: frozenset[str] = frozenset({"mesh"})
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

FORBIDDEN_EDGES: dict[str, frozenset[str]] = {
    "numerics": MESH_PACKAGES | L2_PACKAGES | L3_PACKAGES,
    "geometry": MESH_PACKAGES | L2_PACKAGES | L3_PACKAGES,
    "data": MESH_PACKAGES | L2_PACKAGES | L3_PACKAGES,
    "mesh": L2_PACKAGES | L3_PACKAGES,
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
        # RETIRE_IN_P3_FOLLOWUP — sood_registry lazy-imports CPParams
        # inside a function body to avoid CP transitive deps at import time.
        ("derivations/continuous/sood_registry/builders.py", "cp"),
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
                results.append((module, in_tc))
        elif isinstance(node, ast.Import):
            for alias in node.names:
                results.append((alias.name, in_tc))
        for child in ast.iter_child_nodes(node):
            _visit(child, in_tc)

    _visit(tree, False)
    return results


def _check_source(rel: pathlib.PurePosixPath, src: str) -> list[str]:
    """The layer violations of one module's source, ``rel`` relative to ``orpheus/``."""
    src_pkg = _top_level_package(rel)
    if src_pkg not in FORBIDDEN_EDGES:
        return []
    forbidden = FORBIDDEN_EDGES[src_pkg]
    rel_str = rel.as_posix()
    violations: list[str] = []
    for module_name, is_tc in _imports_of_source(rel, src):
        if not module_name.startswith("orpheus."):
            continue
        tgt_pkg = module_name.split(".")[1]
        if tgt_pkg not in forbidden:
            continue
        # TYPE_CHECKING tolerance: L1/L2 importing L3 names for typing.
        if is_tc and src_pkg in (L1_PACKAGES | L2_PACKAGES):
            continue
        if (rel_str, tgt_pkg) in WHITELIST:
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


# ---------------------------------------------------------------------------
# The numerics -> geometry back-edge, and the discipline that makes it safe
# ---------------------------------------------------------------------------
#
# `numerics/__init__.py` imports `symmetry`, and `symmetry` imports
# `geometry.transformation` (the rigid-motion core). That is a genuine package
# CYCLE: importing `orpheus.numerics` runs `orpheus.geometry.__init__`, which
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

_INPUT_PACKAGE_ROOTS = sorted(INPUT_PACKAGES | MESH_PACKAGES)


@pytest.mark.foundation
@pytest.mark.parametrize("package", _INPUT_PACKAGE_ROOTS)
def test_input_layer_imports_numerics_only_by_submodule(package: str) -> None:
    """An input-layer package may not import the `numerics` PACKAGE itself.

    `from orpheus.numerics.measure import DiscreteMeasure` is fine;
    `from orpheus.numerics import DiscreteMeasure` is not, because the second
    form is an attribute lookup that fails while `numerics/__init__` is still
    executing — which is exactly the state it is in when `symmetry` reaches
    across to `geometry.transformation`.
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
        "orpheus.numerics.symmetry",
        # R2 of #434 (2026-09-03): the four entries the V1 cycle killed that
        # the list above could not see — `[M]` with `manifold -> symmetry` at
        # module scope AND symmetry still importing the axis table from
        # manifold, 6 of 9 entry points died while `import orpheus` stayed
        # green (order-dependent; plan-authoring §6d).
        "orpheus.numerics.manifold",
        "orpheus.numerics.measure",
        "orpheus.numerics.invariance",
        "orpheus.numerics.quadrature.registry",
        "orpheus.geometry",
        "orpheus.geometry.transformation",
        # P1 step 1 of #405 (2026-09-25): the mesh became its own package,
        # importing geometry (whose boundary package now holds `BC`) and
        # numerics; each moved module, and the boundary package, cold.
        "orpheus.geometry.boundary",
        "orpheus.mesh",
        "orpheus.mesh.structured",
        "orpheus.mesh.factories",
        "orpheus.mesh.axis",
        "orpheus.sn.solver",
        # step 2 of the consumers campaign (2026-09-13): the Strategy value
        # module — solver -> splitting -> operators/loss_representation, no
        # runtime edge back to the hub (plan-authoring §6d).
        "orpheus.sn.splitting",
        # step 3 U1 (2026-09-17): the Solution-tier types — the gauges (the
        # section that picked the representative) and the kind-typed outcomes
        # + the exit certificate; numerics-tier, importing posing/pencil/
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
        "orpheus.numerics.scalars",
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
]

_ADMITTED_LEGS = [
    ("mesh/x.py", "from orpheus.geometry.boundary import BC\n"),
    ("mesh/x.py", "from orpheus.numerics.measure import DiscreteMeasure\n"),
    ("mesh/x.py", "from .structured import Mesh1D\n"),
    ("transport/x.py", "from orpheus.mesh import Mesh1D\n"),
    ("derivations/x.py", "from orpheus.mesh import Mesh1D\n"),
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
