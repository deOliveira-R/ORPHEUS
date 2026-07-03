"""Import-linter: enforce the L0 / L1 / (input) / L2 / L3 layer contract.

Layers (per plan §P3.0):
  L0  (derivations/)   — Branch-1 references (SymPy / mpmath);
                         BELOW L1 in the import hierarchy.
  L1  (numerics/)      — math primitives; knows no neutrons.
  (input) (geometry/, data/) — geometry + nuclear data.
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
from collections.abc import Iterable

import pytest

ORPHEUS_ROOT = pathlib.Path(__file__).parent.parent / "orpheus"

# Layer assignment.
L0_PACKAGES: frozenset[str] = frozenset({"derivations"})
L1_PACKAGES: frozenset[str] = frozenset({"numerics"})
INPUT_PACKAGES: frozenset[str] = frozenset({"geometry", "data"})
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
    "numerics": L2_PACKAGES | L3_PACKAGES,
    "geometry": L2_PACKAGES | L3_PACKAGES,
    "data": L2_PACKAGES | L3_PACKAGES,
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
        # RETIRE_IN_P3_FOLLOWUP — inline-import benchmark cross-check
        # (derivations/cases/diffusion.py uses the production solver
        # as a black-box reference inside a function body).
        ("derivations/continuous/cases/diffusion.py", "diffusion"),
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


def _imports_of(module_path: pathlib.Path) -> list[tuple[str, bool]]:
    """Parse imports, marking TYPE_CHECKING-guarded ones."""
    src = module_path.read_text()
    tree = ast.parse(src, filename=str(module_path))
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
        if isinstance(node, ast.ImportFrom) and node.module:
            results.append((node.module, in_tc))
        elif isinstance(node, ast.Import):
            for alias in node.names:
                results.append((alias.name, in_tc))
        for child in ast.iter_child_nodes(node):
            _visit(child, in_tc)

    _visit(tree, False)
    return results


def _check_module(module_path: pathlib.Path) -> list[str]:
    rel = module_path.relative_to(ORPHEUS_ROOT)
    src_pkg = _top_level_package(rel)
    if src_pkg not in FORBIDDEN_EDGES:
        return []
    forbidden = FORBIDDEN_EDGES[src_pkg]
    rel_str = str(rel).replace("\\", "/")
    violations: list[str] = []
    for module_name, is_tc in _imports_of(module_path):
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


_ALL_MODULES = sorted(_iter_python_modules(ORPHEUS_ROOT))


@pytest.mark.foundation
@pytest.mark.parametrize("module_path", _ALL_MODULES, ids=str)
def test_no_forbidden_imports(module_path: pathlib.Path) -> None:
    violations = _check_module(module_path)
    assert not violations, "\n".join(violations)
