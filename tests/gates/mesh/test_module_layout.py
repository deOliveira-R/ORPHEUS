"""The mesh is its own package, with no shim left at its old homes.

P1 step 1 of #405 (``.claude/plans/reference_cache.md``, "P1, the carve
order"): the mesh is the discretisation overlay on the geometry, so it lives
in :mod:`orpheus.mesh`, which imports the geometry and is never imported by
it. A re-export left at an old home would let a consumer keep the old spelling
and hide the layer from the linter, so each old home is checked empty; the
boundary tag :class:`~orpheus.geometry.boundary.BC` moved into the geometry's
boundary package and is ONE definition however it is reached.
"""
from __future__ import annotations

import ast
import importlib
import importlib.util
import pathlib

import pytest

import orpheus.mesh

pytestmark = pytest.mark.foundation


@pytest.mark.parametrize(
    "old_module",
    ["orpheus.geometry.mesh", "orpheus.geometry.factories", "orpheus.transport.mesh.axis"],
)
def test_no_module_remains_at_an_old_home(old_module: str) -> None:
    assert importlib.util.find_spec(old_module) is None, (
        f"{old_module} still exists; the mesh modules live in orpheus.mesh "
        f"and no shim is kept at their old homes"
    )


@pytest.mark.parametrize(
    ("package", "names"),
    [
        ("orpheus.geometry", ("Mesh1D", "Mesh2D", "RegionMesh", "pwr_pin_2d")),
        (
            "orpheus.transport.mesh",
            ("AxisMesh", "RadialAxisMesh", "Axis1D", "AxisCoord", "FaceLabel", "AXIS_NAMES"),
        ),
    ],
)
def test_no_old_package_re_exports_a_mesh_name(package: str, names: tuple[str, ...]) -> None:
    module = importlib.import_module(package)
    leaked = [name for name in names if hasattr(module, name)]
    assert not leaked, f"{package} still exports {leaked}; import them from orpheus.mesh"


def test_the_boundary_tag_is_one_definition_in_the_geometry_boundary_package() -> None:
    import orpheus.geometry
    import orpheus.geometry.boundary

    assert orpheus.geometry.boundary.BC.__module__ == "orpheus.geometry.boundary._tag"
    assert orpheus.geometry.BC is orpheus.geometry.boundary.BC


def test_the_mesh_package_exports_the_moved_names() -> None:
    import orpheus.mesh

    for name in ("Mesh1D", "Mesh2D", "Mesher", "pwr_pin_2d", "AxisMesh", "RadialAxisMesh",
                 "Axis1D", "AxisCoord", "FaceLabel"):
        assert getattr(orpheus.mesh, name).__module__.startswith("orpheus.mesh."), name


# ── The mesh package's module graph is a DAG ──────────────────────────────
#
# A module of the mesh package imports only modules whose vocabulary is
# smaller than its own, so the intra-package import graph is acyclic at every
# import kind: top level, inside a function (lazy) and under
# ``if TYPE_CHECKING:``. A cycle means a concept sits in the wrong module (the
# layer-DAG campaign, ``.claude/plans/layer_dag.md``). Its first red was the
# legacy-mesh adapter living in the axis primitive's module:
# ``axis -lazy-> structured -> face_laws -> axis``.

_MESH_ROOT = pathlib.Path(orpheus.mesh.__file__).parent


def _imported_mesh_modules(source: str, module: str) -> set[str]:
    """Every ``orpheus.mesh.<m>`` a source imports, at every import kind."""
    package = "orpheus.mesh"
    targets: set[str] = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Import):
            names = [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom):
            base = node.module or ""
            if node.level:
                base = ".".join([package, *([base] if base else [])])
            names = [base, *(f"{base}.{alias.name}" for alias in node.names)]
        else:
            continue
        for name in names:
            parts = name.split(".")
            if parts[:2] == ["orpheus", "mesh"] and len(parts) >= 3:
                targets.add(parts[2])
    targets.discard(module)
    return targets


def _mesh_module_graph() -> dict[str, set[str]]:
    modules = {path.stem for path in _MESH_ROOT.glob("*.py")} - {"__init__"}
    return {
        module: _imported_mesh_modules(
            (_MESH_ROOT / f"{module}.py").read_text(encoding="utf-8"), module
        )
        & modules
        for module in sorted(modules)
    }


def _a_cycle(graph: dict[str, set[str]]) -> list[str] | None:
    """One import cycle of ``graph`` as a closed path, or ``None``."""
    state: dict[str, int] = {}  # 1 on the DFS stack, 2 finished
    stack: list[str] = []

    def visit(node: str) -> list[str] | None:
        state[node] = 1
        stack.append(node)
        for nxt in sorted(graph.get(node, ())):
            if state.get(nxt) == 1:
                return [*stack[stack.index(nxt):], nxt]
            if nxt not in state and (cycle := visit(nxt)) is not None:
                return cycle
        stack.pop()
        state[node] = 2
        return None

    for node in sorted(graph):
        if node not in state and (cycle := visit(node)) is not None:
            return cycle
    return None


def test_the_cycle_finder_finds_the_adapter_cycle() -> None:
    """Positive control: the finder reports the cycle the adapter made."""
    graph = {"axis": {"structured"}, "structured": {"face_laws"}, "face_laws": {"axis"}}
    assert _a_cycle(graph) == ["axis", "structured", "face_laws", "axis"]
    lazy = "def f():\n    from orpheus.mesh.structured import Mesh1D\n"
    assert _imported_mesh_modules(lazy, "axis") == {"structured"}


def test_the_mesh_module_graph_is_acyclic_at_every_import_kind() -> None:
    graph = _mesh_module_graph()
    assert "axis" in graph and "structured" in graph, sorted(graph)
    cycle = _a_cycle(graph)
    assert cycle is None, (
        "the mesh package's module imports form a cycle "
        f"({' -> '.join(cycle or [])}): a concept sits in the wrong module"
    )
