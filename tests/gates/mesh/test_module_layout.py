"""The mesh is its own package, with no shim left at its old homes.

P1 step 1 of #405 (``.claude/plans/reference_cache.md``, "P1, the carve
order"): the mesh is the discretisation overlay on the geometry, so it lives
in :mod:`orpheus.mesh`, which imports the geometry and is never imported by
it. A re-export left at an old home would let a consumer keep the old spelling
and hide the layer from the linter, so each old home is checked empty; the
boundary tag ``BC`` and its package ``orpheus.geometry.boundary`` are gone (wave 2: the
geometry facade exports no response vocabulary, and the package was deleted in WF2.4).
"""
from __future__ import annotations

import importlib.util
import types

import pytest

import orpheus.geometry
import orpheus.transport.mesh


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
        (orpheus.geometry, ("Mesh1D", "Mesh2D", "RegionMesh", "pwr_pin_2d")),
        (
            orpheus.transport.mesh,
            ("AxisMesh", "RadialAxisMesh", "Axis1D", "AxisCoord", "FaceLabel", "AXIS_NAMES"),
        ),
    ],
    ids=["orpheus.geometry", "orpheus.transport.mesh"],
)
def test_no_old_package_re_exports_a_mesh_name(package: types.ModuleType, names: tuple[str, ...]) -> None:
    leaked = [name for name in names if hasattr(package, name)]
    assert not leaked, f"{package.__name__} still exports {leaked}; import them from orpheus.mesh"


def test_the_boundary_tag_and_its_package_are_gone() -> None:
    assert importlib.util.find_spec("orpheus.geometry.boundary") is None
    assert not hasattr(orpheus.geometry, "BC") and "BC" not in orpheus.geometry.__all__
    assert importlib.util.find_spec("orpheus.geometry.deck") is not None  # the positive control: the successor home resolves


def test_the_mesh_package_exports_the_moved_names() -> None:
    import orpheus.mesh

    for name in ("Mesh1D", "Mesh2D", "Mesher", "pwr_pin_2d", "AxisMesh", "RadialAxisMesh",
                 "Axis1D", "AxisCoord", "FaceLabel"):
        assert getattr(orpheus.mesh, name).__module__.startswith("orpheus.mesh."), name
