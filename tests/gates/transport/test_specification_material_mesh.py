"""A two-interval specification assembles a ``MaterialMesh`` (#405 P1 step 8, S8.4 (i)).

Specified by the test-architect (2026-10-02, ``.claude/plans/reference_p1_spec.md`` §1.8, S8.4) as one of the two rows of the
specification's assembly rung, an integration rung green on arrival (no first red of its own): it shows that a specification
carries what a method needs, through the test's own lift, since posing a method from a specification is P4's. Its teeth are the
reference side's activation legs and the battery arms named in the spec. Moved from
``tests/gates/specification/test_specification_assembly.py`` by the layer-DAG re-filing (Phase 1b, 2026-10-26), body unchanged:
the claim is about the ``MaterialMesh`` built from the specification's geometry and materials, so the row lives beside
``MaterialMesh``'s gates. The other row of the rung, the infinite-medium k against the dense pencil, is
``tests/gates/cross_method/test_specification_infinite_medium_k.py``.
"""

from __future__ import annotations

import pytest

from orpheus.data.materials import Materials
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.question import Eigen
from orpheus.specification import GeometrySpecification
from tests.gates._content_identity_helpers import require
from tests.gates.specification._fixtures import (
    every_fission as _k,
    moderator,
    slab2,
    two_materials as _two_materials,
    vacuum_exterior,
)

pytestmark = pytest.mark.l1


def test_s8_4_a_two_interval_specification_assembles_a_material_mesh() -> None:
    from orpheus.transport.mesh import MaterialMesh

    spec = GeometrySpecification(Materials({**_two_materials(), 9: moderator()}), slab2(), Eigen(_k()), exterior=vacuum_exterior(slab2()))
    require(sorted(spec.materials.ids) == [0, 1], f"the spectator 9 was kept: {sorted(spec.materials.ids)}")
    mesh = Mesher(spec.geometry).partition(CellsByCount.uniform_width(3)).mesh
    mm = MaterialMesh(mesh, spec.materials)
    require(mm.ng == spec.materials.uniform_group_count() == 2, f"ng {mm.ng}")
    require(mm.materials == spec.materials, "the mesh holds another declaration")
    occupied = sorted(m for m, idx in mm.cells_by_material.items() if idx[0].size)
    require(occupied == [0, 1], f"occupied materials {occupied}")
