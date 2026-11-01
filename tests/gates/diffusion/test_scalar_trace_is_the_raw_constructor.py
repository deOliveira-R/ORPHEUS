r"""Two owners, one space: ``DiffusionMesh.scalar_trace`` and ``full_field_space`` equal the raw constructors on the
mesh's own data (layer DAG Phase 4). DEMOTED with respect to the rule, as the SN twin
(``tests/gates/sn/mesh/test_the_mints_are_the_raw_constructors.py``): the mint is one call to
``phase_spaces.scalar_trace_space``; what this row tests is that the hub reads it and that its material-mesh data
equals a standalone ``MaterialMesh``'s. First red: ``diffusion_trace_reads_a_one_group_mesh`` (the mint re-spelled
with unit areas: red on the sphere and the cylinder, null on the slab, whose two faces have area 1).
"""
from __future__ import annotations

import pytest

from orpheus.diffusion import DiffusionMesh
from orpheus.numerics.spaces.full_field_space import FullFieldSpace
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.phase_spaces import scalar_trace_space
from tests.gates.transport import _phase_space_cases as cases
from tests.gates.transport._materials import CURVED_REFLECTIVE

pytestmark = pytest.mark.foundation


@pytest.mark.parametrize("geometry", ["slab", "sphere", "hollow_sphere", "cylinder", "hollow_cylinder"])
@pytest.mark.rests_on("tests/gates/transport/test_phase_spaces.py::test_the_scalar_trace_is_the_numerics_trace_over_the_hand_faces_and_areas")
def test_the_scalar_trace_and_the_composite_are_the_constructors_on_the_meshs_data(geometry: str) -> None:
    material_mesh = cases.GEOMETRIES[geometry]()
    conditions = {} if geometry == "slab" else CURVED_REFLECTIVE
    mesh = cases.legacy_mesh(material_mesh)
    hub = DiffusionMesh(mesh, material_mesh.materials, boundary_conditions=conditions)
    standalone = MaterialMesh(mesh, material_mesh.materials)
    trace = scalar_trace_space(standalone)
    if hub.scalar_trace != trace or hash(hub.scalar_trace) != hash(trace):
        pytest.fail(f"{geometry}: DiffusionMesh.scalar_trace differs from scalar_trace_space on the same mesh")
    composite = FullFieldSpace.from_blocks(standalone.bulk_space, trace)
    if hub.full_field_space != composite:
        pytest.fail(f"{geometry}: DiffusionMesh.full_field_space differs from the composite of the constructors")
