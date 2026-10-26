"""Diffusion refuses a hollow body's declared inner condition unless it is the unit specular wall (#511).

The radial axis has one trace slot, the outer surface, so the declared inner condition of a hollow cylinder or sphere has no
operator to land in, and diffusion, built on the same axes as SN, would compute a reflective cavity whatever the tag held. The
refusal is a declared scope boundary (``SCOPE-BOUNDARY[guard]``), stated by the shared ``hollow_inner_refusal``.

**Diffusion's leg was a strict ``xfail`` in the S3 build, and is a refusal since the S3 integration.** The section 4.2 draft
admitted diffusion at every site, so diffusion computed the same reflective cavity whatever the tag held (``[M]`` 2026-10-23,
``solve_diffusion_1d`` on a hollow body, material A 2G, 8 cells on [0.5, 2.0], outer tag vacuum: the eigenvalue was bit-identical for
an inner ``vacuum()``, ``specular(1.0)`` and ``white(0.5)``, 3 of 3 conditions on each of 2 of 2 coordinate systems). The stencil
assembles its boundary closures from the trace inventory, which omits the inner face, so that face's flow slot stays zero (zero net
current, :math:`J^- = J^+`): diffusion's table now admits only ``specular(1.0)`` there and refuses the rest through the shared
``hollow_inner_refusal`` (#511), the refusal SN's table states.

Moved from ``tests/gates/mesh/test_hollow_inner_law.py`` by the layer-DAG re-filing (Phase 1b, 2026-10-26): the refusal is
diffusion's, so the rows live with diffusion's gates; the bodies are unchanged. The hollow body, its tags and the refused conditions shared
with SN are in ``tests/gates/transport/_hollow_inner_helpers.py``; the ``None``-role refusal of the mesh stayed in
``tests/gates/mesh/test_hollow_inner_law.py``.

SN's refusal of the same condition, and the void-cavity verification of the admitted one, are
``tests/gates/sn/architecture/test_sn_hollow_inner_law.py``."""
from __future__ import annotations

import pytest

from orpheus.diffusion.augmented_mesh import DiffusionMesh
from orpheus.geometry import CoordSystem
from orpheus.transport.boundary_conditions import BoundaryRefused, WallCondition
from tests.gates.transport._hollow_inner_helpers import _REFUSED_INNER_LAWS, _conditions, _hollow_mesh, _materials

pytestmark = pytest.mark.foundation


@pytest.mark.catches("ERR-120")
@pytest.mark.parametrize("law", _REFUSED_INNER_LAWS, ids=str)
def test_diffusion_refuses_it_too(law):
    """Diffusion reads the same axes, so it would drop the same condition: its table refuses it, citing #511 as SN's does."""
    mesh = _hollow_mesh(CoordSystem.SPHERICAL)
    with pytest.raises(BoundaryRefused) as info:
        DiffusionMesh(mesh, _materials(), boundary_conditions=_conditions(_REFUSED_INNER_LAWS[law]))
    assert info.value.admission.issue == 511


def test_diffusion_admits_the_unit_specular_cavity_at_the_inner_site():
    """The positive control of the leg above: the one condition both tables admit at the inner site constructs (so the refusal is
    the inner site's and not a construction that is broken)."""
    mesh = _hollow_mesh(CoordSystem.SPHERICAL)
    assert DiffusionMesh(mesh, _materials(), boundary_conditions=_conditions(WallCondition.specular(1.0))).bc.keys() == {"xmax"}
