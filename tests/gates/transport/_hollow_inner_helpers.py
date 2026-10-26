"""The hollow body shared by SN's and diffusion's refusals of a declared inner condition (#511; layer-DAG re-filing, Phase 1b).

SN and diffusion both build on the radial axis, which has one trace slot (the outer surface), so both refuse a declared inner
condition through the one ``hollow_inner_refusal`` of ``orpheus/transport/boundary_conditions/admission.py``. Their gates,
``tests/gates/sn/architecture/test_sn_hollow_inner_law.py`` and ``tests/gates/diffusion/test_diffusion_hollow_inner_law.py``, were
one file, ``tests/gates/mesh/test_hollow_inner_law.py``, until the re-filing; the two methods may not import each other, so the
body, its tags and the table of refused conditions they share live here, at the layer of the shared refusal. The definitions are
that file's, unchanged.
"""
from __future__ import annotations

from orpheus.data.macro_xs.synthetic import get_mixture
from orpheus.geometry import CoordSystem, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesh1D, Mesher
from orpheus.transport.boundary_conditions import WallCondition

_CURVILINEAR = (CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL)

#: The inner face of a hollow body is the tag 1 and the outer face the tag 0 (rule R6: an ``inner=`` declaration is a tag, face-inventory order is
#: inner first); the outer wall is vacuum.
_INNER, _OUTER = 1, 0

_REFUSED_INNER_LAWS = {
    "vacuum": WallCondition.vacuum(),
    "white": WallCondition.white(0.5),
    "specular_half": WallCondition.specular(0.5),
}


def _conditions(inner) -> dict:
    return {_INNER: inner, _OUTER: WallCondition.vacuum()}


def _materials():
    return {0: get_mixture("B", "2g")}


def _hollow_mesh(coord: CoordSystem) -> Mesh1D:
    """A hollow body [0.5, 2.0] declared on a geometry, 8 cells; the inner face carries tag 1, the outer face tag 0."""
    geometry = StructuredGeometry(
        coord=coord,
        breakpoints=(0.5, 2.0),
        mat_ids=(0,),
        face_roles=(_INNER, _OUTER),
    )
    return Mesher(geometry).partition(CellsByCount.uniform_volume(8)).mesh
