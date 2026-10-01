"""The body CP computes, for tests parametrised over the coordinate system.

CP's slab kernel puts a mirror at the left face whatever is declared there, so
the slab CP computes is half of a symmetric slab: the left face is reflective
and only the right face carries the law under test. A solid cylinder or sphere
has one face, the outer one. ``cp_body`` declares exactly that, so a test never
states a left law CP would replace (``orpheus.cp.solver._refuse_a_law_cp_drops``,
#513).
"""

from __future__ import annotations

from orpheus.geometry import BC, CoordSystem, StructuredGeometry


def cp_body(coord: CoordSystem, breakpoints, mat_ids, outer) -> StructuredGeometry:
    """A solid body with ``outer`` on its outer face, and a mirror at a slab's left face."""
    if coord is CoordSystem.CARTESIAN:
        return StructuredGeometry.slab(breakpoints, mat_ids, left=BC.reflective, right=outer)
    return StructuredGeometry.uniform_boundary(coord, breakpoints, mat_ids, outer)
