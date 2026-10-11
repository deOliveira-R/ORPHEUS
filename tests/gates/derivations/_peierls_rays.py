"""One-ray readings of the Peierls ray layer, for the hand-algebra rows of the direct geometry tests (#578, P2).

Until #578 the Peierls primitives ``rho_max``, ``rho_inner_intersections``,
``optical_depth_along_ray`` and ``source_position`` answered one ray at a time;
they were retired onto :meth:`CurvilinearGeometry.rays_from`, which returns
the geometric kernel's :class:`~orpheus.geometry.chord.Rays` from one observer. The
direct tests that pinned the retired primitives keep their hand-algebra
expected values and read the successor through these one-element batches, so
each helper is a reading of the kernel's ``Rays``, never a second spelling of
the geometry.
"""
from __future__ import annotations

import numpy as np

#: An outer radius beyond every observer the cavity rows place, for the readings that do not depend on it (the walls
#: of a hollow body's cavity lie ahead of the observer before any outer surface).
_FAR = 10.0


def _rays(geometry, r_obs: float, cos_omega: float, radii):
    return geometry.rays_from(float(r_obs), np.array([float(cos_omega)]), np.asarray(radii, dtype=float))


def depth(geometry, r_obs: float, cos_omega: float, rho: float, radii, sig_t) -> float:
    """The optical depth from the observer to the distance ``rho`` along the ray: the kernel's
    :meth:`~orpheus.geometry.chord.Chord.line_integral` over :meth:`~orpheus.geometry.chord.Rays.lengths_to`; the
    cavity and the outside are void."""
    rays = _rays(geometry, r_obs, cos_omega, radii)
    return float(rays.chord.line_integral(np.asarray(sig_t, dtype=float), rays.lengths_to(np.array([float(rho)])))[0])


def exit_distance(geometry, r_obs: float, cos_omega: float, R: float) -> float:
    """The distance to the outer surface ``R`` (:attr:`~orpheus.geometry.chord.Rays.ahead`), a hollow body's cavity
    crossed."""
    return float(_rays(geometry, r_obs, cos_omega, [R]).ahead[0].sum())


def position(geometry, r_obs: float, rho: float, cos_omega: float) -> float:
    """The radius at the distance ``rho`` along the ray (:meth:`~orpheus.geometry.chord.Rays.orbit_coordinate_at`)."""
    return float(_rays(geometry, r_obs, cos_omega, [_FAR]).orbit_coordinate_at(np.array([[float(rho)]]))[0, 0])


def cavity_walls(geometry, r_obs: float, cos_omega: float) -> tuple[float | None, float | None]:
    """The distances to the near and far walls of the cavity ahead of the observer, ``None`` where absent.

    The near wall closes the run :attr:`~orpheus.geometry.chord.Rays.to_wall` where the first wall ahead is the
    cavity's (``wall == 0``); the far wall adds the cavity slot's length ahead. An observer inside the cavity has no
    near wall, and its far wall is the cavity slot ahead of it. A ray whose first wall is the outer surface, or a solid
    body or a slab, has neither.
    """
    rays = _rays(geometry, r_obs, cos_omega, [_FAR])
    partition = rays.chord.partition
    if not (partition.chart.acts_on_kept_space and partition.breakpoints[0] > 0.0):
        return None, None
    cavity = rays.chord.slot_region[0] == partition.inner_exterior
    ahead_in_cavity = float(np.sum(rays.ahead[0][cavity]))
    if int(partition.region_containing(np.array(float(r_obs)))) == partition.inner_exterior:
        return None, ahead_in_cavity
    if int(rays.wall[0]) != 0:
        return None, None
    near = float(rays.to_wall[0].sum())
    return near, near + ahead_in_cavity
