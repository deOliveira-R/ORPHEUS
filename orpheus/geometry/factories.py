"""Mesh construction factories — 2-D Cartesian only.

Phase F retired the 1-D ``Zone`` / ``mesh1d_from_zones`` /
``pwr_pin_equivalent`` / ``pwr_slab_half_cell`` / ``homogeneous_1d``
/ ``slab_fuel_moderator`` factories. The 1-D path is now
:class:`~orpheus.geometry.structured_geometry.StructuredGeometry` →
:meth:`~orpheus.geometry.mesh.Mesh1D.from_geometry`, with
:meth:`StructuredGeometry.wigner_seitz_pin_cell` and
:meth:`StructuredGeometry.pwr_slab_half_cell` for the conventional
PWR pin-cell shapes.

What survives here:

* :func:`pwr_pin_2d` — 2-D Cartesian factory. There is no 2-D
  :class:`StructuredGeometry` yet, so this stays as a standalone
  helper.
* :func:`_subdivide_zone` — private equal-volume subdivision helper
  used internally by :meth:`Mesh1D.from_geometry`. Kept as the single
  algebraic invariant for "equal-volume cells" across all three
  coordinate systems (catches ERR-020 round-trip drift).
"""

from __future__ import annotations

import numpy as np

from .coord import CoordSystem
from .mesh import Mesh1D, Mesh2D


# ── Equal-volume subdivision (private helper for Mesh1D.from_geometry) ─

def _subdivide_zone(
    inner: float,
    outer: float,
    n: int,
    coord: CoordSystem,
) -> tuple[np.ndarray, np.ndarray]:
    """Return *n + 1* edge positions AND *n* exact cell volumes for a region.

    Subdivision guarantees equal-volume cells in each coordinate
    system — and this function returns volumes computed **directly
    from the algebraic invariant**, not re-derived from the edges
    after the fact. Deriving from edges via
    :func:`~orpheus.geometry.coord.compute_volumes_1d` loses ~1 ULP
    per cell because ``cbrt(x)**3 != x`` exactly (and likewise
    ``sqrt(x)**2``), which breaks the "equal-volume region" property
    at ``rtol=1e-14``.

    * Cartesian:   ``x_k = inner + k/n * (outer - inner)``,
      ``V_cell = (outer - inner) / n``.
    * Cylindrical: ``r_k = sqrt(inner^2 + k/n * (outer^2 - inner^2))``,
      ``V_cell = π (outer^2 - inner^2) / n``.
    * Spherical:   ``r_k = cbrt(inner^3 + k/n * (outer^3 - inner^3))``,
      ``V_cell = (4/3) π (outer^3 - inner^3) / n``.

    Each cell gets the same ``V_cell`` (a scalar broadcast), so every
    cell in the region is **bit-identical** by construction.
    """
    fracs = np.linspace(0.0, 1.0, n + 1)
    match coord:
        case CoordSystem.CARTESIAN:
            edges = inner + fracs * (outer - inner)
            v_cell = (outer - inner) / n
        case CoordSystem.CYLINDRICAL:
            edges = np.sqrt(inner**2 + fracs * (outer**2 - inner**2))
            v_cell = np.pi * (outer**2 - inner**2) / n
        case CoordSystem.SPHERICAL:
            edges = np.cbrt(inner**3 + fracs * (outer**3 - inner**3))
            v_cell = (4.0 / 3.0) * np.pi * (outer**3 - inner**3) / n
        case _:
            raise ValueError(f"Unknown coordinate system: {coord}")
    volumes = np.full(n, v_cell, dtype=float)
    return edges, volumes


# ── 2-D Cartesian PWR pin factory ────────────────────────────────────

def pwr_pin_2d(
    radii: list[float] | None = None,
    mat_ids: list[int] | None = None,
    pitch: float = 3.6,
    n_cells: int = 10,
) -> Mesh2D:
    """2-D Cartesian mesh from concentric annular regions.

    Each cell in the uniform (n_cells x n_cells) grid is assigned a
    material ID based on its distance from the pin centre (pitch / 2).

    Parameters
    ----------
    radii : list[float], optional
        Outer radii of each annular region.  Default: [0.9, 1.1]
        (fuel, clad; everything beyond is coolant).
    mat_ids : list[int], optional
        Material ID for each annulus, plus one for the region beyond
        the outermost radius.  Default: [2, 1, 0].
    pitch : float
        Unit cell side length (cm).
    n_cells : int
        Number of mesh cells per side.
    """
    if radii is None:
        radii = [0.9, 1.1]
    if mat_ids is None:
        mat_ids = [2, 1, 0]

    if len(mat_ids) != len(radii) + 1:
        raise ValueError(
            f"len(mat_ids)={len(mat_ids)} must equal len(radii)+1={len(radii) + 1}"
        )

    delta = pitch / n_cells
    edges = np.linspace(0.0, pitch, n_cells + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    cx, cy = np.meshgrid(centres, centres, indexing="ij")
    r = np.sqrt((cx - pitch / 2) ** 2 + (cy - pitch / 2) ** 2)

    mat_map = np.full((n_cells, n_cells), mat_ids[-1], dtype=int)
    for k in range(len(radii) - 1, -1, -1):
        mat_map[r <= radii[k]] = mat_ids[k]

    return Mesh2D(edges_x=edges, edges_y=edges, mat_map=mat_map)


__all__ = [
    "pwr_pin_2d",
]
