"""Mesh construction factories — 2-D Cartesian only.

The 1-D path is :class:`~orpheus.mesh.mesher.Mesher` over a
:class:`~orpheus.geometry.structured_geometry.StructuredGeometry`, with
:meth:`StructuredGeometry.wigner_seitz_pin_cell` and
:meth:`StructuredGeometry.pwr_slab_half_cell` for the conventional PWR
pin-cell shapes. What survives here is :func:`pwr_pin_2d`, the 2-D Cartesian
factory: there is no 2-D :class:`StructuredGeometry` yet.
"""

from __future__ import annotations

import numpy as np

from orpheus.geometry.boundary import BC, BoundaryTraceLaw
from orpheus.mesh.structured import Mesh2D


# ── 2-D Cartesian PWR pin factory ────────────────────────────────────

def pwr_pin_2d(
    radii: list[float] | None = None,
    mat_ids: list[int] | None = None,
    pitch: float = 3.6,
    n_cells: int = 10,
    *,
    law: "BC | BoundaryTraceLaw",
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
    law : BC or BoundaryTraceLaw
        The law on all four faces of the cell; required, because a factory
        does not choose a physical boundary on the caller's behalf (a lattice
        cell is reflective or periodic, an isolated one vacuum).
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

    return Mesh2D(
        edges, edges, mat_map,
        face_laws={face: law for face in ("xmin", "xmax", "ymin", "ymax")},
    )


__all__ = [
    "pwr_pin_2d",
]
