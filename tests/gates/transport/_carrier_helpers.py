r"""Shared carrier fixture — a GENUINE unit-width one-cell ``Mesh1D`` carrier.

Until the CS4c coda (2026-09-08) the tree carried a factory —
``MaterialMesh.from_materials`` — that fabricated a mesh-less one-cell
carrier for the infinite homogeneous medium: ``[0, 1]`` edges, a node at
``0.5`` and a Cartesian chart that `[M]` nothing consumed. The homogeneous
problem now poses on its own space (``HomogeneousProblem.space``) and
builds no carrier, and the factory retired with it. Tests that used the
factory as a CONVENIENCE fixture (a facade over one cell, a two-material
carrier with a spectator entry, an axis-built ``bulk_space``) migrate here.

`[M]` like-for-like on the pre-retirement tree (2026-09-08; four
declarations: A-2g, A-4g, a spectator pair, a spectator pair carrying an
energy grid): ``volumes`` ``[1.0]``, ``mat_map``, ``spatial_shape``,
``ndim``, ``nx``, ``ng``, ``coord``, ``axes[0].edges``, the declaration's
keys, ``cells_by_material`` (spectator retention preserved —
``{0: ([0],), 1: ([],)}``), ``bulk_space`` (``==`` the fabricated
carrier's AND ``==`` the problem's pose), and EVERY facade array
(``total_cross_section``, ``fission_production``, ``sig_s_legendre``,
``n2n_matrix``, ``chi_per_material``, ``fission_production_per_material``,
the three cross-section fields' values, the fields' space ``==`` the pose)
are bit-identical to the fabricated carrier's. The ONE difference is
``volume_measure.nodes.shape`` — ``(1,)`` here (``Mesh1D``'s measure) vs
the fabricated ``(1, 1)`` (an axis-native coordinate tuple) — which no
migrated site reads; the values (``[0.5]``, weight ``[1.0]``) agree.
"""

from __future__ import annotations

from collections.abc import Mapping

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry import BC, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.transport.mesh.material_mesh import MaterialMesh


def unit_cell_carrier(materials: Mapping[int, Mixture]) -> MaterialMesh:
    """A genuine unit-width one-cell Cartesian ``Mesh1D`` carrier (reflective
    faces, material id ``0`` in the cell) over ``materials`` — which may
    declare SPECTATOR ids the cell does not reference."""
    mesh = Mesher(StructuredGeometry.from_homogeneous(1.0, BC("reflective"))).partition(
        CellsByCount.uniform_width(1),
    ).mesh
    return MaterialMesh(mesh, materials)
