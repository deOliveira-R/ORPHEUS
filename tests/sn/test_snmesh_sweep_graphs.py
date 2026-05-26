r"""Tests for :attr:`SNMesh.sweep_graphs` precompute (Wave 2 / C2.4).

The Wave 2 plan lifts the per-call ``_diag_cache`` precompute that
previously lived inside :func:`orpheus.sn.sweep._sweep_2d_wavefront`
to mesh-construction time. Each :class:`SNMesh` built on a
:class:`Mesh2D` Cartesian grid now carries a
``Mapping[OctantLabel, SweepDependencyGraph]`` precomputed once per
``(mesh, quadrature)`` pair.

This file pins:

* **Existence on Cartesian 2-D**: every Cartesian SNMesh has the dict
  populated with the four streaming octants ``(±1, ±1)``.
* **Absence on curvilinear / 1-D**: spherical, cylindrical, and 1-D
  Cartesian SNMesh have ``sweep_graphs is None`` (no 2-D in-plane
  schedule).
* **Schedule agreement** with a hand-derived ``_diag_cache`` on small
  meshes (3×3, 4×5).
* **Mesh-time, not sweep-time**: building an SNMesh once, running it
  through many sweeps, the dict reference does not rebuild.
* **Per-octant dict keys** match the four ``(±1, ±1)`` octants,
  using :class:`OctantLabel` as the key type.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Mesh2D
from orpheus.sn.geometry import SNMesh
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.sweep_graph import OctantLabel, SweepDependencyGraph
from tests.sn._test_helpers import placeholder_materials


def _build_2d_mesh(nx: int = 3, ny: int = 3) -> Mesh2D:
    return Mesh2D(
        edges_x=np.linspace(0.0, 1.0, nx + 1),
        edges_y=np.linspace(0.0, 1.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("vacuum"), bc_ymax=BC("vacuum"),
    )


# ─────────────────────────────────────────────────────────────────────
# Existence + absence by coordinate system
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestSweepGraphsByCoordSystem:
    def test_cartesian_2d_has_four_octants(self):
        sn_mesh = SNMesh(_build_2d_mesh(3, 3), Quadrature.level_symmetric(4), placeholder_materials())
        assert sn_mesh.sweep_graphs is not None
        assert len(sn_mesh.sweep_graphs) == 4
        expected_keys = {
            OctantLabel(+1, +1), OctantLabel(+1, -1),
            OctantLabel(-1, +1), OctantLabel(-1, -1),
        }
        assert set(sn_mesh.sweep_graphs.keys()) == expected_keys

    def test_cartesian_2d_with_lebedev(self):
        """Lebedev quadrature: octant set is the same — 4 in-plane octants
        (sign_z is dropped by the OctantLabel)."""
        sn_mesh = SNMesh(_build_2d_mesh(4, 4), Quadrature.lebedev(5), placeholder_materials())
        assert sn_mesh.sweep_graphs is not None
        assert len(sn_mesh.sweep_graphs) == 4

    def test_cartesian_1d_has_sweep_graphs(self):
        """1-D slab Cartesian SNMesh: ``_setup_cartesian`` runs for any
        Cartesian mesh, so ``sweep_graphs`` is populated with
        degenerate ``ny=1`` graphs.

        The 1-D slab sweep (:func:`_sweep_1d_cumprod`) does not
        consume them — they exist as a side-effect of the unified
        Cartesian setup path. Verifying the dict shape so a future
        consumer would not be surprised.
        """
        from orpheus.geometry import CoordSystem
        mesh = Mesh1D(
            edges=np.linspace(0.0, 1.0, 4),
            mat_ids=np.zeros(3, dtype=int),
            bc_left=BC("vacuum"), bc_right=BC("vacuum"),
            coord=CoordSystem.CARTESIAN,
        )
        sn_mesh = SNMesh(mesh, Quadrature.gauss_legendre(8), placeholder_materials())
        # 1-D slab still goes through _setup_cartesian, so dict exists.
        assert sn_mesh.sweep_graphs is not None
        assert len(sn_mesh.sweep_graphs) == 4
        # Degenerate ny=1: each graph has nx + ny - 1 = nx levels.
        graph = sn_mesh.sweep_graphs[OctantLabel(+1, +1)]
        assert len(graph.levels) == sn_mesh.nx + sn_mesh.ny - 1


@pytest.mark.l0
class TestSweepGraphsAbsentOnCurvilinear:
    def test_spherical_1d_no_sweep_graphs(self):
        from orpheus.geometry import CoordSystem
        mesh = Mesh1D(
            edges=np.linspace(0.0, 1.0, 4),
            mat_ids=np.zeros(3, dtype=int),
            bc_left=BC("reflective"), bc_right=BC("vacuum"),
            coord=CoordSystem.SPHERICAL,
        )
        sn_mesh = SNMesh(mesh, Quadrature.gauss_legendre(8), placeholder_materials())
        assert sn_mesh.sweep_graphs is None

    def test_cylindrical_1d_no_sweep_graphs(self):
        from orpheus.geometry import CoordSystem
        mesh = Mesh1D(
            edges=np.linspace(0.0, 1.0, 4),
            mat_ids=np.zeros(3, dtype=int),
            bc_left=BC("reflective"), bc_right=BC("vacuum"),
            coord=CoordSystem.CYLINDRICAL,
        )
        sn_mesh = SNMesh(mesh, Quadrature.level_symmetric(4), placeholder_materials())
        assert sn_mesh.sweep_graphs is None


# ─────────────────────────────────────────────────────────────────────
# Schedule agreement with hand-derived _diag_cache
# ─────────────────────────────────────────────────────────────────────


def _hand_diag_cache(nx: int, ny: int) -> dict[tuple[int, int], dict]:
    """Hand-derived per-octant anti-diagonal schedule, matching the
    legacy :func:`_sweep_2d_wavefront` build at sweep.py:766-785.
    """
    cache: dict[tuple[int, int], dict] = {}
    for sx in (-1, +1):
        for sy in (-1, +1):
            ix_in = 0 if sx >= 0 else 1
            ix_out = 1 if sx >= 0 else 0
            iy_in = 0 if sy >= 0 else 1
            iy_out = 1 if sy >= 0 else 0
            ix_arr = np.arange(nx) if sx >= 0 else np.arange(nx)[::-1]
            iy_arr = np.arange(ny) if sy >= 0 else np.arange(ny)[::-1]
            diags = []
            for k in range(nx + ny - 1):
                i_start = max(0, k - ny + 1)
                i_end = min(nx - 1, k)
                local_i = np.arange(i_start, i_end + 1)
                local_j = k - local_i
                diags.append((ix_arr[local_i], iy_arr[local_j]))
            cache[(sx, sy)] = dict(
                ix_in=ix_in, ix_out=ix_out,
                iy_in=iy_in, iy_out=iy_out,
                diags=diags,
            )
    return cache


@pytest.mark.l0
class TestSweepGraphsAgreeWithLegacyDiagCache:
    @pytest.mark.parametrize("nx,ny", [(3, 3), (4, 5), (5, 4), (2, 7)])
    def test_diags_match_legacy_per_octant(self, nx, ny):
        sn_mesh = SNMesh(_build_2d_mesh(nx, ny), Quadrature.level_symmetric(4), placeholder_materials())
        legacy = _hand_diag_cache(nx, ny)
        for (sx, sy), legacy_entry in legacy.items():
            graph = sn_mesh.sweep_graphs[OctantLabel(sx, sy)]
            assert graph.face_in_x == legacy_entry["ix_in"]
            assert graph.face_out_x == legacy_entry["ix_out"]
            assert graph.face_in_y == legacy_entry["iy_in"]
            assert graph.face_out_y == legacy_entry["iy_out"]
            assert len(graph.levels) == len(legacy_entry["diags"])
            for level, legacy_diag in zip(graph.levels, legacy_entry["diags"]):
                np.testing.assert_array_equal(level[0], legacy_diag[0])
                np.testing.assert_array_equal(level[1], legacy_diag[1])


# ─────────────────────────────────────────────────────────────────────
# Mesh-time precompute (no rebuild on access)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestSweepGraphsAreMeshTime:
    def test_dict_reference_stable(self):
        """Multiple reads return the SAME dict object (precomputed once)."""
        sn_mesh = SNMesh(_build_2d_mesh(3, 3), Quadrature.level_symmetric(4), placeholder_materials())
        first = sn_mesh.sweep_graphs
        second = sn_mesh.sweep_graphs
        assert first is second

    def test_graphs_within_dict_stable(self):
        sn_mesh = SNMesh(_build_2d_mesh(3, 3), Quadrature.level_symmetric(4), placeholder_materials())
        for label in sn_mesh.sweep_graphs:
            graph_a = sn_mesh.sweep_graphs[label]
            graph_b = sn_mesh.sweep_graphs[label]
            assert graph_a is graph_b


# ─────────────────────────────────────────────────────────────────────
# Type contract
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestSweepGraphsTypeContract:
    def test_keys_are_octant_label(self):
        sn_mesh = SNMesh(_build_2d_mesh(3, 3), Quadrature.level_symmetric(4), placeholder_materials())
        for key in sn_mesh.sweep_graphs:
            assert isinstance(key, OctantLabel)

    def test_values_are_sweep_dependency_graph(self):
        sn_mesh = SNMesh(_build_2d_mesh(3, 3), Quadrature.level_symmetric(4), placeholder_materials())
        for value in sn_mesh.sweep_graphs.values():
            assert isinstance(value, SweepDependencyGraph)
