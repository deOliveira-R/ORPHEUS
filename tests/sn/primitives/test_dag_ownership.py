r"""DAG ownership + the per-shape graph cache (S6.4(c), #222).

The per-octant :class:`SweepDependencyGraph` family is OWNED by the
``_DAGWavefront`` representation family through the per-shape cache
:meth:`SweepDependencyGraph.for_shape` — NOT by the mesh.  The mesh is pure
geometry; DAG-free representations (``CumprodScan``, ``ScanMarch``) never
mention the substrate; the historical curvilinear
``mesh.sweep_graphs = None`` slot (an illegal state) is gone.

History: this file is the migrated successor of
``test_snmesh_sweep_graphs.py`` (Wave 2 / C2.4 pinned the graphs as a
mesh-construction precompute; S6.4(c) relocated ownership — retirement =
test migration).  The graph-CONTENT pins (the four in-plane octants, the
genuine d=1 chain, the hand-derived legacy ``_diag_cache`` golden, the
type contract) carry over verbatim, now driven through ``for_shape``; the
"absent on curvilinear" pins became the structural
mesh-does-not-expose-the-DAG assertions.

``-O``-safe (vv Mode 8): ``pytest.fail`` / ``np.testing`` only — the
migration also retired this file's pre-existing bare-``assert`` gap.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.loss_representation import (
    CumprodScan,
    FullFieldWavefront,
    MovingFrontierWindow,
    ScanMarch,
)
from orpheus.sn.sweep_graph import OctantLabel, SweepDependencyGraph
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation


def _build_2d_mesh(nx: int = 3, ny: int = 3) -> Mesh2D:
    return Mesh2D(
        edges_x=np.linspace(0.0, 1.0, nx + 1),
        edges_y=np.linspace(0.0, 1.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("vacuum"), bc_ymax=BC("vacuum"),
    )


def _build_mesh(coord: str) -> SNMesh:
    """The standard per-coordinate-system SNMesh builder."""
    if coord == "cart2d":
        return SNMesh(
            _build_2d_mesh(3, 4), Quadrature.level_symmetric(4),
            placeholder_materials(),
        )
    coord_sys = {
        "slab": CoordSystem.CARTESIAN,
        "sphere": CoordSystem.SPHERICAL,
        "cyl": CoordSystem.CYLINDRICAL,
    }[coord]
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 4),
        mat_ids=np.zeros(3, dtype=int),
        bc_left=BC("reflective"), bc_right=BC("vacuum"),
        coord=coord_sys,
    )
    quad = (
        Quadrature.gauss_legendre(8)
        if coord in ("slab", "sphere") else Quadrature.level_symmetric(4)
    )
    return SNMesh(mesh, quad, placeholder_materials())


# ─────────────────────────────────────────────────────────────────────
# Ownership — the mesh is pure geometry (illegal states unrepresentable)
# ─────────────────────────────────────────────────────────────────────


class TestDagOwnership:
    @pytest.mark.parametrize("coord", ["slab", "cart2d", "sphere", "cyl"])
    def test_mesh_does_not_expose_sweep_graphs(self, coord):
        """[L0 structural] The DAG is OWNED by _DAGWavefront, NOT the mesh.

        Post-S6.4(c) the mesh is pure geometry — a mesh attribute would be
        the illegal-state None-slot the curvilinear arms used to carry.
        """
        sn = _build_mesh(coord)
        if hasattr(sn, "sweep_graphs"):
            pytest.fail(
                f"{coord}: SNMesh still exposes `sweep_graphs` — the DAG "
                "must be owned by the _DAGWavefront family (cached per "
                "mesh shape), not the mesh (the curvilinear None-slot was "
                "the illegal-state smell S6.4(c) closed)."
            )

    def test_geometry_module_is_dag_free(self):
        """[L0 structural] geometry.py no longer mentions the DAG substrate.

        Catches both the curvilinear ``= None`` slots and the Cartesian
        build site — any surviving mention means the ownership move is
        incomplete.
        """
        import inspect

        from orpheus.sn import geometry

        if "sweep_graphs" in inspect.getsource(geometry):
            pytest.fail(
                "geometry.py still mentions `sweep_graphs` — the DAG is "
                "owned by _DAGWavefront via SweepDependencyGraph.for_shape; "
                "the mesh must stay pure geometry."
            )

    def test_dag_free_representations_never_reference_the_substrate(self):
        """[L0 structural] CumprodScan / ScanMarch never touch the DAG.

        DAG-free representations are constructed on curvilinear / 1-D /
        row-march meshes where no anti-hyperplane DAG is consumed; a
        ``sweep_graphs`` reference would re-introduce the substrate
        coupling the carve removed.
        """
        import inspect

        for rep_cls in (CumprodScan, ScanMarch):
            for name in (
                "sweep", "loss_action", "loss_action_transpose",
                "_sweep_interior", "_loss_action_interior",
            ):
                method = getattr(rep_cls, name, None)
                if method is None:
                    continue
                if "sweep_graphs" in inspect.getsource(method):
                    pytest.fail(
                        f"{rep_cls.__name__}.{name} references "
                        "`sweep_graphs` — a DAG-free representation must "
                        "never mention the DAG substrate."
                    )

    @pytest.mark.parametrize("rep_cls", [MovingFrontierWindow, FullFieldWavefront])
    def test_dag_family_exposes_the_accessor(self, rep_cls):
        """[L0 structural] The _DAGWavefront family OWNS the accessor: the
        property resolves to the per-shape graphs on a real mesh."""
        sn = _build_mesh("cart2d")
        graphs = rep_cls(sn).sweep_graphs
        if set(graphs.keys()) != {
            OctantLabel((+1, +1)), OctantLabel((+1, -1)),
            OctantLabel((-1, +1)), OctantLabel((-1, -1)),
        }:
            pytest.fail(
                f"{rep_cls.__name__}.sweep_graphs did not resolve the four "
                f"in-plane octants: {set(graphs.keys())}"
            )


# ─────────────────────────────────────────────────────────────────────
# The per-shape cache — relocation, not rebuild (byte-identical graphs)
# ─────────────────────────────────────────────────────────────────────


class TestPerShapeCache:
    def test_same_shape_meshes_share_the_graphs(self):
        """[L0] Two meshes of the SAME spatial shape share the SAME cached
        family (the mesh-time precompute contract carries over: repeated
        access does not rebuild)."""
        rep_a = MovingFrontierWindow(_build_mesh("cart2d"))
        rep_b = MovingFrontierWindow(_build_mesh("cart2d"))
        if rep_a.sweep_graphs is not rep_b.sweep_graphs:
            pytest.fail(
                "same-shape meshes did not share the cached DAG family — "
                "the per-shape cache is not keying on shape."
            )
        for label in rep_a.sweep_graphs:
            if rep_a.sweep_graphs[label] is not rep_b.sweep_graphs[label]:
                pytest.fail(f"octant {label}: graphs not shared")

    def test_cache_equals_fresh_from_cartesian_build(self):
        """[L0] The cached graph equals a freshly-built one — the ownership
        move is a relocation, not a rebuild (byte-identical level structure,
        face offsets, octant key set)."""
        shape = (3, 4)
        cached = SweepDependencyGraph.for_shape(shape)
        for label, graph in cached.items():
            fresh = SweepDependencyGraph.from_cartesian(shape, label=label)
            np.testing.assert_array_equal(
                len(graph.levels), len(fresh.levels),
                err_msg=f"{label}: level count",
            )
            for k, (lvl_c, lvl_f) in enumerate(zip(graph.levels, fresh.levels)):
                for a, (idx_c, idx_f) in enumerate(zip(lvl_c, lvl_f)):
                    np.testing.assert_array_equal(
                        idx_c, idx_f,
                        err_msg=f"{label}: level {k} axis {a} cell indices",
                    )
            np.testing.assert_array_equal(
                graph.face_in, fresh.face_in, err_msg=f"{label}: face_in",
            )
            np.testing.assert_array_equal(
                graph.face_out, fresh.face_out, err_msg=f"{label}: face_out",
            )

    def test_d1_shape_yields_the_chain_graphs(self):
        """[L0] ``for_shape((nx,))`` is the genuine d=1 family: 2 streaming
        octants (``±x``), 1-tuple sign labels, the pure-chain DAG (one cell
        per level) — phantom-axis elimination (C3.3) preserved through the
        ownership move."""
        nx = 3
        graphs = SweepDependencyGraph.for_shape((nx,))
        if set(graphs.keys()) != {OctantLabel((-1,)), OctantLabel((+1,))}:
            pytest.fail(
                f"d=1 octant labels {set(graphs.keys())} != "
                "{OctantLabel((-1,)), OctantLabel((+1,))}"
            )
        graph = graphs[OctantLabel((+1,))]
        np.testing.assert_array_equal(
            len(graph.levels), nx,
            err_msg="d=1 chain must have nx levels (one cell per level)",
        )
        for k, level in enumerate(graph.levels):
            np.testing.assert_array_equal(
                len(level), 1,
                err_msg=f"level {k}: d=1 level must be a 1-axis tuple",
            )
            np.testing.assert_array_equal(
                len(level[0]), 1,
                err_msg=f"level {k}: chain has one cell per level",
            )


# ─────────────────────────────────────────────────────────────────────
# Schedule agreement with the hand-derived legacy _diag_cache (golden)
# ─────────────────────────────────────────────────────────────────────


def _hand_diag_cache(nx: int, ny: int) -> dict[tuple[int, int], dict]:
    """Hand-derived per-octant anti-diagonal schedule, matching the
    legacy ``_sweep_jacobi`` per-call build (pre-Wave-2)."""
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


class TestGraphsAgreeWithLegacyDiagCache:
    @pytest.mark.parametrize("nx,ny", [(3, 3), (4, 5), (5, 4), (2, 7)])
    def test_diags_match_legacy_per_octant(self, nx, ny):
        graphs = SweepDependencyGraph.for_shape((nx, ny))
        legacy = _hand_diag_cache(nx, ny)
        for (sx, sy), legacy_entry in legacy.items():
            graph = graphs[OctantLabel((sx, sy))]
            np.testing.assert_array_equal(
                [graph.face_in_x, graph.face_out_x,
                 graph.face_in_y, graph.face_out_y],
                [legacy_entry["ix_in"], legacy_entry["ix_out"],
                 legacy_entry["iy_in"], legacy_entry["iy_out"]],
                err_msg=f"octant ({sx},{sy}): face offsets",
            )
            np.testing.assert_array_equal(
                len(graph.levels), len(legacy_entry["diags"]),
                err_msg=f"octant ({sx},{sy}): level count",
            )
            for level, legacy_diag in zip(graph.levels, legacy_entry["diags"]):
                np.testing.assert_array_equal(level[0], legacy_diag[0])
                np.testing.assert_array_equal(level[1], legacy_diag[1])


# ─────────────────────────────────────────────────────────────────────
# Type contract
# ─────────────────────────────────────────────────────────────────────


class TestForShapeTypeContract:
    def test_keys_are_octant_label(self):
        for key in SweepDependencyGraph.for_shape((3, 3)):
            if not isinstance(key, OctantLabel):
                pytest.fail(f"for_shape key {key!r} is not an OctantLabel")

    def test_values_are_sweep_dependency_graph(self):
        for value in SweepDependencyGraph.for_shape((3, 3)).values():
            if not isinstance(value, SweepDependencyGraph):
                pytest.fail(
                    f"for_shape value {type(value).__name__} is not a "
                    "SweepDependencyGraph"
                )

    def test_lebedev_quadrature_same_octant_set(self):
        """Quadrature-independence: the DAG family depends only on the
        shape — a Lebedev-quadrature mesh consumes the SAME 4 in-plane
        octant graphs (sign_z is dropped by the in-plane projection)."""
        sn = SNMesh(
            _build_2d_mesh(4, 4), Quadrature.lebedev(5), placeholder_materials(),
        )
        graphs = MovingFrontierWindow(sn).sweep_graphs
        np.testing.assert_array_equal(
            len(graphs), 4, err_msg="lebedev mesh: 4 in-plane octants",
        )
