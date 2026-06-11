r"""End-to-end ``windowed ≡ full-field`` 2-D equivalence — the storage-B oracle.

Phase 5b replaced the full per-axis interior face cochain with a rolling
``_MovingFrontier`` window for BOTH the 2-D sweep (`_sweep_2d_wavefront`) and
its matvec twin (`MovingFrontierWindow.loss_action`). The retained
full-field paths — `orpheus.sn.sweep._sweep_full_field` and
`FullFieldWavefront.loss_action` (the d-generic spine) — are the VERIFICATION
ORACLES: the fuller view (the full interior face field, carried as a typed
``WavefrontFlux`` with whole-trace ι_*/ι* seed/absorb) that the moving frontier
is cross-checked against END-TO-END.

This is STRONGER than `test_sweep_graph_window_equivalence.py` (which pins only
the per-octant WALK at the graph level): it pins the whole ``transport_sweep``
/ matvec ORCHESTRATOR — the per-octant `boundary_flux.face_view` inflow/outflow
of the window vs the whole-trace `WavefrontFlux.seed`/`absorb` of the full
field — bit-for-bit. Both share the cell kernel
(`DiamondDifference.cell_kernel_batch`/`residual_kernel_batch`), so the MATH
cannot drift; only the storage walk + the boundary bookkeeping differ, which is
exactly what this pins. It also keeps the typed ``WavefrontFlux`` cochain (the
fuller-view type the carve had orphaned) genuinely exercised.

``foundation`` — a software invariant (window ≡ full-field), no theory label.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.loss_representation import FullFieldWavefront, MovingFrontierWindow
from orpheus.sn.operator import StreamingOperator
from orpheus.sn.sweep import _sweep_2d_wavefront, _sweep_full_field
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField

pytestmark = pytest.mark.foundation

# (nx, ny, level, ng, bc)
CASES = [
    (8, 8, 4, 2, "reflective"),
    (8, 8, 4, 2, "vacuum"),
    (12, 7, 6, 2, "reflective"),
    (5, 9, 4, 4, "reflective"),
]


def _build_mesh(nx, ny, lvl, ng, bc):
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 2.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_xmin=BC(bc), bc_xmax=BC(bc), bc_ymin=BC(bc), bc_ymax=BC(bc),
    )
    return SNMesh(mesh, Quadrature.level_symmetric(lvl), {0: get_mixture("A", f"{ng}g")})


def _random_sig_t(rng, ng, nx, ny):
    """Heterogeneous, strictly positive Σ_t (Cardinal Rule ≥2G + het)."""
    return rng.uniform(0.3, 3.0, size=(ng, nx, ny))


def _seed_random_inflow(rng, boundary):
    """Fill every boundary face slab with a random inflow trace (so the
    orchestrator's boundary handling is genuinely exercised)."""
    for face in boundary.layout.faces:
        fv = boundary.face_view(face)
        fv[...] = rng.uniform(0.0, 1.0, size=fv.shape)


@pytest.mark.parametrize("nx,ny,lvl,ng,bc", CASES)
def test_sweep_window_equals_full_field_end_to_end(nx, ny, lvl, ng, bc):
    """transport_sweep (windowed) ≡ _sweep_full_field (oracle), bit-for-bit:
    angular flux, scalar flux, AND the post-sweep boundary trace."""
    rng = np.random.default_rng(abs(hash((nx, ny, lvl, ng, bc))) % (2**32))
    sn_mesh = _build_mesh(nx, ny, lvl, ng, bc)
    N = sn_mesh.quad.N
    sig_t = _random_sig_t(rng, ng, nx, ny)
    Q = rng.uniform(0.0, 2.0, size=(N, ng, nx, ny))

    bf_win = BoundaryFlux.zeros_on(sn_mesh)
    _seed_random_inflow(rng, bf_win)
    bf_full = BoundaryFlux.from_mesh(bf_win.values.copy(), sn_mesh)

    ang_w, scal_w = _sweep_2d_wavefront(Q, sig_t, sn_mesh, bf_win)
    ang_f, scal_f = _sweep_full_field(Q, sig_t, sn_mesh, bf_full)

    np.testing.assert_array_equal(ang_w, ang_f, err_msg="angular flux")
    np.testing.assert_array_equal(scal_w, scal_f, err_msg="scalar flux")
    np.testing.assert_array_equal(
        bf_win.values, bf_full.values, err_msg="post-sweep boundary trace",
    )


@pytest.mark.parametrize("nx,ny,lvl,ng,bc", CASES)
def test_matvec_window_equals_full_field_end_to_end(nx, ny, lvl, ng, bc):
    """MovingFrontierWindow.loss_action (windowed) ≡ FullFieldWavefront.loss_action
    (oracle), bit-for-bit: bulk (L+C)ψ AND boundary-block residual.

    S6.3: both representations' ``loss_action`` return ``(L+C)ψ`` (the operator
    does the ``−C`` glue in :meth:`apply`); the window≡full comparison is
    unchanged — two storage policies over ONE walk."""
    rng = np.random.default_rng((abs(hash((nx, ny, lvl, ng, bc))) % (2**32)) ^ 9)
    sn_mesh = _build_mesh(nx, ny, lvl, ng, bc)
    N = sn_mesh.quad.N
    sig_t = _random_sig_t(rng, ng, nx, ny)
    L = StreamingOperator(sn_mesh, sig_t)

    state = TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh)
    state.bulk.values[...] = rng.uniform(-1.0, 1.0, size=state.bulk.values.shape)
    _seed_random_inflow(rng, state.boundary)

    out_win = MovingFrontierWindow(sn_mesh).loss_action(L, state)
    out_full = FullFieldWavefront(sn_mesh).loss_action(L, state)

    np.testing.assert_array_equal(
        out_win.bulk.values, out_full.bulk.values, err_msg="bulk residual",
    )
    np.testing.assert_array_equal(
        out_win.boundary.values, out_full.boundary.values,
        err_msg="boundary-block residual",
    )
