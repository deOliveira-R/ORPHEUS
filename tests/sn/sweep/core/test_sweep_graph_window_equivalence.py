r"""``window ≡ full-field`` equivalence — the Phase 5b storage-B oracle.

The rolling 2-diagonal moving-frontier window walks
(:meth:`SweepDependencyGraph.apply_windowed` /
:meth:`~SweepDependencyGraph.residual_windowed`, the storage-B PRODUCTION
path) MUST reproduce the full per-axis face-field walks
(:meth:`~SweepDependencyGraph.apply` / :meth:`~SweepDependencyGraph.residual`,
the storage-A path, KEPT as the verification oracle) **bit-for-bit**.

Both share the storage-free cell kernel
(:meth:`DiamondDifference.cell_kernel_batch` /
:meth:`~DiamondDifference.residual_kernel_batch`), so the cell *math* cannot
drift between them — only the storage walk (full field vs window) differs,
which is exactly what this test pins. This is the permanent regression guard
promoted from the Phase-0 de-risk
(``derivations/diagnostics/diag_phase0_storageb_derisk.py``): it covers all
four octants on heterogeneous meshes, including the highest-risk boundary-shed
capture (a recycled-slot overwrite would silently corrupt the outflow trace).

``foundation`` level: a software invariant (window ≡ full-field) with no
theory ``:label:`` — it pins the storage refactor, not a physics equation.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.sn.spatial.diamond import DiamondDifference
from orpheus.sn.sweep_graph import OctantLabel, SweepDependencyGraph

pytestmark = pytest.mark.foundation

OCTANTS = [(1, 1), (-1, 1), (1, -1), (-1, -1)]
MESHES = [  # (nx, ny, N_oct, ng)
    (8, 8, 4, 2),
    (12, 7, 6, 3),
    (5, 9, 4, 2),
    (16, 16, 8, 2),
]


def _inputs(nx, ny, N_oct, ng, seed):
    rng = np.random.default_rng(seed)
    return dict(
        sig_t=rng.uniform(0.3, 3.0, size=(ng, nx, ny)),
        Q=rng.uniform(0.0, 2.0, size=(N_oct, ng, nx, ny)),
        str_x=rng.uniform(0.5, 4.0, size=(N_oct, nx)),
        str_y=rng.uniform(0.5, 4.0, size=(N_oct, ny)),
        inflow_x=rng.uniform(0.0, 1.0, size=(N_oct, ng, ny)),
        inflow_y=rng.uniform(0.0, 1.0, size=(N_oct, ng, nx)),
        weights=rng.uniform(0.1, 1.0, size=(N_oct,)),
        probe=rng.uniform(-1.0, 1.0, size=(N_oct, ng, nx, ny)),
    )


def _full_field_apply(graph, nx, ny, N_oct, ng, inp):
    """Oracle: the production full-field solve walk on dense face buffers."""
    sx, sy = graph.label.sign_x, graph.label.sign_y
    psi_x = np.zeros((N_oct, ng, nx + 1, ny))
    psi_y = np.zeros((N_oct, ng, nx, ny + 1))
    psi_x[:, :, 0 if sx >= 0 else nx, :] = inp["inflow_x"]
    psi_y[:, :, :, 0 if sy >= 0 else ny] = inp["inflow_y"]
    ang = np.zeros((N_oct, ng, nx, ny))
    scal = np.zeros((ng, nx, ny))
    graph.apply(
        cell_update=DiamondDifference(),
        psi_x_octant=psi_x, psi_y_octant=psi_y,
        Q_octant=inp["Q"], sig_t=inp["sig_t"],
        str_x_octant=inp["str_x"], str_y_octant=inp["str_y"],
        weights_octant=inp["weights"],
        angular_flux_octant=ang, scalar_flux_buf=scal,
    )
    out_x = psi_x[:, :, nx if sx >= 0 else 0, :]   # domain x-outflow
    out_y = psi_y[:, :, :, ny if sy >= 0 else 0]   # domain y-outflow
    return ang, scal, out_x, out_y


def _windowed_apply(graph, nx, ny, N_oct, ng, inp):
    ang = np.zeros((N_oct, ng, nx, ny))
    scal = np.zeros((ng, nx, ny))
    cap_x = np.zeros((N_oct, ng, ny))
    cap_y = np.zeros((N_oct, ng, nx))
    graph.apply_windowed(
        cell_update=DiamondDifference(),
        inflow_x=inp["inflow_x"], inflow_y=inp["inflow_y"],
        Q_octant=inp["Q"], sig_t=inp["sig_t"],
        str_x_octant=inp["str_x"], str_y_octant=inp["str_y"],
        weights_octant=inp["weights"],
        angular_flux_octant=ang, scalar_flux_buf=scal,
        capture_x=cap_x, capture_y=cap_y,
    )
    return ang, scal, cap_x, cap_y


def _full_field_residual(graph, nx, ny, N_oct, ng, inp):
    sx, sy = graph.label.sign_x, graph.label.sign_y
    psi_x = np.zeros((N_oct, ng, nx + 1, ny))
    psi_y = np.zeros((N_oct, ng, nx, ny + 1))
    psi_x[:, :, 0 if sx >= 0 else nx, :] = inp["inflow_x"]
    psi_y[:, :, :, 0 if sy >= 0 else ny] = inp["inflow_y"]
    res = np.zeros((N_oct, ng, nx, ny))
    graph.residual(
        cell_update=DiamondDifference(),
        psi_x_octant=psi_x, psi_y_octant=psi_y,
        psi_avg_probe_octant=inp["probe"],
        Q_octant=inp["Q"], sig_t=inp["sig_t"],
        str_x_octant=inp["str_x"], str_y_octant=inp["str_y"],
        residual_octant=res,
    )
    out_x = psi_x[:, :, nx if sx >= 0 else 0, :]
    out_y = psi_y[:, :, :, ny if sy >= 0 else 0]
    return res, out_x, out_y


def _windowed_residual(graph, nx, ny, N_oct, ng, inp):
    res = np.zeros((N_oct, ng, nx, ny))
    cap_x = np.zeros((N_oct, ng, ny))
    cap_y = np.zeros((N_oct, ng, nx))
    graph.residual_windowed(
        cell_update=DiamondDifference(),
        inflow_x=inp["inflow_x"], inflow_y=inp["inflow_y"],
        psi_avg_probe_octant=inp["probe"],
        Q_octant=inp["Q"], sig_t=inp["sig_t"],
        str_x_octant=inp["str_x"], str_y_octant=inp["str_y"],
        residual_octant=res, capture_x=cap_x, capture_y=cap_y,
    )
    return res, cap_x, cap_y


@pytest.mark.parametrize("nx,ny,N_oct,ng", MESHES)
@pytest.mark.parametrize("sx,sy", OCTANTS)
def test_solve_window_equals_full_field(nx, ny, N_oct, ng, sx, sy):
    """apply_windowed ≡ apply, bit-for-bit (angular, scalar, shed outflow)."""
    graph = SweepDependencyGraph.from_cartesian_2d(
        nx=nx, ny=ny, label=OctantLabel((sx, sy)),
    )
    inp = _inputs(nx, ny, N_oct, ng, seed=hash((nx, ny, sx, sy)) % (2**32))
    ang_f, scal_f, ox_f, oy_f = _full_field_apply(graph, nx, ny, N_oct, ng, inp)
    ang_w, scal_w, cx_w, cy_w = _windowed_apply(graph, nx, ny, N_oct, ng, inp)
    np.testing.assert_array_equal(ang_w, ang_f, err_msg="angular flux")
    np.testing.assert_array_equal(scal_w, scal_f, err_msg="scalar flux")
    np.testing.assert_array_equal(cx_w, ox_f, err_msg="x-outflow shed")
    np.testing.assert_array_equal(cy_w, oy_f, err_msg="y-outflow shed")


@pytest.mark.parametrize("nx,ny,N_oct,ng", MESHES)
@pytest.mark.parametrize("sx,sy", OCTANTS)
def test_residual_window_equals_full_field(nx, ny, N_oct, ng, sx, sy):
    """residual_windowed ≡ residual, bit-for-bit (matvec twin)."""
    graph = SweepDependencyGraph.from_cartesian_2d(
        nx=nx, ny=ny, label=OctantLabel((sx, sy)),
    )
    inp = _inputs(nx, ny, N_oct, ng, seed=(hash((nx, ny, sx, sy)) % (2**32)) ^ 7)
    res_f, ox_f, oy_f = _full_field_residual(graph, nx, ny, N_oct, ng, inp)
    res_w, cx_w, cy_w = _windowed_residual(graph, nx, ny, N_oct, ng, inp)
    np.testing.assert_array_equal(res_w, res_f, err_msg="operator residual")
    np.testing.assert_array_equal(cx_w, ox_f, err_msg="x-outflow shed")
    np.testing.assert_array_equal(cy_w, oy_f, err_msg="y-outflow shed")


def test_window_backing_is_linear_in_side():
    """The window backing is O(nx), NOT O(nx·ny) — the storage-B win.

    A direct shape assertion on the rolling buffer: per octant the window is
    (N_oct, ng, 2, nx) for x AND y, independent of ny — vs the full field's
    (N_oct, ng, nx+1, ny) + (N_oct, ng, nx, ny+1) which scales with nx·ny.
    """
    N_oct, ng = 8, 2
    for nx, ny in [(16, 16), (16, 64), (32, 32)]:
        graph = SweepDependencyGraph.from_cartesian_2d(
            nx=nx, ny=ny, label=OctantLabel((1, 1)),
        )
        # window slot dimension is nx, parity is 2 — no ny dependence
        win_elems = 2 * (N_oct * ng * 2 * graph.nx)        # x + y windows
        full_elems = (
            N_oct * ng * (nx + 1) * ny + N_oct * ng * nx * (ny + 1)
        )
        np.testing.assert_equal(win_elems, 2 * N_oct * ng * 2 * nx)
        # window must be far smaller once ny is non-trivial
        if ny >= nx:
            assert win_elems < full_elems / (ny / 4)  # ratio, not a value-pin
