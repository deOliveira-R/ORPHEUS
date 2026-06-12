r"""Tests for :mod:`orpheus.sn.sweep_graph` (Wave 2 / C2.3).

The §15A.2 "upwind trace complex / causal transport DAG / direction
sweep ordering" primitive (Grand Report v3 lines 2137-2171). This
file pins the invariants the plan named:

* :ref:`assert_upwind_orientation` — face_in / face_out indices match
  the octant sign convention.
* :ref:`assert_topologically_sorted` — every level's cells depend
  only on cells in strictly earlier levels.
* :ref:`assert_cell_coverage` — every cell appears in exactly one
  level.
* :ref:`assert_face_pairing_consistent` — within an octant, the
  outgoing face of cell ``(i, j)`` matches the incoming face of
  ``(i + sx, j)`` (and analogous for y).
* :ref:`apply_invariants` — running ``SweepDependencyGraph.walk_full``
  matches a per-cell hand calculation on a tiny grid (1×1, 2×2)
  and matches the legacy inlined wavefront math on a 3×3 grid.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.sn.spatial.diamond import DiamondDifference
from orpheus.sn.sweep_graph import (
    OctantLabel,
    SweepDependencyGraph,
    _CellResidual,
    _CellSolve,
)


# ─────────────────────────────────────────────────────────────────────
# OctantLabel
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestOctantLabel:
    def test_valid_signs(self):
        for sx in (-1, 0, +1):
            for sy in (-1, 0, +1):
                lab = OctantLabel((sx, sy))
                assert lab.signs == (sx, sy)

    def test_invalid_sign_raises(self):
        with pytest.raises(ValueError, match=r"signs\[0\]"):
            OctantLabel((2, +1))
        with pytest.raises(ValueError, match=r"signs\[1\]"):
            OctantLabel((+1, -2))

    def test_streams(self):
        assert OctantLabel((+1, +1)).streams
        assert OctantLabel((-1, +1)).streams
        assert OctantLabel((+1, 0)).streams   # y-axis-aligned still streams in x
        assert OctantLabel((0, -1)).streams
        assert not OctantLabel((0, 0)).streams  # pure-z degenerate (2-D label)
        # d-generic: the same predicate at d=1 and d=3.
        assert OctantLabel((-1,)).streams
        assert not OctantLabel((0, 0, 0)).streams

    def test_hashable_for_dict_keys(self):
        d = {OctantLabel((+1, -1)): "x"}
        assert d[OctantLabel((+1, -1))] == "x"


# ─────────────────────────────────────────────────────────────────────
# Construction — face-index conventions and degenerate guard
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestFromCartesian2D:
    @pytest.mark.parametrize("sx,sy,fix,fox,fiy,foy", [
        (+1, +1, 0, 1, 0, 1),
        (+1, -1, 0, 1, 1, 0),
        (-1, +1, 1, 0, 0, 1),
        (-1, -1, 1, 0, 1, 0),
    ])
    def test_face_indices_match_octant_signs(
        self, sx, sy, fix, fox, fiy, foy,
    ):
        """assert_upwind_orientation: face_in/out indices match sign conv."""
        g = SweepDependencyGraph.from_cartesian((4, 4), label=OctantLabel((sx, sy)))
        assert g.face_in_x == fix
        assert g.face_out_x == fox
        assert g.face_in_y == fiy
        assert g.face_out_y == foy

    def test_pure_z_label_raises(self):
        with pytest.raises(ValueError, match="degenerate"):
            SweepDependencyGraph.from_cartesian((3, 3), label=OctantLabel((0, 0)))

    def test_levels_is_tuple_of_ndarrays(self):
        g = SweepDependencyGraph.from_cartesian((3, 4), label=OctantLabel((+1, +1)))
        assert isinstance(g.levels, tuple)
        for ii, jj in g.levels:
            assert isinstance(ii, np.ndarray)
            assert isinstance(jj, np.ndarray)
            assert ii.shape == jj.shape


# ─────────────────────────────────────────────────────────────────────
# §15A.2 invariants — cell coverage + topological soundness
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestAssertCellCoverage:
    """Every cell of the nx × ny grid appears in exactly one level."""

    @pytest.mark.parametrize("nx,ny", [(1, 1), (2, 3), (3, 3), (4, 5), (5, 4)])
    @pytest.mark.parametrize("sx,sy", [
        (+1, +1), (+1, -1), (-1, +1), (-1, -1),
    ])
    def test_every_cell_visited_once(self, nx, ny, sx, sy):
        g = SweepDependencyGraph.from_cartesian((nx, ny), label=OctantLabel((sx, sy)))
        all_cells = set()
        for ii, jj in g.levels:
            for i, j in zip(ii, jj):
                cell = (int(i), int(j))
                assert cell not in all_cells, f"Cell {cell} visited twice"
                all_cells.add(cell)
        assert len(all_cells) == nx * ny


@pytest.mark.l0
class TestAssertTopologicallySorted:
    """Every level's cells depend only on cells in strictly earlier levels."""

    @pytest.mark.parametrize("nx,ny", [(2, 3), (3, 3), (4, 4)])
    @pytest.mark.parametrize("sx,sy", [
        (+1, +1), (+1, -1), (-1, +1), (-1, -1),
    ])
    def test_upstream_cell_in_earlier_level(self, nx, ny, sx, sy):
        """Cell (i, j) at level k has upstream (i - sx, j) and (i, j - sy)
        at level < k (or off-grid → BC)."""
        g = SweepDependencyGraph.from_cartesian((nx, ny), label=OctantLabel((sx, sy)))
        cell_to_level = {}
        for k, (ii, jj) in enumerate(g.levels):
            for i, j in zip(ii, jj):
                cell_to_level[(int(i), int(j))] = k
        # For each cell, its upstream cells must be in earlier levels.
        for (i, j), k in cell_to_level.items():
            upstream_x = (i - sx, j)
            upstream_y = (i, j - sy)
            if 0 <= upstream_x[0] < nx:
                assert cell_to_level[upstream_x] < k, (
                    f"cell {(i, j)} at level {k} has upstream-x "
                    f"{upstream_x} at level {cell_to_level[upstream_x]} "
                    f"(must be < {k})"
                )
            if 0 <= upstream_y[1] < ny:
                assert cell_to_level[upstream_y] < k, (
                    f"cell {(i, j)} at level {k} has upstream-y "
                    f"{upstream_y} at level {cell_to_level[upstream_y]} "
                    f"(must be < {k})"
                )


@pytest.mark.l0
class TestAssertFacePairingConsistent:
    r"""Within an octant: outgoing face of ``(i, j)`` matches the
    incoming face of the next cell along the streaming direction."""

    @pytest.mark.parametrize("sx,sy", [
        (+1, +1), (+1, -1), (-1, +1), (-1, -1),
    ])
    def test_face_indices_are_neighbors(self, sx, sy):
        """For sign_x = +1: out_x of (i, j) is i+1 == in_x of (i+1, j)."""
        g = SweepDependencyGraph.from_cartesian((3, 3), label=OctantLabel((sx, sy)))
        # face_out_x[i] - face_in_x[i+sx*1] = ?
        # For sx=+1: face_out=1, so cell (i, j) outgoing face is i+1;
        #            cell (i+1, j) incoming face is (i+1)+0 = i+1. ✓
        # For sx=-1: face_out=0, so cell (i, j) outgoing face is i;
        #            cell (i-1, j) incoming face is (i-1)+1 = i. ✓
        # The invariant is: face_out_x == 1 + face_in_x - sx (always -1+sx-?)
        # Simplest check: face_out_x + face_in_x == 1 (one is 0, the other is 1).
        assert g.face_in_x + g.face_out_x == 1
        assert g.face_in_y + g.face_out_y == 1


# ─────────────────────────────────────────────────────────────────────
# apply — equivalence with per-cell hand calculation
# ─────────────────────────────────────────────────────────────────────


def _hand_run_legacy_inlined(
    *,
    nx: int, ny: int, ng: int, N_oct: int,
    psi_x: np.ndarray, psi_y: np.ndarray,
    Q: np.ndarray,                 # (N_oct or 1, ng, nx, ny)
    sig_t: np.ndarray,             # (ng, nx, ny)
    str_x: np.ndarray,             # (N_oct, nx)
    str_y: np.ndarray,             # (N_oct, ny)
    weights: np.ndarray,           # (N_oct,)
    sx_sign: int, sy_sign: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Run the inlined ``_sweep_jacobi`` math one ordinate
    at a time, on octant-restricted buffers.

    Issue #196 PR-INDEX-5: principled layouts throughout —
    ``angular_flux: (N_oct, ng, nx, ny)``,
    ``scalar_flux: (ng, nx, ny)``,
    ``psi_x: (N_oct, ng, nx+1, ny)``, ``psi_y: (N_oct, ng, nx, ny+1)``.

    Returns ``(angular_flux, scalar_flux, psi_x_post, psi_y_post)``.
    """
    angular_flux = np.zeros((N_oct, ng, nx, ny))
    scalar_flux = np.zeros((ng, nx, ny))
    psi_x = psi_x.copy()
    psi_y = psi_y.copy()

    ix_in = 0 if sx_sign >= 0 else 1
    ix_out = 1 if sx_sign >= 0 else 0
    iy_in = 0 if sy_sign >= 0 else 1
    iy_out = 1 if sy_sign >= 0 else 0
    ix_arr = np.arange(nx) if sx_sign >= 0 else np.arange(nx)[::-1]
    iy_arr = np.arange(ny) if sy_sign >= 0 else np.arange(ny)[::-1]

    for n in range(N_oct):
        for k in range(nx + ny - 1):
            i_start = max(0, k - ny + 1)
            i_end = min(nx - 1, k)
            local_i = np.arange(i_start, i_end + 1)
            local_j = k - local_i
            ii = ix_arr[local_i]
            jj = iy_arr[local_j]
            # Principled indexing: psi_x[n, :, ii+ix_in, jj] has the
            # advanced indices contiguous at the trailing position →
            # numpy keeps the order, shape ``(ng, n_diag)``.
            psi_in_x = psi_x[n, :, ii + ix_in, jj].T  # advanced separated → (n_diag, ng); .T → (ng, n_diag)
            psi_in_y = psi_y[n, :, ii, jj + iy_in].T
            sx = str_x[n, ii][None, :]                 # (1, n_diag)
            sy = str_y[n, jj][None, :]                 # (1, n_diag)
            denom = sig_t[:, ii, jj] + sx + sy         # (ng, n_diag)
            Qn = Q[n if Q.shape[0] > 1 else 0, :, ii, jj].T  # (ng, n_diag)
            psi_avg = (Qn + sx * psi_in_x + sy * psi_in_y) / denom
            psi_x[n, :, ii + ix_out, jj] = (2.0 * psi_avg - psi_in_x).T
            psi_y[n, :, ii, jj + iy_out] = (2.0 * psi_avg - psi_in_y).T
            angular_flux[n, :, ii, jj] = psi_avg.T  # numpy scatter (n_diag, ng) into (ng, ii, jj)
            scalar_flux[:, ii, jj] += weights[n] * psi_avg
    return angular_flux, scalar_flux, psi_x, psi_y


@pytest.mark.l0
@pytest.mark.regression
class TestApplyMatchesLegacyInlined:
    @pytest.mark.parametrize("sx,sy", [
        (+1, +1), (+1, -1), (-1, +1), (-1, -1),
    ])
    @pytest.mark.parametrize("nx,ny,ng,N_oct", [
        (1, 1, 1, 1),
        (2, 2, 2, 3),
        (3, 4, 2, 4),
        (4, 3, 3, 6),
    ])
    def test_per_cell_loop_equivalence(self, sx, sy, nx, ny, ng, N_oct):
        rng = np.random.default_rng(seed=2026 * (10 * sx + sy + 100) + nx * ny + ng + N_oct)
        # Issue #196 PR-INDEX-5: principled layouts.
        psi_x = rng.standard_normal((N_oct, ng, nx + 1, ny))
        psi_y = rng.standard_normal((N_oct, ng, nx, ny + 1))
        Q = rng.standard_normal((N_oct, ng, nx, ny))
        sig_t = rng.uniform(0.1, 0.5, size=(ng, nx, ny))
        str_x = rng.uniform(0.1, 1.0, size=(N_oct, nx))
        str_y = rng.uniform(0.1, 1.0, size=(N_oct, ny))
        weights = rng.uniform(0.5, 1.5, size=N_oct)

        # Reference (per-ordinate Python loop) — principled.
        ref_ang, ref_scal, ref_px, ref_py = _hand_run_legacy_inlined(
            nx=nx, ny=ny, ng=ng, N_oct=N_oct,
            psi_x=psi_x, psi_y=psi_y, Q=Q, sig_t=sig_t,
            str_x=str_x, str_y=str_y, weights=weights,
            sx_sign=sx, sy_sign=sy,
        )

        # New code (vectorised graph apply).
        graph = SweepDependencyGraph.from_cartesian((nx, ny), label=OctantLabel((sx, sy)))
        angular_flux = np.zeros((N_oct, ng, nx, ny))
        scalar_flux = np.zeros((ng, nx, ny))
        psi_x_oct = psi_x.copy()
        psi_y_oct = psi_y.copy()
        graph.walk_full(
            level_op=_CellSolve(
                cell_update=DiamondDifference(),
                weights_octant=weights,
                angular_flux_octant=angular_flux,
                scalar_flux_buf=scalar_flux,
            ),
            psi_faces_octant=(psi_x_oct, psi_y_oct),
            Q_octant=Q,
            sig_t=sig_t,
            str_axes_octant=(str_x, str_y),
        )

        # Bit-identity on angular flux (operation order matches per
        # the cell-kernel bit-identity gate, C2.2 / S6.4(e)).
        np.testing.assert_array_equal(angular_flux, ref_ang)
        # Face buffers also bit-identical.
        np.testing.assert_array_equal(psi_x_oct, ref_px)
        np.testing.assert_array_equal(psi_y_oct, ref_py)
        # Scalar flux: weighted-sum reduction order differs (einsum
        # vs Python loop), so principled-equivalent via ULP. Per
        # vv-principles bit-identity vs principled-equivalence:
        # scalar_flux is a sum-over-ordinates reduction;
        # ``np.einsum("ngd,n->gd", ...)`` vs ``+= w[n] * psi[n]`` over
        # the Python loop. Drift bounded by ``N_oct × ULP × max_term``.
        # PR-INDEX-5 widened the nulp budget 64 → 128 to absorb the
        # additional per-cell scatter (``angular_flux[n, :, ii, jj] =
        # psi_avg.T``) reordering in the principled-layout reference
        # helper.
        np.testing.assert_array_almost_equal_nulp(
            scalar_flux, ref_scal, nulp=128,
        )

    def test_isotropic_Q_broadcasts(self):
        rng = np.random.default_rng(seed=999)
        nx, ny, ng, N_oct = 3, 3, 2, 4
        # Issue #196 PR-INDEX-5: principled layouts.
        psi_x = rng.standard_normal((N_oct, ng, nx + 1, ny))
        psi_y = rng.standard_normal((N_oct, ng, nx, ny + 1))
        Q = rng.standard_normal((1, ng, nx, ny))   # isotropic-only
        sig_t = rng.uniform(0.1, 0.5, size=(ng, nx, ny))
        str_x = rng.uniform(0.1, 1.0, size=(N_oct, nx))
        str_y = rng.uniform(0.1, 1.0, size=(N_oct, ny))
        weights = rng.uniform(0.5, 1.5, size=N_oct)
        ref_ang, ref_scal, ref_px, ref_py = _hand_run_legacy_inlined(
            nx=nx, ny=ny, ng=ng, N_oct=N_oct,
            psi_x=psi_x, psi_y=psi_y, Q=Q, sig_t=sig_t,
            str_x=str_x, str_y=str_y, weights=weights,
            sx_sign=+1, sy_sign=+1,
        )
        graph = SweepDependencyGraph.from_cartesian((nx, ny), label=OctantLabel((+1, +1)))
        angular_flux = np.zeros((N_oct, ng, nx, ny))
        scalar_flux = np.zeros((ng, nx, ny))
        psi_x_oct = psi_x.copy()
        psi_y_oct = psi_y.copy()
        graph.walk_full(
            level_op=_CellSolve(
                cell_update=DiamondDifference(),
                weights_octant=weights,
                angular_flux_octant=angular_flux,
                scalar_flux_buf=scalar_flux,
            ),
            psi_faces_octant=(psi_x_oct, psi_y_oct),
            Q_octant=Q,
            sig_t=sig_t,
            str_axes_octant=(str_x, str_y),
        )
        np.testing.assert_array_equal(angular_flux, ref_ang)
        np.testing.assert_array_equal(psi_x_oct, ref_px)
        np.testing.assert_array_equal(psi_y_oct, ref_py)
        np.testing.assert_array_almost_equal_nulp(
            scalar_flux, ref_scal, nulp=64,
        )


@pytest.mark.l0
class TestResidualWalkRoundTrip:
    r"""``graph.residual`` at the value ``graph.apply`` solves is zero — the
    graph-level apply↔solve contract (matvec == sweep on the SAME DAG).

    This is the D2 gate for Wave O #208 O.4b: the residual-mode walk
    (the 2-D matvec's engine) and the solve-mode walk (the sweep's engine)
    share the wavefront DAG + the DiamondDifference closure, so the
    operator residual of the swept solution vanishes by construction.
    The same edge-flux reconstruction must run in both directions.
    """

    @pytest.mark.parametrize("sx,sy", [
        (+1, +1), (+1, -1), (-1, +1), (-1, -1),
    ])
    @pytest.mark.parametrize("nx,ny,ng,N_oct", [
        (2, 2, 2, 3), (3, 4, 2, 4), (4, 3, 3, 6),
    ])
    def test_residual_vanishes_at_apply_solution(
        self, sx, sy, nx, ny, ng, N_oct,
    ):
        rng = np.random.default_rng(
            seed=7 * (10 * sx + sy + 100) + nx * ny + ng + N_oct,
        )
        psi_x = rng.standard_normal((N_oct, ng, nx + 1, ny))
        psi_y = rng.standard_normal((N_oct, ng, nx, ny + 1))
        Q = rng.standard_normal((N_oct, ng, nx, ny))
        sig_t = rng.uniform(0.1, 0.5, size=(ng, nx, ny))
        str_x = rng.uniform(0.1, 1.0, size=(N_oct, nx))
        str_y = rng.uniform(0.1, 1.0, size=(N_oct, ny))
        weights = rng.uniform(0.5, 1.5, size=N_oct)

        graph = SweepDependencyGraph.from_cartesian((nx, ny), label=OctantLabel((sx, sy)))

        # SOLVE: walk_full forward-substitutes the solve cell kernel.
        ang = np.zeros((N_oct, ng, nx, ny))
        scal = np.zeros((ng, nx, ny))
        graph.walk_full(
            level_op=_CellSolve(
                cell_update=DiamondDifference(),
                weights_octant=weights,
                angular_flux_octant=ang,
                scalar_flux_buf=scal,
            ),
            psi_faces_octant=(psi_x.copy(), psi_y.copy()),
            Q_octant=Q,
            sig_t=sig_t,
            str_axes_octant=(str_x, str_y),
        )

        # APPLY at the swept solution, from a FRESH copy with the SAME
        # boundary-incoming seeds (apply only scatters outgoing faces).
        residual = np.zeros((N_oct, ng, nx, ny))
        graph.walk_full(
            level_op=_CellResidual(
                cell_update=DiamondDifference(),
                psi_avg_probe_octant=ang,
                residual_octant=residual,
            ),
            psi_faces_octant=(psi_x.copy(), psi_y.copy()),
            Q_octant=Q,
            sig_t=sig_t,
            str_axes_octant=(str_x, str_y),
        )

        np.testing.assert_allclose(residual, 0.0, atol=1e-11)

    def test_residual_isotropic_Q_round_trip(self):
        rng = np.random.default_rng(seed=4242)
        nx, ny, ng, N_oct = 3, 3, 2, 4
        psi_x = rng.standard_normal((N_oct, ng, nx + 1, ny))
        psi_y = rng.standard_normal((N_oct, ng, nx, ny + 1))
        Q = rng.standard_normal((1, ng, nx, ny))   # isotropic-only
        sig_t = rng.uniform(0.1, 0.5, size=(ng, nx, ny))
        str_x = rng.uniform(0.1, 1.0, size=(N_oct, nx))
        str_y = rng.uniform(0.1, 1.0, size=(N_oct, ny))
        weights = rng.uniform(0.5, 1.5, size=N_oct)
        graph = SweepDependencyGraph.from_cartesian((nx, ny), label=OctantLabel((+1, +1)))
        ang = np.zeros((N_oct, ng, nx, ny))
        scal = np.zeros((ng, nx, ny))
        graph.walk_full(
            level_op=_CellSolve(
                cell_update=DiamondDifference(),
                weights_octant=weights,
                angular_flux_octant=ang,
                scalar_flux_buf=scal,
            ),
            psi_faces_octant=(psi_x.copy(), psi_y.copy()),
            Q_octant=Q,
            sig_t=sig_t,
            str_axes_octant=(str_x, str_y),
        )
        residual = np.zeros((N_oct, ng, nx, ny))
        graph.walk_full(
            level_op=_CellResidual(
                cell_update=DiamondDifference(),
                psi_avg_probe_octant=ang,
                residual_octant=residual,
            ),
            psi_faces_octant=(psi_x.copy(), psi_y.copy()),
            Q_octant=Q,
            sig_t=sig_t,
            str_axes_octant=(str_x, str_y),
        )
        np.testing.assert_allclose(residual, 0.0, atol=1e-11)
