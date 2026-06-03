r"""Tests for :meth:`DiamondDifference.update_batch` and the
:class:`SweepCellSlice` packet (Wave 2 / C2.2).

The Wave 2 plan introduces a vectorised per-level update for the 2-D
Cartesian wavefront sweep. ``update_batch`` consumes a
:class:`SweepCellSlice` packet and returns ``psi_avg`` of shape
``(N_oct, n_diag, ng)``, scattering the outgoing face fluxes back into
the persistent BC buffers ``psi_x`` / ``psi_y`` in place.

This file pins three properties:

1. **Bit-identity to the legacy inlined wavefront math** at
   ``orpheus.sn.sweep._sweep_2d_wavefront`` lines 847-871. The
   refactor's load-bearing constraint per the explorer surface map
   is that ``denom = sig_t + sx + sy``,
   ``psi_avg = (Q + sx*psi_in_x + sy*psi_in_y) / denom``, and
   ``psi_out = 2*psi_avg - psi_in`` retain their operation order
   exactly.
2. **Face-flux scatter** — the WDD spatial closure values land at
   ``face_out_x_idx`` / ``face_out_y_idx`` and nowhere else.
3. **Default NotImplementedError** on the base class — strategies
   that don't override ``update_batch`` fail loudly rather than
   silently.

V&V tags: ``L0`` (term-level verification of the closed-form
algebraic update against the WDD balance formula) plus a
``regression``-style bit-identity check against the legacy inlined
math.
"""

from __future__ import annotations

from dataclasses import replace
from typing import ClassVar

import numpy as np
import pytest

from orpheus.sn.spatial.cell_update import (
    CellResult,
    CellUpdateBase,
    CellVisit,
    SweepCellSlice,
    UpstreamState,
)
from orpheus.sn.spatial.diamond import DiamondDifference


# ─────────────────────────────────────────────────────────────────────
# Helper builders
# ─────────────────────────────────────────────────────────────────────


def _build_slice_kwargs(
    *,
    nx: int,
    ny: int,
    N_oct: int,
    ng: int,
    diag_cells: list[tuple[int, int]],
    sx_sign: int,
    sy_sign: int,
    seed: int = 0,
    Q_shape_leading: int | None = None,  # None -> N_oct; 1 -> isotropic
) -> tuple[SweepCellSlice, np.ndarray, np.ndarray]:
    """Build a SweepCellSlice + return references to psi_x and psi_y.

    Issue #196 PR-INDEX-5: principled layout — ``psi_x: (N_oct, ng,
    nx+1, ny)``, ``psi_y: (N_oct, ng, nx, ny+1)``,
    ``Q: (N_oct or 1, ng, nx, ny)``, ``sig_t: (ng, nx, ny)``.
    """
    rng = np.random.default_rng(seed)
    psi_x = rng.standard_normal((N_oct, ng, nx + 1, ny))
    psi_y = rng.standard_normal((N_oct, ng, nx, ny + 1))
    Q_lead = N_oct if Q_shape_leading is None else Q_shape_leading
    Q = rng.standard_normal((Q_lead, ng, nx, ny))
    sig_t = rng.uniform(0.1, 0.5, size=(ng, nx, ny))
    str_x = rng.uniform(0.1, 1.0, size=(N_oct, nx))
    str_y = rng.uniform(0.1, 1.0, size=(N_oct, ny))

    ix_in = 0 if sx_sign >= 0 else 1
    ix_out = 1 if sx_sign >= 0 else 0
    iy_in = 0 if sy_sign >= 0 else 1
    iy_out = 1 if sy_sign >= 0 else 0

    ii = np.array([c[0] for c in diag_cells], dtype=int)
    jj = np.array([c[1] for c in diag_cells], dtype=int)

    return (
        SweepCellSlice(
            ii=ii, jj=jj,
            face_in_x_idx=ii + ix_in,
            face_out_x_idx=ii + ix_out,
            face_in_y_idx=jj + iy_in,
            face_out_y_idx=jj + iy_out,
            psi_x=psi_x, psi_y=psi_y,
            Q=Q, sig_t=sig_t, str_x=str_x, str_y=str_y,
        ),
        psi_x, psi_y,
    )


def _legacy_inlined_psi_avg(
    slice_args: SweepCellSlice, psi_x_pre: np.ndarray, psi_y_pre: np.ndarray,
) -> np.ndarray:
    """Reference implementation: per-ordinate Python loop in principled
    ``(N_oct, ng, n_diag)`` layout.

    Issue #196 PR-INDEX-5: ``psi_x`` is ``(N_oct, ng, nx+1, ny)`` etc.;
    output principled ``(N_oct, ng, n_diag)``.
    """
    s = slice_args
    N_oct = s.psi_x.shape[0]
    ng = s.psi_x.shape[1]
    n_diag = len(s.ii)
    psi_avg_ref = np.empty((N_oct, ng, n_diag))
    for n in range(N_oct):
        psi_in_x = psi_x_pre[n, :, s.face_in_x_idx, s.jj]   # (n_diag, ng) — advanced at end
        psi_in_y = psi_y_pre[n, :, s.ii, s.face_in_y_idx]   # (n_diag, ng)
        # Note: numpy gives (n_diag, ng) when advanced indices trail two
        # basic slices; transpose to (ng, n_diag) for principled output.
        psi_in_x = psi_in_x.T                                # (ng, n_diag)
        psi_in_y = psi_in_y.T                                # (ng, n_diag)
        sx_ii = s.str_x[n, s.ii][None, :]                   # (1, n_diag)
        sy_jj = s.str_y[n, s.jj][None, :]                   # (1, n_diag)
        denom = s.sig_t[:, s.ii, s.jj] + sx_ii + sy_jj      # (ng, n_diag)
        Q_n = s.Q[n if s.Q.shape[0] > 1 else 0, :, s.ii, s.jj].T  # (ng, n_diag)
        psi_avg_ref[n] = (
            Q_n + sx_ii * psi_in_x + sy_jj * psi_in_y
        ) / denom
    return psi_avg_ref


# ─────────────────────────────────────────────────────────────────────
# L0: Closed-form check on a 1-cell, 1-ordinate, 1-group batch
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestSingleCellClosedForm:
    """Smallest possible batch — analytical hand calculation."""

    def test_single_cell_psi_avg_matches_balance_formula(self):
        """psi_avg = (Q + sx*psi_in_x + sy*psi_in_y) / (Σ_t + sx + sy).

        Issue #196 PR-INDEX-5: psi_x/psi_y/Q/sig_t principled layout.
        """
        # 1 ordinate, 1 cell at (i=0, j=0), 1 group, sx > 0, sy > 0.
        # psi_x principled: (N_oct=1, ng=1, nx+1=2, ny=1).
        psi_x = np.zeros((1, 1, 2, 1))
        psi_y = np.zeros((1, 1, 1, 2))
        psi_x[0, 0, 0, 0] = 4.0   # face_in_x at i=0 (sx>=0 -> ix_in=0)
        psi_y[0, 0, 0, 0] = 8.0   # face_in_y at j=0
        # Q principled: (N_oct or 1, ng=1, nx=1, ny=1).
        Q = np.array([[[[16.0]]]])              # (1, 1, 1, 1)
        sig_t = np.array([[[2.0]]])             # (ng=1, nx=1, ny=1)
        str_x = np.array([[3.0]])               # (N_oct, nx)
        str_y = np.array([[5.0]])               # (N_oct, ny)
        slice_args = SweepCellSlice(
            ii=np.array([0]), jj=np.array([0]),
            face_in_x_idx=np.array([0]),  face_out_x_idx=np.array([1]),
            face_in_y_idx=np.array([0]),  face_out_y_idx=np.array([1]),
            psi_x=psi_x, psi_y=psi_y,
            Q=Q, sig_t=sig_t, str_x=str_x, str_y=str_y,
        )
        psi_avg = DiamondDifference().update_batch(slice_args)
        # (16 + 3*4 + 5*8) / (2 + 3 + 5) = (16 + 12 + 40) / 10 = 6.8
        np.testing.assert_array_equal(psi_avg, 6.8)
        # Outgoing face-flux: 2*6.8 - 4 = 9.6 (x), 2*6.8 - 8 = 5.6 (y)
        np.testing.assert_array_equal(psi_x[0, 0, 1, 0], 9.6)  # face_out_x at ix=1
        np.testing.assert_array_equal(psi_y[0, 0, 0, 1], 5.6)  # face_out_y at iy=1
        # Incoming face values must be untouched.
        np.testing.assert_array_equal(psi_x[0, 0, 0, 0], 4.0)
        np.testing.assert_array_equal(psi_y[0, 0, 0, 0], 8.0)

    def test_single_cell_negative_sign_octant(self):
        """sx < 0, sy < 0: ix_in=1, ix_out=0; iy_in=1, iy_out=0."""
        psi_x = np.zeros((1, 1, 2, 1))
        psi_y = np.zeros((1, 1, 1, 2))
        psi_x[0, 0, 1, 0] = 4.0   # face_in_x at i=1 (sx<0 -> ix_in=1)
        psi_y[0, 0, 0, 1] = 8.0   # face_in_y at j=1
        Q = np.array([[[[16.0]]]])
        sig_t = np.array([[[2.0]]])
        str_x = np.array([[3.0]])
        str_y = np.array([[5.0]])
        slice_args = SweepCellSlice(
            ii=np.array([0]), jj=np.array([0]),
            face_in_x_idx=np.array([1]),  face_out_x_idx=np.array([0]),
            face_in_y_idx=np.array([1]),  face_out_y_idx=np.array([0]),
            psi_x=psi_x, psi_y=psi_y,
            Q=Q, sig_t=sig_t, str_x=str_x, str_y=str_y,
        )
        psi_avg = DiamondDifference().update_batch(slice_args)
        np.testing.assert_array_equal(psi_avg, 6.8)
        # face_out_x at i=0; face_out_y at j=0.
        np.testing.assert_array_equal(psi_x[0, 0, 0, 0], 9.6)
        np.testing.assert_array_equal(psi_y[0, 0, 0, 0], 5.6)


# ─────────────────────────────────────────────────────────────────────
# L0 / regression: Bit-identity to the legacy inlined sweep math
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
@pytest.mark.regression
class TestBitIdenticalToLegacyInlinedMath:
    r"""Per-element bit-equality against ``_sweep_2d_wavefront`` math.

    The test runs the legacy per-ordinate inlined math (a Python
    loop over ``n``) and compares to ``update_batch``'s output. They
    must agree at IEEE-754 ULP — the operation order is identical
    by construction.
    """

    @pytest.mark.parametrize("sx_sign,sy_sign", [
        (+1, +1), (+1, -1), (-1, +1), (-1, -1),
    ])
    def test_bit_identical_3x3_4ord_2g(self, sx_sign, sy_sign):
        """3×3 grid, 4 ordinates per octant, 2 groups, 3-cell anti-diag."""
        # Anti-diagonal of length 3 on a 3×3 grid (e.g. cells [(0,2),(1,1),(2,0)]
        # for an octant traversal — the actual indices depend on direction
        # but for this unit test we just need a valid (ii, jj) layout.)
        diag_cells = [(0, 0), (1, 1), (2, 2)]   # diagonal-of-3
        slice_args, psi_x, psi_y = _build_slice_kwargs(
            nx=3, ny=3, N_oct=4, ng=2,
            diag_cells=diag_cells,
            sx_sign=sx_sign, sy_sign=sy_sign,
            seed=12,
        )
        psi_x_pre = psi_x.copy()
        psi_y_pre = psi_y.copy()
        legacy_psi_avg = _legacy_inlined_psi_avg(
            slice_args, psi_x_pre, psi_y_pre,
        )
        # Now reset psi_x, psi_y and call update_batch.
        psi_x[...] = psi_x_pre
        psi_y[...] = psi_y_pre
        psi_avg = DiamondDifference().update_batch(slice_args)
        # Bit-equality: must match per-ordinate Python-loop reference.
        np.testing.assert_array_equal(psi_avg, legacy_psi_avg)

    def test_isotropic_Q_broadcasts_correctly(self):
        """Q with leading dim 1 (isotropic source) must broadcast cleanly."""
        slice_args, psi_x, psi_y = _build_slice_kwargs(
            nx=3, ny=3, N_oct=4, ng=2,
            diag_cells=[(0, 0), (1, 1), (2, 2)],
            sx_sign=+1, sy_sign=+1, seed=21,
            Q_shape_leading=1,        # isotropic-shaped Q
        )
        psi_x_pre = psi_x.copy()
        psi_y_pre = psi_y.copy()
        legacy_psi_avg = _legacy_inlined_psi_avg(
            slice_args, psi_x_pre, psi_y_pre,
        )
        psi_x[...] = psi_x_pre
        psi_y[...] = psi_y_pre
        psi_avg = DiamondDifference().update_batch(slice_args)
        np.testing.assert_array_equal(psi_avg, legacy_psi_avg)


# ─────────────────────────────────────────────────────────────────────
# L0: Face-flux scatter writes the right indices and only those
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestFaceFluxScatter:
    """Outgoing face fluxes land at face_out_x_idx / face_out_y_idx only."""

    def test_only_face_out_indices_change(self):
        """psi_x[face_out_x_idx, jj] changes; other entries unchanged.

        Issue #196 PR-INDEX-5: psi_x principled (N_oct, ng, nx+1, ny).
        """
        slice_args, psi_x, psi_y = _build_slice_kwargs(
            nx=3, ny=3, N_oct=4, ng=2,
            diag_cells=[(0, 0), (1, 1), (2, 2)],
            sx_sign=+1, sy_sign=+1, seed=33,
        )
        psi_x_pre = psi_x.copy()
        psi_y_pre = psi_y.copy()
        DiamondDifference().update_batch(slice_args)

        # Face-out indices: ii + 1, jj + 1 (sx,sy > 0).
        face_out_x_idx = slice_args.ii + 1   # [1, 2, 3]
        face_out_y_idx = slice_args.jj + 1   # [1, 2, 3]

        # Build a mask of expected-changed locations.
        # psi_x: (N_oct, ng, nx+1, ny); change at (:, :, face_out_x[k], jj[k]).
        x_changed_mask = np.zeros_like(psi_x, dtype=bool)
        for k, (i, j) in enumerate(zip(slice_args.ii, slice_args.jj)):
            x_changed_mask[:, :, face_out_x_idx[k], j] = True
        y_changed_mask = np.zeros_like(psi_y, dtype=bool)
        for k, (i, j) in enumerate(zip(slice_args.ii, slice_args.jj)):
            y_changed_mask[:, :, i, face_out_y_idx[k]] = True

        # Outside the masked region nothing changed.
        np.testing.assert_array_equal(
            psi_x[~x_changed_mask], psi_x_pre[~x_changed_mask],
        )
        np.testing.assert_array_equal(
            psi_y[~y_changed_mask], psi_y_pre[~y_changed_mask],
        )

    def test_face_out_value_matches_wdd_closure(self):
        """psi_x[face_out_x_idx] == 2*psi_avg - psi_in_x.

        Issue #196 PR-INDEX-5: principled indexing.
        """
        slice_args, psi_x, psi_y = _build_slice_kwargs(
            nx=3, ny=3, N_oct=4, ng=2,
            diag_cells=[(0, 0), (1, 1), (2, 2)],
            sx_sign=+1, sy_sign=+1, seed=44,
        )
        psi_x_pre = psi_x.copy()
        psi_y_pre = psi_y.copy()
        psi_avg = DiamondDifference().update_batch(slice_args)
        # Reconstruct expected outgoing face values from psi_avg + pre-buffer.
        # psi_x principled: index as [:, :, face_idx, jj].
        psi_in_x = psi_x_pre[:, :, slice_args.face_in_x_idx, slice_args.jj]
        psi_in_y = psi_y_pre[:, :, slice_args.ii, slice_args.face_in_y_idx]
        expected_psi_out_x = 2.0 * psi_avg - psi_in_x
        expected_psi_out_y = 2.0 * psi_avg - psi_in_y
        actual_psi_out_x = psi_x[
            :, :, slice_args.face_out_x_idx, slice_args.jj,
        ]
        actual_psi_out_y = psi_y[
            :, :, slice_args.ii, slice_args.face_out_y_idx,
        ]
        np.testing.assert_array_equal(actual_psi_out_x, expected_psi_out_x)
        np.testing.assert_array_equal(actual_psi_out_y, expected_psi_out_y)


# ─────────────────────────────────────────────────────────────────────
# L0: NotImplementedError for strategies that don't override
# ─────────────────────────────────────────────────────────────────────


class _NoBatchStrategy(CellUpdateBase, key="_no_batch_strategy_test"):
    """Test stub: only overrides ``update`` and ``residual``, not ``update_batch``."""

    is_linear: ClassVar[bool] = True
    is_positivity_preserving: ClassVar[bool] = False

    def update(
        self,
        visit: CellVisit,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> CellResult:
        return CellResult(cell_average_flux=source / total_xs)

    def residual(
        self,
        cell_avg: np.ndarray,
        visit: CellVisit,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> np.ndarray:
        del visit, upstream_state
        return total_xs * cell_avg - source


@pytest.mark.l0
class TestDefaultRaisesNotImplemented:
    """Strategies that don't override update_batch fail loudly."""

    def test_default_raises(self):
        slice_args, _, _ = _build_slice_kwargs(
            nx=2, ny=2, N_oct=1, ng=1,
            diag_cells=[(0, 0)],
            sx_sign=+1, sy_sign=+1, seed=55,
        )
        strat = _NoBatchStrategy()
        with pytest.raises(NotImplementedError, match="update_batch"):
            strat.update_batch(slice_args)

    def test_residual_batch_default_raises(self):
        """Strategies that don't override residual_batch fail loudly too."""
        slice_args, _, _ = _build_slice_kwargs(
            nx=2, ny=2, N_oct=1, ng=1,
            diag_cells=[(0, 0)],
            sx_sign=+1, sy_sign=+1, seed=55,
        )
        strat = _NoBatchStrategy()
        with pytest.raises(NotImplementedError, match="residual_batch"):
            strat.residual_batch(slice_args)


# ─────────────────────────────────────────────────────────────────────
# Wave O #208 O.4b — residual_batch (the batched apply direction)
# ─────────────────────────────────────────────────────────────────────


def _probe_field_from_level(
    slice_args: SweepCellSlice, psi_avg_level: np.ndarray,
) -> np.ndarray:
    """Scatter a per-level ``(N_oct, ng, n_diag)`` value into a full
    ``(N_oct, ng, nx, ny)`` probe field (zeros elsewhere)."""
    N_oct, ng = slice_args.psi_x.shape[0], slice_args.psi_x.shape[1]
    nx, ny = slice_args.sig_t.shape[1], slice_args.sig_t.shape[2]
    probe = np.zeros((N_oct, ng, nx, ny))
    probe[:, :, slice_args.ii, slice_args.jj] = psi_avg_level
    return probe


@pytest.mark.l0
class TestResidualBatchClosedForm:
    """Single-cell closed-form: residual = denom·ψ̄ − (Q + sx·ψ_in_x + sy·ψ_in_y)."""

    def test_residual_at_solution_is_zero(self):
        """At the swept ψ̄ (= 6.8 from the update_batch closed-form test),
        the residual vanishes — the per-cell round-trip contract."""
        psi_x = np.zeros((1, 1, 2, 1))
        psi_y = np.zeros((1, 1, 1, 2))
        psi_x[0, 0, 0, 0] = 4.0
        psi_y[0, 0, 0, 0] = 8.0
        Q = np.array([[[[16.0]]]])
        sig_t = np.array([[[2.0]]])
        str_x = np.array([[3.0]])
        str_y = np.array([[5.0]])
        probe = np.array([[[[6.8]]]])           # the solved ψ̄ (denom=10)
        slice_args = SweepCellSlice(
            ii=np.array([0]), jj=np.array([0]),
            face_in_x_idx=np.array([0]), face_out_x_idx=np.array([1]),
            face_in_y_idx=np.array([0]), face_out_y_idx=np.array([1]),
            psi_x=psi_x, psi_y=psi_y,
            Q=Q, sig_t=sig_t, str_x=str_x, str_y=str_y,
            psi_avg_probe=probe,
        )
        residual = DiamondDifference().residual_batch(slice_args)
        # 10*6.8 - (16 + 3*4 + 5*8) = 68 - 68 = 0.
        np.testing.assert_allclose(residual, 0.0, atol=1e-13)
        # Diamond closure still scatters the outgoing faces with the probe.
        np.testing.assert_array_equal(psi_x[0, 0, 1, 0], 2 * 6.8 - 4.0)  # 9.6
        np.testing.assert_array_equal(psi_y[0, 0, 0, 1], 2 * 6.8 - 8.0)  # 5.6

    def test_residual_off_solution_is_affine(self):
        """A probe shifted by δ from the solution shifts the residual by
        denom·δ — the residual is linear in ψ̄."""
        psi_x = np.zeros((1, 1, 2, 1)); psi_x[0, 0, 0, 0] = 4.0
        psi_y = np.zeros((1, 1, 1, 2)); psi_y[0, 0, 0, 0] = 8.0
        Q = np.array([[[[16.0]]]]); sig_t = np.array([[[2.0]]])
        str_x = np.array([[3.0]]); str_y = np.array([[5.0]])
        probe = np.array([[[[7.0]]]])           # 0.2 above the solution 6.8
        slice_args = SweepCellSlice(
            ii=np.array([0]), jj=np.array([0]),
            face_in_x_idx=np.array([0]), face_out_x_idx=np.array([1]),
            face_in_y_idx=np.array([0]), face_out_y_idx=np.array([1]),
            psi_x=psi_x, psi_y=psi_y,
            Q=Q, sig_t=sig_t, str_x=str_x, str_y=str_y,
            psi_avg_probe=probe,
        )
        residual = DiamondDifference().residual_batch(slice_args)
        # 10*7.0 - 68 = 2.0  (== denom · δ = 10 · 0.2).
        np.testing.assert_allclose(residual, 2.0, atol=1e-13)


@pytest.mark.l0
class TestResidualBatchRoundTrip:
    r"""The batched apply↔solve contract: residual_batch at the value
    update_batch returns is zero (the analogue of the per-cell
    DiamondDifference.residual ↔ .update round-trip)."""

    @pytest.mark.parametrize("sx_sign,sy_sign", [
        (+1, +1), (+1, -1), (-1, +1), (-1, -1),
    ])
    def test_residual_vanishes_at_update_batch_solution(self, sx_sign, sy_sign):
        slice_args, psi_x, psi_y = _build_slice_kwargs(
            nx=3, ny=3, N_oct=4, ng=2,
            diag_cells=[(0, 0), (1, 1), (2, 2)],
            sx_sign=sx_sign, sy_sign=sy_sign, seed=77,
        )
        psi_x_pre = psi_x.copy()
        psi_y_pre = psi_y.copy()
        # SOLVE: update_batch returns ψ̄ and scatters outgoing faces.
        psi_avg = DiamondDifference().update_batch(slice_args)
        probe = _probe_field_from_level(slice_args, psi_avg)
        # Reset the buffers so residual_batch sees the SAME incoming faces.
        psi_x[...] = psi_x_pre
        psi_y[...] = psi_y_pre
        slice_apply = replace(slice_args, psi_avg_probe=probe)
        # APPLY at the swept ψ̄ — residual must vanish.
        residual = DiamondDifference().residual_batch(slice_apply)
        np.testing.assert_allclose(residual, 0.0, atol=1e-13)

    def test_isotropic_Q_round_trip(self):
        """Round-trip holds with an isotropic-shaped (leading-1) Q too."""
        slice_args, psi_x, psi_y = _build_slice_kwargs(
            nx=3, ny=3, N_oct=4, ng=2,
            diag_cells=[(0, 0), (1, 1), (2, 2)],
            sx_sign=+1, sy_sign=+1, seed=88,
            Q_shape_leading=1,
        )
        psi_x_pre = psi_x.copy()
        psi_y_pre = psi_y.copy()
        psi_avg = DiamondDifference().update_batch(slice_args)
        probe = _probe_field_from_level(slice_args, psi_avg)
        psi_x[...] = psi_x_pre
        psi_y[...] = psi_y_pre
        slice_apply = replace(slice_args, psi_avg_probe=probe)
        residual = DiamondDifference().residual_batch(slice_apply)
        np.testing.assert_allclose(residual, 0.0, atol=1e-13)
