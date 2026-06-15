r"""GATE 1 (#158) — LinearDiscontinuous per-cell occupant, isolated.

Verifies the slab/Cartesian LD cell update in ISOLATION (no sweep wiring), per
the staged plan (`.claude/plans/issue_158_spatial_cellupdate_carve.md`) and the
test-architect spec (`agent-memory/test-architect/issue_158_ld_spatial_verification.md`).

The centrepiece is the **linear-exactness oracle**: LD represents the in-cell
flux linearly, so it must reproduce a linear-in-x manufactured flux
(:math:`\psi=a+bx`) to machine precision — recovering the cell-average AND the
outflow face exactly.  This is a *structurally-independent* correctness gate
(a physical property of LD, not a re-statement of the 2×2 algebra), so it
catches the slope-row sign trap that the round-trip alone cannot (the trap is
invisible to ``residual(update(...))=0`` if the sign error is consistent
between the solve and apply directions — LM-1989 memo §1.4/§6).  The 2×2 and
its linear-exactness were first verified symbolically with SymPy
(`.claude/jobs/.../tmp/ld_slab_derivation.py`); this gate verifies the *code*
reproduces it.

Companion gates: the per-cell round-trip (lockstep), residual linearity, the
class traits, and the slab-only / source-shape guards.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, slab_streaming
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.spatial import LinearDiscontinuous, UpstreamState
from orpheus.sn.spatial.affine_closure import cell_average, source_emission
from orpheus.sn.spatial.cell_update import CellResult, CellUpdateBase, CellVisit


def _slab_mesh(nx: int = 5, length: float = 1.0) -> Mesh1D:
    return Mesh1D(
        edges=np.linspace(0.0, length, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )


def _slab_visit(cell_idx: int = 2, n_ord: int = 4, *, positive: bool = True):
    """Build a slab CellVisit + its |μ|, h, and cell edges for cell_idx.

    Returns ``(visit, st, x0, h)`` — the most-positive ordinate (forward
    sweep, inflow = left face) when ``positive``, else the most-negative.
    """
    mesh = _slab_mesh(nx=5, length=1.0)
    quad = Quadrature.gauss_legendre(n_ord)
    op = slab_streaming(mesh, quad)
    direction_idx = quad.N - 1 if positive else 0
    st = op.streaming_terms(cell_idx, direction_idx)
    x0 = float(mesh.edges[cell_idx])
    h = float(mesh.edges[cell_idx + 1] - mesh.edges[cell_idx])
    visit = CellVisit(
        cell_idx=cell_idx, streaming_terms=st, face_area_downstream=1.0,
    )
    return visit, st, x0, h


def _moment_source(s_bar: np.ndarray, s_hat: np.ndarray) -> np.ndarray:
    """Stack the (average, slope) moment source into the LD ``(2, ng)`` shape."""
    return np.stack([s_bar, s_hat], axis=0)


# ═══════════════════════════════════════════════════════════════════════
# THE linear-exactness oracle (structurally-independent correctness gate)
# ═══════════════════════════════════════════════════════════════════════

class TestLDLinearExactness:
    r"""LD reproduces a linear-in-x flux exactly (cell-average AND outflow).

    Manufacture :math:`\psi_g(x) = a_g + b_g x`.  The slab SN equation for
    ordinate :math:`\mu` gives the per-ordinate source
    :math:`Q_g(x) = \mu b_g + \Sigma_{t,g}(a_g + b_g x)` (linear), with
    cell-integrated moments

    .. math::

        \bar Q_g = \mu b_g + \Sigma_{t,g}(a_g + b_g x_m), \qquad
        \hat Q_g = \Sigma_{t,g} b_g h/2 ,

    and prepared source moments ``source[0] = Q̄·h``, ``source[1] = Q̂·h``.
    Expected: :math:`\bar\psi_g = a_g + b_g x_m`,
    :math:`\psi_{\rm out,g} = a_g + b_g(x_0+h)`.
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize("n_groups", [1, 2, 4])
    def test_linear_flux_recovered_exactly(self, n_groups: int) -> None:
        visit, st, x0, h = _slab_visit(cell_idx=2, n_ord=4, positive=True)
        assert st.abs_mu is not None          # slab populates it; narrow type
        mu = float(st.abs_mu)
        xm = x0 + h / 2.0

        rng = np.random.default_rng(158)
        a = rng.uniform(0.5, 2.0, n_groups)            # intercepts (a_g > 0)
        b = rng.uniform(-1.0, 1.0, n_groups)           # slopes (mixed signs)
        sig = rng.uniform(0.6, 1.8, n_groups)          # heterogeneous Σ_t

        psi_in = a + b * x0                            # ψ_exact(left face)
        Q_bar = mu * b + sig * (a + b * xm)            # cell-average source
        Q_hat = sig * b * h / 2.0                      # slope-moment source
        source = _moment_source(Q_bar * h, Q_hat * h)
        upstream = UpstreamState(spatial_upstream=psi_in, angular_upstream=None)

        result = LinearDiscontinuous().update(visit, sig, source, upstream)

        psi_bar_exact = a + b * xm
        psi_out_exact = a + b * (x0 + h)
        np.testing.assert_allclose(
            result.cell_average_flux, psi_bar_exact, rtol=1e-12, atol=1e-13,
        )
        np.testing.assert_allclose(
            result.outgoing_spatial_flux, psi_out_exact, rtol=1e-12, atol=1e-13,
        )
        # The slope is recoverable from (average, outflow): ψ̂ = ψ_out − ψ̄.
        psi_hat = result.outgoing_spatial_flux - result.cell_average_flux
        np.testing.assert_allclose(psi_hat, b * h / 2.0, rtol=1e-12, atol=1e-13)
        assert result.outgoing_angular_state is None   # slab: no redistribution


# ═══════════════════════════════════════════════════════════════════════
# Round-trip (lockstep) — update / residual describe ONE system
# ═══════════════════════════════════════════════════════════════════════

class TestLDRoundTrip:
    """``residual(update(q).cell_average_flux, q) == 0`` to FP noise.

    Term activation: streaming + collision (slab has NO angular redistribution
    — the angular column is nulled BY GEOMETRY, not by ansatz).  ``n_groups``
    in {1,2,4} with heterogeneous Σ_t: 1G is the degenerate control, 2G/4G the
    real gate (L2).
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize("n_groups", [1, 2, 4])
    def test_residual_zero_at_solved_cell_avg(self, n_groups: int) -> None:
        visit, _, _, _ = _slab_visit(cell_idx=2, n_ord=4)
        rng = np.random.default_rng(7 + n_groups)
        sig = rng.uniform(0.6, 1.5, n_groups)
        source = _moment_source(
            rng.uniform(0.2, 2.4, n_groups), rng.uniform(-0.6, 0.6, n_groups),
        )
        psi_in = rng.uniform(0.05, 0.4, n_groups)
        upstream = UpstreamState(spatial_upstream=psi_in, angular_upstream=None)

        strat = LinearDiscontinuous()
        result = strat.update(visit, sig, source, upstream)
        residual = strat.residual(
            result.cell_average_flux, visit, sig, source, upstream,
        )
        np.testing.assert_allclose(residual, 0.0, atol=1e-12)

    @pytest.mark.foundation
    def test_residual_linear_in_cell_avg(self) -> None:
        """LD closure is linear → residual is linear in the probe (mode #2)."""
        visit, _, _, _ = _slab_visit(cell_idx=2, n_ord=4)
        sig = np.array([1.2, 0.7])
        source = _moment_source(np.array([1.5, 0.4]), np.array([0.3, -0.2]))
        upstream = UpstreamState(
            spatial_upstream=np.array([0.2, 0.1]), angular_upstream=None,
        )
        strat = LinearDiscontinuous()
        rng = np.random.default_rng(42)
        pa, pb, lam = rng.normal(1.0, 0.5, 2), rng.normal(2.0, 0.5, 2), 0.37
        r_a = strat.residual(pa, visit, sig, source, upstream)
        r_b = strat.residual(pb, visit, sig, source, upstream)
        r_mix = strat.residual(
            lam * pa + (1 - lam) * pb, visit, sig, source, upstream,
        )
        np.testing.assert_allclose(
            r_mix, lam * r_a + (1 - lam) * r_b, rtol=1e-12, atol=1e-13,
        )

    @pytest.mark.foundation
    def test_residual_affine_in_average_source(self) -> None:
        """∂r/∂(average source) = -1 (mode #3 missing-factor / convention)."""
        visit, _, _, _ = _slab_visit(cell_idx=2, n_ord=4)
        sig = np.array([1.2, 0.7])
        source = _moment_source(np.array([1.5, 0.4]), np.array([0.3, -0.2]))
        upstream = UpstreamState(
            spatial_upstream=np.array([0.2, 0.1]), angular_upstream=None,
        )
        strat = LinearDiscontinuous()
        rng = np.random.default_rng(7)
        probe = rng.normal(1.0, 0.5, 2)
        d_avg = rng.normal(0.0, 0.3, 2)
        shifted = _moment_source(source[0] + d_avg, source[1])
        r0 = strat.residual(probe, visit, sig, source, upstream)
        r1 = strat.residual(probe, visit, sig, shifted, upstream)
        np.testing.assert_allclose(r1, r0 - d_avg, rtol=1e-12, atol=1e-13)


# ═══════════════════════════════════════════════════════════════════════
# Traits + guards
# ═══════════════════════════════════════════════════════════════════════

class TestLDTraits:
    @pytest.mark.foundation
    def test_is_affine_scannable_true(self) -> None:
        """LD's slope is eliminated by the Schur complement → the single-upstream
        affine recurrence holds → LD IS affine-scannable (#158 Increment B) and
        rides CumprodScan/ScanMarch via the coefficient model."""
        assert LinearDiscontinuous.is_affine_scannable is True

    @pytest.mark.foundation
    def test_is_linear_true_not_positivity_preserving(self) -> None:
        assert LinearDiscontinuous.is_linear is True
        assert LinearDiscontinuous.is_positivity_preserving is False

    @pytest.mark.foundation
    def test_registered_under_key(self) -> None:
        assert CellUpdateBase.registry["linear_discontinuous"] is LinearDiscontinuous

    @pytest.mark.foundation
    def test_create_from_registry(self) -> None:
        assert isinstance(
            CellUpdateBase.create("linear_discontinuous"), LinearDiscontinuous,
        )


class TestLDGuards:
    @pytest.mark.foundation
    def test_curvilinear_visit_raises(self) -> None:
        """Curvilinear LD is unpublished/deferred — a curvilinear visit
        (angular_upstream not None) must raise, not silently mis-solve."""
        visit, _, _, _ = _slab_visit(cell_idx=2, n_ord=4)
        sig = np.array([1.0, 0.8])
        source = _moment_source(np.array([1.0, 0.5]), np.array([0.1, 0.0]))
        curvilinear = UpstreamState(
            spatial_upstream=np.array([0.2, 0.1]),
            angular_upstream=np.array([0.15, 0.05]),   # present → curvilinear
        )
        with pytest.raises(NotImplementedError, match="curvilinear"):
            LinearDiscontinuous().update(visit, sig, source, curvilinear)

    @pytest.mark.foundation
    def test_diamond_shaped_source_raises(self) -> None:
        """A DD-shaped (ng,) source drops the LD slope → must fail loudly."""
        visit, _, _, _ = _slab_visit(cell_idx=2, n_ord=4)
        sig = np.array([1.0, 0.8])
        dd_source = np.array([1.0, 0.5])               # (ng,) — wrong for LD
        upstream = UpstreamState(
            spatial_upstream=np.array([0.2, 0.1]), angular_upstream=None,
        )
        with pytest.raises(ValueError, match="two-moment source"):
            LinearDiscontinuous().update(visit, sig, dd_source, upstream)


class TestLDKernel:
    r"""Group-2 batched kernel (the polymorphic DAG-wavefront path).

    LD's ``cell_kernel_batch`` / ``residual_kernel_batch`` are the contract
    :class:`FullFieldWavefront` (the any-d DAG oracle) calls — exactly as it
    calls DD.  These pin the batched round-trip and that the batched kernel
    (÷V) agrees with the per-cell ``update`` (×V) — the SAME LD, two forms
    (the project precedent: DD keeps its per-cell and batched residuals
    separate, pinned consistent here).
    """

    @pytest.mark.foundation
    def test_batched_round_trip(self) -> None:
        """residual_kernel_batch at the solved psi_avg vanishes (FP noise)."""
        _, st, _, _ = _slab_visit(cell_idx=2, n_ord=4)
        assert st.abs_mu is not None and st.volume is not None
        mu, h = float(st.abs_mu), float(st.volume)
        sig = np.array([1.2, 0.7])
        s_axes = (np.array([[[2.0 * mu / h]]]),)        # (1,1,1)
        sigt = sig[:, None]                              # (2,1)
        Q = np.array([2.0, 0.5])[None, :, None]          # (1,2,1)
        pin = (np.array([0.3, 0.1])[None, :, None],)     # (1,2,1)
        strat = LinearDiscontinuous()
        psi_avg, (psi_out,) = strat.cell_kernel_batch(
            psi_in=pin, s_axes=s_axes, sigt_cells=sigt, Q_cells=Q,
        )
        resid, (psi_out2,) = strat.residual_kernel_batch(
            psi_bar=psi_avg, psi_in=pin, s_axes=s_axes, sigt_cells=sigt, Q_cells=Q,
        )
        np.testing.assert_allclose(resid, 0.0, atol=1e-12)
        np.testing.assert_allclose(psi_out, psi_out2)

    @pytest.mark.foundation
    def test_group1_equals_group2_flat(self) -> None:
        """Per-cell update (×V) ≡ batched kernel (÷V), flat source (Q̂=0).

        The two capability groups implement the SAME LD; this pins them
        consistent (the per-cell ``×V`` contract scaling vs the kernel ``÷V``
        Cartesian convention)."""
        visit, st, _, _ = _slab_visit(cell_idx=2, n_ord=4)
        assert st.abs_mu is not None and st.volume is not None
        mu, h = float(st.abs_mu), float(st.volume)
        sig = np.array([1.2, 0.7])
        q_bar = np.array([2.0, 0.5])
        psi_in = np.array([0.3, 0.1])
        strat = LinearDiscontinuous()
        # group 2 (÷V kernel): s = 2|μ|/h, Q_cells = Q̄, flat (Q̂=0)
        psi_avg2, (psi_out2,) = strat.cell_kernel_batch(
            psi_in=(psi_in[None, :, None],),
            s_axes=(np.array([[[2.0 * mu / h]]]),),
            sigt_cells=sig[:, None], Q_cells=q_bar[None, :, None],
        )
        # group 1 (×V per-cell): source = (Q̄·h, 0) → flat
        res1 = strat.update(
            visit, sig, _moment_source(q_bar * h, np.zeros(2)),
            UpstreamState(spatial_upstream=psi_in, angular_upstream=None),
        )
        np.testing.assert_allclose(
            res1.cell_average_flux, psi_avg2.ravel(), rtol=1e-12, atol=1e-13,
        )
        np.testing.assert_allclose(
            res1.outgoing_spatial_flux, psi_out2.ravel(), rtol=1e-12, atol=1e-13,
        )

    @pytest.mark.foundation
    def test_group3_equals_group2_scan_flat(self) -> None:
        r"""Affine-scan coefficients (×V, group 3) ≡ batched kernel (÷V, group 2).

        #158 Inc B sign-trap catcher: LD's ``affine_scan_coefficients`` builds
        the ×V Schur ``(a, inverse_denom, w)`` — a *genuinely different formula*
        from the ÷V ``cell_kernel_batch`` (the LM-1989 §1.4/§6 slope-row sign
        trap lives in the ×V transmission ``a``).  Applying the generic
        coefficient-model ops (``ψ_out = a·ψ_in + b``, ``b = QV·inv/w``;
        ``ψ̄ = (1−w)ψ_in + w·ψ_out``) must reproduce the trusted Increment-A
        kernel cell-for-cell.  Single cell, flat source (Q̂ = 0), 2G.
        """
        _, st, _, _ = _slab_visit(cell_idx=2, n_ord=4)
        assert st.abs_mu is not None and st.volume is not None
        mu, h = float(st.abs_mu), float(st.volume)
        sig = np.array([1.2, 0.7])
        q_bar = np.array([2.0, 0.5])
        psi_in = np.array([0.3, 0.1])
        strat = LinearDiscontinuous()
        # group 2 (÷V kernel) — the trusted Increment-A reference path.
        psi_avg2, (psi_out2,) = strat.cell_kernel_batch(
            psi_in=(psi_in[None, :, None],),
            s_axes=(np.array([[[2.0 * mu / h]]]),),
            sigt_cells=sig[:, None], Q_cells=q_bar[None, :, None],
        )
        # group 3 (×V coefficients) + the generic affine_closure ops.
        a, inv, w = strat.affine_scan_coefficients(
            abs_mu=np.array([mu]), A_down=np.array([[1.0]]),
            A_total=np.array([[2.0]]), dA_w=np.array([[0.0]]),
            c_out=np.array([[0.0]]), V=np.array([[h]]),
            sig_t=sig[None, :, None],
        )                                                    # each (1, 2, 1)
        psi_in_b = psi_in[None, :, None]                     # (1, 2, 1)
        qv = (q_bar * h)[None, :, None]                      # ×V volumetric source
        psi_out3 = a * psi_in_b + source_emission(qv, inv, w)
        psi_avg3 = cell_average(psi_in_b, psi_out3, w)
        np.testing.assert_allclose(psi_avg3, psi_avg2, rtol=1e-12, atol=1e-13)
        np.testing.assert_allclose(psi_out3, psi_out2, rtol=1e-12, atol=1e-13)

    @pytest.mark.foundation
    def test_cell_kernel_batch_rejects_multi_d(self) -> None:
        """A multi-D LD is bilinear (deferred #158 Inc. D) → d≥2 raises."""
        strat = LinearDiscontinuous()
        two = np.zeros((1, 1, 1))
        with pytest.raises(NotImplementedError, match="d=1"):
            strat.cell_kernel_batch(
                psi_in=(two, two), s_axes=(two, two),
                sigt_cells=np.ones((1, 1)), Q_cells=np.zeros((1, 1, 1)),
            )


def test_ld_returns_cell_result() -> None:
    """update returns a CellResult (slab: outgoing_angular_state is None)."""
    visit, _, _, _ = _slab_visit(cell_idx=1, n_ord=4)
    sig = np.array([1.1])
    source = _moment_source(np.array([0.9]), np.array([0.2]))
    upstream = UpstreamState(
        spatial_upstream=np.array([0.3]), angular_upstream=None,
    )
    result = LinearDiscontinuous().update(visit, sig, source, upstream)
    assert isinstance(result, CellResult)
    assert result.outgoing_angular_state is None
