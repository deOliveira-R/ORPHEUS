r"""Branch-2 gate — the production (numpy) d-generic UBLD primitive + d=1 link.

Pins the numpy production primitive
(:mod:`orpheus.sn.spatial._ubld`) — sub-step **D5b-S1 Branch 2** of issues
#240 / #38 / #37 — against the SymPy algebra-of-record
(:mod:`orpheus.derivations.discrete.sn.ld_ubld`, the Branch-1 reference) and
against the LIVE production
:class:`~orpheus.sn.spatial.LinearDiscontinuous` scheme.

The three deliverable gates (the brief):

* **Primitive unit tests** — the numpy ``assemble_ubld`` / ``per_cell_solve``
  reproduce the symbolic dense system at ``d = 1`` (matrices + solved moments)
  AND are exact on a bilinear flux at ``d = 2`` (the ``xy`` cross-coupling
  exercised — the multi-D analog of the 1-D linear-exactness oracle).  This is
  the structural-independence check: the numpy primitive is a SEPARATE
  implementation, not a re-typed copy of the symbolic one (they share only
  ``numpy`` / ``sympy`` primitives, below the trusted-library line).
* **The shared d=1 closed form == the dense primitive's d=1 reduction** — the
  vectorized ``d1_closed_form`` (the fast path) reproduces the batched dense
  ``per_cell_solve`` at ``d = 1`` to machine ε, in all three views (÷V kernel,
  ×V scan, ×V per-cell Schur).  This is the elegance-enforcer's Branch-1
  CONCERN closed IN CODE: the symbolic oracle proved the three views equal the
  reduction symbolically; this proves it numerically, end to end.
* **The link proof** — the LIVE production
  :meth:`~orpheus.sn.spatial.LinearDiscontinuous.update` /
  :meth:`~orpheus.sn.spatial.LinearDiscontinuous.cell_kernel_batch` /
  :meth:`~orpheus.sn.spatial.LinearDiscontinuous.affine_scan_coefficients`
  (which now single-source through ``d1_closed_form``) reproduce the dense
  primitive's ``d = 1`` solve — so the production scheme's three views are
  anchored to the d-generic primitive, not only to each other.

``-O``-safe (Mode 8): the load-bearing checks are FUNCTION CALLS
(``np.testing`` / ``pytest.fail``), never bare ``assert`` (ORPHEUS runs
``python -O``, which strips bare asserts to no-ops).
"""

from __future__ import annotations

import numpy as np
import pytest
import sympy as sp

from orpheus.derivations.discrete.sn import ld_ubld as sym
from orpheus.geometry import BC, CoordSystem, Mesh1D, slab_streaming
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.spatial import LinearDiscontinuous, UpstreamState
from orpheus.sn.spatial._ubld import (
    assemble_inflow_axis,
    assemble_ubld,
    d1_closed_form,
    per_cell_solve,
)
from orpheus.sn.spatial.scheme import CellVisit, DiscretizationSchemeBase

THETA = 1.0 / 3.0

cell_average = DiscretizationSchemeBase.cell_average
source_emission = DiscretizationSchemeBase.source_emission
outgoing_face_from_average = DiscretizationSchemeBase.outgoing_face_from_average


# ─────────────────────────────────────────────────────────────────────────
# -O-safe gate helper (Mode 8)
# ─────────────────────────────────────────────────────────────────────────


def _require(condition: bool, message: str) -> None:
    """Fail the test if ``condition`` is falsy (fires under ``python -O``)."""
    if not condition:
        pytest.fail(message)


# ═══════════════════════════════════════════════════════════════════════
# Primitive unit tests — numpy assemble/solve == symbolic, exact-on-bilinear
# ═══════════════════════════════════════════════════════════════════════


class TestPrimitiveMatchesSymbolic:
    r"""The numpy d-generic primitive reproduces the SymPy algebra-of-record."""

    @pytest.mark.foundation
    def test_d1_assembled_matrices_match_symbolic(self) -> None:
        """The numpy d=1 A / M / G / F_out equal the symbolic ones (small input)."""
        mu_v, h_v, sig_v = 0.7, 0.5, 1.3
        mu, h, sig_t = sp.symbols("mu h Sigma_t", positive=True)
        sasm = sym.assemble_ubld([h], [mu], sig_t, sym.THETA)
        subs = {mu: mu_v, h: h_v, sig_t: sig_v, sym.THETA: sp.Rational(1, 3)}
        nasm = assemble_ubld([np.array(h_v)], [np.array(mu_v)], np.array(sig_v), THETA)
        for key in ("A", "M", "G", "F_out"):
            sym_mat = np.array(sasm[key].subs(subs)).astype(float)
            np.testing.assert_allclose(
                nasm[key], sym_mat, rtol=1e-13, atol=1e-14,
                err_msg=f"numpy {key} != symbolic {key}",
            )

    @pytest.mark.foundation
    def test_d1_solved_moments_match_symbolic(self) -> None:
        """The numpy d=1 per_cell_solve reproduces the symbolic dense solve."""
        mu_v, h_v, sig_v, Q_v, pin_v = 0.7, 0.5, 1.3, 1.1, 0.25
        # Symbolic dense (flat Q̂ = 0): R = M·[Q̄,0] + |μ|·[1,-1]·ψ_in.
        mu, h, sig_t = sp.symbols("mu h Sigma_t", positive=True)
        Qbar, psi_in = sp.symbols("Qbar psi_in", real=True)
        sasm = sym.assemble_ubld([h], [mu], sig_t, sym.THETA)
        R = sasm["M"] * sp.Matrix([Qbar, 0]) + mu * sym.fin_trace_weight() * psi_in
        psi = sym.per_cell_solve(sasm, R)
        subs = {
            mu: mu_v, h: h_v, sig_t: sig_v, Qbar: Q_v, psi_in: pin_v,
            sym.THETA: sp.Rational(1, 3),
        }
        psi_bar_sym = float(psi[0].subs(subs))
        psi_out_sym = float((psi[0] + psi[1]).subs(subs))

        # numpy dense (same flat-source moment convention).
        nasm = assemble_ubld(
            [np.array(h_v)], [np.array(mu_v)], np.array(sig_v), THETA,
        )
        Svec = np.array([Q_v, 0.0])            # per-volume source moments
        R_n = np.einsum("...ij,...j->...i", nasm["M"], Svec)
        R_n = R_n + mu_v * np.array([1.0, -1.0]) * pin_v
        psi_n = per_cell_solve(nasm, R_n)
        np.testing.assert_allclose(psi_n[0], psi_bar_sym, rtol=1e-13, atol=1e-14)
        np.testing.assert_allclose(
            psi_n[0] + psi_n[1], psi_out_sym, rtol=1e-13, atol=1e-14,
        )

    @pytest.mark.foundation
    def test_d2_assembled_matrices_match_symbolic(self) -> None:
        r"""The numpy d=2 A / M / G / F_out equal the symbolic ones ENTRY-WISE.

        D5b-S2 carry-forward (the S1 hindsight audit): at d=1 the matrix-equality
        pin (``test_d1_assembled_matrices_match_symbolic``) fires, but the d≥2
        numpy CELL assembly was pinned ONLY by the physics exact-on-bilinear
        oracle.  An ENTRY-WISE ``A == A`` (+ M / G / F_out) at d=2 is the
        structural pin that catches a dropped / mis-scaled factor in the numpy
        Kronecker CELL assembly directly (e.g. the streaming factor in ``G`` or a
        θ-weight in ``M``), not only through a solved-flux coincidence.

        Scope (why NO ``catches("ERR-060")`` marker): ``assemble_ubld`` builds
        the CELL matrices only — it does NOT carry the upwind-inflow ``|μ_axis|``
        factor (that lives in ``assemble_inflow_axis``).  A dropped ``|μ_axis|``
        inflow factor (ERR-060) leaves A/M/G/F_out unchanged → this pin stays
        GREEN; the inflow factor is caught by ``test_d2_exact_on_bilinear``
        (which routes through the inflow assembler).  Marking this pin
        ``catches("ERR-060")`` would be a coverage overclaim."""
        mx, my, hx_v, hy_v, sig_v = 0.6, 0.4, 0.5, 0.7, 1.2
        mu_x, mu_y, hx, hy, sig_t = sp.symbols(
            "mu_x mu_y h_x h_y Sigma_t", positive=True,
        )
        sasm = sym.assemble_ubld([hx, hy], [mu_x, mu_y], sig_t, sym.THETA)
        subs = {
            mu_x: mx, mu_y: my, hx: hx_v, hy: hy_v, sig_t: sig_v,
            sym.THETA: sp.Rational(1, 3),
        }
        nasm = assemble_ubld(
            [np.array(hx_v), np.array(hy_v)],
            [np.array(mx), np.array(my)],
            np.array(sig_v), THETA,
        )
        for key in ("A", "M", "G", "F_out"):
            sym_mat = np.array(sasm[key].subs(subs)).astype(float)
            np.testing.assert_allclose(
                nasm[key], sym_mat, rtol=1e-13, atol=1e-14,
                err_msg=f"numpy d=2 {key} != symbolic {key} (entry-wise)",
            )

    @pytest.mark.foundation
    def test_d1_assemble_is_batched(self) -> None:
        """The d=1 assembler/solve vectorize over a (N_oct, ng, n_diag) stack."""
        rng = np.random.default_rng(240)
        shape = (3, 2, 4)
        mu = np.broadcast_to(rng.uniform(0.1, 0.9, (3,))[:, None, None], shape)
        h = np.broadcast_to(rng.uniform(0.2, 1.5, (4,))[None, None, :], shape)
        sig = rng.uniform(0.3, 2.0, shape)
        Qbar = rng.uniform(0.2, 2.0, shape)
        psi_in = rng.uniform(0.05, 0.6, shape)
        asm = assemble_ubld([h], [mu], sig, THETA)
        _require(asm["A"].shape == shape + (2, 2), f"batched A wrong shape: {asm['A'].shape}")
        Svec = np.stack([Qbar, np.zeros_like(Qbar)], axis=-1)
        R = np.einsum("...ij,...j->...i", asm["M"], Svec)
        R = R + mu[..., None] * np.array([1.0, -1.0]) * psi_in[..., None]
        psi = per_cell_solve(asm, R)
        _require(psi.shape == shape + (2,), f"batched solve wrong shape: {psi.shape}")

    @pytest.mark.foundation
    @pytest.mark.catches("ERR-060")
    def test_d2_exact_on_bilinear(self) -> None:
        r"""The numpy d=2 UBLD recovers ψ = a + bx + cy + dxy exactly (xy on).

        The multi-D analog of the 1-D linear-exactness oracle, in numpy: feed
        the DG-exact upstream-x / upstream-y face moments + projected source
        moments for a bilinear ψ; the solved ``2² = 4`` tensor-Legendre moments
        must equal the EXACT projections.  This is the structurally-independent
        correctness gate for the d-generic primitive (it catches a missing
        ``|μ_axis|`` streaming factor on the multi-axis inflow — ERR-060 — which
        the d=1 paths are blind to).
        """
        a, b, c, d = 0.7, 1.3, -0.4, 0.9
        mu_x, mu_y, sig = 0.6, 0.5, 1.2
        hx, hy = 0.5, 0.7
        x0, y0 = 0.3, 0.2
        xmx, xmy = x0 + hx / 2.0, y0 + hy / 2.0

        # EXACT tensor-Legendre moments [bar, ŷ, x̂, x̂y] (kron x-outer, y-inner).
        pbar = a + b * xmx + c * xmy + d * xmx * xmy
        phat_x = b * hx / 2.0 + d * (hx / 2.0) * xmy
        phat_y = c * hy / 2.0 + d * xmx * (hy / 2.0)
        phat_xy = d * (hx / 2.0) * (hy / 2.0)
        psi_exact = np.array([pbar, phat_y, phat_x, phat_xy])

        # Manufactured source Q = Ω·∇ψ + Σ_t ψ; closed-form moment projection.
        xx, yy = sp.symbols("x y", real=True)
        Qfield = (
            mu_x * (b + d * yy)
            + mu_y * (c + d * xx)
            + sig * (a + b * xx + c * yy + d * xx * yy)
        )

        def moment(expr: sp.Expr, ox: int, oy: int) -> float:
            Lx = 2 * (xx - xmx) / hx
            Ly = 2 * (yy - xmy) / hy
            basis = (Lx**ox) * (Ly**oy)
            norm = sp.Rational(1, 1)
            if ox == 1:
                norm *= sp.Rational(1, 3)
            if oy == 1:
                norm *= sp.Rational(1, 3)
            integ = sp.integrate(
                sp.integrate(expr * basis, (xx, x0, x0 + hx)), (yy, y0, y0 + hy)
            ) / (hx * hy)
            return float(sp.simplify(integ / norm))

        S_moments = np.array([
            moment(Qfield, 0, 0), moment(Qfield, 0, 1),
            moment(Qfield, 1, 0), moment(Qfield, 1, 1),
        ])

        def y_face(expr_y: sp.Expr) -> np.ndarray:
            m0 = float(sp.integrate(expr_y, (yy, y0, y0 + hy)) / hy)
            m1 = float(
                sp.integrate(expr_y * (2 * (yy - xmy) / hy), (yy, y0, y0 + hy))
                / hy / sp.Rational(1, 3)
            )
            return np.array([m0, m1])

        def x_face(expr_x: sp.Expr) -> np.ndarray:
            m0 = float(sp.integrate(expr_x, (xx, x0, x0 + hx)) / hx)
            m1 = float(
                sp.integrate(expr_x * (2 * (xx - xmx) / hx), (xx, x0, x0 + hx))
                / hx / sp.Rational(1, 3)
            )
            return np.array([m0, m1])

        trace_x = a + b * x0 + c * yy + d * x0 * yy   # ψ(x0, y) upstream-x
        trace_y = a + b * xx + c * y0 + d * xx * y0   # ψ(x, y0) upstream-y
        hs = [np.array(hx), np.array(hy)]
        mus = [np.array(mu_x), np.array(mu_y)]
        asm = assemble_ubld(hs, mus, np.array(sig), THETA)
        R_moment = np.einsum("...ij,...j->...i", asm["M"], S_moments)
        Fin_x = assemble_inflow_axis(hs, mus, 0, y_face(trace_x), THETA)
        Fin_y = assemble_inflow_axis(hs, mus, 1, x_face(trace_y), THETA)
        psi_solved = per_cell_solve(asm, R_moment + Fin_x + Fin_y)
        np.testing.assert_allclose(psi_solved, psi_exact, rtol=1e-12, atol=1e-13)


# ═══════════════════════════════════════════════════════════════════════
# The shared d=1 closed form == the dense primitive's d=1 reduction
# (the elegance-enforcer's Branch-1 CONCERN, closed in code)
# ═══════════════════════════════════════════════════════════════════════


def _dense_d1_reference(
    mu: np.ndarray, h: np.ndarray, sig: np.ndarray, Qbar: np.ndarray, psi_in: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Batched dense d=1 solve (flat Q̂ = 0) → (ψ̄, ψ_out) — the reference."""
    asm = assemble_ubld([h], [mu], sig, THETA)
    Svec = np.stack([Qbar, np.zeros_like(Qbar)], axis=-1)
    R = np.einsum("...ij,...j->...i", asm["M"], Svec)
    R = R + mu[..., None] * np.array([1.0, -1.0]) * psi_in[..., None]
    psi = per_cell_solve(asm, R)
    return psi[..., 0], psi[..., 0] + psi[..., 1]


class TestClosedFormEqualsDenseReduction:
    r"""``d1_closed_form`` (the fast path) == the dense primitive's d=1 solve."""

    @pytest.mark.foundation
    def test_divV_kernel_view_equals_dense(self) -> None:
        """The ÷V kernel view (eff_denom / kernel_rhs / w) == the dense solve."""
        rng = np.random.default_rng(1)
        shape = (3, 2, 4)
        mu = np.broadcast_to(rng.uniform(0.1, 0.9, (3,))[:, None, None], shape)
        h = np.broadcast_to(rng.uniform(0.2, 1.5, (4,))[None, None, :], shape)
        sig = rng.uniform(0.3, 2.0, shape)
        Qbar = rng.uniform(0.2, 2.0, shape)
        psi_in = rng.uniform(0.05, 0.6, shape)
        pbar_d, pout_d = _dense_d1_reference(mu, h, sig, Qbar, psi_in)

        cf = d1_closed_form(mu / h, sig, THETA)         # g = |μ|A_down/V = μ/h
        pbar = cf.kernel_rhs(Qbar, psi_in) / cf.eff_denom
        pout = outgoing_face_from_average(pbar, psi_in, cf.w)
        np.testing.assert_allclose(pbar, pbar_d, rtol=1e-12, atol=1e-13)
        np.testing.assert_allclose(pout, pout_d, rtol=1e-12, atol=1e-13)

    @pytest.mark.foundation
    def test_timesV_scan_view_equals_dense(self) -> None:
        """The ×V scan view (a, inverse_denom, w) == the dense solve."""
        rng = np.random.default_rng(2)
        shape = (3, 2, 4)
        mu = np.broadcast_to(rng.uniform(0.1, 0.9, (3,))[:, None, None], shape)
        h = np.broadcast_to(rng.uniform(0.2, 1.5, (4,))[None, None, :], shape)
        sig = rng.uniform(0.3, 2.0, shape)
        Qbar = rng.uniform(0.2, 2.0, shape)
        psi_in = rng.uniform(0.05, 0.6, shape)
        pbar_d, pout_d = _dense_d1_reference(mu, h, sig, Qbar, psi_in)

        cf = d1_closed_form(mu / h, sig, THETA)
        a, inv, w = cf.scan_xV(h)
        pout = a * psi_in + source_emission(Qbar * h, inv, w)   # ψ_out = a·ψ_in + b
        pbar = cell_average(psi_in, pout, w)
        np.testing.assert_allclose(pbar, pbar_d, rtol=1e-12, atol=1e-13)
        np.testing.assert_allclose(pout, pout_d, rtol=1e-12, atol=1e-13)
        _require(
            bool(np.all(np.isfinite(a))),
            "scan transmission a must be finite (source-independent)",
        )

    @pytest.mark.foundation
    def test_xV_schur_view_equals_dense(self) -> None:
        """The ×V per-cell Schur view (schur_xV) == the dense solve, incl. slope."""
        rng = np.random.default_rng(3)
        shape = (2, 3)
        mu = rng.uniform(0.1, 0.9, shape)
        h = rng.uniform(0.2, 1.5, shape)
        sig = rng.uniform(0.3, 2.0, shape)
        Qbar = rng.uniform(0.2, 2.0, shape)
        psi_in = rng.uniform(0.05, 0.6, shape)
        pbar_d, pout_d = _dense_d1_reference(mu, h, sig, Qbar, psi_in)

        cf = d1_closed_form(mu / h, sig, THETA)
        # ×V per-cell contract: s_bar = Q̄·h, s_hat = 0 (flat).
        S, eff_source, eff_numer, slope_source, mu_Adown, d2p = cf.schur_xV(
            h, Qbar * h, np.zeros_like(Qbar), psi_in,
        )
        psi_bar = (eff_source + eff_numer) / S
        psi_hat = (mu_Adown * psi_bar + slope_source - mu_Adown * psi_in) / d2p
        psi_out = psi_bar + psi_hat
        np.testing.assert_allclose(psi_bar, pbar_d, rtol=1e-12, atol=1e-13)
        np.testing.assert_allclose(psi_out, pout_d, rtol=1e-12, atol=1e-13)


# ═══════════════════════════════════════════════════════════════════════
# The link proof — the LIVE production scheme == the dense primitive's d=1
# ═══════════════════════════════════════════════════════════════════════


class TestProductionViewsAnchoredToPrimitive:
    r"""The production scheme's three views (single-sourced through the helper)
    reproduce the dense primitive's d=1 solve, end to end."""

    @staticmethod
    def _slab_cell():
        """A concrete slab cell + ordinate from the production factory."""
        mesh = Mesh1D(
            edges=np.linspace(0.0, 1.0, 6),
            mat_ids=np.zeros(5, dtype=int),
            coord=CoordSystem.CARTESIAN,
            bc_left=BC("vacuum"),
            bc_right=BC("vacuum"),
        )
        quad = Quadrature.gauss_legendre(4)
        op = slab_streaming(mesh, quad)
        cell_idx, direction_idx = 2, quad.N - 1     # most-positive ordinate
        st = op.streaming_terms(cell_idx, direction_idx)
        visit = CellVisit(
            cell_idx=cell_idx, streaming_terms=st, face_area_downstream=1.0,
        )
        return visit, float(st.abs_mu), float(st.volume)

    @pytest.mark.foundation
    def test_production_update_equals_dense(self) -> None:
        """LinearDiscontinuous.update (×V per-cell) == the dense primitive d=1."""
        visit, mu_v, h_v = self._slab_cell()
        sig = np.array([1.2, 0.7])              # 2G heterogeneous (multi-group teeth)
        q_bar = np.array([2.0, 0.5])
        psi_in = np.array([0.3, 0.1])

        # Dense primitive d=1 reference at these numbers.
        pbar_d, pout_d = _dense_d1_reference(
            np.full(2, mu_v), np.full(2, h_v), sig, q_bar, psi_in,
        )
        # Production update (×V contract: source = (Q̄·h, 0) → flat).
        source = np.stack([q_bar * h_v, np.zeros(2)], axis=0)
        res = LinearDiscontinuous().update(
            visit, sig, source, UpstreamState(spatial_upstream=psi_in),
        )
        np.testing.assert_allclose(res.cell_average_flux, pbar_d, rtol=1e-12, atol=1e-13)
        np.testing.assert_allclose(
            res.outgoing_spatial_flux, pout_d, rtol=1e-12, atol=1e-13,
        )

    @pytest.mark.foundation
    def test_production_kernel_equals_dense(self) -> None:
        """LinearDiscontinuous.cell_kernel_batch (÷V) == the dense primitive d=1."""
        _, mu_v, h_v = self._slab_cell()
        sig = np.array([1.2, 0.7])
        q_bar = np.array([2.0, 0.5])
        psi_in = np.array([0.3, 0.1])
        pbar_d, pout_d = _dense_d1_reference(
            np.full(2, mu_v), np.full(2, h_v), sig, q_bar, psi_in,
        )
        psi_avg, (psi_out,) = LinearDiscontinuous().cell_kernel_batch(
            psi_in=(psi_in[None, :, None],),
            s_axes=(np.array([[[mu_v / h_v]]]),),
            sigt_cells=sig[:, None], Q_cells=q_bar[None, :, None],
        )
        np.testing.assert_allclose(psi_avg.ravel(), pbar_d, rtol=1e-12, atol=1e-13)
        np.testing.assert_allclose(psi_out.ravel(), pout_d, rtol=1e-12, atol=1e-13)

    @pytest.mark.foundation
    def test_production_scan_equals_dense(self) -> None:
        """affine_scan_coefficients (×V scan) + base reconstruction == dense d=1."""
        _, mu_v, h_v = self._slab_cell()
        sig = np.array([1.2, 0.7])
        q_bar = np.array([2.0, 0.5])
        psi_in = np.array([0.3, 0.1])
        pbar_d, pout_d = _dense_d1_reference(
            np.full(2, mu_v), np.full(2, h_v), sig, q_bar, psi_in,
        )
        a, inv, w = LinearDiscontinuous().affine_scan_coefficients(
            abs_mu=np.array([mu_v]), A_down=np.array([[1.0]]),
            A_total=np.array([[2.0]]), dA_w=np.array([[0.0]]),
            c_out=np.array([[0.0]]), V=np.array([[h_v]]),
            sig_t=sig[None, :, None],
        )                                                    # each (1, 2, 1)
        psi_in_b = psi_in[None, :, None]
        qv = (q_bar * h_v)[None, :, None]
        psi_out = a * psi_in_b + source_emission(qv, inv, w)
        psi_bar = cell_average(psi_in_b, psi_out, w)
        np.testing.assert_allclose(psi_bar.ravel(), pbar_d, rtol=1e-12, atol=1e-13)
        np.testing.assert_allclose(psi_out.ravel(), pout_d, rtol=1e-12, atol=1e-13)
