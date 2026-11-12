r"""Branch-1 (SymPy) foundation gate — the d-generic UBLD symbolic reference.

Pins the symbolic algebra-of-record for the dimension-generic Unlumped
BiLinear Discontinuous (UBLD) per-cell Galerkin system
(:mod:`orpheus.derivations.discrete.sn.ld_ubld`) — sub-step D5b-S1 Branch 1
of issues #240 / #38 / #37.  One ``@pytest.mark.foundation`` test per
``derive_*()`` claim (the test count equals the V_n claim count); NO
``verifies(...)`` (these are software/algebra invariants, not L0/L1
solver-equation claims — ``vv-principles`` "foundation" level).

The two oracles (per the brief):

* **Oracle (i) — d=1 reduction.**  The assembled d=1 UBLD ``A`` / ``R`` /
  Schur ``S`` / slope ``ψ̂`` / ``D₂'`` reduce EXACTLY to the production
  slab LD closed forms, AND the production ÷V (``_kernel_terms``) and ×V
  (``affine_scan_coefficients``) views equal that reduction (the
  "single-source the math" proof for Branch 2).  The d=1 test ALSO probes
  the production :class:`~orpheus.transport.spatial.LinearDiscontinuous` ``update``
  numerically at the SymPy-derived flat-source reduction, so the gate ties
  the symbolic primitive to the live production algebra (not just to a
  re-typed copy of its closed form).
* **Oracle (ii) — exact-on-bilinear.**  The d=2 UBLD recovers ANY bilinear
  flux ``ψ = a + bx + cy + dxy`` exactly (the multi-D analog of the 1-D
  "exact on linear-in-x" oracle), the ``xy`` cross moment exercised.

This is the same epistemic role as ``test_linear_discontinuous.py``'s
``TestLDLinearExactness`` (the 1-D linear-exactness oracle), lifted to the
d-generic symbolic primitive that Branch 2 will be verified against.
"""

from __future__ import annotations

import numpy as np
import pytest
import sympy as sp

from orpheus.derivations.discrete.sn.ld_ubld import (
    THETA,
    assemble_ubld,
    derive_d1_kernel_view_equals,
    derive_curvilinear_ld_cell_gram,
    derive_d1_coordinate_pairing_is_the_integral,
    derive_d1_integrated_slope_source_is_the_galerkin_moment,
    derive_d1_reduction_to_production,
    derive_d1_scan_view_equals,
    derive_d1_transpose_equals_At_Minv,
    derive_d2_exact_on_bilinear,
    derive_d3_assembles,
    derive_interval_legendre_gram,
    derive_octant_frame_sign_is_involution,
    fin_trace_weight,
    per_cell_solve,
)
from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from tests.gates.transport._phase_space_cases import slab_streaming_terms
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.spatial import LinearDiscontinuous, UpstreamState
from orpheus.transport.spatial._ubld import d1_closed_form
from orpheus.transport.spatial.scheme import CellVisit


# ─────────────────────────────────────────────────────────────────────────
# -O-safe gate helpers (Mode 8): ORPHEUS runs ``python -O``, which strips
# bare ``assert`` to a no-op.  These are FUNCTION CALLS (``pytest.fail`` /
# ``np.testing``), so the verification fires under ``-O`` — never a bare
# ``assert`` on the load-bearing symbolic-zero / boolean claims.
# ─────────────────────────────────────────────────────────────────────────


def _require(condition: bool, message: str) -> None:
    """Fail the test if ``condition`` is falsy (fires under ``python -O``)."""
    if not condition:
        pytest.fail(message)


def _require_zero(expr: sp.Expr, name: str) -> None:
    """Fail unless ``expr`` is the symbolic zero (fires under ``python -O``)."""
    if sp.simplify(expr) != 0:
        pytest.fail(f"{name} is not symbolically zero: {expr}")


def _require_zero_matrix(mat: sp.Matrix, name: str) -> None:
    """Fail unless ``mat`` is the symbolic zero matrix (fires under ``-O``)."""
    if not mat.is_zero_matrix:
        pytest.fail(f"{name} is not the zero matrix: {mat}")


# ═══════════════════════════════════════════════════════════════════════
# Oracle (i) — d=1 reduction to the production slab LD
# ═══════════════════════════════════════════════════════════════════════


class TestOracleId1Reduction:
    r"""The assembled d=1 UBLD == the production slab LD (algebra + views)."""

    @pytest.mark.foundation
    def test_d1_reduction_to_production_schur(self) -> None:
        """V_d1: A / R / S / ψ̂ / D₂' reduce to the production closed forms."""
        result = derive_d1_reduction_to_production()
        _require(result["pass"], f"V_d1 failed: {result}")
        # Each reported difference must be the symbolic zero.
        _require_zero_matrix(result["diff_A"], "diff_A")
        _require_zero_matrix(result["diff_R"], "diff_R")
        _require_zero(result["diff_face"], "diff_face")
        _require_zero(result["diff_psi_bar"], "diff_psi_bar")
        _require_zero(result["diff_psi_hat"], "diff_psi_hat")

    @pytest.mark.foundation
    def test_d1_divV_kernel_view_equals_reduction(self) -> None:
        """V_d1_kernel: the production ÷V _kernel_terms form == d=1 (flat Q̂=0)."""
        result = derive_d1_kernel_view_equals()
        _require(result["pass"], f"V_d1_kernel failed: {result}")
        _require_zero(result["diff_psi_bar"], "diff_psi_bar")
        _require_zero(result["diff_psi_out"], "diff_psi_out")

    @pytest.mark.foundation
    def test_d1_timesV_scan_view_equals_reduction(self) -> None:
        """V_d1_scan: the production ×V affine_scan_coefficients form == d=1."""
        result = derive_d1_scan_view_equals()
        _require(result["pass"], f"V_d1_scan failed: {result}")
        _require_zero(result["diff_psi_bar"], "diff_psi_bar")
        _require_zero(result["diff_psi_out"], "diff_psi_out")
        _require(result["a_source_independent"], "scan transmission a is not source-independent")

    @pytest.mark.foundation
    def test_d1_symbolic_primitive_matches_production_update(self) -> None:
        r"""Tie the symbolic d=1 primitive to the LIVE production LD.

        The symbolic oracle proves the primitive equals the production
        *closed form*.  This test closes the loop by evaluating the symbolic
        d=1 ``ψ̄`` / ``ψ_out`` (flat ``Q̂ = 0``) at concrete numbers and
        asserting the production
        :meth:`~orpheus.transport.spatial.LinearDiscontinuous.update` reproduces
        them — so the gate is anchored to the running algebra, not only to a
        re-typed copy of its formula (structural-independence of the cross-check).
        """
        # Concrete slab cell + ordinate: the slab terms by hand (the SN operator feeds the same; see slab_streaming_terms).
        mesh = Mesher(StructuredGeometry.from_homogeneous(1.0, 0)).partition(
            CellsByCount.uniform_width(5),
        ).mesh
        quad = Quadrature.gauss_legendre(4)
        cell_idx, direction_idx = 2, quad.n_ordinates - 1     # most-positive ordinate
        st = slab_streaming_terms(mesh, quad, cell_idx, direction_idx)
        _require(
            st.abs_mu is not None and st.volume is not None,
            "slab StreamingTerms must populate abs_mu and volume",
        )
        mu_v = float(st.abs_mu)
        h_v = float(st.volume)
        visit = CellVisit(
            cell_idx=cell_idx, streaming_terms=st, face_area_downstream=1.0,
        )

        # Concrete numbers (2-group, heterogeneous Σ_t — multi-group teeth).
        sig_v = np.array([1.2, 0.7])
        q_bar_v = np.array([2.0, 0.5])
        psi_in_v = np.array([0.3, 0.1])

        # Symbolic d=1 primitive, FLAT source (Q̂ = 0), evaluated at the numbers.
        mu, h, sig_t = sp.symbols("mu h Sigma_t", positive=True)
        Qbar, psi_in = sp.symbols("Qbar psi_in", real=True)
        asm = assemble_ubld([h], [mu], sig_t, THETA)
        R = asm["M"] * sp.Matrix([Qbar, 0]) + mu * fin_trace_weight() * psi_in
        psi = per_cell_solve(asm, R)
        subs_common = {mu: mu_v, h: h_v, THETA: sp.Rational(1, 3)}
        psi_bar_fn = sp.lambdify((sig_t, Qbar, psi_in), psi[0].subs(subs_common), "numpy")
        psi_out_expr = (psi[0] + psi[1]).subs(subs_common)
        psi_out_fn = sp.lambdify((sig_t, Qbar, psi_in), psi_out_expr, "numpy")
        psi_bar_sym = psi_bar_fn(sig_v, q_bar_v, psi_in_v)
        psi_out_sym = psi_out_fn(sig_v, q_bar_v, psi_in_v)

        # Production update (×V contract: source = (Q̄·h, 0) → flat).
        source = np.stack([q_bar_v * h_v, np.zeros(2)], axis=0)
        upstream = UpstreamState(spatial_upstream=psi_in_v)
        res = LinearDiscontinuous().update(visit, sig_v, source, upstream)

        np.testing.assert_allclose(
            res.cell_average_flux, psi_bar_sym, rtol=1e-12, atol=1e-13,
        )
        np.testing.assert_allclose(
            res.outgoing_spatial_flux, psi_out_sym, rtol=1e-12, atol=1e-13,
        )


# ═══════════════════════════════════════════════════════════════════════
# Oracle (ii) — exact on a bilinear flux (d=2, xy coupling exercised)
# ═══════════════════════════════════════════════════════════════════════


class TestOracleIIBilinearExactness:
    r"""The d=2 UBLD recovers ψ = a + bx + cy + dxy exactly (xy coupling on)."""

    @pytest.mark.foundation
    @pytest.mark.catches("ERR-060")
    def test_d2_exact_on_bilinear(self) -> None:
        result = derive_d2_exact_on_bilinear()
        _require(result["pass"], f"V_d2_bilinear failed: diff = {result['diff']}")
        _require_zero_matrix(result["diff"], "d2 bilinear residual")


# ═══════════════════════════════════════════════════════════════════════
# d=3 structural readiness
# ═══════════════════════════════════════════════════════════════════════


class TestD3StructuralReadiness:
    r"""The d-generic assembler scales to the d=3 trilinear (8×8, θ³ weight)."""

    @pytest.mark.foundation
    def test_d3_assembles_8x8_with_theta_cubed(self) -> None:
        result = derive_d3_assembles()
        _require(result["pass"], f"V_d3 failed: {result}")
        _require(result["size"] == (8, 8), f"d=3 system is not 8×8: {result['size']}")
        # The xyz triple-cross moment carries the θ³ collision diagonal weight.
        _require_zero(result["xyz_collision_weight"] - THETA**3, "xyz θ³ collision weight")


# ═══════════════════════════════════════════════════════════════════════
# The transpose oracles (#310 C2 — spec §3.1, the R1b keystones)
# ═══════════════════════════════════════════════════════════════════════


class TestTransposeOracles:
    r"""The LD cell VJP algebra-of-record (#310 C2 R1b, spec §3.1)."""

    @pytest.mark.foundation
    def test_d1_transpose_equals_At_Minv(self) -> None:
        """V_d1_T: VJP == Aᵀ·M⁻¹ (mass-inverse FIRST) + inflow/face pullbacks.

        The order discriminant must be NONZERO — had ``Aᵀ M⁻¹`` and
        ``M⁻¹ Aᵀ`` commuted, the M-R1b-MASSORDER mutation class would be
        undetectable through this oracle (a toothless gate).
        """
        result = derive_d1_transpose_equals_At_Minv()
        _require(result["pass"], f"V_d1_T failed: {result}")
        _require_zero_matrix(result["diff_mass_order"], "diff_mass_order")
        _require(
            result["order_discriminant_nonzero"],
            "Aᵀ·M⁻¹ == M⁻¹·Aᵀ at generic symbols — the mass-order gate has "
            "no teeth",
        )
        _require_zero_matrix(result["diff_inflow_pullback"], "diff_inflow_pullback")
        _require_zero_matrix(result["diff_face_pullback"], "diff_face_pullback")

    @pytest.mark.foundation
    def test_interval_legendre_gram(self) -> None:
        """V_gram: the Gram 1/(2k+1), the endpoint values (±1)^k, the parity, θ = 1/3 (n = 5)."""
        result = derive_interval_legendre_gram()
        _require(result["pass"], f"V_gram failed: {result}")

    @pytest.mark.foundation
    def test_octant_frame_sign_is_involution(self) -> None:
        """V_frame_T: D_s² = I and (D_s·A·D_s)ᵀ = D_s·Aᵀ·D_s (d=1, d=2 rows)."""
        result = derive_octant_frame_sign_is_involution()
        _require(result["pass"], f"V_frame_T failed: {result}")
        for label, ok in result["checks"].items():
            _require(ok, f"frame-sign law failed at {label}")


# ═══════════════════════════════════════════════════════════════════════
# The data side (#707 Q-D): the integrated chart of the moment rows
# ═══════════════════════════════════════════════════════════════════════

_THIS = "tests/gates/transport/spatial/test_ld_ubld_symbolic.py::TestDataSideIntegratedChart"


class TestDataSideIntegratedChart:
    r"""The integrated (V\*) value of an LD moment row is the cell mass on the coefficients (#707 Q-D)."""

    @pytest.mark.foundation
    def test_integrated_slope_source_is_the_galerkin_moment(self) -> None:
        """V_qint: ∫q b_i = M_cell·q_coef (slab, d=2 V·θ^|o|); weak slope row == scan_reconstruct; V-only negative."""
        result = derive_d1_integrated_slope_source_is_the_galerkin_moment()
        _require(result["pass"], f"V_qint failed: {result}")
        _require_zero_matrix(result["diff_mass_slab"], "diff_mass_slab")
        _require_zero_matrix(result["diff_mass_d2"], "diff_mass_d2")
        _require_zero(result["diff_slope_row_vs_scan_reconstruct"], "diff_slope_row_vs_scan_reconstruct")
        _require(sp.simplify(result["residual_V_only"]) != 0, "the V-only negative leg reads zero: no teeth")

    @pytest.mark.foundation
    def test_coordinate_pairing_is_the_integral(self) -> None:
        """V_pair: ψ_coef·q_int == ∫ψq for any q; the V chart misses 2hψ̂Q̂/3."""
        result = derive_d1_coordinate_pairing_is_the_integral()
        _require(result["pass"], f"V_pair failed: {result}")
        _require_zero(result["diff_any_source"], "diff_any_source")
        _require_zero(result["diff_mass_chart"], "diff_mass_chart")
        _require(sp.simplify(result["miss_V_chart"]) != 0, "the V-chart negative leg reads zero: no teeth")

    @pytest.mark.foundation
    def test_curvilinear_ld_cell_gram(self) -> None:
        """V_gram_curv (forward-looking): sphere pole cell M01/V = 1/2, M11/V = 2/5; full Gram."""
        result = derive_curvilinear_ld_cell_gram()
        _require(result["pass"], f"V_gram_curv failed: {result}")
        _require(result["sphere"]["pole"] == (sp.Rational(1, 2), sp.Rational(2, 5)), "sphere pole-cell Gram moved")
        _require(not LinearDiscontinuous.supports_curvilinear,
                 "LD now admits curvilinear meshes: promote V_gram_curv from forward-looking to a production pin")

    @pytest.mark.l1
    @pytest.mark.rests_on(
        f"{_THIS}::test_integrated_slope_source_is_the_galerkin_moment",
        f"{_THIS}::test_coordinate_pairing_is_the_integral",
    )
    def test_slab_mms_order_in_the_integrated_chart(self) -> None:
        r"""Slab LD MMS, forward and adjoint: the V·θ chart keeps the orders, the V chart loses one.

        Fixture: slab [0, 1], one ordinate μ = 0.6, Σ_t = 1.3, no scattering;
        ψ = 1 + 0.5 sin 3x + 0.3 e^{−x} (inflow ψ(0)); adjoint −μψ†' + Σψ† = r,
        ψ† = sin(2(1−x)) e^{x/2}, ψ†(1) = 0; n ∈ {16, 32, 64}. The global
        integrated Galerkin system is assembled from the SYMBOLIC cell
        (:func:`assemble_ubld`, θ = 1/3, :func:`fin_trace_weight`) and gated
        against the PRODUCTION closed-form march (``d1_closed_form(...).scan_reconstruct``).
        The V chart is the production closure fed ``s_hat = VQ̂/θ``.

        [M] 2026-11-12 orders per halving: V·θ forward/adjoint mean 7.9–8.0 (h³),
        slopes ~4; V chart forward mean ~3.9 (h²), adjoint slope ~2 (h).
        """
        mu, sig, theta = 0.6, 1.3, 1.0 / 3.0
        xg, wg = np.polynomial.legendre.leggauss(12)
        psi_ex = lambda x: 1.0 + 0.5 * np.sin(3.0 * x) + 0.3 * np.exp(-x)
        q_ex = lambda x: mu * (1.5 * np.cos(3.0 * x) - 0.3 * np.exp(-x)) + sig * psi_ex(x)
        ad_ex = lambda x: np.sin(2.0 * (1.0 - x)) * np.exp(0.5 * x)
        r_ex = lambda x: -mu * ((-2.0 * np.cos(2.0 * (1.0 - x)) + 0.5 * np.sin(2.0 * (1.0 - x))) * np.exp(0.5 * x)) + sig * ad_ex(x)

        hs = sp.Symbol("h", positive=True)
        cell = assemble_ubld([hs], [sp.Float(mu)], sp.Float(sig), sp.Rational(1, 3))
        A_cell = sp.lambdify(hs, cell["A"], "numpy")
        M_cell = sp.lambdify(hs, cell["M"], "numpy")
        f_in = np.array(fin_trace_weight(), dtype=float).ravel() * mu  # μ·B(−1)
        b_out = np.ones(2)                                             # B(+1)

        def coefs(f, n):
            h = 1.0 / n
            X = (np.arange(n)[:, None] + 0.5) * h + xg[None, :] * h / 2
            F = f(X)
            return np.stack([(F * wg).sum(1) / 2, (F * wg * xg).sum(1) / 2 / theta], 1), h

        def march(n, Q):
            h = 1.0 / n
            cf = d1_closed_form(mu / h, sig, theta)
            out, pin = np.zeros((n, 2)), psi_ex(0.0)
            for i in range(n):
                out[i] = cf.scan_reconstruct(h, Q[i, 0] * h, Q[i, 1] * h, pin)
                pin = out[i].sum()
            return out

        errs = {k: [] for k in ("fb_t", "ah_t", "ab_t", "fb_v", "ah_v")}
        for n in (16, 32, 64):
            Q, h = coefs(q_ex, n)
            E, _ = coefs(psi_ex, n)
            R, _ = coefs(r_ex, n)
            Ead, _ = coefs(ad_ex, n)
            A = np.zeros((2 * n, 2 * n))
            Md = np.diag(M_cell(h)).astype(float)
            for i in range(n):
                A[2 * i:2 * i + 2, 2 * i:2 * i + 2] = A_cell(h)
                if i:
                    A[2 * i:2 * i + 2, 2 * i - 2:2 * i] = -np.outer(f_in, b_out)
            M = np.tile(Md, n)
            b_in = np.zeros(2 * n); b_in[:2] = f_in * psi_ex(0.0)
            fwd = np.linalg.solve(A, M * Q.ravel() + b_in).reshape(n, 2)
            prod = march(n, Q)
            np.testing.assert_allclose(fwd, prod, rtol=0, atol=1e-12)   # symbolic system == production march
            prod_v = march(n, Q / np.array([1.0, theta]))                # V chart: integrated slope = V·Q̂
            adj = np.linalg.solve(A.T, M * R.ravel()).reshape(n, 2)
            D_V = np.full(2 * n, h)
            A_V = (D_V / M)[:, None] * A
            adj_v = np.linalg.solve(A_V.T, D_V * R.ravel()).reshape(n, 2)
            errs["fb_t"].append(np.abs(prod[:, 0] - E[:, 0]).max())
            errs["ab_t"].append(np.abs(adj[:, 0] - Ead[:, 0]).max())
            errs["ah_t"].append(np.abs(adj[:, 1] - Ead[:, 1]).max())
            errs["fb_v"].append(np.abs(prod_v[:, 0] - E[:, 0]).max())
            errs["ah_v"].append(np.abs(adj_v[:, 1] - Ead[:, 1]).max())
        ratio = {k: np.array(v[:-1]) / np.array(v[1:]) for k, v in errs.items()}
        _require(np.all(ratio["fb_t"] > 7.5), f"V·θ forward mean not h³: {ratio['fb_t']}")
        _require(np.all(ratio["ab_t"] > 7.5), f"V·θ adjoint mean not h³: {ratio['ab_t']}")
        _require(np.all(ratio["ah_t"] > 3.7), f"V·θ adjoint slope not h²: {ratio['ah_t']}")
        # negative legs: the V chart loses an order (the mutation bit, X1)
        _require(np.all(ratio["fb_v"] < 4.5), f"V-chart forward mean kept h³: {ratio['fb_v']}")
        _require(np.all(ratio["ah_v"] < 2.5), f"V-chart adjoint slope kept h²: {ratio['ah_v']}")
