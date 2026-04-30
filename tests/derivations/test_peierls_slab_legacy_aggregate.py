"""Per-face vs single-surface-aggregate parity for slab Peierls primitives.

ERR-033 (Issue #131 follow-up, audit 2026-04-23) — the legacy single-
surface aggregates ``compute_P_esc`` / ``compute_G_bc`` retained the
finite-N Gauss-Legendre branch on the ``len(radii) > 1`` path, while
their per-face counterparts (``compute_P_esc_outer/inner`` and
``compute_G_bc_outer/inner``) shipped the closed-form
``½·E_2(τ_face)`` / ``2·E_2(τ_face)`` identities. The slab angular
integral has the form

    ∫₀¹ f(µ) · exp(-τ/µ) dµ

with τ piecewise-constant in σ_t (and therefore µ-independent), so it
IS closed-form as ``E_n(τ)`` regardless of ``n_regions``. Finite-N GL
leaves a ~4e-3 → 6e-5 relative error decaying only as O(N⁻²).

The per-face decomposition

    compute_P_esc = compute_P_esc_outer + compute_P_esc_inner
    compute_G_bc  = compute_G_bc_outer  + compute_G_bc_inner

is structural — the aggregate must equal the sum to machine precision
once both legs share the same closed-form route. After the ERR-033
fix, this test should pass at ``rel_err < 1e-12``. Without the fix it
would FAIL with the ~4e-3 rel err Issue #131 signature replicated on
the legacy aggregates.

Promoted from
``derivations/diagnostics/diag_l131a_single_surface_slab_multi_region.py``
(2026-04-30 triage).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.peierls_geometry import (
    SLAB_POLAR_1D,
    composite_gl_r,
    compute_G_bc,
    compute_G_bc_inner,
    compute_G_bc_outer,
    compute_P_esc,
    compute_P_esc_inner,
    compute_P_esc_outer,
)


@pytest.mark.foundation
@pytest.mark.catches("ERR-033")
class TestSlabSingleSurfaceVsPerFaceAggregate:
    """Multi-region slab legacy aggregate must equal per-face closed-form sum.

    The per-face sum identity holds by construction (each face gets
    ``½·E_2(τ_face)``), and the aggregate has the same closed-form
    representation for any piecewise-constant σ_t. This is a
    multi-region branch reduction invariant — independent of the
    quadrature order — and is the structural gate for ERR-033.
    """

    @pytest.mark.parametrize(
        "radii, sig_t",
        [
            (np.array([1.0, 2.0]), np.array([0.5, 2.0])),
            (np.array([0.3, 1.0, 1.5]), np.array([1.5, 0.2, 3.0])),
        ],
        ids=["2region_thinThick", "3region_mixed"],
    )
    def test_compute_P_esc_matches_per_face_sum(self, radii, sig_t):
        """``compute_P_esc`` must equal the sum of per-face primitives.

        Before the ERR-033 fix: finite-N GL fall-through gives
        ~4e-3 rel err at N=24. After the fix: closed-form sum,
        ~1e-15 rel err.
        """
        r_nodes, _, _ = composite_gl_r(radii, 2, 4, dps=25)
        P_tot = compute_P_esc(
            SLAB_POLAR_1D, r_nodes, radii, sig_t, n_angular=24,
        )
        P_out = compute_P_esc_outer(
            SLAB_POLAR_1D, r_nodes, radii, sig_t,
        )
        P_in = compute_P_esc_inner(
            SLAB_POLAR_1D, r_nodes, radii, sig_t,
        )
        expected = P_out + P_in
        rel_err = np.max(
            np.abs(P_tot - expected) / np.maximum(np.abs(expected), 1e-30)
        )
        assert rel_err < 1e-12, (
            f"Multi-region slab compute_P_esc deviates from per-face "
            f"closed-form sum: rel_err = {rel_err:.3e}. This is the "
            f"ERR-033 / Issue #131 signature replicated on the legacy "
            f"single-surface aggregate."
        )

    @pytest.mark.parametrize(
        "radii, sig_t",
        [
            (np.array([1.0, 2.0]), np.array([0.5, 2.0])),
            (np.array([0.3, 1.0, 1.5]), np.array([1.5, 0.2, 3.0])),
        ],
        ids=["2region_thinThick", "3region_mixed"],
    )
    def test_compute_G_bc_matches_per_face_sum(self, radii, sig_t):
        """``compute_G_bc`` must equal the sum of per-face primitives.

        Mirror of the P_esc gate using the ``2·E_2(τ_face)`` form.
        """
        r_nodes, _, _ = composite_gl_r(radii, 2, 4, dps=25)
        G_tot = compute_G_bc(
            SLAB_POLAR_1D, r_nodes, radii, sig_t, n_surf_quad=24,
        )
        G_out = compute_G_bc_outer(
            SLAB_POLAR_1D, r_nodes, radii, sig_t,
        )
        G_in = compute_G_bc_inner(
            SLAB_POLAR_1D, r_nodes, radii, sig_t,
        )
        expected = G_out + G_in
        rel_err = np.max(
            np.abs(G_tot - expected) / np.maximum(np.abs(expected), 1e-30)
        )
        assert rel_err < 1e-12, (
            f"Multi-region slab compute_G_bc deviates from per-face "
            f"closed-form sum: rel_err = {rel_err:.3e}. This is the "
            f"ERR-033 / Issue #131 signature replicated on the legacy "
            f"single-surface aggregate."
        )
