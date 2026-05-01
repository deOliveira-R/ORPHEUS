r"""Paired symbolic-vs-production contract test for the slab per-face
specular primitives.

Math-origin pattern: the SymPy derivation in
:mod:`orpheus.derivations.peierls_specular.slab` is the **source of
truth** for the per-face slab primitives and the block-diagonal
reflection operator
:math:`R_{\rm slab} = {\rm blockdiag}(R_{\rm spec}, R_{\rm spec})`.

This test consumes the SymPy origin as a contract and pins:

1. **Mode-0 P primitive** equals :math:`(1/2)\,E_2(\tau)` numerically
   for slab outer/inner faces — i.e. matches
   :func:`compute_P_esc_outer` / :func:`compute_P_esc_inner` (slab
   branch) bit-exactly.
2. **Mode-0 G primitive** equals :math:`2\,E_2(\tau)` numerically for
   slab outer/inner faces — matches :func:`compute_G_bc_outer` /
   :func:`compute_G_bc_inner` (slab branch) bit-exactly.
3. **Block-diagonal partial-current contract**:
   :math:`2\,M_{\rm slab}\,R_{\rm slab} = I_{2N}` per face, i.e. each
   :math:`N \times N` block satisfies the curvilinear specular
   identity from :mod:`.r_matrix`.

The verified equation is :eq:`peierls-slab-Pesc-mode` (slab per-face
mode-:math:`N` primitive, no-:math:`\mu`-weight basis).

The canonical derivation source is
``orpheus/derivations/peierls_specular/slab.py:1`` (SymPy module
docstring with the basis-choice resolution). The shipped slab
production primitives are at
``orpheus/derivations/peierls_geometry.py:2906`` (P_esc_outer) and
``orpheus/derivations/peierls_geometry.py:3218`` (G_bc_outer).
"""
from __future__ import annotations

import numpy as np
import pytest
import sympy as sp

from orpheus.derivations.continuous.peierls.geometry import (
    SLAB_POLAR_1D,
    compute_G_bc_inner,
    compute_G_bc_outer,
    compute_P_esc_inner,
    compute_P_esc_outer,
)
from orpheus.derivations.continuous.peierls.origins.specular import (
    build_M_symbolic,
    build_R_slab_blockdiag,
    derive_slab_g_outer_n0,
    derive_slab_p_outer_n0,
)


# ═══════════════════════════════════════════════════════════════════════
# Mode-0 reductions: ½ E_2 and 2 E_2 match production slab branch
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-slab-Pesc-mode")
def test_slab_p_primitive_n0_matches_E2():
    r"""For :math:`n = 0`, the slab outer P primitive
    :func:`derive_slab_p_outer_n0` equals :math:`(1/2)\,E_2(\tau)`
    symbolically; production :func:`compute_P_esc_outer` (slab
    branch, ``orpheus/derivations/peierls_geometry.py:2933``) uses
    the same closed form. Numerical agreement at machine precision is
    required.
    """
    p_n0 = derive_slab_p_outer_n0()
    tau = sp.Symbol("tau", positive=True)
    from scipy.special import expn
    p_lambda = sp.lambdify(
        tau, p_n0, modules=[{"expint": lambda n, x: float(expn(int(n), float(x)))},
                            "numpy"],
    )

    radii = np.array([2.0])
    sig_t = np.array([1.5])
    x_nodes = np.array([0.1, 0.7, 1.5, 1.9])

    p_prod = compute_P_esc_outer(
        SLAB_POLAR_1D, x_nodes, radii, sig_t, n_angular=32, dps=25,
    )
    # τ_outer(x_i) = sig_t · (L - x_i) for homogeneous slab
    tau_vals = float(sig_t[0]) * (float(radii[0]) - x_nodes)
    p_ref = np.array([float(p_lambda(t)) for t in tau_vals])

    np.testing.assert_allclose(
        p_prod, p_ref, rtol=1e-12, atol=1e-13,
        err_msg=(
            "compute_P_esc_outer (slab) drifted from symbolic "
            "(1/2) E_2(τ) closed form."
        ),
    )


@pytest.mark.l1
@pytest.mark.verifies("peierls-slab-Gbc-mode")
def test_slab_g_primitive_n0_matches_2E2():
    r"""For :math:`n = 0`, the slab outer G primitive
    :func:`derive_slab_g_outer_n0` equals :math:`2\,E_2(\tau)`
    symbolically; production :func:`compute_G_bc_outer` (slab branch,
    ``orpheus/derivations/peierls_geometry.py:3242``) uses the same
    closed form. Numerical agreement at machine precision is required.
    """
    g_n0 = derive_slab_g_outer_n0()
    tau = sp.Symbol("tau", positive=True)
    from scipy.special import expn
    g_lambda = sp.lambdify(
        tau, g_n0, modules=[{"expint": lambda n, x: float(expn(int(n), float(x)))},
                            "numpy"],
    )

    radii = np.array([2.0])
    sig_t = np.array([1.5])
    x_nodes = np.array([0.1, 0.7, 1.5, 1.9])

    g_prod = compute_G_bc_outer(
        SLAB_POLAR_1D, x_nodes, radii, sig_t, n_surf_quad=32, dps=25,
    )
    tau_vals = float(sig_t[0]) * (float(radii[0]) - x_nodes)
    g_ref = np.array([float(g_lambda(t)) for t in tau_vals])

    np.testing.assert_allclose(
        g_prod, g_ref, rtol=1e-12, atol=1e-13,
        err_msg=(
            "compute_G_bc_outer (slab) drifted from symbolic "
            "2 E_2(τ) closed form."
        ),
    )


@pytest.mark.l1
@pytest.mark.verifies("peierls-slab-Pesc-mode")
def test_slab_p_primitive_n0_inner_matches_E2():
    r"""Inner-face counterpart of :func:`test_slab_p_primitive_n0_matches_E2`.

    Production :func:`compute_P_esc_inner` (slab branch,
    ``orpheus/derivations/peierls_geometry.py:3020``) uses
    :math:`(1/2)\,E_2(\Sigma_t\,x_i)` for the perpendicular optical
    depth from :math:`x_i` to the inner face at :math:`x = 0`.
    """
    p_n0 = derive_slab_p_outer_n0()
    tau = sp.Symbol("tau", positive=True)
    from scipy.special import expn
    p_lambda = sp.lambdify(
        tau, p_n0, modules=[{"expint": lambda n, x: float(expn(int(n), float(x)))},
                            "numpy"],
    )

    radii = np.array([2.0])
    sig_t = np.array([1.5])
    x_nodes = np.array([0.1, 0.7, 1.5, 1.9])

    p_prod = compute_P_esc_inner(
        SLAB_POLAR_1D, x_nodes, radii, sig_t, n_angular=32, dps=25,
    )
    tau_vals = float(sig_t[0]) * x_nodes  # τ from x_i to inner face
    p_ref = np.array([float(p_lambda(t)) for t in tau_vals])

    np.testing.assert_allclose(
        p_prod, p_ref, rtol=1e-12, atol=1e-13,
        err_msg=(
            "compute_P_esc_inner (slab) drifted from symbolic "
            "(1/2) E_2(τ) closed form."
        ),
    )


@pytest.mark.l1
@pytest.mark.verifies("peierls-slab-Gbc-mode")
def test_slab_g_primitive_n0_inner_matches_2E2():
    r"""Inner-face counterpart of :func:`test_slab_g_primitive_n0_matches_2E2`.

    Production :func:`compute_G_bc_inner` (slab branch) uses
    :math:`2\,E_2(\Sigma_t\,x_i)`.
    """
    g_n0 = derive_slab_g_outer_n0()
    tau = sp.Symbol("tau", positive=True)
    from scipy.special import expn
    g_lambda = sp.lambdify(
        tau, g_n0, modules=[{"expint": lambda n, x: float(expn(int(n), float(x)))},
                            "numpy"],
    )

    radii = np.array([2.0])
    sig_t = np.array([1.5])
    x_nodes = np.array([0.1, 0.7, 1.5, 1.9])

    g_prod = compute_G_bc_inner(
        SLAB_POLAR_1D, x_nodes, radii, sig_t, n_surf_quad=32, dps=25,
    )
    tau_vals = float(sig_t[0]) * x_nodes
    g_ref = np.array([float(g_lambda(t)) for t in tau_vals])

    np.testing.assert_allclose(
        g_prod, g_ref, rtol=1e-12, atol=1e-13,
        err_msg=(
            "compute_G_bc_inner (slab) drifted from symbolic "
            "2 E_2(τ) closed form."
        ),
    )


# ═══════════════════════════════════════════════════════════════════════
# Block-diagonal 2 M R = I per face — non-trivial new content from lift
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-specular-bc-defn")
@pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6])
def test_slab_blockdiag_2MR_eq_I_per_face(n):
    r"""The slab per-face partial-current contract:

    .. math::

       2\,{\rm blockdiag}(M, M)\,R_{\rm slab} \;=\; I_{2N}.

    This is the **non-trivial algebraic identity** the slab
    decomposition lifts above the curvilinear case — the per-face
    block structure is locked in by the requirement that the
    rank-:math:`N` partial-current identity
    :math:`J^{-}_m = J^{+}_m` holds **separately at each face**, not
    coupled across faces. Each :math:`N \times N` block of
    :math:`R_{\rm slab}` is identical to the curvilinear
    :math:`R_{\rm spec}(N) = (1/2)\,M^{-1}`.
    """
    M = build_M_symbolic(n)
    M_block = sp.zeros(2 * n, 2 * n)
    M_block[:n, :n] = M
    M_block[n:, n:] = M

    R_slab = build_R_slab_blockdiag(n)

    diff = sp.simplify(2 * M_block * R_slab - sp.eye(2 * n))
    assert diff == sp.zeros(2 * n, 2 * n), (
        f"2 M_block R_slab != I at N={n}: residual = {diff}"
    )
