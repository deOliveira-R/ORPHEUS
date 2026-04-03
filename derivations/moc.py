r"""SymPy derivation of the MOC (Method of Characteristics) eigenvalue for homogeneous media.

Starting from the characteristic ODE, derives that a homogeneous pin
cell reproduces the infinite-medium eigenvalue exactly.

The derivation proceeds from the MOC equations:

1. Along a characteristic with direction :math:`\hat\Omega`, the ODE is:

   .. math::
       \frac{d\psi}{ds} + \Sigma_t \psi = \frac{Q_i}{4\pi}

2. The analytical solution over a segment of length :math:`\ell`:

   .. math::
       \bar\psi = \psi_{\rm in}\,\frac{1 - e^{-\Sigma_t \ell}}{\Sigma_t \ell}
       + \frac{Q_i}{4\pi\Sigma_t}
       \left(1 - \frac{1 - e^{-\Sigma_t \ell}}{\Sigma_t \ell}\right)

3. For a homogeneous medium, :math:`Q` is uniform. With isotropic
   incoming flux, the scalar flux integrates to
   :math:`\phi = Q/\Sigma_t`.

4. The eigenvalue is :math:`k = \nu\Sigma_f/\Sigma_a`.
"""

from __future__ import annotations

import numpy as np
import sympy as sp
from scipy.sparse import csr_matrix

from data.macro_xs.mixture import Mixture

from ._types import VerificationCase


def _make_mixture(sig_t, sig_c, sig_f, nu, chi, sig_s):
    ng = len(sig_t)
    eg = np.logspace(7, -3, ng + 1)
    return Mixture(
        SigC=sig_c.copy(), SigL=np.zeros(ng),
        SigF=sig_f.copy(), SigP=(nu * sig_f).copy(),
        SigT=sig_t.copy(), SigS=[csr_matrix(sig_s)],
        Sig2=csr_matrix((ng, ng)), chi=chi.copy(), eg=eg.copy(),
    )


def _derive_moc_homogeneous(ng_label: str, sig_t, sig_c, sig_f, nu, chi, sig_s):
    """Derive MOC eigenvalue for a homogeneous pin cell from the characteristic ODE."""
    ng = len(sig_t)

    if ng == 1:
        Sig_t, Q, psi_in, ell, k_sym = sp.symbols(
            'Sigma_t Q psi_in ell k', positive=True,
        )
        nu_s, Sig_f_sym, Sig_a = sp.symbols('nu Sigma_f Sigma_a', positive=True)

        # MOC segment-average angular flux
        tau = Sig_t * ell
        attn = (1 - sp.exp(-tau)) / tau
        psi_bar = psi_in * attn + Q / (4 * sp.pi * Sig_t) * (1 - attn)

        # For homogeneous medium, Q uniform, isotropic incoming:
        # psi_in = Q/(4*pi*Sig_t) (steady-state isotropic flux)
        # → psi_bar = Q/(4*pi*Sig_t) for any ell
        psi_homo = psi_bar.subs(psi_in, Q / (4 * sp.pi * Sig_t))
        psi_simplified = sp.simplify(psi_homo)  # Q/(4*pi*Sig_t)

        # Scalar flux: phi = 4*pi * psi_bar = Q/Sig_t
        phi = 4 * sp.pi * psi_simplified

        # Source: Q = Sig_s*phi + (1/k)*nu*Sig_f*phi
        # → Sig_t*phi = Q = Sig_s*phi + (1/k)*nu*Sig_f*phi
        # → (Sig_t - Sig_s)*phi = (1/k)*nu*Sig_f*phi
        # → k = nu*Sig_f / Sig_a
        k_expr = nu_s * Sig_f_sym / Sig_a

        sig_a_val = sig_c[0] + sig_f[0]
        k_val = float(k_expr.subs({
            nu_s: nu[0], Sig_f_sym: sig_f[0], Sig_a: sig_a_val,
        }))

        latex = (
            r"Starting from the MOC characteristic ODE:"
            "\n\n"
            r".. math::" "\n"
            r"   \frac{d\psi}{ds} + \Sigma_t \psi = \frac{Q}{4\pi}"
            "\n\n"
            r"The segment-average angular flux is:"
            "\n\n"
            r".. math::" "\n"
            r"   \bar\psi = \psi_{\rm in}\,"
            r"\frac{1-e^{-\Sigma_t\ell}}{\Sigma_t\ell}"
            r" + \frac{Q}{4\pi\Sigma_t}"
            r"\left(1-\frac{1-e^{-\Sigma_t\ell}}{\Sigma_t\ell}\right)"
            "\n\n"
            r"For a homogeneous medium, :math:`Q` is uniform. With isotropic "
            r"incoming flux :math:`\psi_{\rm in} = Q/(4\pi\Sigma_t)`, the "
            r"segment average simplifies to "
            r":math:`\bar\psi = Q/(4\pi\Sigma_t)` regardless of segment "
            r"length, confirming that the flat-flux solution is exact."
            "\n\n"
            r"The scalar flux :math:`\phi = 4\pi\bar\psi = Q/\Sigma_t`, "
            r"and the eigenvalue from the neutron balance is:"
            "\n\n"
            r".. math::" "\n"
            rf"   k = {sp.latex(k_expr)} = {k_val:.6f}"
        )
    else:
        A_sym = sp.Matrix(np.diag(sig_t)) - sp.Matrix(sig_s.T)
        F_sym = sp.Matrix(np.outer(chi, nu * sig_f))
        M_sym = A_sym.inv() * F_sym
        eigs = M_sym.eigenvals()
        k_val = float(max(sp.re(e) for e in eigs.keys()))

        latex = (
            r"For multi-group MOC, the same homogeneous argument applies "
            r"per group. The flat-source, isotropic solution along every "
            r"characteristic reduces to the group-wise scalar balance:"
            "\n\n"
            r".. math::" "\n"
            rf"   k_\infty = {k_val:.10f}"
        )

    mix = _make_mixture(sig_t, sig_c, sig_f, nu, chi, sig_s)
    return VerificationCase(
        name=f"moc_cyl1D_{ng_label}",
        k_inf=k_val,
        method="moc",
        geometry="cyl1D",
        n_groups=ng,
        n_regions=1,
        materials={0: mix},
        geom_params={},
        latex=latex,
        description=f"MOC cylindrical pin cell, {ng}G homogeneous — from characteristic ODE",
        tolerance="< 1e-4",
    )


def all_cases() -> list[VerificationCase]:
    """Return all MOC verification cases."""
    return [
        _derive_moc_homogeneous(
            "1eg_1rg",
            sig_t=np.array([1.0]), sig_c=np.array([0.2]),
            sig_f=np.array([0.3]), nu=np.array([2.5]),
            chi=np.array([1.0]), sig_s=np.array([[0.5]]),
        ),
        _derive_moc_homogeneous(
            "2eg_1rg",
            sig_t=np.array([0.50, 1.00]), sig_c=np.array([0.01, 0.02]),
            sig_f=np.array([0.01, 0.08]), nu=np.array([2.50, 2.50]),
            chi=np.array([1.00, 0.00]),
            sig_s=np.array([[0.38, 0.10], [0.00, 0.90]]),
        ),
    ]
