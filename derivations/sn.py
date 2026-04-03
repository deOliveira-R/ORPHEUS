r"""SymPy derivation of the SN (Discrete Ordinates) eigenvalue for homogeneous media.

Starting from the 1D SN transport equation, derives that a homogeneous
slab with reflective boundary conditions reproduces the infinite-medium
eigenvalue exactly, independent of angular quadrature and spatial mesh.

The derivation proceeds from the SN equations:

1. The 1D SN transport equation for discrete direction :math:`\mu_m`:

   .. math::
       \mu_m \frac{\partial \psi_m}{\partial x}
       + \Sigma_t \psi_m = \frac{Q}{2}

2. For a homogeneous medium with reflective BCs, the flux is spatially
   flat (:math:`\partial\psi_m/\partial x = 0`) and isotropic.

3. Therefore :math:`\psi_m = Q/(2\Sigma_t)` for every direction.

4. Integrating over directions (GL weights sum to 2 on [-1,1]):
   :math:`\phi = \sum_m w_m \psi_m = Q/\Sigma_t`.

5. Substituting the fission + scattering source
   :math:`Q = \Sigma_s \phi + (1/k)\,\nu\Sigma_f \phi`
   yields :math:`k = \nu\Sigma_f / \Sigma_a`.

For multi-group, the same argument applies per group, giving the
matrix eigenvalue :math:`k = \lambda_{\max}(\mathbf{A}^{-1}\mathbf{F})`.
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


def _derive_sn_homogeneous(ng_label: str, sig_t, sig_c, sig_f, nu, chi, sig_s):
    """Derive SN eigenvalue for a homogeneous slab from the transport equation."""
    ng = len(sig_t)

    # ── Symbolic derivation ──────────────────────────────────────────
    if ng == 1:
        Sig_t, Sig_s_sym, Sig_a, nu_s, Sig_f_sym = sp.symbols(
            'Sigma_t Sigma_s Sigma_a nu Sigma_f', positive=True,
        )
        psi, Q, phi, k = sp.symbols('psi Q phi k')
        mu = sp.Symbol('mu_m')

        # Step 1-3: SN equation with ∂ψ/∂x = 0
        eq_transport = sp.Eq(Sig_t * psi, Q / 2)
        psi_sol = sp.solve(eq_transport, psi)[0]  # Q/(2*Sig_t)

        # Step 4: scalar flux (weights sum to 2)
        phi_expr = Q / Sig_t

        # Step 5: source = scattering + fission
        source = Sig_s_sym * phi + (1 / k) * nu_s * Sig_f_sym * phi

        # Step 6: Sig_t * phi = source  →  solve for k
        balance = sp.Eq(Sig_t * phi, source.subs(phi, phi))
        k_sol = sp.solve(balance, k)[0]

        sig_a_val = sig_c[0] + sig_f[0]
        k_val = float(k_sol.subs({
            nu_s: nu[0], Sig_f_sym: sig_f[0], Sig_a: sig_a_val,
            Sig_t: sig_t[0], Sig_s_sym: sig_s[0, 0],
        }))

        latex = (
            r"Starting from the 1D S\ :sub:`N` transport equation for direction "
            r":math:`\mu_m`:"
            "\n\n"
            r".. math::" "\n"
            r"   \mu_m \frac{\partial\psi_m}{\partial x} + \Sigma_t \psi_m "
            r"= \frac{Q}{2}"
            "\n\n"
            r"For a homogeneous medium with reflective BCs, "
            r":math:`\partial\psi_m/\partial x = 0`, so "
            r":math:`\psi_m = Q/(2\Sigma_t)` for all directions. "
            r"Integrating with Gauss-Legendre weights "
            r"(:math:`\sum w_m = 2`):"
            "\n\n"
            r".. math::" "\n"
            r"   \phi = \sum_m w_m \psi_m = \frac{Q}{\Sigma_t}"
            "\n\n"
            r"Substituting the source "
            r":math:`Q = \Sigma_s\phi + \frac{1}{k}\nu\Sigma_f\phi` "
            r"and cancelling :math:`\phi`:"
            "\n\n"
            r".. math::" "\n"
            rf"   k = {sp.latex(k_sol)} = "
            rf"\frac{{{nu[0]} \times {sig_f[0]}}}{{{sig_a_val}}} "
            rf"= {k_val:.6f}"
        )
    else:
        # Multi-group: matrix eigenvalue from SN balance per group
        A_sym = sp.Matrix(np.diag(sig_t)) - sp.Matrix(sig_s.T)
        F_sym = sp.Matrix(np.outer(chi, nu * sig_f))
        M_sym = A_sym.inv() * F_sym
        eigs = M_sym.eigenvals()
        k_val = float(max(sp.re(e) for e in eigs.keys()))

        latex = (
            r"For multi-group S\ :sub:`N` in a homogeneous medium, the same "
            r"argument applies per group. The spatially-flat, isotropic "
            r"solution reduces each group's transport equation to the "
            r"scalar balance, yielding the matrix eigenvalue problem:"
            "\n\n"
            r".. math::" "\n"
            r"   \mathbf{A}\phi = \frac{1}{k}\mathbf{F}\phi, \quad "
            r"\mathbf{A} = \text{diag}(\Sigma_t) - \Sigma_s^T, \quad "
            r"\mathbf{F} = \chi \otimes (\nu\Sigma_f)"
            "\n\n"
            r".. math::" "\n"
            rf"   k_\infty = {k_val:.10f}"
        )

    mix = _make_mixture(sig_t, sig_c, sig_f, nu, chi, sig_s)
    return VerificationCase(
        name=f"sn_slab_{ng_label}",
        k_inf=k_val,
        method="sn",
        geometry="slab",
        n_groups=ng,
        n_regions=1,
        materials={0: mix},
        geom_params={},
        latex=latex,
        description=f"SN 1D reflective slab, {ng}G homogeneous — from transport equation",
        tolerance="< 1e-8",
    )


def all_cases() -> list[VerificationCase]:
    """Return all SN verification cases."""
    return [
        _derive_sn_homogeneous(
            "1eg_1rg",
            sig_t=np.array([1.0]), sig_c=np.array([0.2]),
            sig_f=np.array([0.3]), nu=np.array([2.5]),
            chi=np.array([1.0]), sig_s=np.array([[0.5]]),
        ),
        _derive_sn_homogeneous(
            "2eg_1rg",
            sig_t=np.array([0.50, 1.00]), sig_c=np.array([0.01, 0.02]),
            sig_f=np.array([0.01, 0.08]), nu=np.array([2.50, 2.50]),
            chi=np.array([1.00, 0.00]),
            sig_s=np.array([[0.38, 0.10], [0.00, 0.90]]),
        ),
        _derive_sn_homogeneous(
            "4eg_1rg",
            sig_t=np.array([0.01 + 0.005 + 0.28 + 0.08 + 0.02 + 0.005,
                            0.02 + 0.01 + 0.40 + 0.12 + 0.06,
                            0.03 + 0.05 + 0.55 + 0.22,
                            0.05 + 0.10 + 0.90]),
            sig_c=np.array([0.01, 0.02, 0.03, 0.05]),
            sig_f=np.array([0.005, 0.01, 0.05, 0.10]),
            nu=np.array([2.80, 2.60, 2.50, 2.45]),
            chi=np.array([0.60, 0.35, 0.05, 0.00]),
            sig_s=np.array([
                [0.28, 0.08, 0.02, 0.005],
                [0.00, 0.40, 0.12, 0.06],
                [0.00, 0.00, 0.55, 0.22],
                [0.00, 0.00, 0.00, 0.90],
            ]),
        ),
    ]
