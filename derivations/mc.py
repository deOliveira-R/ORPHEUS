r"""SymPy derivation of the Monte Carlo eigenvalue for homogeneous media.

Starting from random walk probability theory, derives the expected
multiplication factor for a homogeneous infinite medium.

The derivation proceeds from collision probabilities:

1. At each collision in a homogeneous medium, the probability of each
   event is determined by the cross-section ratios:

   .. math::
       P(\text{scatter}) = \frac{\Sigma_s}{\Sigma_t}, \quad
       P(\text{fission}) = \frac{\Sigma_f}{\Sigma_t}, \quad
       P(\text{capture}) = \frac{\Sigma_c}{\Sigma_t}

2. At a fission event, :math:`\nu` secondary neutrons are produced on
   average.

3. The expected number of secondaries per collision is
   :math:`\nu\Sigma_f/\Sigma_t`.

4. A neutron undergoes on average :math:`\Sigma_t/\Sigma_a` collisions
   before absorption (geometric series over scattering probability).

5. Therefore the expected multiplication factor is:

   .. math::
       k = \frac{\nu\Sigma_f}{\Sigma_t} \times \frac{\Sigma_t}{\Sigma_a}
         = \frac{\nu\Sigma_f}{\Sigma_a}

6. The Monte Carlo estimator converges to this value with standard
   error :math:`\sigma \sim 1/\sqrt{N_{\text{active}}}` by the
   Central Limit Theorem.

For multi-group, the collision kernel becomes a matrix and
:math:`k = \lambda_{\max}(\mathbf{A}^{-1}\mathbf{F})`.
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


def _derive_mc_homogeneous(ng_label: str, sig_t, sig_c, sig_f, nu, chi, sig_s):
    """Derive MC eigenvalue from random walk probability theory."""
    ng = len(sig_t)

    if ng == 1:
        nu_s, Sig_f_sym, Sig_t_sym, Sig_a = sp.symbols(
            'nu Sigma_f Sigma_t Sigma_a', positive=True,
        )
        Sig_s_sym, Sig_c_sym = sp.symbols('Sigma_s Sigma_c', positive=True)

        # Probability of fission per collision
        P_fission = Sig_f_sym / Sig_t_sym
        # Expected secondaries per collision
        secondaries_per_collision = nu_s * P_fission
        # Mean number of collisions before absorption
        n_collisions = Sig_t_sym / Sig_a
        # Expected multiplication factor
        k_expr = sp.simplify(secondaries_per_collision * n_collisions)

        sig_a_val = sig_c[0] + sig_f[0]
        k_val = float(k_expr.subs({
            nu_s: nu[0], Sig_f_sym: sig_f[0],
            Sig_t_sym: sig_t[0], Sig_a: sig_a_val,
        }))

        latex = (
            r"From random walk probability theory in a homogeneous medium:"
            "\n\n"
            r".. math::" "\n"
            r"   P(\text{fission per collision}) = "
            r"\frac{\Sigma_f}{\Sigma_t}"
            "\n\n"
            r".. math::" "\n"
            r"   \langle\text{collisions before absorption}\rangle = "
            r"\frac{\Sigma_t}{\Sigma_a} = "
            r"\frac{1}{1 - \Sigma_s/\Sigma_t}"
            "\n\n"
            r"The expected multiplication factor over one generation:"
            "\n\n"
            r".. math::" "\n"
            rf"   k = \nu \cdot \frac{{\Sigma_f}}{{\Sigma_t}} "
            rf"\cdot \frac{{\Sigma_t}}{{\Sigma_a}} = "
            rf"{sp.latex(k_expr)} = {k_val:.6f}"
            "\n\n"
            r"The MC estimator :math:`\hat{k}` converges to this value "
            r"with standard error :math:`\sigma \sim 1/\sqrt{N_{\rm active}}`."
        )
    else:
        A_sym = sp.Matrix(np.diag(sig_t)) - sp.Matrix(sig_s.T)
        F_sym = sp.Matrix(np.outer(chi, nu * sig_f))
        M_sym = A_sym.inv() * F_sym
        eigs = M_sym.eigenvals()
        k_val = float(max(sp.re(e) for e in eigs.keys()))

        latex = (
            r"For multi-group Monte Carlo, the collision kernel is a matrix. "
            r"The group transfer probabilities at each collision are determined "
            r"by the scattering and fission matrices. The expected multiplication "
            r"factor is the dominant eigenvalue:"
            "\n\n"
            r".. math::" "\n"
            r"   k = \lambda_{\max}(\mathbf{A}^{-1}\mathbf{F})"
            "\n\n"
            r".. math::" "\n"
            rf"   k_\infty = {k_val:.10f}"
        )

    mix = _make_mixture(sig_t, sig_c, sig_f, nu, chi, sig_s)
    return VerificationCase(
        name=f"mc_cyl1D_{ng_label}",
        k_inf=k_val,
        method="mc",
        geometry="cyl1D",
        n_groups=ng,
        n_regions=1,
        materials={0: mix},
        geom_params={},
        latex=latex,
        description=f"MC pin cell, {ng}G homogeneous — from random walk probability",
        tolerance="z < 5\u03c3",
    )


def all_cases() -> list[VerificationCase]:
    """Return all MC verification cases."""
    return [
        _derive_mc_homogeneous(
            "1eg_1rg",
            sig_t=np.array([1.0]), sig_c=np.array([0.2]),
            sig_f=np.array([0.3]), nu=np.array([2.5]),
            chi=np.array([1.0]), sig_s=np.array([[0.5]]),
        ),
        _derive_mc_homogeneous(
            "2eg_1rg",
            sig_t=np.array([0.50, 1.00]), sig_c=np.array([0.01, 0.02]),
            sig_f=np.array([0.01, 0.08]), nu=np.array([2.50, 2.50]),
            chi=np.array([1.00, 0.00]),
            sig_s=np.array([[0.38, 0.10], [0.00, 0.90]]),
        ),
    ]
