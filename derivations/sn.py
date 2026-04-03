r"""SN (Discrete Ordinates) eigenvalue derivation for homogeneous media.

Starting from the 1D SN transport equation, derives that a homogeneous
slab with reflective boundary conditions reproduces the infinite-medium
eigenvalue exactly, independent of angular quadrature and spatial mesh.

For multi-group, the matrix eigenvalue :math:`k = \lambda_{\max}(\mathbf{A}^{-1}\mathbf{F})`
is derived from the SN equations with spatially-flat, isotropic flux.
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from ._types import VerificationCase
from ._xs_library import get_xs, get_mixture


def _derive_sn_homogeneous(ng_key: str) -> VerificationCase:
    """Derive SN eigenvalue for a homogeneous slab from the transport equation."""
    xs = get_xs("A", ng_key)
    ng = len(xs["sig_t"])

    if ng == 1:
        nu_s, Sig_f_sym, Sig_a = sp.symbols(
            'nu Sigma_f Sigma_a', positive=True,
        )
        k_expr = nu_s * Sig_f_sym / Sig_a

        sig_a_val = xs["sig_c"][0] + xs["sig_f"][0]
        k_val = float(k_expr.subs({
            nu_s: xs["nu"][0], Sig_f_sym: xs["sig_f"][0], Sig_a: sig_a_val,
        }))

        latex = (
            r"From the 1D S\ :sub:`N` equation with "
            r":math:`\partial\psi_m/\partial x = 0` (homogeneous, reflective BCs):"
            "\n\n"
            r".. math::" "\n"
            rf"   k = {sp.latex(k_expr)} = {k_val:.6f}"
        )
    else:
        sig_t = xs["sig_t"]
        sig_s = xs["sig_s"]
        chi = xs["chi"]
        nu_sig_f = xs["nu"] * xs["sig_f"]

        A_sym = sp.Matrix(np.diag(sig_t)) - sp.Matrix(sig_s.T)
        F_sym = sp.Matrix(np.outer(chi, nu_sig_f))
        M_sym = A_sym.inv() * F_sym
        eigs = M_sym.eigenvals()
        k_val = float(max(sp.re(e) for e in eigs.keys()))

        latex = (
            r"Multi-group S\ :sub:`N` homogeneous: flat flux reduces to "
            r":math:`k = \lambda_{\max}(\mathbf{A}^{-1}\mathbf{F})`."
            "\n\n"
            r".. math::" "\n"
            rf"   k_\infty = {k_val:.10f}"
        )

    mix = get_mixture("A", ng_key)
    return VerificationCase(
        name=f"sn_slab_{ng}eg_1rg",
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
    return [_derive_sn_homogeneous(ng_key) for ng_key in ["1g", "2g", "4g"]]
