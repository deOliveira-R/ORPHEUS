r"""MOC (Method of Characteristics) eigenvalue derivation for homogeneous media.

Starting from the characteristic ODE, derives that a homogeneous pin
cell reproduces the infinite-medium eigenvalue exactly.
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from ._types import VerificationCase
from ._xs_library import get_xs, get_mixture


def _derive_moc_homogeneous(ng_key: str) -> VerificationCase:
    """Derive MOC eigenvalue for a homogeneous pin cell."""
    xs = get_xs("A", ng_key)
    ng = len(xs["sig_t"])

    if ng == 1:
        nu_s, Sig_f_sym, Sig_a = sp.symbols('nu Sigma_f Sigma_a', positive=True)
        k_expr = nu_s * Sig_f_sym / Sig_a

        sig_a_val = xs["sig_c"][0] + xs["sig_f"][0]
        k_val = float(k_expr.subs({
            nu_s: xs["nu"][0], Sig_f_sym: xs["sig_f"][0], Sig_a: sig_a_val,
        }))

        latex = (
            r"From the MOC characteristic ODE with homogeneous, isotropic "
            r"source: :math:`\bar\psi = Q/(4\pi\Sigma_t)`, giving"
            "\n\n"
            r".. math::" "\n"
            rf"   k = {sp.latex(k_expr)} = {k_val:.6f}"
        )
    else:
        A_sym = sp.Matrix(np.diag(xs["sig_t"])) - sp.Matrix(xs["sig_s"].T)
        F_sym = sp.Matrix(np.outer(xs["chi"], xs["nu"] * xs["sig_f"]))
        M_sym = A_sym.inv() * F_sym
        eigs = M_sym.eigenvals()
        k_val = float(max(sp.re(e) for e in eigs.keys()))

        latex = (
            r"Multi-group MOC homogeneous: flat-source along every "
            r"characteristic reduces to the matrix eigenvalue."
            "\n\n"
            r".. math::" "\n"
            rf"   k_\infty = {k_val:.10f}"
        )

    mix = get_mixture("A", ng_key)
    return VerificationCase(
        name=f"moc_cyl1D_{ng}eg_1rg",
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
    return [_derive_moc_homogeneous(ng_key) for ng_key in ["1g", "2g", "4g"]]
