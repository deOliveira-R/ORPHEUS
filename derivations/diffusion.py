"""SymPy derivation for 2-group diffusion buckling eigenvalue.

Derives the analytical k_eff for a bare homogeneous slab with
vacuum boundary conditions using the buckling B² = (π/H)².
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from ._types import VerificationCase


def derive_2g_diffusion(fuel_height: float = 50.0) -> VerificationCase:
    r"""2-group diffusion eigenvalue for a bare slab.

    The removal matrix with buckling :math:`B^2 = (\pi/H)^2`:

    .. math::
        \mathbf{A} = \text{diag}(D B^2 + \Sigma_a + \Sigma_{s,\text{out}})
        - \Sigma_s^T

    and the fission production matrix:

    .. math::
        \mathbf{F} = \chi \otimes (\nu\Sigma_f)

    The eigenvalue is :math:`k = \lambda_{\max}(\mathbf{A}^{-1}\mathbf{F})`.
    """
    # Cross sections from the default diffusion problem
    transport = np.array([0.2181, 0.7850])
    absorption = np.array([0.0096, 0.0959])
    fission = np.array([0.0024, 0.0489])
    production = np.array([0.0061, 0.1211])
    chi = np.array([1.0, 0.0])
    scattering = np.array([0.0160, 0.0])  # downscatter only

    # Symbolic derivation
    D1, D2 = sp.symbols('D_1 D_2', positive=True)
    B2_sym = sp.Symbol('B2', positive=True)
    Sa1, Sa2 = sp.symbols('Sigma_a1 Sigma_a2', positive=True)
    Ss12 = sp.Symbol('Sigma_s12', positive=True)
    nSf1, nSf2 = sp.symbols('nuSigma_f1 nuSigma_f2', positive=True)
    chi1 = sp.Symbol('chi_1')

    D_coeff = 1.0 / (3.0 * transport)
    B2_val = (np.pi / fuel_height) ** 2

    # Symbolic A matrix (removal + leakage - inscatter)
    A_sym = sp.Matrix([
        [D1 * B2_sym + Sa1 + Ss12, 0],
        [-Ss12, D2 * B2_sym + Sa2],
    ])
    F_sym = sp.Matrix([
        [chi1 * nSf1, chi1 * nSf2],
        [0, 0],
    ])

    # Substitute numeric values
    subs = {
        D1: D_coeff[0], D2: D_coeff[1],
        B2_sym: B2_val,
        Sa1: absorption[0], Sa2: absorption[1],
        Ss12: scattering[0],
        nSf1: production[0], nSf2: production[1],
        chi1: chi[0],
    }

    A_num = A_sym.subs(subs)
    F_num = F_sym.subs(subs)
    M_num = A_num.inv() * F_num
    eigs = M_num.eigenvals()
    k_val = float(max(sp.re(e) for e in eigs.keys()))

    # Also compute via numpy for cross-check
    A_np = np.diag(D_coeff * B2_val + absorption + scattering) \
        - np.array([[0.0, 0.0], [scattering[0], 0.0]])
    F_np = np.outer(chi, production)
    M_np = np.linalg.solve(A_np, F_np)
    k_numpy = float(np.max(np.real(np.linalg.eigvals(M_np))))
    assert abs(k_val - k_numpy) < 1e-10, \
        f"SymPy/numpy mismatch: {k_val} vs {k_numpy}"

    latex = (
        rf"For a bare homogeneous slab of height :math:`H = {fuel_height}` cm "
        r"with vacuum boundary conditions:"
        "\n\n"
        r".. math::" "\n"
        rf"   B^2 = \left(\frac{{\pi}}{{{fuel_height}}}\right)^2 "
        rf"= {B2_val:.6e}"
        "\n\n"
        r".. math::" "\n"
        rf"   \mathbf{{A}} = {sp.latex(A_sym)}"
        "\n\n"
        r".. math::" "\n"
        rf"   \mathbf{{F}} = {sp.latex(F_sym)}"
        "\n\n"
        r".. math::" "\n"
        rf"   k_{{\text{{eff}}}} = \lambda_{{\max}}"
        rf"(\mathbf{{A}}^{{-1}}\mathbf{{F}}) = {k_val:.10f}"
    )

    # Store XS as a dict (diffusion solver uses TwoGroupXS, not Mixture)
    xs_dict = dict(
        transport=transport, absorption=absorption,
        fission=fission, production=production,
        chi=chi, scattering=scattering,
    )

    return VerificationCase(
        name="dif_slab_2eg_1rg",
        k_inf=k_val,
        method="dif",
        geometry="slab",
        n_groups=2,
        n_regions=1,
        materials=xs_dict,  # type: ignore[arg-type]
        geom_params=dict(fuel_height=fuel_height),
        latex=latex,
        description=(
            f"2-group diffusion bare slab (H={fuel_height} cm, vacuum BCs)"
        ),
        tolerance="O(h²)",
    )


def all_cases() -> list[VerificationCase]:
    """Return all diffusion verification cases."""
    return [derive_2g_diffusion()]
