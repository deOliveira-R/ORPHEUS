"""Homogeneous infinite-medium reactor eigenvalue solver.

Solves for the neutron spectrum and k-infinity in an infinite homogeneous
medium.  All spatial and angular dependence integrates out; the transport
equation reduces to the pure energy-balance eigenvalue problem

    A φ = (1/k) F φ,
        A = diag(Σ_t) − Σ_s0ᵀ − 2·Σ₂ᵀ,
        F = χ ⊗ νΣ_f,

with k_inf = λ_max(A⁻¹F).

The loss matrix A is assembled as the transport-operator composition
``C − K_iso`` over a MESHLESS single-cell MaterialMesh
(``MaterialMesh.from_materials``): the collision diagonal C = diag(Σ_t)
minus the model-shared isotropic energy operators ``IsotropicScattering``
(Σ_s0ᵀ) and ``IsotropicN2N`` (2·Σ₂ᵀ).  Streaming L is identically zero in
an infinite medium and is dropped — so the whole infinite-medium spectrum
runs through the SAME operator algebra the meshed SN solver uses, not a
bespoke matrix (#276).

(n,2n) convention: the (n,2n) reaction is a loss-side multiplicity-2
transfer.  It lives ONLY in A (as 2·Σ₂ᵀ), NOT in the fission production
F — the two emitted neutrons are redistributed by 2·Σ₂, they are not
produced with the fission spectrum χ.  Production is νΣ_f only.

.. seealso:: :ref:`theory-homogeneous` — Key Facts, eigenvalue equations, scattering convention.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.operators.isotropic_scattering import (
    IsotropicN2N,
    IsotropicScattering,
)


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class HomogeneousResult:
    """Result of a homogeneous infinite reactor calculation.

    The energy-grid fields (``eg_mid``, ``de``, ``du``) are ``None`` when
    the underlying :class:`~orpheus.data.macro_xs.mixture.Mixture` has no
    physical energy grid (synthetic / Sood-style XS, post-Phase-E). In
    that case ``flux_per_energy`` / ``flux_per_lethargy`` raise — the
    quantities are not defined without a grid.
    """

    k_inf: float
    flux: np.ndarray  # (NG,) — group fluxes normalised to 100 n/cm³/s production
    eg_mid: np.ndarray | None  # (NG,) — mid-group energies (eV); None if no grid
    de: np.ndarray | None  # (NG,) — energy bin widths (eV); None if no grid
    du: np.ndarray | None  # (NG,) — lethargy bin widths; None if no grid
    sig_prod: float  # one-group production XS (1/cm)
    sig_abs: float  # one-group absorption XS (1/cm)
    mixture: Mixture

    @property
    def flux_per_energy(self) -> np.ndarray:
        if self.de is None:
            raise ValueError(
                "flux_per_energy is undefined for synthetic XS without an "
                "energy grid (mixture.eg is None). Build the Mixture from "
                "an Isotope library or pass an explicit eg= to make_mixture."
            )
        return self.flux / self.de

    @property
    def flux_per_lethargy(self) -> np.ndarray:
        if self.du is None:
            raise ValueError(
                "flux_per_lethargy is undefined for synthetic XS without "
                "an energy grid (mixture.eg is None). Build the Mixture "
                "from an Isotope library or pass an explicit eg= to "
                "make_mixture."
            )
        return self.flux / self.du


# ---------------------------------------------------------------------------
# Solver — k∞ = λ_max(A⁻¹F) over the transport operator algebra
# ---------------------------------------------------------------------------

def solve_homogeneous_infinite(mix: Mixture) -> HomogeneousResult:
    r"""Solve the infinite-medium eigenvalue problem for a homogeneous mixture.

    Assembles the loss matrix :math:`\mathbf{A} = C - K_\mathrm{iso} =
    \operatorname{diag}(\Sigma_t) - \Sigma_{s0}^{T} - 2\Sigma_2^{T}` from
    the model-shared transport operators over a meshless single-cell
    :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh`, and the
    fission production dyad :math:`\mathbf{F} = \chi \otimes \nu\Sigma_f`,
    then returns the dominant eigenpair of
    :math:`\mathbf{A}^{-1}\mathbf{F}`: :math:`k_\infty = \lambda_{\max}`
    and the flux spectrum :math:`\varphi` (the corresponding right
    eigenvector), normalised so the fission production rate
    :math:`\nu\Sigma_f \cdot \varphi = 100` n/cm³/s.

    (n,2n) enters ONLY through :math:`\mathbf{A}` (as :math:`2\Sigma_2^T`),
    never the production :math:`\mathbf{F}` — see the module docstring.

    Parameters
    ----------
    mix : Mixture
        Macroscopic cross sections for the homogeneous medium.

    Returns
    -------
    HomogeneousResult
    """
    # The meshless phase space: one cell, one region, no streaming.
    mat_xs = MaterialMesh.from_materials({0: mix}).material_xs_field()

    # Loss matrix A = C − K_iso, assembled from the transport operators:
    #   C     = diag(Σ_t)        (collision)
    #   K_iso = Σ_s0ᵀ + 2·Σ₂ᵀ    (isotropic energy transfer: scatter + n2n)
    # Streaming L is identically zero in an infinite medium and is dropped.
    sig_t = mat_xs.total_cross_section[:, 0]
    iso = IsotropicScattering(mat_xs).dense_per_material()[0]
    n2n = IsotropicN2N(mat_xs).dense_per_material()[0]
    A = np.diag(sig_t) - iso - n2n

    # Production dyad F = χ ⊗ νΣ_f  (the FissionOperator rank-1 form).
    chi = mat_xs.emission_spectrum[:, 0]
    nu_sig_f = mat_xs.fission_production[:, 0]
    F = np.outer(chi, nu_sig_f)

    # k∞ = λ_max(A⁻¹F); the flux spectrum is the dominant right eigenvector.
    eigvals, eigvecs = np.linalg.eig(np.linalg.solve(A, F))
    dominant = int(np.argmax(np.real(eigvals)))
    k_inf = float(np.real(eigvals[dominant]))
    phi = np.real(eigvecs[:, dominant])
    if phi.sum() < 0:  # sign-normalise to a physical (non-negative) spectrum
        phi = -phi

    # Normalise the flux so the fission production rate νΣ_f·φ = 100 n/cm³/s.
    phi = phi * (100.0 / float(nu_sig_f @ phi))

    prod_rate = float(nu_sig_f @ phi)
    abs_rate = float(mix.absorption_xs @ phi)
    total_flux = float(phi.sum())

    ng = mix.ng
    if mix.eg is None:
        # Synthetic XS — no physical energy grid, so lethargy / per-energy
        # diagnostics are not defined.  k_inf and the flux spectrum still
        # carry meaningful information; only the per-energy plotting path
        # is unavailable.
        eg_mid = None
        de = None
        du = None
    else:
        eg = mix.eg
        eg_mid = 0.5 * (eg[:ng] + eg[1:ng + 1])
        de = np.abs(eg[1:ng + 1] - eg[:ng])
        du = np.abs(np.log(eg[1:ng + 1] / eg[:ng]))

    return HomogeneousResult(
        k_inf=k_inf,
        flux=phi,
        eg_mid=eg_mid,
        de=de,
        du=du,
        sig_prod=prod_rate / total_flux,
        sig_abs=abs_rate / total_flux,
        mixture=mix,
    )
