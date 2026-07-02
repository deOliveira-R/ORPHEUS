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
from typing import TYPE_CHECKING

import numpy as np

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.numerics.eigenvalue import direct_eigenvalue
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.operators.fission import FissionOperator
from orpheus.transport.operators.isotropic_scattering import (
    IsotropicN2N,
    IsotropicScattering,
)
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.transport.reaction_rate_functional import IntegratedReactionRate

if TYPE_CHECKING:
    from orpheus.transport.mesh.material_xs_field import MaterialXSField


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class HomogeneousResult:
    """Result of a homogeneous infinite reactor calculation.

    The energy-grid diagnostics (``representative_energy``, ``energy_widths``,
    ``lethargy_widths``) are ``None`` when the underlying
    :class:`~orpheus.data.macro_xs.mixture.Mixture` has no physical energy grid
    (synthetic / Sood-style XS, post-Phase-E). In that case ``flux_per_energy`` /
    ``flux_per_lethargy`` raise — the quantities are not defined without a grid.

    These three are the :class:`~orpheus.data.energy_grid.EnergyGrid` value
    object's own properties — ``representative_energy`` is the **geometric** group
    centre :math:`\\sqrt{E_{\\rm up}E_{\\rm lo}}` (the natural abscissa on the log
    energy axis, NOT the arithmetic midpoint) — read from ``mixture.energy_grid``
    rather than re-deriving the group geometry here.
    """

    k_inf: float
    flux: np.ndarray  # (NG,) — group fluxes normalised to 100 n/cm³/s production
    representative_energy: np.ndarray | None  # (NG,) — geometric group-centre energies (eV); None if no grid
    energy_widths: np.ndarray | None  # (NG,) — ΔE group widths (eV); None if no grid
    lethargy_widths: np.ndarray | None  # (NG,) — Δu lethargy widths; None if no grid
    sig_prod: float  # one-group production XS (1/cm)
    sig_abs: float  # one-group absorption XS (1/cm)
    mixture: Mixture

    @property
    def flux_per_energy(self) -> np.ndarray:
        if self.energy_widths is None:
            raise ValueError(
                "flux_per_energy is undefined for synthetic XS without an "
                "energy grid (mixture.eg is None). Build the Mixture from "
                "an Isotope library or pass an explicit eg= to make_mixture."
            )
        return self.flux / self.energy_widths

    @property
    def flux_per_lethargy(self) -> np.ndarray:
        if self.lethargy_widths is None:
            raise ValueError(
                "flux_per_lethargy is undefined for synthetic XS without "
                "an energy grid (mixture.eg is None). Build the Mixture "
                "from an Isotope library or pass an explicit eg= to "
                "make_mixture."
            )
        return self.flux / self.lethargy_widths


# ---------------------------------------------------------------------------
# Solver — k∞ = λ_max(A⁻¹F) over the transport operator algebra
# ---------------------------------------------------------------------------

def _assemble_loss_matrix(mat_xs: "MaterialXSField") -> np.ndarray:
    r"""The dense loss matrix :math:`A = C - K_\mathrm{iso}` for an infinite medium.

    Composes the model-shared transport operators on the meshless single-cell
    carrier — collision :math:`C = M[\Sigma_t]`
    (:class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`)
    minus the isotropic energy transfer :math:`K_\mathrm{iso} = \Sigma_{s0}^T
    + 2\Sigma_2^T`
    (:class:`~orpheus.transport.operators.isotropic_scattering.IsotropicScattering`
    + :class:`~orpheus.transport.operators.isotropic_scattering.IsotropicN2N`).
    Streaming :math:`L` is identically zero in an infinite medium and dropped.
    The operator algebra ``C - K_iso`` is an
    :class:`~orpheus.numerics.operator.OperatorSum`, realized as a dense
    ``(ng, ng)`` matrix via the operator's own
    :meth:`~orpheus.numerics.operator.LinearOperator.as_matrix` — the
    apply-to-basis materialization this module's retired ``_as_dense``
    prototyped (promoted at taxonomy §12 step 5).  ``basis_shape=(ng, 1)``
    is the operators' group-leading ``(ng, *spatial)`` bare-ndarray
    contract on the meshless single cell — passed explicitly because the
    meshless operators carry no :attr:`~orpheus.numerics.operator.LinearOperator.domain`
    space to derive it from.
    """
    collision = MultiplicationOperator.from_mesh(
        mat_xs.total_cross_section_field, mat_xs.mesh,
    )
    k_iso = IsotropicScattering(mat_xs) + IsotropicN2N(mat_xs)
    loss = collision - k_iso
    return loss.as_matrix(basis_shape=(mat_xs.mesh.ng, 1))


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
    ng = mix.ng

    # Loss matrix A = C − K_iso and production dyad F = χ ⊗ νΣ_f, both realized
    # from the transport operators on the meshless single cell (the SAME
    # operators the meshed SN solver uses; #276) — the dense (ng,ng) matrices
    # come from the operators' own as_matrix() apply-to-basis materialization,
    # NOT hand-rolled np.diag/np.outer.
    A = _assemble_loss_matrix(mat_xs)
    F = FissionOperator.from_solver_data(mat_xs=mat_xs).as_matrix(
        basis_shape=(ng, 1),
    )

    # k∞ and the flux spectrum φ are the EXACT dominant eigenpair of the dense
    # resolvent A⁻¹F — the direct (non-iterative) realization of the eigenvalue
    # boundary (:func:`~orpheus.numerics.eigenvalue.direct_eigenvalue`), the
    # sibling of the meshed ``power_iteration``.  The 0-D infinite-medium
    # spectrum is exactly solvable, so the dense engine is the right tool, not
    # an iterative approximation.
    k_inf, phi = direct_eigenvalue(A, F)

    # The reaction rates are the §5.6 reaction-rate functional ∫⟨Σx, φ⟩dV
    # (:class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`)
    # on the meshless unit-volume cell — production (νΣf) and absorption (Σa),
    # each the φ†=1 degenerate of the homogenization PG bilinear ⟨φ†, M[Σx]φ⟩.
    # Production is νΣf ONLY: the (n,2n) reaction is a loss-side transfer folded
    # into A as 2Σ₂ᵀ, never a production channel.
    production = IntegratedReactionRate(mat_xs.fission_production_field)
    absorption = IntegratedReactionRate(mat_xs.absorption_cross_section_field)

    # Normalise the flux so the fission production rate νΣf·φ = 100 n/cm³/s.
    phi = phi * (100.0 / production.evaluate(phi.reshape(ng, 1)))

    prod_rate = production.evaluate(phi.reshape(ng, 1))
    abs_rate = absorption.evaluate(phi.reshape(ng, 1))
    total_flux = float(phi.sum())

    if mix.eg is None:
        # Synthetic XS — no physical energy grid, so lethargy / per-energy
        # diagnostics are not defined.  k_inf and the flux spectrum still
        # carry meaningful information; only the per-energy plotting path
        # is unavailable.
        representative_energy = None
        energy_widths = None
        lethargy_widths = None
    else:
        # The group geometry is the EnergyGrid value object's own — the
        # GEOMETRIC group centre (the natural log-axis abscissa) and the
        # energy / lethargy widths — read off ``mixture.energy_grid``, NOT
        # re-derived here (single source of the group structure).
        eg = mix.energy_grid
        representative_energy = eg.representative_energy
        energy_widths = eg.energy_widths
        lethargy_widths = eg.lethargy_widths

    return HomogeneousResult(
        k_inf=k_inf,
        flux=phi,
        representative_energy=representative_energy,
        energy_widths=energy_widths,
        lethargy_widths=lethargy_widths,
        sig_prod=prod_rate / total_flux,
        sig_abs=abs_rate / total_flux,
        mixture=mix,
    )
