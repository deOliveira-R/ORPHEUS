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
    from orpheus.numerics.operator import LinearOperator
    from orpheus.transport.mesh.material_xs_field import MaterialXSField


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

def _as_dense(op: "LinearOperator", ng: int) -> np.ndarray:
    r"""Realize a meshless group operator as a dense ``(ng, ng)`` matrix.

    For the single-cell (0-D) infinite-medium phase space an energy operator
    acts on the group vector as an ``ng × ng`` matrix; column ``i`` is the
    operator applied to the ``i``-th group basis vector.  The basis columns
    are shaped ``(ng, 1)`` — matching the operators' ``(ng, *spatial)`` bare-
    ndarray contract — so a squeezed ``(ng,)`` column (which would broadcast
    against a ``(ng, 1)`` coefficient into a spurious ``(ng, ng)``) cannot
    occur.  Works on ANY meshless-safe operator, including the composition
    ``C - K_iso`` (an :class:`~orpheus.numerics.operator.OperatorSum`).
    """
    cols = []
    for i in range(ng):
        e_i = np.zeros((ng, 1))
        e_i[i, 0] = 1.0
        cols.append(np.asarray(op.apply(e_i)).ravel())
    return np.column_stack(cols)


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
    ``(ng, ng)`` matrix via :func:`_as_dense`.
    """
    collision = MultiplicationOperator.from_mesh(
        mat_xs.total_cross_section_field, mat_xs.mesh,
    )
    k_iso = IsotropicScattering(mat_xs) + IsotropicN2N(mat_xs)
    loss = collision - k_iso
    return _as_dense(loss, mat_xs.mesh.ng)


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
    # come from apply-to-basis-vectors, NOT hand-rolled np.diag/np.outer.
    A = _assemble_loss_matrix(mat_xs)
    F = _as_dense(FissionOperator.from_solver_data(mat_xs=mat_xs), ng)

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
