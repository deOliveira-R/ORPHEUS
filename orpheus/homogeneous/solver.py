"""Homogeneous infinite-medium reactor eigenvalue solver.

Solves for the neutron spectrum and k-infinity in an infinite homogeneous
medium.  All spatial and angular dependence integrates out; the transport
equation reduces to the pure energy-balance eigenvalue problem

    A φ = (1/k) F φ,
        A = diag(Σ_t) − Σ_s0ᵀ − 2·Σ₂ᵀ,
        F = χ ⊗ νΣ_f,

with k_inf = λ_max(A⁻¹F).

The eigenproblem is spelled in the operator algebra itself (taxonomy step
5b): the loss operator ``A = C − K_iso`` is composed over a MESHLESS
single-cell MaterialMesh (``MaterialMesh.from_materials``) — the collision
diagonal C = diag(Σ_t) minus the model-shared isotropic energy operators
``IsotropicScattering`` (Σ_s0ᵀ) and ``IsotropicN2N`` (2·Σ₂ᵀ); streaming L
is identically zero in an infinite medium and is dropped — and the
multiplication operator is the composition
``K = MatrixInverseOperator(A) @ F`` (one eager LU factorization at
construction; the first production consumer of the dense direct inverse),
whose materialization feeds the shared Perron–Frobenius extraction
:func:`~orpheus.numerics.eigenvalue.dominant_eigenpair`.  The whole
infinite-medium spectrum runs through the SAME operator algebra the
meshed SN solver uses, not a bespoke matrix (#276).

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
from orpheus.numerics.axis import Axis, BasisKind, EnergyAxis
from orpheus.numerics.eigenvalue import dominant_eigenpair
from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator
from orpheus.numerics.space import FunctionSpace
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.operators.fission import FissionOperator
from orpheus.transport.operators.isotropic_scattering import (
    IsotropicN2N,
    IsotropicScattering,
)
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator

if TYPE_CHECKING:
    from orpheus.numerics.operator import OperatorSum
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

def _pose_space(mix: Mixture) -> FunctionSpace:
    r"""The space the infinite-medium problem poses on: Energy ⊗ the quotient point.

    Minted from the MIXTURE — the problem's own physics — never read off
    a carrier (CS4a K2): the energy axis goes through the ONE energy-arm
    rule (:meth:`~orpheus.numerics.axis.EnergyAxis.from_materials`, the
    same rule ``MaterialMesh.bulk_space`` routes through, so the two
    spellings cannot diverge), and the spatial factor is the explicit
    quotient point with the COUNTING weight (``weights=None`` — the
    normalized per-unit-volume density convention; the counting-measure
    premise the rate pairing below rests on).

    The degenerate carrier's ``bulk_space`` mints an ``==`` space (the
    identity bridge gate pins it) but is no longer what production
    consumes — the carrier supplies cross sections, not the posing.
    """
    return FunctionSpace.of_axes(
        EnergyAxis.from_materials([mix]),
        Axis("spatial", (1,), kind=BasisKind.NODAL),
    )


def _assemble_loss_operator(
    mat_xs: "MaterialXSField", space: FunctionSpace
) -> "OperatorSum":
    r"""The loss operator :math:`A = C - K_\mathrm{iso}` for an infinite medium.

    Composes the model-shared transport operators on the meshless single-cell
    carrier — collision :math:`C = M[\Sigma_t]`
    (:class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`)
    minus the isotropic energy transfer :math:`K_\mathrm{iso} = \Sigma_{s0}^T
    + 2\Sigma_2^T`
    (:class:`~orpheus.transport.operators.isotropic_scattering.IsotropicScattering`
    + :class:`~orpheus.transport.operators.isotropic_scattering.IsotropicN2N`).
    Streaming :math:`L` is identically zero in an infinite medium and dropped.

    Returned UN-materialized (an
    :class:`~orpheus.numerics.operator.OperatorSum`) — the consumer chooses
    the realization.  :func:`solve_homogeneous_infinite` hands it to
    :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`,
    whose constructor materializes it through the operator's own
    :meth:`~orpheus.numerics.operator.LinearOperator.as_matrix` (the
    apply-to-basis materialization this module's retired ``_as_dense``
    prototyped, promoted at taxonomy §12 step 5) and LU-factors it once
    (the early ``.as_matrix()`` this function performed until taxonomy
    step 5b moved into that constructor).

    Since campaign 1 (CS1) the meshless operators pose on a REAL space,
    and since CS4a (K2) that space is the CALLER's — the mixture-minted
    Energy ⊗ point from :func:`_pose_space`, threaded explicitly into
    all three arms (``C`` by direct construction, no ``from_mesh``
    default chain), so ``C − (IsoS + IsoN2N)`` agrees on one space and
    the ``OperatorSum`` guard VALIDATES the sum instead of skipping it.
    Consumers no longer pass ``basis_shape=(ng, 1)`` by hand:
    ``as_matrix``/``MatrixInverseOperator`` derive it from the threaded
    domain (the pre-CS1 idiom existed only because these operators
    carried no
    :attr:`~orpheus.numerics.operator.LinearOperator.domain` to derive
    it from). (Until K2 the space was read off the carrier's
    ``bulk_space``; the carrier now supplies cross sections only.)
    """
    collision = MultiplicationOperator(
        coefficient=mat_xs.total_cross_section_field, space=space,
    )
    k_iso = IsotropicScattering(mat_xs, space=space) + IsotropicN2N(
        mat_xs, space=space
    )
    return collision - k_iso


def solve_homogeneous_infinite(mix: Mixture) -> HomogeneousResult:
    r"""Solve the infinite-medium eigenvalue problem for a homogeneous mixture.

    Spells the multiplication operator :math:`\mathbf{K} =
    \mathbf{A}^{-1}\mathbf{F}` in the operator algebra itself —
    ``K = MatrixInverseOperator(loss) @ production`` — from the loss
    operator :math:`\mathbf{A} = C - K_\mathrm{iso} =
    \operatorname{diag}(\Sigma_t) - \Sigma_{s0}^{T} - 2\Sigma_2^{T}`
    (model-shared transport operators over a meshless single-cell
    :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh`) and the
    fission production dyad :math:`\mathbf{F} = \chi \otimes \nu\Sigma_f`,
    then returns the dominant eigenpair of the materialized
    :math:`\mathbf{K}`: :math:`k_\infty = \lambda_{\max}` and the flux
    spectrum :math:`\varphi` (the corresponding right eigenvector),
    normalised so the fission production rate
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
    # The posing: the mixture-minted Energy ⊗ point space (the problem's
    # own physics names its space); the meshless carrier supplies the
    # cross sections — one cell, one region, no streaming.
    space = _pose_space(mix)
    mat_xs = MaterialMesh.from_materials({0: mix}).material_xs_field()
    ng = mix.ng

    # The multiplication operator K = A⁻¹F, spelled in the operator algebra
    # itself from the transport operators on the meshless single cell (the
    # SAME operators the meshed SN solver uses; #276): the loss operator
    # A = C − K_iso, inverted DIRECTLY by MatrixInverseOperator — one eager
    # LU factorization at construction, the realization the exactly-solvable
    # 0-D problem earns (the structure-keyed ``loss.inverse()`` would return
    # the ITERATIVE GreenOperator splitting; constructing the matrix inverse
    # explicitly IS the strategy choice) — composed with the fission
    # production dyad F = χ ⊗ νΣ_f.
    loss = _assemble_loss_operator(mat_xs, space)
    production = FissionOperator.from_solver_data(mat_xs=mat_xs, space=space)
    K = MatrixInverseOperator(loss) @ production

    # k∞ and the flux spectrum φ are the EXACT dominant eigenpair of the
    # materialized K, extracted by the shared Perron–Frobenius primitive
    # (:func:`~orpheus.numerics.eigenvalue.dominant_eigenpair`, the one home
    # of the complex-rejection + sign convention).  The 0-D infinite-medium
    # spectrum is exactly solvable, so the dense direct engine is the right
    # tool, not an iterative approximation.
    k_inf, phi = dominant_eigenpair(K.as_matrix())

    # The reaction rates are the SPACE's pairing ⟨Σx, φ⟩ — production
    # (νΣf) and absorption (Σa), each the φ†=1 degenerate of the
    # homogenization PG bilinear ⟨φ†, M[Σx]φ⟩. On the quotient point the
    # measure is COUNTING (the normalized per-unit-volume density
    # convention _pose_space states), so the pairing IS the bare group
    # contraction — and it reads the POSING's measure, never a carrier's
    # (CS4a K2: the carrier's volume_measure is un-wired from this path).
    # Production is νΣf ONLY: the (n,2n) reaction is a loss-side transfer
    # folded into A as 2Σ₂ᵀ, never a production channel.
    nu_sig_f = np.asarray(mat_xs.fission_production_field.values)
    sig_a = np.asarray(mat_xs.absorption_cross_section_field.values)

    # Normalise the flux so the fission production rate νΣf·φ = 100 n/cm³/s.
    phi = phi * (100.0 / space.inner_product(nu_sig_f, phi.reshape(ng, 1)))

    prod_rate = space.inner_product(nu_sig_f, phi.reshape(ng, 1))
    abs_rate = space.inner_product(sig_a, phi.reshape(ng, 1))
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
