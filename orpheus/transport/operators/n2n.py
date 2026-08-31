r"""The :math:`(n,2n)` source operator on the angular composite — first-class.

**Why its own operator (CS4c §14.1, ruled 2026-08-30).** The
:math:`(n,2n)` channel is scattering-like (a group transfer, in principle
with its own anisotropy) AND production-like (it carries the multiplicity
:math:`\nu_{2n}`): its bundling is CONTEXT-dependent — with :math:`S` when
scattering anisotropy is the interesting axis, with :math:`F` when
production accounting is — and a context-dependent bundling must not be
decided at the operator level. So the within-group algebra spells it
explicitly,

.. math::

    A \;=\; L + C - S - N_{2n} - B,

and any bundling is a solver-side :class:`~orpheus.numerics.operator.OperatorSum`
grouping. (Before the extraction the term hid inside
:class:`~orpheus.transport.operators.scattering.ScatteringOperator`'s iso
accumulator — the operator-level commitment this module retires.)

**What it is, algebraically.** :math:`(n,2n)` emission is isotropic (the
kernel is a single :math:`\ell=0` transfer matrix), so the composite
action is the isotropic lift of the energy operator: with
:math:`K = \nu_{2n}\,\Sigma_{2n}^{T}` per cell and :math:`W` the angular
measure's total weight,

.. math::

    N_{2n}\,\psi \;=\; \frac{1}{W}\, E\, K \int_{S^2} \psi \, d\Omega,

— the frame's :math:`\ell=0` conjugation, realized on the reaction-rate
fast path (no moment tensor; the producer-side :math:`/W` combine is the
shared :func:`~orpheus.transport.operators._per_ordinate.assemble_per_ordinate_isotropic`).
Its Euclidean transpose is the lift's reversal,
:math:`N_{2n}^{T}\chi\big|_{n} = \tfrac{w_n}{W}\, K^{T} \sum_m \chi_m` —
the same :math:`w`-weighted embedding the harmonic frame's analysis face
carries at :math:`\ell=0`.

Future anisotropy is the KERNEL's growth (an ℓ-stack on
:class:`~orpheus.transport.kernels.N2NKernel`), never an S entanglement.

Carrier arms mirror :class:`ScatteringOperator`'s until step 5's arm
deletion: composite ``FullField`` (bulk-only; zero trace), per-ordinate
:class:`~orpheus.transport.fields.angular_flux.AngularFlux`, and the
windowed :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
(whose :math:`\ell=0` moment IS the scalar flux). A bare
:class:`~orpheus.transport.fields.scalar_flux.ScalarFlux` consumer wants
the ENERGY binding (:class:`~orpheus.transport.operators.isotropic_scattering.IsotropicN2N`)
— this operator is the ANGULAR binding of the same datum and refuses
scalar carriers loudly.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import singledispatchmethod
from typing import TYPE_CHECKING, Any, cast, overload

import numpy as np

from orpheus.numerics.operator import BlockRole
from orpheus.numerics.space import FunctionSpace
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.operators._per_ordinate import (
    assemble_per_ordinate_isotropic,
)
from orpheus.transport.operators.bound_operator import BoundOperator
from orpheus.transport.source_sinks import (
    AngularBoundarySourceSink,
    AngularSourceSink,
    ScalarSourceSink,
)

if TYPE_CHECKING:
    from orpheus.transport.mesh.material_xs_field import MaterialXSField
    from orpheus.transport.operators.isotropic_scattering import IsotropicN2N

__all__ = ["N2NOperator"]


@dataclass(eq=False)
class N2NOperator(BoundOperator["FullField"]):
    r"""The :math:`(n,2n)` source operator :math:`N_{2n}` on the angular
    composite — the isotropic lift of the energy binding (module
    docstring).

    Parameters
    ----------
    energy : IsotropicN2N
        The scalar-space energy binding :math:`\nu_{2n}\,\Sigma_{2n}^{T}`
        of the same kernel datum — the ONE arithmetic home (its
        :class:`~orpheus.transport.material_field.N2NMaterialField` verbs
        carry the per-material dispatch and the multiplicity).
    weights : np.ndarray
        The per-ordinate angular weights ``(N,)`` of the binding measure
        — the transpose's :math:`w`-embedding; :attr:`total_weight`
        (:math:`W`) derives from it (single source).
    domain, codomain : FunctionSpace
        The two mandatory ends (kw-only, write-once): the composite
        full-field space, both — the SAME instance the within-group
        siblings carry, so the ``(L+C) − S − N_{2n} − B`` OperatorSum
        guard validates the arm natively.
    """

    energy: "IsotropicN2N"
    weights: np.ndarray

    # (n,2n) emission is volumetric — bulk only, no face-trace action.
    # Class-level constant (unannotated: not a dataclass field).
    block_role = BlockRole.BULK

    def __post_init__(self) -> None:
        w = np.asarray(self.weights, dtype=float) + 0.0
        if w.ndim != 1 or w.size == 0:
            raise ValueError(
                f"N2NOperator requires per-ordinate weights (N,); got "
                f"shape {w.shape}"
            )
        w.setflags(write=False)
        # dataclass(eq=False) is not frozen; plain re-assignment is the
        # canonicalizing write (the ends stay write-once via the base).
        self.weights = w
        self._assert_energy_extent_both_ends(
            self.energy.n2n.ng, operator="N2NOperator",
        )

    @classmethod
    def from_solver_data(
        cls, *, mat_xs: "MaterialXSField", space: "FunctionSpace",
    ) -> "N2NOperator":
        r"""Tier-2 extract-and-mint: the energy binding from the facade's
        :math:`(n,2n)` channel (on the composite's scalar interior), the
        weights from the angular axis's generator (the CS5 channel), the
        endomorphic composite ends from one ``space=``."""
        from orpheus.numerics.quadrature.directional import Quadrature
        from orpheus.numerics.spaces.full_field_space import FullFieldSpace
        from orpheus.transport.material_field import N2NMaterialField
        from orpheus.transport.operators.isotropic_scattering import (
            IsotropicN2N,
        )

        interior = (
            space.interior_space
            if isinstance(space, FullFieldSpace)
            else None
        )
        if interior is None or interior.axes is None:
            raise TypeError(
                "N2NOperator.from_solver_data requires an axis-built "
                "composite FullFieldSpace (the angular interior names the "
                "quadrature the lift embeds over); got "
                f"{type(space).__name__}."
            )
        scalar_space = FunctionSpace.of_axes(*interior.axes[1:])
        energy = IsotropicN2N(
            N2NMaterialField.from_material_xs(mat_xs),
            domain=scalar_space,
            codomain=scalar_space,
        )
        weights = interior.axes[0].generator_as(
            Quadrature, consumer="N2NOperator.from_solver_data",
        ).weights
        return cls(energy, weights, domain=space, codomain=space)

    @property
    def total_weight(self) -> float:
        r""":math:`W = \sum_n w_n` — the binding measure's total angular
        weight (derived; the :math:`/W` of the lift)."""
        return float(self.weights.sum())

    @property
    def is_adjointable(self) -> bool:
        # The Euclidean transpose is spelled (apply_transpose), so the
        # metric-aware .H is free; is_invertible inherits base False —
        # an emission source operator is not invertible.
        return True

    # ── the action (carrier arms mirror S's, until step 5) ────────────

    @singledispatchmethod
    def _apply_impl(self, psi) -> "Any":
        raise TypeError(
            f"N2NOperator.apply: unsupported input type "
            f"{type(psi).__name__}; expected FullField, AngularFlux, or "
            f"HarmonicMomentFlux. (A bare ScalarFlux consumer wants the "
            f"ENERGY binding — IsotropicN2N — not the angular lift.)"
        )

    @_apply_impl.register
    def _(self, psi: FullField) -> "FullField":
        # Bulk-only: the interior dispatches to the per-ordinate arms;
        # the boundary is the implicit zero (volumetric emission).
        combined = self.apply(
            cast("AngularFlux | HarmonicMomentFlux", psi.interior),
        )
        return FullField(
            interior=combined,
            boundary=AngularBoundarySourceSink.zeros(psi.boundary.space),
        )

    @_apply_impl.register
    def _(self, psi: AngularFlux) -> "AngularSourceSink":
        phi = psi.integrate_angular()
        iso = ScalarSourceSink(
            values=cast(np.ndarray, self.energy.apply(phi.values)),
            space=phi.space,
        )
        return assemble_per_ordinate_isotropic(
            iso, None, psi.space, self.total_weight,
        )

    @_apply_impl.register
    def _(self, phi_moments: HarmonicMomentFlux) -> "AngularSourceSink":
        # The ℓ=0 moment IS the scalar flux (Y_0^0 = 1); the per-ordinate
        # target is the posed composite's interior — same route as S's
        # windowed arm.
        interior = self._interior_space()
        if interior.axes is None:
            raise TypeError(
                "N2NOperator windowed arm: the composite interior must be "
                "axis-built to name the scalar sub-space."
            )
        phi = phi_moments.scalar_flux(
            space=FunctionSpace.of_axes(*interior.axes[1:]),
        )
        iso = ScalarSourceSink(
            values=cast(np.ndarray, self.energy.apply(phi.values)),
            space=phi.space,
        )
        return assemble_per_ordinate_isotropic(
            iso, None, interior, self.total_weight,
        )

    @_apply_impl.register
    def _(self, phi: ScalarFlux) -> "Any":
        raise TypeError(
            "N2NOperator.apply: a ScalarFlux consumer wants the ENERGY "
            "binding (IsotropicN2N) — this operator is the ANGULAR "
            "binding of the (n,2n) datum."
        )

    def _interior_space(self) -> "FunctionSpace":
        from orpheus.numerics.spaces.full_field_space import FullFieldSpace

        domain = self.domain
        if (
            not isinstance(domain, FullFieldSpace)
            or domain.interior_space is None
        ):
            raise TypeError(
                "N2NOperator: the bound composite domain carries no "
                "interior space to size the per-ordinate target."
            )
        return domain.interior_space

    if TYPE_CHECKING:
        @overload
        def apply(self, psi: FullField, /) -> "FullField": ...
        @overload
        def apply(self, psi: AngularFlux, /) -> "AngularSourceSink": ...
        @overload
        def apply(
            self, phi_moments: HarmonicMomentFlux, /,
        ) -> "AngularSourceSink": ...
        def apply(self, x: Any, /) -> Any: ...
    else:
        apply = _apply_impl

    @overload
    def apply_transpose(self, chi: "FullField", /) -> "FullField": ...
    @overload
    def apply_transpose(self, chi: np.ndarray, /) -> np.ndarray: ...
    def apply_transpose(self, chi: "Any") -> "Any":
        r"""The Euclidean transpose of the lift:
        :math:`(N_{2n}^{T}\chi)_n = \tfrac{w_n}{W}\,K^{T}\sum_m \chi_m`
        (module docstring) — the :math:`(n,2n)` channel of the adjoint
        within-group loss ``(L+C−S−N_{2n}−B).H``. Bulk-only on the
        composite (the pullback trace is the implicit zero), bare
        ``(N, ng, *spatial)`` ndarray otherwise.
        """
        if isinstance(chi, FullField):
            bulk_bar = self.apply_transpose(np.asarray(chi.interior.values))
            return FullField(
                interior=AngularSourceSink(
                    values=bulk_bar, space=chi.interior.space,
                ),
                boundary=AngularBoundarySourceSink.zeros(chi.boundary.space),
            )
        chi_values = np.asarray(getattr(chi, "values", chi))
        summed = chi_values.sum(axis=0)                     # (ng, *spatial)
        core = self.energy.apply_transpose(summed)          # K^T Σχ
        w = self.weights.reshape(
            (self.weights.size,) + (1,) * core.ndim,
        )
        return (w * core[np.newaxis]) / self.total_weight
