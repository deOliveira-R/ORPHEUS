r"""The fission source operator on the angular composite — the frame's :math:`\ell=0` conjugation of the fission dyad.

This module owns the ANGULAR binding of the fission channel in the
operator-algebra view of the Boltzmann transport equation,

.. math::

    (L + C - S - N_{2n} - B)\,\psi = q
    \qquad\text{(fixed source)}

.. math::

    (L + C - S - N_{2n} - B)\,\psi = \tfrac{1}{k}\,F\,\psi
    \qquad\text{(eigenvalue)}

where :math:`L` is the :math:`\sigma`-free **streaming** leaf, :math:`C =
M[\sigma_t]` the **collision** diagonal, :math:`S` the scattering
operator, :math:`N_{2n}` the first-class :math:`(n,2n)` source
(:mod:`orpheus.transport.operators.n2n`, the CS4c §14.1 extraction),
:math:`B` the boundary gain, and :math:`F` the **fission emission
operator**

.. math::

    (F\,\phi)_g(\vec r) = \chi_g\,
    \sum_{g'} \nu\Sigma_{f,g'}(\vec r)\,\phi_{g'}(\vec r).

**The two bindings of one datum (CS4c step 4, §16.2).** The fission
channel's representation-free datum is the per-material factor pair
:class:`~orpheus.transport.kernels.FissionKernel` :math:`(\chi,
\nu\Sigma_f)`, held as a
:class:`~orpheus.transport.material_field.FissionMaterialField`. It is
bound at exactly two space kinds:

* the **ENERGY binding**
  :class:`~orpheus.transport.operators.isotropic_scattering.IsotropicFission`
  — the rank-1 dyad :math:`|\chi\rangle\langle\nu\Sigma_f|` on the scalar
  flux; what the k-eigenvalue outer iteration, the homogeneous
  :math:`K = A^{-1}F`, and the diffusion scalar composite consume;
* the **ANGULAR binding** — this module's :class:`FissionOperator`, the
  frame's :math:`\ell=0` conjugation of that dyad on the posed angular
  composite. It retains the energy binding as its middle factor (the
  :class:`~orpheus.transport.operators.n2n.N2NOperator` /
  :class:`~orpheus.transport.operators.isotropic_scattering.IsotropicN2N`
  relation, at rank 1).

**Fission is the** :math:`\ell=0` **rank-1 degenerate of the scattering
binding.** Scattering is :math:`S = R\,\Lambda_{\ell\le L}\,M/W`; fission
is :math:`F = R_0\,(|\chi\rangle\langle\nu\Sigma_f|)\,M_0/W` — the same
faces from the same hub-interned frame
(:meth:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame.for_space`,
so S and F posed on one space share ONE metric), one moment instead of
:math:`(L{+}1)^2`, an energy factor of rank 1 instead of a general
transfer stack. That single sentence fixes both action spellings:

* **forward** — the reaction-rate fast path (the S idiom): reduce
  :math:`\phi = \int\psi\,d\Omega`, apply the dyad, embed through the
  shared producer-side combine
  :func:`~orpheus.transport.operators._per_ordinate.assemble_per_ordinate_isotropic`
  (no moment tensor on the hot path);
* **Euclidean transpose** — the frame form: the cached
  :attr:`~FissionOperator.full_fission_kernel` product
  ``frame.conjugate(FissionMomentOperator)`` reversed by
  :meth:`~orpheus.numerics.operator.OperatorProduct.apply_transpose`,
  then the producer :math:`/W` — factor reversal, ZERO fission-side
  :math:`w`-arithmetic. (The pre-step-4 hand-rolled
  ``np.multiply.outer(w, ·)/W`` transpose died here; it divided by
  :math:`W` before the dyad where the product divides after — a pure
  IEEE-754 order change, principled-equivalent, gated at tolerance.)

**The adjoint is the machinery's, end to end (§16.2 Riesz maximal use).**
``F.apply_transpose`` is factor reversal; the metric Hilbert adjoint
``F.H`` composes :math:`\sharp_V \circ F^{\mathsf T} \circ \flat_W` from
the bound spaces' own first-class Riesz legs
(:attr:`~orpheus.numerics.space.FunctionSpace.riesz_raise` /
:attr:`~orpheus.numerics.space.FunctionSpace.riesz_lower`) — nothing
fission-specific anywhere on the adjoint path. The kernel level stays
dagger-free by ruling (§6): ``FissionKernel(chi=νΣf, nu_sig_f=χ)`` is
refused by its own simplex guard — the adjoint exists only where the two
metrics do, i.e. on the BOUND operator, by theorem.

**Pencil / resolvent readiness (§16.1).** With mandatory composite ends
:math:`F` is a peer of the within-group loss on ONE space: the
k-eigenvalue pencil pairs :math:`(A, F)` under one discipline, and an
α-resolvent that shifts :math:`F` to the loss side composes it into the
``OperatorSum`` under the same ends guard — no adapter, no special case.
The rank-1 structure is why power iteration converges geometrically
(dominance ratio :math:`|k_1/k_0|`): :math:`F` has one energy mode per
cell because :math:`\chi` is a rank-1 spectrum.

The §5.6 Kernel reading — fission as a rank-1 integral kernel
=============================================================

In the grand-report §5.6 suffix law fission is a **Kernel**: a nonlocal
operator integrating the flux against a measure on the group axis.
:attr:`FissionOperator.kernel` exposes that integral structure (the
2-factor :class:`~orpheus.numerics.operator.TensorProductOperator`
``outer(χ, production_rate) & Identity``, delegated to the energy
binding — ONE arithmetic home) and satisfies the
:class:`~orpheus.transport.operators.integral_kernel_operator.IntegralKernelOperator`
Protocol. Its transpose is the **dual dyad**
:math:`|\nu\Sigma_f\rangle\langle\chi|` — the χ↔νΣf role swap, a theorem
of the rank-1 ``outer`` primitive, never re-derived here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property, singledispatchmethod
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
    from orpheus.numerics.operator import (
        OperatorProduct,
        TensorProductOperator,
    )
    from orpheus.transport.frames.harmonic_frame import HarmonicFrame
    from orpheus.transport.mesh.material_xs_field import MaterialXSField
    from orpheus.transport.operators.isotropic_scattering import (
        IsotropicFission,
    )
    from orpheus.transport.reaction_rate_functional import (
        ReactionRateFunctional,
    )

__all__ = ["FissionMomentOperator", "FissionOperator"]


@dataclass(eq=False)
class FissionMomentOperator(BoundOperator):
    r"""The :math:`\ell=0` moment-space fission factor — the dyad on the
    ``[0, 0]`` harmonic block.

    The fission sibling of
    :class:`~orpheus.transport.operators.scattering.LegendreMomentScattering`
    (:math:`\Lambda`) and
    :class:`~orpheus.transport.operators.scattering.N2NMomentOperator`:
    endomorphic on the :math:`L=0` spherical-harmonic coefficient space,
    applying the energy dyad to the single :math:`Y_0^0` block (fission
    emission is isotropic by construction — every :math:`\ell\ge 1`
    block is zero, and here there are none).

    **Single consumer, stated honestly:** the middle factor of
    :attr:`FissionOperator.full_fission_kernel` (the transpose product).
    The forward hot path rides the reaction-rate fast path and never
    builds a moment tensor. Bare-``ndarray`` legs only — the
    :class:`~orpheus.numerics.operator.OperatorProduct` chain composes
    on raw values; a typed-carrier surface lands with a typed consumer
    (defer-until-consumer).

    The arithmetic routes through the retained energy binding's cached
    :attr:`~orpheus.transport.operators.isotropic_scattering.IsotropicFission.kernel`
    — ONE dyad home; this class is pure moment-layout plumbing.
    """

    energy: "IsotropicFission"

    @property
    def is_adjointable(self) -> bool:
        # The transpose leg is spelled (the dual dyad on the [0,0]
        # block); is_invertible inherits base False (rank 1).
        return True

    def _admit(self, moments: np.ndarray) -> None:
        if moments.ndim < 3 or moments.shape[0] != 1 or moments.shape[1] != 1:
            raise ValueError(
                f"FissionMomentOperator acts on L=0 moment tensors "
                f"(1, 1, ng, *spatial); got shape {moments.shape}"
            )

    def apply(self, moments: np.ndarray) -> np.ndarray:
        r"""Dyad on the :math:`\ell=0` block: ``out[0,0] = F₀ m[0,0]``."""
        m = np.asarray(moments)
        self._admit(m)
        out = np.zeros_like(m)
        out[0, 0] = self.energy.kernel.apply(m[0, 0])
        return out

    def apply_transpose(self, moments: np.ndarray) -> np.ndarray:
        r"""Dual dyad on the :math:`\ell=0` block (the factor swap)."""
        m = np.asarray(moments)
        self._admit(m)
        out = np.zeros_like(m)
        out[0, 0] = self.energy.kernel.apply_transpose(m[0, 0])
        return out


@dataclass(eq=False)
class FissionOperator(BoundOperator["FullField"]):
    r"""The fission source operator :math:`F` on the angular composite —
    the frame's :math:`\ell=0` conjugation of the energy dyad (module
    docstring).

    Parameters
    ----------
    energy : IsotropicFission
        The scalar-space energy binding of the same kernel datum — the
        ONE arithmetic home of the dyad (its cached
        :attr:`~orpheus.transport.operators.isotropic_scattering.IsotropicFission.kernel`
        serves every arm and the moment factor).
    frame : HarmonicFrame
        The :math:`L=0` hub-interned frame of the posed angular space —
        OPERATIVE state, not provenance: the transpose product
        (:attr:`full_fission_kernel`) and the binding measure
        (:attr:`total_weight`) read it. Minted by tier 2 through the
        blessed chain, so it is the SAME frame object every sibling
        binding on this space reaches.
    domain, codomain : FunctionSpace
        The two mandatory ends (kw-only, write-once): the composite
        full-field space, both — the SAME instance the within-group
        siblings carry, so the eigen-pencil pairing and any α-resolvent
        ``OperatorSum`` grouping validate natively.
    """

    energy: "IsotropicFission"
    frame: "HarmonicFrame" = field(kw_only=True)

    # Fission emission is volumetric — bulk only, no face-trace action
    # (χ lives on cell-centred volumes). Issue #208 / Wave O; encoded as
    # the class-level BlockRole constant (unannotated: not a field).
    block_role = BlockRole.BULK

    def __post_init__(self) -> None:
        self._assert_energy_extent_both_ends(
            self.energy.fission.ng, operator="FissionOperator",
        )
        # Frame/space agreement: the frame's measure must be the posed
        # angular axis's (N ordinates) — a frame minted elsewhere would
        # make the transpose product and the forward combine silently
        # inconsistent (the S face-binding guard, one level up).
        interior = self._interior_space()
        n_ordinates = int(interior.shape[0])
        n_frame = int(np.asarray(self.frame.measure.weights).shape[0])
        if n_frame != n_ordinates:
            raise TypeError(
                f"FissionOperator: the frame carries {n_frame} ordinates "
                f"but the bound interior has {n_ordinates} — mint the "
                f"frame from the SAME posed space (tier 2 does)."
            )
        # Energy-binding agreement: the middle factor must be bound on
        # exactly this composite's scalar sub-space shape.
        scalar_shape = tuple(interior.shape[1:])  # (ng, *spatial)
        energy_shape = tuple(self.energy._bulk_scalar_space().shape)
        if energy_shape != scalar_shape:
            raise TypeError(
                f"FissionOperator: the energy binding rides shape "
                f"{energy_shape} but this composite's scalar sub-space "
                f"is {scalar_shape} — bind both from one space= (tier 2 "
                f"does)."
            )

    @classmethod
    def from_solver_data(
        cls, *, mat_xs: "MaterialXSField", space: "FunctionSpace",
    ) -> "FissionOperator":
        r"""Tier-2 extract-and-mint: the energy binding from the facade's
        fission channel (on the composite's scalar interior), the
        :math:`L=0` frame through the blessed chain (the CS5 generator
        channel inside
        :meth:`HarmonicFrame.for_space
        <orpheus.transport.frames.harmonic_frame.HarmonicFrame.for_space>`),
        the endomorphic composite ends from one ``space=``."""
        from orpheus.numerics.spaces.full_field_space import FullFieldSpace
        from orpheus.transport.frames.harmonic_frame import HarmonicFrame
        from orpheus.transport.material_field import FissionMaterialField
        from orpheus.transport.operators.isotropic_scattering import (
            IsotropicFission,
        )

        interior = (
            space.interior_space
            if isinstance(space, FullFieldSpace)
            else None
        )
        if interior is None or interior.axes is None:
            raise TypeError(
                "FissionOperator.from_solver_data requires an axis-built "
                "composite FullFieldSpace (the angular interior names the "
                "quadrature the lift embeds over); got "
                f"{type(space).__name__}. A scalar consumer wants the "
                "ENERGY binding — IsotropicFission.from_material_xs."
            )
        # Frame FIRST: for_space demands a Quadrature-generated leading
        # axis, so a SCALAR composite (diffusion's — leading axis is
        # energy) refuses here with the generator's own message before
        # any energy-binding shape confusion can occur.
        frame = HarmonicFrame.for_space(interior, 0)
        scalar_space = FunctionSpace.of_axes(*interior.axes[1:])
        energy = IsotropicFission(
            FissionMaterialField.from_material_xs(mat_xs),
            domain=scalar_space,
            codomain=scalar_space,
        )
        return cls(energy, frame=frame, domain=space, codomain=space)

    # ── derived structure (single sources) ───────────────────────────

    @property
    def total_weight(self) -> float:
        r""":math:`W = \sum_n w_n` — the binding measure's total angular
        weight (read off the retained frame's measure; the producer-side
        :math:`/W`)."""
        return float(np.asarray(self.frame.measure.weights).sum())

    @property
    def is_adjointable(self) -> bool:
        # The Euclidean transpose is the reversed face product (below);
        # ``.H`` composes the Riesz legs around it. is_invertible
        # inherits base False — a rank-1 production operator is singular.
        return True

    @property
    def kernel(self) -> "TensorProductOperator":
        r"""The §5.6 rank-1 TP kernel — DELEGATED to the energy binding
        (one dyad home; consumed by
        :class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicEmission`
        and the integral-kernel Protocol gates)."""
        return self.energy.kernel

    @property
    def production_rate(self) -> "ReactionRateFunctional":
        r"""The production-rate co-vector — delegated to the energy
        binding (the Pattern-3 criticality diagnostic; see
        :attr:`IsotropicFission.production_rate
        <orpheus.transport.operators.isotropic_scattering.IsotropicFission.production_rate>`)."""
        return self.energy.production_rate

    @cached_property
    def full_fission_kernel(self) -> "OperatorProduct":
        r"""The frame form :math:`R_0\,F_0\,M_0` — fission's whole action
        as ONE conjugated product (the S
        :attr:`~orpheus.transport.operators.scattering.ScatteringOperator.full_scatter_kernel`
        sibling at :math:`\ell=0`). The per-ordinate source is
        ``(1/W)·full_fission_kernel.apply(ψ)``; its
        :meth:`~orpheus.numerics.operator.OperatorProduct.apply_transpose`
        is what makes :meth:`apply_transpose` factor reversal instead of
        arithmetic. Cached at first access (the satellite ruling — the
        kernel field is immutable, so the cache cannot go stale). The
        forward production path keeps the reaction-rate fast path (the S
        ruling: the frame form is the validated NON-hot-path spelling)."""
        # The ℓ=0 ends are the bound frame's BASIS space — read, never
        # minted from the integer 0 (#429 tracker 2.5): which family spans
        # the moments is the quadrature's decision.
        ends = self.frame.basis.space
        return self.frame.conjugate(
            FissionMomentOperator(self.energy, domain=ends, codomain=ends),
        )

    def _interior_space(self) -> "FunctionSpace":
        from orpheus.numerics.spaces.full_field_space import FullFieldSpace

        domain = self.domain
        if (
            not isinstance(domain, FullFieldSpace)
            or domain.interior_space is None
        ):
            raise TypeError(
                "FissionOperator: the bound composite domain carries no "
                "interior space to size the per-ordinate target."
            )
        return domain.interior_space

    # ── the action (carrier arms mirror N2N's, until step 5) ──────────

    @singledispatchmethod
    def _apply_impl(self, psi) -> "Any":
        raise TypeError(
            f"FissionOperator.apply: unsupported input type "
            f"{type(psi).__name__}; expected FullField, AngularFlux, or "
            f"HarmonicMomentFlux. (A scalar consumer — the k-outer, "
            f"diffusion, homogeneous — wants the ENERGY binding, "
            f"IsotropicFission, not the angular lift.)"
        )

    @_apply_impl.register
    def _(self, psi: FullField) -> "FullField":
        # Bulk-only: the interior dispatches to the per-ordinate arms;
        # the boundary is the implicit zero (volumetric emission). The
        # eigenvalue posing stacks this row over the ray fold — see
        # sn/solver's pose site.
        combined = self.apply(
            cast("AngularFlux | HarmonicMomentFlux", psi.interior),
        )
        return FullField(
            interior=combined,
            boundary=AngularBoundarySourceSink.zeros(psi.boundary.space),
        )

    @_apply_impl.register
    def _(self, psi: AngularFlux) -> "AngularSourceSink":
        # φ = ∫ψ dΩ, dyad, producer-side /W — the reaction-rate fast
        # path (bit-identical to the pre-step-4 composite arm; the frame
        # form full_fission_kernel is the transpose's spelling).
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
        # target is the posed composite's interior — same route as the
        # N2N windowed arm.
        interior = self._interior_space()
        if interior.axes is None:
            raise TypeError(
                "FissionOperator windowed arm: the composite interior "
                "must be axis-built to name the scalar sub-space."
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
            "FissionOperator.apply: a ScalarFlux consumer wants the "
            "ENERGY binding (IsotropicFission) — this operator is the "
            "ANGULAR binding of the fission datum."
        )

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
        r"""Adjoint fission source — factor reversal of the frame form:
        :math:`F^{\mathsf T}\psi^* = (R_0 F_0 M_0)^{\mathsf T}\psi^*/W
        = M_0^{\mathsf T} F_0^{\mathsf T} R_0^{\mathsf T}\,\psi^*/W`,
        i.e. the :math:`w`-weighted embedding of the dual dyad applied
        to the ordinate sum — but SPELLED as
        ``full_fission_kernel.apply_transpose(ψ*)/W``, so the identity
        is the :class:`~orpheus.numerics.operator.OperatorProduct`
        chain's, not this module's.

        Physically: high importance :math:`\psi^*` in the emission
        groups :math:`\chi` makes a cell a strong adjoint source
        weighted by its production :math:`\nu\Sigma_f`. This is the
        Euclidean transpose (L12); the metric Hilbert adjoint ``.H``
        composes the spaces' Riesz legs around it (§16.2 — nothing
        hand-rolled on the adjoint path).

        The COMPOSITE (``FullField``) arm mirrors the forward: fission
        is PURE BULK, so the transpose reads only the bulk cotangent
        and emits the implicit-zero trace, letting the daggered
        eigen-pencil compose through ``OperatorSum.apply_transpose``.
        Bare ``(N, ng, *spatial)`` ndarray otherwise.
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
        return (
            np.asarray(self.full_fission_kernel.apply_transpose(chi_values))
            / self.total_weight
        )

