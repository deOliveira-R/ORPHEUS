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
  :class:`~orpheus.transport.operators.isotropic_transfer.IsotropicFission`
  — the rank-1 dyad :math:`|\chi\rangle\langle\nu\Sigma_f|` on the scalar
  flux; what the k-eigenvalue outer iteration, the homogeneous
  :math:`K = A^{-1}F`, and the diffusion scalar composite consume;
* the **ANGULAR binding** — this module's :class:`FissionOperator`, the
  frame's :math:`\ell=0` conjugation of that dyad on the posed angular
  composite: the shared lift
  :class:`~orpheus.transport.operators.angular_lift.AngularLift` over the
  fission datum (CS4c step 5, R-1 — ``{S, N₂ₙ} | {F}`` share the ℓ = 0
  base). It DERIVES the energy binding as its middle factor
  (:attr:`~orpheus.transport.operators.angular_lift.AngularLift.isotropic_energy`
  — the :class:`~orpheus.transport.operators.n2n.N2NOperator` /
  :class:`~orpheus.transport.operators.isotropic_transfer.IsotropicN2N`
  relation, at rank 1) from the ONE datum field it carries (#450: the
  operator holds what its ends select, not a second binding).

**Fission is the** :math:`\ell=0` **rank-1 degenerate of the scattering
binding.** Scattering is :math:`S = R\,\Lambda_{\ell\le L}\,M/W`; fission
is :math:`F = R_0\,(|\chi\rangle\langle\nu\Sigma_f|)\,M_0/W` — the same
faces from the same hub-interned frame
(:meth:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame.for_space`,
so S and F posed on one space share ONE metric), one moment instead of
:math:`(L{+}1)^2`, an energy factor of rank 1 instead of a general
transfer stack. That single sentence fixes both action spellings:

* **forward** — the lift base's reaction-rate fast path (the S idiom):
  reduce :math:`\phi` (the angular integral, or the ℓ = 0 slot of a
  moment operand — selected by the binding's ends at construction),
  apply the dyad, embed through the base's producer-side combine (no
  moment tensor on the hot path);
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

from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np

from orpheus.numerics.space import FunctionSpace
from orpheus.transport.operators.angular_lift import AngularLift
from orpheus.transport.operators.bound_operator import BoundOperator
from orpheus.transport.operators.isotropic_transfer import IsotropicFission
from orpheus.transport.operators.lift import interior_space_of

if TYPE_CHECKING:
    from orpheus.numerics.operator import OperatorProduct, TensorProductOperator
    from orpheus.numerics.spaces.moment_head import MomentHead
    from orpheus.transport.material_field import FissionMaterialField
    from orpheus.transport.mesh.material_xs_field import MaterialXSField
    from orpheus.transport.reaction_rate_functional import ReactionRateFunctional

__all__ = ["FissionMomentOperator", "FissionOperator"]


@dataclass(eq=False)
class FissionMomentOperator(BoundOperator):
    r"""The :math:`\ell=0` moment-space fission factor — the dyad on the
    ``[0, 0]`` harmonic block.

    The fission sibling of
    :class:`~orpheus.transport.operators.transfer.LegendreMomentTransfer`
    (:math:`\Lambda`, the per-ℓ factor both transfer terms share):
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

    The arithmetic routes through the energy binding's OWN verbs
    (:meth:`IsotropicFission.apply
    <orpheus.transport.operators.isotropic_transfer.IsotropicFission.apply>`
    / ``apply_transpose`` — the operator, not its kernel: one-level
    delegation, so the energy binding is the ONE dyad home every route
    reaches); this class is pure moment-layout plumbing.
    """

    energy: "IsotropicFission"

    @property
    def is_adjointable(self) -> bool:
        # The transpose leg is spelled (the dual dyad on the [0,0]
        # block); is_invertible inherits base False (rank 1).
        return True

    @property
    def _head(self) -> "MomentHead":
        r"""The angular HEAD (this operator's domain) — which index tuple is the :math:`\ell = 0` slot is ITS to say (#429)."""
        from orpheus.numerics.spaces.moment_head import MomentHead

        if not isinstance(self.domain, MomentHead):
            raise TypeError(
                f"FissionMomentOperator: the moment ends must be an angular "
                f"HEAD space; got {type(self.domain).__name__}."
            )
        return self.domain

    def _admit(self, moments: np.ndarray) -> None:
        head = self._head
        rank = len(head.shape)
        if moments.ndim < rank + 1 or moments.shape[:rank] != tuple(head.shape):
            raise ValueError(
                f"FissionMomentOperator acts on L=0 moment tensors whose "
                f"leading axes are the head {tuple(head.shape)} of "
                f"{head.name!r} (then ng, *spatial); got shape {moments.shape}"
            )

    def apply(self, moments: np.ndarray) -> np.ndarray:
        r"""Dyad on the :math:`\ell=0` block: ``out[iso] = F₀ m[iso]``, ``iso`` the head's isotropic slot."""
        m = np.asarray(moments)
        self._admit(m)
        iso = self._head.isotropic_slot
        out = np.zeros_like(m)
        out[iso] = self.energy.apply(m[iso])
        return out

    def apply_transpose(self, moments: np.ndarray) -> np.ndarray:
        r"""Dual dyad on the :math:`\ell=0` block (the factor swap)."""
        m = np.asarray(moments)
        self._admit(m)
        iso = self._head.isotropic_slot
        out = np.zeros_like(m)
        out[iso] = self.energy.apply_transpose(m[iso])
        return out


@dataclass(eq=False)
class FissionOperator(AngularLift[IsotropicFission]):
    r"""The fission source operator :math:`F` on the angular composite —
    the frame's :math:`\ell=0` conjugation of the energy dyad (module
    docstring): the shared lift
    :class:`~orpheus.transport.operators.angular_lift.AngularLift` over
    the fission datum, and nothing else (#450 — the exact ctor holds the
    datum, the faces, and the ends; the energy binding is DERIVED).

    Parameters
    ----------
    fission : FissionMaterialField
        The fission channel's kernel field — validated per-material
        :math:`(\chi, \nu\Sigma_f)` pairs over the mesh layout. The
        binding's whole datum; :attr:`isotropic_energy` is its energy
        binding on the emitted interior's scalar sub-space (the ONE dyad
        home every arm and the moment factor reach), cached once.
    flux_analysis, source_reconstruction : the two minted :math:`L=0` faces
        (the base's fields) — the hub-interned frame's, so S and F posed
        on one space share ONE metric.
    domain, codomain : FunctionSpace
        The two mandatory ends (kw-only, write-once): the composite
        full-field space, both — the SAME instance the within-group
        siblings carry, so the eigen-pencil pairing and any α-resolvent
        ``OperatorSum`` grouping validate natively. A moment-domain
        sibling (:meth:`~AngularLift.on_moment_domain`) is admissible by
        the base — machinery-first (§16.1): no consumer today, none
        special-cased.
    """

    fission: "FissionMaterialField"

    # ── the lift's subclass contract ─────────────────────────────────

    @property
    def data_ng(self) -> int:
        return self.fission.ng

    def _bind_energy(self, scalar_space: FunctionSpace) -> IsotropicFission:
        return IsotropicFission(
            self.fission, domain=scalar_space, codomain=scalar_space,
        )

    def _frame_form(self) -> "OperatorProduct":
        return self.full_fission_kernel

    @classmethod
    def from_solver_data(
        cls, *, mat_xs: "MaterialXSField", space: "FunctionSpace",
    ) -> "FissionOperator":
        r"""Tier-2 extract-and-mint: the facade's fission channel as the
        datum field, the :math:`L=0` faces through the blessed frame chain
        (the CS5 generator channel inside
        :meth:`HarmonicFrame.for_space
        <orpheus.transport.frames.harmonic_frame.HarmonicFrame.for_space>`
        — which refuses a SCALAR composite, diffusion's, with the
        generator's own message: that consumer wants the ENERGY binding,
        ``IsotropicFission.from_material_xs``), the endomorphic composite
        ends from one ``space=``."""
        from orpheus.transport.frames.harmonic_frame import HarmonicFrame
        from orpheus.transport.material_field import FissionMaterialField

        interior = interior_space_of(space, owner=f"{cls.__name__}.from_solver_data")
        if interior.axes is None:
            raise TypeError(
                "FissionOperator.from_solver_data requires an axis-built "
                "composite FullFieldSpace (the angular interior names the "
                "quadrature the lift embeds over). A scalar consumer wants "
                "the ENERGY binding — IsotropicFission.from_material_xs."
            )
        frame = HarmonicFrame.for_space(interior, 0)
        return cls(
            FissionMaterialField.from_material_xs(mat_xs),
            flux_analysis=frame.flux_analysis_on(interior),
            source_reconstruction=frame.source_reconstruction_on(interior),
            domain=space,
            codomain=space,
        )

    # ── derived structure (single sources, delegated to the energy binding) ──

    @property
    def kernel(self) -> "TensorProductOperator":
        r"""The §5.6 rank-1 TP kernel — DELEGATED to the energy binding
        (one dyad home; consumed by the integral-kernel Protocol gates)."""
        return self.isotropic_energy.kernel

    @property
    def production_rate(self) -> "ReactionRateFunctional":
        r"""The production-rate co-vector — delegated to the energy
        binding (the Pattern-3 criticality diagnostic; see
        :attr:`IsotropicFission.production_rate
        <orpheus.transport.operators.isotropic_transfer.IsotropicFission.production_rate>`)."""
        return self.isotropic_energy.production_rate

    @cached_property
    def full_fission_kernel(self) -> "OperatorProduct":
        r"""The frame form :math:`R_0\,F_0\,M_0` (:math:`R_0\,F_0` on the
        moment end) — fission's whole action as ONE conjugated product
        (the transfer family's
        :attr:`~orpheus.transport.operators.transfer.TransferOperator.full_transfer_kernel`
        sibling at :math:`\ell=0`). The per-ordinate source is
        ``(1/W)·full_fission_kernel.apply(ψ)``; its
        :meth:`~orpheus.numerics.operator.OperatorProduct.apply_transpose`
        is what makes the lift's ``apply_transpose`` factor reversal
        instead of arithmetic. Cached at first access (the satellite
        ruling — the datum is immutable, so the cache cannot go stale).
        The forward production path keeps the reaction-rate fast path
        (the S ruling: the frame form is the validated NON-hot-path
        spelling)."""
        # The ℓ=0 ends are the bound frame's BASIS space — read, never
        # minted from the integer 0 (#429 tracker 2.5): which family spans
        # the moments is the quadrature's decision.
        ends = self._moment_space
        return self._end.conjugate(
            self,
            FissionMomentOperator(self.isotropic_energy, domain=ends, codomain=ends),
        )
