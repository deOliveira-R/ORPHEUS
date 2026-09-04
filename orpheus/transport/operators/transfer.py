r"""The transfer family's shared cores — the moment factor :math:`\Lambda` and the angular binding.

A **transfer channel** is a Legendre stack with a yield —
:class:`~orpheus.transport.kernels.TransferKernel`
:math:`(\{\Sigma_{c,\ell}\}, y_c)` — and the two collision-gain terms of
the within-group algebra

.. math::

    A \;=\; L + C - S - N_{2n} - B

are two instances of ONE object: the scattering gain :math:`S`
(:math:`y = 1`) and the first-class :math:`(n,2n)` gain :math:`N_{2n}`
(:math:`y = 2`), each the angular binding of its channel's field at the
solve's Legendre order. This module owns the bindings' shared arithmetic
(#426 step 2, ruled 2026-09-03 — the F2/F3 rulings in
``.claude/plans/n2n_anisotropy_kept.md`` §2b):

* :class:`LegendreMomentTransfer` — the per-ℓ block-diagonal moment
  factor :math:`\Lambda = y \sum_\ell \mathbf{P}_\ell \otimes
  \Sigma_{c,\ell}` on the frame's coefficient space;
* :class:`TransferOperator` — the angular binding :math:`T = R\,\Lambda\,M
  / W` on the posed composite, realized on the reaction-rate fast path
  for :math:`\ell = 0` (the scalar energy binding
  :class:`~orpheus.transport.operators.isotropic_scattering.IsotropicTransfer`,
  lifted through the shared producer-side combine) plus the
  frame-conjugated :math:`\ell \ge 1` redistribution.

**The kernel tier names the mathematical object; the operator tier names
the TERM.** The two terms of the algebra are thin role subclasses —
:class:`~orpheus.transport.operators.scattering.ScatteringOperator` and
:class:`~orpheus.transport.operators.n2n.N2NOperator` — whose only
content is the extraction classmethod (which ``Mixture`` channel the
tier-2 mint reads) and the role name; every verb lives HERE, once. An AST
gate (``tests/transport/test_transfer_roles.py``) asserts the roles define
nothing else, so the twin path the carve removed cannot regrow one
override at a time.

Until 2026-09-04 the :math:`(n,2n)` binding re-spelled this class's arms
with ``aniso = None`` and a frame minted at :math:`L = 0`: the operator
tier imposed a :math:`P_0` model on a channel whose tape stores the same
seven Legendre moments as elastic — `[M]` −414 Δk·10⁵ on a Be-reflected
fast slab (issue #426; ``docs/theory/methods/sn/adjoint.rst``
§sn-n2n-p0-truncation). The measured size is what retired the twin.

The theory lives in the book — one concept, one home:

* P0 in-scatter, the (n,2n) source, the Pℓ reconstruction, and the
  :math:`1/W` normalisation chain —
  ``docs/theory/methods/sn/slab_multigroup.rst §mg-inscatter-source``,
  §pn-scatter, §n2n-source, §pn-scatter-rlm.
* The no-prefactor :math:`Y_\ell^m` convention and the Funk–Hecke
  eigenbasis (why :math:`T = R\circ\Lambda\circ M`) —
  ``docs/theory/foundations/spherical_harmonics.rst
  §spherical-harmonics-eigenbasis``.
* The §5.6 integral-kernel reading and the apply-only capability
  surface — ``docs/theory/foundations/operator_algebra.rst
  §integral-kernel-category``, §capability-set-semantics.
* The Euclidean adjoint :math:`T^{T}` (forward fast-path vs adjoint
  frame-form) — ``docs/theory/methods/sn/adjoint.rst
  §sn-scattering-adjoint-source``, §sn-n2n-adjoint-source.

Capability surface
==================

``apply`` + ``apply_transpose`` (``is_adjointable=True``); **no**
``solve`` (``is_invertible=False``). A transfer gain is rank
:math:`O(N_{\text{cells}}\cdot N_{\text{groups}})` with no tractable
inverse — it is *applied*, never *inverted*, and the missing ``solve``
is structural method-absence (a composer refuses to build :math:`T^{-1}`
at construction time), not an advertising flag. The adjoint :math:`T^{T}`
rides the harmonic-frame :attr:`~TransferOperator.full_transfer_kernel`
(closes `#118 <https://github.com/deOliveira-R/ORPHEUS/issues/118>`_):
see :meth:`~TransferOperator.apply_transpose`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property, singledispatchmethod
from typing import TYPE_CHECKING, Any, ClassVar, Self, cast, overload

import numpy as np

from orpheus.numerics.space import FunctionSpace
from orpheus.numerics.basis.base import TruncatedBasis
from orpheus.numerics.spaces.moment_head import MomentHead
from orpheus.numerics.spaces.full_field_space import FullFieldSpace
from orpheus.transport.frames import (
    HarmonicAnalysisOperator,
    HarmonicFrame,
    HarmonicReconstructionOperator,
)
from orpheus.numerics.operator import (
    BlockRole,
    LinearOperator,
    OperatorProduct,
)

# Runtime imports of the flux / source types — required at module load
# time because :func:`singledispatchmethod.register` dispatches on the
# runtime class.  These modules form a leaf in the transport package
# dependency graph (they do not import the operators), so the imports
# are circular-import-safe.
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from orpheus.transport.source_sinks import (
    ScalarSourceSink,
    AngularSourceSink,
    AngularBoundarySourceSink,
    HarmonicMomentSourceSink,
)
from orpheus.transport.full_field import FullField
from orpheus.transport.kernels import TransferKernel
from orpheus.transport.material_field import TransferMaterialField
from orpheus.transport.operators._per_ordinate import (
    assemble_per_ordinate_isotropic,
)
from orpheus.transport.operators.bound_operator import BoundOperator
from orpheus.transport.operators.isotropic_scattering import IsotropicTransfer

if TYPE_CHECKING:
    from orpheus.transport.mesh.material_xs_field import MaterialXSField


__all__ = ["LegendreMomentTransfer", "TransferOperator"]


def _moment_head_of(space: "FunctionSpace | None", owner: str) -> MomentHead:
    """Narrow a moment operator's end to the MomentHead surface, loudly."""
    if not isinstance(space, MomentHead):
        raise TypeError(
            f"{owner}: the moment ends must be an angular HEAD space "
            f"(SphericalHarmonicSpace or LegendreSpace); got "
            f"{type(space).__name__}."
        )
    return space


@dataclass(eq=False)
class LegendreMomentTransfer(BoundOperator):
    r"""Per-ℓ block-diagonal transfer :math:`\Lambda` on the moment space.

    The diagonal spectrum of the sum-of-tensor-products form
    :math:`\Lambda = y\sum_{\ell=0}^{L} \mathbf{P}_\ell \otimes
    \Sigma_{c,\ell}` (:math:`\mathbf{P}_\ell` selects the :math:`\ell`-th
    harmonic block, :math:`\Sigma_{c,\ell}` the per-material per-ℓ
    Legendre matmul on the group axis, :math:`y` the channel's yield) —
    see ``docs/theory/foundations/operator_algebra.rst
    §scattering-as-tensor-product-sum``.

    For an input moment field :math:`\phi_\ell^m(\vec r)` (the head's
    axes leading, then ``g``, then the spatial axes), the action is
    per-group-transfer per :math:`\ell` block,

    .. math::

        (\Lambda \phi)_\ell^m(\vec r)\bigg|_{g}
        \;=\; y \sum_{g'} \Sigma_{c,\ell}^{m(\vec r)}(g' \to g)\,
              \phi_\ell^m(\vec r)\bigg|_{g'},

    with :math:`m(\vec r)` the material id at cell :math:`\vec r` (the
    per-material structure IS the bound datum — see below).

    **The CS4c rebind (design record §14), generic in the channel (#426
    step 2):** the operator holds the representation-free datum — a
    :class:`~orpheus.transport.material_field.TransferMaterialField`
    already :meth:`at_order
    <orpheus.transport.material_field.TransferMaterialField.at_order>`
    the binding's order — plus its two mandatory ends (the
    :class:`~orpheus.transport.operators.bound_operator.BoundOperator`
    base; :math:`\Lambda` is endomorphic on the coefficient space of its
    order). ``L`` is DERIVED from the field (the order IS the field's —
    single source); the per-material dispatch and the yield live on the
    field's :meth:`~orpheus.transport.material_field.TransferMaterialField.moment_source`
    verb, whose shape guard refuses a moments tensor at any other order.
    The :math:`(n,2n)` instance is this SAME class over the channel's
    field: until 2026-09-04 a twin (``N2NMomentOperator``) spelled its
    :math:`\ell = 0` block alone.

    ``skip_l0`` (default ``True``) skips the :math:`\ell = 0` block, which
    the project's P0 emission handles on a separate reaction-rate fast
    path. Set ``False`` for the full :math:`R\Lambda M\psi` composition on
    the LinearOperator surface — an ℓ-range selector inside the ONE datum,
    not a path switch.

    Capability set: ``{apply, apply_transpose}``; no efficient ``solve``
    (rank-deficient on the :math:`\ell = 0` block by design).
    :math:`\Lambda^{T}` (:meth:`apply_transpose`) is the per-ℓ group-axis
    transpose — the ONLY group-asymmetric factor of the kernel
    :math:`R\circ\Lambda\circ M`, so
    :math:`(R\circ\Lambda\circ M)^{T} = M^{T}\circ\Lambda^{T}\circ R^{T}`
    (``docs/theory/methods/sn/adjoint.rst §sn-scattering-adjoint-source``).

    Parameters
    ----------
    transfer : TransferMaterialField
        The per-material Legendre stacks (with their yield) over the mesh
        layout, at this binding's order.
    skip_l0 : bool, default ``True``
        Skip the :math:`\ell = 0` block (handled by the P0 fast path). Set
        ``False`` for the full :math:`R \Lambda M \psi` composition.
    domain, codomain : FunctionSpace
        The two mandatory ends (kw-only, write-once — the base) — the
        coefficient space of the field's order, both, for the shipped
        endomorphic binding.
    """

    transfer: "TransferMaterialField"
    skip_l0: bool = True

    @classmethod
    def from_field(
        cls,
        transfer: "TransferMaterialField",
        basis: TruncatedBasis,
        *,
        skip_l0: bool = True,
    ) -> "LegendreMomentTransfer":
        r"""Tier-2 mint: bring a channel's field to the basis's order and
        bind the endomorphic ends on the BASIS's coefficient space
        (``basis.space`` supplying both — the endomorphism sugar lives
        HERE, never on the exact ctor). The basis is the single source of
        the moment space (#429 tracker 2.5): an integer cannot say which
        family — full harmonics on a sphere rule, Legendre on a 1-D rule
        — so the caller hands the basis its frame bound, and the ends are
        that basis's, never re-minted."""
        ends = basis.space
        return cls(
            transfer.at_order(basis.L),
            skip_l0=skip_l0,
            domain=ends,
            codomain=ends,
        )

    @property
    def _head(self) -> MomentHead:
        r"""The angular HEAD the moment verbs read the layout from — this
        operator's own domain (the bound basis's coefficient space)."""
        return _moment_head_of(self.domain, type(self).__name__)

    @property
    def L(self) -> int:
        r"""The Legendre truncation :math:`L` — DERIVED from the bound
        field (the order is the field's; single source)."""
        return self.transfer.order

    @property
    def is_adjointable(self) -> bool:
        # Λ exposes its group-transpose Σ_{c,ℓ}^T (apply_transpose), so the
        # metric-aware .H is free. is_invertible
        # inherits base False — a per-ℓ source map is not invertible.
        return True

    @overload
    def apply(self, moments: "HarmonicMomentFlux") -> "HarmonicMomentSourceSink": ...

    @overload
    def apply(self, moments: np.ndarray) -> np.ndarray: ...

    def apply(
        self, moments: "np.ndarray | HarmonicMomentFlux",
    ) -> "np.ndarray | HarmonicMomentSourceSink":
        r"""Apply :math:`\Lambda` to a moment field — the **role-changing** edge.

        :math:`\Lambda` maps a flux moment to the emission **source**
        moment it produces (flux → source); the typed arm makes that role
        change explicit in the signature. (Why the role change lives on
        the operator and not on the frame:
        ``docs/theory/foundations/operator_algebra.rst
        §integral-kernel-category``.)

        Parameters
        ----------
        moments : np.ndarray or HarmonicMomentFlux
            Flux moment field with the head's axes leading, then the group
            axis, then the spatial axes. On the rectangular harmonic head
            the :math:`m`-axis is the addition-theorem-shifted index where
            slot ``l + m`` holds the :math:`(\ell, m)` entry; entries
            outside :math:`|m| \le \ell` are conventionally zero.

            Typed
            :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
            → typed
            :class:`~orpheus.transport.source_sinks.harmonic_moment_source_sink.HarmonicMomentSourceSink`
            (the emitted moment SOURCE) with matching ``L`` / ``mesh`` /
            spatial-moment width.  Bare ndarray → ndarray (the endomorphic
            moment-space view the :math:`R\circ\Lambda\circ M`
            kernel ``OperatorProduct`` composes on).

        Returns
        -------
        np.ndarray or HarmonicMomentSourceSink
            The emitted moment source, same shape as ``moments``.  The
            :math:`\ell = 0` block is zero when ``skip_l0`` is ``True``;
            otherwise the P0 contribution is included.  Typed in →
            typed source out; ndarray in → ndarray out.

        Notes
        -----
        Both arms route through the field's single per-material verb
        :meth:`~orpheus.transport.material_field.TransferMaterialField.moment_source`;
        they differ only in the carrier wrap.
        """
        if isinstance(moments, HarmonicMomentFlux):
            out_values = self.transfer.moment_source(
                moments.values, skip_l0=self.skip_l0, head=self._head,
            )
            # flux moment → source moment: the explicit role change
            # (CS4b S4 — same space, new class; role is class identity).
            return HarmonicMomentSourceSink(
                values=out_values, space=moments.space, L=moments.L,
                spatial_moments=moments.spatial_moments,
            )
        return self.transfer.moment_source(
            moments, skip_l0=self.skip_l0, head=self._head,
        )

    def apply_transpose(
        self, moments: "np.ndarray | HarmonicMomentSourceSink",
    ) -> "np.ndarray | HarmonicMomentFlux":
        r"""Apply :math:`\Lambda^{T}` — the per-ℓ group-transpose (the role-REVERSING edge).

        The bare Euclidean transpose of :meth:`apply`: :math:`\Lambda^{T}` maps a
        source moment back into the flux-moment space it came from (source →
        flux), transposing the per-ℓ group-transfer
        :math:`\Sigma_{c,\ell}(g'\to g) \mapsto (g\to g')`.  Routes through the
        field's transpose verb
        :meth:`~orpheus.transport.material_field.TransferMaterialField.moment_source_transpose`.

        This is the Euclidean transpose, **not** the metric Hilbert adjoint
        :math:`\Lambda^{\dagger} = G^{-1}\Lambda^{T}G` (the
        :attr:`~orpheus.numerics.operator.LinearOperator.H` wrapper's job). As the
        only group-asymmetric factor of the kernel :math:`R\circ\Lambda\circ M`,
        it is what lets the whole kernel transpose fall out for free — see
        ``docs/theory/methods/sn/adjoint.rst §sn-scattering-adjoint-source``.

        Typed
        :class:`~orpheus.transport.source_sinks.harmonic_moment_source_sink.HarmonicMomentSourceSink`
        → :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
        (the explicit role reversal); bare ndarray → ndarray (the endomorphic
        moment-space view the ``OperatorProduct.apply_transpose`` chain composes on).
        """
        if isinstance(moments, HarmonicMomentSourceSink):
            out_values = self.transfer.moment_source_transpose(
                moments.values, skip_l0=self.skip_l0, head=self._head,
            )
            return HarmonicMomentFlux(
                values=out_values, space=moments.space, L=moments.L,
                spatial_moments=moments.spatial_moments,
            )
        return self.transfer.moment_source_transpose(
            moments, skip_l0=self.skip_l0, head=self._head,
        )


@dataclass(eq=False)
class TransferOperator(BoundOperator["FullField"]):
    r"""The angular binding of a transfer channel: :math:`T = R\,\Lambda\,M/W` (P0 + Pℓ).

    **The CS4c step-3 binding (design record §14), generic in the channel
    (#426 step 2):** the exact ctor retains the representation-free datum
    and the minted products, and nothing richer —

    * :attr:`transfer` — the
      :class:`~orpheus.transport.material_field.TransferMaterialField`
      (per-material Legendre stacks with their yield, over the mesh
      layout), already at this binding's order (the order IS
      :attr:`scattering_order` — single source);
    * :attr:`flux_analysis` / :attr:`source_reconstruction` — the two
      typed faces minted from the HUB-interned
      :class:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame`
      (tier 2 mints them and forgets the frame; the :attr:`frame`
      accessor is PROVENANCE, riding on the faces);
    * the two mandatory ends (kw-only, write-once —
      :class:`~orpheus.transport.operators.bound_operator.BoundOperator`):
      the composite full-field space, both — the SAME instance
      ``L``/``C``/``B`` carry, so the within-group
      ``(L + C) − S − N₂ₙ − B`` OperatorSum guard validates each gain arm
      natively. (⚠ the 2-D windowed operand carries a MOMENT interior —
      the shipped non-endomorphism the step-0 census measured; the
      carrier dispatch serves it until step 5's arm deletion.)

    **Two terms, one binding.** The scattering gain and the
    :math:`(n,2n)` gain are the role subclasses
    :class:`~orpheus.transport.operators.scattering.ScatteringOperator`
    and :class:`~orpheus.transport.operators.n2n.N2NOperator`, built on
    the SAME posed space at the SAME solve order through the same
    interned frame (their ``from_solver_data`` mints); the within-group
    algebra spells ``− S − N₂ₙ`` explicitly (§14.1) and bundles them as
    the solver sees fit. Every verb below is channel-agnostic: the yield
    enters the P0 fast path and the moment factor from the field's own
    datum, and the production accounting a :math:`y > 1` channel adds is
    arithmetic that vanishes for :math:`y = 1`.

    Capability surface: ``{apply, apply_transpose}`` — no efficient
    ``solve``; the adjoint :math:`T^{T}` is free via the harmonic-frame
    :attr:`full_transfer_kernel` (see :meth:`apply_transpose`).
    """

    transfer: "TransferMaterialField"
    #: The minted FLUX analysis face :math:`M \otimes I` on the posed
    #: interior (``AngularFlux → HarmonicMomentFlux``) — bound at tier 2,
    #: retained. Consumers: the windowed bulk projection and the S6
    #: adjoint gates.
    flux_analysis: "HarmonicAnalysisOperator[AngularFlux, HarmonicMomentFlux]" = field(
        kw_only=True,
    )
    #: The minted SOURCE reconstruction face :math:`R \otimes I` landing
    #: on the posed interior (``HarmonicMomentSourceSink →
    #: AngularSourceSink``) — the windowed arm's typed R.
    source_reconstruction: "HarmonicReconstructionOperator[HarmonicMomentSourceSink, AngularSourceSink]" = field(
        kw_only=True,
    )

    #: The P0 ENERGY binding this term lifts — the role subclass names its
    #: own (``IsotropicScattering`` for :math:`S`, ``IsotropicN2N`` for
    #: :math:`N_{2n}`); the core mints the shared
    #: :class:`~orpheus.transport.operators.isotropic_scattering.IsotropicTransfer`.
    #: A ClassVar, not a field: it is the role's identity, not a datum.
    isotropic_binding: ClassVar[type[IsotropicTransfer]] = IsotropicTransfer

    # A transfer gain is a BULK operator — the moment-folding
    # `Σ_c · ⟨P_ℓ, ψ⟩` reads and writes the bulk flux only (A_bb), no
    # boundary action. Class-level constant (unannotated so the dataclass
    # does not treat it as a field).
    block_role = BlockRole.BULK

    def __post_init__(self) -> None:
        # Per-END energy admission (CS4c §1).
        self._assert_energy_extent_both_ends(
            self.transfer.ng, operator=type(self).__name__,
        )
        # Face-binding agreement: both minted faces must be bound to THIS
        # binding's interior — a face bound elsewhere would make the
        # windowed arm and the per-ordinate combine silently inconsistent.
        # (The XD-1 wrong-EMBEDDING controls stay spellable: a doctored
        # face carries the RIGHT spaces with the WRONG measure.)
        interior = self._interior_space
        if self.flux_analysis.domain != interior:
            raise TypeError(
                f"{type(self).__name__}: the flux-analysis face is bound to "
                f"a different angular space than this binding's interior "
                f"— mint the faces from the SAME posed space (tier 2 does)."
            )

    @classmethod
    def from_field(
        cls,
        transfer: "TransferMaterialField",
        *,
        scattering_order: int,
        space: "FunctionSpace",
    ) -> Self:
        r"""Tier-2 mint (CS4c §14): bring a channel's field to the solve's
        order, mint the two faces from the HUB-interned frame, and bind
        the endomorphic composite ends from one ``space=``. The role
        subclasses' ``from_solver_data`` extract the channel and call
        this.

        ``space`` is the composite
        :attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.full_field_space`
        the solver threads — MANDATORY since the flip (the ends are
        write-once fields; the OperatorSum guard validates every build).
        The quadrature is reached through the space's angular axis (the
        CS5 generator channel inside
        :meth:`HarmonicFrame.for_space
        <orpheus.transport.frames.harmonic_frame.HarmonicFrame.for_space>`)
        — no ``quadrature=`` parameter survives, so a frame/space metric
        mismatch is unspellable. Both terms of the algebra mint at the
        SAME ``(rule, L)`` and therefore share ONE interned frame.
        """
        interior = (
            space.interior_space
            if isinstance(space, FullFieldSpace)
            else None
        )
        if interior is None:
            raise TypeError(
                f"{cls.__name__}.from_solver_data requires the posed "
                f"composite FullFieldSpace (its interior is the angular "
                f"space the faces bind); got "
                f"{type(space).__name__}."
            )
        frame = HarmonicFrame.for_space(interior, scattering_order)
        return cls(
            transfer.at_order(scattering_order),
            flux_analysis=frame.flux_analysis_on(interior),
            source_reconstruction=frame.source_reconstruction_on(interior),
            domain=space,
            codomain=space,
        )

    @property
    def is_adjointable(self) -> bool:
        # T = R∘Λ∘M exposes its Euclidean transpose T^T via
        # :attr:`full_transfer_kernel` (the OperatorProduct adjoint chain);
        # is_invertible inherits base False —
        # an emission source operator is not invertible.
        return True

    @property
    def scattering_order(self) -> int:
        r"""Maximum Legendre order :math:`L` retained — DERIVED from the
        bound field (the order IS the field's; single source). ``0``
        means P0 only. The solve's ``scattering_order`` for both terms:
        the :math:`(n,2n)` stack is brought to the scattering stack's
        order (ruling O-1), so the two bindings share one frame."""
        return self.transfer.order

    @property
    def is_isotropic(self) -> bool:
        r"""``True`` iff this binding's :math:`\Lambda_{\ell\ge 1}` is the
        zero operator — order 0, or every moment above :math:`\ell = 0`
        exactly zero (an absent section, an ``NL = 1`` evaluation, a stack
        padded to the solve's order). The anisotropic arms then emit
        nothing and are skipped: the same statement ``scattering_order ==
        0`` used to make from the SHAPE, now made from the VALUES — the
        result is bit-identical (an all-zero :math:`\Lambda` reconstructs
        to exact zeros) and the :math:`R\Lambda M` product is not run."""
        return self.transfer.is_isotropic

    @property
    def total_weight(self) -> float:
        r""":math:`W = \int_{S^2} d\Omega` — the binding measure's total
        angular weight (the producer-side :math:`/W`). Read off the
        retained faces' frame MEASURE: the measure is the binding's
        metric (operative data — unlike the :attr:`frame` accessor,
        which is provenance)."""
        return float(self.flux_analysis.frame.measure.weights.sum())

    @property
    def _moment_space(self) -> "FunctionSpace":
        r"""The coefficient space of the bound frame's BASIS — the
        endomorphic ends of the internally-minted moment factors.

        READ off the retained faces' frame (``frame.basis.space``), never
        minted from :attr:`scattering_order`: which family spans the
        moments is the quadrature's decision (full harmonics on a sphere
        rule, Legendre on a 1-D rule), and the frame already carries it.
        The continuum-metric space (the basis's own), not the frame's
        Parseval-dressed ``basis_space``: `[M]` #429 tracker 2.5 the two are
        ``(name, shape)``-equal and metric-DIFFERENT on 33 of 33 shipped
        (rule, L) rows (the per-:math:`\ell` ratio is exactly
        :math:`[(2\ell+1)/4\pi]^2`), and under the continuum end the
        factor's Hilbert adjoint is its transpose EXACTLY
        (:math:`\Lambda^* = \Lambda^{\mathsf T}`, 0.0 on 33/33) while the
        dressed end would move it on 10 of 33 rows (every 1-D and folded
        rule at :math:`\ell \ge 1`) — binding the basis's own space is what
        keeps every number and every ``.H`` bit-identical to the
        ``from_L(L)`` mint it replaces. (Equality is ``(name, shape)`` and
        cannot see the fork; the gate asserts the metric ARRAY.)"""
        return self.flux_analysis.frame.basis.space

    def _moment_transfer(self, *, skip_l0: bool) -> LegendreMomentTransfer:
        r"""Mint the moment-space :math:`\Lambda` factor on this binding's
        datum + moment ends — the ONE internal spelling (three consumers:
        the §5.6 kernel, the full conjugation, and the aniso moment
        route; the windowed arm consumes the cached :attr:`kernel`
        factors)."""
        ends = self._moment_space
        return LegendreMomentTransfer(
            self.transfer, skip_l0=skip_l0, domain=ends, codomain=ends,
        )

    @property
    def frame(self) -> HarmonicFrame:
        r"""PROVENANCE accessor (design record §2 — retirement-tracked):
        the HUB-interned frame the retained faces were minted from,
        riding on :attr:`flux_analysis` (zero extra state). Production
        reads the FACES and the kernel-field datum, never this; it stays
        for provenance and prototyping until proven removable.
        """
        return self.flux_analysis.frame

    @property
    def _interior_space(self) -> "FunctionSpace":
        r"""The posed composite's interior — the angular field space of
        the binding (the mints' input at tier 2; the per-ordinate target
        of the windowed arm). Mandatory since the CS4c flip: the ends
        always carry it, and a non-composite binding refuses loudly.
        """
        domain = self.domain
        if (
            not isinstance(domain, FullFieldSpace)
            or domain.interior_space is None
        ):
            raise TypeError(
                f"{type(self).__name__} binds the composite full-field "
                f"space — this instance's domain carries no interior to "
                f"size the angular arms."
            )
        return domain.interior_space

    @cached_property
    def isotropic_energy(self) -> IsotropicTransfer:
        r"""The P0 ENERGY binding of this operator's own datum — the
        scalar-space :attr:`isotropic_binding` (the role's own:
        :class:`~orpheus.transport.operators.isotropic_scattering.IsotropicScattering`
        under :math:`S`,
        :class:`~orpheus.transport.operators.isotropic_scattering.IsotropicN2N`
        under :math:`N_{2n}`) the per-ordinate fast path lifts — and the
        solver's K_iso assembly consumes: ``K_iso = S.isotropic_energy +
        N2N.isotropic_energy``, composed at the ONE within-group
        construction site (§14.1; :math:`\ell = 0` by physics — the ray
        seed is the scalar flux's).

        Bound to the composite interior's scalar sub-space; holds the
        field at order 0 (the P0 head — the datum it reads and nothing
        richer; :math:`y` travels in the kernels). Cached at
        construction-time semantics (the kernel field is immutable).
        """
        interior = self._interior_space
        if interior.axes is None:
            raise TypeError(
                f"{type(self).__name__}.isotropic_energy: the composite "
                f"interior must be axis-built to name the scalar "
                f"sub-space."
            )
        scalar_space = FunctionSpace.of_axes(*interior.axes[1:])
        return type(self).isotropic_binding(
            self.transfer.at_order(0),
            domain=scalar_space,
            codomain=scalar_space,
        )

    # ── Anisotropic emission: the moment→source map R·Λ_{ℓ≥1} ──────────

    def _aniso_source_from_moment_values(
        self, moment_values: np.ndarray,
    ) -> np.ndarray:
        r"""Reconstruct the per-ordinate :math:`\ell\ge 1` emission from a
        flux-moment tensor — the **moment → source** half
        :math:`R\,\Lambda_{\ell\ge 1}` of the anisotropic composition
        :math:`T_{\rm aniso} = \tfrac{1}{W}\,R\,\Lambda\,M`.

        Built as ``frame.reconstruct_after(Λ)`` — the §5.6 :attr:`kernel`
        sub-operator for inputs ALREADY in moment space, with
        :math:`\Lambda_{\ell\ge 1}` = :class:`LegendreMomentTransfer`
        (``skip_l0=True``) and :math:`R` the :attr:`frame`'s
        :attr:`~orpheus.numerics.frame.FrameBase.reconstruction` face. The
        trailing :math:`1/W` is **not** applied here — it is the producer-side
        normalisation at the :meth:`apply` boundary. Its sole caller is the
        windowed moment-iterate :meth:`apply` arm (:math:`M` already done). The
        two aniso realisations (this :math:`R\Lambda` vs the full-angular
        :attr:`kernel` :math:`R\Lambda M`) share the frame's :math:`R` and
        :math:`\Lambda` and are pinned at 0 ULP by
        ``tests/sn/operators/test_scattering_kernel_crosscheck.py`` — see
        ``docs/theory/foundations/operator_algebra.rst §integral-kernel-category``.

        Parameters
        ----------
        moment_values : np.ndarray
            Flux-moment tensor (head axes, ``ng``, spatial) — typically
            ``M.apply(psi)`` or the windowed iterate's
            :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
            ``.values``.

        Returns
        -------
        np.ndarray
            ``(N, ng, *spatial)`` per-ordinate :math:`\ell\ge 1` emission,
            **pre** :math:`1/W`.
        """
        redistribution = self._moment_transfer(skip_l0=True)
        # R∘Λ as ONE typed operator (M already applied — the moments ARE M·ψ).
        return self.frame.reconstruct_after(redistribution).apply(moment_values)

    @cached_property
    def kernel(self) -> LinearOperator:
        r"""The §5.6 integral kernel :math:`R \circ \Lambda_{\ell\ge 1} \circ M`.

        The anisotropic Legendre redistribution as a typed
        :class:`~orpheus.numerics.operator.OperatorProduct`
        ``frame.conjugate(Λ)`` ``= OperatorProduct(R, OperatorProduct(Λ, M))``,
        so ``kernel.apply(ψ.values) = R(Λ(M ψ))``. The factors:

        * :math:`M` = ``frame.analysis`` (the :attr:`frame`'s
          :attr:`~orpheus.numerics.frame.FrameBase.analysis` face);
        * :math:`\Lambda_{\ell\ge 1}` = :class:`LegendreMomentTransfer`
          over this binding's field (``skip_l0=True``);
        * :math:`R` = ``frame.reconstruction`` (the :attr:`frame`'s
          :attr:`~orpheus.numerics.frame.FrameBase.reconstruction` face).

        This is the production :math:`\ell\ge 1` map (not a parallel semantic
        reading): :meth:`build_aniso_source` is ``(1/W)·kernel`` and the windowed
        arm consumes the sub-operator ``frame.reconstruct_after(Λ)``. The
        producer-side :math:`1/W` lives OUTSIDE the kernel (at the :meth:`apply`
        boundary), so ``kernel.apply`` returns the source **pre**-:math:`1/W`.
        With it, :class:`TransferOperator` satisfies the
        :class:`~orpheus.transport.operators.integral_kernel_operator.IntegralKernelOperator`
        Protocol — the theory (why scattering IS a nonlocal integral kernel,
        why P0 is the local component) is in
        ``docs/theory/foundations/operator_algebra.rst §integral-kernel-category``.
        CACHED at first access (the kernel field is immutable).

        Raises
        ------
        ValueError
            If ``scattering_order == 0`` — a P0-only operator has no anisotropic
            kernel (the :math:`\ell\ge 1` redistribution is empty); the P0
            emission is the LOCAL component handled by :meth:`add_iso_source`.

        Returns
        -------
        LinearOperator
            The typed ``R∘Λ∘M`` anisotropic redistribution.
        """
        if self.scattering_order == 0:
            raise ValueError(
                f"{type(self).__name__}.kernel requires scattering_order >= 1: "
                f"an isotropic (P0-only) operator has no anisotropic integral "
                f"kernel (the R∘Λ∘M angular redistribution is empty). The P0 "
                f"emission is the LOCAL component, handled by add_iso_source."
            )
        # R ∘ Λ ∘ M as ONE typed operator: the frame conjugates Λ (per-ℓ
        # moment-space transfer) between its analysis (M) and reconstruction (R)
        # faces. Λ carries real spaces (== frame.basis_space), so the
        # OperatorProduct composability guard validates R∘Λ∘M natively — NO cast.
        return self.frame.conjugate(
            self._moment_transfer(skip_l0=True)
        )

    @cached_property
    def full_transfer_kernel(self) -> OperatorProduct:
        r"""The FULL emission kernel :math:`R\circ\Lambda_{\ell\ge 0}\circ M`.

        The COMPLETE P0 + anisotropic emission as ONE frame-conjugated
        operator: the isotropic ℓ=0 transfer and the anisotropic ℓ≥1
        redistribution (one :class:`LegendreMomentTransfer`,
        ``skip_l0=False``) conjugated by the frame. The per-ordinate
        source is ``(1/W)·full_transfer_kernel.apply(ψ)``; its transpose
        ``(1/W)·full_transfer_kernel.apply_transpose(ψ*)`` is the adjoint
        :math:`T^{T}` (:meth:`apply_transpose`). Riding the same frame
        conjugation for iso and aniso is what lets the whole transpose
        fall out for free —
        ``docs/theory/methods/sn/adjoint.rst §sn-scattering-adjoint-source``.
        Under :math:`N_{2n}` this is the product whose reversal the
        §sn-n2n-adjoint-source equation states (its :math:`\ell = 0`
        block is the lift's reversal; the :math:`\ell \ge 1` blocks are
        the anisotropy's, since #426 step 2).

        Distinct from :attr:`kernel` (the §5.6 ℓ≥1 ANISOTROPIC subcomponent + the
        0-ULP crosscheck oracle): this is the FULL ℓ≥0 emission.
        CACHED at first access (CS4c §14.7 — the satellite mint drops
        from once-per-apply to once-per-construction; the kernel field is
        immutable, so the cache cannot go stale).
        """
        return self.frame.conjugate(
            self._moment_transfer(skip_l0=False)
        )

    def _assemble_per_ordinate_source(
        self,
        phi: "ScalarFlux",
        aniso_or_none: "AngularSourceSink | None",
        angular_space: "FunctionSpace",
    ) -> "AngularSourceSink":
        r"""Combine the P0 emission (from the scalar flux :math:`\phi_0`)
        with the pre-:math:`/W` :math:`\ell\ge 1` aniso emission into the
        per-ordinate source :math:`(\text{iso}/W) + \text{aniso}`.

        The **single source of truth for the producer-side** :math:`1/W`
        **combine**: every ``apply`` arm that emits a per-ordinate
        :class:`~orpheus.transport.source_sinks.AngularSourceSink` (the
        full-angular :class:`AngularFlux` arm and the windowed
        :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
        arm) routes here — the :math:`/W` convention lives HERE, once. The
        normalisation chain is in ``docs/theory/methods/sn/slab_multigroup.rst``
        (the normalization-chain section).

        Parameters
        ----------
        phi : ScalarFlux
            Scalar flux :math:`\phi_0` (iso magnitude) driving the P0 emission.
        aniso_or_none : AngularSourceSink or None
            Per-ordinate :math:`\ell\ge 1` source ALREADY in per-ordinate
            magnitude (post-:math:`/W`), or ``None`` for an isotropic binding.
        angular_space : FunctionSpace
            The per-ordinate target space (CS4b S4 — sizes the zero
            accumulator; the caller holding the pose supplies it: the
            operand's own space on the full-angular arm, the posed
            composite's interior on the windowed moment arm).
        """
        # The isotropic P0 energy emission through this binding's own
        # scalar-space energy operator; the iso rides the driving flux's
        # own space (width travels in it), and the producer-side /W
        # combine is the ONE shared home
        # (:func:`~orpheus.transport.operators._per_ordinate.assemble_per_ordinate_isotropic`
        # — FissionOperator rides the same primitive).
        iso: ScalarSourceSink = ScalarSourceSink(
            values=cast(np.ndarray, self.isotropic_energy.apply(phi)),
            space=phi.space,
        )
        return assemble_per_ordinate_isotropic(
            iso, aniso_or_none, angular_space, self.total_weight,
        )

    # ── In-place helpers (preserve bit-identity vs SNSolver pre-Wave-D) ─

    @overload
    def add_iso_source(
        self, Q: "ScalarSourceSink", phi: "np.ndarray | ScalarFlux",
    ) -> "ScalarSourceSink": ...

    @overload
    def add_iso_source(
        self, Q: np.ndarray, phi: "np.ndarray | ScalarFlux",
    ) -> None: ...

    def add_iso_source(
        self,
        Q: "np.ndarray | ScalarSourceSink",
        phi: "np.ndarray | ScalarFlux",
    ) -> "ScalarSourceSink | None":
        r"""Add the P0 emission :math:`y\,\Sigma_{c,0}^T\phi` to :math:`Q`.

        Vectorised by material: per cell ``c`` of material ``mid``,
        ``Q[:, ic, jc] += y · p0[mid].T @ phi[:, ic, jc]`` where
        ``p0[mid]`` is ``(ng, ng)`` indexed ``[g_from, g_to]``.

        Typed-action overload:

        * **Raw-in (legacy)** — ``Q: np.ndarray`` is mutated in place and
          ``None`` is returned.
        * **Typed-in (return-new)** — a frozen ``Q: ScalarSourceSink`` stays
          immutable; a NEW :class:`ScalarSourceSink` carrying
          ``Q.values + emission`` is returned (the caller spells the algebra
          ``Q = transfer.add_iso_source(Q, phi)``).

        Parameters
        ----------
        Q : np.ndarray or ScalarSourceSink
            Isotropic source carrier.  Raw ``(ng, nx, ny)`` ndarray is
            mutated in place; typed :class:`ScalarSourceSink` returns a
            new instance.
        phi : np.ndarray or ScalarFlux
            Scalar flux.  Either form is accepted; the underlying
            values are unwrapped before the per-material dispatch.

        Returns
        -------
        np.ndarray or ScalarSourceSink or None
            Raw-in returns ``None`` (legacy in-place); typed-in returns
            a fresh :class:`ScalarSourceSink`.

        Notes
        -----
        The per-material dispatch lives inside
        :meth:`~orpheus.transport.material_field.TransferMaterialField.add_p0_source`.
        """
        phi_values = phi.values if isinstance(phi, ScalarFlux) else phi
        if isinstance(Q, ScalarSourceSink):
            Q_values = Q.values.copy()
            self.transfer.add_p0_source(Q_values, phi_values)
            # Preserve Q's spatial-moment width (#240 D5b-S3 — the φ̂ accumulator
            # carries the trailing 2^d axis; re-wrapping without it would lose
            # the typed factor and raise on the shape).
            return ScalarSourceSink(values=Q_values, space=Q.space)
        self.transfer.add_p0_source(Q, phi_values)
        return None

    def build_aniso_source(
        self,
        angular_flux: "AngularFlux | None",
    ) -> "AngularSourceSink | None":
        r"""Build the per-ordinate Pℓ (:math:`\ell \ge 1`) emission.

        Implements the Galerkin reconstruction :eq:`pn-scatter` from the
        angular-flux moments :eq:`flux-moments` as the literal operator
        composition :math:`Q^{\rm aniso}_n = \tfrac{1}{W}\,(R\,\Lambda\,M\,
        \psi)_n` — i.e. ``(1/W)·`` :attr:`kernel`. The trailing :math:`1/W`
        is the producer-side per-ordinate projection (the source enters the
        sweep already in per-ordinate magnitude, so the sweep does NOT apply
        ``/W`` again). The full derivation — M/Λ/R faces, the addition-theorem
        :math:`(2\ell+1)` factor, and the :math:`1/W` normalisation chain —
        is in ``docs/theory/methods/sn/slab_multigroup.rst §pn-scatter-rlm``
        and ``docs/theory/foundations/spherical_harmonics.rst``.

        Parameters
        ----------
        angular_flux : AngularFlux or None
            Angular flux shape ``(N, ng, nx, ny)`` from the most recent sweep,
            or ``None`` on the first source iteration before any sweep has run.
            Only the typed :class:`AngularFlux` is accepted (it carries the mesh
            the output :class:`AngularSourceSink` requires).

        Returns
        -------
        AngularSourceSink or None
            ``(N, ng, nx, ny)`` per-ordinate :math:`\ell \ge 1` contribution in
            **per-ordinate magnitude** (the trailing ``/W`` is applied here).
            Returns ``None`` when the binding :attr:`is_isotropic` (order 0,
            or exactly zero above :math:`\ell = 0`) or ``angular_flux is None``.
        """
        if self.is_isotropic or angular_flux is None:
            return None
        # T_aniso = (1/W)·kernel: the §5.6 :attr:`kernel` is the R∘Λ∘M
        # redistribution; the producer-side /W lives OUTSIDE it (applied here).
        out_values = self.kernel.apply(angular_flux.values) / self.total_weight
        # RΛM is spatial-moment-axis-agnostic (#240 D5b-S3): the source
        # rides the iterate's own space (CS4b S4 — same space, new role),
        # so the moment factor travels with it.
        return AngularSourceSink(values=out_values, space=angular_flux.space)

    # ── Foldable / residual split ─────────────────────────────────────
    #
    # T = T_foldable + T_residual, additive at rtol=1e-14: T_foldable is the
    # P0 within-group self-transfer (diagonal y·Σ_{c,0}^{g→g}, foldable into
    # the removal cross-section σ_r = σ_t − y·Σ_{c,0}^{g→g}); T_residual
    # carries everything else (cross-group P0, all Pℓ≥1). Generic in the
    # yield (#426 step 2, ruling §4.4): the SI splitting reads S's; folding
    # the (n,2n) within-group block is a later, measured decision. Data API
    # only — no solver/sweep/iteration consumes these methods yet; the
    # intended consumer is a consistent DSA preconditioner (#2). Theory (why
    # each residual piece is unfoldable):
    # docs/theory/methods/sn/loss_representation.rst §loss-rep-removal-sigma.
    #
    # ⚠ LATENT CORRECTNESS TRAP (#215): do NOT wire the σ_r-SWEEP as the
    # within-group A_wg.inverse(). The σ_r-sweep inverts a DIAGONAL-in-angle
    # removal, but T_foldable = Σ_c0·P_iso is the ISOTROPIC-PROJECTION
    # self-transfer — the two coincide ONLY for isotropic flux, so the wiring
    # ships 46–56 % silent flux errors on anisotropic problems. Use consistent
    # DSA (#2) or Krylov (splitting-invariant, already production). Any
    # within-group accelerator MUST be gated on an ANISOTROPIC config — the
    # isotropic box cannot see this error. Full failure table:
    # docs/theory/methods/sn/slab_one_group.rst §si-sigma-r-fold-mismatch.

    def _sibling(
        self, field_: "TransferMaterialField",
    ) -> Self:
        r"""A sibling binding of ``field_`` on the SAME ends and of the SAME
        role — faces re-minted from the HUB at the sibling field's own
        order (the interned frame chain, so an order-0 sibling gets
        order-0 faces).
        """
        interior = self._interior_space
        frame = HarmonicFrame.for_space(interior, field_.order)
        return type(self)(
            field_,
            flux_analysis=frame.flux_analysis_on(interior),
            source_reconstruction=frame.source_reconstruction_on(interior),
            domain=self.domain,
            codomain=self.codomain,
        )

    def foldable_part(self) -> Self:
        r"""Return the P0 within-group self-transfer sibling of :math:`T`.

        Carries only the diagonal of each material's :math:`\Sigma_{c,0}`
        — the within-group self-transfer cross-section
        :math:`\Sigma_{c,0}^{g\to g}` per energy group (the yield travels
        in the sibling's kernels, so the sibling emits
        :math:`y\,\Sigma_{c,0}^{g\to g}`). All other channels (cross-group
        P0, every :math:`P_\ell \ge 1`) live in :meth:`residual_part`.

        Returns
        -------
        Self
            An order-0 sibling (of this operator's role) whose
            per-material kernel is the diagonal P0 head, on the same
            bound ends.
        """
        fold = TransferMaterialField(
            per_material={
                mid: TransferKernel(
                    moments=(np.diag(np.diag(k.p0)),),
                    multiplicity=k.multiplicity,
                )
                for mid, k in self.transfer.per_material.items()
            },
            cells_by_material=self.transfer.cells_by_material,
        )
        return self._sibling(fold)

    def residual_part(self) -> Self:
        r"""Return the non-foldable sibling of :math:`T`.

        Carries everything :meth:`foldable_part` does not: the
        off-diagonal of each material's :math:`\Sigma_{c,0}` (cross-group
        P0) and every :math:`P_\ell \ge 1` block verbatim.

        Notes
        -----
        Algebraic contract:
        ``T.apply(ψ) ≈ T.foldable_part().apply(ψ) +
        T.residual_part().apply(ψ)`` at ``rtol=1e-14``.
        """
        residual = TransferMaterialField(
            per_material={
                mid: TransferKernel(
                    moments=(
                        k.p0 - np.diag(np.diag(k.p0)),
                        *k.moments[1:],
                    ),
                    multiplicity=k.multiplicity,
                )
                for mid, k in self.transfer.per_material.items()
            },
            cells_by_material=self.transfer.cells_by_material,
        )
        return self._sibling(residual)

    def is_foldable_into_sigma_r(self) -> bool:
        r"""``True`` iff this operator is structurally the
        :meth:`foldable_part` of some parent :math:`T` — order 0 with
        every material's P0 diagonal.
        """
        if self.scattering_order != 0:
            return False
        return all(
            np.allclose(k.p0, np.diag(np.diag(k.p0)))
            for k in self.transfer.per_material.values()
        )

    def foldable_sigma(self) -> dict[int, np.ndarray]:
        r"""The per-material foldable cross-section
        :math:`(y\,\Sigma_{c,0}^{g\to g})_g` — ``{mid: (ng,)}``, fresh
        copies, read off the bound kernel field (the removal fold
        :math:`\Sigma_r = \Sigma_t - y\,\Sigma_{c,0}^{g\to g}`; the yield
        is 1 for scattering, so the shipped fold is unchanged)."""
        return {
            mid: k.multiplicity * np.diag(k.p0)
            for mid, k in self.transfer.per_material.items()
        }

    # ── LinearOperator surface ─────────────────────────────────────────

    @singledispatchmethod
    def _apply_impl(self, psi) -> "Any":
        r"""Runtime dispatch for :meth:`apply` — see the typed overloads.

        Applies the full emission :math:`T\,\psi`, dispatched on the
        input *carrier* type via :func:`functools.singledispatchmethod`. The
        public :meth:`apply` aliases this dispatcher and carries the per-carrier
        ``@overload`` surface, so callers statically see the output type per
        input carrier:

        * :class:`~orpheus.transport.timed_full_field.TimedFullField` /
          :class:`~orpheus.transport.full_field.FullField`
          → :class:`~orpheus.transport.full_field.FullField` — composite bulk +
          boundary variant. Bulk = the full :math:`P_\ell` path; boundary =
          implicit-zero (a transfer gain is volumetric, ``block_role = BULK``).
        * :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`
          → :class:`~orpheus.transport.source_sinks.ScalarSourceSink` —
          :math:`P_0` only, in **iso scalar magnitude** (no
          :math:`P_\ell`; no :math:`1/W`).
        * :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
          → :class:`~orpheus.transport.source_sinks.AngularSourceSink` —
          full :math:`P_\ell` Galerkin in **per-ordinate magnitude** (the
          trailing :math:`1/W` lives at this producer boundary).
        * :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
          → :class:`~orpheus.transport.source_sinks.AngularSourceSink` — the
          angular-windowing path: :math:`T` consumes flux MOMENTS (already
          :math:`M\psi`, so :math:`M` is skipped), bit-identical to the
          :class:`AngularFlux` arm for :math:`\phi = M\psi`.

        The internal helpers :meth:`add_iso_source` and
        :meth:`build_aniso_source` remain available for callers that need
        the iso / aniso pieces separately.
        """
        raise TypeError(
            f"{type(self).__name__}.apply: unsupported input type "
            f"{type(psi).__name__}; expected TimedFullField, ScalarFlux, "
            f"AngularFlux, or HarmonicMomentFlux.  Dispatch table is "
            f"registered via @singledispatchmethod."
        )

    @_apply_impl.register
    def _(self, psi: FullField) -> "FullField":
        r"""Composite :class:`FullField` variant — bulk-only emission.

        Registered on the timeless :class:`FullField`: a
        :class:`TimedFullField` iterate dispatches here via MRO (it IS a
        ``FullField``), and a bare ``FullField`` dispatches correctly; the arm
        reads only ``psi.interior`` (history-blind). Bulk follows the same math
        as the :class:`AngularFlux` arm; the boundary is the **implicit-zero**
        :class:`~orpheus.transport.source_sinks.AngularBoundarySourceSink` —
        a transfer gain is volumetric (``block_role = BlockRole.BULK``), no
        face-trace contribution, so ``(L + C - S - N₂ₙ - B)`` composes under
        :meth:`TimedFullField.__sub__`.
        """
        # Delegate the bulk source to the bulk-type dispatch arm and wrap with
        # the implicit-zero boundary. ``psi.interior`` is either the full-angular
        # AngularFlux or the windowed HarmonicMomentFlux (both return an
        # AngularSourceSink); the cast names that runtime truth so the typed
        # ``apply`` overloads resolve.
        combined = self.apply(cast("AngularFlux | HarmonicMomentFlux", psi.interior))
        # T is PURE BULK (#282 route (a)): the (ray, bulk) closed-μ emission
        # lives on the SN coupling operator
        # :class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicEmission`,
        # not here — see docs/theory/foundations/coupled_block_operator.rst.
        return FullField(
            interior=combined,
            boundary=AngularBoundarySourceSink.zeros(psi.boundary.space),
        )

    @_apply_impl.register
    def _(self, phi: ScalarFlux) -> "ScalarSourceSink":
        r"""Typed ScalarFlux variant — iso scalar magnitude output (P0 only).

        :math:`Q_g = y\,\Sigma_{c,0}(g'\to g)\,\phi_{g'}`. No :math:`P_\ell`
        (scalar flux lacks angular info); no :math:`1/W` (scalar
        consumers — diffusion / CP / kinetics — do not project to
        per-ordinate).

        **Deliberately retained — a named-future-consumer surface, NOT dead
        weight.** This arm has no current production caller (the within-group
        SI/Krylov path feeds the composite / :class:`AngularFlux` /
        :class:`HarmonicMomentFlux` arms, never a bare :class:`ScalarFlux`, and
        no sibling arm delegates here). It is kept as the typed entry-point for
        the scalar-carrier cross-method consumers the cross-method field
        architecture (`#205
        <https://github.com/deOliveira-R/ORPHEUS/issues/205>`_) will wire;
        the keep-vs-retire call is recorded here rather than left a silent
        orphan.
        """
        # CS4b S4 — the accumulator rides the operand's space (role-blind
        # shared mint; role is class identity).
        iso: ScalarSourceSink = ScalarSourceSink.zeros(phi.space)
        iso = self.add_iso_source(iso, phi)
        return iso

    @_apply_impl.register
    def _(self, psi: AngularFlux) -> "AngularSourceSink":
        r"""Typed :class:`AngularFlux` variant — per-ordinate magnitude output.

        Reduce ``psi`` angular → scalar, build the iso :math:`P_0`
        emission and the per-ordinate :math:`P_\ell\ge 1` Galerkin contribution
        (:meth:`build_aniso_source`), then combine via the producer-side
        :math:`1/W` in :meth:`_assemble_per_ordinate_source`.
        """
        # φ = ∫ψ dΩ (scalar), aniso = (1/W) RΛM ψ (per-ordinate), then the
        # shared producer-side assembly. The iso P0 keeps the cheap
        # reaction-rate fast path (NO moment tensor) — a load-bearing PERF
        # optimisation on the SI-sweep hot path; routing it through the frame
        # regresses LD/P0 badly. The frame form lives on as
        # :attr:`full_transfer_kernel` for the (non-hot-path) adjoint transpose
        # only — see docs/theory/methods/sn/adjoint.rst §sn-scattering-adjoint-source.
        return self._assemble_per_ordinate_source(
            psi.integrate_angular(), self.build_aniso_source(psi), psi.space,
        )

    @_apply_impl.register
    def _(self, phi_moments: HarmonicMomentFlux) -> "AngularSourceSink":
        r"""Windowed moment-iterate variant — :math:`T` consumes flux MOMENTS.

        The angular-windowing path: when the within-group SI iterate is stored
        as harmonic moments :math:`\phi_\ell^m` (the 2-D Cartesian windowed
        iterate) instead of the full per-ordinate :class:`AngularFlux`,
        :math:`T` consumes the moments WITHOUT the :math:`M` projection — the
        moments ARE :math:`M\psi`, so projecting again is redundant.

        Structurally parallel to the :class:`AngularFlux` arm and bit-identical
        to it for :math:`\phi = M\psi`: the :math:`\ell=0` moment IS the scalar
        flux (:math:`Y_0^0 = 1`, so :meth:`HarmonicMomentFlux.scalar_flux`
        equals :meth:`AngularFlux.integrate_angular`), feeding the identical
        P0 fast path; the :math:`\ell\ge 1` aniso takes the explicit
        typed grid path (:math:`\Lambda` then the frame's :math:`R`), equal to
        the full-angular path's :math:`R\Lambda` after its :math:`M`. Both arms
        end at the shared :meth:`_assemble_per_ordinate_source` (per-ordinate
        :class:`AngularSourceSink`, producer-side :math:`1/W`). The
        explicit-typed vs fused-kernel choice is in
        ``docs/theory/foundations/operator_algebra.rst §integral-kernel-category``.
        """
        # F-1: the per-ordinate target is the posed composite's interior —
        # the same space the minted faces are bound to (their codomain /
        # domain read). Reading it refuses loudly on a space-less operator
        # (the A1 refusal, relocated to the mint seam), so no local guard
        # remains here.
        angular_target = self._interior_space
        if self.is_isotropic:
            aniso = None
        else:
            # Explicit typed grid path: Λ maps flux moments → source moments
            # (HarmonicMomentFlux → HarmonicMomentSourceSink, the role-changing
            # edge), the minted source-reconstruction FACE synthesises the
            # per-ordinate AngularSourceSink (riding its bound angular
            # codomain), then the producer-side 1/W. Numerically equals the
            # kernel's ndarray reconstruct_after(Λ) reference.
            emitted = self._moment_transfer(skip_l0=True).apply(
                phi_moments,
            )
            aniso = self.source_reconstruction.apply(
                emitted,
            ) / self.total_weight
        # ℓ=0 moment IS the scalar flux (Y_0^0 = 1) — the typed accessor
        # carries that convention (== integrate_angular bit-exactly).
        assert angular_target.axes is not None  # type-narrowing (axis-built)
        return self._assemble_per_ordinate_source(
            phi_moments.scalar_flux(
                space=FunctionSpace.of_axes(*angular_target.axes[1:]),
            ),
            aniso,
            angular_target,
        )

    if TYPE_CHECKING:
        # Honest per-carrier typing surface (#257 S8c).  ``TransferOperator``
        # is NOT an endomorphism ``V -> V`` (the mixin's nominal contract):
        # it maps each input carrier to a DISTINCT output carrier.  These
        # ``@overload`` stubs exist only for the type checker; the public
        # ``apply`` IS the runtime dispatcher (``apply = _apply_impl`` below),
        # so callers statically see e.g. ``T.apply(ScalarFlux) ->
        # ScalarSourceSink`` instead of the dispatcher's untyped fallback.
        @overload
        def apply(self, psi: FullField, /) -> "FullField": ...
        @overload
        def apply(self, phi: ScalarFlux, /) -> "ScalarSourceSink": ...
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
        r"""The adjoint emission :math:`T^{T}\chi =
        (1/W)\,\mathrm{full\_transfer\_kernel}^{T}\chi` (closes
        `#118 <https://github.com/deOliveira-R/ORPHEUS/issues/118>`_).

        :math:`T^{T}` is the group-and-angle transpose the adjoint transport
        equation :math:`(L+C-S-N_{2n})^{T}\psi^{*}=q^{*}` needs. The production
        FORWARD keeps the scalar fast-path (:attr:`isotropic_energy` +
        :meth:`build_aniso_source`) for SI-sweep performance, so the adjoint —
        NOT the hot path — rides the validated harmonic-frame
        :attr:`full_transfer_kernel`, whose transpose falls out of
        :meth:`~orpheus.numerics.operator.OperatorProduct.apply_transpose` in
        ONE expression (iso :math:`\ell=0` + aniso :math:`\ell\ge1`).

        This is the **Euclidean** transpose (L12) — NOT the metric Hilbert
        adjoint ``.H`` (which would carry the angular Gram). The
        forward/adjoint structural asymmetry (which makes the reciprocity
        :math:`\langle T\psi,\chi\rangle=\langle\psi,T^{T}\chi\rangle` a genuine
        cross-check, not a tautology) and the proof are in
        ``docs/theory/methods/sn/adjoint.rst §sn-scattering-adjoint-source``.

        The COMPOSITE (``FullField``) arm mirrors the forward lift: a transfer
        gain is volumetric, so :math:`T^{T}` reads ONLY the bulk cotangent
        (``chi.interior``) and emits the implicit-zero trace, letting the full
        within-group loss ``(L + C - S - N₂ₙ - B).H`` compose through
        ``OperatorSum.apply_transpose``.

        Parameters
        ----------
        chi : np.ndarray, carrier, or FullField
            The daggered per-ordinate field :math:`\chi` of shape
            ``(N, ng, *spatial)`` — a typed carrier's ``.values`` are unwrapped
            and a bare ``ndarray`` returns; or the composite ``FullField``,
            returning a composite with bulk-only content.
        """
        if isinstance(chi, FullField):
            bulk_bar = self.apply_transpose(np.asarray(chi.interior.values))
            # T is PURE BULK, so T^T is too (#282 route (a)): the seed-cotangent
            # pullback lives on RadialCharacteristicEmission.apply_transpose —
            # see docs/theory/foundations/coupled_block_operator.rst.
            # CS4b S4 — the space route: output blocks ride the operand's.
            return FullField(
                interior=AngularSourceSink(
                    values=bulk_bar, space=chi.interior.space,
                ),
                boundary=AngularBoundarySourceSink.zeros(chi.boundary.space),
            )
        chi_values = np.asarray(getattr(chi, "values", chi))
        return (
            self.full_transfer_kernel.apply_transpose(chi_values)
            / self.total_weight
        )
