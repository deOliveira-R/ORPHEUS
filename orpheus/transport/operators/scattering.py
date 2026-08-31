r"""Multigroup scattering source operator :math:`S` as a :class:`LinearOperator`.

This module owns the **scattering gain** :math:`S` of the honest
within-group operator algebra :math:`A = L + C - S - B` (streaming
:math:`L`, collision :math:`C`, scattering gain :math:`S`, boundary
:math:`B`; the eigenvalue problem is :math:`K = A^{-1}F`). The operator
letters are pinned in ``docs/theory/conventions/notation.rst
§notation-symbol-table``; the multigroup posing in
``docs/theory/methods/sn/slab_multigroup.rst
§sn-scattering-fission-operators``.

:math:`S` aggregates **all secondary-emission channels that depend on
the in-cell flux**:

* **P0 isotropic in-scatter** :math:`\Sigma_s^0(g'\to g)\,\phi_{g'}`;
* **Pℓ Galerkin reconstruction** on real spherical harmonics
  :math:`Y_\ell^m` (:math:`\ell\ge 1`) from the angular-flux moments;
* **(n,2n) doubling** :math:`2\,\Sigma_{2n}(g'\to g)\,\phi_{g'}`.

The theory lives in the book — one concept, one home:

* P0 in-scatter, the (n,2n) fold, the Pℓ reconstruction, and the
  :math:`1/W` normalisation chain —
  ``docs/theory/methods/sn/slab_multigroup.rst §mg-inscatter-source``,
  §pn-scatter, §n2n-source, §pn-scatter-rlm.
* The no-prefactor :math:`Y_\ell^m` convention and the Funk–Hecke
  eigenbasis (why :math:`S = R\circ\Lambda\circ M`) —
  ``docs/theory/foundations/spherical_harmonics.rst
  §spherical-harmonics-eigenbasis``.
* The §5.6 integral-kernel reading and the apply-only capability
  surface — ``docs/theory/foundations/operator_algebra.rst
  §integral-kernel-category``, §capability-set-semantics.
* The Euclidean adjoint :math:`S^{T}` (forward fast-path vs adjoint
  frame-form) — ``docs/theory/methods/sn/adjoint.rst
  §sn-scattering-adjoint-source``.

Capability surface
==================

``apply`` + ``apply_transpose`` (``is_adjointable=True``); **no**
``solve`` (``is_invertible=False``). :math:`S` is rank
:math:`O(N_{\text{cells}}\cdot N_{\text{groups}})` with no tractable
inverse — it is *applied*, never *inverted*, and the missing ``solve``
is structural method-absence (a composer refuses to build :math:`S^{-1}`
at construction time), not an advertising flag. The adjoint :math:`S^{T}`
rides the harmonic-frame :attr:`~ScatteringOperator.full_scatter_kernel`
(closes `#118 <https://github.com/deOliveira-R/ORPHEUS/issues/118>`_):
see :meth:`~ScatteringOperator.apply_transpose`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property, singledispatchmethod
from typing import TYPE_CHECKING, Any, cast, overload

import numpy as np

from orpheus.numerics.frame import GalerkinFrame
from orpheus.numerics.space import FunctionSpace
from orpheus.numerics.spaces import SphericalHarmonicSpace
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
    OperatorSum,
)

# Runtime imports of the flux / source types — required at module load
# time because :func:`singledispatchmethod.register` dispatches on the
# runtime class.  These three modules form a leaf in the SN package
# dependency graph (they do not import scattering.py), so the imports
# are circular-import-safe.
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.source_sinks import (
    ScalarSourceSink,
    AngularSourceSink,
    AngularBoundarySourceSink,
    HarmonicMomentSourceSink,
)
from orpheus.transport.full_field import FullField
from orpheus.transport.kernels import ScatteringKernel
from orpheus.transport.material_field import (
    N2NMaterialField,
    ScatteringMaterialField,
)
from orpheus.transport.operators._per_ordinate import (
    assemble_per_ordinate_isotropic,
)
from orpheus.transport.operators.bound_operator import BoundOperator
from orpheus.transport.timed_full_field import TimedFullField
# Runtime (not TYPE_CHECKING-only) — the windowed moment-iterate ``apply``
# arm registers on this type via ``@apply.register``, which needs the
# concrete class at class-definition time.  Circular-import-safe (the
# transport field types do not import scattering.py; same as ScalarFlux /
# AngularFlux above).
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux

if TYPE_CHECKING:
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.transport.mesh.material_xs_field import MaterialXSField
    from orpheus.transport.operators.isotropic_scattering import (
        IsotropicScattering,
    )


__all__ = ["LegendreMomentScattering", "ScatteringOperator"]


@dataclass(eq=False)
class LegendreMomentScattering(BoundOperator):
    r"""Per-ℓ block-diagonal scattering :math:`\Lambda` on harmonic-moment space.

    The diagonal spectrum of the sum-of-tensor-products form
    :math:`\Lambda = \sum_{\ell=1}^{L} \mathbf{P}_\ell \otimes
    \Sigma_{s,\ell}` (:math:`\mathbf{P}_\ell` selects the :math:`\ell`-th
    harmonic block, :math:`\Sigma_{s,\ell}` the per-material per-ℓ Legendre
    matmul on the group axis) — see
    ``docs/theory/foundations/operator_algebra.rst
    §scattering-as-tensor-product-sum``.

    For an input moment field
    :math:`\phi_\ell^m(\vec r) \in \mathbb{R}^{(L+1) \times (2L+1) \times N_g \times *\text{spatial}}`
    (``g`` leading the spatial axes), the action is per-group-transfer
    per :math:`\ell` block,

    .. math::

        (\Lambda \phi)_\ell^m(\vec r)\bigg|_{g}
        \;=\; \sum_{g'} \Sigma_{s,\ell}^{m(\vec r)}(g' \to g)\,
              \phi_\ell^m(\vec r)\bigg|_{g'},

    with :math:`m(\vec r)` the material id at cell :math:`\vec r` (the
    per-material structure IS the bound datum — see below).

    **The CS4c rebind (design record §14):** the operator holds the
    representation-free datum — a
    :class:`~orpheus.transport.material_field.ScatteringMaterialField`
    already :meth:`truncated
    <orpheus.transport.material_field.ScatteringMaterialField.truncated>`
    to the binding's order — plus its two mandatory ends (the
    :class:`~orpheus.transport.operators.bound_operator.BoundOperator`
    base; :math:`\Lambda` is endomorphic on the SH coefficient space of
    its order). ``L`` is DERIVED from the field (the truncation IS the
    order — single source); the per-material dispatch lives on the
    field's :meth:`~orpheus.transport.material_field.ScatteringMaterialField.moment_source`
    verb, whose shape guard refuses a moments tensor at any other order.

    ``skip_l0`` (default ``True``) skips the :math:`\ell = 0` block, which
    the project's P0 in-scatter handles on a separate reaction-rate fast
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
    scattering : ScatteringMaterialField
        The per-material Legendre transfer stacks over the mesh layout,
        truncated to this binding's order.
    skip_l0 : bool, default ``True``
        Skip the :math:`\ell = 0` block (handled by P0 in-scatter). Set
        ``False`` for the full :math:`R \Lambda M \psi` composition.
    domain, codomain : FunctionSpace
        The two mandatory ends (kw-only, write-once — the base) — the SH
        coefficient space of the field's order, both, for the shipped
        endomorphic binding.
    """

    scattering: "ScatteringMaterialField"
    skip_l0: bool = True

    @classmethod
    def from_material_xs(
        cls,
        mat_xs: "MaterialXSField",
        L: int,
        *,
        skip_l0: bool = True,
    ) -> "LegendreMomentScattering":
        r"""Tier-2 extract-and-mint: pull the scattering channel of a
        :class:`~orpheus.transport.mesh.material_xs_field.MaterialXSField`
        facade, truncate to ``L``, and bind the endomorphic SH ends
        (``SphericalHarmonicSpace.from_L(L)`` supplying both — the
        endomorphism sugar lives HERE, never on the exact ctor)."""
        sh_space = SphericalHarmonicSpace.from_L(L)
        return cls(
            ScatteringMaterialField.from_material_xs(mat_xs).truncated(L),
            skip_l0=skip_l0,
            domain=sh_space,
            codomain=sh_space,
        )

    @property
    def L(self) -> int:
        r"""The Legendre truncation :math:`L` — DERIVED from the bound
        field (the truncation is the order; single source)."""
        return self.scattering.order

    @property
    def is_adjointable(self) -> bool:
        # Λ exposes its group-transpose Σ_{s,ℓ}^T (apply_transpose), so the
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

        :math:`\Lambda` maps a flux moment to the in-scatter **source** moment
        it emits (flux → source); the typed arm makes that role change explicit
        in the signature. (Why the role change lives on the operator and not
        on the frame: ``docs/theory/foundations/operator_algebra.rst
        §integral-kernel-category``.)

        Parameters
        ----------
        moments : np.ndarray or HarmonicMomentFlux
            Flux moment field of shape ``(L+1, 2L+1, ng, *spatial)``.  The
            :math:`m`-axis is the addition-theorem-shifted index where slot
            ``l + m`` holds the :math:`(\ell, m)` entry; entries outside
            :math:`|m| \le \ell` are conventionally zero.

            Typed
            :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
            → typed
            :class:`~orpheus.transport.source_sinks.harmonic_moment_source_sink.HarmonicMomentSourceSink`
            (the scattered moment SOURCE) with matching ``L`` / ``mesh`` /
            spatial-moment width.  Bare ndarray → ndarray (the endomorphic
            moment-space view the :math:`R\circ\Lambda\circ M`
            kernel ``OperatorProduct`` composes on).

        Returns
        -------
        np.ndarray or HarmonicMomentSourceSink
            The scattered moment source, same shape as ``moments``.  The
            :math:`\ell = 0` block is zero when ``skip_l0`` is ``True``;
            otherwise the P0 in-scatter contribution is included.  Typed in →
            typed source out; ndarray in → ndarray out.

        Notes
        -----
        Both arms route through the field's single per-material verb
        :meth:`~orpheus.transport.material_field.ScatteringMaterialField.moment_source`;
        they differ only in the carrier wrap.
        """
        from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
        if isinstance(moments, HarmonicMomentFlux):
            out_values = self.scattering.moment_source(
                moments.values, skip_l0=self.skip_l0,
            )
            # flux moment → source moment: the explicit role change
            # (CS4b S4 — same space, new class; role is class identity).
            return HarmonicMomentSourceSink(
                values=out_values, space=moments.space, L=moments.L,
                spatial_moments=moments.spatial_moments,
            )
        return self.scattering.moment_source(moments, skip_l0=self.skip_l0)

    def apply_transpose(
        self, moments: "np.ndarray | HarmonicMomentSourceSink",
    ) -> "np.ndarray | HarmonicMomentFlux":
        r"""Apply :math:`\Lambda^{T}` — the per-ℓ group-transpose (the role-REVERSING edge).

        The bare Euclidean transpose of :meth:`apply`: :math:`\Lambda^{T}` maps a
        source moment back into the flux-moment space it scattered from (source →
        flux), transposing the per-ℓ group-transfer
        :math:`\Sigma_{s,\ell}(g'\to g) \mapsto (g\to g')`.  Routes through the
        field's transpose verb
        :meth:`~orpheus.transport.material_field.ScatteringMaterialField.moment_source_transpose`.

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
        from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux

        if isinstance(moments, HarmonicMomentSourceSink):
            out_values = self.scattering.moment_source_transpose(
                moments.values, skip_l0=self.skip_l0,
            )
            return HarmonicMomentFlux(
                values=out_values, space=moments.space, L=moments.L,
                spatial_moments=moments.spatial_moments,
            )
        return self.scattering.moment_source_transpose(
            moments, skip_l0=self.skip_l0,
        )


@dataclass(eq=False)
class N2NMomentOperator(BoundOperator):
    r"""The (n,2n) isotropic :math:`\ell=0` moment operator :math:`\nu_{2n}\,\Sigma_{2n}`.

    The (n,2n) reaction is a DISTINCT isotropic (:math:`\ell=0`) group transfer —
    a *multiplication* channel (each event emits two neutrons), NOT scattering —
    so it is its own named operator, summed with :math:`\Lambda` in moment
    space (an :class:`~orpheus.numerics.operator.OperatorSum`) where a full
    in-scatter conjugation wants both. Keeping the
    multiplication reaction a visible distinct operator, rather than hidden in
    the scattering matmul, is the physics-faithful choice — see
    ``docs/theory/methods/sn/adjoint.rst §sn-scattering-adjoint-source``.

    **The CS4c rebind (design record §14):** the bound datum is a
    :class:`~orpheus.transport.material_field.N2NMaterialField` — the raw
    per-material reaction matrices over the mesh layout, with the
    multiplicity :math:`\nu_{2n}` entering from its one home
    (:attr:`~orpheus.transport.kernels.N2NKernel.multiplicity`) inside the
    field's verbs, never as a literal here. The two mandatory ends are the
    SH coefficient space it composes on (it reads/writes ONLY the
    :math:`\ell=0` block, so it joins an ``OperatorSum`` with
    :math:`\Lambda` on the same space).

    Parameters
    ----------
    n2n : N2NMaterialField
        The per-material :math:`(n,2n)` reaction matrices over the layout.
    domain, codomain : FunctionSpace
        The two mandatory ends (kw-only, write-once — the base): the SH
        coefficient space of the composing order, both.
    """

    n2n: "N2NMaterialField"

    @classmethod
    def from_material_xs(
        cls, mat_xs: "MaterialXSField", L: int,
    ) -> "N2NMomentOperator":
        r"""Tier-2 extract-and-mint: pull the :math:`(n,2n)` channel of a
        facade and bind the endomorphic SH ends at order ``L``."""
        sh_space = SphericalHarmonicSpace.from_L(L)
        return cls(
            N2NMaterialField.from_material_xs(mat_xs),
            domain=sh_space,
            codomain=sh_space,
        )

    @property
    def is_adjointable(self) -> bool:
        # ν·Σ_{2n}^T (apply_transpose) is the ℓ=0 group-transpose; caps ⊇
        # apply_transpose. is_invertible inherits base False.
        return True

    def apply(self, moments: np.ndarray) -> np.ndarray:
        r""":math:`\nu_{2n}\,\Sigma_{2n}` applied to the :math:`\ell=0` moment (ℓ≥1 zero).

        Bare ndarray (the endomorphic moment-space view a ``frame.conjugate``
        ``OperatorProduct`` chain composes on).
        """
        return self.n2n.moment_emission(moments)

    def apply_transpose(self, moments: np.ndarray) -> np.ndarray:
        r"""The :math:`\ell=0` group-transpose :math:`(\nu_{2n}\,\Sigma_{2n})^{T}` (bare ndarray)."""
        return self.n2n.moment_emission_transpose(moments)


@dataclass(eq=False)
class ScatteringOperator(BoundOperator["FullField"]):
    r"""Scattering source operator :math:`S` (P0 + Pℓ).

    **The CS4c step-3 binding (design record §14):** the exact ctor
    retains the representation-free datum and the minted products, and
    nothing richer —

    * :attr:`scattering` — the
      :class:`~orpheus.transport.material_field.ScatteringMaterialField`
      (per-material Legendre stacks over the mesh layout), already
      truncated to this binding's order (the truncation IS
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
      ``(L + C) − S − N₂ₙ − B`` OperatorSum guard validates the ``− S``
      arm natively. (⚠ the 2-D windowed operand carries a MOMENT
      interior — the shipped non-endomorphism the step-0 census
      measured; the carrier dispatch serves it until step 5's arm
      deletion.)

    **The (n,2n) channel is NOT here** (§14.1, ruled 2026-08-30): it is
    the first-class
    :class:`~orpheus.transport.operators.n2n.N2NOperator`, and the
    within-group algebra spells ``− S − N₂ₙ`` explicitly. Use
    :meth:`from_solver_data` to build instances from a
    :class:`MaterialXSField` + the posed composite space (the quadrature
    is reached through the space's angular axis — the CS5 generator
    channel; no ``quadrature`` field survives).

    Capability surface: ``{apply, apply_transpose}`` — no efficient
    ``solve``; the adjoint :math:`S^{T}` is free via the harmonic-frame
    :attr:`full_scatter_kernel` (see :meth:`apply_transpose`).
    """

    scattering: "ScatteringMaterialField"
    #: The minted FLUX analysis face :math:`M \otimes I` on the posed
    #: interior (``AngularFlux → HarmonicMomentFlux``) — bound at tier 2,
    #: retained. Consumers: the windowed bulk projection and the S6
    #: adjoint gates.
    flux_analysis: "HarmonicAnalysisOperator[AngularFlux, HarmonicMomentFlux]" = field(
        kw_only=True,
    )
    #: The minted SOURCE reconstruction face :math:`R \otimes I` landing
    #: on the posed interior (``HarmonicMomentSourceSink →
    #: AngularSourceSink``) — the windowed in-scatter arm's typed R.
    source_reconstruction: "HarmonicReconstructionOperator[HarmonicMomentSourceSink, AngularSourceSink]" = field(
        kw_only=True,
    )

    # Scattering is a BULK operator — the moment-folding `Σ_s · ⟨P_ℓ, ψ⟩`
    # reads and writes the bulk flux only (A_bb), no boundary action.
    # Class-level constant (unannotated so the dataclass does not treat it
    # as a field).
    block_role = BlockRole.BULK

    def __post_init__(self) -> None:
        # Per-END energy admission (CS4c §1; NEW wiring — S carried no
        # energy-conformity guard before this flip: the reach census in
        # tests/transport/test_kernels.py G2.10 row 3 re-derives).
        self._assert_energy_extent_both_ends(
            self.scattering.ng, operator="ScatteringOperator",
        )
        # Face-binding agreement: both minted faces must be bound to THIS
        # binding's interior — a face bound elsewhere would make the
        # windowed arm and the per-ordinate combine silently inconsistent.
        # (The XD-1 wrong-EMBEDDING controls stay spellable: a doctored
        # face carries the RIGHT spaces with the WRONG measure.)
        interior = self._interior_space
        if self.flux_analysis.domain != interior:
            raise TypeError(
                "ScatteringOperator: the flux-analysis face is bound to "
                "a different angular space than this binding's interior "
                "— mint the faces from the SAME posed space (tier 2 does)."
            )

    @property
    def is_adjointable(self) -> bool:
        # S = R∘Λ∘M exposes its Euclidean transpose S^T via
        # :attr:`full_scatter_kernel` (the OperatorProduct adjoint chain);
        # is_invertible inherits base False —
        # a scattering source operator is not invertible.
        return True

    @property
    def scattering_order(self) -> int:
        r"""Maximum Legendre order :math:`L` retained — DERIVED from the
        bound field (the truncation IS the order; single source). ``0``
        means P0 only."""
        return self.scattering.order

    # (The legacy n_ordinates / weights / Y / ng / spatial_shape /
    # sig_s / sig2 / sig_s0 / cells_by_mat read-throughs retired with
    # the CS4c step-3 rebind — closing #306. The data lives on
    # :attr:`scattering`; the harmonics table on the faces' shared
    # frame; the (n,2n) channel on N2NOperator.)

    @property
    def total_weight(self) -> float:
        r""":math:`W = \int_{S^2} d\Omega` — the binding measure's total
        angular weight (the producer-side :math:`/W`). Read off the
        retained faces' frame MEASURE: the measure is the binding's
        metric (operative data — unlike the :attr:`frame` accessor,
        which is provenance)."""
        return float(self.flux_analysis.frame.measure.weights.sum())

    @property
    def _sh_space(self) -> "FunctionSpace":
        r"""The SH coefficient space of this binding's order — the
        endomorphic ends of the internally-minted moment factors."""
        return SphericalHarmonicSpace.from_L(self.scattering_order)

    def _moment_scattering(self, *, skip_l0: bool) -> LegendreMomentScattering:
        r"""Mint the moment-space :math:`\Lambda` factor on this binding's
        datum + SH ends — the ONE internal spelling (three consumers:
        the §5.6 kernel, the full conjugation, and the aniso moment
        route; the windowed arm consumes the cached :attr:`kernel`
        factors)."""
        sh = self._sh_space
        return LegendreMomentScattering(
            self.scattering, skip_l0=skip_l0, domain=sh, codomain=sh,
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
                "ScatteringOperator binds the composite full-field space"
                " — this instance's domain carries no interior to size"
                " the angular arms."
            )
        return domain.interior_space

    @cached_property
    def isotropic_energy(self) -> "IsotropicScattering":
        r"""The P0 ENERGY binding of this operator's own datum — the
        scalar-space :class:`~orpheus.transport.operators.isotropic_scattering.IsotropicScattering`
        the per-ordinate fast path lifts (and the solver's K_iso
        assembly consumes: ``K_iso = S.isotropic_energy + N2N.energy``,
        composed at the ONE within-group construction site since the
        §14.1 extraction retired ``S.isotropic_kernel``).

        Bound to the composite interior's scalar sub-space; shares the
        p0 arrays with :attr:`scattering` (truncation shares storage).
        Cached at construction-time semantics (the kernel field is
        immutable — snapshot identical to the retired facade caches).
        """
        from orpheus.transport.operators.isotropic_scattering import (
            IsotropicScattering,
        )

        interior = self._interior_space
        if interior.axes is None:
            raise TypeError(
                "ScatteringOperator.isotropic_energy: the composite "
                "interior must be axis-built to name the scalar "
                "sub-space."
            )
        scalar_space = FunctionSpace.of_axes(*interior.axes[1:])
        return IsotropicScattering(
            self.scattering.truncated(0),
            domain=scalar_space,
            codomain=scalar_space,
        )

    # ── Anisotropic in-scatter: the moment→source map R·Λ_{ℓ≥1} ────────

    def _aniso_source_from_moment_values(
        self, moment_values: np.ndarray,
    ) -> np.ndarray:
        r"""Reconstruct the per-ordinate :math:`\ell\ge 1` in-scatter source
        from a flux-moment tensor — the **moment → source** half
        :math:`R\,\Lambda_{\ell\ge 1}` of the anisotropic composition
        :math:`S_{\rm aniso} = \tfrac{1}{W}\,R\,\Lambda\,M`.

        Built as ``frame.reconstruct_after(Λ)`` — the §5.6 :attr:`kernel`
        sub-operator for inputs ALREADY in moment space, with
        :math:`\Lambda_{\ell\ge 1}` = :class:`LegendreMomentScattering`
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
            Flux-moment tensor ``(L+1, 2L+1, ng, nx, ny)`` — typically
            ``M.apply(psi)`` or the windowed iterate's
            :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
            ``.values``.

        Returns
        -------
        np.ndarray
            ``(N, ng, nx, ny)`` per-ordinate :math:`\ell\ge 1` in-scatter
            source, **pre** :math:`1/W`.
        """
        scatter = self._moment_scattering(skip_l0=True)
        # R∘Λ as ONE typed operator (M already applied — the moments ARE M·ψ).
        return self.frame.reconstruct_after(scatter).apply(moment_values)

    @cached_property
    def kernel(self) -> LinearOperator:
        r"""The §5.6 integral kernel :math:`R \circ \Lambda_{\ell\ge 1} \circ M`.

        The anisotropic Legendre redistribution as a typed
        :class:`~orpheus.numerics.operator.OperatorProduct`
        ``frame.conjugate(Λ)`` ``= OperatorProduct(R, OperatorProduct(Λ, M))``,
        so ``kernel.apply(ψ.values) = R(Λ(M ψ))``. The factors:

        * :math:`M` = ``frame.analysis`` (the :attr:`frame`'s
          :attr:`~orpheus.numerics.frame.FrameBase.analysis` face);
        * :math:`\Lambda_{\ell\ge 1}` = :class:`LegendreMomentScattering`
          ``(mat_xs, L=scattering_order, skip_l0=True)``;
        * :math:`R` = ``frame.reconstruction`` (the :attr:`frame`'s
          :attr:`~orpheus.numerics.frame.FrameBase.reconstruction` face).

        This is the production :math:`\ell\ge 1` map (not a parallel semantic
        reading): :meth:`build_aniso_source` is ``(1/W)·kernel`` and the windowed
        arm consumes the sub-operator ``frame.reconstruct_after(Λ)``. The
        producer-side :math:`1/W` lives OUTSIDE the kernel (at the :meth:`apply`
        boundary), so ``kernel.apply`` returns the source **pre**-:math:`1/W`.
        With it, :class:`ScatteringOperator` satisfies the
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
            in-scatter is the LOCAL component handled by :meth:`add_iso_source`.

        Returns
        -------
        LinearOperator
            The typed ``R∘Λ∘M`` anisotropic redistribution.
        """
        if self.scattering_order == 0:
            raise ValueError(
                "ScatteringOperator.kernel requires scattering_order >= 1: "
                "an isotropic (P0-only) operator has no anisotropic integral "
                "kernel (the R∘Λ∘M angular redistribution is empty). The P0 "
                "in-scatter is the LOCAL component, handled by add_iso_source."
            )
        # R ∘ Λ ∘ M as ONE typed operator: the frame conjugates Λ (per-ℓ
        # moment-space scattering) between its analysis (M) and reconstruction (R)
        # faces. Λ now carries real spaces (== frame.basis_space), so the
        # OperatorProduct composability guard validates R∘Λ∘M natively — NO cast.
        return self.frame.conjugate(
            self._moment_scattering(skip_l0=True)
        )

    @cached_property
    def full_scatter_kernel(self) -> OperatorProduct:
        r"""The FULL in-scatter source kernel :math:`R\circ(\Lambda_{\ell\ge 0} + N_{2n})\circ M`.

        The COMPLETE P0 + anisotropic + (n,2n) in-scatter as ONE frame-conjugated
        operator: the isotropic ℓ=0 scattering, the anisotropic ℓ≥1
        redistribution (one :class:`LegendreMomentScattering`,
        ``skip_l0=False``) conjugated by the frame. The per-ordinate
        source is ``(1/W)·full_scatter_kernel.apply(ψ)``; its transpose
        ``(1/W)·full_scatter_kernel.apply_transpose(ψ*)`` is the adjoint
        :math:`S^{T}` (:meth:`apply_transpose`). Riding the same frame
        conjugation for iso and aniso is what lets the whole transpose
        fall out for free —
        ``docs/theory/methods/sn/adjoint.rst §sn-scattering-adjoint-source``.
        (The :math:`(n,2n)` term left this composite with the §14.1
        extraction — :class:`~orpheus.transport.operators.n2n.N2NOperator`
        carries its own lift and transpose.)

        Distinct from :attr:`kernel` (the §5.6 ℓ≥1 ANISOTROPIC subcomponent + the
        0-ULP crosscheck oracle): this is the FULL ℓ≥0 scattering source.
        CACHED at first access (CS4c §14.7 — the satellite mint drops
        from once-per-apply to once-per-construction; the kernel field is
        immutable, so the cache cannot go stale).
        """
        return self.frame.conjugate(
            self._moment_scattering(skip_l0=False)
        )

    # (``isotropic_kernel`` — the IsoS + IsoN2N sum — retired with the
    # §14.1 extraction: the K_iso composition is the SOLVER's grouping
    # now, assembled at the one within-group construction site from
    # ``S.isotropic_energy + N2N.energy``.)

    def _assemble_per_ordinate_source(
        self,
        phi: "ScalarFlux",
        aniso_or_none: "AngularSourceSink | None",
        angular_space: "FunctionSpace",
    ) -> "AngularSourceSink":
        r"""Combine the P0 + (n,2n) iso source (from the scalar flux
        :math:`\phi_0`) with the pre-:math:`/W` :math:`\ell\ge 1` aniso source
        into the per-ordinate scattering source :math:`(\text{iso}/W) +
        \text{aniso}`.

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
            Scalar flux :math:`\phi_0` (iso magnitude) driving P0 + (n,2n).
        aniso_or_none : AngularSourceSink or None
            Per-ordinate :math:`\ell\ge 1` source ALREADY in per-ordinate
            magnitude (post-:math:`/W`), or ``None`` for
            ``scattering_order == 0``.
        angular_space : FunctionSpace
            The per-ordinate target space (CS4b S4 — sizes the zero
            accumulator; the caller holding the pose supplies it: the
            operand's own space on the full-angular arm, the posed
            composite's interior on the windowed moment arm).
        """
        # The isotropic P0 energy source through this binding's own
        # scalar-space energy operator (the (n,2n) term left with §14.1);
        # the iso rides the driving flux's own space (width travels in
        # it), and the producer-side /W combine is the ONE shared home
        # (:func:`~orpheus.transport.operators._per_ordinate.assemble_per_ordinate_isotropic`
        # — N2NOperator rides the same primitive).
        iso: ScalarSourceSink = ScalarSourceSink(
            values=cast(np.ndarray, self.isotropic_energy.apply(phi)),
            space=phi.space,
        )
        return assemble_per_ordinate_isotropic(
            iso, aniso_or_none, angular_space, self.total_weight,
        )

    @classmethod
    def from_solver_data(
        cls,
        *,
        mat_xs: "MaterialXSField",
        scattering_order: int,
        space: "FunctionSpace",
    ) -> "ScatteringOperator":
        r"""Tier-2 extract-and-mint (CS4c §14): extract the scattering
        channel of the facade, truncate to the binding's order, mint the
        two faces from the HUB-interned frame, and bind the endomorphic
        composite ends from one ``space=``.

        ``space`` is the composite
        :attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.full_field_space`
        the solver threads — MANDATORY since the flip (the ends are
        write-once fields; the OperatorSum guard validates every build).
        The quadrature is reached through the space's angular axis (the
        CS5 generator channel inside
        :meth:`HarmonicFrame.for_space
        <orpheus.transport.frames.harmonic_frame.HarmonicFrame.for_space>`)
        — no ``quadrature=`` parameter survives, so a frame/space metric
        mismatch is unspellable.
        """
        from orpheus.numerics.spaces.full_field_space import FullFieldSpace

        interior = (
            space.interior_space
            if isinstance(space, FullFieldSpace)
            else None
        )
        if interior is None:
            raise TypeError(
                "ScatteringOperator.from_solver_data requires the posed "
                "composite FullFieldSpace (its interior is the angular "
                "space the faces bind); got "
                f"{type(space).__name__}."
            )
        frame = HarmonicFrame.for_space(interior, scattering_order)
        return cls(
            ScatteringMaterialField.from_material_xs(mat_xs).truncated(
                scattering_order,
            ),
            flux_analysis=frame.flux_analysis_on(interior),
            source_reconstruction=frame.source_reconstruction_on(interior),
            domain=space,
            codomain=space,
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
        r"""Add P0 in-scatter source to :math:`Q`.

        Vectorised by material: per cell ``c`` of material ``mid``,
        ``Q[:, ic, jc] += sig_s0[mid].T @ phi[:, ic, jc]`` where
        ``sig_s0[mid]`` is ``(ng, ng)`` indexed ``[g_from, g_to]``.

        Typed-action overload:

        * **Raw-in (legacy)** — ``Q: np.ndarray`` is mutated in place and
          ``None`` is returned.
        * **Typed-in (return-new)** — a frozen ``Q: ScalarSourceSink`` stays
          immutable; a NEW :class:`ScalarSourceSink` carrying
          ``Q.values + in_scatter`` is returned (the caller spells the algebra
          ``Q = scattering.add_iso_source(Q, phi)``).

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
        :meth:`~orpheus.transport.material_field.ScatteringMaterialField.add_p0_source`.
        """
        from orpheus.transport.fields.scalar_flux import ScalarFlux
        from orpheus.transport.source_sinks import ScalarSourceSink
        phi_values = phi.values if isinstance(phi, ScalarFlux) else phi
        if isinstance(Q, ScalarSourceSink):
            Q_values = Q.values.copy()
            self.scattering.add_p0_source(Q_values, phi_values)
            # Preserve Q's spatial-moment width (#240 D5b-S3 — the φ̂ accumulator
            # carries the trailing 2^d axis; re-wrapping without it would lose
            # the typed factor and raise on the shape).
            return ScalarSourceSink(values=Q_values, space=Q.space)
        self.scattering.add_p0_source(Q, phi_values)
        return None

    def build_aniso_source(
        self,
        angular_flux: "AngularFlux | None",
    ) -> "AngularSourceSink | None":
        r"""Build per-ordinate Pℓ (:math:`\ell \ge 1`) scattering source.

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
            Returns ``None`` when ``scattering_order == 0`` or
            ``angular_flux is None``.
        """
        if self.scattering_order == 0 or angular_flux is None:
            return None
        # S_aniso = (1/W)·kernel: the §5.6 :attr:`kernel` is the R∘Λ∘M
        # redistribution; the producer-side /W lives OUTSIDE it (applied here).
        out_values = self.kernel.apply(angular_flux.values) / self.total_weight
        # RΛM is spatial-moment-axis-agnostic (#240 D5b-S3): the source
        # rides the iterate's own space (CS4b S4 — same space, new role),
        # so the moment factor travels with it.
        return AngularSourceSink(values=out_values, space=angular_flux.space)

    # ── Foldable / residual split ─────────────────────────────────────
    #
    # S = S_foldable + S_residual, additive at rtol=1e-14: S_foldable is the
    # P0 within-group self-scatter (diagonal Σ_s0^{g→g}, foldable into the
    # removal cross-section σ_r = σ_t − Σ_s0^{g→g}); S_residual carries
    # everything else (cross-group P0, all Pℓ≥1, (n,2n)). Data API only — no
    # solver/sweep/iteration consumes these methods yet; the intended consumer
    # is a consistent DSA preconditioner (#2). Theory (why each residual piece
    # is unfoldable): docs/theory/methods/sn/loss_representation.rst
    # §loss-rep-removal-sigma.
    #
    # ⚠ LATENT CORRECTNESS TRAP (#215): do NOT wire the σ_r-SWEEP as the
    # within-group A_wg.inverse(). The σ_r-sweep inverts a DIAGONAL-in-angle
    # removal, but S_foldable = Σ_s0·P_iso is the ISOTROPIC-PROJECTION
    # self-scatter — the two coincide ONLY for isotropic flux, so the wiring
    # ships 46–56 % silent flux errors on anisotropic problems. Use consistent
    # DSA (#2) or Krylov (splitting-invariant, already production). Any
    # within-group accelerator MUST be gated on an ANISOTROPIC config — the
    # isotropic box cannot see this error. Full failure table:
    # docs/theory/methods/sn/slab_one_group.rst §si-sigma-r-fold-mismatch.

    def _sibling(
        self, field_: "ScatteringMaterialField",
    ) -> "ScatteringOperator":
        r"""A sibling binding of ``field_`` on the SAME ends — faces
        re-minted from the HUB at the sibling field's own order (the
        interned frame chain, so an order-0 sibling gets order-0 faces).
        """
        interior = self._interior_space
        frame = HarmonicFrame.for_space(interior, field_.order)
        return ScatteringOperator(
            field_,
            flux_analysis=frame.flux_analysis_on(interior),
            source_reconstruction=frame.source_reconstruction_on(interior),
            domain=self.domain,
            codomain=self.codomain,
        )

    def foldable_part(self) -> "ScatteringOperator":
        r"""Return the P0 within-group self-scatter sibling of :math:`S`.

        Carries only the diagonal of each material's :math:`\Sigma_{s,0}`
        — the within-group self-scatter cross-section
        :math:`\Sigma_{s,0}^{g\to g}` per energy group. All other
        scattering channels (cross-group P0, every :math:`P_\ell \ge 1`)
        live in :meth:`residual_part`. (The :math:`(n,2n)` channel is not
        this operator's since the §14.1 extraction — the split is a pure
        scattering-kernel split now.)

        Returns
        -------
        ScatteringOperator
            An order-0 sibling whose per-material kernel is the diagonal
            P0 head, on the same bound ends.
        """
        fold = ScatteringMaterialField(
            per_material={
                mid: ScatteringKernel(
                    moments=(np.diag(np.diag(k.p0)),),
                )
                for mid, k in self.scattering.per_material.items()
            },
            cells_by_material=self.scattering.cells_by_material,
        )
        return self._sibling(fold)

    def residual_part(self) -> "ScatteringOperator":
        r"""Return the non-foldable sibling of :math:`S`.

        Carries everything :meth:`foldable_part` does not: the
        off-diagonal of each material's :math:`\Sigma_{s,0}` (cross-group
        P0) and every :math:`P_\ell \ge 1` block verbatim.

        Notes
        -----
        Algebraic contract:
        ``S.apply(ψ) ≈ S.foldable_part().apply(ψ) +
        S.residual_part().apply(ψ)`` at ``rtol=1e-14``.
        """
        residual = ScatteringMaterialField(
            per_material={
                mid: ScatteringKernel(
                    moments=(
                        k.p0 - np.diag(np.diag(k.p0)),
                        *k.moments[1:],
                    ),
                )
                for mid, k in self.scattering.per_material.items()
            },
            cells_by_material=self.scattering.cells_by_material,
        )
        return self._sibling(residual)

    def is_foldable_into_sigma_r(self) -> bool:
        r"""``True`` iff this operator is structurally the
        :meth:`foldable_part` of some parent :math:`S` — order 0 with
        every material's P0 diagonal. (The pre-§14.1 zero-(n,2n) clause
        dissolved with the extraction: S carries no (n,2n) channel.)
        """
        if self.scattering_order != 0:
            return False
        return all(
            np.allclose(k.p0, np.diag(np.diag(k.p0)))
            for k in self.scattering.per_material.values()
        )

    def foldable_sigma(self) -> dict[int, np.ndarray]:
        r"""The per-material foldable cross-section
        :math:`(\sigma_{s,0}^{g\to g})_g` — ``{mid: (ng,)}``, fresh
        copies, read off the bound kernel field."""
        return {
            mid: np.array(np.diag(k.p0))
            for mid, k in self.scattering.per_material.items()
        }

    # ── LinearOperator surface ─────────────────────────────────────────

    @singledispatchmethod
    def _apply_impl(self, psi) -> "Any":
        r"""Runtime dispatch for :meth:`apply` — see the typed overloads.

        Applies the full scattering source :math:`S\,\psi`, dispatched on the
        input *carrier* type via :func:`functools.singledispatchmethod`. The
        public :meth:`apply` aliases this dispatcher and carries the per-carrier
        ``@overload`` surface, so callers statically see the output type per
        input carrier:

        * :class:`~orpheus.transport.timed_full_field.TimedFullField` /
          :class:`~orpheus.transport.full_field.FullField`
          → :class:`~orpheus.transport.full_field.FullField` — composite bulk +
          boundary variant. Bulk = the full :math:`P_\ell` path; boundary =
          implicit-zero (scattering is volumetric, ``block_role = BULK``).
        * :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`
          → :class:`~orpheus.transport.source_sinks.ScalarSourceSink` —
          :math:`P_0` + :math:`(n,2n)` only, in **iso scalar magnitude** (no
          :math:`P_\ell`; no :math:`1/W` — scalar consumers do not project to
          per-ordinate).
        * :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
          → :class:`~orpheus.transport.source_sinks.AngularSourceSink` —
          full :math:`P_\ell` Galerkin in **per-ordinate magnitude** (the
          trailing :math:`1/W` lives at this producer boundary).
        * :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
          → :class:`~orpheus.transport.source_sinks.AngularSourceSink` — the
          angular-windowing path: :math:`S` consumes flux MOMENTS (already
          :math:`M\psi`, so :math:`M` is skipped), bit-identical to the
          :class:`AngularFlux` arm for :math:`\phi = M\psi`.

        The internal helpers :meth:`add_iso_source` and
        :meth:`build_aniso_source` remain available for callers that need
        the iso / aniso pieces separately (the (n,2n) verbs live on
        :class:`~orpheus.transport.operators.n2n.N2NOperator` since §14.1).
        """
        raise TypeError(
            f"ScatteringOperator.apply: unsupported input type "
            f"{type(psi).__name__}; expected TimedFullField, ScalarFlux, "
            f"AngularFlux, or HarmonicMomentFlux.  Dispatch table is "
            f"registered via @singledispatchmethod."
        )

    @_apply_impl.register
    def _(self, psi: FullField) -> "FullField":
        r"""Composite :class:`FullField` variant — bulk-only scattering source.

        Registered on the timeless :class:`FullField`: a
        :class:`TimedFullField` iterate dispatches here via MRO (it IS a
        ``FullField``), and a bare ``FullField`` dispatches correctly; the arm
        reads only ``psi.interior`` (history-blind). Bulk follows the same math
        as the :class:`AngularFlux` arm; the boundary is the **implicit-zero**
        :class:`~orpheus.transport.source_sinks.AngularBoundarySourceSink` —
        scattering is volumetric (``block_role = BlockRole.BULK``), no
        face-trace contribution, so ``(L + C - S - B)`` composes under
        :meth:`TimedFullField.__sub__`.
        """
        # Delegate the bulk source to the bulk-type dispatch arm and wrap with
        # the implicit-zero boundary. ``psi.interior`` is either the full-angular
        # AngularFlux or the windowed HarmonicMomentFlux (both return an
        # AngularSourceSink); the cast names that runtime truth so the typed
        # ``apply`` overloads resolve.
        combined = self.apply(cast("AngularFlux | HarmonicMomentFlux", psi.interior))
        # S is PURE BULK (#282 route (a)): the (ray, bulk) closed-μ emission
        # lives on the SN coupling operator
        # :class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicEmission`,
        # not here — see docs/theory/foundations/coupled_block_operator.rst.
        return FullField(
            interior=combined,
            boundary=AngularBoundarySourceSink.zeros(psi.boundary.space),
        )

    @_apply_impl.register
    def _(self, phi: ScalarFlux) -> "ScalarSourceSink":
        r"""Typed ScalarFlux variant — iso scalar magnitude output (P0 + n2n only).

        :math:`Q_g = \Sigma_{s,0}(g'\to g)\,\phi_{g'} + 2\,\Sigma_{2n}(g'\to g)
        \,\phi_{g'}`. No :math:`P_\ell` (scalar flux lacks angular info); no
        :math:`1/W` (scalar consumers — diffusion / CP / kinetics — do not
        project to per-ordinate).

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
        # shared mint; role is class identity). P0 only since §14.1 — the
        # (n,2n) term is N2NOperator's.
        iso: ScalarSourceSink = ScalarSourceSink.zeros(phi.space)
        iso = self.add_iso_source(iso, phi)
        return iso

    @_apply_impl.register
    def _(self, psi: AngularFlux) -> "AngularSourceSink":
        r"""Typed :class:`AngularFlux` variant — per-ordinate magnitude output.

        Reduce ``psi`` angular → scalar, build the iso :math:`P_0 + (n,2n)`
        source and the per-ordinate :math:`P_\ell\ge 1` Galerkin contribution
        (:meth:`build_aniso_source`), then combine via the producer-side
        :math:`1/W` in :meth:`_assemble_per_ordinate_source`.
        """
        # φ = ∫ψ dΩ (scalar), aniso = (1/W) RΛM ψ (per-ordinate), then the
        # shared producer-side assembly. The iso (P0 + n2n) keeps the cheap
        # reaction-rate fast path (NO moment tensor) — a load-bearing PERF
        # optimisation on the SI-sweep hot path; routing it through the frame
        # regresses LD/P0 badly. The frame form lives on as
        # :attr:`full_scatter_kernel` for the (non-hot-path) adjoint transpose
        # only — see docs/theory/methods/sn/adjoint.rst §sn-scattering-adjoint-source.
        return self._assemble_per_ordinate_source(
            psi.integrate_angular(), self.build_aniso_source(psi), psi.space,
        )

    @_apply_impl.register
    def _(self, phi_moments: HarmonicMomentFlux) -> "AngularSourceSink":
        r"""Windowed moment-iterate variant — :math:`S` consumes flux MOMENTS.

        The angular-windowing path: when the within-group SI iterate is stored
        as harmonic moments :math:`\phi_\ell^m` (the 2-D Cartesian windowed
        iterate) instead of the full per-ordinate :class:`AngularFlux`,
        :math:`S` consumes the moments WITHOUT the :math:`M` projection — the
        moments ARE :math:`M\psi`, so projecting again is redundant.

        Structurally parallel to the :class:`AngularFlux` arm and bit-identical
        to it for :math:`\phi = M\psi`: the :math:`\ell=0` moment IS the scalar
        flux (:math:`Y_0^0 = 1`, so :meth:`HarmonicMomentFlux.scalar_flux`
        equals :meth:`AngularFlux.integrate_angular`), feeding the identical
        P0 + (n,2n) fast path; the :math:`\ell\ge 1` aniso takes the explicit
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
        if self.scattering_order == 0:
            aniso = None
        else:
            # Explicit typed grid path: Λ scatters flux moments → source moments
            # (HarmonicMomentFlux → HarmonicMomentSourceSink, the role-changing
            # edge), the minted source-reconstruction FACE synthesises the
            # per-ordinate AngularSourceSink (riding its bound angular
            # codomain), then the producer-side 1/W. Numerically equals the
            # kernel's ndarray reconstruct_after(Λ) reference.
            scattered = self._moment_scattering(skip_l0=True).apply(
                phi_moments,
            )
            aniso = self.source_reconstruction.apply(
                scattered,
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
        # Honest per-carrier typing surface (#257 S8c).  ``ScatteringOperator``
        # is NOT an endomorphism ``V -> V`` (the mixin's nominal contract):
        # it maps each input carrier to a DISTINCT output carrier.  These
        # ``@overload`` stubs exist only for the type checker; the public
        # ``apply`` IS the runtime dispatcher (``apply = _apply_impl`` below),
        # so callers statically see e.g. ``S.apply(ScalarFlux) ->
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
        r"""The adjoint scattering source :math:`S^{T}\chi =
        (1/W)\,\mathrm{full\_scatter\_kernel}^{T}\chi` (closes
        `#118 <https://github.com/deOliveira-R/ORPHEUS/issues/118>`_).

        :math:`S^{T}` is the group-and-angle transpose the adjoint transport
        equation :math:`(L+C-S)^{T}\psi^{*}=q^{*}` needs. The production FORWARD
        keeps the scalar fast-path (:attr:`isotropic_energy` +
        :meth:`build_aniso_source`) for SI-sweep performance, so the adjoint —
        NOT the hot path — rides the validated harmonic-frame
        :attr:`full_scatter_kernel`, whose transpose falls out of
        :meth:`~orpheus.numerics.operator.OperatorProduct.apply_transpose` in
        ONE expression (iso :math:`\ell=0` + aniso :math:`\ell\ge1` + (n,2n)).

        This is the **Euclidean** transpose (L12) — NOT the metric Hilbert
        adjoint ``.H`` (which would carry the angular Gram). The
        forward/adjoint structural asymmetry (which makes the reciprocity
        :math:`\langle S\psi,\chi\rangle=\langle\psi,S^{T}\chi\rangle` a genuine
        cross-check, not a tautology) and the proof are in
        ``docs/theory/methods/sn/adjoint.rst §sn-scattering-adjoint-source``.

        The COMPOSITE (``FullField``) arm mirrors the forward lift: scattering
        is volumetric, so :math:`S^{T}` reads ONLY the bulk cotangent
        (``chi.interior``) and emits the implicit-zero trace, letting the full
        within-group loss ``(L + C - S - B).H`` compose through
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
            # S is PURE BULK, so S^T is too (#282 route (a)): the seed-cotangent
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
            self.full_scatter_kernel.apply_transpose(chi_values)
            / self.total_weight
        )
