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
from orpheus.numerics.spaces.full_field_space import FullFieldSpace
from orpheus.transport.frames import HarmonicFrame
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
    from orpheus.numerics.quadrature import Quadrature


__all__ = ["LegendreMomentScattering", "ScatteringOperator"]


@dataclass(frozen=True)
class LegendreMomentScattering(LinearOperator):
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

    with :math:`m(\vec r)` the material id at cell :math:`\vec r`
    (per-material structure folded into the cell axis via ``cells_by_mat``).

    ``skip_l0`` (default ``True``) skips the :math:`\ell = 0` block, which
    the project's P0 in-scatter handles on a separate reaction-rate fast
    path (:meth:`ScatteringOperator.add_iso_source`). Set ``False`` for the
    full :math:`R\Lambda M\psi` composition on the LinearOperator surface.

    Capability set: ``{apply, apply_transpose}``; no efficient ``solve``
    (rank-deficient on the :math:`\ell = 0` block by design).
    :math:`\Lambda^{T}` (:meth:`apply_transpose`) is the per-ℓ group-axis
    transpose — the ONLY group-asymmetric factor of the kernel
    :math:`R\circ\Lambda\circ M`, so
    :math:`(R\circ\Lambda\circ M)^{T} = M^{T}\circ\Lambda^{T}\circ R^{T}`
    (``docs/theory/methods/sn/adjoint.rst §sn-scattering-adjoint-source``).

    Parameters
    ----------
    mat_xs : MaterialXSField
        Macroscopic XS field carrying the per-material Legendre scattering
        matrices and the cell-to-material map. The per-material dispatch
        lives inside :meth:`MaterialXSField.apply_legendre_scattering_moments`.
    L : int
        Maximum Legendre order :math:`L` retained.
    skip_l0 : bool, default ``True``
        Skip the :math:`\ell = 0` block (handled by P0 in-scatter). Set
        ``False`` for the full :math:`R\Lambda M\psi` composition.
    """

    mat_xs: "MaterialXSField"
    L: int
    skip_l0: bool = True

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
            :attr:`~ScatteringOperator.kernel` ``OperatorProduct`` composes on).

        Returns
        -------
        np.ndarray or HarmonicMomentSourceSink
            The scattered moment source, same shape as ``moments``.  The
            :math:`\ell = 0` block is zero when ``skip_l0`` is ``True``;
            otherwise the P0 in-scatter contribution is included.  Typed in →
            typed source out; ndarray in → ndarray out.

        Notes
        -----
        Both arms route through the single per-material kernel
        :meth:`MaterialXSField.apply_legendre_scattering_moments`; they differ
        only in the carrier wrap.
        """
        from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
        if isinstance(moments, HarmonicMomentFlux):
            out_values = self.mat_xs.apply_legendre_scattering_moments(
                moments.values, L=self.L, skip_l0=self.skip_l0,
            )
            # flux moment → source moment: the explicit role change
            # (CS4b S4 — same space, new class; role is class identity).
            return HarmonicMomentSourceSink(
                values=out_values, space=moments.space, L=moments.L,
                spatial_moments=moments.spatial_moments,
            )
        return self.mat_xs.apply_legendre_scattering_moments(
            moments, L=self.L, skip_l0=self.skip_l0,
        )

    def apply_transpose(
        self, moments: "np.ndarray | HarmonicMomentSourceSink",
    ) -> "np.ndarray | HarmonicMomentFlux":
        r"""Apply :math:`\Lambda^{T}` — the per-ℓ group-transpose (the role-REVERSING edge).

        The bare Euclidean transpose of :meth:`apply`: :math:`\Lambda^{T}` maps a
        source moment back into the flux-moment space it scattered from (source →
        flux), transposing the per-ℓ group-transfer
        :math:`\Sigma_{s,\ell}(g'\to g) \mapsto (g\to g')`.  Routes through the
        transpose twin
        :meth:`~orpheus.transport.mesh.material_xs_field.MaterialXSField.apply_legendre_scattering_moments_transpose`.

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
            out_values = self.mat_xs.apply_legendre_scattering_moments_transpose(
                moments.values, L=self.L, skip_l0=self.skip_l0,
            )
            return HarmonicMomentFlux(
                values=out_values, space=moments.space, L=moments.L,
                spatial_moments=moments.spatial_moments,
            )
        return self.mat_xs.apply_legendre_scattering_moments_transpose(
            moments, L=self.L, skip_l0=self.skip_l0,
        )

    @property
    def domain(self) -> "FunctionSpace":
        r"""The spherical-harmonic coefficient space — :math:`\Lambda` is endomorphic.

        :math:`\Lambda` acts block-diagonally per :math:`\ell` ON moment space, so
        domain and codomain are BOTH the SH coefficient space
        :class:`~orpheus.numerics.spaces.SphericalHarmonicSpace` of order :attr:`L`
        (``== frame.basis_space``). Carrying real spaces lets
        :class:`~orpheus.numerics.operator.OperatorProduct` admit the ``R∘Λ∘M``
        composition natively (the composability guard validates it).
        """
        from orpheus.numerics.spaces import SphericalHarmonicSpace

        return SphericalHarmonicSpace.from_L(self.L)

    @property
    def codomain(self) -> "FunctionSpace":
        r"""Endomorphic on the SH coefficient space — codomain ``==`` :attr:`domain`."""
        return self.domain


@dataclass(frozen=True)
class N2NMomentOperator(LinearOperator):
    r"""The (n,2n) isotropic :math:`\ell=0` moment operator :math:`2\,\Sigma_{2n}`.

    The (n,2n) reaction is a DISTINCT isotropic (:math:`\ell=0`) group transfer —
    a *multiplication* channel (each event emits two neutrons), NOT scattering —
    so it is its own named operator, summed with :math:`\Lambda` in moment space
    (an :class:`~orpheus.numerics.operator.OperatorSum`) and frame-conjugated
    with it into the one full in-scatter source :math:`R\circ(\Lambda + N_{2n})
    \circ M` (:attr:`ScatteringOperator.full_scatter_kernel`). Keeping the
    multiplication reaction a visible distinct operator, rather than hidden in
    the scattering matmul, is the physics-faithful choice — see
    ``docs/theory/methods/sn/adjoint.rst §sn-scattering-adjoint-source``.

    Endomorphic on the spherical-harmonic coefficient space of order :attr:`L`
    (it reads/writes ONLY the :math:`\ell=0` block, so it composes in an
    ``OperatorSum`` with :math:`\Lambda` on the same space); per-material dispatch
    lives in :meth:`MaterialXSField.apply_n2n_moments`.
    """

    mat_xs: "MaterialXSField"
    L: int

    @property
    def is_adjointable(self) -> bool:
        # 2Σ_{2n}^T (apply_transpose) is the ℓ=0 group-transpose; caps ⊇
        # apply_transpose. is_invertible inherits base False.
        return True

    def apply(self, moments: np.ndarray) -> np.ndarray:
        r""":math:`2\,\Sigma_{2n}` applied to the :math:`\ell=0` moment (ℓ≥1 zero).

        Bare ndarray (the endomorphic moment-space view the ``frame.conjugate``
        ``OperatorProduct`` chain composes on — both for the forward
        :attr:`~ScatteringOperator.full_scatter_kernel` and its
        :meth:`apply_transpose`).
        """
        return self.mat_xs.apply_n2n_moments(moments)

    def apply_transpose(self, moments: np.ndarray) -> np.ndarray:
        r"""The :math:`\ell=0` group-transpose :math:`(2\,\Sigma_{2n})^{T}` (bare ndarray)."""
        return self.mat_xs.apply_n2n_moments_transpose(moments)

    @property
    def domain(self) -> "FunctionSpace":
        r"""The SH coefficient space of order :attr:`L` — :math:`N_{2n}` is endomorphic."""
        from orpheus.numerics.spaces import SphericalHarmonicSpace

        return SphericalHarmonicSpace.from_L(self.L)

    @property
    def codomain(self) -> "FunctionSpace":
        r"""Endomorphic — codomain ``==`` :attr:`domain`."""
        return self.domain


@dataclass
class ScatteringOperator(LinearOperator):
    r"""Scattering source operator :math:`S` (P0 + Pℓ + (n,2n)).

    Holds the precomputed per-material Legendre scattering matrices
    :math:`\Sigma_{s,\ell}`, the (n,2n) matrices :math:`\Sigma_{2n}`,
    the precomputed real spherical harmonics :math:`Y_\ell^m` evaluated
    at every quadrature direction, and a per-material cell-index map
    for vectorised assembly. Use :meth:`from_solver_data` to build
    instances from a :class:`MaterialXSField` + quadrature.

    Attributes
    ----------
    mat_xs : MaterialXSField
        Macroscopic XS field — the single source of truth for both
        per-material scattering / (n,2n) data AND the cell-to-material
        topology.  Every per-material loop routes through a typed verb
        (``mat_xs.apply_p0_in_scatter``, ``apply_n2n``, …) which
        encapsulates the dispatch.
    quadrature : Quadrature
        The angular quadrature.  Carries ``N``, ``weights``, and
        :meth:`spherical_harmonics`.
    scattering_order : int
        Maximum Legendre order :math:`L` retained. ``0`` means P0 only.

    Capability surface: ``{apply, apply_transpose}`` — no efficient
    ``solve``; the adjoint :math:`S^{T}` is free via the harmonic-frame
    :attr:`full_scatter_kernel` (see :meth:`apply_transpose`).
    """

    mat_xs: "MaterialXSField"
    quadrature: "Quadrature"
    scattering_order: int

    # Scattering is a BULK operator — the moment-folding `Σ_s · ⟨P_ℓ, ψ⟩`
    # reads and writes the bulk flux only (A_bb), no boundary action.
    # Class-level constant (unannotated so the dataclass does not treat it
    # as a field).
    block_role = BlockRole.BULK

    # Lazy cache for the precomputed spherical harmonics — only
    # populated when ``scattering_order > 0`` (avoids paying the
    # cost on P0-only problems).
    _Y_cached: np.ndarray | None = field(
        default=None, init=False, repr=False,
    )

    #: The space this operator acts on (renamed from ``full_field_space`` in
    #: campaign 1 CS1 — the slot names the ROLE; typed ``FunctionSpace |
    #: None``, the operator family's uniform slot type, so the whole family
    #: greps as one pattern. What actually flows today is the SN composite:
    #: threaded from the solver's ``sn_mesh.full_field_space`` via
    #: :meth:`from_solver_data`; ``None`` for the bare/test constructor
    #: (then ``domain``/``codomain`` report ``None`` and the composition guard
    #: skips — backward-compatible). When present it is the SAME instance
    #: ``L``/``C``/``B`` carry, so the within-group ``(L+C) - S``
    #: :class:`~orpheus.numerics.operator.OperatorSum` guard validates the
    #: ``- S`` arm natively (the load-bearing guard is instance AGREEMENT,
    #: not the annotation's family). ``S`` depends on this numerics
    #: ``FunctionSpace``, NOT on an SN mesh (``S`` scatters in every method).
    space: "FunctionSpace | None" = field(
        default=None, repr=False, compare=False,
    )

    # ── Operator-algebra space metadata (P4.5 W-D) ───────────────────
    @property
    def domain(self) -> "FunctionSpace | None":
        r"""The composite full-field space, or ``None`` if unthreaded.

        Although :math:`S` is a BULK operator, it composes into the within-group
        loss ``(L+C) - S`` as a composite-field operator, so it advertises the
        SAME :attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.full_field_space`
        instance ``L``/``C``/``B`` carry (threaded via :meth:`from_solver_data`),
        which lets the residual
        :class:`~orpheus.numerics.operator.OperatorSum` guard VALIDATE the
        ``- S`` arm instead of silently skipping it. The ``None``-spaced default
        keeps bare/test construction backward-compatible. :math:`S` reads and
        writes the bulk block only, so domain == codomain on the composite.
        """
        return self.space

    @property
    def codomain(self) -> "FunctionSpace | None":
        # Endomorphic on the composite full-field space (see :meth:`domain`).
        return self.space

    @property
    def is_adjointable(self) -> bool:
        # S = R∘(Λ+N2N)∘M exposes its Euclidean transpose S^T via
        # :attr:`full_scatter_kernel` (the OperatorProduct adjoint chain);
        # is_invertible inherits base False —
        # a scattering source operator is not invertible.
        return True

    # ── Convenience read-throughs (the legacy n_ordinates / nx / ny
    # / ng / weights / Y / cells_by_mat / sig_s / sig2 / sig_s0
    # attributes the test suite + internal call sites consumed).
    # Issue #197 PR-TYPED-1: these become read-throughs onto
    # mat_xs + quadrature so consumers don't break.  Marked TRANSIENT
    # in source comments where appropriate; full retirement deferred
    # to a follow-up PR once consumers all read mat_xs.* directly.

    @property
    def n_ordinates(self) -> int:
        """Number of angular ordinates :math:`N`."""
        return self.quadrature.N

    @property
    def spatial_shape(self) -> tuple[int, ...]:
        """Per-axis cell counts — read-through from :attr:`mat_xs`.

        C5.2 (#225): replaces the retired ``nx``/``ny`` pair (rank-
        generic; honest at any spatial dimension).
        """
        return self.mat_xs.spatial_shape

    @property
    def ng(self) -> int:
        """Energy group count."""
        return self.mat_xs.ng

    @property
    def weights(self) -> np.ndarray:
        """Quadrature weights ``(N,)``."""
        return self.quadrature.weights

    @property
    def Y(self) -> np.ndarray | None:
        r"""Real spherical harmonics :math:`Y_\ell^m(\Omega_n)`,
        shape ``(N, L+1, 2L+1)``, or ``None`` when ``L == 0``."""
        if self.scattering_order == 0:
            return None
        if self._Y_cached is None:
            self._Y_cached = self.quadrature.spherical_harmonics(
                self.scattering_order,
            )
        return self._Y_cached

    @cached_property
    def frame(self) -> HarmonicFrame:
        r"""The angular discrete frame — :class:`SphericalHarmonicBasis`
        (order :math:`L=` :attr:`scattering_order`) bound to the quadrature's
        :math:`S^2` measure.

        The single source of the analysis (:math:`M`) and reconstruction
        (:math:`R`) faces that realise the anisotropic Legendre redistribution
        :math:`R\circ\Lambda\circ M` (:attr:`kernel`); both the §5.6
        :attr:`kernel` and the angular-windowing in-sweep moment accumulation
        read THIS frame, so the projection table is shared term-for-term.

        A :class:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame` — the
        carrier-typed specialization of the generic
        :class:`~orpheus.numerics.frame.GalerkinFrame`
        ``quadrature.angular_frame(L)``, adding the typed verbs
        :meth:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame.analyse`
        (:math:`M`) and
        :meth:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame.reconstruct`
        (:math:`R`) the windowed in-scatter arm uses. ``frame.table`` is the
        :math:`Y_\ell^m(\hat\Omega_n)` tabulation (equal to :attr:`Y`). The
        HarmonicFrame-IS-A-GalerkinFrame relation and the shared-table
        rationale: ``docs/theory/foundations/operator_algebra.rst
        §integral-kernel-category``.
        """
        return HarmonicFrame.from_galerkin(
            self.quadrature.angular_frame(self.scattering_order),
        )

    @property
    def sig_s(self) -> dict[int, list[np.ndarray]]:
        """TRANSIENT — per-material dense Legendre scattering dict.

        Read-through onto :meth:`MaterialXSField.sig_s_legendre`, pending
        retirement — tracked in issue #306 (its original ``_build_rhs_*``
        consumers in :mod:`orpheus.sn.solver` are gone; kept until remaining
        callers read ``mat_xs.*`` directly).
        """
        return {
            mid: self.mat_xs.sig_s_legendre(mid)
            for mid in self.mat_xs.materials
        }

    @property
    def sig2(self) -> dict[int, np.ndarray]:
        """TRANSIENT — per-material dense (n,2n) dict.  See :attr:`sig_s`."""
        return {
            mid: self.mat_xs.n2n_matrix(mid)
            for mid in self.mat_xs.materials
        }

    @property
    def sig_s0(self) -> dict[int, np.ndarray]:
        """TRANSIENT — per-material P0 scattering matrix dict.
        See :attr:`sig_s`."""
        return {
            mid: self.mat_xs.sig_s_legendre(mid)[0]
            for mid in self.mat_xs.materials
        }

    @property
    def cells_by_mat(self) -> dict[int, tuple[np.ndarray, ...]]:
        """TRANSIENT — per-material cell-index dict (one index array per
        mesh axis).  See :attr:`sig_s`."""
        return self.mat_xs.cells_by_material

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
        scatter = LegendreMomentScattering(
            mat_xs=self.mat_xs, L=self.scattering_order, skip_l0=True,
        )
        # R∘Λ as ONE typed operator (M already applied — the moments ARE M·ψ).
        return self.frame.reconstruct_after(scatter).apply(moment_values)

    @property
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
        why P0 + (n,2n) are the local components) is in
        ``docs/theory/foundations/operator_algebra.rst §integral-kernel-category``.
        Built fresh per access (read-through to :attr:`mat_xs` / :attr:`Y`).

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
            LegendreMomentScattering(
                mat_xs=self.mat_xs, L=self.scattering_order, skip_l0=True,
            )
        )

    @property
    def full_scatter_kernel(self) -> OperatorProduct:
        r"""The FULL in-scatter source kernel :math:`R\circ(\Lambda_{\ell\ge 0} + N_{2n})\circ M`.

        The COMPLETE P0 + anisotropic + (n,2n) in-scatter as ONE frame-conjugated
        operator: the isotropic ℓ=0 scattering, the anisotropic ℓ≥1
        redistribution (both :class:`LegendreMomentScattering`, ``skip_l0=False``),
        and the (n,2n) multiplication (the DISTINCT :class:`N2NMomentOperator`) are
        summed in moment space and conjugated by the frame TOGETHER. The
        per-ordinate source is ``(1/W)·full_scatter_kernel.apply(ψ)``; its
        transpose ``(1/W)·full_scatter_kernel.apply_transpose(ψ*)`` is the adjoint
        :math:`S^{T}` (:meth:`apply_transpose`). Riding the same frame conjugation
        for iso and aniso is what lets the whole transpose fall out for free —
        ``docs/theory/methods/sn/adjoint.rst §sn-scattering-adjoint-source``.

        Distinct from :attr:`kernel` (the §5.6 ℓ≥1 ANISOTROPIC subcomponent + the
        0-ULP crosscheck oracle): this is the FULL ℓ≥0 + (n,2n) source. Built
        fresh per access (read-through to :attr:`mat_xs`).
        """
        return self.frame.conjugate(
            LegendreMomentScattering(
                mat_xs=self.mat_xs, L=self.scattering_order, skip_l0=False,
            )
            + N2NMomentOperator(mat_xs=self.mat_xs, L=self.scattering_order)
        )

    @cached_property
    def isotropic_kernel(self) -> "OperatorSum":
        r"""The model-shared isotropic in-scatter energy operator
        :math:`\Sigma_{s,0} + 2\,\Sigma_{2n}` on the scalar flux.

        The :math:`\ell=0` energy source — P0 in-scatter plus the :math:`(n,2n)`
        doubling — as the cross-model
        :class:`~orpheus.transport.operators.isotropic_scattering.IsotropicScattering`
        ``+``
        :class:`~orpheus.transport.operators.isotropic_scattering.IsotropicN2N`
        sum. The same energy operation is shared by every transport model
        (CP / MoC / diffusion / homogeneous / MC); the model-shared K_iso
        narrative is in ``docs/theory/methods/sn/adjoint.rst
        §sn-scattering-adjoint-source`` and
        ``docs/theory/foundations/infinite_medium.rst``.

        Annotated as the concrete :class:`~orpheus.numerics.operator.OperatorSum`
        (not the bare ``LinearOperator`` erasure) so the sum surface's
        ``apply_transpose`` (both leaves implement it) stays visible to the
        checker — the starting-direction seed arm's transpose consumes it.
        Space-anonymous (``space=None``: no composition-guard space at this
        producer-side use). Cached, but reads the XS *through* :attr:`mat_xs`
        at ``apply`` time, so depletion / feedback updates show up immediately.
        """
        from orpheus.transport.operators.isotropic_scattering import (
            IsotropicN2N,
            IsotropicScattering,
        )
        return IsotropicScattering(self.mat_xs) + IsotropicN2N(self.mat_xs)

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
        # Route the isotropic (P0 + (n,2n)) energy source through the
        # model-shared K_iso operators (isotropic_kernel = Σ_s0 + 2Σ_2n);
        # the iso rides the driving flux's own space (width travels in it).
        iso: ScalarSourceSink = ScalarSourceSink(
            values=self.isotropic_kernel.apply(phi), space=phi.space,
        )
        aniso = (
            aniso_or_none
            if aniso_or_none is not None
            else AngularSourceSink.zeros(angular_space)
        )
        sum_w = float(self.quadrature.weights.sum())
        # The containment dunder's cross-class arm returns the LARGER
        # (angular) class — the #288 principled LSP exception the static
        # union cannot carry.
        return cast("AngularSourceSink", (iso / sum_w) + aniso)

    @classmethod
    def from_solver_data(
        cls,
        *,
        mat_xs: "MaterialXSField",
        quadrature: "Quadrature",
        scattering_order: int,
        space: "FunctionSpace | None" = None,
    ) -> "ScatteringOperator":
        """Construct from a :class:`MaterialXSField` + quadrature.

        ``space`` (P4.5 W-D) is the composite
        :attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.full_field_space` the solver
        threads so ``S``'s ``domain``/``codomain`` match ``L``/``C``/``B``
        and the within-group ``(L+C) - S`` composition guard validates the
        ``- S`` arm. ``None`` (the default) leaves ``S`` space-less — the
        guard skips it, preserving the legacy contract for direct callers.

        The :class:`MaterialXSField` carries everything per-material plus the
        spatial topology; the :class:`Quadrature` carries ``N`` / ``weights`` /
        harmonics.
        """
        return cls(
            mat_xs=mat_xs,
            quadrature=quadrature,
            scattering_order=scattering_order,
            space=space,
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
        :meth:`MaterialXSField.apply_p0_in_scatter`.
        """
        from orpheus.transport.fields.scalar_flux import ScalarFlux
        from orpheus.transport.source_sinks import ScalarSourceSink
        phi_values = phi.values if isinstance(phi, ScalarFlux) else phi
        if isinstance(Q, ScalarSourceSink):
            Q_values = Q.values.copy()
            self.mat_xs.apply_p0_in_scatter(Q_values, phi_values)
            # Preserve Q's spatial-moment width (#240 D5b-S3 — the φ̂ accumulator
            # carries the trailing 2^d axis; re-wrapping without it would lose
            # the typed factor and raise on the shape).
            return ScalarSourceSink(values=Q_values, space=Q.space)
        self.mat_xs.apply_p0_in_scatter(Q, phi_values)
        return None

    @overload
    def add_n2n_source(
        self, Q: "ScalarSourceSink", phi: "np.ndarray | ScalarFlux",
    ) -> "ScalarSourceSink": ...

    @overload
    def add_n2n_source(
        self, Q: np.ndarray, phi: "np.ndarray | ScalarFlux",
    ) -> None: ...

    def add_n2n_source(
        self,
        Q: "np.ndarray | ScalarSourceSink",
        phi: "np.ndarray | ScalarFlux",
    ) -> "ScalarSourceSink | None":
        r"""Add (n,2n) source to :math:`Q`.

        Vectorised by material: per cell ``c`` of material ``mid``,
        ``Q[:, ic, jc] += 2 · sig2[mid].T @ phi[:, ic, jc]``.

        Same typed-action overload as :meth:`add_iso_source` — raw-in mutates
        in place and returns ``None``; typed-in returns a fresh
        :class:`ScalarSourceSink`.

        Parameters
        ----------
        Q : np.ndarray or ScalarSourceSink
            Isotropic source carrier.
        phi : np.ndarray or ScalarFlux
            Scalar flux.

        Returns
        -------
        np.ndarray or ScalarSourceSink or None
            Raw-in returns ``None`` (legacy in-place); typed-in returns
            a fresh :class:`ScalarSourceSink`.

        Notes
        -----
        The per-material dispatch lives inside :meth:`MaterialXSField.apply_n2n`.
        """
        from orpheus.transport.fields.scalar_flux import ScalarFlux
        from orpheus.transport.source_sinks import ScalarSourceSink
        phi_values = phi.values if isinstance(phi, ScalarFlux) else phi
        if isinstance(Q, ScalarSourceSink):
            Q_values = Q.values.copy()
            self.mat_xs.apply_n2n(Q_values, phi_values)
            # Preserve Q's spatial-moment width (#240 D5b-S3, as add_iso_source).
            return ScalarSourceSink(values=Q_values, space=Q.space)
        self.mat_xs.apply_n2n(Q, phi_values)
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
        sum_w = float(self.weights.sum())
        out_values = self.kernel.apply(angular_flux.values) / sum_w
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

    def foldable_part(self) -> "ScatteringOperator":
        r"""Return the P0 within-group self-scatter sibling of :math:`S`.

        Carries only the diagonal of ``sig_s[mid][0]`` per material —
        the within-group self-scatter cross-section
        :math:`\Sigma_{s,0}^{g\to g}` per energy group. All other
        scattering channels (cross-group P0, all :math:`P_\ell \ge 1`,
        and (n,2n)) live in :meth:`residual_part`.

        Returns
        -------
        ScatteringOperator
            A sibling with ``scattering_order = 0``, diagonal-only
            ``sig_s[mid][0]``, zero ``sig2``, and a derived
            :class:`MaterialXSField` carrying the foldable overrides.

        Notes
        -----
        The per-material loop that builds the ``sig_s_foldable`` /
        ``sig2_foldable`` dicts lives inside
        :meth:`MaterialXSField.foldable_sig_s`.
        """
        sig_s_foldable = self.mat_xs.foldable_sig_s()
        sig2_foldable = {
            mid: np.zeros_like(self.mat_xs.n2n_matrix(mid))
            for mid in self.mat_xs.materials
        }
        derived_mat_xs = self.mat_xs.with_overridden_sig_s_and_n2n(
            sig_s_dense=sig_s_foldable,
            n2n_dense=sig2_foldable,
        )
        return ScatteringOperator(
            mat_xs=derived_mat_xs,
            quadrature=self.quadrature,
            scattering_order=0,
            space=self.space,
        )

    def residual_part(self) -> "ScatteringOperator":
        r"""Return the non-foldable sibling of :math:`S`.

        Carries everything :meth:`foldable_part` does not: the
        off-diagonal of ``sig_s[mid][0]`` (cross-group P0), every
        :math:`P_\ell \ge 1` block verbatim, and ``sig2[mid]``
        verbatim.

        Returns
        -------
        ScatteringOperator
            A sibling with ``scattering_order ==
            self.scattering_order``, zero-diagonal P0 matrix, and
            unchanged :math:`P_\ell \ge 1` + (n,2n) data.

        Notes
        -----
        Algebraic contract:
        ``S.apply(\psi) \approx S.foldable_part().apply(\psi) +
        S.residual_part().apply(\psi)`` at ``rtol=1e-14``.

        The per-material loop that builds the ``sig_s_residual`` dict lives
        inside :meth:`MaterialXSField.residual_sig_s`.
        """
        sig_s_residual = self.mat_xs.residual_sig_s()
        # (n,2n) carried verbatim — pull from mat_xs.
        sig2_residual = {
            mid: self.mat_xs.n2n_matrix(mid)
            for mid in self.mat_xs.materials
        }
        derived_mat_xs = self.mat_xs.with_overridden_sig_s_and_n2n(
            sig_s_dense=sig_s_residual,
            n2n_dense=sig2_residual,
        )
        return ScatteringOperator(
            mat_xs=derived_mat_xs,
            quadrature=self.quadrature,
            scattering_order=self.scattering_order,
            space=self.space,
        )

    def is_foldable_into_sigma_r(self) -> bool:
        r"""Return ``True`` iff this operator is structurally the
        :meth:`foldable_part` of some parent :math:`S`.

        Returns
        -------
        bool
            ``True`` iff ``scattering_order == 0`` AND every material's
            P0 is diagonal AND every material's (n,2n) is zero.

        Notes
        -----
        The per-material allclose checks live inside
        :meth:`MaterialXSField.is_p0_diagonal_with_zero_n2n`.
        """
        if self.scattering_order != 0:
            return False
        return self.mat_xs.is_p0_diagonal_with_zero_n2n()

    def foldable_sigma(self) -> dict[int, np.ndarray]:
        r"""Return the per-material foldable cross-section :math:`(\sigma_{s,0}^{g\to g})_g`.

        Returns
        -------
        dict[int, np.ndarray]
            ``{mid: (ng,) array}`` keyed by material id. Each array
            is a fresh copy.

        Notes
        -----
        Delegates to :meth:`MaterialXSField.foldable_sigma`.
        """
        return self.mat_xs.foldable_sigma()

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

        The internal helpers :meth:`add_iso_source`, :meth:`add_n2n_source`,
        and :meth:`build_aniso_source` remain available for callers that need
        the iso / aniso pieces separately.
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
        # shared mint; role is class identity).
        iso: ScalarSourceSink = ScalarSourceSink.zeros(phi.space)
        iso = self.add_iso_source(iso, phi)
        iso = self.add_n2n_source(iso, phi)
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
        # CS4b S4: the per-ordinate TARGET space comes from this
        # operator's own posed composite space (its interior — the
        # carrier's angular mint, width-widened on LD), because neither
        # the moment operand nor the (basis, measure) frame carries the
        # quadrature/scheme axes the target needs.
        if (
            not isinstance(self.space, FullFieldSpace)
            or self.space.interior_space is None
        ):
            raise TypeError(
                "ScatteringOperator's windowed moment arm needs the posed"
                " composite space (space=) to name its per-ordinate"
                " reconstruction target; this operator was built"
                " space-less."
            )
        angular_target = self.space.interior_space
        if self.scattering_order == 0:
            aniso = None
        else:
            # Explicit typed grid path: Λ scatters flux moments → source moments
            # (HarmonicMomentFlux → HarmonicMomentSourceSink, the role-changing
            # edge), the frame's R reconstructs the per-ordinate
            # AngularSourceSink, then the producer-side 1/W. Numerically equals
            # the kernel's ndarray reconstruct_after(Λ) reference.
            scattered = LegendreMomentScattering(
                mat_xs=self.mat_xs, L=self.scattering_order, skip_l0=True,
            ).apply(phi_moments)
            aniso = self.frame.reconstruct(
                scattered, space=angular_target,
            ) / float(self.weights.sum())
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
        keeps the scalar fast-path (:attr:`isotropic_kernel` +
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
        return self.full_scatter_kernel.apply_transpose(chi_values) / float(
            self.weights.sum()
        )
