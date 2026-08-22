r"""The :class:`HarmonicFrame` — the carrier-typed angular spherical-harmonic frame.

A thin transport-layer specialization of the numerics
:class:`~orpheus.numerics.frame.GalerkinFrame` that adds **carrier-typed** verbs
on top of the generic, ``np.ndarray``-valued analysis/reconstruction faces. The
generic Galerkin frame is the angular spherical-harmonic projection
``quadrature.angular_frame(L)``; this wrapper lets a consumer write the projection
as ``moments = frame.analyse(psi)`` — typed :class:`Field` in, typed :class:`Field`
out — instead of hand-unwrapping ``.values`` and re-wrapping the result at every
call site (the Frame campaign's "define Frame, analyse, done" payoff).

Why a transport-layer specialization (the layering)
===================================================

The two carriers the angular frame maps between —
:class:`~orpheus.transport.fields.angular_flux.AngularFlux` and
:class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux` (and
their source/sink siblings) — share their deepest primitive,
:class:`~orpheus.numerics.field.Field`, in **numerics**. But the part that makes
them *castable* — the ``mesh`` binding and the ``from_mesh`` / ``from_mesh_and_L``
factories — lives in the transport
:class:`~orpheus.transport.fields._bases.BulkField` base. And
``Quadrature.angular_frame`` (which builds the generic
:class:`~orpheus.numerics.frame.GalerkinFrame`) is in **numerics**, which cannot
import the transport carriers. So a generic numerics face *cannot* return a typed
``HarmonicMomentFlux`` without inverting the layer order. The clean home for the
casting is therefore one layer up, in transport: ``HarmonicFrame`` IS-A
``GalerkinFrame`` (Liskov — the angular SH projection, the canonical pure-Galerkin
frame), constructed from the generic frame's basis + measure via
:meth:`from_galerkin` (zero duplication of the basis/measure build), adding only
the typed verbs. The generic numerics faces stay untouched (carrier-agnostic,
shared with the P3 indicator-homogenisation frame, 0-ULP-safe).

Role-polymorphic verbs (the carrier grid)
=========================================

:meth:`analyse` (:math:`M`) and :meth:`reconstruct` (:math:`R`) **preserve the
role** (flux ↔ flux, source ↔ source) and change the **axis** (angular ↔ moment),
dispatching on the input carrier type over an enumerated isinstance chain
under typed ``@overload`` faces (the grid is architecturally CLOSED, so
callers get the precise per-carrier return and an off-grid carrier is a
call-site type error) — the two vertical edges of the
``(angular ⊗ moment) × (flux ⊗ source)`` grid::

    analyse:      AngularFlux ─────▶ HarmonicMomentFlux
                  AngularSourceSink ▶ HarmonicMomentSourceSink
    reconstruct:  HarmonicMomentFlux ─────▶ AngularFlux
                  HarmonicMomentSourceSink ▶ AngularSourceSink

The **role-changing** edge :math:`\Lambda` (scatter,
``HarmonicMomentFlux → HarmonicMomentSourceSink``) is NOT a frame verb — it is the
scattering operator's job (:class:`~orpheus.transport.operators.scattering.LegendreMomentScattering`),
where the physics puts the role change. The hot anisotropic-scatter kernel
:math:`S = (1/W)\,R\,\Lambda\,M` stays the single composed
:meth:`~orpheus.numerics.frame.FrameBase.conjugate` operator (one ``np.ndarray``
chain, the 0-ULP canary); these typed verbs serve the consumers that apply
:math:`M` or :math:`R` in isolation.

The frame binds only the angular basis + measure (mesh-agnostic). The output
carrier's SPACE is derived from the operand for analysis and supplied by the
caller for reconstruction (CS4b S4 — the frame square's structural asymmetry;
see the verbs); the truncation order ``L`` comes from the frame's own basis,
and the optional within-cell spatial-moment width is threaded through as a
typed factor (read off the input carrier's space, #240 D5b-S3).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, overload

from orpheus.numerics.basis.spherical_harmonic_basis import SphericalHarmonicBasis
from orpheus.numerics.frame import GalerkinFrame

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from orpheus.transport.source_sinks.angular_source_sink import AngularSourceSink
from orpheus.transport.source_sinks.harmonic_moment_source_sink import (
    HarmonicMomentSourceSink,
)


__all__ = ["HarmonicFrame"]


class HarmonicFrame(GalerkinFrame):
    r"""The angular spherical-harmonic :class:`GalerkinFrame` with typed verbs.

    Constructed from a generic ``quadrature.angular_frame(L)`` via
    :meth:`from_galerkin` (or directly from ``(basis, measure)`` — the
    constructor narrows the inherited signature to require a
    :class:`~orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis`:
    a harmonic frame over any other trial basis is an illegal state, and the
    truncation order ``L`` the typed verbs read exists only on the SH basis).
    Adds the role-polymorphic :meth:`analyse` (:math:`M`) and
    :meth:`reconstruct` (:math:`R`) carrier verbs; the generic ``np.ndarray``
    faces (:attr:`~orpheus.numerics.frame.FrameBase.analysis`,
    :attr:`~orpheus.numerics.frame.FrameBase.reconstruction`,
    :meth:`~orpheus.numerics.frame.FrameBase.conjugate`) are inherited
    unchanged.
    """

    # Covariant narrowing of the inherited frozen (read-only) field: the
    # constructor below guarantees it, and the typed verbs read ``basis.L``.
    basis: SphericalHarmonicBasis

    def __init__(
        self, basis: SphericalHarmonicBasis, measure: DiscreteMeasure,
    ) -> None:
        super().__init__(basis, measure)

    @classmethod
    def from_galerkin(cls, frame: GalerkinFrame) -> "HarmonicFrame":
        r"""Upgrade a generic angular :class:`GalerkinFrame` to a typed
        :class:`HarmonicFrame`, reusing its basis + measure (no rebuild).

        ``frame`` is the ``quadrature.angular_frame(L)`` Galerkin frame; the
        returned :class:`HarmonicFrame` shares its
        :class:`~orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis`
        trial basis and its angular :class:`~orpheus.numerics.measure.DiscreteMeasure`,
        so the table / spaces / faces are bit-identical — only the typed verbs
        are added. A frame over any other trial basis (e.g. an indicator
        basis) is rejected here, at the upgrade boundary — not later, when a
        typed verb first reads the SH-only truncation order ``L``.
        """
        basis = frame.basis
        if not isinstance(basis, SphericalHarmonicBasis):
            raise TypeError(
                f"HarmonicFrame.from_galerkin requires a spherical-harmonic "
                f"trial basis; got {type(basis).__name__}."
            )
        return cls(basis, frame.measure)

    # ── M : project (role-preserving, angular → moment) ──────────────
    #
    # The carrier grid is architecturally CLOSED (the two vertical edges of
    # the (angular ⊗ moment) × (flux ⊗ source) grid — module docstring), so
    # the verbs dispatch over an enumerated isinstance chain under typed
    # @overload faces: callers get the precise per-carrier return type, and
    # an off-grid carrier is a call-site type error before it is a runtime
    # TypeError. (An open singledispatch registry modeled extensibility the
    # grid deliberately does not have.)

    @overload
    def analyse(self, field: AngularFlux) -> HarmonicMomentFlux: ...

    @overload
    def analyse(self, field: AngularSourceSink) -> HarmonicMomentSourceSink: ...

    def analyse(
        self, field: AngularFlux | AngularSourceSink,
    ) -> HarmonicMomentFlux | HarmonicMomentSourceSink:
        r"""Project a per-ordinate carrier to its spherical-harmonic moments
        (:math:`M`), preserving the role.

        :class:`AngularFlux` → :class:`HarmonicMomentFlux`,
        :class:`AngularSourceSink` → :class:`HarmonicMomentSourceSink`. The
        output's space is DERIVED from the operand's (CS4b S4): the moment
        target is ``SH(L) * cell_group[* SMS]`` where the cell group is
        the operand's non-angular, non-moment axes — analysis's target is
        determined by its operand + the frame's ``L``, the structural
        half of the frame square that self-derives (reconstruction's
        angular target needs carrier knowledge and takes ``space=``).
        """
        if isinstance(field, (AngularFlux, AngularSourceSink)):
            space = self._moment_space_for(field)
            values = self.analysis.apply(field.values)
            per_axis = field.spatial_moments_per_axis
            if isinstance(field, AngularFlux):
                return HarmonicMomentFlux(
                    values=values, space=space, L=self.basis.L,
                    spatial_moments=per_axis,
                )
            return HarmonicMomentSourceSink(
                values=values, space=space, L=self.basis.L,
                spatial_moments=per_axis,
            )
        raise TypeError(
            f"HarmonicFrame.analyse: unsupported carrier "
            f"{type(field).__name__}; expected AngularFlux or "
            f"AngularSourceSink (a per-ordinate carrier)."
        )

    def _moment_space_for(
        self, field: AngularFlux | AngularSourceSink,
    ) -> "FunctionSpace":
        r"""The moment target space derived from a per-ordinate operand.

        ``SH.from_L(L) * of_axes(<non-angular, non-moment axes>)`` — the
        cell group is the operand's own energy/spatial axes (the same
        instances the carrier's mints share, so the product content-equals
        ``MomentField._space_for_mesh_and_L``'s) — with the densified
        ``SpatialMomentSpace`` factor appended for a widened operand,
        mirroring ``_compose_spatial_moments``'s densified arm.
        """
        from orpheus.numerics.moment_layout import SPATIAL_MOMENT_AXIS_LABEL
        from orpheus.numerics.space import FunctionSpace
        from orpheus.numerics.spaces.spatial_moment_space import (
            SpatialMomentSpace,
        )
        from orpheus.numerics.spaces.spherical_harmonic_space import (
            SphericalHarmonicSpace,
        )

        axes = field.space.axes
        if axes is None:  # unreachable for shipped angular spaces (S2)
            raise TypeError(
                "HarmonicFrame.analyse: the per-ordinate operand must ride "
                "an axis-built angular space."
            )
        cell_axes = [
            ax for ax in axes[1:] if ax.label != SPATIAL_MOMENT_AXIS_LABEL
        ]
        base = SphericalHarmonicSpace.from_L(self.basis.L) * (
            FunctionSpace.of_axes(*cell_axes)
        )
        per_axis = field.spatial_moments_per_axis
        if per_axis == 1:
            return base
        ndim = next(
            len(ax.shape) for ax in cell_axes if ax.label == "spatial"
        )
        return base * SpatialMomentSpace.from_per_axis(per_axis, ndim)

    # ── R : reconstruct (role-preserving, moment → angular) ──────────

    @overload
    def reconstruct(
        self, moment: HarmonicMomentFlux, *, space: "FunctionSpace",
    ) -> AngularFlux: ...

    @overload
    def reconstruct(
        self, moment: HarmonicMomentSourceSink, *, space: "FunctionSpace",
    ) -> AngularSourceSink: ...

    def reconstruct(
        self,
        moment: HarmonicMomentFlux | HarmonicMomentSourceSink,
        *,
        space: "FunctionSpace",
    ) -> AngularFlux | AngularSourceSink:
        r"""Reconstruct a per-ordinate carrier from its spherical-harmonic
        moments (:math:`R`), preserving the role.

        :class:`HarmonicMomentFlux` → :class:`AngularFlux`,
        :class:`HarmonicMomentSourceSink` → :class:`AngularSourceSink`.

        ``space`` is the per-ordinate TARGET space, supplied by the caller
        (CS4b S4): the angular target carries the quadrature axis and, on
        a widened iterate, the scheme's mass-bearing moment axis — carrier
        knowledge neither the moment operand nor this (basis, measure)
        frame holds. The caller that poses the problem holds it (the
        production caller passes its composite space's ``interior_space``).
        This is the frame square's structural asymmetry: analysis
        self-derives, reconstruction is told where to land.
        """
        if isinstance(moment, HarmonicMomentFlux):
            return AngularFlux(
                values=self.reconstruction.apply(moment.values), space=space,
            )
        if isinstance(moment, HarmonicMomentSourceSink):
            return AngularSourceSink(
                values=self.reconstruction.apply(moment.values), space=space,
            )
        raise TypeError(
            f"HarmonicFrame.reconstruct: unsupported carrier "
            f"{type(moment).__name__}; expected HarmonicMomentFlux or "
            f"HarmonicMomentSourceSink (a moment carrier)."
        )
