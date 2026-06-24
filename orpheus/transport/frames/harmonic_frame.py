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
dispatching on the input carrier type
(:func:`~functools.singledispatchmethod`) — the two vertical edges of the
``(angular ⊗ moment) × (flux ⊗ source)`` grid::

    analyse:      AngularFlux ─────▶ HarmonicMomentFlux
                  AngularSourceSink ▶ HarmonicMomentSourceSink
    reconstruct:  HarmonicMomentFlux ─────▶ AngularFlux
                  HarmonicMomentSourceSink ▶ AngularSourceSink

The **role-changing** edge :math:`\Lambda` (scatter,
``HarmonicMomentFlux → HarmonicMomentSourceSink``) is NOT a frame verb — it is the
scattering operator's job (:class:`~orpheus.sn.scattering.LegendreMomentScattering`),
where the physics puts the role change. The hot anisotropic-scatter kernel
:math:`S = (1/W)\,R\,\Lambda\,M` stays the single composed
:meth:`~orpheus.numerics.frame.FrameBase.conjugate` operator (one ``np.ndarray``
chain, the 0-ULP canary); these typed verbs serve the consumers that apply
:math:`M` or :math:`R` in isolation.

The output carrier inherits the input's ``mesh`` (the frame is mesh-agnostic — it
binds only the angular basis + measure) and the truncation order ``L`` from the
frame's own basis; the optional within-cell spatial-moment width is threaded
through as a typed factor (read off the input carrier's space, #240 D5b-S3).
"""

from __future__ import annotations

from functools import singledispatchmethod

from orpheus.numerics.frame import GalerkinFrame
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
    :meth:`from_galerkin` (or directly from ``(basis, measure)`` — the inherited
    :class:`~orpheus.numerics.frame.GalerkinFrame` constructor). Adds the
    role-polymorphic :meth:`analyse` (:math:`M`) and :meth:`reconstruct`
    (:math:`R`) carrier verbs; the generic ``np.ndarray`` faces
    (:attr:`~orpheus.numerics.frame.FrameBase.analysis`,
    :attr:`~orpheus.numerics.frame.FrameBase.reconstruction`,
    :meth:`~orpheus.numerics.frame.FrameBase.conjugate`) are inherited unchanged.
    """

    @classmethod
    def from_galerkin(cls, frame: GalerkinFrame) -> "HarmonicFrame":
        r"""Upgrade a generic angular :class:`GalerkinFrame` to a typed
        :class:`HarmonicFrame`, reusing its basis + measure (no rebuild).

        ``frame`` is the ``quadrature.angular_frame(L)`` Galerkin frame; the
        returned :class:`HarmonicFrame` shares its
        :class:`~orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis`
        trial basis and its angular :class:`~orpheus.numerics.measure.DiscreteMeasure`,
        so the table / spaces / faces are bit-identical — only the typed verbs
        are added.
        """
        return cls(frame.basis, frame.measure)

    # ── M : project (role-preserving, angular → moment) ──────────────

    @singledispatchmethod
    def analyse(self, field):
        r"""Project a per-ordinate carrier to its spherical-harmonic moments
        (:math:`M`), preserving the role.

        Dispatches on the carrier type: :class:`AngularFlux` →
        :class:`HarmonicMomentFlux`, :class:`AngularSourceSink` →
        :class:`HarmonicMomentSourceSink`. The output carries the input's
        ``mesh`` and the frame's truncation order ``L``.
        """
        raise TypeError(
            f"HarmonicFrame.analyse: unsupported carrier "
            f"{type(field).__name__}; expected AngularFlux or "
            f"AngularSourceSink (a per-ordinate carrier)."
        )

    @analyse.register
    def _(self, psi: AngularFlux) -> HarmonicMomentFlux:
        return HarmonicMomentFlux.from_mesh_and_L(
            self.analysis.apply(psi.values), psi.mesh, self.basis.L,
            spatial_moments=psi.spatial_moments_per_axis,
        )

    @analyse.register
    def _(self, q: AngularSourceSink) -> HarmonicMomentSourceSink:
        return HarmonicMomentSourceSink.from_mesh_and_L(
            self.analysis.apply(q.values), q.mesh, self.basis.L,
            spatial_moments=q.spatial_moments_per_axis,
        )

    # ── R : reconstruct (role-preserving, moment → angular) ──────────

    @singledispatchmethod
    def reconstruct(self, moment):
        r"""Reconstruct a per-ordinate carrier from its spherical-harmonic
        moments (:math:`R`), preserving the role.

        Dispatches on the carrier type: :class:`HarmonicMomentFlux` →
        :class:`AngularFlux`, :class:`HarmonicMomentSourceSink` →
        :class:`AngularSourceSink`. The output carries the input's ``mesh``.
        """
        raise TypeError(
            f"HarmonicFrame.reconstruct: unsupported carrier "
            f"{type(moment).__name__}; expected HarmonicMomentFlux or "
            f"HarmonicMomentSourceSink (a moment carrier)."
        )

    @reconstruct.register
    def _(self, phi: HarmonicMomentFlux) -> AngularFlux:
        return AngularFlux.from_mesh(
            self.reconstruction.apply(phi.values), phi.mesh,
            spatial_moments=phi.spatial_moments_per_axis,
        )

    @reconstruct.register
    def _(self, q: HarmonicMomentSourceSink) -> AngularSourceSink:
        return AngularSourceSink.from_mesh(
            self.reconstruction.apply(q.values), q.mesh,
            spatial_moments=q.spatial_moments_per_axis,
        )
