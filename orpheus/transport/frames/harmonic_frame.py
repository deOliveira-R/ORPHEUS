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
them *castable* — the carrier classes themselves — lives in the transport
field layer. And ``Quadrature.angular_frame`` (which builds the generic
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

The typed lift is BOUND at construction (the S4-amendment)
==========================================================

The typed verbs realise the lift :math:`M \otimes I_{\rm cells}` — a pair of
formal operators between two full field spaces — and an operator is not an
operator without its domain and codomain. The frame therefore binds BOTH at
construction: the constructor takes the per-ordinate **angular field space**
(:attr:`angular_space`) and derives + stores the **moment space**
(:attr:`moment_space`) once. The derivation direction is a constructor-input
fact: the moment space is a function of the angular space and the truncation
order ``L`` — never the reverse, because the angular space carries the
quadrature axis and, on a widened iterate, the scheme's mass-bearing
spatial-moment axis, which no moment operand determines. (An earlier revision
read this as a "structural asymmetry" of the frame square and made
``reconstruct`` take a per-call ``space=`` — that parameter was an unbound
operator's missing codomain leaking into the apply signature; the user
diagnosed it and the binding repaired it, 2026-08-22.)

Both verbs ADMIT by content equality against the bound spaces and their
outputs ride the bound spaces; the truncation order ``L`` comes from the
frame's own basis, and the within-cell spatial-moment width is read off the
bound angular space (#240 D5b-S3, hoisted to
:meth:`~orpheus.transport.fields._bases.BulkField._spatial_moments_per_axis_of`).
"""

from __future__ import annotations

from dataclasses import dataclass
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


def _derive_moment_space(
    angular_space: "FunctionSpace", L: int,
) -> "FunctionSpace":
    r"""The moment codomain derived from the bound angular domain (+ ``L``).

    ``SH.from_L(L) * of_axes(<non-angular, non-moment axes>)`` — the cell
    group is the angular space's own energy/spatial axes (the same instances
    the carrier's mints share, so the product content-equals
    ``MomentField._space_for_mesh_and_L``'s) — with the densified
    ``SpatialMomentSpace`` factor appended for a widened angular space,
    mirroring ``_compose_spatial_moments``'s densified arm. Runs ONCE, at
    frame construction (the S4-amendment): the derivation direction is
    moment = f(angular, L), never the reverse.
    """
    from orpheus.numerics.moment_layout import SPATIAL_MOMENT_AXIS_LABEL
    from orpheus.numerics.space import FunctionSpace
    from orpheus.numerics.spaces.spatial_moment_space import (
        SpatialMomentSpace,
    )
    from orpheus.numerics.spaces.spherical_harmonic_space import (
        SphericalHarmonicSpace,
    )
    from orpheus.transport.fields._bases import BulkField

    axes = angular_space.axes
    if axes is None:
        raise TypeError(
            "HarmonicFrame: the bound angular space must be axis-built (an "
            "S2 per-ordinate space); a shape-only space cannot name the "
            "cell axes the moment codomain is derived from."
        )
    cell_axes = [
        ax for ax in axes[1:] if ax.label != SPATIAL_MOMENT_AXIS_LABEL
    ]
    base = SphericalHarmonicSpace.from_L(L) * (
        FunctionSpace.of_axes(*cell_axes)
    )
    per_axis = BulkField._spatial_moments_per_axis_of(angular_space)
    if per_axis == 1:
        return base
    ndim = next(
        len(ax.shape) for ax in cell_axes if ax.label == "spatial"
    )
    return base * SpatialMomentSpace.from_per_axis(per_axis, ndim)


@dataclass(frozen=True, init=False)
class HarmonicFrame(GalerkinFrame):
    r"""The angular spherical-harmonic :class:`GalerkinFrame` with typed verbs,
    BOUND to its two full field spaces.

    Constructed from a generic ``quadrature.angular_frame(L)`` via
    :meth:`from_galerkin` (or directly from ``(basis, measure)`` — the
    constructor narrows the inherited signature to require a
    :class:`~orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis`:
    a harmonic frame over any other trial basis is an illegal state, and the
    truncation order ``L`` the typed verbs read exists only on the SH basis)
    — plus the **angular field space** the typed lift acts on. The moment
    space is derived and stored at construction (:func:`_derive_moment_space`;
    moment = f(angular, L), never the reverse). Adds the role-polymorphic
    :meth:`analyse` (:math:`M`, ``angular_space → moment_space``) and
    :meth:`reconstruct` (:math:`R`, ``moment_space → angular_space``) carrier
    verbs; the generic ``np.ndarray`` faces
    (:attr:`~orpheus.numerics.frame.FrameBase.analysis`,
    :attr:`~orpheus.numerics.frame.FrameBase.reconstruction`,
    :meth:`~orpheus.numerics.frame.FrameBase.conjugate`) are inherited
    unchanged and remain bound at the numerics level (basis space / measure
    space).
    """

    # Covariant narrowing of the inherited frozen (read-only) field: the
    # constructor below guarantees it, and the typed verbs read ``basis.L``.
    basis: SphericalHarmonicBasis

    #: The per-ordinate angular field space — :meth:`analyse`'s domain and
    #: :meth:`reconstruct`'s codomain. In production this is the posed
    #: composite's iterate-width ``interior_space``
    #: (:attr:`~orpheus.transport.operators.scattering.ScatteringOperator.frame`).
    angular_space: "FunctionSpace"

    #: The spherical-harmonic moment space — :meth:`analyse`'s codomain and
    #: :meth:`reconstruct`'s domain. DERIVED from :attr:`angular_space` +
    #: the basis ``L`` at construction; never a constructor input.
    moment_space: "FunctionSpace"

    def __init__(
        self,
        basis: SphericalHarmonicBasis,
        measure: DiscreteMeasure,
        *,
        angular_space: "FunctionSpace",
    ) -> None:
        super().__init__(basis, measure)
        object.__setattr__(self, "angular_space", angular_space)
        object.__setattr__(
            self, "moment_space", _derive_moment_space(angular_space, basis.L),
        )

    @classmethod
    def from_galerkin(
        cls, frame: GalerkinFrame, *, angular_space: "FunctionSpace",
    ) -> "HarmonicFrame":
        r"""Upgrade a generic angular :class:`GalerkinFrame` to a typed, bound
        :class:`HarmonicFrame`, reusing its basis + measure (no rebuild).

        ``frame`` is the ``quadrature.angular_frame(L)`` Galerkin frame; the
        returned :class:`HarmonicFrame` shares its
        :class:`~orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis`
        trial basis and its angular :class:`~orpheus.numerics.measure.DiscreteMeasure`,
        so the table / numerics spaces / faces are bit-identical — only the
        typed verbs and the field-space binding (``angular_space`` in, moment
        space derived) are added. A frame over any other trial basis (e.g. an
        indicator basis) is rejected here, at the upgrade boundary — not
        later, when a typed verb first reads the SH-only truncation order
        ``L``.
        """
        basis = frame.basis
        if not isinstance(basis, SphericalHarmonicBasis):
            raise TypeError(
                f"HarmonicFrame.from_galerkin requires a spherical-harmonic "
                f"trial basis; got {type(basis).__name__}."
            )
        return cls(basis, frame.measure, angular_space=angular_space)

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
        :class:`AngularSourceSink` → :class:`HarmonicMomentSourceSink`. A
        bound operator's verb: the operand is ADMITTED by content equality
        against the bound :attr:`angular_space`, and the output rides the
        bound :attr:`moment_space` (derived once, at construction — the
        S4-amendment; no per-operand space derivation, no per-call
        ``space=``).
        """
        if isinstance(field, (AngularFlux, AngularSourceSink)):
            if field.space != self.angular_space:
                raise TypeError(
                    f"HarmonicFrame.analyse: operand rides space "
                    f"{field.space.name!r} but this frame is bound to "
                    f"angular domain {self.angular_space.name!r} — a bound "
                    f"frame verb admits only elements of its bound spaces."
                )
            values = self.analysis.apply(field.values)
            per_axis = field.spatial_moments_per_axis
            if isinstance(field, AngularFlux):
                return HarmonicMomentFlux(
                    values=values, space=self.moment_space, L=self.basis.L,
                    spatial_moments=per_axis,
                )
            return HarmonicMomentSourceSink(
                values=values, space=self.moment_space, L=self.basis.L,
                spatial_moments=per_axis,
            )
        raise TypeError(
            f"HarmonicFrame.analyse: unsupported carrier "
            f"{type(field).__name__}; expected AngularFlux or "
            f"AngularSourceSink (a per-ordinate carrier)."
        )

    # ── R : reconstruct (role-preserving, moment → angular) ──────────

    @overload
    def reconstruct(self, moment: HarmonicMomentFlux) -> AngularFlux: ...

    @overload
    def reconstruct(
        self, moment: HarmonicMomentSourceSink,
    ) -> AngularSourceSink: ...

    def reconstruct(
        self, moment: HarmonicMomentFlux | HarmonicMomentSourceSink,
    ) -> AngularFlux | AngularSourceSink:
        r"""Reconstruct a per-ordinate carrier from its spherical-harmonic
        moments (:math:`R`), preserving the role.

        :class:`HarmonicMomentFlux` → :class:`AngularFlux`,
        :class:`HarmonicMomentSourceSink` → :class:`AngularSourceSink`. A
        bound operator's verb: the operand is ADMITTED by content equality
        against the bound :attr:`moment_space`, and the output rides the
        bound :attr:`angular_space`. (The pre-amendment per-call ``space=``
        parameter was an unbound operator's missing codomain — the angular
        target is carrier knowledge the CONSTRUCTOR takes, once, from the
        caller that poses the problem; see the module docstring.)
        """
        if isinstance(moment, (HarmonicMomentFlux, HarmonicMomentSourceSink)):
            if moment.space != self.moment_space:
                raise TypeError(
                    f"HarmonicFrame.reconstruct: operand rides space "
                    f"{moment.space.name!r} but this frame is bound to "
                    f"moment domain {self.moment_space.name!r} — a bound "
                    f"frame verb admits only elements of its bound spaces."
                )
            values = self.reconstruction.apply(moment.values)
            if isinstance(moment, HarmonicMomentFlux):
                return AngularFlux(values=values, space=self.angular_space)
            return AngularSourceSink(values=values, space=self.angular_space)
        raise TypeError(
            f"HarmonicFrame.reconstruct: unsupported carrier "
            f"{type(moment).__name__}; expected HarmonicMomentFlux or "
            f"HarmonicMomentSourceSink (a moment carrier)."
        )
