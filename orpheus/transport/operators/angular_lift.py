r"""``AngularLift`` — an ENERGY binding lifted onto the angular composite by
the frame's :math:`\ell = 0` conjugation, acting through the body its ends
select.

Every collision gain of the within-group algebra

.. math::

    A \;=\; L + C - S - N_{2n} - B, \qquad A\,\psi = \tfrac{1}{k}\,F\,\psi

is the same shape: an energy operator :math:`E` on the scalar flux
(:math:`y\,\Sigma_{c,0}^{\mathsf T}` for a transfer channel, the dyad
:math:`|\chi\rangle\langle\nu\Sigma_f|` for fission), lifted onto the
per-ordinate composite by the frame's :math:`\ell = 0` conjugation
:math:`R_0\,E\,M_0 / W` — realised on the reaction-rate fast path
(:math:`\phi = \int\psi\,d\Omega`, then :math:`E\phi`, then the
producer-side :math:`/W` broadcast; no moment tensor on the hot loop) —
plus, for a transfer channel, the :math:`\ell \ge 1` redistribution
:math:`R\,\Lambda_{\ell\ge1}\,M / W`. This class is that lift, ONCE
(CS4c step 5, ruling R-1: ``{S, N₂ₙ} | {F}`` share the ℓ = 0 base and
differ by the datum, the energy binding it derives, and whether an
:math:`\ell \ge 1` part exists — never by a second spelling of the lift).

**Each binding acts through the body its ends select** (the step-5
outcome). The retained analysis face :math:`M \otimes I` has two ends —
the per-ordinate space it reads and the moment space it writes — and the
binding's DOMAIN interior is one of them:

* ``domain.interior == flux_analysis.domain`` — the **angular** end: the
  operand is the per-ordinate flux; :math:`\phi` is its angular
  integral; a coefficient-space operator acts through
  ``frame.conjugate`` (:math:`R\,\Lambda\,M`); the transpose's cotangent
  is a per-ordinate source/sink;
* ``domain.interior == flux_analysis.codomain`` — the **moment** end:
  the operand is ALREADY :math:`M\psi` (the 2-D Cartesian windowed
  iterate, which never materialises the per-ordinate flux between
  sweeps); :math:`\phi` is its :math:`\ell = 0` slot; a coefficient-space
  operator acts through ``frame.reconstruct_after`` (:math:`R\,\Lambda`,
  :math:`M` skipped — re-projecting would double-project); the cotangent
  is a moment source/sink.

The selection happens ONCE, at construction; a third interior is refused.
Until step 5 the moment operand was handed to the ANGULAR-bound operator,
which dispatched on the carrier's class per call (`[M]` 2026-09-04: 143
such feeds per windowed solve, on a bit-exact frozen snapshot) — the
shipped non-endomorphism the step-0 census measured. Now the windowed
driver binds the gains on the moment composite (:meth:`on_moment_domain`)
and every operand rides its operator's own domain
(:func:`~orpheus.transport.operators.lift.admit_composite`).

Why the scalar sub-space is read off the CODOMAIN (F-1 of the step-5
verification plan): the moment composite's interior is a tensor-product
space with no axes, so it cannot name its own energy ⊗ spatial factor;
the codomain is the angular composite for BOTH ends and its interior is
axis-built. The energy binding therefore lives on the codomain's scalar
sub-space — which is also where the emission lands.

The subclass contract is three members and no verbs:

* :attr:`data_ng` — the bound datum's group count (the per-end energy
  admission's operand);
* :meth:`_bind_energy` — the datum's ENERGY binding on a given scalar
  space (the ℓ = 0 middle factor; cached once as :attr:`isotropic_energy`);
* :meth:`_frame_form` — the WHOLE action as one frame-conjugated product
  (the transpose's spelling: factor reversal, no arithmetic of its own).

A subclass with an :math:`\ell \ge 1` part overrides
:meth:`_interior_action` to add it (the transfer core does; fission does
not — :math:`\chi` carries no angle).
"""

from __future__ import annotations

from abc import abstractmethod
from dataclasses import dataclass, field, replace
from functools import cached_property
from typing import TYPE_CHECKING, Generic, Self, TypeVar, cast

import numpy as np

from orpheus.numerics.operator import BlockRole
from orpheus.numerics.space import FunctionSpace
from orpheus.numerics.spaces.full_field_space import FullFieldSpace
from orpheus.transport.fields._bases import BulkField, FieldRole
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.operators.bound_operator import BoundOperator
from orpheus.transport.operators.lift import (
    admit_composite,
    interior_space_of,
    lift_bulk_action,
)
from orpheus.transport.source_sinks import (
    AngularSourceSink,
    HarmonicMomentSourceSink,
    ScalarSourceSink,
)

if TYPE_CHECKING:
    from orpheus.numerics.operator import LinearOperator, OperatorProduct
    from orpheus.transport.frames.harmonic_frame import (
        HarmonicAnalysisOperator,
        HarmonicFrame,
        HarmonicReconstructionOperator,
    )

__all__ = ["AngularLift"]

#: The ENERGY binding a lift derives from its datum (the ℓ = 0 middle factor).
EnergyT = TypeVar("EnergyT", bound=BoundOperator)


# ═══════════════════════════════════════════════════════════════════════
# The two ends of the analysis face — what the operand IS
# ═══════════════════════════════════════════════════════════════════════


class _AngularEnd:
    r"""The domain's interior is the analysis face's DOMAIN: the operand is
    the per-ordinate flux :math:`\psi`."""

    @staticmethod
    def scalar_flux(lift: "AngularLift", bulk: BulkField) -> ScalarFlux:
        # φ = ∫ψ dΩ — the ONE reduction body, on the space's memoised
        # retraction (bit-identical to the pre-step-5 arm).
        return cast(AngularFlux, bulk).integrate_angular()

    @staticmethod
    def conjugate(lift: "AngularLift", operator: "LinearOperator") -> "OperatorProduct":
        # R ∘ A ∘ M — the frame's production composition.
        return lift.frame.conjugate(operator)

    @staticmethod
    def cotangent(lift: "AngularLift", values: np.ndarray) -> BulkField:
        return AngularSourceSink(values=values, space=lift._domain_interior)


class _MomentEnd:
    r"""The domain's interior is the analysis face's CODOMAIN: the operand
    is already :math:`M\psi` (the windowed moment iterate)."""

    @staticmethod
    def scalar_flux(lift: "AngularLift", bulk: BulkField) -> ScalarFlux:
        # The ℓ=0 moment IS the scalar flux (Y_0^0 = 1) — the typed
        # accessor carries the convention; the target is the codomain's
        # scalar sub-space (a widened iterate cannot self-derive it).
        return cast(HarmonicMomentFlux, bulk).scalar_flux(
            space=lift._scalar_interior_space,
        )

    @staticmethod
    def conjugate(lift: "AngularLift", operator: "LinearOperator") -> "OperatorProduct":
        # R ∘ A — M already applied; conjugating would double-project.
        return cast("OperatorProduct", lift.frame.reconstruct_after(operator))

    @staticmethod
    def cotangent(lift: "AngularLift", values: np.ndarray) -> BulkField:
        domain_interior = lift._domain_interior
        return HarmonicMomentSourceSink(
            values=values,
            space=domain_interior,
            L=lift.frame.truncation_order,
            spatial_moments=BulkField._spatial_moments_per_axis_of(domain_interior),
        )


# ═══════════════════════════════════════════════════════════════════════
# The lift
# ═══════════════════════════════════════════════════════════════════════


@dataclass(eq=False)
class AngularLift(BoundOperator["FullField", "FullField"], Generic[EnergyT]):
    r"""An energy binding lifted onto the angular composite — see the module
    docstring.

    Parameters
    ----------
    flux_analysis : HarmonicAnalysisOperator
        The minted FLUX analysis face :math:`M \otimes I`
        (``AngularFlux → HarmonicMomentFlux``) bound on the angular space
        this binding EMITS on (the codomain's interior). Its two ends are
        the two admissible domain interiors; its frame is the binding's
        (:attr:`frame` rides on it — provenance, zero extra state).
    source_reconstruction : HarmonicReconstructionOperator
        The minted SOURCE reconstruction face :math:`R \otimes I`
        (``HarmonicMomentSourceSink → AngularSourceSink``) landing on the
        same angular space — the typed :math:`R` of the moment end's
        :math:`\ell \ge 1` route (the transfer core's).
    domain, codomain : FunctionSpace
        The two mandatory ends (kw-only, write-once —
        :class:`~orpheus.transport.operators.bound_operator.BoundOperator`):
        composite full-field spaces. The codomain's interior is the angular
        space the faces are bound on; the domain's interior is one of the
        analysis face's two ends and SELECTS the body.
    """

    flux_analysis: "HarmonicAnalysisOperator[AngularFlux, HarmonicMomentFlux]" = field(
        kw_only=True,
    )
    source_reconstruction: "HarmonicReconstructionOperator[HarmonicMomentSourceSink, AngularSourceSink]" = field(
        kw_only=True,
    )

    # A collision gain is volumetric — bulk only, no face-trace action; the
    # lift enters the composite by extension-by-zero on the trace. A
    # class-level default of the base's ``block_role`` instance attribute,
    # deliberately unannotated (a ClassVar annotation would override the
    # base's instance variable; a plain annotation would make it a field).
    block_role = BlockRole.BULK

    #: The selected end — set once in :meth:`__post_init__` (not a ctor
    #: argument: it is DERIVED from the ends, never declared).
    _end: "type[_AngularEnd] | type[_MomentEnd]" = field(init=False, repr=False)

    # ── the subclass contract ────────────────────────────────────────

    @property
    @abstractmethod
    def data_ng(self) -> int:
        """The bound datum's group count (the per-end energy admission's operand)."""

    @abstractmethod
    def _bind_energy(self, scalar_space: FunctionSpace) -> EnergyT:
        """The datum's ENERGY binding, endomorphic on ``scalar_space``."""

    @abstractmethod
    def _frame_form(self) -> "OperatorProduct":
        r"""The whole action as ONE frame-conjugated product (pre-:math:`/W`)
        — the transpose's spelling."""

    # ── construction: admission, then the SELECTION ──────────────────

    def __post_init__(self) -> None:
        owner = type(self).__name__
        # Per-END energy admission (CS4c §1). On an axes-less composite
        # the guard is declaredly inert (its own docstring) — the moment
        # composite's domain leg is one such, exactly as every shipped
        # composite binding's already is.
        self._assert_energy_extent_both_ends(self.data_ng, operator=owner)
        # Face admission against the angular space this binding EMITS on:
        # both faces are minted on the codomain's interior.
        emitted = self._codomain_interior
        if self.flux_analysis.domain != emitted:
            raise TypeError(
                f"{owner}: the flux-analysis face is bound to a different "
                f"angular space than this binding's interior — mint the "
                f"faces from the SAME posed space (tier 2 does)."
            )
        if self.source_reconstruction.codomain != emitted:
            raise TypeError(
                f"{owner}: the source-reconstruction face lands on a "
                f"different angular space than this binding's interior — "
                f"mint the faces from the SAME posed space (tier 2 does)."
            )
        # The selection: which END of the analysis face the domain's
        # interior is decides the body, once.
        consumed = self._domain_interior
        if consumed == self.flux_analysis.domain:
            self._end = _AngularEnd
        elif consumed == self.flux_analysis.codomain:
            self._end = _MomentEnd
        else:
            raise TypeError(
                f"{owner}: the domain's interior {consumed!r} is neither "
                f"end of the analysis face (angular "
                f"{self.flux_analysis.domain!r} / moment "
                f"{self.flux_analysis.codomain!r}) — a binding acts through "
                f"the body its ends select, and no body reads this space."
            )
        # Bind the energy operator EAGERLY: its plain-scalar admission on
        # the axis-built scalar sub-space is the effective ng guard (the
        # composite legs above are declaredly inert on an axes-less end),
        # and the F-1 requirement — the emitted interior names a scalar
        # sub-space — is checked at construction, not at first apply.
        self.isotropic_energy

    # ── derived structure (single sources) ───────────────────────────

    @property
    def frame(self) -> "HarmonicFrame":
        r"""The HUB-interned frame the faces were minted from, riding on
        :attr:`flux_analysis` (provenance, zero extra state) — the
        conjugation properties read it to compose the frame forms."""
        return self.flux_analysis.frame

    @property
    def total_weight(self) -> float:
        r""":math:`W = \int_{S^2} d\Omega` — the binding measure's total
        angular weight (the producer-side :math:`/W`), read off the
        frame's MEASURE (operative data)."""
        return float(np.asarray(self.frame.measure.weights).sum())

    @property
    def _domain_interior(self) -> FunctionSpace:
        """The domain's bulk block — the operand's space."""
        return interior_space_of(self.domain, owner=type(self).__name__)

    @property
    def _codomain_interior(self) -> FunctionSpace:
        """The codomain's bulk block — the angular space the lift emits on."""
        return interior_space_of(self.codomain, owner=type(self).__name__)

    @property
    def _moment_space(self) -> FunctionSpace:
        r"""The coefficient space of the bound frame's BASIS — the
        endomorphic ends of the internally-minted moment factors.

        READ off the frame (``frame.basis.space``), never minted from an
        order: which family spans the moments is the quadrature's
        decision. The continuum-metric space (the basis's own), not the
        frame's Parseval-dressed ``basis_space``: `[M]` #429 tracker 2.5
        the two are ``(name, shape)``-equal and metric-DIFFERENT on 33 of
        33 shipped (rule, L) rows, and under the continuum end the
        factor's Hilbert adjoint is its transpose EXACTLY
        (:math:`\Lambda^* = \Lambda^{\mathsf T}`, 0.0 on 33/33) while the
        dressed end would move it on 10 of 33."""
        return self.frame.basis.space

    @property
    def _scalar_interior_space(self) -> FunctionSpace:
        r"""The emitted interior's scalar ``(ng, *spatial)`` sub-space —
        the energy binding's ends and the moment end's scalar-flux target,
        minted here once. Read off the CODOMAIN (axis-built for both
        ends — see the module docstring)."""
        interior = self._codomain_interior
        if interior.axes is None:
            raise TypeError(
                f"{type(self).__name__}: the composite interior must be "
                f"axis-built to name the scalar sub-space."
            )
        return FunctionSpace.of_axes(*interior.axes[1:])

    @cached_property
    def isotropic_energy(self) -> EnergyT:
        r"""The :math:`\ell = 0` ENERGY binding of this operator's own
        datum, on the emitted interior's scalar sub-space — the middle
        factor the fast path lifts, and what the scalar consumers read
        (the solver's ``K_iso = S.isotropic_energy + N2N.isotropic_energy``
        at the ONE within-group construction site; the eigen-posing's
        ray seed). Cached at construction-time semantics (the datum is
        immutable)."""
        return self._bind_energy(self._scalar_interior_space)

    @property
    def is_adjointable(self) -> bool:
        # The Euclidean transpose is the frame form's factor reversal;
        # is_invertible inherits base False — an emission is not
        # invertible.
        return True

    # ── re-binding on the other end ──────────────────────────────────

    def on_moment_domain(self) -> Self:
        r"""This binding re-bound to CONSUME the moment representation —
        the same datum, the same faces, the same codomain; the domain's
        interior becomes the analysis face's codomain (:math:`M\psi`).

        The windowed SI driver's gain: the 2-D Cartesian iterate is held
        as harmonic moments, so the gains that read it are bound here
        (``S.on_moment_domain()``, ``N2N.on_moment_domain()``). A moment
        binding is not an endomorphism — its domain is the moment
        composite, its codomain the angular composite — and the
        ``OperatorSum`` guard never sees it: the windowed loop consumes
        the gains one by one. Built through :func:`dataclasses.replace`,
        so every admission re-runs and the selection lands on the moment
        end.
        """
        domain = self.domain
        if not isinstance(domain, FullFieldSpace) or domain.trace_space is None:
            raise TypeError(
                f"{type(self).__name__}.on_moment_domain: the bound domain "
                f"carries no trace block to pair with the moment interior."
            )
        return replace(
            self,
            domain=FullFieldSpace.from_blocks(
                self.flux_analysis.codomain, domain.trace_space,
            ),
        )

    # ── the action ───────────────────────────────────────────────────

    def apply(self, x: FullField, /) -> FullField:
        r"""The lifted emission :math:`T\psi` on the composite — the
        interior body selected at construction, the zero source/sink on
        the trace (a collision gain is volumetric)."""
        psi = admit_composite(self, x, end="domain")
        return lift_bulk_action(
            psi, self._interior_action, trace_role=FieldRole.SOURCE_SINK,
        )

    def _interior_action(self, bulk: BulkField) -> AngularSourceSink:
        r"""The bulk emission — the :math:`\ell = 0` lift alone here; a
        subclass with an :math:`\ell \ge 1` part overrides this to add it
        through :meth:`_combine`."""
        return self._combine(self._isotropic_source(bulk), None)

    def _isotropic_source(self, bulk: BulkField) -> ScalarSourceSink:
        r"""The :math:`\ell = 0` emission in iso scalar magnitude —
        :math:`E\,\phi` through the energy binding, riding the scalar
        flux's own space (the reaction-rate fast path: no moment tensor)."""
        phi = self._end.scalar_flux(self, bulk)
        return ScalarSourceSink(
            values=np.asarray(self.isotropic_energy.apply(phi.values)),
            space=phi.space,
        )

    def _combine(
        self, iso: ScalarSourceSink, aniso: AngularSourceSink | None,
    ) -> AngularSourceSink:
        r"""The producer-side combine :math:`(\text{iso}/W) + \text{aniso}`
        — the :math:`1/W` convention's ONE home (its normalisation chain:
        ``docs/theory/methods/sn/slab_multigroup.rst``). The zero
        accumulator of a purely isotropic lift rides the emitted
        interior; the containment dunder's cross-class arm returns the
        LARGER (angular) class (the #288 principled LSP exception)."""
        aniso_part = (
            aniso if aniso is not None
            else AngularSourceSink.zeros(self._codomain_interior)
        )
        return cast(
            AngularSourceSink, (iso / self.total_weight) + aniso_part,
        )

    def apply_transpose(self, x: FullField, /) -> FullField:
        r"""The Euclidean transpose :math:`T^{\mathsf T}\chi` on the
        composite — the frame form reversed by the
        :class:`~orpheus.numerics.operator.OperatorProduct` chain, then
        the producer :math:`/W`; the cotangent lands on the DOMAIN's
        interior in the end's own source/sink class; the trace is the
        implicit zero (a volumetric gain's transpose is volumetric).

        This is the Euclidean transpose (L12), not the metric Hilbert
        adjoint ``.H`` (which composes the spaces' Riesz legs around it).
        """
        chi = admit_composite(self, x, end="codomain")
        return lift_bulk_action(
            chi, self._interior_transpose, trace_role=FieldRole.SOURCE_SINK,
        )

    def _interior_transpose(self, bulk: BulkField) -> BulkField:
        values = (
            np.asarray(self._frame_form().apply_transpose(np.asarray(bulk.values)))
            / self.total_weight
        )
        return self._end.cotangent(self, values)
