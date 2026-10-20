r"""Angular windowing as a typed composition — ``P @ A.inverse()``
(#226 taxonomy step 2, §17 W1).

The Phase 5a/5c "moment-windowed resolvent" was a duck-typed wrapper whose
``solve`` returned harmonic MOMENTS — an output-mode config that silently
changed the operator's codomain.  A codomain change is not a config of one
morphism; it is a DIFFERENT morphism, a composition (§1's two-layer law).
This module reifies the two factors:

.. math::

    \text{windowed} \;=\; P \circ A^{-1}, \qquad
    P \;=\; \underbrace{M_{\rm frame}}_{\text{analysis on the bulk}}
            \oplus \underbrace{\mathrm{Id}}_{\text{on the trace}},

where :math:`M_{\rm frame}` is the scattering operator's MINTED flux-analysis
face (F-1 — ``S.flux_analysis``, an ``AngularFlux → HarmonicMomentFlux``
bound :class:`~orpheus.numerics.projection.AnalysisOperator` member whose
kernel is the frame's angular→moment reduction ``φ_ℓ^m = Σ_n w_n Y_ℓ^m ψ_n``)
and the trace passes through un-reduced (windowing is interior-bulk-only —
the reflective coupling needs the full per-ordinate face trace).  :math:`P` is a **block coisometry**
(``analysis ∘ reconstruction = 4π·I`` under the no-prefactor SH convention —
pinned by ``test_pi_R_is_4pi_identity_through_the_frame``), NOT invertible —
so the composite honestly makes no round-trip promise (``is_invertible`` is
``False``, by the product closure laws).

**Fusion is an evaluation strategy, not a new operator.**
:class:`WindowedSweep` is an atom whose :meth:`~WindowedSweep.expansion` is the
plain word ``P @ A⁻¹`` and whose :meth:`~WindowedSweep.apply` DEFORESTS the
intermediate: instead of materializing the full per-ordinate angular field and
then projecting (``P.apply(A_inv.apply(rhs))`` — the expansion's own
``apply``, the fuller view the windowed-≡-post-projection gate spells
independently), it runs the sweep in the
substrate's MOMENT-emit mode (``_SweepEmit`` accumulating per anti-diagonal),
so the ``(N, ng, nx, ny)`` intermediate is never built (the ~3× linear
peak-memory win).  The two evaluations differ only in the ordinate-sum
reduction order ⇒ principled-equivalence within the scale-relative
``4·N·eps`` bound (the windowed-≡-post-projection gate).

The fused object is built EXPLICITLY, ``WindowedSweep(P, A.inverse())``, by
its one owner, the solver's ``_maybe_window`` factory; ``P @ A.inverse()`` is
the plain word, because ``@`` never fuses (numerics_untangling.md ruling 52:
a fusion chosen by a type test inside a dunder is a realisation chosen by a
spelling, and the realisation is the Strategy's). The former public
``solve_moments`` cross-reach (three sites) retired into this one
composition.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Final, Optional, cast

from orpheus.numerics.operator import FusedAtom, LinearOperator

from .sweep_operator import SweepOperator

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace
    from orpheus.numerics.spaces.full_field_space import FullFieldSpace
    from orpheus.transport.fields.angular_flux import AngularFlux
    from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
    from orpheus.transport.frames import HarmonicAnalysisOperator
    from orpheus.transport.full_field import FullField
    from orpheus.transport.timed_full_field import TimedFullField


__all__ = ["BulkAnalysisOperator", "WindowedSweep"]


class BulkAnalysisOperator(LinearOperator["FullField", "TimedFullField"]):
    r"""The minted flux-analysis face on the BULK ⊕ identity on the trace —
    the spelled word :math:`P = \mathrm{lift}(M) \oplus \mathrm{Id}_{\Gamma}`.

    Maps the composite carrier's per-ordinate bulk to its harmonic moments
    through the MINTED face (``AngularFlux → HarmonicMomentFlux`` — a bound
    :class:`~orpheus.numerics.projection.AnalysisOperator` member) while the
    boundary trace passes through VERBATIM (the biproduct keeps the trace a
    distinct summand — the reflective ``B`` coupling reads full per-ordinate
    face values).  A block coisometry, not invertible. Its transpose is the
    block transpose :math:`P^{\mathsf T} = M^{\mathsf T} \oplus \mathrm{Id}_{\Gamma}`:
    the face's own transpose on the moment bulk, the trace passed through.
    Until 2026-10-18 it advertised no transpose although the face has one
    (``tests/gates/numerics/test_adjointability_is_honest.py``).

    ``face`` must be the SCATTERING operator's own minted
    :attr:`~orpheus.transport.operators.scattering.ScatteringOperator.flux_analysis`
    (F-1): its kernel IS ``S.apply``'s internal projection — the same frame,
    the same cached table, single source — so the emitted moments equal what
    ``S`` would project term-for-term, and the bulk admission (the
    space-content invariant) is the face's own bound-domain admission.
    """

    def __init__(
        self,
        face: "HarmonicAnalysisOperator[AngularFlux, HarmonicMomentFlux]",
        full_field_space: "FullFieldSpace",
    ) -> None:
        from orpheus.numerics.spaces.full_field_space import FullFieldSpace

        #: The minted flux-analysis face whose kernel reduces the bulk.
        self.face = face
        #: The composite carrier (bulk ⊕ trace) this operator acts on —
        #: the one thing the mesh used to contribute (un-weld arc O-1).
        self._domain = full_field_space
        # F-1: the moment-bulk composite is posed from the face's bound
        # codomain ⊕ the angular composite's trace block (the trace passes
        # through untouched). The pre-F-1 ``None`` ("not yet typed") debt
        # died with the mint — composition guards now check this end.
        assert full_field_space.trace_space is not None  # SN composite; narrowing only
        self._codomain = FullFieldSpace.from_blocks(
            self.face.codomain, full_field_space.trace_space,
        )

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        return self._domain

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        return self._codomain

    def apply(self, field: "FullField") -> "TimedFullField":
        r"""``P·ψ`` — face-project the bulk, pass the trace through.

        The bulk admission is the FACE's (exactly ``AngularFlux`` on its
        bound domain — the space-content invariant, one spelling); the wrap
        rides the face's bound moment codomain (no per-apply re-mint).
        """
        from orpheus.transport.fields.angular_flux import AngularFlux
        from orpheus.transport.timed_full_field import TimedFullField

        interior = field.interior
        assert isinstance(interior, AngularFlux)  # narrowing only — the
        # REAL admission (carrier + bound space) is the face's, next line.
        bulk = self.face.apply(interior)
        return TimedFullField(
            interior=bulk,
            boundary=field.boundary,
            _history=(),
            history_depth=(
                field.history_depth
                if isinstance(field, TimedFullField)
                else 0
            ),
        )

    def apply_transpose(self, field: "FullField") -> "FullField":
        r"""``Pᵀ·φ = Mᵀ φ_bulk ⊕ φ_Γ`` — the face's transpose on the moment bulk, the trace passed through."""
        from orpheus.transport.fields.angular_flux import AngularFlux
        from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
        from orpheus.transport.full_field import FullField

        # The admission is the face's, as in :meth:`apply`; a typed moment
        # maps to a typed flux.
        bulk = self.face.apply_transpose(cast(HarmonicMomentFlux, field.interior))
        return FullField(interior=cast(AngularFlux, bulk), boundary=field.boundary)

    @property
    def is_adjointable(self) -> bool:
        # The block transpose exists wherever the face's does (the identity
        # on the trace always transposes).
        return self.face.is_adjointable

    def __repr__(self) -> str:
        return f"BulkAnalysisOperator({self.face!r})"


class WindowedSweep(FusedAtom["FullField", "TimedFullField"]):
    r"""The fused evaluation of ``P @ A.inverse()`` — same morphism, deforested.

    An atom with two named legs, :attr:`p` and :attr:`sweep`, whose algebraic
    content is :meth:`expansion`, the plain word ``P @ A⁻¹`` (ruling 52).
    :meth:`apply` is the substrate moment emit: the sweep accumulates
    ``w·Y·ψ`` per anti-diagonal (``_SweepEmit`` MOMENT mode), so the full
    per-ordinate intermediate is never materialized. The expansion's own
    ``apply`` — ``P.apply(A⁻¹.apply(rhs))`` — is the fuller view; the
    windowed-≡-post-projection gate
    (``tests/gates/sn/solve/test_2d_anisotropic_windowing.py``) spells that
    oracle independently, as the flat analysis of the full-angular sweep. The spaces, the transpose and the
    honest ``is_invertible = False`` (the coisometry factor) are the
    expansion's.

    2-D Cartesian only, by construction: the moment emit is a 2-D kernel and
    the solver's ``_maybe_window`` gate (``is_cartesian and ndim == 2``) is
    the single production factory.
    """

    def __init__(
        self, p: "BulkAnalysisOperator", sweep: "SweepOperator",
    ) -> None:
        if not isinstance(p, BulkAnalysisOperator):
            raise TypeError(
                f"WindowedSweep: left factor must be a BulkAnalysisOperator; "
                f"got {type(p).__name__}."
            )
        if not isinstance(sweep, SweepOperator):
            raise TypeError(
                f"WindowedSweep: right factor must be a SweepOperator; "
                f"got {type(sweep).__name__}."
            )
        if p.domain != sweep.inner.domain:
            raise ValueError(
                "WindowedSweep: the analysis factor and the sweep factor "
                "must act on the same composite space (space-content "
                "invariant)."
            )
        #: The analysis factor ``P``.
        self.p: Final = p
        #: The inverse factor ``A⁻¹``.
        self.sweep: Final = sweep

    def expansion(self) -> "LinearOperator[FullField, TimedFullField]":
        r"""The plain word ``P @ A⁻¹`` this atom evaluates fused (ruling 52)."""
        return self.p @ self.sweep

    def apply(
        self, rhs: "FullField", *,
        initial_guess: "FullField | None" = None,
    ) -> "TimedFullField":
        r"""``(P ∘ A⁻¹)·rhs`` via the substrate moment emit (fused).

        Routes to the forward operator's uniform private solve surface with
        ``moment_frame`` set — the ONE moment-emitting entry (the former
        ``solve_moments`` cross-reach retired here).  Principled-equivalent
        (≤ scale-relative ``4·N·eps``) to the expansion's ``apply``, the
        fuller view the windowed-≡-post-projection gate spells independently.

        ``initial_guess`` is the inverse family's canonical driver
        signature, accepted and unused here: the multi-D Cartesian
        wavefront walk has no bulk-seed consumer (the Carlson coupled-pole
        seed is 1-D curvilinear), and the previous iterate's boundary lag
        rides the driver's ``B`` / ``B_upper`` gain — never this kwarg.
        The windowed moment iterate could not seed the angular walk anyway
        (wrong representation); accepting-and-ignoring is the honest
        contract, stated here rather than behind an adapter.
        """
        del initial_guess  # no bulk-seed consumer in the multi-D walk
        # F-1: the substrate's moment emit consumes the FRAME (its analysis
        # kernel — table + weights); the minted face carries it, one hop in.
        return self.sweep.inner._solve_timed_full_field(
            rhs, moment_frame=self.p.face.frame,
        )

    def __repr__(self) -> str:
        return f"WindowedSweep({self.p!r} @ {self.sweep!r})"
