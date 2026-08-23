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
:class:`WindowedSweep` is the composition itself (an
:class:`~orpheus.numerics.operator.OperatorProduct` subclass) whose
:meth:`~WindowedSweep.apply` DEFORESTS the intermediate: instead of
materializing the full per-ordinate angular field and then projecting
(``P.apply(A_inv.apply(rhs))`` — the inherited base body, retained verbatim
as the fuller-view verification oracle), it runs the sweep in the
substrate's MOMENT-emit mode (``_SweepEmit`` accumulating per anti-diagonal),
so the ``(N, ng, nx, ny)`` intermediate is never built (the ~3× linear
peak-memory win).  The two evaluations differ only in the ordinate-sum
reduction order ⇒ principled-equivalence within the scale-relative
``4·N·eps`` bound (the windowed-≡-post-projection gate).

The ``@`` spelling builds the fused object:
:meth:`BulkAnalysisOperator.__matmul__` dispatches on a
:class:`~orpheus.sn.operators.sweep_operator.SweepOperator` right factor —
mirroring the ``L + C`` fusion precedent — so ``P @ A.inverse()`` IS the
windowed path.  The solver's ``_maybe_window`` is the factory; the former
public ``solve_moments`` cross-reach (three sites) retired into this one
composition.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional, overload

from orpheus.numerics.operator import (
    D2,
    LinearOperator,
    OperatorProduct,
)

from .sweep_operator import SweepOperator

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace
    from orpheus.transport.fields.angular_flux import AngularFlux
    from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
    from orpheus.transport.frames import HarmonicAnalysisOperator
    from orpheus.transport.full_field import FullField
    from orpheus.transport.timed_full_field import TimedFullField
    from ..mesh.augmented_mesh import SNMesh


__all__ = ["BulkAnalysisOperator", "WindowedSweep"]


class BulkAnalysisOperator(LinearOperator["FullField", "TimedFullField"]):
    r"""The minted flux-analysis face on the BULK ⊕ identity on the trace —
    the spelled word :math:`P = \mathrm{lift}(M) \oplus \mathrm{Id}_{\Gamma}`.

    Maps the composite carrier's per-ordinate bulk to its harmonic moments
    through the MINTED face (``AngularFlux → HarmonicMomentFlux`` — a bound
    :class:`~orpheus.numerics.projection.AnalysisOperator` member) while the
    boundary trace passes through VERBATIM (the biproduct keeps the trace a
    distinct summand — the reflective ``B`` coupling reads full per-ordinate
    face values).  A block coisometry: not invertible, no transpose
    advertised (mint the synthesis-side lift when a consumer arrives), so
    apply-only, and the two-axis contract holds by the base defaults.

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
        sn_mesh: "SNMesh",
    ) -> None:
        from orpheus.numerics.spaces.full_field_space import FullFieldSpace

        #: The minted flux-analysis face whose kernel reduces the bulk.
        self.face = face
        #: The mesh carrying the composite-carrier geometry for the wrap.
        self.sn_mesh = sn_mesh
        # F-1: the moment-bulk composite is posed from the face's bound
        # codomain ⊕ the angular composite's trace block (the trace passes
        # through untouched). The pre-F-1 ``None`` ("not yet typed") debt
        # died with the mint — composition guards now check this end.
        full = sn_mesh.full_field_space
        assert full.trace_space is not None  # SN composite; narrowing only
        self._codomain = FullFieldSpace.from_blocks(
            self.face.codomain, full.trace_space,
        )

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        return self.sn_mesh.full_field_space

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

    @overload
    def __matmul__(self, other: "SweepOperator") -> "WindowedSweep": ...
    @overload
    def __matmul__(
        self, other: "LinearOperator[D2, FullField]",
    ) -> "OperatorProduct[D2, TimedFullField]": ...
    def __matmul__(self, other):
        r"""``P @ A.inverse()`` — the windowed composition, FUSED.

        Dispatches a :class:`SweepOperator` right factor to
        :class:`WindowedSweep` (whose ``apply`` evaluates through the
        substrate's moment emit); any other operand falls through to the
        generic :class:`~orpheus.numerics.operator.OperatorProduct`.
        The typed ``@overload`` spells the fusion rule (C4, as
        ``StreamingOperator.__add__``). Mirrors the ``L + C`` fusion
        precedent (#261: one-directional dispatch on the specific leaf).
        """
        if isinstance(other, SweepOperator):
            return WindowedSweep(self, other)
        return super().__matmul__(other)

    def __repr__(self) -> str:
        return f"BulkAnalysisOperator({self.face!r})"


class WindowedSweep(
    OperatorProduct[
        "FullField", "TimedFullField", BulkAnalysisOperator, SweepOperator,
    ],
):
    r"""The fused evaluation of ``P @ A.inverse()`` — same morphism, deforested.

    :meth:`apply` overrides the inherited factor-by-factor body with the
    substrate moment emit: the sweep accumulates ``w·Y·ψ`` per anti-diagonal
    (``_SweepEmit`` MOMENT mode), so the full per-ordinate intermediate is
    never materialized.  The inherited
    :meth:`OperatorProduct.apply <orpheus.numerics.operator.OperatorProduct.apply>`
    body — ``P.apply(A⁻¹.apply(rhs))`` — is the fuller-view verification
    ORACLE (the windowed-≡-post-projection gate calls it unbound), retained
    per the aggressive-retirement oracle exception.  Everything else
    (spaces, the honest ``is_invertible = False`` from the
    coisometry factor) is the generic product's.

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
        if p.sn_mesh is not sweep.inner.sn_mesh:
            raise ValueError(
                "WindowedSweep: the analysis factor and the sweep factor "
                "must act on the same mesh instance (mesh-identity "
                "invariant)."
            )
        super().__init__(p, sweep)

    @property
    def p(self) -> "BulkAnalysisOperator":
        """The analysis factor ``P`` (alias for ``self.a``)."""
        return self.a

    @property
    def sweep(self) -> "SweepOperator":
        """The inverse factor ``A⁻¹`` (alias for ``self.b``)."""
        return self.b

    def apply(
        self, rhs: "FullField", *,
        initial_guess: "FullField | None" = None,
    ) -> "TimedFullField":
        r"""``(P ∘ A⁻¹)·rhs`` via the substrate moment emit (fused).

        Routes to the forward operator's uniform private solve surface with
        ``moment_frame`` set — the ONE moment-emitting entry (the former
        ``solve_moments`` cross-reach retired here).  Principled-equivalent
        (≤ scale-relative ``4·N·eps``) to the inherited deforested body,
        which remains the verification oracle.

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
        return f"WindowedSweep({self.a!r} @ {self.b!r})"
