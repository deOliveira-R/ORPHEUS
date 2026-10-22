r"""System B as an independent composite: ``Composite[interior ⊕ boundary]`` of ψ½.

The coupled-block campaign's Phase B poses the curvilinear ψ½ ray as **System B** —
its own ``interior ⊕ boundary`` composite, exactly parallel to System A's
:class:`~orpheus.transport.full_field.FullField`
(``Composite[AngularFlux, AngularBoundaryFlux]``). Here that composite is
realized:

.. code-block:: python

    RadialCharacteristicField
        = Composite[RadialCharacteristicInteriorField,    # the marched cells (A_BB)
                    RadialCharacteristicBoundaryField]     # the r = R corner (B_b)

It is a **trivial** subclass — no hook overrides — because System B adds no third
block (unlike the historical 3-block ``FullField``): the generic
:class:`~orpheus.transport.full_field.Composite` base (its slot bounds relaxed to
``Field`` in Phase B) holds the entire 2-block vector-space algebra, and System B
inherits it whole. The subclass adds only the concrete-locus ``__post_init__``
guard and the role-keyed birth factories (:meth:`flux_zeros` /
:meth:`source_zeros` — both presence-gated at the member-space parse, CS4b
S5). The per-ordinate source fold onto a q½ composite is SN's algorithm (the
level structure of the SN quadrature and Hébert's Eq. 3.432 fold), so it lives
in the SN package:
:func:`orpheus.sn.operators.radial_characteristic.radial_characteristic_source_from_angular`
(layer-DAG P3).

**Role-erased slots (B.2b DP2, the FullField precedent).** The static parameters
bind the locus FIELD BASES, not the flux leaves — exactly as ``FullField`` binds
``BulkField``/``BoundaryField`` — so ONE composite class carries the flux state
(the iterate, whose same-class differences are ordinary signed ψ½ composites —
campaign 1 CS3) and the source emission (an operator ``.apply`` output:
``A_BA`` / ``B_b``). Role identity lives
on the MEMBERS (the Field class-identity gate rejects cross-role sums); a
consumer that needs a specific role parses it off the member (the #289-F2
discipline).

Since Phase C (4e) this composite is the NATIVE ψ½ representation end to end:
the fused ``(L+C)`` walk marches its members directly, so the historical
unified ``RadialCharacteristicField`` leaf family, its space, and the
``from_unified``/``to_unified`` bridge are RETIRED (the demotion licence —
``to_unified ∘ from_unified == id`` — was discharged by the walk-slot rewrite
landing bit-identically on the frozen baselines).

References
----------

* ``.claude/plans/archive/coupled_block_operator_campaign.md`` — Phase B (pose System B
  as an independent composite); B.1a (base relaxation), B.1b (split spaces),
  B.1c (split leaves), B.1d (this composite), Phase C 4e (the walk-slot
  rewrite + unified retirement).
* ``coding-elegance`` Pattern 2 (the 2-block algebra lives ONCE on ``Composite``;
  System B is a no-op extension) + Pattern 4 (the leaf-type guard makes a wrong
  pairing unconstructable).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from orpheus.transport.fields._bases import (
    RadialCharacteristicBoundaryField,
    RadialCharacteristicInteriorField,
)
from orpheus.transport.fields.radial_characteristic_boundary_flux import (
    RadialCharacteristicBoundaryFlux,
)
from orpheus.transport.fields.radial_characteristic_interior_flux import (
    RadialCharacteristicInteriorFlux,
)
from orpheus.transport.full_field import Composite
from orpheus.transport.source_sinks.radial_characteristic_boundary_source_sink import (
    RadialCharacteristicBoundarySourceSink,
)
from orpheus.transport.source_sinks.radial_characteristic_interior_source_sink import (
    RadialCharacteristicInteriorSourceSink,
)

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace
    from orpheus.numerics.spaces.full_field_space import FullFieldSpace

__all__ = ["RadialCharacteristicField"]

@dataclass(frozen=True, kw_only=True)
class RadialCharacteristicField(
    Composite[RadialCharacteristicInteriorField, RadialCharacteristicBoundaryField],
):
    r"""System B: the ψ½ ray as an ``interior ⊕ boundary`` composite.

    ``interior`` is a
    :class:`~orpheus.transport.fields._bases.RadialCharacteristicInteriorField`
    (the marched cells — flux state or source emission; a difference of two
    flux composites is the same flux class carrying signed values, campaign 1
    CS3); ``boundary`` a
    :class:`~orpheus.transport.fields._bases.RadialCharacteristicBoundaryField`
    (the r = R corner). Role-erased slots (the FullField precedent — see the
    module docstring): role identity lives on the members. The 2-block algebra
    (``±``, scalar ``·``, ``to_flat`` / ``from_flat``, ``copy`` — plain V
    arithmetic per block since campaign 1 CS3) is inherited whole from
    :class:`~orpheus.transport.full_field.Composite`.
    """

    def __post_init__(self) -> None:
        # System B narrows the generic slots to the ψ½ split loci (the concrete
        # guard belongs with the concrete specialization, as FullField guards
        # BulkField / BoundaryField). Guard the FIELD BASE, not a role leaf, so
        # both flux and source (operator-emission) composites are admitted —
        # role-erased slots, B.2b DP2.
        if not isinstance(self.interior, RadialCharacteristicInteriorField):
            raise TypeError(
                f"{type(self).__name__}: interior must be a "
                f"RadialCharacteristicInteriorField (its flux / source "
                f"leaf); got {type(self.interior).__name__}"
            )
        if not isinstance(self.boundary, RadialCharacteristicBoundaryField):
            raise TypeError(
                f"{type(self).__name__}: boundary must be a "
                f"RadialCharacteristicBoundaryField (its flux / source "
                f"leaf); got {type(self.boundary).__name__}"
            )
        super().__post_init__()  # base: Field + mesh-identity.

    # ── Construction ─────────────────────────────────────────────────

    @classmethod
    def flux_zeros(
        cls, space: "FullFieldSpace | None",
    ) -> "RadialCharacteristicField":
        r"""A zero ψ½ FLUX composite on System B's member space (presence-gated).

        Space-keyed since CS4b S5 (the mesh-keyed ``from_mesh`` retired with
        the sugar tier): the caller passes the carrier's cached
        :attr:`~orpheus.sn.problem.SNProblem.radial_characteristic_field_space`,
        whose blocks ARE the split ψ½ spaces (``is``-shared with
        ``radial_characteristic_interior_space`` /
        ``radial_characteristic_boundary_space``). ``None`` — what a
        NON-carrying mesh's cached property returns (R12a) — is REFUSED with
        the presence diagnosis, so System B stays unconstructable exactly
        where it does not exist (presence = block existence), now diagnosed
        at the composite seam instead of inside a leaf factory.
        """
        cells, corner = cls._require_member_space(
            space, context=f"{cls.__name__}.flux_zeros",
        )
        return cls(
            interior=RadialCharacteristicInteriorFlux.zeros(cells),
            boundary=RadialCharacteristicBoundaryFlux.zeros(corner),
        )

    @classmethod
    def _require_member_space(
        cls, space: "FullFieldSpace | None", *, context: str,
    ) -> "tuple[FunctionSpace, FunctionSpace]":
        r"""Parse the member-space argument — the shared presence gate.

        The ONE spelling of the R12a absence diagnosis for the composite
        allocators (:meth:`flux_zeros` / :meth:`source_zeros` — Pattern 2).
        Returns the narrowed ``(cells, corner)`` block pair so the leaf
        allocations type-check without re-narrowing.
        """
        if space is None:
            raise ValueError(
                f"{context}: System B is absent on this carrier — "
                f"radial_characteristic_field_space is None (R12a: no "
                f"μ-level consumes independent starting-direction state; "
                f"Cartesian and the production cylinder rules land here). "
                f"System B's blocks are spelled absent, never zero-DOF "
                f"fields."
            )
        return space.interior_space, space.trace_space

    @classmethod
    def require_member(
        cls, x: object, *, space: "FullFieldSpace", context: str,
    ) -> "RadialCharacteristicField":
        r"""Parse ``x`` as a System-B member carrier on ``space`` — the shared
        block-boundary guard.

        Every System-B block boundary (``A_BB`` / ``A_AB`` / ``B_b``) receives
        a :class:`RadialCharacteristicField` and must refuse (i) a foreign
        carrier type and (ii) block spaces disagreeing with the operator's
        ray member spaces in CONTENT (the space-content invariant, per block
        since CS4b S3 — the member's legs, the operator's coefficients, and
        the radial widths must not desync). One parse, three consumers (coding-elegance
        Pattern 2; parse-don't-validate at the boundary). ``space`` is System
        B's member composite (interior ⊕ boundary corner — the operator's
        bound space; un-weld arc O-1: the SPACE is the contract, not the
        mesh that minted it). ``context`` names the refusing surface
        (``"RadialCharacteristicOperator.apply"``) so the error reads at the
        call site.
        """
        if not isinstance(x, cls):
            raise TypeError(
                f"{context}: expected a RadialCharacteristicField "
                f"(System B's member carrier); got {type(x).__name__}."
            )
        if (
            x.interior.space != space.interior_space
            or x.boundary.space != space.trace_space
        ):
            raise ValueError(
                f"{context}: the input field's block spaces must agree with "
                f"the operator's ray member spaces in content (space-content "
                f"invariant, per block — the composite's own ``==`` is "
                f"name+shape and cannot see blocks); got interior "
                f"{x.interior.space!r} vs {space.interior_space!r}."
            )
        return x

    # ── The source-role birth factories ───────────────────────────────

    @classmethod
    def source_zeros(
        cls, space: "FullFieldSpace | None",
    ) -> "RadialCharacteristicField":
        r"""A zero q½ SOURCE composite on System B's member space (presence-gated).

        The source-role sibling of :meth:`flux_zeros` — the buffer the joint
        matvec/transpose surfaces fill with System B's emitted rows (an
        operator ``.apply`` output is source-role, never flux-role). Same
        space-keyed contract and R12a ``None`` refusal (CS4b S5; the
        mesh-keyed ``source_zeros_on`` retired with the sugar tier).
        """
        cells, corner = cls._require_member_space(
            space, context=f"{cls.__name__}.source_zeros",
        )
        return cls(
            interior=RadialCharacteristicInteriorSourceSink.zeros(cells),
            boundary=RadialCharacteristicBoundarySourceSink.zeros(corner),
        )
