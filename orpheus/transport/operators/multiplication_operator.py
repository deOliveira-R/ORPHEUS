r"""The multiplication operator :math:`M[f]` — a coefficient field promoted.

This is the **promotion** of the grand report §5.5–5.7: a
:class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`
is the *symbol* of a zeroth-order operator, and
:class:`MultiplicationOperator` is that field *promoted* to a
:class:`~orpheus.numerics.operator.LinearOperator`.

The multiplier algebra (the §5.7 embedding)
============================================

For a (essentially bounded) coefficient field :math:`f \in L^\infty`,
the **multiplication operator** :math:`M[f]` acts on a state by
pointwise multiplication:

.. math::

    (M[f]\,\psi)(\vec r, \hat\Omega, g) \;=\;
        f(\vec r, g)\,\psi(\vec r, \hat\Omega, g).

The map :math:`M : L^\infty \to B(L^2)`, :math:`f \mapsto M[f]`, is a
**faithful unital ``*``-homomorphism** onto the diagonal subalgebra of
the bounded operators — the multiplier-algebra embedding. It satisfies
the law-suite that makes "a coefficient IS a (diagonal) operator"
literally true:

* **unital** — :math:`M[1] = I` (the constant-one field promotes to
  the identity);
* **zero** — :math:`M[0] = 0` (the zero coefficient promotes to the
  :class:`~orpheus.numerics.operator.ZeroOperator`, codomain-aware: a
  *zero source*, never a zero flux);
* **linear** — :math:`M[a f + b g] = a\,M[f] + b\,M[g]` (the promotion
  is a linear map on the coefficient vector space);
* **multiplicative (homomorphism)** — :math:`M[f]\,M[g] = M[f\,g]`
  (composition of two multipliers is the multiplier of the product;
  the algebra of coefficients embeds as the algebra of diagonal
  operators);
* **``*`` (self-adjoint for real f)** — :math:`M[f]^{*} = M[\bar f] =
  M[f]` for a real-valued coefficient, so ``M.H == M``;
* **spectrum** — :math:`\mathrm{spec}(M[f]) = \mathrm{ess\,range}(f)`,
  so :math:`M[f]` is invertible **iff** :math:`f` is
  bounded away from zero (:math:`\min|f| > 0`); its inverse is
  :math:`M[1/f]`.

Collision as a named instance
=============================

The §5.7 named instance is the collision multiplier
:math:`C = M[\sigma_t]` — a plain :class:`MultiplicationOperator` carrying
the total cross-section field as its coefficient. (#261 retired the former
``CollisionOperator`` thin subclass: it added nothing the base lacked once
the base gained the optional :attr:`space` for the W-D composition guard;
the ``L + C → StreamingCollisionOperator`` sweep dispatch lives on the SN-specific
:class:`~orpheus.sn.operators.streaming.StreamingOperator`, keyed on this base type —
a transport multiplier cannot dispatch back onto an ``sn`` operator.)
The collision rate :math:`\sigma_t\,\psi` turns a flux into a *source*
(a collision-rate density), so :meth:`apply` emits an
:class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`
codomain (units gain
:data:`~orpheus.numerics.units.ANGULAR_FLUX_UNITS` ×
:data:`~orpheus.numerics.units.CROSS_SECTION_UNITS` =
:data:`~orpheus.numerics.units.ANGULAR_RATE_UNITS`), while :meth:`solve`
returns to a flux.

Delegation to the numerics broadcast engine
============================================

The raw pointwise multiply is delegated to the N-D
:class:`~orpheus.numerics.operator.DiagonalOperator` broadcast engine
(#257 S3a): for a per-cell-per-group coefficient on a per-ordinate
carrier ``(N, ng, *spatial)``, the engine
``DiagonalOperator(f.values, broadcast_axes=(0,))`` broadcasts the
coefficient over the leading ordinate axis, so
``engine.apply(psi.interior.values) == f.values[None] * psi.interior.values``.
This operator owns only the *typed codomain* (the field wrapping); the
arithmetic lives once, in the engine (``coding-elegance`` Pattern 2).

The ends select the body (CS4c step 5)
======================================

:class:`MultiplicationOperator` stores ONLY the coefficient field and its
two ends. Which body ``apply`` / ``solve`` run is decided ONCE, at
construction, from the DOMAIN: a composite ``FullFieldSpace`` selects the
lifted body (:func:`~orpheus.transport.operators.lift.lift_bulk_action` —
the multiply on the bulk block, the zero source/sink — or, for ``solve``,
the zero flux — on the trace; the output leaf is the operand's role
partner, so one body serves the angular and the scalar families without
a carrier parse), a plain bulk space selects the bare-array body (the
cross-model escape hatch: the homogeneous / diffusion / depletion outer
loops feed bare arrays of the bound shape). Every off-binding carrier is
a typed refusal naming the operator
(:func:`~orpheus.transport.operators.lift.admit_composite` /
:func:`~orpheus.transport.operators.lift.admit_array`); the
``singledispatchmethod`` table and its per-family ``isinstance`` arms
retired with the selection. Both bodies run the ONE engine: a scalar
operand whose rank equals the coefficient's rides the engine's degenerate
one-ordinate lift (bit-identical to the direct multiply), an angular
operand the engine's leading-axis broadcast.

References
----------

* Grand Report v3 §5.5–5.7 — the field hierarchy and the operator
  promotion ``C = M[σ_t]``.
* ``.claude/agent-memory/cross-domain-attacker/coefficient_field_promotion_frames.md``
  — Frame 1 (the multiplier-algebra embedding).
* ``.claude/plans/issue_257_coefficient_field_promotion.md`` — S3.
* Reed, M. & Simon, B. (1980). *Methods of Modern Mathematical
  Physics I: Functional Analysis*. §VII (multiplication operators,
  spec = ess-range).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Callable, overload

import numpy as np

from orpheus.transport.operators.bound_operator import BoundOperator
from orpheus.transport.operators.lift import (
    admit_array,
    admit_composite,
    embed_bulk_assembly,
    lift_bulk_action,
)
from orpheus.numerics.operator import (
    BlockRole,
    DiagonalOperator,
    InverseOperator,
    LinearOperator,
    NotInvertible,
)
from orpheus.numerics.spaces.full_field_space import FullFieldSpace
from orpheus.transport.fields._bases import BulkField, FieldRole
from orpheus.transport.full_field import FullField

if TYPE_CHECKING:
    from orpheus.numerics.assembled_operator import SparseAssembledOperator
    from orpheus.numerics.field import Field
    from orpheus.numerics.space import FunctionSpace
    from orpheus.transport.fields.cross_section_field import CrossSectionField


__all__ = ["MultiplicationOperator"]

#: The engine verb a body runs — ``engine.apply`` or ``engine.solve``.
_EngineVerb = Callable[[np.ndarray], np.ndarray]


@dataclass(eq=False)
class MultiplicationOperator(BoundOperator):
    r"""The promotion :math:`M[f]` of a coefficient field to a diagonal operator.

    Stores ONLY the coefficient and its ends; the raw pointwise multiply
    is delegated to the N-D
    :class:`~orpheus.numerics.operator.DiagonalOperator` broadcast engine,
    and the body (composite vs bare array) is selected from the domain
    at construction (module docstring).

    See the module docstring for the multiplier-algebra law-suite
    (:math:`M[1]=I`, :math:`M[0]=0`, linearity, homomorphism,
    self-adjointness, :math:`\mathrm{spec} = \mathrm{ess\,range}`).

    Parameters
    ----------
    coefficient : CrossSectionField
        The coefficient field :math:`f` (units ``1/cm`` for a cross
        section). Its ``.values`` is the ``(ng, *spatial)`` per-cell
        per-group array the broadcast engine consumes.

    Structural axes
    ---------------

    Always adjointable (self-adjoint); invertible
    **iff** :math:`\min|f| > 0` (the spectrum law: invertible iff the
    coefficient is bounded away from zero). This is the honest
    value-dependent gate (``coding-elegance`` Pattern 4); the legacy
    collision operator solved unconditionally and produced silent IEEE
    NaN on a zero entry — :class:`MultiplicationOperator` refuses
    eagerly (:class:`~orpheus.numerics.operator.NotInvertible` from
    ``inverse()``/``solve``) instead, so a downstream composer never
    hits a broken inverse at call time. Always assemblable: the bulk
    diagonal on the plain ends, embedded in the composite flat layout on
    composite ends.
    """

    coefficient: "CrossSectionField"

    # The two spaces are the INHERITED kw-only mandatory ends
    # (:class:`~orpheus.transport.operators.bound_operator.BoundOperator`,
    # CS4c step 2): the multiplier is space-endomorphic — flux block →
    # source block, same shape — so every binding passes the SAME space
    # to both, and the tier-2 :meth:`from_mesh` sugar spells exactly
    # that. With the ends mandatory, the OperatorSum/OperatorProduct
    # composition guard validates every build, and ``.H``'s Riesz legs
    # can never see a ``None``.

    #: The N-D broadcast engine (#257 S3a) the multiply delegates to,
    #: built ONCE in :meth:`__post_init__` over the immutable
    #: ``coefficient.values`` (``coding-elegance`` Pattern 2: the
    #: arithmetic AND the spectrum gate live in ONE engine object that
    #: both the capability freeze and the action read; the coefficient is
    #: immutable, so the stored engine cannot go stale).
    engine: "DiagonalOperator" = field(init=False, repr=False)

    # Multiplication by a per-cell per-group coefficient is a BULK
    # operator — diagonal in (cell, group, ordinate), no boundary action
    # (A_bb only; Wave O / Issue #208). Class-level constant (unannotated
    # so the dataclass does not treat it as a field).
    block_role = BlockRole.BULK

    def __post_init__(self) -> None:
        # Build the broadcast engine ONCE over the coefficient. The N-D
        # engine (DiagonalOperator(f.values, broadcast_axes=(0,))) broadcasts
        # the (ng, *spatial) coefficient over the leading ordinate axis of a
        # (N, ng, *spatial) carrier, so engine.apply(x) == f.values[None]*x.
        # It already encodes "invertible iff every entry != 0", so inheriting
        # its capability set single-sources the spectrum law
        # spec(M[f]) = ess-range(f) (Pattern 4 — the transport operator and
        # the numerics engine agree by construction, not by a copied
        # predicate). Stored, not rebuilt per call: a fresh engine would
        # re-run the O(size) `np.all(coeff != 0)` spectrum scan on every
        # matvec and discard a capability already frozen here.
        self.engine = DiagonalOperator(
            self.coefficient.values, broadcast_axes=(0,),
        )
        # CS4a K2: the binding refuses a space whose EnergyAxis states a
        # different group count than the coefficient (guard reach + the
        # deliberate axes-less inertness: _energy_conformity docstring).
        # ``values.shape[0]`` IS ng by CrossSectionField's own declared
        # ``(ng, *spatial)`` layout (#196) — deliberately NOT
        # ``coefficient.ng``, which is a MESH read-through (fields/_bases)
        # that CS4b retires: the guard compares the space against the
        # DATA, never against a second metadata source (CS4a-R EE-3).
        self._assert_energy_extent_both_ends(
            self.coefficient.values.shape[0],
            operator="MultiplicationOperator",
        )
        # The SELECTION (CS4c step 5): the domain names the body. This is
        # a parse of the SPACE at construction, not of a carrier per call.
        domain = self.domain
        bulk_shape = (
            tuple(domain.interior_space.shape)
            if isinstance(domain, FullFieldSpace) and domain.interior_space is not None
            else tuple(domain.shape)
        )
        self._composite = isinstance(domain, FullFieldSpace)
        # A bulk whose rank equals the coefficient's has NO ordinate axis
        # (the scalar families): the ONE engine (and its frozen spectrum
        # gate) serves it through a degenerate one-ordinate lift, dropped
        # after — bit-identical to the direct ``coefficient.values * x``.
        # An angular bulk rides the engine's leading-axis broadcast.
        self._degenerate_lift = len(bulk_shape) == self.coefficient.values.ndim

    # ── the ONE multiply, on either bulk rank ───────────────────────

    def _run(self, verb: _EngineVerb, values: np.ndarray) -> np.ndarray:
        if self._degenerate_lift:
            return verb(values[None])[0]
        return verb(values)

    def _lift(
        self, x: object, verb: _EngineVerb, *, role: FieldRole,
    ) -> FullField:
        r"""The composite body: the engine on the bulk block, the zero
        field of ``role`` on the trace, the output leaf the operand's
        role partner (one body for the angular and the scalar families)."""
        psi = admit_composite(self, x, end="domain")

        def act(bulk: BulkField) -> "Field":
            return bulk.into_role(role, self._run(verb, bulk.values))

        return lift_bulk_action(psi, act, trace_role=role)

    @property
    def is_invertible(self) -> bool:
        # M[f] is invertible iff every coefficient entry is nonzero (the
        # spectrum law min|f| > 0). DELEGATE to the broadcast engine (a
        # DiagonalOperator) — the same single-source that already carries the
        # capability set, so the predicate cannot drift from the spectrum.
        return self.engine.is_invertible

    def inverse(self) -> "InverseOperator":
        r"""Return :math:`M[f]^{-1}` as an :class:`InverseOperator` over this leaf.

        Delegation, NOT a reciprocal-field twin (#226 taxonomy step 1): the
        returned object's ``apply`` IS :meth:`solve` — the typed division
        :math:`q \mapsto q/f` with the flux-role codomain — bit-identical to
        the gated call. Materializing :math:`M[1/f]` instead would mint a
        units-dishonest "reciprocal cross-section" coefficient
        (:math:`1/\Sigma` is a mean free path, a different named quantity)
        and flip the flux/source carrier roles this class hard-binds; the
        division realization carries the inverse semantics without either.
        """
        if not self.is_invertible:
            raise NotInvertible(
                "MultiplicationOperator.inverse requires min|f| > 0 (the "
                "spectrum law); this coefficient has a zero entry."
            )
        return InverseOperator(self)

    @property
    def is_adjointable(self) -> bool:
        # Diagonal multiplication is its own structural transpose; delegate to
        # the engine (always adjointable).
        return self.engine.is_adjointable

    @property
    def is_metric_free_adjoint(self) -> bool:
        r"""``True`` — M[f] is a real multiplier, so its Hilbert adjoint is
        metric-free (the module docstring's self-adjointness law: pointwise
        multiplication commutes with every diagonal metric). A BOUND class
        can carry the metric-free predicate: bound-ness says where it acts,
        metric-freeness says its adjoint needs no sandwich. (Until CS4c
        step 2 this also exempted the then-legal SPACE-LESS ``M`` from the
        S4-amendment's unbound-``.H`` refusal; with the ends mandatory
        that exemption arm is unreachable from this class — the Riesz
        legs always execute, and commute.)"""
        return True

    # ── The assembly mode (stencil-assembly 2b) ────────────────────────

    @property
    def is_assemblable(self) -> bool:
        r"""``True`` — the bulk diagonal :math:`\mathrm{diag}(f)` is
        emittable on either binding: on the plain bulk space directly, on
        a composite end embedded in the ``[bulk C-ravel | trace]`` flat
        layout (zero trace rows — a multiplier has no face action)."""
        return True

    def _bulk_assembly(self, bulk_space: "FunctionSpace") -> "SparseAssembledOperator":
        from scipy import sparse

        from orpheus.numerics.assembled_operator import SparseAssembledOperator

        bulk_shape = tuple(bulk_space.shape)
        n_bulk = int(np.prod(bulk_shape))
        # One-source discipline: the diagonal IS ``self.coefficient.values``
        # broadcast over the bulk shape — the SAME array (and the same
        # prepend-broadcast semantics) the engine multiplies at apply
        # time, with NO arithmetic performed on it, so ``assembled @ x``
        # reproduces ``apply``'s per-entry multiply bit-for-bit.
        # Family-blind: the scalar ``(ng, nx)`` bulk broadcasts identically
        # and an angular ``(N, ng, …)`` bulk gains the leading ordinate
        # axis, exactly as the engine's ``broadcast_axes=(0,)``.
        diagonal = np.ascontiguousarray(
            np.broadcast_to(self.coefficient.values, bulk_shape)
        ).ravel()
        nonzero = diagonal != 0.0
        indices = np.arange(n_bulk)[nonzero]
        matrix = sparse.coo_array(
            (diagonal[nonzero], (indices, indices)), shape=(n_bulk, n_bulk),
        )
        return SparseAssembledOperator(matrix, domain=bulk_space, codomain=bulk_space)

    def assemble(self) -> "SparseAssembledOperator":
        r"""Emit the bulk diagonal :math:`\mathrm{diag}(f)` — on the plain
        ends directly, on composite ends in the composite flat layout
        (:func:`~orpheus.transport.operators.lift.embed_bulk_assembly`)."""
        domain = self.domain
        if isinstance(domain, FullFieldSpace) and domain.interior_space is not None:
            return embed_bulk_assembly(
                self._bulk_assembly(domain.interior_space),
                domain=domain, codomain=self.codomain,
            )
        return self._bulk_assembly(domain)

    @classmethod
    def from_mesh(
        cls,
        sigma: "np.ndarray | CrossSectionField",
        mesh,
        *,
        space: "FunctionSpace | None" = None,
    ) -> "MultiplicationOperator":
        r"""Construct from a coefficient that may be a bare ndarray.

        ``sigma`` is either a bare ``(ng, *spatial)`` :class:`numpy.ndarray`
        — wrapped into a
        :class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`
        on ``mesh.bulk_space``, exactly as the production
        ``mat_xs.total_cross_section_field`` accessor builds it — or an
        already-typed :class:`CrossSectionField` (passed straight through).
        ``space`` is the optional composite
        :class:`~orpheus.numerics.space.FunctionSpace` for the composition
        guard (e.g. ``sn_mesh.full_field_space``).

        #261: folded up from the retired ``CollisionOperator(sn_mesh, sigma)``
        constructor — the legacy / test-caller convenience that accepts a raw
        array. Production passes a :class:`CrossSectionField` directly to the
        dataclass; this classmethod is the bare-array entry point.
        """
        from orpheus.transport.fields.cross_section_field import CrossSectionField

        coefficient = (
            sigma
            if isinstance(sigma, CrossSectionField)
            else CrossSectionField(values=np.asarray(sigma), space=mesh.bulk_space)
        )
        # The space defaults to the mesh's own, most-structured-first (CS1):
        # a method mesh's composite ``full_field_space`` wins (SNMesh /
        # DiffusionMesh short-circuit here, untouched), else the carrier's
        # axis-built ``bulk_space`` (every MaterialMesh carries one). Since
        # CS4a K2 this bulk_space arm has NO production caller — the
        # homogeneous solver mints its own Energy ⊗ point space via
        # ``_pose_space`` and constructs C directly, so the resolution
        # chain below serves bare-array/test callers only (CS4a-R CEN-2).
        # Pass ``space=`` to override. The first arm
        # makes ``from_mesh(σ, sn_mesh)`` a faithful drop-in for the retired
        # ``CollisionOperator(sn_mesh, σ)``, which reached the same space via
        # ``sn_mesh.full_field_space`` (#261).
        if space is None:
            space = getattr(mesh, "full_field_space", None)
        if space is None:
            space = getattr(mesh, "bulk_space", None)
        if space is None:
            raise TypeError(
                "MultiplicationOperator.from_mesh: the mesh carries neither "
                "full_field_space nor bulk_space and no space= was passed — "
                "a binding needs its ends (CS4c: the ctor's domain/codomain "
                "are mandatory; this classmethod is the endomorphism sugar "
                "that supplies both from ONE space)."
            )
        return cls(coefficient=coefficient, domain=space, codomain=space)

    # ── The §5.7 promotion: f ↦ M[f], through the body the ends select ──

    @overload
    def apply(self, psi: FullField, /) -> FullField: ...
    @overload
    def apply(self, psi: np.ndarray, /) -> np.ndarray: ...
    def apply(self, psi: "FullField | np.ndarray", /) -> "FullField | np.ndarray":
        r"""Forward action :math:`M[f]\,\psi = f \cdot \psi`.

        On a composite binding: the pointwise multiply on the bulk block
        (per cell, per group, broadcast over any leading ordinate axis)
        repackaged as a SOURCE composite — multiplying a flux by a cross
        section yields a collision-rate density — with the implicit-zero
        source/sink of the operand's family on the trace (a multiplier
        has no face-trace action; the cell-balance :math:`f\,\psi` term
        is a CELL quantity). On a plain binding: the bare
        ``(ng, *spatial)`` multiply, mesh-free (the cross-model escape
        hatch, #276).
        """
        if self._composite:
            return self._lift(psi, self.engine.apply, role=FieldRole.SOURCE_SINK)
        return self._run(self.engine.apply, admit_array(self, psi, end="domain"))

    @overload
    def solve(self, q: FullField, /) -> FullField: ...
    @overload
    def solve(self, q: np.ndarray, /) -> np.ndarray: ...
    def solve(self, q: "FullField | np.ndarray", /) -> "FullField | np.ndarray":
        r"""Inverse action :math:`M[f]^{-1}\,q = q / f = M[1/f]\,q`.

        Requires invertibility (the spectrum law: :math:`\min|f| > 0`);
        the engine raises
        :class:`~orpheus.numerics.operator.NotInvertible` otherwise.
        The codomain returns to a flux (the inverse of "flux → source" is
        "source → flux") in the operand's family — the same body as
        :meth:`apply` with the FLUX role on the output leaf and the trace.
        """
        if self._composite:
            return self._lift(q, self.engine.solve, role=FieldRole.FLUX)
        return self._run(self.engine.solve, admit_array(self, q, end="codomain"))

    @overload
    def apply_transpose(self, psi: FullField, /) -> FullField: ...
    @overload
    def apply_transpose(self, psi: np.ndarray, /) -> np.ndarray: ...
    def apply_transpose(self, psi: "FullField | np.ndarray", /) -> "FullField | np.ndarray":
        r"""Adjoint action :math:`M[f]^{*}\,\psi = M[\bar f]\,\psi = M[f]\,\psi`.

        Equal to :meth:`apply` — a real-valued multiplier is self-adjoint
        (:math:`M[f]^{*} = M[\bar f] = M[f]`), so ``M.H == M`` (the
        metric-blind Euclidean transpose).
        """
        return self.apply(psi)
