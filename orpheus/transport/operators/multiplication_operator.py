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

Mesh-free (the carrier carries the mesh)
========================================

:class:`MultiplicationOperator` stores ONLY the coefficient field; the
output spaces are read off the operand's blocks at apply time,
faithful to the structure the legacy collision operator used. The
coefficient field is itself mesh-bound (it was built ``from_mesh``), so
the operator's domain is implicit in the field it carries.

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
from functools import singledispatchmethod
from typing import TYPE_CHECKING, Any, overload

import numpy as np

from orpheus.transport.operators._energy_conformity import (
    assert_energy_extent_conforms,
)
from orpheus.numerics.operator import (
    BlockRole,
    DiagonalOperator,
    InverseOperator,
    LinearOperator,
    NotInvertible,
)
# Runtime import for ``singledispatchmethod.register`` (mirrors fission.py):
# ``FullField`` is a leaf in the SN dependency graph (it imports no operators),
# so this module-level import is cycle-free.
from orpheus.transport.full_field import FullField

if TYPE_CHECKING:
    from orpheus.numerics.assembled_operator import SparseAssembledOperator
    from orpheus.numerics.space import FunctionSpace
    from orpheus.transport.fields.cross_section_field import CrossSectionField


__all__ = ["MultiplicationOperator"]


@dataclass(eq=False)
class MultiplicationOperator(LinearOperator["FullField"]):
    r"""The promotion :math:`M[f]` of a coefficient field to a diagonal operator.

    Stores ONLY the coefficient; the mesh is read off the carrier at
    apply time. The raw pointwise multiply is delegated to the N-D
    :class:`~orpheus.numerics.operator.DiagonalOperator` broadcast engine
    on ``psi.interior.values``; this operator owns the typed codomain.

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
    hits a broken inverse at call time.
    """

    coefficient: "CrossSectionField"

    #: The composite :class:`~orpheus.numerics.space.FunctionSpace` this
    #: multiplier acts on (optional). When supplied, :attr:`domain` /
    #: :attr:`codomain` report it, so the operator joins the
    #: :class:`~orpheus.numerics.operator.OperatorSum` / ``OperatorProduct``
    #: composition guard (W-D) — e.g. the SN ``full_field_space`` makes the
    #: ``L + C`` build VALIDATE on every within-group solve. When ``None``
    #: (the default) the operator is space-anonymous: the guard skips it and
    #: the carrier still supplies the mesh at apply time (the mesh-free
    #: contract). #261: folded up from the retired ``CollisionOperator``,
    #: which reached this space through ``sn_mesh.full_field_space``.
    space: "FunctionSpace | None" = field(default=None)

    #: The N-D broadcast engine (#257 S3a) the multiply delegates to,
    #: built ONCE in :meth:`__post_init__` over the immutable
    #: ``coefficient.values`` (``coding-elegance`` Pattern 2: the
    #: arithmetic AND the spectrum gate live in ONE engine object that
    #: both the capability freeze and the action read; the coefficient is
    #: immutable, so the stored engine cannot go stale).
    engine: "DiagonalOperator" = field(init=False, repr=False)

    #: Inherited from the engine's spectrum gate in :meth:`__post_init__`.

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
        assert_energy_extent_conforms(
            self.space, self.coefficient.values.shape[0],
            operator="MultiplicationOperator",
        )

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
        metric-freeness says its adjoint needs no sandwich — so the
        S4-amendment's unbound-``.H`` refusal correctly exempts a
        space-less ``M`` (its Euclidean ``.H`` IS the Hilbert adjoint)."""
        return True

    # ── The assembly mode (stencil-assembly 2b) ────────────────────────

    @property
    def is_assemblable(self) -> bool:
        r"""``True`` iff the composite flat layout is known — a block-bearing
        :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`
        was threaded at construction. A multiplier without one honestly
        refuses: a bare/space-less multiplier has no layout at all, and a
        plain bulk space (e.g. the carrier's axis-built ``bulk_space``,
        CS1) carries no bulk ⊕ trace composite flat layout — there is no
        global DOF numbering to emit into."""
        from orpheus.numerics.spaces.full_field_space import FullFieldSpace

        return (
            isinstance(self.space, FullFieldSpace)
            and self.space.interior_space is not None
        )

    def assemble(self) -> "SparseAssembledOperator":
        r"""Emit the bulk diagonal :math:`\mathrm{diag}(f)` in the composite
        flat layout (zero trace rows — a multiplier has no face action).

        One-source discipline: the diagonal IS
        ``self.coefficient.values`` broadcast over the bulk shape — the
        SAME array (and the same prepend-broadcast semantics) the
        engine multiplies at apply time, with NO arithmetic performed
        on it, so ``assembled @ x`` reproduces ``apply``'s per-entry
        multiply bit-for-bit. Family-blind: the scalar ``(ng, nx)``
        bulk broadcasts identically and an angular ``(N, ng, …)`` bulk
        gains the leading ordinate axis, exactly as the engine's
        ``broadcast_axes=(0,)``.
        """
        from orpheus.numerics.assembled_operator import SparseAssembledOperator
        from orpheus.numerics.operator import MissingAssembly
        from orpheus.numerics.spaces.full_field_space import FullFieldSpace
        from scipy import sparse

        space = self.space
        if not isinstance(space, FullFieldSpace) or space.interior_space is None:
            raise MissingAssembly(
                f"MultiplicationOperator.assemble requires a block-bearing "
                f"FullFieldSpace (the composite flat layout); this "
                f"multiplier's space is "
                f"{'None' if space is None else type(space).__name__}, "
                f"which carries no composite flat layout to emit into."
            )
        interior_shape = tuple(space.interior_space.shape)
        n_interior = int(np.prod(interior_shape))
        n_total = int(space.shape[0])
        diagonal = np.ascontiguousarray(
            np.broadcast_to(self.coefficient.values, interior_shape)
        ).ravel()
        nonzero = diagonal != 0.0
        indices = np.arange(n_interior)[nonzero]
        matrix = sparse.coo_array(
            (diagonal[nonzero], (indices, indices)),
            shape=(n_total, n_total),
        )
        return SparseAssembledOperator(matrix, domain=space, codomain=space)

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
        on ``mesh`` via the same ``from_mesh`` factory the production
        ``mat_xs.total_cross_section_field`` accessor uses — or an
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
        return cls(coefficient=coefficient, space=space)

    # ── Operator-algebra space metadata (W-D / #261) ─────────────────────
    @property
    def domain(self) -> "FunctionSpace | None":
        r"""The composite space the multiplier acts on, or ``None``.

        Returns the optional :attr:`space`: when set (e.g. the SN
        ``full_field_space`` threaded at construction), the ``L + C``
        :class:`~orpheus.numerics.operator.OperatorSum` composition guard
        VALIDATES the build (equal domains AND codomains) on every
        within-group solve, instead of silently skipping a ``None``-spaced
        operand; when ``None``, the operator stays space-anonymous (the
        mesh-free default — the carrier supplies the mesh at apply time).
        The multiplier is space-endomorphic — flux block → source block,
        same shape — so ``codomain == domain``.
        """
        return self.space

    @property
    def codomain(self) -> "FunctionSpace | None":
        # Endomorphic on the composite space (see :meth:`domain`).
        return self.space

    # ── The §5.7 promotion: f ↦ M[f], dispatched on the input carrier ────
    #
    # Mirrors :class:`~orpheus.transport.operators.fission.FissionOperator`:
    # ``_apply_impl`` is the runtime ``singledispatchmethod``; the public
    # ``apply`` name aliases it (below) with the per-carrier ``@overload``
    # typing surface. Two arms, ONE multiply (``self.coefficient.values``):
    #
    # * :class:`FullField` — the SN per-ordinate carrier (the ~206-caller
    #   path): the engine's leading-axis broadcast, repackaged as a source
    #   composite (flux → collision-rate source).
    # * bare :class:`numpy.ndarray` — the MESHLESS ``(ng, *spatial)`` scalar/
    #   group block (the cross-model escape hatch: homogeneous / diffusion /
    #   depletion outer loops feed bare arrays). The pointwise multiply needs
    #   no mesh, so ``M[σ]`` composes on a mesh-free ``MaterialMesh`` (#276).

    @singledispatchmethod
    def _apply_impl(self, psi) -> "Any":
        raise TypeError(
            f"MultiplicationOperator.apply: unsupported input type "
            f"{type(psi).__name__}; expected FullField or numpy.ndarray. "
            f"Dispatch table is registered via @singledispatchmethod."
        )

    @_apply_impl.register
    def _(self, psi: FullField) -> "FullField":
        r"""Forward action :math:`M[f]\,\psi = f \cdot \psi` on the composite.

        The pointwise multiply :math:`f\,\psi` is per-cell per-group.
        The codomain is a *source* — multiplying a flux by a cross
        section yields a collision-rate density — and the boundary is
        the implicit-zero source/sink of the input's family (a
        multiplier has no face-trace action — the cell-balance
        :math:`f\,\psi` term is a CELL quantity). Two family arms
        (parsed loudly at the seam, #289 discipline):

        * **angular** bulk (SN, the ~206-caller path) — the engine's
          leading-ordinate-axis broadcast; out =
          :class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`
          ⊕ zero
          :class:`~orpheus.transport.source_sinks.angular_boundary_source_sink.AngularBoundarySourceSink`.
        * **scalar** bulk (diffusion / CP, #290 P4) — the direct
          per-cell-per-group multiply (no ordinate axis; the engine's
          broadcast axis is degenerate over the leading ``ng`` axis, so
          the coefficient is applied directly); out =
          :class:`~orpheus.transport.source_sinks.scalar_source_sink.ScalarSourceSink`
          ⊕ zero
          :class:`~orpheus.transport.source_sinks.scalar_boundary_source_sink.ScalarBoundarySourceSink`.
        """
        from orpheus.transport.fields._bases import AngularField, ScalarField
        from orpheus.transport.source_sinks import (
            AngularSourceSink,
            AngularBoundarySourceSink,
            ScalarSourceSink,
            ScalarBoundarySourceSink,
        )

        # Parse the family loudly at the seam (#289 discipline). The mesh
        # is read off the PARSED bulk, whose family declaration carries
        # the mesh type the widened composite surfaces erase (#290 P2).
        bulk = psi.interior
        if isinstance(bulk, AngularField):
            # CS4b S4 — the space route: every output block rides its
            # OPERAND block's space (role transition = new class, same
            # space; the zero trace block rides the input's trace space).
            out_bulk = self.engine.apply(bulk.values)
            return FullField(
                interior=AngularSourceSink(values=out_bulk, space=bulk.space),
                boundary=AngularBoundarySourceSink.zeros(psi.boundary.space),
            )
        if isinstance(bulk, ScalarField):
            # Scalar arm (#290 P4): the (ng, *spatial) bulk has NO
            # ordinate axis — lift it onto a degenerate 1-ordinate
            # leading axis so the ONE broadcast engine (and its frozen
            # spectrum gate) serves both families, then drop the axis.
            # Bit-identical to the direct ``coefficient.values * values``.
            out_bulk = self.engine.apply(bulk.values[None])[0]
            return FullField(
                interior=ScalarSourceSink(values=out_bulk, space=bulk.space),
                boundary=ScalarBoundarySourceSink.zeros(psi.boundary.space),
            )
        raise TypeError(
            f"MultiplicationOperator composite apply: angular- or "
            f"scalar-family bulk required; got {type(bulk).__name__}."
        )

    @_apply_impl.register
    def _(self, phi: np.ndarray) -> np.ndarray:
        r"""Meshless bare-:class:`numpy.ndarray` arm — the pointwise multiply.

        :math:`(M[f]\,\phi)_g(\vec r) = f_g(\vec r)\,\phi_g(\vec r)` on a bare
        ``(ng, *spatial)`` scalar/group block — the diagonal action with NO
        ordinate axis and NO mesh. This is the SAME per-block multiply the
        :class:`FullField` arm broadcasts over ordinates (single source of
        truth: ``self.coefficient.values``), so the two arms agree on every
        ordinate (pinned by the cross-arm consistency gate). Preserved for the
        cross-model outer-iteration consumers (homogeneous / diffusion /
        depletion) that feed bare arrays.
        """
        return self.coefficient.values * np.asarray(phi)

    if TYPE_CHECKING:
        # Honest per-carrier typing surface (mirrors FissionOperator, #257 S8c):
        # the public ``apply`` IS the runtime dispatcher (``apply = _apply_impl``
        # below), so callers statically see the exact output type per carrier.
        @overload
        def apply(self, psi: FullField, /) -> "FullField": ...
        @overload
        def apply(self, phi: np.ndarray, /) -> np.ndarray: ...
        def apply(self, psi: Any, /) -> Any: ...
    else:
        apply = _apply_impl

    def solve(self, q: "FullField") -> "FullField":
        r"""Inverse action :math:`M[f]^{-1}\,q = q / f = M[1/f]\,q`.

        Requires invertibility (the spectrum law: :math:`\min|f| > 0`);
        the engine raises
        :class:`~orpheus.numerics.operator.NotInvertible` otherwise.
        The codomain returns to a flux (the inverse of "flux → source" is
        "source → flux") in the input's family (the same two-arm parse as
        :meth:`apply`): an angular bulk returns
        :class:`~orpheus.transport.fields.angular_flux.AngularFlux` ⊕ zero
        :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`;
        a scalar bulk (#290 P4) returns
        :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux` ⊕ zero
        :class:`~orpheus.transport.fields.scalar_boundary_flux.ScalarBoundaryFlux`.
        """
        from orpheus.transport.fields._bases import AngularField, ScalarField
        from orpheus.transport.fields.angular_flux import AngularFlux
        from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
        from orpheus.transport.fields.scalar_flux import ScalarFlux
        from orpheus.transport.fields.scalar_boundary_flux import ScalarBoundaryFlux

        # Same #289 seam parse as :meth:`apply`; the mesh comes off the
        # parsed bulk's family-typed declaration.
        bulk = q.interior
        if isinstance(bulk, AngularField):
            # CS4b S4 — the space route (see apply).
            out_bulk = self.engine.solve(bulk.values)
            return FullField(
                interior=AngularFlux(values=out_bulk, space=bulk.space),
                boundary=AngularBoundaryFlux.zeros(q.boundary.space),
            )
        if isinstance(bulk, ScalarField):
            # Scalar arm (#290 P4): the typed division q/f through the
            # ONE engine (degenerate 1-ordinate lift — see apply), which
            # gates the spectrum law exactly as on the angular arm.
            out_bulk = self.engine.solve(bulk.values[None])[0]
            return FullField(
                interior=ScalarFlux(values=out_bulk, space=bulk.space),
                boundary=ScalarBoundaryFlux.zeros(q.boundary.space),
            )
        raise TypeError(
            f"MultiplicationOperator composite solve: angular- or "
            f"scalar-family bulk required; got {type(bulk).__name__}."
        )

    @overload
    def apply_transpose(self, psi: FullField, /) -> "FullField": ...
    @overload
    def apply_transpose(self, phi: np.ndarray, /) -> np.ndarray: ...
    def apply_transpose(self, psi: "Any") -> "Any":
        r"""Adjoint action :math:`M[f]^{*}\,\psi = M[\bar f]\,\psi = M[f]\,\psi`.

        Equal to :meth:`apply` — a real-valued multiplier is self-adjoint
        (:math:`M[f]^{*} = M[\bar f] = M[f]`), so ``M.H == M`` (the
        metric-blind Euclidean transpose; domain ``None``). Dispatches on the
        carrier exactly like :meth:`apply` (:class:`FullField` or bare ndarray).
        """
        return self.apply(psi)
