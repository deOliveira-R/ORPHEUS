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
  so :math:`M[f]` is invertible (``CAP_SOLVE``) **iff** :math:`f` is
  bounded away from zero (:math:`\min|f| > 0`); its inverse is
  :math:`M[1/f]`.

Collision as a named instance
=============================

The §5.7 named instance is the collision multiplier
:math:`C = M[\sigma_t]` — a plain :class:`MultiplicationOperator` carrying
the total cross-section field as its coefficient. (#261 retired the former
``CollisionOperator`` thin subclass: it added nothing the base lacked once
the base gained the optional :attr:`space` for the W-D composition guard;
the ``L + C → InvertibleOperator`` sweep dispatch lives on the SN-specific
:class:`~orpheus.sn.operator.StreamingOperator`, keyed on this base type —
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
``engine.apply(psi.bulk.values) == f.values[None] * psi.bulk.values``.
This operator owns only the *typed codomain* (the field wrapping); the
arithmetic lives once, in the engine (``coding-elegance`` Pattern 2).

Mesh-free (the carrier carries the mesh)
========================================

:class:`MultiplicationOperator` stores ONLY the coefficient field; the
mesh is read off the carrier at apply time (``mesh = psi.bulk.mesh``),
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
from typing import TYPE_CHECKING

import numpy as np

from orpheus.numerics.operator import (
    BlockRole,
    DiagonalOperator,
    LinearOperator,
)

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace
    from orpheus.transport.fields.cross_section_field import CrossSectionField
    from orpheus.transport.full_field import FullField


__all__ = ["MultiplicationOperator"]


@dataclass(eq=False)
class MultiplicationOperator(LinearOperator["FullField"]):
    r"""The promotion :math:`M[f]` of a coefficient field to a diagonal operator.

    Stores ONLY the coefficient; the mesh is read off the carrier at
    apply time. The raw pointwise multiply is delegated to the N-D
    :class:`~orpheus.numerics.operator.DiagonalOperator` broadcast engine
    on ``psi.bulk.values``; this operator owns the typed codomain.

    See the module docstring for the multiplier-algebra law-suite
    (:math:`M[1]=I`, :math:`M[0]=0`, linearity, homomorphism,
    self-adjointness, :math:`\mathrm{spec} = \mathrm{ess\,range}`).

    Parameters
    ----------
    coefficient : CrossSectionField
        The coefficient field :math:`f` (units ``1/cm`` for a cross
        section). Its ``.values`` is the ``(ng, *spatial)`` per-cell
        per-group array the broadcast engine consumes.

    Capabilities
    ------------

    ``{CAP_APPLY, CAP_APPLY_TRANSPOSE}`` always; adds ``CAP_SOLVE``
    **iff** :math:`\min|f| > 0` (the spectrum law: invertible iff the
    coefficient is bounded away from zero). This is the honest capability
    gate (``coding-elegance`` Pattern 4); the legacy collision operator
    advertised ``CAP_SOLVE`` unconditionally and produced silent IEEE
    NaN on a zero entry — :class:`MultiplicationOperator` revokes
    ``CAP_SOLVE`` instead, so a downstream composer never hits a broken
    inverse at call time.
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
    capabilities: frozenset[str] = field(default=frozenset(), init=False)

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
        # It already encodes "CAP_SOLVE iff every entry != 0", so inheriting
        # its capability set single-sources the spectrum law
        # spec(M[f]) = ess-range(f) (Pattern 4 — the transport operator and
        # the numerics engine agree by construction, not by a copied
        # predicate). Stored, not rebuilt per call: a fresh engine would
        # re-run the O(size) `np.all(coeff != 0)` spectrum scan on every
        # matvec and discard a capability already frozen here.
        self.engine = DiagonalOperator(
            self.coefficient.values, broadcast_axes=(0,),
        )
        self.capabilities = self.engine.capabilities

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
            else CrossSectionField.from_mesh(np.asarray(sigma), mesh)
        )
        # The composite space defaults to the mesh's own (an SNMesh carries
        # ``full_field_space``); pass ``space=`` to override. This makes
        # ``from_mesh(σ, sn_mesh)`` a faithful drop-in for the retired
        # ``CollisionOperator(sn_mesh, σ)``, which reached the same space via
        # ``sn_mesh.full_field_space`` (#261).
        if space is None:
            space = getattr(mesh, "full_field_space", None)
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

    # ── The §5.7 promotion: f ↦ M[f] on the leading ordinate axis ────────

    def apply(self, psi: "FullField") -> "FullField":
        r"""Forward action :math:`M[f]\,\psi = f \cdot \psi` on the composite.

        The pointwise multiply :math:`f\,\psi` is per-cell per-group,
        broadcast across every ordinate (the engine's leading-axis
        broadcast). The codomain is a *source* — multiplying a flux by a
        cross section yields a collision-rate density (units gain
        :data:`~orpheus.numerics.units.ANGULAR_RATE_UNITS`), so the bulk
        is an
        :class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`;
        the boundary is the implicit-zero
        :class:`~orpheus.transport.source_sinks.boundary_source_sink.BoundarySourceSink`
        (a multiplier has no face-trace action — the cell-balance
        :math:`f\,\psi` term is a CELL quantity).
        """
        from orpheus.transport.full_field import FullField
        from orpheus.transport.source_sinks import (
            AngularSourceSink,
            BoundarySourceSink,
        )

        mesh = psi.bulk.mesh
        out_bulk = self.engine.apply(psi.bulk.values)
        return FullField(
            bulk=AngularSourceSink.from_mesh(out_bulk, mesh),
            boundary=BoundarySourceSink.zeros_on(mesh),
        )

    def solve(self, q: "FullField") -> "FullField":
        r"""Inverse action :math:`M[f]^{-1}\,q = q / f = M[1/f]\,q`.

        Requires ``CAP_SOLVE`` (the spectrum law: :math:`\min|f| > 0`);
        the engine raises
        :class:`~orpheus.numerics.operator.MissingCapability` otherwise.
        The codomain returns to a flux (the inverse of "flux → source" is
        "source → flux"), so the bulk is an
        :class:`~orpheus.transport.fields.angular_flux.AngularFlux` and the
        boundary the implicit-zero
        :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`.
        """
        from orpheus.transport.fields.angular_flux import AngularFlux
        from orpheus.transport.fields.boundary_flux import BoundaryFlux
        from orpheus.transport.full_field import FullField

        mesh = q.bulk.mesh
        out_bulk = self.engine.solve(q.bulk.values)
        return FullField(
            bulk=AngularFlux.from_mesh(out_bulk, mesh),
            boundary=BoundaryFlux.zeros_on(mesh),
        )

    def apply_transpose(self, psi: "FullField") -> "FullField":
        r"""Adjoint action :math:`M[f]^{*}\,\psi = M[\bar f]\,\psi = M[f]\,\psi`.

        Equal to :meth:`apply` — a real-valued multiplier is self-adjoint
        (:math:`M[f]^{*} = M[\bar f] = M[f]`), so ``M.H == M`` (the
        metric-blind Euclidean transpose; domain ``None``).
        """
        return self.apply(psi)
