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
:meth:`source_zeros` / :meth:`source_from_angular` — all presence-gated at
the member-space parse, CS4b S5).

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

import numpy as np

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
    from numpy.typing import NDArray

    from orpheus.numerics.space import FunctionSpace
    from orpheus.numerics.spaces.full_field_space import FullFieldSpace
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.transport.source_sinks.angular_boundary_source_sink import (
        AngularBoundarySourceSink,
    )

__all__ = ["RadialCharacteristicField"]

#: The ψ½ direction legs, in canonical (DAG) order: inward first, then the
#: pole-continued outward leg (the same order the split spaces lay out).
_SIGNS: tuple[int, int] = (-1, +1)


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
        :attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.radial_characteristic_field_space`,
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
        if space.interior_space is None or space.trace_space is None:
            raise ValueError(
                f"{context}: the ψ½ member space must carry BOTH blocks "
                f"(interior_space / trace_space); got {space!r}."
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

    @classmethod
    def source_from_angular(
        cls, angular_source_values: "NDArray", mesh: "SNMesh",
        *, boundary_trace: "AngularBoundarySourceSink | None" = None,
    ) -> "RadialCharacteristicField | None":
        r"""Fold a per-ordinate volumetric source into its q½ composite.

        The ONE source-side birth factory of #282 route (a) (Pattern 2 —
        the solver cold-starts, the fixed-source rhs, and the eigen-finalize
        reconstruction all route
        through here): ``None`` on a non-carrying mesh; on a carrying mesh
        (1-D curvilinear, R12a) the starting-direction legs receive the
        value of the source AT the level's start ray, synthesized from
        the level's own orthogonal analysis — one shape, two families:

        * **SPHERE** (the polar interval): the Hébert Eq. 3.432 Legendre
          fold at :math:`\mu = \pm 1` (the R14 full :math:`(-1)^\ell`
          fold),

          .. math::

             \bar q_{1/2}(\mu = \pm 1)
               \;=\; \sum_\ell \tfrac{2\ell+1}{2}\,q_\ell\,(\pm 1)^\ell,
             \qquad
             q_\ell(r) \;=\; \sum_n w_n\,P_\ell(\mu_n)\,q_n(r),

          via the R14 helper
          (:func:`~orpheus.numerics.spaces.radial_characteristic_space.fold_moments_to_radial_characteristic`).
        * **FOLDED CYLINDER** (the arc, Q5.6): the same analysis-synthesis
          in the arc's OWN Gauss family — the staggered arc is
          Gauss–Chebyshev-1 in :math:`x = \cos\omega = \eta/\sin\theta_p`
          (T25), so the endpoint synthesis uses
          :math:`T_k(\pm 1) = (\pm 1)^k` with the weight-free GC1
          discrete-orthogonality analysis
          :math:`c_k = \tfrac{2 - \delta_{k0}}{M}\sum_n T_k(x_n)\,q_n`.
          `[M]` the sphere's Legendre fold run on these levels produced
          seed values of −0.58 / +3.72 on a flat unit problem (an 82 %
          solution error, Mode-12-masked in every two-spelling test);
          the arc fold restores the flat L0 to 2.8e-13.
        The full fold is REQUIRED for an anisotropic source: even an
        isotropic trial flux :math:`\psi = A(r)` streams to a
        :math:`\mu`-linear source :math:`q = \mu A'(r) + \sigma_t A(r)`,
        whose value at :math:`\mu = -1` is :math:`\sigma_t A - A'` — the
        :math:`\ell = 1` term carries the :math:`-A'` that an
        :math:`\ell = 0`-only fold drops (which floored the anisotropic
        curvilinear MMS; #282 route (a)).  For an isotropic source the
        higher moments vanish and the fold collapses to
        :math:`\tfrac12 q_0` bit-exactly (so the isotropic eigenvalue /
        fixed-source paths are unchanged).

        **The inflow-corner law has THREE arms**, one per BC family of the
        outer face — each datum enters through its own system channel:

        * **vacuum** ⇒ 0 (this factory's zero corners, no trace supplied);
        * **reflective** ⇒ the ``B_b`` corner GAIN arm — iterate-dependent,
          added into the SI rhs by the splitting's :math:`N` (never here);
        * **prescribed inflow** ⇒ the SOURCE's own given datum, delivered
          HERE via ``boundary_trace``: the seed equation's r = R Dirichlet
          is :math:`\bar\psi_{1/2}(R) = \psi_{\rm in}(\mu=-1, R)`, and the
          composite rhs's boundary member carries :math:`\psi_{\rm in}` at
          the quadrature nodes, so each carrying level's inflow corner
          receives the trace's MOST-INWARD-ordinate ``xmax`` value — the
          nearest-node proxy for :math:`\mu = -1` (the pre-#282-route-(a)
          ``bc_outer_value = inflow[most_inward]`` datum, restored through
          the source channel; its :math:`O(\Delta\mu)` error sits below the
          #229 angular floor that governs every curvilinear gate; a
          half-range Legendre fold to :math:`\mu = -1` is the named
          upgrade seam if a sub-floor consumer ever appears).

        Omitting ``boundary_trace`` (or passing a zero trace — every
        vacuum/reflective rhs) leaves the corners zero, byte-identical to
        the pre-widening factory. The regression this arm closes: the d3
        carve dropped the prescribed-corner datum when it retired the
        lagged ``bc_outer_value`` read, and the sphere prescribed-inflow
        MMS converged to a wrong fixed point (L2 ≈ 0.21 vs the 2.4e-3
        floor) until step 7 restored it.

        Parameters
        ----------
        angular_source_values : NDArray
            The per-ordinate source in principled 1-D ``(N, ng, nx)``
            layout (carrying meshes are 1-D curvilinear).
        mesh : SNMesh
            The phase-space carrier (its
            ``radial_characteristic_field_space`` is the R12a presence
            predicate; its ``pole_angular_closure.level_indices`` give each
            level's ordinate bundle for the per-level moment integration).
        boundary_trace : trace source-sink, optional
            The SAME composite rhs's boundary member (an
            ``AngularBoundarySourceSink`` — duck-typed on ``face_view``):
            the prescribed-inflow given data whose ``xmax`` face feeds the
            inflow corners as above. ``None`` (the default) for callers
            whose source carries no boundary member (cold starts, the
            eigen fission seed).
        """
        from orpheus.geometry import CoordSystem
        from orpheus.numerics.spaces.radial_characteristic_space import (
            fold_moments_to_radial_characteristic,
        )

        if mesh.radial_characteristic_field_space is None:
            return None
        vals = np.asarray(angular_source_values)
        if vals.ndim != 3:
            raise ValueError(
                "RadialCharacteristicField.source_from_angular expects "
                f"the principled 1-D (N, ng, nx) per-ordinate layout; got "
                f"shape {vals.shape} (carrying meshes are 1-D curvilinear, "
                f"R12a)."
            )
        mu = mesh.quad.mu_x
        weights = mesh.quad.weights
        level_indices = mesh.pole_angular_closure.level_indices
        assert mesh.reduced is not None  # carrying ⇒ curvilinear ⇒ reduced
        arc_family = mesh.reduced.coord is CoordSystem.CYLINDRICAL
        seed = cls.source_zeros(mesh.radial_characteristic_field_space)
        for p in seed.interior.space.levels:
            ords = np.asarray(level_indices[p])
            mu_p = mu[ords]
            q_p = vals[ords]                                  # (M_p, ng, nx)
            if arc_family:
                # The ARC family (Q5.6): a folded level's nodes are the
                # Gauss-Chebyshev-1 angles in x = cos ω = η/sinθ_p (T25 —
                # the staggered arc IS the GC1 rule), so the source's value
                # at the start ray ω = π (x = −1) is the level's
                # trigonometric interpolant synthesized at the endpoint —
                # the exact arc analogue of the sphere's Legendre fold,
                # with T_k(±1) = (±1)^k in place of P_ℓ(±1) = (±1)^ℓ and
                # the GC1 discrete orthogonality
                # Σ_n T_j(x_n) T_k(x_n) = (M/2)(1+δ_{k0}) δ_{jk}
                # making the analysis weight-free:
                # c_k = (2 − δ_{k0})/M · Σ_n T_k(x_n) q_n.  For a level-
                # constant source c_0 = q and the fold returns q exactly.
                mu_z0 = float(mesh.quad.mu_z[ords[0]])
                sin_theta = float(np.sqrt(1.0 - mu_z0 * mu_z0))
                x_p = mu_p / sin_theta
                cheb = np.polynomial.chebyshev.chebvander(
                    x_p, ords.size - 1
                )
                coeffs = np.einsum("nk,ngx->kgx", cheb, q_p) * (
                    2.0 / ords.size
                )
                coeffs[0] *= 0.5
                k = np.arange(ords.size)
                for sign in _SIGNS:
                    seed.interior.cells(p, sign)[...] = np.einsum(
                        "k,kgx->gx", np.float64(sign) ** k, coeffs
                    )
            else:
                w_p = weights[ords]
                # Legendre moments of the source over the level's μ-nodes:
                # q_ℓ = Σ_n w_n P_ℓ(μ_n) q_n, ℓ = 0 … M_p−1 (the full
                # angular content the level resolves; the fold reconstructs
                # q(μ=±1)).
                legendre = np.polynomial.legendre.legvander(
                    mu_p, ords.size - 1
                )
                moments = np.einsum("n,nl,ngx->lgx", w_p, legendre, q_p)
                for sign in _SIGNS:
                    seed.interior.cells(p, sign)[...] = (
                        fold_moments_to_radial_characteristic(moments, sign)
                    )
        if boundary_trace is not None:
            # The prescribed-inflow corner arm (three-arm law above): the
            # r = R inward-leg Dirichlet datum per carrying level = the
            # trace's most-inward (most negative μ) ordinate on the outer
            # face. Carrying meshes are 1-D curvilinear (R12a), whose one
            # boundary face is "xmax" (the pole r = 0 is not a face).
            inflow_at_outer = np.asarray(boundary_trace.face_view("xmax"))
            for p in seed.interior.space.levels:
                ords = np.asarray(level_indices[p])
                # ords is stored eta-ascending (the level's SORT
                # CONTRACT), so the most-inward member IS ords[0] — the
                # one tie-break spelling. (An argmin re-derivation here
                # was #326's "second independent tie-break": equivalent
                # only by argmin's first-min convention.)
                most_inward = ords[0]
                seed.boundary.corner(p, -1)[...] = inflow_at_outer[most_inward, :]
        return seed
