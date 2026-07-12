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
guard and the role-keyed birth factories (:meth:`from_mesh` /
:meth:`source_zeros_on` / :meth:`source_from_angular` — all presence-gated by
the leaf factories).

**Role-erased slots (B.2b DP2, the FullField precedent).** The static parameters
bind the locus FIELD BASES, not the flux leaves — exactly as ``FullField`` binds
``BulkField``/``BoundaryField`` — so ONE composite class carries the flux state
(the iterate), the source emission (an operator ``.apply`` output: ``A_BA`` /
``B_b``), and the displacement (minted per block by ``⊖``). Role identity lives
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

* ``.claude/plans/coupled_block_operator_campaign.md`` — Phase B (pose System B
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

    from orpheus.sn.mesh.augmented_mesh import SNMesh

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
    (the marched cells — flux state, source emission, or, after
    ``composite ⊖ composite``, an interior displacement); ``boundary`` a
    :class:`~orpheus.transport.fields._bases.RadialCharacteristicBoundaryField`
    (the r = R corner). Role-erased slots (the FullField precedent — see the
    module docstring): role identity lives on the members. The 2-block algebra
    (``±``, scalar ``·``, ``to_flat`` / ``from_flat``, ``copy``, the affine
    torsor propagated to both blocks) is inherited whole from
    :class:`~orpheus.transport.full_field.Composite`.
    """

    def __post_init__(self) -> None:
        # System B narrows the generic slots to the ψ½ split loci (the concrete
        # guard belongs with the concrete specialization, as FullField guards
        # BulkField / BoundaryField). Guard the FIELD BASE, not a role leaf, so
        # flux, source (operator emissions), and displacement (from ⊖)
        # composites are all admitted — role-erased slots, B.2b DP2.
        if not isinstance(self.interior, RadialCharacteristicInteriorField):
            raise TypeError(
                f"{type(self).__name__}: interior must be a "
                f"RadialCharacteristicInteriorField (its flux / source / "
                f"displacement leaf); got {type(self.interior).__name__}"
            )
        if not isinstance(self.boundary, RadialCharacteristicBoundaryField):
            raise TypeError(
                f"{type(self).__name__}: boundary must be a "
                f"RadialCharacteristicBoundaryField (its flux / source / "
                f"displacement leaf); got {type(self.boundary).__name__}"
            )
        super().__post_init__()  # base: Field + mesh-identity.

    # ── Construction ─────────────────────────────────────────────────

    @classmethod
    def from_mesh(cls, mesh: "SNMesh") -> "RadialCharacteristicField":
        r"""A zero ψ½ flux composite on ``mesh`` (presence-gated).

        Builds the two zero flux leaves from ``mesh``'s split spaces — so on a
        NON-carrying mesh (no R12a seed level) the leaf factories raise, i.e.
        System B is unconstructable exactly where it does not exist (presence =
        block existence).
        """
        return cls(
            interior=RadialCharacteristicInteriorFlux.zeros_on(mesh),
            boundary=RadialCharacteristicBoundaryFlux.zeros_on(mesh),
        )

    @classmethod
    def require_member(
        cls, x: object, *, mesh: "SNMesh", context: str,
    ) -> "RadialCharacteristicField":
        r"""Parse ``x`` as a System-B member carrier on ``mesh`` — the shared
        block-boundary guard.

        Every System-B block boundary (``A_BB`` / ``A_AB`` / ``B_b``) receives
        a :class:`RadialCharacteristicField` and must refuse (i) a foreign
        carrier type and (ii) a foreign mesh (the mesh-identity invariant —
        the member's legs, the operator's coefficients, and the radial widths
        must not desync). One parse, three consumers (coding-elegance
        Pattern 2; parse-don't-validate at the boundary). ``context`` names
        the refusing surface (``"RadialCharacteristicOperator.apply"``) so the
        error reads at the call site.
        """
        if not isinstance(x, cls):
            raise TypeError(
                f"{context}: expected a RadialCharacteristicField "
                f"(System B's member carrier); got {type(x).__name__}."
            )
        if x.mesh is not mesh:
            raise ValueError(
                f"{context}: the input field and the operator must share the "
                f"same SNMesh instance (mesh-identity invariant); got field "
                f"mesh {x.mesh!r} vs operator mesh {mesh!r}."
            )
        return x

    # ── The source-role birth factories ───────────────────────────────

    @classmethod
    def source_zeros_on(cls, mesh: "SNMesh") -> "RadialCharacteristicField":
        r"""A zero q½ SOURCE composite on ``mesh`` (presence-gated).

        The source-role sibling of :meth:`from_mesh` — the buffer the joint
        matvec/transpose surfaces fill with System B's emitted rows (an
        operator ``.apply`` output is source-role, never flux-role).  On a
        non-carrying mesh the leaf factories raise (presence = block
        existence).
        """
        return cls(
            interior=RadialCharacteristicInteriorSourceSink.zeros_on(mesh),
            boundary=RadialCharacteristicBoundarySourceSink.zeros_on(mesh),
        )

    @classmethod
    def source_from_angular(
        cls, angular_source_values: "NDArray", mesh: "SNMesh",
    ) -> "RadialCharacteristicField | None":
        r"""Fold a per-ordinate volumetric source into its q½ composite.

        The ONE source-side birth factory of #282 route (a) (Pattern 2 —
        the solver cold-starts, the fixed-source rhs, and the operator-free
        :func:`~orpheus.sn.loss_representation.transport_sweep` all route
        through here): ``None`` on a non-carrying mesh; on a carrying mesh
        (1-D curvilinear, R12a) the starting-direction legs receive the
        value of the source at the starting direction :math:`\mu = \pm 1`,
        reconstructed from ALL its Legendre moments (Hébert Eq. 3.432, the
        R14 full :math:`(-1)^\ell` fold):

        .. math::

           \bar q_{1/2}(\mu = \pm 1)
             \;=\; \sum_\ell \tfrac{2\ell+1}{2}\,q_\ell\,(\pm 1)^\ell,
           \qquad
           q_\ell(r) \;=\; \sum_n w_n\,P_\ell(\mu_n)\,q_n(r),

        via the R14 helper
        (:func:`~orpheus.numerics.spaces.radial_characteristic_space.fold_moments_to_radial_characteristic`).
        The full fold is REQUIRED for an anisotropic source: even an
        isotropic trial flux :math:`\psi = A(r)` streams to a
        :math:`\mu`-linear source :math:`q = \mu A'(r) + \sigma_t A(r)`,
        whose value at :math:`\mu = -1` is :math:`\sigma_t A - A'` — the
        :math:`\ell = 1` term carries the :math:`-A'` that an
        :math:`\ell = 0`-only fold drops (which floored the anisotropic
        curvilinear MMS; #282 route (a)).  For an isotropic source the
        higher moments vanish and the fold collapses to
        :math:`\tfrac12 q_0` bit-exactly (so the isotropic eigenvalue /
        fixed-source paths are unchanged).  The boundary corners stay zero:
        the inflow-corner datum is the BOUNDARY's job (vacuum ⇒ 0;
        reflective ⇒ the ``B`` corner arm into the SI rhs).

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
        """
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
        seed = cls.source_zeros_on(mesh)
        for p in seed.interior.space.levels:
            ords = np.asarray(level_indices[p])
            mu_p = mu[ords]
            w_p = weights[ords]
            q_p = vals[ords]                                  # (M_p, ng, nx)
            # Legendre moments of the source over the level's μ-nodes:
            # q_ℓ = Σ_n w_n P_ℓ(μ_n) q_n, ℓ = 0 … M_p−1 (the full angular
            # content the level resolves; the fold reconstructs q(μ=±1)).
            legendre = np.polynomial.legendre.legvander(mu_p, ords.size - 1)
            moments = np.einsum("n,nl,ngx->lgx", w_p, legendre, q_p)
            for sign in _SIGNS:
                seed.interior.cells(p, sign)[...] = (
                    fold_moments_to_radial_characteristic(moments, sign)
                )
        return seed
