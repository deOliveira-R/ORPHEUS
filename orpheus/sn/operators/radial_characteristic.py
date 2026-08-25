r"""The ψ½ block operators — all four blocks of the 2×2 coupled system.

This module hosts the FOUR curvilinear-SN ψ½ block operators of the augmented
within-group system, re-partitioned as two systems::

    [ A_AA   A_AB ] [ transport ]
    [ A_BA   A_BB ] [ ray       ]

* ``A_BB`` = :class:`RadialCharacteristicOperator` — the System-B self-block,
  the radial straight-characteristic march (documented below);
* ``A_AB`` = :class:`RadialCharacteristicSeeding` — the ray→bulk seed injection
  (cell-local angular; see its class docstring);
* ``A_BA`` = :class:`RadialCharacteristicEmission` — the bulk→ray coupling, the
  emission :math:`\mathrm{Fold} \circ K_{\rm iso} \circ \mathrm{integrate}` a
  within-group gain lags (see its class docstring);
* the shared factor ``Fold`` = :class:`RadialCharacteristicReconstruction` —
  reconstructs a bulk moment source at the closed :math:`\mu = \pm 1` rays.

``A_AA`` (the transport bulk⊕trace self-block) has no explicit object — it is
the driver-level composite ``L + C − S − B``.  The 2×2 posing, the augmented
composite, and the direct-seed route (a) are the record at
``docs/theory/methods/sn/curvilinear_one_group.rst §sn-direct-seed-solve``; the
remainder of THIS docstring documents ``A_BB``'s contract.

**What it is.** :class:`RadialCharacteristicOperator` is the banded radial
transport operator :math:`\mu\,\partial_r + \sigma_t` acting on the
starting-direction flux :math:`\psi_{1/2}` at the closed :math:`\mu = \pm 1`
rays of a curvilinear μ-level (Hébert §3.9.4).  At those rays the angular
redistribution vanishes (:math:`\alpha_{1/2}=0`), so the balance decouples
from the α-cascade to a plain DD recurrence in radius — a *straight
characteristic* (the physics:
``docs/theory/methods/sn/curvilinear_one_group.rst §sn-direct-seed-pole-straight-characteristic``).
Domain = codomain = System B's member space
(the composite the carrying mesh mints as
``radial_characteristic_field_space``, BOUND at construction — un-weld
arc O-1; the carrier of
:class:`~orpheus.transport.radial_characteristic_field.RadialCharacteristicField`);
it exists ONLY on a seed-carrying mesh (the sphere, R12a).

**How it is posed and solved.** System B carries TWO boundary conditions —
r = R Dirichlet (the outer-face :math:`\mu=-1` inflow corner is GIVEN data:
0 for vacuum, the reflected outflow corner for reflective, fed by the
boundary block ``B_b``
:class:`~orpheus.sn.operators.boundary.RadialCharacteristicBoundaryOperator`)
and r = 0 pole continuation (:math:`\psi_{1/2}^{+}(0)=\psi_{1/2}^{-}(0)`,
internal to the march) — a well-posed two-point radial BVP.  Because the
recurrence is triangular in radius (the #284 forward-substitution
certificate), the two-leg Carlson march
(:func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`,
Hébert 3.434–3.435) is the EXACT direct inverse :math:`A_{BB}^{-1}` — no
iteration.  :meth:`solve` is that march (a thin WRAP of the standalone
engine — single source, Cardinal Rule 2 — the SAME engine the in-sweep block
:mod:`orpheus.sn.loss_representation` runs), and :meth:`solve_transpose` its
Euclidean reverse-mode adjoint
(:func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_transpose`).
The BVP posing and the block-triangular certificate are the record at
``§sn-direct-seed-solve`` / ``§sn-direct-seed-block-triangular`` in
``docs/theory/methods/sn/curvilinear_one_group.rst``.

Realized surface — the full operator
====================================

This operator realizes BOTH directions of the complete ``A_BB`` machinery:
the **forward** :math:`A_{BB}\,\psi_{1/2}` as :meth:`apply` and its Euclidean
transpose :math:`A_{BB}^{\mathsf T}` as :meth:`apply_transpose`; the
**resolvent** :math:`A_{BB}^{-1}` as :meth:`solve` (the direct Carlson march)
and its adjoint as :meth:`solve_transpose`; and the operator-returning
:meth:`inverse` (a generic
:class:`~orpheus.numerics.operator.InverseOperator` whose ``apply`` IS
:meth:`solve`), so ``is_invertible`` and ``is_adjointable`` both read
``True``.

**Single source (Cardinal Rule 2 — no twin).** The forward and its transpose
wrap ONE shared kernel each
(:func:`~orpheus.sn.sweep.psi_half_angle_seed.radial_characteristic_forward_residual`
/ ``…_transpose``, the PURE ``A_BB`` transpose — the ``A_AB`` seed→bulk
coupling's transpose is the explicit
:meth:`RadialCharacteristicSeeding.apply_transpose`, not part of ``A_BB``);
and the two-leg solve ORCHESTRATION (source views → inward leg →
pole-continue → reversed outward leg → flux views) lives ONLY here — the
production ``(L+C)`` walk routes its ray solve THROUGH :meth:`solve` /
:meth:`solve_transpose`, the former in-sweep inline block retired.  Do NOT
re-implement the march at a call site; construct this operator over the local
:math:`\sigma_t` and call it.  The routing is pinned by the Mode-11
engine-execution sentinels (the S1 driver gates); the un-weave record is at
``docs/theory/methods/sn/curvilinear_one_group.rst §sn-direct-seed-solve``.

Sourcing
========

The assembly supplies the ray carrier (``field_space``), the radial widths
``dr``, and the per-level start cosines (:func:`march_start_cosines`) — all
read off the carrying mesh AT the assembly, bound as values here (un-weld
arc O-1); ``total_cross_section`` is the :math:`\sigma_t`
:class:`~orpheus.transport.fields.cross_section_field.CrossSectionField` — the
collision coefficient :math:`C_{\rm ray}` the march reads (``.values`` are
``(ng, nx)``). It is a TYPED, mesh-bound field (not a bare array) so the
mesh-identity invariant — σ_t and the operator march the SAME radial widths —
is enforced at construction (Pattern 4; the sibling
:class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` guards its
diagonal coefficient the same way). Step 4 sources it from the solver's
macroscopic XS as ``mat_xs.total_cross_section_field``.

References
==========

* Hébert, A. (2009). *Applied Reactor Physics*. §3.9.4, Eqs. (3.432)–(3.435).
* GH #282 (the direct ψ½ seed / route (a)), #280 (the walk unification),
  #284 (the forward-substitution / triangular certificate).
* ``.claude/plans/archive/coupled_block_operator_campaign.md`` — the 2×2 posing,
  rulings P1/P2, and this step's design reconnaissance.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional, Protocol

import numpy as np

from orpheus.numerics.operator import (
    InverseOperator,
    LinearOperator,
    NotInvertible,
    SystemRole,
)
from orpheus.numerics.spaces.radial_characteristic_space import (
    RadialCharacteristicInteriorSpace,
    fold_moments_to_radial_characteristic,
    fold_moments_to_radial_characteristic_transpose,
)
from orpheus.sn.sweep.psi_half_angle_seed import (
    carlson_inward_sweep_from_source,
    carlson_inward_sweep_transpose,
    radial_characteristic_forward_residual,
    radial_characteristic_forward_residual_transpose,
)

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping

    from orpheus.geometry import CoordSystem
    from orpheus.geometry.reduced_operator import ReducedStreamingOperator
    from orpheus.numerics.quadrature.directional import Quadrature
    from orpheus.numerics.space import FunctionSpace
    from orpheus.numerics.spaces.full_field_space import FullFieldSpace
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.transport.fields.cross_section_field import CrossSectionField
    from orpheus.transport.full_field import FullField
    from orpheus.transport.radial_characteristic_field import (
        RadialCharacteristicField,
    )


def march_start_cosines(
    reduced: "ReducedStreamingOperator", levels: "Iterable[int]",
) -> dict[int, float]:
    r"""``|η_start|`` per carried level — the march ray's radial cosine.

    The DD engine
    (:func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`)
    marches in PATH length; every consumer divides the radial widths by
    this cosine so the optical step is :math:`\Delta r\,\sigma/|\eta|` —
    the ray at polar angle θ traverses :math:`\Delta r` of radius over
    :math:`\Delta r/\sin\theta` of path.  Sphere: the diameter ray,
    :math:`|\mu| = 1` (the historical hard-coded case — dividing by 1.0
    is bit-exact, so every sphere path is byte-identical to the
    pre-Q5.6 spelling).  Folded cylinder: level p's ξ = 0 start ray at
    :math:`|\eta| = \sin\theta_p`.  Single-sourced from the reduced
    operator's own start-direction fields (``mu_start`` /
    ``mu_start_per_level``), never recomputed from the quadrature.
    Public since the un-weld arc (O-1): the assembly reads
    ``(sn_mesh.reduced, sn_mesh.radial_characteristic_levels)`` and hands
    the computed map to :class:`RadialCharacteristicOperator` — the
    operator binds the VALUES; this producer owns the provenance ruling.
    """
    lvls = tuple(levels)
    if reduced.mu_start_per_level is not None:
        return {p: abs(float(reduced.mu_start_per_level[p])) for p in lvls}
    assert reduced.mu_start is not None  # spherical arm always sets it
    diameter_ray = abs(float(reduced.mu_start))
    return {p: diameter_ray for p in lvls}


class _EmissionKernel(Protocol):
    r"""The isotropic :math:`\ell = 0` emission kernel — ``A_BA``'s factor between
    the angular integral and the fold: an ``ndarray → ndarray`` operator carrying
    a Euclidean transpose.

    A **structural** contract (no nominal base): the bare
    :class:`~orpheus.numerics.operator.LinearOperator` Protocol does NOT surface
    ``apply_transpose`` — that is a runtime capability of *adjointable* operators
    only (#276 P3; the same reason
    :attr:`~orpheus.transport.operators.scattering.ScatteringOperator.isotropic_kernel`
    is annotated as the concrete ``OperatorSum``, not ``LinearOperator``, so the
    checker sees the transpose). ``A_BA`` needs BOTH directions, so it types its
    kernel by the capability it consumes. Satisfied by the scattering
    ``isotropic_kernel`` (``K_iso`` = an ``OperatorSum``) and the fission
    ``kernel`` (the rank-1 ``χ ⊗ νΣf`` dyad).
    """

    def apply(self, phi: "np.ndarray", /) -> "np.ndarray": ...

    def apply_transpose(self, chi: "np.ndarray", /) -> "np.ndarray": ...


__all__ = [
    "RadialCharacteristicOperator",
    "RadialCharacteristicSeeding",
    "RadialCharacteristicReconstruction",
    "RadialCharacteristicEmission",
]


class RadialCharacteristicOperator(LinearOperator["RadialCharacteristicField"]):
    r"""``A_BB`` — the radial straight-characteristic transport operator on ψ½.

    System B's self-block of the 2×2 coupled block operator: the banded
    radial DD recurrence :math:`\mu\,\partial_r + \sigma_t` at the closed
    :math:`\mu = \pm 1` rays (Hébert §3.9.4). Endomorphic on System B's
    member carrier
    :class:`~orpheus.transport.radial_characteristic_field.RadialCharacteristicField`
    (domain = codomain = the bound System-B member composite; the
    ``CoupledOperator`` grid places it at the ``(B, B)`` slot).  All four
    action surfaces parse the composite at the block boundary and march the
    split member views directly.

    See the module docstring's "Realized surface" section for the
    operator-algebra posing (the two-point radial BVP) and the full realized
    surface — the forward :meth:`apply` / :meth:`apply_transpose`, the
    resolvent :meth:`solve` / :meth:`solve_transpose` (the direct Carlson
    march IS :math:`A_{BB}^{-1}`), and the operator-returning :meth:`inverse`
    (so ``is_invertible`` and ``is_adjointable`` both read ``True``).

    Parameters
    ----------
    field_space : FullFieldSpace or None
        System B's member composite (interior ⊕ boundary corner) — the
        endomorphic domain/codomain. ``None`` (a seedless pose — Cartesian /
        non-carrying cylinder, R12a) is REFUSED: no System B, no A_BB.
    total_cross_section : CrossSectionField
        The total cross-section :math:`\sigma_t` as a typed
        :class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`
        (``.values`` are ``(ng, nx)``) — the collision coefficient
        :math:`C_{\rm ray}` the march reads. Must be posed on ``bulk_space``
        (the pose-consistency invariant — refused otherwise) and strictly
        positive everywhere (the DD denominator :math:`\Delta r\,\sigma + 2`
        is then well-defined and the march never emits NaN).
    bulk_space : FunctionSpace
        The scalar bulk σ_t must be posed on (the admission reference —
        consulted at construction, not stored).
    dr : np.ndarray
        The radial cell widths Δr — the geometry commitment's spatial
        measure along the ray axis.
    start_cosines : Mapping[int, float]
        ``|η_start|`` per carried level — produce with
        :func:`march_start_cosines` (single-sourced from the reduced
        operator's start-direction fields; never recompute from the
        quadrature).
    """

    def __init__(
        self,
        field_space: "FullFieldSpace | None",
        total_cross_section: "CrossSectionField",
        *,
        bulk_space: "FunctionSpace",
        dr: np.ndarray,
        start_cosines: "Mapping[int, float]",
    ) -> None:
        if field_space is None:
            raise ValueError(
                "RadialCharacteristicOperator: the pose carries no "
                "starting-direction ray (the split ψ½ spaces are None) — a "
                "seedless mesh (Cartesian, or a cylinder whose levels start on "
                "an edge node or an η-tie, R12a) has no System B. A_BB "
                "exists only on a seed-carrying mesh — the GL sphere, the "
                "σ_y-folded cylinder (Q5.6)."
            )
        # Pose-consistency invariant (Pattern 4 — make the illegal state
        # unrepresentable): the march reads the radial widths ``dr`` against
        # σ_t, so a σ_t posed on a foreign bulk — even one with an identical
        # (ng, nx) — would march the wrong Δr. The typed CrossSectionField
        # carries its space, so we refuse the mismatch at construction (the
        # sibling StreamingCollisionOperator guards σ↔pose the same way,
        # streaming.py). The field's own space invariant then GUARANTEES
        # values.shape == (ng, nx) for this pose — so no separate shape check
        # is needed (an explicit one would be redundant ceremony).
        if total_cross_section.space != bulk_space:
            raise ValueError(
                "RadialCharacteristicOperator: total_cross_section must be a "
                "CrossSectionField agreeing with the operator's scalar "
                "bulk in content (space-content invariant); got field space "
                f"{total_cross_section.space!r} vs "
                f"{bulk_space!r}."
            )
        sigma = total_cross_section.values
        if np.any(sigma <= 0.0):
            raise ValueError(
                f"RadialCharacteristicOperator: total_cross_section must be "
                f"strictly positive everywhere for the DD march denominator "
                f"(Δr·σ + 2) to be well-defined; got "
                f"min(σ_t) = {float(np.min(sigma)):.3e}."
            )
        #: System B's member composite (interior ⊕ boundary corner) —
        #: domain AND codomain (endomorphic; the one space the mesh used to
        #: contribute, un-weld arc O-1).
        self._field_space = field_space
        #: The total cross-section :math:`\sigma_t` — the ``C_ray`` collision
        #: coefficient, a typed :class:`CrossSectionField` on the scalar bulk
        #: (``.values`` are ``(ng, nx)``; the pose-consistency guard above
        #: makes a foreign-pose σ_t unconstructable). The march reads ``.values``.
        self.total_cross_section = total_cross_section
        #: The radial cell widths Δr the march divides by |η_start| — the
        #: geometry commitment's spatial measure along the ray axis.
        self._dr = np.asarray(dr, dtype=float)
        #: |η_start| per carried level (see :func:`march_start_cosines`,
        #: the single-sourced producer the assembly calls).
        self._start_cosines = dict(start_cosines)
        interior_space = field_space.interior_space
        assert isinstance(interior_space, RadialCharacteristicInteriorSpace)  # System-B mint; narrowing only
        #: The ψ½ interior split space (levels/ng/nx metadata for repr;
        #: non-None by the ctor guard — presence-coextensive with the
        #: composite space).
        self._ray_space = interior_space

    # ── Spaces ────────────────────────────────────────────────────────
    #
    # Endomorphic on System B's member space (B.2c): both the source (input
    # to solve) and the flux (its output) live on the SAME composite space —
    # the role (source / flux) lives on the member field types, not the
    # space. The unified `_ray_space` stays the ENGINE's carrier (slot views
    # inside the march), never the declared block boundary.
    # block_role stays None (the base default): A_BB is System B's self-block
    # on the ray space, not a bulk/full/boundary block of the FullField
    # composite. Its two-system membership rides the SystemRole axis (campaign
    # step 4a) — A_BB acts within System B alone.
    system_role = SystemRole.B

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        # Non-None by the ctor guard (the composite space is presence-
        # coextensive with the unified one — both None exactly off-R12a).
        return self._field_space

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        return self._field_space

    # ── Predicates — the forward closes the involution web (step 4b) ──

    @property
    def is_invertible(self) -> bool:
        # A_BB's ``solve`` IS the exact direct inverse (the two-leg Carlson
        # march); with the forward ``apply`` realized (step 4b) the involution
        # web ``inverse().solve == apply`` closes, so ``inverse()`` is a
        # faithful capability (a carrying-mesh A_BB is always invertible).
        return True

    @property
    def is_adjointable(self) -> bool:
        # Both transposes are realized: the forward ``apply_transpose``
        # (Euclidean :math:`A_{BB}^{\mathsf T}`) and the resolvent
        # ``solve_transpose`` (:math:`(A_{BB}^{-1})^{\mathsf T}`).
        return True

    # ── Forward action — A_BB ψ = (μ ∂_r + σ_t)ψ (campaign step 4b) ────

    def apply(
        self, x: "RadialCharacteristicField", /,
    ) -> "RadialCharacteristicField":
        r"""The forward matvec :math:`A_{BB}\,\psi_{1/2} = (\mu\,\partial_r
        + \sigma_t)\,\psi_{1/2}` — the exact algebraic inverse of :meth:`solve`.

        A thin WRAP of the single-sourced
        :func:`~orpheus.sn.sweep.psi_half_angle_seed.radial_characteristic_forward_residual`
        (Cardinal Rule 2: one forward, no twin).  Reads the flux state ψ½
        (cells + corners) off the member composite; returns the residual
        source :math:`q_{1/2} = A_{BB}\,\psi_{1/2}` as a source-member
        composite.  The ``apply ∘ solve`` outflow-corner defect closes to
        ``0.0`` exactly; the cell round-trip is principled-equiv at ~FP ULP
        (the forward's :math:`2/\Delta r` and the march's :math:`\Delta r\,
        \sigma + 2` reassociate).

        Parameters
        ----------
        x : RadialCharacteristicField
            The ψ½ state as System B's member composite (member roles erased
            — the bridge preserves whatever role the members carry). Must
            share this operator's mesh.

        Returns
        -------
        RadialCharacteristicField
            The residual source :math:`A_{BB}\,\psi_{1/2}` (source members).
        """
        from orpheus.transport.radial_characteristic_field import (
            RadialCharacteristicField,
        )
        from orpheus.transport.source_sinks.radial_characteristic_boundary_source_sink import (
            RadialCharacteristicBoundarySourceSink,
        )
        from orpheus.transport.source_sinks.radial_characteristic_interior_source_sink import (
            RadialCharacteristicInteriorSourceSink,
        )

        comp = self._require_member_composite(x, "apply")
        interior_out, boundary_out = radial_characteristic_forward_residual(
            comp.interior.values, comp.boundary.values,
            comp.interior.space, comp.boundary.space,
            self.total_cross_section.values,
            self._dr,
            start_cosines=self._start_cosines,
        )
        return RadialCharacteristicField(
            interior=RadialCharacteristicInteriorSourceSink(
                values=interior_out, space=comp.interior.space,
            ),
            boundary=RadialCharacteristicBoundarySourceSink(
                values=boundary_out, space=comp.boundary.space,
            ),
        )

    def apply_transpose(
        self, y: "RadialCharacteristicField", /,
    ) -> "RadialCharacteristicField":
        r"""The Euclidean transpose of the forward matvec —
        :math:`A_{BB}^{\mathsf T}`.

        A thin WRAP of
        :func:`~orpheus.sn.sweep.psi_half_angle_seed.radial_characteristic_forward_residual_transpose`
        — the PURE ``A_BB`` transpose. The ``A_AB`` seed→bulk recurrence
        coupling's transpose is the explicit grid block
        (:meth:`RadialCharacteristicSeeding.apply_transpose`), NOT part of
        ``A_BB`` in isolation. This is the flat
        EUCLIDEAN adjoint :math:`A_{BB}^{\mathsf T}`; the metric Hilbert
        adjoint :math:`G^{-1}A_{BB}^{\mathsf T}G` (``.H``) is realized once at
        the composite (L19). Composite in, source-member composite out — the
        same block-boundary bridge as :meth:`apply`.
        """
        from orpheus.transport.radial_characteristic_field import (
            RadialCharacteristicField,
        )
        from orpheus.transport.source_sinks.radial_characteristic_boundary_source_sink import (
            RadialCharacteristicBoundarySourceSink,
        )
        from orpheus.transport.source_sinks.radial_characteristic_interior_source_sink import (
            RadialCharacteristicInteriorSourceSink,
        )

        comp = self._require_member_composite(y, "apply_transpose")
        interior_out, boundary_out = radial_characteristic_forward_residual_transpose(
            comp.interior.values, comp.boundary.values,
            comp.interior.space, comp.boundary.space,
            self.total_cross_section.values,
            self._dr,
            start_cosines=self._start_cosines,
        )
        return RadialCharacteristicField(
            interior=RadialCharacteristicInteriorSourceSink(
                values=interior_out, space=comp.interior.space,
            ),
            boundary=RadialCharacteristicBoundarySourceSink(
                values=boundary_out, space=comp.boundary.space,
            ),
        )

    def inverse(self) -> "InverseOperator":
        r"""Return :math:`A_{BB}^{-1}` as an
        :class:`~orpheus.numerics.operator.InverseOperator`.

        The generic solve-backed inverse (#226 taxonomy §13): the returned
        operator's ``apply`` IS :meth:`solve` (the direct Carlson march) and
        its ``solve`` is this operator's forward :meth:`apply` — the
        involution web, no reciprocal twin. ``A_BB`` earns only the generic
        ``InverseOperator`` (round-trip alone); the distinguishing triangular
        invariant lives in
        :class:`~orpheus.sn.operators.sweep_operator.SweepOperator`.
        """
        if not self.is_invertible:  # pragma: no cover — always invertible
            raise NotInvertible(
                "RadialCharacteristicOperator.inverse requires a working "
                "solve (the Carlson march); a carrying-mesh A_BB always has "
                "one."
            )
        return InverseOperator(self)

    # ── Resolvent action — the direct Carlson march IS A_BB⁻¹ ─────────

    def solve(
        self, source: "RadialCharacteristicField",
    ) -> "RadialCharacteristicField":
        r"""Solve :math:`A_{BB}\,\psi_{1/2} = q_{1/2}` by the two-leg Carlson march.

        The EXACT direct inverse :math:`A_{BB}^{-1}` (no iteration): per
        seed-carrying level, the inward :math:`\mu=-1` leg marches from the
        given r = R inflow corner to the pole, then the outward
        :math:`\mu=+1` leg rides the SAME engine on reversed cell data
        (orientation is carried by the DATA, never a flag) from the
        pole-continued face out to the r = R outflow corner. A thin WRAP of
        :func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`
        (Hébert 3.434–3.435) — the single-sourced DD engine.  The production
        ``(L+C)`` walk routes its ray solve through this method (the former
        in-sweep inline twin is retired — see the module docstring "Single
        source").

        The per-level slot key is the carrier's own ``space.levels`` member
        (the level POSITION, ``p_idx`` in the in-sweep) — the coordinate
        that keys the space slots, pinned by
        :mod:`tests.sn.mesh.test_radial_characteristic_slot_coordination`.

        Parameters
        ----------
        source : RadialCharacteristicField
            The q½ source as System B's member composite — the interior
            member's cells legs carry the folded starting-direction source
            and the boundary member's :math:`\mu=-1` corner the r = R inflow
            (Dirichlet) datum. Must carry SOURCE members (the role parse at
            the block boundary — the resolvent inverts :math:`A_{BB}\psi = q`,
            so a flux composite here is a caller error worth
            raising loudly) and share this operator's mesh.

        Returns
        -------
        RadialCharacteristicField
            The ψ½ state :math:`\psi_{1/2}` satisfying the two-point radial
            BVP (flux members) — cells + both corner legs filled per carried
            level.
        """
        from orpheus.transport.radial_characteristic_field import (
            RadialCharacteristicField,
        )
        from orpheus.transport.source_sinks.radial_characteristic_interior_source_sink import (
            RadialCharacteristicInteriorSourceSink,
        )

        comp = self._require_member_composite(source, "solve")
        if not isinstance(comp.interior, RadialCharacteristicInteriorSourceSink):
            raise TypeError(
                f"RadialCharacteristicOperator.solve: the source composite "
                f"must carry q½ SOURCE members (the role parse at the block "
                f"boundary, #289-F2 — the resolvent inverts A_BB·ψ = q); got "
                f"the {type(comp.interior).__name__} interior role."
            )
        sigma = self.total_cross_section.values
        dr = self._dr
        start_cosines = self._start_cosines

        flux = RadialCharacteristicField.flux_zeros(self._field_space)
        for level in comp.interior.space.levels:
            # The engine marches in PATH length: the start ray at radial
            # cosine |η_start| traverses Δr of radius over Δr/|η_start| of
            # path (sphere: |μ| = 1, unchanged bit-for-bit; folded-cylinder
            # level p rides its ξ = 0 ray at sinθ_p — Q5.6).
            dr_path = dr / start_cosines[level]
            q_minus = comp.interior.cells(level, -1)
            q_plus = comp.interior.cells(level, +1)
            corner_in = comp.boundary.corner(level, -1)
            # inward starting-direction leg: enter at the r=R inflow corner,
            # exit at the pole.
            cells_minus, pole_face = carlson_inward_sweep_from_source(
                q_minus, sigma, dr_path, corner_in,
            )
            # outward leg: the SAME engine on reversed data, entering at
            # the pole-continued face (ψ½⁺(0) = ψ½⁻(0)) and exiting at r=R.
            cells_plus_rev, corner_out = carlson_inward_sweep_from_source(
                q_plus[:, ::-1], sigma[:, ::-1], dr_path[::-1], pole_face,
            )
            flux.interior.cells(level, -1)[...] = cells_minus
            flux.boundary.corner(level, -1)[...] = corner_in
            flux.interior.cells(level, +1)[...] = cells_plus_rev[:, ::-1]
            # The outflow corner is the free-sign defect row (forward:
            # ``streamed − stored = q_out``), so the exact inverse emits
            # ``ψ_out = streamed − q_out`` — the same completion the
            # System-A sweep's outflow restore got (ERR-071; this is its
            # System-B twin, ERR-078). Every PHYSICAL rhs carries
            # q_out = 0 (the defect row reads 0 = 0 at any fixed point),
            # so production values are bit-unchanged; dropping the term
            # made the coupled composite's inverse silently non-inverse
            # on the per-level outflow slots.
            flux.boundary.corner(level, +1)[...] = (
                corner_out - comp.boundary.corner(level, +1)
            )
        return flux

    def solve_transpose(
        self, cotangent: "RadialCharacteristicField",
    ) -> "RadialCharacteristicField":
        r"""The Euclidean adjoint of :meth:`solve` — :math:`(A_{BB}^{-1})^{\mathsf T}`.

        The reverse-mode adjoint of the two-leg march: given a cotangent on
        the flux (the solve's codomain), return the cotangent on the source
        (its domain). Per level, the OUTWARD leg is reversed first (its exit
        corner feeds the pole-face cotangent), then the INWARD leg (the pole
        cotangent is its exit), threading the running face cotangent back to
        the r = R inflow corner — the transpose of :meth:`solve`'s leg chain,
        via
        :func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_transpose`.

        This is the ISOLATED ray-block transpose :math:`(A_{BB}^{-1})^{\mathsf T}`
        — the pure resolvent adjoint. The full ``(L+C)`` reverse-scan adds a
        seed→bulk term (the M-M thread cotangent) on the inward cells; that
        term is the ``A_AB`` coupling's transpose (campaign steps 2–3), NOT
        part of ``A_BB`` in isolation.

        **Duality typing (#276 A4):** the input is the dual of the solve's
        codomain — dual-of-flux, i.e. an adjoint ray SOURCE (a cotangent
        legitimately rides any role's carrier — member roles erased on
        input); the output is the dual of its domain — dual-of-source
        under the G-pairing, i.e. the adjoint ray FLUX (flux members).
        This is the ray-leg sibling of
        :meth:`~orpheus.sn.operators.streaming.StreamingCollisionOperator.solve_transpose`'s
        A4 re-classing: the daggered coupled fixed-point iteration closes
        over the SAME class pattern as the forward (flux iterate → source
        gains → transposed substitution → flux iterate), and the typed
        cross-class guard caught the pre-A4 source-family wrap on the
        daggered sphere's first contact.

        Parameters
        ----------
        cotangent : RadialCharacteristicField
            A cotangent on the flux composite (member roles erased — a
            cotangent legitimately rides any role's carrier). Must share this
            operator's mesh.

        Returns
        -------
        RadialCharacteristicField
            The domain cotangent = the adjoint ray FLUX (flux members).
            The :math:`\mu=+1` source corner carries
            :math:`-\bar y_{\rm out}` — the transpose of the defect
            row's :math:`-I` coupling (``ψ_out = streamed − q_out``,
            ERR-078); the q½ fold's own writes remain cells + the
            :math:`\mu=-1` corner.
        """
        from orpheus.transport.radial_characteristic_field import (
            RadialCharacteristicField,
        )

        comp = self._require_member_composite(cotangent, "solve_transpose")
        sigma = self.total_cross_section.values
        dr = self._dr
        start_cosines = self._start_cosines

        # Duality typing (#276 A4, docstring above): dual-of-source = the
        # adjoint ray FLUX — the flux-role zeros buffer.
        src_bar = RadialCharacteristicField.flux_zeros(self._field_space)
        for level in comp.interior.space.levels:
            # The transpose of the PATH-length march uses the same per-level
            # path widths as the forward (sphere: ÷1.0, byte-identical).
            dr_path = dr / start_cosines[level]
            cells_minus_bar = comp.interior.cells(level, -1)
            cells_plus_bar = comp.interior.cells(level, +1)
            corner_in_bar = comp.boundary.corner(level, -1).copy()
            corner_out_bar = comp.boundary.corner(level, +1)
            # reverse the OUTWARD (+1) leg — marched on reversed data; its exit
            # corner cotangent seeds the pole-face cotangent for the inward leg.
            q_plus_rev_bar, pole_face_bar = carlson_inward_sweep_transpose(
                cells_plus_bar[:, ::-1], corner_out_bar, sigma[:, ::-1],
                dr_path[::-1],
            )
            q_plus_bar = q_plus_rev_bar[:, ::-1]
            # reverse the INWARD (−1) leg — the pole face is its exit.
            q_minus_bar, corner_in_from_minus = carlson_inward_sweep_transpose(
                cells_minus_bar, pole_face_bar, sigma, dr_path,
            )
            # the r=R inflow corner both passes through to the flux corner AND
            # enters the inward leg — its cotangent is the sum of both paths.
            corner_in_bar = corner_in_bar + corner_in_from_minus
            src_bar.interior.cells(level, -1)[...] = q_minus_bar
            src_bar.interior.cells(level, +1)[...] = q_plus_bar
            src_bar.boundary.corner(level, -1)[...] = corner_in_bar
            # The defect row's −I coupling (ψ_out = streamed − q_out,
            # ERR-078): the outflow-flux cotangent feeds −1 into the
            # outflow-SOURCE cotangent, alongside its march-threading
            # path above.
            src_bar.boundary.corner(level, +1)[...] = -corner_out_bar
        return src_bar

    # ── Internals ─────────────────────────────────────────────────────

    def _require_member_composite(
        self, x: "RadialCharacteristicField", method: str,
    ) -> "RadialCharacteristicField":
        r"""Parse the member composite at the block boundary.

        Single source for the four action-surface guards: the shared
        :meth:`~orpheus.transport.radial_characteristic_field.RadialCharacteristicField.require_member`
        parse (carrier class + mesh-identity — the field's legs, this
        operator's :math:`\sigma_t`, and ``axis_widths`` cannot desync).
        Since 4e the engines march the split members directly — the
        parse returns the composite itself (the unified bridge retired).
        """
        from orpheus.transport.radial_characteristic_field import (
            RadialCharacteristicField,
        )

        return RadialCharacteristicField.require_member(
            x,
            space=self._field_space,
            context=f"RadialCharacteristicOperator.{method}",
        )

    def __repr__(self) -> str:
        return (
            f"RadialCharacteristicOperator(levels={self._ray_space.levels}, "
            f"ng={self._ray_space.ng}, nx={self._ray_space.nx})"
        )


class RadialCharacteristicSeeding(
    LinearOperator["RadialCharacteristicField", "FullField"],
):
    r"""``A_AB`` — the ray→bulk seed injection (the Morel–Montry angular seed).

    The off-diagonal ``(transport, ray)`` block of the 2×2 coupled block
    operator: the ψ½ starting-direction ray seeds the bulk angular recurrence.
    Domain = System B's member space
    ``sn_mesh.radial_characteristic_field_space`` (the operator reads the
    inward :math:`\mu=-1` ``cells(p, -1)`` leg); codomain = System A's
    ``sn_mesh.full_field_space`` — :meth:`apply` emits the seed's bulk
    contribution as the interior member of a
    :class:`~orpheus.transport.full_field.FullField` over a zero trace.  It
    exists ONLY on a seed-carrying mesh (the sphere, R12a).

    **A CELL-LOCAL ANGULAR coupling.** At each radial cell :math:`i` the ray
    value :math:`\psi_{1/2}(i)` is the seed of the Morel–Montry **weighted**
    angular recurrence (Morel & Montry 1984; implemented form
    Bailey--Morel--Chang 2010 Eqs. (42)/(43) — Hébert §3.9.4 is the
    authority for the Carlson march that PRODUCES this seed, Eqs.
    3.432-3.435, not for the recurrence it feeds), run over ORDINATES at a
    FIXED cell; the
    upstream half-flux :math:`\psi_{m-1/2,\,i}` then enters that cell's balance
    as the angular numerator :math:`(\Delta A/w)\,c_{\rm in}\,\psi_{m-1/2,\,i}`.
    So the seed at cell :math:`i` couples ONLY to the bulk ordinates at the
    SAME cell — there is **no spatial cell-cell coupling** (the radial march of
    :math:`\psi_{1/2}` is ``A_BB``'s job).  This is what separates ``A_AB``
    (cell-local angular, both directions realized HERE) from ``A_BB``
    (spatially woven into ``(L+C).apply``).  The seeding formula and the
    augmented posing:
    ``docs/theory/methods/sn/curvilinear_one_group.rst §sn-direct-seed-solve``.

    **Realized as WRAPs of the shared M-M closure** — the SAME
    :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
    machinery the ``(L+C)`` matvec runs (single source — Cardinal Rule 2).
    ``A_AB`` isolates its own block by ZEROING the bulk (forward: the
    recurrence is linear in :math:`(\psi_{\rm bulk}, \psi_{1/2})` jointly, so a
    zero bulk leaves only the ``A_AB`` part) and DISCARDING the bulk cotangent
    (transpose).  The forward places the angular numerator as :math:`-(\Delta
    A/w)\,c_{\rm in}\,\psi_{m-1/2}/V`; the transpose reverses it with the local
    gather :math:`\bar n_{p,\,m,\,i} = -\bar o_{m,\,i}/V_i` then
    ``angular_adjoint`` to the seed cotangent on the ``cells(p, -1)`` leg — the
    Euclidean transpose :math:`A_{AB}^{\mathsf T}`.  Both directions realized ⇒
    ``is_adjointable = True``; ``is_invertible`` inherits ``False``
    (rectangular ray → bulk).

    The shared kernel serves BOTH ``A_AA`` (bulk redistribution, ``psi_view ≠
    0``) and ``A_AB`` (the seed, ``psi_view = 0``) — never a twin.  A tracked
    transient twin remains at a DIFFERENT production entry (the in-sweep
    placement fused into the ``(L+C)`` matvec,
    :mod:`orpheus.sn.loss_representation`), retired when steps 4/5 route the
    ``(L+C)`` bulk rows through the ``CoupledOperator`` block matvec; until then
    both sides are behaviour-pinned (the in-sweep by the regression floor +
    sweep suite, this by the bit-identity gates ``TestA_AB_SeedInjection``).

    Parameters
    ----------
    sn_mesh : SNMesh
        The augmented geometry — seed-carrying (1-D curvilinear, R12a). Supplies
        the ray carrier (the domain), the M-M closure ``pole_angular_closure``
        (the single-sourced kernel), the cell volumes ``volumes``, and the
        quadrature ``quad``. A seedless mesh (Cartesian, or a non-carrying
        cylinder) has NO ray→bulk coupling: constructing over one is rejected. Unlike ``A_BB``,
        ``A_AB`` needs NO :math:`\sigma_t` — with the bulk zeroed the
        collision/streaming terms drop out and only the σ-independent angular
        numerator survives.
    """

    # A_AB maps System B (the ray seed) → System A (the bulk residual): an
    # off-diagonal block, so it spans both systems (campaign step 4a).
    system_role = SystemRole.COUPLED

    def __init__(self, sn_mesh: "SNMesh") -> None:
        space = sn_mesh.radial_characteristic_interior_space
        if space is None:
            raise ValueError(
                "RadialCharacteristicSeeding: the mesh carries no "
                "starting-direction ray (the split ψ½ spaces are None) "
                "— a seedless mesh (Cartesian, or a non-carrying cylinder, "
                "R12a) has no System B, hence no ray→bulk seed coupling to "
                "inject. A_AB exists only on a seed-carrying mesh — the GL "
                "sphere, the σ_y-folded cylinder (Q5.6)."
            )
        #: The augmented geometry (ray carrier + the M-M closure + volumes).
        self.sn_mesh = sn_mesh
        #: The ψ½ interior split space (level metadata for the gather loops).
        #: The declared DOMAIN is the member composite space (B.2c) — see
        #: :attr:`domain`.
        self._ray_space = space

    # ── Predicates / spaces ───────────────────────────────────────────

    @property
    def is_adjointable(self) -> bool:
        # Both directions are realized: :meth:`apply` (seed → bulk residual
        # injection) and :meth:`apply_transpose` (bulk cotangent → ray seed
        # cotangent, the seed_cells_bar term). A_AB is cell-local angular (no
        # spatial weaving), so — unlike A_BB's forward matvec — neither
        # direction defers. is_invertible inherits False: a rectangular
        # coupling (ray → bulk) is not square.
        return True

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        # System B's member space — the ψ½ seed composite (the input to
        # :meth:`apply`). Non-None by the ctor guard (presence-coextensive
        # with the unified engine carrier).
        return self.sn_mesh.radial_characteristic_field_space

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        # System A's composite carrier (:meth:`apply` emits the seed's bulk
        # contribution as a FullField interior member — B.2c).
        return self.sn_mesh.full_field_space

    # ── Forward — seed → bulk angular-numerator contribution ──────────

    def apply(self, seed: "RadialCharacteristicField", /) -> "FullField":
        r"""Inject the ψ½ ray seed into the bulk angular recurrence.

        Build the seed-only Morel–Montry half-angle grid
        (``precompute_psi_state`` with an all-zero ``psi_view`` — the zero
        bulk isolates ``A_AB`` from ``A_AA`` by linearity), then per carrying
        level and cell place the angular numerator :math:`-(\Delta A/w)\,
        c_{\rm in}\,\psi_{m-1/2}/V` (the seed's term in :math:`m =
        (\mathrm{denom}\cdot\psi - \mathrm{numer})/V` with :math:`\psi_{\rm
        cell}=0`).  Non-carrying-level ordinates stay zero.  Bit-identical to
        the in-sweep injection; the bulk contribution lands as the FullField
        INTERIOR member over a zero trace.

        Parameters
        ----------
        seed : RadialCharacteristicField
            The ψ½ seed as System B's member composite (member roles erased —
            the bridge preserves them). Only the interior member's inward
            ``cells(p, -1)`` leg is read (the recurrence seed); the ``+1``
            leg and corners are ``A_BB``-internal. Must share this operator's
            mesh.

        Returns
        -------
        FullField
            The seed's contribution to System A's residual row — the
            per-ordinate bulk term ``(N, ng, nx)`` as the interior member
            over a zero trace.
        """
        from orpheus.transport.full_field import FullField
        from orpheus.transport.radial_characteristic_field import (
            RadialCharacteristicField,
        )
        from orpheus.transport.source_sinks import (
            AngularBoundarySourceSink,
        )
        from orpheus.transport.source_sinks.angular_source_sink import (
            AngularSourceSink,
        )

        mesh = self.sn_mesh
        rc_space = mesh.radial_characteristic_field_space
        assert rc_space is not None  # ctor guard: carrying mesh; narrowing only
        comp = RadialCharacteristicField.require_member(
            seed, space=rc_space, context="RadialCharacteristicSeeding.apply",
        )
        closure = mesh.pole_angular_closure
        space = self._ray_space
        N = mesh.quad.N
        ng = space.ng
        nx = mesh.nx
        V = mesh.volumes
        level_indices = closure.level_indices

        # Bulk zeroed → the recurrence propagates ONLY the seed (linearity).
        # The closure reads the INTERIOR member (its ``.cells(p, -1)`` seed).
        psi_state = closure.precompute_psi_state(
            np.zeros((N, ng, nx)), radial_characteristic=comp.interior,
        )
        out_g_first = np.zeros((ng, N, nx))
        for p in space.levels:
            ordinates = np.asarray(level_indices[p])
            within = np.arange(ordinates.size)
            for i in range(nx):
                _, upstream_numer = closure.cell_contribution(
                    psi_state, i, p, within,
                )
                out_g_first[:, ordinates, i] = -upstream_numer / V[i]
        return FullField(
            interior=AngularSourceSink(values=out_g_first.swapaxes(0, 1), space=mesh.angular_bulk_space),
            boundary=AngularBoundarySourceSink.zeros(mesh.angular_trace),
        )

    # ── Euclidean transpose — bulk cotangent → ray seed cotangent ─────

    def apply_transpose(
        self, cotangent: "FullField", /,
    ) -> "RadialCharacteristicField":
        r"""Euclidean transpose :math:`A_{AB}^{\mathsf T}` — System A → ray seed.

        The adjoint of :meth:`apply`: given a cotangent on System A's row
        space (the codomain), return the cotangent on the ray seed composite
        (the domain). The forward writes ONLY the interior member, so the
        transpose reads ONLY ``cotangent.interior`` — the trace part is
        annihilated structurally (discarded, not zeroed). Reverse the
        forward's :math:`-\,\mathrm{numer}/V` placement with the local gather
        :math:`\bar n_{p,\,m,\,i} = -\bar o_{m,\,i}/V_i`, then
        ``angular_adjoint`` reverses the M-M recurrence to the per-carrying-
        level seed cotangent ``seed_cells_bar[p]``, written on the inward
        ``cells(p, -1)`` leg. The bulk-redistribution cotangent
        (``psi_ang_bar``, ``A_AA``'s share of the shared kernel) is discarded.
        The ``+1`` leg and corners stay zero (the forward writes only the
        inward leg).

        Parameters
        ----------
        cotangent : FullField
            A cotangent on System A's composite (member roles erased — only
            the interior member ``(N, ng, nx)`` is read). Must share this
            operator's mesh.

        Returns
        -------
        RadialCharacteristicField
            The ray seed cotangent (source members) — the interior member's
            ``cells(p, -1)`` filled per carrying level, everything else zero.
        """
        from orpheus.transport.full_field import FullField
        from orpheus.transport.radial_characteristic_field import (
            RadialCharacteristicField,
        )

        mesh = self.sn_mesh
        if not isinstance(cotangent, FullField):
            raise TypeError(
                f"RadialCharacteristicSeeding.apply_transpose: expected a "
                f"FullField (System A's composite — the codomain cotangent); "
                f"got {type(cotangent).__name__}."
            )
        self._check_mesh(cotangent, "apply_transpose")
        closure = mesh.pole_angular_closure
        V = mesh.volumes
        level_indices = closure.level_indices

        ob_g_first = cotangent.interior.values.swapaxes(0, 1)  # (ng, N, nx)
        # numer_bar for EVERY level (angular_adjoint needs the full tuple);
        # the local −ō/V gather is the exact transpose of the forward −·/V.
        numer_bar = tuple(
            -ob_g_first[:, np.asarray(level_idx), :] / V[None, None, :]
            for level_idx in level_indices
        )
        _, seed_cells_bar = closure.angular_adjoint(numer_bar)

        src_bar = RadialCharacteristicField.source_zeros(mesh.radial_characteristic_field_space)
        for p, cells_bar in seed_cells_bar.items():
            src_bar.interior.cells(p, -1)[...] = cells_bar
        return src_bar

    # ── Internals ─────────────────────────────────────────────────────

    def _check_mesh(self, field: "FullField", method: str) -> None:
        r"""Enforce the space-content invariant (Pattern 4; CS4b S3).

        The input composite, this operator's ``pole_angular_closure``, and
        the ``volumes`` must all agree with ONE
        :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh`'s content, so the
        seed legs, the M-M coefficients, and the ``/V`` scaling cannot
        desync. Compared per BLOCK (the composite's own ``==`` is
        name+shape and cannot see blocks).
        """
        if field.interior.space != field.interior.space_on(self.sn_mesh):
            raise ValueError(
                f"RadialCharacteristicSeeding.{method}: the input composite "
                f"must agree with the operator mesh in space content "
                f"(space-content invariant); got interior space "
                f"{field.interior.space!r} on operator mesh {self.sn_mesh!r}."
            )

    def __repr__(self) -> str:
        return (
            f"RadialCharacteristicSeeding(levels={self._ray_space.levels}, "
            f"ng={self._ray_space.ng}, nx={self._ray_space.nx})"
        )


class RadialCharacteristicReconstruction(LinearOperator):
    r"""Reconstruct a bulk angular source at the closed μ=±1 rays — the ``A_BA`` fold.

    The shared factor of the ``A_BA`` bulk→ray coupling block
    (``A_BA = RadialCharacteristicReconstruction ∘ Emission``, where the full
    coupling is :class:`RadialCharacteristicEmission`): given a bulk within-group
    angular source in its Legendre-moment representation, it evaluates that
    source at the two closed radial-characteristic rays :math:`\mu = \pm 1` and
    writes the result as the q½ ray source on every carried μ-level.  It lives
    here, beside its ``A_BB`` / ``A_BA`` siblings, because all are blocks of the
    one ψ½ coupled operator (the fold *math* kernel
    :func:`~orpheus.numerics.spaces.radial_characteristic_space.fold_moments_to_radial_characteristic`
    stays at the numerics data layer).

    **What it is.** The 1-D angular reconstruction :math:`\mathcal R`
    (Hébert Eq. 3.432) SAMPLED at the closed rays,
    :math:`\bar q_{1/2}(\pm 1) = \sum_\ell (2\ell+1)/2\,q_\ell\,(\pm 1)^\ell`
    (:math:`P_\ell(\pm 1) = (\pm 1)^\ell`).  Its Euclidean transpose injects a
    ray cotangent back into moment space
    (:func:`~orpheus.numerics.spaces.radial_characteristic_space.fold_moments_to_radial_characteristic_transpose`)
    — so ``is_adjointable = True``; ``is_invertible`` inherits ``False`` (a
    reconstruction from ``n_moments`` to the two rays is not square).  Why the
    fold keeps ALL Legendre moments (streaming manufactures angular structure
    the flux lacks): ``docs/theory/methods/sn/curvilinear_one_group.rst
    §sn-direct-seed-source-fold``.

    **Single source of the fold.** Both the forward reconstruction and its
    transpose route through the ONE math kernel
    :func:`~orpheus.numerics.spaces.radial_characteristic_space.fold_moments_to_radial_characteristic`
    (Cardinal Rule 2 — the :math:`P_1(-1)` sign is spelled once).  The
    production emitter (:class:`RadialCharacteristicEmission`) feeds the
    isotropic :math:`\ell = 0` emission (``n_moments = 1``), and the fold
    accepts any order so the anisotropic reach is testable before it is needed.
    (The distinct
    :meth:`~orpheus.transport.radial_characteristic_field.RadialCharacteristicField.source_from_angular`
    data-factory folds a *per-ordinate* source through the SAME kernel after a
    Legendre projection — a different typed input, not a twin.)

    **Broadcast across levels.** The same moment source is folded onto every
    carried level (exact for the isotropic emission — a level-independent
    value); the GL sphere carries ONE level, a σ_y-folded cylinder SEVERAL
    (R12a / Q5.6), and the transpose sums the per-level, per-sign
    cotangents either way.  Corners stay zero: the fold writes only the cells legs (the
    inflow-corner datum is the boundary block ``B_b``'s job).

    Parameters
    ----------
    field_space : FullFieldSpace or None
        System B's member composite — the codomain (its interior IS the
        ψ½ split space with the cells-leg layout). ``None`` (a seedless
        pose) is REFUSED: no bulk→ray coupling to fold.
    coord : CoordSystem
        The geometry commitment — selects the fold family (CYLINDRICAL =
        the arc family's ℓ = 0 reproducing weight; otherwise the sphere's
        Legendre kernel path).
    quadrature : Quadrature
        The angular rule whose total weight Σw spells the arc family's
        constant-reproducing fold weight ``1/Σw``.
    n_moments : int, default 1
        The operator's domain dimension — the number of angular Legendre
        moments it reconstructs from. ``1`` is the isotropic production reach
        (:math:`\ell = 0`); a larger order is the manufactured anisotropic
        case. Fixed at construction so the transpose has a well-defined
        codomain ``(n_moments, ng, nx)``.
    """

    # A_BA maps System A (the bulk emission) → System B (the ray source): an
    # off-diagonal block, so it spans both systems (campaign step 4a).
    system_role = SystemRole.COUPLED

    def __init__(
        self,
        field_space: "FullFieldSpace | None",
        *,
        coord: "CoordSystem",
        quadrature: "Quadrature",
        n_moments: int = 1,
    ) -> None:
        from orpheus.geometry import CoordSystem

        if field_space is None:
            raise ValueError(
                "RadialCharacteristicReconstruction: the pose carries no "
                "radial-characteristic ray (the split ψ½ spaces are None) "
                "— a seedless mesh (Cartesian, or a cylinder whose levels "
                "start on an edge node or an η-tie, R12a) has no bulk→ray "
                "coupling to fold."
            )
        space = field_space.interior_space
        assert isinstance(space, RadialCharacteristicInteriorSpace)  # System-B mint; narrowing only
        if n_moments < 1:
            raise ValueError(
                f"RadialCharacteristicReconstruction: n_moments must be ≥ 1 "
                f"(at least the ℓ = 0 moment); got {n_moments!r}."
            )
        arc_family = coord is CoordSystem.CYLINDRICAL
        if arc_family and n_moments > 1:
            raise NotImplementedError(
                f"RadialCharacteristicReconstruction: an arc-family "
                f"(folded-cylinder) mesh reconstructs from the ℓ = 0 "
                f"moment only (n_moments = 1); got {n_moments}. The "
                f"Legendre reconstruction weights are the polar "
                f"interval's — a folded arc's higher moments live in the "
                f"Chebyshev basis (T25), a seam no consumer needs yet."
            )
        #: System B's member composite — the declared codomain (the one
        #: space the mesh used to contribute, un-weld arc O-1).
        self._field_space = field_space
        #: The angular Legendre-moment order of the domain (``1`` = isotropic,
        #: the production reach; larger is the manufactured anisotropic case).
        self.n_moments = n_moments
        #: The ψ½ interior split space (level metadata + fold-loop keys;
        #: non-None by the ctor guard). The declared codomain is the member
        #: composite space.
        self._ray_space = space
        #: The arc family's ℓ = 0 fold weight — the CONSTANT function's
        #: reproducing weight ``1/Σw`` (a level-constant moment's value at
        #: ANY direction is q₀/∫dμ; the sphere's Legendre (2·0+1)/2 = ½ is
        #: the same identity at Σw = 2).  ``None`` selects the sphere's
        #: Legendre kernel path.  Q5.6: the sphere's hard-coded ½ applied
        #: to a folded cylinder (Σw = 4π) over-injected the ray scattering
        #: source 2π-fold — `[M]` the c = 0.4 flat-flux equilibrium read
        #: 158 % off before this weight, machine-exact after.
        self._ell0_scale: float | None = (
            1.0 / float(quadrature.weights.sum()) if arc_family else None
        )

    # ── Predicates / spaces ───────────────────────────────────────────

    @property
    def is_adjointable(self) -> bool:
        # Both the forward fold and its Euclidean transpose are realized
        # (:meth:`apply` / :meth:`apply_transpose`) — the same shape as the
        # shared IsotropicScattering kernel. is_invertible inherits False (a
        # reconstruction from n_moments to the two rays is not square).
        return True

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        # The bulk angular-moment source — an untyped ``(n_moments, ng, nx)``
        # intermediate (like K_iso's ndarray domain); no minted moment space.
        return None

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        # System B's member space — the fold emits the q½ composite (4e).
        return self._field_space

    # ── Forward fold — reconstruct at μ=±1 ────────────────────────────

    def apply(self, moments: np.ndarray, /) -> "RadialCharacteristicField":
        r"""Reconstruct the bulk moment source at μ=±1 → the q½ ray source.

        Folds the Legendre-moment source ``moments`` (shape
        ``(n_moments, ng, nx)``) onto every carried level's cells legs at both
        closed rays via
        :func:`~orpheus.numerics.spaces.radial_characteristic_space.fold_moments_to_radial_characteristic`
        (the single-source fold). The corners stay zero — the fold is
        volumetric (the inflow-corner datum is ``B_b``'s job).

        Parameters
        ----------
        moments : np.ndarray
            The bulk within-group source in ℓ-leading Legendre moments,
            shape ``(n_moments, ng, nx)``. The production emitter passes the
            isotropic emission with a unit ℓ axis (``emission[None]``,
            ``n_moments=1``).

        Returns
        -------
        RadialCharacteristicField
            The q½ ray source composite — interior cells legs folded,
            boundary corners zero.
        """
        from orpheus.transport.radial_characteristic_field import (
            RadialCharacteristicField,
        )

        arr = np.asarray(moments)
        expected = (self.n_moments, self._ray_space.ng, self._ray_space.nx)
        if arr.shape != expected:
            raise ValueError(
                f"RadialCharacteristicReconstruction.apply expects a bulk "
                f"Legendre-moment source of shape (n_moments, ng, nx) = "
                f"{expected}; got {arr.shape}."
            )
        seed = RadialCharacteristicField.source_zeros(self._field_space)
        for level in self._ray_space.levels:
            for sign in (-1, +1):
                if self._ell0_scale is not None:
                    # Arc family: the ℓ = 0 constant's reproducing weight
                    # 1/Σw — sign-independent (a constant's value at any
                    # direction), guarded to n_moments = 1 at the ctor.
                    seed.interior.cells(level, sign)[...] = (
                        arr[0] * self._ell0_scale
                    )
                else:
                    seed.interior.cells(level, sign)[...] = (
                        fold_moments_to_radial_characteristic(arr, sign)
                    )
        return seed

    # ── Euclidean transpose — inject a ray cotangent into moment space ─

    def apply_transpose(
        self, cotangent: "RadialCharacteristicField", /,
    ) -> np.ndarray:
        r"""Euclidean transpose — a ray cotangent → the bulk moment cotangent.

        The adjoint of :meth:`apply` with respect to the moments: sum the
        per-level, per-sign ray-cells cotangents expanded back onto the
        moments via
        :func:`~orpheus.numerics.spaces.radial_characteristic_space.fold_moments_to_radial_characteristic_transpose`
        (the SAME reconstruction weights the forward contracts — the sign is
        spelled once). This is the single source of the scattering seed
        adjoint (``∂S/∂ψ½`` cotangent → the ℓ = 0 bulk moment, which
        ``K_isoᵀ`` then scatters into the bulk) — the transpose factor of
        :meth:`RadialCharacteristicEmission.apply_transpose`.

        Parameters
        ----------
        cotangent : RadialCharacteristicField
            A cotangent on the ray source composite (member roles erased).
            Must share this operator's mesh. Only the interior cells legs
            are read (the forward writes only cells); the boundary corner
            cotangents are ignored.

        Returns
        -------
        np.ndarray
            The bulk moment cotangent, shape ``(n_moments, ng, nx)``.
        """
        if cotangent.interior.space != self._ray_space:
            raise ValueError(
                f"RadialCharacteristicReconstruction.apply_transpose: the "
                f"cotangent composite must agree with the operator's ray "
                f"interior in space content (space-content invariant); got "
                f"interior space {cotangent.interior.space!r} vs "
                f"{self._ray_space!r}."
            )
        moment_bar = np.zeros(
            (self.n_moments, self._ray_space.ng, self._ray_space.nx),
            dtype=float,
        )
        for level in self._ray_space.levels:
            for sign in (-1, +1):
                cells_bar = cotangent.interior.cells(level, sign)
                if self._ell0_scale is not None:
                    # Arc family: the exact transpose of the ℓ = 0
                    # reproducing-weight forward (n_moments = 1 by ctor).
                    moment_bar[0] += cells_bar * self._ell0_scale
                else:
                    moment_bar += (
                        fold_moments_to_radial_characteristic_transpose(
                            cells_bar, sign, self.n_moments,
                        )
                    )
        return moment_bar

    def __repr__(self) -> str:
        return (
            f"RadialCharacteristicReconstruction("
            f"n_moments={self.n_moments}, levels={self._ray_space.levels}, "
            f"ng={self._ray_space.ng}, nx={self._ray_space.nx})"
        )


class RadialCharacteristicEmission(LinearOperator):
    r"""``A_BA`` — the bulk→ray coupling (the ψ½ emission ``Fold ∘ K ∘ integrate``).

    The off-diagonal ``(ray, transport)`` block of the 2×2 coupled block
    operator: the bulk within-group flux emits an isotropic source that seeds
    the ψ½ starting-direction ray. It composes THREE factors —

    .. math::

        A_{BA} \;=\; \underbrace{\mathrm{Fold}}_{\text{Reconstruction}}
                 \;\circ\; \underbrace{K}_{\text{emission kernel}}
                 \;\circ\; \underbrace{\textstyle\int\! d\mu}_{\text{integrate}} ,

    the angular integral of the bulk flux to :math:`\phi_0`, the operator's
    isotropic emission kernel :math:`K` (the scattering
    ``K_{\rm iso} = \Sigma_{s0} + 2\Sigma_{2n}`` — the production instantiation),
    and the reconstruction of that :math:`\ell = 0` moment at the closed
    :math:`\mu = \pm 1` rays (:class:`RadialCharacteristicReconstruction`). It
    exists ONLY on a seed-carrying mesh (the sphere, R12a).

    **Generic over the emission kernel.** What distinguishes a bulk→ray
    emission coupling is only its emission kernel :math:`K` — an ``ndarray →
    ndarray`` operator (``(ng, nx) → (ng, nx)``) with ``apply`` /
    ``apply_transpose``.  The SCATTERING coupling passes
    :attr:`~orpheus.transport.operators.scattering.ScatteringOperator.isotropic_kernel`
    (``K_iso`` — the production consumer, a within-group lagged gain); the
    fission dyad
    :attr:`~orpheus.transport.operators.fission.FissionOperator.kernel`
    (``χ ⊗ νΣf``) is a smoke-verified SECOND kernel that exercises the same
    machinery.

    .. warning::

       **HAZARD — do NOT route fission's PRODUCTION seed through this operator.**
       Fission is the eigenvalue OUTER source: its ``K ∘ integrate`` is
       pre-computed as the fission source, so the ray seed is a DIRECT
       :class:`RadialCharacteristicReconstruction` fold at the ``q_ext`` seam.
       Routing ``F.kernel`` through this operator's full
       ``Fold ∘ K ∘ integrate`` would DOUBLE-apply ``K ∘ integrate``.  S is a
       within-group gain, F the outer source — different seams.  The
       kernel-genericity is a clean dependency injection, NOT a claim that
       fission wires through here.

    **A within-group lagged gain, a true System A → System B block.** The SI
    driver lags this coupling as its own gain (the Wave-O #208 pattern that
    separated ``B`` from ``S``): the within-group scattering gain rides ``(S,
    A_BA, B)``, while fission's ``A_BA`` rides the OUTER ``q_ext`` (within-group
    fission is zero).  ``apply`` is typed ``FullField → RadialCharacteristic
    Field`` — it reads System A and returns System B's OWN carrier with SOURCE
    members, no present-zero bulk padding, so the old "A_BA writes into the
    bulk" double-count is UNSPELLABLE (Pattern 4).  Since B.2d the production
    SI/Krylov driver consumes this block natively at the ``(B, A)`` slot of the
    within-group gain grid
    (:func:`orpheus.sn.coupled_system.build_within_group_system`).  The lift
    that moved the ψ½ seed off the model-generic ``S``/``F`` and onto this
    first-class coupling: ``docs/theory/methods/sn/curvilinear_one_group.rst
    §sn-direct-seed-solve``.

    **The transpose IS the S-adjoint bulk pullback (single source).**
    :meth:`apply_transpose` reads a ray cotangent and pulls it back into the
    bulk as :math:`\int\!d\mu^{\mathsf T}\,K^{\mathsf T}\,\mathrm{Fold}^{\mathsf
    T}` — exactly the ``w · K_isoᵀ(Reconstructionᵀ χ_seed)`` term the S-adjoint
    carries.  ``S.apply_transpose`` is pure-bulk and the composite ``(L+C − S −
    A_BA − B).H`` reconstructs the monolithic adjoint.  There is NO adjoint SI
    driver in production, so the pullback is reached ONLY by the ``.H``
    reciprocity gates — which is why :meth:`apply_transpose` is realized NOW and
    pinned by a NONZERO-seed-cotangent gate (a present-zero seed would hide a
    lost pullback).  ``is_adjointable = True``; ``is_invertible`` inherits
    ``False`` (rectangular bulk→ray).

    Parameters
    ----------
    emission_kernel : LinearOperator
        The operator's isotropic :math:`\ell = 0` emission kernel :math:`K`,
        an ``ndarray → ndarray`` map ``(ng, nx) → (ng, nx)`` with
        ``apply``/``apply_transpose``. In production pass
        ``scattering_op.isotropic_kernel`` (``K_iso``) — the SAME shared object
        the bulk scattering gain uses, so the emission is single-sourced
        (computed once in ``S``'s bulk arm and once here — one shared kernel
        call, not a twin). ``fission_op.kernel`` is accepted (the machinery is
        kernel-generic) but fission's production ray seed rides the outer
        ``q_ext`` seam as a direct ``Fold``, NOT through this operator (see
        the class docstring — routing it here would double-apply
        ``K ∘ integrate``).
    field_space : FullFieldSpace or None
        System B's member composite — the codomain. ``None`` (a seedless
        pose) is REFUSED: no System B, no bulk→ray emission coupling.
    full_field_space : FullFieldSpace
        System A's composite carrier — the domain (bulk ⊕ trace).
    angular_bulk_space, angular_trace : FunctionSpace
        The transpose's output spaces (per-ordinate bulk pullback over a
        zero trace).
    quadrature : Quadrature
        The angular rule whose weights spell the ``(∫dμ)ᵀ`` broadcast.
    coord : CoordSystem
        The geometry commitment, forwarded to the internal fold factor.
    """

    # A_BA maps System A (bulk) → System B (ray): an off-diagonal block spanning
    # both systems (campaign step 4a).
    system_role = SystemRole.COUPLED

    def __init__(
        self,
        emission_kernel: "_EmissionKernel",
        *,
        field_space: "FullFieldSpace | None",
        full_field_space: "FullFieldSpace",
        angular_bulk_space: "FunctionSpace",
        angular_trace: "FunctionSpace",
        quadrature: "Quadrature",
        coord: "CoordSystem",
    ) -> None:
        if field_space is None:
            raise ValueError(
                "RadialCharacteristicEmission: the pose carries no "
                "radial-characteristic ray (the split ψ½ spaces are None) "
                "— a seedless mesh (Cartesian, or a non-carrying cylinder, "
                "R12a) has no System B, hence no bulk→ray emission "
                "coupling. A_BA exists only on a seed-carrying mesh — the "
                "GL sphere, the σ_y-folded cylinder (Q5.6)."
            )
        #: The isotropic ℓ=0 emission kernel (K_iso for S, the fission dyad
        #: for F) — the SHARED ``ndarray → ndarray`` object the bulk gain uses.
        self.emission_kernel = emission_kernel
        #: System B's member composite — the declared codomain.
        self._field_space = field_space
        #: System A's composite carrier — the declared domain (bulk ⊕ trace).
        self._full_field_space = full_field_space
        #: The transpose's output spaces (per-ordinate bulk + trace) and the
        #: quadrature whose weights spell (∫dμ)ᵀ — bound at construction
        #: (un-weld arc O-1: the operator binds VALUES and SPACES, no mesh).
        self._angular_bulk_space = angular_bulk_space
        self._angular_trace = angular_trace
        self._quadrature = quadrature
        #: The fold factor (the migrated reconstruction, ℓ=0 production reach).
        self._fold = RadialCharacteristicReconstruction(
            field_space, coord=coord, quadrature=quadrature,
        )
        space = field_space.interior_space
        assert isinstance(space, RadialCharacteristicInteriorSpace)  # System-B mint; narrowing only
        #: The ψ½ interior split space (level metadata for repr; non-None by
        #: the ctor guard). The declared codomain is the composite space.
        self._ray_space = space

    # ── Predicates / spaces ───────────────────────────────────────────

    @property
    def is_adjointable(self) -> bool:
        # Both directions realized: :meth:`apply` (bulk → ray emission) and
        # :meth:`apply_transpose` (ray cotangent → bulk pullback, the
        # S-adjoint term). is_invertible inherits False: a rectangular
        # coupling (bulk → ray) is not square.
        return True

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        # System A — the FullField composite carrier (B.2b re-type; the
        # 2-blockification of System A lands with the B.2d ray eviction).
        return self._full_field_space

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        # System B — the ψ½ composite member space (B.2b DP1; non-None by
        # the ctor guard: A_BA exists only where System B does).
        return self._field_space

    # ── Forward — bulk flux → ψ½ ray emission ─────────────────────────

    def apply(self, psi: "FullField", /) -> "RadialCharacteristicField":
        r"""Emit the bulk within-group source onto the ψ½ ray.

        Integrate the bulk flux to :math:`\phi_0`, apply the emission kernel
        :math:`K` (the isotropic ℓ=0 source), and fold that moment at the
        closed :math:`\mu = \pm 1` rays. Returns **System B's own carrier**: a
        :class:`~orpheus.transport.radial_characteristic_field.RadialCharacteristicField`
        with SOURCE members — the folded emission in the interior cells, a
        zero boundary corner (the fold writes cells only). No present-zero
        bulk/boundary padding: the codomain has no such slots, so the old
        disjointness with ``S``'s bulk is structural (Pattern 4).

        Parameters
        ----------
        psi : FullField
            The within-group System-A iterate. Its ``interior`` must be a
            full-angular
            :class:`~orpheus.transport.fields.angular_flux.AngularFlux` (1-D
            curvilinear); the seed emission needs the per-ordinate flux to
            integrate.

        Returns
        -------
        RadialCharacteristicField
            The ψ½ ray source emission (interior = folded cells, boundary =
            zero corner; both SOURCE role).
        """
        from orpheus.transport.fields.angular_flux import AngularFlux
        from orpheus.transport.radial_characteristic_field import (
            RadialCharacteristicField,
        )

        bulk = psi.interior
        if not isinstance(bulk, AngularFlux):
            raise TypeError(
                "RadialCharacteristicEmission.apply: a seed-carrying composite "
                "must have a full-angular AngularFlux bulk (1-D curvilinear); "
                f"got {type(bulk).__name__}."
            )
        phi0 = bulk.integrate_angular().values              # (ng, nx)
        q0 = np.asarray(self.emission_kernel.apply(phi0))   # (ng, nx) iso ℓ=0
        return self._fold.apply(q0[None])                   # the q½ composite

    # ── Euclidean transpose — ray cotangent → bulk pullback ───────────

    def apply_transpose(
        self, cotangent: "RadialCharacteristicField", /,
    ) -> "FullField":
        r"""Euclidean transpose :math:`A_{BA}^{\mathsf T}` — ray → bulk pullback.

        The adjoint of :meth:`apply`: pull a System-B cotangent back into the
        bulk as
        :math:`\int\!d\mu^{\mathsf T}\,K^{\mathsf T}\,\mathrm{Fold}^{\mathsf T}`.
        Reconstructionᵀ lifts the ray cotangent to the :math:`\ell = 0` bulk
        moment, :math:`K^{\mathsf T}` scatters it in energy, and the
        ``integrate`` transpose broadcasts it per ordinate with the quadrature
        weight :math:`w_n` (:math:`(\int d\mu)^{\mathsf T}\,x = w_n\,x`). This is
        exactly the ``w · K_isoᵀ(Reconstructionᵀ χ_seed)`` bulk pullback the
        S-adjoint carries (single source).  The fold transpose reads the
        composite's interior member directly.

        Parameters
        ----------
        cotangent : RadialCharacteristicField
            A cotangent on System B (the ray source codomain). The old "the
            cotangent carries no ψ½ block" failure is now UNSPELLABLE — the
            composite always carries both members (Pattern 4).

        Returns
        -------
        FullField
            The System-A cotangent: ``interior`` = the bulk per-ordinate
            pullback over a zero trace (System A's honest 2-block row
            space, post-B.2d).
        """
        from orpheus.transport.full_field import FullField
        from orpheus.transport.radial_characteristic_field import (
            RadialCharacteristicField,
        )
        from orpheus.transport.source_sinks import (
            AngularBoundarySourceSink,
            AngularSourceSink,
        )

        if not isinstance(cotangent, RadialCharacteristicField):
            raise TypeError(
                "RadialCharacteristicEmission.apply_transpose: the cotangent "
                "must be a RadialCharacteristicField (System B's member "
                f"carrier); got {type(cotangent).__name__}."
            )
        m_bar = self._fold.apply_transpose(cotangent)        # (1, ng, nx)
        phi0_bar = np.asarray(
            self.emission_kernel.apply_transpose(m_bar[0]),  # (ng, nx)
        )
        # (∫dμ)ᵀ: broadcast the bulk moment cotangent per ordinate with the
        # quadrature weight w_n (the exact transpose of ``integrate_angular``'s
        # Σ_n w_n ψ_n) — the same w the S-adjoint pullback used.
        w = np.asarray(self._quadrature.weights, dtype=float)  # (N,)
        bulk_bar = w.reshape((w.size,) + (1,) * phi0_bar.ndim) * phi0_bar[None]
        return FullField(
            interior=AngularSourceSink(values=bulk_bar, space=self._angular_bulk_space),
            boundary=AngularBoundarySourceSink.zeros(self._angular_trace),
        )

    def __repr__(self) -> str:
        return (
            f"RadialCharacteristicEmission("
            f"emission_kernel={type(self.emission_kernel).__name__}, "
            f"levels={self._ray_space.levels}, ng={self._ray_space.ng}, "
            f"nx={self._ray_space.nx})"
        )


# ``A_BA`` sits block-native in the (B, A) slot of the within-group gain grid
# (:func:`orpheus.sn.coupled_system.build_within_group_system`).
