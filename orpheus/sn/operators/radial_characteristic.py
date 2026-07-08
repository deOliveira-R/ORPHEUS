r"""``A_BB`` — System B's radial straight-characteristic transport operator.

The self-block of the ψ½ ray in the 2×2 coupled block operator (campaign
:mod:`coupled block operator <orpheus.sn>` — the augmented within-group
system re-partitioned as two systems)::

    [ A_AA   A_AB ] [ transport ]      A_AA = L + C − S − B   (System A)
    [ A_BA   A_BB ] [ ray       ]      A_BB = this operator   (System B)

**What it is (operator algebra).** :class:`RadialCharacteristicOperator`
is the banded radial transport operator :math:`\mu\,\partial_r + \sigma_t`
acting on the starting-direction flux :math:`\psi_{1/2}` at the closed
:math:`\mu = \pm 1` rays of a curvilinear μ-level (Hébert §3.9.4). At those
rays the angular redistribution coefficient :math:`1-\mu^2` vanishes
(:math:`\alpha_{1/2}=0`, Hébert 3.423), so the streaming–collision balance
DECOUPLES from the α-cascade and reduces to a plain diamond-difference (DD)
recurrence in radius — a *straight characteristic*. Its domain and codomain
are the ψ½ carrier
:class:`~orpheus.numerics.spaces.radial_characteristic_space.RadialCharacteristicSpace`
(``sn_mesh.radial_characteristic_space``); it exists ONLY on a seed-carrying
mesh (the sphere, R12a — see that space's module docstring).

**How it is posed — a two-point radial BVP.** System B carries TWO boundary
conditions, so the ray solve is a well-posed two-point boundary-value
problem closed by its own boundary block ``B_b``
(:class:`~orpheus.sn.operators.boundary.RadialCharacteristicBoundaryOperator`):

* **r = R Dirichlet** — the outer-face inflow corner (:math:`\mu=-1`) is
  GIVEN data (0 for vacuum; the reflected outflow corner for reflective —
  ``B_b``'s specular swap feeds it into the source);
* **r = 0 pole continuation** — :math:`\psi_{1/2}^{+}(0)=\psi_{1/2}^{-}(0)`;
  the inward leg's exit face IS the outward leg's entry face (internal to
  the march, not a free datum).

**How it is solved — the direct march IS the resolvent.** Because the
recurrence is triangular in radius (:math:`\rho=0` — the #284
forward-substitution certificate), the two-leg Carlson march
(:func:`~orpheus.sn.spatial.psi_half_angle_seed.carlson_inward_sweep_from_source`,
Hébert 3.434–3.435) is the EXACT direct inverse :math:`A_{BB}^{-1}` — no
iteration, no previous-iterate seed. :meth:`solve` is that march; it is a
thin WRAP of the standalone engine (single source — Cardinal Rule 2 — the
SAME engine the in-sweep direct-seed block
:mod:`orpheus.sn.loss_representation` runs), and :meth:`solve_transpose` is
its Euclidean reverse-mode adjoint
(:func:`~orpheus.sn.spatial.psi_half_angle_seed.carlson_inward_sweep_transpose`).

Scope of this realization (campaign step 1c) — resolvent action only
====================================================================

This operator realizes the **resolvent action** :math:`A_{BB}^{-1}` as
:meth:`solve` (the march) and its adjoint as :meth:`solve_transpose`. The
**forward** action :math:`A_{BB}\,\psi = (\mu\,\partial_r + \sigma_t)\psi`
(:meth:`apply`) is DEFERRED to campaign step 4, and with it the
operator-returning :meth:`inverse` and ``is_invertible = True``. This is a
principled coherence deferral, not an omission:

* The forward matvec is currently woven into the composite ``(L+C).apply``
  (the ray residual on the :math:`\mu=+1` outflow corner, ruling R13).
  Spelling it standalone HERE would either duplicate that DD recurrence (a
  Pattern-2 twin path that could drift by a rounding) or force the
  hot-path extraction out of ``(L+C).apply`` prematurely — the campaign's
  highest-risk carve, gated at step 4 (the ``CoupledOperator`` assembly,
  which is the forward matvec's first real consumer).

**Tracked transient twin (Cardinal Rule 2 — retired at step 4/5).** The
:math:`\sigma_t`-driven DD *engine*
(:func:`~orpheus.sn.spatial.psi_half_angle_seed.carlson_inward_sweep_from_source`)
is single-sourced, but the two-leg ORCHESTRATION that :meth:`solve` runs
(read the source views → inward leg → pole-continue → reversed outward leg
→ write the flux views) currently mirrors the in-sweep production block
``orpheus/sn/loss_representation:4104-4119`` — the direct-seed march that
``(L+C).solve`` runs inline. This transient twin is the *reason* step 1c
poses the named operator FIRST: steps 4/5 route the production ``(L+C)``
ray solve THROUGH this operator (the coupled block-triangular resolvent),
retiring the inline block so the orchestration lives in ONE place. Until
then both sides are behaviour-pinned (the in-sweep by the regression floor
+ sweep suite; this by the Mode-11 WRAP bit-identity gate
``test_psi_half_coupling.py::TestA_BB_RadialBVP.test_wrap_executes_engine_bit_identical``),
and the campaign plan's step-4/5 retirement list carries the inline block.
* The inverse-family involution web is two-directional
  (``inverse().solve`` == the forward ``apply``); shipping :meth:`inverse`
  while ``apply`` is deferred would yield a resolvent whose ``solve``
  raises — a false capability. So the forward ``apply``, :meth:`inverse`,
  and ``is_invertible = True`` land TOGETHER in step 4, keeping the
  faithfulness contract (``is_invertible=True`` ⟺ a working ``solve``)
  honest here.

``is_invertible`` therefore reads ``False`` (the base default) in this
step — it advertises the ``inverse()`` OPERATOR-returning act, which does
not yet exist — while the mathematical invertibility is expressed by the
presence of a working :meth:`solve`. Likewise ``is_adjointable`` stays
``False`` (the FORWARD transpose :math:`A_{BB}^{\mathsf T}` defers with the
forward ``apply``; the *solve*-transpose :meth:`solve_transpose` is a
distinct realized verb, the adjoint of the inverse action).

Sourcing
========

``sn_mesh`` supplies the ray carrier and the radial widths
(:math:`\Delta r = ` ``sn_mesh.axis_widths[0]``); ``total_cross_section``
is the :math:`\sigma_t`
:class:`~orpheus.transport.fields.cross_section_field.CrossSectionField` — the
collision coefficient :math:`C_{\rm ray}` the march reads (``.values`` are
``(ng, nx)``). It is a TYPED, mesh-bound field (not a bare array) so the
mesh-identity invariant — σ_t and the operator march the SAME radial widths —
is enforced at construction (Pattern 4; the sibling
:class:`~orpheus.sn.operators.streaming.InvertibleOperator` guards its
diagonal coefficient the same way). Step 4 sources it from the solver's
macroscopic XS as ``mat_xs.total_cross_section_field``.

References
==========

* Hébert, A. (2009). *Applied Reactor Physics*. §3.9.4, Eqs. (3.432)–(3.435).
* GH #282 (the direct ψ½ seed / route (a)), #280 (the walk unification),
  #284 (the forward-substitution / triangular certificate).
* ``.claude/plans/coupled_block_operator_campaign.md`` — the 2×2 posing,
  rulings P1/P2, and this step's design reconnaissance.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

import numpy as np

from orpheus.numerics.operator import LinearOperator
from orpheus.sn.spatial.psi_half_angle_seed import (
    carlson_inward_sweep_from_source,
    carlson_inward_sweep_transpose,
)

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace
    from orpheus.numerics.spaces.radial_characteristic_space import (
        RadialCharacteristicSpace,
    )
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.transport.fields._bases import RadialCharacteristicField
    from orpheus.transport.fields.cross_section_field import CrossSectionField
    from orpheus.transport.fields.radial_characteristic_flux import (
        RadialCharacteristicFlux,
    )
    from orpheus.transport.source_sinks import RadialCharacteristicSourceSink


__all__ = ["RadialCharacteristicOperator"]


class RadialCharacteristicOperator(LinearOperator["RadialCharacteristicField"]):
    r"""``A_BB`` — the radial straight-characteristic transport operator on ψ½.

    System B's self-block of the 2×2 coupled block operator: the banded
    radial DD recurrence :math:`\mu\,\partial_r + \sigma_t` at the closed
    :math:`\mu = \pm 1` rays (Hébert §3.9.4). Endomorphic on the ψ½ carrier
    :class:`~orpheus.numerics.spaces.radial_characteristic_space.RadialCharacteristicSpace`
    (``sn_mesh.radial_characteristic_space``).

    See the module docstring for the operator-algebra posing (two-point
    radial BVP), the solve (the direct Carlson march IS the resolvent
    :math:`A_{BB}^{-1}`), and the step-1c realization scope (resolvent
    action :meth:`solve` + adjoint :meth:`solve_transpose`; the forward
    :meth:`apply` / :meth:`inverse` / ``is_invertible`` land in step 4).

    Parameters
    ----------
    sn_mesh : SNMesh
        The augmented geometry — seed-carrying (1-D curvilinear, R12a).
        Supplies the ray carrier ``radial_characteristic_space`` and the
        radial widths ``axis_widths[0]``. A seedless mesh (Cartesian /
        cylinder) has NO System B: constructing over one is rejected.
    total_cross_section : CrossSectionField
        The total cross-section :math:`\sigma_t` as a typed, mesh-bound
        :class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`
        on ``sn_mesh`` (``.values`` are ``(ng, nx)``) — the collision
        coefficient :math:`C_{\rm ray}` the march reads. Must be on the SAME
        mesh as the operator (the mesh-identity invariant — refused otherwise)
        and strictly positive everywhere (the DD denominator
        :math:`\Delta r\,\sigma + 2` is then well-defined and the march never
        emits NaN).
    """

    def __init__(
        self, sn_mesh: "SNMesh", total_cross_section: "CrossSectionField",
    ) -> None:
        space = sn_mesh.radial_characteristic_space
        if space is None:
            raise ValueError(
                "RadialCharacteristicOperator: the mesh carries no "
                "starting-direction ray (radial_characteristic_space is "
                "None) — a seedless mesh (Cartesian / cylinder, R12a) has "
                "no System B. A_BB exists only on a seed-carrying mesh "
                "(the sphere)."
            )
        # Mesh-identity invariant (Pattern 4 — make the illegal state
        # unrepresentable): the march reads THIS mesh's radial widths
        # (axis_widths[0]) against σ_t, so a σ_t from a foreign mesh — even one
        # with an identical (ng, nx) — would march the wrong Δr. The typed
        # CrossSectionField carries its mesh, so we can refuse the mismatch at
        # construction (the sibling InvertibleOperator guards σ↔mesh the same
        # way, streaming.py). The field's own space invariant then GUARANTEES
        # values.shape == (ng, nx) for this mesh — so no separate shape check is
        # needed (an explicit one would be redundant ceremony).
        if total_cross_section.mesh is not sn_mesh:
            raise ValueError(
                "RadialCharacteristicOperator: total_cross_section must be a "
                "CrossSectionField on the SAME SNMesh as the operator "
                "(mesh-identity invariant); got field mesh "
                f"{total_cross_section.mesh!r} vs operator mesh {sn_mesh!r}."
            )
        sigma = total_cross_section.values
        if np.any(sigma <= 0.0):
            raise ValueError(
                f"RadialCharacteristicOperator: total_cross_section must be "
                f"strictly positive everywhere for the DD march denominator "
                f"(Δr·σ + 2) to be well-defined; got "
                f"min(σ_t) = {float(np.min(sigma)):.3e}."
            )
        #: The augmented geometry (ray carrier + radial widths).
        self.sn_mesh = sn_mesh
        #: The total cross-section :math:`\sigma_t` — the ``C_ray`` collision
        #: coefficient, a typed :class:`CrossSectionField` on :attr:`sn_mesh`
        #: (``.values`` are ``(ng, nx)``; the mesh-identity guard above makes a
        #: foreign-mesh σ_t unconstructable). The march reads ``.values``.
        self.total_cross_section = total_cross_section
        #: The ψ½ carrier — domain and codomain (non-None by the ctor guard).
        self._ray_space: "RadialCharacteristicSpace" = space

    # ── Spaces ────────────────────────────────────────────────────────
    #
    # Endomorphic on the ψ½ carrier: both the source (input to solve) and
    # the flux (its output) live on the SAME RadialCharacteristicSpace —
    # the role (source / flux) lives on the field types, not the space.
    # block_role stays None (the base default): A_BB is System B's self-
    # block on the ray space, not a bulk/full/boundary block of the
    # FullField composite; the SystemRole {A, B, COUPLED} lattice widening
    # is campaign step 4.

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        return self._ray_space

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        return self._ray_space

    # ── Forward action — DEFERRED to campaign step 4 ──────────────────

    def apply(self, x: "RadialCharacteristicField", /) -> "RadialCharacteristicField":
        r"""The forward matvec :math:`A_{BB}\,\psi = (\mu\,\partial_r + \sigma_t)\psi`.

        DEFERRED to campaign step 4 (loud, never a silent wrong answer).
        The forward ray action is currently woven into the composite
        ``(L+C).apply`` (the ray residual on the :math:`\mu=+1` outflow
        corner, ruling R13); extracting it standalone here would duplicate
        that DD recurrence (a twin path) or force the hot-path carve
        prematurely. It lands in step 4 with the ``CoupledOperator``
        assembly — the forward matvec's first real consumer — as the single
        source (see the module docstring "Scope of this realization").
        """
        raise NotImplementedError(
            "RadialCharacteristicOperator.apply (the forward matvec "
            "A_BB ψ = (μ ∂_r + σ_t) ψ) is deferred to campaign step 4: it "
            "is currently woven into (L+C).apply (ruling R13), and is "
            "extracted there — with inverse() and is_invertible=True — as "
            "the CoupledOperator's (B,B) block. Step 1c realizes the "
            "resolvent action via .solve (the direct Carlson march)."
        )

    # ── Resolvent action — the direct Carlson march IS A_BB⁻¹ ─────────

    def solve(
        self, source: "RadialCharacteristicSourceSink",
    ) -> "RadialCharacteristicFlux":
        r"""Solve :math:`A_{BB}\,\psi_{1/2} = q_{1/2}` by the two-leg Carlson march.

        The EXACT direct inverse :math:`A_{BB}^{-1}` (no iteration): per
        seed-carrying level, the inward :math:`\mu=-1` leg marches from the
        given r = R inflow corner to the pole, then the outward
        :math:`\mu=+1` leg rides the SAME engine on reversed cell data
        (orientation is carried by the DATA, never a flag) from the
        pole-continued face out to the r = R outflow corner. A thin WRAP of
        :func:`~orpheus.sn.spatial.psi_half_angle_seed.carlson_inward_sweep_from_source`
        (Hébert 3.434–3.435) — the SAME single-sourced DD engine the
        in-sweep direct-seed block runs, so this is bit-identical to that
        block (pinned by the Mode-11 WRAP gate). The two-leg orchestration
        HERE is a tracked transient twin of that block, retired at step 4/5
        — see the module docstring "Tracked transient twin".

        The per-level slot key is the carrier's own ``space.levels`` member
        (the level POSITION, ``p_idx`` in the in-sweep) — the coordinate
        that keys the space slots, pinned by
        :mod:`tests.sn.mesh.test_radial_characteristic_slot_coordination`.

        Parameters
        ----------
        source : RadialCharacteristicSourceSink
            The q½ source block — its cells legs carry the folded
            starting-direction source and its :math:`\mu=-1` corner the
            r = R inflow (Dirichlet) datum. Must share this operator's mesh.

        Returns
        -------
        RadialCharacteristicFlux
            The ψ½ state :math:`\psi_{1/2}` satisfying the two-point radial
            BVP — cells + both corner legs filled per carried level.
        """
        from orpheus.transport.fields.radial_characteristic_flux import (
            RadialCharacteristicFlux,
        )

        mesh = self.sn_mesh
        space = self._require_same_carrier(source, "solve")
        sigma = self.total_cross_section.values
        dr = mesh.axis_widths[0]
        src_vals = source.values

        flux = RadialCharacteristicFlux.zeros_on(mesh)
        buf = flux.values
        for level in space.levels:
            q_minus = space.cells_view(src_vals, level, -1)
            q_plus = space.cells_view(src_vals, level, +1)
            corner_in = space.corner_view(src_vals, level, -1)
            # inward μ=−1 leg: enter at the r=R inflow corner, exit at the pole.
            cells_minus, pole_face = carlson_inward_sweep_from_source(
                q_minus, sigma, dr, corner_in,
            )
            # outward μ=+1 leg: the SAME engine on reversed data, entering at
            # the pole-continued face (ψ½⁺(0) = ψ½⁻(0)) and exiting at r=R.
            cells_plus_rev, corner_out = carlson_inward_sweep_from_source(
                q_plus[:, ::-1], sigma[:, ::-1], dr[::-1], pole_face,
            )
            space.cells_view(buf, level, -1)[...] = cells_minus
            space.corner_view(buf, level, -1)[...] = corner_in
            space.cells_view(buf, level, +1)[...] = cells_plus_rev[:, ::-1]
            space.corner_view(buf, level, +1)[...] = corner_out
        return flux

    def solve_transpose(
        self, cotangent: "RadialCharacteristicField",
    ) -> "RadialCharacteristicSourceSink":
        r"""The Euclidean adjoint of :meth:`solve` — :math:`(A_{BB}^{-1})^{\mathsf T}`.

        The reverse-mode adjoint of the two-leg march: given a cotangent on
        the flux (the solve's codomain), return the cotangent on the source
        (its domain). Per level, the OUTWARD leg is reversed first (its exit
        corner feeds the pole-face cotangent), then the INWARD leg (the pole
        cotangent is its exit), threading the running face cotangent back to
        the r = R inflow corner — the transpose of :meth:`solve`'s leg chain,
        via
        :func:`~orpheus.sn.spatial.psi_half_angle_seed.carlson_inward_sweep_transpose`.

        This is the ISOLATED ray-block transpose :math:`(A_{BB}^{-1})^{\mathsf T}`
        — the pure resolvent adjoint. The full ``(L+C)`` reverse-scan adds a
        seed→bulk term (the M-M thread cotangent) on the inward cells; that
        term is the ``A_AB`` coupling's transpose (campaign steps 2–3), NOT
        part of ``A_BB`` in isolation.

        Parameters
        ----------
        cotangent : RadialCharacteristicField
            A cotangent on the flux space (role-erased: a
            :class:`~orpheus.transport.fields.radial_characteristic_flux.RadialCharacteristicFlux`
            cotangent). Must share this operator's mesh.

        Returns
        -------
        RadialCharacteristicSourceSink
            The cotangent on the q½ source. The :math:`\mu=+1` source corner
            is unused by the march (the q½ fold writes only cells + the
            :math:`\mu=-1` corner), so it stays zero.
        """
        from orpheus.transport.source_sinks import RadialCharacteristicSourceSink

        mesh = self.sn_mesh
        space = self._require_same_carrier(cotangent, "solve_transpose")
        sigma = self.total_cross_section.values
        dr = mesh.axis_widths[0]
        cot_vals = cotangent.values

        src_bar = np.zeros_like(cot_vals)
        for level in space.levels:
            cells_minus_bar = space.cells_view(cot_vals, level, -1)
            cells_plus_bar = space.cells_view(cot_vals, level, +1)
            corner_in_bar = space.corner_view(cot_vals, level, -1).copy()
            corner_out_bar = space.corner_view(cot_vals, level, +1)
            # reverse the OUTWARD (+1) leg — marched on reversed data; its exit
            # corner cotangent seeds the pole-face cotangent for the inward leg.
            q_plus_rev_bar, pole_face_bar = carlson_inward_sweep_transpose(
                cells_plus_bar[:, ::-1], corner_out_bar, sigma[:, ::-1], dr[::-1],
            )
            q_plus_bar = q_plus_rev_bar[:, ::-1]
            # reverse the INWARD (−1) leg — the pole face is its exit.
            q_minus_bar, corner_in_from_minus = carlson_inward_sweep_transpose(
                cells_minus_bar, pole_face_bar, sigma, dr,
            )
            # the r=R inflow corner both passes through to the flux corner AND
            # enters the inward leg — its cotangent is the sum of both paths.
            corner_in_bar = corner_in_bar + corner_in_from_minus
            space.cells_view(src_bar, level, -1)[...] = q_minus_bar
            space.cells_view(src_bar, level, +1)[...] = q_plus_bar
            space.corner_view(src_bar, level, -1)[...] = corner_in_bar
            # source corner(+1) unused by the q½ fold → stays zero.
        return RadialCharacteristicSourceSink(
            values=src_bar, space=space, mesh=mesh,
        )

    # ── Internals ─────────────────────────────────────────────────────

    def _require_same_carrier(
        self, field: "RadialCharacteristicField", method: str,
    ) -> "RadialCharacteristicSpace":
        r"""Enforce the mesh-identity invariant and return the shared carrier.

        Single source for the ``solve`` / ``solve_transpose`` guard: the
        input field and this operator must share the SAME
        :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` instance (so the
        field's carrier, this operator's :math:`\sigma_t`, and
        ``axis_widths`` cannot desync). Returns the field's own carrier —
        the same object as :attr:`_ray_space` under the invariant.
        """
        if field.mesh is not self.sn_mesh:
            raise ValueError(
                f"RadialCharacteristicOperator.{method}: the input field and "
                f"the operator must share the same SNMesh instance "
                f"(mesh-identity invariant); got field mesh {field.mesh!r} "
                f"vs operator mesh {self.sn_mesh!r}."
            )
        return field.space

    def __repr__(self) -> str:
        return (
            f"RadialCharacteristicOperator(levels={self._ray_space.levels}, "
            f"ng={self._ray_space.ng}, nx={self._ray_space.nx})"
        )
