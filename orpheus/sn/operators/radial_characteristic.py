r"""The ψ½ System-B block operators — ``A_BB`` (self-block) + ``A_AB`` (ray→bulk).

The two System-B blocks of the 2×2 coupled block operator (campaign
:mod:`coupled block operator <orpheus.sn>` — the augmented within-group
system re-partitioned as two systems)::

    [ A_AA   A_AB ] [ transport ]   A_AB = RadialCharacteristicSeeding (ray→bulk seed)
    [ A_BA   A_BB ] [ ray       ]   A_BB = RadialCharacteristicOperator (System B self-block)

This module hosts BOTH: :class:`RadialCharacteristicOperator` (``A_BB``, the
radial straight-characteristic self-block, documented in depth below) and
:class:`RadialCharacteristicSeeding` (``A_AB``, the cell-local angular ray→bulk
seed injection — see its class docstring for the operator-algebra posing). The
remainder of this module docstring documents ``A_BB``.

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

Scope of this realization (campaign step 4b) — the full operator
================================================================

This operator realizes BOTH directions of the complete ``A_BB`` machinery:

* the **forward** action :math:`A_{BB}\,\psi_{1/2} = (\mu\,\partial_r +
  \sigma_t)\,\psi_{1/2}` as :meth:`apply`, and its Euclidean transpose
  :math:`A_{BB}^{\mathsf T}` as :meth:`apply_transpose`;
* the **resolvent** action :math:`A_{BB}^{-1}` as :meth:`solve` (the direct
  Carlson march) and its adjoint :math:`(A_{BB}^{-1})^{\mathsf T}` as
  :meth:`solve_transpose`;
* the operator-returning :meth:`inverse` (a generic
  :class:`~orpheus.numerics.operator.InverseOperator` whose ``apply`` IS
  :meth:`solve`), so ``is_invertible`` and ``is_adjointable`` both read
  ``True`` — the involution web ``inverse().solve == apply`` closes.

**Single-sourced forward (campaign step 4b — extract, not twin).** The
forward matvec is the exact algebraic inverse of the march; it lives in ONE
place —
:func:`~orpheus.sn.spatial.psi_half_angle_seed.radial_characteristic_forward_residual`
— which BOTH :meth:`apply` here AND the fused ``(L+C)`` walk's ψ½ rows
(:meth:`~orpheus.sn.loss_representation._OneDimScanWalk._seed_rows_forward`)
call. There is NO forward twin: step 4b reduced the walk method to a thin
wrapper over the shared kernel (the user ruled "extract the shared kernel
now" over a tracked twin, closing the "could drift by a rounding" risk by
construction — Cardinal Rule 2). The transpose is single-sourced the same
way
(:func:`~orpheus.sn.spatial.psi_half_angle_seed.radial_characteristic_forward_residual_transpose`,
the PURE ``A_BB`` transpose; the walk's ``_seed_rows_transpose`` wraps it and
ADDS the ``A_AB`` seed→bulk coupling, which is not part of ``A_BB`` in
isolation).

**Tracked transient twin — the SOLVE orchestration (Cardinal Rule 2, retired
at step 4e).** The :math:`\sigma_t`-driven DD *engine*
(:func:`~orpheus.sn.spatial.psi_half_angle_seed.carlson_inward_sweep_from_source`)
is single-sourced, but the two-leg ORCHESTRATION that :meth:`solve` runs
(read the source views → inward leg → pole-continue → reversed outward leg
→ write the flux views) still mirrors the in-sweep production block
``orpheus/sn/loss_representation:4104-4131`` — the direct-seed march that
``(L+C).solve`` runs inline. Step 4e routes the production ``(L+C)`` ray
solve THROUGH this operator (the coupled block-triangular resolvent),
retiring the inline block so the orchestration lives in ONE place. Until
then both sides are behaviour-pinned (the in-sweep by the regression floor +
sweep suite; this by the Mode-11 WRAP bit-identity gate) and the campaign
plan's step-4e retirement list carries the inline block.

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

from orpheus.numerics.operator import (
    InverseOperator,
    LinearOperator,
    NotInvertible,
    SystemRole,
)
from orpheus.sn.spatial.psi_half_angle_seed import (
    carlson_inward_sweep_from_source,
    carlson_inward_sweep_transpose,
    radial_characteristic_forward_residual,
    radial_characteristic_forward_residual_transpose,
)

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace
    from orpheus.numerics.spaces.radial_characteristic_space import (
        RadialCharacteristicSpace,
    )
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.transport.fields._bases import (
        AngularField,
        RadialCharacteristicField,
    )
    from orpheus.transport.fields.cross_section_field import CrossSectionField
    from orpheus.transport.fields.radial_characteristic_flux import (
        RadialCharacteristicFlux,
    )
    from orpheus.transport.source_sinks import RadialCharacteristicSourceSink
    from orpheus.transport.source_sinks.angular_source_sink import (
        AngularSourceSink,
    )


__all__ = ["RadialCharacteristicOperator", "RadialCharacteristicSeeding"]


class RadialCharacteristicOperator(LinearOperator["RadialCharacteristicField"]):
    r"""``A_BB`` — the radial straight-characteristic transport operator on ψ½.

    System B's self-block of the 2×2 coupled block operator: the banded
    radial DD recurrence :math:`\mu\,\partial_r + \sigma_t` at the closed
    :math:`\mu = \pm 1` rays (Hébert §3.9.4). Endomorphic on the ψ½ carrier
    :class:`~orpheus.numerics.spaces.radial_characteristic_space.RadialCharacteristicSpace`
    (``sn_mesh.radial_characteristic_space``).

    See the module docstring's "Scope of this realization" section for the
    operator-algebra posing (two-point radial BVP) and the full realized
    surface — the forward :meth:`apply` / :meth:`apply_transpose`, the
    resolvent :meth:`solve` / :meth:`solve_transpose` (the direct Carlson
    march IS :math:`A_{BB}^{-1}`), and the operator-returning :meth:`inverse`
    (so ``is_invertible`` and ``is_adjointable`` both read ``True``; step 4b
    completed the forward, single-sourced with the ``(L+C)`` walk).

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
    # block_role stays None (the base default): A_BB is System B's self-block
    # on the ray space, not a bulk/full/boundary block of the FullField
    # composite. Its two-system membership rides the SystemRole axis (campaign
    # step 4a) — A_BB acts within System B alone.
    system_role = SystemRole.B

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        return self._ray_space

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        return self._ray_space

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
    ) -> "RadialCharacteristicSourceSink":
        r"""The forward matvec :math:`A_{BB}\,\psi_{1/2} = (\mu\,\partial_r
        + \sigma_t)\,\psi_{1/2}` — the exact algebraic inverse of :meth:`solve`.

        A thin WRAP of the single-sourced
        :func:`~orpheus.sn.spatial.psi_half_angle_seed.radial_characteristic_forward_residual`
        — the SAME body the fused ``(L+C)`` walk runs for the ψ½ rows
        (``_OneDimScanWalk._seed_rows_forward``), so ``A_BB.apply`` is
        bit-identical to that production forward (Cardinal Rule 2: one
        forward, no twin). Reads the flux state ψ½ (cells + corners); returns
        the residual source :math:`q_{1/2} = A_{BB}\,\psi_{1/2}`. The
        ``apply ∘ solve`` outflow-corner defect closes to ``0.0`` exactly; the
        cell round-trip is principled-equiv at ~FP ULP — the forward's
        :math:`2/\Delta r` and the march's :math:`\Delta r\,\sigma + 2`
        reassociate.

        Parameters
        ----------
        x : RadialCharacteristicField
            The ψ½ flux state (role-erased). Must share this operator's mesh.

        Returns
        -------
        RadialCharacteristicSourceSink
            The residual source :math:`A_{BB}\,\psi_{1/2}`.
        """
        from orpheus.transport.source_sinks import RadialCharacteristicSourceSink

        space = self._require_same_carrier(x, "apply")
        out = radial_characteristic_forward_residual(
            x.values, space, self.total_cross_section.values,
            self.sn_mesh.axis_widths[0],
        )
        return RadialCharacteristicSourceSink(values=out, space=space, mesh=self.sn_mesh)

    def apply_transpose(
        self, y: "RadialCharacteristicField", /,
    ) -> "RadialCharacteristicSourceSink":
        r"""The Euclidean transpose of the forward matvec —
        :math:`A_{BB}^{\mathsf T}`.

        A thin WRAP of
        :func:`~orpheus.sn.spatial.psi_half_angle_seed.radial_characteristic_forward_residual_transpose`
        — the PURE ``A_BB`` transpose. The walk's ``_seed_rows_transpose``
        wraps the SAME body then ADDS the ``A_AB`` seed→bulk recurrence
        coupling, which is NOT part of ``A_BB`` in isolation. This is the flat
        EUCLIDEAN adjoint :math:`A_{BB}^{\mathsf T}`; the metric Hilbert
        adjoint :math:`G^{-1}A_{BB}^{\mathsf T}G` (``.H``) is realized once at
        the composite (L19).
        """
        from orpheus.transport.source_sinks import RadialCharacteristicSourceSink

        space = self._require_same_carrier(y, "apply_transpose")
        out = radial_characteristic_forward_residual_transpose(
            y.values, space, self.total_cross_section.values,
            self.sn_mesh.axis_widths[0],
        )
        return RadialCharacteristicSourceSink(values=out, space=space, mesh=self.sn_mesh)

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


class RadialCharacteristicSeeding(
    LinearOperator["RadialCharacteristicField", "AngularField"],
):
    r"""``A_AB`` — the ray→bulk seed injection (the Morel–Montry angular seed).

    The off-diagonal ``(transport, ray)`` block of the 2×2 coupled block
    operator: the ψ½ starting-direction ray seeds the bulk angular recurrence.
    Domain = the ψ½ carrier
    :class:`~orpheus.numerics.spaces.radial_characteristic_space.RadialCharacteristicSpace`
    (``sn_mesh.radial_characteristic_space`` — the operator reads the inward
    :math:`\mu=-1` ``cells(p, -1)`` leg of a seed); codomain = the bulk
    per-ordinate residual
    :class:`~orpheus.transport.source_sinks.AngularSourceSink`. It exists ONLY
    on a seed-carrying mesh (the sphere, R12a).

    **What it is (operator algebra) — a CELL-LOCAL ANGULAR coupling.** At each
    radial cell :math:`i`, the ray value :math:`\psi_{1/2}(i)` is the seed of
    the Morel–Montry angular recurrence (Hébert §3.9.4)

    .. math::

        \psi_{m+1/2,\,i}
          = \frac{\psi_{m,\,i} - (1-\tau_m)\,\psi_{m-1/2,\,i}}{\tau_m},
        \qquad \psi_{-1/2,\,i} \equiv \psi_{1/2}(i),

    run over ORDINATES :math:`m` at a FIXED cell :math:`i`. The upstream
    half-flux :math:`\psi_{m-1/2,\,i}` then enters that cell's balance as the
    angular numerator :math:`(\Delta A/w)\,c_{\rm in}\,\psi_{m-1/2,\,i}`. So the
    seed at cell :math:`i` couples ONLY to the bulk ordinates at the SAME cell
    :math:`i` — there is **no spatial cell-cell coupling** (the radial march of
    :math:`\psi_{1/2}` is ``A_BB``'s job; the bulk's spatial DD face march is
    seed-independent). This is exactly what separates ``A_AB`` from ``A_BB``:
    ``A_BB``'s forward matvec is spatially woven into ``(L+C).apply`` and
    DEFERS to step 4, whereas ``A_AB`` is cell-local angular and realizes BOTH
    directions HERE as thin WRAPs of the single-sourced closure methods.

    **How it is realized — WRAPs of the shared M-M closure.** The seed→bulk map
    is the SAME
    :class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`
    machinery the ``(L+C)`` matvec runs (single source — Cardinal Rule 2);
    ``A_AB`` isolates its own block by ZEROING the bulk (forward) and
    DISCARDING the bulk cotangent (transpose):

    * :meth:`apply` (ray → bulk) — ``precompute_psi_state`` with an all-zero
      ``psi_view`` builds the seed-only half-angle grid. Because the recurrence
      is linear in :math:`(\psi_{\rm bulk}, \psi_{1/2})` jointly, the zero bulk
      isolates the ``A_AB`` part exactly (``A_AA``'s angular redistribution
      contributes nothing). Then per (carrying level, cell) ``cell_contribution``
      gives the angular numerator, placed as :math:`-(\Delta A/w)\,c_{\rm
      in}\,\psi_{m-1/2}/V` — the seed's contribution to the bulk residual
      :math:`m = (\mathrm{denom}\cdot\psi - \mathrm{numer})/V` with
      :math:`\psi_{\rm cell}=0` (the σ-diagonal and streaming terms drop out).
    * :meth:`apply_transpose` (bulk cotangent → ray seed cotangent) — the local
      gather :math:`\bar n_{p,\,m,\,i} = -\bar o_{m,\,i}/V_i` (the exact
      transpose of the forward :math:`-\cdot/V` placement), then
      ``angular_adjoint`` reverses the M-M recurrence to the seed cotangent
      ``seed_cells_bar[p]``, written on the ray ``cells(p, -1)`` leg. This is
      the :math:`\bar{c}_{-}\mathrel{+}=\mathrm{seed\_cells\_bar}[p]` term the
      in-sweep reverse (:mod:`orpheus.sn.loss_representation`) adds — the
      Euclidean transpose :math:`A_{AB}^{\mathsf T}`.

    Because both directions are realized, ``is_adjointable = True``.
    ``is_invertible`` inherits ``False`` — a rectangular coupling (ray → bulk)
    has no inverse.

    **The shared kernel (a step-4 note).** ``precompute_psi_state`` /
    ``cell_contribution`` / ``angular_adjoint`` each serve BOTH ``A_AA`` (the
    bulk angular redistribution, ``psi_view ≠ 0``) and ``A_AB`` (the seed).
    ``A_AB`` projects out its block by zeroing / discarding; the step-4
    ``CoupledOperator`` calls the ONE angular kernel and routes it into the
    ``A_AA`` and ``A_AB`` blocks — never a twin kernel.

    **Tracked transient twin (Cardinal Rule 2 — retired at step 4/5).** The M-M
    *kernel* is single-sourced, but the thin :math:`\mp\,\mathrm{numer}/V`
    orchestration that :meth:`apply` / :meth:`apply_transpose` run around it
    mirrors the in-sweep placement fused into the ``(L+C)`` matvec
    (:mod:`orpheus.sn.loss_representation` — the seed's ``angular_numer`` term in
    ``m = (denom·ψ − numer)/V`` on the forward, and the ``numer_bar →
    angular_adjoint → seed_cells_bar`` term on the reverse). This is the SAME
    kind of transient twin ``A_BB.solve`` carries, at a DIFFERENT production
    entry point (the ``apply`` walk, not the ``solve`` march). Steps 4/5 route
    the production ``(L+C)`` bulk rows through ``A_AA + A_AB`` (the
    ``CoupledOperator`` block matvec), retiring the inline placement so the
    orchestration lives in ONE place (campaign retirement-list entries 3–4).
    Until then both sides are behaviour-pinned — the in-sweep by the regression
    floor + sweep suite, this by the bit-identity gates ``TestA_AB_SeedInjection``.

    Parameters
    ----------
    sn_mesh : SNMesh
        The augmented geometry — seed-carrying (1-D curvilinear, R12a). Supplies
        the ray carrier ``radial_characteristic_space`` (the domain), the M-M
        closure ``pole_angular_closure`` (the single-sourced kernel), the cell
        volumes ``volumes``, and the quadrature ``quad``. A seedless mesh
        (Cartesian / cylinder) has NO ray→bulk coupling: constructing over one
        is rejected. Unlike ``A_BB``, ``A_AB`` needs NO :math:`\sigma_t` — with
        the bulk zeroed the collision/streaming terms drop out and only the
        σ-independent angular numerator survives (so the coupling is a pure
        function of the mesh geometry + quadrature).
    """

    # A_AB maps System B (the ray seed) → System A (the bulk residual): an
    # off-diagonal block, so it spans both systems (campaign step 4a).
    system_role = SystemRole.COUPLED

    def __init__(self, sn_mesh: "SNMesh") -> None:
        space = sn_mesh.radial_characteristic_space
        if space is None:
            raise ValueError(
                "RadialCharacteristicSeeding: the mesh carries no "
                "starting-direction ray (radial_characteristic_space is None) "
                "— a seedless mesh (Cartesian / cylinder, R12a) has no System "
                "B, hence no ray→bulk seed coupling to inject. A_AB exists "
                "only on a seed-carrying mesh (the sphere)."
            )
        from orpheus.transport.source_sinks.angular_source_sink import (
            AngularSourceSink,
        )
        #: The augmented geometry (ray carrier + the M-M closure + volumes).
        self.sn_mesh = sn_mesh
        #: The ψ½ carrier — the DOMAIN. The operator reads the inward
        #: ``cells(p, -1)`` leg of a seed on this space.
        self._ray_space: "RadialCharacteristicSpace" = space
        #: The bulk per-ordinate residual space — the CODOMAIN (the seed's
        #: contribution to the ``(L+C)`` bulk rows). Derived from the single
        #: source of the ``AngularSourceSink`` space identity.
        self._bulk_space: "FunctionSpace" = AngularSourceSink._space_for_mesh(
            sn_mesh,
        )

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
        # The ψ½ ray seed (the input to :meth:`apply`).
        return self._ray_space

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        # The bulk per-ordinate residual contribution (:meth:`apply` output).
        return self._bulk_space

    # ── Forward — seed → bulk angular-numerator contribution ──────────

    def apply(self, seed: "RadialCharacteristicField", /) -> "AngularSourceSink":
        r"""Inject the ψ½ ray seed into the bulk angular recurrence.

        The seed's contribution to the ``(L+C)`` bulk residual: build the
        seed-only Morel–Montry half-angle grid (``precompute_psi_state`` with
        an all-zero ``psi_view`` — the zero bulk isolates ``A_AB`` from
        ``A_AA``'s redistribution by linearity), then per carrying level and
        cell take the angular numerator :math:`(\Delta A/w)\,c_{\rm in}\,
        \psi_{m-1/2}` (``cell_contribution``) and place :math:`-\,\cdot\,/V` —
        the seed's term in :math:`m = (\mathrm{denom}\cdot\psi -
        \mathrm{numer})/V` with :math:`\psi_{\rm cell}=0`. Non-carrying-level
        ordinates stay zero (they have no ray seed). Bit-identical to the
        in-sweep injection (same single-sourced closure methods).

        Parameters
        ----------
        seed : RadialCharacteristicField
            The ψ½ ray flux. Only its inward ``cells(p, -1)`` leg is read
            (the recurrence seed); the ``+1`` leg and corners are ``A_BB``
            -internal. Must share this operator's mesh.

        Returns
        -------
        AngularSourceSink
            The seed's per-ordinate contribution to the bulk residual,
            shape ``(N, ng, nx)``.
        """
        from orpheus.transport.source_sinks.angular_source_sink import (
            AngularSourceSink,
        )

        mesh = self.sn_mesh
        self._check_mesh(seed, "apply")
        closure = mesh.pole_angular_closure
        space = self._ray_space
        N = mesh.quad.N
        ng = space.ng
        nx = mesh.nx
        V = mesh.volumes
        level_indices = closure.level_indices

        # Bulk zeroed → the recurrence propagates ONLY the seed (linearity).
        psi_state = closure.precompute_psi_state(
            np.zeros((N, ng, nx)), radial_characteristic=seed,
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
        return AngularSourceSink.from_mesh(out_g_first.swapaxes(0, 1), mesh)

    # ── Euclidean transpose — bulk cotangent → ray seed cotangent ─────

    def apply_transpose(
        self, cotangent: "AngularField", /,
    ) -> "RadialCharacteristicSourceSink":
        r"""Euclidean transpose :math:`A_{AB}^{\mathsf T}` — bulk → ray seed.

        The adjoint of :meth:`apply`: given a cotangent on the bulk residual
        (the codomain), return the cotangent on the ray seed (the domain).
        Reverse the forward's :math:`-\,\mathrm{numer}/V` placement with the
        local gather :math:`\bar n_{p,\,m,\,i} = -\bar o_{m,\,i}/V_i`, then
        ``angular_adjoint`` reverses the M-M recurrence to the per-carrying-
        level seed cotangent ``seed_cells_bar[p]``, written on the inward
        ``cells(p, -1)`` leg. The bulk-redistribution cotangent
        (``psi_ang_bar``, ``A_AA``'s share of the shared kernel) is discarded.
        The ``+1`` leg and corners stay zero (the forward writes only the
        inward leg).

        Parameters
        ----------
        cotangent : AngularField
            A cotangent on the bulk residual space (role-erased — any
            ``(N, ng, nx)`` bulk field). Must share this operator's mesh.

        Returns
        -------
        RadialCharacteristicSourceSink
            The ray seed cotangent — ``cells(p, -1)`` filled per carrying
            level, everything else zero.
        """
        from orpheus.transport.source_sinks.radial_characteristic_source_sink import (
            RadialCharacteristicSourceSink,
        )

        mesh = self.sn_mesh
        self._check_mesh(cotangent, "apply_transpose")
        closure = mesh.pole_angular_closure
        space = self._ray_space
        V = mesh.volumes
        level_indices = closure.level_indices

        ob_g_first = cotangent.values.swapaxes(0, 1)          # (ng, N, nx)
        # numer_bar for EVERY level (angular_adjoint needs the full tuple);
        # the local −ō/V gather is the exact transpose of the forward −·/V.
        numer_bar = tuple(
            -ob_g_first[:, np.asarray(level_idx), :] / V[None, None, :]
            for level_idx in level_indices
        )
        _, seed_cells_bar = closure.angular_adjoint(numer_bar)

        src_bar = RadialCharacteristicSourceSink.zeros_on(mesh)
        for p, cells_bar in seed_cells_bar.items():
            space.cells_view(src_bar.values, p, -1)[...] = cells_bar
        return src_bar

    # ── Internals ─────────────────────────────────────────────────────

    def _check_mesh(self, field: "RadialCharacteristicField | AngularField", method: str) -> None:
        r"""Enforce the mesh-identity invariant (Pattern 4).

        The input field, this operator's ``pole_angular_closure``, and the
        ``volumes`` must all be on the SAME
        :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` instance, so the seed
        legs, the M-M coefficients, and the ``/V`` scaling cannot desync.
        """
        if field.mesh is not self.sn_mesh:
            raise ValueError(
                f"RadialCharacteristicSeeding.{method}: the input field and "
                f"the operator must share the same SNMesh instance "
                f"(mesh-identity invariant); got field mesh {field.mesh!r} "
                f"vs operator mesh {self.sn_mesh!r}."
            )

    def __repr__(self) -> str:
        return (
            f"RadialCharacteristicSeeding(levels={self._ray_space.levels}, "
            f"ng={self._ray_space.ng}, nx={self._ray_space.nx})"
        )
