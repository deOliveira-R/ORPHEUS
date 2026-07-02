"""SN operator algebra leaves — typed :class:`TimedFullField` contract.

Provides the four-operator algebra leaves consumed by the within-group
equation :math:`A_{\\rm wg} = L + C - S_{\\rm foldable}`:

* :class:`StreamingOperator` — :math:`L = \\Omega\\cdot\\nabla
  + \\text{angular redistribution}` (the curvilinear pole term lives
  here for sphere / cylinder).
* the collision multiplier :math:`C = M[\\sigma_t]` (a plain
  :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
  — diagonal in position, group, and direction; #261 retired the former
  ``CollisionOperator`` thin subclass).
* :class:`InvertibleOperator` — the sweep-invertible specialisation
  :math:`(L + C)` returned by ``L + C``; ``is_invertible=True``
  via the WDD sweep.

All three operators consume and emit
:class:`~orpheus.transport.timed_full_field.TimedFullField` — the
typed composite carrier (bulk = :class:`~orpheus.transport.fields.angular_flux.AngularFlux`,
boundary = :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`).
Producer-side normalisation (Pattern 7): the typed contract is
enforced at every operator entry; no bare-ndarray packed-vector
adapter.  The geometry-agnostic 1-D matvec kernel is
:meth:`~orpheus.sn.loss_representation._OneDimScanWalk._apply_walk`
(1-D slab / sphere / cylinder; the fused ``(L+C)ψ`` single-emission
body, the apply-direction twin of the sweep; #206 Phase C moved it off
this operator INTO the representation).  The curvilinear Morel–Montry
angular redistribution is computed IN-SWEEP there (it is not a separable
typed leaf — the #238-retired ``M_spatial`` / ``M_angular_redist`` split
had no production consumer).  The multi-D Cartesian matvec is the
representation's ``loss_action`` walk that S6.3 moved OFF this operator
(``orpheus.sn.loss_representation``; production default ``ScanMarch``
since the S6.9 Fork-B2 flip, with ``MovingFrontierWindow`` a selectable
peer).

Three geometries are supported:

* **Cartesian 2D** — ``L = μ_x ∂/∂x + μ_y ∂/∂y + Σ_t``
* **Spherical 1D** — ``L = μ (A ∂/∂r)/V + (α ∂/∂μ)/V + Σ_t``
* **Cylindrical 1D** — per-level azimuthal redistribution

History
-------

* Wave E retired the standalone ``build_rhs*`` and
  ``build_transport_linear_operator*`` helpers along with the
  ``angular_flux_to_scalar`` aggregator.
* D-H (2026-05) — typed :class:`TimedFullField` composite carrier
  replaced the legacy packed-vector convention; operators bridged
  the bare-ndarray↔typed boundary internally.
* D-I (2026-05-29) — the bare-ndarray adapter retired from every
  operator leaf (L / C / S / F).  Operators consume only
  :class:`TimedFullField`.
* D-J (2026-05-30) — the supporting packed-vector codec family
  retired: :class:`EquationMap`, :func:`build_equation_map_*`,
  :func:`solution_to_angular_flux*`, :func:`pack_with_traces`.  The
  legacy 2-D Cartesian FD matvec ``transport_operator_matvec``
  (and its cartesian / spherical / cylindrical predecessors) had
  retired in D-H.2-C4e.6 (commit ``a614610``).

.. note:: Symmetric-closure invariant

   The operator :math:`L` uses the **symmetric** closure that makes
   the Krylov path agree with analytical references:

   * **Cartesian**: upwind cell-center finite differences for the
     streaming gradient — first-order accurate and consistent with DD
     on uniform meshes, first-order inconsistent on non-uniform meshes.
   * **Curvilinear**: arithmetic averages for spatial face fluxes and
     τ-weighted interpolation for angular face fluxes — distinct from
     the WDD asymmetric closure :math:`\\psi_{\\rm out} =
     (\\overline{\\psi} - (1 - \\tau)\\,\\psi_{\\rm in})/\\tau` used
     by the sweeps.

   On uniform meshes the symmetric-closure operator and the WDD
   sweep converge to the same physics in the fine-mesh limit; on
   curvilinear, the sweep's WDD closure has the ERR-026 closure-bias-
   driven self-consistent fixed point that the Krylov-on-:meth:`apply`
   path bypasses.

.. note:: Boundary-condition handling (Wave O step O.4a.2, Issue #208)

   The realized boundary law ``B`` is a **first-class sibling
   operator** of ``L``, NOT re-applied inside this matvec.  The
   canonical SN loss is ``(L_full + C - S - F - B)`` on the direct-sum
   state ``V = V_bulk ⊕ V_inflow ⊕ V_outflow``.  For the **1-D** path
   (slab / sphere / cylinder), the representation's ``loss_action``
   (:meth:`~orpheus.sn.loss_representation._OneDimScanWalk.loss_action`,
   which #206 Phase C moved off this operator)
   reads ``psi.boundary.inflow`` as a GIVEN, keeps the outflow
   self-consistency defect ``psi.outflow - streamed`` on the outflow
   trace row, and adds the inflow identity ``I·psi.inflow`` on the
   inflow row — with NO ``bc.apply``.  The reflective coupling
   ``psi.inflow = B·psi.outflow`` is delivered by the sibling ``-B``
   (:class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`), and the
   outer Krylov / SI loop drives the boundary consistency residual
   ``psi.inflow - B·psi.outflow - q.inflow → 0``.  See
   :ref:`bc-extraction` for the full block-matrix derivation, the three
   design corrections, the two ``-B`` delivery routes, and the O.2
   forcing function.

   The **multi-D Cartesian** path (the representation's ``loss_action``,
   which S6.3 moved off this operator — ``ScanMarch`` default since S6.9,
   ``MovingFrontierWindow`` peer) is ALSO bare (O.4b
   Phase E landed): it seeds the octant-incoming face slots from the
   GIVEN ``psi.boundary.inflow`` via the typed ``wavefront.seed`` (ι_*)
   with NO ``bc.apply``, walks the same per-octant
   :class:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph` (the apply-direction
   level operation → the diamond-difference ``DiscretizationScheme`` closure) the multi-D *sweep*
   ``_sweep_jacobi`` uses — so matvec ≡ sweep in 2-D by
   construction (L21, one discretization) — and emits the boundary
   consistency residual (outflow defect ``streamed − given`` + inflow
   identity) as a :class:`BoundarySourceSink`.  The reflective coupling
   is the sibling ``-B`` exactly as in 1-D.  The pre-extraction Phase C
   insight that the BC must consume the WDD-propagated outflow face
   vector (not cell centres) is preserved and strengthened: post-O.4a.2
   the outflow trace is the explicit solved unknown ``psi.outflow`` that
   ``-B`` reads, closing ERR-026 by construction for the 1-D curvilinear
   path.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Optional, cast

import numpy as np

from functools import cached_property

from orpheus.numerics.operator import (
    BlockRole,
    LinearOperator,
    OperatorSum,
)

from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator

if TYPE_CHECKING:
    from collections.abc import Callable

    from orpheus.transport.fields.angular_flux import AngularFlux
    from orpheus.transport.fields.boundary_flux import BoundaryFlux
    from orpheus.transport.fields.cross_section_field import CrossSectionField
    from orpheus.transport.full_field import FullField
    from orpheus.transport.timed_full_field import TimedFullField
    from orpheus.numerics.space import FunctionSpace
    from ..mesh.augmented_mesh import SNMesh
    from orpheus.numerics.frame import FrameBase
    from orpheus.transport.source_sinks import ScalarSourceSink, AngularSourceSink
    from ..loss_representation.sweep_schedule import SweepSchedule
    from ..spatial.pole_angular_closure import PoleAngularClosure
    from ..loss_representation import LossRepresentation
    # Type-only (the runtime constructions are late imports inside ``inverse``
    # / ``__sub__`` to break the operator <-> composite import cycles).
    from .boundary import SNMaskedBoundaryOperator
    from .scheduled_invertible import ScheduledInvertibleOperator
    from .sweep_operator import SweepOperator

__all__ = [
    "StreamingOperator",
]
# #206 Phase C / #238: the 1-D matvec walk lives in `_OneDimScanWalk`
# (orpheus.sn.loss_representation) — the fused `(L+C)ψ` (`loss_action`)
# rides its `_apply_walk` core (Cardinal Rule 2, one source). The
# production hot path is `StreamingOperator.apply` → `loss_action`
# (single emission). External consumers call `(L + C).apply(state)` via
# the public operator-algebra path. The former separately-applicable
# `(M_spatial, M_angular_redist)` operator-leaf split (`_SpatialSweepDirection`
# / `_MSpatialOperatorSum` / `AngularRedistributionOperator`) was retired in
# #238 — it had no production consumer; the curvilinear Morel–Montry angular
# redistribution is computed IN-SWEEP inside `_apply_walk`, not as a separable
# typed leaf, and is verified by the anisotropic curvilinear MMS
# (`tests/sn/verification/mms/test_curvilinear_aniso_convergence.py`).


def _require_typed_composite(
    method: str, sn_mesh: "SNMesh", field: "FullField",
) -> None:
    r"""The shared matvec input contract — timeless composite + mesh identity.

    The two guards :meth:`StreamingOperator.apply` introduced (D-I.3d typed
    contract + the mesh-identity invariant) are now consumed by EVERY SN
    matvec entry that takes a :class:`FullField`:
    :meth:`StreamingOperator.apply` / :meth:`apply_transpose` AND the #240
    Phase 2 Step B :meth:`InvertibleOperator.apply` / :meth:`apply_transpose`
    overrides.  Single source of the contract (``coding-elegance`` Pattern 2 /
    Pattern 4 — illegal inputs unrepresentable at one place, not re-validated
    per leaf).

    #257 S8a — the input contract is the TIMELESS
    :class:`~orpheus.transport.full_field.FullField` (the matvec leaves are
    base arrows ``FullField -> FullField``; only the iteration driver carries
    the history-bearing :class:`TimedFullField` comonad).  A
    :class:`TimedFullField` IS a :class:`FullField`, so the SI / Krylov drivers
    that pass a timed iterate still satisfy the guard; a bare ndarray does not.

    Parameters
    ----------
    method : str
        Qualified method name for the error message (e.g.
        ``"StreamingOperator.apply"``).
    sn_mesh : SNMesh
        The operator's mesh — ``field.bulk.mesh`` must be the SAME instance.
    field : FullField
        The matvec input (``psi`` for apply, ``phi`` for the transpose).
        A timeless :class:`FullField` or its timed subclass.
    """
    from orpheus.transport.full_field import FullField

    if not isinstance(field, FullField):
        raise TypeError(
            f"{method}: expected FullField, got "
            f"{type(field).__name__}.  D-I.3d (2026-05-29) retired the "
            "bare-ndarray packed-vector contract; construct a timeless "
            "composite via ``FullField(bulk=AngularFlux(...), "
            "boundary=BoundaryFlux(...))`` (or the timed "
            "``TimedFullField(bulk=..., boundary=...)`` for an iterate)."
        )
    if sn_mesh is not field.bulk.mesh:
        raise ValueError(
            f"{method}: operator and composite must share the SAME "
            "SNMesh instance (mesh-identity invariant)."
        )


@dataclass
class StreamingOperator(LinearOperator["FullField"]):
    r"""Pure streaming + angular-redistribution operator :math:`L` as a
    :class:`~orpheus.numerics.operator.LinearOperator` leaf.

    The "L" of the Phase G four-operator algebra
    :math:`A_{\rm wg} = L + C - S_{\rm foldable}`. Carries the spatial
    streaming math plus the curvilinear angular redistribution — the
    cell-collision term lives in the separate collision multiplier
    :math:`C = M[\sigma_t]` (a
    :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`).
    The split lets the within-group operator be pure algebra
    (``L + C - S.foldable_part()``); no ``WithinGroupOperator`` wrapper.

    Pure σ-free streaming ``apply`` (#257 S8b)
    ------------------------------------------

    .. math::

        L.{\rm apply}(\psi) \;:=\; \Omega\cdot\nabla\psi \;=\;
            \text{streaming\_action}(\psi)

    :math:`L` computes pure streaming **directly, with no
    :math:`\sigma`**: :meth:`apply` calls the
    :attr:`loss_representation`'s named σ-free
    :meth:`~orpheus.sn.loss_representation.LossRepresentation.streaming_action`
    leaf (the ONE streaming discretization, single-sourced through
    ``loss_action`` at :math:`\sigma = 0`).  The collision diagonal
    :math:`C = M[\sigma_t]` is the separate shared collision multiplier
    (a :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`),
    and the composition

    .. math::

        (L + C).{\rm apply}(\psi) \;=\; \text{streaming\_action}(\psi)
            \;+\; \sigma_t \odot \psi \;=\; M(\psi;\;\sigma_t)

    recovers the full within-group loss (the WDD matvec is affine in
    :math:`\sigma` in the forward direction — see
    :meth:`~orpheus.sn.loss_representation.LossRepresentation.streaming_action`).

    Why :math:`L` carries NO σ (Pattern 4 — #257 S8b)
    -------------------------------------------------

    The DISCRETE curvilinear matvec threads :math:`\sigma_t` through the
    Carlson coupled-pole seed (Hébert §3.9.4 Eqs. 3.432-3.435), but that
    :math:`\sigma`-dependence is *exactly the collision diagonal* the seed
    injects — it cancels into the :math:`\sigma\cdot\psi` term and belongs
    to :math:`C`, not :math:`L`.  Probe-confirmed σ-freedom: ``streaming\_action``
    is byte-stable (to ≤ a few ULP) across wildly different :math:`\sigma`
    fields (#257 S8b).  The CONTINUOUS :math:`L`
    (:math:`\Omega\cdot\nabla\psi + (1-\mu^2)/r\,\partial\psi/\partial\mu`)
    is σ-independent, and so now is the discrete leaf.  Pattern 4 (illegal
    states unrepresentable): ``StreamingOperator(sn_mesh)`` takes ONLY the
    mesh — a σ on :math:`L` would be a parameter the leaf never reads.

    Capability set
    --------------

    Pure streaming alone is **not
    invertible** (the streaming operator is rank-deficient without a
    collision term to make the within-group cell balance non-singular).
    The ``solve`` capability appears at the
    :class:`~orpheus.numerics.operator.OperatorSum` level: ``(L + C
    - S_foldable).solve(q)`` would route to the within-group sweep via a
    σ_r fusion hook — but ⚠ that σ_r-sweep is exact ONLY for isotropic flux
    (it removes the diagonal-in-angle ``Σ_s0·I``, not the isotropic-projection
    ``Σ_s0·P_iso``); wiring it as a within-group accelerator ships 46–56 %
    errors on anisotropic problems (issue #215; the stable+correct fold is
    consistent DSA #2 / Krylov). ``apply_transpose`` IS available
    (Wave O / O.2b) — the analytic reverse-direction adjoint matvec
    :math:`L^{\mathsf T}` (see :meth:`apply_transpose`), so the operator
    carries a working ``apply_transpose`` and ``L.H`` is the physical G-adjoint.

    Parameters
    ----------
    sn_mesh : SNMesh
        The augmented geometry carrying quadrature, BCs (the
        face-name-keyed ``sn_mesh.bc`` dict), pole closure, and (for
        curvilinear) the precomputed connection coefficients — no
        ``boundary`` constructor parameter.  The SOLE parameter: pure
        :math:`L` reads no :math:`\sigma` (#257 S8b).

    Notes
    -----
    Depth B D-I.3d (2026-05-29) — :meth:`apply` accepts ONLY
    :class:`~orpheus.transport.timed_full_field.TimedFullField`.
    Depth B D-J (2026-05-30) — the underlying packed-vector codec
    family (``EquationMap`` / ``build_equation_map_*`` /
    ``solution_to_angular_flux_*`` / ``pack_with_traces``) retired
    too.  Producer-side normalisation (Pattern 7) at the operator:
    the entire matvec consumes / emits the typed direct-sum carrier,
    never any legacy packed flat layout.
    """

    sn_mesh: "SNMesh"

    # Streaming is the sole FULL operator — it couples bulk ↔ boundary
    # (reads the inflow trace to seed the sweep, writes the outflow
    # trace). Issue #208 / Wave O. Class-level constant (unannotated so
    # the dataclass does not treat it as a field).
    # apply_transpose (Wave O / O.2b): the analytic reverse-direction
    # adjoint matvec landed — see :meth:`apply_transpose`.
    block_role = BlockRole.FULL

    # D-J (2026-05-30): ``_eq_map`` / ``_ensure_eq_map`` / ``n_unknowns``
    # retired alongside the :class:`EquationMap` codec family — the
    # typed :class:`~orpheus.transport.timed_full_field.TimedFullField`
    # contract has no need for the legacy packed-vector slot map.

    @property
    def is_adjointable(self) -> bool:
        # The analytic reverse-direction adjoint matvec L^T landed (Wave O /
        # O.2b — see :meth:`apply_transpose`).
        # is_invertible inherits base False — pure streaming L is not
        # sweep-invertible; only the (L+C) InvertibleOperator is.
        return True

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        r"""The composite carrier :math:`V_{\rm bulk}\oplus V_{\rm trace}` (Wave O / O.2b).

        :math:`L` is the sole FULL operator — it couples bulk :math:`\leftrightarrow`
        boundary (seeds the sweep from the inflow trace, emits the outflow
        trace). Advertising :attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.full_field_space`
        is what lets :class:`~orpheus.numerics.operator._AdjointOperator`
        read the **block-diagonal G-adjoint metric** (bulk :math:`V\,w_n`
        :math:`\oplus` trace :math:`|\Omega\cdot\hat n|\,w_n`) for ``L.H`` —
        without it the adjoint silently reduces to the metric-blind Euclidean
        transpose (Issue #208 risk R5).

        Since P4.5 W-D, ``C`` / ``S`` / ``F`` advertise the SAME composite
        domain (the cross-method composition-guard close — see
        :attr:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator.domain`),
        so the within-group
        ``(L + C) - S - B`` :class:`~orpheus.numerics.operator.OperatorSum`
        guard VALIDATES the build rather than silently skipping the formerly
        ``None``-spaced summands. ``L``'s composite domain therefore agrees
        with every summand, and the **transpose-closed** sub-sums whose ``.H``
        is actually reachable — ``(L + C)`` and ``(L + C - B)`` — G-conjugate
        every bulk leaf for free via the op-level
        :math:`G^{-1}(\sum \text{leaf}^{\mathsf T})G`. The full prompt loss
        ``(L + C - S - F - B).H`` stays **intentionally unreachable**, and W-D
        does NOT change that — the unreachability is enforced by the CAPABILITY
        lattice, not by ``None`` spaces: ``S`` / ``F`` advertise no
        ``apply_transpose``, so ``OperatorSum`` does not propagate
        a working ``apply_transpose`` and
        :class:`~orpheus.numerics.operator._AdjointOperator` raises
        :class:`MissingAdjoint` (fails loud, never silently Euclidean —
        the capability lattice makes the metric-blind state unrepresentable).
        The foldable scattering / fission contributions are handled at the
        eigenvalue / DSA outer layer, not via this within-group adjoint.
        """
        return self.sn_mesh.full_field_space

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        # Endomorphism on the composite (see :meth:`domain`).
        return self.sn_mesh.full_field_space

    def apply(self, psi: "FullField") -> "FullField":
        r"""Pure σ-free forward streaming :math:`L\,\psi = \Omega\cdot\nabla\psi`.

        Computes pure streaming DIRECTLY via the :attr:`loss_representation`'s
        named σ-free
        :meth:`~orpheus.sn.loss_representation.LossRepresentation.streaming_action`
        leaf (#257 S8b — the ``(L+C)−C`` fold is retired).  The streaming
        discretization is single-sourced through ``loss_action`` at
        :math:`\sigma = 0` (the WDD matvec is affine in :math:`\sigma`);
        :math:`L` reads NO :math:`\sigma`.  The collision diagonal
        :math:`C = M[\sigma_t]` is the separate shared collision multiplier
        (a :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`);
        the composition
        :math:`(L + C).{\rm apply}(\psi) = \text{streaming\_action}(\psi)
        + \sigma_t\odot\psi = M(\psi;\sigma_t)` recovers the full loss.

        ``L`` is the ONLY operator that emits a non-zero face
        residual on its output ``.boundary`` — the collision multiplier
        :math:`C = M[\sigma_t]`
        (a :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`),
        :class:`~orpheus.transport.operators.scattering.ScatteringOperator`, and
        :class:`~orpheus.transport.operators.fission.FissionOperator` all leave the
        output boundary at the auto-allocated zero.

        Parameters
        ----------
        psi : FullField
            Composite carrier with bulk
            (:class:`~orpheus.transport.fields.angular_flux.AngularFlux`)
            and boundary
            (:class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`).
            Operator and ``psi.bulk.mesh`` MUST be the same
            :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` instance.

        Returns
        -------
        FullField
            ``L·ψ`` as a timeless composite — bulk carries the pure
            streaming cell action, boundary carries the face residual at
            the layout-assigned trace slots (non-zero at outer face for
            curvilinear; non-zero at outer + inner faces for slab).
            History-free (#257 S8a — the matvec leaf is a base arrow
            ``FullField -> FullField``; the comonad lives on the driver).
        """
        _require_typed_composite("StreamingOperator.apply", self.sn_mesh, psi)
        return self.loss_representation.streaming_action(psi)

    def apply_transpose(self, phi: "FullField") -> "FullField":
        r"""Hilbert transpose :math:`L^{\mathsf T}\,\phi` (Wave O / O.2b, #208).

        The σ-free adjoint streaming leaf (#257 S8b): :math:`L^{\mathsf T}\phi`
        via the :attr:`loss_representation`'s named
        :meth:`~orpheus.sn.loss_representation.LossRepresentation.streaming_action_transpose`
        (single-sourced through ``loss_action_transpose`` at :math:`\sigma = 0`;
        the curvilinear reverse walk carries the angular second triangular
        factor, the multi-D Cartesian adjoint stays a deferral raise — never a
        silent wrong answer).  Since :math:`C = \sigma_t\odot` is a self-adjoint
        diagonal, the full adjoint loss factors as
        :math:`(L + C)^{\mathsf T} = L^{\mathsf T} + C` — the collision
        multiplier :math:`C = M[\sigma_t]`
        (a :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`)
        is self-adjoint, so its transpose is the same multiplier.

        This returns the **plain Euclidean transpose** :math:`L^{\mathsf T}`.
        The metric conjugation :math:`G^{-1}\!\cdot^{\mathsf T}\!\cdot G` of the
        physical **G-adjoint** ``L.H`` is applied AROUND this by
        :class:`~orpheus.numerics.operator._AdjointOperator`, which reads the
        ``domain`` / ``codomain`` ``inner_product_weights`` (bulk volume on the
        cell block, the ``|Ω·n|·w`` partial-current metric on the trace block).

        Verified by the G-adjoint reciprocity gate
        ``test_g_adjoint_reciprocity`` (slab / sphere / cylinder, -O-firing) +
        its L11 wrong-trace-metric negative control.
        """
        _require_typed_composite(
            "StreamingOperator.apply_transpose", self.sn_mesh, phi,
        )
        return self.loss_representation.streaming_action_transpose(phi)

    # ── LossRepresentation carve (S2) — the polymorphic matvec dispatch ─────

    @cached_property
    def loss_representation(self) -> "LossRepresentation":
        r"""THE loss-operator representation for this operator's mesh (S6.5).

        The ONE first-class ``LossRepresentation``
        (``orpheus.sn.loss_representation``) carrying BOTH actions of
        :math:`(L+C)`: :meth:`apply` routes through
        ``representation.loss_action`` / ``loss_action_transpose`` (the
        matvec), and :meth:`InvertibleOperator.solve` runs the forward
        substitution on the SAME object via
        :attr:`InvertibleOperator.loss_representation` — L21 ("matvec ≡
        sweep") as a type fact.  Selection is by geometry
        (``default_for``): 1-D → ``CumprodScan``; multi-D Cartesian →
        ``ScanMarch`` (the S6.9 Fork-B2 default).  ``cached_property`` because
        the selection is fixed by the mesh, stable across the operator's
        lifetime; the lazy import breaks the operator ↔ loss_representation
        module cycle.
        """
        from ..loss_representation import default_for

        return default_for(self.sn_mesh)

    # ── Algebra dispatch — sweep-invertible composite (R-1 Step C) ────

    def __add__(self, other):
        r"""Compose :math:`L + X`.

        When ``X`` is a
        :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
        (the collision diagonal :math:`C = M[\sigma_t]`), returns the
        sweep-invertible specialisation :class:`InvertibleOperator`
        carrying the algebraic identity :math:`(L + C)^{-1} \approx
        \text{WDD sweep}`.  Otherwise falls through to the generic
        :class:`OperatorSum` via the mixin.

        #261: ``L + C`` is the canonical (and only) ordering — the dispatch
        lives here on the SN-specific streaming leaf, because the
        transport-level multiplier cannot dispatch back onto ``StreamingOperator``
        (that would be a ``transport → sn`` upward import).
        """
        if isinstance(other, MultiplicationOperator):
            return InvertibleOperator(self, other)
        return super().__add__(other)


# ─────────────────────────────────────────────────────────────────────────
# InvertibleOperator — sweep-invertible composite (L + C)
# ─────────────────────────────────────────────────────────────────────────


class InvertibleOperator(OperatorSum["FullField"]):
    r"""Sweep-invertible composite :math:`L + C` carrying ``.solve`` = WDD sweep.

    R-1 Step C (2026-05-19) — the SN-specific algebraic identity

    .. math::

        (L_{\rm streaming} + C_{\rm diagonal})^{-1} \;\approx\;
        \text{WDD sweep}

    has no generic ``(A+B)^{-1}`` formula — a plain
    :class:`OperatorSum` can only invert ITERATIVELY (the
    preconditioned-splitting
    :class:`~orpheus.numerics.green_operator.GreenOperator` its generic
    ``.inverse()`` returns, taxonomy §12 step 4).  The WDD sweep IS the
    DIRECT inverse algorithm for this specific composite — that's the
    algebraic foundation of the entire SN method (Lewis & Miller §3.2;
    Adams & Larsen 2002 §III) — and this subclass's ``.inverse()``
    override (→ :class:`~orpheus.sn.operators.sweep_operator.SweepOperator`)
    shadows the generic Green by MRO: the type IS the structure
    (taxonomy §11.1).  :class:`InvertibleOperator` is the specialisation
    that carries the identity at the type level: it OWNS its full action
    algebra.  :meth:`apply` / :meth:`apply_transpose` OVERRIDE the
    :class:`OperatorSum` leaf-sum to return the within-group loss
    :math:`(L+C)\psi = M(\sigma)\psi` (and its transpose) DIRECTLY via
    :attr:`loss_representation`, single-sourcing :math:`\sigma` from the
    diagonal — the SAME :math:`\sigma` ``solve`` threads into the WDD
    sweep (#240 Phase 2 Step B).  ``solve`` is the forward substitution
    on that SAME
    :class:`~orpheus.sn.loss_representation.LossRepresentation` instance
    (S6.5, #222), so matvec, adjoint, and sweep are three actions of ONE
    operator (L21).

    Construction
    ============

    Two equivalent paths:

    * **Operator algebra dispatch** — ``L + C`` where ``L`` is a
      :class:`StreamingOperator` and ``C`` is the collision multiplier
      :math:`M[\sigma_t]`
      (a :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`)
      returns this class automatically.  ``L + C`` is the canonical (and
      only) ordering: the dispatch lives one-directionally on
      :meth:`StreamingOperator.__add__`, because a transport-level
      multiplier cannot dispatch back onto an ``sn`` operator (#261).  The
      composite reads as math.
    * **Explicit construction** — ``InvertibleOperator(L, C)``.  Useful
      when composing variants such as
      ``InvertibleOperator(L_leaf, MultiplicationOperator.from_mesh(σ_r, mesh))``
      where
      ``σ_r = σ_t - Σ_{s,0}^{g→g}`` is the removal cross-section that
      lets one fold the within-group self-scatter into the diagonal
      collision term (Adams & Larsen 2002 §III; tracked by issue
      `#200 <https://github.com/deOliveira-R/ORPHEUS/issues/200>`_).

    The two paths produce structurally identical objects — the choice
    only changes the call-site readability.

    Capability set
    ==============

    ``is_invertible=True`` — adds
    ``solve`` (the WDD sweep) to the parent :class:`OperatorSum`'s set;
    ``apply_transpose`` propagates by the :class:`OperatorSum` closure
    law (both :math:`L` and :math:`C` advertise it) and is OVERRIDDEN to
    the composite's own :math:`M(\sigma)^{\mathsf T}` action (Wave O #208 /
    #240 Step B).  The multi-D Cartesian adjoint raises (the
    representation's deferral contract — never a silent wrong answer).

    The ``.solve`` API
    ==================

    The ``rhs`` parameter is the timeless composite
    :class:`~orpheus.transport.full_field.FullField` (W-C, P4.5).  The
    history-bearing :class:`~orpheus.transport.timed_full_field.TimedFullField`
    iterate passes through by inheritance (it IS a ``FullField``), and a
    bare ``FullField`` is admitted as the ``history_depth = 0``
    degenerate.  ``rhs`` carries:

    * ``rhs.bulk.values`` — per-ordinate source ``(N, ng, nx, ny)``.
      This is treated as the per-ordinate anisotropic source
      :math:`Q^{\rm aniso}` that the sweep consumes (the isotropic
      source is zero).
    * ``rhs.boundary`` — face source / BC inflow trace.  Typically
      zero for volumetric SI/Krylov sources (which carry no face
      contribution); the persistent reflective-BC state lives on the
      :class:`SNMesh` and is handled inside the sweep.  It seeds the
      sweep's mutable boundary buffer (per-face copy) and falls back
      as the inflow seed when no ``initial_guess`` is supplied.

    The previous iterate :math:`\psi^{(k-1)}` is passed via the explicit
    ``initial_guess`` keyword (NOT piggy-backed on a lag-1 frame of
    ``rhs``).  The curvilinear sweep reads the Carlson coupled-pole seed
    from its ``.bulk.values`` (container-agnostically — see
    :func:`~orpheus.sn.loss_representation._initial_guess_values`).
    ``None`` (cold start) → the sweep falls back to its
    in-iteration-source default.

    Parameters
    ----------
    streaming : StreamingOperator
        :math:`L = \Omega\cdot\nabla + \text{angular redistribution}`.
        Resolution A subtractive form:
        ``L.apply(ψ) = M(ψ; σ_t) - σ_t ⊙ ψ_cell``.
    diagonal : MultiplicationOperator
        :math:`C = M[\sigma]`.  Its ``.coefficient.values`` is the
        per-cell per-group coefficient used by the sweep (canonically
        ``σ_t``; can be ``σ_r`` for the foldable variant).

    Notes
    -----
    The validation ``σ > 0`` everywhere guards against the
    ``σ_r < 0`` case that can arise when within-group self-scatter
    exceeds total cross-section (rare; not physically meaningful but
    mathematically possible for ill-conditioned multi-group sets).
    The sweep would emit NaN at those cells — surfacing the
    inconsistency at construction is friendlier.
    """

    def __init__(
        self,
        streaming: "StreamingOperator",
        diagonal: "MultiplicationOperator",
    ) -> None:
        if not isinstance(streaming, StreamingOperator):
            raise TypeError(
                f"InvertibleOperator: 'streaming' must be a "
                f"StreamingOperator; got {type(streaming).__name__}."
            )
        if not isinstance(diagonal, MultiplicationOperator):
            raise TypeError(
                f"InvertibleOperator: 'diagonal' must be a "
                f"MultiplicationOperator; got {type(diagonal).__name__}."
            )
        # Mesh-identity invariant (#261 re-anchor): the WDD sweep threads the
        # diagonal's σ against the STREAMING geometry, so the two must act on
        # the SAME mesh object — geometric consistency, strictly stronger than
        # the (name, shape) shape-equality the OperatorSum composition guard
        # checks. The diagonal multiplier is mesh-free; its mesh is carried by
        # its CrossSectionField coefficient.
        if streaming.sn_mesh is not diagonal.coefficient.mesh:
            raise ValueError(
                "InvertibleOperator: the streaming operator and the diagonal "
                "multiplier must act on the same mesh instance — the "
                "mesh-identity invariant (streaming.sn_mesh is "
                "diagonal.coefficient.mesh): the WDD sweep pairs the "
                "diagonal's σ with the streaming geometry."
            )
        if not np.all(diagonal.coefficient.values > 0):
            min_sigma = float(np.min(diagonal.coefficient.values))
            raise ValueError(
                f"InvertibleOperator: diagonal coefficient must be "
                f"strictly positive everywhere for the WDD sweep to be "
                f"well-defined; got min(sigma) = {min_sigma:.3e}.  If "
                f"sigma_r = sigma_t - Sigma_(s,0)^(g->g) is dipping "
                f"negative, the multi-group cross-section set is "
                f"physically inconsistent."
            )
        super().__init__(streaming, diagonal)
        # block_role is now DERIVED by OperatorSum.__init__ (Wave O / O.2b 4.5):
        # join(L=FULL, C=BULK) = FULL. The former hand-stamp here was the
        # twin-path retired in 4.5 — the role is carried by construction.

    @property
    def is_invertible(self) -> bool:
        # (L+C) is sweep-invertible: the WDD forward-substitution sweep IS its
        # inverse operator (in
        # __init__). This is the SOLE invertible OperatorSum — the base
        # OperatorSum.is_invertible is False. is_adjointable inherits the
        # OperatorSum a∧b law (both L and C advertise apply_transpose).
        return True

    # ── Convenience accessors ─────────────────────────────────────────

    @property
    def streaming(self) -> "StreamingOperator":
        """The streaming operand (alias for ``self.a``)."""
        return cast("StreamingOperator", self.a)

    @property
    def diagonal(self) -> "MultiplicationOperator":
        """The diagonal-multiplier operand :math:`C = M[\\sigma]` (alias for ``self.b``)."""
        return cast("MultiplicationOperator", self.b)

    @property
    def loss_representation(self) -> "LossRepresentation":
        r"""The ONE :class:`LossRepresentation` for this operator (S6.5, #222).

        Delegates to the streaming leaf's cached instance — the SAME
        object :meth:`StreamingOperator.apply` consumes for the matvec
        :math:`(L+C)\psi`.  :meth:`solve` runs the forward substitution
        :math:`(L+C)^{-1}q` on it, so "matvec ≡ sweep — two actions of
        ONE operator" (L21) is a type fact enforced by construction,
        not a coincidence of two ``default_for`` calls agreeing.
        """
        return self.streaming.loss_representation

    @property
    def sn_mesh(self) -> "SNMesh":
        """The shared :class:`SNMesh` (validated mesh-identity at init)."""
        return self.streaming.sn_mesh

    @property
    def sigma(self) -> np.ndarray:
        r"""The diagonal coefficient used by ``solve`` (σ_t or σ_r).

        Single-sources :math:`\sigma` off the diagonal multiplier's
        :class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`
        coefficient (``coding-elegance`` Pattern 2 — no duplicate storage).
        """
        return self.diagonal.coefficient.values

    # ── apply / apply_transpose: the composite's OWN matvec (#240 Step B) ──

    def apply(self, psi: "FullField") -> "FullField":
        r"""Matvec :math:`(L+C)\,\psi = M(\sigma)\,\psi` — the composite OWNS it.

        #240 Phase 2 Step B.  Both the matvec and the sweep are actions of the
        ONE :math:`(L+C)` operator (L21 "matvec ≡ sweep"), realised with THIS
        composite's diagonal :math:`\sigma` (``self.sigma`` = the collision
        leaf's :math:`\sigma` — the SAME array :meth:`solve` threads into the
        WDD sweep).  The representation's :meth:`loss_action` returns the FULL
        within-group loss :math:`(L+C)\psi = M(\sigma)\psi` directly.

        This OVERRIDES the inherited :meth:`OperatorSum.apply` (``L.apply +
        C.apply``).  The leaf sum is value-equal *only by coincidence*: in the
        forward direction the WDD matvec is AFFINE in :math:`\sigma`
        (:math:`M(\sigma)\psi = \text{streaming\_action}(\psi) + \sigma\cdot\psi`),
        so ``L.apply(σ_t) + C.apply(σ_r)`` collapses to
        :math:`\text{streaming\_action}(\psi) + \sigma_r\cdot\psi = M(\sigma_r)\psi`
        — the right value, but sourcing :math:`\sigma` from ``L.sigma_t`` (the
        streaming leaf) while :meth:`solve` sources it from ``C``.  Two sources
        that agree only because production builds :math:`L` and :math:`C` from
        the same :math:`\sigma_t`.  The override single-sources :math:`\sigma`
        from the diagonal (``coding-elegance`` Pattern 2: one ``loss_action``,
        one source of :math:`\sigma`), removing the latent affine-in-:math:`\sigma`
        coupling — the composite never asks the leaf for a :math:`\sigma`-bearing
        action it must then undo.
        """
        _require_typed_composite("InvertibleOperator.apply", self.sn_mesh, psi)
        return self.loss_representation.loss_action(self.sigma, psi)

    def apply_transpose(self, phi: "FullField") -> "FullField":
        r"""Adjoint matvec :math:`(L+C)^{\mathsf T}\,\phi = M(\sigma)^{\mathsf T}\,\phi`.

        The adjoint sibling of :meth:`apply` (#240 Phase 2 Step B): the
        representation's :meth:`loss_action_transpose` realises
        :math:`(L+C)^{\mathsf T}\phi = M(\sigma)^{\mathsf T}\phi` directly with
        THIS composite's diagonal :math:`\sigma`, overriding the inherited
        :meth:`OperatorSum.apply_transpose` leaf sum.  Multi-D Cartesian raises
        (the representation's deferral contract — never a silent wrong answer).
        The plain Euclidean transpose; the metric conjugation of the physical
        G-adjoint ``.H`` is applied AROUND this by
        :class:`~orpheus.numerics.operator._AdjointOperator` (pinned by
        ``test_g_adjoint_reciprocity``).
        """
        _require_typed_composite(
            "InvertibleOperator.apply_transpose", self.sn_mesh, phi,
        )
        return self.loss_representation.loss_action_transpose(self.sigma, phi)

    # ── Algebra dispatch — schedule-folded composite (#226 step 2) ────

    def __sub__(self, other):
        r"""Compose :math:`(L+C) - X`.

        When ``X`` is the strictly-lower boundary half ``B_lower`` (an
        :class:`~orpheus.sn.operators.boundary.SNMaskedBoundaryOperator`
        from :meth:`~orpheus.sn.operators.boundary.SNBoundaryOperator.split`),
        returns the sweep-invertible schedule-folded specialisation
        :class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`
        — the reified splitting matrix :math:`M = (L+C-B_{\rm lower})` whose
        ``solve`` is the octant-group forward substitution (§17 W2).
        Otherwise falls through to the generic difference via the mixin.

        ``(L+C) - B_lower`` is the canonical spelling — the dispatch lives
        here on the SN composite, mirroring :meth:`StreamingOperator.__add__`
        (#261: one-directional, the operand cannot dispatch back).
        """
        from .boundary import SNMaskedBoundaryOperator

        if isinstance(other, SNMaskedBoundaryOperator):
            from .scheduled_invertible import ScheduledInvertibleOperator

            return ScheduledInvertibleOperator(self, other)
        return super().__sub__(other)

    # ── solve: WDD sweep ─────────────────────────────────────────────

    def inverse(self) -> "SweepOperator":
        r"""Return the inverse OPERATOR :math:`(L+C)^{-1}` (the carve, #226).

        ``A.inverse().apply(b)`` is the WDD forward-substitution sweep,
        BIT-IDENTICAL to ``A.solve(b)`` — the returned
        :class:`~orpheus.sn.operators.sweep_operator.SweepOperator` delegates to
        :meth:`solve`. This is the operator normal form
        ``K = A_loss.inverse() @ F`` (Grand Report v3 §1): the forward view
        :meth:`apply` and the inverse view ``inverse().apply`` are the two views
        of ONE operator, the way ``A`` and ``A.H`` are. Coexists with
        :meth:`solve` through Phase 2–3; ``solve`` retires at Phase 4.

        **Forward-side back-half twin (collapse trigger).**  The
        ``is_invertible``/``inverse``/``solve`` back-half here is deliberately
        coextensive with
        :class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`'s
        (the schedule-folded sibling, #226 step 2) — two witnesses, kept
        twinned per defer-until-≥2.  TRIGGER: at the 3rd sweep-invertible
        FORWARD composite, extract a shared mixin; do not hand-re-derive it.
        (Distinct from the INVERSE-side twin noted on ``SweepOperator`` —
        Green/Matrix inverses grow that shape, not this one.)
        """
        from orpheus.sn.operators.sweep_operator import SweepOperator

        return SweepOperator(self)

    def solve(
        self,
        rhs: "FullField",
        *,
        initial_guess: "FullField | None" = None,
    ) -> "TimedFullField":
        r"""Invert :math:`(L + C)\,\psi = \text{rhs}` via the WDD sweep.

        The cell-balance equation
        :math:`(\Omega\cdot\nabla + \sigma)\,\psi = Q` is integrated
        cell-by-cell in inflow-to-outflow order; the angular closure
        (Cartesian → identity, curvilinear → Morel-Montry) is bound
        on the mesh.

        Parameters
        ----------
        rhs : FullField
            The timeless composite source (W-C, P4.5).  The
            history-bearing :class:`TimedFullField` iterate passes
            through by inheritance (it IS a ``FullField``), and a bare
            :class:`~orpheus.transport.full_field.FullField` is admitted
            as the ``history_depth = 0`` degenerate.  Carries:

            * ``bulk.values`` — per-ordinate source
              :math:`Q^{\rm aniso}`, shape ``(N, ng, nx, ny)``.
            * ``boundary`` — BC inflow trace (typically zero for
              SI/Krylov volumetric sources; falls back as the seed when
              no ``initial_guess`` is supplied).

            The legacy :class:`AngularFlux` arm is retired (it is NOT a
            ``FullField``, so the guard rejects it).
        initial_guess : FullField or None, optional
            Previous iterate :math:`\psi^{(k-1)}` for the curvilinear
            Carlson coupled-pole seed (M-M reads its ``.bulk.values`` as
            the level's angular flux to derive :math:`\bar Q` at
            :math:`\mu = -1`).  ``None`` (default) → cold start: M-M's
            Carlson seed degenerates to the zero-input result (same
            as ``ZeroSeed`` ablation).

            Explicit kwarg (post-Phase-1.2) — the seed is no longer
            piggy-backed on a lag-1 frame of ``rhs`` (history threading
            is the driver iterate's concern, an unrelated matter).
            Outer iterators (:class:`SourceIteration`,
            :class:`KEigenvalue`) pass the previous iterate
            explicitly; GMRES residual calls pass ``None`` so the
            preconditioner doesn't silently route through stateful
            seed state (closes the R-1 Step D silent-fallback bug
            class).

        Returns
        -------
        TimedFullField
            The angular flux satisfying :math:`(L + C)\,\psi =
            \text{rhs}`, with the sweep's outflow face state in
            ``.boundary``.  The WDD sweep emits a
            :class:`TimedFullField` iterate (the genuine driver-side
            comonad carrier); its ``history_depth`` matches
            ``rhs.history_depth`` (0 when ``rhs`` is a bare
            ``FullField``), and ``_history`` is empty — the outer
            SI / Krylov loop owns history threading.
        """
        return self._solve_timed_full_field(
            rhs, initial_guess=initial_guess,
        )

    # ``solve_moments`` (Phase 5c) retired in #226 step 2 (§17 W1): a public
    # method whose output-mode argument silently changed the operator's
    # codomain was a composition wearing a config.  The moment-emitting entry
    # is now the typed windowed product ``P @ A.inverse()``
    # (:class:`~orpheus.sn.operators.windowing.WindowedSweep`), whose fused
    # ``apply`` calls :meth:`_solve_timed_full_field` with ``moment_frame`` —
    # the ONE application-context body, private.

    def _solve_timed_full_field(
        self,
        rhs: "FullField",
        *,
        initial_guess: "FullField | None" = None,
        moment_frame: "FrameBase | None" = None,
        schedule: "SweepSchedule | None" = None,
        reflect: "Callable[[BoundaryFlux, tuple[str, ...]], None] | None" = None,
    ) -> "TimedFullField":
        r"""Composite :class:`TimedFullField` body of :meth:`solve` (D-H.1c stage 1).

        Runs the WDD forward substitution on
        :attr:`loss_representation` — the operator's ONE
        :class:`~orpheus.sn.loss_representation.LossRepresentation`
        instance (S6.5, #222) — and handles the L2 field plumbing at
        the public-entry boundary.

        The boundary plumbing seeds the sweep's mutable write-through
        buffer from the source trace: ``boundary_buf =
        BoundaryFlux.zeros_on(sn_mesh)`` is filled per-face from
        ``rhs.boundary`` via the L2 ``face_view`` copy
        (``boundary_buf.face_view(face_name)[:] =
        seed_boundary.face_view(face_name)`` for each shared face —
        works for slab, curvilinear, and 2-D Cartesian).  The sweep
        mutates ``boundary_buf`` in place; the result is re-wrapped as
        an L2 composite at the end.

        Parameters
        ----------
        rhs : FullField
            Per-ordinate source on the composite carrier (the timed
            iterate passes via inheritance; a bare ``FullField`` is the
            ``history_depth = 0`` degenerate).
        initial_guess : FullField or None
            Previous iterate (carries the curvilinear Carlson seed via
            its ``.bulk.values``).  ``None`` → cold start.

        Returns
        -------
        TimedFullField
            Solve output with ``bulk`` = ``(L + C)^{-1} rhs.bulk`` and
            ``boundary`` = the sweep's outflow face state.
            ``history_depth`` matches ``rhs.history_depth``; ``_history``
            is empty (solver outputs carry no iteration history — the
            outer SI / Krylov loop owns history threading).
        """
        from orpheus.transport.fields.angular_flux import (
            AngularFlux,
        )
        from orpheus.transport.fields.boundary_flux import (
            BoundaryFlux,
        )
        from orpheus.transport.fields.harmonic_moment_flux import (
            HarmonicMomentFlux,
        )
        from orpheus.transport.full_field import FullField
        from orpheus.transport.timed_full_field import TimedFullField

        # W-C (P4.5): the operator boundary speaks the timeless
        # :class:`FullField` composite; the timed iterate passes via
        # inheritance (``TimedFullField`` IS a ``FullField``), and a bare
        # ``FullField`` is admitted as the ``history_depth=0`` degenerate.
        # Legacy :class:`AngularFlux` stays retired — it is NOT a
        # ``FullField``, so the guard still rejects it.  Single guard site
        # for both :meth:`solve` and :meth:`solve_moments`.
        if not isinstance(rhs, FullField):
            raise TypeError(
                f"InvertibleOperator: 'rhs' must be FullField; "
                f"got {type(rhs).__name__}.  Legacy AngularFlux retired "
                f"in D-H.2-C3."
            )
        if initial_guess is not None and not isinstance(
            initial_guess, FullField,
        ):
            raise TypeError(
                f"InvertibleOperator: 'initial_guess' must be "
                f"FullField or None; got "
                f"{type(initial_guess).__name__}."
            )

        sn_mesh = self.sn_mesh
        if rhs.bulk.mesh is not sn_mesh:
            raise ValueError(
                "InvertibleOperator.solve(FullField): rhs and "
                "operator must share the same SNMesh instance "
                "(mesh-identity invariant)."
            )
        if initial_guess is not None and initial_guess.bulk.mesh is not sn_mesh:
            raise ValueError(
                "InvertibleOperator.solve(FullField): initial_guess "
                "and operator must share the same SNMesh instance "
                "(mesh-identity invariant)."
            )

        # ── L2 boundary buffer for the sweep (D-H.2-C2) ───────────────
        #
        # The sweep mutates ``boundary_buf`` (the L2 mutable
        # write-through; ``frozen=True`` freezes field rebinding but
        # the underlying flat ndarray remains writable through
        # :meth:`face_view`).
        #
        # Wave O (#208) O.4a.2 — BARE SWEEP: the inflow seed is the
        # boundary SOURCE ``rhs.boundary`` (the inflow slots carry
        # ``q.boundary + B·ψ.outflow`` — the SI driver folds ``S + B`` so
        # the ``Bψ`` reflective inflow rides in ``rhs.boundary``).  This
        # REPLACES the pre-extraction seed-from-``initial_guess.boundary``:
        # the bare sweep no longer re-applies ``bc`` to the iterate's
        # outflow, so the iterate's boundary is NOT the inflow seed.  The
        # iterate (``initial_guess``) still threads the BULK Carlson /
        # angular warm-start through the representation sweep below —
        # that path reads ``initial_guess.bulk``, not its boundary.
        boundary_buf = BoundaryFlux.zeros_on(sn_mesh)  # L2 after C2
        seed_boundary = rhs.boundary
        # Per-face copy via L2 face_view — works for slab (xmin, xmax),
        # curvilinear (xmax only), and 2-D Cartesian (all 4).
        for face_name in boundary_buf.layout.faces:
            if face_name in seed_boundary.layout.faces:
                boundary_buf.face_view(face_name)[:] = (
                    seed_boundary.face_view(face_name)
                )

        # ── Sweep on the operator's ONE representation (S6.5, #222) —
        # the SAME :class:`LossRepresentation` instance the matvec
        # (:meth:`StreamingOperator.apply`) consumes, so L21 ("matvec ≡
        # sweep") is a type fact, not two ``default_for`` calls agreeing.
        # ``rhs.bulk.values`` IS the per-ordinate source by producer
        # contract (R-1 Step 4 A1) — typed at the ``rhs`` guard above, so
        # no wrap-unwrap round trip through :class:`AngularSourceSink`
        # (the module-level :func:`transport_sweep` keeps that typed
        # boundary for operator-free callers).
        #
        # The composite ``initial_guess`` passes straight through —
        # D-H.1c stage 4: the sweep kernels read ``.bulk.values`` via the
        # container-agnostic :func:`_initial_guess_values` extractor.
        #
        # Phase 5c: ONE sweep for BOTH output modes — the moment
        # frame rides as an optional kwarg (the 1-D representation
        # raises on it, since moment output is 2-D Cartesian only;
        # ``moment_*`` and ``initial_guess`` are mutually 2-D-vs-1-D, so
        # the 2-D branch harmlessly drops the unused seed).  Only the
        # OUTPUT WRAP differs: full angular field vs harmonic moments.
        bulk_values, _scalar = self.loss_representation.sweep(
            rhs.bulk.values,
            self.sigma,
            boundary_buf,
            initial_guess=initial_guess,
            moment_frame=moment_frame,
            schedule=schedule,
            reflect=reflect,
        )
        # The sweep output carries the trailing 2^d spatial-moment axis at a
        # multi-moment closure (the φ̂ iterate, #240 D5b-S3); the typed wrap
        # selects the SpatialMomentSpace factor so the iterate is a legal typed
        # state.  DD/Step (per_axis == 1) → no factor, byte-identical.
        per_axis = sn_mesh.scheme.spatial_basis_per_axis
        if moment_frame is None:
            bulk = AngularFlux.from_mesh(
                bulk_values, sn_mesh, spatial_moments=per_axis,
            )
        else:
            # In moment mode the sweep returns the (L+1, 2L+1, ...) moment
            # tensor, so its own leading axis fixes L (no basis-specific read).
            bulk = HarmonicMomentFlux.from_mesh_and_L(
                bulk_values, sn_mesh, bulk_values.shape[0] - 1,
                spatial_moments=per_axis,
            )

        # ── L2 direct return — no adapter needed (D-H.2-C2). ───────────
        return TimedFullField(
            bulk=bulk,
            boundary=boundary_buf,
            _history=(),
            history_depth=(
                rhs.history_depth if isinstance(rhs, TimedFullField) else 0
            ),
        )
