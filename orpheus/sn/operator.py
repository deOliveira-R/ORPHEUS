"""SN operator algebra leaves — typed :class:`TimedFullField` contract.

Provides the four-operator algebra leaves consumed by the within-group
equation :math:`A_{\\rm wg} = L + C - S_{\\rm foldable}`:

* :class:`StreamingOperator` — :math:`L = \\Omega\\cdot\\nabla
  + \\text{angular redistribution}` (the curvilinear pole term lives
  here for sphere / cylinder).
* :class:`CollisionOperator` — :math:`C = \\sigma\\cdot` (diagonal
  in position, group, and direction).
* :class:`InvertibleOperator` — the sweep-invertible specialisation
  :math:`(L + C)` returned by ``L + C``; advertises ``CAP_SOLVE``
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
   (:class:`~orpheus.sn.boundary_operator.SNBoundaryOperator`), and the
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
   :class:`~orpheus.sn.sweep_graph.SweepDependencyGraph` (the apply-direction
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
from typing import TYPE_CHECKING, Optional

import numpy as np

from functools import cached_property

from orpheus.numerics.operator import (
    BlockRole,
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    CAP_SOLVE,
    LinearOperatorMixin,
    OperatorSum,
)

from orpheus.numerics.quadrature import Quadrature

if TYPE_CHECKING:
    from orpheus.transport.fields.angular_flux import AngularFlux
    from orpheus.transport.fields.boundary_flux import BoundaryFlux
    from orpheus.transport.timed_full_field import TimedFullField
    from orpheus.numerics.space import FunctionSpace
    from .geometry import SNMesh
    from orpheus.numerics.projection import MomentProjection
    from orpheus.transport.source_sinks import ScalarSourceSink, AngularSourceSink
    from .spatial.pole_angular_closure import PoleAngularClosure
    from .loss_representation import LossRepresentation

__all__ = [
    "StreamingOperator",
    "CollisionOperator",
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
    method: str, sn_mesh: "SNMesh", field: "TimedFullField",
) -> None:
    r"""The shared matvec input contract — typed composite + mesh identity.

    The two guards :meth:`StreamingOperator.apply` introduced (D-I.3d typed
    contract + the mesh-identity invariant) are now consumed by EVERY SN
    matvec entry that takes a :class:`TimedFullField`:
    :meth:`StreamingOperator.apply` / :meth:`apply_transpose` AND the #240
    Phase 2 Step B :meth:`InvertibleOperator.apply` / :meth:`apply_transpose`
    overrides.  Single source of the contract (``coding-elegance`` Pattern 2 /
    Pattern 4 — illegal inputs unrepresentable at one place, not re-validated
    per leaf).

    Parameters
    ----------
    method : str
        Qualified method name for the error message (e.g.
        ``"StreamingOperator.apply"``).
    sn_mesh : SNMesh
        The operator's mesh — ``field.bulk.mesh`` must be the SAME instance.
    field : TimedFullField
        The matvec input (``psi`` for apply, ``phi`` for the transpose).
    """
    from orpheus.transport.timed_full_field import TimedFullField

    if not isinstance(field, TimedFullField):
        raise TypeError(
            f"{method}: expected TimedFullField, got "
            f"{type(field).__name__}.  D-I.3d (2026-05-29) retired the "
            "bare-ndarray packed-vector contract; construct a typed "
            "composite via ``TimedFullField.zeros(bulk=AngularFlux, "
            "boundary=BoundaryFlux, mesh=sn_mesh)`` or explicit "
            "``TimedFullField(bulk=..., boundary=...)``."
        )
    if sn_mesh is not field.bulk.mesh:
        raise ValueError(
            f"{method}: operator and composite must share the SAME "
            "SNMesh instance (mesh-identity invariant)."
        )


@dataclass
class StreamingOperator(LinearOperatorMixin["TimedFullField"]):
    r"""Pure streaming + angular-redistribution operator :math:`L` as a
    :class:`~orpheus.numerics.operator.LinearOperator` leaf.

    The "L" of the Phase G four-operator algebra
    :math:`A_{\rm wg} = L + C - S_{\rm foldable}`. Carries the spatial
    streaming math plus the curvilinear angular redistribution — the
    cell-collision term lives in the separate :class:`CollisionOperator`
    leaf. The split lets the within-group operator be pure algebra
    (``L + C - S.foldable_part()``); no ``WithinGroupOperator`` wrapper.

    Resolution A — subtractive ``apply``
    ------------------------------------

    .. math::

        L.{\rm apply}(\psi) \;:=\; M(\psi;\;\sigma_t) \;-\;
                                  \sigma_t \odot \psi_{\rm packed}

    The matvec primitive is called at the user's σ_t (constructor-stored),
    then the cell-collision term :math:`\sigma_t \odot \psi` is
    subtracted at the packed-vector level. Combined with
    :math:`C.{\rm apply}(\psi) = \sigma_t \odot \psi`, this gives
    :math:`(L + C).{\rm apply}(\psi) = M(\psi;\;\sigma_t)` bit-exact.

    Why σ_t is a CONSTRUCTOR parameter (Pattern 4)
    ----------------------------------------------

    The discrete curvilinear matvec is RATIONAL in σ_t through the
    Carlson coupled-pole seed (Hébert §3.9.4 Eqs. 3.432-3.435):

    .. math::

        \bar\phi_i \;=\; \frac{\Delta r_i\,\bar Q_i + 2\,\bar\phi_{i+1/2}}
                              {\Delta r_i\,\Sigma_t(i) + 2}

    where :math:`\bar Q_i = \Sigma_t(i)\,\phi_0(i) / W` at ℓ=0. The
    discrete L (with M-M angular closure) is therefore intrinsically
    σ_t-coupled — analogous to the DD coefficient :math:`\alpha_{DD}
    (\sigma_t\,\Delta x)` in characteristic-line methods, or the
    exponential characteristic :math:`\exp(-\sigma_t\,s)` in MoC/CP.
    Discrete streaming operators routinely carry σ_t through their
    closure choices.

    The CONTINUOUS L (:math:`\Omega\cdot\nabla\psi + (1-\mu^2)/r\,
    \partial\psi/\partial\mu`) is σ_t-independent. The DISCRETE L's
    σ_t parametrisation is a property of the discretisation. Pattern 4
    (illegal states unrepresentable) — ``StreamingOperator(sn_mesh)``
    without σ_t is not a meaningful object for curvilinear SN; the
    constructor refusing to construct without σ_t encodes that fact.

    Capability set
    --------------

    ``frozenset({CAP_APPLY})`` — pure streaming alone is **not
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
    advertises ``CAP_APPLY_TRANSPOSE`` and ``L.H`` is the physical G-adjoint.

    Parameters
    ----------
    sn_mesh : SNMesh
        The augmented geometry carrying quadrature, BCs (the
        face-name-keyed ``sn_mesh.bc`` dict), pole closure, and (for
        curvilinear) the precomputed connection coefficients — no
        ``boundary`` constructor parameter.
    sigma_t : np.ndarray
        Total cross-section, shape ``(ng, nx, ny)`` (Issue #196 PR-INDEX-3).
        Carried at constructor time per Resolution A's subtractive
        definition.

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
    sigma_t: np.ndarray

    capabilities: frozenset[str] = field(
        default_factory=lambda: frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})
    )
    # Streaming is the sole FULL operator — it couples bulk ↔ boundary
    # (reads the inflow trace to seed the sweep, writes the outflow
    # trace). Issue #208 / Wave O. Class-level constant (unannotated so
    # the dataclass does not treat it as a field).
    # CAP_APPLY_TRANSPOSE (Wave O / O.2b): the analytic reverse-direction
    # adjoint matvec landed — see :meth:`apply_transpose`.
    block_role = BlockRole.FULL

    # D-J (2026-05-30): ``_eq_map`` / ``_ensure_eq_map`` / ``n_unknowns``
    # retired alongside the :class:`EquationMap` codec family — the
    # typed :class:`~orpheus.transport.timed_full_field.TimedFullField`
    # contract has no need for the legacy packed-vector slot map.

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        r"""The composite carrier :math:`V_{\rm bulk}\oplus V_{\rm trace}` (Wave O / O.2b).

        :math:`L` is the sole FULL operator — it couples bulk :math:`\leftrightarrow`
        boundary (seeds the sweep from the inflow trace, emits the outflow
        trace). Advertising :attr:`~orpheus.sn.geometry.SNMesh.full_field_space`
        is what lets :class:`~orpheus.numerics.operator._AdjointOperator`
        read the **block-diagonal G-adjoint metric** (bulk :math:`V\,w_n`
        :math:`\oplus` trace :math:`|\Omega\cdot\hat n|\,w_n`) for ``L.H`` —
        without it the adjoint silently reduces to the metric-blind Euclidean
        transpose (Issue #208 risk R5).

        ``C`` / ``S`` / ``F`` report ``None`` (skipped by the composition
        guard), so ``L``'s composite domain propagates through ``OperatorSum``
        to the **transpose-closed** sub-sums whose ``.H`` is actually reachable
        — ``(L + C)`` and ``(L + C - B)`` — and every bulk leaf in those is
        G-conjugated for free by the op-level :math:`G^{-1}(\sum \text{leaf}^{\mathsf T})G`.
        The full prompt loss ``(L + C - S - F - B).H`` is **intentionally
        unreachable**: ``S`` / ``F`` advertise no ``apply_transpose``, so
        ``OperatorSum`` does not propagate ``CAP_APPLY_TRANSPOSE`` and
        :class:`~orpheus.numerics.operator._AdjointOperator` raises
        :class:`MissingCapability` (fails loud, never silently Euclidean —
        the capability lattice makes the metric-blind state unrepresentable).
        The foldable scattering / fission contributions are handled at the
        eigenvalue / DSA outer layer, not via this within-group adjoint.
        """
        return self.sn_mesh.full_field_space

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        # Endomorphism on the composite (see :meth:`domain`).
        return self.sn_mesh.full_field_space

    def apply(self, psi: "TimedFullField") -> "TimedFullField":
        r"""Subtractive forward action :math:`L\,\psi = M(\psi;\sigma_t)
        - \sigma_t \odot \psi.bulk`.

        The within-group loss action :math:`(L+C)\psi` is the
        :attr:`loss_representation`'s walk (the matvec twin of the forward
        sweep; S6.3 moved it OFF this operator) — the selection fact lives
        on that property, single-sourced.  This leaf then applies the ONLY
        algebra glue — the Resolution-A collision subtraction below.

        Resolution A — :math:`(L + C).{\rm apply}(\psi) \equiv
        M(\psi;\sigma_t)` bit-exact: the per-cell σ_t·ψ cancellation
        lives at the algebra layer; this leaf returns the subtractive
        :math:`L`-only action on the typed direct-sum carrier.

        D-I.3d (2026-05-29) collapsed the bare-ndarray packed-vector
        adapter; D-J (2026-05-30) retired the supporting codec family.
        Producer-side normalisation (Pattern 7) — the operator
        consumes / emits ONLY
        :class:`~orpheus.transport.timed_full_field.TimedFullField`;
        the typed↔packed boundary collapses at the producer site.

        ``L`` is the ONLY operator that emits a non-zero face
        residual on its output ``.boundary`` —
        :class:`CollisionOperator`,
        :class:`~orpheus.sn.scattering.ScatteringOperator`, and
        :class:`~orpheus.sn.fission.FissionOperator` all leave the
        output boundary at the auto-allocated zero.

        Parameters
        ----------
        psi : TimedFullField
            Composite carrier with bulk
            (:class:`~orpheus.transport.fields.angular_flux.AngularFlux`)
            and boundary
            (:class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`).
            Operator and ``psi.bulk.mesh`` MUST be the same
            :class:`~orpheus.sn.geometry.SNMesh` instance.

        Returns
        -------
        TimedFullField
            ``L·ψ`` as a typed composite — bulk carries the cell
            action, boundary carries the face residual at the
            layout-assigned trace slots (non-zero at outer face for
            curvilinear; non-zero at outer + inner faces for slab).
        """
        from orpheus.transport.source_sinks import AngularSourceSink
        from orpheus.transport.timed_full_field import TimedFullField

        _require_typed_composite("StreamingOperator.apply", self.sn_mesh, psi)

        # L = (L+C) − C (Resolution A): the representation returns the FULL
        # within-group loss (L+C)ψ (its walk — L21: matvec and sweep are two
        # actions of the SAME (L+C) operator; S6.3 moved the walk OFF this
        # operator INTO the representation). The operator's ONLY remaining glue
        # is the collision-diagonal subtraction C = σ_t⊙, applied ONCE here (it
        # was 5× duplicated across the retired _apply_*). C has no boundary
        # action, so the (L+C)ψ boundary residual passes through unchanged.
        lpc = self.loss_representation.loss_action(self.sigma_t, psi)
        out_bulk = lpc.bulk.values - self.sigma_t[None] * psi.bulk.values
        return TimedFullField(
            bulk=AngularSourceSink.from_mesh(out_bulk, self.sn_mesh),
            boundary=lpc.boundary,
            _history=(),
            history_depth=psi.history_depth,
        )

    def apply_transpose(self, phi: "TimedFullField") -> "TimedFullField":
        r"""Hilbert transpose :math:`L^{\mathsf T}\,\phi` (Wave O / O.2b, #208).

        Resolution A: :math:`L = (L + C) - C` with :math:`C = \sigma_t\odot`
        a self-adjoint diagonal, so :math:`L^{\mathsf T} = (L+C)^{\mathsf T} - C`.
        The :math:`(L+C)^{\mathsf T}` core is the representation's
        ``loss_action_transpose``
        (:meth:`~orpheus.sn.loss_representation._OneDimScanWalk.loss_action_transpose`,
        the reverse-mode adjoint of the forward sweep; #206 Phase C moved it off
        this operator); the :math:`-\sigma_t\odot\phi.bulk` subtraction mirrors
        :meth:`apply`'s Resolution-A factoring (``C`` has no boundary action, so
        the transpose's boundary block is :math:`(L+C)^{\mathsf T}`'s).

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
        from orpheus.transport.source_sinks import AngularSourceSink
        from orpheus.transport.timed_full_field import TimedFullField

        _require_typed_composite(
            "StreamingOperator.apply_transpose", self.sn_mesh, phi,
        )
        # Lᵀ = (L+C)ᵀ − C (Resolution A; C = σ_t⊙ is a self-adjoint diagonal).
        # The representation returns (L+C)ᵀφ (its reverse walk — CumprodScan
        # carries the curvilinear angular second triangular factor; the multi-D
        # Cartesian adjoint stays a deferral raise, never a silent wrong
        # answer). The operator subtracts C here, ONCE, mirroring :meth:`apply`.
        lpct = self.loss_representation.loss_action_transpose(self.sigma_t, phi)
        out_bulk = lpct.bulk.values - self.sigma_t[None] * phi.bulk.values
        return TimedFullField(
            bulk=AngularSourceSink.from_mesh(out_bulk, self.sn_mesh),
            boundary=lpct.boundary,
            _history=(),
            history_depth=phi.history_depth,
        )

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
        from .loss_representation import default_for

        return default_for(self.sn_mesh)

    # ── Algebra dispatch — sweep-invertible composite (R-1 Step C) ────

    def __add__(self, other):
        r"""Compose :math:`L + X`.

        When ``X`` is a :class:`CollisionOperator`, returns the
        sweep-invertible specialisation :class:`InvertibleOperator`
        carrying the algebraic identity :math:`(L + C)^{-1} \approx
        \text{WDD sweep}`.  Otherwise falls through to the generic
        :class:`OperatorSum` via the mixin.
        """
        if isinstance(other, CollisionOperator):
            return InvertibleOperator(self, other)
        return super().__add__(other)


@dataclass
class CollisionOperator(LinearOperatorMixin["TimedFullField"]):
    r"""Pure collision operator :math:`C = \sigma\cdot` as a
    :class:`~orpheus.numerics.operator.LinearOperator` leaf.

    The "C" of the Phase G four-operator algebra
    :math:`A_{\rm wg} = L + C - S_{\rm foldable}`. Diagonal in position,
    group, and direction:

    .. math::

        (C\psi)_{n,i,j,g} \;=\; \sigma(i,j,g)\,\psi_{n,i,j,g}

    The supplied ``sigma`` is per-cell per-group only — the operator
    broadcasts identically across every ordinate ``n`` because the
    collision cross-section is direction-independent.

    Convention-agnostic σ
    ---------------------

    The same operator class supports the full total cross-section
    :math:`\sigma_t` AND the within-group removal cross-section
    :math:`\sigma_r = \sigma_t - \Sigma_{s,0}^{g\to g}`. The operator
    does NOT carry an interpretation flag — both quantities are
    ``(ng, nx, ny)`` arrays applied as :math:`\sigma\cdot\psi` (see
    :ref:`theory-sn-index-convention` for the canonical layout). The
    fusion hook (substep 3+4.b.ii) builds a lazy :math:`\sigma_r`
    :class:`CollisionOperator` from
    :meth:`~orpheus.sn.scattering.ScatteringOperator.foldable_sigma` at
    the moment it detects ``L + C - S.foldable_part()`` in an
    :class:`~orpheus.numerics.operator.OperatorSum`. (Pattern 4 —
    illegal states unrepresentable via the type, not a runtime flag
    the operator never inspects.)

    ⚠ A :math:`\sigma_r`-sweep is NOT a correct within-group self-scatter
    accelerator for anisotropic flux — it inverts the *diagonal-in-angle*
    removal :math:`\sigma_r\cdot I`, whereas the foldable self-scatter is the
    *isotropic-projection* operator :math:`\Sigma_{s,0}\,P_{\rm iso}`; the two
    agree only for isotropic ψ. Wiring ``A_wg.solve(S_residual)`` with this
    sweep ships 46–56 % errors on vacuum / heterogeneous problems (issue
    #215). The stable+correct fold is consistent DSA (#2) or Krylov. See the
    foldable/residual split note in :mod:`orpheus.sn.scattering`.

    Capability set
    --------------

    ``frozenset({CAP_APPLY, CAP_SOLVE, CAP_APPLY_TRANSPOSE})`` —
    collision is the simplest sort of operator. ``solve(q) = q / σ``
    is element-wise division; ``apply_transpose == apply`` because the
    operator is self-adjoint (diagonal in every basis the SN method
    operates in). All three capabilities are analytic.

    Parameters
    ----------
    sn_mesh : SNMesh
        The augmented geometry — carries mesh + quadrature + boundary
        operators (same as :class:`StreamingOperator`).
    sigma : np.ndarray
        Per-cell per-group cross-section, shape ``(ng, nx, ny)``
        (Issue #196 PR-INDEX-3). May be σ_t (full collision) or σ_r
        (removal — within-group self-scatter folded). The operator's
        action is identical.
    """

    sn_mesh: "SNMesh"
    sigma: np.ndarray

    capabilities: frozenset[str] = field(
        default_factory=lambda: frozenset(
            {CAP_APPLY, CAP_SOLVE, CAP_APPLY_TRANSPOSE}
        )
    )
    # Collision is a BULK operator — diagonal in (cell, group, ordinate),
    # no boundary action (A_bb only). Issue #208 / Wave O. Class-level
    # constant (unannotated so the dataclass does not treat it as a field).
    block_role = BlockRole.BULK

    # D-I.1 (2026-05-29): lazy ``_eq_map`` cache + ``_ensure_eq_map`` +
    # ``_sigma_at_unknowns`` retired together with the bare-ndarray
    # packed-vector apply / solve arms.  All consumers route through
    # :class:`TimedFullField`; no packed-vector decoder needed.

    def apply(self, psi: "TimedFullField") -> "TimedFullField":
        r"""Forward action :math:`C\,\psi = \sigma\cdot\psi` on the composite carrier.

        The diagonal action :math:`\sigma\cdot\psi` is per-cell per-group,
        broadcast across every ordinate.  Bulk receives
        ``σ ⊙ ψ.bulk.values``; boundary is the implicit-zero
        :class:`BoundaryFlux` — collision has no face-trace contribution
        (the cell-balance :math:`\sigma\cdot\psi` term is a CELL quantity;
        the boundary residual is a TRACE equation).  Option β3 / Issue #208
        will encode this bulk-only nature in the type via
        :class:`BulkOperator`.

        D-I.1 (2026-05-29) retired the legacy bare-ndarray packed-vector
        arm.  :class:`TimedFullField` is the sole accepted carrier.
        """
        from orpheus.transport.timed_full_field import TimedFullField
        from orpheus.transport.fields.angular_flux import AngularFlux
        from orpheus.transport.source_sinks import AngularSourceSink, BoundarySourceSink
        mesh = psi.bulk.mesh
        return TimedFullField(
            bulk=AngularSourceSink.from_mesh(
                self.sigma[None] * psi.bulk.values, mesh,
            ),
            boundary=BoundarySourceSink.zeros_on(mesh),
            _history=(),
            history_depth=psi.history_depth,
        )

    def solve(self, q: "TimedFullField") -> "TimedFullField":
        r"""Inverse action :math:`C^{-1}\,q = q/\sigma` on the composite carrier.

        Trivially invertible on the bulk: collision is diagonal, so the
        inverse is per-slot reciprocal scaling.  Returns NaN / Inf at
        slots where ``σ == 0`` per the IEEE-754 division contract —
        consumers constructing :math:`\sigma_r = \sigma_t -
        \Sigma_{s,0}^{g\to g}` must guarantee positivity (the operator
        does not check).

        Boundary is the implicit-zero :class:`BoundaryFlux` (Option β3,
        formal pseudoinverse on the rank-deficient face block — face
        slots of ``q`` are NOT inverted because collision contributes
        no volumetric term on the trace).

        D-I.1 retired the legacy bare-ndarray packed-vector arm.
        :class:`TimedFullField` is the sole accepted carrier.
        """
        from orpheus.transport.timed_full_field import TimedFullField
        from orpheus.transport.fields.angular_flux import AngularFlux
        from orpheus.transport.fields.boundary_flux import BoundaryFlux
        mesh = q.bulk.mesh
        return TimedFullField(
            bulk=AngularFlux.from_mesh(
                q.bulk.values / self.sigma[None], mesh,
            ),
            boundary=BoundaryFlux.zeros_on(mesh),
            _history=(),
            history_depth=q.history_depth,
        )

    def apply_transpose(self, psi: "TimedFullField") -> "TimedFullField":
        r"""Adjoint action :math:`C^*\,\psi = \sigma\cdot\psi`.

        Equal to :meth:`apply` — collision is self-adjoint (diagonal
        operator). Returned bit-equal to ``apply(psi)``.
        """
        return self.apply(psi)

    # ── Algebra dispatch — sweep-invertible composite (R-1 Step C) ────

    def __add__(self, other):
        r"""Compose :math:`C + X`.

        When ``X`` is a :class:`StreamingOperator`, returns the
        sweep-invertible specialisation :class:`InvertibleOperator`
        with the streaming operator placed first (the canonical
        ``L + C`` ordering for the algebraic identity).  Otherwise
        falls through to the generic :class:`OperatorSum` via the
        mixin.
        """
        if isinstance(other, StreamingOperator):
            return InvertibleOperator(other, self)
        return super().__add__(other)


# ─────────────────────────────────────────────────────────────────────────
# InvertibleOperator — sweep-invertible composite (L + C)
# ─────────────────────────────────────────────────────────────────────────


class InvertibleOperator(OperatorSum["TimedFullField"]):
    r"""Sweep-invertible composite :math:`L + C` carrying ``.solve`` = WDD sweep.

    R-1 Step C (2026-05-19) — the SN-specific algebraic identity

    .. math::

        (L_{\rm streaming} + C_{\rm diagonal})^{-1} \;\approx\;
        \text{WDD sweep}

    has no generic ``(A+B)^{-1}`` formula — :class:`OperatorSum` by
    itself cannot ``solve``.  The WDD sweep IS the inverse algorithm
    for this specific composite — that's the algebraic foundation of
    the entire SN method (Lewis & Miller §3.2; Adams & Larsen 2002
    §III).  :class:`InvertibleOperator` is the specialisation that
    carries the identity at the type level: it OWNS its full action
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
      :class:`StreamingOperator` and ``C`` is a
      :class:`CollisionOperator` returns this class automatically (see
      :meth:`StreamingOperator.__add__` and
      :meth:`CollisionOperator.__add__`).  The composite reads as math.
    * **Explicit construction** — ``InvertibleOperator(L, C)``.  Useful
      when composing variants such as
      ``InvertibleOperator(L_leaf, CollisionOperator(σ_r))`` where
      ``σ_r = σ_t - Σ_{s,0}^{g→g}`` is the removal cross-section that
      lets one fold the within-group self-scatter into the diagonal
      collision term (Adams & Larsen 2002 §III; tracked by issue
      `#200 <https://github.com/deOliveira-R/ORPHEUS/issues/200>`_).

    The two paths produce structurally identical objects — the choice
    only changes the call-site readability.

    Capability set
    ==============

    ``frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE, CAP_SOLVE})`` — adds
    ``solve`` (the WDD sweep) to the parent :class:`OperatorSum`'s set;
    ``apply_transpose`` propagates by the :class:`OperatorSum` closure
    law (both :math:`L` and :math:`C` advertise it) and is OVERRIDDEN to
    the composite's own :math:`M(\sigma)^{\mathsf T}` action (Wave O #208 /
    #240 Step B).  The multi-D Cartesian adjoint raises (the
    representation's deferral contract — never a silent wrong answer).

    The ``.solve`` API
    ==================

    The ``rhs`` parameter is a typed
    :class:`~orpheus.sn.angular_flux.AngularFlux` carrying:

    * ``rhs.values`` — per-ordinate source ``(N, ng, nx, ny)``.  This
      is treated as the per-ordinate anisotropic source
      :math:`Q^{\rm aniso}` that the sweep consumes (the isotropic
      source is zero).
    * ``rhs.boundary`` — face source / BC inflow trace.  Typically
      zero for volumetric SI/Krylov sources (which carry no face
      contribution); the persistent reflective-BC state lives on the
      :class:`SNMesh` and is handled inside the sweep.
    * ``rhs(1)`` (the lag-1 frame, when ``rhs.history_depth >= 2``) —
      the previous iterate that GENERATED this source.  Threaded as
      :pydata:`initial_guess` to ``loss_representation.sweep`` so the
      curvilinear sweep can read the Carlson coupled-pole seed from
      it.  ``None`` (cold start) → the sweep falls back to its
      in-iteration-source default.

    Parameters
    ----------
    streaming : StreamingOperator
        :math:`L = \Omega\cdot\nabla + \text{angular redistribution}`.
        Resolution A subtractive form:
        ``L.apply(ψ) = M(ψ; σ_t) - σ_t ⊙ ψ_cell``.
    diagonal : CollisionOperator
        :math:`C = \sigma\cdot`.  Its ``.sigma`` attribute is the
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
        diagonal: "CollisionOperator",
    ) -> None:
        if not isinstance(streaming, StreamingOperator):
            raise TypeError(
                f"InvertibleOperator: 'streaming' must be a "
                f"StreamingOperator; got {type(streaming).__name__}."
            )
        if not isinstance(diagonal, CollisionOperator):
            raise TypeError(
                f"InvertibleOperator: 'diagonal' must be a "
                f"CollisionOperator; got {type(diagonal).__name__}."
            )
        if streaming.sn_mesh is not diagonal.sn_mesh:
            raise ValueError(
                "InvertibleOperator: streaming and diagonal operators "
                "must share the same SNMesh instance "
                "(mesh-identity invariant)."
            )
        if not np.all(diagonal.sigma > 0):
            min_sigma = float(np.min(diagonal.sigma))
            raise ValueError(
                f"InvertibleOperator: diagonal coefficient must be "
                f"strictly positive everywhere for the WDD sweep to be "
                f"well-defined; got min(sigma) = {min_sigma:.3e}.  If "
                f"sigma_r = sigma_t - Sigma_(s,0)^(g->g) is dipping "
                f"negative, the multi-group cross-section set is "
                f"physically inconsistent."
            )
        super().__init__(streaming, diagonal)
        # OperatorSum.__init__ set capabilities = {CAP_APPLY, ...};
        # we add CAP_SOLVE because this composite IS sweep-invertible.
        self.capabilities = self.capabilities | frozenset({CAP_SOLVE})
        # block_role is now DERIVED by OperatorSum.__init__ (Wave O / O.2b 4.5):
        # join(L=FULL, C=BULK) = FULL. The former hand-stamp here was the
        # twin-path retired in 4.5 — the role is carried by construction.

    # ── Convenience accessors ─────────────────────────────────────────

    @property
    def streaming(self) -> "StreamingOperator":
        """The streaming operand (alias for ``self.a``)."""
        return self.a  # type: ignore[return-value]

    @property
    def diagonal(self) -> "CollisionOperator":
        """The diagonal-collision operand (alias for ``self.b``)."""
        return self.b  # type: ignore[return-value]

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
        r"""The diagonal coefficient used by ``solve`` (σ_t or σ_r)."""
        return self.diagonal.sigma

    # ── apply / apply_transpose: the composite's OWN matvec (#240 Step B) ──

    def apply(self, psi: "TimedFullField") -> "TimedFullField":
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

    def apply_transpose(self, phi: "TimedFullField") -> "TimedFullField":
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

    # ── solve: WDD sweep ─────────────────────────────────────────────

    def solve(
        self,
        rhs: "AngularFlux | TimedFullField",
        *,
        initial_guess: "AngularFlux | TimedFullField | None" = None,
    ) -> "AngularFlux | TimedFullField":
        r"""Invert :math:`(L + C)\,\psi = \text{rhs}` via the WDD sweep.

        The cell-balance equation
        :math:`(\Omega\cdot\nabla + \sigma)\,\psi = Q` is integrated
        cell-by-cell in inflow-to-outflow order; the angular closure
        (Cartesian → identity, curvilinear → Morel-Montry) is bound
        on the mesh.

        Parameters
        ----------
        rhs : AngularFlux or TimedFullField
            Legacy :class:`AngularFlux` or composite
            :class:`TimedFullField` carrying:

            * ``values`` / ``bulk.values`` — per-ordinate source
              :math:`Q^{\rm aniso}`, shape ``(N, ng, nx, ny)``.
            * ``boundary`` / ``bulk.boundary`` — BC inflow trace
              (typically zero for SI/Krylov volumetric sources;
              falls back as the seed when no ``initial_guess`` is
              supplied).

            Both branches dispatch via runtime isinstance; the
            composite branch (D-H.1c stage 1) bridges through the
            legacy sweep path internally and returns
            :class:`TimedFullField`.  Output type matches input type.
        initial_guess : AngularFlux, TimedFullField, or None, optional
            Previous iterate :math:`\psi^{(k-1)}` for the curvilinear
            Carlson coupled-pole seed (M-M reads it as the level's
            angular flux to derive :math:`\bar Q` at :math:`\mu = -1`).
            Its ``.boundary`` also seeds the reflective-BC partner-flux
            trace when present (the previous outflow IS the partner-
            ordinate inflow).  ``None`` (default) → cold start: M-M's
            Carlson seed degenerates to the zero-input result (same
            as ``ZeroSeed`` ablation); BC trace falls back to
            ``rhs.boundary``.

            Type must match ``rhs`` (both legacy or both composite);
            mixed inputs raise :class:`TypeError`.

            Explicit kwarg (post-Phase-1.2) — the seed is no longer
            piggy-backed on ``rhs(1)`` (the AngularFlux lag-1 history
            is reserved for time-derivative tracking, an unrelated
            concern).  Outer iterators (:class:`SourceIteration`,
            :class:`KEigenvalue`) pass the previous iterate
            explicitly; GMRES residual calls pass ``None`` so the
            preconditioner doesn't silently route through stateful
            seed state (closes the R-1 Step D silent-fallback bug
            class).

        Returns
        -------
        AngularFlux or TimedFullField
            The angular flux satisfying :math:`(L + C)\,\psi =
            \text{rhs}`, with the sweep's outflow face state in
            ``.boundary`` and ``history_depth`` inherited from
            ``rhs.history_depth``.  Return type matches ``rhs`` input
            type.
        """
        return self._solve_timed_full_field(
            rhs, initial_guess=initial_guess,
        )

    def solve_moments(
        self,
        rhs: "TimedFullField",
        projection: "MomentProjection",
        *,
        initial_guess: "TimedFullField | None" = None,
    ) -> "TimedFullField":
        r"""Invert :math:`(L + C)\,\psi = \text{rhs}` and return the harmonic
        MOMENTS of :math:`\psi`, projected IN-SWEEP per anti-diagonal — the full
        per-ordinate angular field is never materialized.

        The Phase 5c moment-emitting sibling of :meth:`solve`: the SAME WDD
        sweep + boundary handling, but the bulk of the returned
        :class:`TimedFullField` is a
        :class:`~orpheus.transport.fields.harmonic_moment_field.HarmonicMomentField`
        ``(L+1, 2L+1, ng, nx, ny)`` rather than an
        :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
        ``(N, ng, nx, ny)``.  ``projection`` is the scattering operator's
        :class:`~orpheus.numerics.projection.MomentProjection` (its harmonics +
        weights), so the in-sweep moments equal ``S``'s internal projection
        term-for-term; the cross-octant accumulation reorders the ordinate sum
        vs the flat post-sweep projection ⇒ principled-equivalence, NOT
        bit-identity.  2-D Cartesian ONLY (the windowed-SI path; the windowing
        gate — the genuine ``is_cartesian and ndim == 2`` condition in
        ``_maybe_window`` since C5.4 (#225), replacing the ``reduced is
        None`` proxy that was ALSO true at d=3 Cartesian — guarantees it).
        ``solve`` followed by the flat :meth:`MomentProjection.apply` is the
        fuller-view verification oracle (``vv-principles``; the
        aggressive-retirement "verification oracle" exception).
        """
        return self._solve_timed_full_field(
            rhs, initial_guess=initial_guess, moment_projection=projection,
        )

    def _solve_timed_full_field(
        self,
        rhs: "TimedFullField",
        *,
        initial_guess: "TimedFullField | None" = None,
        moment_projection: "MomentProjection | None" = None,
    ) -> "TimedFullField":
        r"""Composite :class:`TimedFullField` body of :meth:`solve` (D-H.1c stage 1).

        Runs the WDD forward substitution on
        :attr:`loss_representation` — the operator's ONE
        :class:`~orpheus.sn.loss_representation.LossRepresentation`
        instance (S6.5, #222) — and handles the L2 field plumbing at
        the public-entry boundary.

        The boundary plumbing for the reflective-BC partner-flux state
        is preserved: ``initial_guess.boundary`` (composite) → legacy
        BoundaryFlux via :meth:`TimedFullField.to_legacy_angular_flux`
        → ``boundary_buf`` (the sweep's mutable write-through buffer)
        via :func:`_copy_boundary_face_state`.  The sweep mutates
        ``boundary_buf`` in place; the result is re-wrapped as an L2
        composite at the end.

        Parameters
        ----------
        rhs : TimedFullField
            Per-ordinate source on the composite carrier.
        initial_guess : TimedFullField or None
            Previous iterate (carries the partner-flux trace and
            curvilinear Carlson seed).  ``None`` → cold start.

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
        from orpheus.transport.fields.harmonic_moment_field import (
            HarmonicMomentField,
        )
        from orpheus.transport.timed_full_field import TimedFullField

        # D-H.2-C3: only the :class:`TimedFullField` composite branch remains;
        # legacy :class:`AngularFlux` retired.  ``rhs`` and ``initial_guess``
        # MUST be :class:`TimedFullField` (or ``None`` for ``initial_guess``).
        # Single guard site for both :meth:`solve` and :meth:`solve_moments`.
        if not isinstance(rhs, TimedFullField):
            raise TypeError(
                f"InvertibleOperator: 'rhs' must be TimedFullField; "
                f"got {type(rhs).__name__}.  Legacy AngularFlux retired "
                f"in D-H.2-C3."
            )
        if initial_guess is not None and not isinstance(
            initial_guess, TimedFullField,
        ):
            raise TypeError(
                f"InvertibleOperator: 'initial_guess' must be "
                f"TimedFullField or None; got "
                f"{type(initial_guess).__name__}."
            )

        sn_mesh = self.sn_mesh
        if rhs.bulk.mesh is not sn_mesh:
            raise ValueError(
                "InvertibleOperator.solve(TimedFullField): rhs and "
                "operator must share the same SNMesh instance "
                "(mesh-identity invariant)."
            )
        if initial_guess is not None and initial_guess.bulk.mesh is not sn_mesh:
            raise ValueError(
                "InvertibleOperator.solve(TimedFullField): initial_guess "
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
        # projection rides as an optional kwarg (the 1-D representation
        # raises on it, since moment output is 2-D Cartesian only;
        # ``moment_*`` and ``initial_guess`` are mutually 2-D-vs-1-D, so
        # the 2-D branch harmlessly drops the unused seed).  Only the
        # OUTPUT WRAP differs: full angular field vs harmonic moments.
        bulk_values, _scalar = self.loss_representation.sweep(
            rhs.bulk.values,
            self.sigma,
            boundary_buf,
            initial_guess=initial_guess,
            moment_projection=moment_projection,
        )
        # The sweep output carries the trailing 2^d spatial-moment axis at a
        # multi-moment closure (the φ̂ iterate, #240 D5b-S3); the typed wrap
        # selects the SpatialMomentSpace factor so the iterate is a legal typed
        # state.  DD/Step (per_axis == 1) → no factor, byte-identical.
        per_axis = sn_mesh.scheme.spatial_basis_per_axis
        if moment_projection is None:
            bulk = AngularFlux.from_mesh(
                bulk_values, sn_mesh, spatial_moments=per_axis,
            )
        else:
            bulk = HarmonicMomentField.from_mesh_and_L(
                bulk_values, sn_mesh, moment_projection.L,
                spatial_moments=per_axis,
            )

        # ── L2 direct return — no adapter needed (D-H.2-C2). ───────────
        return TimedFullField(
            bulk=bulk,
            boundary=boundary_buf,
            _history=(),
            history_depth=rhs.history_depth,
        )
