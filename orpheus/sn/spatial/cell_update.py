r"""Per-cell update strategy contract for the SN sweep.

Why this abstraction exists
===========================

The SN sweep walks the spatial mesh once per ordinate, and for each
cell it must solve the same shape of algebraic problem: combine the
local volumetric source, the upstream face flux (and, in curvilinear
geometry, the upstream angular half-flux) with the local total cross
section and a *spatial closure* — Diamond Difference (DD), Linear
Discontinuous (LD), Exponential Characteristic (EC), Step — to
produce the cell-average flux plus the downstream face flux and (for
sphere/cylinder) the downstream angular half-flux.

Per Cardinal Rule 2 (architecture), the cell-update math is **the
same algebra** in slab, sphere, and cylindrical 1-D — only the
populated fields of :class:`~orpheus.geometry.reduced_operator.StreamingTerms`
change.  A monolithic sweep that inlines the closure is the
historical pattern (see ``orpheus.sn.sweep``); lifting the closure to
a :class:`CellUpdate` strategy makes it possible to:

* unit-test the cell-update math in isolation (without spinning up a
  full SNMesh + iteration loop),
* swap closures (DD → LD → EC → Step) without rewriting the sweep
  driver,
* keep the sweep orchestration thin and the closure-specific algebra
  local to a strategy module.

This is Round 1 of Wave C of the SN reshape campaign — see Issues
#157 (this contract) and #158 (the bit-identical
:class:`DiamondDifference` extraction).  Wave D Issue 12 (#159)
rewrites the sweep itself around ``cell_update.update(...)`` so the
sweep dispatches to whichever strategy was selected; until then, the
existing sweep still inlines the math and the strategies live as a
parallel, testable abstraction.

What each dataclass holds
=========================

:class:`UpstreamState` — input to a single cell update::

    spatial_upstream:  (ng,) face flux entering the cell
                       (face area :math:`A_{i-1/2}` for
                       sphere/cylinder; the upstream face for slab).
    angular_upstream:  (ng,) :math:`\psi_{n-1/2,\,i}` half-flux at the
                       upstream half-angle, for spherical /
                       cylindrical geometry only.  Slab / Cartesian
                       set this to ``None`` — there is no angular
                       redistribution.

:class:`CellResult` — output of a single cell update::

    cell_average_flux:        (ng,) :math:`\overline{\psi}_i`,
                              the per-group average flux on the
                              cell.
    outgoing_spatial_flux:    (ng,) :math:`\psi_{i+1/2}`, the
                              downstream face flux via the closure
                              relation (e.g. DD's
                              :math:`\psi_{i+1/2} = 2\overline{\psi}_i
                              - \psi_{i-1/2}`).  ``None`` for the
                              cylindrical pure-azimuthal degenerate
                              case (``streaming_terms.abs_mu <
                              1e-15``) where no spatial face flow
                              exists in the cell.
    outgoing_angular_state:   (ng,) :math:`\psi_{n+1/2,\,i}`, the
                              downstream half-angle flux via the
                              Morel--Montry closure
                              :math:`\psi_{n+1/2} =
                              (\overline{\psi} - (1-\tau_{mm})\,
                              \psi_{n-1/2})/\tau_{mm}`.  ``None``
                              for slab — slab geometry has no
                              angular redistribution.

Geometry-as-data — single cell-balance algebra (Step 2.5)
========================================================

Issue #196 Phase G Step 2.5 retired the historical
``alpha_in is None`` / ``abs_mu < 1e-15`` runtime branches in
favour of geometry-blind data populated by the
:class:`StreamingTerms` factories:

* **Slab / Cartesian** — ``StreamingTerms`` carries neutral
  curvature: ``face_area_inner = face_area_outer = 1.0``,
  ``delta_A_over_w = 0.0``, ``alpha_in = alpha_out = 0.0``,
  ``tau_mm = 1.0`` (synthetic neutral element of the M-M
  closure — the closure is the identity for slab).  The
  ``CellVisit`` carries ``face_area_downstream = 1.0``.
  ``upstream_state.angular_upstream`` is ``None`` (slab has no
  half-angles); the strategy returns
  ``CellResult(outgoing_angular_state=None, ...)``.
* **Sphere / cylinder** — curvature fields physically populated.
  ``CellVisit.face_area_downstream`` is the sweep-direction-
  resolved outgoing face area (outer for outward, inner for inward).
  ``upstream_state.angular_upstream`` carries
  :math:`\psi_{n-1/2,\,i}`.  The strategy returns the M-M-closed
  ``outgoing_angular_state``.
* **Cylindrical pure-azimuthal degenerate** — the cell has no
  radial face flow on this ordinate; the iterator emits
  ``CellVisit.face_area_downstream = 0.0`` (geometric truth).
  The strategy's spatial closure outputs ``None`` when
  ``face_area_downstream == 0.0`` (no downstream face exists),
  while the M-M angular closure remains active.

The "is this output exists" checks remaining inside the strategy
(``face_area_downstream > 0.0`` for the spatial closure;
``angular_upstream is not None`` for the angular closure) are NOT
geometry dispatch — they test the structural presence of a
direction, not the geometry kind.

Where downstream consumers will call this
=========================================

In Wave D, the sweep at ``orpheus.sn.loss_representation`` (the dissolved ``sweep.py``) dispatches via
:meth:`SNMesh.dag_walk` — the per-visit packets pre-resolve
the sweep direction so the strategy sees no sign-of-:math:`\mu`
branching::

    for visit in sn_mesh.dag_walk(ordinate_idx=n,
                                  mu_level_idx=p):
        upstream = UpstreamState(
            spatial_upstream=psi_face,
            angular_upstream=psi_angle[visit.cell_idx]
                if reduced_op.requires_upstream_angular_state
                else None,
        )
        result = cell_update.update(visit, total_xs, source, upstream)
        psi_avg[visit.cell_idx] = result.cell_average_flux
        if result.outgoing_spatial_flux is not None:
            psi_face = result.outgoing_spatial_flux
        if result.outgoing_angular_state is not None:
            psi_angle[visit.cell_idx] = result.outgoing_angular_state

Round 1 shipped the contract only; Round 2 shipped
``DiamondDifference`` bit-identically; the SN-sweep-DAG refactor
adds :class:`CellVisit` and migrates the strategy contract to
consume it.

References
==========

* Lewis, E. E., & Miller, W. F. (1984). *Computational Methods of
  Neutron Transport*.  §5.3 — cell-update closures (Diamond
  Difference, weighted-DD, Linear Discontinuous); §4.5 — the
  Morel--Montry angular closure used for
  :math:`\psi_{n+1/2,\,i}`.
* Bailey, T. S., Adams, M. L., Yang, B., & Zika, M. R. (2009).
  *A piecewise linear finite element discretization of the
  diffusion equation for arbitrary polyhedral grids*.
  JCP 227, 3738--3757.  Eq. 50 (dome recursion) and Eq. 74
  (Morel--Montry) feed the curvilinear cell update.
* See also :doc:`/theory/discrete_ordinates`, "Cell update
  strategies (the strategy contract)".

See also
========

* :class:`~orpheus.geometry.reduced_operator.StreamingTerms` —
  the per-(cell, direction) packet a strategy receives.
* :class:`~orpheus.geometry.reduced_operator.ReducedStreamingOperator` —
  builds the streaming terms; its ``streaming_terms()`` method is
  the canonical extraction site.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import ClassVar, Protocol, runtime_checkable

import numpy as np

from orpheus.geometry.reduced_operator import StreamingTerms
from orpheus.numerics.registry import RegistryMixin


# ═══════════════════════════════════════════════════════════════════════
# CellVisit — one visit to one cell during an SN sweep
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True, slots=True)
class CellVisit:
    r"""One visit to one cell during an SN sweep.

    The SN sweep is a topological sort of the **directed cell graph**
    for a given ordinate, where edges are oriented by
    :math:`\mathrm{sign}(\Omega \cdot \hat n_{\text{face}})`.  This
    dataclass is the per-visit packet a
    :class:`CellUpdate` strategy receives — all sweep-direction data
    is **pre-resolved** here.  The strategy sees only the
    upstream/downstream view of this cell, never the geometric
    inner/outer view.

    SN-specific by design.  Produced by
    :meth:`orpheus.sn.geometry.SNMesh.dag_walk`.  MoC will
    define its own analog (per-ray traversal) — different DAG shape,
    different mathematical structure (fiber bundles + solution
    sheaves rather than a topological sort over a cell graph).
    Premature abstraction across SN/MoC is avoided per Cardinal Rule
    2 — there is no shared ``SweepGraph`` Protocol because there is
    no shared structure.

    Attributes
    ----------
    cell_idx : int
        Spatial cell index in the SNMesh.  The cell update reads
        ``total_xs[cell_idx]`` and ``source[cell_idx]`` at the call
        site (the strategy itself does not see ``cell_idx``).
    streaming_terms : StreamingTerms
        Pure geometric primitive from
        :class:`~orpheus.geometry.reduced_operator.ReducedStreamingOperator`.
        Carries cell volume, face areas (inner / outer — geometric
        labels), connection coefficients (:math:`\alpha`,
        :math:`\Delta A / w`, :math:`\tau_{mm}`), and signed +
        absolute primary direction cosine.
    face_area_downstream : float
        **Sweep-direction-resolved**: which of the two cell faces
        (inner or outer) is the downstream face for this visit.
        Issue #196 Phase G Step 2.5 made this float-always (was
        ``float | None``) so the unified cell-balance helper can
        consume one number regardless of geometry.

        * For slab / Cartesian: ``1.0`` — neutral value; the slab
          denominator's ``2|\mu|`` term reads ``2 * abs_mu *
          face_area_downstream = 2|μ|·1``.
        * For curvilinear with :math:`\mu \ge 0` (outward sweep,
          centre → boundary): equals
          ``streaming_terms.face_area_outer``.
        * For curvilinear with :math:`\mu < 0` (inward sweep,
          boundary → centre): equals
          ``streaming_terms.face_area_inner``.
        * For the cylindrical pure-azimuthal degenerate case
          (``abs_mu < 1e-15``): ``0.0`` — the cell has no spatial
          face flow, so the ``2|\mu| A_{\rm down}`` term vanishes
          via ``A_down = 0`` (geometric truth) rather than via the
          numerical threshold ``abs_mu < 1e-15``.

    Notes
    -----
    The companion :attr:`UpstreamState.spatial_upstream` is **also**
    sweep-direction-resolved by the orchestrator (it is the flux
    flowing INTO the cell from the previously-visited cell along
    the topological ordering).  Together,
    :attr:`face_area_downstream` and
    :attr:`UpstreamState.spatial_upstream` give the cell-update
    strategy a fully resolved view of "what flows into me, what
    flows out of me" — no sign-of-:math:`\mu` branching inside the
    strategy.

    Frozen + slotted: instances are immutable and lightweight.
    """

    cell_idx: int
    streaming_terms: StreamingTerms
    face_area_downstream: float = 0.0


# ═══════════════════════════════════════════════════════════════════════
# Per-cell input state
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True, slots=True)
class UpstreamState:
    r"""Per-cell input state that a :class:`CellUpdate` consumes.

    Attributes
    ----------
    spatial_upstream :
        Shape ``(ng,)``.  Face flux entering this cell from the
        upstream radial face — :math:`\psi_{i-1/2}` for sphere /
        cylinder; the upstream face flux for slab.  Always populated
        (the sweep always carries a face flux into the cell update).
    angular_upstream :
        Shape ``(ng,)`` for sphere / cylinder; ``None`` for slab.
        :math:`\psi_{n-1/2,\,i}` — the half-flux at the upstream
        half-angle on this cell, used by the Morel--Montry closure
        and the :math:`\Delta A/w \cdot \alpha_{n-1/2}` source
        contribution in curvilinear geometry.

    Frozen + slotted: instances are immutable and lightweight.
    """

    spatial_upstream: np.ndarray
    angular_upstream: np.ndarray | None = None


# ═══════════════════════════════════════════════════════════════════════
# Per-cell output state
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True, slots=True)
class CellResult:
    r"""Per-cell output state produced by a :class:`CellUpdate`.

    Attributes
    ----------
    cell_average_flux :
        Shape ``(ng,)``.  The cell-average flux
        :math:`\overline{\psi}_i = \mathrm{numer}/\mathrm{denom}` from
        the closure's algebraic solve.
    outgoing_spatial_flux :
        Shape ``(ng,)`` for the typical case; ``None`` for the
        cylindrical pure-azimuthal degenerate case
        (``streaming_terms.abs_mu < 1e-15``) where the cell has no
        radial face flow.  For Diamond Difference,
        :math:`\psi_{i+1/2} = 2\overline{\psi}_i - \psi_{i-1/2}`.
    outgoing_angular_state :
        Shape ``(ng,)`` for sphere / cylinder; ``None`` for slab.
        :math:`\psi_{n+1/2,\,i}` from the Morel--Montry closure
        :math:`\psi_{n+1/2} = (\overline{\psi} -
        (1-\tau_{mm})\,\psi_{n-1/2})/\tau_{mm}`.  Slab geometry has
        no angular redistribution and returns ``None``.

    Frozen + slotted: instances are immutable and lightweight.
    """

    cell_average_flux: np.ndarray
    outgoing_spatial_flux: np.ndarray | None = None
    outgoing_angular_state: np.ndarray | None = None


# ═══════════════════════════════════════════════════════════════════════
# CellUpdate Protocol
# ═══════════════════════════════════════════════════════════════════════

@runtime_checkable
class CellUpdate(Protocol):
    r"""Strategy contract for a per-cell SN sweep update.

    Concrete implementations (e.g. :class:`DiamondDifference` shipping
    in Round 2 of Wave C) are ``@dataclass`` instances exposing the
    two class-level traits below and the :meth:`update` method.  The
    contract follows the project's pattern (see
    :class:`~orpheus.numerics.operator.LinearOperator`,
    :class:`~orpheus.geometry.boundary.BoundaryTraceLaw`,
    :class:`~orpheus.geometry.reduced_operator.AngularMeasure`):
    ``@runtime_checkable Protocol``, satisfied by structural
    typing — concrete strategies do **not** need to inherit.

    Attributes
    ----------
    is_linear : bool
        Whether the closure is linear in the inputs.  Diamond
        Difference is linear; weighted-DD with a flux-dependent
        weight, Step's positivity-fixup, and exponential
        characteristic with a clipped argument are not.  Read-only
        class attribute.
    is_positivity_preserving : bool
        Whether non-negative inputs guarantee non-negative outputs.
        Diamond Difference is **not** positivity preserving (Lewis
        & Miller §5.3); Step is.  Read-only class attribute.
    is_affine_scannable : bool
        Whether the closure admits the closed-form affine recurrence
        :math:`\psi_{\rm out} = a\,\psi_{\rm in} + b` (Blelloch §1.5) so
        the DAG-free *scan* schedules (``CumprodScan``, ``ScanMarch``)
        can consume it.  ``True`` requires that the cell-average is an
        affine function of a **single** upstream face flux — Diamond
        Difference qualifies; a Linear-Discontinuous closure (two coupled
        face moments) does **not**.  ``False`` is the safe default: a
        non-affine-scannable scheme is routed to the DAG wavefront
        schedule by :meth:`LossRepresentation.default_for` (the scan
        strategies gate their ``supports`` on this trait).  Read-only
        class attribute.

    Notes
    -----
    All three traits are class-level so they can be inspected without
    instantiating the strategy.  Code that selects a closure based on
    cell-thickness or stiffness criteria reads ``is_linear``;
    diagnostics that gate on whether negative fluxes can appear read
    ``is_positivity_preserving``; the sweep-schedule selection reads
    ``is_affine_scannable``.
    """

    is_linear: bool
    is_positivity_preserving: bool
    is_affine_scannable: bool

    def update(
        self,
        visit: CellVisit,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> CellResult:
        r"""Compute the per-cell average flux and downstream states.

        Parameters
        ----------
        visit :
            One visit to this cell during the sweep.  Contains the
            pure-geometric :class:`StreamingTerms` packet plus the
            sweep-direction-resolved :attr:`face_area_downstream`.
            Geometry is data, not control-flow (Issue #196 Phase G
            Step 2.5): ``visit.face_area_downstream > 0.0`` signals
            "this cell has a downstream spatial face" (any geometry);
            ``visit.face_area_downstream == 0.0`` signals "no
            spatial face flow" (cylindrical pure-azimuthal
            degenerate).
        total_xs :
            Shape ``(ng,)``.  Per-group total cross section
            :math:`\Sigma_t` on this cell.
        source :
            Shape ``(ng,)``.  Per-group volumetric source
            :math:`Q_i` on this cell, **already weight-normalized**
            by the sweep (the sweep applies its
            :math:`1/\sum_n w_n` factor at the call site so the
            strategy receives the full per-ordinate volumetric
            source ready to plug into the numerator).
        upstream_state :
            Per-cell input state.  See :class:`UpstreamState`.
            ``upstream_state.spatial_upstream`` is the flux flowing
            **into** this cell from the previously-visited cell —
            already sweep-direction-resolved by the orchestrator.

        Returns
        -------
        CellResult
            The cell-average flux + downstream states.  Slab
            consumers ignore ``outgoing_angular_state``;
            cylindrical-degenerate consumers handle
            ``outgoing_spatial_flux is None`` by skipping the
            face-flux update.
        """
        ...

    def residual(
        self,
        cell_avg: np.ndarray,
        visit: CellVisit,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> np.ndarray:
        r"""Per-cell operator residual :math:`L_{\text{cell}}\,\bar\psi - q`.

        Companion to :meth:`update`: where ``update`` answers the
        **solve direction** ("given the source and upstream state,
        find the cell-average flux that satisfies the cell balance"),
        :meth:`residual` answers the **apply direction** ("given a
        probe value for the cell-average flux, evaluate the discrete
        per-cell operator action minus the right-hand side").

        Round-trip invariant
        --------------------

        For any ``visit``, ``total_xs``, ``source``, ``upstream_state``::

            result   = strategy.update(visit, total_xs, source,
                                       upstream_state)
            residual = strategy.residual(result.cell_average_flux,
                                         visit, total_xs, source,
                                         upstream_state)
            assert np.allclose(residual, 0.0, atol=1e-13)

        i.e. evaluating the operator residual at the cell-average
        flux returned by ``update`` yields zero to floating-point
        rounding.  This is the strategy-layer apply-vs-solve
        consistency contract.

        Use case
        --------

        The residual is the discrete per-cell operator action used by
        the matvec form ``(L + C).apply(\psi)`` — needed for explicit
        Krylov inner solvers, GMRES preconditioning, and operator-
        algebra composition above the cell-update layer.  At any
        probe point ``cell_avg`` (not necessarily the solved value),
        the residual is the per-group residual of the discretised
        cell-balance equation.

        Parameters
        ----------
        cell_avg :
            Shape ``(ng,)``.  Per-group cell-average flux probe
            point at which the operator residual is evaluated.
        visit :
            Same as :meth:`update` — pre-resolved sweep-direction
            view of the cell.
        total_xs :
            Shape ``(ng,)``.  Per-group total cross section
            :math:`\Sigma_t` on this cell.
        source :
            Shape ``(ng,)``.  Per-group volumetric source on this
            cell, **on the same already-weight-normalised convention
            that** :meth:`update` **consumes**.  The residual is
            affine in this argument (a shift ``source += δ`` shifts
            the residual by ``-δ``).
        upstream_state :
            Per-cell upstream state — same convention as
            :meth:`update`.

        Returns
        -------
        np.ndarray
            Shape ``(ng,)``.  Per-group residual
            :math:`L_{\text{cell}}\,\bar\psi - q`.  Linear in
            ``cell_avg`` (the cell-update closures shipped today are
            all linear; non-linear closures may shadow this with a
            non-linear ``residual``).
        """
        ...


# ═══════════════════════════════════════════════════════════════════════
# CellUpdateBase — concrete ABC carrying RegistryMixin
# ═══════════════════════════════════════════════════════════════════════


class CellUpdateBase(RegistryMixin, ABC):
    r"""Concrete abstract base for self-registering cell-update strategies.

    Round 1 of Wave C shipped :class:`CellUpdate` as a
    ``runtime_checkable`` Protocol so concrete strategies could be
    matched by structural typing without inheritance. Issue 9.6 adds
    :class:`CellUpdateBase` as a parallel **concrete ABC** layered on
    top of :class:`RegistryMixin`: strategies that inherit it
    self-register under their ``key=`` class-creation kwarg, gaining
    name-keyed lookup via :meth:`create`. The Protocol stays — third-
    party / synthetic strategies (e.g. test mocks) that do not want to
    inherit can still satisfy the contract by providing the right
    methods.

    Subclasses MUST declare:

    * ``is_linear: ClassVar[bool]`` — whether the closure is linear
      in the inputs (Diamond Difference is linear; weighted-DD with a
      flux-dependent weight is not).
    * ``is_positivity_preserving: ClassVar[bool]`` — whether
      non-negative inputs guarantee non-negative outputs (Step is
      positivity-preserving; DD is not).
    * :meth:`update` — the per-cell algorithm.

    Notes
    -----

    Adding a new strategy is now a one-line edit::

        class Step(CellUpdateBase, key="step"):
            is_linear: ClassVar[bool] = True
            is_positivity_preserving: ClassVar[bool] = True

            def update(self, visit, total_xs, source, upstream_state):
                ...

    No registry insert; ``CellUpdateBase.create("step")`` is
    immediately callable.

    Reconstruction ops — the generic affine cell algebra (#158 Inc B / #240 D2)
    ===========================================================================

    A *consistent* affine spatial discretization (Diamond Difference, Linear
    Discontinuous, Step, ...) is fully characterized, per cell, by **three
    :math:`\Sigma_t`-epoch coefficients** returned by
    :meth:`affine_scan_coefficients`:

    * ``a`` — the face transmission multiplier (``ψ_out = a·ψ_in + b``),
    * ``inverse_denom`` — the reciprocal of the cell-balance diagonal ``1/S``,
    * ``w`` — the **cell-average blend weight** (``ψ̄ = (1−w)ψ_in + w·ψ_out``).

    The DAG-free **scan SOLVE** (forward substitution) applies these
    coefficients to *factored* unknowns — it has no concrete inflow until the
    prefix scan has run — so it consumes the three generic staticmethods below
    (:meth:`source_emission`, :meth:`cell_average`,
    :meth:`outgoing_face_from_average`), parameterized by the coefficients and
    carrying NO scheme-specific constant.  These ops are pure functions of the
    face fluxes / source and the blend weight ``w`` — no instance state — so
    they live ON the scheme **base** as ``@staticmethod``\\s: they are the
    GENERIC advection–reaction reconstruction (Step ↔ first-order upwind, DD ↔
    central / Keller box, LD ↔ DG-P1-upwind; the coefficient triple is generic
    in wave-speed + reaction-rate), shared by every affine scheme and inherited
    by the future diffusion scheme.  The discretization math therefore lives in
    exactly one place (the scheme's :meth:`affine_scan_coefficients` for the
    coefficients, these staticmethods for the reconstruction), never duplicated
    in a sweep body (Cardinal Rule 2 / the coefficient-model litmus: *if an
    explicit-matrix representation would have to re-derive a calculation the
    sweep does, that calculation belongs on the scheme*).

    Why these forms are universal
    -----------------------------

    Exactness on a spatially-constant flux (``ψ_in = ψ_out = ψ̄`` under a
    matched source) forces the cell-average to be a **convex** face blend — the
    two weights sum to one — so ``ψ̄ = (1−w)ψ_in + w·ψ_out`` holds for *any*
    consistent affine scheme, with ``w`` the only free per-cell number.  The
    source emission then follows algebraically (SymPy-verified, all schemes):

    .. math::

        b   &= QV \cdot \mathrm{inverse\_denom} / w
           \qquad(\text{DD: } w=\tfrac12 \Rightarrow b = 2\,QV\cdot\mathrm{inverse\_denom})\\
        \bar\psi &= (1-w)\,\psi_{\rm in} + w\,\psi_{\rm out}
           \qquad(\text{DD: } w=\tfrac12 \Rightarrow \bar\psi = \tfrac12(\psi_{\rm in}+\psi_{\rm out}))

    These two ops are in the ×V "denom" convention (``inverse_denom = 1/S``
    with ``S`` the ×V cell-balance diagonal) — the convention the scan cache
    (:class:`~orpheus.sn.spatial.sweep_cache.CollisionCache`) uses.  The
    group-3 ≡ group-2 gate (``test_group3_equals_group2_scan_flat``) pins the
    SOLVE direction (``ψ̄``/``ψ_out``) against the trusted Increment-A kernel.

    The apply direction (matvec) is NOT a generic op here: with a CONCRETE
    probe ``ψ̄`` it rides the scheme's group-2 :meth:`residual_kernel_batch`
    (the ÷V ``g=|μ|/Δ`` form — it returns the density residual AND the outgoing
    face in one call, the natural apply twin of the scan solve).  EVERY affine
    scheme routes its Cartesian matvec through this one kernel UNIFORMLY
    (#158/#240 — no per-scheme matvec branch, no capability flag): Diamond
    Difference reproduces its diamond march ``ψ_out = 2ψ̄ − ψ_in``,
    Linear-Discontinuous its Schur residual (its kernel halves the
    scheme-agnostic ``s = 2|μ|/Δ`` to ``g = |μ|/Δ`` internally).  A future
    ``ExplicitMatrix`` representation would assemble the **×V** coefficients
    into matrix entries (diagonal ``S = 1/inverse_denom``, upstream coupling
    ``(1−w+w·a)/inverse_denom`` — which equals DD's ``2|μ|`` and LD's
    ``|μ|(1+k)``); that is the ×V matrix convention, NOT the ÷V
    ``residual_kernel_batch`` density form (they differ by the cell volume
    ``V``).

    Bit-identity vs principled-equivalence
    --------------------------------------

    **The scan SOLVE ops** (:meth:`source_emission` / :meth:`cell_average`)
    are, for ``w=½`` (Diamond Difference), **byte-identical** to the
    pre-coefficient-model closures: ``0.5*in + 0.5*out`` equals
    ``0.5*(in+out)`` and ``QV·inv/0.5`` equals ``2·QV·inv`` bit-for-bit,
    because multiply/divide by ½ is an exact power-of-2 scaling that commutes
    with IEEE rounding.  So DD's SCAN snapshots stay byte-identical (the
    ``tests/sn/sweep/core tests/sn/solve -W error::DriftWarning`` gate is
    green).

    **The matvec APPLY is a deliberate principled-equivalence re-baseline.**
    DD's Cartesian matvec moved off the ×V ``cell_balance`` density path onto
    the ÷V ``residual_kernel_batch`` kernel (#240, retiring the
    bit-identity-only ``matvec_via_kernel`` special case), which re-associates
    ~1 ULP on non-power-of-2 cell widths.  This is sanctioned by
    ``vv-principles`` §"Bit-identity vs principled-equivalence": bit-identity
    is never a design constraint — the architecture (one uniform matvec kernel,
    no scheme flag) is the compounding asset; a regenerated ~1-ULP golden is
    the negligible cost.  The same principled-equivalence (~1 ULP) governs LD's
    ×V scan vs ÷V kernel two-paths agreement.
    """

    registry: ClassVar[dict[str, type["CellUpdateBase"]]] = {}

    is_linear: ClassVar[bool]
    is_positivity_preserving: ClassVar[bool]
    is_affine_scannable: ClassVar[bool] = False
    """Opt-out default (``False``): a scheme is assumed NOT affine-scannable
    unless it explicitly declares otherwise (mirroring how
    :meth:`cell_kernel_batch` raises by default until a scheme opts into the
    batched wavefront walks).  An affine scheme overrides this to ``True``
    AND supplies the per-cell coefficient triple ``(a, inverse_denom, w)`` via
    :meth:`affine_scan_coefficients` — that is the scheme's *entire* group-3
    SOLVE contribution (#158 the coefficient model).  The generic scan-solve
    *operations* (source emission, cell-average) are NOT per-scheme: they live
    once on the base as the :meth:`source_emission` / :meth:`cell_average`
    staticmethods, parameterized by those coefficients.  The matvec APPLY (a
    concrete probe ψ̄) rides the scheme's
    group-2 :meth:`residual_kernel_batch` UNIFORMLY for every affine scheme
    (Cartesian) — DD reproduces its diamond march, LD its Schur residual, with
    no scheme branch (#158/#240 the coefficient model — there is no per-scheme
    matvec capability flag; bit-identity is never a design constraint, only a
    free bonus when the kernel re-association happens to be exact)."""

    @classmethod
    def _registry_base(cls) -> type:
        return CellUpdateBase

    @abstractmethod
    def update(
        self,
        visit: CellVisit,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> CellResult:
        ...

    @abstractmethod
    def residual(
        self,
        cell_avg: np.ndarray,
        visit: CellVisit,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> np.ndarray:
        r"""Per-cell operator residual :math:`L_{\text{cell}}\,\bar\psi - q`.

        See :meth:`CellUpdate.residual` for the full contract.
        Subclasses MUST implement this in lockstep with :meth:`update`
        — the two methods describe the same per-cell linear system
        in solve direction (``update``) and apply direction
        (``residual``).
        """
        ...

    def cell_kernel_batch(
        self,
        *,
        psi_in: tuple[np.ndarray, ...],
        s_axes: tuple[np.ndarray, ...],
        sigt_cells: np.ndarray,
        Q_cells: np.ndarray,
    ) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
        r"""Pure batched SOLVE cell kernel — the level-vectorised extension point.

        Solve the strategy's cell closure for the cell-average flux on one
        topological level, vectorised over ``(N_oct, ng, n_diag)``: given the
        per-axis incoming face fluxes ``psi_in`` + streaming coefficients
        ``s_axes`` (positional-by-axis, ``d = 1, 2, 3``), the level's
        ``sigt_cells`` ``(ng, n_diag)`` and source ``Q_cells``, return
        ``(psi_avg, psi_out)`` with ``psi_out`` the d-tuple of outgoing face
        fluxes.

        STORAGE-FREE by contract (S6.4(e)): gathering the inputs from and
        scattering ``psi_out`` back into the face cochain is the WALK's job
        (:meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.walk_full` /
        :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.walk_windowed`
        via the ``_CellSolve`` level operation) — a strategy supplies ONLY
        its cell algebra.  Default raises :exc:`NotImplementedError`;
        :class:`DiamondDifference` overrides; future closures (Step / LD /
        EC) override this pair to join the batched wavefront walks (the
        per-cell :meth:`update` stays the canonical reference contract).
        """
        raise NotImplementedError(
            f"{type(self).__name__} does not implement cell_kernel_batch. "
            "Override the kernel pair (cell_kernel_batch / "
            "residual_kernel_batch) to enable the batched wavefront walks."
        )

    def residual_kernel_batch(
        self,
        *,
        psi_bar: np.ndarray,
        psi_in: tuple[np.ndarray, ...],
        s_axes: tuple[np.ndarray, ...],
        sigt_cells: np.ndarray,
        Q_cells: np.ndarray,
    ) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
        r"""Pure batched APPLY cell kernel — the level-vectorised residual.

        The apply-direction companion of :meth:`cell_kernel_batch`, exactly
        as the per-cell :meth:`residual` is the companion of :meth:`update`:
        evaluate the strategy's per-cell operator residual at the GIVEN
        probe cell-average ``psi_bar`` ``(N_oct, ng, n_diag)`` and
        reconstruct the outgoing faces from it.  Returns
        ``(residual, psi_out)``.

        Round-trip contract: at ``psi_bar`` equal to the value
        :meth:`cell_kernel_batch` solves for (same inputs), the residual
        vanishes to FP noise — the batched analogue of the per-cell
        :meth:`residual`/:meth:`update` contract.  Implement the pair in
        LOCKSTEP (one per-level linear system, two directions).  Default
        raises :exc:`NotImplementedError`.
        """
        raise NotImplementedError(
            f"{type(self).__name__} does not implement residual_kernel_batch. "
            "Override the kernel pair (cell_kernel_batch / "
            "residual_kernel_batch) to enable the batched wavefront walks."
        )

    # ── Scan-family capability (Issue #236 §2 — the DAG-free schedules) ──
    #
    # The DAG-free *scan* schedules (CumprodScan, ScanMarch) traverse cells
    # coupled along the sweep direction (cell i inflow = cell i-1 outflow),
    # so they cannot use the independent-batch ``cell_kernel_batch`` (which
    # solves a level of cells whose inflows are ALL already known).  They
    # need instead the closed-form affine recurrence ``ψ_out = a·ψ_in + b``
    # plus the two solve-directions of the closure relation.  A scheme that
    # is ``is_affine_scannable`` overrides this triple; the defaults raise so
    # a non-affine scheme is never silently scanned (the scan strategies'
    # ``supports`` gate keeps it off the scan path, and these raises are the
    # belt-and-braces backstop).

    def affine_scan_coefficients(
        self,
        *,
        abs_mu: np.ndarray,
        A_down: np.ndarray,
        A_total: np.ndarray,
        dA_w: np.ndarray,
        c_out: np.ndarray,
        V: np.ndarray,
        sig_t: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        r""":math:`\Sigma_t`-epoch affine-scan coefficients ``(a, inverse_denom, w)``.

        The scheme's **entire** group-3 contribution (#158 the coefficient
        model): the three per-ordinate-per-cell-per-group coefficients that
        characterize a consistent affine cell, all source-INDEPENDENT (one
        constant-:math:`\Sigma_t` epoch):

        * ``a`` — the closed-form scan (Blelloch §1.5) transmission multiplier
          of the first-order recurrence :math:`\psi^{s}[i+1]=a[i]\psi^{s}[i]+b[i]`;
        * ``inverse_denom`` — the reciprocal cell-balance diagonal :math:`1/S`;
        * ``w`` — the **cell-average blend weight**
          (:math:`\bar\psi=(1-w)\psi_{\rm in}+w\,\psi_{\rm out}`); DD's ``w=½``,
          LD's ``w=1/(1+k)``.

        The source-dependent emission ``b`` and the order-dependent
        ``cumprod_a`` are NOT here — ``cumprod_a`` is a scan-schedule transform
        (:class:`~orpheus.sn.spatial.sweep_cache.CollisionCache`), and ``b``
        plus the cell-average / outgoing-face / matvec-residual operations are
        the representation's job, applied generically from these coefficients by
        the base reconstruction staticmethods (:meth:`source_emission`,
        :meth:`cell_average`, :meth:`outgoing_face_from_average`; the same three
        coefficients feed the scan, the matvec, and a future explicit-matrix
        assembly — no per-scheme closure method, no cell math in any sweep
        body).

        Default raises :exc:`NotImplementedError`; an ``is_affine_scannable``
        scheme (:class:`DiamondDifference`, :class:`LinearDiscontinuous`)
        overrides.  A non-affine scheme leaves it raising — and the scan
        strategies' ``supports`` gate on ``is_affine_scannable`` keeps it
        off the scan path in the first place.
        """
        raise NotImplementedError(
            f"{type(self).__name__} does not implement affine_scan_coefficients "
            "(is_affine_scannable is False).  Only affine-scannable closures "
            "supply the (a, inverse_denom, w) coefficient triple consumed by "
            "CumprodScan / ScanMarch; non-affine closures run on the DAG "
            "wavefront schedule."
        )

    # ── Generic affine reconstruction ops (#158 Inc B / #240 D2) ─────────────
    #
    # Pure functions of the face fluxes / source and the blend weight ``w`` —
    # no instance state, no scheme-specific constant.  The GENERIC
    # advection–reaction reconstruction shared by every affine scheme (and
    # inherited by the future diffusion scheme).  See the class docstring
    # §"Reconstruction ops" for the universality argument and the
    # bit-identity-vs-principled-equivalence note.

    @staticmethod
    def source_emission(
        QV: np.ndarray, inverse_denom: np.ndarray, w: "float | np.ndarray",
    ) -> np.ndarray:
        r"""Affine-scan additive source coefficient ``b = QV · inverse_denom / w``.

        The per-cell source term of the recurrence ``ψ_out = a·ψ_in + b``.  ``QV``
        is the volumetric source on the cell (already weight-normalised per
        ordinate, times cell volume; for a curvilinear sweep it carries the
        Morel–Montry angular contribution folded in by the caller).  DD's
        historical ``2·QV·inverse_denom`` is the ``w=½`` special case.
        """
        return QV * inverse_denom / w

    @staticmethod
    def cell_average(
        face_in: np.ndarray, face_out: np.ndarray, w: "float | np.ndarray",
    ) -> np.ndarray:
        r"""Cell-average from the face pair: ``ψ̄ = (1−w)·ψ_in + w·ψ_out``.

        The universal convex face blend (exactness-on-constants forces the
        weights to sum to one).  DD's ``½(ψ_in + ψ_out)`` is the ``w=½`` special
        case.
        """
        return (1.0 - w) * face_in + w * face_out

    @staticmethod
    def outgoing_face_from_average(
        psi_bar: np.ndarray, face_in: np.ndarray, w: "float | np.ndarray",
    ) -> np.ndarray:
        r"""Downstream face flux from the cell average: ``ψ_out = (ψ̄ − (1−w)·ψ_in)/w``.

        The inverse of :meth:`cell_average` (the convex face blend
        ``ψ̄ = (1−w)·ψ_in + w·ψ_out``).  The universal affine-scheme outflow
        reconstruction: DD's diamond mean ``2ψ̄ − ψ_in`` is the ``w=½`` case;
        LD's ``(1+k)ψ̄ − k·ψ_in`` is the ``w=1/(1+k)`` case.

        .. note::

           For ``w=½`` (Diamond Difference) this is **byte-identical** to the
           inlined ``2·ψ̄ − ψ_in``: ``÷0.5`` is the exact power-of-2 ``×2`` and
           round-to-nearest commutes with exact doubling, so
           ``fl(2·(ψ̄ − 0.5·ψ_in)) == fl(2ψ̄ − ψ_in)`` bit-for-bit.  For LD's
           ``w=1/(1+k)`` it is algebraically equal to ``ψ̄ + (g/θ)(ψ̄ − ψ_in)/D₂``
           but takes a DIFFERENT reduction tree (compute ``w`` then ``÷w`` vs the
           direct ``ψ̄ + k·(…)``), so LD reconstruction is a principled
           ~1-ULP re-baseline (``vv-principles`` §"Bit-identity vs
           principled-equivalence"), not a byte-identical move.
        """
        return (psi_bar - (1.0 - w) * face_in) / w


__all__ = [
    "CellResult",
    "CellUpdate",
    "CellUpdateBase",
    "CellVisit",
    "UpstreamState",
]
