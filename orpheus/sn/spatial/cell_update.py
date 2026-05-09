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

Slab vs curvilinear discrimination
==================================

A strategy distinguishes slab from curvilinear geometry by inspecting
the :class:`~orpheus.geometry.reduced_operator.StreamingTerms`
instance that arrives:

* **Slab / Cartesian** — ``streaming_terms.alpha_in is None``
  (and ``alpha_out``, ``delta_A_over_w``, ``tau_mm``, ``face_area_*``
  are all ``None`` too).  ``upstream_state.angular_upstream`` is
  ``None``.  The strategy returns
  ``CellResult(outgoing_angular_state=None, ...)``.
* **Sphere / cylinder** — ``streaming_terms.alpha_in is not None``
  (and the full curvature bundle is populated).
  ``upstream_state.angular_upstream`` carries
  :math:`\psi_{n-1/2,\,i}`.  The strategy returns the M-M-closed
  ``outgoing_angular_state``.
* **Cylindrical pure-azimuthal degenerate** — when
  ``streaming_terms.abs_mu < 1e-15`` (the level has axial direction
  cosine :math:`|\eta| \to 0`), the cell has no radial face flow.
  The strategy must return ``outgoing_spatial_flux=None`` while still
  applying the M-M angular closure (the redistribution physics
  remains active).

The discrimination convention is locked into Round 1 by the
contract's protocol tests — concrete strategies in Round 2 (and the
sweep rewrite in Wave D) honour it.

Where downstream consumers will call this
=========================================

In Wave D, the sweep at :mod:`orpheus.sn.sweep` dispatches via
:meth:`SNMesh.iter_cell_visits` — the per-visit packets pre-resolve
the sweep direction so the strategy sees no sign-of-:math:`\mu`
branching::

    for visit in sn_mesh.iter_cell_visits(ordinate_idx=n,
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

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np

from orpheus.geometry.reduced_operator import StreamingTerms


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
    :meth:`orpheus.sn.geometry.SNMesh.iter_cell_visits`.  MoC will
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
    face_area_downstream : float | None
        **Sweep-direction-resolved**: which of the two cell faces
        (inner or outer) is the downstream face for this visit.

        * For slab / Cartesian, ``None`` — slab DD is built around
          the chord-length / :math:`|\mu|` recurrence and does not
          read face areas.
        * For curvilinear with :math:`\mu \ge 0` (outward sweep,
          centre → boundary): equals
          ``streaming_terms.face_area_outer``.
        * For curvilinear with :math:`\mu < 0` (inward sweep,
          boundary → centre): equals
          ``streaming_terms.face_area_inner``.
        * For the cylindrical pure-azimuthal degenerate case
          (``abs_mu < 1e-15``), ``None`` — the cell has no spatial
          face flow.

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
    face_area_downstream: float | None = None


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
    :class:`~orpheus.geometry.boundary.ResolvedBC`,
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

    Notes
    -----
    Both traits are class-level so they can be inspected without
    instantiating the strategy.  Code that selects a closure based on
    cell-thickness or stiffness criteria reads ``is_linear``;
    diagnostics that gate on whether negative fluxes can appear read
    ``is_positivity_preserving``.
    """

    is_linear: bool
    is_positivity_preserving: bool

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
            ``visit.streaming_terms.alpha_in is None`` discriminates
            slab from curvilinear geometry; for cylindrical
            geometry, ``visit.streaming_terms.abs_mu < 1e-15`` flags
            the pure-azimuthal degenerate case (in which
            ``visit.face_area_downstream is None``).
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


__all__ = [
    "CellResult",
    "CellUpdate",
    "CellVisit",
    "UpstreamState",
]
