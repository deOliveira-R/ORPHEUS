r"""Per-cell update strategy contract for a transport sweep.

The sweep walks the spatial mesh once per ordinate, and for each cell
it solves the same shape of algebraic problem: combine the local
volumetric source, the upstream face flux (and, in curvilinear
geometry, the upstream angular half-flux) with the local total cross
section and a *spatial closure* — Diamond Difference (DD), Linear
Discontinuous (LD), Exponential Characteristic (EC), Step — to
produce the cell-average flux plus the downstream face flux and (for
sphere/cylinder) the downstream angular half-flux.

The cell-update math is **the same algebra** in slab, sphere, and
cylindrical 1-D — only the populated fields of
:class:`~orpheus.geometry.reduced_operator.StreamingTerms` change; the
remaining ``if`` checks in a strategy test the *structural presence* of a
direction (``face_area_downstream > 0.0``; ``angular_upstream is not
None``), NOT the geometry kind (the geometry-as-data collapse is
``docs/theory/foundations/discretization.rst §discretization-space-angle``).
Lifting the closure out of the sweep body into a
:class:`DiscretizationScheme` strategy (Cardinal Rule 2) lets the
cell-update math be unit-tested in isolation, lets closures be swapped
(DD → LD → EC → Step) without touching the sweep driver, and keeps the
closure-specific algebra local to a strategy module.  Concrete strategies
live in sibling modules
(:mod:`~orpheus.transport.spatial.diamond`,
:mod:`~orpheus.transport.spatial.linear_discontinuous`); the sweep consumes
whichever one the mesh selected, via per-visit :class:`CellVisit` packets
from :meth:`SNMesh.dag_walk` that pre-resolve the sweep direction so the
strategy sees no sign-of-:math:`\mu` branching.  The strategy-pattern
rationale and the geometry-blind closure algebra are documented in
``docs/theory/methods/sn/index.rst`` and
``docs/theory/foundations/discretization.rst``.

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
* See also :doc:`/theory/methods/sn/index`, "Cell update
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
from enum import Enum
from functools import cache
from typing import TYPE_CHECKING, Any, ClassVar, Protocol, runtime_checkable

import numpy as np

from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.reduced_operator import StreamingTerms
from orpheus.numerics.moment_layout import cell_moment_count
from orpheus.numerics.registry import RegistryMixin

if TYPE_CHECKING:  # pragma: no cover
    from ._ubld import D1ClosedForm
    from orpheus.numerics.axis import Axis


# ═══════════════════════════════════════════════════════════════════════
# CellVisit — one visit to one cell during an SN sweep
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True, slots=True)
class CellVisit:
    r"""One visit to one cell during a transport sweep.

    A cell-graph transport sweep (S\ :sub:`N` is the producer today) is
    a topological sort of the **directed cell graph**
    for a given ordinate, where edges are oriented by
    :math:`\mathrm{sign}(\Omega \cdot \hat n_{\text{face}})`.  This
    dataclass is the per-visit packet a
    :class:`DiscretizationScheme` strategy receives — all sweep-direction data
    is **pre-resolved** here.  The strategy sees only the
    upstream/downstream view of this cell, never the geometric
    inner/outer view.

    A **transport-layer** per-cell visit packet — method-generic, not
    SN-only.  It lives in :mod:`orpheus.transport.spatial`, and its
    sweep-resolved angular-closure fields are plain floats with
    **Cartesian-neutral defaults** (:attr:`c_in` / :attr:`c_out`
    ``= 0.0``, :attr:`tau` ``= 1.0``), so a cell-graph sweep in
    Cartesian geometry populates only ``cell_idx`` +
    ``streaming_terms`` and inherits the neutral rest.  SN's mesh is
    the producer today
    (:meth:`orpheus.sn.mesh.augmented_mesh.SNMesh.dag_walk` stamps
    it); any other cell-graph transport method may stamp the same
    packet.  MoC is deliberately **not** such a consumer (its per-ray
    traversal has a different mathematical structure — fiber bundles /
    solution sheaves, not a topological sort over a cell graph), so
    there is no shared ``SweepGraph`` Protocol.

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
    c_in : float
        **Angular-closure-owned**: the Morel--Montry weighted-diamond
        upstream-numerator constant
        :math:`c_{\rm in} = (1-\tau_m)/\tau_m\,\alpha_{m+1/2}
        + \alpha_{m-1/2}` for THIS visit's ordinate, stamped on the visit
        by the mesh from its bound ``pole_angular_closure``
        (:class:`~orpheus.sn.angular.closure.PoleAngularClosureBase`).
        ``0.0`` for slab / Cartesian (the identity closure carries no
        angular redistribution).  The spatial scheme consumes it as DATA
        on the visit; it must NOT be rebuilt from ``streaming_terms``
        (that would re-fuse the spatial and angular closures).
    c_out : float
        **Angular-closure-owned**: the Morel--Montry weighted-diamond
        denominator constant :math:`c_{\rm out} = \alpha_{m+1/2}/\tau_m`
        for THIS visit's ordinate, stamped by the mesh (see :attr:`c_in`).
        ``0.0`` for slab / Cartesian.
    tau : float
        **Angular-closure-owned**: the Morel--Montry angular weight
        :math:`\tau_m` for THIS visit's ordinate — the FUNDAMENTAL weight
        (Bailey--Morel--Chang 2010 Eq. 43) from which :attr:`c_in` /
        :attr:`c_out` are derived, stamped by the mesh (see :attr:`c_in`).
        The DEFAULT is ``1.0`` (NOT ``0.0``): :math:`\tau = 1` is the
        neutral M-M weight the Cartesian identity closure supplies, making
        the angular recurrence
        :math:`(\bar\psi - (1-\tau)\psi^{\theta}_{\rm in})/\tau` the
        identity and keeping :math:`c_{\rm out} = \alpha_{\rm out}/\tau`
        well-defined.  ⚠ Do NOT default :attr:`tau` to ``0.0`` (like
        :attr:`c_in` / :attr:`c_out`) — it is a divide-by-zero landmine in
        :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.update`'s
        angular thread.

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
    c_in: float = 0.0
    c_out: float = 0.0
    tau: float = 1.0


# ═══════════════════════════════════════════════════════════════════════
# Per-cell input state
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True, slots=True)
class UpstreamState:
    r"""Per-cell input state that a :class:`DiscretizationScheme` consumes.

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
    r"""Per-cell output state produced by a :class:`DiscretizationScheme`.

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
# Face-transmission damping — the gauge-freedom half that belongs to the
# CLOSURE (the other half is the boundary set; see
# ``SNMesh.reflective_axis_pairs``)
# ═══════════════════════════════════════════════════════════════════════

class FaceModeDamping(Enum):
    r"""Does this closure damp every mode of its own face transmission?

    The three states are genuinely distinct and the third is not a
    failure to be swallowed: ``UNDETERMINED`` means *we could not ask*,
    which is neither "there is freedom" nor "there is none".
    """

    DAMPED = "damped"
    r"""Every face mode decays.  A reflective loop cannot sustain one."""

    UNDAMPED = "undamped"
    r"""Some face mode survives transmission unattenuated (:math:`|\lambda| = 1`).
    Around a closed reflective loop it returns to itself, so the assembled
    loss operator acquires a null space — a GAUGE FREEDOM."""

    UNDETERMINED = "undetermined"
    r"""The closure could not be driven at this dimension.  Callers MUST
    treat this as "unknown", never as :attr:`DAMPED` — the honest response
    is to say so and leave the answer ungauged."""


@dataclass(frozen=True, slots=True)
class FaceTransmissionSpectrum:
    r"""The measured verdict of :meth:`DiscretizationSchemeBase.face_transmission_spectrum`.

    Attributes
    ----------
    damping :
        The three-state verdict.
    spectral_radius :
        :math:`\rho(\Sigma)` — the measured number the verdict rests on.
        ``None`` iff :attr:`FaceModeDamping.UNDETERMINED` (nothing was
        measured; per the ``balance_defect`` discipline, ``None`` means
        *not measured*, never *measured and small*).
    undetermined_because :
        Why the closure could not be driven — carried so the caller's
        warning can name the real obstruction instead of guessing.
        ``None`` unless ``UNDETERMINED``.
    """

    damping: FaceModeDamping
    spectral_radius: float | None = None
    undetermined_because: str | None = None


# ``rho`` is measured EXACTLY 1.0 for an undamped closure (the sawtooth
# eigenvalue is -1 by construction, not by cancellation), so the tolerance
# only absorbs the eigensolver's own noise.
_UNDAMPED_ATOL = 1e-9

# Two deliberately DIFFERENT probe cells.  The classification is a property
# of the CLOSURE, not of a cell — so the two must agree, and a scheme whose
# verdict is cell-dependent is reported UNDETERMINED rather than classified
# from a single lucky draw (`vv-principles` #13: one sample is not a survey).
# `[M]` 2026-08-14 diamond_difference reads rho = 1.0000000000 on BOTH;
# linear_discontinuous reads 0.6373008435 / 0.8607021518 — the margin moves
# with the cell, the verdict does not.
_PROBE_CELLS: tuple[tuple[tuple[float, float, float], float], ...] = (
    ((0.7, 1.3, 0.4), 0.9),   # asymmetric streaming, thick
    ((1.0, 1.0, 1.0), 0.3),   # symmetric streaming, thin
)


def _face_transmission_matrix(
    scheme: "DiscretizationSchemeBase",
    ndim: int,
    w: tuple[float, ...],
    sigma_t_volume: float,
) -> np.ndarray:
    r"""Assemble :math:`\Sigma` by driving the closure one unit inflow at a time.

    Column :math:`j` is the outgoing face vector produced by a
    source-free cell whose inflow is the :math:`j`-th face-moment unit
    vector — i.e. the closure's own kernel IS the definition, read
    column by column.  Square of side ``ndim * moments_per_face``.
    """
    moments_per_face = scheme.spatial_basis_per_axis ** (ndim - 1)
    face_shape = (
        (1, 1, 1, moments_per_face) if moments_per_face > 1 else (1, 1, 1)
    )
    s_axes = tuple(np.full((1, 1, 1), float(w[a])) for a in range(ndim))
    reaction_xs = np.full((1, 1), float(sigma_t_volume))
    q_cells = np.zeros((1, 1))

    columns = []
    for j in range(ndim * moments_per_face):
        inflow_axis, moment = divmod(j, moments_per_face)
        psi_in = []
        for axis in range(ndim):
            face = np.zeros(face_shape)
            if axis == inflow_axis:
                face.reshape(-1)[moment] = 1.0
            psi_in.append(face)
        _avg, psi_out = scheme.cell_kernel_batch(
            psi_in=tuple(psi_in), s_axes=s_axes,
            reaction_xs=reaction_xs, Q_cells=q_cells,
        )
        columns.append(np.concatenate([np.ravel(o) for o in psi_out]))
    return np.asarray(columns).T


@cache
def _face_transmission_spectrum(
    scheme_type: "type[DiscretizationSchemeBase]", ndim: int,
) -> FaceTransmissionSpectrum:
    r"""Cached per ``(closure, ndim)`` — see :meth:`face_transmission_spectrum`.

    Keyed on the TYPE, not the instance: the verdict is a property of the
    closure and `[M]` does not move with the probe cell, so one
    eigendecomposition serves every mesh, group and iterate.
    """
    scheme = scheme_type()
    radii: list[float] = []
    for w, sigma_t_volume in _PROBE_CELLS:
        try:
            sigma = _face_transmission_matrix(
                scheme, ndim, w[:ndim], sigma_t_volume,
            )
        except Exception as exc:  # the closure cannot be driven here
            return FaceTransmissionSpectrum(
                damping=FaceModeDamping.UNDETERMINED,
                undetermined_because=(
                    f"{scheme_type.__name__} could not be driven at "
                    f"ndim={ndim}: {type(exc).__name__}: {exc}"
                ),
            )
        radii.append(float(np.abs(np.linalg.eigvals(sigma)).max()))

    verdicts = {r >= 1.0 - _UNDAMPED_ATOL for r in radii}
    if len(verdicts) != 1:
        # A cell-DEPENDENT verdict means this predicate does not describe
        # the closure, only the draw — refuse rather than pick one.
        return FaceTransmissionSpectrum(
            damping=FaceModeDamping.UNDETERMINED,
            undetermined_because=(
                f"{scheme_type.__name__} at ndim={ndim} classifies "
                f"differently on the two probe cells (rho = "
                f"{', '.join(f'{r:.10f}' for r in radii)}), so the damping "
                f"is a property of the cell rather than of the closure."
            ),
        )
    return FaceTransmissionSpectrum(
        damping=(
            FaceModeDamping.UNDAMPED if verdicts.pop()
            else FaceModeDamping.DAMPED
        ),
        spectral_radius=max(radii),
    )


# ═══════════════════════════════════════════════════════════════════════
# DiscretizationScheme Protocol
# ═══════════════════════════════════════════════════════════════════════

@runtime_checkable
class DiscretizationScheme(Protocol):
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
        :math:`\psi_{\rm out} = a\,\psi_{\rm in} + b` (Blelloch §1.5), the
        cell-average being an affine function of a **single** upstream face
        flux, so the DAG-free *scan* schedules (``CumprodScan``,
        ``ScanMarch``) can consume it.  ``True`` for Diamond Difference;
        ``True`` for Linear-Discontinuous too — its two face moments couple,
        but the local slope is eliminated by the Schur complement, leaving a
        single-upstream recurrence.  ``False`` is the safe default (the scan
        strategies gate their ``supports`` on this trait; a non-scannable
        scheme is routed to the DAG wavefront).  A **single-axis** (1-D
        prefix-scannability) statement, distinct from
        ``transverse_coupling_is_facewise`` (the cross-axis statement).
        Read-only class attribute.
    transverse_coupling_is_facewise : bool
        Whether, in :math:`d \ge 2`, a NON-swept axis couples through a
        **0th-order face value** (a Dirichlet trace) rather than a
        **1st-order slope moment** — equivalently, whether the d-D cell
        closure is **tensor-product separable** into independent per-axis
        1-D scans chained by scalar face traces.  ``True`` for the slopeless
        cell-average closures (Diamond Difference / Step): the transverse
        term folds into the scan's affine source, so a row-march
        (``ScanMarch``) admits them.  ``False`` (the default) for the
        bilinear DG-P1 Linear-Discontinuous closure: its transverse coupling
        is a per-axis slope moment that CANNOT reduce to a single face value,
        so the d-D LD closure is irreducibly axis-coupled and rides the DAG
        wavefront.  Distinct from the single-axis ``is_affine_scannable``;
        the :math:`d \ge 2` scan-march gates on THIS trait
        (``docs/theory/methods/sn/loss_representation.rst §loss-rep-scanmarch-facewise``).
        Read-only class attribute.
    spatial_basis_per_axis : int
        The number of spatial moments the cell closure carries **per spatial
        axis** (its 1-D moment-basis size).  ``1`` for the slopeless
        cell-average closures (Diamond Difference / Step — basis ``{1}``);
        ``2`` for Linear-Discontinuous (Legendre ``{1, P₁}`` per axis).  The
        tensor-product (Kronecker) UBLD structure makes the per-cell
        (:math:`(\text{per\_axis})^d`) and per-face
        (:math:`(\text{per\_axis})^{d-1}`) moment counts dimension-dependent
        derivations of this single per-axis number (per-AXIS, not per-cell:
        LD's per-cell count is dimension-dependent).  The boolean
        ``per_axis > 1`` is the GATE on every multi-moment wiring site;
        DD/Step at ``per_axis == 1`` keep the scalar face + ``psi_avg``
        byte-identical (no trailing length-1 axis appended — #240 D5b).
        Derivation:
        ``docs/theory/methods/sn/cartesian_multid.rst §spatial-moment-space``.
        Read-only class attribute.
    diffusion_limit_consistent : bool
        Whether the scheme's THICK-DIFFUSION limit is a *consistent*
        diffusion discretization for the leading-order scalar flux — the
        SPATIAL half of the (spatial ⊗ angular) diffusion-limit condition.
        ``True`` for Diamond Difference (Larsen–Morel–Miller 1987
        Eq. (4.24)) and full Linear-Discontinuous (Larsen–Morel 1989 Part II
        Eq. (4.16), which requires the slope SOURCE moment :math:`\hat Q`
        threaded); ``False`` (the default) for Step (LMM-1987 Eq. (5.20)).
        This is the SPATIAL axis ONLY; the ANGULAR first-order condition
        lives on the redistribution closure as
        ``beta_first_order_consistent``, and a (scheme, closure) PAIR is
        valid iff both hold
        (:func:`~orpheus.sn.sweep.pairing.pair_diffusion_limit_consistent`).
        ⚠ Do NOT conflate the spatial DD verdict (``True``) with
        DD-in-ANGLE's first-order :math:`\beta`-failure (the curvilinear flux
        dip) — that is an angular-axis result.  Read-only class attribute.
    supports_curvilinear : bool
        Whether the scheme has a CURVILINEAR (sphere/cylinder) cell closure —
        i.e. whether :meth:`update` / :meth:`residual` handle the Morel–Montry
        angular-redistribution thread (a non-``None`` ``angular_upstream``).
        ``True`` for Diamond Difference; ``False`` (the default) for
        Linear-Discontinuous (the curvilinear LD closure is not yet implemented —
        Issue #158 curvilinear arm / #6; published: Adams-Martin 1992 App. A).  The 1-D sweep-strategy selection
        reads this so a curvilinear mesh paired with a slab-only scheme is
        rejected at SELECTION (a clear reason), rather than passing
        ``supports()`` on the geometry-blind ``is_affine_scannable`` and
        raising mid-sweep.  Read-only class attribute.

    Notes
    -----
    All seven traits are class-level so they can be inspected without
    instantiating the strategy: ``is_linear`` gates stiffness-based closure
    selection; ``is_positivity_preserving`` is the cone-preservation flag
    (honest values, no production reader yet — #390 tracks wiring the
    element predicate ``Field.cone_violations`` as its first consumer);
    ``is_affine_scannable`` gates the 1-D sweep-schedule; and
    ``transverse_coupling_is_facewise`` gates the multi-D scan-march (and a
    future diffusion ADI / line-SOR preconditioner — #240's confirmed next
    consumer).
    """

    is_linear: bool
    is_positivity_preserving: bool
    is_affine_scannable: bool
    transverse_coupling_is_facewise: bool
    spatial_basis_per_axis: int
    diffusion_limit_consistent: bool
    supports_curvilinear: bool

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

    # NOTE: the scan-family capability methods — ``affine_scan_coefficients``
    # (1-D / curvilinear chain), ``cartesian_scan_coefficients`` +
    # ``reflect_scan_coefficients`` (the multi-D row-march), and the batched
    # kernel pair ``cell_kernel_batch`` / ``residual_kernel_batch`` — are NOT
    # part of this Protocol.  They are OPT-IN capabilities (guarded by the
    # ``is_affine_scannable`` / ``transverse_coupling_is_facewise`` traits and
    # the strategy ``supports`` gates), carried on
    # :class:`DiscretizationSchemeBase` with raising defaults.  The Protocol
    # declares ONLY the universal per-cell contract (``update`` / ``residual``)
    # every scheme MUST provide, so a minimal synthetic scheme passes
    # ``@runtime_checkable`` ``isinstance`` on the two required members alone.


# ═══════════════════════════════════════════════════════════════════════
# DiscretizationSchemeBase — concrete ABC carrying RegistryMixin
# ═══════════════════════════════════════════════════════════════════════


class DiscretizationSchemeBase(RegistryMixin, ABC):
    r"""Concrete abstract base for self-registering cell-update strategies.

    Layered on :class:`RegistryMixin`: strategies that inherit it self-register
    under their ``key=`` class-creation kwarg, gaining name-keyed lookup via
    :meth:`create`.  The :class:`DiscretizationScheme` Protocol stays alongside
    — third-party / synthetic strategies (e.g. test mocks) that do not want to
    inherit can still satisfy the contract by providing the right methods.

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

    Adding a new strategy is a one-line edit::

        class Step(DiscretizationSchemeBase, key="step"):
            is_linear: ClassVar[bool] = True
            is_positivity_preserving: ClassVar[bool] = True

            def update(self, visit, total_xs, source, upstream_state):
                ...

    No registry insert; ``DiscretizationSchemeBase.create("step")`` is
    immediately callable.

    Σ-stateless advection–reaction contract
    =======================================

    The coefficient / kernel methods (:meth:`affine_scan_coefficients`,
    :meth:`cell_kernel_batch`, :meth:`residual_kernel_batch`) hold **NO**
    cross-section state — a scheme instance carries no :math:`\Sigma`, no
    materials, no mesh.  The reaction-rate is an EXPLICIT argument
    (``reaction_xs`` on the coefficient + kernel methods; ``total_xs`` on the
    per-cell :meth:`update` / :meth:`residual`), so any consumer may pass an
    arbitrary reaction-rate: SN's :math:`\Sigma_t`, an off-diagonal removal
    :math:`\Sigma_r = \Sigma_t - \Sigma_{s0}`, a diffusion-removal, a
    pure-advection :math:`\Sigma \to 0`.  This is the **model-agnostic
    advection–reaction** spatial discretization the diffusion solver will
    consume (standalone AND as the consistent-DSA preconditioner, #2).
    ``streaming`` / ``s_axes`` are kept distinct from ``reaction_xs`` — the one
    genuine SN/CFD frame conflict (diffusion has no advective :math:`\mu`
    coefficient).  The Step ↔ upwind / DD ↔ central-box / LD ↔ DG-P1-upwind
    reading and the CFD Péclet blend are
    ``docs/theory/foundations/discretization.rst §discretization-closures``;
    the Σ-statelessness teeth are
    :mod:`tests.transport.spatial.test_scheme_reaction_rate_contract`.

    Generic affine reconstruction ops
    =================================

    A consistent affine scheme is fully characterized, per cell, by the three
    :math:`\Sigma_t`-epoch coefficients ``(a, inverse_denom, w)`` from
    :meth:`affine_scan_coefficients` (``a`` the face transmission multiplier,
    ``inverse_denom`` the reciprocal cell-balance diagonal, ``w`` the
    cell-average blend weight — DD's ``½``, LD's ``1/(1+k)``).  Both the scan
    SOLVE and the matvec APPLY consume those coefficients through the generic
    base staticmethods (:meth:`source_emission`, :meth:`cell_average`,
    :meth:`outgoing_face_from_average`; the reverse walks ride the matching
    VJP :meth:`outgoing_face_from_average_transpose`, #311) — pure functions
    of the faces / source and ``w``, no instance state — so the
    discretization math lives in exactly
    one place per scheme (its :meth:`affine_scan_coefficients` + these shared
    staticmethods), never duplicated in a sweep body (Cardinal Rule 2).  The
    universality of the convex face blend
    :math:`\bar\psi = (1-w)\psi_{\rm in} + w\,\psi_{\rm out}` (exactness on a
    constant flux forces the two weights to sum to one, leaving ``w`` the only
    free per-cell number) and its outflow inverse are
    ``docs/theory/foundations/discretization.rst §discretization-closures``.

    The coefficients are in the ×V "denom" convention (``inverse_denom = 1/S``
    with ``S`` the ×V cell-balance diagonal) — the convention the scan cache
    (:class:`~orpheus.sn.sweep.cache.CollisionCache`) uses; the matvec APPLY
    instead rides the scheme's ÷V ``g=|μ|/Δ`` :meth:`residual_kernel_batch`
    UNIFORMLY (every affine scheme, no per-scheme matvec branch — #158/#240).

    At ``w=½`` (Diamond Difference) the scan SOLVE ops are **byte-identical**
    to the pre-coefficient-model closures (``÷½`` is an exact power-of-2
    ``×2``), so DD's scan snapshots stay strict (the
    ``tests/sn/sweep/core -W error::DriftWarning`` gate).  DD's matvec APPLY is
    a deliberate principled ~1-ULP re-baseline onto the ÷V kernel (one uniform
    matvec kernel, no per-scheme flag), sanctioned by ``vv-principles``
    §"Bit-identity vs principled-equivalence"; the same governs LD's ×V-scan vs
    ÷V-kernel agreement.
    """

    registry: ClassVar[dict[str, type["DiscretizationSchemeBase"]]] = {}

    is_linear: ClassVar[bool]
    is_positivity_preserving: ClassVar[bool]
    is_affine_scannable: ClassVar[bool] = False
    """Opt-out default (``False``): a scheme is assumed NOT affine-scannable
    unless it declares otherwise.  An affine scheme overrides to ``True`` AND
    supplies the per-cell coefficient triple ``(a, inverse_denom, w)`` via
    :meth:`affine_scan_coefficients` — its *entire* scan SOLVE contribution
    (#158 the coefficient model); the generic reconstruction staticmethods on
    this base do the rest (see the class docstring §"Generic affine
    reconstruction ops" and the Protocol trait contract for what qualifies)."""

    transverse_coupling_is_facewise: ClassVar[bool] = False
    r"""Opt-in (``False`` default): does the multi-D cell closure couple a
    NON-swept axis only through a **0th-order face value** (so the :math:`d`-D
    update is tensor-product separable into per-axis 1-D scans chained by scalar
    traces)?  ``True`` for Diamond Difference / Step (facewise); LD opts OUT
    (its transverse coupling is a 1st-order slope moment — bilinear,
    non-separable — so it rides the DAG wavefront).  ⚠ This is the CROSS-AXIS
    claim, distinct from the single-axis ``is_affine_scannable``: overloading
    the 1-D trait to license the multi-D scan-march silently misroutes a 2-D LD
    mesh into inline DD (#240 D5-0), which minting this trait makes
    unrepresentable.  A diffusion ADI / line-SOR preconditioner (#240's next
    consumer) reads this SAME predicate — the trait is named for the scheme
    property, not the ``ScanMarch`` strategy.  Derivation:
    ``docs/theory/methods/sn/loss_representation.rst §loss-rep-scanmarch-facewise``.
    Read-only class attribute."""

    spatial_basis_per_axis: ClassVar[int] = 1
    r"""The cell closure's number of spatial moments **per spatial axis** (its
    1-D moment-basis size); default ``1`` for the slopeless cell-average
    closures (Diamond Difference / Step), ``2`` for Linear-Discontinuous.  The
    tensor-product (Kronecker) UBLD structure makes the per-cell
    (:math:`(\text{per\_axis})^d`) and per-face
    (:math:`(\text{per\_axis})^{d-1}`) moment counts dimension-dependent
    derivations of this single per-axis number (per-AXIS, not per-cell: LD's
    per-cell count is dimension-dependent).  The boolean ``per_axis > 1`` is the
    GATE on every multi-moment wiring site (#240 D5b); DD/Step at
    ``per_axis == 1`` keep the scalar face + rank-3 ``psi_avg`` byte-identical
    (no trailing length-1 axis appended — a uniform length-1 axis would
    re-associate DD's emit einsum and break bit-identity).  Derivation (the
    Kronecker layout + the append policy):
    ``docs/theory/methods/sn/cartesian_multid.rst §spatial-moment-space``; the
    d=1 reduction is :mod:`orpheus.transport.spatial._ubld` (``d1_closed_form``).
    Read-only class attribute."""

    diffusion_limit_consistent: ClassVar[bool] = False
    r"""Whether the scheme's thick-diffusion limit is a consistent diffusion
    discretization for the leading-order scalar flux — the SPATIAL half of the
    (spatial ⊗ angular) diffusion-limit condition.  Opt-in (``False`` default).
    ``True`` for Diamond Difference (Larsen–Morel–Miller 1987 Eq. (4.24)) and
    full Linear-Discontinuous (Larsen–Morel 1989 Part II Eq. (4.16) — requires
    the slope SOURCE moment); ``False`` for Step (LMM-1987 Eq. (5.20)).  The
    ANGULAR first-order condition is the redistribution closure's
    ``beta_first_order_consistent`` (Bailey–Morel–Chang 2010 Eq. (42)); the
    PAIR's validity is their conjunction
    (:func:`~orpheus.sn.sweep.pairing.pair_diffusion_limit_consistent`).  ⚠ The
    spatial DD verdict (``True``) is NOT the DD-in-angle :math:`\beta`-failure
    (the curvilinear flux dip is an ANGULAR artefact — BMC 2010, not LMM 1987).
    Read-only class attribute."""

    supports_curvilinear: ClassVar[bool] = False
    r"""Whether the scheme has a curvilinear (sphere/cylinder) cell closure
    (handles the Morel–Montry angular-redistribution thread).  Opt-in
    (``False`` default — slab/Cartesian only until a curvilinear closure is
    declared).  The 1-D sweep-strategy selection reads this so a curvilinear
    mesh with a slab-only scheme is rejected at SELECTION (a clear
    ``Compatibility`` reason), not passed on the geometry-blind
    ``is_affine_scannable`` trait and raised mid-sweep.  ``True`` for Diamond
    Difference; ``False`` for Linear-Discontinuous (the curvilinear LD closure
    is not yet implemented — #158/#6).  Read-only class attribute."""

    has_transpose_kernel: ClassVar[bool] = False
    r"""Whether this scheme REGISTERS its transpose kernels — DERIVED, never declared.

    ``True`` iff the subclass overrides
    :meth:`residual_kernel_batch_transpose` (the Cartesian batch VJP — the
    reverse walks' scheme-uniform kernel) AND, if it claims
    :attr:`supports_curvilinear`, also :meth:`streaming_cell_transpose` (the
    1-D curvilinear cell-balance VJP).  Set by ``__init_subclass__`` from the
    overrides themselves, so a declared capability with no kernel behind it
    (the pre-2.5a "predicate lie") is unrepresentable — the registration IS
    the flag (#310 ruling 2), and the registration must COVER the scheme's
    own forward span (#310 C2).  ``True`` for Diamond Difference (both
    kernels); ``False`` for Linear-Discontinuous until its UBLD
    Schur-residual batch VJP registers (#310 C2 — LD is slab-only, so the
    batch kernel alone covers its span).
    Consumers:
    :attr:`~orpheus.sn.operators.streaming.StreamingOperator.is_adjointable`
    (an eager ``.H`` on a non-transposable scheme raises ``MissingAdjoint`` at
    construction, Pattern 4) and the reverse walk's entry guard (a typed
    ``NotImplementedError``, never a silent scalar-buffer broadcast against
    moment-tailed cotangents).  The discrete Euclidean reverse-DAG transpose is
    ``docs/theory/methods/sn/loss_representation.rst §loss-rep-orientation-two-frames``."""

    def __init_subclass__(cls, **kwargs: Any) -> None:
        super().__init_subclass__(**kwargs)
        # Registration-coupled trait (#310 ruling 2): the flag DERIVES from
        # the transpose-kernel overrides, so flag and kernels cannot drift
        # apart — declaring the capability without registering the kernels
        # (or vice versa) is unrepresentable.  The registration must cover
        # the scheme's own forward span (#310 C2): the Cartesian batch VJP
        # always, plus the curvilinear cell-balance VJP iff the scheme
        # claims curvilinear support.
        cls.has_transpose_kernel = (
            cls.residual_kernel_batch_transpose
            is not DiscretizationSchemeBase.residual_kernel_batch_transpose
        ) and (
            not cls.supports_curvilinear
            or cls.streaming_cell_transpose
            is not DiscretizationSchemeBase.streaming_cell_transpose
        )

    @property
    def is_multi_moment(self) -> bool:
        r"""``True`` iff this closure carries within-cell spatial moments.

        The named, scheme-level form of the ``per_axis > 1`` gate the
        multi-moment wiring keys on (#240 D5b): a closure is multi-moment iff
        it tracks a within-cell slope (``ψ̂``) beyond the cell-average — LD →
        ``True``; DD / Step → ``False``.  This is the UNCONDITIONAL scheme
        truth (it reads the closure's basis size, not any field provenance), so
        the inner-walk moment-buffer checks route through it rather than
        re-spelling ``spatial_basis_per_axis > 1`` (single source of truth).
        Equivalent to
        :func:`~orpheus.transport.spatial._ubld.octant_moment_frame_signs`
        ``is not None`` at a streaming octant and to the 1-D scan's
        ``moment_tail != ()`` width-presence test."""
        return self.spatial_basis_per_axis > 1

    def face_transmission_spectrum(self, ndim: int) -> FaceTransmissionSpectrum:
        r"""Does this closure leave any face mode UNDAMPED? — ASKED, never declared.

        Drives the scheme's own :meth:`cell_kernel_batch` on a
        source-free probe cell, one unit inflow face-moment at a time, to
        assemble the face-to-face transmission :math:`\Sigma` and return
        :math:`\rho(\Sigma)`.  Nothing is tabulated: a future closure
        answers for itself with no edit here (the *ask-don't-tabulate*
        ruling, Q6-E — the sibling of :attr:`has_transpose_kernel`'s
        "DERIVED, never declared", one level up: that one reads which
        methods exist, this one reads what they compute).

        **Why the spectral radius is the whole question.**  For a
        single-moment closure the transmission is
        :math:`\Sigma = (2/D)\,\mathbf{1}w^{\mathsf T} - I` with
        :math:`w_a = 2|\mu_a|A_a`, :math:`D = \Sigma_t V + \sum_b w_b`,
        whose spectrum is

        =========================  ==============================  =============
        eigenvalue                 eigenvector                     multiplicity
        =========================  ==============================  =============
        :math:`1 - 2\Sigma_t V/D`  :math:`\mathbf 1`                          1
        :math:`-1`                 :math:`\{v: w^{\mathsf T}v = 0\}`  :math:`d-1`
        =========================  ==============================  =============

        The first is the physical, **absorption-damped** mode.  The
        second is a **cell-average-blind sawtooth**: it drives
        :math:`\psi_c = 0`, so the absorption term :math:`\Sigma_t V\psi_c`
        never sees it and it transmits at unit magnitude.  So
        :math:`\rho(\Sigma) \ge 1` **iff** some mode is undamped — no
        subspace restriction is needed, which is what lets this work
        unchanged on a moment-tailed face (LD's face carries
        ``spatial_basis_per_axis ** (ndim - 1)`` moments) and makes
        ``ndim = 1`` fall out on its own rather than as a special case.

        Around a **closed reflective loop** two :math:`-1`\ s compose to
        :math:`+1`, so an undamped mode returns to itself and the
        assembled loss operator :math:`A = L + C - S - B` acquires a null
        space.  That — together with the boundary set
        (:attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.reflective_axis_pairs`)
        — is the gauge freedom the :class:`LossKernelGauge` fixes.

        `[M]` 2026-08-14, on both probe cells:

        ==================== ============ ============ ==================
        scheme               ``ndim=1``   ``ndim=2``   ``ndim=3``
        ==================== ============ ============ ==================
        diamond_difference   damped       **1.0**      **1.0**
        linear_discontinuous damped       0.637/0.861  ``UNDETERMINED``
        ==================== ============ ============ ==================

        Diamond reads **exactly** ``1.0000000000`` on both probes; LD's
        margin moves with the cell but its verdict does not — the
        classification is a property of the CLOSURE, which is why two
        probe cells are run and required to agree.

        ⚠ ``UNDETERMINED`` is a third state, not a synonym for
        :attr:`FaceModeDamping.DAMPED`.  LD at ``ndim = 3`` raises from
        its own ``assemble_inflow_axis`` (``axis in {0, d-1}`` only) — a
        pre-existing LD limitation.  A caller that cannot classify MUST
        say so and leave the answer ungauged; silently assuming "damped"
        is the blindness this whole surface exists to remove.
        """
        return _face_transmission_spectrum(type(self), int(ndim))

    @classmethod
    def _registry_base(cls) -> type:
        return DiscretizationSchemeBase

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

        See :meth:`DiscretizationScheme.residual` for the full contract.
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
        reaction_xs: np.ndarray,
        Q_cells: np.ndarray,
    ) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
        r"""Pure batched SOLVE cell kernel — the level-vectorised extension point.

        Solve the strategy's cell closure for the cell-average flux on one
        topological level, vectorised over ``(N_oct, ng, n_diag)``: given the
        per-axis incoming face fluxes ``psi_in`` + streaming coefficients
        ``s_axes`` (positional-by-axis, ``d = 1, 2, 3``), the level's
        ``reaction_xs`` ``(ng, n_diag)`` and source ``Q_cells``, return
        ``(psi_avg, psi_out)`` with ``psi_out`` the d-tuple of outgoing face
        fluxes.

        STORAGE-FREE by contract (S6.4(e)): gathering the inputs from and
        scattering ``psi_out`` back into the face cochain is the WALK's job
        (:meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_full` /
        :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed`
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
        reaction_xs: np.ndarray,
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

    def residual_kernel_batch_transpose(
        self,
        *,
        res_bar: np.ndarray,
        psi_out_bar: tuple[np.ndarray, ...],
        s_axes: tuple[np.ndarray, ...],
        reaction_xs: np.ndarray,
    ) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
        r"""Pure batched VJP of :meth:`residual_kernel_batch` — the reverse-mode pair.

        Given the residual cotangent ``res_bar`` and the d-tuple of
        downstream-face cotangents ``psi_out_bar`` (the cotangents of the
        forward's ``(residual, psi_out)`` outputs — same batch layout, same
        positional-by-axis convention), return ``(psi_bar_cot, psi_in_cots)``
        — the cotangents of the forward's ``(psi_bar, psi_in)`` inputs —
        computed from the SAME :math:`\Sigma_t`-epoch data ``(s_axes,
        reaction_xs)`` the forward consumed (Pattern 2 — no twin
        coefficients).  The matvec transpose is source-free: the forward's
        ``Q_cells`` enters ``−R`` and its cotangent belongs to the solve
        path, so no source cotangent is returned.

        The scheme-uniform CARTESIAN reverse kernel: the 1-D reverse loop
        walk's Cartesian arm (#310 C2) and the multi-D reverse wavefront's
        ``_CellResidualTranspose`` level op (#310 C3) both bottom here — the
        exact mirror of the forward arms bottoming in
        :meth:`residual_kernel_batch`.  Overriding it is one conjunct of the
        derived :attr:`has_transpose_kernel` registration; a scheme claiming
        :attr:`supports_curvilinear` must ALSO register
        :meth:`streaming_cell_transpose` (the curvilinear cell-balance arm's
        kernel).  Default raises.
        """
        raise NotImplementedError(
            f"{type(self).__name__} registers no batch transpose kernel "
            "(residual_kernel_batch_transpose) — the adjoint matvec on this "
            "scheme is a typed deferral (#310)."
        )

    def streaming_cell_transpose(
        self,
        *,
        res_bar: np.ndarray,
        psi_out_bar: np.ndarray,
        denom: np.ndarray,
        abs_mu_A_total: np.ndarray,
        volume: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""Per-cell VJP of the streamed cell relation — the CURVILINEAR reverse arm's kernel.

        The transpose extension point of the 1-D loop-walk cell relation
        {residual :math:`m = (\text{denom}\,\bar\psi - |\mu|A_{\rm tot}\,
        \psi_{\rm in} - \text{numer})/V`, face chain
        :math:`\psi_{\rm out} = f(\bar\psi, \psi_{\rm in})`}: given the
        residual cotangent ``res_bar`` ``(ng, n)`` and the downstream face
        cotangent ``psi_out_bar`` ``(ng, n)``, return
        ``(psi_bar_cot, psi_in_bar)`` — the cell-average cotangent (the walk
        scatter-adds it; each cell slot is visited by exactly one leg, so the
        single scatter is bit-identical to in-place accumulation) and the
        upstream face cotangent (the reversed march's carry).

        SPATIAL relation ONLY (#310 ruling 1): the Morel–Montry
        angular-numerator cotangent stays with the walk and
        ``PoleAngularClosure.angular_adjoint`` — a kernel folding the angular
        thread would re-fuse the spatial/angular closure separation the
        codebase deliberately keeps.

        STORAGE-FREE by contract, like the batch pair: the walk owns
        gather/scatter and the ψ-independent coefficient recomputation
        (``cell_balance_for_streaming`` — Pattern 2, the SAME coefficients as
        the forward).  Since #310 C2 this kernel serves the 1-D CURVILINEAR
        reverse arm only (single-occupant DD geometry — the mirror of the
        forward's ``cell_balance`` arm); the Cartesian reverse rides the
        scheme-uniform :meth:`residual_kernel_batch_transpose`.  Overriding
        it is the curvilinear conjunct of the derived
        :attr:`has_transpose_kernel` (required iff
        :attr:`supports_curvilinear`).
        """
        raise NotImplementedError(
            f"{type(self).__name__} registers no transpose kernel "
            "(streaming_cell_transpose) — the curvilinear adjoint "
            "matvec/solve on this scheme is a typed deferral (#310)."
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
        reaction_xs: np.ndarray,
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
        (:class:`~orpheus.sn.sweep.cache.CollisionCache`), and ``b``
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

    def cartesian_scan_coefficients(
        self,
        *,
        s_scan: np.ndarray,
        s_transverse: tuple[np.ndarray, ...],
        reaction_xs: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[np.ndarray, ...]]:
        r"""Cartesian row-march scan coefficients ``(a, inverse_denom, w, transverse_couplings)``.

        The multi-D **scan-march** analogue of :meth:`affine_scan_coefficients`
        (which is the 1-D / curvilinear chain producer): given the RAW down-face
        streaming ``g = |μ|/Δ`` on the **scanned** axis (``s_scan``) and on the
        already-known **transverse** axes (``s_transverse``, one entry per
        non-swept Cartesian axis — the 0th-order face traces a row-march
        absorbs into its affine source), plus the local :math:`\Sigma_t`,
        return the four :math:`\Sigma_t`-epoch coefficients that drive the
        in-row first-order x-scan ``ψ_out = a·ψ_in + b``:

        * ``a`` — the scan transmission multiplier (source-INDEPENDENT);
        * ``inverse_denom`` — the reciprocal cell-balance diagonal
          :math:`1/S` with :math:`S = \Sigma_t + (\text{scan diag}) +
          \sum (\text{transverse diag})`;
        * ``w`` — the cell-average blend weight
          (:math:`\bar\psi = (1-w)\psi_{\rm in} + w\,\psi_{\rm out}`);
        * ``transverse_couplings`` — one coefficient per transverse axis, the
          scheme-OWNED coupling on the known upstream transverse face flux that
          enters BOTH the diagonal (above) AND the affine source (the caller
          folds :math:`\sum c_\perp\,\psi_\perp^{\rm in}` into ``QV`` before
          :meth:`source_emission`).  For Diamond Difference each is
          :math:`2 g_\perp` (its diamond :math:`2 = 1/w_{\rm DD}`).

        The caller (``ScanMarch._sweep_interior``) composes the recurrence
        generically via the base staticmethods — folding
        :math:`\sum c_\perp\,\psi_\perp^{\rm in}` into ``QV`` for
        :meth:`source_emission`, then :meth:`cell_average` /
        :meth:`outgoing_face_from_average` — so the diamond ``2`` and the blend
        ``w`` live HERE, never inline in the row body (#239 the 2-D
        coefficient-model lift).

        Only a ``transverse_coupling_is_facewise`` scheme overrides this — the
        separable facewise transverse coupling is exactly the precondition the
        d-D row-march needs (a bilinear LD closure carries its transverse term
        as a slope moment, not a face value, and rides the DAG wavefront
        instead — its ``False`` trait keeps it off the scan-march, and this
        default raise is the belt-and-braces backstop).

        Default raises :exc:`NotImplementedError`.
        """
        raise NotImplementedError(
            f"{type(self).__name__} does not implement "
            "cartesian_scan_coefficients (transverse_coupling_is_facewise is "
            "False).  Only a facewise/separable closure supplies the row-march "
            "scan coefficients consumed by the multi-D ScanMarch; a bilinear "
            "closure rides the DAG wavefront schedule."
        )

    def reflect_scan_coefficients(
        self, psi_bar: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""Apply-direction x-face reflection scan coefficients ``(α, β)``.

        The matvec twin of :meth:`cartesian_scan_coefficients`: in the APPLY
        direction the cell-average ``psi_bar`` is KNOWN (the probe), so the
        downstream x-face is the pure reconstruction
        :math:`\psi_{\rm out} = (\bar\psi - (1-w)\psi_{\rm in})/w =
        \beta + \alpha\,\psi_{\rm in}` — a first-order recurrence in the faces
        with ``α = −(1−w)/w`` and ``β = ψ̄/w`` (the recurrence form of
        :meth:`outgoing_face_from_average`, so the blend ``w`` lives in the
        scheme, never inline).  For Diamond Difference ``w = ½`` gives ``α = −1``
        (a pure reflection, :math:`|α| = 1` so no scan underflow) and
        ``β = 2ψ̄``.

        Only a ``transverse_coupling_is_facewise`` scheme overrides it (the
        scan-march apply twin); the default raises.
        """
        raise NotImplementedError(
            f"{type(self).__name__} does not implement reflect_scan_coefficients "
            "(transverse_coupling_is_facewise is False)."
        )

    def reflect_scan_coefficients_transpose(
        self, res_bar: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""Apply-TRANSPOSE reflection scan coefficients ``(α, β_pullback)``.

        The reverse-direction sibling of :meth:`reflect_scan_coefficients`
        (#310 C4 — the row-march reverse): the forward apply reconstructs the
        scanned-axis faces from the KNOWN probe via ``ψ_out = α·ψ_in + β(ψ̄)``;
        its reverse-mode needs the ψ̄-INDEPENDENT coefficient pair —

        * ``α = −(1−w)/w`` — the transpose of a first-order linear chain is a
          first-order linear chain with the SAME multiplier, run in the
          opposite direction (the ``ordinate_scan_transpose`` chain);
        * ``β_pullback = ∂β/∂ψ̄ = 1/w`` — the diagonal pullback of the
          β-leg onto the probe cotangent (``ψ̄† += β̄·β_pullback``).

        Both are the unit application of the #311 generic VJP
        :meth:`outgoing_face_from_average_transpose` at the scheme's blend
        weight (``(β_pullback, α) = VJP(1, w)``) — the reflect scan IS
        :meth:`outgoing_face_from_average` in recurrence form, so its
        transpose coefficients are that op's, single-sourced.  ``res_bar``
        is a SHAPE template only (the coefficients broadcast against the
        row batch).  Only a ``transverse_coupling_is_facewise`` scheme
        overrides it; the default raises.
        """
        raise NotImplementedError(
            f"{type(self).__name__} does not implement "
            "reflect_scan_coefficients_transpose "
            "(transverse_coupling_is_facewise is False)."
        )

    def moment_scan_closure(
        self,
        *,
        abs_mu: np.ndarray,
        A_down: np.ndarray,
        V: np.ndarray,
        reaction_xs: np.ndarray,
    ) -> "D1ClosedForm":
        r"""Per-cell closed form for the 1-D multi-moment SCAN (#240 D5b-S3).

        The slope-source-aware companion of
        :meth:`affine_scan_coefficients`, consumed by the scan-march's
        multi-moment branch (gated on :attr:`is_multi_moment`) to fold the
        threaded scattering-slope source ``Σ_s·φ̂`` into the face-chain
        affine source and reconstruct the ``(ψ̄, ψ̂)`` cell moments.  Only a
        multi-moment scheme overrides it (LD today; see
        :meth:`~orpheus.transport.spatial.linear_discontinuous.LinearDiscontinuous.moment_scan_closure`);
        the default raises — the same capability-method idiom as
        :meth:`affine_scan_coefficients` / :meth:`reflect_scan_coefficients`.
        """
        raise NotImplementedError(
            f"{type(self).__name__} does not implement moment_scan_closure "
            "(is_multi_moment is False)."
        )

    def _assert_moment_mass_is_expressible(self, coord: CoordSystem) -> None:
        r"""Refuse a moment mass this codebase cannot yet SPELL.

        .. admonition:: ``coord`` is TRANSITIONAL — do not inherit it as a
           decision
           :class: warning

           This method (and the two producers that call it) take the chart as
           a **tag**, and a tag is the wrong thing to take.  The mass is
           :math:`M_{kj} = \\int b_k b_j\\,\\mathrm{d}V`: what it needs is the
           **measure**, and ``coord`` is a proxy the callee has to
           re-interpret.  Discriminating on a tag to recover a measure is the
           missing-type smell ``nexus discriminations`` exists to find.

           The tell that it is a proxy is the predicate below: it is not a
           three-way dispatch but ``is not CARTESIAN`` — i.e. *"is the volume
           element constant across the cell?"*, which is a property of the
           MEASURE's behaviour, not of the chart's identity, and is exactly
           what makes the unit-width shortcut valid on a slab.

           ⭐ The right shape is already derived, and it is a sibling
           relationship rather than a new abstraction::

               M_kj = ∫ b_k b_j            dV     # Galerkin matrix of 1
               R_kj = ∫ b_k b_j (∇·ê_r)    dV     # Galerkin matrix of ∇·ê_r

           the SAME construction against different integrands, so the moment
           mass and the redistribution pairing are two products of ONE
           bilinear form and should be minted together by the object that
           owns the measure — ``scheme.mint(chart)``, §5 of
           ``.claude/plans/streaming_path_says_what_it_is.md``, phase P4.

           ⚠ Why the tag ships in the meantime, stated so the trade is
           visible rather than implied: the guard is a CORRECTNESS repair —
           `[M]` the slab's mass was being installed on the sphere AND the
           cylinder, bit-identical, silently — and it could not exist at all
           while the producer was never told its chart.  One tag-dispatch
           was the smallest thing that makes the wrong value unspellable
           today.  P4 dissolves it; nothing here should be built upon.


        The cell mass is :math:`M_{kj} = \int b_k b_j \, \mathrm{d}V`, so
        it depends on the basis (the scheme's) **and** on the measure (the
        chart's).  On a Cartesian chart :math:`M/V` is width-independent,
        which is why the tensor-Legendre diagonal evaluated at unit width is
        exact there.  On a curvilinear chart it is neither:

        * **non-diagonal** — `[M]` the true :math:`M/V` at a spherical pole
          cell is ``[[1, 0.5], [0.5, 0.4]]``;
        * **cell-dependent** — it varies with :math:`r_i`, i.e. along a
          DIFFERENT axis.

        A per-axis :class:`~orpheus.numerics.axis.Axis` weight vector can
        express neither, so there is no honest value to return — and the
        value that WAS returned is the slab's, silently
        (`[M]` 2026-08-26: sphere and cylinder both installed
        ``[1, 0.3333]``, bit-identical to a slab).  Refusing is the honest
        answer until the machinery exists.

        Blocked on two independent things, deliberately named because they
        need opposite repairs: the **machinery** (a non-Hadamard space
        metric — ORPHEUS #409) and the **value** (which :math:`G` is right
        is pinned by physical functionals, and needs #158's cell solve to
        have a consumer).

        Single-moment schemes are unaffected on every chart: the cell
        average's mass is :math:`V/V = 1` whatever the measure.
        """
        if self.is_multi_moment and coord is not CoordSystem.CARTESIAN:
            raise NotImplementedError(
                f"{type(self).__name__} has no moment mass on a "
                f"{coord.value} chart: the true M/V there is cell-dependent "
                "AND non-diagonal (a spherical pole cell wants "
                "[[1, 0.5], [0.5, 0.4]]), which a per-axis Axis weight "
                "vector cannot express.  Returning the Cartesian diagonal "
                "would install the SLAB's mass on a curved chart.  See "
                "ORPHEUS #409 (the non-Hadamard metric) and #158 (the cell "
                "solve that gives the value a consumer)."
            )

    def moment_mass_diagonal(self, ndim: int, coord: CoordSystem) -> np.ndarray:
        r"""Unit-volume Legendre mass diagonal over the ``2^d`` cell moments.

        The per-moment mass weight ``M_ii/V = ∏_a θ^{o_a}`` of the scheme's
        tensor-Legendre cell basis — ``(1,)`` for the slopeless closures
        (the trivial average-only mass), ``[1, θ]`` for LD at d=1,
        ``[1, θ, θ, θ²]`` at d=2.  The Hilbert-metric consumer:
        the bulk measure of a moment-tailed field is
        ``G_bulk = V·w_n ⊗ moment_mass_diagonal`` (#310 C2 ruling 3 — an
        average-only ``V·w_n`` mis-weights the slope DOF, making ``.H`` a
        WRONG adjoint on the slope rows AND reciprocity Mode-12-blind to a
        slope-row transpose).  Single-sourced from the scheme because the
        mass IS the scheme's (the same diagonal
        :meth:`residual_kernel_batch` normalises by, ÷V).
        """
        self._assert_moment_mass_is_expressible(coord)
        return np.ones(cell_moment_count(self.spatial_basis_per_axis, ndim))

    def moment_axis(self, ndim: int, coord: CoordSystem) -> "Axis":
        r"""The within-cell spatial-moment factor as a typed space AXIS.

        The :class:`~orpheus.numerics.axis.Axis` form of this scheme's
        tensor-Legendre cell basis — label ``"spatial_moment"``, shape
        ``(2^d,)``, ``weights = ``:meth:`moment_mass_diagonal`, kind
        **MODAL** (expansion coefficients: a positive function has legally
        signed slope coefficients, so per-component positivity is not a
        meaningful predicate — ``has_coordinate_cone`` reads ``False``
        through this factor, which is the honest answer campaign 1 CS4b
        routes :meth:`~orpheus.numerics.field.Field.cone_violations`
        through).

        Scheme-owned because the mass IS the scheme's (``θ`` enters
        :meth:`moment_mass_diagonal` — basis ↔ mass single-sources here);
        every axis-built bulk space that carries the moment tail composes
        THIS object, so the field layer and the composite interior cannot
        spell the mass twice (CS4b crosswalk B5,
        ``.claude/plans/cs4b_crosswalk.md``).

        **Refuses on a slopeless closure** (``is_multi_moment`` is False —
        the same gate :meth:`moment_scan_closure` keys on): the moment-tail
        policy appends NO factor at width 1 (the byte-identity invariant,
        #240 D5b-S3-A0), so a trivial ``(1,)`` moment axis is a state no
        composed space admits — minting one would invite it (Pattern 4).
        """
        if not self.is_multi_moment:
            raise NotImplementedError(
                f"{type(self).__name__} has no moment axis "
                "(is_multi_moment is False): the moment-tail policy "
                "appends no factor at width 1, so there is no moment axis "
                "to mint."
            )
        from orpheus.numerics.axis import Axis, BasisKind
        from orpheus.numerics.moment_layout import SPATIAL_MOMENT_AXIS_LABEL

        return Axis(
            SPATIAL_MOMENT_AXIS_LABEL,
            (cell_moment_count(self.spatial_basis_per_axis, ndim),),
            weights=self.moment_mass_diagonal(ndim, coord),
            kind=BasisKind.MODAL,
        )

    # ── Generic affine reconstruction ops (#158 Inc B / #240 D2) ─────────────
    #
    # Pure functions of the face fluxes / source and the blend weight ``w`` —
    # no instance state, no scheme-specific constant.  The GENERIC
    # advection–reaction reconstruction shared by every affine scheme (and
    # inherited by the future diffusion scheme).  See the class docstring
    # §"Generic affine reconstruction ops" for the universality argument and
    # the bit-identity-vs-principled-equivalence note.

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
        reconstruction: DD's diamond mean ``2ψ̄ − ψ_in`` is the ``w=½`` case
        (byte-identical: ``÷½`` is an exact power-of-2 ``×2``); LD's
        ``(1+k)ψ̄ − k·ψ_in`` is the ``w=1/(1+k)`` case (a principled ~1-ULP
        re-baseline off a different reduction tree — see the class docstring
        §"Generic affine reconstruction ops").
        """
        return (psi_bar - (1.0 - w) * face_in) / w

    @staticmethod
    def outgoing_face_from_average_transpose(
        face_out_bar: np.ndarray, w: "float | np.ndarray",
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""VJP of :meth:`outgoing_face_from_average`: ``(ψ̄†, ψ_in†) = (f̄/w, −((1−w)/w)·f̄)``.

        The w-generic reverse-mode pair of the affine outflow reconstruction —
        each returned cotangent is the partial derivative times the
        outgoing-face cotangent ``f̄``: ``∂ψ_out/∂ψ̄ = 1/w`` and
        ``∂ψ_out/∂ψ_in = −(1−w)/w``.  The single source of the reverse walks'
        face-chain pair (#311): DD's hand-transposed diamond chain
        ``(2f̄, −f̄)`` is the ``w=½`` case (byte-identical — ``÷½`` is an
        exact power-of-2 ``×2`` and ``(1−½)/½ = 1`` exactly, so ``×(−1)`` is
        the exact sign flip); LD's ``((1+k)·f̄, −k·f̄)`` is the ``w=1/(1+k)``
        case (#310 C2).
        """
        avg_cot = face_out_bar / w
        in_cot = (-(1.0 - w) / w) * face_out_bar
        return avg_cot, in_cot


__all__ = [
    "CellResult",
    "DiscretizationScheme",
    "DiscretizationSchemeBase",
    "CellVisit",
    "UpstreamState",
]
