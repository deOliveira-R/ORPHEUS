r"""Unified S\ :sub:`N` transport sweep — fold over the per-ordinate DAG.

Issue #196 Phase G Step 2.5 collapses three 1-D sweep bodies
(:func:`_sweep_1d_cumprod`, :func:`_sweep_1d_spherical`,
:func:`_sweep_1d_cylindrical`) into ONE skeleton that folds over the
DAG-ordered cell visits emitted by
:meth:`~orpheus.sn.geometry.SNMesh.dag_walk`.  The per-cell algebra
lives in ``sn_mesh.cell_update.update(visit, total_xs, source,
upstream)`` (Wave C strategy contract; the Step-2.5 unified
:class:`DiamondDifference` body in
:mod:`orpheus.sn.spatial.diamond`).  Geometry is data carried by
:class:`StreamingTerms` / :class:`CellVisit`; the sweep skeleton has
no internal geometry dispatch.

The structural claim
====================

The SN sweep is **forward substitution** on the block-triangular
(streaming + collision) operator under the DAG ordering induced by
the ordinate's direction.  Forward substitution is a fold over a
topologically-ordered work-unit stream.  This module makes that
visible:

.. code-block:: python

    psi_face_in = bc.inflow(direction_sign)
    for visit in sn_mesh.dag_walk(direction_sign, mu_level_idx=p):
        upstream = UpstreamState(spatial_upstream=psi_face_in,
                                 angular_upstream=psi_angle[visit.cell_idx])
        result = cell_update.update(visit, total_xs[i], source[i], upstream)
        psi_avg[i] = result.cell_average_flux
        if result.outgoing_spatial_flux is not None:
            psi_face_in = result.outgoing_spatial_flux
        if result.outgoing_angular_state is not None:
            psi_angle[i] = result.outgoing_angular_state

The fold accumulator is ``(psi_face_in, psi_angle[:])`` — one scalar
per direction (spatial face flux threaded across cells along the
DAG) and one per-cell array (angular state threaded across ordinates
within the μ-level).  For slab, ``angular_upstream`` is ``None``
(no angular redistribution) and the angular loop is trivial.  For
sphere, the μ-level loop is trivial (single level).  For cylinder,
both loops run.

The 2-D wavefront sweep (:func:`_sweep_2d_wavefront`) retains its
own batched anti-diagonal scheduling — that is a fold over LEVELS,
each level a :class:`SweepCellSlice` consumed by
:meth:`CellUpdate.update_batch`.  The 1-D and 2-D folds are the
same operator-algebra concept (forward substitution on block-
triangular L+C); the only difference is the work-unit shape
(per-cell vs per-level).

What this module retired (Step 2.5)
===================================

* ``_sweep_1d_cumprod`` (slab) — replaced by the unified fold over
  :meth:`SNMesh.dag_walk(direction_sign)`.  The cumprod operation
  order was migration-phase scaffolding for the pre-Wave-C bit-
  identity contract; Phase G is the natural endpoint.
* ``_sweep_1d_spherical`` and ``_sweep_1d_cylindrical`` — merged
  into :func:`_sweep_1d_curvilinear` which handles both via
  ``sn_mesh.reduced.coord``.

The slab Cartesian regression snapshots may drift at IEEE-754 ULP
due to the operation-order change (cumprod vectorised across cells
vs per-cell fold) — bounded by ``iteration_count × cell_count ×
ULP`` and well within the existing ``rtol=1e-12`` regression
contract.  Per ``vv-principles`` §"Bit-identity vs principled-
equivalence" all three criteria hold (principled at every step,
structurally-independent reference via Phase E Variant α
attestation, drift bounded by FP-non-associativity).

References
==========

* Hébert, A. (2009). *Applied Reactor Physics*.  Ch. 3 §3.9.4
  (pp. 141-144) — primary source for the curvilinear S\ :sub:`N`
  cell-balance + DD difference relations the sweep implements.
* Bailey, T. S., Morel, J. E., & Chang, J. H. (2010).  *Asymptotic
  Diffusion-Limit Accuracy of Sn Angular Differencing Schemes*.
  NSE 165(2):149-169.  Auxiliary justification for the M-M clamp.
* Lewis, E. E., & Miller, W. F. (1984).  *Computational Methods of
  Neutron Transport.*  §4.5 (curvilinear DD); §5.3 (DD, weighted-DD,
  Step, LD); §6.4 (sweep ordering).

See also
========

* :class:`~orpheus.sn.spatial.cell_update.CellUpdate` — the
  per-cell strategy contract.
* :class:`~orpheus.sn.spatial.diamond.DiamondDifference` — the
  default cell-update strategy (Step 2.5 unified body).
* :meth:`~orpheus.sn.geometry.SNMesh.dag_walk` — the unified
  iteration primitive.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from orpheus.geometry import CoordSystem

from .spatial.cell_update import UpstreamState
from .spatial.diamond import DiamondDifference
from .spatial.psi_half_angle_seed import carlson_inward_sweep_from_source
from .sweep_graph import OctantLabel

if TYPE_CHECKING:
    from .geometry import SNMesh


# ═══════════════════════════════════════════════════════════════════════
# Top-level dispatch
# ═══════════════════════════════════════════════════════════════════════

def transport_sweep(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: SNMesh,
    psi_bc: dict,
    Q_aniso: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Perform one full transport sweep.

    Boundary conditions are read from ``sn_mesh`` (resolved at
    construction time from the geometry mesh's :class:`BC` declarations).
    The cell-update strategy is read from ``sn_mesh.cell_update``
    (defaults to :class:`DiamondDifference`).

    Dispatch:

    * ``reduced.requires_upstream_angular_state == True`` →
      :func:`_sweep_1d_curvilinear` (sphere / cylinder; the per-
      μ-level loop subsumes both via ``reduced.coord``).
    * 2-D Cartesian → :func:`_sweep_2d_wavefront`.
    * 1-D Cartesian → :func:`_sweep_1d_cartesian`.

    Parameters
    ----------
    Q : (nx, ny, ng) isotropic source density.
    sig_t : (nx, ny, ng) total macroscopic cross section.
    sn_mesh : SNMesh — augmented geometry with precomputed stencil,
        resolved boundary conditions, and a ``cell_update`` strategy.
    psi_bc : mutable dict storing persistent boundary fluxes
        between outer iterations.
    Q_aniso : (N, nx, ny, ng) per-ordinate anisotropic source (P1+
        scattering).  ``None`` for isotropic-only (P0).

    Returns
    -------
    angular_flux : (N, nx, ny, ng) angular flux per ordinate.
    scalar_flux : (nx, ny, ng) = Σ_n w_n ψ_n.
    """
    reduced = sn_mesh.reduced
    if reduced is not None and reduced.requires_upstream_angular_state:
        return _sweep_1d_curvilinear(Q, sig_t, sn_mesh, psi_bc, Q_aniso)
    # 1-D Cartesian sweep requires ``reduced`` populated (slab_streaming).
    # ``Mesh2D`` with ``ny == 1`` is genuinely 2-D (no ``reduced``) and
    # must take the wavefront path — its ``iter_cell_visits`` is unavailable.
    if reduced is not None:
        return _sweep_1d_cartesian(Q, sig_t, sn_mesh, psi_bc, Q_aniso)
    return _sweep_2d_wavefront(Q, sig_t, sn_mesh, psi_bc, Q_aniso)


# ═══════════════════════════════════════════════════════════════════════
# 1-D Cartesian sweep — unified fold (replaces _sweep_1d_cumprod)
# ═══════════════════════════════════════════════════════════════════════

def _sweep_1d_cartesian(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: SNMesh,
    psi_bc: dict,
    Q_aniso: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Slab 1-D sweep — fold over :meth:`SNMesh.dag_walk` (Step 2.5).

    Replaces ``_sweep_1d_cumprod`` (retired in Issue #196 Phase G
    Step 2.5) with the same fold-over-DAG-visits skeleton that the
    curvilinear sweep uses.  Slab has no angular redistribution
    (``upstream_state.angular_upstream = None``) so the M-M closure
    in :class:`DiamondDifference` returns ``outgoing_angular_state =
    None`` and the angular-state thread does not run.

    Bit-identity caveat
    -------------------

    The pre-Step-2.5 cumprod path vectorised the per-cell DD
    recurrence as a single cumulative product across all cells in a
    row.  The Step-2.5 fold computes one cell at a time via
    :meth:`DiamondDifference.update`, which uses the unified
    ``(source + numer_upstream)/denom`` operation order (not the
    cumprod's ``a·ψ_in + 2q/denom`` order).  The two are
    algebraically identical but ULP-different at IEEE-754; slab
    regression snapshots may drift by ``iteration_count × ULP``
    (well within the existing ``rtol=1e-12`` contract).
    """
    quad = sn_mesh.quad
    N = quad.N
    nx = sn_mesh.nx
    ng = Q.shape[2]
    weights = quad.weights
    mu = quad.mu_x

    Q_1d = Q[:, 0, :]          # (nx, ng)
    sig_t_1d = sig_t[:, 0, :]  # (nx, ng)
    V = sn_mesh.volumes[:, 0]  # (nx,) — for slab, V == dx
    cell_update = sn_mesh.cell_update

    weight_norm = 1.0 / weights.sum()
    QV_iso = Q_1d * V[:, None] * weight_norm  # (nx, ng)

    bc_left_obj = sn_mesh.bc_left
    bc_right_obj = sn_mesh.bc_right

    # Persistent per-ordinate face-flux buffers (carry reflective BC
    # state across iterations).  Stored at the full N axis so the BC
    # operator's ``apply`` can read partner ordinates.
    if "bc_1d_left_face" not in psi_bc:
        psi_bc["bc_1d_left_face"] = np.zeros((N, ng))
        psi_bc["bc_1d_right_face"] = np.zeros((N, ng))
    bc_left_face = psi_bc["bc_1d_left_face"]   # (N, ng) — outgoing at x=0
    bc_right_face = psi_bc["bc_1d_right_face"]  # (N, ng) — outgoing at x=L

    angular_flux = np.zeros((N, nx, 1, ng))
    scalar_flux = np.zeros((nx, ng))

    has_aniso = Q_aniso is not None
    if has_aniso:
        Q_aniso_1d = Q_aniso[:, :, 0, :] * weight_norm  # (N, nx, ng)

    # BC inflow at both faces (resolved once per sweep call).  For
    # vacuum: zeros.  For specular reflective: the previous-sweep
    # outgoing flux at the partner ordinate.  Linear realisation —
    # commutes with per-ordinate extraction.
    inflow_left = bc_left_obj.apply(bc_left_face)    # (N, ng)
    inflow_right = bc_right_obj.apply(bc_right_face)  # (N, ng)

    for n in range(N):
        mu_n = mu[n]
        w_n = weights[n]
        direction_sign = +1 if mu_n >= 0 else -1

        QV = QV_iso
        if has_aniso:
            QV = QV_iso + Q_aniso_1d[n] * V[:, None]

        # Initial spatial-upstream face flux for this ordinate's
        # sweep.  Outward: read left BC.  Inward: read right BC.
        if direction_sign == +1:
            psi_spatial_in = inflow_left[n]
        else:
            psi_spatial_in = inflow_right[n]

        # Fold over the DAG-ordered cell visits.  Slab has no angular
        # redistribution, so we pass ``angular_upstream = None``.
        # ``iter_cell_visits(ordinate_idx=n)`` (NOT ``dag_walk``) is
        # required because slab DD reads per-ordinate quantities from
        # ``streaming_terms`` (``abs_mu``, ``mu``) that depend on the
        # specific ordinate, not just the direction sign.
        for visit in sn_mesh.iter_cell_visits(ordinate_idx=n):
            i = visit.cell_idx
            upstream = UpstreamState(
                spatial_upstream=psi_spatial_in,
                angular_upstream=None,
            )
            result = cell_update.update(
                visit=visit,
                total_xs=sig_t_1d[i],
                source=QV[i],
                upstream_state=upstream,
            )
            psi = result.cell_average_flux
            angular_flux[n, i, 0, :] = psi
            scalar_flux[i] += w_n * psi
            # Spatial closure always populated for slab
            # (visit.face_area_downstream == 1.0).
            psi_spatial_in = result.outgoing_spatial_flux

        # Store the outgoing flux at the downstream boundary face for
        # reflective BCs to read on the next iteration.
        if direction_sign == +1:
            bc_right_face[n] = psi_spatial_in
        else:
            bc_left_face[n] = psi_spatial_in

    return angular_flux, scalar_flux[:, None, :]


# ═══════════════════════════════════════════════════════════════════════
# 1-D Curvilinear sweep — unified fold over μ-levels (sphere + cylinder)
# ═══════════════════════════════════════════════════════════════════════

def _sweep_1d_curvilinear(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: SNMesh,
    psi_bc: dict,
    Q_aniso: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    r"""1-D curvilinear sweep — fold over μ-levels and DAG visits.

    Issue #196 Phase G Step 2.5: replaces ``_sweep_1d_spherical``
    and ``_sweep_1d_cylindrical`` with one body that handles both
    via ``sn_mesh.reduced.coord``.  Sphere is the single-level case
    (``mu_level_idx = None``); cylinder loops over
    ``quad.level_indices``.

    Per-level structure
    -------------------

    For each μ-level :math:`p`:

    1. Build the Carlson coupled-pole half-angle seed for the level
       (Hébert §3.9.4 Eqs. 3.432-3.435) from the previous-iteration
       scalar flux φ_0 (or cold-start from the source Q).  Replaces
       the legacy ``psi_angle = 0`` seed (ERR-026 manifestation #6
       fix; Phase F).
    2. Iterate ordinates within the level (sphere: global ordinates;
       cylinder: within-level azimuthal ordinates).  For each
       ordinate, fold over :meth:`SNMesh.dag_walk` to walk the
       cells in topological order.
    3. Within the fold, the cell update reads the upstream angular
       state from ``psi_angle[i]`` and writes the downstream into
       the same slot for the NEXT ordinate's fold to consume — the
       M-M angular DAG threads orthogonally to the spatial DAG.

    Pole-face initial condition (Phase G Step 2 Path C)
    ---------------------------------------------------

    On outward sweeps the pole-face ψ inflow is set to the cell-
    centre value from the previous outer iteration (Lewis-Miller
    §4.5 Carlson starting-direction convention).  This makes the SI
    Picard fixed point coincide with the apply-matvec fixed point on
    the homogeneous reflective sphere streaming-equilibrium test.
    Cold-start: pre-iter-1 there is no ψ_pole cache; fall back to
    the legacy zero IC.

    Coord switch
    ------------

    The only coord-dependent decisions are:

    * **Level enumeration**: sphere = single level
      (``mu_level_idx=None``); cylinder = ``quad.level_indices``.
    * **psi_bc keys**: separate keys per geometry to avoid
      collisions in mixed test fixtures.
    """
    nx = sn_mesh.nx
    ng = Q.shape[2]
    quad = sn_mesh.quad
    N = quad.N
    mu = quad.mu_x
    weights = quad.weights

    Q_1d = Q[:, 0, :]          # (nx, ng)
    sig_t_1d = sig_t[:, 0, :]  # (nx, ng)
    V = sn_mesh.volumes[:, 0]  # (nx,) cell volumes
    cell_update = sn_mesh.cell_update

    coord = sn_mesh.reduced.coord
    is_sphere = coord is CoordSystem.SPHERICAL

    # Persistent boundary-flux buffer (per ordinate, indexed by
    # global ordinate).
    bc_key = "bc_sph" if is_sphere else "bc_cyl"
    pole_key = "psi_pole" if is_sphere else "psi_pole_cyl"
    phi_prev_key = "phi_0_prev" if is_sphere else "phi_0_prev_cyl"

    bc_outer_obj = sn_mesh.bc_right
    if bc_key not in psi_bc:
        psi_bc[bc_key] = np.zeros((N, ng))
    bc_outer = psi_bc[bc_key]

    # Resolved outer-face inflow (single BC apply; linear, commutes
    # with per-level / per-ordinate extraction).
    inflow_full = bc_outer_obj.apply(bc_outer)  # (N, ng)

    angular_flux = np.zeros((N, nx, 1, ng))
    scalar_flux = np.zeros((nx, ng))

    # Source pre-scale.
    weight_norm = 1.0 / weights.sum()
    QV_iso = Q_1d * V[:, None] * weight_norm  # (nx, ng)

    has_aniso = Q_aniso is not None
    if has_aniso:
        Q_aniso_1d = Q_aniso[:, :, 0, :] * weight_norm  # (N, nx, ng)

    # Carlson seed source (Phase G Step 2 Path C):
    # Q̄ = Σ_t · φ_0_prev / Σw_full.  The geometry-general
    # normalisation (Σw_full = quad.weights.sum()) replaces the
    # legacy sphere-only ``0.5`` literal.  Cold-start: no
    # phi_0_prev cache → fall back to Q_1d / Σw (P_0(−1)=1).
    sigma_t_gx = sig_t_1d.T  # (ng, nx)
    dr = sn_mesh.dx
    Sw_full = weights.sum()
    if phi_prev_key in psi_bc:
        phi_0_prev = psi_bc[phi_prev_key]
        Q_bar_iso = sigma_t_gx * phi_0_prev.T / Sw_full
    else:
        Q_bar_iso = Q_1d.T / Sw_full

    # Level enumeration: sphere = single virtual level; cylinder =
    # quad.level_indices.
    if is_sphere:
        levels = [None]   # mu_level_idx=None for sphere
        # Sphere's "level" contains every global ordinate; ordinate
        # index is the global ordinate itself.
        level_ordinate_lists: list[list[int]] = [list(range(N))]
    else:
        level_indices = quad.level_indices
        levels = list(range(len(level_indices)))
        level_ordinate_lists = [list(level_indices[p]) for p in levels]

    for p_idx, level in enumerate(levels):
        ordinates_in_level = level_ordinate_lists[p_idx]

        # Carlson coupled-pole seed for this level's angular DAG.
        # ``bc_outer_value``: outer-face inflow at the most-inward
        # ordinate of this level (μ = −1 for sphere; η = −sin θ for
        # cylinder).  Derived from the BC-realised outflow via the
        # already-computed ``inflow_full``.
        ords_arr = np.asarray(ordinates_in_level)
        mu_in_level = mu[ords_arr]
        most_inward_global = int(ords_arr[int(np.argmin(mu_in_level))])
        bc_outer_value = inflow_full[most_inward_global, :]  # (ng,)

        phi_aux = carlson_inward_sweep_from_source(
            Q_bar=Q_bar_iso,
            sigma_t=sigma_t_gx,
            dr=dr,
            bc_outer_value=bc_outer_value,
        )  # (ng, nx)
        psi_angle = phi_aux.T.copy()  # (nx, ng)

        # Iterate ordinates within this level.  For sphere, the
        # within-level index IS the global ordinate index.  For
        # cylinder, we resolve the global ordinate via
        # ``level_indices[p]``.
        for m_local, global_n in enumerate(ordinates_in_level):
            mu_n = mu[global_n]
            w_n = weights[global_n]

            QV = QV_iso
            if has_aniso:
                QV = QV_iso + Q_aniso_1d[global_n] * V[:, None]

            # Spatial-upstream face flux at the entry cell.
            if mu_n < 0:
                # Inward: read outer-face inflow.
                psi_spatial_in = inflow_full[global_n]
            else:
                # Outward (or degenerate): pole-face IC at i=0.
                # Phase G Step 2 Path C: cell-centre ψ from the
                # previous outer iteration (Carlson starting
                # direction).  Cold-start: zeros.
                if pole_key in psi_bc:
                    psi_spatial_in = psi_bc[pole_key][global_n]
                else:
                    psi_spatial_in = np.zeros(ng)

            # Fold over the DAG-ordered cell visits.  We use
            # ``iter_cell_visits(ordinate_idx=...)`` (NOT
            # ``dag_walk``) because the per-cell ``streaming_terms``
            # bundle is direction-IDX-specific: ``alpha_in``,
            # ``alpha_out``, ``tau_mm`` are per-ordinate (sphere) or
            # per-(level, within-level-idx) (cylinder), and ``abs_mu``
            # depends on the exact ordinate.  ``dag_walk`` picks a
            # representative ordinate for the direction sign and is
            # appropriate only when the consumer recomputes its own
            # per-ordinate quantities (the sweep-frame matvec at
            # operator.py).
            #
            # For sphere: ordinate_idx == global ordinate index;
            # mu_level_idx is None.
            # For cylinder: ordinate_idx is the within-level index;
            # mu_level_idx selects the level.
            mu_level_idx = level   # None for sphere; p for cylinder
            ordinate_idx = global_n if is_sphere else m_local
            for visit in sn_mesh.iter_cell_visits(
                ordinate_idx=ordinate_idx,
                mu_level_idx=mu_level_idx,
            ):
                i = visit.cell_idx
                upstream = UpstreamState(
                    spatial_upstream=psi_spatial_in,
                    angular_upstream=psi_angle[i],
                )
                result = cell_update.update(
                    visit=visit,
                    total_xs=sig_t_1d[i],
                    source=QV[i],
                    upstream_state=upstream,
                )
                psi = result.cell_average_flux
                # M-M angular state threads orthogonally across
                # ordinates within the level.
                psi_angle[i] = result.outgoing_angular_state
                angular_flux[global_n, i, 0, :] = psi
                scalar_flux[i] += w_n * psi
                # Spatial closure: ``None`` for cylindrical pure-
                # azimuthal degenerate (visit.face_area_downstream
                # == 0.0); skip the face-flux thread on that ordinate.
                if result.outgoing_spatial_flux is not None:
                    psi_spatial_in = result.outgoing_spatial_flux

            # Store outflow at the outer boundary for reflective BC
            # — only on non-degenerate outward sweeps.  Inward sweeps
            # don't write to the outer boundary (they READ from it);
            # degenerate ordinates have no spatial face flow.
            abs_mu_n = abs(mu_n)
            if mu_n >= 0 and abs_mu_n >= sn_mesh._DEGENERATE_ABS_ETA_THRESHOLD:
                bc_outer[global_n] = psi_spatial_in

    # Cache previous-iter state for next sweep's Carlson seed source
    # + pole-face IC fixes.  Phase G Step 2 Path C.
    psi_bc[pole_key] = angular_flux[:, 0, 0, :].copy()
    psi_bc[phi_prev_key] = scalar_flux.copy()

    return angular_flux, scalar_flux[:, None, :]


# ═══════════════════════════════════════════════════════════════════════
# 2D wavefront path (vectorized along anti-diagonals; inlined DD math)
# ═══════════════════════════════════════════════════════════════════════
#
# The 2-D wavefront sweep retains its inlined DD math.  Parameterising
# anti-diagonal-vectorised cell updates by the
# :class:`~orpheus.sn.spatial.cell_update.CellUpdate` Protocol while
# preserving 1-ULP equality is the open architectural problem for
# Wave C-extension's LD / EC / Step rollout — the per-cell call shape
# would need to accept ``(n_diag, ng)`` slice arguments rather than
# the contract's ``(ng,)`` per-cell shape.  Wave D's gating contract
# is bit-identity for ``cell_update = DiamondDifference``; the inlined
# DD math here meets that contract by construction.

def _sweep_2d_wavefront(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: SNMesh,
    psi_bc: dict,
    Q_aniso: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    r"""2-D wavefront sweep — per-octant batched (Wave 2 / C2.6).

    Iterates over angular **octants** (lexicographic order from
    :attr:`AngularQuadrature.octants`). For each in-plane octant
    :math:`\sigma = (\mathrm{sign}\,\mu_x, \mathrm{sign}\,\mu_y)`:

    1. **BC apply once** on the octant-incoming face(s) — replaces
       the per-ordinate ``bc.apply(...)[n]`` calls of the legacy
       implementation (saves ``N`` invocations per sweep).
    2. **Dispatch to ``SNMesh.sweep_graphs[OctantLabel(σ)]``** —
       the per-octant ``SweepDependencyGraph`` precomputed once at
       mesh construction (Wave 2 / C2.4).
    3. The graph's ``apply`` walks topological levels (anti-diagonals)
       and dispatches each level to
       ``sn_mesh.cell_update.update_batch(slice_args)`` — vectorised
       over ``(N_oct, n_diag, ng)`` simultaneously.

    The smoking gun ``for n in range(N)`` is gone: the outer loop is
    now ``for octant in quad.octants`` (4-8 iterations, structural),
    and within each octant the work is one BC apply + one
    ``SweepDependencyGraph.apply`` call. The ordinate axis is
    INTERNAL to every numpy operation.

    Bit-identity to legacy
    ----------------------

    For LS-family quadratures whose ordinate ordering is
    octant-grouped in lexicographic order (``LevelSymmetricSN``,
    ``ProductQuadrature``), this implementation is bit-identical to
    the legacy per-ordinate loop on every snapshot — verified by
    ``tests/sn/test_2d_octant_sweep_equivalence.py`` (the C2.5
    TESTS-FIRST harness). The argument:

    * BC apply for octant σ reads partners (in the x-reflected
      octant). For LS4, lex order matches legacy n-order at the
      partner-state granularity, so the same iteration's value is
      observed at the same point.
    * Within an octant, the per-ordinate cell sweeps are independent
      (different rows of ``psi_x`` / ``psi_y``); batching is
      bit-identical to any per-ordinate sequencing of the same set.

    For quadratures whose ordering is NOT lex (e.g. Lebedev), the
    converged answer is the same Gauss-Seidel fixed point but the
    iter-to-iter values may differ. Per ``vv-principles`` Bit-identity
    vs principled-equivalence: principled at every step (octant-batched
    BC apply IS the §15A.2 form), structurally-independent grounding
    (closed-form L1 anchor + MMS regression suite both pass),
    convergence to the same answer (verified empirically on regression
    snapshots).
    """
    nx, ny, ng = Q.shape
    quad = sn_mesh.quad
    N = quad.N
    weights = quad.weights

    angular_flux = np.zeros((N, nx, ny, ng))
    scalar_flux = np.zeros((nx, ny, ng))

    # Persistent boundary-flux buffers for reflective BCs (mutated in
    # place across sweep calls — partner reads need the previous
    # iteration's outgoing-face writes).
    if "bc_2d_x" not in psi_bc:
        psi_bc["bc_2d_x"] = np.zeros((N, nx + 1, ny, ng))
        psi_bc["bc_2d_y"] = np.zeros((N, nx, ny + 1, ng))
    psi_x = psi_bc["bc_2d_x"]   # (N, nx+1, ny, ng)
    psi_y = psi_bc["bc_2d_y"]   # (N, nx, ny+1, ng)

    # Source pre-scale once (was inside the ordinate loop in legacy
    # for a no-op O(1) load — kept here for clarity).
    weight_norm = 1.0 / weights.sum()
    Q_scaled = Q * weight_norm
    has_aniso = Q_aniso is not None
    if has_aniso:
        Q_aniso_scaled = Q_aniso * weight_norm  # (N, nx, ny, ng)

    str_x = sn_mesh.streaming_x   # (N, nx)
    str_y = sn_mesh.streaming_y   # (N, ny)
    cell_update = sn_mesh.cell_update

    for octant in quad.octants:
        label_tuple = octant.label   # e.g. (-1, +1, +1) or (-1,) for 1-D
        oct_idx = octant.indices     # (N_oct,) int into N
        N_oct = oct_idx.size
        # In-plane signs — drop sign_z if the label is 3-D.
        sx = label_tuple[0] if len(label_tuple) >= 1 else +1
        sy = label_tuple[1] if len(label_tuple) >= 2 else 0
        # Pure-z degenerate octant: no in-plane streaming. The
        # angular flux is the volumetric balance ψ = Q_n / Σ_t and
        # the scalar flux gets a weighted contribution.
        if sx == 0 and sy == 0:
            Q_pure_z = Q_scaled[None, :, :, :]    # (1, nx, ny, ng)
            if has_aniso:
                Q_pure_z = Q_pure_z + Q_aniso_scaled[oct_idx]
            psi_avg_pure_z = Q_pure_z / sig_t     # broadcasts to (N_oct, nx, ny, ng)
            angular_flux[oct_idx] = psi_avg_pure_z
            scalar_flux += np.einsum(
                "nijg,n->ijg", psi_avg_pure_z, weights[oct_idx],
            )
            continue

        # Effective in-plane sign for sweep-graph lookup. Match
        # legacy's ``key = (1 if mx >= 0 else -1, ...)`` mapping:
        # ordinates with ``mx == 0`` are treated as ``+1`` (the BC
        # apply uses xmin, the streaming coefficient is zero, and
        # the WDD result is identical regardless of sign choice).
        sx_eff = +1 if sx == 0 else sx
        sy_eff = +1 if sy == 0 else sy
        sweep_graph = sn_mesh.sweep_graphs[OctantLabel(sx_eff, sy_eff)]

        # ── BC apply once per octant ───────────────────────────────
        #
        # The octant's incoming-x face is index 0 if sx >= 0 else nx.
        # Boundary operators take the FULL (N, ny, ng) face buffer
        # (so reflective BCs can read partner-octant rows) and return
        # an updated full buffer; we scatter only this octant's rows
        # back into psi_x / psi_y.
        if sx_eff >= 0:
            full_face_x = sn_mesh.bc_xmin.apply(psi_x[:, 0, :, :])
            psi_x[oct_idx, 0, :, :] = full_face_x[oct_idx]
        else:
            full_face_x = sn_mesh.bc_xmax.apply(psi_x[:, nx, :, :])
            psi_x[oct_idx, nx, :, :] = full_face_x[oct_idx]

        if sy_eff >= 0:
            full_face_y = sn_mesh.bc_ymin.apply(psi_y[:, :, 0, :])
            psi_y[oct_idx, :, 0, :] = full_face_y[oct_idx]
        else:
            full_face_y = sn_mesh.bc_ymax.apply(psi_y[:, :, ny, :])
            psi_y[oct_idx, :, ny, :] = full_face_y[oct_idx]

        # ── Per-octant buffers for the graph apply ────────────────
        #
        # Fancy indexing creates copies (not views), so we extract
        # mutable per-octant buffers, run the graph, and scatter back.
        # The cost (~2 × N_oct × nx × ny × ng × 8 bytes per octant)
        # is negligible vs the loop savings.
        psi_x_oct = psi_x[oct_idx].copy()    # (N_oct, nx+1, ny, ng)
        psi_y_oct = psi_y[oct_idx].copy()    # (N_oct, nx, ny+1, ng)
        Q_octant = Q_scaled[None, :, :, :]   # (1, nx, ny, ng) broadcastable
        if has_aniso:
            Q_octant = Q_octant + Q_aniso_scaled[oct_idx]   # (N_oct, nx, ny, ng)
        angular_flux_oct = np.zeros((N_oct, nx, ny, ng))

        sweep_graph.apply(
            cell_update=cell_update,
            psi_x_octant=psi_x_oct,
            psi_y_octant=psi_y_oct,
            Q_octant=Q_octant,
            sig_t=sig_t,
            str_x_octant=str_x[oct_idx],
            str_y_octant=str_y[oct_idx],
            weights_octant=weights[oct_idx],
            angular_flux_octant=angular_flux_oct,
            scalar_flux_buf=scalar_flux,
        )

        # Scatter back the persistent BC buffers + angular flux.
        psi_x[oct_idx] = psi_x_oct
        psi_y[oct_idx] = psi_y_oct
        angular_flux[oct_idx] = angular_flux_oct

    return angular_flux, scalar_flux
