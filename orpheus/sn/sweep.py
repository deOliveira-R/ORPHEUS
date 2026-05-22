r"""Unified S\ :sub:`N` transport sweep — fold over the per-ordinate DAG.

Issue #196 Phase G Step 2.5c collapses the two 1-D sweep bodies
(``_sweep_1d_cartesian`` and ``_sweep_1d_curvilinear``) into ONE
:func:`_sweep_1d_unified` that consumes the two-stratum precomputed
cache (:class:`~orpheus.sn.spatial.sweep_cache.GeometryCoefficients`
+ :class:`~orpheus.sn.spatial.sweep_cache.CollisionCache`).  The
per-cell algebra lives in
:class:`~orpheus.sn.spatial.diamond.DiamondDifference`'s ``update``
method (slow per-cell reference; used for degenerate cylindrical
ordinates and the L1 dual-view validator); the hot path is THREE
numpy tensor ops per ordinate — ``b``-build, scan, average — driven
entirely off the cache.

The structural claim
====================

The SN sweep is **forward substitution** on the block-triangular
(streaming + collision) operator under the DAG ordering induced by
the ordinate's direction.  Forward substitution is a fold over a
topologically-ordered work-unit stream.  Step 2.5c's cache
precomputes the per-ordinate transmission-emission pair :math:`(a, b)`
(Lewis-Miller §5.3) at solver construction; the sweep then evaluates
the Blelloch §1.5 closed form

.. math::

    \psi^{s}[n, i, g]
        \;=\; \mathrm{cumprod\_a}[n, i, g]\,
              \bigl(\psi^{s}_0[n, g] + \mathrm{cumsum}_{k \le i}
                    (b[n, k, g] / \mathrm{cumprod\_a}[n, k, g])\bigr)

via :func:`~orpheus.sn.spatial.scan.ordinate_scan`.

The 2-D wavefront sweep (:func:`_sweep_2d_wavefront`) retains its
own batched anti-diagonal scheduling — Step 2.6 Q2 explicitly
defers a unification of the work-unit shape; 2D's anti-diagonal
level scan and 1D's per-ordinate chain scan are genuinely different
folds.

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
* Blelloch, G. E. (1990).  *Prefix Sums and Their Applications*.
  CMU-CS-90-190 §1.5.  First-order linear recurrence closed form.

See also
========

* :class:`~orpheus.sn.spatial.cell_update.CellUpdate` — the per-cell
  strategy contract (Pattern 2 dual-view: ``update`` is the
  human-legible reference; the cache populator + scan is the fast
  path).
* :class:`~orpheus.sn.spatial.diamond.DiamondDifference` — the
  default cell-update strategy.
* :class:`~orpheus.sn.spatial.sweep_cache.GeometryCoefficients` /
  :class:`~orpheus.sn.spatial.sweep_cache.CollisionCache` — the
  two-stratum precomputed cache.
* :func:`~orpheus.sn.spatial.scan.ordinate_scan` — the Blelloch
  scan primitive.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from orpheus.geometry import CoordSystem

from .spatial.cell_update import UpstreamState
from .spatial.diamond import DiamondDifference
from .spatial.psi_half_angle_seed import CarlsonSweepContext
from .spatial.scan import ordinate_scan
from .spatial.sweep_cache import CollisionCache, GeometryCoefficients
from .sweep_graph import OctantLabel

if TYPE_CHECKING:
    from .angular_flux import AngularFlux
    from .boundary_flux import BoundaryFlux
    from .geometry import SNMesh
    from .sources import IsotropicSource, PerOrdinateSource


# ═══════════════════════════════════════════════════════════════════════
# Top-level dispatch
# ═══════════════════════════════════════════════════════════════════════


def transport_sweep(
    source: "PerOrdinateSource",
    sig_t: np.ndarray,
    sn_mesh: "SNMesh",
    boundary_flux: "BoundaryFlux",
    *,
    initial_guess: "AngularFlux | None" = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Perform one full transport sweep.

    Boundary conditions are read from ``sn_mesh`` (resolved at
    construction time from the geometry mesh's :class:`BC` declarations).
    The cell-update strategy is read from ``sn_mesh.cell_update``
    (defaults to :class:`DiamondDifference`).

    Single-source contract (R-1 Step 4 A1)
    --------------------------------------

    The sweep consumes ONE :class:`PerOrdinateSource` carrying the
    combined per-ordinate source magnitude
    :math:`Q_n(\\vec r, g)` — whatever combination of iso (P0 + n2n +
    fission) and aniso (P_ℓ ≥ 1) the caller produced.  The producer
    side has already applied the :math:`1/W` projection (Pattern 7
    per ``coding-elegance`` SKILL.md §"Convention crosswalk template",
    lesson L18); the sweep does NOT apply ``/W`` internally to ANY
    part of the source.

    External iso scalar sources :math:`Q(\\vec r, g)` (e.g. user-
    supplied fixed-source problems) project to per-ordinate at
    construction time via
    :meth:`~orpheus.sn.sources.PerOrdinateSource.from_isotropic`.
    Scattering-generated sources project at the producer boundary
    via the singledispatched
    :meth:`~orpheus.sn.scattering.ScatteringOperator.apply`.  Fission-
    generated sources project at the producer boundary via
    :meth:`~orpheus.sn.fission.FissionOperator.apply`.

    The legacy two-parameter convention (``iso_source: IsotropicSource``
    + ``aniso_source: PerOrdinateSource | None`` with sweep-internal
    ``/W``) is GONE.  See `#205
    <https://github.com/deOliveira-R/ORPHEUS/issues/205>`_ for the
    cross-method field architecture that will further refine the
    typed contract.

    History (chronological)
    -----------------------

    * Issue #196 PR-INDEX-5: PUBLIC contract is principled
      ``(ng, nx, ny)`` / ``(N, ng, nx, ny)``.
    * Issue #197 PR-TYPED-2: ``psi_bc: dict`` retired in favour of
      typed :class:`BoundaryFlux`.
    * Issue #197 PR-TYPED-3: typed :class:`IsotropicSource` /
      :class:`PerOrdinateSource` inputs.
    * Issue #197 PR-TYPED-4: bare-``np.ndarray`` overload retired.
    * R-1 Step 0 (2026-05-19): curvilinear Carlson seed derives from
      ``initial_guess`` (= previous iterate; ``None`` fallback uses
      the in-iteration source angular average).
    * R-1 Step 4 A1 (2026-05-21): ``iso_source`` parameter retired;
      sweep takes one :class:`PerOrdinateSource` carrying the combined
      per-ordinate source magnitude.  Producer-side ``/W`` projection
      everywhere.

    Parameters
    ----------
    source : PerOrdinateSource
        Per-ordinate volumetric source, shape ``(N, ng, nx, ny)``.
        Convention: **per-ordinate magnitude** (the producer has
        already applied any required ``/W`` projection).
    sig_t : np.ndarray
        Total cross-section, shape ``(ng, nx, ny)``.
    sn_mesh : SNMesh
        :class:`SNMesh` carrying geometry, BCs, quadrature, cell-update
        strategy.
    boundary_flux : BoundaryFlux
        Persistent :class:`BoundaryFlux` (mutated in place).  Build a
        zero-initialised instance via ``sn_mesh.zeros_boundary_flux()``.
    initial_guess : AngularFlux or None, optional
        Previous-iteration angular flux estimate, used for the
        curvilinear Carlson coupled-pole seed and the per-ordinate
        spatial-upstream inflow at the pole cell.  Ignored on slab.
        ``None`` (default) selects the in-iteration fallback seed.

    Returns
    -------
    angular_flux
        Shape ``(N, ng, nx, ny)``.
    scalar_flux
        Shape ``(ng, nx, ny)`` — :math:`\\sum_n w_n \\psi_n`.

    Dispatch:

    * 1-D meshes (``sn_mesh.reduced is not None``) →
      :func:`_sweep_1d_unified` (slab, sphere, cylinder — one body via
      the two-stratum cache).
    * 2-D Cartesian → :func:`_sweep_2d_wavefront` (anti-diagonal
      scheduling; Step 2.6 Q2 deferred).
    """
    Q = _unwrap_source(source)
    reduced = sn_mesh.reduced
    if reduced is not None:
        return _sweep_1d_unified(
            Q, sig_t, sn_mesh, boundary_flux,
            initial_guess=initial_guess,
        )
    return _sweep_2d_wavefront(
        Q, sig_t, sn_mesh, boundary_flux,
    )


def _unwrap_source(source: "PerOrdinateSource") -> np.ndarray:
    """Unwrap typed :class:`PerOrdinateSource` to bare ndarray.

    Issue #197 PR-TYPED-4 — strict typed input.  R-1 Step 4 A1 collapsed
    the iso / aniso parameter pair into a single per-ordinate source.
    The internal hot path consumes bare ndarray; this helper performs
    the unwrap once at the public boundary.
    """
    from .sources import PerOrdinateSource
    if not isinstance(source, PerOrdinateSource):
        raise TypeError(
            f"transport_sweep: source must be "
            f"PerOrdinateSource (R-1 Step 4 A1); got "
            f"{type(source).__name__}"
        )
    return source.values


# ═══════════════════════════════════════════════════════════════════════
# 1-D unified sweep — slab + sphere + cylinder via the cache
# ═══════════════════════════════════════════════════════════════════════


def _sweep_1d_unified(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: "SNMesh",
    boundary_flux: "BoundaryFlux",
    *,
    initial_guess: "AngularFlux | None" = None,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Geometry-blind 1-D SN sweep — three numpy tensor ops per ordinate.

    Replaces ``_sweep_1d_cartesian`` and ``_sweep_1d_curvilinear`` with
    one body driven by the two-stratum precomputed cache.  Slab,
    sphere, and cylinder share THE SAME scan expression; the
    cache abstracts the geometry difference (slab carries neutral
    curvature values; the M-M angular thread and Carlson seed run only
    when ``geom.level_ordinates is not None``).

    Per-ordinate hot path
    ---------------------

    1. ``b = 2 · (QV_chain + ang_contrib) · coll.inverse_denom[n]``
       — per-cell (in chain order) affine additive coefficient.
    2. ``psi_face = ordinate_scan(coll.a_attenuation[n], b, psi_in)``
       — the Blelloch closed form, three numpy ops internally.
    3. ``psi_avg = 0.5 · (psi_in_chain + psi_face)``
       — DD spatial closure.

    For the rare degenerate cylindrical pure-azimuthal ordinate
    (``geom.is_degenerate[n] == True``, ``|η| < 10^{-15}``), the scan
    is meaningless and the slow per-cell ``cell_update.update`` path
    runs instead.

    Cache provenance
    ----------------

    The cache is stashed on ``sn_mesh`` by :class:`SNSolver.__init__`.
    If the sweep is invoked outside the solver (e.g. ad-hoc tests),
    the cache is built lazily on first call and held on the mesh.

    Bit-identity contract
    ---------------------

    The cache-driven path produces algebraically the SAME values as the
    per-cell ``cell_update.update`` reference iteration (the Pattern 2
    dual-view contract).  The cache's ``a_attenuation`` field IS the
    per-ordinate sequence of transmission coefficients that
    Step 2.5b's ``affine_coefficients`` builder produced — but
    precomputed once at solver construction rather than rebuilt every
    sweep.
    The Pattern 2 dual-view test
    (``tests/sn/spatial/test_sweep_cache.py``) pins this at
    ``rtol=1e-13`` across the parametrised geometry × ng × source
    grid.  Slab regression snapshots stay bit-identical at
    ``rtol=1e-12``.
    """
    geom = _ensure_geom_cache(sn_mesh)
    coll = _ensure_coll_cache(sn_mesh, sig_t, geom)
    return _run_1d_sweep(
        Q, sig_t, sn_mesh, boundary_flux, geom, coll,
        initial_guess=initial_guess,
    )


def _ensure_geom_cache(sn_mesh: "SNMesh") -> GeometryCoefficients:
    """Return the geometry cache, building it on first use if absent."""
    cache = getattr(sn_mesh, "_geom_cache", None)
    if cache is None:
        cache = GeometryCoefficients.from_mesh_and_quad(sn_mesh)
        sn_mesh._geom_cache = cache  # type: ignore[attr-defined]
    return cache


def _ensure_coll_cache(
    sn_mesh: "SNMesh",
    sig_t: np.ndarray,
    geom: GeometryCoefficients,
) -> CollisionCache:
    """Return the collision cache, building it on first use if absent.

    The expected invariant (per cache-invariance test #4) is that the
    cache is constructed by :class:`SNSolver.__init__` and consumed by
    every sweep without rebuild.  Ad-hoc test callers may bypass the
    solver — in that case the cache is built lazily here.

    No bridge needed under PR-INDEX-3: ``sig_t`` arrives as principled
    ``(ng, nx, ny=1)`` and the cache consumes ``(ng, nx)`` — a single
    slice on the degenerate ``ny`` axis suffices.
    """
    cache = getattr(sn_mesh, "_coll_cache", None)
    if cache is None:
        # 1-D meshes: sig_t shape is (ng, nx, 1).  Slice the trailing
        # degenerate axis — the principled (ng, nx) layout the cache
        # expects natively.
        sig_t_1d = sig_t[:, :, 0]  # (ng, nx)
        cache = CollisionCache.from_geometry(geom, sig_t_1d)
        sn_mesh._coll_cache = cache  # type: ignore[attr-defined]
    return cache


def _run_1d_sweep(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: "SNMesh",
    boundary_flux: "BoundaryFlux",
    geom: GeometryCoefficients,
    coll: CollisionCache,
    *,
    initial_guess: "AngularFlux | None" = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Inner body of the unified 1-D sweep.

    Issue #196 PR-INDEX-1 through PR-INDEX-5: internal arrays AND the
    public ``transport_sweep`` signature both carry the principled
    ``(N, ng, nx, ny=1)`` layout (energy ``g`` is the *second* axis,
    NOT trailing; see :ref:`theory-sn-index-convention`).  No
    entry/exit transposes are required at the public boundary —
    caller-side principled-layout inputs flow directly through the
    sweep body.

    Issue #196 PR-INDEX-2: :class:`CollisionCache` fields carry the
    principled ``(N, ng, nx)`` layout natively; the bridge transposes at
    the cache-access sites are gone.  :class:`GeometryCoefficients` stays
    on ``(N, nx)`` / ``(N,)`` shapes — no group axis, no flip needed.

    Splits cleanly into setup (BC inflow, source pre-scale, Carlson
    seed when curvilinear) and a per-direction or per-ordinate scan:

    * **SLAB** (joint-batch): ordinates within a chain direction are
      independent (no M-M angular thread), so one
      :func:`ordinate_scan` call per chain handles the entire chain's
      ordinates at once with shape ``(nx, K, ng)`` where ``K`` is the
      number of ordinates in that direction (``N/2`` for symmetric GL).
      Exactly 2 scan calls per sweep regardless of ``N`` or ``ng``.

    * **CURVILINEAR** (sphere/cylinder, per-ordinate): the M-M angular
      thread couples ordinates sequentially within a μ-level (the
      Hébert §3.9.4 Eqs. 3.437/3.439 recurrence reads
      ``psi_angle[chain]`` updated by the *previous* ordinate in the
      level).  One ``ordinate_scan`` per ordinate per level — unchanged
      from PR-INDEX-1's pre-state.  A future parallel-prefix
      reformulation of the M-M recurrence could unlock joint-batch for
      curvilinear too (research-level; deferred per plan §7).
    """
    quad = sn_mesh.quad
    N = quad.N
    nx = sn_mesh.nx
    ng = Q.shape[1]                                          # (N, ng, nx, ny=1)
    weights = quad.weights
    mu = quad.mu_x

    # ── Entry layout — PR-INDEX-5: public contract is principled ──────
    # Drop the degenerate trailing ``ny=1`` axis to obtain the
    # (N, ng, nx) working layout.  No transpose required.
    Q_per_ord = Q[:, :, :, 0]                                # (N, ng, nx)
    sig_t_p = sig_t[:, :, 0]                                 # (ng, nx)
    V = sn_mesh.volumes[:, 0]                                # (nx,) — no group axis
    cell_update = sn_mesh.cell_update

    coord = sn_mesh.reduced.coord  # type: ignore[union-attr]
    is_slab = coord is CoordSystem.CARTESIAN
    is_sphere = coord is CoordSystem.SPHERICAL

    # ── Common pre-scale ──────────────────────────────────────────────
    # R-1 Step 4 A1 — single per-ordinate source.  The producer applied
    # ``1/W`` already; the sweep multiplies by cell volume V only.
    # No iso/aniso distinction internally — every WDD recurrence
    # consumes the same ``QV_per_ord``.
    QV_per_ord = Q_per_ord * V[None, None, :]                # (N, ng, nx)

    # Internal principled layout — angular flux (N, ng, nx, 1),
    # scalar flux (ng, nx) working buffer (ny added at return).
    angular_flux = np.zeros((N, ng, nx, 1))
    scalar_flux = np.zeros((ng, nx))

    # ── BC inflow + per-level Carlson seed (curvilinear only) ─────────
    if is_slab:
        bc_left_obj = sn_mesh.bc_left
        bc_right_obj = sn_mesh.bc_right
        # Issue #197 PR-TYPED-2: BoundaryFlux carries named face buffers.
        # ``zeros_boundary_flux`` pre-allocates them; the
        # ``if ... is None`` guard handles the rare case where a caller
        # builds a bare BoundaryFlux without zeros (e.g. a synthetic test).
        if boundary_flux.xmin_face is None:
            boundary_flux.xmin_face = np.zeros((N, ng))
        if boundary_flux.xmax_face is None:
            boundary_flux.xmax_face = np.zeros((N, ng))
        bc_left_face = boundary_flux.xmin_face
        bc_right_face = boundary_flux.xmax_face
        inflow_left = bc_left_obj.apply(bc_left_face)    # (N, ng)
        inflow_right = bc_right_obj.apply(bc_right_face)  # (N, ng)
        levels = [None]
        level_ordinates_list = [list(range(N))]
        bc_outer = None
        sigma_t_gx = None
        dr = None
        inflow_full = None
    else:
        bc_outer_obj = sn_mesh.bc_right
        # Issue #197 PR-TYPED-2: curvilinear outer radial face routes
        # through :attr:`BoundaryFlux.xmax_face` (the unified naming for
        # 1-D outer face, retiring the ``"bc_sph"`` / ``"bc_cyl"`` dict
        # keys that previously bifurcated spherical vs cylindrical).
        if boundary_flux.xmax_face is None:
            boundary_flux.xmax_face = np.zeros((N, ng))
        bc_outer = boundary_flux.xmax_face
        inflow_full = bc_outer_obj.apply(bc_outer)  # (N, ng)

        # Per-level Carlson coupled-pole seed delegates to the M-M
        # closure's ``psi_half_seed`` strategy — the SAME strategy the
        # matvec uses (Pattern 2 single source of truth).  Pre-amble
        # only stashes the (σ_t, Δr) bundle the per-level context
        # needs; the level loop builds the per-level CarlsonSweepContext
        # and calls ``closure.psi_half_seed(psi_level, context)``.
        # See ``MorelMontryAngularSweep.precompute_psi_state`` (the
        # matvec entry point) for the symmetric routing.
        sigma_t_gx = sig_t_p                                  # (ng, nx)
        dr = sn_mesh.dx

        if is_sphere:
            levels = [None]
            level_ordinates_list = [list(range(N))]
        else:
            level_indices = quad.level_indices  # type: ignore[attr-defined]
            levels = list(range(len(level_indices)))
            level_ordinates_list = [list(level_indices[p]) for p in levels]

        inflow_left = inflow_right = None
        bc_left_face = bc_right_face = None

    # ── SLAB joint-batch fast path ────────────────────────────────────
    #
    # Slab has no M-M angular thread, no degenerate ordinates, and one
    # chain per direction.  Group ordinates by chain direction and run
    # ONE ordinate_scan per chain.  Exactly 2 scan calls per sweep.
    if is_slab:
        # Partition ordinates by direction sign (μ ≥ 0 → forward chain).
        forward_mask = mu >= 0
        forward_ords = np.where(forward_mask)[0]
        backward_ords = np.where(~forward_mask)[0]

        for direction_sign, ords in ((+1, forward_ords), (-1, backward_ords)):
            if ords.size == 0:
                continue
            K = ords.size

            # Chain order is identical across ordinates in one direction
            # for slab — pick from the first ordinate.
            chain = geom.chain_idx[ords[0]]                   # (nx,)
            inv = geom.chain_idx_inv[ords[0]]                 # (nx,)

            # Per-ordinate inflow (cells degenerate, group axis full).
            psi_in_chain = (
                inflow_left[ords] if direction_sign == +1
                else inflow_right[ords]
            )                                                  # (K, ng)

            # Per-ordinate source in chain order — R-1 Step 4 A1's
            # single-source convention: ``QV_per_ord`` already encodes
            # per-ordinate magnitude × cell volume.  Slice the K
            # ordinates and reorder along the chain axis.
            QV_full_chain = QV_per_ord[ords][:, :, chain]      # (K, ng, nx)

            # Cache fields are (N, ng, nx) natively under PR-INDEX-2.
            # Indexed slice [ords] yields (K, ng, nx) — no transpose.
            inv_denom_chain = coll.inverse_denom[ords]         # (K, ng, nx)
            a_atten_chain = coll.a_attenuation[ords]           # (K, ng, nx)

            # b shape needed for ordinate_scan: (nx, K, ng) with cells
            # on axis 0 (scan axis).  Build (K, ng, nx) first, then
            # transpose.
            b_chain = 2.0 * QV_full_chain * inv_denom_chain   # (K, ng, nx)
            # Scan-input layout: (nx, K, ng).
            a_scan = np.transpose(a_atten_chain, (2, 0, 1))   # (nx, K, ng)
            b_scan = np.transpose(b_chain, (2, 0, 1))         # (nx, K, ng)

            # ONE scan call per chain — joint-batched over (K, ng).
            psi_face_chain_scan = ordinate_scan(
                a_scan, b_scan, psi_in_chain,
            )                                                  # (nx, K, ng)

            # DD spatial closure — face-in shifts upstream by 1.
            psi_face_in_chain = np.empty_like(psi_face_chain_scan)
            psi_face_in_chain[0] = psi_in_chain
            psi_face_in_chain[1:] = psi_face_chain_scan[:-1]
            psi_avg_scan = 0.5 * (psi_face_in_chain + psi_face_chain_scan)
            # (nx, K, ng) → per-ordinate (ng, nx) via reorder.
            psi_avg_per_ord = np.transpose(psi_avg_scan, (1, 2, 0))  # (K, ng, nx)

            # Scatter back to cell-index order + write angular_flux,
            # accumulate scalar_flux.
            psi_avg_cell_order = psi_avg_per_ord[:, :, inv]   # (K, ng, nx)
            angular_flux[ords, :, :, 0] = psi_avg_cell_order

            # scalar_flux += Σ_n w_n · ψ_n  (broadcast over K).
            w_ords = weights[ords]                            # (K,)
            scalar_flux += np.einsum(
                "k,kgx->gx", w_ords, psi_avg_cell_order,
            )

            # Persist outflow at the appropriate boundary face — the
            # last chain output is the outgoing-face flux on that side.
            if direction_sign == +1:
                bc_right_face[ords] = psi_face_chain_scan[-1]  # (K, ng)
            else:
                bc_left_face[ords] = psi_face_chain_scan[-1]   # (K, ng)

    # ── CURVILINEAR per-ordinate path ─────────────────────────────────
    #
    # M-M angular thread couples ordinates sequentially within a level
    # (psi_angle[chain] is updated by ordinate m → consumed by m+1).
    # Joint-batch over ordinates is blocked; loop stays per-ordinate.
    else:
        # M-M closure owns the per-level Carlson coupled-pole seed
        # (Pattern 2 single source of truth — the matvec consumes the
        # SAME ``psi_half_seed`` strategy via
        # :meth:`MorelMontryAngularSweep.precompute_psi_state`).  The
        # strategy reads ψ at the level's ordinates (the previous
        # iterate when ``initial_guess`` is supplied; zeros on cold
        # start) and emits the cell-centred half-angle face flux
        # ``φ̄_{1/2,i,g}`` per Hébert §3.9.4 Eqs. (3.432)-(3.435).
        closure = sn_mesh.pole_angular_closure
        if initial_guess is not None:
            # (N, ng, nx, 1) → (ng, N, nx)
            psi_g_first = initial_guess.values.transpose(1, 0, 2, 3)[..., 0]
        else:
            psi_g_first = None

        for p_idx, level in enumerate(levels):
            ordinates_in_level = level_ordinates_list[p_idx]
            ords_arr = np.asarray(ordinates_in_level)
            mu_in_level = mu[ords_arr]
            most_inward_global = int(ords_arr[int(np.argmin(mu_in_level))])
            bc_outer_value = inflow_full[most_inward_global, :]  # type: ignore[index]
            level_weights = weights[ords_arr]
            if psi_g_first is not None:
                psi_level = psi_g_first[:, ords_arr, :]          # (ng, M_p, nx)
            else:
                psi_level = np.zeros((ng, ords_arr.size, nx))
            carlson_ctx = CarlsonSweepContext(
                sigma_t=sigma_t_gx,
                dr=dr,
                mu_quad=mu_in_level.copy(),
                weights=level_weights.copy(),
                bc_outer_value=bc_outer_value,
            )
            phi_aux = closure.psi_half_seed(psi_level, carlson_ctx)  # (ng, nx)
            psi_angle = phi_aux.copy()                            # (ng, nx) — principled

            for m_local, global_n in enumerate(ordinates_in_level):
                mu_n = mu[global_n]
                w_n = weights[global_n]
                chain = geom.chain_idx[global_n]

                # Per-ordinate source assembly (R-1 Step 4 A1):
                # ``QV_per_ord[global_n]`` is the per-ordinate source ×
                # cell volume for ordinate ``global_n``, shape (ng, nx).
                QV_full = QV_per_ord[global_n]                  # (ng, nx)
                QV_chain = QV_full[:, chain]                    # (ng, nx)

                # Per-ordinate spatial-upstream inflow (ng,).
                if mu_n < 0:
                    psi_in = inflow_full[global_n]               # type: ignore[index]
                else:
                    # R-1 Step 0: pole spatial-upstream derived from the
                    # caller's ``initial_guess`` — same trace-space logic
                    # as the matvec (read the input flux's pole cell).
                    if initial_guess is not None:
                        psi_in = initial_guess.values[global_n, :, 0, 0]
                    else:
                        psi_in = np.zeros(ng)

                # Degenerate cyl-axis ordinate: slow per-cell path.
                if geom.is_degenerate[global_n]:
                    ordinate_idx = global_n if is_sphere else m_local
                    visits = list(sn_mesh.dag_walk(
                        ordinate_idx=ordinate_idx,
                        mu_level_idx=level,
                    ))
                    for visit in visits:
                        i = visit.cell_idx
                        upstream = UpstreamState(
                            spatial_upstream=psi_in,
                            angular_upstream=psi_angle[:, i],
                        )
                        # cell_update.update expects per-cell (ng,)
                        # arrays — sig_t / source slice on the cell axis.
                        result = cell_update.update(
                            visit=visit,
                            total_xs=sig_t_p[:, i],
                            source=QV_full[:, i],
                            upstream_state=upstream,
                        )
                        psi = result.cell_average_flux           # (ng,)
                        psi_angle[:, i] = result.outgoing_angular_state
                        angular_flux[global_n, :, i, 0] = psi
                        scalar_flux[:, i] += w_n * psi
                    continue

                # Non-degenerate fast path: per-ordinate scan (ng, nx).
                # psi_angle on (ng, nx); chain reorders the nx axis.
                psi_a_in_chain = psi_angle[:, chain].copy()      # (ng, nx)
                ang_contrib = (
                    geom.dA_w[global_n] * geom.c_in[global_n]
                )[None, :] * psi_a_in_chain                       # (ng, nx)

                # Cache fields are (N, ng, nx) natively under PR-INDEX-2.
                # Indexed slice [global_n] yields (ng, nx) — no transpose.
                inv_denom_p = coll.inverse_denom[global_n]       # (ng, nx)
                a_atten_p = coll.a_attenuation[global_n]         # (ng, nx)
                b = 2.0 * (QV_chain + ang_contrib) * inv_denom_p  # (ng, nx)

                # ordinate_scan: leading axis is the scan/cell axis.
                # Pass (nx, ng) — transpose from (ng, nx).
                psi_face_chain = ordinate_scan(
                    a_atten_p.T, b.T, psi_in,
                )                                                 # (nx, ng)

                # DD spatial closure — vectorised cell-average.
                psi_face_in_chain = np.empty_like(psi_face_chain)
                psi_face_in_chain[0] = psi_in
                psi_face_in_chain[1:] = psi_face_chain[:-1]
                psi_avg_chain = 0.5 * (psi_face_in_chain + psi_face_chain)
                # Principled view: (ng, nx).
                psi_avg_chain_p = psi_avg_chain.T                # (ng, nx)

                # M-M angular thread (curvilinear-only).
                psi_angle_out_chain_p = (
                    geom.tau_inv[global_n] * psi_avg_chain_p
                    - geom.mm_a_in_coeff[global_n] * psi_a_in_chain
                )                                                 # (ng, nx)
                psi_angle[:, chain] = psi_angle_out_chain_p

                # Scatter back to cell-index order + writes.
                inv = geom.chain_idx_inv[global_n]
                psi_avg_p = psi_avg_chain_p[:, inv]              # (ng, nx)
                angular_flux[global_n, :, :, 0] = psi_avg_p
                scalar_flux += w_n * psi_avg_p

                # Persist outflow at the outer face for outward ordinates.
                if mu_n >= 0 and abs(mu_n) >= sn_mesh._DEGENERATE_ABS_ETA_THRESHOLD:
                    bc_outer[global_n] = psi_face_chain[-1]      # (ng,)

    # ── Exit — PR-INDEX-5: caller consumes principled layout ──────────
    # R-1 Step 0: NO iteration-cache write-back.  The caller threads the
    # returned ``angular_flux`` as ``initial_guess`` to the NEXT sweep —
    # that's all the "history" needed.  The matvec already operates this
    # way; the sweep now mirrors it.
    # angular_flux is already (N, ng, nx, 1); scalar_flux needs the ny=1
    # axis added back to satisfy the (ng, nx, ny) public contract.
    return angular_flux, scalar_flux[:, :, None]


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
    sn_mesh: "SNMesh",
    boundary_flux: "BoundaryFlux",
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

    R-1 Step 4 A1: single per-ordinate source ``Q`` shape
    ``(N, ng, nx, ny)`` carries the producer-side-projected magnitude.
    Sweep does NOT apply ``/W`` internally.  The legacy iso/aniso
    parameter pair is GONE.

    Issue #196 PR-INDEX-5: fully principled.  ``Q`` is consumed in
    principled ``(N, ng, nx, ny)``; ``sig_t`` is principled
    ``(ng, nx, ny)``.  The persistent ``psi_x`` / ``psi_y`` BC
    buffers are principled ``(N, ng, nx+1, ny)`` / ``(N, ng, nx,
    ny+1)``; ``angular_flux`` is returned principled
    ``(N, ng, nx, ny)``; ``scalar_flux`` is principled
    ``(ng, nx, ny)``.  The PR-INDEX-4 ``BRIDGE_pure_z_to_legacy``
    transpose is GONE — the pure-z degenerate branch writes directly
    into the principled angular_flux buffer.

    Bit-identity to legacy
    ----------------------

    The PR-INDEX-1..4 sequence preserved bit-identity at every step;
    PR-INDEX-5 flips the persistent BC buffers (so the BC apply face
    slices reorder), which IS principled-equivalent to the legacy
    operation per ``vv-principles`` § "Bit-identity vs principled-
    equivalence" — the values written / read at each face slot are
    the same, only the memory layout reorders.  Snapshots regenerate
    under the principled layout (the snapshot generator stores the
    final values in principled order).
    """
    ng, nx, ny = sig_t.shape
    quad = sn_mesh.quad
    N = quad.N
    weights = quad.weights

    angular_flux = np.zeros((N, ng, nx, ny))
    scalar_flux = np.zeros((ng, nx, ny))

    # Persistent boundary-flux buffers for reflective BCs (mutated in
    # place across sweep calls — partner reads need the previous
    # iteration's outgoing-face writes).  PR-INDEX-5: principled
    # ``(N, ng, nx+1, ny)`` / ``(N, ng, nx, ny+1)``.
    # Issue #197 PR-TYPED-2: the two persistent buffers live on the
    # typed :class:`BoundaryFlux` as named attributes; per-face slices
    # (``boundary_flux.xmin`` / ``xmax`` / ``ymin`` / ``ymax``) are
    # auto-derived views.
    if boundary_flux.xmin_xmax_buf is None:
        boundary_flux.xmin_xmax_buf = np.zeros((N, ng, nx + 1, ny))
    if boundary_flux.ymin_ymax_buf is None:
        boundary_flux.ymin_ymax_buf = np.zeros((N, ng, nx, ny + 1))
    psi_x = boundary_flux.xmin_xmax_buf  # (N, ng, nx+1, ny) — principled
    psi_y = boundary_flux.ymin_ymax_buf  # (N, ng, nx, ny+1) — principled

    # R-1 Step 4 A1: ``Q`` is per-ordinate magnitude (N, ng, nx, ny).
    # No ``/W`` applied here — the producer normalised at the apply
    # boundary (Pattern 7).
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
            Q_pure_z = Q[oct_idx]                         # (N_oct, ng, nx, ny)
            # sig_t (ng, nx, ny) broadcasts against (N_oct, ng, nx, ny).
            psi_avg_pure_z = Q_pure_z / sig_t              # (N_oct, ng, nx, ny)
            angular_flux[oct_idx] = psi_avg_pure_z
            scalar_flux += np.einsum(
                "ngij,n->gij",
                psi_avg_pure_z, weights[oct_idx],
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
        # Under principled layout the BC operator consumes a
        # ``(N, ng_axis_2)`` face block.  Slice the face slot, then
        # transpose to ``(N, ng_axis_2)`` for the BC's ordinate-leading
        # contract.  ``bc.apply`` returns the same shape; transpose back.
        if sx_eff >= 0:
            full_face_x = sn_mesh.bc_xmin.apply(psi_x[:, :, 0, :])
            psi_x[oct_idx, :, 0, :] = full_face_x[oct_idx]
        else:
            full_face_x = sn_mesh.bc_xmax.apply(psi_x[:, :, nx, :])
            psi_x[oct_idx, :, nx, :] = full_face_x[oct_idx]

        if sy_eff >= 0:
            full_face_y = sn_mesh.bc_ymin.apply(psi_y[:, :, :, 0])
            psi_y[oct_idx, :, :, 0] = full_face_y[oct_idx]
        else:
            full_face_y = sn_mesh.bc_ymax.apply(psi_y[:, :, :, ny])
            psi_y[oct_idx, :, :, ny] = full_face_y[oct_idx]

        # ── Per-octant buffers for the graph apply ────────────────
        psi_x_oct = psi_x[oct_idx].copy()    # (N_oct, ng, nx+1, ny)
        psi_y_oct = psi_y[oct_idx].copy()    # (N_oct, ng, nx, ny+1)
        # R-1 Step 4 A1 — slice the per-ordinate source for this octant.
        Q_octant = Q[oct_idx]                # (N_oct, ng, nx, ny)
        angular_flux_oct = np.zeros((N_oct, ng, nx, ny))

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


__all__ = [
    "transport_sweep",
]
