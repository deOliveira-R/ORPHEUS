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
from .spatial.psi_half_angle_seed import carlson_inward_sweep_from_source
from .spatial.scan import ordinate_scan
from .spatial.sweep_cache import CollisionCache, GeometryCoefficients
from .sweep_graph import OctantLabel

if TYPE_CHECKING:
    from .geometry import SNMesh


# ═══════════════════════════════════════════════════════════════════════
# Top-level dispatch
# ═══════════════════════════════════════════════════════════════════════


def transport_sweep(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: "SNMesh",
    psi_bc: dict,
    Q_aniso: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Perform one full transport sweep.

    Boundary conditions are read from ``sn_mesh`` (resolved at
    construction time from the geometry mesh's :class:`BC` declarations).
    The cell-update strategy is read from ``sn_mesh.cell_update``
    (defaults to :class:`DiamondDifference`).

    Dispatch:

    * 1-D meshes (``sn_mesh.reduced is not None``) →
      :func:`_sweep_1d_unified` (slab, sphere, cylinder — one body via
      the two-stratum cache).
    * 2-D Cartesian → :func:`_sweep_2d_wavefront` (anti-diagonal
      scheduling; Step 2.6 Q2 deferred).
    """
    reduced = sn_mesh.reduced
    if reduced is not None:
        return _sweep_1d_unified(Q, sig_t, sn_mesh, psi_bc, Q_aniso)
    return _sweep_2d_wavefront(Q, sig_t, sn_mesh, psi_bc, Q_aniso)


# ═══════════════════════════════════════════════════════════════════════
# 1-D unified sweep — slab + sphere + cylinder via the cache
# ═══════════════════════════════════════════════════════════════════════


def apply_sweep_1d(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: "SNMesh",
    psi_bc: dict,
    Q_aniso: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Public alias for :func:`_sweep_1d_unified`.

    Provided so the JAX ``lax.scan(f, init, xs)`` vocabulary
    (cross-domain-attacker Frame 3) is reachable by name — the cache
    plays ``xs``, the BC inflow plays ``init``, the returned ψ is the
    scan trajectory.  No code change; vocabulary only.
    """
    return _sweep_1d_unified(Q, sig_t, sn_mesh, psi_bc, Q_aniso)


def _sweep_1d_unified(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: "SNMesh",
    psi_bc: dict,
    Q_aniso: np.ndarray | None = None,
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
    return _run_1d_sweep(Q, sig_t, sn_mesh, psi_bc, Q_aniso, geom, coll)


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
    """
    cache = getattr(sn_mesh, "_coll_cache", None)
    if cache is None:
        # 1-D meshes: sig_t shape is (nx, 1, ng); reduce to (nx, ng).
        sig_t_1d = sig_t[:, 0, :]
        cache = CollisionCache.from_geometry(geom, sig_t_1d)
        sn_mesh._coll_cache = cache  # type: ignore[attr-defined]
    return cache


def _run_1d_sweep(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: "SNMesh",
    psi_bc: dict,
    Q_aniso: np.ndarray | None,
    geom: GeometryCoefficients,
    coll: CollisionCache,
) -> tuple[np.ndarray, np.ndarray]:
    """Inner body of the unified 1-D sweep.

    Splits cleanly into setup (BC inflow, source pre-scale, Carlson
    seed when curvilinear) and a per-ordinate loop whose body is three
    tensor ops (slab) or four (curvilinear, with the M-M angular
    update + Carlson-seeded angular thread).
    """
    quad = sn_mesh.quad
    N = quad.N
    nx = sn_mesh.nx
    ng = Q.shape[2]
    weights = quad.weights
    mu = quad.mu_x

    Q_1d = Q[:, 0, :]           # (nx, ng)
    sig_t_1d = sig_t[:, 0, :]   # (nx, ng)
    V = sn_mesh.volumes[:, 0]   # (nx,)
    cell_update = sn_mesh.cell_update

    coord = sn_mesh.reduced.coord  # type: ignore[union-attr]
    is_slab = coord is CoordSystem.CARTESIAN
    is_sphere = coord is CoordSystem.SPHERICAL

    # ── Common pre-scale ──────────────────────────────────────────────
    weight_norm = 1.0 / weights.sum()
    QV_iso = Q_1d * V[:, None] * weight_norm   # (nx, ng)
    has_aniso = Q_aniso is not None
    if has_aniso:
        Q_aniso_1d = Q_aniso[:, :, 0, :] * weight_norm  # (N, nx, ng)
    else:
        Q_aniso_1d = None

    angular_flux = np.zeros((N, nx, 1, ng))
    scalar_flux = np.zeros((nx, ng))

    # ── BC inflow + per-level Carlson seed (curvilinear only) ─────────
    if is_slab:
        bc_left_obj = sn_mesh.bc_left
        bc_right_obj = sn_mesh.bc_right
        if "bc_1d_left_face" not in psi_bc:
            psi_bc["bc_1d_left_face"] = np.zeros((N, ng))
            psi_bc["bc_1d_right_face"] = np.zeros((N, ng))
        bc_left_face = psi_bc["bc_1d_left_face"]
        bc_right_face = psi_bc["bc_1d_right_face"]
        inflow_left = bc_left_obj.apply(bc_left_face)    # (N, ng)
        inflow_right = bc_right_obj.apply(bc_right_face)  # (N, ng)
        levels = [None]
        level_ordinates_list = [list(range(N))]
        bc_outer = None
        pole_key = None
        phi_prev_key = None
        Q_bar_iso = None
        sigma_t_gx = None
        dr = None
        inflow_full = None
    else:
        bc_key = "bc_sph" if is_sphere else "bc_cyl"
        pole_key = "psi_pole" if is_sphere else "psi_pole_cyl"
        phi_prev_key = "phi_0_prev" if is_sphere else "phi_0_prev_cyl"

        bc_outer_obj = sn_mesh.bc_right
        if bc_key not in psi_bc:
            psi_bc[bc_key] = np.zeros((N, ng))
        bc_outer = psi_bc[bc_key]
        inflow_full = bc_outer_obj.apply(bc_outer)  # (N, ng)

        # Carlson seed source (Phase G Step 2 Path C).
        sigma_t_gx = sig_t_1d.T  # (ng, nx)
        dr = sn_mesh.dx
        if phi_prev_key in psi_bc:
            phi_0_prev = psi_bc[phi_prev_key]
            Q_bar_iso = sigma_t_gx * phi_0_prev.T / weights.sum()
        else:
            Q_bar_iso = Q_1d.T / weights.sum()

        if is_sphere:
            levels = [None]
            level_ordinates_list = [list(range(N))]
        else:
            level_indices = quad.level_indices  # type: ignore[attr-defined]
            levels = list(range(len(level_indices)))
            level_ordinates_list = [list(level_indices[p]) for p in levels]

        inflow_left = inflow_right = None
        bc_left_face = bc_right_face = None

    # ── Per-level / per-ordinate loop ─────────────────────────────────
    for p_idx, level in enumerate(levels):
        ordinates_in_level = level_ordinates_list[p_idx]

        # Carlson coupled-pole seed for this level's angular DAG.
        if is_slab:
            psi_angle = None
        else:
            ords_arr = np.asarray(ordinates_in_level)
            mu_in_level = mu[ords_arr]
            most_inward_global = int(ords_arr[int(np.argmin(mu_in_level))])
            bc_outer_value = inflow_full[most_inward_global, :]  # type: ignore[index]
            phi_aux = carlson_inward_sweep_from_source(
                Q_bar=Q_bar_iso,
                sigma_t=sigma_t_gx,
                dr=dr,
                bc_outer_value=bc_outer_value,
            )                                                    # (ng, nx)
            psi_angle = phi_aux.T.copy()                          # (nx, ng)

        for m_local, global_n in enumerate(ordinates_in_level):
            mu_n = mu[global_n]
            w_n = weights[global_n]
            chain = geom.chain_idx[global_n]

            # ── Per-ordinate source assembly ─────────────────────────
            if has_aniso:
                QV_full = QV_iso + Q_aniso_1d[global_n] * V[:, None]  # type: ignore[index]
            else:
                QV_full = QV_iso
            QV_chain = QV_full[chain]                                # (nx, ng)

            # ── Per-ordinate spatial-upstream inflow ─────────────────
            if is_slab:
                direction_sign = +1 if mu_n >= 0 else -1
                psi_in = (
                    inflow_left[global_n] if direction_sign == +1     # type: ignore[index]
                    else inflow_right[global_n]                       # type: ignore[index]
                )
            else:
                if mu_n < 0:
                    psi_in = inflow_full[global_n]                    # type: ignore[index]
                else:
                    if pole_key in psi_bc:
                        psi_in = psi_bc[pole_key][global_n]
                    else:
                        psi_in = np.zeros(ng)

            # ── Degenerate cyl-axis ordinate: slow per-cell path ─────
            if geom.is_degenerate[global_n]:
                # Per-cell update — no spatial chain.
                ordinate_idx = global_n if is_sphere else m_local
                visits = list(sn_mesh.dag_walk(
                    ordinate_idx=ordinate_idx,
                    mu_level_idx=level,
                ))
                for visit in visits:
                    i = visit.cell_idx
                    upstream = UpstreamState(
                        spatial_upstream=psi_in,
                        angular_upstream=psi_angle[i] if psi_angle is not None else None,
                    )
                    result = cell_update.update(
                        visit=visit,
                        total_xs=sig_t_1d[i],
                        source=QV_full[i],
                        upstream_state=upstream,
                    )
                    psi = result.cell_average_flux
                    if psi_angle is not None:
                        psi_angle[i] = result.outgoing_angular_state
                    angular_flux[global_n, i, 0, :] = psi
                    scalar_flux[i] += w_n * psi
                continue

            # ── Non-degenerate fast path: three tensor ops ───────────
            if psi_angle is not None:
                psi_a_in_chain = psi_angle[chain].copy()
                ang_contrib = (
                    geom.dA_w[global_n] * geom.c_in[global_n]
                )[:, None] * psi_a_in_chain
                b = 2.0 * (QV_chain + ang_contrib) * coll.inverse_denom[global_n]
            else:
                psi_a_in_chain = None
                b = 2.0 * QV_chain * coll.inverse_denom[global_n]

            psi_face_chain = ordinate_scan(
                coll.a_attenuation[global_n], b, psi_in,
            )                                                          # (nx, ng)

            # DD spatial closure — vectorised cell-average.
            psi_face_in_chain = np.empty_like(psi_face_chain)
            psi_face_in_chain[0] = psi_in
            psi_face_in_chain[1:] = psi_face_chain[:-1]
            psi_avg_chain = 0.5 * (psi_face_in_chain + psi_face_chain)

            # M-M angular thread for curvilinear (slab: no angular state).
            if psi_angle is not None:
                psi_angle_out_chain = (
                    geom.tau_inv[global_n] * psi_avg_chain
                    - geom.mm_a_in_coeff[global_n] * psi_a_in_chain
                )
                psi_angle[chain] = psi_angle_out_chain

            # ── Scatter back to cell-index order ─────────────────────
            inv = geom.chain_idx_inv[global_n]
            psi_avg = psi_avg_chain[inv]
            angular_flux[global_n, :, 0, :] = psi_avg
            scalar_flux += w_n * psi_avg

            # ── Persist outflow at the appropriate boundary face ─────
            if is_slab:
                if direction_sign == +1:
                    bc_right_face[global_n] = psi_face_chain[-1]  # type: ignore[index]
                else:
                    bc_left_face[global_n] = psi_face_chain[-1]   # type: ignore[index]
            else:
                # Curvilinear: outward sweeps' last chain output IS the
                # outer-face outflow; inward sweeps don't write the
                # outer face.
                if mu_n >= 0 and abs(mu_n) >= sn_mesh._DEGENERATE_ABS_ETA_THRESHOLD:
                    bc_outer[global_n] = psi_face_chain[-1]      # type: ignore[index]

    # ── Cache previous-iteration state for next sweep's seeds ─────────
    if not is_slab:
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
    sn_mesh: "SNMesh",
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


__all__ = [
    "apply_sweep_1d",
    "transport_sweep",
]
