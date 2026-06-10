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
from .sweep_schedule import SweepSchedule

if TYPE_CHECKING:
    from collections.abc import Callable

    from orpheus.transport.fields.boundary_flux import BoundaryFlux
    from orpheus.numerics.projection import MomentProjection
    from .geometry import SNMesh
    from .sweep_schedule import OctantSweepGroup
    from orpheus.transport.source_sinks import ScalarSourceSink, AngularSourceSink
    from orpheus.transport.timed_full_field import TimedFullField


# ═══════════════════════════════════════════════════════════════════════
# Top-level dispatch
# ═══════════════════════════════════════════════════════════════════════


def transport_sweep(
    source: "AngularSourceSink",
    sig_t: np.ndarray,
    sn_mesh: "SNMesh",
    boundary_flux: "BoundaryFlux",
    *,
    initial_guess: "AngularFlux | TimedFullField | None" = None,
    moment_projection: "MomentProjection | None" = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Perform one full transport sweep.

    Boundary conditions are read from ``sn_mesh`` (resolved at
    construction time from the geometry mesh's :class:`BC` declarations).
    The cell-update strategy is read from ``sn_mesh.cell_update``
    (defaults to :class:`DiamondDifference`).

    Single-source contract (R-1 Step 4 A1)
    --------------------------------------

    The sweep consumes ONE :class:`AngularSourceSink` carrying the
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
    :meth:`~orpheus.sn.sources.AngularSourceSink.from_isotropic`.
    Scattering-generated sources project at the producer boundary
    via the singledispatched
    :meth:`~orpheus.sn.scattering.ScatteringOperator.apply`.  Fission-
    generated sources project at the producer boundary via
    :meth:`~orpheus.sn.fission.FissionOperator.apply`.

    The legacy two-parameter convention (``iso_source: ScalarSourceSink``
    + ``aniso_source: AngularSourceSink | None`` with sweep-internal
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
    * Issue #197 PR-TYPED-3: typed :class:`ScalarSourceSink` /
      :class:`AngularSourceSink` inputs.
    * Issue #197 PR-TYPED-4: bare-``np.ndarray`` overload retired.
    * R-1 Step 0 (2026-05-19): curvilinear Carlson seed derives from
      ``initial_guess`` (= previous iterate; ``None`` fallback uses
      the in-iteration source angular average).
    * R-1 Step 4 A1 (2026-05-21): ``iso_source`` parameter retired;
      sweep takes one :class:`AngularSourceSink` carrying the combined
      per-ordinate source magnitude.  Producer-side ``/W`` projection
      everywhere.

    Parameters
    ----------
    source : AngularSourceSink
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
        zero-initialised instance via ``BoundaryFlux.zeros_on(sn_mesh)``.
    initial_guess : AngularFlux, TimedFullField, or None, optional
        Previous-iteration angular flux estimate, used for the
        curvilinear Carlson coupled-pole seed and the per-ordinate
        spatial-upstream inflow at the pole cell.  Ignored on slab.
        ``None`` (default) selects the in-iteration fallback seed.

        Accepts both the legacy
        :class:`~orpheus.sn.angular_flux.AngularFlux` (reads
        ``.values``) and the composite
        :class:`~orpheus.transport.timed_full_field.TimedFullField`
        (reads ``.bulk.values``) via D-H.1c stage 4's
        :func:`_initial_guess_values` extractor — the kernel reads
        only the per-ordinate bulk ndarray, container-agnostic.
    moment_projection : MomentProjection or None, optional
        Phase 5c moment OUTPUT mode (2-D Cartesian ONLY — raises on a 1-D
        mesh).  ``None`` (default) returns the full per-ordinate angular flux
        (every full-angular consumer).  When given (the windowed-SI path), the
        2-D sweep accumulates the harmonic moment tensor per anti-diagonal and
        returns it INSTEAD of the angular field — the full ``(N, ng, nx, ny)``
        field is never materialized (the ~3× peak-memory win).  See
        :func:`_sweep_2d_scheduled`.

    Returns
    -------
    bulk
        ``moment_projection is None`` → ``angular_flux`` ``(N, ng, nx, ny)``.
        Given → the harmonic moment tensor ``(L+1, 2L+1, ng, nx, ny)``.
    scalar_flux
        ``(ng, nx, ny)`` :math:`\\sum_n w_n \\psi_n` in angular mode; ``None`` in
        moment mode (the scalar IS :math:`\\phi_0^0` = ``moments[0, 0]``).

    Dispatch:

    The sweep algorithm is a first-class, selectable sweep strategy
    (``orpheus.sn.sweep_strategy.SweepStrategy``; ``default_for`` picks the
    default for the mesh).  This replaces the historical scattered branch —
    the ``reduced is not None`` test here and the five ``not is_1d`` gates in
    the operator algebra — with one polymorphic selection.  ``default_for``
    reproduces the legacy choice exactly:

    * 1-D meshes → ``CumprodScan`` (:func:`_sweep_1d_unified`; slab, sphere,
      cylinder — one body via the two-stratum cache).
    * 2-D Cartesian → ``MovingFrontierWindow`` (:func:`_sweep_2d_wavefront`;
      anti-diagonal scheduling).

    The ``moment_projection`` guard (moment output is 2-D Cartesian only)
    now lives in ``CumprodScan.sweep`` — the strategy that cannot produce it
    carries its own refusal.

    (The strategy architecture's rendered API + theory page land in
    Sweep-Strategy carve phase S5; this docstring names the symbols as
    literals until then.)
    """
    Q = _unwrap_source(source)
    # Lazy import breaks the sweep ↔ sweep_strategy cycle: the strategy
    # module wraps the ``_sweep_*`` functions defined here, so it imports
    # ``sweep`` at module load; ``transport_sweep`` reaches back for the
    # selector only at call time.
    from .sweep_strategy import default_for

    return default_for(sn_mesh).sweep(
        Q, sig_t, boundary_flux,
        initial_guess=initial_guess,
        moment_projection=moment_projection,
    )


def _unwrap_source(source: "AngularSourceSink") -> np.ndarray:
    """Unwrap typed :class:`AngularSourceSink` to bare ndarray.

    Issue #197 PR-TYPED-4 — strict typed input.  R-1 Step 4 A1 collapsed
    the iso / aniso parameter pair into a single per-ordinate source.
    The internal hot path consumes bare ndarray; this helper performs
    the unwrap once at the public boundary.
    """
    from orpheus.transport.source_sinks import AngularSourceSink
    if not isinstance(source, AngularSourceSink):
        raise TypeError(
            f"transport_sweep: source must be "
            f"AngularSourceSink (R-1 Step 4 A1); got "
            f"{type(source).__name__}"
        )
    return source.values


def _initial_guess_values(
    initial_guess: "AngularFlux | TimedFullField | None",
) -> "np.ndarray | None":
    """Extract the per-ordinate bulk ndarray from either container type.

    D-H.1c stage 4 — the sweep kernel reads only the per-ordinate bulk
    values (shape ``(N, ng, nx, ny)``) at two sites: the M-M Carlson
    coupled-pole seed (transposed slice at level p) and the pole-cell
    spatial-upstream inflow (single-cell slice at ordinate global_n).

    Both legacy :class:`~orpheus.sn.angular_flux.AngularFlux` and
    composite :class:`~orpheus.transport.timed_full_field.TimedFullField`
    carry the same ndarray under different attribute paths.  This
    helper centralises the access so the kernel stays
    container-agnostic.

    Parameters
    ----------
    initial_guess : AngularFlux, TimedFullField, or None
        Container carrying the previous iterate, OR ``None`` for
        cold-start.

    Returns
    -------
    np.ndarray or None
        The per-ordinate ``(N, ng, nx, ny)`` ndarray, or ``None``
        when ``initial_guess`` is ``None``.

    Notes
    -----
    Uses duck-typing on ``.bulk`` to detect the composite — avoids a
    runtime import of :class:`TimedFullField` (which would create a
    circular-dependency risk through transport↔sn).
    """
    if initial_guess is None:
        return None
    # Composite container exposes ``.bulk`` (an L2 AngularFlux); the
    # legacy bundle does not.
    bulk = getattr(initial_guess, "bulk", None)
    if bulk is not None:
        return bulk.values
    return initial_guess.values  # type: ignore[union-attr]


# ═══════════════════════════════════════════════════════════════════════
# 1-D unified sweep — slab + sphere + cylinder via the cache
# ═══════════════════════════════════════════════════════════════════════


def _sweep_1d_unified(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: "SNMesh",
    boundary_flux: "BoundaryFlux",
    *,
    initial_guess: "AngularFlux | TimedFullField | None" = None,
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
        # 1-D meshes: sig_t is the principled (ng, nx) layout the cache
        # expects natively (rank-d (N, ng, *spatial); no phantom ny axis).
        sig_t_1d = sig_t  # (ng, nx)
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
    # NOTE: ``initial_guess`` typing widens to also accept
    # :class:`TimedFullField` after D-H.1c stage 4; the container-
    # agnostic extractor :func:`_initial_guess_values` centralises the
    # read.
    *,
    initial_guess: "AngularFlux | TimedFullField | None" = None,
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

    # ── Entry layout — the public contract is the principled
    # (N, ng, *spatial) = (N, ng, nx) for 1-D (no phantom ny axis).
    Q_per_ord = Q                                            # (N, ng, nx)
    sig_t_p = sig_t                                          # (ng, nx)
    V = sn_mesh.volumes                                      # (nx,) — no group axis
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
    angular_flux = np.zeros((N, ng, nx))
    scalar_flux = np.zeros((ng, nx))

    # ── BC inflow + per-level Carlson seed (curvilinear only) ─────────
    #
    # Wave O (#208) O.4a.2 — BARE SWEEP: the entry ``bc_*.apply`` is GONE.
    # The reflective coupling ``ψ.inflow = B·ψ.outflow`` is no longer
    # re-applied inside the sweep; it is supplied by the CALLER as the
    # ``−B`` source term (the SI driver folds ``S + B`` into the source;
    # the direct fixed-source loops + the final reconstruction reflect the
    # persisted outflow into the inflow slots via ``SNBoundaryOperator``
    # before each sweep — see ``solver.py``). The sweep now reads the
    # SEEDED inflow trace DIRECTLY: the incoming-ordinate slots of the
    # face view ARE the inflow seed, and the outgoing-ordinate slots are
    # persisted in place after the sweep. Reading the inflow ords (before)
    # and writing the outflow ords (after) touch DISJOINT ordinate sets,
    # so aliasing the face view is safe.
    if is_slab:
        # D-H.2-C2: L2 :class:`BoundaryFlux` provides writable per-face
        # views via :meth:`face_view`.  Slab layout has both ``xmin``
        # and ``xmax`` slots (shape ``(N, ng)`` each); writes through
        # the view propagate to the flat backing buffer.  Per-cell-call
        # outflow persistence below (``bc_left_face[ords] = ...``)
        # mutates these views in place.
        bc_left_face = boundary_flux.face_view("xmin")   # (N, ng)
        bc_right_face = boundary_flux.face_view("xmax")  # (N, ng)
        inflow_left = bc_left_face    # incoming-ord slots = seeded inflow
        inflow_right = bc_right_face  # incoming-ord slots = seeded inflow
        levels = [None]
        level_ordinates_list = [list(range(N))]
        bc_outer = None
        sigma_t_gx = None
        dr = None
        inflow_full = None
    else:
        # D-H.2-C2: 1-D curvilinear layout has only the outer radial
        # ``xmax`` face (the geometric pole at r=0 is a regularity
        # condition, not a BC face).  Writable view into the L2 flat
        # backing buffer.
        bc_outer = boundary_flux.face_view("xmax")  # (N, ng)
        inflow_full = bc_outer  # incoming-ord slots = seeded inflow (bare sweep)

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
            angular_flux[ords, :, :] = psi_avg_cell_order

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
        # D-H.1c stage 4 — container-agnostic bulk extraction (works for
        # legacy AngularFlux ``.values`` and composite TimedFullField
        # ``.bulk.values`` identically).
        ig_values = _initial_guess_values(initial_guess)
        if ig_values is not None:
            # (N, ng, nx) → (ng, N, nx)
            psi_g_first = ig_values.swapaxes(0, 1)
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
                    # D-H.1c stage 4 — ``ig_values`` carries the bulk
                    # ndarray regardless of container type (legacy
                    # AngularFlux or composite TimedFullField).
                    if ig_values is not None:
                        psi_in = ig_values[global_n, :, 0]
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
                        angular_flux[global_n, :, i] = psi
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
                angular_flux[global_n, :, :] = psi_avg_p
                scalar_flux += w_n * psi_avg_p

                # Persist outflow at the outer face for outward ordinates.
                if mu_n >= 0 and abs(mu_n) >= sn_mesh._DEGENERATE_ABS_ETA_THRESHOLD:
                    bc_outer[global_n] = psi_face_chain[-1]      # (ng,)

    # ── Exit — PR-INDEX-5: caller consumes principled layout ──────────
    # R-1 Step 0: NO iteration-cache write-back.  The caller threads the
    # returned ``angular_flux`` as ``initial_guess`` to the NEXT sweep —
    # that's all the "history" needed.  The matvec already operates this
    # way; the sweep now mirrors it.
    # angular_flux is (N, ng, nx); scalar_flux is (ng, nx) — the principled
    # (N, ng, *spatial) / (ng, *spatial) public contract (no phantom ny).
    return angular_flux, scalar_flux


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

def sweep_octant_group(
    group: "OctantSweepGroup",
    *,
    boundary_flux: "BoundaryFlux",
    Q: np.ndarray,
    sig_t: np.ndarray,
    str_x: np.ndarray,
    str_y: np.ndarray,
    weights: np.ndarray,
    sn_mesh: "SNMesh",
    angular_flux: np.ndarray | None = None,
    scalar_flux: np.ndarray | None = None,
    moment_buf: np.ndarray | None = None,
    Y: np.ndarray | None = None,
) -> None:
    r"""Sweep one octant GROUP's ordinates on the rolling moving-frontier
    window, in place (Phase 3 sub-step 3b primitive; Phase 5b storage-B).

    For each :class:`~orpheus.sn.sweep_schedule.OctantSweep` in ``group``,
    dispatch to the precomputed ``SNMesh.sweep_graphs[OctantLabel(σ)]`` and
    drive :meth:`SweepDependencyGraph.apply_windowed` — the rolling 2-diagonal
    window walk that carries the interior face cochain in an O(N·ng·nx) buffer
    instead of the full O(N·ng·nx·ny) per-axis field (Phase 5b storage-B). The
    octant's domain-edge INFLOW is read from ``boundary_flux`` (its incoming
    faces); the OUTFLOW is shed back into ``boundary_flux`` (its outgoing
    faces) as the frontier advances.

    Boundary coupling via ``boundary_flux`` (replacing the storage-A
    whole-trace seed/absorb on a persistent ``WavefrontFlux``): each octant
    reads ``boundary_flux.face_view(<incoming face>)[oct_idx]`` as its inflow
    and writes ``boundary_flux.face_view(<outgoing face>)[oct_idx]`` with its
    shed outflow. Because distinct octants own DISJOINT ordinate slices of a
    face, an octant's outflow write never clobbers another octant's inflow —
    so the Jacobi single-group call is bit-identical to the legacy per-octant
    loop, and the inflowing ordinates retain their seed value exactly as the
    storage-A ``absorb`` left them. The Gauss-Seidel schedule calls this once
    per in-plane octant group, with the resolvent reflecting the just-shed
    outflow between groups so a later group reads the fresh current-iterate
    inflow off the SAME ``boundary_flux`` (the ``(L+C−B_lower)⁻¹`` forward
    substitution).

    Parameters
    ----------
    group : OctantSweepGroup
        The octants to sweep (each an :class:`OctantSweep`: in-plane
        :class:`~orpheus.sn.sweep_graph.OctantLabel` + ordinate indices).
        ``group.reflect_faces`` is consumed by the resolvent, NOT here — this
        function is the bare sweep, blind to the boundary coupling.
    boundary_flux : BoundaryFlux
        The per-ordinate boundary trace. Carries each octant's inflow on its
        incoming faces (seeded by the caller / freshly reflected by an earlier
        G-S group) and receives the shed outflow on its outgoing faces.
    Q, sig_t, str_x, str_y, weights : np.ndarray
        Per-ordinate source ``(N, ng, nx, ny)``, total XS ``(ng, nx, ny)``,
        the streaming stencils ``(N, nx)`` / ``(N, ny)``, quadrature weights
        ``(N,)`` — sliced per octant by the ordinate indices.
    sn_mesh : SNMesh
        Carries the per-octant ``sweep_graphs`` and the ``cell_update``
        strategy.
    angular_flux, scalar_flux : np.ndarray, optional
        ANGULAR-mode output buffers ``(N, ng, nx, ny)`` / ``(ng, nx, ny)`` —
        mutated in place (the angular flux scattered per octant, the scalar flux
        accumulated ``Σ_n w_n ψ_n``). Given iff ``moment_buf`` is ``None``.
    moment_buf, Y : np.ndarray, optional
        MOMENT-mode output (Phase 5c): the harmonic moment accumulator
        ``(L+1, 2L+1, ng, nx, ny)`` and the full harmonic table
        ``(N, L+1, 2L+1)`` (octant-sliced ``Y[oct_idx]`` per octant). Given iff
        ``angular_flux`` is ``None``. The windowed SI iterate is moments (5a), so
        the full per-ordinate angular field is never materialized — the walk
        projects each anti-diagonal directly into ``moment_buf``
        (:math:`\phi_\ell^m \mathrel{+}= \sum_n w_n Y_\ell^m \psi_n`). The
        cross-octant ``+=`` (octants share the moment buffer) reorders the
        ordinate sum vs the post-sweep flat projection ⇒ principled-equivalence,
        NOT bit-identity. The scalar is subsumed (``moment_buf[0, 0]``).
    """
    ng, nx, ny = sig_t.shape
    cell_update = sn_mesh.cell_update

    for sweep in group.sweeps:
        # The schedule already projected the quadrature octant to its in-plane
        # OctantLabel (sign_z dropped, mirroring the legacy ``sx = label[0];
        # sy = label[1] if len>=2 else 0``); convert the ordinate-index tuple
        # back to an array for value-based fancy indexing (a bare tuple would
        # be read as a multi-axis index — NOT what we want).
        oct_idx = np.asarray(sweep.indices)   # (N_oct,) int into N
        sx = sweep.label.sign_x
        sy = sweep.label.sign_y

        # Pure-z degenerate octant: no in-plane streaming. The angular flux is
        # the volumetric balance ψ = Q_n / Σ_t and the scalar flux gets a
        # weighted contribution. No faces, no boundary interaction.
        if sx == 0 and sy == 0:
            Q_pure_z = Q[oct_idx]                         # (N_oct, ng, nx, ny)
            # sig_t (ng, nx, ny) broadcasts against (N_oct, ng, nx, ny).
            psi_avg_pure_z = Q_pure_z / sig_t              # (N_oct, ng, nx, ny)
            if moment_buf is None:
                angular_flux[oct_idx] = psi_avg_pure_z
                scalar_flux += np.einsum(
                    "ngij,n->gij",
                    psi_avg_pure_z, weights[oct_idx],
                )
            else:
                moment_buf += np.einsum(
                    "nlm,ngij,n->lmgij",
                    Y[oct_idx], psi_avg_pure_z, weights[oct_idx],
                )
            continue

        # Effective in-plane sign for sweep-graph lookup. Match legacy's
        # ``key = (1 if mx >= 0 else -1, ...)`` mapping: ordinates with
        # ``mx == 0`` are treated as ``+1`` (the streaming coefficient is
        # zero, and the WDD result is identical regardless of sign choice).
        sx_eff = +1 if sx == 0 else sx
        sy_eff = +1 if sy == 0 else sy
        sweep_graph = sn_mesh.sweep_graphs[OctantLabel((sx_eff, sy_eff))]
        N_oct = oct_idx.size

        # ── Domain-edge inflow/outflow faces for this octant ──────────
        # The incoming face is the LOW face of the sweep direction (xmin for
        # +x, xmax for −x); the outgoing face is the opposite. The octant
        # streams from inflow → outflow; ``boundary_flux`` carries both.
        x_in_face = "xmin" if sx_eff >= 0 else "xmax"
        x_out_face = "xmax" if sx_eff >= 0 else "xmin"
        y_in_face = "ymin" if sy_eff >= 0 else "ymax"
        y_out_face = "ymax" if sy_eff >= 0 else "ymin"

        # Read this octant's inflow (fancy-index → a (N_oct, ng, ·) copy).
        inflow_x = boundary_flux.face_view(x_in_face)[oct_idx]   # (N_oct, ng, ny)
        inflow_y = boundary_flux.face_view(y_in_face)[oct_idx]   # (N_oct, ng, nx)

        capture_x = np.empty((N_oct, ng, ny))   # shed domain x-outflow
        capture_y = np.empty((N_oct, ng, nx))   # shed domain y-outflow

        # Angular mode allocates a per-octant angular buffer (scattered into the
        # global field below); moment mode accumulates directly into the shared
        # moment buffer per anti-diagonal, so NO per-octant angular field is
        # materialized (the Phase 5c peak-memory win).
        angular_flux_oct = (
            np.zeros((N_oct, ng, nx, ny)) if moment_buf is None else None
        )
        sweep_graph.apply_windowed(
            cell_update=cell_update,
            inflow_x=inflow_x,
            inflow_y=inflow_y,
            Q_octant=Q[oct_idx],
            sig_t=sig_t,
            str_x_octant=str_x[oct_idx],
            str_y_octant=str_y[oct_idx],
            weights_octant=weights[oct_idx],
            capture_x=capture_x,
            capture_y=capture_y,
            angular_flux_octant=angular_flux_oct,
            scalar_flux_buf=scalar_flux,
            moment_buf=moment_buf,
            Y_octant=None if Y is None else Y[oct_idx],
        )

        # Shed the outflow into the boundary trace. The outflow write touches
        # only this octant's ordinates (disjoint from other octants' inflow
        # ordinates), so it is ι*-faithful. Angular mode also scatters the
        # per-octant angular field into the global buffer.
        if moment_buf is None:
            angular_flux[oct_idx] = angular_flux_oct
        boundary_flux.face_view(x_out_face)[oct_idx] = capture_x
        boundary_flux.face_view(y_out_face)[oct_idx] = capture_y


def _sweep_2d_scheduled(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: "SNMesh",
    boundary_flux: "BoundaryFlux",
    *,
    schedule: "SweepSchedule",
    reflect: "Callable[[BoundaryFlux, tuple[str, ...]], None] | None" = None,
    moment_projection: "MomentProjection | None" = None,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Polymorphic schedule-driven 2-D wavefront sweep (Phase 3 sub-step 3c).

    ONE uniform sweep-and-reflect loop parameterized by ``schedule`` — the
    realization of the polymorphic Jacobi / Gauss-Seidel design (there is NO
    ``if jacobi/gs`` branch; the splitting IS the schedule):

    1. The GIVEN inflow ``boundary_flux`` (carrying ``B·ψₙ`` — the lagged
       reflection of the previous iterate, prepared by the caller) is read
       per-octant by the window walk; there is no separate whole-trace seed.
    2. ``for group in schedule.groups``: :func:`sweep_octant_group` the group's
       octants on the rolling window, sheding each octant's outflow into
       ``boundary_flux`` (the ι* absorb, now per-octant during the walk). Then
       if ``reflect`` is given AND the group has reflective outgoing faces,
       apply ``reflect`` (the ``−B`` outflow→inflow reflection, in place,
       face-restricted) so a LATER group reads the fresh current-iterate inflow
       directly off ``boundary_flux`` (the ``(L+C−B_lower)⁻¹`` forward
       substitution) — no re-seed needed, the next walk reads the trace fresh.

    * **Jacobi** (``reflect=None``, one all-octants group) — every octant reads
      the frozen seed; the inter-group reflect never fires. This is exactly the
      bare sweep :func:`_sweep_2d_wavefront` passes.
    * **Gauss-Seidel** (``reflect`` = the face-restricted ``−B``, one group per
      in-plane octant) — later groups see earlier groups' fresh reflected
      outflow. The SI scheduled resolvent supplies both: its ``.solve`` seeds
      ``B·ψₙ`` onto ``boundary_flux`` then calls this with the G-S schedule +
      the reflect closure. The walk's per-octant shed populates the fresh
      outflow that ``reflect`` then maps to inflow (replacing the storage-A
      ``absorb``-before-``reflect`` step).

    The converged fixed point is INVARIANT under ``schedule`` (any consistent
    splitting of ``(L+C−S−B)ψ=q`` shares ψ\*); only the SI spectral rate
    changes. NOTE (Phase 3 spike, issue #2/#215): this folds the BOUNDARY
    coupling ``B`` only — a modest reflective-SI rate gain. The dominant
    within-group SCATTERING ``c``-mode is NOT folded here (it cannot be folded
    into a directional sweep); that is consistent DSA / Krylov territory.

    Phase 5b storage-B: the interior face cochain is a rolling 2-diagonal
    moving-frontier window (O(N·ng·nx)), carried inside
    :meth:`SweepDependencyGraph.apply_windowed` per octant — NOT the full
    O(N·ng·nx·ny) per-axis field. The full-field walk
    (:meth:`SweepDependencyGraph.apply` on a ``WavefrontFlux``) is retained as
    the bit-identity verification oracle (see the ``window ≡ full-field``
    test); the converged solution is unchanged.

    Phase 5c moment output: when ``moment_projection`` is given (the 2-D
    Cartesian windowed-SI path), the walk accumulates the harmonic moment tensor
    ``(L+1, 2L+1, ng, nx, ny)`` per anti-diagonal directly — the full
    per-ordinate angular OUTPUT ``(N, ng, nx, ny)`` is never materialized (the
    ~3× linear peak-memory win; the persistent SI iterate is already moments,
    5a).  Returns ``(moment_buf, None)`` — the scalar flux is :math:`\phi_0^0`
    = ``moment_buf[0, 0]`` (``Y_0^0 = 1``), read off the tensor, NOT returned
    separately (the angular-mode scalar is an independent array; ``None`` keeps
    the modes' second slot from being mistaken).  Principled-equivalence, NOT
    bit-identity: the cross-octant ``+=`` reorders the ordinate sum vs the
    post-sweep flat :class:`~orpheus.numerics.projection.MomentProjection`
    reduce (≤ 4 ULP de-risk).  ``moment_projection is None`` (every full-angular
    consumer — reconstruction, Krylov, 1-D) returns ``(angular_flux,
    scalar_flux)`` exactly as before.
    """
    ng, nx, ny = sig_t.shape
    N = sn_mesh.quad.N
    if moment_projection is None:
        angular_flux = np.zeros((N, ng, nx, ny))
        scalar_flux = np.zeros((ng, nx, ny))
        moment_buf = None
        Y = None
    else:
        L = moment_projection.L
        moment_buf = np.zeros((L + 1, 2 * L + 1, ng, nx, ny))
        Y = moment_projection.Y
        angular_flux = None
        scalar_flux = None

    str_x = sn_mesh.streaming_x   # (N, nx)
    str_y = sn_mesh.streaming_y   # (N, ny)
    weights = sn_mesh.quad.weights

    for group in schedule.groups:
        sweep_octant_group(
            group,
            boundary_flux=boundary_flux,
            Q=Q,
            sig_t=sig_t,
            str_x=str_x,
            str_y=str_y,
            weights=weights,
            sn_mesh=sn_mesh,
            angular_flux=angular_flux,
            scalar_flux=scalar_flux,
            moment_buf=moment_buf,
            Y=Y,
        )
        if reflect is not None and group.reflect_faces:
            # G-S inter-group reflect (a no-op for the Jacobi schedule, whose
            # sole group carries no reflect_faces): the walk already shed this
            # group's fresh outflow into ``boundary_flux``; reflect ``−B``
            # (outflow→inflow, in place, face-restricted) so the NEXT group
            # reads the fresh current-iterate reflected inflow off the trace.
            reflect(boundary_flux, group.reflect_faces)

    if moment_projection is None:
        return angular_flux, scalar_flux
    # Moment mode: (moments, None).  The scalar IS φ_0^0 = ``moments[0, 0]``
    # (Y_0^0 = 1), read off the tensor by the caller — NOT returned separately
    # (returning the live ``moment_buf[0, 0]`` view invites aliasing; the
    # angular-mode scalar is an independent array, so a None here keeps the two
    # modes' second slot from being mistaken for the same kind of value).
    return moment_buf, None


def _sweep_2d_wavefront(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: "SNMesh",
    boundary_flux: "BoundaryFlux",
    *,
    moment_projection: "MomentProjection | None" = None,
) -> tuple[np.ndarray, np.ndarray]:
    r"""2-D wavefront sweep — per-octant batched (Wave 2 / C2.6).

    Iterates over angular **octants** (lexicographic order from
    :attr:`AngularQuadrature.octants`). For each in-plane octant
    :math:`\sigma = (\mathrm{sign}\,\mu_x, \mathrm{sign}\,\mu_y)`:

    1. **BARE inflow seed** (Wave O #208 O.4b Phase E1) — the
       octant-incoming face slot is seeded from the GIVEN inflow trace
       ``boundary_flux.face_view(...)``; there is NO ``bc.apply``. The
       reflective coupling ``ψ.inflow = B·ψ.outflow`` is delivered
       externally by ``_reflect_outflow_into_inflow`` / the sibling
       ``-B`` between sweeps (mirroring the 1-D O.4a.2 bare sweep), so
       the sweep is the pure bulk solve ``ψ = (L+C)^{-1} q`` reading the
       inflow as a fixed boundary datum.
    2. **Dispatch to ``SNMesh.sweep_graphs[OctantLabel(σ)]``** —
       the per-octant ``SweepDependencyGraph`` precomputed once at
       mesh construction (Wave 2 / C2.4).
    3. The graph's ``apply`` walks topological levels (anti-diagonals)
       and dispatches each level to
       ``sn_mesh.cell_update.update_batch(slice_args)`` — vectorised
       over ``(N_oct, n_diag, ng)`` simultaneously.

    The smoking gun ``for n in range(N)`` is gone: the outer loop is now
    ``for group in SweepSchedule.jacobi(sn_mesh).groups`` (ONE all-octants
    group), and the per-octant work is delegated to
    :func:`sweep_octant_group` — one ``SweepDependencyGraph.apply`` per octant
    over the typed interior cochain
    :class:`~orpheus.transport.fields.wavefront_flux.WavefrontFlux`
    (Wave O #205); the boundary trace is exchanged once per sweep via the
    typed :math:`\iota_*` (``seed``) / :math:`\iota^*` (``absorb``). The
    ordinate axis is INTERNAL to every numpy operation.

    Phase 3 sub-step 3b — this body IS the **Jacobi** octant schedule (one
    group, all octants, no inter-group reflect). The carve into
    :func:`sweep_octant_group` + the polymorphic
    :class:`~orpheus.sn.sweep_schedule.SweepSchedule` lets the SI resolvent's
    **Gauss-Seidel** schedule (sub-step 3c) re-use the SAME per-octant sweep at
    finer group granularity, recovering the intra-sweep reflective coupling
    Wave O externalised — while THIS default sweep stays bit-identical Jacobi.

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
    # The bare 2-D sweep IS the JACOBI octant schedule (Phase 3 sub-step 3c):
    # ONE group (all octants), NO inter-group reflect. Delegates to the
    # polymorphic :func:`_sweep_2d_scheduled` with ``reflect=None`` — all
    # octants read the same frozen inflow seed (``boundary_flux``, carrying
    # ``B·ψₙ`` for the SI driver via ``rhs.boundary``); the inter-group reflect
    # never fires (the Jacobi group carries no ``reflect_faces``). The
    # Gauss-Seidel SI resolvent calls the SAME orchestrator with the
    # per-in-plane-octant schedule + a ``−B`` reflect closure, recovering the
    # boundary reflective coupling Wave O externalised. Bit-identical to the
    # pre-3c ``for octant in quad.octants`` loop.
    return _sweep_2d_scheduled(
        Q, sig_t, sn_mesh, boundary_flux,
        schedule=SweepSchedule.jacobi(sn_mesh),
        reflect=None,
        moment_projection=moment_projection,
    )


# ═══════════════════════════════════════════════════════════════════════
# VERIFICATION ORACLE — the full-field DAG-walk sweep (d-generic; NOT production)
# ═══════════════════════════════════════════════════════════════════════


def _sweep_full_field(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: "SNMesh",
    boundary_flux: "BoundaryFlux",
) -> tuple[np.ndarray, np.ndarray]:
    r"""Full-field DAG-walk Jacobi sweep — the dimension-generic VERIFICATION
    ORACLE (:class:`~orpheus.sn.sweep_strategy.FullFieldWavefront`).

    The dimension-agnostic SPINE: ONE body for d = 1 (slab) and d = 2
    (Cartesian), and the verification reference the d-specific production
    optimizations are cross-checked against — the 1-D
    :class:`~orpheus.sn.sweep_strategy.CumprodScan` (principled-equivalence at
    nulp) and the 2-D :class:`~orpheus.sn.sweep_strategy.MovingFrontierWindow`
    (``window ≡ full`` bit-identity). It carries the FULL interior face cochain
    (the storage-A path Phase 5b superseded in production); its sole purpose is
    verification. Three reasons it is a legitimate kept reference (not a
    twin-path smell):

    1. It carries the FULL interior face cochain as the typed
       :class:`~orpheus.transport.fields.wavefront_flux.WavefrontFlux` (the
       fuller view — all cells' face flux, not just the rolling frontier),
       with the :math:`\iota_*` (``seed``) / :math:`\iota^*` (``absorb``)
       whole-trace boundary algebra.
    2. It walks via :meth:`SweepDependencyGraph.apply` — the full-field walk
       that shares the SAME cell kernel
       (:meth:`DiamondDifference.cell_kernel_batch`) as the windowed
       :meth:`~SweepDependencyGraph.apply_windowed`. So the cell MATH cannot
       drift between oracle and production; only the storage walk differs —
       which is exactly what the ``window ≡ full`` equivalence test pins
       (``np.array_equal`` at d = 2).
    3. Jacobi only (one all-octants group, no inter-group reflect): the default
       :func:`_sweep_2d_wavefront` path. The window≡full bit-identity is
       per-octant (schedule-independent), so the Jacobi oracle suffices; the
       Gauss-Seidel windowed path is exercised by the eigenvalue tier.

    Dimension-generic via the wavefront DAG (B7-verified d = 1 / d = 3): the
    per-axis tuples (faces, streaming, octant signs) are built over
    ``range(ndim)`` — at d = 2 byte-for-byte the legacy ``(sign_x, sign_y)`` /
    ``(streaming_x, streaming_y)`` path (the window anchor pins this); at d = 1
    a single-axis chain walk. The schedule's octant labels live in the full
    angular space (``signs`` may carry a phantom out-of-plane component, e.g.
    ``(-1, 0)`` on a 1-D mesh), so they are projected to ``signs[:ndim]`` —
    which at d = 2 is the in-plane ``(sign_x, sign_y)`` the 2-D graphs are keyed
    by, and at d = 1 the single in-plane sign the chain graphs are keyed by.
    """
    from orpheus.transport.fields.wavefront_flux import WavefrontFlux

    ng = sig_t.shape[0]
    spatial = sig_t.shape[1:]                            # (nx,) at d=1; (nx, ny) at d=2
    ndim = sn_mesh.ndim
    N = sn_mesh.quad.N
    angular_flux = np.zeros((N, ng, *spatial))
    scalar_flux = np.zeros((ng, *spatial))

    # The FULL interior cochain (typed) — ι_*-seeded whole-trace from the given
    # inflow, ι*-absorbed whole-trace at the end. This is the fuller view the
    # window replaces with a rolling (d-1)-frontier in production.
    wavefront = WavefrontFlux.zeros_on(sn_mesh)
    wavefront.seed(boundary_flux)
    # Typed axis→buffer map: the per-axis face fields are taken via
    # ``WavefrontFlux.face(a)`` over ``WavefrontFlux.axes``; the streaming
    # tuple is the axis-keyed ``sn_mesh.streaming(a)`` map over the SAME
    # ``range(ndim)`` — so ``str_axes[a]`` cannot drift from ``psi_faces[a]``
    # (axis a's streaming coeff pairs with axis a's face in the cell kernel).
    psi_faces = tuple(wavefront.face(a) for a in wavefront.axes)
    str_axes = tuple(sn_mesh.streaming(a) for a in range(ndim))
    weights = sn_mesh.quad.weights
    cell_update = sn_mesh.cell_update

    for group in SweepSchedule.jacobi(sn_mesh).groups:
        for sweep in group.sweeps:
            oct_idx = np.asarray(sweep.indices)
            signs = sweep.label.signs[:ndim]             # in-plane projection
            if not any(signs):                           # pure-z degenerate (d≥2 only)
                psi_avg_pure_z = Q[oct_idx] / sig_t
                angular_flux[oct_idx] = psi_avg_pure_z
                scalar_flux += np.einsum(
                    "ng...,n->g...", psi_avg_pure_z, weights[oct_idx],
                )
                continue
            signs_eff = tuple(+1 if s == 0 else s for s in signs)
            graph = sn_mesh.sweep_graphs[OctantLabel(signs_eff)]
            # FULL per-octant face fields — octant-restricted working copies.
            psi_faces_oct = tuple(pf[oct_idx].copy() for pf in psi_faces)
            angular_flux_oct = np.zeros((oct_idx.size, ng, *spatial))
            graph.apply(                                 # the full-field walk
                cell_update=cell_update,
                psi_faces_octant=psi_faces_oct,
                Q_octant=Q[oct_idx], sig_t=sig_t,
                str_axes_octant=tuple(s[oct_idx] for s in str_axes),
                weights_octant=weights[oct_idx],
                angular_flux_octant=angular_flux_oct, scalar_flux_buf=scalar_flux,
            )
            for a in wavefront.axes:
                psi_faces[a][oct_idx] = psi_faces_oct[a]
            angular_flux[oct_idx] = angular_flux_oct

    wavefront.absorb(boundary_flux)                      # ι* whole-trace
    return angular_flux, scalar_flux


__all__ = [
    "transport_sweep",
    "sweep_octant_group",
]
