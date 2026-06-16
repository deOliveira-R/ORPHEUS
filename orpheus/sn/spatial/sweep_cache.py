r"""Two-stratum precomputed cache for 1-D SN sweeps (Issue #196 Phase G Step 2.5c).

Issue #196 Phase G Step 2.5c factors the per-ordinate-per-cell DD coefficients
into TWO frozen dataclasses by **mutation cadence**, eliminating per-sweep
:func:`np.fromiter` plus per-cell ``StreamingTerms`` dataclass construction
from the hot path:

* :class:`GeometryCoefficients` — Stratum 1, geometry+quadrature only.  Built
  once at :class:`~orpheus.sn.solver.SNSolver` construction.  No ``ng`` axis;
  invariant under cross-section rebinds, BC changes, every outer / inner /
  Picard iteration.  Lifetime = ``SNMesh`` × ``AngularQuadrature``.

* :class:`CollisionCache` — Stratum 2, geometry × :math:`\Sigma_t`.  Built
  when :math:`\Sigma_t` is bound.  Lifetime = constant-:math:`\Sigma_t` epoch
  (one fixed-source solve, one eigenvalue solve, one depletion micro-step).
  Rebuilt by :meth:`~orpheus.sn.solver.SNSolver.rebind_cross_sections`.

The strata are deliberately separate (cross-domain-attacker Smell #16): mixing
geometry-only and :math:`\Sigma_t`-dependent fields into one tensor ties the
cache lifetime to its most-mutable member.  Depletion / thermal-feedback
micro-steps that change :math:`\Sigma_t` then force a full geometry-cache
rebuild for zero gain.

The math notation
=================

For one ordinate :math:`n`, one cell :math:`i` (in chain order), one group
:math:`g`.  Storage layout: ``(N, ng, nx)`` — energy second, cell trailing
(Issue #196 PR-INDEX-2 principled layout).

.. math::

   \mathrm{denom}[n, g, i]
       &\;=\; 2|\mu_n|\,A_{\rm down}[n, i]
              + \frac{\Delta A[n, i]}{w_n}\,c_{\rm out}[n]
              + \Sigma_t[g, i]\,V[i],\\
   a[n, g, i]
       &\;=\; \frac{2|\mu_n|\,A_{\rm total}[i]}{\mathrm{denom}[n, g, i]} - 1,\\
   b[n, g, i]
       &\;=\; \frac{2\,\bigl(q[n, g, i] + (\Delta A[n, i]/w_n)\,c_{\rm in}[n]
                              \,\psi^{a}_{\rm in}[n, g, i]\bigr)}
                   {\mathrm{denom}[n, g, i]},\\
   \psi^{s}[n, g, i+1]
       &\;=\; a[n, g, i]\,\psi^{s}[n, g, i] + b[n, g, i].

Closed-form scan (Blelloch §1.5; :func:`~orpheus.sn.spatial.scan.ordinate_scan`):

.. math::

   \psi^{s}[n, g, i]
       \;=\; \mathrm{cumprod\_a}[n, g, i]\,
             \bigl(\psi^{s}_0[n, g] + \mathrm{cumsum}_{k \le i}
                   (b[n, g, k] / \mathrm{cumprod\_a}[n, g, k])\bigr).

The :class:`GeometryCoefficients` populator hoists the geometry tensors
(``A_down``, ``A_total``, ``dA_w``, ``V``, ``c_in``, ``c_out``, chain
indices, M-M closure constants ``mm_a_in_coeff`` and ``tau_inv``) ONCE at
solver construction.  The :class:`CollisionCache` populator combines them
with :math:`\Sigma_t` to produce ``inverse_denom``, ``a_attenuation``, and
``cumprod_a`` — three numpy ops per ordinate.

Slab degeneracy (one cache, all geometries)
===========================================

The cache fields are populated for slab with neutral curvature values that
make the curvilinear formula degenerate to the slab form:

* ``A_total = 2`` (``face_area_inner + face_area_outer = 1 + 1``);
* ``A_down = 1`` (the neutral downstream face area for slab);
* ``dA_w = 0`` (no curvature redistribution);
* ``c_in = c_out = 0``, ``tau_inv = 1``, ``mm_a_in_coeff = 0``.

With these values the curvilinear ``denom`` collapses to ``2|μ| + σ_t·V``
(the slab form), ``a = 2·2|μ|/denom - 1 = (2|μ| - σ_t·V)/(2|μ| + σ_t·V)``,
and ``b = 2·q/denom``.  The same three tensor ops in
:meth:`~orpheus.sn.loss_representation._OneDimScanWalk.sweep` drive slab, sphere, and cylinder.

Pattern 2 anchor — single source of truth
=========================================

The cache populator derives from :func:`~orpheus.sn.spatial.cell_balance.cell_balance_terms`
indirectly via the per-cell :class:`~orpheus.geometry.reduced_operator.StreamingTerms`
factory.  The L1 dual-view validator
(``test_cache_populator_matches_cell_balance_terms``) asserts that for any
cell, the cache's ``(a, denom)`` agree with the per-cell ``cell_balance_terms``
output to ``rtol=1e-14``.  This is the Pattern 2 contract that keeps the
fast cache-driven path and the slow per-cell-update reference path consistent.

References
==========

* Lewis, E. E., & Miller, W. F. (1984). *Computational Methods of Neutron
  Transport*.  §5.3 — the transmission-emission pair ``(a, b)`` for DD.
* Blelloch, G. E. (1990). *Prefix Sums and Their Applications*. CMU-CS-90-190.
  §1.5 — first-order linear recurrence closed form.
* Hébert, A. (2009). *Applied Reactor Physics*. §3.9.4 — curvilinear DD with
  Morel-Montry angular closure.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from orpheus.sn.geometry import SNMesh
    from orpheus.sn.spatial.scheme import DiscretizationSchemeBase


# ═══════════════════════════════════════════════════════════════════════
# Stratum 1 — Geometry coefficients (built once at solver construction)
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class GeometryCoefficients:
    r"""Stratum 1: geometry + quadrature only.  Built ONCE per ``SNMesh``.

    Carries every per-ordinate-per-cell geometric quantity that enters the DD
    coefficients :math:`(a, b)` of the per-ordinate spatial scan.  Survives
    :math:`\Sigma_t` rebinds (depletion, thermal feedback), BC changes, every
    outer Picard / power iteration, every inner SI / Krylov iteration.

    No ``ng`` axis — Stratum 1 has no cross-section dependence.

    Shape contract (slab :math:`\eta\equiv` cyl-degenerate excluded — see
    ``is_degenerate``):

    +-----------------+----------------+---------------------------------------+
    | Field           | Shape          | Meaning                               |
    +=================+================+=======================================+
    | ``chain_idx``   | ``(N, nx)`` int| Cell index in chain order per ordinate|
    +-----------------+----------------+---------------------------------------+
    | ``chain_idx_inv``|``(N, nx)`` int| Inverse: cell index → chain position  |
    +-----------------+----------------+---------------------------------------+
    | ``abs_mu``      | ``(N,)``       | :math:`|\mu_n|`                       |
    +-----------------+----------------+---------------------------------------+
    | ``A_down``      | ``(N, nx)``    | Chain-ordered downstream face area    |
    |                 |                | (slab: ``1``; sphere/cyl: outer / inner |
    |                 |                | per ordinate sign; cyl-deg: ``0``)    |
    +-----------------+----------------+---------------------------------------+
    | ``A_total``     | ``(N, nx)``    | Chain-ordered ``A_inner + A_outer``   |
    |                 |                | (slab: ``2``; curvilinear: physical)  |
    +-----------------+----------------+---------------------------------------+
    | ``dA_w``        | ``(N, nx)``    | Chain-ordered :math:`\Delta A / w_n`  |
    |                 |                | (slab: ``0``; curvilinear: physical)  |
    +-----------------+----------------+---------------------------------------+
    | ``V``           | ``(N, nx)``    | Chain-ordered cell volume per ordinate|
    +-----------------+----------------+---------------------------------------+
    | ``c_in``        | ``(N,)``       | M-M closure constant                  |
    |                 |                | :math:`(1-\tau)/\tau\,\alpha_{\rm out}|
    |                 |                | + \alpha_{\rm in}` (slab: ``0``)      |
    +-----------------+----------------+---------------------------------------+
    | ``c_out``       | ``(N,)``       | :math:`\alpha_{\rm out}/\tau`         |
    |                 |                | (slab: ``0``)                         |
    +-----------------+----------------+---------------------------------------+
    | ``tau_inv``     | ``(N,)``       | :math:`1/\tau_n` (slab: ``1``)        |
    +-----------------+----------------+---------------------------------------+
    | ``mm_a_in_coeff``| ``(N,)``      | :math:`(1-\tau_n)/\tau_n` (slab: ``0``)|
    +-----------------+----------------+---------------------------------------+
    | ``is_degenerate``| ``(N,)`` bool | Cylindrical pure-azimuthal ordinate   |
    |                 |                | (rare; routes to slow per-cell path)  |
    +-----------------+----------------+---------------------------------------+
    | ``level_ordinates``| see notes   | Curvilinear only: per-:math:`\mu`-    |
    |                 |                | level global-ordinate lists.          |
    +-----------------+----------------+---------------------------------------+

    The ``level_ordinates`` field is ``None`` for slab geometry and a list of
    1-D ``int`` arrays for sphere / cylinder.  For sphere it is a single-level
    list ``[range(N)]``; for cylinder it is one array per :math:`\mu`-level
    (cylindrical ordinate grouping).

    Per ``vv-principles`` Smell #16 (cross-domain-attacker memo): the strata
    are deliberately separated; do NOT pack :math:`\Sigma_t`-dependent fields
    here.  Those belong on :class:`CollisionCache`.
    """

    chain_idx: np.ndarray              # (N, nx) int
    chain_idx_inv: np.ndarray          # (N, nx) int — inverse permutation
    abs_mu: np.ndarray                 # (N,)
    A_down: np.ndarray                 # (N, nx) — chain-ordered
    A_total: np.ndarray                # (N, nx) — chain-ordered
    dA_w: np.ndarray                   # (N, nx) — chain-ordered
    V: np.ndarray                      # (N, nx) — chain-ordered
    c_in: np.ndarray                   # (N,)
    c_out: np.ndarray                  # (N,)
    tau_inv: np.ndarray                # (N,)
    mm_a_in_coeff: np.ndarray          # (N,)
    mu_start: np.ndarray               # (N,) — level's starting-direction edge
    is_degenerate: np.ndarray          # (N,) bool
    level_ordinates: tuple[np.ndarray, ...] | None = None
    r"""Curvilinear only: per-:math:`\mu`-level ordinate index lists.

    For sphere: ``(np.arange(N),)`` (single virtual level — every ordinate).
    For cylinder: one ``np.ndarray`` per :math:`\mu`-level, each carrying the
    global ordinate indices of its within-level azimuthal ordinates.
    ``None`` for slab.
    """

    @classmethod
    def from_mesh_and_quad(cls, sn_mesh: "SNMesh") -> "GeometryCoefficients":
        r"""Populate Stratum 1 from one :class:`SNMesh` + its quadrature.

        Iterates ``sn_mesh.dag_walk(ordinate_idx=...)`` (slow Python path —
        but ONLY ONCE per solver lifetime; cost amortised across every
        subsequent sweep).  The per-cell
        :class:`~orpheus.geometry.reduced_operator.StreamingTerms`
        dataclass is unpacked into chain-ordered numpy arrays once and never
        rematerialised.
        """
        from orpheus.geometry import CoordSystem  # local import: cyclic risk

        quad = sn_mesh.quad
        N = quad.N
        nx = sn_mesh.nx
        reduced = sn_mesh.reduced
        assert reduced is not None, (
            "GeometryCoefficients requires a ReducedStreamingOperator "
            "(1-D Cartesian / spherical / cylindrical).  2-D Cartesian "
            "wavefront uses anti-diagonal scheduling, not the chain scan."
        )
        coord = reduced.coord

        # ── Per-ordinate scalars (slab carries neutral M-M constants) ─
        abs_mu = np.abs(np.asarray(quad.mu_x, dtype=np.float64))  # (N,)

        # ── Per-ordinate-per-cell chain-ordered tensors ───────────────
        chain_idx = np.empty((N, nx), dtype=np.int64)
        A_down = np.empty((N, nx), dtype=np.float64)
        A_total = np.empty((N, nx), dtype=np.float64)
        dA_w = np.empty((N, nx), dtype=np.float64)
        V = np.empty((N, nx), dtype=np.float64)
        tau = np.empty(N, dtype=np.float64)
        alpha_in = np.empty(N, dtype=np.float64)
        alpha_out = np.empty(N, dtype=np.float64)
        mu_start = np.empty(N, dtype=np.float64)
        is_degenerate = np.zeros(N, dtype=bool)

        # ── Level enumeration (sphere = single virtual level) ─────────
        level_ordinates: tuple[np.ndarray, ...] | None
        if coord is CoordSystem.CARTESIAN:
            level_ordinates = None
            level_visits_iter = [(None, n, n) for n in range(N)]
        elif coord is CoordSystem.SPHERICAL:
            level_ordinates = (np.arange(N, dtype=np.int64),)
            level_visits_iter = [(None, n, n) for n in range(N)]
        elif coord is CoordSystem.CYLINDRICAL:
            level_indices = quad.level_indices  # type: ignore[attr-defined]
            level_ordinates = tuple(
                np.asarray(lvl, dtype=np.int64) for lvl in level_indices
            )
            level_visits_iter = []
            for p_idx, lvl in enumerate(level_indices):
                for m_local, global_n in enumerate(lvl):
                    level_visits_iter.append((p_idx, m_local, int(global_n)))
        else:
            raise ValueError(  # pragma: no cover — exhaustive match
                f"Unknown coord system: {coord!r}"
            )

        for mu_level_idx, ordinate_idx, global_n in level_visits_iter:
            visits = list(sn_mesh.dag_walk(
                ordinate_idx=ordinate_idx,
                mu_level_idx=mu_level_idx,
            ))
            # Cell indices in chain (sweep-direction-resolved) order.
            chain = np.fromiter(
                (v.cell_idx for v in visits),
                dtype=np.int64, count=nx,
            )
            chain_idx[global_n] = chain
            # Unpack streaming-terms fields once per visit; chain-order
            # by construction (the visit iterator yields them in chain
            # order already).
            A_down[global_n] = np.fromiter(
                (v.face_area_downstream for v in visits),
                dtype=np.float64, count=nx,
            )
            A_total[global_n] = np.fromiter(
                (v.streaming_terms.face_area_inner
                 + v.streaming_terms.face_area_outer for v in visits),
                dtype=np.float64, count=nx,
            )
            dA_w[global_n] = np.fromiter(
                (v.streaming_terms.delta_A_over_w for v in visits),
                dtype=np.float64, count=nx,
            )
            V[global_n] = np.fromiter(
                (v.streaming_terms.volume for v in visits),
                dtype=np.float64, count=nx,
            )
            # The M-M closure scalars are ordinate-only; pick the first
            # visit's value (all visits within one ordinate carry the
            # same alpha/tau pair).
            st0 = visits[0].streaming_terms
            tau[global_n] = st0.tau_mm
            alpha_in[global_n] = st0.alpha_in
            alpha_out[global_n] = st0.alpha_out
            mu_start[global_n] = st0.mu_start
            # Cylindrical pure-azimuthal degenerate: visit carries
            # face_area_downstream == 0.0 (geometric truth).  The slow
            # per-cell path handles these ordinates.
            if visits[0].face_area_downstream == 0.0 and abs_mu[global_n] < 1e-15:
                is_degenerate[global_n] = True

        # ── M-M closure constants (slab: alpha=0, tau=1 → c=0) ────────
        c_out = alpha_out / tau
        c_in = (1.0 - tau) / tau * alpha_out + alpha_in
        tau_inv = 1.0 / tau
        mm_a_in_coeff = (1.0 - tau) / tau

        # ── Inverse chain index for scatter-back ──────────────────────
        chain_idx_inv = np.empty_like(chain_idx)
        cols = np.arange(nx, dtype=np.int64)[None, :]
        np.put_along_axis(chain_idx_inv, chain_idx, cols, axis=1)

        return cls(
            chain_idx=chain_idx,
            chain_idx_inv=chain_idx_inv,
            abs_mu=abs_mu,
            A_down=A_down,
            A_total=A_total,
            dA_w=dA_w,
            V=V,
            c_in=c_in,
            c_out=c_out,
            tau_inv=tau_inv,
            mm_a_in_coeff=mm_a_in_coeff,
            mu_start=mu_start,
            is_degenerate=is_degenerate,
            level_ordinates=level_ordinates,
        )


# ═══════════════════════════════════════════════════════════════════════
# Stratum 2 — Collision cache (built when σ_t binds; instrumented counter)
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class CollisionCache:
    r"""Stratum 2: geometry :math:`\times \Sigma_t`.  Built once per :math:`\Sigma_t`.

    Combines :class:`GeometryCoefficients` with the per-cell-per-group
    :math:`\Sigma_t` to produce the per-ordinate DD scan coefficients that
    survive every source iteration within one constant-:math:`\Sigma_t`
    epoch:

    * ``inverse_denom[n, g, i] = 1 / (2|\mu_n|\,A_{\rm down}[n,i] + (\Delta A/w)\,c_{\rm out}[n] + \Sigma_t[g,i]\,V[i])``
    * ``a_attenuation[n, g, i] = 2|\mu_n|\,A_{\rm total}[i] \cdot \mathrm{inverse\_denom}[n, g, i] - 1``
    * ``cumprod_a[n, g, i] = \prod_{k \le i} a[n, g, k]``  (along chain / cell axis)

    Index convention (Issue #196 PR-INDEX-2).  The principled storage layout
    is ``(N, ng, nx)`` — energy ``g`` is the *second* axis, NOT trailing.
    See ``.claude/plans/principled_index_migration.md`` §1 derivation table.
    Under this layout the chain / cell axis is axis 2; the ``cumprod_a``
    cumulative product runs along ``axis=2`` (cell-wise).

    Cache invariance contract — the load-bearing performance gate.  Within
    one constant-:math:`\Sigma_t` epoch (one fixed-source solve, one
    eigenvalue solve, ...), :meth:`from_geometry` is called EXACTLY ONCE.
    The :attr:`_build_count` class variable instruments this — pinning the
    invariant in
    ``tests/sn/spatial/test_sweep_cache.py::test_collision_cache_invariance_under_source_iteration``.

    Storage: ``(N, ng, nx)`` per field × 3 fields × 8 bytes.  Canonical
    ``(N=16, ng=2, nx=160)`` problem ≈ 240 kB.
    """

    inverse_denom: np.ndarray        # (N, ng, nx) — chain-ordered along nx (axis 2)
    a_attenuation: np.ndarray        # (N, ng, nx) — chain-ordered along nx (axis 2)
    cumprod_a: np.ndarray            # (N, ng, nx) — chain-ordered along nx (axis 2)
    cell_average_weight: np.ndarray  # (N, ng, nx) — the scheme's blend weight w
    r"""The per-cell cell-average blend weight ``w`` (#158 coefficient model):
    :math:`\bar\psi=(1-w)\psi_{\rm in}+w\,\psi_{\rm out}`.  DD is ``w=½``
    everywhere; LD is ``w=1/(1+k)``.  Stored chain-ordered alongside
    ``a_attenuation`` / ``inverse_denom`` so the scan body and the matvec apply
    the generic base reconstruction staticmethods
    (:meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.source_emission` /
    :meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.cell_average` /
    :meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.outgoing_face_from_average`)
    without re-deriving any per-scheme cell math."""

    _build_count: ClassVar[int] = 0
    """Class-level counter incremented on every :meth:`from_geometry`.

    Instrumentation for the cardinal cache-invariance test (#4 in plan
    §"Test catalog"): under a converged Picard loop, this counter must
    advance by exactly 1.
    """

    @classmethod
    def reset_build_count(cls) -> None:
        """Reset the class-level build counter to zero.

        Used by the cache-invariance test to establish a clean baseline
        before a measured run.
        """
        cls._build_count = 0

    @classmethod
    def from_geometry(
        cls,
        geom: GeometryCoefficients,
        sig_t: np.ndarray,
        scheme: "DiscretizationSchemeBase",
    ) -> "CollisionCache":
        r"""Populate Stratum 2 from Stratum 1 + per-cell :math:`\Sigma_t`.

        The :math:`\Sigma_t`-epoch DD scan coefficients
        ``(a_attenuation, inverse_denom)`` are owned by the cell-update
        scheme (Issue #236 §2): this method delegates their three numpy
        ops to :meth:`scheme.affine_scan_coefficients
        <orpheus.sn.spatial.scheme.DiscretizationScheme.affine_scan_coefficients>`
        so the cache reflects whichever spatial closure ``SNMesh`` selected
        — the cache keeps storage + lifetime; the scheme owns the math
        (Cardinal Rule 2 / Pattern 2, single source of truth).  This cache
        feeds the DAG-free scan schedules (``CumprodScan`` / ``ScanMarch``),
        so ``scheme`` is always an ``is_affine_scannable`` scheme here
        (the scan strategies' ``supports`` gate guarantees it).

        Pure broadcasting — three numpy ops in
        ``affine_scan_coefficients`` plus the order-dependent ``cumprod_a``
        here, no Python per-cell loop.

        Parameters
        ----------
        geom : GeometryCoefficients
            The Stratum 1 cache.
        sig_t : ndarray, shape ``(ng, nx)``.
            Per-group per-cell total cross section.  Shape matches the
            principled 1-D sweep contract (``ng``, ``sn_mesh.nx``) — see
            Issue #196 PR-INDEX-2 in
            ``.claude/plans/principled_index_migration.md``.
        scheme : DiscretizationSchemeBase
            The selected spatial closure scheme (e.g.
            :class:`~orpheus.sn.spatial.diamond.DiamondDifference`).  Must
            be ``is_affine_scannable``; supplies the closed-form recurrence
            coefficients via :meth:`affine_scan_coefficients`.

        Notes
        -----
        Output layout is ``(N, ng, nx)`` for every field.  The cumulative
        product ``cumprod_a`` runs along ``axis=2`` (the trailing cell
        axis); under the principled layout the cell axis is NOT axis 1.
        The ``cumprod_a`` stays HERE (not in the cell-update) because it is
        a *scan-schedule* transform — a prefix product along the chain
        order — not closure math.
        """
        cls._build_count += 1

        # ── σ_t chain-ordered per ordinate: (N, ng, nx) ───────────────
        # sig_t is (ng, nx); reorder the cell axis (axis 1 of sig_t)
        # by geom.chain_idx (N, nx).  Result has shape (ng, N, nx);
        # transpose to (N, ng, nx) to match the principled layout.  This
        # chain reorder is a scan-schedule data-prep step, so it stays in
        # the cache builder (the cell-update is ordering-agnostic).
        sig_t_chain = sig_t[:, geom.chain_idx].transpose(1, 0, 2)        # (N, ng, nx)

        # ── (a, 1/denom) delegated to the cell-update scheme ──────────
        # The σ_t-epoch, source-independent recurrence coefficients
        # a = 2|μ|·A_total/denom − 1, denom = 2|μ|·A_down + dA_w·c_out
        # + Σ_t·V.  The scheme owns the closure math; the cache keeps the
        # storage and the (order-dependent) cumprod.
        a_attenuation, inverse_denom, cell_average_weight = (
            scheme.affine_scan_coefficients(
                abs_mu=geom.abs_mu,
                A_down=geom.A_down,
                A_total=geom.A_total,
                dA_w=geom.dA_w,
                c_out=geom.c_out,
                V=geom.V,
                sig_t=sig_t_chain,
            )
        )                                                                # all (N, ng, nx)

        # ── cumprod along the cell axis (axis 2 in principled layout) ─
        cumprod_a = np.cumprod(a_attenuation, axis=2)                    # (N, ng, nx)

        return cls(
            inverse_denom=inverse_denom,
            a_attenuation=a_attenuation,
            cumprod_a=cumprod_a,
            cell_average_weight=cell_average_weight,
        )


__all__ = ["CollisionCache", "GeometryCoefficients"]
