r"""Two-stratum precomputed cache for 1-D SN sweeps (Issue #196 Phase G Step 2.5c).

Issue #196 Phase G Step 2.5c factors the per-ordinate-per-cell DD coefficients
into TWO frozen dataclasses by **mutation cadence**, eliminating per-sweep
:func:`np.fromiter` plus per-cell ``StreamingTerms`` dataclass construction
from the hot path:

* :class:`StreamingCoefficientCache` — Stratum 1, geometry+quadrature only.  Built
  once at :class:`~orpheus.sn.solver.SNSolver` construction.  No ``ng`` axis;
  invariant under cross-section rebinds, BC changes, every outer / inner /
  Picard iteration.  Lifetime = ``SNMesh`` × ``Quadrature``.

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

Closed-form scan (Blelloch §1.5; :func:`~orpheus.sn.sweep.scan.ordinate_scan`):

.. math::

   \psi^{s}[n, g, i]
       \;=\; \mathrm{cumprod\_a}[n, g, i]\,
             \bigl(\psi^{s}_0[n, g] + \mathrm{cumsum}_{k \le i}
                   (b[n, g, k] / \mathrm{cumprod\_a}[n, g, k])\bigr).

The :class:`StreamingCoefficientCache` populator hoists the geometry tensors
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

The cache populator computes its ``(a, denom)`` from its own vectorized
expression over the per-cell
:class:`~orpheus.transport.spatial.scheme.StreamingTerms` factory's data —
independently of
:func:`~orpheus.transport.spatial.cell_balance.cell_balance_for_streaming`,
the per-cell balance's single algebra source (P4.9a retired its scalar
twin ``cell_balance_terms``).  The L1 dual-view validator
(``test_cache_populator_matches_cell_balance_for_streaming``) asserts that
for any cell the two independent implementations agree to ``rtol=1e-14``.
This is the Pattern 2 contract that keeps the fast cache-driven path and
the slow per-cell-update reference path consistent.

References
==========

* Lewis, E. E., & Miller, W. F. (1984). *Computational Methods of Neutron
  Transport*.  §5.3 — the transmission-emission pair ``(a, b)`` for DD.
* Blelloch, G. E. (1990). *Prefix Sums and Their Applications*. CMU-CS-90-190.
  §1.5 — first-order linear recurrence closed form.
* Hébert, A. (2009). *Applied Reactor Physics*. §3.9.3 (cylinder) / §3.9.4
  (sphere) — the curvilinear DD cell balance and sweep mechanics.  The
  Morel-Montry weighted angular closure the DD update composes with is
  NOT his: see Morel & Montry (1984), TTSP 13(5):615-633, and
  Bailey-Morel-Chang (2010), NSE 165(2):149-169 Eqs. (42)/(43).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from ..angular.closure import AngularClosureBase
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.transport.spatial.scheme import DiscretizationSchemeBase


# ═══════════════════════════════════════════════════════════════════════
# Stratum 1 — Geometry coefficients (built once at solver construction)
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class StreamingCoefficientCache:
    r"""The Σ\ :sub:`t`\ -free half of the DD scan coefficients, built ONCE per
    ``SNMesh`` × ``Quadrature``.

    Carries every per-ordinate-per-cell quantity that enters the DD
    coefficients :math:`(a, b)` of the per-ordinate spatial scan.  Survives
    :math:`\Sigma_t` rebinds (depletion, thermal feedback), BC changes, every
    outer Picard / power iteration, every inner SI / Krylov iteration.

    ⛔ **The name is deliberately stratum-AGNOSTIC, because the contents are
    not one stratum.**  This class was ``GeometryCoefficients`` until
    2026-08-26; that name was measurably false — `[M]` **0 of 13** fields are
    un-permuted chart data.  Its fields span **three** invalidation classes,
    measured by building three meshes against one quadrature (uniform ``nx=6``,
    uniform ``nx=20``, and a GRADED ``nx=6``):

    * **S0 — mesh-free** (`[M]` bit-identical on all three): ``abs_mu``,
      ``c_in``, ``c_out``, ``tau_inv``, ``mm_a_in_coeff``, ``is_degenerate``,
      ``level_ordinates`` — **7 fields**, invalidated by a new quadrature or
      chart ONLY.
    * **S1 — chart × basis** (`[M]` differ on every re-mesh): ``A_down``,
      ``A_total``, ``dA_w``, ``V`` — **4 fields**.
    * **S3 — traversal**: ``chain_idx``, ``chain_idx_inv`` — **2 fields**,
      `[M]` identical between the uniform and GRADED ``nx=6`` meshes, so they
      turn on ordinate sign and cell COUNT, not on edge positions.

    ⚠ So the whole object rebuilds on any re-mesh, including the 7 fields that
    provably cannot change — the F3 weld.  A name asserting any one stratum
    would bless it: ``ChainScanCoefficients`` (the un-weld plan's first
    proposal) names the **S3** half, 2 of 13, which is *worse* than the
    ``Geometry`` it replaces.

    ⛔ **Do not propose ``SweepCoefficientCache`` — the name is TAKEN by a
    REFUTED design.**  It was the un-weld plan's second candidate and was
    rejected on discovery: ``.claude/lessons.md`` **L15** records
    ``SweepCoefficientCache (N, nx, ng)``, a single cache mixing geometry and
    :math:`\Sigma_t` fields, as *"the wrong shape"* — the very monolith whose
    split produced this class and :class:`CollisionCache`.  `[M]` no class by
    that name ever existed (`git log -S "class SweepCoefficientCache"` is
    empty); it was a proposal in a cross-domain-attacker memo.  Reusing it
    would make L15 read as condemning the shipped design, in a file loaded at
    every session start.

    ⟹ The name that survived says what is true of all three strata AND
    realizes the algebra L15 itself states: *"L (streaming + curvature) lives
    in [this class]; C joins via* :math:`1/(g_{\rm streaming} + \Sigma_t
    g_{\rm volume})` *to form* :class:`CollisionCache`\ *"*.  This is **L**'s
    coefficient cache; its sibling is **L + C**'s.  The stratum split is
    scheduled as phase **P4b**, not forgotten — its measurement and its
    done-when are at ``.claude/plans/streaming_path_says_what_it_is.md``
    **§4bis** (§4 states the strata; §4bis measures THIS class against them).

    ⚠ **And L15's own lesson applies one level down, to this class.** L15 is
    *"cache shape that mixes immutability strata hurts twice"*, and it holds
    THIS class up as the right shape — while `[M]` it mixes three.  The 2026
    measurement above is that lesson re-applied to its own worked example.

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
    is_degenerate: np.ndarray          # (N,) bool
    level_ordinates: tuple[np.ndarray, ...] | None = None
    r"""Curvilinear only: per-:math:`\mu`-level ordinate index lists.

    For sphere: ``(np.arange(N),)`` (single virtual level — every ordinate).
    For cylinder: one ``np.ndarray`` per :math:`\mu`-level, each carrying the
    global ordinate indices of its within-level azimuthal ordinates.
    ``None`` for slab.
    """

    @classmethod
    def from_mesh_and_quad(
        cls,
        sn_mesh: "SNMesh",
        angular_closure: "AngularClosureBase",
    ) -> "StreamingCoefficientCache":
        r"""Populate Stratum 1 from one :class:`SNMesh` + its quadrature.

        Iterates ``sn_mesh.dag_walk(ordinate_idx=...)`` (slow Python path —
        but ONLY ONCE per solver lifetime; cost amortised across every
        subsequent sweep).  The per-cell
        :class:`~orpheus.transport.spatial.scheme.StreamingTerms`
        dataclass is unpacked into chain-ordered numpy arrays once and never
        rematerialised.
        """
        from orpheus.geometry import CoordSystem  # local import: cyclic risk

        quad = sn_mesh.quad
        N = quad.N
        nx = sn_mesh.nx
        if not sn_mesh.is_1d:
            # A domain/admission contract, not type-narrowing — a bare
            # ``assert`` here is a NO-OP under the canonical ``python -O``
            # runner (coding-standards, the bare-assert clause; converted
            # at P4.9b step 2c with its first witness).  Keyed on the
            # honest predicate since P4.5: the chain scan is a 1-D
            # construct, and ``reduced`` presence is its ctor-guaranteed
            # realization (populated iff ``is_1d``).
            raise TypeError(
                "StreamingCoefficientCache requires a ReducedStreamingOperator "
                "(1-D Cartesian / spherical / cylindrical).  2-D Cartesian "
                "wavefront uses anti-diagonal scheduling, not the chain scan."
            )
        coord = sn_mesh.coord

        # ── Per-ordinate scalars (slab carries neutral M-M constants) ─
        abs_mu = np.abs(np.asarray(quad.mu_x, dtype=np.float64))  # (N,)

        # ── Per-ordinate-per-cell chain-ordered tensors ───────────────
        chain_idx = np.empty((N, nx), dtype=np.int64)
        A_down = np.empty((N, nx), dtype=np.float64)
        A_total = np.empty((N, nx), dtype=np.float64)
        dA_w = np.empty((N, nx), dtype=np.float64)
        V = np.empty((N, nx), dtype=np.float64)
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
            # Cylindrical pure-azimuthal degenerate: visit carries
            # face_area_downstream == 0.0 (geometric truth).  The slow
            # per-cell path handles these ordinates.
            if visits[0].face_area_downstream == 0.0 and abs_mu[global_n] < 1e-15:
                is_degenerate[global_n] = True

        # ── M-M angular weight τ + closure constants (Issue #236) ──
        # τ, ``c_out = α_{m+1/2}/τ``, ``c_in = (1−τ)/τ·α_{m+1/2} +
        # α_{m−1/2}``, AND the scan-normal march constants ``1/τ`` /
        # ``(1−τ)/τ`` are ANGULAR-closure properties — the pole-angular
        # closure is their canonical owner.  Read the ``(N,)``
        # global-ordinate views directly instead of rebuilding any formula
        # here (Cardinal Rule 2; Pattern 7 — normalise at the definition
        # site).  The closure dispatches by TYPE: MorelMontry
        # (sphere/cylinder) returns its precomputed values; the Cartesian
        # IdentityAngularClosure returns the neutral ones (α=0, τ=1 ⇒ c=0,
        # tau_inv=1, a_in_coeff=0) — so this populator stays geometry-blind
        # by data.
        #
        # ⛔ REVISED at P4.9a (2026-08-28, user-ruled), overturning the
        # 2026-08-26 ruling that stood here — "the closure exposes only the
        # PRIMITIVE τ (Pattern 5); the cache derives 1/τ and (1−τ)/τ as a
        # legit L16 perf hoist".  The revision: ``(1−τ)/τ`` is not a
        # convenience product, it is the ψ_in coefficient of the closure's
        # OWN update in scan-normal form — the PAIRING is relation
        # knowledge, and the cache deriving it was the cache spelling a
        # fragment of the Morel-Montry relation (the same smell, one notch
        # down, as the DiamondDifference twin P4.9a removed).  The closure
        # now MINTS both scan constants (tau_inv_per_ordinate /
        # march_a_in_coeff_per_ordinate, closure.py); L16's hoisting half
        # SURVIVES here — this cache is still where the values are stored
        # once per solve for the per-iteration scan.  Consumed at the
        # ``loss_representation`` scan fast path
        # ``geom.tau_inv·ψ_avg − geom.mm_a_in_coeff·ψ_in`` — the closure's
        # own scan-normal representation, welded to its march by gate
        # (test_angular_closure::TestMintedScanConstants; the M7
        # mutation arm prices the realistic 1-2 ULP respelling).
        # Bit-identical to the pre-handing derivation: the accessors
        # compute the same expressions on the same cached τ array
        # (test_cache's closure-algebra field gate pins it array_equal).
        closure = angular_closure
        c_out = closure.c_out_per_ordinate
        c_in = closure.c_in_per_ordinate
        tau_inv = closure.tau_inv_per_ordinate
        mm_a_in_coeff = closure.march_a_in_coeff_per_ordinate

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
            is_degenerate=is_degenerate,
            level_ordinates=level_ordinates,
        )


# ═══════════════════════════════════════════════════════════════════════
# Stratum 2 — Collision cache (built when σ_t binds; instrumented counter)
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class CollisionCache:
    r"""Stratum 2: geometry :math:`\times \Sigma_t`.  Built once per :math:`\Sigma_t`.

    Combines :class:`StreamingCoefficientCache` with the per-cell-per-group
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
    ``tests/sn/sweep/core/test_cache.py::test_collision_cache_invariance_under_source_iteration``.

    Storage: ``(N, ng, nx)`` per field × 3 fields × 8 bytes.  Canonical
    ``(N=16, ng=2, nx=160)`` problem ≈ 240 kB.
    """

    inverse_denom: np.ndarray        # (N, ng, nx) — chain-ordered along nx (axis 2)
    a_attenuation: np.ndarray        # (N, ng, nx) — chain-ordered along nx (axis 2)
    cumprod_a: np.ndarray            # (N, ng, nx) — chain-ordered along nx (axis 2)
    face_blend_weight: np.ndarray  # (N, ng, nx) — the scheme's blend weight w
    r"""The per-cell cell-average blend weight ``w`` (#158 coefficient model):
    :math:`\bar\psi=(1-w)\psi_{\rm in}+w\,\psi_{\rm out}`.  DD is ``w=½``
    everywhere; LD is ``w=1/(1+k)``.  Stored chain-ordered alongside
    ``a_attenuation`` / ``inverse_denom`` so the scan body and the matvec apply
    the generic base reconstruction staticmethods
    (:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.source_emission` /
    :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.cell_average` /
    :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.outgoing_face_from_average`)
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
        geom: StreamingCoefficientCache,
        sig_t: np.ndarray,
        scheme: "DiscretizationSchemeBase",
    ) -> "CollisionCache":
        r"""Populate Stratum 2 from Stratum 1 + per-cell :math:`\Sigma_t`.

        The :math:`\Sigma_t`-epoch DD scan coefficients
        ``(a_attenuation, inverse_denom)`` are owned by the cell-update
        scheme (Issue #236 §2): this method delegates their three numpy
        ops to :meth:`scheme.affine_scan_coefficients
        <orpheus.transport.spatial.scheme.DiscretizationSchemeBase.affine_scan_coefficients>`
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
        geom : StreamingCoefficientCache
            The Stratum 1 cache.
        sig_t : ndarray, shape ``(ng, nx)``.
            Per-group per-cell total cross section.  Shape matches the
            principled 1-D sweep contract (``ng``, ``sn_mesh.nx``) — see
            Issue #196 PR-INDEX-2 in
            ``.claude/plans/principled_index_migration.md``.
        scheme : DiscretizationSchemeBase
            The selected spatial closure scheme (e.g.
            :class:`~orpheus.transport.spatial.diamond.DiamondDifference`).  Must
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
        # P4.9a row 3b: the caller ASSEMBLES the closure's denominator
        # contribution ((ΔA/w)·c_out — same expression DD used to build
        # in-scheme, bit-identical) so the scheme family sees no
        # closure-named constant.
        a_attenuation, inverse_denom, face_blend_weight = (
            scheme.affine_scan_coefficients(
                abs_mu=geom.abs_mu,
                A_down=geom.A_down,
                A_total=geom.A_total,
                angular_denom_term=geom.dA_w * geom.c_out[:, None],
                V=geom.V,
                reaction_xs=sig_t_chain,
            )
        )                                                                # all (N, ng, nx)

        # ── cumprod along the cell axis (axis 2 in principled layout) ─
        cumprod_a = np.cumprod(a_attenuation, axis=2)                    # (N, ng, nx)

        return cls(
            inverse_denom=inverse_denom,
            a_attenuation=a_attenuation,
            cumprod_a=cumprod_a,
            face_blend_weight=face_blend_weight,
        )


__all__ = ["CollisionCache", "StreamingCoefficientCache"]
