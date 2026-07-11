r"""The SN curvilinear starting-direction (ψ½) block Hilbert metric ``G_sd``.

Foundation gates pinning the state metric of the starting-direction block of
the augmented loss composite :math:`A = L + C` acting on
:math:`V = V_{\rm bulk}(\text{AngularFlux}) \oplus V_{\rm trace}
(\text{AngularBoundaryFlux}) \oplus V_{\rm sd}(\text{RadialCharacteristicFlux})`
(#282 route (a) / #280 phase 2.5d).

The derivation of record
(``.claude/agent-memory/numerics-investigator/radial_characteristic_metric_gauge_derivation.md``)
established, from the adjoint constraint :math:`A^{\mathsf T} G = G A^\dagger`:

* ``apply_transpose`` is the **EXACT Euclidean transpose** :math:`A^{\mathsf T}`
  (``‖T − Aᵀ‖ ≈ 3e-16``, every block including :math:`T_{sb} = A_{bs}^{\mathsf T}`),
  so ``A.H = G⁺·apply_transpose·G`` is the honest metric adjoint
  :math:`G^{-1} A^{\mathsf T} G` for ANY invertible ``G``;
* ⟹ reciprocity :math:`\langle A\psi,\varphi\rangle_G =
  \langle\psi, A^\dagger\varphi\rangle_G` is **GAUGE-FREE** — it holds for every
  SPD ``G_sd``; the determining equation pins ``G_sd`` only up to non-degeneracy;
* ``G_sd = 0`` (the retired "ghost metric") is the **single forbidden value** —
  the boundary of the SPD cone. It puts the seed in ``ker G``, so ``A.H`` severs
  the seed→bulk coupling ``A_bs ≠ 0`` (the Morel–Montry recurrence) and is a
  WRONG adjoint the instant ψ½ carries data — green ONLY because every reciprocity
  gate historically fed a present-but-ZERO seed (**vv-principles Mode 12**,
  ERR-067).

The shipped fix (``RadialCharacteristicSpace.for_levels`` builds ``G_sd = V_cell``,
the radial cell volume — mirroring ``G_bulk = V_cell·w_n``; the angular ``w`` is a
free gauge, a single :math:`\mu=\pm 1` ray carries no canonical quadrature weight)
moves the seed rows OUT of the reciprocity functional's invariance group, closing
Mode 12. ψ½ is a first-class radial STATE field (its self-block ``A_ss`` is a
banded radial transport operator :math:`\mu\partial_r + \sigma_t`, NOT a
face-restriction), so its Hilbert metric is the bulk's spatial measure restricted
to the ray.

These are algebraic / adjoint-structure invariants (exact by construction), so the
suite is ``@pytest.mark.foundation`` — no ``verifies()`` (no theory equation; the
claim is a software/operator-algebra invariant). Promoted from
``derivations/diagnostics/diag_gsd_0{1..5}_*.py`` (2026-07-06) after the
``G_sd = V_cell`` fix landed; the pre-fix RED gates were inverted to the post-fix
truth.

``-O``-safe (vv Mode 8): every gate is a ``pytest.fail`` / ``np.testing`` function
call — a bare ``assert`` is stripped to a NO-OP under the canonical ``python -O``.

References
==========

* GH #282 (the spherical seed lag), #280 (the walk unification), #229 (the
  cylinder τ clamp fact); ERR-067 (the ghost-metric Mode-12 blindness).
* Hébert, A. (2009). *Applied Reactor Physics*. §3.9.4 (the starting-direction
  equation, Eqs. 3.432–3.435).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.full_field import FullField
from orpheus.derivations.common.xs_library import make_mixture

pytestmark = pytest.mark.foundation


# ═══════════════════════════════════════════════════════════════════════
# Mode-8-safe scalar bound helpers (pytest.fail fires under ``python -O``;
# a bare ``assert`` in a helper would be stripped to a NO-OP)
# ═══════════════════════════════════════════════════════════════════════


def _fail_unless_below(value, bound, msg):
    """Fail (under ``-O``) unless ``float(value) < bound``."""
    v = float(value)
    if not (v < bound):
        pytest.fail(f"{msg}: {v:.3e} is not < {bound:.1e}")


def _fail_unless_above(value, bound, msg):
    """Fail (under ``-O``) unless ``float(value) > bound``."""
    v = float(value)
    if not (v > bound):
        pytest.fail(f"{msg}: {v:.3e} is not > {bound:.1e}")


# ═══════════════════════════════════════════════════════════════════════
# Builders + dense-probe helpers (single source — the five diagnostics each
# carried a near-copy; consolidation collapses them to one set)
# ═══════════════════════════════════════════════════════════════════════


def _mixture(sig_t: float, sig_s: float, ng: int):
    """Non-degenerate per-group mixture (only used to construct the SNMesh)."""
    st = np.array([sig_t * (1.0 + 0.4 * g) for g in range(ng)])
    ss = np.diag([sig_s * (1.0 + 0.4 * g) for g in range(ng)])
    return make_mixture(
        sig_t=st, sig_c=st - ss.sum(axis=0), sig_f=np.zeros(ng),
        nu=np.zeros(ng), chi=np.zeros(ng), sig_s=ss,
    )


def _build_sphere(nx: int, ng: int, sigma: float):
    """A seed-carrying sphere SNMesh (GL S4) + the loss composite ``A = L + C``.

    The collision σ_t (``sigma·(1+0.3g)``) differs per group from the mixture's
    (``sigma·(1+0.4g)``) so ``A`` is non-degenerate across the group axis — the
    adjoint-structure claims here are geometry/algebra facts, insensitive to the
    specific σ values.
    """
    mesh = Mesh1D(
        edges=np.linspace(0.0, 4.0, nx + 1), mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL, bc_right=BC("vacuum"),
    )
    sn = SNMesh(mesh, Quadrature.gauss_legendre(4), {0: _mixture(sigma, 0.4 * sigma, ng)})
    sig_t = np.stack(
        [np.full(sn.spatial_shape, sigma * (1.0 + 0.3 * g)) for g in range(ng)], axis=0)
    return sn, StreamingOperator(sn) + MultiplicationOperator.from_mesh(sig_t, sn)


def _joint_M(sn, LC):
    """The joint ``M`` — ``LC``'s fused walk on the coupled pair (B.2d: the
    ψ½ legs ride the explicit leaf kwargs; the pair is the honest carrier)."""
    from orpheus.numerics.coupled_system import CoupledSpace
    from orpheus.sn.coupled_system import CoupledInvertibleOperator

    space = CoupledSpace.from_systems(
        (sn.full_field_space, sn.radial_characteristic_composite_space),
    )
    return CoupledInvertibleOperator(LC, space=space, sn_mesh=sn), space


def _composite(sn, *, bulk: bool, trace: bool, seed: bool, rng):
    """The coupled pair with each block either random (True) or zero (False).

    Draws in block order bulk → trace → seed, only for the active blocks
    (the SAME order as the pre-eviction 3-block builder).
    """
    from orpheus.numerics.coupled_system import CoupledField
    from orpheus.transport.radial_characteristic_composite import (
        RadialCharacteristicComposite,
    )

    N, nx, ng = sn.quad.N, sn.nx, sn.ng
    n_tr = int(sn.angular_trace.layout.total_size)
    psi_a = FullField(
        interior=AngularFlux.from_mesh(
            rng.standard_normal((N, ng, nx)) if bulk else np.zeros((N, ng, nx)), sn),
        boundary=AngularBoundaryFlux(
            values=rng.standard_normal(n_tr) if trace else np.zeros(n_tr),
            space=sn.angular_trace, mesh=sn),
    )
    # System B is the native split composite (4e); a random seed fills its
    # to_flat (interior ⊕ boundary) — the reciprocity identity is layout-agnostic
    # (psi and phi share this convention).
    system_b = RadialCharacteristicComposite.from_mesh(sn)
    if seed:
        n_sd = sn.radial_characteristic_composite_space.shape[0]
        system_b = RadialCharacteristicComposite.from_flat(
            rng.standard_normal(n_sd), system_b,
        )
    return CoupledField(systems=(psi_a, system_b))


def _template(sn):
    """An all-zero coupled pair (a shape template)."""
    return _composite(sn, bulk=False, trace=False, seed=False, rng=np.random.default_rng(0))


def _dense(fn, tpl):
    """Assemble a dense matrix by unit-vector probing over ``to_flat`` coords."""
    n = tpl.to_flat().size
    M = np.zeros((n, n))
    for j in range(n):
        e = np.zeros(n)
        e[j] = 1.0
        M[:, j] = fn(type(tpl).from_flat(e, tpl)).to_flat()
    return M


def _blocks(sn):
    """``(bulk, trace, System-B)`` index slices in the COUPLED ``to_flat``
    order + their sizes (System B is the trailing member; its internal
    member order differs from the old unified interleave — every consumer
    here treats it as a BLOCK, so the analysis is permutation-invariant
    within the slice)."""
    N, nx, ng = sn.quad.N, sn.nx, sn.ng
    nb = N * ng * nx
    nt = int(sn.angular_trace.layout.total_size)
    ns = sn.radial_characteristic_composite_space.shape[0]
    return slice(0, nb), slice(nb, nb + nt), slice(nb + nt, nb + nt + ns), (nb, nt, ns)


def _seed_to_coupled_layout(sn, seed_values):
    """Identity passthrough (4e).

    Pre-4e this re-labelled a UNIFIED-layout seed diagonal onto the coupled
    System-B member layout via the ``from_unified`` gather. Since 4e the
    composite IS the native representation, so ``_v_cell_seed`` / ``_seed_vw``
    build directly in the composite ``to_flat`` order — the re-label is now
    the identity (its referent, the unified layout, is retired)."""
    return np.asarray(seed_values, dtype=float)


def _prod_metric(sn, tpl, space):
    """The SHIPPED production metric ``G`` as a flat diagonal (the coupled
    member-wise metric — System A's ``full_field_space`` ⊕ System B's
    composite instance)."""
    ones = type(tpl).from_flat(np.ones(tpl.to_flat().size), tpl)
    return space.apply_metric(ones).to_flat()


def _v_cell_seed(sn):
    """The ``G_sd = V_cell`` seed diagonal in the composite ``to_flat``
    (interior ⊕ boundary) order — per (level, sign): cells (ng, nx) = radial
    cell volume group-broadcast, corner (ng,) = outer-cell volume (gauge).

    Built by HAND from ``sn.volumes`` through the SPLIT spaces' own
    ``slot_view`` — the SAME layout ``for_levels`` populates — so it reproduces
    the shipped seed block bit-for-bit (a structurally-independent reference:
    the diagonal comes from ``sn.volumes``, not from the space's stored metric).
    """
    interior_space = sn.radial_characteristic_interior_space
    boundary_space = sn.radial_characteristic_boundary_space
    V = np.asarray(sn.volumes, dtype=float).ravel()
    interior = np.zeros(interior_space.shape[0])
    boundary = np.zeros(boundary_space.shape[0])
    for p in sn.radial_characteristic_levels:
        for sign in (-1, +1):
            interior_space.slot_view(interior, p, sign)[:] = V[None, :]
            boundary_space.slot_view(boundary, p, sign)[:] = V[-1]
    return np.concatenate([interior, boundary])


def _seed_vw(sn):
    """A ``G_sd = V_cell · w_start`` seed diagonal (composite ``to_flat`` order) —
    the "fold an angular weight in" gauge variant. Also SPD, also valid (the
    ``w`` is a free gauge)."""
    interior_space = sn.radial_characteristic_interior_space
    boundary_space = sn.radial_characteristic_boundary_space
    V = np.asarray(sn.volumes, dtype=float).ravel()
    w0 = float(np.asarray(sn.quad.weights, dtype=float)[0])
    interior = np.zeros(interior_space.shape[0])
    boundary = np.zeros(boundary_space.shape[0])
    for p in sn.radial_characteristic_levels:
        for sign in (-1, +1):
            interior_space.slot_view(interior, p, sign)[:] = V[None, :] * w0
            boundary_space.slot_view(boundary, p, sign)[:] = V[-1] * w0
    return np.concatenate([interior, boundary])


def _recip_defect(Afwd, T, g, s, *, seed_scale, flip_seed_self, rng):
    r"""Reciprocity defect :math:`|\langle A\psi,\varphi\rangle_G -
    \langle\psi, A.H\varphi\rangle_G| / \text{norm}` using dense operators and the
    diagonal metric ``g`` (faithful to production ``A.H = G⁺TG``).

    ``seed_scale`` : 0.0 → zero seed data in ψ,φ; 1.0 → random seed data.
    ``flip_seed_self`` : negate the forward ``A_ss`` block (the seed-row transpose
    inconsistency) WITHOUT flipping the transpose ``T`` (mirrors the production
    monkeypatch that a Mode-12-blind gate cannot see).
    """
    n = Afwd.shape[0]
    A = Afwd.copy()
    if flip_seed_self:
        A[s, s] *= -1.0
    gpos = g > 0
    ginv = np.where(gpos, 1.0 / np.where(gpos, g, 1.0), 0.0)
    G, Gpinv = np.diag(g), np.diag(ginv)
    H = Gpinv @ T @ G  # production adjoint A.H = G⁺TG

    def vec():
        v = rng.standard_normal(n)
        v[s] *= seed_scale
        return v

    psi, phi = vec(), vec()
    lhs = (A @ psi) @ (G @ phi)
    rhs = psi @ (G @ (H @ phi))
    return abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)


# ═══════════════════════════════════════════════════════════════════════
# (R1) The foundation the gauge-freedom rests on: T == Aᵀ EXACTLY
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("nx,ng,sigma", [
    (3, 2, 0.5), (5, 2, 1.0), (8, 2, 0.3), (4, 3, 0.7), (6, 3, 2.0),
])
def test_r1_transpose_is_euclidean_across_configs(nx, ng, sigma):
    r"""``apply_transpose`` is the EXACT Euclidean transpose ``Aᵀ`` — checked
    against numpy's ``Afwd.T`` (a structurally-INDEPENDENT ground, NOT the
    operator's own metric machinery) — at every (nx, ng, σ), not just the minimal
    case. This is why ``A.H = G⁺·apply_transpose·G`` is the honest metric adjoint
    for ANY invertible ``G``, i.e. reciprocity is gauge-free in ``G_sd``.

    Mutation tooth: a sign/index error in ``apply_transpose`` (e.g. dropping the
    ``T_sb = A_bsᵀ`` seed↔bulk coupling) makes ``‖T − Aᵀ‖rel`` jump O(1) → RED.
    """
    sn, A = _build_sphere(nx, ng, sigma)
    A, cspace = _joint_M(sn, A)
    tpl = _template(sn)
    Afwd = _dense(A.apply, tpl)
    T = _dense(A.apply_transpose, tpl)
    rel = np.max(np.abs(T - Afwd.T)) / (np.max(np.abs(Afwd)) + 1e-300)
    print(f"  nx={nx} ng={ng} σ={sigma}: ‖T−Aᵀ‖rel = {rel:.3e}, dim={Afwd.shape[0]}")
    _fail_unless_below(rel, 1e-12, f"apply_transpose ≠ exact Euclidean transpose (nx={nx},ng={ng})")


# ═══════════════════════════════════════════════════════════════════════
# The canonical Mode-12-closure gate (ERR-067)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.catches("ERR-067")
def test_derive_gsd_and_close_mode12():
    r"""THE canonical Mode-12-closure gate for ``G_sd`` (ERR-067).

    A dense, structurally-independent ``G⁺TG`` sweep over
    ``{0, identity, V_cell, V·w_start} × {zero-seed, random-seed, seed-row-flip}``,
    faithfulness-checked against the production ``A.H``. It proves:

    * **(a) zero-seed** — reciprocity holds for ALL candidates (the OLD gate's
      blind regime: with a zero seed the coupling term ``ψ½ᵀ A_bsᵀ G_b φ_b`` is
      identically zero regardless of ``G_sd``);
    * **(b) CONTROL LEG — unmutated random seed** — ``G_sd = 0`` BREAKS
      (unmatched ``A_bs`` coupling, defect ~1e-1); every SPD candidate HOLDS
      ``< 1e-12``. This is the leg that makes a (c) RED *attributable to the flip*:
      it pins the honest baseline, so a broken ``G_sd=0`` baseline (also ``> tol``
      on a random seed) cannot masquerade as "caught";
    * **(c) Mode-12 CLOSURE — seed-row (A_ss) sign flip** — REDS every SPD
      candidate (defect ~1e-2 to 1); ``G_sd = 0`` stays BLIND (the flip does not
      move the already-broken defect — flip ≡ no-flip).

    Mutation teeth: (i) reverting ``G_sd → 0`` collapses leg (b)/(c) to the blind
    ghost row (the closure is lost); (ii) the seed-row flip in leg (c) is the
    live tooth that REDS every SPD candidate. Structural independence: the
    reference adjoint is ``G⁺·apply_transpose·G`` with ``T`` cross-checked against
    numpy ``Afwd.T`` in ``test_r1_...`` (an independent ground).
    """
    sn, A = _build_sphere(4, 2, 0.5)
    A, cspace = _joint_M(sn, A)
    tpl = _template(sn)
    b, t, s, (nb, nt, ns) = _blocks(sn)
    Afwd = _dense(A.apply, tpl)
    T = _dense(A.apply_transpose, tpl)
    Hprod = _dense(A.H.apply, tpl)

    g_cur = _prod_metric(sn, tpl, cspace)

    def make_g(seed_diag):
        g = g_cur.copy()
        g[s] = seed_diag
        return g

    candidates = {
        "g_sd = 0 (ghost)":    make_g(np.zeros(ns)),
        "g_sd = 1 (identity)": make_g(np.ones(ns)),
        "g_sd = V_cell":       make_g(_seed_to_coupled_layout(sn, _v_cell_seed(sn))),
        "g_sd = V·w_start":    make_g(_seed_to_coupled_layout(sn, _seed_vw(sn))),
    }
    _spd = ("g_sd = 1 (identity)", "g_sd = V_cell", "g_sd = V·w_start")

    # Faithfulness: dense G⁺TG reproduces production A.H ⟹ the dense analysis IS
    # the production adjoint (not a parallel re-derivation).
    gpos = g_cur > 0
    ginv = np.where(gpos, 1.0 / np.where(gpos, g_cur, 1.0), 0.0)
    H_recon = np.diag(ginv) @ T @ np.diag(g_cur)
    faith = np.max(np.abs(H_recon - Hprod))
    print(f"\n=== Faithfulness ‖G⁺TG (dense) − A.H (production)‖ = {faith:.3e} ===")
    _fail_unless_below(faith, 1e-10, "dense G⁺TG does not reproduce production A.H")

    print(f"{'candidate':22s} {'zero-seed':>12s} {'random-seed':>12s} {'flip+random':>12s}")
    results = {}
    for name, g in candidates.items():
        rng = np.random.default_rng(2026)
        d_zero = _recip_defect(Afwd, T, g, s, seed_scale=0.0, flip_seed_self=False, rng=rng)
        rng = np.random.default_rng(2026)
        d_rand = _recip_defect(Afwd, T, g, s, seed_scale=1.0, flip_seed_self=False, rng=rng)
        rng = np.random.default_rng(2026)
        d_flip = _recip_defect(Afwd, T, g, s, seed_scale=1.0, flip_seed_self=True, rng=rng)
        results[name] = (d_zero, d_rand, d_flip)
        print(f"{name:22s} {d_zero:12.3e} {d_rand:12.3e} {d_flip:12.3e}")

    z0, r0, f0 = results["g_sd = 0 (ghost)"]

    # (a) zero-seed: reciprocity holds for ALL candidates (the blind regime).
    for name, (z, _, _) in results.items():
        _fail_unless_below(z, 1e-12, f"{name}: zero-seed reciprocity should hold")

    # (b) CONTROL LEG — unmutated random seed: ghost BREAKS; every SPD holds < 1e-12.
    _fail_unless_above(r0, 1e-3,
                       "ghost g_sd=0 random-seed reciprocity should BREAK (unmatched A_bs coupling)")
    for name in _spd:
        _fail_unless_below(results[name][1], 1e-12,
                           f"{name}: CONTROL LEG — unmutated random-seed reciprocity must hold")

    # (c) Mode-12 CLOSURE — seed-row flip REDS every SPD; the ghost stays BLIND.
    for name in _spd:
        _fail_unless_above(results[name][2], 1e-3,
                           f"{name}: seed-row flip should RED (Mode-12 closes)")
    _fail_unless_below(abs(f0 - r0) / (r0 + 1e-300), 1e-6,
                       "ghost g_sd=0 flip changed the defect — expected BLIND (no change)")

    print("=== VERDICT: T = Aᵀ exactly ⟹ reciprocity is GAUGE-FREE in g_sd; any SPD "
          "closes Mode-12 (flip reds) + is the honest adjoint for all seed data; "
          "g_sd=0 is blind AND a wrong adjoint. ===")


# ═══════════════════════════════════════════════════════════════════════
# Production-path corroborators (INVERTED from the pre-fix RED diagnostics —
# the shipped metric is now V_cell, so the seed-carrying regimes HOLD)
# ═══════════════════════════════════════════════════════════════════════


def test_production_path_vcell_honest_adjoint_for_nonzero_seed():
    r"""Production-path corroboration: with the SHIPPED ``G_sd = V_cell``, the REAL
    ``A.H`` (via the production ``full_field_space`` metric) is the honest adjoint
    even when the seed carries data — reciprocity HOLDS for a random-seed composite,
    not merely a zero-seed one.

    INVERTED from the pre-fix diagnostic (which asserted the ghost ``G_sd=0`` made
    the random-seed defect BREAK ``> 1e-3``). Post-fix the random-seed defect is
    ``~8e-16``.

    Mutation tooth: reverting production ``for_levels`` to the ``np.zeros`` ghost
    makes the random-seed defect jump to ``~1e-1`` (unmatched ``A_bs`` coupling) →
    RED. The zero-seed leg is the control (holds under both metrics).
    """
    sn, A = _build_sphere(4, 2, 0.5)
    A, cspace = _joint_M(sn, A)
    space = cspace
    rng = np.random.default_rng(7)

    def defect(seed):
        psi = _composite(sn, bulk=True, trace=True, seed=seed, rng=rng)
        phi = _composite(sn, bulk=True, trace=True, seed=seed, rng=rng)
        lhs = space.inner_product(A.apply(psi), phi)
        rhs = space.inner_product(psi, A.H.apply(phi))
        return abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)

    d_zero = defect(False)
    d_rand = defect(True)
    print(f"\n=== Production A.H reciprocity (shipped G_sd=V_cell) ===")
    print(f"  zero-seed   defect = {d_zero:.3e}  (control — holds under any metric)")
    print(f"  random-seed defect = {d_rand:.3e}  (HOLDS — V_cell is the honest adjoint)")
    _fail_unless_below(d_zero, 1e-12, "zero-seed reciprocity should hold")
    _fail_unless_below(d_rand, 1e-12,
                       "random-seed reciprocity must HOLD — shipped V_cell is the honest adjoint")


def test_shipped_metric_block_values():
    r"""The SHIPPED ``full_field_space`` metric: bulk ``V·w`` and trace ``|Ω·n|·w``
    are nonzero, and the seed block is now ``V_cell`` (SPD, strictly positive) — the
    metric is FULLY non-singular.

    INVERTED from the pre-fix diagnostic (which asserted the seed block was the
    identically-zero ghost and ``(g>0).sum() == nb+nt``). Post-fix the seed block
    is ``V_cell ∈ [~4, ~155]`` and ``(g>0).sum() == nb+nt+ns``.

    Mutation tooth: reverting production ``for_levels`` to ``np.zeros`` makes the
    seed block all-zero → the positivity check and the count ``nb+nt+ns`` both RED.
    """
    sn, A = _build_sphere(4, 2, 0.5)
    A, cspace = _joint_M(sn, A)
    tpl = _template(sn)
    b, t, s, (nb, nt, ns) = _blocks(sn)
    g = _prod_metric(sn, tpl, cspace)
    g_bulk, g_trace, g_seed = g[b], g[t], g[s]
    print(f"\n=== Shipped metric block ranges ===")
    print(f"  bulk  V·w   : [{g_bulk.min():.3e}, {g_bulk.max():.3e}]")
    print(f"  trace |Ω·n|w: [{g_trace.min():.3e}, {g_trace.max():.3e}]")
    print(f"  seed  V_cell: [{g_seed.min():.3e}, {g_seed.max():.3e}]")

    _fail_unless_above(np.min(g_bulk), 0.0, "bulk metric V·w must be strictly positive")
    _fail_unless_above(np.min(g_trace), 0.0, "trace metric |Ω·n|·w must be strictly positive")
    # The shipped seed metric is V_cell (SPD, nonzero) — NOT the zero ghost.
    _fail_unless_above(np.min(g_seed), 0.0,
                       "shipped seed metric must be V_cell (SPD), not the zero ghost")
    np.testing.assert_allclose(
        g_seed, _seed_to_coupled_layout(sn, _v_cell_seed(sn)), rtol=0.0, atol=0.0,
        err_msg="shipped seed metric block is not exactly V_cell")
    # ...so G is FULLY non-singular: every block contributes to the norm.
    n_pos = int((g > 0).sum())
    if n_pos != nb + nt + ns:
        pytest.fail(f"G should be fully non-singular: (g>0).sum()={n_pos}, expected {nb+nt+ns}")


def test_trace_and_pole_faces_both_hold_reciprocity():
    r"""Under the SHIPPED metric, activating the bulk, the SPATIAL-trace face
    (``|Ω·n|·w``), OR the ANGULAR-pole/seed face (``V_cell``) all leave G-adjoint
    reciprocity intact — the seed face is now a valid SPD metric, not a degeneracy.

    INVERTED from the pre-fix diagnostic (which asserted activating the pole/seed
    face BROKE reciprocity ``> 1e-3`` under the zero ghost, with a ``> 1e6``
    trace/pole asymmetry ratio — that ratio assertion is DELETED). Post-fix the
    seed-active defect is ``~6.6e-16``.

    Mutation tooth: reverting production ``for_levels`` to ``np.zeros`` makes the
    ``bulk + SEED`` leg jump to ``~1e-1`` → RED, while bulk-only and bulk+trace
    stay green (the seed metric is the only thing that changed).
    """
    sn, A = _build_sphere(4, 2, 0.5)
    A, cspace = _joint_M(sn, A)
    space = cspace

    def defect(*, bulk, trace, seed):
        rng = np.random.default_rng(2026)
        psi = _composite(sn, bulk=bulk, trace=trace, seed=seed, rng=rng)
        rng = np.random.default_rng(4052)
        phi = _composite(sn, bulk=bulk, trace=trace, seed=seed, rng=rng)
        lhs = space.inner_product(A.apply(psi), phi)
        rhs = space.inner_product(psi, A.H.apply(phi))
        return abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)

    d_bulk = defect(bulk=True, trace=False, seed=False)
    d_trace = defect(bulk=True, trace=True, seed=False)
    d_seed = defect(bulk=True, trace=False, seed=True)
    print("\n=== Reciprocity under the shipped metric ===")
    print(f"  bulk only           : {d_bulk:.3e}  (holds)")
    print(f"  bulk + TRACE active : {d_trace:.3e}  (holds — |Ω·n|·w)")
    print(f"  bulk + SEED  active : {d_seed:.3e}  (holds — V_cell is the SPD pole metric)")

    _fail_unless_below(d_bulk, 1e-11, "bulk-only reciprocity should hold")
    _fail_unless_below(d_trace, 1e-11, "trace-active reciprocity should hold (|Ω·n|·w correct)")
    # INVERTED: with V_cell shipped, activating the pole/seed face now HOLDS too.
    _fail_unless_below(d_seed, 1e-11,
                       "seed-active reciprocity now HOLDS (V_cell is the SPD pole metric)")


# ═══════════════════════════════════════════════════════════════════════
# The structural "why V_cell not zero" gate (metric-independent — it reads
# the FORWARD operator's block structure, unaffected by the fix)
# ═══════════════════════════════════════════════════════════════════════


def test_seed_selfblock_is_transport_not_reflection():
    r"""ψ½'s self-block ``A_ss`` is a BANDED radial-transport operator
    (:math:`\mu\partial_r + \sigma_t` — large off-diagonal cell-to-cell coupling),
    NOT a reflection/restriction map like the trace's ``A_tt`` (diag in ``[-1,1]``,
    tiny off-diagonal). This is WHY the seed's Hilbert metric is the radial volume
    ``V_cell`` (bulk-like), not the grazing angular-face measure 0.

    Metric-independent (it reads the forward operator ``A``, not ``G``) — green
    both before and after the fix. Also the #284 triangularity certificate:
    ``A_sb = A_st = 0`` (seed rows self-contained), ``A_bs ≠ 0`` (seed feeds bulk).

    Mutation tooth: were ``A_ss`` a face-restriction (a reflection map like the
    trace) its off-diagonal norm would collapse below 10 and the
    ``A_ss/A_tt`` interior-coupling ratio below 5 → RED. A seed row that read the
    bulk/trace (``A_sb`` or ``A_st`` ≠ 0) also REDS the triangularity pins.
    """
    sn, A = _build_sphere(6, 2, 0.5)
    A, cspace = _joint_M(sn, A)
    tpl = _template(sn)
    M = _dense(A.apply, tpl)
    b, t, s, (nb, nt, ns) = _blocks(sn)

    def offdiag(X):
        return float(np.linalg.norm(X - np.diag(np.diag(X))))

    A_tt, A_ss = M[t, t], M[s, s]
    print("\n=== Diagonal self-block character ===")
    for name, X in (("trace A_tt", A_tt), ("seed  A_ss", A_ss)):
        dr = (float(np.diag(X).min()), float(np.diag(X).max()))
        print(f"  {name}: ‖X‖={np.linalg.norm(X):7.3f}  diag∈[{dr[0]:6.3f},{dr[1]:6.3f}]  "
              f"offdiag‖·‖={offdiag(X):7.3f}")

    # The trace self-block is a reflection/restriction map: |diag| ≤ 1, tiny offdiag.
    _fail_unless_below(np.max(np.abs(np.diag(A_tt))), 1.0 + 1e-9,
                       "trace self-block should be a |·|≤1 map")
    _fail_unless_below(offdiag(A_tt), 5.0,
                       "trace self-block should have tiny off-diagonal (a restriction)")
    # The seed self-block is a banded radial-transport operator (interior dynamics).
    _fail_unless_above(offdiag(A_ss), 10.0,
                       "seed self-block should carry large radial-streaming off-diagonal")
    ratio = offdiag(A_ss) / max(offdiag(A_tt), 1e-300)
    print(f"  interior-coupling ratio offdiag(A_ss)/offdiag(A_tt) = {ratio:.2f}")
    _fail_unless_above(ratio, 5.0, "ψ½ carries far more interior dynamics than the trace")

    # Structural triangularity (#284 certificate): seed rows are self-contained.
    _fail_unless_below(np.linalg.norm(M[s, b]), 1e-12, "seed rows should be self-contained (A_sb=0)")
    _fail_unless_below(np.linalg.norm(M[s, t]), 1e-12, "seed rows should be self-contained (A_st=0)")
    # ...but the seed FEEDS the bulk (the M-M recurrence seed→ordinate-0).
    _fail_unless_above(np.linalg.norm(M[b, s]), 1e-3, "seed should feed the bulk (A_bs ≠ 0)")


# ═══════════════════════════════════════════════════════════════════════
# Gauge & install properties (R2, R3, R4)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("nx,ng,sigma", [(3, 2, 0.5), (5, 2, 1.0), (6, 3, 2.0)])
def test_r2_adjoint_bulk_is_gauge_invariant(nx, ng, sigma):
    r"""(R2 — observable-level gauge inertness) The physically-observable adjoint
    (bulk & trace rows of ``A.H``) is BITWISE-invariant under the seed gauge
    ``G_sd``. ``A.H = G⁻¹AᵀG`` is block-upper-triangular with the seed at the TOP
    (adjoint seed reads bulk; adjoint bulk reads nothing in the seed), so changing
    ``G_sd`` (identity vs V_cell vs 10·V_cell) leaves the bulk/trace rows exact
    while ONLY the internal seed rows move. This is the structural proof that
    ``G_sd`` is a genuine gauge — the ONLY physical requirement is non-degeneracy.

    Mutation tooth (gauge-inertness, note (c)): a corner-gauge change
    ``V[-1]→2·V[-1]`` (a subset of the identity↔V_cell↔10·V_cell sweep here)
    leaves the bulk/trace ``A.H`` rows GREEN (Δ = 0). Conversely a real bug that
    let the bulk adjoint READ the seed metric would make ``d_bulk`` move O(1) → RED.
    """
    sn, A = _build_sphere(nx, ng, sigma)
    A, cspace = _joint_M(sn, A)
    tpl = _template(sn)
    b, t, s, (nb, nt, ns) = _blocks(sn)
    T = _dense(A.apply_transpose, tpl)
    g0 = _prod_metric(sn, tpl, cspace)

    def adjoint_matrix(seed_diag):
        g = g0.copy()
        g[s] = seed_diag
        gpos = g > 0
        ginv = np.where(gpos, 1.0 / np.where(gpos, g, 1.0), 0.0)
        return np.diag(ginv) @ T @ np.diag(g)

    V = _v_cell_seed(sn)
    H_id = adjoint_matrix(np.ones(ns))
    H_V = adjoint_matrix(V)
    H_10V = adjoint_matrix(10.0 * V)

    bt = np.r_[np.arange(nb), np.arange(nb, nb + nt)]
    d_bulk_V = np.max(np.abs(H_id[bt, :] - H_V[bt, :]))
    d_bulk_10 = np.max(np.abs(H_id[bt, :] - H_10V[bt, :]))
    d_seed = np.max(np.abs(H_id[s, :] - H_V[s, :]))
    print(f"  nx={nx} ng={ng}: bulk⊕trace A.H gauge-Δ = {max(d_bulk_V, d_bulk_10):.3e} "
          f"| seed rows gauge-Δ = {d_seed:.3e}")
    _fail_unless_below(d_bulk_V, 1e-11, f"bulk A.H moved under gauge (V_cell): {d_bulk_V:.2e}")
    _fail_unless_below(d_bulk_10, 1e-11, f"bulk A.H moved under gauge (10·V_cell): {d_bulk_10:.2e}")
    _fail_unless_above(d_seed, 1e-6,
                       "seed A.H did NOT move under gauge — the gauge should be non-vacuous")


@pytest.mark.parametrize("nx,ng,sigma", [(3, 2, 0.5), (5, 2, 1.0), (8, 2, 0.3), (6, 3, 2.0)])
def test_r3_vcell_closes_mode12_across_configs(nx, ng, sigma):
    r"""(R3 — across-config closure) The recommended ``G_sd = V_cell`` closes
    Mode-12 at every config: with V_cell installed and RANDOM seed data, (a)
    reciprocity holds, and (b) a seed-row (``A_ss``) sign flip REDS it.

    Mutation tooth: the seed-row flip is the live tooth (``d_flip`` must exceed
    ``1e-3``); reverting to ``G_sd=0`` would make the flip blind (``d_flip`` would
    not exceed the pre-flip defect) — the across-config generalisation of the
    canonical Mode-12 gate's leg (c).
    """
    sn, A = _build_sphere(nx, ng, sigma)
    A, cspace = _joint_M(sn, A)
    tpl = _template(sn)
    b, t, s, (nb, nt, ns) = _blocks(sn)
    Afwd = _dense(A.apply, tpl)
    T = _dense(A.apply_transpose, tpl)
    g = _prod_metric(sn, tpl, cspace)
    g[s] = _v_cell_seed(sn)
    gpos = g > 0
    ginv = np.where(gpos, 1.0 / np.where(gpos, g, 1.0), 0.0)
    G, Gpinv = np.diag(g), np.diag(ginv)
    H = Gpinv @ T @ G
    rng = np.random.default_rng(99)
    n = Afwd.shape[0]
    psi, phi = rng.standard_normal(n), rng.standard_normal(n)

    def defect(Amat):
        lhs = (Amat @ psi) @ (G @ phi)
        rhs = psi @ (G @ (H @ phi))
        return abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)

    d_ok = defect(Afwd)
    A_flip = Afwd.copy()
    A_flip[s, s] *= -1.0
    d_flip = defect(A_flip)
    print(f"  nx={nx} ng={ng} σ={sigma}: recip(V_cell)={d_ok:.2e}  flip-reds={d_flip:.2e}")
    _fail_unless_below(d_ok, 1e-11, f"V_cell reciprocity should hold ({d_ok:.2e})")
    _fail_unless_above(d_flip, 1e-3, f"seed-row flip should RED under V_cell ({d_flip:.2e})")


@pytest.mark.parametrize("nx,ng,sigma", [(4, 2, 0.5), (6, 3, 1.0)])
def test_r4_gauge_perturbing_gsd_leaves_forward_bit_identical(nx, ng, sigma):
    r"""(R4 — install/gauge de-risk) The seed metric is read ONLY by the adjoint
    path (``A.H``) and ``inner_product``, NEVER the forward sweep/matvec — so
    gauge-perturbing ``G_sd`` leaves EVERY forward result BYTE-for-byte unchanged
    (the same guarantee the O.2b trace-metric install carried, #208).

    Proven by monkeypatching the cached seed space's ``inner_product_weights``
    with an explicit non-vacuous gauge metric (``2·V_cell`` — a uniform gauge
    scaling that includes the corner ``V[-1]→2·V[-1]`` tooth (c)) and confirming
    ``A.apply`` on a random composite is bit-identical. The perturbation is a
    freshly-computed ``2·V_cell`` (NOT ``2×`` the shipped value), so the gate is a
    pure forward-blindness probe, decoupled from the shipped metric's value.

    Mutation tooth: if the forward path READ the metric, ``y_after`` would differ
    from ``y_before`` → the ``assert_array_equal`` REDS. The gauge change here is
    the observable-inert side of note (c): the forward stays GREEN under it.
    """
    sn, A = _build_sphere(nx, ng, sigma)
    A, cspace = _joint_M(sn, A)
    tpl = _template(sn)
    n = tpl.to_flat().size
    rng = np.random.default_rng(5)
    x = type(tpl).from_flat(rng.standard_normal(n), tpl)

    y_before = A.apply(x).to_flat()

    # A non-vacuous GAUGE change on the (cached) seed metric: an explicit 2·V_cell
    # (includes the corner V[-1]→2·V[-1]). The forward path must never read the
    # metric, so the result is byte-for-byte unchanged.
    interior_space = sn.radial_characteristic_interior_space
    boundary_space = sn.radial_characteristic_boundary_space
    ni = interior_space.shape[0]
    perturbed = 2.0 * _v_cell_seed(sn)
    object.__setattr__(interior_space, "inner_product_weights", perturbed[:ni])
    object.__setattr__(boundary_space, "inner_product_weights", perturbed[ni:])

    y_after = A.apply(x).to_flat()
    np.testing.assert_array_equal(
        y_before, y_after,
        err_msg="forward A.apply changed when the seed metric was gauge-perturbed — "
                "the forward path MUST NOT read the metric")
    # sanity: the metric really did change (non-vacuous perturbation).
    _fail_unless_above(np.max(np.abs(perturbed)), 0.0,
                       "perturbed seed metric should be nonzero (non-vacuous gauge change)")
