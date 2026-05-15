r"""Tests for the Step 2.5c two-stratum sweep cache (Issue #196 Phase G).

Twelve tests in five thematic groups (per plan §"Test catalog"):

* **Cache structure** (#1-3): :class:`GeometryCoefficients` + :class:`CollisionCache`
  populate the expected fields with the expected shapes; the two strata are
  separate by ``ng`` axis.
* **Cache-invariance** (#4-5 — the CARDINAL tests): the collision cache is
  built EXACTLY ONCE across a 5+ iteration Picard fixed-point;
  ``rebind_cross_sections`` invalidates only :class:`CollisionCache`.
* **Dual-view consistency** (#6-7): the cache-driven sweep result matches the
  per-cell ``cell_update.update`` iteration to ``rtol=1e-13`` (Pattern 2).
* **Performance gates** (#8-9): slab benchmark ≤ 1.5 ms; full SN suite < 5 min.
* **Production gates** (#10-12): L0 streaming-equilibrium, regression
  snapshots, Step 2.5b's pair-monoid associativity all stay green.

Per ``vv-principles`` §"Bit-identity vs principled-equivalence": the
cache-driven sweep produces algebraically the same value as the per-cell
``update`` reference; bit-identity is preserved on slab (Blelloch §1.5 closed
form is the SAME numpy expression both paths use) and the dual-view contract
is the L1 cross-check that catches any reduction drift.
"""

from __future__ import annotations

import time
from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn.geometry import SNMesh
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.spatial.cell_balance import cell_balance_terms
from orpheus.sn.spatial.cell_update import UpstreamState
from orpheus.sn.spatial.diamond import DiamondDifference
from orpheus.sn.spatial.scan import ordinate_scan
from orpheus.sn.spatial.sweep_cache import CollisionCache, GeometryCoefficients


# ═══════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════


def _make_slab(nx: int = 10, N: int = 8) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = GaussLegendre1D.create(N)
    return SNMesh(mesh, quad)


def _make_sphere(nx: int = 10, N: int = 8) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = GaussLegendre1D.create(N)
    return SNMesh(mesh, quad)


# ═══════════════════════════════════════════════════════════════════════
# Group 1 — Cache structure (tests #1-3)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l0
def test_geometry_coefficients_built_at_construction() -> None:
    """Test #1 — :class:`GeometryCoefficients` populates every field.

    All Stratum-1 fields present; shapes match the ``(N, nx)`` contract; the
    frozen dataclass refuses post-construction mutation.
    """
    sn_mesh = _make_slab(nx=10, N=8)
    geom = GeometryCoefficients.from_mesh_and_quad(sn_mesh)
    N, nx = 8, 10
    assert geom.chain_idx.shape == (N, nx)
    assert geom.chain_idx_inv.shape == (N, nx)
    assert geom.abs_mu.shape == (N,)
    assert geom.A_down.shape == (N, nx)
    assert geom.A_total.shape == (N, nx)
    assert geom.dA_w.shape == (N, nx)
    assert geom.V.shape == (N, nx)
    assert geom.c_in.shape == (N,)
    assert geom.c_out.shape == (N,)
    assert geom.tau_inv.shape == (N,)
    assert geom.mm_a_in_coeff.shape == (N,)
    assert geom.is_degenerate.shape == (N,)
    assert geom.is_degenerate.dtype == bool
    # Frozen dataclass — refuses re-binding any field
    with pytest.raises(FrozenInstanceError):
        geom.A_down = np.zeros_like(geom.A_down)  # type: ignore[misc]


@pytest.mark.l0
def test_collision_cache_built_at_sigma_t_bind() -> None:
    """Test #2 — :class:`CollisionCache` shape consistent with ``geom``; values match the formula.

    Hand-evaluate the cache against the canonical formula
    ``inverse_denom = 1 / (2|μ|·A_down + dA_w·c_out + σ_t·V)``
    cell-by-cell for one ordinate on a slab fixture (where dA_w=0, A_down=1,
    A_total=2) — the closed-form ``a = (2|μ|·2 − σ_t·V)/(2|μ|·1 + σ_t·V)``.

    Cache storage layout is ``(N, ng, nx)`` under Issue #196 PR-INDEX-2.
    """
    sn_mesh = _make_slab(nx=4, N=4)
    geom = GeometryCoefficients.from_mesh_and_quad(sn_mesh)
    # sig_t is (ng, nx) under PR-INDEX-2.  Two groups × four cells,
    # uniform per group: group 0 has σ_t=1.0, group 1 has σ_t=2.0.
    sig_t = np.array([[1.0] * 4, [2.0] * 4])  # (ng=2, nx=4)
    coll = CollisionCache.from_geometry(geom, sig_t)

    # (N, ng, nx) — N=4 ordinates, ng=2 groups, nx=4 cells.
    assert coll.inverse_denom.shape == (4, 2, 4)
    assert coll.a_attenuation.shape == (4, 2, 4)
    assert coll.cumprod_a.shape == (4, 2, 4)

    # Hand-eval for ordinate n=0 (most-inward μ), group 0, cell 0
    # (in chain order):  |μ|=|mu_x[0]|, V=mesh width = 0.25, σ_t=1.0.
    quad = sn_mesh.quad
    mu = abs(float(quad.mu_x[0]))
    V = 0.25  # uniform mesh
    sig = 1.0
    denom_expected = 2.0 * mu * 1.0 + 0.0 + sig * V  # dA_w=0 slab
    a_expected = 2.0 * mu * 2.0 / denom_expected - 1.0
    # Indexing semantics: [n=0, g=0, i=0] under (N, ng, nx) layout.
    assert np.isclose(coll.inverse_denom[0, 0, 0], 1.0 / denom_expected, rtol=1e-14)
    assert np.isclose(coll.a_attenuation[0, 0, 0], a_expected, rtol=1e-14)


@pytest.mark.l0
def test_two_strata_independence_by_ng_axis() -> None:
    """Test #3 — Stratum 1 has NO ``ng`` axis; Stratum 2 has ``(N, ng, nx)``.

    Proves the strata are deliberately separated by mutation cadence — a
    Smell #16 instance would surface as Stratum-1 carrying an ``ng`` axis.
    Cache storage layout is ``(N, ng, nx)`` under Issue #196 PR-INDEX-2.
    """
    sn_mesh = _make_slab(nx=5, N=4)
    geom = GeometryCoefficients.from_mesh_and_quad(sn_mesh)
    # Stratum 1 — no ng axis on ANY field.
    for name in ("A_down", "A_total", "dA_w", "V"):
        field_arr = getattr(geom, name)
        assert field_arr.ndim == 2, f"{name} should be (N, nx); got shape {field_arr.shape}"
    for name in ("abs_mu", "c_in", "c_out", "tau_inv", "mm_a_in_coeff", "is_degenerate"):
        field_arr = getattr(geom, name)
        assert field_arr.ndim == 1, f"{name} should be (N,); got shape {field_arr.shape}"

    # Stratum 2 — every tensor has the (N, ng, nx) shape (PR-INDEX-2).
    sig_t = np.ones((3, 5))  # (ng=3, nx=5) under PR-INDEX-2
    coll = CollisionCache.from_geometry(geom, sig_t)
    for name in ("inverse_denom", "a_attenuation", "cumprod_a"):
        field_arr = getattr(coll, name)
        assert field_arr.shape == (4, 3, 5), f"{name} shape {field_arr.shape} != (4, 3, 5)"


# ═══════════════════════════════════════════════════════════════════════
# Group 2 — Cache invariance (tests #4-5)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l0
def test_collision_cache_invariance_under_source_iteration() -> None:
    """Test #4 (CARDINAL) — :meth:`CollisionCache.from_geometry` called ONCE per σ_t epoch.

    Run a 5-iteration Picard fixed-point on a homogeneous slab eigenvalue
    problem.  The inner SI updates ``Q`` each iteration; the outer power
    iteration updates the fission source.  Across both loops, σ_t is bound
    once at ``SNSolver.__init__`` and the cache MUST be built exactly once.

    This is the cardinal invariance gate: it proves the cache placement
    on :class:`SNSolver` survives Picard, and that no sweep path is
    secretly rebuilding the cache on every iteration.
    """
    from orpheus.data.macro_xs.mixture import Mixture
    from orpheus.sn.solver import solve_sn
    from scipy.sparse import csr_matrix

    # Single-group multiplying mixture: c = σ_s/σ_t = 0.5, νΣ_f = 0.3.
    mix = Mixture(
        SigT=np.array([1.0]),
        SigC=np.array([0.5]),                       # capture
        SigL=np.array([0.0]),                       # leakage / (n, 2n) loss
        SigF=np.array([0.3]),
        SigP=np.array([0.3]),                       # ν · Σ_f (production)
        SigS=[csr_matrix(np.array([[0.5]]))],
        Sig2=csr_matrix(np.array([[0.0]])),
        chi=np.array([1.0]),
    )
    materials = {0: mix}
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 11),
        mat_ids=np.zeros(10, dtype=int),
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    quad = GaussLegendre1D.create(4)

    # Reset counter, then run a converged eigenvalue (≥ 5 outer × N inner).
    CollisionCache.reset_build_count()
    result = solve_sn(
        materials=materials,
        mesh=mesh,
        quadrature=quad,
        inner_solver="source_iteration",
        max_outer=20,
        max_inner=50,
        inner_tol=1e-8,
        keff_tol=1e-6,
        flux_tol=1e-5,
    )
    assert len(result.keff_history) >= 5, (
        "Test fixture is too trivial — converged in fewer than 5 outer "
        "iterations; raise max_outer to exercise the Picard loop."
    )
    # ``solve_sn`` constructs ONE fresh SNSolver, which builds the cache
    # exactly once at ``__init__``.  Subsequent sweeps (every outer ×
    # every inner ≥ 5 × 5 = 25 sweep calls) MUST NOT rebuild it.
    assert CollisionCache._build_count == 1, (
        f"CollisionCache rebuilt {CollisionCache._build_count} times — "
        f"expected exactly 1 across {len(result.keff_history)} outer × "
        f"~tens of inner iterations.  Some sweep path is re-instantiating "
        f"the cache on every iteration."
    )


@pytest.mark.l0
def test_geometry_coefficients_invariance_under_sigma_t_change() -> None:
    """Test #5 — :class:`GeometryCoefficients` survives ``rebind_cross_sections``.

    After ``solver.rebind_cross_sections(new_sig_t)``, the geometry cache
    is the SAME object (identity check).  Only :class:`CollisionCache`
    rebuilds.
    """
    from orpheus.data.macro_xs.mixture import Mixture
    from orpheus.sn.solver import SNSolver
    from scipy.sparse import csr_matrix

    mix = Mixture(
        SigT=np.array([1.0]),
        SigC=np.array([0.5]),
        SigL=np.array([0.0]),
        SigF=np.array([0.0]),
        SigP=np.array([0.0]),
        SigS=[csr_matrix(np.array([[0.5]]))],
        Sig2=csr_matrix(np.array([[0.0]])),
        chi=np.array([1.0]),
    )
    materials = {0: mix}
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 6),
        mat_ids=np.zeros(5, dtype=int),
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = GaussLegendre1D.create(4)
    sn_mesh = SNMesh(mesh, quad)
    solver = SNSolver(materials=materials, sn_mesh=sn_mesh)

    geom_before = solver.geom_cache
    coll_before = solver.coll_cache
    new_sig_t = solver.sig_t * 2.0
    solver.rebind_cross_sections(new_sig_t)
    geom_after = solver.geom_cache
    coll_after = solver.coll_cache

    assert geom_after is geom_before, (
        "GeometryCoefficients should be invariant under σ_t rebinds; "
        "rebind_cross_sections accidentally invalidated Stratum 1."
    )
    assert coll_after is not coll_before, (
        "CollisionCache should be rebuilt on σ_t rebind."
    )


# ═══════════════════════════════════════════════════════════════════════
# Group 3 — Dual-view consistency (tests #6-7)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l0
@pytest.mark.parametrize("geometry", ["slab", "sphere"])
@pytest.mark.parametrize("ng", [1, 2, 3])
@pytest.mark.parametrize("source_kind", ["uniform", "linear", "gaussian"])
def test_cache_driven_sweep_matches_per_cell_update(
    geometry: str, ng: int, source_kind: str,
) -> None:
    """Test #6 — cache-driven sweep equals per-cell ``cell_update.update``.

    For each (geometry × ng × source) combination, run ONE ordinate's
    spatial scan through:

    1. The cache-driven path: ``b = 2·QV·inverse_denom`` (slab) /
       ``b = 2·(QV + ang_contrib)·inverse_denom`` (curvilinear), then
       ``psi_face = ordinate_scan(a, b, psi_in)``;
    2. The per-cell reference path: ``update(visit, ...)`` cell-by-cell.

    The two MUST agree at ``rtol=1e-13`` (Pattern 2 dual-view contract).
    Disagreement signals either a cache-populator algebra bug or a
    ``cell_balance_terms`` non-affine term.  Per plan decision-point
    checkpoints: STOP and trace algebraically; do NOT widen rtol.
    """
    nx = 8
    N = 4
    if geometry == "slab":
        sn_mesh = _make_slab(nx=nx, N=N)
    else:
        sn_mesh = _make_sphere(nx=nx, N=N)

    geom = GeometryCoefficients.from_mesh_and_quad(sn_mesh)
    rng = np.random.default_rng(42)
    # Issue #196 PR-INDEX-2: cache consumes σ_t as (ng, nx).
    sig_t = 1.0 + 0.5 * rng.random((ng, nx))                  # (ng, nx)
    coll = CollisionCache.from_geometry(geom, sig_t)

    # Build a representative source in (ng, nx) — principled layout.
    if source_kind == "uniform":
        Q = np.ones((ng, nx))
    elif source_kind == "linear":
        Q = np.ones((ng, 1)) * (np.arange(nx)[None, :] + 1.0)
    else:
        x = np.linspace(0, 1, nx)
        Q = np.ones((ng, 1)) * np.exp(-((x - 0.5) ** 2) * 10)[None, :]

    # Pick the most-inward ordinate (μ<0) so chain order is reverse.
    quad = sn_mesh.quad
    mu = quad.mu_x
    n = int(np.argmin(mu))

    # ψ_a_in is an INPUT to one ordinate's spatial scan (it's the
    # frozen output of the previous-ordinate's M-M angular thread).
    # Both paths consume the SAME pre-computed (ng, nx) chain-ordered
    # array — for the dual-view contract, the cache and the reference
    # MUST be driven by the same input ψ_a_in.
    chain = geom.chain_idx[n]
    V_chain = geom.V[n]                                       # (nx,)
    # Q is (ng, nx); chain reorders the nx axis to chain order.
    QV_chain = Q[:, chain] * V_chain[chain] / quad.weights.sum()  # (ng, nx)
    psi_in = np.zeros(ng)
    if geometry == "sphere":
        # Frozen ψ_a_in input — a representative bump so the curvature
        # term is genuinely exercised (zeros would null dA_w·c_in·ψ_a_in).
        rng2 = np.random.default_rng(7)
        psi_a_in_chain = 0.1 * rng2.random((ng, nx))          # (ng, nx)
        ang_contrib = (geom.dA_w[n] * geom.c_in[n])[None, :] * psi_a_in_chain  # (ng, nx)
        b = 2.0 * (QV_chain + ang_contrib) * coll.inverse_denom[n]  # (ng, nx)
    else:
        psi_a_in_chain = None
        b = 2.0 * QV_chain * coll.inverse_denom[n]            # (ng, nx)
    # ordinate_scan expects scan axis (cell axis) leading — transpose
    # the principled (ng, nx) views to (nx, ng) at the call.
    psi_face_chain_fast = ordinate_scan(
        coll.a_attenuation[n].T, b.T, psi_in,
    )                                                          # (nx, ng)

    # Reference path — per-cell update over the same chain, feeding the
    # SAME frozen ψ_a_in_chain per cell (NOT the M-M output).  This is
    # the dual-view contract: at fixed input state, the cache-driven
    # scan equals the per-cell update.
    dd = DiamondDifference()
    psi_face_chain_ref = np.full_like(psi_face_chain_fast, np.nan)
    has_downstream = np.zeros(nx, dtype=bool)
    psi_face_in = psi_in
    visits_full = list(sn_mesh.dag_walk(ordinate_idx=n))
    for k_chain in range(nx):
        cell_i = int(chain[k_chain])
        visit = visits_full[k_chain]
        psi_a_in_cell = (
            psi_a_in_chain[:, k_chain] if psi_a_in_chain is not None else None
        )
        upstream = UpstreamState(
            spatial_upstream=psi_face_in,
            angular_upstream=psi_a_in_cell,
        )
        result = dd.update(
            visit=visit,
            total_xs=sig_t[:, cell_i],                         # (ng,) — group axis at axis 0
            source=QV_chain[:, k_chain],                       # (ng,)
            upstream_state=upstream,
        )
        if result.outgoing_spatial_flux is not None:
            psi_face_chain_ref[k_chain] = result.outgoing_spatial_flux
            psi_face_in = result.outgoing_spatial_flux
            has_downstream[k_chain] = True

    # Dual-view contract — only cells with a downstream spatial face
    # are part of the chain's spatial scan.  The last sphere-inward
    # cell (at r=0) has ``A_down == 0`` and contributes no scan output;
    # ``update`` returns ``None`` there and the cache's formal scan
    # output is undefined.  Compare only the well-defined cells.
    np.testing.assert_allclose(
        psi_face_chain_fast[has_downstream],
        psi_face_chain_ref[has_downstream],
        rtol=1e-13,
        err_msg=(
            f"Dual-view disagreement on {geometry}/ng={ng}/{source_kind} "
            f"at ordinate n={n}: cache-driven and per-cell update "
            f"diverge.  Per plan decision-point: trace algebraically."
        ),
    )
    # Must have well-defined chain on all but at most one cell (the
    # innermost sphere cell on an inward sweep; slab has none).
    assert has_downstream.sum() >= nx - 1


@pytest.mark.l0
def test_cache_populator_matches_cell_balance_terms() -> None:
    """Test #7 — cache ``(a, 1/denom)`` agrees with :func:`cell_balance_terms`.

    Pattern 2 anchor at ``rtol=1e-14``.  For any cell, the cache's per-cell
    ``(a, inverse_denom)`` MUST equal what ``cell_balance_terms`` would
    produce — the two paths derive from the same algebra.
    """
    sn_mesh = _make_sphere(nx=8, N=4)
    geom = GeometryCoefficients.from_mesh_and_quad(sn_mesh)
    ng = 2
    # Issue #196 PR-INDEX-2: cache consumes σ_t as (ng, nx).
    # Build (nx, ng) first via outer product for readability, then transpose.
    sig_t_xg = np.linspace(0.5, 1.5, 8)[:, None] * np.array([[1.0, 2.0]])  # (8, 2)
    sig_t = sig_t_xg.T                                                     # (ng=2, nx=8)
    coll = CollisionCache.from_geometry(geom, sig_t)

    # Sample two ordinates × two cells (chain positions).
    quad = sn_mesh.quad
    for n in (int(np.argmin(quad.mu_x)), int(np.argmax(quad.mu_x))):
        visits = list(sn_mesh.dag_walk(ordinate_idx=n))
        chain = geom.chain_idx[n]
        for k_chain in (0, 3, 7):
            visit = visits[k_chain]
            cell_i = int(chain[k_chain])
            # Use a zero upstream + zero angular probe — the cell_balance_terms
            # output's ``denom`` is independent of upstream, so this is fine.
            upstream = UpstreamState(
                spatial_upstream=np.zeros(ng),
                angular_upstream=np.zeros(ng),
            )
            terms = cell_balance_terms(
                visit.streaming_terms,
                visit.face_area_downstream,
                sig_t[:, cell_i],                                          # (ng,)
                upstream,
            )
            # Cache layout (N, ng, nx) — fix n, fix k_chain ⇒ (ng,) vector.
            # Indexing pattern updated from [n, k_chain] (legacy axis 1 was
            # cell) to [n, :, k_chain] (PR-INDEX-2 axis 1 is group, axis 2
            # is cell).  The semantic intent is the same: "per-cell denom
            # vector across groups."
            denom_cached = 1.0 / coll.inverse_denom[n, :, k_chain]         # (ng,)
            a_cached = coll.a_attenuation[n, :, k_chain]                   # (ng,)
            # Algebraic a from terms: a = 2|μ|·A_total / denom − 1.
            A_total = (
                visit.streaming_terms.face_area_inner
                + visit.streaming_terms.face_area_outer
            )
            a_expected = (
                2.0 * visit.streaming_terms.abs_mu * A_total / terms.denom - 1.0
            )
            np.testing.assert_allclose(
                denom_cached, terms.denom, rtol=1e-14,
                err_msg=f"denom mismatch n={n} k={k_chain}",
            )
            np.testing.assert_allclose(
                a_cached, a_expected, rtol=1e-14,
                err_msg=f"a mismatch n={n} k={k_chain}",
            )


# ═══════════════════════════════════════════════════════════════════════
# Group 4 — Performance gates (tests #8-9)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l0
@pytest.mark.slow
def test_slab_sweep_benchmark_under_2ms() -> None:
    """Test #8 — slab sweep ``nx=160 N=16 ng=4`` runs in ≤ 2 ms.

    Step 2.5b baseline: 15.43 ms/sweep.  Target: ≤ 1.5 ms (10× speedup);
    acceptance gate: ≤ 2.0 ms (a slim safety margin for CI machine noise).
    Marked ``@slow`` — skipped by default but runs in CI.
    """
    from orpheus.sn.sweep import transport_sweep

    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 161),
        mat_ids=np.zeros(160, dtype=int),
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = GaussLegendre1D.create(16)
    sn_mesh = SNMesh(mesh, quad)
    Q = np.ones((160, 1, 4))
    sig_t = np.ones((4, 160, 1))  # (ng, nx, ny) — PR-INDEX-3
    psi_bc: dict = {}

    # Warm-up — first call also caches inside SNMesh.
    for _ in range(3):
        transport_sweep(Q, sig_t, sn_mesh, psi_bc)

    # Measured wall clock over 100 sweeps.
    n_iters = 100
    t0 = time.perf_counter()
    for _ in range(n_iters):
        transport_sweep(Q, sig_t, sn_mesh, psi_bc)
    elapsed_per_sweep_ms = (time.perf_counter() - t0) / n_iters * 1000.0
    print(f"\nSlab sweep nx=160 N=16 ng=4: {elapsed_per_sweep_ms:.3f} ms/sweep")
    assert elapsed_per_sweep_ms < 2.0, (
        f"Slab sweep at {elapsed_per_sweep_ms:.3f} ms/sweep exceeds "
        f"the 2.0 ms gate (target ≤ 1.5 ms; Step 2.5b baseline 15.43 ms). "
        "Profile and find the per-cell hot path."
    )


@pytest.mark.l0
@pytest.mark.slow
def test_full_sn_suite_under_5min() -> None:
    """Test #9 — full ``tests/sn`` suite runs in < 5 min.

    Placeholder marker — the gate is exercised at the closeout layer via
    ``time pytest tests/sn/ -q``.  This test is a tag carrier so CI can
    enable / disable the performance gate via marker selection.
    """
    pytest.skip("Performance gate is measured externally via `time pytest tests/sn/ -q`.")


# ═══════════════════════════════════════════════════════════════════════
# Group 5 — Production gates (tests #10-12)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l0
def test_l0_streaming_equilibrium_preserved_after_2_5c() -> None:
    """Test #10 — homogeneous reflective slab streaming-equilibrium φ → Q/σ_t.

    The canonical L0 invariant: in a homogeneous reflective medium with
    uniform source ``Q`` and uniform total cross section ``σ_t``, the
    scalar flux converges to ``Q/σ_t`` to machine precision regardless of
    quadrature family or geometry.  Step 2.5c is structural ONLY — this
    must stay green.
    """
    from orpheus.data.macro_xs.mixture import Mixture
    from orpheus.sn.solver import solve_sn_fixed_source
    from scipy.sparse import csr_matrix

    sigma_t = 2.0
    mix = Mixture(
        SigT=np.array([sigma_t]),
        SigC=np.array([sigma_t]),                   # all absorption is capture
        SigL=np.array([0.0]),
        SigF=np.array([0.0]),
        SigP=np.array([0.0]),
        SigS=[csr_matrix(np.array([[0.0]]))],
        Sig2=csr_matrix(np.array([[0.0]])),
        chi=np.array([1.0]),
    )
    materials = {0: mix}
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 11),
        mat_ids=np.zeros(10, dtype=int),
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    quad = GaussLegendre1D.create(4)
    # external_source has shape (N, nx, ny, ng); uniform Q=1.0 per ordinate.
    Q_external = np.ones((quad.N, 10, 1, 1))

    result = solve_sn_fixed_source(
        materials=materials,
        mesh=mesh,
        quadrature=quad,
        external_source=Q_external,
        boundary_condition="reflective",
        max_inner=500,
        inner_tol=1e-12,
    )
    # In a homogeneous reflective medium, φ = Q/σ_t where Q is the
    # per-cell isotropic source.  external_source carries Q_n=1.0 per
    # ordinate; the source-iteration form is q = ∑_n w_n · Q_n / W
    # → φ → 1.0 / σ_t.
    expected = 1.0 / sigma_t
    np.testing.assert_allclose(result.scalar_flux[:, 0, 0], expected, rtol=1e-10)


@pytest.mark.l0
def test_pair_monoid_associativity_still_passes() -> None:
    """Test #12 — Step 2.5b's pair-monoid associativity invariant survives.

    The Blelloch §1.5 closed form rests on associativity of the
    transmission-emission pair-monoid.  This test re-runs the algebraic
    anchor; if it fails, ``ordinate_scan`` is mis-wired.
    """
    rng = np.random.default_rng(0)
    nx = 100
    ng = 1
    a = 0.5 + 0.5 * rng.random((nx, ng))
    b = rng.random((nx, ng))
    psi_0 = np.zeros(ng)

    # Direct ordinate_scan
    psi_direct = ordinate_scan(a, b, psi_0)

    # Pair-monoid composition: ψ_{n+1} = a_n · ψ_n + b_n
    psi_loop = np.empty_like(psi_direct)
    psi = psi_0.copy()
    for i in range(nx):
        psi = a[i] * psi + b[i]
        psi_loop[i] = psi

    np.testing.assert_allclose(psi_direct, psi_loop, rtol=1e-13)


# Test #11 (regression snapshots bit-identical) is exercised externally
# via ``pytest tests/sn/regression/`` — listed in the closeout test-pin
# block for verbatim paste-back.
