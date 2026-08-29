r"""Tests for the Step 2.5c two-stratum sweep cache (Issue #196 Phase G).

Twelve tests in five thematic groups (per plan §"Test catalog"):

* **Cache structure** (#1-3): :class:`StreamingCoefficientCache` + :class:`CollisionCache`
  populate the expected fields with the expected shapes; the two strata are
  separate by ``ng`` axis.
* **Cache-invariance** (#4-5 — the CARDINAL tests): the collision cache is
  built EXACTLY ONCE across a 5+ iteration Picard fixed-point;
  ``rebind_cross_sections`` invalidates only :class:`CollisionCache`.
* **Dual-view consistency** (#6-7): the cache-driven sweep result matches the
  per-cell ``scheme.update`` iteration to ``rtol=1e-13`` (Pattern 2).
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
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.spatial.cell_balance import cell_balance_for_streaming
from orpheus.transport.spatial.scheme import UpstreamState
from orpheus.transport.spatial.diamond import DiamondDifference
from orpheus.sn.sweep.scan import ordinate_scan
from orpheus.sn.sweep.cache import CollisionCache, StreamingCoefficientCache
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux


# ═══════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════


def _trivial_materials(ng: int = 1) -> dict:
    """Build a minimal materials dict for ng-group geometry-only tests.

    Issue #197 PR-TYPED-0: SNMesh requires ``materials``.  This helper
    provides a placeholder mixture for tests that exercise pure-
    geometry/cache structure and don't consume cross-section values.
    """
    from orpheus.data.macro_xs.mixture import Mixture
    from scipy.sparse import csr_matrix
    z = np.zeros(ng)
    z_mat = csr_matrix(np.zeros((ng, ng)))
    return {0: Mixture(
        SigC=z.copy(), SigL=z.copy(), SigF=z.copy(),
        SigP=z.copy(), SigT=np.ones(ng),
        SigS=[z_mat], Sig2=z_mat, chi=z.copy(),
    )}


def _make_slab(nx: int = 10, N: int = 8) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(N)
    return SNMesh(mesh, quad, _trivial_materials(ng=1))


def _make_sphere(nx: int = 10, N: int = 8) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(N)
    return SNMesh(mesh, quad, _trivial_materials(ng=1))


# ═══════════════════════════════════════════════════════════════════════
# Group 1 — Cache structure (tests #1-3)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l0
def test_geometry_coefficients_built_at_construction() -> None:
    """Test #1 — :class:`StreamingCoefficientCache` populates every field.

    All Stratum-1 fields present; shapes match the ``(N, nx)`` contract; the
    frozen dataclass refuses post-construction mutation.
    """
    sn_mesh = _make_slab(nx=10, N=8)
    geom = StreamingCoefficientCache.from_mesh_and_quad(sn_mesh, sn_mesh.angular_closure)
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
    geom = StreamingCoefficientCache.from_mesh_and_quad(sn_mesh, sn_mesh.angular_closure)
    # sig_t is (ng, nx) under PR-INDEX-2.  Two groups × four cells,
    # uniform per group: group 0 has σ_t=1.0, group 1 has σ_t=2.0.
    sig_t = np.array([[1.0] * 4, [2.0] * 4])  # (ng=2, nx=4)
    coll = CollisionCache.from_geometry(geom, sig_t, sn_mesh.scheme)

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
    geom = StreamingCoefficientCache.from_mesh_and_quad(sn_mesh, sn_mesh.angular_closure)
    # Stratum 1 — no ng axis on ANY field.
    for name in ("A_down", "A_total", "dA_w", "V"):
        field_arr = getattr(geom, name)
        assert field_arr.ndim == 2, f"{name} should be (N, nx); got shape {field_arr.shape}"
    for name in ("abs_mu", "c_in", "c_out", "tau_inv", "mm_a_in_coeff", "is_degenerate"):
        field_arr = getattr(geom, name)
        assert field_arr.ndim == 1, f"{name} should be (N,); got shape {field_arr.shape}"

    # Stratum 2 — every tensor has the (N, ng, nx) shape (PR-INDEX-2).
    sig_t = np.ones((3, 5))  # (ng=3, nx=5) under PR-INDEX-2
    coll = CollisionCache.from_geometry(geom, sig_t, sn_mesh.scheme)
    for name in ("inverse_denom", "a_attenuation", "cumprod_a"):
        field_arr = getattr(coll, name)
        assert field_arr.shape == (4, 3, 5), f"{name} shape {field_arr.shape} != (4, 3, 5)"


# ═══════════════════════════════════════════════════════════════════════
# Group 2 — Cache invariance (tests #4-5)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l0
@pytest.mark.filterwarnings(
    "ignore::orpheus.numerics.convergence.ConvergenceWarning"
)
def test_collision_cache_invariance_under_source_iteration() -> None:
    """Test #4 (CARDINAL) — :meth:`CollisionCache.from_geometry` called ONCE per σ_t epoch.

    Run a multi-iteration Picard fixed-point on a **heterogeneous**
    fuel|moderator|fuel 2-group slab eigenvalue problem.  The inner SI
    updates ``Q`` each iteration; the outer power iteration updates the
    fission source.  Across both loops, σ_t is bound once at
    ``SNSolver.__init__`` and the cache MUST be built exactly once.

    This is the cardinal invariance gate: it proves the cache placement
    on :class:`SNSolver` survives Picard, and that no sweep path is
    secretly rebuilding the cache on every iteration.

    A heterogeneous 2G case is used deliberately: a homogeneous-reflective
    slab is a degenerate eigenvalue problem (flat flux, k = νΣ_f/Σ_a is
    shape-independent), so the fuel|moderator|fuel slab is what gives the
    flux a genuine spatial + spectral shape for the Picard loop to settle.

    ⭐ **The truncation is the FIXTURE, and it is declared** (#340 R2).
    ``max_inner=50`` at ``inner_tol=1e-8`` is far short of what this slab
    costs — `[M]` 2026-08-10 all 7 inners hit the cap at ρ ≈ 0.958 — and
    **the >= 5 outer floor below depends on that starvation.**  A starved
    inner suppresses the outer increments (the #340 F2 mechanism), which
    is what stretches this solve to 7 outer steps; `[M]` at
    ``max_inner=425`` the SAME heterogeneous slab converges in **3**
    outers and the floor assertion FAILS.  The cardinal invariant
    (``_build_count == 1``) is unaffected either way — σ_t is bound once
    at ``SNSolver.__init__`` — so the truncation costs the gate nothing
    and buys it the long Picard history it asserts against.

    ⛔ Until 2026-08-10 this paragraph credited the ~7 outer steps to the
    fixture's heterogeneity.  They do not come from it; they come from the
    starved inner.  The distinction says which knob is load-bearing: raise
    ``max_inner`` and this test goes red on its own non-degeneracy floor,
    not on its cardinal assertion.

    ⚠ **The floor may be measuring the wrong quantity — filed as #351.**
    The cardinal claim is that no sweep path rebuilds the cache, which is
    a *sweep*-count claim, not an outer-count one.  `[M]` at
    ``max_inner=425`` this solve performs **972** sweeps against **350**
    today, so the invariant is exercised ~2.8× HARDER at the budget that
    breaks the floor.  Re-founding the non-degeneracy floor on sweeps (or
    on the cache call count) would let the budget rise while the gate gets
    stronger — but that changes the assertion, not just the budget, so it
    is deliberately out of scope for an R2 declaration pass.
    """
    from orpheus.derivations.common.xs_library import get_mixture
    from orpheus.geometry import Region, RegionMesh, StructuredGeometry
    from orpheus.sn.solver import solve_sn

    fuel = get_mixture("A", "2g")
    mod = get_mixture("B", "2g")
    materials = {0: fuel, 1: mod}
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(
            Region(mat_id=0, outer_thickness_cm=0.5),
            Region(mat_id=1, outer_thickness_cm=1.0),
            Region(mat_id=0, outer_thickness_cm=0.5),
        ),
        bcs=(BC("reflective"), BC("reflective")),
    )
    mesh = Mesh1D.from_geometry(
        geom,
        region_meshes=(
            RegionMesh(n_cells=2),
            RegionMesh(n_cells=4),
            RegionMesh(n_cells=2),
        ),
    )
    quad = Quadrature.gauss_legendre(4)

    # Reset counter, then run a converged eigenvalue (~7 outer × N inner).
    CollisionCache.reset_build_count()
    result = solve_sn(
        materials=materials,
        mesh=mesh,
        quadrature=quad,
        inner_solver="source_iteration",
        max_outer=50,
        max_inner=50,
        inner_tol=1e-8,
        keff_tol=1e-6,
        flux_tol=1e-5,
    )
    assert len(result.keff_history) >= 5, (
        "Test fixture is too trivial — converged in fewer than 5 outer "
        "iterations; the heterogeneous slab should need ~7. Raise "
        "max_outer or tighten tolerances to exercise the Picard loop."
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
    """Test #5 — :class:`StreamingCoefficientCache` survives ``rebind_cross_sections``.

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
        chi=np.zeros(1),  # non-fissile ⇒ null spectrum (S10a __post_init__ guard)
    )
    materials = {0: mix}
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 6),
        mat_ids=np.zeros(5, dtype=int),
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(4)
    sn_mesh = SNMesh(mesh, quad, materials)
    solver = SNSolver(sn_mesh=sn_mesh)

    geom_before = solver.geom_cache
    coll_before = solver.coll_cache
    new_sig_t = solver.mat_xs.total_cross_section * 2.0
    solver.rebind_cross_sections(new_sig_t)
    geom_after = solver.geom_cache
    coll_after = solver.coll_cache

    assert geom_after is geom_before, (
        "StreamingCoefficientCache should be invariant under σ_t rebinds; "
        "rebind_cross_sections accidentally invalidated Stratum 1."
    )
    assert coll_after is not coll_before, (
        "CollisionCache should be rebuilt on σ_t rebind."
    )
    # P4.9b step 2c — the re-posed halves of the two-stratum contract:
    # (1) Stratum 1 now lives in the strategy layer's INTERN; the mesh-attr
    #     memo is RETIRED (the walk resolves through geometry_cache_for).
    from orpheus.sn.loss_representation import geometry_cache_for

    assert not hasattr(sn_mesh, "_geom_cache"), (
        "the mesh-attr _geom_cache memo was retired at P4.9b step 2c — "
        "the strategy layer's intern is the one home (Q1 ruling)"
    )
    assert geometry_cache_for(
        sn_mesh, sn_mesh.angular_closure,
    ) is geom_before, "the solver's Stratum 1 IS the interned instance"
    # (2) STALENESS — the σ stratum the WALK sees is the fresh one: the
    #     rebind re-stamped the mesh memo, and its values carry the NEW σ
    #     (pre-carve this was structural; post-carve it is asserted).
    #     _coll_cache / _pole_mirror_cache deliberately SURVIVE as mesh
    #     memos — the σ-stratum posing is Campaign 2's consumer-side
    #     territory (the design memo §9 records the ruling).
    assert sn_mesh._coll_cache is coll_after
    assert not np.array_equal(
        coll_after.inverse_denom, coll_before.inverse_denom,
    ), "the rebuilt σ stratum must carry the NEW σ_t (the staleness leg)"


# ═══════════════════════════════════════════════════════════════════════
# Group 3 — Dual-view consistency (tests #6-7)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l0
@pytest.mark.parametrize("geometry", ["slab", "sphere"])
@pytest.mark.parametrize("ng", [1, 2, 3])
@pytest.mark.parametrize("source_kind", ["uniform", "linear", "gaussian"])
def test_cache_driven_sweep_matches_per_cell_scheme_update(
    geometry: str, ng: int, source_kind: str,
) -> None:
    """Test #6 — cache-driven sweep equals per-cell ``scheme.update``.

    For each (geometry × ng × source) combination, run ONE ordinate's
    spatial scan through:

    1. The cache-driven path: ``b = 2·QV·inverse_denom`` (slab) /
       ``b = 2·(QV + ang_contrib)·inverse_denom`` (curvilinear), then
       ``psi_face = ordinate_scan(a, b, psi_in)``;
    2. The per-cell reference path: ``update(visit, ...)`` cell-by-cell.

    The two MUST agree at ``rtol=1e-13`` (Pattern 2 dual-view contract).
    Disagreement signals either a cache-populator algebra bug or a
    non-affine term in the per-cell balance.  Per plan decision-point
    checkpoints: STOP and trace algebraically; do NOT widen rtol.
    """
    nx = 8
    N = 4
    if geometry == "slab":
        sn_mesh = _make_slab(nx=nx, N=N)
    else:
        sn_mesh = _make_sphere(nx=nx, N=N)

    geom = StreamingCoefficientCache.from_mesh_and_quad(sn_mesh, sn_mesh.angular_closure)
    rng = np.random.default_rng(42)
    # Issue #196 PR-INDEX-2: cache consumes σ_t as (ng, nx).
    sig_t = 1.0 + 0.5 * rng.random((ng, nx))                  # (ng, nx)
    coll = CollisionCache.from_geometry(geom, sig_t, sn_mesh.scheme)

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
    # P4.9a: the caller assembles the closure contributions from the
    # closure's own per-ordinate constants (the walk's production idiom).
    closure = sn_mesh.angular_closure
    c_in_n = closure.c_in_per_ordinate[n]
    c_out_n = closure.c_out_per_ordinate[n]
    for k_chain in range(nx):
        cell_i = int(chain[k_chain])
        visit = visits_full[k_chain]
        st_v = visit.streaming_terms
        psi_a_in_cell = (
            psi_a_in_chain[:, k_chain] if psi_a_in_chain is not None else None
        )
        result = dd.update(
            visit=visit,
            total_xs=sig_t[:, cell_i],                         # (ng,) — group axis at axis 0
            source=QV_chain[:, k_chain],                       # (ng,)
            upstream_state=UpstreamState(spatial_upstream=psi_face_in),
            angular_denom_term=st_v.delta_A_over_w * c_out_n,
            angular_numer_upstream=(
                None if psi_a_in_cell is None
                else st_v.delta_A_over_w * c_in_n * psi_a_in_cell
            ),
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
def test_cache_populator_matches_cell_balance_for_streaming() -> None:
    """Test #7 — cache ``(a, 1/denom)`` agrees with the per-cell balance helper.

    Pattern 2 anchor at ``rtol=1e-14``.  For any cell, the cache's per-cell
    ``(a, inverse_denom)`` MUST equal what :func:`cell_balance_for_streaming`
    (at ``n_mask=1``, angular contributions assembled from the visit's
    closure data) would produce — [M] ``affine_scan_coefficients`` computes
    its ``denom`` from its OWN expression and calls neither helper, so
    cache-populator vs balance-helper remain two independent
    implementations and this gate keeps its cross-implementation claim
    (P4.9a rewire of ``…_matches_cell_balance_terms``; the scalar twin
    retired).
    """
    sn_mesh = _make_sphere(nx=8, N=4)
    geom = StreamingCoefficientCache.from_mesh_and_quad(sn_mesh, sn_mesh.angular_closure)
    ng = 2
    # Issue #196 PR-INDEX-2: cache consumes σ_t as (ng, nx).
    # Build (nx, ng) first via outer product for readability, then transpose.
    sig_t_xg = np.linspace(0.5, 1.5, 8)[:, None] * np.array([[1.0, 2.0]])  # (8, 2)
    sig_t = sig_t_xg.T                                                     # (ng=2, nx=8)
    coll = CollisionCache.from_geometry(geom, sig_t, sn_mesh.scheme)

    # Sample two ordinates × two cells (chain positions).
    quad = sn_mesh.quad
    for n in (int(np.argmin(quad.mu_x)), int(np.argmax(quad.mu_x))):
        visits = list(sn_mesh.dag_walk(ordinate_idx=n))
        chain = geom.chain_idx[n]
        for k_chain in (0, 3, 7):
            visit = visits[k_chain]
            cell_i = int(chain[k_chain])
            # Zero upstream + zero angular-numer probe — ``denom`` is
            # independent of upstream, so this is fine.
            st = visit.streaming_terms
            denom_ref, _numer_ref = cell_balance_for_streaming(
                abs_mu=np.array([st.abs_mu]),
                A_downstream=np.array([visit.face_area_downstream]),
                A_total=np.array(
                    [st.face_area_inner + st.face_area_outer],
                ),
                total_xs=sig_t[:, cell_i],                                 # (ng,)
                volume=st.volume,
                psi_face_in=np.zeros((ng, 1)),
                # Angular contribution assembled from the closure's own
                # per-ordinate constant (P4.9a — the visit is purely
                # spatial; the closure owns the c-map).
                angular_denom_term=np.array(
                    [st.delta_A_over_w
                     * sn_mesh.angular_closure.c_out_per_ordinate[n]],
                ),
                angular_numer_upstream=np.zeros((ng, 1)),
            )
            denom_ref = denom_ref[:, 0]                                   # (ng,)
            # Cache layout (N, ng, nx) — fix n, fix k_chain ⇒ (ng,) vector.
            # Indexing pattern updated from [n, k_chain] (legacy axis 1 was
            # cell) to [n, :, k_chain] (PR-INDEX-2 axis 1 is group, axis 2
            # is cell).  The semantic intent is the same: "per-cell denom
            # vector across groups."
            denom_cached = 1.0 / coll.inverse_denom[n, :, k_chain]         # (ng,)
            a_cached = coll.a_attenuation[n, :, k_chain]                   # (ng,)
            # Algebraic a from the reference denom: a = 2|μ|·A_total/denom − 1.
            A_total = (
                st.face_area_inner + st.face_area_outer
            )
            a_expected = (
                2.0 * st.abs_mu * A_total / denom_ref - 1.0
            )
            np.testing.assert_allclose(
                denom_cached, denom_ref, rtol=1e-14,
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
    from orpheus.transport.source_sinks import AngularSourceSink
    from tests.sn._test_helpers import sweep_once

    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 161),
        mat_ids=np.zeros(160, dtype=int),
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(16)
    sn_mesh = SNMesh(mesh, quad, _trivial_materials(ng=4))
    # Issue #196 PR-INDEX-5: Q principled.
    # R-1 Step 4 A1: single per-ordinate source carrier.
    Q = AngularSourceSink.from_isotropic(np.ones((4, *sn_mesh.spatial_shape)), sn_mesh)
    sig_t = np.ones((4, *sn_mesh.spatial_shape))  # (ng, *spatial)
    # Issue #197 PR-TYPED-2: typed boundary state replaces dict.
    boundary_flux = AngularBoundaryFlux.zeros(sn_mesh.angular_trace)

    # Warm-up — first call also caches inside SNMesh.
    for _ in range(3):
        sweep_once(Q, sig_t, sn_mesh, boundary_flux)

    # Measured wall clock over 100 sweeps.
    n_iters = 100
    t0 = time.perf_counter()
    for _ in range(n_iters):
        sweep_once(Q, sig_t, sn_mesh, boundary_flux)
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
        chi=np.zeros(1),  # non-fissile ⇒ null spectrum (S10a __post_init__ guard)
    )
    materials = {0: mix}
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 11),
        mat_ids=np.zeros(10, dtype=int),
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    quad = Quadrature.gauss_legendre(4)
    sn_mesh = SNMesh(mesh, quad, materials)
    # R-1 Step 4 A1 — ``external_source`` is per-ordinate density
    # (already ``/sum_w``).  Iso scalar magnitude 1 ⇒ per-ord ``1/sum_w``.
    sum_w = float(quad.weights.sum())
    Q_external = np.full((quad.N, 1, *sn_mesh.spatial_shape), 1.0 / sum_w)

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
    # rank-d: scalar_flux principled (ng=1, nx=10) — radial slice at g=0.
    np.testing.assert_allclose(result.scalar_flux.values[0, :], expected, rtol=1e-10)


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


@pytest.mark.foundation
def test_closure_algebra_fields_are_the_closures_minted_values() -> None:
    """The cache's ``(N,)`` closure-algebra block IS the closure's mint.

    P4.9a handing gate (``scratch/p4_9a_verification_plan.md`` §5): the
    closure owns the derivation of the scan-normal march constants; the
    cache stores them.  ``array_equal`` — NO tolerance: [M] the realistic
    defect is the algebraically-equal respelling ``tau_inv − 1.0`` inside
    the closure, which sits 1–2 ULP away, so any tolerance ≥ 1e-15 makes
    this gate a non-catcher (the M7 mutation arm is its proof).

    ⚠ What this gate CANNOT see, stated per coding-standards: after the
    handing the right-hand sides are the closure's own accessors, so this
    is a *storage-fidelity* pin (cache == mint) and a *spelling* pin (via
    the closure-side discrimination leg in
    ``test_angular_closure::TestMintedScanConstants``), NOT an
    independent-value pin of the constants themselves — that anchor is
    the closure-vs-surrogate contract in ``test_closure_constant_map.py``.
    """
    sn_mesh = _make_sphere(nx=8, N=4)
    geom = StreamingCoefficientCache.from_mesh_and_quad(sn_mesh, sn_mesh.angular_closure)
    closure = sn_mesh.angular_closure
    tau = closure.tau_per_ordinate
    np.testing.assert_array_equal(geom.tau_inv, 1.0 / tau)
    np.testing.assert_array_equal(geom.mm_a_in_coeff, (1.0 - tau) / tau)
    np.testing.assert_array_equal(geom.c_in, closure.c_in_per_ordinate)
    np.testing.assert_array_equal(geom.c_out, closure.c_out_per_ordinate)


@pytest.mark.foundation
def test_closure_algebra_fields_slab_neutral_element() -> None:
    """Negative leg: the slab cache carries the exact neutral element.

    ``tau_inv == 1``, ``mm_a_in_coeff == 0``, ``c_in == c_out == 0`` —
    bit-exact.  Without a structurally-different input the gate above
    has no reading that could fail for a wrong-geometry reason.
    """
    sn_mesh = _make_slab(nx=8, N=4)
    geom = StreamingCoefficientCache.from_mesh_and_quad(sn_mesh, sn_mesh.angular_closure)
    N = sn_mesh.quad.N
    np.testing.assert_array_equal(geom.tau_inv, np.ones(N))
    np.testing.assert_array_equal(geom.mm_a_in_coeff, np.zeros(N))
    np.testing.assert_array_equal(geom.c_in, np.zeros(N))
    np.testing.assert_array_equal(geom.c_out, np.zeros(N))


def test_geometry_cache_builds_exactly_once_per_mesh() -> None:
    """[foundation] The COUNT gate — Stratum 1 is built ONCE per (mesh, closure).

    P4.9b step 2c (verification plan F2/§2.3(c-ii)): [M] the operator is
    constructed 6-10 times per solve, so a per-operator memo would build
    the table 6-10× where today's count is 1 — up to 24.65 % of a slab
    solve (GL16/nx=200).  The ruled home (Q1) is the strategy layer's
    hub-interned lazy resolve, whose lifetime spans the MESH: the count
    is 1 across a whole solve AND across a second solve on the same
    mesh.  A wall-clock assertion is forbidden (flaky proxy); the count
    is exact and is the only instrument that can see a memo-scoping
    regression.
    """
    import orpheus.sn.loss_representation  # ensure the module is bound
    from orpheus.sn.solver import solve_sn_fixed_source
    from orpheus.sn.sweep.cache import StreamingCoefficientCache
    from tests.sn._test_helpers import placeholder_materials

    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 5),
        mat_ids=np.zeros(4, dtype=int),
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(4)
    materials = placeholder_materials(ng=2)
    sn_mesh = SNMesh(mesh, quad, materials)

    counts = {"builds": 0}
    real = StreamingCoefficientCache.from_mesh_and_quad.__func__

    def counting(cls, m, closure):
        counts["builds"] += 1
        return real(cls, m, closure)

    q_ext = np.ones((quad.N, sn_mesh.ng, 4))
    try:
        StreamingCoefficientCache.from_mesh_and_quad = classmethod(counting)
        # Leg 1 — one FULL solve builds exactly once (the entry builds its
        # own hub from the raw geometry; [M] 6-10 operators live inside).
        solve_sn_fixed_source(materials, mesh, quad, q_ext)
        assert counts["builds"] == 1, (
            f"Stratum 1 built {counts['builds']}x in one solve — the "
            "interned lazy resolve must build exactly once (F2: a "
            "per-operator memo costs up to 24.65 % of a slab solve)"
        )
        # Leg 2 — the intern's LIFETIME is the hub's: two independently
        # posed operators over ONE hub share one build.
        counts["builds"] = 0
        from orpheus.sn.operators.streaming import StreamingOperator
        from orpheus.transport.operators.multiplication_operator import (
            MultiplicationOperator,
        )
        from orpheus.transport.fields.angular_flux import AngularFlux
        from orpheus.transport.fields.angular_boundary_flux import (
            AngularBoundaryFlux,
        )
        from orpheus.transport.timed_full_field import TimedFullField

        sig = np.ones((sn_mesh.ng, *sn_mesh.spatial_shape))
        rhs = TimedFullField.zeros(
            interior=AngularFlux,
            boundary=AngularBoundaryFlux,
            space=sn_mesh.full_field_space,
        )
        rhs.interior.values[...] = 1.0
        for _ in range(2):
            L = StreamingOperator.pose(sn_mesh)
            (L + MultiplicationOperator.from_mesh(sig, sn_mesh)).solve(rhs)
        assert counts["builds"] == 1, (
            f"two posed operators over ONE hub built Stratum 1 "
            f"{counts['builds']}x — the intern's lifetime is the hub's, "
            "not the operator's"
        )
    finally:
        StreamingCoefficientCache.from_mesh_and_quad = classmethod(real)


def test_cache_builder_refuses_a_meshless_chain_under_dash_O() -> None:
    """[foundation] The admission contract raises — not a stripped assert.

    P4.9b step 2c ride-along: the guard was a bare ``assert`` (a NO-OP
    under the canonical ``python -O`` runner) with [M] ZERO witnesses
    tree-wide — this test is NET-NEW coverage, not a migration.  A 2-D
    Cartesian mesh (``reduced is None``) must be REFUSED with the typed
    message, in optimized and debug mode alike.
    """
    from orpheus.sn.sweep.cache import StreamingCoefficientCache

    from orpheus.geometry import Mesh2D
    from tests.sn._test_helpers import placeholder_materials

    mesh2d = Mesh2D(
        edges_x=np.linspace(0.0, 1.0, 4),
        edges_y=np.linspace(0.0, 1.0, 4),
        mat_map=np.zeros((3, 3), dtype=int),
    )
    sn2d = SNMesh(
        mesh2d, Quadrature.level_symmetric(sn_order=4), placeholder_materials(ng=2),
    )
    assert sn2d.reduced is None  # the witness's own premise
    with pytest.raises(TypeError, match="requires a ReducedStreamingOperator"):
        StreamingCoefficientCache.from_mesh_and_quad(
            sn2d, sn2d.angular_closure,
        )
