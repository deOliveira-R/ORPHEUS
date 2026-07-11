r"""§16.C — the #282 route-(a) FIXED-POINT acceptance gates (2.5d).

The decisive classifiers for the #282 direct starting-direction fix
(#280 Phase 2.5d, ruling R10).  Route (a) retired the lagged
Morel–Montry ψ½ pole seed (extrapolate-from-the-iterate — a walk-order
BACK EDGE that made the spherical SOLVE a non-direct inverse) and
replaced it with a first-class ψ½ STATE block whose legs the sweep
marches DIRECTLY from the true q½ source.  The keystone acceptance
number is the sphere **cold residual 5.18e5 → < 1e-11** — the single
measurement that certifies the lag is dead (solve ≡ apply on the cold
start, a genuine single-pass exact inverse).

Promoted from the diagnostics ``diag_curvilinear_seed_sensitivity`` /
``diag_sphere_fixedpoint_consistency`` probes.  Every gate measures the
FULL augmented field (bulk ⊕ trace ⊕ seed) — a bulk-only norm would be
seed-blind (vv-principles Mode 12 (b)).  The seed COEFFICIENTS are
certified by §16.B (``carlson_inward_sweep_from_source`` convergence)
and the 2.5b Euclidean Mᵀ oracle; THIS file certifies the lag death
(solve↔apply consistency + seed-insensitivity + end-to-end physicality).

References
==========
* Roadmap ``.claude/plans/stencil_assembly_dsa_roadmap.md`` §C2.4b (R10);
  gate spec ``a3_solve_transpose_verification.md`` §16.C.
* Issue #282 (the spherical starting-direction seed lag).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.sn.solver import solve_sn_fixed_source
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.radial_characteristic_flux import RadialCharacteristicFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.transport.source_sinks import (
    AngularBoundarySourceSink,
    AngularSourceSink,
    RadialCharacteristicSourceSink,
)

pytestmark = [pytest.mark.regression]


# ── builders ──────────────────────────────────────────────────────────


def _mixture(sig_t: float, sig_s: float, ng: int = 2) -> Mixture:
    """A homogeneous ≥2G mixture with a chosen scattering ratio ``c``.

    Per-group σ_t varies (non-degenerate group axis); the P0 in-group
    scatter carries ``c = sig_s/sig_t`` (0 ⇒ pure absorber, C(iv))."""
    st = np.array([sig_t * (1.0 + 0.4 * g) for g in range(ng)])
    ss = np.diag([sig_s * (1.0 + 0.4 * g) for g in range(ng)])   # (ng, ng) P0
    sig_c = st - ss.sum(axis=0)                                   # absorption
    return make_mixture(
        sig_t=st, sig_c=sig_c, sig_f=np.zeros(ng), nu=np.zeros(ng),
        chi=np.zeros(ng), sig_s=ss,
    )


def _operator(coord: CoordSystem, nx: int, *, sigma: float, ng: int = 2):
    """``A = L + C`` on a homogeneous curvilinear/slab mesh + its SNMesh."""
    kw = (
        dict(bc_left=BC("vacuum"), bc_right=BC("vacuum"))
        if coord is CoordSystem.CARTESIAN else dict(bc_right=BC("vacuum"))
    )
    mesh = Mesh1D(
        edges=np.linspace(0.0, 4.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int), coord=coord, **kw,
    )
    quad = (
        Quadrature.level_symmetric(4) if coord is CoordSystem.CYLINDRICAL
        else Quadrature.gauss_legendre(4)
    )
    sn = SNMesh(mesh, quad, {0: _mixture(sigma, 0.4 * sigma, ng=ng)})
    sig_t = np.stack(
        [np.full(sn.spatial_shape, sigma * (1.0 + 0.3 * g)) for g in range(ng)],
        axis=0,
    )
    return sn, StreamingOperator(sn) + MultiplicationOperator.from_mesh(sig_t, sn)


def _random_source(sn, rng, ng: int = 2):
    r"""A 2-block volumetric source with ZEROED inflow trace + the q½ leaf
    folded from the bulk (the joint rhs legs on a carrying mesh; the q½
    leaf is ``None`` on slab/cyl — B.2d: the legs are EXPLICIT kwargs)."""
    N, nx = sn.quad.N, sn.nx
    bvals = rng.standard_normal((N, ng, nx))
    q_half = (
        RadialCharacteristicSourceSink.from_angular_source(bvals, sn)
        if sn.radial_characteristic_space is not None else None
    )
    b = FullField(
        interior=AngularSourceSink.from_mesh(bvals, sn),
        boundary=AngularBoundarySourceSink.zeros_on(sn),
    )
    return bvals, b, q_half


def _random_iterate(sn, rng, ng: int = 2):
    """A random 2-block flux composite (an ``initial_guess``)."""
    N, nx = sn.quad.N, sn.nx
    return FullField(
        interior=AngularFlux.from_mesh(rng.standard_normal((N, ng, nx)), sn),
        boundary=AngularBoundaryFlux.zeros_on(sn),
    )


def _joint_solve(A, sn, b, q_half, *, initial_guess=None):
    """``A.solve`` through the B.2d explicit legs; returns ``(sol, flux)``
    where ``flux`` is the marched ψ½ carrier (``None`` seedless)."""
    if q_half is None:
        return A.solve(b, initial_guess=initial_guess), None
    buf = RadialCharacteristicFlux.zeros_on(sn)
    sol = A.solve(
        b, initial_guess=initial_guess,
        radial_characteristic_source=q_half,
        radial_characteristic_flux=buf,
    )
    return sol, buf


def _joint_apply(A, sn, psi, flux):
    """``A.apply`` through the B.2d explicit legs; returns ``(out, rows)``
    where ``rows`` is the emitted ray-block buffer (``None`` seedless)."""
    if flux is None:
        return A.apply(psi), None
    rows = RadialCharacteristicSourceSink.zeros_on(sn)
    out = A.apply(
        psi,
        radial_characteristic_flux=flux,
        radial_characteristic_source=rows,
    )
    return out, rows


def _full_residual_inf(A, sn, sol, flux, q_half, bvals) -> float:
    """‖A·solve(b) − b‖_∞ / ‖b‖_∞ over the FULL joint field (bulk + trace
    + the emitted ray rows vs the q½ source)."""
    r, rows = _joint_apply(A, sn, sol, flux)
    num = np.max(np.abs(r.interior.values - bvals))
    if rows is not None and q_half is not None:
        num = max(num, np.max(np.abs(rows.values - q_half.values)))
    # boundary source is zeroed; r.boundary is the trace defect (≈0 too).
    num = max(num, np.max(np.abs(r.boundary.values)))
    return float(num / np.max(np.abs(bvals)))


_COORDS = [
    pytest.param(CoordSystem.SPHERICAL, id="sphere"),
    pytest.param(CoordSystem.CARTESIAN, id="slab"),
    pytest.param(CoordSystem.CYLINDRICAL, id="cyl"),
]


# ── C(i) — the cold-residual acceptance (lag death) ────────────────────


@pytest.mark.foundation
@pytest.mark.parametrize("coord", _COORDS)
def test_ci_cold_residual_is_machine_zero(coord):
    r"""C(i) — ``‖A·solve(b) − b‖_∞ < 1e-11`` on a COLD start, over the
    FULL augmented field.  The keystone: pre-fix the sphere sat at
    5.18e5 (the ψ½ seed lag made the cold solve a non-inverse); route
    (a) marches the seed directly, so the cold solve is a single-pass
    exact inverse (slab/cyl were already exact — they must STAY)."""
    sn, A = _operator(coord, nx=4, sigma=0.5)
    rng = np.random.default_rng(282)
    bvals, b, q_half = _random_source(sn, rng)
    sol, flux = _joint_solve(A, sn, b, q_half)
    r = _full_residual_inf(A, sn, sol, flux, q_half, bvals)
    assert r < 1e-11, f"[{coord}] cold residual {r:.3e} ≥ 1e-11 (the seed lag)"


# ── C(ii) — seed-insensitivity + the Probe-6 corroboration ─────────────


@pytest.mark.foundation
def test_cii_sphere_solve_is_seed_insensitive_bitwise():
    r"""C(ii) — the lag SIGNATURE dies: two random ``initial_guess`` seeds
    give a BITWISE-identical sphere solve (pre-fix Δ = 4.57e-2; the
    direct march does not read the guess for the ψ½ rows)."""
    sn, A = _operator(CoordSystem.SPHERICAL, nx=4, sigma=0.5)
    rng = np.random.default_rng(6)
    _bvals, b, q_half = _random_source(sn, rng)
    s1, f1 = _joint_solve(A, sn, b, q_half, initial_guess=_random_iterate(sn, rng))
    s2, f2 = _joint_solve(A, sn, b, q_half, initial_guess=_random_iterate(sn, rng))
    np.testing.assert_array_equal(
        s1.interior.values, s2.interior.values,
        err_msg="sphere solve depends on the initial guess — the seed lag lives",
    )
    np.testing.assert_array_equal(
        f1.values, f2.values,
        err_msg="the marched ψ½ carrier depends on the initial guess",
    )


@pytest.mark.foundation
@pytest.mark.parametrize("coord", _COORDS)
def test_cii_probe6_cold_solve_recovers_preimage(coord):
    r"""C(ii) Probe-6 corroboration: ``ψ₀`` arbitrary → ``b = A·ψ₀`` →
    the COLD ``A.solve(b)`` recovers ``ψ₀`` to ``rtol=1e-11`` (pre-fix
    only the WARM sphere solve recovered it).  This ALSO pins that the
    fix did not move the fixed point away from the correct one."""
    sn, A = _operator(coord, nx=4, sigma=0.5)
    rng = np.random.default_rng(60)
    psi0 = _random_iterate(sn, rng)
    flux0 = None
    if sn.radial_characteristic_space is not None:
        flux0 = RadialCharacteristicFlux.zeros_on(sn)
        flux0.values[...] = rng.standard_normal(flux0.values.shape)
    b, rows = _joint_apply(A, sn, psi0, flux0)
    sol, _flux = _joint_solve(A, sn, b, rows)
    np.testing.assert_allclose(
        sol.interior.values, psi0.interior.values, rtol=1e-11, atol=1e-12,
        err_msg=f"[{coord}] cold solve did not recover the pre-image ψ₀",
    )


# ── C(iii) — end-to-end coarse fixed-source physicality ────────────────


@pytest.mark.l1
def test_ciii_coarse_sphere_fixed_source_finite_positive():
    r"""C(iii) — the #282-comment companion: a coarse 16-cell sphere
    fixed-source solve is FINITE + NON-NEGATIVE on BOTH inner drivers
    (pre-fix: SI → NaN, Krylov → negative flux).  A physicality /
    robustness gate (flux-shape layer), not a precision claim."""
    nx = 16
    mesh = Mesh1D(
        edges=np.linspace(0.0, 8.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int), coord=CoordSystem.SPHERICAL,
        bc_right=BC("vacuum"),
    )
    materials = {0: _mixture(1.0, 0.5, ng=2)}
    source = np.ones((Quadrature.gauss_legendre(8).N, 2, nx))
    for driver in ("source_iteration", "krylov"):
        sol = solve_sn_fixed_source(
            materials, mesh, Quadrature.gauss_legendre(8), source,
            boundary_condition="vacuum", inner_solver=driver,
        )
        flux = np.asarray(sol.scalar_flux.values)
        assert np.all(np.isfinite(flux)), f"[{driver}] non-finite flux"
        assert np.all(flux >= 0.0), f"[{driver}] negative flux (min {flux.min():.3e})"


# ── C(iv) — pure-absorber c=0 (the no-outer-iteration degenerate) ──────


@pytest.mark.foundation
def test_civ_pure_absorber_sphere_cold_solve_exact():
    r"""C(iv) — a c = 0 pure-absorber sphere has NO scattering outer loop,
    so the cold solve IS the answer.  Pre-fix it NaN'd (the seed lag with
    no SI loop to mask it); post-fix the direct solve is a genuine
    single-pass exact inverse (C(i) < 1e-11, finite + positive)."""
    kw = dict(bc_right=BC("vacuum"))
    mesh = Mesh1D(
        edges=np.linspace(0.0, 4.0, 5), mat_ids=np.zeros(4, dtype=int),
        coord=CoordSystem.SPHERICAL, **kw,
    )
    sn = SNMesh(mesh, Quadrature.gauss_legendre(4), {0: _mixture(0.8, 0.0)})
    sig_t = np.stack(
        [np.full(sn.spatial_shape, 0.8 * (1.0 + 0.3 * g)) for g in range(2)],
        axis=0,
    )
    A = StreamingOperator(sn) + MultiplicationOperator.from_mesh(sig_t, sn)
    rng = np.random.default_rng(4)
    bvals, b, q_half = _random_source(sn, rng)
    sol, flux = _joint_solve(A, sn, b, q_half)
    r = _full_residual_inf(A, sn, sol, flux, q_half, bvals)
    assert r < 1e-11, f"pure-absorber cold residual {r:.3e} ≥ 1e-11"
    # A single-pass exact inverse of a physical source is finite.
    assert np.all(np.isfinite(sol.interior.values))


# ═══════════════════════════════════════════════════════════════════════
# §16.F — Mode-11 (rewired-path execution) + Mode-12 (invariant-functional)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_mode11_direct_solver_executes_on_sphere_not_slab():
    r"""Mode 11 — the rewired production line IS on the solve's call graph.

    Wrap-sentinel the direct starting-direction solver
    ``carlson_inward_sweep_from_source`` (the name the 1-D scan calls) and
    confirm the sphere COLD solve executes it (counter > 0) while the slab
    does NOT (counter == 0, no carrying levels).  A green solve that never
    reached the new reader would be a vacuous gate (Mode 11)."""
    import orpheus.sn.loss_representation as lr

    calls = {"n": 0}
    orig = lr.carlson_inward_sweep_from_source

    def _counting(*a, **k):
        calls["n"] += 1
        return orig(*a, **k)

    for coord, expect_calls in (
        (CoordSystem.SPHERICAL, True), (CoordSystem.CARTESIAN, False),
    ):
        sn, A = _operator(coord, nx=4, sigma=0.5)
        _bvals, b, q_half = _random_source(sn, np.random.default_rng(11))
        calls["n"] = 0
        with pytest.MonkeyPatch.context() as mp:
            mp.setattr(lr, "carlson_inward_sweep_from_source", _counting)
            _joint_solve(A, sn, b, q_half)
        if expect_calls and calls["n"] == 0:
            pytest.fail(
                f"[{coord}] the direct ψ½ solver was NEVER called by the "
                "sphere cold solve — the rewired path is off the call graph."
            )
        if not expect_calls and calls["n"] != 0:
            pytest.fail(
                f"[{coord}] the slab solve called the ψ½ direct solver "
                f"{calls['n']}× — it has no carrying levels."
            )


def _coupled_space(sn):
    """The mesh's coupled ψ½ space (System A ⊕ System B, member-wise metric)."""
    from orpheus.numerics.coupled_system import CoupledSpace

    return CoupledSpace.from_systems(
        (sn.full_field_space, sn.radial_characteristic_composite_space),
    )


@pytest.mark.foundation
def test_mode12_g_reciprocity_catches_a_seed_row_flip():
    r"""Mode 12 (a) — CLOSED (INVERTED).  With the state metric ``G_sd =
    V_cell`` (SPD) and a NONZERO ψ½ seed, a seed-row (``A_ss``) sign flip
    makes the metric-weighted G-reciprocity RED — the seed rows now carry
    metric weight, so they lie OUTSIDE the functional's invariance group.
    Pre-fix the ghost ``G_sd = 0`` put them INSIDE it and the flip was
    invisible (the standing false-green this gate replaces).

    Mechanism (diag_gsd_02/03): the mutation flips ``_seed_rows_forward``
    (the forward ``A_ss``) but NOT ``_seed_rows_transpose`` (``A.H``'s
    independently-coded reverse mode), creating an apply-vs-transpose
    inconsistency.  With ``G_sd = V_cell`` the seed contributes
    ``Σ V_cell·(A_ss ψ_seed)·φ_seed`` to ⟨Aψ,φ⟩_G that ``A.H``'s unflipped
    transpose cannot match → RED (diag measured 1.000).

    Two legs, BOTH required (this is the subtlety that keeps the gate
    honest): the CONTROL leg — unmutated nonzero-seed reciprocity HOLDS
    (< 1e-12) — proves the baseline adjoint is the honest V_cell one, so the
    mutated RED is attributable to the FLIP.  Without the control leg a
    reverted ghost ``G_sd = 0`` would ALSO give a mutated defect ~0.107 > tol
    (a broken baseline mimicking "caught") and fool the gate."""
    from orpheus.numerics.coupled_system import CoupledField
    from orpheus.sn.coupled_system import CoupledInvertibleOperator
    from orpheus.sn.loss_representation import _OneDimScanWalk
    from orpheus.transport.radial_characteristic_composite import (
        RadialCharacteristicComposite,
    )

    sn, A = _operator(CoordSystem.SPHERICAL, nx=4, sigma=0.5)
    rng = np.random.default_rng(12)
    n_seed = sn.radial_characteristic_space.shape[0]
    n_trace = int(sn.angular_trace.layout.total_size)
    space = _coupled_space(sn)
    # B.2d: the joint operator is M on the coupled pair (the fused walk
    # behind the explicit legs) — its .H is the member-wise G-adjoint over
    # System A's (V·w ⊕ |Ω·n|·w) and System B's (V_cell ⊕ corner-gauge)
    # metrics, so the seed rows carry metric weight OUTSIDE the functional's
    # invariance group.
    M_op = CoupledInvertibleOperator(A, space=space, sn_mesh=sn)

    def _random_pair():
        # NONZERO ψ_B — activates the V_cell metric AND the A_BB self-rows.
        psi_a = FullField(
            interior=AngularFlux.from_mesh(
                rng.standard_normal((sn.quad.N, sn.ng, sn.nx)), sn,
            ),
            boundary=AngularBoundaryFlux(
                values=rng.standard_normal(n_trace),
                space=sn.angular_trace, mesh=sn,
            ),
        )
        seed = RadialCharacteristicFlux(
            values=rng.standard_normal(n_seed),
            space=sn.radial_characteristic_space, mesh=sn,
        )
        return CoupledField(systems=(
            psi_a, RadialCharacteristicComposite.from_unified(seed),
        ))

    psi, phi = _random_pair(), _random_pair()

    # ── CONTROL leg: UNMUTATED nonzero-seed reciprocity HOLDS (V_cell is the
    #    honest adjoint for seed data — pre-fix the ghost broke this ~1.3e-2). ──
    lhs0 = float(space.inner_product(M_op.apply(psi), phi))
    rhs0 = float(space.inner_product(psi, M_op.H.apply(phi)))
    rel0 = abs(lhs0 - rhs0) / (abs(lhs0) + abs(rhs0) + 1e-300)
    if not rel0 < 1e-12:
        pytest.fail(
            f"UNMUTATED nonzero-seed reciprocity broke (rel={rel0:.2e}) — "
            "G_sd=V_cell is not the honest adjoint for seed data (is the ghost "
            "G_sd=0 back?); the mutated RED below would be un-attributable."
        )

    # ── Mode-12 CLOSURE: the seed-row flip now REDS reciprocity. ──
    orig_rows = _OneDimScanWalk._seed_rows_forward

    def _flipped(self, sigma, seed):
        return -orig_rows(self, sigma, seed)   # flip the forward A_BB rows only

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(_OneDimScanWalk, "_seed_rows_forward", _flipped)
        lhs = float(space.inner_product(M_op.apply(psi), phi))
        rhs = float(space.inner_product(psi, M_op.H.apply(phi)))
    recip_rel = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)
    if not recip_rel > 1e-6:
        pytest.fail(
            f"G-reciprocity did NOT move under a seed-row flip "
            f"(rel={recip_rel:.2e}) — Mode-12 is NOT closed; the seed rows are "
            "still in the functional's invariance group (is G_sd the ghost 0?)."
        )


@pytest.mark.foundation
def test_mode10_seed_source_activation_q_half_moves_the_sphere_solve():
    r"""§16.E teeth #3 / Mode 10 — the ACTIVATION proof: zeroing the q½
    source block MOVES the sphere solve.  If it did NOT, the carrier
    augmentation would be inert (a Mode-11-adjacent vacuity — the whole
    route-(a) machinery un-exercised)."""
    sn, A = _operator(CoordSystem.SPHERICAL, nx=4, sigma=0.5)
    rng = np.random.default_rng(103)
    bvals, b, q_half = _random_source(sn, rng)
    sol_with, _f1 = _joint_solve(A, sn, b, q_half)
    # Zero the q½ leg only (keep bulk + trace); re-solve.
    sol_without, _f2 = _joint_solve(
        A, sn, b, RadialCharacteristicSourceSink.zeros_on(sn),
    )
    delta = np.max(np.abs(sol_with.interior.values - sol_without.interior.values))
    if delta <= 1e-12:
        pytest.fail(
            "the q½ source block is DEAD: zeroing it left the sphere solve "
            f"unchanged (Δ={delta:.3e}) — the seed carrier is inert."
        )
