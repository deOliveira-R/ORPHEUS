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
* Roadmap ``.claude/plans/archive/stencil_assembly_dsa_roadmap.md`` §C2.4b (R10);
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
from orpheus.transport.full_field import FullField
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.transport.radial_characteristic_field import (
    RadialCharacteristicField,
)
from orpheus.transport.source_sinks import (
    AngularBoundarySourceSink,
    AngularSourceSink,
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
        Quadrature.folded_product(n_mu=4, n_phi=8)
        if coord is CoordSystem.CYLINDRICAL
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
    leaf is ``None`` on the slab — since Q5.6.3 the folded cylinder
    CARRIES, so slab is the only seedless geometry here — B.2d: the legs
    are EXPLICIT kwargs)."""
    N, nx = sn.quad.N, sn.nx
    bvals = rng.standard_normal((N, ng, nx))
    # source_from_angular returns None on a non-carrying mesh (slab).
    q_half = RadialCharacteristicField.source_from_angular(bvals, sn)
    b = FullField(
        interior=AngularSourceSink(values=bvals, space=sn.angular_bulk_space),
        boundary=AngularBoundarySourceSink.zeros(sn.angular_trace),
    )
    return bvals, b, q_half


def _random_iterate(sn, rng, ng: int = 2):
    """A random 2-block flux composite (an ``initial_guess``)."""
    N, nx = sn.quad.N, sn.nx
    return FullField(
        interior=AngularFlux(values=rng.standard_normal((N, ng, nx)), space=sn.angular_bulk_space),
        boundary=AngularBoundaryFlux.zeros(sn.angular_trace),
    )


def _joint_solve(A, sn, b, q_half, *, initial_guess=None):
    """The joint solve through THE GRID (step 6 — presence structural);
    returns ``(sol, flux)`` where ``flux`` is the marched ψ½ member
    (``None`` seedless)."""
    if q_half is None:
        return A.solve(b, initial_guess=initial_guess), None
    from orpheus.numerics.coupled_system import CoupledField

    from tests.sn._test_helpers import joint_m_grid

    grid, _space = joint_m_grid(sn, A)
    state = grid.solve(CoupledField(systems=(b, q_half)))
    return state.systems[0], state.systems[1]


def _joint_apply(A, sn, psi, flux):
    """The joint matvec through THE GRID (step 6); returns ``(out, rows)``
    where ``rows`` is the emitted ray-block member (``None`` seedless)."""
    if flux is None:
        return A.apply(psi), None
    from orpheus.numerics.coupled_system import CoupledField

    from tests.sn._test_helpers import joint_m_grid

    grid, _space = joint_m_grid(sn, A)
    out = grid.apply(CoupledField(systems=(psi, flux)))
    return out.systems[0], out.systems[1]


def _full_residual_inf(A, sn, sol, flux, q_half, bvals) -> float:
    """‖A·solve(b) − b‖_∞ / ‖b‖_∞ over the FULL joint field (bulk + trace
    + the emitted ray rows vs the q½ source)."""
    r, rows = _joint_apply(A, sn, sol, flux)
    num = np.max(np.abs(r.interior.values - bvals))
    if rows is not None and q_half is not None:
        num = max(num, np.max(np.abs(rows.to_flat() - q_half.to_flat())))
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
    exact inverse.  At the time slab and cyl were seedless-exact (the
    LS4 cylinder was non-carrying) and had to STAY; since Q5.6.3 the
    folded cylinder rides the COUPLED branch like the sphere, so its
    row now exercises the joint grid — including the ψ½ outflow row
    ERR-078 repaired."""
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
        f1.to_flat(), f2.to_flat(),
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
    if sn.radial_characteristic_field_space is not None:
        flux0 = RadialCharacteristicField.flux_zeros(sn.radial_characteristic_field_space)
        for level in sn.radial_characteristic_levels:
            for sign in (-1, +1):
                cells = flux0.interior.cells(level, sign)
                cells[...] = rng.standard_normal(cells.shape)
                corner = flux0.boundary.corner(level, sign)
                corner[...] = rng.standard_normal(corner.shape)
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


# The two Mode-11 routing spies (A_BB.solve / .solve_transpose is on the
# joint solve's call graph) RETIRED at step 6: the joint solve IS the
# grid's block substitution, whose (B, B) leg calls A_BB by
# CONSTRUCTION — the routed-around fused channel the spies guarded is
# deleted, so the claim is a type fact, not an execution fact.


@pytest.mark.foundation
def test_4e_unweave_walk_source_has_no_carlson_reference():
    r"""S2 — the un-weave completion tripwire: the fused ``(L+C)`` walk
    (:mod:`orpheus.sn.loss_representation`) carries ZERO textual references
    to the ``carlson_inward_sweep_*`` engines — System B's marches live ONLY
    behind ``RadialCharacteristicOperator.solve`` / ``solve_transpose``. An
    inline march creeping back (an import, a call, even a stale comment
    naming the engine) REDS here before any value gate could go blind."""
    from pathlib import Path

    import orpheus.sn.loss_representation as lr

    src = Path(lr.__file__).read_text(encoding="utf-8")
    n = src.count("carlson")
    if n != 0:
        pytest.fail(
            f"loss_representation source contains {n} 'carlson' "
            "reference(s) — the 4e un-weave is regressing (the walk must "
            "route System B through RadialCharacteristicOperator, never "
            "inline the engine)."
        )


def _coupled_space(sn):
    """The mesh's coupled ψ½ space (System A ⊕ System B, member-wise metric)."""
    from orpheus.numerics.coupled_system import CoupledSpace

    return CoupledSpace.from_systems(
        (sn.full_field_space, sn.radial_characteristic_field_space),
    )


@pytest.mark.foundation
def test_mode12_g_reciprocity_catches_a_seed_row_flip():
    r"""Mode 12 (a) — CLOSED (INVERTED).  With the state metric ``G_sd =
    V_cell`` (SPD) and a NONZERO ψ½ seed, a seed-row (``A_ss``) sign flip
    makes the metric-weighted G-reciprocity RED — the seed rows now carry
    metric weight, so they lie OUTSIDE the functional's invariance group.
    Pre-fix the ghost ``G_sd = 0`` put them INSIDE it and the flip was
    invisible (the standing false-green this gate replaces).

    Mechanism (diag_gsd_02/03, re-aimed at 5d onto the NAMED block): the
    mutation flips ``RadialCharacteristicSeeding.apply`` (the grid M's
    forward seed placement) but NOT its transpose (``A.H``'s
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
    from orpheus.sn.loss_representation import _OneDimScanWalk
    from orpheus.transport.radial_characteristic_field import (
        RadialCharacteristicField,
    )

    sn, A = _operator(CoordSystem.SPHERICAL, nx=4, sigma=0.5)
    rng = np.random.default_rng(12)
    n_seed = sn.radial_characteristic_field_space.shape[0]
    n_trace = int(sn.angular_trace.layout.total_size)
    space = _coupled_space(sn)
    # Step 5: the joint operator is the honest triangular M grid on the
    # coupled pair — its .H is the member-wise G-adjoint over System A's
    # (V·w ⊕ |Ω·n|·w) and System B's (V_cell ⊕ corner-gauge) metrics, so
    # the seed rows carry metric weight OUTSIDE the functional's
    # invariance group.
    from tests.sn._test_helpers import joint_m_grid
    M_op, _ = joint_m_grid(sn, A)

    def _random_pair():
        # NONZERO ψ_B — activates the V_cell metric AND the A_BB self-rows.
        psi_a = FullField(
            interior=AngularFlux(values=rng.standard_normal((sn.quad.N, sn.ng, sn.nx)), space=sn.angular_bulk_space),
            boundary=AngularBoundaryFlux(
                values=rng.standard_normal(n_trace),
                space=sn.angular_trace,
            ),
        )
        seed = RadialCharacteristicField.from_flat(
            rng.standard_normal(n_seed),
            RadialCharacteristicField.flux_zeros(sn.radial_characteristic_field_space),
        )
        return CoupledField(systems=(psi_a, seed))

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

    # ── Mode-12 CLOSURE: the seed-placement flip now REDS reciprocity. ──
    # Step 5 → step 6: the grid's forward seed placement IS the named
    # RadialCharacteristicSeeding block (the walk's fused ``_seed_rows_*``
    # went production-dead at 5b and were DELETED at step 6 — the block is
    # the one spelling); flip the FORWARD block only (its transpose stays
    # honest) so the SPD G_sd=V_cell reciprocity sees the inconsistency.
    from orpheus.sn.operators.radial_characteristic import (
        RadialCharacteristicSeeding,
    )

    orig_apply = RadialCharacteristicSeeding.apply

    def _flipped(self, x, /):
        return -orig_apply(self, x)

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(RadialCharacteristicSeeding, "apply", _flipped)
        lhs = float(space.inner_product(M_op.apply(psi), phi))
        rhs = float(space.inner_product(psi, M_op.H.apply(phi)))
    recip_rel = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)
    if not recip_rel > 1e-6:
        pytest.fail(
            f"G-reciprocity did NOT move under the forward seed-placement "
            f"flip (rel={recip_rel:.2e}) — Mode-12 is NOT closed; the seed "
            "coupling sits in the functional's invariance group (is G_sd "
            "the ghost 0, or is the mutation off M's call graph?)."
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
        A, sn, b, RadialCharacteristicField.source_zeros(sn.radial_characteristic_field_space),
    )
    delta = np.max(np.abs(sol_with.interior.values - sol_without.interior.values))
    if delta <= 1e-12:
        pytest.fail(
            "the q½ source block is DEAD: zeroing it left the sphere solve "
            f"unchanged (Δ={delta:.3e}) — the seed carrier is inert."
        )
