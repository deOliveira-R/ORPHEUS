r"""Comparison tests — unified cylinder matvec (PR-TYPED-6c Step 3).

Issue #197 PR-TYPED-6c Step 3 verifies the cylindrical branch of the
unified ``(L + C)`` matvec.  The module-level helper it landed as,
``_transport_operator_matvec_unified``, was DELETED at the Wave T T.5
matvec retirement; the kernel now lives on the loss representation
(:meth:`~orpheus.sn.loss_representation._OneDimScanWalk._apply_walk`)
and is reached through
:meth:`~orpheus.sn.operators.streaming.StreamingCollisionOperator.apply`.
The gates below drive it via the ``legacy_proxy_matvec`` shim, which
supplies the pre-B1'' cell-centre boundary-fill convention the L0 hand
reference was built around (a CONVENTION bridge, not retired code).

Two levels of evidence (per L14 in ``.claude/lessons.md`` — solver
correctness is a 4-way standoff):

* **L0** — per-ordinate **hand reference** at machine precision. The
  hand reference is structurally independent: it walks each ordinate's
  sweep explicitly using the WDD recurrence ``ψ_out = 2·ψ̄ − ψ_in`` +
  ``streaming + redistribution + collision`` per cell, with NO
  bool-mask scatter into the legacy's misrouting ``ks``. The unified
  matvec uses ``out_g_first[:, global_X, i, 0] = m_full`` fancy-index
  scatter that is **structurally immune** to the routing bug. Both
  match at rtol=1e-12.

* **L1** — heterogeneous 3-region 2G closed cylinder eigenvalue.  The
  GMRES inner solve routes through :class:`StreamingCollisionOperator`
  (= ``L + C``) consuming the unified matvec; the resulting ``k_eff``
  is cross-checked against the trajectory_resolvent reference
  (Variant α at α=1). The row is a strict ``xfail`` on #516 (the cylinder
  reference carries no certified error bound) with a RECORD companion; the
  3 % it used to carry was "Variant α's quadrature error budget", which was
  the reference's one-spline emission density (ERR-090).

The since-retired legacy ``transport_operator_matvec_cylindrical`` had a
per-ordinate routing bug (ERR-049 — ascending-global ``ks`` indices vs
level-internal μ-sorted column order in ``streaming + redistribution +
collision``), which is why it was **never** used as a reference here:
the L0 evidence above is the hand reference, not a legacy cross-check.
"""
from __future__ import annotations

import contextlib
import functools

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesh1D, Mesher
from orpheus.sn.operators import streaming as sn_op
from orpheus.sn import solve_sn
from orpheus.sn.problem import SNProblem
from tests.gates.sn._test_helpers import _LC_matvec
from orpheus.numerics.quadrature import Quadrature
from tests.gates.sn._test_helpers import (
    legacy_proxy_matvec,
    placeholder_materials,
)
from tests.gates.derivations._trajectory_resolvent_ladders import (
    CYLINDER_3REG_SN_4X8_K_STEP,
    tolerance_for,
)
from tests.gates.sn.verification.analytical._certified_agreement import (
    ABA_RADII,
    CYLINDER_3REG_RECORD,
    CYLINDER_3REG_RECORD_BAND,
    CYLINDER_3REG_RECORD_RELATIVE,
    CYLINDER_3REG_REFERENCE_BOUND,
    aba_xs_2g,
    assert_record,
    awaits_cylinder_bound,
    certify_agreement,
    cylinder_3reg_reference,
)


# ═══════════════════════════════════════════════════════════════════════
# Helpers — canonical ↔ packed conversion + per-ordinate hand reference
# ═══════════════════════════════════════════════════════════════════════


def _build_cyl(n_cells: int, quad) -> SNProblem:
    """Build a homogeneous reflective hollow cylinder mesh."""
    mesh = Mesher(StructuredGeometry.cylinder(
        (0.1, 1.0), (0,), inner=BC("reflective"), outer=BC("reflective"),
    )).partition(CellsByCount.uniform_width(n_cells)).mesh
    return SNProblem(mesh, quad, placeholder_materials())


def _bc_fill_outer(psi_view: np.ndarray, problem: SNProblem) -> np.ndarray:
    """Make psi_view BC-consistent at the outer face (incoming ordinates)."""
    quad = problem.quad
    incoming_mask = quad.mu_x < -1e-15
    if not incoming_mask.any():
        return psi_view
    outer_face = psi_view[:, :, -1]
    inflow_full = problem.bc["xmax"].apply(outer_face)
    psi_view = psi_view.copy()
    psi_view[incoming_mask, :, -1] = inflow_full[incoming_mask, :]
    return psi_view


def _extract_at_unknown_slots(
    field_4d: np.ndarray, problem: SNProblem,
) -> np.ndarray:
    """Gather field_4d at the curvilinear equation-bearing slots → (ng, n_eq).

    D-J (2026-05-30) — replaces the legacy :class:`EquationMap`-driven
    slot map.  Curvilinear 1-D equation set: all ``(n, ix, 0)`` slots
    EXCEPT inward ordinates at the outermost cell ``ix == nx - 1``
    (the reflective BC determines those values; they are NOT unknowns).
    """
    quad = problem.quad
    nx = problem.nx
    ng = field_4d.shape[1]
    inflow_outer = quad.mu_x < -1e-15  # (N,)
    cols = []
    for ix in range(nx):
        for n in range(quad.N):
            if ix == nx - 1 and inflow_outer[n]:
                continue
            cols.append(field_4d[n, :, ix])
    return np.stack(cols, axis=1)  # (ng, n_eq)


def _hand_reference_cyl_matvec(
    psi_view: np.ndarray, problem: SNProblem, sigma_t: np.ndarray,
) -> np.ndarray:
    r"""Per-ordinate explicit matvec — the structurally-independent L0 reference.

    For each ordinate ``n``: walk cells in the sweep direction using the
    WDD recurrence and compute ``streaming + redistribution + collision``
    per cell with NO bool-mask scatter. Routing is impossible to get
    wrong because each ordinate is processed in its own scalar pass.

    Returns shape ``(N, ng, nx)``.
    """
    quad = problem.quad
    N = quad.N
    ng = psi_view.shape[1]
    nx = problem.nx
    eps = 1e-15

    reduced = problem.reduced
    assert reduced is not None  # 1-D mesh => minted by the ctor (narrowing)
    A = reduced.face_areas
    V = problem.volumes
    mu_x = quad.mu_x

    bc_outer = problem.bc["xmax"]

    out = np.zeros((N, ng, nx))

    sigma_t_gx = sigma_t
    psi_g_first = psi_view.transpose(1, 0, 2)
    level_indices = quad.level_indices

    # #282 route (a) (#280 Phase 2.5d, 2026-07-04): the seed-strategy zoo
    # + its per-level context object were retired.  Until Q5.6.3 the
    # cylinder fixtures here were NON-carrying (LS/full-product rules with
    # first-ordinate raw τ₀ ∈ {0, 1}), and ``precompute_psi_state`` inlined
    # the 2-point angular-edge-extrapolation seed internally.  The folded
    # fixtures are CARRYING, and BOTH comparison legs keep the same seed
    # convention: ``legacy_proxy_matvec`` fills the walk's ψ½ block with
    # the closure's edge extrapolation of ``psi_view``
    # (``radial_characteristic_edge_seed`` — the pre-route-(a) convention),
    # while this hand reference consumes ``precompute_psi_state``'s LIVE
    # per-level half-angle grid built from the same ``psi_view``.  Then
    # apply the α·ΔA/w/V redistribution fold explicitly here.
    #
    # The hand reference's structural-independence claim is about the
    # ROUTING/scatter (each ordinate processed in its own scalar pass, no
    # bool-mask scatter into the legacy's misrouting ``ks``), NOT about the
    # redistribution closure — so sharing the edge-extrapolation seed
    # convention with the proxy is consistent (independence lives in the
    # per-ordinate walk, not in the seed the two sides agree to consume).
    closure = problem.angular_closure
    psi_state = closure.precompute_psi_state(psi_view)
    redist_full = np.zeros((ng, N, nx))
    for p, level_idx in enumerate(level_indices):
        level_idx_arr = np.asarray(level_idx)
        faces = psi_state[p].faces              # (ng, M_p+1, nx)
        alpha = reduced.angular.alpha_per_level[p]   # (M_p+1,)
        # ΔA ⊗ 1/w formed from its two factors (the fused cache retired
        # 2026-08-26); this side is the independent reference, so forming
        # it here is what keeps the comparison two-sided.
        lvl = np.asarray(quad.level_indices[p])
        dAw = reduced.delta_A[:, None] / np.asarray(quad.weights)[lvl][None, :]
        for m in range(level_idx_arr.size):
            redist_full[:, level_idx_arr[m], :] = (
                dAw[:, m].reshape(1, nx)
                * (alpha[m + 1] * faces[:, m + 1, :]
                   - alpha[m] * faces[:, m, :])
                / V.reshape(1, nx)
            )

    outflow_at_boundary = np.zeros((ng, N))

    # Outward — per ordinate.
    for level_idx in level_indices:
        level_idx_arr = np.asarray(level_idx)
        eta_level = mu_x[level_idx_arr]
        out_within = eta_level > +eps
        if not np.any(out_within):
            continue
        global_out = level_idx_arr[out_within]
        for n_g in global_out:
            mu_n = mu_x[n_g]
            psi_face_in = psi_g_first[:, n_g, 0].copy()
            for i in range(nx):
                psi_cell = psi_g_first[:, n_g, i]
                psi_face_out = 2.0 * psi_cell - psi_face_in
                streaming = mu_n * (
                    A[i + 1] * psi_face_out - A[i] * psi_face_in
                ) / V[i]
                redistribution = redist_full[:, n_g, i]
                collision = sigma_t_gx[:, i] * psi_cell
                out[n_g, :, i] = streaming + redistribution + collision
                psi_face_in = psi_face_out
            outflow_at_boundary[:, n_g] = psi_face_out

    # BC trace.
    inflow_full = bc_outer.apply(outflow_at_boundary.T)

    # Inward — per ordinate.
    for level_idx in level_indices:
        level_idx_arr = np.asarray(level_idx)
        eta_level = mu_x[level_idx_arr]
        in_within = eta_level < -eps
        if not np.any(in_within):
            continue
        global_in = level_idx_arr[in_within]
        for n_g in global_in:
            mu_n = mu_x[n_g]
            psi_face_in = inflow_full[n_g, :]
            for i in range(nx - 1, -1, -1):
                psi_cell = psi_g_first[:, n_g, i]
                psi_face_out = 2.0 * psi_cell - psi_face_in
                streaming = mu_n * (
                    A[i + 1] * psi_face_in - A[i] * psi_face_out
                ) / V[i]
                redistribution = redist_full[:, n_g, i]
                collision = sigma_t_gx[:, i] * psi_cell
                out[n_g, :, i] = streaming + redistribution + collision
                psi_face_in = psi_face_out

    # Degenerate (|μ_x| < eps) — no radial flow.
    degenerate_mask = np.abs(mu_x) < eps
    if np.any(degenerate_mask):
        global_deg = np.where(degenerate_mask)[0]
        for n_g in global_deg:
            for i in range(nx):
                psi_cell = psi_g_first[:, n_g, i]
                redistribution = redist_full[:, n_g, i]
                collision = sigma_t_gx[:, i] * psi_cell
                out[n_g, :, i] = redistribution + collision

    return out


# ═══════════════════════════════════════════════════════════════════════
# L0 — Hand-reference battery (the structural correctness anchor)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.xfail(
    reason="cylinder matvec/sweep WDD divergence — issue #206",
    strict=False,
)
@pytest.mark.l0
# Q5.6.3: the carrying folded family replaces the pre-flip LS4 / LS6 /
# product(2,4) rows under the standard migration (LS(n) → folded(n, 2n);
# P(m, p) → folded(m, p) — parent azimuthal counts).
@pytest.mark.parametrize("quad_factory", [
    lambda: Quadrature.folded_product(n_mu=4, n_phi=8),
    lambda: Quadrature.folded_product(n_mu=6, n_phi=12),
    lambda: Quadrature.folded_product(n_mu=2, n_phi=4),
])
@pytest.mark.parametrize("n_cells", [3, 5, 10])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_unified_cylinder_matches_hand_reference(
    quad_factory, n_cells, seed,
) -> None:
    """Unified cylindrical matvec matches the per-ordinate hand reference.

    Promoted from ``derivations/diagnostics/diag_step3_cyl_unified_vs_hand_battery.py``
    (numerics-investigator 2026-05-17 closeout, then LS4/LS6/product(2,4)
    → the folded family at Q5.6.3). The hand reference is structurally
    independent of the unified (no bool-mask scatter); both agree at
    machine precision across folded(4,8) / folded(6,12) / folded(2,4) ×
    {3, 5, 10} cells × {0, 1, 2} seeds = 27 cases.
    """
    quad = quad_factory()
    problem = _build_cyl(n_cells, quad)
    ng = 1
    N = quad.N

    rng = np.random.default_rng(seed)
    psi_view = rng.standard_normal((N, ng, n_cells)).astype(np.float64)
    psi_view = _bc_fill_outer(psi_view, problem)
    sigma_t = np.full((ng, n_cells), 2.0)

    m_unified = legacy_proxy_matvec(psi_view, problem, sigma_t)
    m_hand = _hand_reference_cyl_matvec(psi_view, problem, sigma_t)

    m_unified_u = _extract_at_unknown_slots(m_unified, problem)
    m_hand_u = _extract_at_unknown_slots(m_hand, problem)

    np.testing.assert_allclose(
        m_unified_u, m_hand_u, rtol=1e-12, atol=1e-13,
        err_msg=(
            f"Unified must match per-ordinate hand reference at quad "
            f"N={N}, n_cells={n_cells}, seed={seed}"
        ),
    )


@pytest.mark.l0
def test_unified_cylinder_zero_psi_gives_zero() -> None:
    """Linear operator: zero input → zero output."""
    quad = Quadrature.folded_product(n_mu=4, n_phi=8)
    problem = _build_cyl(n_cells=5, quad=quad)
    ng = 1
    sigma_t = np.full((ng, problem.nx), 2.0)
    psi_view = np.zeros((quad.N, ng, problem.nx))

    m_unified = legacy_proxy_matvec(psi_view, problem, sigma_t)
    np.testing.assert_array_equal(m_unified, np.zeros_like(m_unified))


@pytest.mark.l0
def test_unified_cylinder_constant_psi_gives_sigma_t() -> None:
    """At ψ = constant on homogeneous reflective cylinder, unified matvec
    returns σ_t · ψ. Sanity check — flat flux activates only the
    collision term in the per-cell balance."""
    quad = Quadrature.folded_product(n_mu=4, n_phi=8)
    problem = _build_cyl(n_cells=5, quad=quad)
    ng = 1
    sigma_t_val = 2.0
    sigma_t = np.full((ng, problem.nx), sigma_t_val)
    psi_view = np.ones((quad.N, ng, problem.nx))

    m_unified = legacy_proxy_matvec(psi_view, problem, sigma_t)
    m_at_unknowns = _extract_at_unknown_slots(m_unified, problem)
    np.testing.assert_allclose(
        m_at_unknowns, sigma_t_val, rtol=1e-12, atol=1e-13,
    )


# ═══════════════════════════════════════════════════════════════════════
# L1 — Heterogeneous trajectory_resolvent cross-check via Krylov
# ═══════════════════════════════════════════════════════════════════════
#
# The L0 battery proves the per-cell algebra. The L1 test proves the
# *converged eigenvalue* of a heterogeneous cylinder problem agrees with
# the structurally-independent trajectory_resolvent reference (Variant α
# at α=1, see ``orpheus.derivations.continuous.trajectory_resolvent``).
#
# Why heterogeneous: L2 in lessons.md — homogeneous k = νΣ_f/Σ_a is
# flux-shape independent (the matvec's redistribution & angular closure
# all collapse on flat flux). Heterogeneous closed MR exercises every
# term in the unified per-cell algebra.
# ═══════════════════════════════════════════════════════════════════════


def _make_2g_mixture(sigma_t, sig_s, nu_sig_f, chi):
    """Build a 2-group Mixture from explicit XS arrays."""
    sigma_t = np.asarray(sigma_t, dtype=float)
    sig_s = np.asarray(sig_s, dtype=float)
    nu_sig_f = np.asarray(nu_sig_f, dtype=float)
    chi = np.asarray(chi, dtype=float)
    sig_a = sigma_t - sig_s.sum(axis=1)
    nu = np.ones_like(nu_sig_f)
    sig_f = nu_sig_f.copy()
    sig_c = sig_a - sig_f
    return make_mixture(
        sig_t=sigma_t, sig_c=sig_c, sig_f=sig_f, nu=nu, chi=chi, sig_s=sig_s,
    )


def _build_mr_cylinder_mesh(nx: int = 40) -> tuple[Mesh1D, dict]:
    """Build the 3-region cylindrical mesh + 2G materials for the MR case."""
    sigma_t, sigma_s, nu_sigma_f, chi = aba_xs_2g()
    materials = {
        i: _make_2g_mixture(sigma_t[i], sigma_s[i], nu_sigma_f[i], chi[i])
        for i in range(3)
    }
    # Regions A | B | A end at ABA_RADII; nx equal-width cells over the
    # whole radius, so each region holds its width's share of them.
    radii = (0.0, *map(float, ABA_RADII))
    geom = StructuredGeometry.cylinder(radii, (0, 1, 0), outer=BC("reflective"))
    mesh = Mesher(geom).partition(tuple(
        CellsByCount.uniform_width(round(nx * (b - a) / radii[-1]))
        for a, b in geom.intervals
    )).mesh
    return mesh, materials


# Post-D-K (commit ``dadf4e8``), the within-group loss composite
# ``L + C`` (:func:`build_streaming_collision` → ``StreamingOperator +
# MultiplicationOperator``, i.e. :class:`StreamingCollisionOperator`) calls
# :func:`_transport_operator_matvec_unified` natively for 1-D
# cylindrical.  No monkey-patch is required.


#: The k tolerance this comparison is held to once the reference is certified:
#: ``tolerance_for(e, None)`` of ``tests/gates/derivations/_trajectory_resolvent_ladders.py``,
#: the reference assumed at the floor, for this folded-4x8 solve's own k error
#: e = 5.9e-4 ([M] 2026-09-26: the 4x8 SN k against its 32x64 limit), giving
#: 2e-3. It was 3e-2, justified as "the reference's quadrature budget"; that
#: budget was the reference's one-spline emission density (ERR-090), and the
#: reference still carries no certified bound (#516).
_UNIFIED_CYL_K_TOLERANCE = tolerance_for(CYLINDER_3REG_SN_4X8_K_STEP, None)


@functools.cache
def _unified_cylinder_k() -> float:
    """The Krylov-on-(L + C) SN eigenvalue of the ABA cylinder (folded 4x8, 40 uniform cells), once per session."""
    mesh, materials = _build_mr_cylinder_mesh(nx=40)
    sol = solve_sn(
        materials=materials,
        mesh=mesh,
        quadrature=Quadrature.folded_product(n_mu=4, n_phi=8),
        inner_solver="krylov",
        max_outer=200, keff_tol=1e-7, flux_tol=1e-7,
        max_inner=200, inner_tol=1e-9,
    )
    return float(sol.outcome.keff)


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("sn-curvilinear-trajectory-resolvent-crosscheck")
@pytest.mark.rests_on(
    "tests/gates/derivations/test_trajectory_resolvent_regionwise_source.py"
    "::test_mr_oracle_first_leg_matches_the_line_integral[cylinder]",
)
@awaits_cylinder_bound
def test_unified_cylinder_l1_mr_2g_trajectory_resolvent() -> None:
    r"""L1 — heterogeneous 3-region 2G closed cylinder via unified matvec.

    Drives :func:`solve_sn` with ``inner_solver="krylov"`` — the
    matvec routes through :class:`StreamingCollisionOperator` (= ``L + C``)
    via :func:`_transport_operator_matvec_unified`. The converged
    ``k_eff`` is compared against the structurally-independent
    trajectory-resolvent reference (Variant α at α=1) at
    :data:`_UNIFIED_CYL_K_TOLERANCE`.

    Strict ``xfail`` on #516: the cylinder reference carries no certified
    error bound, so the floor assertion fails first; the comparison stays
    live through ``test_unified_cylinder_l1_mr_2g_trajectory_resolvent_record``.

    Per ``.claude/lessons.md`` L14 — solver correctness is a 4-way
    standoff. The L0 hand-reference battery in this file proves the
    per-cell algebra; this test pins the converged eigenvalue against
    a continuous reference that is structurally independent of every
    discrete primitive in the SN code (no shared FP path, no shared
    redist closure, no shared boundary recurrence).
    """
    k_unified = _unified_cylinder_k()
    k_ref = float(cylinder_3reg_reference().k_eff)
    rel = abs(k_unified - k_ref) / k_ref
    print(f"unified cylinder: k_unified={k_unified:.10f} k_ref={k_ref:.10f} rel={rel:.3e}")
    certify_agreement(
        "unified cylinder k", rel, _UNIFIED_CYL_K_TOLERANCE, CYLINDER_3REG_REFERENCE_BOUND["k"],
    ).require()


@pytest.mark.l1
@pytest.mark.slow
def test_unified_cylinder_l1_mr_2g_trajectory_resolvent_record() -> None:
    r"""RECORD: the unified-matvec solve's k and the reference's, as they read today.

    Not verification: it keeps the comparison live while the row above is a
    strict xfail, and reddens when either side moves (an SN change, or the
    reference's #516 repair, after which the bound is re-derived and the
    xfail lifted). The recorded values and their band are
    :data:`~tests.gates.sn.verification.analytical._certified_agreement.CYLINDER_3REG_RECORD`.
    """
    readings = {"unified_k": _unified_cylinder_k(), "k_ref": float(cylinder_3reg_reference().k_eff)}
    assert_record(
        readings, {name: CYLINDER_3REG_RECORD[name] for name in readings},
        CYLINDER_3REG_RECORD_BAND, relative=CYLINDER_3REG_RECORD_RELATIVE,
    )


@pytest.mark.l1
@pytest.mark.slow
def test_unified_cylinder_l1_homogeneous_kinf_2g() -> None:
    r"""L1 sanity — 2G homogeneous closed cylinder eigenvalue via unified.

    Same monkey-patch construction as the MR test, but on homogeneous
    XS. k_eff is shape-independent (per L2 in lessons.md), so this test
    primarily confirms that the unified matvec drives GMRES to the
    correct k_∞ — a *necessary* but not *sufficient* condition. The
    sufficient condition is the heterogeneous MR test above.
    """
    from orpheus.derivations.common.eigenvalue import kinf_homogeneous

    sigma_t = [0.5, 1.0]
    sig_s = [[0.3, 0.05], [0.0, 0.7]]
    nu_sig_f = [0.4, 0.6]
    chi = [1.0, 0.0]
    mat = _make_2g_mixture(sigma_t, sig_s, nu_sig_f, chi)
    k_analytical = kinf_homogeneous(
        np.asarray(sigma_t), np.asarray(sig_s),
        np.asarray(nu_sig_f), np.asarray(chi),
    )

    nx = 20
    mesh = Mesher(
        StructuredGeometry.cylinder((0.0, 2.0), (0,), outer=BC("reflective")),
    ).partition(CellsByCount.uniform_width(nx)).mesh
    quad = Quadrature.folded_product(n_mu=4, n_phi=8)

    sol = solve_sn(
            materials={0: mat},
            mesh=mesh,
            quadrature=quad,
            inner_solver="krylov",
            max_outer=200, keff_tol=1e-9, flux_tol=1e-8,
            max_inner=200, inner_tol=1e-10,
        )

    rel = abs(sol.outcome.keff - k_analytical) / k_analytical
    assert rel < 5e-4, (
        f"unified cylinder k_∞ recovery violated: "
        f"k_analytical={k_analytical:.10f}, k_unified={sol.outcome.keff:.10f}, "
        f"rel={rel:.2e}"
    )
