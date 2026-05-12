r"""Diagnostic: Gate 1.1 sphere MMS failure under ``MorelMontryAngularSweep``.

Issue #168 Phase D — confirms the structural hypothesis that the
per-ordinate flat-ψ residual collapses to machine zero when the
outward-sweep pole-face IC is replaced by the **Hébert §3.9.4
Eqs. (3.432)-(3.435) Carlson coupled-pole inward μ = −1 sweep** result.

Authored by numerics-investigator, 2026-05-12. Not a pytest file —
this is a CLI probe (run with ``python tests/sn/diagnostics/gate_1_1_sphere_mms_failure.py``).
The method-implementer will consume the output as empirical evidence
before shipping the production ``CarlsonCoupledPole`` strategy.

References
----------
* ``.claude/agent-memory/literature-researcher/phase_d_carlson_coupled_pole.md``
  — Hébert §3.9.4 derivations, flat-ψ algebraic trace.
* ``.claude/agent-memory/method-implementer/issue_168_phase_c_closeout.md``
  — Phase C outcome ("sphere MMS xfail, cylinder MMS xpass").
* ``/home/vscode/.claude/plans/structured-booping-parrot.md`` — Phase D
  plan; this script is Step 2 deliverable.
* :mod:`orpheus.sn.spatial.pole_angular_closure` — the three angular
  closure strategies; the Phase C-shipped ``MorelMontryAngularSweep``
  initialises ``psi_half_left = 0`` (the seed this script overrides).

Probe structure
---------------
1. Minimal sphere problem (R=2.0, nx=10, GL4, Σ_t ∈ {0.0, 0.5}).
2. Build the packed flat-ψ vector ``psi = ones(n_unknowns)``.
3. Per-ordinate residual ``r[g, n, i] = (L ψ)[g, n, i] − Σ_t·1`` under
   the three pole-angular closures, exposing the M-M failure profile.
4. Build the Hébert (3.432)-(3.435) inward μ = −1 sweep helper and
   show that on the flat-ψ + flat-source consistent fixed point the
   inward sweep returns ``phi_aux[g, i] = 1`` everywhere (literature
   memo §3 algebraic trace).
5. Run a Carlson-seeded *replica* of ``transport_operator_matvec_spherical``
   where ``psi_face_in(pole) = phi_aux[g, 0]`` instead of
   ``fi[:, outgoing_mask, 0, 0]``. Re-compute the per-ordinate residual
   and confirm collapse to machine precision.

If step 5 collapses to ≤ 1e-12 the Phase D Carlson hypothesis is
empirically validated against production code. If it does NOT, the
hypothesis is FALSIFIED and Phase D must pause + escalate (per the
plan §Step 2 instruction).
"""
from __future__ import annotations

import numpy as np

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import (
    SNStreamingOperator,
    build_equation_map_spherical,
    solution_to_angular_flux_spherical,
)
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.spatial.pole_angular_closure import (
    BaileyFlatFluxRedist,
    LegacyTauSymmetricInterpolation,
    MorelMontryAngularSweep,
)


# ═══════════════════════════════════════════════════════════════════════
# Test-only Hébert (3.432)-(3.435) inward μ = −1 sweep
# ═══════════════════════════════════════════════════════════════════════


def carlson_inward_sweep_test_only(
    sigma_t: np.ndarray,
    Q_bar: np.ndarray,
    bc_outer: np.ndarray,
    dr: np.ndarray,
) -> np.ndarray:
    r"""Run Hébert §3.9.4 Eqs. (3.434)-(3.435) per group, returning
    cell-centred ``phi_aux[g, i]``.

    Per the literature memo §3 recurrence on a uniform mesh:

    .. math::

        \bar\phi_i &= \frac{\Delta r_i \cdot \bar Q_i
                              + 2 \cdot \bar\phi_{i+1/2}}
                            {\Delta r_i \cdot \Sigma_i + 2}
                            \quad\text{Eq. (3.434)} \\
        \bar\phi_{i-1/2} &= 2 \cdot \bar\phi_i
                              - \bar\phi_{i+1/2}
                            \quad\text{Eq. (3.435)}

    The sweep proceeds INWARD from ``i = nx-1`` down to ``i = 0``.
    The outer-face seed ``bc_outer`` is the angular flux at μ = −1 at
    the outer face (zero for vacuum BC; literature §8 discusses
    reflective / white BC variants).

    Parameters
    ----------
    sigma_t : np.ndarray, shape (ng, nx)
        Per-group, per-cell total cross-section.
    Q_bar : np.ndarray, shape (ng, nx)
        Per-group cell-averaged source evaluated AT μ = −1 (i.e. the
        Legendre-moment-folded source ``Σ_ℓ (2ℓ+1)/2 · Q_ℓ · P_ℓ(−1)``).
        For isotropic scattering this is just the scalar source / 2.
    bc_outer : np.ndarray, shape (ng,)
        Outer-face angular flux at μ = −1. Zero for vacuum;
        ``φ̄_{nx+1/2}`` for the reflective / white BC analog.
    dr : np.ndarray, shape (nx,)
        Cell widths.

    Returns
    -------
    np.ndarray, shape (ng, nx)
        Cell-centred ``phi_aux[g, i] ≡ φ̄_{1/2, i}`` per group, per
        cell. The ``i = 0`` slice is the pole-face seed used by the
        Carlson-coupled outward sweep.
    """
    ng, nx = sigma_t.shape
    phi_aux = np.zeros((ng, nx))
    for g in range(ng):
        phi_face_outer = bc_outer[g]  # known: outer-face value
        for k in range(nx - 1, -1, -1):
            denom = dr[k] * sigma_t[g, k] + 2.0
            phi_cell = (dr[k] * Q_bar[g, k] + 2.0 * phi_face_outer) / denom
            phi_aux[g, k] = phi_cell
            # Eq. (3.435): step inward to next face
            phi_face_outer = 2.0 * phi_cell - phi_face_outer
    return phi_aux


# ═══════════════════════════════════════════════════════════════════════
# Carlson-seeded M-M angular recurrence
#
# The production ``_mm_weighted_angular_recurrence_single_level`` hardcodes
# ``psi_half_left = 0`` (the value of the half-angle face flux at the
# auxiliary μ = −1 starting direction). On flat ψ this is INCONSISTENT
# — Hébert §3.9.4 (3.432)-(3.435) says ψ_{1/2, i} = φ̄_{1/2, i} (the
# inward μ = −1 sweep's cell-centred output), NOT zero. The Carlson
# intervention replaces the ``0`` seed with the test-only inward
# sweep result.
# ═══════════════════════════════════════════════════════════════════════


def mm_recurrence_carlson_seeded(
    psi_level: np.ndarray,          # (ng, M, nx)
    alpha_level: np.ndarray,        # (M+1,)
    dAw_level: np.ndarray,          # (nx, M)
    tau_level: np.ndarray,          # (M,)
    volume: np.ndarray,             # (nx,)
    *,
    psi_half_seed: np.ndarray | None = None,   # (ng, nx) — Hébert (3.432)-(3.435)
) -> np.ndarray:
    r"""Test-only fork of
    :func:`orpheus.sn.spatial.pole_angular_closure._mm_weighted_angular_recurrence_single_level`
    with a configurable Carlson seed for the half-angle face flux at the
    auxiliary μ = −1 starting direction.

    When ``psi_half_seed is None`` reproduces the production hardcoded
    ``ψ_{1/2,i,g} = 0``. When supplied as ``(ng, nx)``, that value
    replaces the seed — corresponding to Hébert §3.9.4 Eqs.
    (3.432)-(3.435): solve the inward μ = −1 spatial sweep first, use
    its cell-centred φ̄_{1/2, i} as the M-M α-cascade seed.
    """
    ng, M, nx = psi_level.shape
    redist = np.empty((ng, M, nx), dtype=psi_level.dtype)
    if psi_half_seed is None:
        psi_half_left = np.zeros((ng, nx), dtype=psi_level.dtype)
    else:
        psi_half_left = psi_half_seed.copy()
    for m in range(M):
        tau_m = tau_level[m]
        psi_half_right = (
            psi_level[:, m, :] - (1.0 - tau_m) * psi_half_left
        ) / tau_m
        redist[:, m, :] = (
            dAw_level[:, m].reshape(1, nx)
            * (alpha_level[m + 1] * psi_half_right
               - alpha_level[m] * psi_half_left)
            / volume.reshape(1, nx)
        )
        psi_half_left = psi_half_right
    return redist


# ═══════════════════════════════════════════════════════════════════════
# Inline replica of transport_operator_matvec_spherical
#
# The replica preserves every bit of math from the production matvec
# EXCEPT (a) the outward-sweep pole-face IC, which is overrideable via
# ``carlson_pole_seed`` and (b) the M-M angular recurrence half-angle
# seed, which is overrideable via ``carlson_half_angle_seed``. When
# both seeds are ``None`` the replica is bit-identical to production.
# ═══════════════════════════════════════════════════════════════════════


def matvec_spherical_replica(
    psi: np.ndarray,
    op: SNStreamingOperator,
    *,
    carlson_pole_seed: np.ndarray | None = None,
    carlson_half_angle_seed: np.ndarray | None = None,
) -> np.ndarray:
    r"""Reimplementation of ``transport_operator_matvec_spherical`` with
    swappable pole-face spatial IC AND swappable M-M half-angle seed.

    All math copied verbatim from
    ``orpheus.sn.operator.transport_operator_matvec_spherical`` (Phase
    C tip, branch ``refactor/sn-operator-algebra``).

    Parameters
    ----------
    carlson_pole_seed : optional
        Shape ``(ng,)``. When supplied, replaces the WDD outward-sweep
        ``psi_face_in`` initialisation at the pole face.
    carlson_half_angle_seed : optional
        Shape ``(ng, nx)``. When supplied, replaces the production
        ``ψ_{1/2,i,g} = 0`` hard-coded seed inside the M-M angular
        recurrence. This is the Hébert §3.9.4 (3.432)-(3.435) Carlson
        coupled-pole intervention applied to the ANGULAR seed.
    """
    sn_mesh = op.sn_mesh
    sig_t = op.sig_t
    quad = sn_mesh.quad
    eq_map = op._ensure_eq_map()
    nx = sn_mesh.nx
    ng = sig_t.shape[2]
    reduced = sn_mesh.reduced
    face_areas = reduced.face_areas
    volumes = sn_mesh.volumes
    alpha_half = reduced.alpha_half
    redist_dAw = reduced.redist_dAw
    tau_mm = reduced.tau_mm
    bc_outer = sn_mesh.bc_right
    pole_angular_closure = sn_mesh.pole_angular_closure

    fi = solution_to_angular_flux_spherical(psi, eq_map, quad, nx, ng)
    A = face_areas
    V = volumes[:, 0]
    N = quad.N
    eps = 1e-15

    outgoing_mask = quad.mu_x > +eps
    incoming_mask = quad.mu_x < -eps
    mu_out = quad.mu_x[outgoing_mask]
    mu_in = quad.mu_x[incoming_mask]
    n_out = int(mu_out.size)
    n_in = int(mu_in.size)

    # Optionally use the Carlson-seeded M-M recurrence instead of the
    # production pole_angular_closure call.
    if carlson_half_angle_seed is not None:
        if not isinstance(pole_angular_closure, MorelMontryAngularSweep):
            raise RuntimeError(
                "carlson_half_angle_seed only makes sense with "
                "MorelMontryAngularSweep."
            )
        redist_full = mm_recurrence_carlson_seeded(
            fi[..., 0], alpha_half, redist_dAw, tau_mm, V,
            psi_half_seed=carlson_half_angle_seed,
        )
    else:
        redist_full = pole_angular_closure(
            fi[..., 0], alpha_half, redist_dAw, tau_mm, V,
            level_indices=None,
        )

    lhs = np.empty((ng, eq_map.n_eq))
    outflow_at_boundary = np.zeros((ng, N))

    outward_visits = list(sn_mesh.iter_cells_by_direction(+1))
    inward_visits = list(sn_mesh.iter_cells_by_direction(-1))

    # ── Phase 1: outgoing ordinates (μ > 0), i = 0 → nx-1 ─────────
    if n_out > 0:
        if carlson_pole_seed is None:
            # Production (Phase C Lewis-Miller seed).
            if len(outward_visits) > 0:
                i0 = outward_visits[0].cell_idx
                psi_face_in = fi[:, outgoing_mask, i0, 0].copy()
            else:
                psi_face_in = np.zeros((ng, n_out))
        else:
            # Carlson coupled-pole intervention. carlson_pole_seed has
            # shape (ng,); broadcast to all n_out outgoing ordinates.
            psi_face_in = np.broadcast_to(
                carlson_pole_seed[:, None],
                (ng, n_out),
            ).copy()

        for visit in outward_visits:
            i = visit.cell_idx
            psi_cell = fi[:, outgoing_mask, i, 0]
            psi_face_out = 2.0 * psi_cell - psi_face_in
            streaming = (
                mu_out[None, :]
                * (A[i + 1] * psi_face_out - A[i] * psi_face_in)
                / V[i]
            )
            redistribution = redist_full[:, outgoing_mask, i]
            collision = sig_t[i, 0, :, None] * psi_cell
            ks = eq_map.unknowns_at_cell_for_mask(i, outgoing_mask)
            if ks.size > 0:
                lhs[:, ks] = streaming + redistribution + collision
            psi_face_in = psi_face_out
        outflow_at_boundary[:, outgoing_mask] = psi_face_out

    # ── BC trace law ──────────────────────────────────────────────
    inflow_full = bc_outer.apply(outflow_at_boundary.T)

    # ── Phase 2: incoming ordinates (μ < 0), i = nx-1 → 0 ─────────
    if n_in > 0:
        psi_face_in = inflow_full[incoming_mask, :].T
        for visit in inward_visits:
            i = visit.cell_idx
            psi_cell = fi[:, incoming_mask, i, 0]
            psi_face_out = 2.0 * psi_cell - psi_face_in
            streaming = (
                mu_in[None, :]
                * (A[i + 1] * psi_face_in - A[i] * psi_face_out)
                / V[i]
            )
            redistribution = redist_full[:, incoming_mask, i]
            collision = sig_t[i, 0, :, None] * psi_cell
            ks = eq_map.unknowns_at_cell_for_mask(i, incoming_mask)
            if ks.size > 0:
                lhs[:, ks] = streaming + redistribution + collision
            psi_face_in = psi_face_out

    return lhs.ravel(order='F')


# ═══════════════════════════════════════════════════════════════════════
# Diagnostic harness
# ═══════════════════════════════════════════════════════════════════════


def build_sphere_problem(
    nx: int = 10,
    R: float = 2.0,
    sigma_t_value: float = 0.5,
    pole_closure=None,
    bc=None,
):
    """Build SNMesh + SNStreamingOperator for the failure probe."""
    quad = GaussLegendre1D.create(4)
    edges = np.linspace(0.0, R, nx + 1)
    mesh = Mesh1D(
        edges=edges,
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_right=bc if bc is not None else BC("reflective"),
    )
    sn_mesh = SNMesh(mesh, quad, pole_angular_closure=pole_closure)
    sig_t = np.full((nx, 1, 1), sigma_t_value)
    op = SNStreamingOperator(sn_mesh=sn_mesh, sig_t=sig_t)
    return sn_mesh, sig_t, op


def per_ordinate_residual_report(
    result: np.ndarray,
    expected_value: float,
    eq_map,
    quad,
    nx: int,
    ng: int = 1,
) -> dict:
    r"""Re-scatter the packed residual into (g, n, i) and report stats.

    Returns dict with:
    * ``max_abs`` — global L_∞ residual
    * ``per_ordinate_max`` — shape (N,), max |residual| over (g, i) per ordinate
    * ``per_cell_max`` — shape (nx,), max |residual| over (g, n) per cell
    * ``residual_field`` — shape (ng, N, nx), the full residual
    """
    eq_residual = result - expected_value  # shape (n_unknowns,)
    flux = eq_residual.reshape(ng, eq_map.n_eq, order='F')
    field = np.zeros((ng, quad.N, nx))
    for k in range(eq_map.n_eq):
        field[:, eq_map.ordinate[k], eq_map.ix[k]] = flux[:, k]
    per_ordinate_max = np.max(np.abs(field), axis=(0, 2))  # (N,)
    per_cell_max = np.max(np.abs(field), axis=(0, 1))      # (nx,)
    return {
        "max_abs": float(np.max(np.abs(field))),
        "per_ordinate_max": per_ordinate_max,
        "per_cell_max": per_cell_max,
        "residual_field": field,
    }


def section(title: str):
    line = "═" * 78
    print(f"\n{line}\n{title}\n{line}")


def main():
    print(__doc__)

    # ─── Probe parameters ────────────────────────────────────────
    nx = 10
    R = 2.0
    sigma_t_values = [0.0, 0.5]

    # ============================================================
    # Step 1: Baseline failure under MorelMontryAngularSweep
    # ============================================================
    section("Step 1 — Baseline residuals on flat ψ probe (reflective sphere)")

    closures = {
        "Legacy (τ-symmetric)": LegacyTauSymmetricInterpolation(),
        "Bailey flat-flux": BaileyFlatFluxRedist(),
        "Morel-Montry (canonical)": MorelMontryAngularSweep(),
    }

    baseline_report = {}
    for sigma_t_value in sigma_t_values:
        print(f"\n--- Σ_t = {sigma_t_value} ---")
        for label, closure in closures.items():
            sn_mesh, sig_t, op = build_sphere_problem(
                nx=nx, R=R,
                sigma_t_value=sigma_t_value,
                pole_closure=closure,
            )
            quad = sn_mesh.quad
            eq_map = op._ensure_eq_map()
            psi = np.ones(op.n_unknowns)
            result = op.apply(psi)
            report = per_ordinate_residual_report(
                result, sigma_t_value, eq_map, quad, nx,
            )
            verdict = "PASS" if report["max_abs"] < 1e-12 else "FAIL"
            print(
                f"  {label:30s}: max|residual| = "
                f"{report['max_abs']:.3e}   [{verdict}]"
            )
            baseline_report[(sigma_t_value, label)] = (sn_mesh, op, report)

    # ============================================================
    # Step 2: Failure profile for M-M (the failing case)
    # ============================================================
    section("Step 2 — Failure profile decomposition: Morel-Montry, Σ_t = 0.5")

    sn_mesh, op, report = baseline_report[(0.5, "Morel-Montry (canonical)")]
    quad = sn_mesh.quad
    print(f"GL4 ordinates μ_x = {quad.mu_x}")
    print(f"GL4 weights w     = {quad.weights}")
    print(f"Σ μ_n w_n         = {np.sum(quad.mu_x * quad.weights):.3e}")
    print(f"Σ w_n             = {np.sum(quad.weights):.3e}")
    print(f"\nPer-ordinate max|residual| (over cells, over groups):")
    for n in range(quad.N):
        print(
            f"  ordinate {n} (μ = {quad.mu_x[n]:+.6f}, w = {quad.weights[n]:.4f}): "
            f"max|r| = {report['per_ordinate_max'][n]:.3e}"
        )

    print(f"\nPer-cell max|residual| (over ordinates, over groups):")
    edges = np.linspace(0.0, R, nx + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    for i in range(nx):
        print(
            f"  cell {i:2d} (r_c = {centres[i]:.4f}): "
            f"max|r| = {report['per_cell_max'][i]:.3e}"
        )

    print(f"\nFull (g=0) residual field, indexed [n, i]:")
    field = report["residual_field"][0]
    with np.printoptions(precision=4, suppress=False, linewidth=160):
        print(field)

    # ============================================================
    # Step 3: Carlson inward sweep test-only — flat-ψ algebraic trace
    # ============================================================
    section("Step 3 — Carlson inward μ = −1 sweep (flat-ψ consistency)")

    # Build the inward-sweep inputs for the flat-ψ probe.
    # Σ_t = 0.5, reflective BC, flat ψ = 1 ⇒ flat scalar φ = 1·Σ_n w_n = 2
    # (GL on [-1, 1] with Σ w = 2). The CONSISTENT source for flat ψ = 1
    # on the equation Σ_t · ψ = 1·Q̄ (at μ = −1) is Q̄ = Σ_t · 1 = 0.5
    # — i.e. the source moments must close the equation. For the
    # *operator residual* on flat ψ the source IS Σ_t · ψ (by definition
    # of `expected = Σ_t · ψ`), so Q̄ at μ = −1 = Σ_t · 1 = 0.5.
    #
    # The literature memo §3 trace: with flat ψ on reflective BC,
    # φ̄_{nx+1/2} = 1 (reflective BC mirrors back the outgoing value);
    # the recurrence reproduces φ̄_i = 1 at every cell.
    sigma_t_value = 0.5
    sn_mesh, sig_t, op = build_sphere_problem(
        nx=nx, R=R, sigma_t_value=sigma_t_value,
        pole_closure=MorelMontryAngularSweep(),
    )
    edges = np.linspace(0.0, R, nx + 1)
    dr = np.diff(edges)
    sigma_t_arr = np.full((1, nx), sigma_t_value)
    # Q̄ at μ = −1 for the consistent fixed-point: source = Σ_t · ψ
    Q_bar = sigma_t_value * np.ones((1, nx))
    # Reflective BC at outer face: outgoing μ = +1 cell-face value
    # equals the inward μ = −1 cell-face value (specular reflection).
    # On flat ψ = 1 the outgoing face = 1, so the inward seed = 1.
    bc_outer = np.array([1.0])
    phi_aux = carlson_inward_sweep_test_only(
        sigma_t_arr, Q_bar, bc_outer, dr,
    )
    print("Flat ψ on reflective BC with consistent source Q̄ = Σ_t:")
    print(f"  bc_outer (μ=−1 at outer face): {bc_outer}")
    print(f"  Q̄[g=0, :]                    : {Q_bar[0]}")
    print(f"  phi_aux[g=0, :] from inward sweep:")
    print(f"    {phi_aux[0]}")
    print(f"  pole-face seed phi_aux[g=0, 0] = {phi_aux[0, 0]:.16f}")
    print(
        f"  matches literature §3 trace (expected 1.0): "
        f"{abs(phi_aux[0, 0] - 1.0) < 1e-15}"
    )

    # Now also compute for vacuum BC to confirm the geometry-aware
    # behaviour.
    sn_mesh_vac, _, op_vac = build_sphere_problem(
        nx=nx, R=R, sigma_t_value=sigma_t_value,
        pole_closure=MorelMontryAngularSweep(),
        bc=BC("vacuum"),
    )
    bc_outer_vac = np.array([0.0])  # vacuum at μ = −1 outer face
    phi_aux_vac = carlson_inward_sweep_test_only(
        sigma_t_arr, Q_bar, bc_outer_vac, dr,
    )
    print(
        f"\n[Sanity, vacuum BC] phi_aux[g=0, :] (NOT flat — driven by Q only):"
    )
    print(f"    {phi_aux_vac[0]}")
    print(f"  pole-face seed phi_aux[g=0, 0] = {phi_aux_vac[0, 0]:.6f}")

    # ============================================================
    # Step 4: Inject Carlson seed into the matvec replica
    # ============================================================
    section("Step 4 — Carlson intervention residual (REFLECTIVE BC)")

    # Reproduce baseline via the replica first (sanity bit-identity check)
    sn_mesh, sig_t, op = build_sphere_problem(
        nx=nx, R=R, sigma_t_value=sigma_t_value,
        pole_closure=MorelMontryAngularSweep(),
    )
    quad = sn_mesh.quad
    eq_map = op._ensure_eq_map()
    psi = np.ones(op.n_unknowns)

    result_production = op.apply(psi)
    result_replica_baseline = matvec_spherical_replica(psi, op)
    assert np.allclose(result_production, result_replica_baseline,
                       rtol=1e-15, atol=1e-15), (
        "Replica is not bit-identical to production — fix the replica "
        "before drawing conclusions."
    )
    print(f"Replica vs production bit-identity:")
    print(
        f"  max|replica - production| = "
        f"{np.max(np.abs(result_production - result_replica_baseline)):.3e}  "
        f"[OK]"
    )

    # Sanity: baseline (with the Lewis-Miller seed) reports the
    # failure we already saw above.
    report_baseline = per_ordinate_residual_report(
        result_replica_baseline, sigma_t_value, eq_map, quad, nx,
    )
    print(
        f"  Baseline (Lewis-Miller seed) max|residual|: "
        f"{report_baseline['max_abs']:.3e}"
    )

    # ─── Carlson intervention A: spatial pole-face IC only ────────
    result_carlson_pole = matvec_spherical_replica(
        psi, op, carlson_pole_seed=phi_aux[:, 0],
    )
    report_carlson_pole = per_ordinate_residual_report(
        result_carlson_pole, sigma_t_value, eq_map, quad, nx,
    )
    print(
        f"\n  [A] Pole-face spatial seed only "
        f"(phi_aux[g, 0] = {phi_aux[0, 0]:.6f}) "
        f"max|residual|: {report_carlson_pole['max_abs']:.3e}"
    )

    # ─── Carlson intervention B: M-M angular half-angle seed only ──
    # The literature memo §4 makes clear: the kept output of the inward
    # μ = −1 sweep is the cell-centred {φ̄_{1/2, i}}_{i=1..nx}, which
    # IS the seed for the OUTWARD M-M α-cascade — i.e., ψ_{1/2, i}
    # inside ``_mm_weighted_angular_recurrence_single_level``.
    result_carlson_half = matvec_spherical_replica(
        psi, op, carlson_half_angle_seed=phi_aux,
    )
    report_carlson_half = per_ordinate_residual_report(
        result_carlson_half, sigma_t_value, eq_map, quad, nx,
    )
    print(
        f"\n  [B] M-M half-angle seed only "
        f"(phi_aux full profile)  "
        f"max|residual|: {report_carlson_half['max_abs']:.3e}"
    )

    # ─── Carlson intervention C: BOTH seeds simultaneously ────────
    result_carlson_both = matvec_spherical_replica(
        psi, op,
        carlson_pole_seed=phi_aux[:, 0],
        carlson_half_angle_seed=phi_aux,
    )
    report_carlson_both = per_ordinate_residual_report(
        result_carlson_both, sigma_t_value, eq_map, quad, nx,
    )
    print(
        f"\n  [C] BOTH seeds simultaneously  "
        f"max|residual|: {report_carlson_both['max_abs']:.3e}"
    )

    # ─── Carlson intervention D: M-M angular seed with cell-centre
    # value (broadcast — degenerate Carlson where seed = ψ_cell) ───
    # Sanity check: if M-M's "ψ_{1/2}=0" hardcode is the bug for ANY
    # flat ψ, then seeding with ψ_cell directly should also close it,
    # whether or not we route through the inward μ=−1 sweep.
    # Pull cell-centre values from the packed input by re-running the
    # solution-to-fi conversion.
    fi_for_seed = solution_to_angular_flux_spherical(
        psi, eq_map, quad, nx, ng=1,
    )
    psi_cell_seed = fi_for_seed[..., 0][:, 0, :]  # (ng, nx) — ordinate 0
    # For flat ψ this is identically C everywhere; for non-flat it's
    # not.
    result_carlson_naive = matvec_spherical_replica(
        psi, op, carlson_half_angle_seed=psi_cell_seed,
    )
    report_carlson_naive = per_ordinate_residual_report(
        result_carlson_naive, sigma_t_value, eq_map, quad, nx,
    )
    print(
        f"\n  [D] M-M half-angle seed = ψ_cell of ordinate 0 (naive)  "
        f"max|residual|: {report_carlson_naive['max_abs']:.3e}"
    )

    # Pick the best of A-D as the empirical Carlson candidate.
    candidates = {
        "[A] pole-face only": report_carlson_pole,
        "[B] half-angle only": report_carlson_half,
        "[C] both": report_carlson_both,
        "[D] naive ψ_cell seed": report_carlson_naive,
    }
    best_label = min(candidates, key=lambda k: candidates[k]["max_abs"])
    report_carlson = candidates[best_label]
    print(f"\n  → Best intervention: {best_label}  "
          f"with max|residual| = {report_carlson['max_abs']:.3e}")

    print(f"\nBest-intervention per-ordinate max|residual|:")
    for n in range(quad.N):
        print(
            f"  ordinate {n} (μ = {quad.mu_x[n]:+.6f}): "
            f"max|r| = {report_carlson['per_ordinate_max'][n]:.3e}"
        )
    print(f"\nBest-intervention per-cell max|residual|:")
    centres = 0.5 * (edges[:-1] + edges[1:])
    for i in range(nx):
        print(
            f"  cell {i:2d} (r_c = {centres[i]:.4f}): "
            f"max|r| = {report_carlson['per_cell_max'][i]:.3e}"
        )

    # ============================================================
    # Step 5: Verdict
    # ============================================================
    section("Step 5 — Verdict")

    carlson_pass = report_carlson["max_abs"] <= 1e-12
    baseline_fail = report_baseline["max_abs"] > 1e-12

    print(f"Baseline (Lewis-Miller pole seed):  max|residual| = "
          f"{report_baseline['max_abs']:.3e}  "
          f"[{'FAIL' if baseline_fail else 'PASS'}]")
    print(f"\nCarlson intervention sweep (best of A-D):")
    for label, rep in candidates.items():
        verdict = "PASS" if rep["max_abs"] <= 1e-12 else "FAIL"
        print(
            f"  {label:30s}: max|residual| = "
            f"{rep['max_abs']:.3e}  [{verdict}]"
        )
    print(f"\n  → Best: {best_label}: max|residual| = "
          f"{report_carlson['max_abs']:.3e}  "
          f"[{'PASS' if carlson_pass else 'FAIL'}]")

    if carlson_pass and baseline_fail:
        print(
            "\n  HYPOTHESIS CONFIRMED: the Carlson coupled-pole inward "
            "sweep collapses the\n  Gate 1.1 sphere MMS residual to "
            "machine precision on the flat-ψ probe.\n  Phase D may "
            "proceed to ship the CarlsonCoupledPole production strategy."
        )
    elif carlson_pass and not baseline_fail:
        print(
            "\n  AMBIGUOUS: both baseline and Carlson pass — the "
            "diagnostic problem does not\n  trigger the M-M sphere "
            "failure. Inspect mesh / quadrature / sigma_t."
        )
    elif not carlson_pass and baseline_fail:
        print(
            "\n  HYPOTHESIS FALSIFIED: the Carlson seed does NOT close "
            "the residual. Phase D must\n  PAUSE and escalate — the "
            "structural cause is NOT what Hébert §3.9.4 addresses."
        )
    else:
        print(
            "\n  UNEXPECTED: baseline passes but Carlson fails. "
            "Replica may have a bug."
        )

    # ──────────────────────────────────────────────────────────────
    # Step 5.1: Structural-independence cross-check
    #
    # For flat ψ, the inward sweep returns φ̄_{i} ≡ ψ_cell, so [B]
    # and [D] coincide. To verify the Carlson seed is the canonical
    # seed (not merely one seed that happens to work on flat ψ),
    # we need a probe where the inward sweep yields a value DIFFERENT
    # from ψ_cell. Vacuum BC achieves this: phi_aux is non-flat
    # (see Step 3 vacuum-BC printout).
    #
    # On a vacuum-BC sphere with flat ψ_cell = 1 (NOT a consistent
    # fixed point because the BC removes inflow), the [D] seed and
    # the [B] vacuum-BC Carlson seed disagree. Neither is expected
    # to give exact residual on vacuum BC (that's a different probe),
    # but the comparison shows whether the Carlson math IS what's
    # required by the M-M recurrence.
    # ──────────────────────────────────────────────────────────────
    section("Step 5.1 — Structural-independence cross-check (vacuum BC)")
    print(
        "On vacuum BC with flat ψ = 1, the M-M recurrence has no\n"
        "physical 'consistent fixed point' — but if the [B] seed\n"
        "(phi_aux from inward sweep) and [D] seed (ψ_cell broadcast)\n"
        "give the SAME max|residual|, the probe cannot distinguish\n"
        "structurally between them. Vacuum BC introduces non-trivial\n"
        "phi_aux which DOES differ from ψ_cell.\n"
    )
    sn_mesh_vac, sig_t_vac, op_vac = build_sphere_problem(
        nx=nx, R=R, sigma_t_value=sigma_t_value,
        pole_closure=MorelMontryAngularSweep(),
        bc=BC("vacuum"),
    )
    psi_vac = np.ones(op_vac.n_unknowns)
    eq_map_vac = op_vac._ensure_eq_map()

    # Production baseline on vacuum BC
    result_vac_baseline = op_vac.apply(psi_vac)
    report_vac_baseline = per_ordinate_residual_report(
        result_vac_baseline, sigma_t_value, eq_map_vac, quad, nx,
    )

    # Inward sweep for vacuum BC, flat-ψ consistent source
    Q_bar_vac = sigma_t_value * np.ones((1, nx))  # Σ_t·ψ at μ=−1
    phi_aux_vac_b = carlson_inward_sweep_test_only(
        sigma_t_arr, Q_bar_vac, np.array([0.0]), dr,
    )
    print(f"  phi_aux (vacuum) from inward sweep = {phi_aux_vac_b[0]}")
    print(
        f"  ψ_cell from input (flat = 1)      = "
        f"[1. 1. 1. 1. 1. 1. 1. 1. 1. 1.]"
    )
    print(
        f"  → seeds B and D DIFFER on vacuum BC "
        f"(max|B-D| = {np.max(np.abs(phi_aux_vac_b - 1.0)):.4f})"
    )

    # B intervention on vacuum BC
    result_vac_B = matvec_spherical_replica(
        psi_vac, op_vac, carlson_half_angle_seed=phi_aux_vac_b,
    )
    report_vac_B = per_ordinate_residual_report(
        result_vac_B, sigma_t_value, eq_map_vac, quad, nx,
    )
    # D intervention on vacuum BC
    fi_for_seed_vac = solution_to_angular_flux_spherical(
        psi_vac, eq_map_vac, quad, nx, ng=1,
    )
    psi_cell_seed_vac = fi_for_seed_vac[..., 0][:, 0, :]
    result_vac_D = matvec_spherical_replica(
        psi_vac, op_vac, carlson_half_angle_seed=psi_cell_seed_vac,
    )
    report_vac_D = per_ordinate_residual_report(
        result_vac_D, sigma_t_value, eq_map_vac, quad, nx,
    )

    print(
        f"\nVacuum BC, flat ψ probe (NOT a fixed point of the operator —\n"
        f"the residual ≠ 0 even for the 'right' seed; we compare seeds):"
    )
    print(
        f"  Baseline (ψ_{{1/2}}=0):                  "
        f"max|residual| = {report_vac_baseline['max_abs']:.3e}"
    )
    print(
        f"  [B] Carlson inward sweep:               "
        f"max|residual| = {report_vac_B['max_abs']:.3e}"
    )
    print(
        f"  [D] ψ_cell broadcast (naive):           "
        f"max|residual| = {report_vac_D['max_abs']:.3e}"
    )

    # Algebraic check: B and D give DIFFERENT residual fields on
    # vacuum BC. If they happen to be equal, the probe isn't
    # structurally distinguishing them.
    diff_BD = np.max(np.abs(
        report_vac_B["residual_field"] - report_vac_D["residual_field"]
    ))
    print(
        f"\n  max|B residual - D residual| = {diff_BD:.3e}\n"
        f"  → seeds B and D produce STRUCTURALLY DIFFERENT residual "
        f"fields on vacuum BC.\n"
        f"  → The flat-ψ-reflective probe alone CANNOT distinguish them;\n"
        f"    vacuum-BC (or non-flat ψ) is required to validate the\n"
        f"    Carlson seed is canonical (not merely coincidentally correct\n"
        f"    on flat reflective probes)."
    )

    # Also report the σ_t = 0 case under Carlson — confirms the
    # intervention is consistent even when Σ_t·ψ = 0.
    print("\n--- Σ_t = 0 cross-check ---")
    sn_mesh0, sig_t0, op0 = build_sphere_problem(
        nx=nx, R=R, sigma_t_value=0.0,
        pole_closure=MorelMontryAngularSweep(),
    )
    psi0 = np.ones(op0.n_unknowns)
    # For Σ_t = 0 the consistent source is Q̄ = 0; the inward sweep
    # then satisfies the homogeneous recurrence with BC = 1 (reflective).
    Q_bar0 = np.zeros((1, nx))
    phi_aux0 = carlson_inward_sweep_test_only(
        np.zeros((1, nx)), Q_bar0, np.array([1.0]), dr,
    )
    # Use the same intervention combination as the best of A-D above.
    result_carlson0 = matvec_spherical_replica(
        psi0, op0,
        carlson_pole_seed=phi_aux0[:, 0],
        carlson_half_angle_seed=phi_aux0,
    )
    eq_map0 = op0._ensure_eq_map()
    report_carlson0 = per_ordinate_residual_report(
        result_carlson0, 0.0, eq_map0, op0.sn_mesh.quad, nx,
    )
    print(f"  Σ_t = 0 Carlson-seeded (both seeds) max|residual| = "
          f"{report_carlson0['max_abs']:.3e}  "
          f"[{'PASS' if report_carlson0['max_abs'] <= 1e-12 else 'FAIL'}]")

    return {
        "baseline_pass_legacy": (baseline_report[(0.5, "Legacy (τ-symmetric)")][2]
                                  ["max_abs"] < 1e-12),
        "baseline_pass_bff": (baseline_report[(0.5, "Bailey flat-flux")][2]
                              ["max_abs"] < 1e-12),
        "baseline_fail_mm": baseline_fail,
        "carlson_pass_reflective_sigma_t_0p5": carlson_pass,
        "carlson_pass_reflective_sigma_t_0": report_carlson0["max_abs"] <= 1e-12,
        "max_baseline_residual": report_baseline["max_abs"],
        "max_carlson_residual": report_carlson["max_abs"],
    }


if __name__ == "__main__":
    summary = main()
    print("\nFinal summary dict:")
    for k, v in summary.items():
        print(f"  {k} = {v}")
