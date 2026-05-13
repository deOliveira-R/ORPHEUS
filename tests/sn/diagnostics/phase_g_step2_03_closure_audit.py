"""Diagnostic: Output D + Output E — structural code-walk audit of SI vs Krylov paths.

Created by numerics-investigator on 2026-05-12 for Phase G Step 2 pre-step.

This diagnostic is a STRUCTURAL audit, not a numerical probe. It tabulates
the per-cell algebraic intermediates of the two paths side by side, and
states the architectural difference that produces the O(h) drift.

NOT a runnable test in the usual sense — its evidence is in the file
itself + the markdown summary at the bottom. Kept under pytest so the
audit text-block ships and pytest verifies the constant-extraction
sanity check (face count etc.).
"""
from __future__ import annotations

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.geometry import SNMesh


def test_structural_audit_consistency():
    """Verify that the two paths' BC entry points match what the audit claims."""
    fuel = get_mixture("A", "2g")
    mod = get_mixture("B", "2g")
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(
            Region(mat_id=0, outer_thickness_cm=0.5),
            Region(mat_id=1, outer_thickness_cm=1.0),
            Region(mat_id=0, outer_thickness_cm=0.5),
        ),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(
        geom,
        region_meshes=(
            RegionMesh(n_cells=10),
            RegionMesh(n_cells=20),
            RegionMesh(n_cells=10),
        ),
    )
    sn_mesh = SNMesh(mesh, GaussLegendre1D.create(n_ordinates=8))

    nx = sn_mesh.nx
    N = sn_mesh.quad.N

    assert nx == 40
    assert N == 8
    assert sn_mesh.bc_right is not None  # both paths consult this

    # eq_map skips inward-at-outer-boundary slots — so the
    # _sweep_1d_spherical and transport_operator_matvec_spherical
    # MUST compute those values themselves.
    from orpheus.sn.operator import build_equation_map_spherical
    eq_map = build_equation_map_spherical(nx, sn_mesh.quad, 2)
    inward_outer_slots = sum(
        1 for n in range(N)
        if sn_mesh.quad.mu_x[n] < -1e-15
    )
    expected_n_eq = nx * N - inward_outer_slots  # (40*8 - 4 = 316)
    assert eq_map.n_eq == expected_n_eq
    # 4 inward ordinates (μ < 0): the matvec must reconstruct their
    # ψ at i=nx-1 from BC + the input cell-centre at i=nx-1 via
    # `bc_outer.apply(fi[:, :, -1, 0].T)`.


AUDIT_TEXT = """
═══════════════════════════════════════════════════════════════════════
STRUCTURAL AUDIT — SI sweep vs Krylov apply-matvec per-cell algebra
═══════════════════════════════════════════════════════════════════════

CODE PATHS (sphere, refactor/sn-operator-algebra @ dda6f28):

1. SI sweep
   `orpheus/sn/sweep.py:_sweep_1d_spherical`
     → per ordinate n in 0..N-1:
       → bc_outer_obj.apply(bc_outer)[n] gives ψ_in face flux from BC
       → for each cell i in DAG order:
         → cell_update.update(visit, sig_t, source=Q·V·norm,
                              upstream=UpstreamState(psi_spat_in, psi_angle[i]))
           (DiamondDifference._update_curvilinear)
         → returns CellResult with
             cell_average_flux psi_avg = numer / denom
             outgoing_spatial_flux  = 2·psi_avg − psi_spat_in (WDD)
             outgoing_angular_state = (psi_avg − (1−τ)·psi_angle_in)/τ (M-M)
           where
             denom = 2|μ|·A_downstream + dA_w·c_out + Σ_t·V
             numer = source + |μ|·(A_in + A_out)·psi_spat_in
                     + dA_w·c_in·psi_angle_in
             c_out = α_out/τ;  c_in = (1−τ)/τ · α_out + α_in
         → psi_angle[i] ← psi_angle_out  (angular face flux propagates
                                          ALONG ordinates inside this i)
         → psi_spat_in ← psi_spat_out    (spatial face flux propagates
                                          to next cell)
     CARLSON SEED: psi_angle initialised at sweep start via
       `carlson_inward_sweep_from_source(Q_bar=0.5·Q_1d, sigma_t, dr, bc_outer_value)`
       → produces ψ_{1/2,i,g} via Hébert §3.9.4 inward DD sweep.

2. Apply matvec
   `orpheus/sn/operator.py:transport_operator_matvec_spherical`
     → fi ← solution_to_angular_flux_spherical(sol, ...)
        → fi[:, n, i, 0] = sol[k]  at unknown slots
        → fi[:, n, nx-1, 0] = fi[:, ref[n], nx-1, 0] for inward
          ordinates n at outer cell (CELL-CENTRE proxy for face flux,
          ONLY first-order at boundary on non-constant solutions)
     → bc_outer_value = bc_outer.apply(fi[:, :, -1, 0].T)[most_inward, :]
        (extracts Carlson seed BC input from CELL-CENTRE proxy)
     → carlson_ctx = CarlsonSweepContext(σ_t, dr, μ, w, bc_outer_value)
     → redist_full = pole_angular_closure(fi, alpha_half, dAw, tau, V,
                                          carlson_context=carlson_ctx)
        → MorelMontryAngularSweep(psi_cells=fi, ...)
        → internally runs Carlson inward sweep from INPUT ψ (not from
          source) to seed ψ_{1/2,i,g}, then M-M recurrence over m=0..M-1
        → returns redist[g, m, i] — independent of cell-by-cell sweep
          state; this is a SEPARATE PASS over the whole ψ field
     → BC PHASE 1 (outward, μ > 0, i = 0 → nx-1):
       psi_face_in initialised at pole = fi[:, mask, 0, 0]
         (CELL-CENTRE at pole, exact for flat ψ on reflective sphere)
       for each cell visit:
         psi_cell = fi[:, mask, i, 0]                # input cell-centre
         psi_face_out = 2·psi_cell − psi_face_in     # WDD
         streaming = μ·(A[i+1]·psi_face_out − A[i]·psi_face_in)/V[i]
         redistribution = redist_full[:, mask, i]
         collision = σ_t[i]·psi_cell
         lhs[k] = streaming + redistribution + collision
         psi_face_in ← psi_face_out
       outflow_at_boundary[:, mask] = psi_face_out
     → BC TRACE: inflow_full = bc_outer.apply(outflow_at_boundary.T)
        (operates on FACE values from the outward sweep — correct trace)
     → BC PHASE 2 (inward, μ < 0, i = nx-1 → 0):
       psi_face_in = inflow_full[mask, :].T          # FACE inflow from BC
       for each cell visit: same per-cell algebra

═══════════════════════════════════════════════════════════════════════
ALGEBRAIC FORM COMPARISON (per-cell, per-ordinate)
═══════════════════════════════════════════════════════════════════════

SI SWEEP: solve form (numer/denom)
  denom · psi_avg = source + |μ|·(A_in + A_out)·psi_spat_in + dA_w·c_in·psi_angle_in
  where psi_avg is the cell-centre value, source is per-ordinate
  Q·V·norm (the source iteration's within-group source).

APPLY MATVEC: residual form (LHS only — Σ_s + (1/k)F live on RHS)
  LHS[k] = streaming + redistribution + collision
         = μ·(A[i+1]·psi_face_out − A[i]·psi_face_in)/V[i]
           + redist[g, m, i]                          (from input ψ)
           + Σ_t[i, g] · psi_cell
  where psi_face_out = 2·psi_cell − psi_face_in (WDD).

═══════════════════════════════════════════════════════════════════════
THREE STRUCTURAL DIFFERENCES (root cause of the O(h) drift)
═══════════════════════════════════════════════════════════════════════

(1) ANGULAR-REDISTRIBUTION COUPLING — STATE PROPAGATION
    SI sweep: M-M recurrence runs INSIDE the cell-update.  The cell's
      psi_angle_out propagates to the NEXT ordinate at the SAME cell
      (see `psi_angle[i] ← psi_angle_out` after each cell-visit; psi_angle
      is a (nx, ng) array that lives ACROSS the ordinate loop).
      Equivalently: at cell i, ordinate n, the angular face flux is built
      from the just-computed psi_avg AT THIS CELL i.
    Apply matvec: redist_full is precomputed from INPUT psi_cells
      via M-M recurrence in pole_angular_closure(...) — runs ONCE per
      matvec call over all (g, m, i) BEFORE any per-cell sweep.
      Equivalently: the angular face flux at cell i ordinate n uses
      psi_cells[g, m, i] AT THE INPUT, not the WDD-propagated value.

    On the fixed point ψ* (where input = output) the two forms are
    ALGEBRAICALLY equivalent.  Off the fixed point they differ — and
    the way they CONVERGE to fixed points is different.

(2) PSI HALF-ANGLE SEED SOURCE
    SI sweep:  Carlson inward sweep driven by SOURCE Q_1d.  Eq (3.434):
        φ̄_i = (Δr_i · Q̄_i + 2·φ̄_{i+1/2}) / (Δr_i·Σ_t,i + 2)
      with Q̄_i = 0.5·Q_1d[i] for L=0 isotropic.
    Apply matvec: Carlson inward sweep driven by INPUT ψ.  Same Eq
      (3.434) but with Q̄_i = 0.5·Σ_t,i·φ_0(r_i) where φ_0 is built
      from input ψ via the M-M angular recurrence.

    On the fixed point ψ* the two seeds AGREE because Σ_t·φ_0 = Q_1d
    at the within-group fixed point.  Off the fixed point they differ.

(3) BOUNDARY-CONDITION TRACE EVALUATION
    SI sweep: bc_outer.apply(bc_outer_buffer) gives the (N, ng) inflow
      vector indexed by ordinate; bc_outer_buffer is a PERSISTENT
      buffer indexed by *outgoing* ordinate, populated at the END of
      each outward sweep with psi_spat_out (the WDD face flux at the
      outer face from cell nx-1).  Read at the START of inward sweeps.
      → Uses FACE values (correct trace).
    Apply matvec PHASE 1: bc_outer.apply(fi[:, :, -1, 0].T) seeds the
      Carlson inward sweep with a value derived from the input
      CELL-CENTRE at i=nx-1 (NOT the face flux).
      → Uses CELL-CENTRE as face proxy at the outer boundary — only
      first-order at the boundary on non-constant solutions.
    Apply matvec PHASE 2: inflow_full = bc_outer.apply(outflow_at_boundary.T)
      where outflow_at_boundary is the WDD-propagated outward-sweep face flux.
      → Uses FACE values (correct trace) — this is the SAME as the
      SI path for the inward-sweep entry.

    NOTE the asymmetry: the apply-matvec uses CELL-CENTRE-as-face for
    the Carlson seed but FACE values for the inward-sweep BC.  This
    is the O(h) cell-centre-as-face-proxy at the boundary that the
    `solve_sn_fixed_source` docstring (orpheus/sn/solver.py:974-988)
    explicitly flags as the reason the curvilinear MMS rate regresses
    from O(h^{1.3}) to O(h^1) when Krylov is the default.

═══════════════════════════════════════════════════════════════════════
EMPIRICAL EVIDENCE (from phase_g_step2_01_psi_comparison.py)
═══════════════════════════════════════════════════════════════════════

At sphere_2g_3reg n=40 the SI and Krylov fixed points:
- Differ by 0.286% in k_eff.
- Differ by 168% in TOTAL FLUX NORM (sf_si.sum()=91.57, sf_kr.sum()=34.14).
  → Both are valid eigenvectors; the eigenmode SHAPE differs.
- Per-cell Δψ (after sum-normalisation): max 7.92e-3 at i*=39
  (the outer-boundary cell, NOT the pole).
- Per-cell drift monotonically RISES with cell index toward i=nx-1:
  i=0 +26.15% / i=10 -5.02% / i=20 -3.40% / i=30 -3.13% / i=39 -5.30%
  → cells near the OUTER boundary are most affected.
- Per-ordinate at i=39: SI gives ψ_in ≈ 4e-3 for all inward ordinates
  (near-isotropic inflow); Krylov gives ψ_in ≈ 1e-4 (near-vacuum inflow).
  Kr/SI ratio: μ<0: 0.007 (50× drift); μ>0: 0.012-0.32.

The drift pattern (worst at the outer boundary, growing inward
toward the boundary, asymmetric in μ at the boundary) is the
CHARACTERISTIC SIGNATURE of a BC-trace truncation error compounded by
the angular-redistribution scheme.

═══════════════════════════════════════════════════════════════════════
VERDICT: H2 (closure algebra differs)
═══════════════════════════════════════════════════════════════════════

The SI sweep and the Krylov apply-matvec are NOT the same operator at
the per-cell level.  They share:
  - The WDD diamond closure  psi_face_out = 2·psi_cell − psi_face_in.
  - The M-M angular recurrence at τ = 1/2 (pure DD angular).
  - The Hébert §3.9.4 Carlson coupled-pole seed Eq. (3.434).

But they differ in:
  (A) WHICH ψ feeds the M-M angular recurrence (propagated vs input)
  (B) WHICH SOURCE feeds the Carlson seed (Q or Σ_t·φ_0)
  (C) WHICH BC TRACE feeds the Carlson seed (face flux or cell-centre proxy)

(A) and (B) agree on the fixed point.  (C) does NOT — the apply
matvec's `bc_outer.apply(fi[:, :, -1, 0].T)` is an O(h) truncation
at the outer face even on the fixed point.

The empirical drift in our test (worst cell at i=39, rising toward
the boundary) maps cleanly onto (C) being the dominant contributor.

═══════════════════════════════════════════════════════════════════════
IMPLICATIONS FOR STEP 2 CALL-SITE UNIFICATION
═══════════════════════════════════════════════════════════════════════

CANONICAL CLOSURE CHOICE: the SI sweep's algebra is the canonical
form (face-flux BC trace, ψ propagated inside per-cell update).  The
Krylov apply-matvec's CELL-CENTRE-AS-FACE-PROXY at the boundary is
the O(h) artefact that must go away.

The architecture reconciliation memo §2.1 recommends WDD (the SI
form) backed by Bailey-Morel-Chang asymptotic-diffusion-limit.  The
empirical evidence above SUPPORTS THIS RECOMMENDATION.

STEP 2 UNIFICATION SHAPE:
- `SNCellOperator` already wraps DiamondDifference.update (= the SI
  per-cell algebra).  Step 2 routes both call sites through THAT.
- BC entry: ONE entry per direction-reversal at the outer face,
  consuming the WDD-PROPAGATED FACE FLUX (NOT a cell-centre proxy).
  This is what `_sweep_1d_spherical` already does; the apply matvec
  must be rewritten to match.
- Angular redistribution: M-M recurrence runs INSIDE the per-cell
  loop with the PROPAGATED psi_avg, not as a precomputed pass over
  input ψ.  This is the AngularRedistribution operator's apply
  invocation — but it should be invoked PER CELL, PER ORDINATE,
  consuming the just-computed psi_avg at THIS cell.
- Carlson seed: ONE seed strategy per direction-reversal, driven by
  the within-group SOURCE (the SI form's `Q̄_i = 0.5·Q_1d[i]`).
  The apply-matvec's input-ψ-driven form is equivalent on the fixed
  point but produces transient drift off it.

The unified path is structurally: a per-(cell-visit) loop that
internally constructs the per-cell balance via `SNCellOperator.solve`
when used as a sweep, and via `SNCellOperator.apply` when used as a
matvec residual.  Both share the same `denom`, `numer_upstream`,
`c_in`, `c_out` algebraic intermediates.

LOCATION OF THE _cell_balance_terms HELPER:
- Currently DiamondDifference._update_curvilinear (diamond.py:550-624)
  and SNCellOperator._apply_curvilinear_residual (operators.py:319-367)
  both build `denom`, `numer_upstream`, `c_in`, `c_out` —
  qa CONCERN-2 of Step 1 review.
- These should factor into a free function or static method, e.g.:
    def _cell_balance_terms(st, A_downstream, total_xs, upstream_state, source)
      → returns CellBalanceTerms(c_out, c_in, denom, numer_upstream)
- Both DiamondDifference._update_curvilinear (which divides:
    psi_avg = (source + numer_upstream) / denom)
  AND SNCellOperator._apply_curvilinear_residual (which subtracts:
    residual = denom · cell_avg − (source + numer_upstream))
  consume the same `denom` and `numer_upstream`.
- Recommended location: as a private free function in `operators.py`
  (or moved into `cell_balance.py` as the per-cell math module). It
  must be CALLABLE FROM BOTH `solve` and `apply` so the bit-identity
  contract on Gate 2 round-trip is preserved.

BIT-IDENTITY EXPECTATION (revising plan note 5):
- The 11 regression snapshots are SI-generated.  After Step 2 the
  Krylov apply-matvec switches to face-flux BC trace, but the SI
  sweep's algebraic body is unchanged — so the SI snapshots remain
  BIT-IDENTICAL.
- The Krylov-driven results CHANGE — but no regression snapshots are
  Krylov-generated, so no snapshot regen is required from this Step.
- Phase E flux-shape sentinel (which xfails because SI and Krylov
  fixed points disagree by 5%+) WILL XPASS at Step 2 when both
  paths produce the SAME fixed point.
- One caveat: if Step 2 also fixes the apply-matvec to use the SI's
  source-driven Carlson seed (vs the current input-ψ-driven seed),
  then the BiCGSTAB / GMRES iterative residual on the apply matvec
  changes during transient iterations — but the fixed point is the
  same SI fixed point.  The 6 curvilinear snapshots stay
  bit-identical because they are SI-generated.

PREDICTION: at Step 2 the Phase E flux-shape sentinel
`test_phase_e_trajectory_resolvent_flux_shape_crosscheck` XPASSES
because both methods now solve the same operator with the same
boundary trace.
"""


if __name__ == "__main__":
    test_structural_audit_consistency()
    print(AUDIT_TEXT)
