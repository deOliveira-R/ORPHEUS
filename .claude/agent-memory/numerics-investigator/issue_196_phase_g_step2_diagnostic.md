---
name: issue-196-phase-g-step2-diagnostic
description: Per-cell SI-vs-Krylov drift characterisation for Phase G Step 2. VERDICT H2 (closure algebra differs at the operator level) with very-high confidence. L0 streaming-equilibrium test reveals SI has structural 22%-at-pole error on homogeneous reflective sphere fixed-source — does not improve with refinement; Krylov gives machine-zero. ERR-026 manifestation #6 still ACTIVE post-Phase-F at the L0 level. Recommendation: Step 2 unifies on the APPLY-MATVEC closure (face-flux BC trace + cell-centre-as-pole-face seed), not on the SI sweep closure. Diagnostic memo includes 6-script suite + structural code-walk audit.
metadata:
  type: project
  branch: refactor/sn-operator-algebra
  step: Phase G Step 2 pre-step
  date: 2026-05-12
---

# Phase G Step 2 pre-step — Per-cell SI-vs-Krylov drift diagnostic

**Branch**: `refactor/sn-operator-algebra` 2026-05-12 (Step 1 tip `dda6f28`).

**Plan**: `.claude/plans/issue_196_phase_g_four_operator_unification.md` Step 2.

**Phase F predecessor**: `.claude/agent-memory/method-implementer/issue_168_phase_f_closeout.md`.

## Headline — H2 with very-high confidence; Krylov IS the canonical closure

**Hypothesis H2 confirmed**: the SI sweep (`_sweep_1d_spherical`) and the
Krylov apply-matvec (`transport_operator_matvec_spherical`) are NOT the
same operator at the per-cell level. The empirical evidence is decisive:

**The L0 streaming-equilibrium test on a homogeneous reflective sphere
with isotropic Q = 1 (analytical answer: ψ = Q/Σ_t·1/Σw = 5.0 for every
ordinate at every cell):**

| n_cells | SI max rel error | SI worst cell | Kr max rel error |
|---------|------------------|---------------|------------------|
| 20      | **22.40%**       | **i = 0 (pole)** | 1.93e-11      |
| 40      | **21.96%**       | **i = 0 (pole)** | 1.93e-11      |
| 80      | **21.55%**       | **i = 0 (pole)** | 1.93e-11      |

**The SI error is NOT a discretisation artefact** — it does not improve
with refinement; it plateaus at ~22% structural bias at the pole. The
Krylov apply-matvec hits the analytical answer to machine precision.

**This is an active L0 V&V failure.** The Phase F Carlson seed fix
improved the bug (reduced from the pre-Phase-F divergent ratio of 0.522
to a stable but still-wrong 0.778 ratio) but did NOT close it. ERR-026
manifestation #6 remains ACTIVE at the L0 verification level.

## Empirical outputs

### Output A — per-cell-per-ordinate ψ comparison (k-eigenvalue)

Script: `tests/sn/diagnostics/phase_g_step2_01_psi_comparison.py`.
Snapshot: `/tmp/phase_g_step2_psi.npz`.

Test problem: sphere_2g_3reg n=40, GL-8, reflective BC. Phase F closeout
numbers reproduced exactly:
- k_eff(SI) = 1.38069561, k_eff(Kr) = 1.38464040, diff = 0.286%.
- sf[0]/sf[1] (SI) = 0.7776, (Kr) = 1.0288.

After sum-normalising both fluxes, the worst-cell ψ drift is at:
- **i\* = 39 (outer boundary cell)**, max |Δψ| = 7.92e-3.
- Per-ordinate breakdown at i=39 g=0: SI gives ψ_in ≈ 4-7e-3 for ALL
  inward ordinates (near-isotropic inflow); Krylov gives ψ_in ≈ 1e-4
  (near-vacuum inflow). Kr/SI ratio 0.007 for inward, 0.012-0.32 for
  outward.
- Per-cell |Δψ|_∞ is monotonically RISING from interior toward i=39:
  i=10 → 7.59e-3, i=20 → 5.64e-3, i=30 → 4.97e-3, i=39 → 7.92e-3.

Unnormalised: total flux norms differ by 168% (sf_si.sum() = 91.57,
sf_kr.sum() = 34.14). The two eigenvectors disagree fundamentally —
their k_eff agreement to 0.3% is incidental.

### Output B — apply-matvec residual at each fixed point

Script: `tests/sn/diagnostics/phase_g_step2_02_sncell_residual.py`.

Not load-bearing — the apply-matvec consumes ψ via
`solution_to_angular_flux_spherical` which fills inward-at-outer slots
with the reflected partner's cell-centre value, so a direct residual
comparison is contaminated by the matvec's own input-extension logic.
The fixed-source attribution test below replaces this.

### Output C — fixed-source SI vs Krylov at MATCHED Q (DECISIVE)

Script: `tests/sn/diagnostics/phase_g_step2_04_fixed_source.py`.
Snapshot: `/tmp/phase_g_step2_fixed_source.npz`.

**This is the cleanest H1 vs H2 discriminator** — Q is externally
imposed and identical for both solvers, so any ψ drift is purely from
the within-group closure form, NOT source feedback.

sphere_2g_3reg n=40, reflective BC, isotropic Q=1:
- Worst cell drift: i\* = **0 (pole)** with max |Δψ| = 10.6 (63%
  relative).
- Per-cell relative drift (Kr - SI / SI on sf g=0):
  - i = 0:  **+36.8%** (pole)
  - i = 10: −0.44%
  - i = 20: −1.23%
  - i = 30: −0.71%
  - i = 39: −1.71% (outer)
- The drift is sharply localised at the pole; interior is within 2%.

**Per-ordinate at i = 0 (pole), g = 0, fixed-source isotropic Q = 1**:

| n | μ_n     | ψ_si    | ψ_kr    | Δψ      |
|---|---------|---------|---------|---------|
| 0 | −0.9603 | 2.58    | 3.13    | +0.55   |
| 1 | −0.7967 | 2.70    | 3.15    | +0.45   |
| 2 | −0.5255 | 2.78    | 3.15    | +0.37   |
| 3 | −0.1834 | 2.82    | 3.15    | +0.34   |
| 4 | +0.1834 | 2.64    | 3.16    | +0.52   |
| 5 | +0.5255 | 2.00    | 3.16    | +1.17   |
| 6 | +0.7967 | 1.08    | 3.17    | +2.10   |
| 7 | +0.9603 | 0.38    | 3.19    | +2.81   |

**SI gives strongly anisotropic ψ at the pole** (cv ≈ 1.0, max/min ≈ 7×)
while **Krylov gives near-isotropic ψ** (cv ≈ 0.005, max/min ≈ 1.018).
Pomraning 1989 structural-singularity result predicts NEAR-ISOTROPY at
r = 0 for the sphere — Krylov agrees.

### Output D — structural code-walk audit

Script + AUDIT_TEXT: `tests/sn/diagnostics/phase_g_step2_03_closure_audit.py`.

Three structural differences identified:

**(1) Angular-redistribution coupling — state propagation**
- SI sweep: M-M recurrence runs INSIDE `DiamondDifference._update_curvilinear`.
  `psi_angle[i]` propagates `psi_angle_out` ALONG the ordinate axis. The
  cell's `psi_avg` (just computed at THIS cell, THIS ordinate) drives
  the recurrence into the NEXT ordinate at the SAME cell.
- Apply matvec: `redist_full` precomputed from INPUT `psi_cells` via
  `MorelMontryAngularSweep(psi_cells, ...)` — ONE pass over all (g, m, i)
  BEFORE any per-cell sweep. The angular face flux at cell i ordinate n
  uses `psi_cells[g, m, i]` at the INPUT, not the WDD-propagated value.
- At the fixed point these are algebraically equivalent. Off the fixed
  point they aren't, and the iteration converges differently.

**(2) Carlson seed source**
- SI sweep: source-driven, `Q̄_i = 0.5 · Q_1d[i]` for L=0 isotropic.
- Apply matvec: input-ψ-driven, `Q̄_i = 0.5 · Σ_t,i · φ_0(r_i)` where
  φ_0 is built from input ψ.
- On the within-group fixed point `Σ_t · φ_0 = Q_1d`, so the seeds AGREE.
  Off the fixed point they differ.

**(3) Boundary-condition trace evaluation — THE DOMINANT DIFFERENCE**
- SI sweep: `bc_outer_obj.apply(bc_outer)` where `bc_outer` is a
  persistent buffer populated with `psi_spat_out` (the WDD face flux at
  the outer face) at the END of each outward sweep. Uses **face values**.
- Apply matvec PHASE 1 (Carlson seed): `bc_outer.apply(fi[:, :, -1, 0].T)`
  — operates on **cell-centre values** at the outer cell. This is only
  first-order at the boundary on non-constant solutions.
- Apply matvec PHASE 2 (inward sweep entry):
  `bc_outer.apply(outflow_at_boundary.T)` where `outflow_at_boundary` is
  the WDD-propagated outward face flux. Uses **face values** — correct.

The asymmetry inside the apply-matvec — cell-centre-as-face for the
Carlson seed but face values for the inward sweep — is the documented
O(h) artefact at the boundary face (orpheus/sn/solver.py:974-988
docstring of `solve_sn_fixed_source`). This is one part of the closure
discrepancy.

### Output E — BC application timing

Already covered by Output D. Both paths apply BC inside the within-group
iteration, not at the outer (multi-group / k-eigenvalue) iteration level.
The DIFFERENT BC traces (face flux vs cell-centre) are the actual
structural divergence.

## L0 streaming-equilibrium failure — the smoking gun

Script: `tests/sn/diagnostics/phase_g_step2_05_homogeneous.py`.

Test problem: **homogeneous mixture B** (Σ_t = 2, Σ_s = 1.9, no
fission) **sphere of 2 cm**, **reflective BC**, **isotropic Q = 1**,
**GL-8**.

Analytical solution: with `c_eff = (Σ_s + νΣ_f) / Σ_t = 0.95`, infinite
reflective sphere has `φ = Q / (Σ_t · (1 − c_eff)) = 10.0` everywhere.

Pomraning 1989: per-ordinate isotropy at r = 0. So
`ψ_n(r=0) = φ(0) / Σw = 10.0 / 2.0 = 5.0` for every ordinate.

Result table:

| n_cells | SI max abs error | SI worst cell | SI rel error % | Kr max abs error | Kr rel error % |
|---------|------------------|---------------|----------------|------------------|----------------|
| 20      | 2.24             | i=0 (pole)    | **22.40%**     | 1.93e-10         | 1.93e-9%       |
| 40      | 2.20             | i=0 (pole)    | **21.96%**     | 1.93e-10         | 1.93e-9%       |
| 80      | 2.16             | i=0 (pole)    | **21.55%**     | 1.93e-10         | 1.93e-9%       |

SI per-ordinate at the pole (n=40), expected ψ = 5.0:

```
  ord  mu       psi_si    error
  0  -0.9603    4.5736    -8.5%
  1  -0.7967    4.7208    -5.6%
  2  -0.5255    4.7000    -6.0%
  3  -0.1834    4.6168    -7.7%
  4  +0.1834    4.2190    -15.6%
  5  +0.5255    3.2422    -35.2%
  6  +0.7967    2.0733    -58.5%
  7  +0.9603    1.3278    -73.4%
```

The SI per-ordinate ψ at the pole shows a clear **outward-bias-toward-zero**
pattern: high-|μ| outward ordinates approach zero (1.33 at μ = +0.96),
while inward ordinates stay near 4.5-4.7. **Pomraning isotropy is
severely violated.** Krylov gives ψ = 5.0 exactly for every ordinate.

This is a regression of **L0-SN-001** (streaming-equilibrium) in the
curvilinear SI sweep, surviving Phase F's Carlson seed fix.

## Canonical closure decision

**Step 2 MUST unify on the APPLY-MATVEC closure form.** The SI sweep's
current closure carries an active L0-failing bug at the pole.

The architecture reconciliation memo §2.1 recommends WDD (the SI form)
backed by Bailey-Morel-Chang asymptotic-diffusion-limit. **This
diagnostic empirically refutes that recommendation in its current SI
implementation**: the SI's WDD-with-M-M-recurrence-inside-cell-update
form is broken at the pole. The Krylov apply-matvec's WDD-with-
M-M-recurrence-as-precomputed-redist-pass form is correct.

The actual canonical form is the apply-matvec's structure:

1. M-M angular recurrence runs as a SEPARATE PASS over the ψ field
   (BEFORE per-cell sweep), driven by the INPUT ψ. This gives the
   redistribution term per (g, m, i) up front.
2. WDD diamond closure `psi_face_out = 2·psi_cell − psi_face_in`
   propagates the spatial face flux cell-by-cell.
3. BC trace law applied ONCE per direction reversal at the outer face,
   consuming the WDD-propagated FACE flux (not cell-centre proxy).
4. Carlson coupled-pole seed evaluated at the start of the inward sweep
   from the input ψ field (or from the source on the fixed point —
   equivalent).

What needs to BE FIXED in the apply-matvec before Step 2 unifies onto it:
- The Carlson seed's `bc_outer.apply(fi[:, :, -1, 0].T)` cell-centre
  proxy → replace with the WDD-propagated face flux from the (just-
  computed) outward sweep. This eliminates the O(h) boundary truncation
  the `solve_sn_fixed_source` docstring already flags.

## Implications for Step 2 unification

### One canonical SNCellOperator, called per cell-visit

The `SNCellOperator.apply` method as currently written
(`orpheus/sn/spatial/operators.py:319-367`) computes
`denom · cell_avg − (source + numer_upstream)`. This is the SI sweep's
per-cell residual form. Step 2 routes both call sites through it.

**BUT** the SI sweep's algebra `denom · psi_avg = source + numer_upstream`
inherently encodes the M-M recurrence INSIDE the per-cell solve (via
`psi_angle_in` → `psi_angle_out` propagation through ordinates at the
SAME cell). The apply-matvec separates these into two passes
(redist + per-cell streaming + collision).

**Two implementation paths for Step 2**:

**(a) Recast as cell-balance + angular-redist composition**:
- `SNCellOperator.apply(cell_avg, *, visit, total_xs, upstream_state,
  source) → residual` — uses the SI's `denom · cell_avg − (source +
  numer_upstream)` form WITHOUT the angular redistribution inside.
- A separate `AngularRedistribution.apply(...)` operator runs the M-M
  recurrence over the cell field.
- `(L + C).apply(ψ) = streaming(ψ) + R(ψ) + Σ_t · ψ` composes these.
- This is closer to the apply-matvec's structure and is L0-correct.

**(b) Bind the angular recurrence into the per-cell solve**:
- Keep the SI's "M-M recurrence inside cell-update" form.
- Adds a `psi_angle_state` to the upstream state that persists ACROSS
  ordinates at the same cell. The cell-visit DAG must order ordinates
  correctly for this state propagation.
- This matches DiamondDifference._update_curvilinear's current
  contract.

**(a) is RECOMMENDED.** The L0 streaming-equilibrium failure of (b)
demonstrates it doesn't give the analytical answer even at the
fixed point. The apply-matvec's structure (a) is mathematically the
canonical Hébert §3.9.4 form.

### Location of `_cell_balance_terms` helper

Currently `DiamondDifference._update_curvilinear` (diamond.py:580-624)
and `SNCellOperator._apply_curvilinear_residual` (operators.py:319-367)
build the same intermediates. The qa CONCERN-2 fix: factor into a
named helper.

**Recommended**: a private free function in
`orpheus/sn/spatial/cell_balance.py` (new module) returning a
`CellBalanceTerms` frozen dataclass:

```python
@dataclass(frozen=True, slots=True)
class CellBalanceTerms:
    c_in: float | np.ndarray
    c_out: float | np.ndarray
    denom: np.ndarray            # (ng,)
    numer_upstream: np.ndarray   # (ng,) — excluding source

def cell_balance_terms(
    st: StreamingTerms,
    A_downstream: float,
    total_xs: np.ndarray,
    upstream_state: UpstreamState,
) -> CellBalanceTerms:
    ...
```

Then:
- `DiamondDifference._update_curvilinear`: `psi_avg = (source + terms.numer_upstream) / terms.denom`.
- `SNCellOperator._apply_curvilinear_residual`: `terms.denom * cell_avg − (source + terms.numer_upstream)`.

Both share the algebra; only one place to maintain.

### BC promotion (Step 4) — partial overlap with Step 2

The BC face-vs-cell-centre asymmetry in the apply-matvec is structurally
fixed by Step 4 (BoundaryRealizer as first-class operator). But Step 2's
correctness depends on it. **Recommendation**: pull forward the Carlson
seed's BC-trace fix into Step 2 — replace
`bc_outer.apply(fi[:, :, -1, 0].T)` with
`bc_outer.apply(outflow_face_from_outward_sweep)`. The remaining
infrastructure work (BoundaryRealizer Protocol) lives in Step 4.

### Bit-identity expectation revised

The plan's note 5 says "Bit-identity break on snapshots is expected at
Step 2". **Confirmed**: Step 2 WILL break the 6 curvilinear regression
snapshots (sphere_2g_homogeneous_dd_n20, sphere_2g_3reg_dd_n40, etc.)
because:

1. Those snapshots were generated under the SI default (with the
   broken pole closure).
2. The Step 2 unification flips the default to the Krylov-correct
   closure.
3. The Krylov fixed point differs from the SI fixed point by ~6-30% at
   the pole.

Per `vv-principles` §"Bit-identity vs principled-equivalence":
- Principled new formulation: YES — the L0 streaming-equilibrium
  failure of the SI closure is the principled reason.
- Structurally-independent reference: YES — Pomraning 1989 isotropy +
  analytical Q/Σ_a (1-D homogeneous reflective sphere has zero leakage,
  so φ = Q/Σ_a exactly at every point) provides closed-form ground.
  Variant α's integral form is a second structurally-independent
  reference.
- FP-non-associativity bound: NO — this is a structural change, not
  FP noise. Snapshot REGEN is required, not contract relaxation.

The 6 curvilinear snapshots must be REGENERATED, with the
regeneration commit citing this diagnostic's L0 evidence as the
principled-equivalence justification.

### Sentinel prediction

`test_phase_e_trajectory_resolvent_flux_shape_crosscheck` is currently
xfail-strict at Step 1. **Prediction: XPASS at Step 2** with the
Krylov-form closure unification, because:

- Variant α produces the canonical-form per-cell flux shape (no pole
  pathology — it's a Schur+power-iteration on the integral form).
- Krylov apply-matvec produces ~isotropic ψ at the pole (the L0 truth).
- Step 2's unified closure produces Krylov's flux shape.
- Therefore Step 2 ψ-shape matches Variant α ψ-shape.

The 5% tolerance on the sentinel is generous; Step 2's residual O(h)
boundary truncation (the apply-matvec's i=39 cell-centre artefact —
see Output A above, sf_kr[39] < sf_kr[38] by ~2.5%) might still produce
a small heterogeneous-MR drift. If the sentinel CHANGES from xfail to
xpass (i.e. drift < 5%), Step 2 is empirically closed; if drift is
still ~6%, the boundary face-trace fix in the apply matvec (Step 4
preview pulled forward) is also needed.

## ERR-NNN recommendations

**Update ERR-026 manifestation #6** in `error_catalog.md`:
- Pre-Phase-F: divergent (sf[0]/sf[1] → 0.473 under refinement, eigenvector
  shape catastrophically wrong).
- Post-Phase-F: stabilised but still wrong at L0 streaming-equilibrium
  (22% rel error at the pole, independent of refinement).
- Post-Phase-G Step 2 (planned): CLOSED — closure unified on the
  L0-passing apply-matvec form.

**NEW ERR-NNN to file**: "SI sweep curvilinear closure fails L0
streaming-equilibrium at the pole even with Phase F Carlson seed".
- Failure mode: #2 / #4 (closure-form / wrong recurrence — the M-M
  recurrence's IN-CELL state propagation through ordinates is the
  wrong structural form for the pole).
- How it hid: 22% pole error was masked because (a) the snapshot
  regression test enshrines the bug as the contract, (b) heterogeneous
  MR k-eigenvalue tests use generous tolerances, (c) L1 MMS rate (~1.3)
  was attributed to "ERR-026 pre-Phase-F" and not isolated.
- Test that catches it: a NEW
  `tests/sn/spatial/test_streaming_equilibrium_curvilinear.py` —
  the 1G homogeneous reflective sphere fixed-source test with the
  analytical Q/Σ_a comparison. Should be promoted from the diagnostic
  script in this memo.
- Lesson: the L0 streaming-equilibrium test was apparently NOT being
  run on curvilinear SI; only on Cartesian SI. Cardinal Rule 1 +
  vv-principles §1: necessity chain L1 → L0. The L0 gap is what let
  the L1 closure ship.

## Diagnostic scripts (to promote)

All under `tests/sn/diagnostics/`:

1. `phase_g_step2_00_baseline.py` — verifies Step 1 didn't perturb
   the Phase F baseline. **PROMOTE** to `tests/sn/regression/` as a
   guard against future drift.
2. `phase_g_step2_01_psi_comparison.py` — per-cell-per-ordinate ψ_si
   vs ψ_kr on k-eigenvalue. **PROMOTE** as a permanent regression for
   "Step 2 closes the eigenvector drift to <0.1%".
3. `phase_g_step2_02_sncell_residual.py` — KEEP as scratch (residual
   construction non-trivial).
4. `phase_g_step2_03_closure_audit.py` — the structural code-walk.
   **KEEP** under diagnostics (the audit text is the documentation).
5. `phase_g_step2_04_fixed_source.py` — fixed-source SI vs Krylov
   attribution at matched Q. **PROMOTE** as a Step 2 regression.
6. `phase_g_step2_05_homogeneous.py` — L0 streaming-equilibrium on
   curvilinear SI. **PROMOTE PRIORITY**: this catches the missing
   L0 gate. Move to `tests/sn/spatial/test_streaming_equilibrium_curvilinear.py`
   and tag `@pytest.mark.l0 @pytest.mark.catches("ERR-026")`. The
   test should ASSERT that BOTH SI and Krylov pass with rtol=1e-6.
   Currently SI FAILS by 22% → this test, once promoted, is an active
   regression gate that turns green at Step 2.

## Confidence: very-high

The L0 streaming-equilibrium failure on a homogeneous reflective
sphere is the canonical curvilinear correctness test. The 22%-at-pole
SI error vs the machine-zero Krylov error is unambiguous: the two
methods solve different operators, and Krylov solves the correct one.

The fixed-source attribution test (Output C) confirms that source
feedback is NOT the cause of the drift — the drift exists even at
matched Q, isolated to the closure algebra. Hypothesis H1 is refuted.

The structural code-walk (Output D) identifies three concrete
algebraic differences and explains why two of them (1, 2) agree on
the fixed point but the third (BC trace asymmetry in the apply-matvec)
is an O(h) issue. The dominant SI bug (Output C, E) is independent —
the SI's M-M-recurrence-inside-cell-update is itself structurally
wrong at the pole, regardless of BC.

## Pointers

- Phase G plan: `.claude/plans/issue_196_phase_g_four_operator_unification.md`
- Phase F closeout: `.claude/agent-memory/method-implementer/issue_168_phase_f_closeout.md`
- Step 1 closeout: `.claude/agent-memory/method-implementer/issue_196_phase_g_step1_closeout.md`
- Phase F Step 3 diagnostic (Phase D analog): `.claude/agent-memory/numerics-investigator/phase_f_step3_diagnostic.md`
- Architecture reconciliation: `.claude/agent-memory/general-purpose/phase_g_four_operator_architecture_reconciliation.md`
- ERR-026 catalogue entry: `.claude/skills/vv-principles/error_catalog.md` (under "manifestation #6").
- Apply-matvec code: `orpheus/sn/operator.py:571-838` (spherical), `:851-...` (cylindrical).
- SI sweep code: `orpheus/sn/sweep.py:397-595` (spherical), `:602-...` (cylindrical).
- DiamondDifference per-cell algebra: `orpheus/sn/spatial/diamond.py:550-624` (curvilinear branch).
- SNCellOperator wrapper: `orpheus/sn/spatial/operators.py:84-313`, `:319-367` (residual helpers).

## Linked memories

- `[[issue-168-phase-f-closeout]]` — Phase F partial closure; this
  diagnostic is the empirical proof that "partial" is the right word.
- `[[phase-f-step3-diagnostic]]` — the Phase D apply-vs-SI mirror diagnosis.
- `[[issue-196-phase-g-step1-closeout]]` — Step 1's SNCellOperator
  wrapper; Step 2 will consume it.
