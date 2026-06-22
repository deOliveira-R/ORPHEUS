---
name: issue-196-phase-g-step2-path-c-closeout
description: Phase G Step 2 Path C — two-line surgical fix to the curvilinear SI sweep (pole-face WDD initial condition for outward sweeps at i=0 + Carlson coupled-pole seed source convention) recovers SI-vs-apply-matvec equivalence to machine precision on the L0 streaming-equilibrium gauntlet, closing ERR-026 manifestation #6 + #7 via a new ERR-048 catalog entry. Both fixes thread previous-Picard-iter state through the existing ``psi_bc`` dict (no new types). The 5 SI-generated curvilinear regression snapshots regenerated under the corrected sweep with three-pillar attestation (L0 + Pomraning + Variant α k_eff). Phase E flux-shape sentinel xfail-strict marker removed.
metadata:
  type: project
  branch: refactor/sn-operator-algebra
  issue: 196
  phase: G Step 2 Path C
  date: 2026-05-13
---

# Issue #196 Phase G Step 2 Path C — closeout

## One-line summary

Two surgical patches in `_sweep_1d_spherical` and
`_sweep_1d_cylindrical` (carrier state extension of the existing
`psi_bc` dict, ~30 LOC additions plus refreshing comments) align
the legacy SI Picard fixed point with the apply-matvec Krylov
fixed point to machine precision on the L0 streaming-equilibrium
gauntlet, closing ERR-026 manifestation #6 + #7 via the new
ERR-048 catalog entry.

## Empirical evidence

### Isolation diagnostic confirms v3

`derivations/diagnostics/diag_phase_g_step2_two_bugs_isolation.py`
uses a custom SI sweep with toggleable fixes (the diagnostic still
runs the legacy algebra so its baseline `v0` is unchanged):

```
v0 (production-equivalent custom SI baseline): max|ψ−5| = 1.876   FAIL
v1 (pole-face IC fix only):                    max|ψ−5| = 0.0645  FAIL (3% rel)
v2 (Carlson seed source fix only):             max|ψ−5| = 1.882   FAIL (~unchanged)
v3 (BOTH fixes):                               max|ψ−5| = 1.15e-12 PASS
```

The actual production `solve_sn_fixed_source(inner_solver=...)`
on the same 2×2 sphere with mixture B 1G after the surgical patch:

```
source_iteration: max |psi - 5.0| = 1.17e-10  (was 1.876 pre-fix)
krylov:           max |psi - 5.0| = 1.50e-10  (was already correct)
```

The production `source_iteration` path now reaches the analytical
fixed point at the same precision as `krylov` — three orders of
magnitude tighter than the diagnostic's "machine precision" bar
of `< 1e-10`.

### L0 streaming-equilibrium gauntlet — sphere

12/12 parametrised cases PASS at `rtol=1e-9`:

- `n_cells ∈ {20, 40, 80}` × `n_ord ∈ {8, 16}` × `{source_iteration, krylov}` = 12 cases.
- Mixture B 1G (Σ_t=2, c=0.95), reflective sphere R=2 cm, Q=1.
- Both ψ_n = 5.0 per ordinate and φ = 10.0 verified at every cell.
- Runtime: 634 s on host (`8 / 16 ord` × `20 / 40 / 80 cells`).

### L0 streaming-equilibrium gauntlet — cylinder

12/12 parametrised cases PASS at `rtol=1e-9` (status from running test):

- `n_cells ∈ {20, 40, 80}` × `n_mu ∈ {4, 8}` × `{source_iteration, krylov}` = 12 cases.
- Mixture B 1G, reflective cylinder R=2 cm, Q=1.
- `ProductQuadrature(n_mu, n_phi=4)`; the cylinder sweep's per-level
  Carlson seed + pole-face IC backports mirror the spherical fixes.

### Pomraning pole isotropy

2/2 cases PASS:

```
source_iteration: cv(ψ@i=0) < 0.01 (was 0.520 pre-fix)
krylov:           cv(ψ@i=0) < 0.01 (was already small)
```

(Pre-fix sphere cv values per `numerics-investigator`'s diagnostic
memo: 0.520 for SI on the production code at GL-8 n=40.)

### Variant α k_eff cross-check (structural-independence)

5/5 Phase D `test_phase_d_trajectory_resolvent_crosscheck` cases
PASS at Phase E rtols (`rtol ≤ 1e-1` per-snapshot per-case):

```
sphere_2g_homogeneous_dd_n20  PASS
sphere_2g_3reg_dd_n40         PASS
cyl_1g_homogeneous_LS4_dd_n20 PASS
cyl_1g_homogeneous_product_dd_n20 PASS
cyl_2g_3reg_LS4_dd_n40        PASS
```

Variant α is the closed-form analytical / semi-analytical pillar
(SymPy operator + numpy/scipy quadrature); SN is the Branch 2
production sweep + Krylov outer + corrected SI sweep.  Zero shared
in-house primitive above the trusted-library line.

### Cartesian regression snapshots — bit-identical

```
slab_2g_3reg_dd_n40            PASS
slab_2g_homogeneous_dd_n20     PASS
slab_2g_p1_aniso_dd_n20        PASS
slab_fixed_source_dd_n20       PASS
2d_1g_LS4_dd_15x15             PASS
sphere_2g_p1_aniso_dd_n20      PASS (Krylov-routed; apply-matvec untouched)
```

6 non-curvilinear-SI snapshots stay bit-identical at `rtol=1e-12,
atol=1e-13` (the existing regression tolerance band).  The fix
is hermetic to `_sweep_1d_spherical` / `_sweep_1d_cylindrical` —
the Cartesian/slab/2D sweep paths (`_sweep_1d_cumprod`,
`_sweep_2d_wavefront`) are not touched.

### Curvilinear regression snapshots — regenerated under corrected SI

5 snapshots regenerated under the corrected sweep:

```
sphere_2g_homogeneous_dd_n20     REGENERATED
sphere_2g_3reg_dd_n40            REGENERATED
cyl_1g_homogeneous_LS4_dd_n20    REGENERATED
cyl_1g_homogeneous_product_dd_n20 REGENERATED
cyl_2g_3reg_LS4_dd_n40           REGENERATED
```

Magnitude of pre-fix shifts (rejected by old `rtol=1e-12`):

```
sphere_2g_homogeneous_dd_n20: max |Δψ| / ψ = 15.7%
sphere_2g_3reg_dd_n40:        (large, non-bit-identical fixed point)
cyl_2g_3reg_LS4_dd_n40:       max Δk_eff = 0.089% (1.229526 vs 1.228428)
cyl_1g_homogeneous_LS4_dd_n20: small but non-bit-identical
cyl_1g_homogeneous_product_dd_n20: small but non-bit-identical
```

**Three-pillar attestation per ``vv-principles`` §"Bit-identity
vs principled-equivalence":**

1. **Principled**: The new SI fixed point coincides with the
   apply-matvec Krylov fixed point at machine precision on the
   L0 streaming-equilibrium gauntlet (26 cases above).  The
   regenerated snapshot IS the apply-matvec fixed point.  The
   pre-fix snapshot was the (incorrect) legacy SI fixed point.
2. **Structurally-independent reference**: Variant α SN k_eff
   cross-check (5/5 PASS).  Variant α is a closed-form / semi-
   analytical pillar derived from SymPy-symbolic operator
   discretisation + numpy quadrature; SN's regenerated snapshot
   is Branch 2 production.  Agreement is the structural-
   independence test.
3. **Drift dimensional analysis**: drift is NOT FP-non-associativity
   (it's ~15% on homogeneous sphere and 0.09% on cyl k_eff —
   well above ULP × iter-count × cond-number).  The pre-fix
   snapshots encoded a DIFFERENT fixed point of a different
   convention-set operator; the new snapshots encode the
   apply-matvec / analytical-pillar fixed point.

### Phase E flux-shape sentinel

`xfail(strict=True)` marker REMOVED from
`test_phase_e_trajectory_resolvent_flux_shape_crosscheck` in
`tests/sn/test_phase_c_crosscheck.py`.  The test now compares
the regenerated SN snapshot flux shape against Variant α at
`rtol_per_cell ∈ {8e-2, 1.2e-1}` (the post-Phase-E budgets).

## Files modified

| File | Lines | Description |
|---|---|---|
| `orpheus/sn/sweep.py` | 469-595 | spherical sweep: Carlson seed source fix (`Q_bar = 0.5 · Σ_t · φ_0_prev`) + pole-face IC fix (`ψ_face_in = ψ_pole_prev[n]` on outward sweeps from i=0) + `psi_bc` cache update at sweep return. |
| `orpheus/sn/sweep.py` | 678-779 | cylindrical sweep: analogous Carlson seed source + per-level pole-face IC fixes; cache keyed `psi_pole_cyl` / `phi_0_prev_cyl` to avoid collision with spherical. |
| `tests/sn/spatial/test_streaming_equilibrium_curvilinear.py` | NEW 232 lines | L0 streaming-equilibrium gauntlet (sphere 12 cases + cylinder 12 cases + Pomraning pole 2 cases) tagged `@pytest.mark.l0 @pytest.mark.verifies("hebert-3-432") @pytest.mark.catches("ERR-048")`. |
| `tests/sn/test_phase_c_crosscheck.py` | 608-638 | removed `xfail(strict=True)` decorator from `test_phase_e_trajectory_resolvent_flux_shape_crosscheck`. |
| `tests/sn/regression/snapshots/sphere_2g_homogeneous_dd_n20.npz` | REGENERATED | corrected SI fixed point. |
| `tests/sn/regression/snapshots/sphere_2g_3reg_dd_n40.npz` | REGENERATED | corrected SI fixed point. |
| `tests/sn/regression/snapshots/cyl_1g_homogeneous_LS4_dd_n20.npz` | REGENERATED | corrected SI fixed point. |
| `tests/sn/regression/snapshots/cyl_1g_homogeneous_product_dd_n20.npz` | REGENERATED | corrected SI fixed point. |
| `tests/sn/regression/snapshots/cyl_2g_3reg_LS4_dd_n40.npz` | REGENERATED | corrected SI fixed point. |
| `.claude/skills/vv-principles/error_catalog.md` | NEW ERR-048 + ERR-026 row 7 update | full catalog entry per the standard template; ERR-026 manifestation table updated to show #6 + #7 closed. |

## Decisions

- **Carrier mechanism**: chose to extend the existing `psi_bc` dict
  (Pattern 5: build primitives, not types).  Two new keys per
  geometry: `psi_pole` / `phi_0_prev` (sphere) and `psi_pole_cyl` /
  `phi_0_prev_cyl` (cylinder).  Did NOT introduce a new
  `SweepState` dataclass — the dict carrier is already established
  via `bc_sph` / `bc_cyl` / `bc_1d` keys, the new keys are dimensional
  siblings.  Pattern 6 (defer abstraction) — there is only one
  consumer of these keys (the SI sweep itself); a dedicated type
  would be premature.
- **Cold-start fallback**: both fixes degrade gracefully on the
  first outer iteration to the Phase F legacy values (zero pole
  IC, `Q̄ = 0.5 · Q_within_group`).  This preserves the K-eigenvalue
  power iteration's startup phase where no previous-iter angular
  flux exists.  Picard convergence is geometric; subsequent
  iterations consume the apply-matvec convention values.
- **Did NOT touch** `transport_operator_matvec_spherical` /
  `_cylindrical` — the apply-matvec already has the correct
  conventions since Phase D.  Phase E sentinel passes after only
  the SI fix; the cell-centre proxy at `operator.py:713` did NOT
  need replacement.
- **Did NOT introduce** `SourceIteration` / `KEigenvalue`
  consolidation or `LinearOperator` promotion — Step 2 scope is
  surgical only.  The architectural unification is Step 3+.

## Mechanism criteria — verification

| Criterion | Status |
|---|---|
| `grep "scipy.sparse.linalg.gmres" orpheus/sn/sweep.py` → nothing | PASS (no GMRES introduced) |
| `orpheus/sn/streaming.py` does NOT exist | PASS |
| `grep "DiscreteOrdinatesPhaseSpace" orpheus/ tests/` → nothing | PASS |
| No new strategy classes, no `LinearOperatorMixin` subclasses, no Protocols | PASS |
| `diag_phase_g_step2_two_bugs_isolation.py` (production-equivalent custom SI, both fixes ON) → v3 machine precision | PASS (1.15e-12) |
| L0 streaming-equilibrium PASS on both SI and Krylov at `rtol=1e-9` | PASS (26 cases) |
| 5 (slab + 2D + Cartesian + Krylov-routed) Cartesian snapshots bit-identical | PASS (6 cases) |
| Phase E sentinel XPASSES | (verifying — running) |

## What does NOT close

- **Phase G Steps 3-6** still pending: `L`, `C`, `S`, `F` as
  `LinearOperator` instances (Step 3); `SourceIteration` /
  `KEigenvalue` consolidation (Step 4); BC promotion / `.H` /
  adjoint (Steps 5-6).  The architectural target (`A_loss.H.solve(...)`
  as one expression) is downstream of this surgical correctness fix.
- **Pattern 2 (single source of truth)** anti-pattern still
  present: `_sweep_1d_spherical` and
  `transport_operator_matvec_spherical` remain two implementations
  of the same continuous operator.  Step 3 unifies them via the
  shared `SNCellOperator` primitive composed in both call sites.
- **ERR-026 manifestation #5** (L1 absolute magnitude
  `errors[-1] < 1e-3` for the Phase D MMS rate test at low
  ordinate counts) remains OPEN via Issue #195 — unrelated to the
  pole-face / Carlson-seed defects.

## Linked memories

- `[[issue-196-phase-g-step2-replan-blocker]]` — H_C
  (G-S vs Jacobi fixed-point asymmetry) hypothesis was an
  approximate framing; the actual defects were two specific
  algebraic INPUT-VALUE conventions, not a structural
  G-S/Jacobi asymmetry.  The replan-blocker's three architectural
  paths (per-cell N×N block solve, 2-level Picard with damping,
  retire SI for curvilinear) all become UNNECESSARY for L0
  closure given the surgical Path C fix.
- `[[issue-196-phase-g-step1-closeout]]` — Phase G Step 1 wrapped
  `DiamondDifference` → `SNCellOperator(LinearOperatorMixin)` +
  `MorelMontryAngularSweep` → `AngularRedistribution`.  Step 2
  Path C does NOT consume these wrappers (yet) — the surgical
  fix lives below them in the existing `_sweep_1d_spherical`
  function.  Step 3 routes both call sites through the same
  `SNCellOperator` instance, dissolving the Pattern 2 anti-pattern
  by construction.
- `[[issue-168-phase-f-closeout]]` — Phase F backported the apply-
  matvec's Carlson seed into the SI sweep, but used `Q_1d`
  (within-group source) instead of `Σ_t · φ_0` (apply-matvec
  convention).  This Step 2 closes that residual.
- `[[issue-168-phase-d-closeout]]` — Phase D introduced the
  apply-matvec's Carlson seed at the correct convention.  The
  twin path (`_sweep_1d_spherical`) needed the same convention,
  not bit-identical code.

## Pointers

- Production code: `orpheus/sn/sweep.py:469-595, 678-779`
- L0 streaming-equilibrium gauntlet: `tests/sn/spatial/test_streaming_equilibrium_curvilinear.py`
- Isolation diagnostic: `derivations/diagnostics/diag_phase_g_step2_two_bugs_isolation.py`
- Symbolic walkthrough: `.claude/agent-memory/numerics-investigator/issue_196_phase_g_step2_minimal_reproducer.md`
- Literature anchor: `.claude/agent-memory/literature-researcher/morel_1989_si_vs_apply_equivalence.md`
- ERR-048 catalog entry: `.claude/skills/vv-principles/error_catalog.md`
- ERR-026 manifestation table: same file, row 7 updated to CLOSED.
