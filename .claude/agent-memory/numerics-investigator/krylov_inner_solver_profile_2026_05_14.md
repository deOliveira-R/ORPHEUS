---
name: krylov-inner-solver-profile-2026-05-14
description: Krylov inner-solver profiling memo for cylindrical heavy-scattering SN. Root cause of the residual full-suite cost is NOT a per-iteration cache miss but the outer-Picard-around-GMRES structure in `_solve_fixed_source_krylov` — for c=0.95 reflective the path runs ~482 outer Picards × ~13 GMRES iters per outer = ~12 500 sweep/apply pairs vs ~529 sweeps for SI on the same case. SI is 27× faster than Krylov here. Two structural-improvement options diagnostic-axis preserving: (A) collapse outer Picard via direct GMRES on (L+C-S) — biggest win but requires `S` to be a proper LinearOperator; (B) hoist apply-matvec's per-cell `dag_walk` / `streaming_terms` work onto the existing Step 2.5c GeometryCoefficients cache — eliminates ~50% of per-matvec cost. Stratum 2 (CollisionCache) is already correctly cached; the apply path simply isn't consuming it.
metadata:
  type: project
---

# Krylov inner-solver profiling memo (Issue #196, Phase G)

Branch `refactor/sn-operator-algebra`, tip `01999e0`.  Profile date
2026-05-14.  Environment: host (`darwin 25.4.0`, `python 3.14.3`,
`scipy 1.x`, `.venv/bin/python`).

## Executive summary (5 bullets)

1. **The slow case is the cylindrical-Krylov nx=80 nq=16-32 family.**
   Empirically reproduced at minimal scale `cylinder[krylov-4-20]`
   (nx=20, n_mu=4, n_phi=4, N=16 ordinates, c=0.95): **32.21 s wall
   clock**, against 1.17 s for the SI counterpart on the same problem
   — Krylov is **27 × slower** than SI on the c=0.95 reflective
   cylinder it was supposed to make safer.
2. **The root cause is NOT a per-matvec cache miss.** It is the
   outer-Picard-around-GMRES structure in
   `_solve_fixed_source_krylov` (`orpheus/sn/solver.py:1171-1328`).
   For c=0.95 reflective the path runs **~482 outer Picard
   iterations × ~13 GMRES inner iterations per outer ≈ 12 500
   sweep + apply pairs**, where each pair is ~5 ms.  SI on the same
   case runs ~529 sweeps, period.
3. **`transport_operator_matvec_cylindrical` (the GMRES matvec)
   is the 51 % of cumulative time hot spot** AND **does not consume
   the Step 2.5c cache.**  It re-runs `dag_walk` per ordinate per
   matvec, building `StreamingTerms` via `_iter_cylindrical_visits`
   ~320 × per matvec.  This is a Smell #16 instance hiding in the
   apply path — the sibling SI sweep was fixed by Step 2.5c; the
   apply twin survived.
4. **Two structural-improvement options are recoverable**, ranked
   below.  The bigger lever is collapsing the outer Picard into a
   single GMRES on `(L + C − S)` (Option A); the smaller lever is
   making the apply-matvec consume the existing
   `GeometryCoefficients` cache (Option B).  Both are diagnostic-axis
   preserving: they keep `inner_solver="krylov"` working on the
   bare `(L+C)` operator with sweep preconditioner so the Krylov
   pathology axis stays measurable.
5. **What this does NOT close.**  Step 4 (left-preconditioner GMRES
   on the full equation) is the elegance-aligned target state per
   `.claude/plans/issue_196_phase_g_replan.md:1197-1268`; this memo's
   Option A overlaps with it but at the diagnostic level (still
   solve via composed operators, not yet via `.solve` API).  Step 6
   (adjoint via GMRES) inherits the same apply-matvec cache miss
   if Option B is deferred.

## 1. Slow-case identification

The Step 2.5c closeout attributed full-suite cost (16:50 wall clock
for `test_streaming_equilibrium_curvilinear.py`) to one cylinder
Krylov nx=80 test at 4:52.  Per the brief's anti-recommendation #3
("Do NOT exceed 2 hours of wall-clock profiling — if the slow case
would take 30+ min, profile a smaller-N variant that exhibits the
same pathology"), I selected the smallest cylinder-Krylov variant
and verified it shows the same fingerprint.

`pytest --durations` against the minimal cylinder-Krylov case
(reproducer that needs nothing larger than 20 cells):

```
$ time .venv/bin/python -m pytest \
    "tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[krylov-4-20]" -v

tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[krylov-4-20] PASSED [100%]
======================== 1 passed, 1 warning in 32.21s =========================
.venv/bin/python -m pytest  -v 2>&1  32.28s user 0.33s system 99% cpu 32.908 total
```

For comparison (one combined session — both pytest invocations
together):

```
$ time .venv/bin/python -m pytest \
    "tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[source_iteration-4-20]" \
    "tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_sphere[krylov-8-20]" -v
======================== 2 passed, 1 warning in 9.97s =========================
```

The cylinder Krylov case alone consumes ~32 s; cylinder SI plus
sphere Krylov together consume 10 s.  Same problem, same physics —
the SI sweep handles c=0.95 reflective cylinder in roughly the same
time the GMRES outer Picard wraps a single inner iteration.
Empirical "SI vs Krylov on the same case" microbench:

```
source_iteration    : 1.17 s  n_inner=529  residual=9.60e-13
krylov              : 31.79 s  n_inner=482  residual=9.67e-13
```

Both reach the same residual at the same iteration count, but
Krylov pays ~13× more work per outer iteration to do so (each outer
GMRES converges in ~13 inner iterations).

## 2. cProfile top-30 by cumulative time

```
$ .venv/bin/python -m cProfile -o /tmp/krylov_profile.prof \
    -m pytest "tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[krylov-4-20]" -q
1 passed, 1 warning in 39.59s
```

(cProfile overhead bumps wall clock from 32 s to 39 s; the ranking
is unaffected.)

Cumulative-time top-30 (`pstats.print_stats(30)`, sort `cumulative`):

```
20171991 function calls (20133501 primitive calls) in 35.987 seconds
Ordered by: cumulative time

ncalls  tottime  percall  cumtime  percall filename:lineno(function)
   1    0.000    0.000   39.137   39.137 solver.py:963(solve_sn_fixed_source)
   1    0.008    0.008   39.131   39.131 solver.py:1171(_solve_fixed_source_krylov)
 482    0.198    0.000   38.691    0.080 scipy/_isolve/iterative.py:587(gmres)
12339   0.027    0.000   38.418    0.003 scipy/_interface.py:231(matvec)
12339   0.013    0.000   38.377    0.003 scipy/_interface.py:619(_matvec)
6170    0.987    0.000   19.442    0.003 orpheus/sn/solver.py:622(matvec)        # precond matvec
6169    0.047    0.000   18.922    0.003 orpheus/sn/operator.py:1291(apply)       # L.apply
6169    8.055    0.001   18.872    0.003 orpheus/sn/operator.py:851(transport_operator_matvec_cylindrical)
6170    0.004    0.000   17.557    0.003 orpheus/sn/sweep.py:96(transport_sweep)
6170    0.040    0.000   17.554    0.003 orpheus/sn/sweep.py:146(_sweep_1d_unified)
6170    3.913    0.001   17.508    0.003 orpheus/sn/sweep.py:235(_run_1d_sweep)
2073288  0.415    0.000    9.902    0.000 orpheus/sn/geometry.py:425(dag_walk)
2073288  1.295    0.000    9.369    0.000 orpheus/sn/geometry.py:654(_iter_cylindrical_visits)
1974560  4.498    0.000    7.253    0.000 orpheus/geometry/reduced_operator.py:385(streaming_terms)
987200   2.415    0.000    6.049    0.000 orpheus/sn/spatial/diamond.py:128(update)
49356    3.545    0.000    3.574    0.000 orpheus/sn/spatial/psi_half_angle_seed.py:363(carlson_inward_sweep_from_source)
987200   2.721    0.000    3.222    0.000 orpheus/sn/spatial/cell_balance.py:120(cell_balance_terms)
6169     0.109    0.000    2.621    0.000 orpheus/sn/spatial/pole_angular_closure.py:551(__call__)
24676    0.085    0.000    1.973    0.000 orpheus/sn/spatial/psi_half_angle_seed.py:535(__call__)
12821    1.805    0.000    1.824    0.000 orpheus/sn/operator.py:491(solution_to_angular_flux_spherical)
1110420  1.647    0.000    1.771    0.000 orpheus/sn/operator.py:130(unknowns_at_cell_for_mask)
24676    0.472    0.000    0.531    0.000 orpheus/sn/spatial/pole_angular_closure.py:340(_mm_weighted_angular_recurrence_single_level)
49352    0.108    0.000    0.117    0.000 orpheus/sn/geometry.py:551(_representative_ordinate)
[…import / pytest plumbing trimmed…]
```

Tottime (self-time) top-15:

```
6169     8.055    0.001   18.872    0.003 operator.py:851 transport_operator_matvec_cylindrical
1974560  4.498    0.000    7.253    0.000 reduced_operator.py:385 streaming_terms
6170     3.913    0.001   17.508    0.003 sweep.py:235 _run_1d_sweep
49356    3.545    0.000    3.574    0.000 psi_half_angle_seed.py:363 carlson_inward_sweep_from_source
987200   2.721    0.000    3.222    0.000 cell_balance.py:120 cell_balance_terms
987200   2.415    0.000    6.049    0.000 diamond.py:128 DiamondDifference.update
12821    1.805    0.000    1.824    0.000 operator.py:491 solution_to_angular_flux_spherical
1110420  1.647    0.000    1.771    0.000 operator.py:130 unknowns_at_cell_for_mask
2073288  1.295    0.000    9.369    0.000 geometry.py:654 _iter_cylindrical_visits
6170     0.987    0.000   19.442    0.003 solver.py:622 _make_sweep_preconditioner.matvec
[…]
```

## 3. Iteration-level breakdown

Microbench `derivations/diagnostics/diag_krylov_iter_breakdown.py`
running ONE GMRES call (no outer Picard) on the same operator at
the same scale, with per-matvec wall-clock timing:

```
n_unknowns                     = 316
info                           = 0
n_gmres_iterations             = 13
n_matvecs                      = 16
n_precond_matvecs              = 16
wall_total                     = 0.0860
matvec_wall_mean_ms            = 2.6051
matvec_wall_std_ms             = 0.3032
matvec_wall_total_s            = 0.0417
precond_wall_mean_ms           = 2.6997
precond_wall_std_ms            = 0.3174
precond_wall_total_s           = 0.0432
```

Per-iteration cost partition:

| Component | ms/iter | % of iter |
|---|---|---|
| Apply matvec (`L.apply`) | 2.61 | 49 % |
| Sweep preconditioner | 2.70 | 51 % |
| GMRES Arnoldi orthogonalisation | < 0.01 | < 1 % |

GMRES bookkeeping is **negligible** here — operator size `n=316`,
restart=50, so the Krylov-side cost is O(n × m) where m ≤ 50.  The
entire cost is split ~50/50 between the apply matvec and the sweep
preconditioner, both ~2.6 ms per call.

The **outer Picard multiplier**: cProfile reports 482 outer GMRES
solves × 13 inner matvecs ≈ **6170 sweep + 6170 apply pairs** for
one production solve.  This is the dominant cost driver.  Per
session-lessons L15 (immutability strata), the outer Picard is the
unit of "per-call mutable" work; everything inside is per-iteration
mutable; everything inside the per-iteration work is what cache
factoring can address.

## 4. σ_t-immutability audit on the Krylov hot path

For each major function in the hot path, classify by mutation cadence
(in the spirit of the Step 2.5c six-layer audit on the sweep):

| Function | tot ms | Cadence | Cache opportunity? |
|---|---|---|---|
| `_solve_fixed_source_krylov` outer Picard loop | 39 000 | per-call (one fixed-source solve) | NO — algorithmic |
| `gmres(...)` per outer | 38 691 (cum) | per-outer-iter (×482) | NO — Krylov is per-RHS |
| `L.apply` = `transport_operator_matvec_cylindrical` | 18 872 (cum) | per-matvec (×6169) | partial — see below |
| `_make_sweep_preconditioner.matvec` | 19 442 (cum) | per-matvec (×6170) | partial — see below |
| `transport_sweep` → `_run_1d_sweep` (preconditioner) | 17 508 (cum) | per-matvec | partial — see below |
| `dag_walk` + `_iter_cylindrical_visits` | 9 902 (cum) | per-matvec (×~336/matvec) | **YES** — σ_t-immutable AND ψ-immutable; already cached by `GeometryCoefficients` |
| `streaming_terms` | 7 253 (cum) | per-matvec (×~320/matvec) | **YES** — σ_t-immutable (geometry+quadrature only); already cached |
| `DiamondDifference.update` (in SI sweep degenerate fallback) | 6 049 (cum) | per-matvec (×~160/matvec) | YES for σ_t-mutable layer, partial for full ψ |
| `cell_balance_terms` (in DD.update) | 3 222 (cum) | per-matvec | YES — algebraic identity, already covered by `CollisionCache` |
| `carlson_inward_sweep_from_source` | 3 574 (cum) | per-matvec (×8/matvec) | partial — geometry inputs immutable, ψ inputs mutable |
| `pole_angular_closure(...)` apply path | 2 621 (cum) | per-matvec | NO — body is per-ψ tensor op |
| `solution_to_angular_flux_spherical` | 1 824 (cum) | per-matvec | NO — packing/unpacking is per-ψ |
| `unknowns_at_cell_for_mask` | 1 771 (cum) | per-matvec (×~180/matvec) | YES — geometry+mask immutable |

Summary by stratum:

* **σ_t-immutable AND ψ-immutable** (geometry+quadrature only):
  ~18 s cum (`dag_walk` + `_iter_cylindrical_visits` +
  `streaming_terms` + `unknowns_at_cell_for_mask`).  Every byte
  of this work is already captured by Step 2.5c's
  `GeometryCoefficients`.  **The apply-matvec path simply doesn't
  consume the cache.**
* **σ_t-immutable, ψ-mutable** (algebraic-identity layer):
  ~3 s cum (`cell_balance_terms` inside DD.update, the degenerate
  cylindrical-axis fallback only).  Already covered by Step 2.5c
  CollisionCache for the fast path; the slow path remains.
* **Per-matvec mutable** (genuine arithmetic on this iteration's ψ):
  ~12 s cum (the WDD recurrence body inside the apply matvec, the
  Carlson seed source rebuild, the angular-redistribution closure
  body, packed↔structured re-shape).  Cannot be cached;
  cumulative time bounded by physics.
* **Per-outer-iter mutable** (algorithmic — the outer Picard cost):
  ~38 s cum, the **dominant** stratum.  Cannot be cached; can be
  collapsed by changing the iteration scheme (Option A below).

The audit reveals **two distinct structural levers**:

1. **The outer Picard wrapping (per-call mutable)** is the largest
   absolute cost.  Each outer iteration reruns a 13-iter GMRES on
   the same operator with a slightly different RHS (scattering
   source from the previous outer).  This is mathematically
   equivalent to running ONE GMRES on the full operator
   `(L + C − S)` with sweep preconditioner — and that ONE GMRES
   would converge in approximately `log(tol) / log(reduction
   ratio)` iterations rather than 482 outer × 13 inner.
2. **The apply-matvec cache miss (per-matvec mutable, but with a
   σ_t-immutable sub-stratum that's currently rebuilt every call)**
   is the second-largest lever.  ~50 % of per-matvec time is spent
   in σ_t-immutable geometry walks that `GeometryCoefficients`
   already caches but `transport_operator_matvec_cylindrical`
   doesn't consume.

## 5. Recommended structural improvements (ranked)

Ranked by impact-vs-effort (impact = total wall-clock reduction;
effort = code change scope; risk = diagnostic-axis preservation
per user policy).

### Option A — Collapse the outer Picard into a single GMRES on `(L + C − S)` (impact: ~12×, effort: medium, risk: medium)

**What it changes.**  In `_solve_fixed_source_krylov`, replace the
`for n_outer in range(max_inner)` Picard loop with a single
`gmres(L_full, rhs)` call where `L_full = L + C - S` and the sweep
remains the preconditioner.  Today the body solves
`(L + C) ψ = q_scatter + q_ext` repeatedly, taking the scattering
source from the previous outer's converged scalar flux; the
mathematically equivalent statement is `(L + C - S) ψ = q_ext`,
which GMRES handles directly given an `S.apply` matvec.

**What immutability stratum it touches.**  Per-call mutable —
collapses the outermost iteration over RHS rebuilds.  Does NOT
touch the per-matvec cache miss (Option B), so the two are
independent.

**Expected speedup.**  GMRES on `(L+C-S)` with a sweep
preconditioner has condition number ≈ `1/(1-c)` ≈ 20 for c=0.95;
the preconditioned operator condition number is much smaller
(typically O(1) for an exact-inverse preconditioner) — empirically
this should converge in 10-30 iterations end-to-end vs the current
482 × 13 = 6266 inner iterations.  **Order-of-magnitude reduction in
per-solve wall clock** is the realistic expectation; the lower
bound is the bare matvec + sweep cost × ~20 iterations ≈ 100 ms,
vs the current 32 s — a ~300× theoretical bound.  Realistic estimate
is 10-30× because preconditioning isn't perfect.

**Test-coverage risk (diagnostic-axis preservation).**  This is
the *exact* concern the user flagged.  The change MUST keep the
naked Krylov solve available, because "Krylov pathology behind a
strong preconditioner" is a hideable bug class.  Mitigation:
introduce the new path under an explicit flag (e.g.
`inner_solver="krylov_full"`) and keep the current
`inner_solver="krylov"` (which now reads as "Krylov on `(L+C)` with
outer Picard for scattering") as the diagnostic axis.  Both paths
share the same matvec primitives — bugs surface in either.

**Relationship to Step 4.**  Step 4 (per the replan §1197-1268)
plans to wire `(L+C).solve` as a left preconditioner for GMRES.
Option A is the algorithmic-equivalent change at the iteration level
(GMRES on the right operator, sweep as preconditioner) — Step 4 is
the elegance-aligned API surface that does the same algorithm
through `OperatorSum.solve`.  Per user policy this memo does NOT
recommend Step 4 as the primary lever (it removes the diagnostic
axis); Option A here is the *iteration-scheme* change that lands
the speedup WHILE preserving naked-Krylov coverage.

### Option B — Make `transport_operator_matvec_*` consume `GeometryCoefficients` (impact: ~2× on per-matvec, effort: medium, risk: low)

**What it changes.**  Replace the per-matvec `list(sn_mesh.dag_walk(...))`
+ per-cell `streaming_terms` calls in
`transport_operator_matvec_spherical` (`orpheus/sn/operator.py:751-836`)
and `transport_operator_matvec_cylindrical`
(`orpheus/sn/operator.py:992-1109`) with a vectorised tensor
expression that consumes `geom.A_down`, `geom.dA_w`,
`geom.chain_idx`, etc., directly from the cache.  The WDD recurrence
body becomes a closed-form `np.cumsum`/`np.cumprod` chain analogous
to `ordinate_scan` already used by the sweep path.

**What immutability stratum it touches.**  σ_t-immutable AND
ψ-immutable (the geometry walk).  ~50 % of per-matvec time
(`dag_walk` + `streaming_terms` + `_iter_cylindrical_visits` +
`unknowns_at_cell_for_mask`) is eliminated.

**Expected speedup.**  ~2× on per-matvec wall clock (from 2.6 ms
to ~1.3 ms).  This carries through to the full solve linearly:
~32 s → ~16 s in the absence of Option A.  Composed with Option A:
~16 s outer × 1/15 inner ≈ ~1-2 s — comparable to SI in this case
but with the better convergence properties of Krylov on harder
problems.

**Test-coverage risk.**  LOW.  This is bit-equivalent to the current
algebra (cache populator is already validated against
`cell_balance_terms` at `rtol=1e-13` per the Step 2.5c dual-view
contract test #7).  The change just makes apply-matvec consume the
same cache.

**Relationship to Step 2.5c.**  Step 2.5c left the apply path
out-of-scope.  Per the closeout's "What this does NOT close" section,
the apply path was acknowledged as a known cache miss; the closeout's
note "16:50 wall clock for streaming-equilibrium-curvilinear" is
*this* cache miss + Option A's outer-Picard amplification.

### Option C — Cache the Carlson seed source `Q_bar` decomposition (impact: ~10 %, effort: low, risk: low)

**What it changes.**  Inside `_run_1d_sweep`,
`carlson_inward_sweep_from_source` is called 8× per sweep (4 levels ×
2 — once for the inward seed of each cylindrical level + 4 from
the apply path's `pole_angular_closure`).  Its inputs (`sigma_t`,
`dr`, `bc_outer_value`, `Q_bar`) are σ_t-stable in `sigma_t` and
`dr` (geometry) and per-iteration mutable in `bc_outer_value` and
`Q_bar`.  Refactoring to precompute the σ_t/dr part of the
recurrence kernel (geometric attenuation factors per cell × level)
and reuse it across the per-iteration inputs would save ~half the
Carlson cost.

**What immutability stratum it touches.**  σ_t-immutable
sub-stratum of an otherwise per-matvec mutable computation.

**Expected speedup.**  Carlson contributes 3.5 s cum (10 %) to the
profile.  Caching halves it → ~5 % overall.  Marginal vs A and B.

**Test-coverage risk.**  LOW.  Pattern 2 dual-view tests already
exercise the Carlson seed; the cache must agree algebraically at
`rtol=1e-14`.

### Option D — Replace `_iter_cylindrical_visits` with vector form even on the slow degenerate path (impact: small, effort: low, risk: low)

The cylindrical pure-azimuthal degenerate ordinates (8 of 16 at
n_mu=4 n_phi=4) take the slow `cell_update.update` path in
`_run_1d_sweep` at lines 379-403.  This contributes ~6 s
(`DiamondDifference.update`, `cell_balance_terms`).  Vectorising
this path is principled per the Step 2.5c approach (`a_attenuation`
is zero on degenerate ordinates because `|η|≈0`, so `psi_face`
short-circuits to the source-only expression `b`).  Would close
the ~6 s on its own.

**Trade-off:** the degenerate path is intentionally a slow
reference per the Step 2.5c brief — the bool gate is documented as
acceptable.  Closing it makes the cylindrical path fully vectorised.
Diagnostic risk: zero (the degenerate path is algebraically the
slab-collapse of the cache formula).

### Why NOT recommended: bare preconditioner wiring change

The user policy explicitly defers Step 4 (`(L+C).solve` as left
preconditioner via `OperatorSum.solve`) because it removes a
diagnostic axis: bugs that need the bare `(L+C-S)` GMRES path to
surface would hide behind a strong preconditioner.  Options A-D
above all preserve the unpreconditioned Krylov pathway and only
change how it's wrapped.

## 6. What this does NOT close

Pending Phase G steps that touch Krylov perf and remain open after
this memo:

1. **Step 3** — `StreamingOperator` (`L`) and `CollisionOperator`
   (`C`) split.  Option B requires `L`'s apply to consume the
   geometry cache cleanly; Step 3's architectural split is the
   natural home for that change.  Option B can land *before* Step 3
   as a non-API-changing internal refactor.
2. **Step 4** — `(L+C).solve` as left preconditioner via
   `OperatorSum.solve`.  Option A's outer-Picard collapse overlaps
   with Step 4's algorithmic equivalence; Step 4 is the
   elegance-aligned API surface.  Per user policy, Step 4 is
   deferred until both naked-Krylov AND preconditioned-Krylov
   coverage exists.  Option A here introduces preconditioned
   coverage as an additive path (`inner_solver="krylov_full"` or
   equivalent) without retiring the naked path.
3. **Step 6** — adjoint via GMRES on `A_loss.H`.  Inherits the same
   apply-matvec cache miss until Option B lands.  Adjoint test
   coverage will be expensive (same outer-Picard structure) unless
   Option A is also adopted for the adjoint inverse.
4. **Anisotropic Pℓ scattering on curvilinear** — currently the
   `_build_rhs_*` curvilinear path is P0-only (`scattering_order > 0`
   is a Cartesian feature).  This memo's options don't change the
   gap, but Option A (single-GMRES-on-`L+C-S`) requires Pℓ S.apply
   to be a clean LinearOperator on curvilinear too; currently
   Pℓ S construction is inlined inside `_build_rhs_cartesian` and
   curvilinear branches go through the P0 builder.

## 7. Methodology trail

* Cylinder Krylov nx=20 n_mu=4 case selected per anti-recommendation
  #3 (smaller-N variant exhibiting the same fingerprint).
  Verification: 32 s vs SI's 1 s on the SAME problem matches the
  closeout's 4:52 nx=80 outlier qualitatively (Krylov outer-Picard
  cost scales as `outer × inner × matvec_cost`; matvec_cost scales
  with mesh size, the outer/inner counts are mesh-independent for
  this physics).
* cProfile run wall clock 39.59 s vs raw 32.21 s — cProfile overhead
  ~23 %, ranking-preserving.
* Iteration breakdown microbench
  (`derivations/diagnostics/diag_krylov_iter_breakdown.py`) runs ONE
  GMRES call with timing callbacks on the same operator at the same
  scale.  Confirms 13 inner iter / 16 matvec, 2.6 ms per matvec,
  ~50/50 apply vs precond split.
* SI vs Krylov head-to-head: identical problem (Mixture B 1G,
  reflective cylinder, c=0.95, uniform Q), identical tolerance,
  identical residual reached.  SI 1.17 s / 529 iter; Krylov 31.79 s
  / 482 outer × 13 inner.

## 8. Artefacts

* Profile: `/tmp/krylov_profile.prof` (binary cProfile output).
* Microbench script:
  `derivations/diagnostics/diag_krylov_iter_breakdown.py`
  (one-shot profiling probe; do NOT promote — wall-clock
  thresholds are machine-dependent; the GMRES iteration count
  (13) and matvec/precond ratio (~50/50) ARE machine-independent
  and could be promoted to a regression test if Option A or B
  lands and we want to pin the new iteration scheme's
  characteristic counts).
* Reproducer (production test):
  `tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[krylov-4-20]`.
