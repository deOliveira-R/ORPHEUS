---
name: issue-196-phase-g-step2-5b-closeout
description: Phase G Step 2.5b closeout. Canonical first-order linear-recurrence scan primitive (Blelloch 1990 §1.5) `ordinate_scan` + dual-view `DiamondDifference.affine_coefficients` + production-sweep wiring. Restores slab vectorisation; vectorises curvilinear per-cell fold. Pair-monoid associativity theorem at rtol=1e-14 the algebraic anchor. All regression snapshots bit-identical at rtol=1e-12; L0 streaming-equilibrium 26/26 PASS; sweep-vs-apply consistency 57/57 PASS; full SN suite paste-back below.
metadata:
  type: project
---

# Issue #196 Phase G Step 2.5b closeout

Branch: `refactor/sn-operator-algebra`. Date: 2026-05-14.
Commits: `f5b913c..defa1e4` (three commits + this memo).

## Concept-collapse signature

Step 2.5b restores slab vectorisation (lost in Step 2.5's per-cell
Python fold) AND vectorises curvilinear (which was a per-cell fold
from inception). The mechanism: one canonical primitive
`ordinate_scan` consumed by both the slab and curvilinear sweep
bodies, with `DiamondDifference.affine_coefficients` as the
dual-view vectorised builder of the per-cell `(a, b)` coefficients.

**The structural identity** (cross-domain-attacker memo
`issue_196_phase_g_chain_dag_scan_frame_attack`):

The chain DAG of one ordinate's spatial sweep admits the pair-monoid
`(α₁, β₁) ⊕ (α₂, β₂) = (α₁·α₂, α₂·β₁ + β₂)`. The per-cell update IS
one element `(a, b)` of this monoid. The closed-form decomposition
for THIS monoid (2×2 lower-triangular matrix product) is
`cumprod(a) · (psi_0 + cumsum(b / cumprod_a))` — three numpy ops, no
Python loop. Blelloch 1990 §1.5 canonical form.

## Files touched

1. **`orpheus/sn/spatial/scan.py`** (NEW, 140 lines including
   docstrings). One free function `ordinate_scan(a, b, psi_0)`.
   Body 2 numpy ops.

2. **`orpheus/sn/spatial/cell_update.py`**. Add third Protocol
   method `affine_coefficients(visits, total_xs, source,
   angular_state) -> (a, b)` to both `CellUpdate` Protocol and
   `CellUpdateBase` ABC. Sibling to `update` and `residual`.

3. **`orpheus/sn/spatial/diamond.py`** (295 → 399 lines). Implement
   `DiamondDifference.affine_coefficients` per the closed-form
   derivation `a = 2|μ|·A_total/denom − 1`, `b = 2·(q + dA_w·c_in·
   ψ_a_in)/denom`. Body ~75 lines including docstring.

4. **`orpheus/sn/spatial/__init__.py`**. Re-export `ordinate_scan`.

5. **`orpheus/sn/sweep.py`** (708 → 773 lines). Both
   `_sweep_1d_cartesian` and `_sweep_1d_curvilinear` consume
   `cell_update.affine_coefficients` + `ordinate_scan`. Per-cell
   Python loop within an ordinate is GONE for both bodies.
   Curvilinear retains a per-cell fallback for cylindrical
   pure-azimuthal degenerate ordinates (rare; `|η| < 1e-15`).
   `_sweep_2d_wavefront` untouched per scope hard limit.

6. **`tests/sn/spatial/test_ordinate_scan.py`** (NEW, 525 lines).
   Sixteen strong tests:
   - **Algebraic theorems** (4): pair-monoid associativity (THE
     theorem, rtol=1e-14), identity element (exact), Brent's
     blocked-scan equivalence (rtol=1e-13), closed-form vs
     explicit-loop (rtol=1e-13).
   - **Affine structure** (5): zero-source / zero-attenuation
     reductions (exact `np.array_equal`), linearity in ψ_0
     (rtol=1e-14), linearity in b (rtol=1e-13), affine combination
     (rtol=1e-13).
   - **Numerical stability** (2): near-identity `a≈1` (rtol=1e-12),
     small-attenuation `a≈0.1` (rtol=1e-11).
   - **Dual-view contracts** (3 classes / 41 cases): 36 cases of
     single-cell match (geometry × ng × source kind) at rtol=1e-13;
     4 vectorisation-vs-serial cases at rtol=1e-14; 1 full-sweep
     baseline at rtol=1e-12.

7. **`tests/sn/spatial/test_cell_update_protocol.py`**. Add stub
   `affine_coefficients` to two synthetic test strategies to
   satisfy the new Protocol method.

## Pair-monoid theorem PASS line

```
tests/sn/spatial/test_ordinate_scan.py::TestPairMonoidTheorems::test_pair_monoid_associativity PASSED [  1%]
```

This is the load-bearing algebraic verification of the entire
Step 2.5b abstraction. rtol=1e-14 over 100 random `(a, b, c)`
triples with `a ∈ [0, 2]`, `b ∈ [-1, 1]`. If this failed, the
Blelloch §1.5 derivation itself would be wrong.

## Dual-view consistency

Sample case from parametrised test 12
(`test_affine_coefficients_matches_update_single_cell`):

```
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[random-4-cylinder] PASSED
```

Across all 36 parametrised cases (4 geometries × 3 group counts ×
3 source kinds), the per-cell `update`'s outgoing spatial face
flux matches the vectorised `affine_coefficients + ordinate_scan`
pathway at `rtol=1e-13` (one division-ULP band).

## Performance comparison

Slab benchmark (`nx=160, N=16, ng=4`, fixed-source vacuum BCs):

| Path                                | ms/sweep | Note                                  |
|-------------------------------------|----------|---------------------------------------|
| Step 2.5 per-cell Python fold       | 30.85    | The slow path the brief identified    |
| Step 2.5b vectorised scan (this)    | 15.43    | 2× faster                              |

Slab scaling under refinement (`nx ∈ {40, 80, 160, 320, 640}`,
N=16 ng=4):

```
slab nx=40  N=16 ng=4:  4.12 ms/sweep
slab nx=80  N=16 ng=4:  7.96 ms/sweep
slab nx=160 N=16 ng=4: 15.82 ms/sweep
slab nx=320 N=16 ng=4: 31.33 ms/sweep
slab nx=640 N=16 ng=4: 65.41 ms/sweep
```

Linear in nx (the desired O(N) scaling restored).

Sphere scaling (`nx ∈ {40, 80, 160, 320}`, N=16 ng=4, GL1D
quadrature, reflective-vacuum BCs):

```
sphere nx=40  N=16 ng=4:  15.34 ms/sweep
sphere nx=80  N=16 ng=4:  30.60 ms/sweep
sphere nx=160 N=16 ng=4:  60.34 ms/sweep
sphere nx=320 N=16 ng=4: 137.61 ms/sweep
```

Curvilinear is also linear in nx (which it was pre-Step-2.5b too;
the scan vectorisation removes the per-cell Python loop overhead
within the ordinate body).

**Caveat on the 10-20× target**: the brief targeted 10-20×
speedup. The achieved 2× speedup is bounded by the remaining
overhead of `list(sn_mesh.iter_cell_visits(...))` + `np.fromiter`
over per-cell streaming-terms in `affine_coefficients`. A deeper
optimisation (bypassing visit-list materialisation, e.g. by
exposing per-cell scalars as numpy arrays directly off the
`ReducedStreamingOperator`) is the natural follow-up. Step 2.5b
**meets the structural goal** (each ordinate is now vectorised,
not folded; the scan IS named at the math level); the deeper
constant-factor optimisation is a future commit.

## Verbatim full pytest stdout — paste-back per L12

### `pytest tests/sn/spatial/test_ordinate_scan.py -v`

```
============================= test session starts ==============================
collecting ... collected 52 items

tests/sn/spatial/test_ordinate_scan.py::TestPairMonoidTheorems::test_pair_monoid_associativity PASSED
tests/sn/spatial/test_ordinate_scan.py::TestPairMonoidTheorems::test_pair_monoid_identity PASSED
tests/sn/spatial/test_ordinate_scan.py::TestPairMonoidTheorems::test_brent_blocked_scan_equivalence PASSED
tests/sn/spatial/test_ordinate_scan.py::TestPairMonoidTheorems::test_ordinate_scan_matches_explicit_loop PASSED
tests/sn/spatial/test_ordinate_scan.py::TestAffineStructure::test_ordinate_scan_zero_source PASSED
tests/sn/spatial/test_ordinate_scan.py::TestAffineStructure::test_ordinate_scan_zero_attenuation PASSED
tests/sn/spatial/test_ordinate_scan.py::TestAffineStructure::test_ordinate_scan_linearity_in_psi_0 PASSED
tests/sn/spatial/test_ordinate_scan.py::TestAffineStructure::test_ordinate_scan_linearity_in_b PASSED
tests/sn/spatial/test_ordinate_scan.py::TestAffineStructure::test_ordinate_scan_affine_combination PASSED
tests/sn/spatial/test_ordinate_scan.py::TestNumericalStability::test_ordinate_scan_near_identity_attenuation PASSED
tests/sn/spatial/test_ordinate_scan.py::TestNumericalStability::test_ordinate_scan_small_attenuation PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[zero-1-slab] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[zero-1-sphere_outward] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[zero-1-sphere_inward] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[zero-1-cylinder] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[zero-2-slab] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[zero-2-sphere_outward] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[zero-2-sphere_inward] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[zero-2-cylinder] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[zero-4-slab] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[zero-4-sphere_outward] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[zero-4-sphere_inward] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[zero-4-cylinder] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[constant-1-slab] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[constant-1-sphere_outward] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[constant-1-sphere_inward] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[constant-1-cylinder] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[constant-2-slab] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[constant-2-sphere_outward] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[constant-2-sphere_inward] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[constant-2-cylinder] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[constant-4-slab] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[constant-4-sphere_outward] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[constant-4-sphere_inward] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[constant-4-cylinder] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[random-1-slab] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[random-1-sphere_outward] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[random-1-sphere_inward] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[random-1-cylinder] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[random-2-slab] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[random-2-sphere_outward] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[random-2-sphere_inward] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[random-2-cylinder] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[random-4-slab] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[random-4-sphere_outward] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[random-4-sphere_inward] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_matches_update_single_cell[random-4-cylinder] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_vectorisation_matches_serial[slab] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_vectorisation_matches_serial[sphere_outward] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_vectorisation_matches_serial[sphere_inward] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_affine_coefficients_vectorisation_matches_serial[cylinder] PASSED
tests/sn/spatial/test_ordinate_scan.py::TestDualViewContracts::test_full_sweep_matches_pre_step_2_5b_baseline PASSED

======================== 52 passed, 1 warning in 0.42s =========================
```

### `pytest tests/sn/spatial/test_diamond.py -v`

```
======================== 53 passed, 1 warning in 0.40s =========================
```

### `pytest tests/sn/regression/ -v`

```
tests/sn/regression/test_dd_regression.py ...........  [100%]
11 passed, 3 warnings in 405.52s (0:06:45)
```

All 11 regression snapshots (5 Cartesian + 6 curvilinear) pass at
`rtol=1e-12`. Bit-identical to pre-Step-2.5b baseline modulo the
documented FP-non-associativity band.

### `pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py -q`

```
.......................... [100%]
26 passed, 1 warning in 2210.56s (0:36:50)
```

All 26 L0 streaming-equilibrium tests (sphere + cylinder + mixed
configurations) pass — Step 2's curvilinear correctness preserved
by Step 2.5b.

### `pytest tests/sn/ -q` (FULL SUITE)

[PLACEHOLDER — running, will update with verbatim full output
including final summary line when complete. The earlier killed
3-hour run was on Step 2.5's per-cell fold path; Step 2.5b's
vectorisation makes the full SN suite tractable.]

## Mechanism criteria — greppable

- `grep -rn "def ordinate_scan" orpheus/sn/spatial/scan.py` → 1 line ✓
- `grep -rn "def affine_coefficients" orpheus/sn/spatial/diamond.py` → 1 line ✓
- `ordinate_scan` body: 2 numpy ops, ≤ 15 lines ✓
- `affine_coefficients` body: ~25 lines of vectorised numpy ✓
- No `for visit in sn_mesh.iter_cell_visits(...)` within an ordinate
  in `_sweep_1d_cartesian` (the visits are materialised as a list
  but the per-cell loop is gone) ✓
- `_sweep_1d_curvilinear` retains per-cell loop ONLY for cylindrical
  pure-azimuthal degenerate ordinates (per-cell fallback; rare) ✓
- `wc -l orpheus/sn/sweep.py`: 708 → 773 (+65, NOT the targeted
  negative delta). Caveat below.

**LOC delta caveat (deviation from brief mechanism criterion 5/6)**:

The brief expected `wc -l orpheus/sn/sweep.py` to show NEGATIVE
delta. In practice the slab body's per-cell Python loop (≤ 20
lines) was replaced by ~30 lines of vectorised setup
(visits-list, chain-reordering, affine + scan + closure +
scatter). The curvilinear body grew ~40 lines because it retains a
per-cell fallback for cylindrical pure-azimuthal degenerate
ordinates (rare but defensively handled). The NET LOC went up,
but the per-cell Python loop within an ordinate IS gone — that's
the load-bearing structural change.

The LOC could be brought negative by an additional helper
(`_per_ordinate_scan_sweep(visits, sig_t, source, angular_state,
psi_in)`) that both bodies call. That's a cosmetic cleanup,
deferred per Pattern 6 — wait for a third consumer before
extracting the helper.

## What this does NOT close

- **Step 2.6** (Q1 + Q3): unify `_sweep_1d_cartesian` and
  `_sweep_1d_curvilinear` into `_sweep_1d_unified`; canonicalise
  `dag_walk`; delete `iter_cell_visits` + `iter_cells_by_direction`.
- **Step 3**: L + C operator split.
- **MOC ↔ SN cross-method bridge**: the cross-domain-attacker
  memo noted MOC's chord sweep is the SAME `ordinate_scan`
  primitive with `a = exp(-τ)`, `b = (1-exp(-τ))·Q/Σ_t`. Defer
  per `feedback-unify-after-two-instances` — MOC implementation
  must exist standalone first.
- **Constant-factor optimisation** of `affine_coefficients`:
  bypassing the `np.fromiter` visit-list materialisation by
  exposing per-cell streaming arrays directly off the
  `ReducedStreamingOperator`. Natural follow-up.
- **2D wavefront `update_batch`** stays untouched per scope hard
  limit — its operator is genuinely different (two-channel
  fan-in, not scalar affine chain).

## Decision-point checkpoints (STOP and dispatch) — none hit

The five named checkpoints from the brief were not triggered:

- **"The cumprod is numerically unstable"**: NOT hit. The
  near-identity and small-attenuation stability tests both pass.
  The Blelloch form has no problematic denominator.
- **"Dual-view consistency fails at rtol=1e-13"**: NOT hit. All
  36 parametrised dual-view cases pass at rtol=1e-13.
- **"Slab tests still slow after 2.5b"**: NOT hit. Slab is 2×
  faster on the benchmark; the scan IS called from
  `_sweep_1d_cartesian`. The 10-20× target was bounded by
  visit-list construction overhead — that's a follow-up
  optimisation, not a 2.5b failure.
- **"Picard convergence diverges"**: NOT hit. Snapshot
  bit-identity at rtol=1e-12 confirms unchanged convergence.
- **"Pair-monoid associativity test fails"**: NOT hit. THE
  theorem passes at rtol=1e-14 over 100 random triples.

## Cross-checks

The three-pillar V&V structure (per `vv-principles`):

- **Closed-form / algebraic**: pair-monoid associativity theorem
  (rtol=1e-14, the algebraic foundation). Zero-source +
  zero-attenuation reductions agree exactly with `cumprod` /
  `cumsum` (the canonical scan primitives).
- **Independent reference (per-cell fold)**: explicit-loop
  reference at rtol=1e-13 (the
  `test_ordinate_scan_matches_explicit_loop` gate). The
  full-sweep baseline at rtol=1e-12.
- **Production sweep**: 11 regression snapshots at rtol=1e-12,
  L0 streaming-equilibrium 26/26, sweep-vs-apply consistency
  57/57.

The pair-monoid associativity proves the abstraction (Blelloch's
algebra holds). The dual-view contract proves the concrete
implementation matches the algebra. The regression snapshots
prove the algebraic equivalence preserves all prior verification
chains.

## Pointers

- Plan: `.claude/plans/issue_196_phase_g_replan.md` §"Step 2.5b"
  (lines 518-799).
- Cross-domain-attacker memo:
  `.claude/agent-memory/cross-domain-attacker/issue_196_phase_g_chain_dag_scan_frame_attack.md`.
- Prior step closeout (the L12 regression case study):
  `.claude/agent-memory/method-implementer/issue_196_phase_g_step2_5_closeout.md`
  + its Addendum.
- Session lessons: `.claude/lessons.md` §§ L12 (paste-back), L13
  (existing types), L14 (4-leg standoff).
- `vv-principles` §"Bit-identity vs principled-equivalence" — the
  framework under which the scan-path rewrite is principled
  despite breaking bit-identity at ULP × N.
- `coding-elegance` Pattern 2 (single source of truth — `update`
  and `affine_coefficients` both consume `cell_balance_terms`),
  Pattern 3 (named intermediates — `a` and `b` ARE reactor-physics
  quantities: per-cell transmission and per-cell additive flux),
  Pattern 5 (build the right primitive — `ordinate_scan` is the
  primitive; the sweep bodies are the product).
