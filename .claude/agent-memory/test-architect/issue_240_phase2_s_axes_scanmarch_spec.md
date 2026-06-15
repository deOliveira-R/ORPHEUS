---
name: issue-240-phase2-s-axes-scanmarch-spec
description: Verification spec for #240 Phase 2 — the s_axes streaming convention change (2|μ|/Δ → g=|μ|·A_down/V, factor-2 internalised in the DD kernel) coupled with #239 (2-D ScanMarch onto the coefficient model). Crosswalk, re-baseline list, positive same-value gate, polymorphic gate.
metadata:
  type: project
---

# #240 Phase 2 verification spec — `s_axes` convention + 2-D ScanMarch coefficient-model

**Status:** PRE-IMPLEMENTATION spec (carve NOT written). Branch
`feature/sn-space-angle-tier2`, HEAD `944d3b3` (Phase 1 `fde76ac` +
records). Host env, canonical `python -O -m pytest`.

**The carve in one line.** Move the diamond factor-2 out of the SHARED
streaming interface (`SNMesh.streaming`) and INTO the DD cell kernel,
so the interface returns the scheme-agnostic raw down-face streaming
`g = |μ|·A_down/V` and every affine scheme applies its own closure
factor. Coupled: the 2-D ScanMarch inline DD (which assumes `s = 2g`)
must be rewritten onto `affine_scan_coefficients` in the SAME change,
because both consume `SNMesh.streaming`.

**It is value-preserving.** Real-arithmetic identical; the only drift
is FP-non-associativity (~1 ULP) from the changed reduction grouping.
**The failure mode this spec is built to catch is a MISSED producer →
a clean 2× error** (one site passes `g` where the consumer still
expects `2g`, or vice versa). A 2× error is LOUD — the gates below
make it impossible to ship silently.

---

## 0. Source-of-truth map (read live @ HEAD, file:line)

This is the producer/consumer inventory the implementer edits. The
graph (`callers`/`impact`) is a cross-check; **the file:line is
authoritative** (Nexus was mid-rebuild this session).

### Producers of `streaming` (`2|μ|/Δ` today)

| Site | File:line | What it builds | Role |
|------|-----------|----------------|------|
| **CANONICAL** | `geometry.py:1436-1440` `_setup_cartesian` → `_streaming_axes[a] = 2·|μ_a|/Δa` | the per-axis stencil tuple | THE single producer |
| accessor | `geometry.py:738-764` `SNMesh.streaming(axis)` returns `_streaming_axes[axis]` | accessor over canonical | accessor (no math) |
| **DUPLICATE (Pattern-2 smell)** | `loss_representation.py:2038` `_OneDimScanWalk._apply_walk` → `(2.0 * abs_mu / V[i])` | re-derives `2|μ|/Δ` for the 1-D matvec | DEDUP TARGET |

### `str_axes` plumbing (feeds the DAG walks)

| Site | File:line | Consumes |
|------|-----------|----------|
| apply operands | `loss_representation.py:633` `str_axes=tuple(sn_mesh.streaming(a) …)` | feeds `_ApplyOperands` |
| solve operands | `loss_representation.py:2924` same | feeds `_SolveOperands` |
| DAG slice (solve) | `sweep_graph.py:759` `s[..., c]` per level | `cell_kernel_batch` |
| DAG slice (apply) | `sweep_graph.py:713` `s[..., cell_idx]` | `residual_kernel_batch` |
| ScanMarch slice | `loss_representation.py:1294, 1401` `s_x, s_y = (s[oct_idx] …)` | inline DD row-march |

### Consumers of the streaming value `s`

| Consumer | File:line | Uses `s` as | Convention assumed |
|----------|-----------|-------------|--------------------|
| **DD solve kernel** | `diamond.py:336-340` `cell_kernel_batch` — `denom += s_a`, `numer += s_a·in_a`, closure `2ψ̄−in` | the streaming-face denom term DIRECTLY | `s = 2g` (denom term IS `2|μ|/Δ`) |
| **DD apply kernel** | `diamond.py:370-377` `residual_kernel_batch` — same fold | same | `s = 2g` |
| **LD solve kernel** | `linear_discontinuous.py:429` `g = 0.5 * s_axes[0]` | HALVES it back to `g` | `s = 2g` (halves to `g`) |
| **2-D ScanMarch solve** | `loss_representation.py:1313-1322` `D_row = sigt+s_x+s_y`, `alpha = 2·s_x/D_row−1`, `beta = 2·(Q+s_y·ψy)/D` | inline DD, `s` IS denom term | `s = 2g` |
| **2-D ScanMarch apply** | `loss_representation.py:1421-1434` `D_row = sigt+s_x+s_y`, `−s_x·in_x − s_y·ψy`, closure `2ψ̄−in` | inline DD | `s = 2g` |
| **scan.py `_scanmarch_row`** | `scan.py:341-342` `ψ̄ = ½(in_x+out_x)`, `out_y = 2ψ̄−ψy_in` | reads `alpha/beta` (already DD-baked) | downstream of `s = 2g` |
| **1-D apply DUPLICATE** | `loss_representation.py:2038` (see above) | feeds `residual_kernel_batch` | `s = 2g` |

### The DD coefficient-model producer (already factor-2-explicit)

| Site | File:line | Form |
|------|-----------|------|
| `affine_scan_coefficients` | `diamond.py:462` `streaming_face_term = 2.0·|μ|·A_down` | builds `2g·A_down` denom term EXPLICITLY (factor-2 is HERE, not hidden in `s`) |
| consumed by | `sweep_cache.py:464-465` `CollisionCache.from_geometry` | the `CumprodScan`/`ScanMarch`-1D scan substrate |
| generic ops | `affine_closure.py` `source_emission(b=QV·inv/w)`, `cell_average((1−w)in+w·out)` | scheme-agnostic; DD `w=½` |

**KEY: the coefficient-model path (`affine_scan_coefficients` →
`CumprodScan`) already carries the factor-2 EXPLICITLY in its own
denom build (`diamond.py:462`), independent of `SNMesh.streaming`.**
It does NOT read `SNMesh.streaming`; it reads geometry (`abs_mu`,
`A_down`, `V`) directly. So the 1-D CumprodScan / curvilinear scan
paths are UNAFFECTED by the `streaming` convention change. Only the
DAG-walk consumers (`cell_kernel_batch`, `residual_kernel_batch`, LD
kernel) and the 2-D ScanMarch inline DD read `SNMesh.streaming`.

---

## 1. THE CROSSWALK (L17) — `s_axes` / `streaming`

The convention crossed: **the diamond factor-2 location** (shared
interface vs scheme kernel). Per `coding-elegance` Pattern 7 the fix
moves the bridge to the DEFINITION SITE — but here the "convention" is
*who owns the factor-2*, and the principled answer is: the producer
emits the scheme-AGNOSTIC `g`, each scheme applies its OWN factor. The
factor-2 is DD's diamond closure (`ψ_out=2ψ̄−ψ_in ⟹ denom = Σ_t+2g`),
provably NOT geometry — it does not belong in a shared geometric
accessor.

### Crosswalk table (BEFORE → AFTER)

| Subsystem | Input | Internal | Output | Bridge (where the factor-2 lives) |
|-----------|-------|----------|--------|-----------------------------------|
| **PRODUCER** `_setup_cartesian` :1436 | geometry `μ,Δ,A_down,V` | BEFORE: `2·|μ|/Δ`  AFTER: `g=|μ|·A_down/V` (slab: `A_down=1,V=Δ ⟹ |μ|/Δ`) | BEFORE `2g`  AFTER `g` | factor-2 LEAVES here |
| **DD solve kernel** `cell_kernel_batch` :336 | BEFORE `s=2g`  AFTER `g` | `denom += s_a` → AFTER `denom += 2·s_a`; `numer += s_a·in` → AFTER `numer += 2·s_a·in` | `ψ̄`, closure `2ψ̄−in` (unchanged) | factor-2 ARRIVES here (internalised) |
| **DD apply kernel** `residual_kernel_batch` :370 | BEFORE `s=2g`  AFTER `g` | same: AFTER `denom += 2·s_a`, `numer += 2·s_a·in` | residual `denom·ψ̄−numer` | factor-2 ARRIVES here |
| **LD kernel** `_kernel_terms` :429 | BEFORE `s=2g`  AFTER `g` | BEFORE `g = 0.5*s_axes[0]` → AFTER `g = s_axes[0]` (DROP the 0.5) | Schur `eff_denom, rhs` | factor-2 GONE from LD (it never wanted it — was halving) |
| **2-D ScanMarch solve** `_sweep_interior` :1313 | BEFORE `s=2g`  AFTER `g` | REWRITE off inline DD onto `affine_scan_coefficients`+`affine_closure` (consumes the cache's `a, inverse_denom, w`; the `2·QV·inv` source emission becomes `source_emission(QV, inv, w)`) | `ψ̄`, `out_y` | factor-2 now in the cache build (`affine_scan_coefficients` already explicit) |
| **2-D ScanMarch apply** `_loss_action_interior` :1421 | BEFORE `s=2g`  AFTER `g` | matvec rides `residual_kernel_batch` (the AFTER kernel that internalises 2), OR explicit `2·s_x` if kept inline | residual | factor-2 internalised same as DD apply |
| **1-D apply dup** `_apply_walk` :2038 | geometry `μ,V` | BEFORE `2.0*abs_mu/V` → AFTER DEDUP: read `g` from a single-source (`SNMesh.streaming` or the scheme); kernel internalises 2 | feeds `residual_kernel_batch` | DEDUP + factor-2 to kernel |

### Column-for-column post-change verification

Every producer→consumer pair must agree on `g` (raw) AFTER the change:

- `_setup_cartesian` emits `g`. `str_axes:633/2924` pass `g` through
  unchanged. DAG slices `sweep_graph.py:713/759` pass `g`. DD kernels
  multiply by 2 internally → denom term `2g` = BEFORE. ✓
- LD kernel receives `g`, drops the `0.5` → uses `g` directly
  (BEFORE it received `2g` and halved to `g` — same value `g`). ✓
- ScanMarch slices `:1294/:1401` receive `g`; the rewritten interior
  routes through the cache (factor-2 in `affine_scan_coefficients`) OR
  the `residual_kernel_batch` kernel (factor-2 internal). ✓
- 1-D apply dup `:2038` deduped to read `g`; kernel internalises 2. ✓

**The one trap (write it on the implementer's wall):** the DD kernel
fold is `denom = denom + s_a` TODAY. AFTER it MUST be `denom = denom +
2.0 * s_a` AND `numer = numer + 2.0 * s_a * in_a` — BOTH the denom AND
the numerator gain the factor-2, because `s` appears in both (`denom
+= s_a`, `numer += s_a·in_a`). Forgetting the numerator's factor-2 is
a NON-uniform 2× bug (denom doubled, numerator not) → wrong `ψ̄` (not
even a clean global 2×). The closure `2ψ̄−in` is UNCHANGED (that `2` is
the diamond mean, already correct).

---

## 2. RE-BASELINE SCOPE — every strict bit-id gate that drifts

### Green baseline @ HEAD (L12, verbatim)

`python -O -m pytest tests/sn/sweep/core tests/sn/operators -q`:

```
7 failed, 928 passed, 5 skipped, 4 xfailed, 14 warnings in 22.73s
```

The 7 failures are ALL PRE-EXISTING reds (clean working tree re: Phase
2 source — confirmed `git diff HEAD orpheus/sn/` empty). They are NOT
Phase 2's job and must be ROUTED AROUND:

| Failing test | Issue | Mechanism |
|--------------|-------|-----------|
| `test_bc_extraction_matvec.py::TestVacuumMatvecBitIdentity::test_vacuum_bulk_bit_identical_1d[{0,1,2}-SPH]` | #195/#209 | curvilinear sphere stale-post-ERR-058 |
| `test_streaming_operator.py::TestT4cPreT4RegressionSnapshotCurvilinear::test_sphere_{1g,2g}_apply_bit_identical` | #195/#209 | curvilinear sphere snapshot |
| `test_boundary_conditions.py::TestSNBCResolution::test_2d_mesh_resolution` | #214 | default quad has no genuine `mu_y` (`ny=1`-class) |
| `test_native_matvec.py::TestTwoDCartesianRaises::test_two_d_cartesian_loss_action_returns_result` | #214 | same `mu_y==0` rank-mismatch |

**Route-around invocation for the green-subset check:**
`-k "not (vacuum_bulk_bit_identical_1d and SPH) and not (sphere_1g_apply_bit_identical or sphere_2g_apply_bit_identical) and not test_2d_mesh_resolution and not two_d_cartesian_loss_action"`.
NEVER run all of `tests/sn` (#212 `continuous_get` hang).

### What DRIFTS (matvec re-association ~1 ULP) and what does NOT

The decisive fact, established this session from `affine_closure.py`
docstring + `diamond.py` reads:

- **The SCAN solve is byte-identical for DD (w=½)** — `0.5*in+0.5*out
  == 0.5*(in+out)` and `QV·inv/0.5 == 2·QV·inv` bit-for-bit (×/÷ by ½
  is an exact power-of-2 scaling commuting with IEEE rounding). So
  every SOLVE/SWEEP snapshot stays STRICT bit-id.
- **The matvec/apply re-associates** — moving the denom build to the
  `2·s_a` fold changes the addition grouping on non-power-of-2 widths
  → ~1 ULP. Apply snapshots need the Phase-1 pattern.

#### 2a. DD apply/matvec snapshots → migrate to `assert_regression`

The Phase-1 PRECEDENT to copy verbatim is
`tests/sn/operators/test_streaming_operator.py::TestT4bPreT4RegressionSnapshot`
(slab arms ALREADY migrated by Phase 1): bulk →
`assert_regression(..., kind="direct", reduction_depth=mesh.nx)`,
boundary stays STRICT `assert_array_equal` (the outflow defect
reconstructs from the same `2ψ̄−in` faces → 0 ULP). Bit-identity kept
as escalatable `DriftWarning` (`-W error::DriftWarning` re-pins
bytes).

| Gate | File | Treatment | Notes |
|------|------|-----------|-------|
| **`test_cell_kernel_batch.py` hand-calc** `TestSolveKernelClosedForm` / `TestResidualKernelClosedForm` | `tests/sn/sweep/core/test_cell_kernel_batch.py:84,164` | **RE-BASELINE the hand-calc** (see §2c — this is the kernel's PUBLIC contract change, NOT a snapshot) | the `s_axes=(3,5)`→`denom=10` arithmetic CHANGES to `denom=2+2·3+2·5=18` IF inputs stay `s=2g`; see §2c for the correct fix |
| `test_cell_kernel_batch.py` left-fold bit-id | `:122 TestBitIdenticalToPerOrdinateLoop` | UPDATE the `_per_ordinate_loop_reference` (:109) to mirror the NEW `2·s_a` fold; KEEP `array_equal` (vectorisation is still exact vs loop) | the reference must gain the same factor-2 |
| 2-D ScanMarch ≡ oracle (apply leg) | `tests/sn/sweep/cartesian_2d/test_scan_march_equivalence.py:179 test_scanmarch_residual_equals_oracle` | already `assert_allclose` (nulp-tolerant) — SURVIVES | NON-bit-id by design across schedules |
| window ≡ full-field (apply leg) | `tests/sn/sweep/core/test_sweep_graph_window_equivalence.py:239 test_residual_window_equals_full_field` | `array_equal` — SURVIVES IFF both walks read the SAME kernel (they do — both `residual_kernel_batch`); the re-association is IDENTICAL on both paths → 0 ULP between them | this gate is path-A≡path-B, not path-vs-snapshot — re-association cancels |
| end-to-end matvec window≡full | `tests/sn/sweep/cartesian_2d/test_2d_full_field_oracle.py:109 test_matvec_window_equals_full_field_end_to_end` | `array_equal` — SURVIVES (same reason: both paths, same kernel) | |
| `bc_extraction_matvec` cart2d | `tests/sn/operators/test_bc_extraction_matvec.py` (cart2d arms, NOT the SPH reds) | if it pins matvec vs a frozen snapshot/closed form → migrate to nulp; if path-A≡path-B → survives | INSPECT: confirm which (see §2d action) |

**NOTE on the cart2d apply snapshots:** `pre_t4_snapshots.npz` carries
`cart2d_{1g_vacuum,1g_specular,2g_specular}_apply_{bulk,boundary}` and
`cart2d_2g_specular_LpC_apply_*` — BUT they are **NOT currently
consumed by any test** (grep confirms only slab/sphere/cyl arms read
snapshots in `test_streaming_operator.py`). The 2-D Cartesian matvec
is regression-covered by the path-vs-path oracles (window≡full,
ScanMarch≡oracle), which SURVIVE the re-association because both sides
move together. **Action:** the implementer SHOULD add a cart2d apply
arm to `TestT4bPreT4RegressionSnapshot` (the snapshots already exist)
so the 2-D matvec is also pinned against a STRUCTURALLY-INDEPENDENT
frozen reference, not only against its own twin path — a path-vs-path
gate is blind to a bug that moves BOTH paths identically (e.g. both
forget the numerator factor-2). This is the §3 positive gate applied
to 2-D.

#### 2b. SCAN/SWEEP snapshots → STAY STRICT (byte-identical)

These MUST remain `array_equal` and the strict gate
`pytest tests/sn/sweep/core tests/sn/solve -W "error::…DriftWarning"`
MUST stay clean (DD scan is power-of-2-exact):

| Gate | File | Why unchanged |
|------|------|---------------|
| sweep regression | `tests/sn/sweep/core/test_sweep_regression.py` | SOLVE path (scan) — byte-identical |
| affine carve baseline | `tests/sn/sweep/core/test_affine_carve_baseline.py` | the Phase-1 strict sha256/bit-id gate — scan-side |
| group-3≡group-2 scan flat | `test_sweep_cache.py` (the `affine_scan_coefficients` gate) | scan coefficients unchanged (factor-2 was always explicit at `diamond.py:462`) |
| sweep window≡full (SOLVE) | `test_sweep_graph_window_equivalence.py:216 test_solve_window_equals_full_field` | scan solve, path≡path AND byte-exact |
| CumprodScan≡wavefront | `test_wavefront_cumprod_equivalence.py` | 1-D scan — does NOT read `SNMesh.streaming` (rides `affine_scan_coefficients`) → UNTOUCHED |
| `test_diamond.py` slab snapshots | `tests/sn/sweep/core/test_diamond.py` | INSPECT: solve-side stays strict; any apply-side row migrates per §2a |

**The Phase-1 strict gate command (must stay 505p/1s/4xf per my prior
note `project_issue_158_ld_dag`):**
`python -O -m pytest tests/sn/sweep/core tests/sn/solve -W "error::orpheus.sn.regression._regression_assert.DriftWarning"`
— if this goes RED on a SWEEP snapshot, the change leaked into the
solve path (a bug — scan must stay byte-identical). If it goes red on
an APPLY snapshot, that snapshot needed §2a migration.

#### 2c. The kernel hand-calc — the PUBLIC CONTRACT decision (CRITICAL)

`test_cell_kernel_batch.py:84` feeds `s_axes=(3.0, 5.0)` and asserts
`denom = Σ_t + s_x + s_y = 2+3+5 = 10`, `ψ̄ = 68/10 = 6.8`. This is
**not a regression snapshot** — it is the kernel's documented
closed-form contract, and it ENCODES the convention. After the change
the kernel computes `denom = Σ_t + 2·s_x + 2·s_y`. There are two
admissible designs; the implementer MUST pick deliberately and the
test re-baseline differs:

- **Design A (recommended) — the kernel's `s_axes` param means raw
  `g`.** Then the hand-calc inputs become the raw streaming `g`, and
  the expected `denom = 2 + 2·3 + 2·5 = 18`, `ψ̄ = (16 + 2·3·4 +
  2·5·8)/18 = 120/18`. Rewrite the docstring math to
  `ψ̄ = (Q + Σ 2g_a·in_a)/(Σ_t + Σ 2g_a)`. This is the honest contract:
  the kernel param IS the raw streaming, the kernel owns the factor-2.
- **Design B (reject) — keep the kernel param meaning `2g` and only
  change the PRODUCER.** Then the kernel is unchanged and the hand-calc
  stays `denom=10` — BUT now `SNMesh.streaming` returns `g` while the
  kernel still wants `2g`, so SOMETHING between them must double it.
  That "something" is a new bridge at every call site (`sweep_graph`
  slices, ScanMarch) → re-introduces the very Pattern-7 scatter the
  carve exists to remove. **REJECT** — it defeats the carve.

→ **Design A.** The hand-calc re-baseline is a CONTRACT change, not a
golden regeneration: the new numbers are HAND-DERIVED from the new
documented formula (NOT copied from a fresh run — `coding-elegance`
"never transcribe values"). Verify by independent arithmetic in the
test docstring.

#### 2d. Pre-write inspection actions (the implementer runs these first)

1. `grep -n "def test_\|array_equal\|nulp\|allclose\|snapshot" tests/sn/operators/test_bc_extraction_matvec.py`
   — classify each cart2d arm: path-vs-path (survives) or
   path-vs-frozen (migrate per §2a).
2. `grep -n "apply\|residual\|matvec\|array_equal\|snapshot" tests/sn/sweep/core/test_diamond.py`
   — find any APPLY-side strict gate; migrate per §2a (solve-side
   stays strict).
3. Run the strict DriftWarning gate (§2b command) BEFORE editing →
   record the exact `Np/Ns/Nxf` baseline so a post-change delta is
   attributable.

---

## 3. THE POSITIVE SAME-VALUE GATE (the value-preserving-change linchpin)

A convention change that re-associates only ~1 ULP MUST be proven to
produce the SAME value as the old convention to nULP — otherwise a
missed producer (passing `g` where `2g` is expected) yields a clean
2× error that a re-baselined snapshot would happily absorb as "the new
golden". The §2 snapshots are re-baselined and therefore CANNOT catch
a 2× error introduced by the same change. The gate below CAN.

**The strongest existing gate: `test_scanmarch_sweep_equals_oracle`**
(`tests/sn/sweep/cartesian_2d/test_scan_march_equivalence.py:104`).
WHY it is the strongest:

- `FullFieldWavefront.sweep` reads `str_axes` → routes through
  `cell_kernel_batch` (the AFTER kernel, factor-2 internalised).
- `ScanMarch.sweep` reads `str_axes` → routes through the REWRITTEN
  `_sweep_interior` (factor-2 via `affine_scan_coefficients`).
- These are TWO INDEPENDENT consumers of the SAME `g`. If the
  implementer changes the producer but mis-handles the factor-2 in
  EITHER consumer (forgets it in one, applies it twice in the other),
  the two paths disagree by a clean factor — `assert_allclose` at
  `_RTOL/_ATOL` (nulp-class) fires LOUD (a 2× error is ~1e15 ULP, not
  ~1 ULP). This is a Mode-9-style two-paths-must-agree gate, but here
  the degeneracy is BROKEN by construction (the two paths use
  genuinely different reduction trees), so it is a true
  same-value cross-check, not a degenerate coincidence.

**But a two-paths gate is BLIND to a bug that moves BOTH paths
identically** (e.g. both `cell_kernel_batch` AND `affine_scan_coefficients`
emit `1·g` instead of `2·g`). To close that hole, PAIR it with a
gate against a value FROZEN UNDER THE OLD CONVENTION:

**The frozen-reference positive gate (NEW, the implementer adds it):**

1. **2-D matvec vs `cart2d` apply snapshots** (the snapshots exist,
   frozen pre-#240). Add a cart2d arm to
   `TestT4bPreT4RegressionSnapshot` mirroring the slab arm: bulk →
   `assert_regression(kind="direct", reduction_depth=mesh.nx)` (the
   ONLY admissible drift is ~1 ULP; a 2× error blows the nulp bound by
   ~1e15 and fires). This pins the 2-D matvec against a value computed
   under the OLD `2|μ|/Δ` convention — so a uniform-both-paths 2× bug
   is caught.
2. **1-D/slab matvec vs slab apply snapshots** — ALREADY in place
   (`TestT4bPreT4RegressionSnapshot` slab arms, Phase 1). These pin
   the 1-D path against the frozen old value. SURVIVES the change
   (drift ~1 ULP); a 2× error fires. The implementer must CONFIRM
   these still pass after the change (do NOT regenerate the slab
   snapshot — it is the frozen old-convention reference).

**The clincher:** the slab+cart2d apply snapshots are FROZEN under the
old convention. After the change the matvec re-associates ~1 ULP but
the VALUE is preserved → `assert_regression` passes within nulp. If
the implementer drops a factor-2 anywhere on the matvec path, the
value is 2× off → the nulp bound is exceeded by ~15 orders → LOUD red
naming the exact arm. This is the single most important gate in the
spec: **a frozen-old-convention reference + a value-preserving change
= a 2× error is unmissable.**

(Boundary traces stay STRICT `array_equal` — 0 ULP — on all apply
arms; a factor-2 bug would also corrupt the boundary defect, giving a
second independent red.)

---

## 4. THE #239 POLYMORPHIC GATE — DD 2-D ≡ FullFieldWavefront oracle through the genericised ScanMarch

The oracle ALREADY exists. The genericised ScanMarch path is verified
by:

| Gate | File:line | Claim |
|------|-----------|-------|
| **headline** | `test_scan_march_equivalence.py:104 test_scanmarch_sweep_equals_oracle` | `ScanMarch.sweep ≡ FullFieldWavefront.sweep` (nulp) — the genericised row-march ≡ the d-generic full-field kernel |
| matvec leg | `test_scan_march_equivalence.py:179 test_scanmarch_residual_equals_oracle` | `ScanMarch.loss_action ≡ FullFieldWavefront.loss_action` (nulp) |
| moment leg | `test_scan_march_equivalence.py:129 test_scanmarch_moment_equals_window` | moment-output ≡ window (transitively pins) |
| end-to-end | `test_2d_full_field_oracle.py:79,109` | window ≡ full-field, sweep AND matvec |

**Parametrisation requirement (Mode-9 degeneracy-break):** these gates
MUST run on a het / 2G-asymmetric / non-flat / non-square config — a
flat-flux box NULLS the streaming-redistribution and a square box
hides x↔y swaps. Confirm the existing parametrize covers `ng=2`,
`bc=specular`, and `nx≠ny` (read `:104` decorator — the file already
parametrises `(nx, ny, lvl, ng, bc)`; verify at least one tuple is
`nx≠ny, ng=2, bc=specular`). If not, ADD one.

**Quadrature requirement:** the 2-D path needs a genuine-`mu_y`
quadrature (`level_symmetric` / `product`) — the default GL has
`mu_y==0` (the #214 red). The existing 2-D gates already build the
right quad (they pass green @ HEAD — confirmed 56p this session);
the implementer must NOT regress that.

**Scope note for #239 (DD-verified + polymorphic-READY):** the full
2-D LD TWO-PATHS gate (LD via ScanMarch ≡ LD via FullFieldWavefront)
requires the bilinear 2-D LD kernel, which is DEFERRED (#158 Increment
D — `linear_discontinuous.py:422` raises `NotImplementedError` for
`d>1`). So #239 lands:
- DD 2-D ScanMarch **VERIFIED** (the gates above, all DD).
- LD 2-D **polymorphic-READY** — the ScanMarch interior consumes
  `affine_scan_coefficients` generically (no inline DD), so when the
  bilinear LD kernel lands (Increment D), LD rides the 2-D row-march
  with NO further ScanMarch change. Add an `xfail(reason="#158 Inc D —
  bilinear 2-D LD kernel")` placeholder test asserting `ScanMarch.sweep`
  with an LD `cell_update` ≡ FullFieldWavefront, so the readiness is
  TRACKED and flips to xpass when Inc D lands.

---

## 5. STRUCTURALLY-INDEPENDENT REFERENCES THAT MUST STAY GREEN (criterion b)

Per `vv-principles` §"Bit-identity vs principled-equivalence" criterion
2: old-vs-new ULP-distance is necessary but NEVER sufficient. The
change must be cross-checked against references from a DIFFERENT
structural angle. These are value-claims (NOT bit-id) and MUST stay
green — they prove the re-associated value is CORRECT, not merely close
to the old value.

| Reference | File | Claim layer | Pillar | Why it MUST stay green |
|-----------|------|-------------|--------|------------------------|
| **DD analytical k∞ (MULTI-group)** | `tests/sn/verification/analytical/test_kinf_homogeneous.py` | eigenvalue | closed-form (`λ_max(A⁻¹F)`) | the ONLY structurally-independent eigenvalue anchor; 1G is DEGENERATE (Cardinal Rule) — confirm ≥2G arms run |
| DD 2-D MMS convergence | `tests/sn/verification/mms/test_mms_2d.py` | convergence-order + flux-shape | MMS | proves the 2-D operator (post-rewrite) is consistent; O(h²) to the imposed solution |
| LD MMS O(h²) slab | `tests/sn/verification/mms/test_mms_ld_slab.py` | convergence-order | MMS | proves the LD kernel still O(h²) after dropping the `0.5` — the LD slope-sign trap catcher |
| DD 1-D MMS (slab + het) | `tests/sn/verification/mms/test_mms_heterogeneous.py` | convergence-order | MMS | heterogeneous → activates redistribution (H2); flat flux would null it |
| Phase-C crosscheck | `tests/sn/verification/analytical/test_phase_c_crosscheck.py` | flux-shape | semi-analytical (trajectory_resolvent) | structurally-independent of the sweep kernel |

**Live confirmation this session:** `test_kinf_homogeneous.py` +
`test_mms_2d.py` + `test_mms_ld_slab.py` + the 2-D oracles =
**56 passed, 3 xfailed** under `-O`. This is the green floor the change
must preserve.

**The necessity chain (why both §3 AND §5):** §3 (frozen-old-value +
nulp) proves the change is VALUE-PRESERVING (caught a 2× error). §5
(structurally-independent k∞/MMS) proves the PRESERVED value is
CORRECT. Neither alone suffices — a change could preserve a wrong old
value (§5 needed), or match k∞ while a 2× bug compensates elsewhere
(§3 needed). Run BOTH.

---

## 6. COMPLETE-POLYMORPHISM AUDIT CHECKLIST

The carve's architectural goal (per `feedback_principled_over_bit_identical`
+ user directive): NO scheme special-casing, equivalent capability
across the strategy polymorphism. After the change, AUDIT:

### Grep patterns (run all; each hit is a smell to justify or remove)

```bash
# inline discretization constants in sweep strategies (the carve removes these)
grep -rn "0\.5 \* s\|0\.5\*s\|2\.0 \* s_\|2\.0\*s_\|2 \* s_x\|2 \* s_y\|2\.0 \* psi_bar\|2\.0 \* psi_avg\|2\.0 \* QV\|2 \* QV" orpheus/sn/loss_representation.py orpheus/sn/spatial/scan.py
# the duplicate streaming producer (Pattern-2 dedup target)
grep -rn "2\.0 \* abs_mu\|2 \* abs_mu\|2\.0\*np.abs\|abs_mu / V\|abs_mu/V" orpheus/sn/
# scheme-name / isinstance dispatch (must be NONE in sweep strategies)
grep -rn "isinstance.*DiamondDifference\|isinstance.*LinearDiscontinuous\|cell_update.__class__\|matvec_via_kernel\|is_diamond\|== 'dd'\|== \"dd\"\|scheme ==" orpheus/sn/
# the retired Phase-1 flag must be GONE
grep -rn "matvec_via_kernel" orpheus/
```

### Strategies to audit for inline discretization / scheme-name / isinstance

| Strategy | File:class | Audit focus |
|----------|------------|-------------|
| `CumprodScan` | `loss_representation.py` | UNAFFECTED (rides `affine_scan_coefficients`, no `SNMesh.streaming`) — confirm no new inline `2*` |
| **`ScanMarch`** | `loss_representation.py:_sweep_interior:1268`, `_loss_action_interior:1380` | THE rewrite target — confirm NO inline `2*s`, `2*psi_bar`, `2*QV`; all via `affine_scan_coefficients`+`affine_closure` |
| `FullFieldWavefront` | `loss_representation.py` + `sweep_graph.py` walks | rides `cell_kernel_batch`/`residual_kernel_batch` — confirm factor-2 ONLY inside the DD kernel, not in the walk |
| `MovingFrontierWindow` | `loss_representation.py` | same walk family — confirm no inline `2*s` in the window seam |
| `_OneDimScanWalk` | `loss_representation.py:1730`, `_apply_walk:1880` | DEDUP `:2038` `2.0*abs_mu/V` → single-source `g`; kernel internalises 2 |
| `_OctantWalk` | `loss_representation.py` (the shared apply frame) | confirm it passes `str_axes` (now `g`) through unchanged to `_loss_action_interior` |
| `scan.py` helpers | `_scanmarch_row:310`, `_x_scan_faces` | `_scanmarch_row:341-342` inline DD `½(in+out)`, `2ψ̄−in` — should become `affine_closure.cell_average`/`outgoing_face_from_average` if ScanMarch routes through the coefficient model; OR documented as the apply-direction `α=−1,β=2ψ̄` reflection scan (which is scheme-agnostic — the `2ψ̄` is the KNOWN-probe reconstruction, not a DD constant) |

### The `_apply_walk:2038` Pattern-2 dedup (explicit deliverable)

`loss_representation.py:2038` re-derives `2.0 * abs_mu / V[i]`
independently of `SNMesh.streaming`. AFTER the change BOTH should be
`g` from ONE source. Options:
- (a) `_OneDimScanWalk` reads `SNMesh.streaming(0)` like the 2-D path
  (preferred — single producer), OR
- (b) reads the scheme's geometry primitive.
Either way the literal `2.0 * abs_mu / V[i]` must be DELETED and the
kernel internalises the 2. **Gate:** after dedup, the 1-D apply
snapshot (`slab_*_apply_bulk`, §3) must still pass at nulp — a dedup
that changes the value beyond ~1 ULP is a bug.

### Audit pass criterion

The grep for inline-`2*s` / `0.5*s` / `matvec_via_kernel` /
`isinstance.*(DiamondDifference|LinearDiscontinuous)` in
`orpheus/sn/loss_representation.py` + `scan.py` MUST return ZERO hits
in sweep-strategy bodies after the change (the only legitimate `2.0*`
survivors: the diamond CLOSURE `2ψ̄−in` and the known-probe
reconstruction `2ψ̄` in the apply direction — both are the diamond MEAN
reconstruction, not the streaming factor, and are scheme-correct).

---

## 7. MODE-8 RUNTIME-MODE HAZARD (per gate)

Canonical invocation is `-O` (strips bare `assert`). The 2-D oracle +
scan-march files emit a `-O` warning this session (bare asserts
present). Per gate:

- **`test_scan_march_equivalence.py`, `test_2d_full_field_oracle.py`,
  `test_sweep_graph_window_equivalence.py`** — use `np.testing.assert_*`
  / `assert_allclose` (FIRE under `-O`). The headline §4 gates are
  `-O`-safe. ✓
- **`test_cell_kernel_batch.py` hand-calc** — uses
  `np.testing.assert_array_equal` (FIRES under `-O`). ✓
- **Any NEW gate the implementer adds** (the cart2d apply arm, the LD
  xfail placeholder) MUST use `np.testing.assert_*` / `assert_regression`
  / `pytest.fail`, NEVER bare `assert` — so it fires under the
  canonical `-O`. Bare-`assert`-bearing gates MUST be run WITHOUT `-O`
  or rewritten (vv Mode 8).
- The `affine_carve_baseline` strict DriftWarning gate uses
  `assert_regression` (fires under `-O`). ✓

---

## 8. DELIVERY CHECKLIST (the implementer's gate sequence)

1. Run §2d inspections + §2b strict-gate baseline (record `Np/Ns/Nxf`).
2. Apply the change (producer → `g`; DD kernel internalises 2 in BOTH
   denom AND numer; LD drops `0.5`; ScanMarch rewrite; `:2038` dedup).
3. **§3 positive gate** — slab apply snapshots pass at nulp (confirm,
   do NOT regenerate); ADD cart2d apply arm, passes at nulp.
4. **§4 polymorphic gate** — `test_scanmarch_sweep_equals_oracle` +
   `_residual_equals_oracle` green on het/2G/non-square; add LD xfail
   placeholder.
5. **§5 references** — `test_kinf_homogeneous` (≥2G) + `test_mms_2d` +
   `test_mms_ld_slab` green (the 56p floor).
6. **§2b strict gate** — `tests/sn/sweep/core tests/sn/solve -W
   "error::…DriftWarning"` stays at the recorded baseline (scan
   byte-identical).
7. **§2a re-baseline** — hand-calc (Design A, hand-derived numbers) +
   any apply snapshot migrated to `assert_regression`; boundary stays
   strict.
8. **§6 audit** — zero inline-`2*s`/scheme-name/isinstance hits in
   sweep-strategy bodies; `:2038` dedup done.
9. Route around the 7 pre-existing reds (§2 invocation); NEVER run all
   of `tests/sn` (#212).

---

## Cross-links

- Extends [[issue_158_ld_spatial_verification]] (the Inc-B affine-scan
  coefficient model — §6 B4 sign-trap gate is the sibling of §3 here;
  the 2-D ScanMarch lift is Inc-B's deferred bonus, landing in #239).
- Precedent [[feedback_regression_tolerance_design]] (`assert_regression`
  kind=direct reduction_depth, DriftWarning tripwire — the Phase-1
  pattern §2a copies).
- Precedent [[issue_206_phase_c_verification]] (density-vs-scan denom
  nULP precedent; the `_apply_walk` 1-D matvec single-source).
- Precedent [[snapshot_migration_when_production_goes_bare]] (vacuum
  bit-id correctness gate; frozen-reference inheritance needs a
  structurally-independent anchor — §3+§5 here).
- `coding-elegance` Pattern 7 (normalise at definition site — the
  factor-2 belongs to the scheme, NOT the shared geometric accessor)
  + the convention-crosswalk template (§1).
- `vv-principles` §"Bit-identity vs principled-equivalence" (the
  three-criteria gate — §3=criterion-3 drift-bound, §5=criterion-2
  structural-independence, §2a=narrowed contract + DriftWarning).
