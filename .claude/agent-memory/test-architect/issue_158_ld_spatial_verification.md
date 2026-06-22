---
name: issue-158-ld-spatial-verification
description: Verification plan for the Linear-Discontinuous (LD) SN cell-update scheme (#158). Increment A (LANDED — DAG-oracle gates) + Increment B (PRE-IMPL — the affine-scan bonus that flips is_affine_scannable→True so LD rides CumprodScan/ScanMarch). Two-paths agreement (CumprodScan-LD ≡ FullFieldWavefront-LD, principled-equiv nulp), routing-test INVERSION, DD strict bit-identity pin, ×V slope-row sign-trap gate.
metadata:
  type: project
---

# Issue #158 — Linear-Discontinuous (LD) spatial scheme + dispatch-seam verification

> **ORIENTATION (2026-06-14).** Increment A has LANDED on
> `feature/sn-space-angle-tier2` (HEAD `3f0cbc2`; LD module
> `orpheus/sn/spatial/linear_discontinuous.py`, gate-1
> `tests/sn/spatial/test_linear_discontinuous.py`, MMS gate
> `tests/sn/verification/mms/test_mms_ld_slab.py`). The §below "INCREMENT A
> SHIPPED" reconciles what landed vs the original plan. The Increment-B
> section (the affine-scan bonus) is the CURRENT pre-impl deliverable — read
> it FIRST if you are gating Increment B. The original pre-impl plan (the
> `OneDimPerCellWalk` seam, 5 gates) is PRESERVED below as history but was
> SUPERSEDED: LD landed on the existing `FullFieldWavefront` DAG oracle
> (group-2 `cell_kernel_batch`), NOT a new `OneDimPerCellWalk` rep — the
> `OneDimPerCellWalk` attempt was reverted per the user. Increment B does NOT
> resurrect it; it generalizes the group-3 (scan) contract so LD ALSO rides
> the fast `CumprodScan`/`ScanMarch`.

---

# ════════════════════════════════════════════════════════════════════
# INCREMENT A — SHIPPED (reconciliation; live evidence 2026-06-14)
# ════════════════════════════════════════════════════════════════════

What landed vs the original 5-gate plan (read before extending):

- **The dispatch seam is NOT `OneDimPerCellWalk`.** LD is
  `is_affine_scannable=False` and routes to `FullFieldWavefront` (the
  dimension-generic DAG oracle, incl. 1-D slab) via its group-2 batched
  kernel `cell_kernel_batch`/`residual_kernel_batch` (the ÷V `g=|μ|/Δ`
  form). The GAP my original plan described ("a 1-D non-scannable mesh is
  admitted by NOTHING") was closed by `FullFieldWavefront.supports` already
  admitting any-d Cartesian incl. 1-D — so no new rep was needed. **Original
  gates 1, 3, 5 shipped; gate 2 (the `OneDimPerCellWalk` seam control) is
  N/A (no such rep); gate 4 (#233 pole) stays deferred (curvilinear LD
  unpublished — `update`/`residual` raise on a curvilinear visit).**
- **The DRIVER GAP is CLOSED.** `solve_sn_fixed_source` now threads
  `cell_update` (`solver.py:107/1868` → `SNMesh(..., cell_update=...)`
  `:134/1972`). All Increment-A + Increment-B tests drive LD END-TO-END
  through the high-level driver — NO direct `SNMesh`/`SNSolver`
  construction needed (my original skeleton's workaround is obsolete).
- **Gate 1 (round-trip + linear-exactness) shipped + STRENGTHENED.** The
  shipped `tests/sn/spatial/test_linear_discontinuous.py` adds a
  `TestLDLinearExactness` oracle I did not spec: LD reproduces ψ=a+bx to
  machine precision (cell-avg AND outflow). This is the
  structurally-independent correctness gate that catches the LM-1989
  §1.4/§6 slope-row SIGN TRAP — stronger than the round-trip alone (the
  round-trip is blind to a sign error consistent between solve+apply).
  KEEP this; Increment B leans on it.
- **Gate 3 (MMS O(h²)) shipped on the SHIPPED `build_1d_slab_mms_case`
  (sin ansatz), NOT my `SNSlabLDStressMMSCase` stress ansatz.** The shipped
  gate asserts `orders[-1]>1.95`, `all>1.85`, value band `1e-9<err<1e-2`.
  **The Mode-7 ansatz-bias OVERRIDE is STILL OWED** — see Increment-B §B5
  (the two-paths gate is where the stress ansatz now earns its keep; a
  flat/sin two-paths gate is Mode-9-degenerate-blind).
- **The `is_affine_scannable=False` → group-2-only fact is the Increment-B
  pivot.** LD ships group 1 (per-cell `update`/`residual`) + group 2
  (batched `cell_kernel_batch`/`residual_kernel_batch`) but NOT group 3
  (scan triple — it raises by default). Increment B ADDS group 3.

**Live baselines (worktree of `feature/sn-space-angle-tier2` @ `3f0cbc2`,
2026-06-14, `.venv/bin/python -O -m pytest`):**
```
tests/sn/spatial/test_linear_discontinuous.py tests/sn/verification/mms/test_mms_ld_slab.py
  → 21 passed, 1 xfailed in 6.54s
  (the 1 xfailed = test_ld_thick_diffusive_limit_xfail — the Increment-C
   diffusion-limit tripwire; stays xfail through Increment B.)

tests/sn/sweep/core tests/sn/solve -W "error::tests.sn.regression._regression_assert.DriftWarning"
  → 505 passed, 1 skipped, 4 xfailed in 81.29s
  (the DD strict bit-identity gate — Increment B's negative control.)
```

---

# ════════════════════════════════════════════════════════════════════
# INCREMENT B — the affine-scan bonus (PRE-IMPL — CURRENT deliverable)
# ════════════════════════════════════════════════════════════════════

**SUT (the carve):** generalize the CellUpdate "group 3" (the DAG-free
affine-scan triple) so LD ALSO rides `CumprodScan`/`ScanMarch`,
polymorphically, exactly as DD does. After B,
`LinearDiscontinuous.is_affine_scannable` flips `False→True`, and a 1-D
slab LD mesh routes to `CumprodScan` (the fast path) instead of
`FullFieldWavefront` (the oracle). LD IS provably affine-scannable:
`ψ_out = a·ψ_in + b` with `a` (transmission) source-independent to ε
(verified from the ÷V `_kernel_terms`: source-free `ψ_out/ψ_in` is a
function of `g, θ, Σ_t` only — no `Q`, no `ψ_in` magnitude).

**The four production touch-points B carves (read from the LD-branch
source):**
1. `LinearDiscontinuous.affine_scan_coefficients(...) -> (a_attenuation,
   inverse_denom)` — NEW. The σ_t-epoch LD transmission, consumed by
   `CollisionCache.from_geometry` (`sweep_cache.py:457`). DD's form is
   `a = 2|μ|·A_total·inv_denom − 1` (`diamond.py:462`); LD's is the
   Schur-reduced outflow-over-inflow ratio — NOT the DD form. **THE
   sign-trap habitat** (the LM-1989 slope-row sign re-expressed in the ×V
   cache convention).
2. `LinearDiscontinuous.cell_average_from_faces(face_in, face_out)` — NEW.
   DD's `½(in+out)` (`diamond.py:480`); LD's is the weighted face mean
   `(1−w)·ψ_in + w·ψ_out` with a σ_t-epoch per-cell weight `w` (DD = the
   `w=½` special case).
3. `LinearDiscontinuous.outgoing_face_from_average(cell_avg, face_in)` —
   NEW. DD's `2ψ̄−ψ_in` (`diamond.py:497`); LD reconstructs the outflow
   from (ψ̄, ψ_in) via the slope row.
4. The source emission `b` — currently INLINED at
   `loss_representation.py:2556` (slab joint-batch,
   `b = 2.0·QV_full_chain·inv_denom_chain`) AND `:2718` (curvilinear
   per-ordinate, `b = 2.0·(QV_chain+ang_contrib)·inv_denom_p`). B moves
   this into a scheme method (DD's factor-2 form = the special case). **NB
   the curvilinear `:2718` site carries the M-M `ang_contrib` — slab LD
   does NOT touch it, but the carve MUST NOT perturb its DD reduction (it
   is inside gate 5's strict-bit-id surface).**

**SOLVER-SIDE SEAM (the non-obvious blast point, `solver.py:915`):** after
B, LD (`is_affine_scannable=True`) takes the `if … is_affine_scannable`
branch → `CollisionCache.from_geometry(geom_cache, sig_t_1d,
cell_update=LD)` builds the LD scan cache (calls
`LD.affine_scan_coefficients` + `cumprod_a`). The source-prep for the scan
must build LD's **two-moment** source (avg + slope) — the scan body
consumes `b` which now needs the slope-moment source `ŝ`. This is the seam
where the Increment-A "flat-source" `Q̂=0` cut meets the scan: B is still
flat (`Q̂=0` — Increment C threads the slope source), so the scan `b` uses
only the average-moment source, and LD's `affine_scan_coefficients`/`b`
encode the Q̂=0 Schur form. **Gate B6 pins exactly this: the scan path
flat-source LD ≡ the DAG-oracle flat-source LD.**

---

## L14 four-leg standoff — how the Increment-B gates instantiate it

The carve is a NEW production path for the SAME LD discretization (scan
schedule instead of DAG schedule). The legs:

| Leg | Claim layer | Increment-B gate |
|-----|-------------|------------------|
| (1) structurally-indep reference | convergence-order (math, L1) | **B5** CumprodScan-LD MMS O(h²) on the STRESS ansatz (the new path independently matches the manufactured ref) |
| (3) scan ≡ DAG-oracle (the headline) | software invariant (foundation) | **B6** two-paths: `solve_sn_fixed_source(LD)`-via-CumprodScan ≡ -via-FullFieldWavefront, principled-equiv nulp |
| (3) scan-matvec ≡ scan-sweep | software invariant (foundation) | **B7** LD Krylov (scan `loss_action`) ≡ LD SI (scan sweep) on the scan path |
| no-regression (negative control) | bit-identity inheritance | **B1** DD strict bit-identity survives the group-3 generalization + b-move |
| routing | selection contract | **B2** routing INVERTS (LD→CumprodScan); **B3** DD still→CumprodScan/ScanMarch |
| (×V sign-trap) | software invariant (foundation) | **B4** group2(÷V) ≡ group3(×V scan coeffs) at Q̂=0 — the analogue of Increment-A's group1≡group2 |
| (4) refinement | built INTO B5 (the ladder) | the convergence ladder |

Claim-layer gate: B5/B6/B7 are convergence-order / software-invariant
claims. **NO eigenvalue gate** — LD is a spatial-closure change; MMS does
NOT prove eigenvalues (vv hierarchical taxonomy). The homogeneous `k_inf`
eigenvalue gates (DD path) are untouched and stay green via B1.

---

## GATE B1 — DD strict bit-identity SURVIVES the carve (negative control)

**Claim:** moving `b` out of the scan body and generalizing
`cell_average_from_faces`/`outgoing_face_from_average`/the scan denom MUST
NOT perturb DD's IEEE-754 reduction order. DD is the `w=½`, factor-2
special case — its `cell_average_from_faces`, `affine_scan_coefficients`,
and source emission must stay byte-for-byte.

**The pin (NO new test file — this IS the invocation, confirmed live
2026-06-14 at the LD-branch HEAD = 505/1/4):**
```
.venv/bin/python -O -m pytest tests/sn/sweep/core tests/sn/solve \
  -W "error::tests.sn.regression._regression_assert.DriftWarning"
```
After the Increment-B landing this MUST still be **505 passed / 1 skipped /
4 xfailed** (re-confirm the 4 xfailed ids are the SAME pre-existing standing
reds — #206 cyl-matvec etc. — NOT newly introduced).

**Which snapshots specifically (the DD-bit-id-bearing tests in that
surface):** `tests/sn/sweep/core/test_sweep_regression.py` (the slab/sphere/
cyl DD scan snapshots, pinned `rtol=1e-12`) — these exercise the EXACT three
sites the carve touches: the `b` emission (the source half of the scan
recurrence), `cell_average_from_faces` (the `½(in+out)` fold at
`loss_representation.py:2570/2735`), and `affine_scan_coefficients` (the
`CollisionCache` σ_t stratum). If the carve regroups DD's ops (e.g. folds the
factor-2 differently, or reorders the weighted-mean), the `DriftWarning`
escalates to error HERE. **Why NOT `tests/sn/regression/`:** its
`conftest.py:pytest_configure` forces `simplefilter("always", DriftWarning)`
+ an `always::…DriftWarning` ini line → shadows `-W error`, INERTING the
escalation there. The strict gate runs on `sweep/core`+`solve` (no override —
FINDING-1, [[issue_206_phase_c_verification]]).

**DISCIPLINE for the implementer (coding-elegance Pattern 7, single source):**
DD's `b` emission factor-2 and LD's general `b` MUST be ONE scheme method
(e.g. `cell_update.scan_source_emission(QV, inverse_denom, ...)` with DD
returning `2·QV·inv_denom`). The scan body at `:2556`/`:2718` calls it; DD's
override reproduces the byte-exact `2.0 * QV_full_chain * inv_denom_chain`
order. **The carve MUST preserve the `2.0 * X * Y` operation order in DD's
override verbatim** (IEEE-754 non-associativity — `2.0*(X*Y)` ≠ `(2.0*X)*Y`
at ULP; read the live line, match it). This is the bit-identity trap; B1 is
its catcher.

---

## GATE B2 — Routing INVERSION (the test that must be UPDATED, not added)

**Claim:** after B, `default_for(LD 1-D slab mesh)` returns `CumprodScan`
(was `FullFieldWavefront` in Increment A).

**⭐ THE INCREMENT-A TEST THAT MUST BE UPDATED (NOT added-alongside):**
`tests/sn/verification/mms/test_mms_ld_slab.py::test_ld_slab_mesh_routes_to_full_field_wavefront`.
It currently asserts (lines 55-59):
```python
if not isinstance(default_for(ld_mesh), FullFieldWavefront):
    pytest.fail("LD slab mesh routed to ... expected FullFieldWavefront")
```
After B this assertion INVERTS — LD now routes to `CumprodScan`. **Rename
the test** `test_ld_slab_mesh_routes_to_full_field_wavefront` →
`test_ld_slab_mesh_routes_to_cumprod_scan` and flip the expected type to
`CumprodScan`. Also drop the now-unused `FullFieldWavefront` import (or keep
it for B6's two-paths gate, which constructs both reps explicitly). **Do NOT
leave the old test asserting the old routing** — it would be a FALSE gate
(asserting a contract the carve deliberately reverses). The
clean-before-extend rule (retire-as-you-go): the routing assertion is a
behavioral contract that the successor (`CumprodScan`) now owns; rewire it,
do not duplicate.

The DD-still-selects half of that same test (lines 60-64,
`default_for(dd_mesh) → CumprodScan`) stays as-is (DD routing unchanged).

**Companion mutual-exclusion (NEW, in the same test — vv anti-pattern #11:
a contract gate needs BOTH a must-select and a must-NOT-select):**
```python
# after B: LD is affine-scannable; FullFieldWavefront must NO LONGER be
# the FIRST supporting rep, but it must still ADMIT the LD mesh (any-d
# Cartesian). The discriminator is default_for ORDER, not supports().
assert CumprodScan.supports(ld_mesh).ok        # NEW: LD now supported by the scan
assert FullFieldWavefront.supports(ld_mesh).ok # still true (the spine admits it)
# the INVERSION proof: CumprodScan is registered FIRST, so it wins.
```
NB: this differs from the Increment-A mutual-exclusion (where
`CumprodScan.supports(ld_mesh)` was `False`). After B it becomes `True` —
that is the whole point. The selection is decided by `LOSS_REPRESENTATIONS`
registry ORDER (CumprodScan first), not by mutual exclusion of `supports`.

---

## GATE B3 — DD routing no-regression (the registry ordering is undisturbed)

**Claim:** generalizing the group-3 contract did not shift DD's selection.
`default_for(DD 1-D slab) == CumprodScan`; `default_for(DD 2-D Cartesian)
== ScanMarch`.

**Pinned by (existing):**
`tests/sn/verification/mms/test_mms_ld_slab.py::test_ld_slab_mesh_routes_to_*`
(the DD half) +
`tests/sn/sweep/core/test_unified_sweep_dispatch.py::TestDispatchSelectsStrategy`
(if present — confirm it still pins DD slab/sphere/cyl→CumprodScan and 2-D→
ScanMarch). **No NEW test needed if those already cover the DD selection** —
B3 is the assertion that they STAY green. The B1 strict-bit-id snapshots
ALSO catch a DD misrouting (a DD mesh that fell off CumprodScan would change
its numerical reduction tree → DriftWarning), so B3 is belt-and-braces at
the selection layer.

---

## GATE B4 — group2(÷V) ≡ group3(×V scan coeffs) at Q̂=0 (the sign-trap gate)

**Claim:** software invariant (foundation). LD's NEW group-3 scan triple
(`affine_scan_coefficients` + `cell_average_from_faces` +
`outgoing_face_from_average`, the ×V `CollisionCache` convention) computes
the SAME LD as the LANDED group-2 batched kernel (`cell_kernel_batch`, the
÷V `g=|μ|/Δ` convention), at the flat `Q̂=0` slice. **This is the EXACT
analogue of Increment-A's `test_group1_equals_group2_flat`
(`test_linear_discontinuous.py:285`) — which cross-checks group1(×V per-cell)
≡ group2(÷V kernel). B4 adds the third leg: group2(÷V) ≡ group3(×V scan).**

**Why this is THE sign-trap catcher (the brief's (e) requirement):** the
LM-1989 §1.4/§6 slope-row SIGN TRAP re-appears every time LD's slope
algebra is re-expressed in a new convention. Increment A caught the per-cell
(×V) and kernel (÷V) instances via linear-exactness + group1≡group2. The
group-3 `affine_scan_coefficients` is a THIRD expression of the slope row
(the Schur transmission `a`), in the ×V cache convention — a fresh habitat
for a sign flip / factor drift. A wrong sign in LD's scan `a` would make the
scan transmission diverge from the kernel transmission. **The Q̂=0 slice is
load-bearing:** the slope SOURCE (Q̂≠0) is Increment C; at Q̂=0 the slope
UNKNOWN ψ̂ is still solved (that is what delivers O(h²)), so the transmission
`a` and the source-free `b` ARE exercised — the sign trap on the UNKNOWN
path is fully visible. (The slope-SOURCE sign — the other half of the trap —
stays untested until Increment C; gate B4 documents this scope.)

**New test:**
`tests/sn/spatial/test_linear_discontinuous.py::TestLDKernel::test_group2_equals_group3_scan_flat`
(co-locate with the existing `test_group1_equals_group2_flat`). For a slab
cell with `(|μ|, h, Σ_t, ψ_in)`:
- group 2 (÷V): `cell_kernel_batch(psi_in, s_axes=(2|μ|/h,), sigt, Q_cells)`
  → `(psi_avg2, psi_out2)`.
- group 3 (×V scan): build the affine recurrence by hand from the triple —
  `a_atten, inv_denom = LD.affine_scan_coefficients(abs_mu, A_down=1,
  A_total=2, dA_w=0, c_out=0, V=h, sig_t)`; `b = LD.scan_source_emission(
  QV=Q̄·h, inverse_denom=inv_denom, ...)` (the moved emission, flat Q̂=0);
  one-cell recurrence `psi_out3 = a_atten·psi_in + b`; recover
  `psi_avg3 = LD.cell_average_from_faces(psi_in, psi_out3)`.
- assert `psi_avg2 ≡ psi_avg3` AND `psi_out2 ≡ psi_out3`.

**TOLERANCE:** `rtol=1e-12, atol=1e-13` (the SAME as the Increment-A
group1≡group2 gate). These are TWO algebraic forms of the same per-cell LD
in real arithmetic; the ULP drift is a handful of divisions deep. NOT
`array_equal` — the ÷V and ×V forms group `Σ_t·h` vs `Σ_t` differently (the
exact density-vs-scan denom-grouping precedent from #206 Phase C). `n_groups
∈ {1,2}` (1G control + 2G real gate — a group-coupling sign drift in the ×V
denom must be observable; vv H1). `@pytest.mark.foundation`.

**Failure mode coverage:** this is the `mode #1 sign flip` + `mode #6
convention drift` catcher for the ×V scan `a`. A sign flip in LD's scan
transmission makes `psi_out3 ≠ psi_out2` (the per-cell magnitudes diverge);
a missing/extra factor in the ×V denom makes `psi_avg3 ≠ psi_avg2`.

---

## GATE B5 — CumprodScan-LD MMS O(h²) on the STRESS ansatz (leg 1)

**Claim:** convergence-order (math, L1) + flux-shape. The NEW fast scan
path independently matches the manufactured reference at O(h²). MMS = the
math/flux-shape pillar (NOT eigenvalue).

**The Increment-A MMS gate
(`test_mms_ld_slab.py::test_sn_1d_slab_ld_mms_converges_second_order`) runs
on `FullFieldWavefront`-LD** (Increment A's only path). After B, the DEFAULT
`solve_sn_fixed_source(LD)` routes to `CumprodScan` — so that EXISTING test,
unchanged, now exercises the SCAN path automatically. **This is the cleanest
B5: the Increment-A MMS gate transparently re-targets to the scan after the
routing flip.** Re-confirm it stays green (orders[-1]>1.95, all>1.85, value
band `1e-9<err<1e-2`) — if the scan `b`/transmission carries a bug, the
order or value band breaks.

**⭐ BUT: the Increment-A MMS uses the SHIPPED `build_1d_slab_mms_case` (sin
ansatz) — the Mode-7 simplification bias my original plan flagged and the
shipped code did NOT override.** For Increment B this is the moment to land
the override, because B6 (two-paths) on a flat/sin problem is Mode-9-
degenerate-blind (a smooth single-scale sin field can make the scan and DAG
reduction trees accidentally close). **NEW test (the stress companion):**
`tests/sn/verification/mms/test_mms_ld_slab.py::test_sn_1d_slab_ld_mms_stress_converges_second_order`.
Markers `@pytest.mark.l1 @pytest.mark.slow @pytest.mark.verifies("ld-cartesian-1d",
"ld-slab")`. Drives `solve_sn_fixed_source(LD)` (→ CumprodScan after B) on
the STRESS ansatz; ladder `n_cells=[20,40,80,160]`; assert `orders[-1]>1.95
∧ all>1.85` + value band.

### THE STRESS ANSATZ (the simplification-bias override — owed since Inc. A)
```
ψ_{n,g}(x) = [ a0_g + a1_g·sin(πx/L) + a2_g·cos(3πx/L)
               + μ_n·( b0_g + b1_g·sin(2πx/L) ) ] / W
```
with **a0_g > 0 load-bearing** (NON-vanishing at both boundaries → real
vacuum-inflow BC that the LD boundary closure must satisfy with a non-trivial
interior; the sin ansatz vanishes at the boundary and leaves BC handling
untested — [[phase4_46_nonvacuum_mms_ansatz]]). Concretely:
- `k₁=1, k₂=3, k₃=2` → **mixed spatial scales** (low + high harmonic) so the
  coarse-mesh O(h²) error is high-frequency-dominated → STRESSES the closure
  + gives the ladder real dynamic range (NOT floor-limited; the shipped sin
  gate is "mildly pre-asymptotic on the coarse pair" per its own docstring —
  the stress ansatz is what makes the coarse pair clear 1.85 honestly).
- **angularly non-trivial** via `μ_n·(b0+b1·sin)` — even though slab has no
  angular redistribution, the per-ordinate streaming `μ·∂ψ/∂x` now carries
  TWO spatial scales per ordinate, so LD's cell-average vs slope moments see
  a genuinely μ-dependent field (NOT the isotropic-flat field sin produces).
- **≥2 groups, HETEROGENEOUS Σ_t(x)** (reuse `_default_hetero_xs_functions`,
  `derivations/continuous/mms/sn.py`, Σ_a>0 guaranteed) + non-trivial
  per-group `a`/`b` coeffs so downscatter coupling is exercised (mode #6).

**Source = SymPy (algebra-of-record Branch 1; NEVER hand-transcribe):** a
new `SNSlabLDStressMMSCase` in `derivations/continuous/mms/sn.py`. Substitute
the ansatz into `μ∂_xψ + Σ_tψ = (1/W)(Σ_s^Tφ + Q^ext)`, solve for `Q^ext`.
`φ_g = ∫ψdμ` → the `μ·(…)` term integrates to zero over symmetric GL → `φ_g
= a0+a1·sin+a2·cos`; the streaming derivative carries the full ansatz. Pin
the symbolic source with a `@pytest.mark.foundation` derive-test, then
consume it. (NB the slope-moment source Q̂ of this manufactured ψ is NON-zero
— but Increment B is flat-source LD, so the MMS must supply Q^ext to the
solver as the per-ordinate average source ONLY, matching the flat-LD
contract; the O(h²) holds because the slope UNKNOWN is solved. This is the
same flat-source posture as the shipped sin gate, which already converges
O(h²) — confirm the stress gate does too. If the stress gate reveals an
order DROP below 2 on the flat cut that the sin ansatz masked, that is a
genuine Increment-B finding → file an ERR + decide whether it gates B or
defers to C.)

**vv §5 (rate ≠ correctness):** assert BOTH the rate AND the value band
against the imposed `phi_exact` (the structurally-independent reference IS
the manufactured solution — that is what makes MMS a flux-shape pillar).

---

## GATE B6 — TWO-PATHS AGREEMENT (the headline; CumprodScan-LD ≡ FullFieldWavefront-LD)

**Claim:** software invariant (foundation), the L14 leg-3 standoff. The NEW
fast scan path and the LANDED, already-trusted DAG oracle are the SAME LD
discretization via two different SCHEDULES — they must agree on the scalar
flux.

**The two paths are NOT bit-identical (justify the tolerance regime — the
brief's (a)+(d)):** the scan path
(`CumprodScan` → `_OneDimScanWalk.sweep` → `ordinate_scan` cumprod/cumsum)
builds `denom` via the ×V `CollisionCache.inverse_denom` convention
(`2|μ|·A_down + Σ_t·V` grouping, `sweep_cache.py:351`) and evaluates the
recurrence as a closed-form parallel-prefix `cumprod·(ψ_0 + cumsum(b/cumprod))`
(Blelloch §1.5, `scan.py:31`). The DAG oracle
(`FullFieldWavefront` → `cell_kernel_batch`) builds `denom` via the ÷V
`g=|μ|/Δ` convention (`linear_discontinuous.py:414-418`) and evaluates a
sequential anti-hyperplane walk (per-level explicit fold). The two compute
the SAME LD value in real arithmetic but accumulate over DIFFERENT IEEE-754
reduction trees (density-vs-scan denom grouping × parallel-prefix-vs-walk
accumulation). Addition/multiplication are non-associative → they disagree at
ULP. Per `vv-principles` §Bit-identity vs principled-equivalence, ALL THREE
criteria hold:
1. **principled** — each intermediate (`inverse_denom`, the per-cell
   `cell_average_flux`, the scan `cumprod_a`) is a named, inspectable
   quantity; the scan denom IS the LD effective denominator, the DAG denom
   IS the same denominator regrouped ÷V.
2. **structurally-independent reference** — the agreement gate is PAIRED
   with B5 (BOTH paths independently match the manufactured MMS flux). The
   two-paths nULP-distance alone is necessary-not-sufficient; B5 supplies
   the structural ground (the imposed solution).
3. **drift dimensionally explained** — single-step solve over the cell
   chain → drift bounded by `(reduction depth) × ULP`. The chain is `nx`
   cells deep AND the SI/Krylov outer iterates ~`O(iters)` times → bound
   `≈ (nx · iteration_count) × condition × ULP`.

**TOLERANCE DECISION (justified):** because the comparison is on the
END-TO-END converged scalar flux (an iterative solve, not a single sweep),
use `np.testing.assert_allclose(rtol=SAFETY × inner_tol, atol=...)` with
`inner_tol=1e-11` and `SAFETY≈100`, i.e. `rtol≈1e-9` — the SAME
"iterative→`SAFETY×conv_tol`" regime [[feedback_regression_tolerance_design]]
and the SAME tolerance the Increment-A `test_sn_1d_slab_ld_mms_krylov_matches_si`
already uses (`rtol=1e-9, atol=1e-11`) for an analogous two-path
(SI≡Krylov) end-to-end agreement. **Do NOT use `array_equal` or `nulp`**
on the END-TO-END flux — the iterative solve makes the FP-drift accumulate
beyond a fixed ULP count (the per-SWEEP scan≡DAG drift IS a fixed nULP, but
the converged-fixed-point drift rides the iteration count). If a
SINGLE-SWEEP agreement is also wanted (tighter, cleaner), add it separately:
one `transport_sweep` through each rep on the SAME `(Q, sig_t, boundary)` →
`assert_array_almost_equal_nulp(nulp=nx)` (the per-sweep fixed-depth drift —
the #206 Phase-C density-vs-scan precedent). **Recommend BOTH:** the
single-sweep nULP gate (sharp, fixed-depth) AND the end-to-end converged
`rtol=1e-9` gate (the production claim).

**New test:**
`tests/sn/verification/mms/test_mms_ld_slab.py::test_ld_two_paths_scan_equals_dag_oracle`.
Markers `@pytest.mark.l1 @pytest.mark.verifies("ld-cartesian-1d")`.

**STRESSING CONFIG (the brief's (d) — Mode-9 degenerate-blindness defense):**
the two-paths gate MUST run on a config that BREAKS any FP-coincidence:
- **heterogeneous Σ_t(x)** (the STRESS ansatz's per-cell-varying σ_t, NOT a
  uniform medium — uniform σ_t can make the ×V and ÷V denoms coincide more
  closely);
- **≥2 groups with asymmetric downscatter** (1G is flux-shape-degenerate —
  a group-coupling drift between the two reduction trees is invisible in 1G;
  vv H1 + Mode-9);
- **non-flat per-ordinate field** (the STRESS ansatz's `μ·(b0+b1·sin)` term
  — a flat/isotropic field nulls the per-ordinate variation that the two
  schedules accumulate differently);
- a **moderate mesh** (`nx≈40`, the B5 stress-ansatz mesh) so the chain is
  deep enough that the cumprod-vs-walk reduction-tree difference is real.

Concretely, REUSE `SNSlabLDStressMMSCase` (the B5 case): drive
`solve_sn_fixed_source(case.materials, mesh, case.quadrature, Q,
cell_update=LD)` TWICE — once forcing `CumprodScan`, once forcing
`FullFieldWavefront` — and compare `result.scalar_flux.values`. **Forcing
the rep:** `solve_sn_fixed_source` selects via `default_for`, which after B
returns `CumprodScan` for LD. To force `FullFieldWavefront` for the oracle
leg, either (a) thread a `loss_representation=` override if the driver
exposes one, or (b) construct `SNMesh` + drive the rep directly via
`transport_sweep`/`SNSolver` with the explicit rep (the Increment-A pattern
— check whether `solve_sn_fixed_source` grew a `loss_representation` kwarg;
if not, the explicit-rep path is the test's vehicle and a `loss_representation`
kwarg is a candidate driver-ergonomics follow-up to FILE). **VERIFY the
forcing actually engaged** (assert `type(default_for(ld_mesh)) is CumprodScan`
for the default leg AND that the oracle leg genuinely ran FullFieldWavefront)
— a two-paths gate where both legs secretly ran the same rep is a silent
false green (the canary-that-cannot-trip, Mode-8-adjacent).

---

## GATE B7 — scan-matvec ≡ scan-sweep (leg 3 on the new path)

**Claim:** software invariant (foundation). On the NEW scan path, LD's
forward matvec (`CumprodScan.loss_action` → `_OneDimScanWalk.loss_action`)
≡ LD's SI sweep (`CumprodScan.sweep`). The Increment-A
`test_sn_1d_slab_ld_mms_krylov_matches_si` pins this for the DAG path; after
B's routing flip, that SAME test re-targets to the scan path automatically
(both `inner_solver="krylov"` and the default SI now route LD through
`CumprodScan`). **Re-confirm it stays green** (`rtol=1e-9, atol=1e-11`).

**⚠ HAZARD — does LD's scan `loss_action` exist?** `_OneDimScanWalk.loss_action`
(`loss_representation.py:742`) is the matvec twin of the sweep. It consumes
the SAME group-3 scan triple (`affine_scan_coefficients`,
`outgoing_face_from_average`). The carve MUST supply LD's group-3 such that
BOTH the sweep AND the matvec scan bodies work — `loss_action` reconstructs
the outflow via `outgoing_face_from_average` (`:2032` in the apply walk).
B7 (= the re-targeted Krylov≡SI test) is the catcher if LD's
`outgoing_face_from_average` is wrong in the matvec direction. If the carve
ships the scan SWEEP but not a working scan MATVEC, B7 fails (Krylov can't
converge or diverges from SI). **This is the matvec≡sweep leg — do NOT skip
it; a scan sweep without a verified scan matvec is half a carve** (the L14
leg-3 standoff demands both directions).

---

## L17 CONVENTION CROSSWALK — Increment B (the ×V scan convention)

The carve crosses: the per-cell ÷V kernel convention (group 2, LANDED) ↔
the ×V `CollisionCache` scan convention (group 3, NEW) ↔ the scan body's
`b` emission ↔ the `_OneDimScanWalk` sweep/matvec. (Template:
coding-elegance Pattern 7.)

| Value | group-2 kernel (÷V, LANDED) | group-3 scan (×V, NEW) | scan body (`_OneDimScanWalk`) | Normalizing side |
|-------|----------------------------|------------------------|-------------------------------|------------------|
| streaming coeff | `g = |μ|/Δ` (`s_axes/2`) | `2|μ|·A_down` folded into `inverse_denom` denom | the cache builds `inverse_denom`; the body reads `coll.inverse_denom[ords]` | **CELL OP** owns both conventions; B4 pins them equal at Q̂=0 |
| transmission `a` | implicit in `psi_out/psi_in` source-free | `a_attenuation` from `affine_scan_coefficients` (×V) — **DD's `2|μ|·A_total·inv−1` is NOT LD's**; LD's is the Schur outflow ratio | `coll.a_attenuation[ords]` → `ordinate_scan(a, b, ψ_in)` | **CELL OP** (`affine_scan_coefficients`) — THE sign-trap site |
| source emission `b` | n/a (kernel folds Q into rhs) | `b = scheme.scan_source_emission(QV, inv_denom)`; DD = `2·QV·inv_denom` | `:2556` slab / `:2718` curv — MOVE the inline `b` into the scheme method | **PRODUCER** (the scheme owns its emission; the body just calls it) — Pattern 7 |
| cell-average from faces | `psi_avg` solved directly | `cell_average_from_faces(in,out)`; DD=`½(in+out)`, LD=`(1−w)in+w·out` | `:2570`/`:2735` consume it | **CELL OP** (the weighted-mean weight `w` is σ_t-epoch per-cell) |
| outflow from average | `2ψ̄−in` (kernel) | `outgoing_face_from_average(ψ̄,in)`; DD=`2ψ̄−in`, LD via slope row | the matvec apply walk `:2032` consumes it | **CELL OP** — B7 (matvec) catches a wrong LD form here |
| `/W` source norm | already `/W` (producer) | SAME — `QV = Q·V`, Q already `/W` | the body multiplies by V only (`:2447`) | **PRODUCER** (unchanged from DD; LD MUST NOT re-divide) |

**The Bridge rows (where Pattern 7 demands action):**
- **`b` emission**: bridge is the SCHEME method (move the inline `b` to
  `scan_source_emission`). DD's override reproduces `2.0·QV·inv_denom`
  byte-exact (B1 catcher). LD's encodes the Q̂=0 Schur emission.
- **transmission `a`**: bridge is `affine_scan_coefficients`. DD's is
  `2|μ|·A_total·inv−1`; LD's is its OWN Schur ratio (NOT a copy of DD's
  with different numbers — a genuinely different formula). **The single
  most likely sign-trap re-introduction site; B4 is its catcher.**
- **`/W`**: bridge already at the producer (Pattern 7 satisfied). LD MUST
  NOT re-apply `/W` in `scan_source_emission` (DD doesn't — the `2·QV·inv`
  has no `/W`; the `/W` was applied upstream by `AngularSourceSink`).

---

## Failure-mode coverage — Increment B (vv §6 + numerical-bug-signatures)

| Failure mode | Increment-B gate that catches it |
|---|---|
| #1 sign flip (LD scan transmission `a` / slope row in ×V) | **B4** group2≡group3 at Q̂=0 (per-cell magnitudes diverge) + **B5** stress MMS diverges under refinement |
| #2 variable swap (`a` vs `inv_denom`, in vs out face) | **B4** + **B6** two-paths (scan ≠ DAG oracle) |
| #3 missing/extra factor (the moved `b` drops the 2, or the LD weight `w`) | **B1** DD bit-id breaks (factor-2 changed) + **B4** (LD `psi_avg` wrong) + **B5** value band |
| #4 wrong recursion (cumprod chain index drift in LD scan) | **B6** two-paths (scan ≠ DAG) + **B7** matvec≡sweep |
| #5 index error (scan chain reorder for LD) | **B5** stress MMS (non-uniform field) + **B6** |
| #6 convention drift (×V vs ÷V denom grouping; `/W` re-applied; `Lψ` vs `(L+C)ψ` in scan matvec) | L17 crosswalk Bridge rows; **B4** (denom) + **B7** (matvec `loss_action`) |
| Mode 7 MMS simplification bias | **B5** stress ansatz override (a0>0, mixed k=1+k=3, μ-nontrivial, 2G het) — OWED since Inc. A, landed here |
| Mode 8 `-O` strips bare assert | ALL Increment-B gates use `np.testing.assert_*` / `pytest.fail` (verified the Increment-A LD tests fire under `-O` — 21 passed live) |
| **Mode 9 splitting/acceleration in degenerate regime** | **B6 stressing config** — a two-paths gate on a flat/homogeneous/1G problem is Mode-9-blind (the two reduction trees can be accidentally FP-coincident on a smooth uniform field). B6 runs on heterogeneous ≥2G non-flat (the STRESS ansatz) → breaks the coincidence. **THE Increment-B-specific Mode-9 hazard** (Increment A had none — it was a single discretization on a single path; B introduces a SECOND path → the FP-invariance-across-paths claim is exactly the kind Mode-9 warns about, even though it's a schedule change not a splitting). |

**NEW failure-mode-table consideration (self-improvement trigger):** Mode 9
historically targets *splittings/accelerators* (G-S, σ_r-fold, DSA). The
Increment-B two-paths agreement is a *schedule* change (scan vs DAG), not a
splitting — but the SAME degenerate-regime-blindness applies: the FP-
invariance-across-schedules claim can pass on a degenerate config and hide a
real per-ordinate/group reduction-tree bug. **This is a NEW sub-case of Mode
9 (schedule-equivalence, not splitting-equivalence).** Per the
self-improvement directive, this generalization of Mode 9 ("an
equivalence-across-implementations FP claim — splitting OR schedule OR
representation — must run on a degeneracy-breaking config") should be noted
when the plan is delivered. It does NOT need a new ERR (no caught bug yet),
but the Mode-9 skill row's scope ("a splitting") is narrower than the
phenomenon — flag to sharpen "splitting/acceleration" → "splitting /
acceleration / schedule / representation equivalence" in the vv-principles
Mode-9 row. (The #206 Phase-C scan-vs-DAG nULP gates are the same class and
already in practice — this just names the principle.)

---

## Summary for the dispatcher (Increment B)

- **Headline gate B6 (two-paths):** `solve_sn_fixed_source(LD)`-via-
  `CumprodScan` ≡ -via-`FullFieldWavefront` on the STRESS ansatz
  (heterogeneous, 2G asymmetric downscatter, μ-non-trivial, nx≈40).
  Tolerance = `rtol=1e-9` (iterative→`SAFETY×inner_tol`) end-to-end, PLUS
  a single-sweep `nulp=nx` gate (fixed-depth). NOT `array_equal` (×V-vs-÷V
  denom grouping + cumprod-vs-walk reduction trees). PAIR with B5 (both
  paths independently match the MMS ref — the structural ground). VERIFY
  each leg actually ran its rep (no silent same-rep false green).
- **Routing test INVERTS, do NOT add-alongside:** rename
  `test_ld_slab_mesh_routes_to_full_field_wavefront` →
  `…_routes_to_cumprod_scan`, flip the expected type. The mutual-exclusion
  flips too: `CumprodScan.supports(ld_mesh)` becomes `True` after B.
- **DD strict bit-identity (B1) is the negative control:** `pytest
  tests/sn/sweep/core tests/sn/solve -W
  "error::tests.sn.regression._regression_assert.DriftWarning"` MUST stay
  505/1/4 (confirmed live 2026-06-14). Snapshots = `test_sweep_regression.py`
  (DD slab/sphere/cyl scan). The moved `b` emission MUST reproduce DD's
  `2.0·QV·inv_denom` op-order byte-exact (single-source scheme method; DD =
  the special case).
- **The sign-trap gate (B4):** group2(÷V kernel) ≡ group3(×V scan coeffs)
  at Q̂=0 — the third leg of the Increment-A group1≡group2 cross-check, in
  the NEW ×V scan convention where the LM-1989 slope-row sign re-appears.
  `rtol=1e-12`, 1G+2G. The single most likely re-introduction site is LD's
  `affine_scan_coefficients` transmission `a` (NOT a copy of DD's
  `2|μ|·A_total·inv−1` — LD's own Schur ratio).
- **The owed Mode-7 override (B5):** land `SNSlabLDStressMMSCase` (a0>0,
  k=1+k=3 mixed scales, μ-nontrivial, 2G het) NOW — the two-paths gate
  (B6) needs a degeneracy-breaking config or it's Mode-9-blind. Source via
  SymPy (Branch 1), foundation derive-test, flat-source posture (Q̂=0,
  O(h²) from the solved slope unknown).
- **The matvec leg (B7):** the re-targeted `test_sn_1d_slab_ld_mms_krylov_matches_si`
  must stay green on the scan path — LD's `outgoing_face_from_average` must
  work in BOTH the scan sweep AND the scan matvec (`_OneDimScanWalk.loss_action`).
  A scan sweep without a verified scan matvec is half a carve.
- **Flat-source caveat UNCHANGED:** B is a schedule bonus, NOT the slope
  source. `test_ld_thick_diffusive_limit_xfail` stays xfail. Do NOT claim B
  fixes the diffusion limit.

---

# ════════════════════════════════════════════════════════════════════
# ORIGINAL PRE-IMPL PLAN (Increment-A era — SUPERSEDED, kept as history)
# The OneDimPerCellWalk seam below was NOT how LD landed (reverted per user;
# LD rides FullFieldWavefront's group-2 kernel). The MMS ansatz override
# (gate 3) and the DD-strict-gate invocation (gate 5) carry forward into
# Increment B above. Read the Increment-B section for the CURRENT contract.
# ════════════════════════════════════════════════════════════════════

Pre-implementation verification plan. The SUT (`LinearDiscontinuous`
cell-update + `OneDimPerCellWalk` loss-representation) does NOT exist yet;
every NEW gate ships as a collectible `@pytest.mark.skip(reason="SUT not
yet implemented")` stub. This is an **L17 carve** (per-cell `update`/`residual`
↔ the `_OneDimScanWalk` frame ↔ the matvec apply path cross subsystem
boundaries) so the proactive-dispatch protocol fired correctly.

**Why:** add a non-Diamond-Difference (`is_affine_scannable=False`) spatial
closure and the routing that lets a non-affine scheme reach the production
sweep. DD continues to select the fast `CumprodScan`/`ScanMarch`; LD is the
ONLY occupant of the new DAG-free per-cell path.

**How to apply:** when the method-implementer builds LD, the 5 gates below
are the acceptance contract. Read this BEFORE writing any LD code.

Related skill notes: [[w2_pole_cell_characterization_gate]] (gate 4 is the
#233 gate I designed — already LD-ready by lower-bound construction);
[[w1_sphere_clamp_removal_verification]] (the sphere pole lift is gated by
the angular floor #235 — out of scope, do NOT make gate 4 hinge on
sphere-pole order≈2).

---

## Live evidence (lesson L12 — claims need verbatim stdout)

Env: `CLAUDE_ENVIRONMENT=` (Host) → `.venv/bin/python -O -m pytest`.

DD round-trip (the gate-1 model) is green:
```
.venv/bin/python -O -m pytest \
  tests/sn/sweep/core/test_diamond.py::TestResidual::test_residual_zero_at_solved_cell_avg \
  tests/sn/sweep/core/test_diamond.py::TestResidual::test_residual_zero_multi_group_heterogeneous -q
→ 20 passed, 1 warning in 0.37s
```
(the warning is the `-O`-strips-bare-assert PytestConfigWarning — vv Mode 8;
the `TestResidual` round-trip uses `np.testing.assert_allclose`, which FIRES
under `-O`. Confirmed.)

DD-no-regression strict gate (gate 5) — the exact invocation, currently GREEN:
```
.venv/bin/python -O -m pytest tests/sn/sweep/core tests/sn/solve -q \
  -W "error::tests.sn.regression._regression_assert.DriftWarning"
→ 505 passed, 1 skipped, 4 xfailed, 2 warnings in 83.15s
```

Why NOT `tests/sn/regression/`: its `conftest.py:pytest_configure` calls
`warnings.simplefilter("always", DriftWarning)` AND adds an
`always::…DriftWarning` ini line → the `always` filter shadows `-W error`,
making the escalation INERT there. The strict bit-identity gate MUST run on
`sweep/core` + `solve`, which carry no such override (verified: their
conftests do not mention DriftWarning). This matches FINDING-1 in
[[issue_206_phase_c_verification]] (`-W error::DriftWarning` escalates ONLY
in sweep/core + solve).

Existing slab MMS + 2G heterogeneous MMS green:
```
.venv/bin/python -O -m pytest \
  tests/sn/verification/mms/test_mms.py::test_sn_1d_slab_mms_converges_second_order \
  tests/sn/verification/mms/test_mms_heterogeneous.py -q
→ 5 passed, 2 warnings in 2.35s
```

---

## SUT seams confirmed by reading production (do NOT invent)

- `CellUpdateBase` ABC: `orpheus/sn/spatial/cell_update.py:509`. Abstract
  `update` (`:565`, → `CellResult`) + `residual` (`:575`, → `(ng,)` array).
  `is_affine_scannable: ClassVar[bool] = False` (`:553`) is the OPT-OUT
  default — LD inherits `False`, DD overrides `True` (`diamond.py:127`).
  Round-trip contract note: `cell_update.py:443-458` ("residual at the
  cell-average returned by update yields zero to FP rounding").
- DD: `orpheus/sn/spatial/diamond.py`. `update` (`:138`), `residual`
  (`:189`, delegates to `cell_balance_for_streaming` at n_mask=1). The shape
  LD mirrors. DD also owns 3 capability groups (per-cell ref pair; batched
  `cell_kernel_batch`/`residual_kernel_batch`; scan-family triple) — LD
  supplies group 1 (mandatory) and, being non-affine, NOT group 3.
- Selection: `orpheus/sn/loss_representation.py`. `CumprodScan.supports`
  (`:701`) = `mesh.is_1d AND mesh.cell_update.is_affine_scannable`.
  `ScanMarch.supports` (`:1217`) = `(is_1d OR (is_cartesian AND ndim==2))
  AND is_affine_scannable`. `LOSS_REPRESENTATIONS` tuple (`:1460`) =
  `(CumprodScan, ScanMarch, MovingFrontierWindow, FullFieldWavefront)`;
  `default_for` (`:1468`) returns the FIRST whose `supports().ok`.
  - **THE GAP THE SEAM FILLS:** a 1-D `is_affine_scannable=False` mesh is
    admitted by NOTHING today — CumprodScan/ScanMarch reject (not
    scannable), the two `_DAGWavefront` reps require Cartesian d≥2.
    `default_for` would raise `IncompatibleRepresentation`. The new
    `OneDimPerCellWalk` (`supports = is_1d AND NOT is_affine_scannable`)
    fills exactly this hole. It MUST be inserted into
    `LOSS_REPRESENTATIONS`; placement vs CumprodScan/ScanMarch is
    irrelevant for correctness (mutually exclusive on the scannable trait),
    but place it adjacent to CumprodScan for readability (both are the
    1-D family).
- `SNMesh(..., cell_update=…)`: `geometry.py:203` constructor kwarg →
  `self.cell_update` defaults to `DiamondDifference()` (`:271`).
- **DRIVER GAP (the method-implementer MUST close):**
  `solve_sn_fixed_source` (`solver.py:1843`) takes a `Mesh1D`, NOT an
  `SNMesh`, and builds `SNMesh(geometry, quadrature, materials)` internally
  at `solver.py:1971` with NO `cell_update` passed → it ALWAYS uses DD.
  To run an LD MMS through the high-level driver, thread a
  `cell_update: CellUpdateBase | None = None` kwarg through
  `solve_sn_fixed_source` → `SNMesh(... , cell_update=cell_update)`.
  Until that lands, the gate-3 skeleton constructs the `SNMesh` directly
  (`cell_update=LinearDiscontinuous()`) and drives `SNSolver(sn_mesh, …)`.

---

## L14 four-leg standoff — how the 5 gates instantiate it for LD

A scheme is correct only when (1) it matches a structurally-independent
reference, (2) the matvec matches a reference, (3) sweep≡matvec, (4) all
hold under refinement.

| Leg | Claim layer | Gate |
|-----|-------------|------|
| (3) sweep≡matvec | software invariant (foundation) | **Gate 1** per-cell round-trip `residual(update(q).ψ̄)=0` |
| (3) sweep frame correct | software invariant (foundation) | **Gate 2** seam-control: DD-via-`OneDimPerCellWalk` ≡ DD-via-`CumprodScan` |
| (1) structurally-indep reference | convergence-order (math, L1) | **Gate 3** slab-Cartesian LD MMS O(h²) (MMS = math/flux-shape pillar; NOT eigenvalue) |
| (1) curvilinear closure char. | characterization (L1, no `verifies`) | **Gate 4** #233 pole-cell (LD-ready by lower-bound; sphere lift deferred #235) |
| no-regression (the negative control) | bit-identity inheritance | **Gate 5** DD bit-identical everywhere |
| (4) refinement | built INTO gate 3 (the ladder) | the convergence ladder |

Note the claim layers (vv-principles hierarchical taxonomy): gate 3 is a
**convergence-order/flux-shape** claim — MMS does NOT and CANNOT prove an
eigenvalue. There is NO eigenvalue gate in this plan because LD is a
spatial-closure change verified at the convergence-order + flux-shape
layers; the existing homogeneous `k_inf` eigenvalue gates (DD path) are
untouched by LD and stay green via gate 5. If a future task claims an LD
eigenvalue improvement, that needs a closed-form reference (homogeneous
`k_inf = λ_max(A⁻¹F)`), NOT this MMS.

---

## GATE 1 — Per-cell round-trip (lockstep), isolated

**Claim:** software invariant (foundation). `update`/`residual` describe the
SAME per-cell linear system: `residual(update(...).cell_average_flux, ...)`
vanishes to FP noise.

**Pinned by (existing, for DD):**
`tests/sn/sweep/core/test_diamond.py::TestResidual::test_residual_zero_at_solved_cell_avg`
(parametrized over 5 geometries) +
`::test_residual_zero_multi_group_heterogeneous` (×{1,2,4}G, heterogeneous XS).
This is exactly the round-trip the brief asks LD to mirror.

**New test (catches LD):**
`tests/sn/spatial/test_linear_discontinuous.py::TestLDRoundTrip` — mirror DD's
`TestResidual` structure for the LD occupant alone. LD is `is_linear` (per
Lewis & Miller §5.3 the LD closure IS linear in source + upstream), so keep
the linearity-in-cell_avg and affine-in-source companions from DD's
`TestResidual` (modes 4 + 5 catchers).

**File path:** `tests/sn/spatial/test_linear_discontinuous.py` (the new
home; `tests/sn/spatial/` exists). Markers: `@pytest.mark.foundation` (no
`verifies` — it pins a software contract, not an equation; mirrors DD's
`TestResidual`).

**Crucial design notes:**
- L2: parametrize `n_groups ∈ {1, 2, 4}` with HETEROGENEOUS per-group
  total_xs (1G is the degenerate control, must hold; 2G/4G are the real
  gate). DD's `test_residual_zero_multi_group_heterogeneous` is the template.
- LD couples TWO face moments (cell-average + slope/second moment). The
  `CellResult` shape may need a slope output; if LD returns extra moments,
  the round-trip must feed ALL of them back to `residual` (a partial
  feed-back would pass spuriously). The round-trip is `residual(update(q)) =
  0` on the FULL solved state, not just the average.
- atol band: a few division ULP — `atol=1e-13` like DD (one division ULP).
  If LD's per-cell solve is a 2×2 (two moments), bump to `atol=1e-12`
  (reduction depth ~ a 2×2 solve) and document the bound.
- Term-activation declaration (vv Mode 7): each test states which terms
  its inputs ACTIVATE. The round-trip activates streaming + collision +
  (slab: no angular redistribution). NO term is nulled.

---

## GATE 2 — Seam-opening control (principled-equiv)

**Claim:** software invariant (foundation). Route DD itself through the NEW
`OneDimPerCellWalk` per-cell path and assert it reproduces DD's fast-scan
answer (`CumprodScan`). Proves the per-cell frame is correct INDEPENDENTLY
of LD (DD is the structurally-independent reference here — its fast-scan
answer is already bit-identity-pinned by gate 5's regression snapshots).

**Pinned by (existing):** the SELECTION half is
`tests/sn/sweep/core/test_unified_sweep_dispatch.py::TestDispatchSelectsStrategy`
(slab/sphere/cyl → `CumprodScan`). NO existing test routes DD through a
per-cell walk and compares to the scan — gate 2 is genuinely NEW.

**New test (catches the seam):**
`tests/sn/sweep/core/test_unified_sweep_dispatch.py::TestOneDimPerCellSeam`
(co-locate with the dispatch contract it extends). Two sub-tests:
1. **Selection pin** — a 1-D mesh with `cell_update=LinearDiscontinuous()`
   makes `default_for(mesh)` return `OneDimPerCellWalk` (and a 1-D mesh with
   DD still returns `CumprodScan`). `pytest.fail` on mismatch (-O-safe).
   ALSO assert `CumprodScan.supports(ld_mesh).ok is False` and
   `OneDimPerCellWalk.supports(dd_mesh).ok is False` — the mutual exclusion
   (vv anti-pattern #11: a contract gate needs BOTH a must-select and a
   must-NOT-select case).
2. **Equivalence (the control)** — sweep a fixed `(Q, sig_t, boundary)`
   through DD twice: once via `CumprodScan(dd_mesh).sweep(...)` (the fast
   scan) and once via `OneDimPerCellWalk(dd_mesh).sweep(...)` (the per-cell
   loop, DD as the occupant). Assert agreement.

**TOLERANCE DECISION (justified):** `assert_array_almost_equal_nulp(nulp=nx)`
— **principled-equivalence, NOT bit-identity.** Reason (vv-principles
§Bit-identity vs principled-equivalence): the two paths compute the SAME
per-cell DD algebra in real arithmetic but accumulate over a DIFFERENT
IEEE-754 reduction tree. `CumprodScan` evaluates the affine recurrence as a
closed-form parallel-prefix `cumprod`/`cumsum` (Blelloch §1.5) — the
streaming attenuation is multiplied in a prefix-product chain. The
`OneDimPerCellWalk` evaluates the SAME recurrence as a sequential cell-by-cell
`update` loop. Addition/multiplication are non-associative, so the chain
differs at ULP. The drift is FP-non-associativity, dimensionally bounded by
`(reduction depth) × ULP = nx × eps` (the chain length is the cell count) —
criterion 3. The reference (DD's scan answer) is itself bit-identity-pinned
to a structurally-independent ground (the regression snapshots + the
homogeneous `k_inf` analytic limit) — criterion 2. Each intermediate (the
per-cell `cell_average_flux`) is a named, inspectable reactor-physics
quantity — criterion 1. All three hold ⟹ relax to `nulp=nx`. **Do NOT use
`np.array_equal`** here (it would be the wrong gate — it demands a reduction
tree the two paths legitimately do not share). This is the EXACT precedent
of the `kind=direct reduction_depth=nx` A-NEW matvec leg in
[[issue_206_phase_c_verification]] (density-vs-scan denom grouping → nULP
`≤nx`).
- Companion: ALSO assert the per-cell-walk path agrees with `CumprodScan`
  on a 2G heterogeneous mesh (L2 — 1G is flux-shape-degenerate; the seam
  must be exercised on a multi-group field so a group-coupling drift in the
  frame would show). `markers: @pytest.mark.foundation`.

---

## GATE 3 — MMS O(h²) convergence on Cartesian SLAB for LD

**Claim:** convergence-order (math, L1) + flux-shape. The interference-free
spatial oracle: slab Cartesian has NO space-angle coupling (no angular
redistribution term — `upstream_state.angular_upstream is None` for slab),
so it isolates the spatial order. MMS = the math/flux-shape pillar (NOT
eigenvalue).

**Pinned by (existing, for DD):**
`tests/sn/verification/mms/test_mms.py::test_sn_1d_slab_mms_converges_second_order`
(1G) + `test_mms_heterogeneous.py::test_sn_heterogeneous_mms_converges_second_order`
(2G heterogeneous continuous-Σ). The LD gate is a SIBLING.

**New test (catches LD spatial order):**
`tests/sn/verification/mms/test_mms_ld_slab.py::test_sn_1d_slab_ld_mms_converges_second_order`.
Markers: `@pytest.mark.l1 @pytest.mark.slow @pytest.mark.verifies(...)` with
LD-specific equation labels (the LD theory page MUST add
`:label: ld-cartesian-1d` / `ld-slab` — flag to the archivist; the
`verifies` here is the gate that pins them). Ladder `n_cells = [10, 20, 40,
80]` fixed quadrature; assert `np.all(orders > 1.9)` + magnitude band
`1e-9 < errors[-1] < 1e-3`.

### THE MMS ANSATZ (the simplification-bias override — load-bearing)

**Every existing slab MMS uses `A(x) = sin(πx/L)`** (verified by reading
`orpheus/derivations/continuous/mms/sn.py:79-503`). That ansatz is the
canonical Mode-7 simplification bias:
1. it VANISHES at both boundaries (`A(0)=A(L)=0`) → vacuum BC is automatic
   but the BC handling is **untested** (MMS rule: "trial MUST be non-trivial
   at boundaries to verify BC handling");
2. it is **isotropic-in-μ** (`ψ_{n,g} = c_g·A(x)/W`) → the only per-ordinate
   variation is the smooth streaming term `μ_n·c_g·A'(x)`;
3. it is "the simplest trig that satisfies the BCs".

For an LD verification this is INSUFFICIENT — LD's whole selling point is
O(h²) + diffusion-limit consistency on a stressed spatial profile. **MUST
override the bias per vv-principles §MMS operational rules.** The
test-architect-chosen ansatz:

```
ψ_{n,g}(x) = [ a0_g + a1_g·sin(k₁ π x/L) + a2_g·cos(k₂ π x/L)
               + μ_n·( b0_g + b1_g·sin(k₃ π x/L) ) ] / W
```

with **a0_g > 0 load-bearing** (the constant term makes ψ NON-vanishing at
both boundaries → vacuum BC is now a real inflow=0 condition that the LD
boundary closure must satisfy with a non-trivial interior field; this is the
non-vacuum-flavored stress that [[phase4_46_nonvacuum_mms_ansatz]] showed is
the untested corner). Concretely:
- `k₁ = 1`, `k₂ = 3`, `k₃ = 2` → **mixed spatial scales** (a low + a high
  harmonic) so the O(h²) error is dominated by the high-frequency component
  at coarse mesh — this STRESSES the closure (the simplification bias would
  pick a single `k=1`).
- **angularly non-trivial** via the `μ_n·(b0 + b1·sin)` term — even though
  slab has no angular redistribution, the per-ordinate streaming `μ·∂ψ/∂x`
  now carries TWO spatial scales per ordinate, and the cell-average vs
  slope moments of LD see a genuinely μ-dependent field (NOT the
  isotropic-flat field the `sin(πx/L)` ansatz produces). L27: the
  manufactured source is per-ordinate by construction here.
- **≥2 groups, HETEROGENEOUS Σ_t(x)** (reuse the smooth continuous-Σ
  profile from `_default_hetero_xs_functions`, `sn.py:506`, which guarantees
  Σ_a>0). Non-trivial group ratio `c`-vector (e.g. via per-group `a`/`b`
  coefficients) so downscatter coupling is exercised (mode 6 / convention
  drift catcher) — exactly as the existing 2G case argues at `sn.py:290`.

**Manufactured source** (derive symbolically — algebra-of-record Branch 1;
do NOT hand-transcribe): substitute the ansatz into
`μ_n ∂_x ψ_{n,g} + Σ_{t,g}(x) ψ_{n,g} = (1/W)(Σ_g' Σ_{s,g'→g}(x) φ_g' +
Q^ext_{n,g})` and solve for `Q^ext`. φ_g = ∫ψ dμ — the `μ_n·(…)` term
integrates to zero over symmetric GL ordinates, so `φ_g(x) = a0_g + a1_g·sin
+ a2_g·cos`; the streaming derivative carries the full ansatz. This is the
load-bearing equation — generate it with SymPy (a new
`SNSlabLDStressMMSCase` in `orpheus/derivations/continuous/mms/sn.py`), pin
the symbolic source with a `@pytest.mark.foundation` derive-test, then
consume it in the L1 convergence test.

**Why O(h²) and not higher:** standard LD is O(h²) in the cell-average flux
(Lewis & Miller §5.3); the gate asserts `orders > 1.9`. Do NOT assert >2.9 —
LD is second-order, not third. **NEVER** treat the rate alone as correctness
(vv §5): the gate ALSO pins the converged VALUE via the magnitude band
`1e-9 < errors[-1] < 1e-3` against the manufactured `phi_exact` (the
structurally-independent reference IS the imposed solution — that is what
makes MMS a flux-shape pillar). The high-frequency `k=3` component means the
coarse-mesh error is genuinely large (so the ladder has real dynamic range,
not floor-dominated) — confirm `errors[0]` is well above the magnitude band
floor at `nx=10`.

**Pole/curvilinear companion NOT required for gate 3** (gate 3 is
slab-only, the interference-free oracle by design). The curvilinear LD lift
is gate 4 + the deferred #235 angular floor.

---

## GATE 4 — #233 pole-cell characterization (LD-ready, sphere lift deferred)

**Claim:** characterization (L1, NO `verifies` — pins a limitation, not a
correctness claim). The existing gate
`tests/sn/verification/mms/test_curvilinear_pole_cell_characterization.py`
pins the DD pole-cell as LOWER-bounded first order with NO upper bound. A
higher-order scheme (LD) should make the pole order rise — and the gate
STAYS GREEN by construction (`order 2.0 > 0.8 lower bound`), DOCUMENTING the
improvement.

**This gate already exists and is already LD-ready — I designed it that way**
([[w2_pole_cell_characterization_gate]]). The four tests:
- `test_sphere_global_L2_second_order_dual_reference` (GUARANTEE, `verifies`)
- `test_cylinder_global_L2_second_order` (GUARANTEE, `verifies`)
- `test_sphere_pole_cell_first_order_and_Linf_dominant` (CHARACTERIZATION,
  no `verifies`, `catches("ERR-059")`, pole order `> 0.8` lower bound only)
- `test_cylinder_pole_first_order_vs_volume_average_masked_by_midpoint`
  (CHARACTERIZATION, lower bound only)

**What the LD task DOES here (NO new test by default):** run this gate with
an LD curvilinear mesh and DOCUMENT the measured pole order in the LD
closeout memo. **CRITICAL — the plan's success MUST NOT hinge on the sphere
pole asserting order≈2 now:** on the sphere the pole lift is GATED by the
ANGULAR floor (issue #235 — the M-M half-angle interpolation, NOT in scope).
LD fixes the SPATIAL central-cell closure but the sphere pole error is
co-limited by the angular floor, so the sphere pole order may stay ~O(h)
even with LD until #235 lands. Per [[w1_sphere_clamp_removal_verification]]
the sphere aniso/pole floor is #229-interpolation-bound and needs S32 to
clean — same family. So:
- The cylinder pole (no angular floor co-limit on the radial closure) is
  the cleaner LD discriminator; if LD lifts the cylinder pole-vs-volume-avg
  from O(h) toward O(h²), the existing `> 0.8` lower bound stays green and
  the closeout records `ord_va` rising.
- **OPTIONAL new test (stage as `xfail(strict=False, reason="#235 angular
  floor co-limits the sphere pole until S32/M-M-interp fix")`):** an
  LD-specific `test_sphere_pole_LD_lifts_toward_second_order` that ASPIRES
  to `pole_orders > 1.5`. Staged xfail so it flips to xpass WHEN #235 lands
  WITH LD — a forward tripwire, not a blocking gate. Do NOT make it strict;
  do NOT let gate-4 pass/fail depend on it.

---

## GATE 5 — DD bit-identical everywhere (no regression)

**Claim:** bit-identity inheritance — adding LD + `OneDimPerCellWalk` must
NOT perturb DD's reduction tree. The negative control for the whole carve.

**The strict live gate (confirmed green above):**
```
.venv/bin/python -O -m pytest tests/sn/sweep/core tests/sn/solve \
  -W "error::tests.sn.regression._regression_assert.DriftWarning"
```
505 passed / 1 skipped / 4 xfailed at HEAD. After the LD landing this MUST
still be 505/… green (the 4 xfailed are pre-existing standing reds, NOT
introduced by LD — re-confirm their ids are unchanged). NO new test file —
gate 5 IS this invocation on the EXISTING snapshots, run in CI as the carve's
acceptance.

**Plus a selection-no-regression assertion (NEW, cheap):** in
`TestOneDimPerCellSeam` (gate 2 file) add `test_dd_mesh_still_selects_scan`:
a default (DD) 1-D mesh after the `LOSS_REPRESENTATIONS` edit STILL returns
`CumprodScan`, and a 2-D DD mesh still returns `ScanMarch`. This pins that
inserting `OneDimPerCellWalk` into the registry tuple did NOT shadow the
scan reps for DD meshes (a registry-ordering regression the bit-identity
snapshots would ALSO catch, but this localizes it to the selection layer).
`pytest.fail` on mismatch (-O-safe).

**Why these two dirs and not `tests/sn/regression/`:** see Live-evidence
§. The `regression/` conftest forces `always::DriftWarning` → `-W error`
inert there; `sweep/core` + `solve` carry no override (FINDING-1).

---

## L17 CONVENTION CROSSWALK (required deliverable)

The carve crosses: per-cell `update` SOLVE ↔ per-cell `residual` APPLY ↔ the
`_OneDimScanWalk` / new `OneDimPerCellWalk` frame's walk inputs ↔ the matvec
apply path. For each value crossing a boundary: meaning / units / layout on
each side, and which side normalizes. (Template: coding-elegance Pattern 7.)

| Value | `update` (SOLVE) side | `residual` (APPLY) side | `OneDimPerCellWalk.sweep/loss_action` frame side | Normalizing side |
|-------|----------------------|--------------------------|---------------------------------------------------|------------------|
| `source` (Q) | `(ng,)` per-group volumetric source, **already weight-normalised** `Q·V·(1/Σw)` (cell_update.py:402-408) | `(ng,)` SAME already-`/W` convention; residual is affine in it (`δsource → −δresidual`, cell_update.py:485) | the frame builds `Q_cell = QV·(1/Σw)` from the `AngularSourceSink` (producer already applied `/W` per Pattern 7, lesson L18; the sweep does NOT re-apply) | **PRODUCER** (`AngularSourceSink.from_isotropic` applies `1/W` at the boundary; `transport_sweep` docstring `loss_representation.py:1535-1545`) |
| `spatial_upstream` (ψ_in) | `(ng,)` face flux flowing INTO the cell, **sweep-direction-resolved** by the orchestrator (cell_update.py:411) | `(ng,)` SAME — the residual reads `upstream_state.spatial_upstream[:,None]` (diamond.py:233) | the frame resolves which face is upstream from the chain order + sign(μ); the per-cell op sees only "what flows in" (NO sign-of-μ branch inside the op) | **FRAME** (`OneDimPerCellWalk` resolves direction; the cell op is direction-blind — geometry-as-data, cell_update.py:111-131) |
| `face_area_downstream` | `float`; slab=1.0 (neutral), curvilinear=outer(μ≥0)/inner(μ<0), degenerate=0.0 (cell_update.py:215-235) | SAME float; residual builds `A_down_arr=[face_area_downstream]` (diamond.py:228) | the frame picks outer/inner from sign(μ) BEFORE issuing the `CellVisit` | **FRAME** (sweep-direction resolution is the orchestrator's job) |
| `cell_average_flux` (ψ̄) | OUTPUT of `update` → `CellResult.cell_average_flux` `(ng,)` | INPUT probe to `residual` `(ng,)`; residual = `denom·ψ̄ − (source+numer_upstream)` (diamond.py:279) | the frame stores ψ̄ into the per-ordinate field and feeds the outgoing face downstream | **lockstep contract** — the round-trip (gate 1) IS this crosswalk's correctness gate: `residual(update(q).ψ̄)=0` |
| `outgoing_spatial_flux` (ψ_out) | OUTPUT `CellResult.outgoing_spatial_flux` `(ng,)` or `None` (degenerate); DD: `2ψ̄−ψ_in` (diamond.py:167) | reconstructed from the probe: `2ψ̄−ψ_in` (diamond.py:374, the matvec twin) | the frame writes ψ_out as the NEXT cell's `spatial_upstream` | **CELL OP** owns the closure formula; frame just chains it. LD couples TWO moments → `outgoing_spatial_flux` may be richer (LD slope); the frame MUST chain ALL emitted moments |
| `outgoing_angular_state` (ψ_{n+½}) | OUTPUT `(ng,)` or `None` (slab); M-M closure (diamond.py:177) | M-M coefficients read inline (diamond.py:250-265) | curvilinear only; slab `None` (no angular redistribution) — **gate 3 is slab so this is `None` throughout gate 3** | **CELL OP** (M-M closure) — but slab nulls it; the LD slab MMS does NOT exercise this column (deliberately — gate 4 + #235 do) |
| matvec result `(L+C)ψ̄` | n/a | `residual` returns `L_cell·ψ̄ − q` per cell `(ng,)`; the operator subtracts `C=σ_t⊙` in `apply` (loss_representation.py:730-742) | `OneDimPerCellWalk.loss_action` chains `residual` per cell → full `(L+C)ψ` (the apply twin of `sweep` — L21 "matvec ≡ sweep") | **OPERATOR** subtracts `C` at the `apply` layer; the loss_action returns the FULL `(L+C)ψ` (NOT bare `Lψ`) — same convention as `CumprodScan.loss_action` (loss_representation.py:730) |

**The Bridge rows (where Pattern 7 demands action):**
- `source` `/W`: bridge is at the PRODUCER (Pattern 7 satisfied — the sweep
  never re-applies `/W`; LD's `update` must consume the already-normalised
  source, same as DD). **LD MUST NOT re-divide by `Σw`.**
- direction resolution: bridge is the FRAME (`OneDimPerCellWalk`), NOT the
  cell op. **LD's `update`/`residual` MUST be direction-blind** (read
  `face_area_downstream` + `spatial_upstream`, never sign(μ)). This is the
  invariant gate 2's equivalence test protects: if LD's frame resolved
  direction differently from DD's `CumprodScan`/`_OneDimScanWalk`, the
  DD-through-the-new-frame control would drift beyond `nulp=nx`.
- `(L+C)` vs bare `L`: bridge is the OPERATOR `apply`/`apply_transpose`
  (subtracts `C`). **LD's `loss_action` MUST return `(L+C)ψ`** (full
  within-group loss), matching `CumprodScan.loss_action` — NOT bare `Lψ`.
  A bare-`Lψ` return would double-count or drop the collision diagonal; the
  round-trip (gate 1) and the matvec leg (the `loss_action` consumed by the
  Krylov path) are the catchers.

---

## Failure-mode coverage (vv-principles §6 + numerical-bug-signatures)

| Failure mode | LD gate that catches it |
|---|---|
| #1 sign flip (LD slope/moment sign) | gate 3 het convergence DIVERGES under refinement |
| #2 variable swap (cell-avg vs slope moment) | gate 1 round-trip fails (the two moments don't round-trip) + gate 3 wrong order |
| #3 missing factor (LD weight/volume) | gate 3 magnitude band fails (wrong converged VALUE) + fixed-source `Q/Σt` flat check (add as a gate-1 companion) |
| #4 wrong recursion (per-cell loop index drift in the frame) | gate 2 seam control — DD-via-frame ≠ DD-via-scan beyond nulp |
| #5 index error (frame chains wrong cell's ψ_out) | gate 2 (non-uniform field) + gate 3 |
| #6 convention drift (`/W` re-applied, or `Lψ` vs `(L+C)ψ`) | L17 crosswalk Bridge rows; gate 1 (round-trip) + the matvec leg |
| Mode 7 MMS simplification bias | gate 3 ansatz override (mixed scales, a0>0, μ-nontrivial, 2G het) |
| Mode 8 `-O` strips bare assert | ALL new gates use `np.testing.assert_*` / `pytest.fail` — verified the DD round-trip fires under `-O` |
| Mode 9 splitting in degenerate regime | N/A — LD is a discretization, not a splitting/acceleration; no FP-coincident-regime risk |

**NO new failure mode introduced** → no skill-table append needed (the
self-improvement trigger does not fire for this plan).

---

## pytest skeletons (NEW gates — collectible, SUT-skipped)

All three NEW test files below collect today (imports of the not-yet-existing
`LinearDiscontinuous` / `OneDimPerCellWalk` are deferred INSIDE the skipped
tests so collection does not ImportError). Each is `@pytest.mark.skip(reason=
"SUT not yet implemented (#158 LinearDiscontinuous + OneDimPerCellWalk)")`.

### `tests/sn/spatial/test_linear_discontinuous.py` (GATE 1)

```python
"""GATE 1 — LinearDiscontinuous per-cell round-trip (lockstep), isolated.

Mirror of tests/sn/sweep/core/test_diamond.py::TestResidual for the LD
occupant: solve update(...) → feed the solved state to residual(...) →
residual vanishes to FP noise. Proves update/residual describe ONE per-cell
linear system (L21). LD is is_affine_scannable=False (couples two face
moments), so the round-trip MUST feed back the FULL solved state, not just
the cell average.

SUT not yet implemented (#158) — every test is skip-collectible.
"""
from __future__ import annotations

import numpy as np
import pytest

_SKIP = pytest.mark.skip(
    reason="SUT not yet implemented (#158 LinearDiscontinuous)"
)


def _ld_slab_visit_inputs(n_groups: int = 2):
    """Heterogeneous slab visit + inputs for the LD round-trip.

    Mirrors tests/sn/sweep/core/test_diamond.py::_slab_visit_inputs but for
    the LD strategy. source is already weight-normalised (Q·chord·1/Σw) per
    the CellUpdate contract.
    """
    from orpheus.geometry import BC, CoordSystem, Mesh1D, slab_streaming
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.spatial import UpstreamState
    from orpheus.sn.spatial.cell_update import CellVisit

    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 6), mat_ids=np.zeros(5, dtype=int),
        coord=CoordSystem.CARTESIAN, bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(4)
    op = slab_streaming(mesh, quad)
    st = op.streaming_terms(cell_idx=2, direction_idx=quad.N - 1)
    total_xs = np.linspace(0.6, 1.5, n_groups)          # heterogeneous (L2)
    Q = np.linspace(0.2, 2.4, n_groups)
    psi_in = np.linspace(0.05, 0.35, n_groups)
    source = Q * st.chord_length / quad.weights.sum()
    upstream = UpstreamState(spatial_upstream=psi_in, angular_upstream=None)
    visit = CellVisit(cell_idx=2, streaming_terms=st, face_area_downstream=1.0)
    return visit, total_xs, source, upstream


class TestLDRoundTrip:
    """update/residual lockstep for LinearDiscontinuous (foundation)."""

    @_SKIP
    @pytest.mark.foundation
    @pytest.mark.parametrize("n_groups", [1, 2, 4])
    def test_residual_zero_at_solved_state(self, n_groups: int) -> None:
        """residual(update(q).solved_state, q) == 0 to FP noise.

        Term activation: streaming + collision (slab has NO angular
        redistribution — angular column nulled BY GEOMETRY, not by ansatz).
        n_groups={1,2,4}: 1G is the degenerate control; 2G/4G heterogeneous
        is the real gate (vv H1/H2). LD couples cell-avg + slope; the
        round-trip feeds the FULL CellResult back to residual — a partial
        feed-back passes spuriously.
        """
        from orpheus.sn.spatial.linear_discontinuous import LinearDiscontinuous

        visit, total_xs, source, upstream = _ld_slab_visit_inputs(n_groups)
        strat = LinearDiscontinuous()
        result = strat.update(visit, total_xs, source, upstream)
        residual = strat.residual(
            result.cell_average_flux, visit, total_xs, source, upstream,
        )
        # Two-moment per-cell solve → reduction depth ~ a 2x2 solve; band
        # 1e-12 (one decade above DD's single-division 1e-13). Document if
        # LD's solve is deeper.
        np.testing.assert_allclose(residual, 0.0, atol=1e-12)

    @_SKIP
    @pytest.mark.foundation
    def test_residual_linear_in_cell_avg(self) -> None:
        """LD closure is linear → residual is linear in the probe (mode #2)."""
        from orpheus.sn.spatial.linear_discontinuous import LinearDiscontinuous

        visit, total_xs, source, upstream = _ld_slab_visit_inputs(2)
        strat = LinearDiscontinuous()
        rng = np.random.default_rng(42)
        a = rng.normal(1.0, 0.5, 2)
        b = rng.normal(2.0, 0.5, 2)
        lam = 0.37
        r_a = strat.residual(a, visit, total_xs, source, upstream)
        r_b = strat.residual(b, visit, total_xs, source, upstream)
        r_mix = strat.residual(
            lam * a + (1 - lam) * b, visit, total_xs, source, upstream,
        )
        np.testing.assert_allclose(
            r_mix, lam * r_a + (1 - lam) * r_b, rtol=1e-12, atol=1e-13,
        )

    @_SKIP
    @pytest.mark.foundation
    def test_residual_affine_in_source(self) -> None:
        """∂r/∂q = -1 (mode #3 missing-factor / convention catcher)."""
        from orpheus.sn.spatial.linear_discontinuous import LinearDiscontinuous

        visit, total_xs, source, upstream = _ld_slab_visit_inputs(2)
        strat = LinearDiscontinuous()
        rng = np.random.default_rng(7)
        probe = rng.normal(1.0, 0.5, 2)
        ds = rng.normal(0.0, 0.3, 2)
        r0 = strat.residual(probe, visit, total_xs, source, upstream)
        r1 = strat.residual(probe, visit, total_xs, source + ds, upstream)
        np.testing.assert_allclose(r1, r0 - ds, rtol=1e-12, atol=1e-13)


class TestLDTraits:
    """LD class-level traits (foundation)."""

    @_SKIP
    @pytest.mark.foundation
    def test_is_affine_scannable_false(self) -> None:
        """LD couples two face moments → NOT affine-scannable (the trait that
        routes it to OneDimPerCellWalk and AWAY from CumprodScan/ScanMarch)."""
        from orpheus.sn.spatial.linear_discontinuous import LinearDiscontinuous
        assert LinearDiscontinuous.is_affine_scannable is False

    @_SKIP
    @pytest.mark.foundation
    def test_registered_under_key(self) -> None:
        """LD self-registers under key='linear_discontinuous'."""
        from orpheus.sn.spatial.cell_update import CellUpdateBase
        from orpheus.sn.spatial.linear_discontinuous import LinearDiscontinuous
        assert CellUpdateBase.registry["linear_discontinuous"] is LinearDiscontinuous
```

### `tests/sn/sweep/core/test_unified_sweep_dispatch.py` ADDITIONS (GATE 2 + the gate-5 selection pin)

```python
# ─── append to tests/sn/sweep/core/test_unified_sweep_dispatch.py ───
# GATE 2 (seam-opening control) + GATE 5 selection-no-regression.

_SKIP_LD = pytest.mark.skip(
    reason="SUT not yet implemented (#158 OneDimPerCellWalk + LinearDiscontinuous)"
)


def _ld_slab_sn_mesh(nx: int = 8, length: float = 1.0):
    """Slab SNMesh carrying an LD cell-update (the non-affine-scannable mesh)."""
    from orpheus.sn.spatial.linear_discontinuous import LinearDiscontinuous
    mesh = Mesh1D(
        edges=np.linspace(0.0, length, nx + 1), mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN, bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    return SNMesh(mesh, quad, placeholder_materials(),
                  cell_update=LinearDiscontinuous())


class TestOneDimPerCellSeam:
    """GATE 2 — the per-cell-walk frame is correct independently of LD.

    Routes DD itself through the NEW OneDimPerCellWalk per-cell path and
    asserts it reproduces DD's fast-scan (CumprodScan) answer. DD is the
    structurally-independent reference (its scan answer is bit-identity-pinned
    by the regression snapshots + the homogeneous k_inf analytic limit).
    """

    @_SKIP_LD
    @pytest.mark.foundation
    def test_ld_mesh_selects_one_dim_per_cell_walk(self):
        """SELECTION + mutual-exclusion (vv anti-pattern #11: BOTH a
        must-select and a must-NOT-select case)."""
        from orpheus.sn.loss_representation import (
            CumprodScan, OneDimPerCellWalk, default_for,
        )
        ld_mesh = _ld_slab_sn_mesh()
        dd_mesh = _slab_sn_mesh()
        if not isinstance(default_for(ld_mesh), OneDimPerCellWalk):
            pytest.fail(
                f"LD 1-D mesh → {type(default_for(ld_mesh)).__name__}, "
                "expected OneDimPerCellWalk"
            )
        # mutual exclusion on the is_affine_scannable trait
        if CumprodScan.supports(ld_mesh).ok:
            pytest.fail("CumprodScan must REJECT a non-affine-scannable mesh")
        if OneDimPerCellWalk.supports(dd_mesh).ok:
            pytest.fail("OneDimPerCellWalk must REJECT an affine-scannable (DD) mesh")

    @_SKIP_LD
    @pytest.mark.foundation
    def test_dd_via_per_cell_walk_equals_cumprod_scan(self):
        """EQUIVALENCE CONTROL — DD swept via OneDimPerCellWalk reproduces
        DD's CumprodScan answer to principled-equivalence (nulp=nx).

        TOLERANCE: assert_array_almost_equal_nulp(nulp=nx), NOT array_equal.
        The two paths compute the SAME DD algebra in real arithmetic over a
        DIFFERENT IEEE-754 reduction tree (closed-form parallel-prefix
        cumprod vs sequential per-cell loop); addition is non-associative.
        Drift bounded by (reduction depth = nx) × ULP — vv §Bit-identity:
        principled (named per-cell ψ̄ intermediate), reference is
        bit-identity-pinned DD scan, drift dimensionally explained. See the
        #206 Phase-C density-vs-scan nULP precedent.

        2G heterogeneous field (L2 — 1G is flux-shape degenerate; a
        group-coupling drift in the frame must be observable).
        """
        from orpheus.sn.loss_representation import CumprodScan, OneDimPerCellWalk
        from orpheus.transport.source_sinks import AngularSourceSink
        from orpheus.transport.fields.boundary_flux import BoundaryFlux

        nx = 16
        dd_mesh = _slab_sn_mesh(nx=nx)            # DD occupant (affine-scannable)
        ng, N = dd_mesh.ng, dd_mesh.quad.N
        spatial = dd_mesh.spatial_shape
        rng = np.random.default_rng(2026)
        # heterogeneous (per-cell) σ_t, non-flat per-ordinate source.
        sig_t = rng.uniform(0.6, 1.5, (ng, *spatial))
        Q = rng.uniform(0.1, 2.0, (N, ng, *spatial))
        src = AngularSourceSink.from_per_ordinate(Q, dd_mesh)   # adjust to real API
        bflux = BoundaryFlux.zeros_on(dd_mesh)

        psi_scan, out_scan = CumprodScan(dd_mesh).sweep(
            src.magnitude, sig_t, bflux,
        )
        psi_walk, out_walk = OneDimPerCellWalk(dd_mesh).sweep(
            src.magnitude, sig_t, bflux,
        )
        np.testing.assert_array_almost_equal_nulp(psi_walk, psi_scan, nulp=nx)
        np.testing.assert_array_almost_equal_nulp(out_walk, out_scan, nulp=nx)

    @_SKIP_LD
    @pytest.mark.foundation
    def test_dd_mesh_still_selects_scan_after_registry_edit(self):
        """GATE 5 selection-no-regression — inserting OneDimPerCellWalk into
        LOSS_REPRESENTATIONS did NOT shadow the scan reps for DD meshes."""
        from orpheus.sn.loss_representation import (
            CumprodScan, ScanMarch, default_for,
        )
        if not isinstance(default_for(_slab_sn_mesh()), CumprodScan):
            pytest.fail("DD 1-D mesh no longer selects CumprodScan (registry regression)")
        if not isinstance(default_for(_2d_sn_mesh()), ScanMarch):
            pytest.fail("DD 2-D mesh no longer selects ScanMarch (registry regression)")
```

### `tests/sn/verification/mms/test_mms_ld_slab.py` (GATE 3)

```python
"""GATE 3 — LinearDiscontinuous slab-Cartesian MMS, O(h²) spatial order.

The interference-free spatial oracle: slab Cartesian has NO space-angle
coupling, so the convergence ladder isolates the LD spatial order. MMS is
the math/flux-shape pillar (NOT eigenvalue — vv hierarchical taxonomy).

ANSATZ OVERRIDE (vv §MMS operational rules, Mode 7): the existing slab MMS
all use A(x)=sin(πx/L) — vanishing-at-boundaries, isotropic-in-μ, the
simplest trig. For LD that under-stresses the closure. This gate's ansatz:
  ψ_{n,g} = [a0_g + a1_g·sin(πx/L) + a2_g·cos(3πx/L)
             + μ_n·(b0_g + b1_g·sin(2πx/L))] / W
with a0_g > 0 (NON-vanishing at boundaries → real vacuum-inflow BC test),
mixed scales (k=1 + k=3 high harmonic), angularly non-trivial (μ·(…)),
≥2 groups heterogeneous Σ_t(x). Source derived symbolically (Branch 1,
SNSlabLDStressMMSCase in orpheus/derivations/continuous/mms/sn.py).

SUT not yet implemented (#158) — skip-collectible.
"""
from __future__ import annotations

import numpy as np
import pytest

_SKIP = pytest.mark.skip(
    reason="SUT not yet implemented (#158 LinearDiscontinuous + SNSlabLDStressMMSCase)"
)


def _l2_error(phi_num, phi_ref, widths):
    diff = phi_num - phi_ref
    return float(np.sqrt(np.sum(widths * diff * diff)))


@_SKIP
@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("ld-cartesian-1d", "ld-slab", "transport-cartesian")
def test_sn_1d_slab_ld_mms_converges_second_order():
    r"""LD slab MMS shows measured O(h²) on the angular-stress ansatz.

    Asserts BOTH the rate (orders > 1.9 — LD is exactly 2nd order, NOT
    higher) AND the converged VALUE (magnitude band against the imposed
    phi_exact — vv §5: rate alone is necessary, never sufficient). The
    structurally-independent reference IS the manufactured solution (that is
    what makes MMS a flux-shape pillar). The k=3 high harmonic gives the
    ladder real dynamic range so the coarse-mesh error is not floor-limited.
    """
    from orpheus.derivations.continuous.mms.sn import build_1d_slab_ld_stress_mms_case
    from orpheus.sn.geometry import SNMesh
    from orpheus.sn.solver import SNSolver  # direct construction until the
    # cell_update kwarg is threaded through solve_sn_fixed_source
    from orpheus.sn.spatial.linear_discontinuous import LinearDiscontinuous

    case = build_1d_slab_ld_stress_mms_case()   # 2G heterogeneous, a0>0, mixed scales
    n_cells = [10, 20, 40, 80]
    errors = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        materials = case.build_materials(mesh)
        Q = case.external_source(mesh)
        sn_mesh = SNMesh(
            mesh, case.quadrature, materials,
            cell_update=LinearDiscontinuous(),
        )
        # If solve_sn_fixed_source(cell_update=...) lands, prefer it:
        #   result = solve_sn_fixed_source(materials, mesh, case.quadrature, Q,
        #       cell_update=LinearDiscontinuous(), max_inner=500, inner_tol=1e-13)
        solver = SNSolver(sn_mesh, inner_solver="source_iteration",
                          scattering_order=0, max_inner=500, inner_tol=1e-13)
        result = solver.solve_fixed_source(Q)   # adjust to the real entry
        phi_num = result.scalar_flux.values[0, :]
        phi_ref = case.phi_exact(mesh.centers, g=0)
        errors.append(_l2_error(phi_num, phi_ref, mesh.widths))

    errors = np.asarray(errors)
    orders = np.log2(errors[:-1] / errors[1:])
    assert np.all(orders > 1.9), (
        f"Expected LD O(h²), got orders={orders} from errors={errors}"
    )
    # Converged VALUE (NOT just rate). Floor well above inner_tol; ceiling
    # below the coarse-mesh error so the band has real dynamic range.
    assert 1e-9 < errors[-1] < 1e-3, (
        f"LD finest-mesh error {errors[-1]:.3e} outside the value band — "
        "rate may be O(h²) to the WRONG limit (vv §5)"
    )


@_SKIP
@pytest.mark.l1
@pytest.mark.verifies("ld-cartesian-1d")
def test_sn_1d_slab_ld_mms_second_group_coupling():
    r"""The g=2 MMS source couples to c_1 via downscatter — a wrong/transposed
    scatter assembly (mode #6) produces an incorrect φ_2 the O(h²) rate (on
    g=2 specifically) catches. Same ansatz, asserted on the g=1 slice."""
    from orpheus.derivations.continuous.mms.sn import build_1d_slab_ld_stress_mms_case
    # ... identical ladder, phi_ref = case.phi_exact(centers, g=1), assert
    # orders>1.9 on the g=1 (second-group) error.
    pytest.skip("body mirrors the g=0 test on the g=1 slice")
```

### GATE 4 — no new mandatory test (run the existing #233 gate with LD)

The existing
`tests/sn/verification/mms/test_curvilinear_pole_cell_characterization.py`
(4 tests, lower-bound `> 0.8` pole order, no upper bound) is ALREADY
LD-ready by construction. The LD task RUNS it on an LD curvilinear mesh and
records the measured pole order in the closeout. OPTIONAL forward tripwire
(staged xfail, NON-strict — does NOT gate the plan):

```python
@pytest.mark.skip(reason="SUT not yet implemented (#158 LinearDiscontinuous)")
@pytest.mark.xfail(strict=False,
    reason="#235 angular floor co-limits the sphere pole until S32/M-M-interp")
@pytest.mark.l1
@pytest.mark.catches("ERR-059")
def test_sphere_pole_LD_lifts_toward_second_order():
    """ASPIRATIONAL forward tripwire — LD should lift the sphere pole order
    above first order; GATED by #235 (out of scope) so staged xfail. Flips to
    xpass when #235 lands WITH LD. Do NOT make strict; do NOT let gate-4
    pass/fail depend on it."""
    ...  # run build_spherical_mms_case with cell_update=LinearDiscontinuous(),
         # ladder [40,80,160,320], assert pole_orders > 1.5
```

---

## Summary for the dispatcher

- **Slab-Cartesian LD MMS design:** override the universal `sin(πx/L)`
  simplification bias. Ansatz `ψ = [a0 + a1·sin(πx/L) + a2·cos(3πx/L) +
  μ·(b0 + b1·sin(2πx/L))]/W`, **a0>0** (non-vanishing-at-boundary → real
  vacuum BC test), mixed scales (k=1+k=3), μ-non-trivial, 2G heterogeneous
  Σ_t(x). Source symbolic (new `SNSlabLDStressMMSCase`). Assert orders>1.9
  AND the converged value band (rate ≠ correctness, vv §5).
- **Seam-control tolerance:** `assert_array_almost_equal_nulp(nulp=nx)` —
  principled-equivalence, NOT `array_equal`. DD-via-`OneDimPerCellWalk` and
  DD-via-`CumprodScan` compute the same algebra over different IEEE-754
  reduction trees (sequential loop vs parallel-prefix cumprod); drift
  bounded by reduction-depth nx × ULP. (#206 Phase-C precedent.)
- **DD-no-regression is already pinned by EXISTING tests:** the strict gate
  is `pytest tests/sn/sweep/core tests/sn/solve -W
  "error::tests.sn.regression._regression_assert.DriftWarning"` (NOT
  `tests/sn/regression/` — its conftest forces `always::DriftWarning`,
  inerting `-W error`). Currently 505 passed / 1 skipped / 4 xfailed.
  Plus a NEW selection-no-regression assertion that DD meshes still pick
  CumprodScan/ScanMarch after the registry edit.
```
