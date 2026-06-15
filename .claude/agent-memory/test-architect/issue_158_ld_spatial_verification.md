---
name: issue-158-ld-spatial-verification
description: Verification plan + pytest skeletons for the Linear-Discontinuous (LD) SN cell-update scheme + the OneDimPerCellWalk dispatch seam (#158). 5 gates, L17 crosswalk, slab-Cartesian LD MMS ansatz (angular-stress override of the simplification bias), seam-control nULP tolerance, DD-no-regression strict invocation.
metadata:
  type: project
---

# Issue #158 — Linear-Discontinuous (LD) spatial scheme + dispatch-seam verification

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
