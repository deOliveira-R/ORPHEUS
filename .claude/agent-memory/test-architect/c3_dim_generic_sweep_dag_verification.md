---
name: c3-dim-generic-sweep-dag-verification
description: C3 carve V&V plan (EXPANDED — spine+optimization layering). Dim-generic SweepDependencyGraph spine = ORACLE; cumprod = d=1 opt; _MovingFrontier = d≥2 opt. wavefront(1-D)≡cumprod(1-D) @ nulp anchored on k_inf=1.875; perf measurement; d=1 real spine pins; x↔y-swap gap.
metadata:
  type: project
---

C3 makes the SN sweep DAG + diamond cell kernel + octant enumeration
DIMENSION-GENERIC (`from_cartesian_2d`→`from_cartesian(d)`; explicit
`sx`/`sy`→`denom=σ_t+Σ_axis s_axis`; octant→`itertools.product((-1,+1),repeat=d)`).

⭐ **ARCHITECTURE EXPANDED (user override — settled, do not re-litigate): SPINE +
PLUGGABLE-OPTIMIZATIONS.** The wavefront DAG full-field walk
(`from_cartesian(shape)` d=1,2,3 + d-generic `cell_kernel_batch`) is the
DIMENSION-AGNOSTIC SPINE and the VERIFICATION ORACLE — ONE path for every d. Two
optimizations plug in WITHOUT changing the spine's value: **(a) cumprod
parallel-prefix scan** (`ordinate_scan`/`_sweep_1d_unified`, Blelloch §1.5
`cumprod_a·(psi_0+cumsum(b/cumprod_a))`) is the **d=1 OPTIMIZATION** (default-selected
for d=1 production); **(b) `_MovingFrontier` storage-B window** is the **d≥2
OPTIMIZATION** (pinned ≡ full-field by `tests/sn/sweep/cartesian_2d/test_2d_full_field_oracle.py`).
So 1-D is NO LONGER "untouched" — it is a REAL spine compute path validated against
its cumprod optimization. 2-D stays BIT-IDENTICAL (refactor); synthetic-3-D = SHAPE
ADMISSION only.

**Why:** N-D layout foundation (eliminate phantom ny=1); C3 is the sweep-DAG arm.

**ORIGINAL sections STILL STAND** (verified live, branch HEAD f994a0c, 2026-06-09):
§A 2-D bit-identity moat (hand-loop oracle `TestApplyMatchesLegacyInlined`), B1-B5
synthetic-3-D pins, C d=2 legacy-equivalence pin, §D gate ladder, §E standing reds.
ADDED below: NEW-1 (wavefront≡cumprod oracle), NEW-2 (perf measurement), NEW-3 (d=1
real spine pins B6), NEW-4 (d≥2 oracle confirmation).

---

## NEW-1..4 (EXPANDED scope — 2 UNCOMMITTED stub files, for main-agent review)

**NEW-1 wavefront(1-D)≡cumprod(1-D) oracle** — file
`tests/sn/sweep/core/test_wavefront_cumprod_equivalence.py` (NEW). DECISION:
**PRINCIPLED-EQUIVALENCE at nulp, NOT bit-identity.** The cumprod scan
(parallel-prefix reduction TREE over all nx cells) and the wavefront
forward-substitution walk (nx sequential cell-by-cell) compute the SAME affine
recurrence `psi_out[i]=2psi_bar−psi_in, psi_bar=(Q+s·psi_in)/(σ_t+s)` in DIFFERENT
FP association orders → NOT bit-identical (IEEE-754 non-associative). All 3 vv
criteria hold: (1) intermediates named (psi_out[i], psi_bar); (2) anchor =
ANALYTICAL `k_inf=λ_max(A⁻¹F)=1.8750000000000009` mixture A 2g (transfer-matrix
`kinf_and_spectrum_homogeneous(SigT,SigS[0],SigP,chi)`, VERIFIED real; reflective
slab `solve_sn` → 1.8750000086 to ~1e-8) — old-vs-new ULP necessary-never-sufficient,
so the cumprod path itself terminates in closed form; (3) single-sweep drift bounded
`(reduction depth=nx)×ULP`. **`_NULP_BOUND=128`** (8×16 safety; carve MEASURES actual
once spine lands, tighten/loosen+RE-JUSTIFY). Config: mixture A 2g (asymmetric
downscatter SigS[0,1]=0.1, ≥2G H1), heterogeneous random Σ_t(ng,nx) + non-uniform
`from_isotropic` source + NON-ZERO BoundaryFlux inflow ⟹ psi_in≠psi_out everywhere
(redistribution out of cancellation, H2). Gate `assert_array_almost_equal_nulp`,
`-O`, foundation. Wavefront-1D driven via the ONLY adapter `_wavefront_1d_sweep`
(builds `from_cartesian((nx,))` + d-generic single-axis kernel; raises
`_SpineNotLanded` pre-carve → row `xfail(strict=False)`). The analytical-anchor row
(`test_cumprod_path_hits_analytical_kinf`) runs TODAY (live cumprod) ⟹ file not inert.
⚠ Bit-identity NOT claimed achievable: the d=1 wavefront is forward-substitution, not
realizable AS the cumprod tree (different reduction), so principled-equiv is correct,
not a fallback.

**NEW-2 perf measurement** — `test_cumprod_faster_than_wavefront_d1_loose_tripwire`,
SAME file. `@pytest.mark.slow` + LOOSE final-gate tripwire (NOT hard correctness gate;
wall-clock noisy — Phase-5b Leg-2 loose-discipline precedent). nx=4096 long chain;
min-over-5-repeats; PRINTS `speedup=X×` (the real deliverable); tripwire `ratio<2.0 →
pytest.fail` (catches "wavefront accidentally became d=1 default, 1-D got slow" — the
2× bound is FAR below expected cumprod O(N/logN)-vs-sequential ratio). xfail pre-carve.

**NEW-3 d=1 real spine pins (B6)** — ADDED to
`tests/sn/sweep/core/test_sweep_graph_nd_admission.py`. d=1 is now a REAL wavefront
compute path: **B6a** octant count 2=2^1 (`[(-1,),(+1,)]`, exact order, PASSES today);
**B6b** `from_cartesian((nx,))` builds nx levels, ONE cell/level (singleton chain;
`Σ(n−1)+1=(nx−1)+1=nx`; (+1,) ⟹ cell i==ℓ; total coverage), nx∈{5,8,12} both octants;
**B6c** d-generic kernel at d=1 → `denom=σ_t+s_x` single-axis, hand oracle (PASSES
today). B6b xfail pre-carve (needs builder).

**NEW-4 d≥2 oracle CONFIRMED EXISTS** — `tests/sn/sweep/cartesian_2d/test_2d_full_field_oracle.py`
(LIVE, foundation). 2 parametrized tests × 4 CASES `[(8,8,4,2,reflective),(8,8,4,2,vacuum),
(12,7,6,2,reflective),(5,9,4,4,reflective)]` (≥2G + het random Σ_t + random inflow + NON-SQUARE
12×7/5×9). Key assertions PASTED:
`test_sweep_window_equals_full_field_end_to_end`: `_sweep_2d_wavefront`(windowed) vs
`_sweep_2d_full_field`(oracle) → `assert_array_equal(ang_w,ang_f)` + `(scal_w,scal_f)` +
`(bf_win.values,bf_full.values)` post-sweep boundary trace, BIT-FOR-BIT.
`test_matvec_window_equals_full_field_end_to_end`: `L._apply_2d_cartesian`(windowed) vs
`L._apply_2d_cartesian_full_field`(oracle) → `assert_array_equal` bulk residual +
boundary-block residual. This REMAINS the gate pinning `_MovingFrontier` (storage-B,
d≥2 opt) against the d-generic spine. Both share `DiamondDifference.cell_kernel_batch`
⟹ math cannot drift; only storage walk + boundary bookkeeping differ.

**STUB STATE (both UNCOMMITTED, -O verified):** admission file = 6 passed + 26 xfailed
(builder-gated via shared `_needs_spine = xfail(strict=False, not hasattr(...,'from_cartesian'))`);
equivalence file = 2 passed (anchor + import guard) + 3 xfailed (spine rows). Whole
`tests/sn/sweep/core/` = **401 passed, 1 skipped, 33 xfailed, 0 failed (6.1s, -O)** —
no regression. xfails FLIP to xpass automatically when `from_cartesian` lands.

⭐ **API SPELLINGS CONFIRMED LIVE (only adapters change if carve differs):** cumprod path
= `transport_sweep(AngularSourceSink, sig_t(ng,nx), SNMesh, BoundaryFlux)` →
`_sweep_1d_unified`→`ordinate_scan` (returns `(ang(N,ng,nx,1), scal(ng,nx,1))`);
`AngularSourceSink.from_isotropic(iso(ng,nx), sn_mesh)`; `BoundaryFlux.zeros_on/from_mesh`
+ `.face_view(face)` seed; `Mesh1D(edges,mat_ids,coord,bc_left,bc_right)`;
`Quadrature.gauss_legendre(n_ordinates=8)`. ASSUMED-NOT-LANDED: `from_cartesian((nx,),
label=<d-tuple>)` (spine builder; only `from_cartesian_2d` exists at HEAD).

---

**How to apply (the moat that already exists, branch `worktree-sn-nd-layout`):**

- **The load-bearing d=2 oracle is `tests/sn/sweep/core/test_sweep_graph.py::
  TestApplyMatchesLegacyInlined::test_per_cell_loop_equivalence`** — a
  STRUCTURALLY-INDEPENDENT hand-written per-ordinate Python loop
  (`_hand_run_legacy_inlined`) vs `graph.apply`, on NON-SQUARE (3,4)/(4,3),
  multi-group, randomized per-ordinate inputs. `assert_array_equal` (angular+faces,
  fires under -O) + `assert_array_almost_equal_nulp(nulp=128)` (scalar). This is
  a hand-loop-vs-vectorized gate, NOT a snapshot — it subsumes the sum-order +
  sign-tuple bit-identity claim at d=2. Plus `TestAssertCellCoverage` /
  `TestAssertTopologicallySorted` / `TestAssertFacePairingConsistent` pin DAG
  structure on (2,3)/(4,5)/(5,4) — BUT those use BARE `assert` ⟹ NO-OP under -O
  (vv Mode 8). Run that file WITHOUT -O or migrate to np.testing.
- **The end-to-end moat:** `test_2d_octant_sweep_equivalence.py` (6 snapshot
  cases + 1 closed-form sentinel anchor, all `nx=ny=3` SQUARE) + `test_dd_regression.py`
  (`assert_regression` = SAFETY×conv_tol). The DD-regression 2-D cases ARE
  non-square (`_cartesian_2d_het_si` 8×4, `_cartesian_2d_p1_aniso_het_si` 8×4,
  2G het reflective, P1 aniso) ⟹ they catch an x↔y swap via converged flux.

**THE GAP (x↔y axis swap, AI failure mode 2):** the octant-snapshot moat is
SQUARE-ONLY (nx=ny=3) ⟹ blind to a swap. Closed by (a) the existing non-square
hand-loop oracle above + the non-square DD-regression cases, and (b) a NEW
non-square d=2 equivalence pin (5,7)/(4,3) in the stub.

**Synthetic-3-D exact numbers (computed, put in asserts):** octant count
`2^d` (d=3→8, distinct); n_levels = `Σ(n_axis−1)+1` (d=2 `nx+ny−1` generalised);
(2,3,2)→5 levels, per-level cell counts `[1,3,4,3,1]`, Σ=12; (3,3,3)→7 levels
`[1,3,6,7,6,3,1]`; (4,2,3)→7; level membership `i+j+k==ℓ`; diamond
`denom=σ_t+Σ_axis s_axis` (3-tuple `(.3,.7,1.1)`+σ_t=.5→2.6); per-axis closure
`psi_out_axis=2·psi_avg−psi_in_axis` independent per axis. 3-D pins are
`@pytest.mark.foundation` (no `:label:`), np.testing/pytest.fail (fire under -O).

**Standing reds (all held as `xfail`, gate stays GREEN; confirmed clean HEAD
f994a0c):** #206 cyl-matvec (27 xfailed in `test_unified_matvec_cylinder.py`),
#195 MMS@160 (`xfail strict=True`), #212 `continuous_get` hang (DESELECT
`test_keff_slab.py::test_heterogeneous_absolute_keff`, NO `continuous_get`).
**#214 is CLOSED** (no remaining skip/xfail). Fast ladder baseline:
`primitives/ sweep/core/ sweep/slab/ sweep/curvilinear/ solve/` =
**862 passed, 1 skipped, 31 xfailed, 0 failed** (858s; xdist deadlocks per
`sn_taxonomy_reorg_mapping` so sequential).

**UNCOMMITTED stub written (for main-agent review):**
`tests/sn/sweep/core/test_sweep_graph_nd_admission.py` — B1 octant count,
B2 d=3 build, B3 level count, B4 hyperplane membership+coverage, B5 3-axis
kernel reduction, C non-square d=2 equivalence vs legacy. Wired to the ASSUMED
C3 API `from_cartesian(shape, *, label=<sign-tuple>)` via 3 adapter helpers
(`_build_nd`/`_levels_of`/`_octant_labels`) — only those update if the final
spelling differs.
