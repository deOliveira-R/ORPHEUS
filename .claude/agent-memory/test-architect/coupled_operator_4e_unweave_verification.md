---
name: coupled-operator-4e-unweave-verification
description: The Phase-C step-4e "WALK UN-WEAVE" gate-spec delta for the Coupled Block Operator campaign (branch refactor/sn-walk-unification @ 9fdbd47, B.2d COMPLETE, wall 3733/0 · 37 xfailed). 4e = the EXECUTIONER: (i) EXTRACT the four welded ray legs of the fused (L+C) walk (A_AB seed feed fwd/transpose, A_BB ray march solve/solve_transpose) into RadialCharacteristicOperator/RadialCharacteristicSeeding (or a shared kernel) + DELETE the welded copies — principled-equiv on the carrying sphere, bit-exact on non-carrying slab/cyl CONTROL; (ii) WALK-SLOT REWRITE split-native — retire to_unified/from_unified + the unified RadialCharacteristicField leaf + unified RadialCharacteristicSpace + the F2 space-identity aliasing. 4e does NOT flip grid.is_invertible (step 5 #41). THE headline finding: the Mode-12 ERR-067 seed-row-flip closure (test_282::test_mode12_g_reciprocity_catches_a_seed_row_flip) has its tooth on _seed_rows_forward — a welded leg 4e DELETES → silently toothless unless re-aimed onto RadialCharacteristicSeeding. Load-bearing survivors: row-6 dense-LU oracle (the durable principled-equiv floor 5.5e-16), test_282 cold-residual <1e-11 (the L18 acceptance number), G-d3.1 φ=Q/Σ_t (redistribution catcher), affine-carve bit-exact (the (A,A)-untouched CONTROL). Sequencing: e1 walk-slot rewrite (bit-id-ABLE, keeps every frozen 0-ULP baseline valid) THEN e2 extraction (principled-equiv, sphere re-captures anchored to row-6). Sibling of step4 / b2c / b2d memos (L13/L16/L18/L19).
metadata:
  type: reference
---

# Verifying step 4e — THE WALK UN-WEAVE (the EXTRACT executioner)

**When this recipe applies.** A coupled-block campaign has (B.2a) landed the
N-general `CoupledField`/`CoupledOperator` machinery, (B.2b/c) re-typed the
off-diagonal blocks + built `build_coupled_system`, and (B.2d) made the driver
block-native + **atomically evicted the ray from FullField** (FullField is 2-block;
the driver iterate is `CoupledField[ψ_A, ψ_B]`; the six walk signatures speak
explicit leaf kwargs). What REMAINS welded: the fused `(L+C)` walk still (i)
**inlines** the four ray legs (the engines `RadialCharacteristicOperator`/`…Seeding`
EXIST and are proven kernel-equivalent, but the PRODUCTION walk routes around them),
and (ii) marches a **unified flat buffer** via `cells_view`/`corner_view`, with the
`CoupledInvertibleOperator` (M) round-tripping `to_unified`/`from_unified` at the
block boundary. **4e is the executioner:** route the legs through the engines +
DELETE the welded copies (EXTRACT), and make the walk march split-native + retire
the unified leaf/space/bridge (WALK-SLOT REWRITE). Governing memos (same dir):
`coupled_operator_step4_verification.md` (the RE-SCOPED 4d.x — the row-6 oracle +
A2 reciprocity + the assemble≡probe floor), `coupled_operator_b2d_driver_eviction_
verification.md` (d1/d2/d3 — the current gate names + the F1–F5 findings).

**Verified against HEAD `9fdbd47`.** Every `file:line` below is read off that tree
(⚠ best-estimate for production internals — the walk methods drift; grep the symbol,
not the line). Canonical `.venv/bin/python -O -m pytest -p no:xdist --timeout=300
-p no:cacheprovider` SERIAL; every tooth in-process `monkeypatch` (NEVER `git checkout`
an uncommitted file — L4); every gate fires under `-O` (`np.testing`/`pytest.fail`,
never bare `assert` outside a pytest-rewritten test body — vv Mode 8).

**The step-5 boundary (do NOT cross).** 4e does NOT change the solve's joint-pair
REQUIRE contract and does NOT flip `grid.is_invertible`. The block-solve spelling +
the A2a inverse-arm flip (`test_a2a_grid_swap_law_inverse_arm`, currently
`xfail(strict=False)` — the wall's 37th xfail) are **step 5 (#41)**. Through 4e that
arm STAYS `xfail`; a 4e that greened it would be out of scope.

## The object (arithmetic-verified against the code)

**What already exists (the extraction is a RE-ROUTE, not a build).** The engines are
built and proven kernel-equivalent to the welded legs:

- **A_BB forward** (`RadialCharacteristicOperator.apply`) routes through the SHARED
  `radial_characteristic_forward_residual`; **the walk's ψ½ forward rows route
  through the SAME kernel** — pinned TODAY by the Mode-11 spy pair
  `test_walk_forward_routes_through_the_shared_kernel` + `…apply_routes_through_the_
  shared_forward_kernel` (psi_half `TestA_BB_Forward` :1285/:1273).
- **A_BB solve march** = `carlson_inward_sweep_from_source`/`…_transpose`
  (`radial_characteristic.py:510/515`), the same recurrence the walk marches.
- **A_AB fwd/transpose** (`RadialCharacteristicSeeding.apply`/`.apply_transpose`) ≡ the
  in-sweep seed injection / `seed_cells_bar` pullback — pinned TODAY BIT-IDENTICALLY
  by `test_apply_matches_the_in_sweep_seed_injection` (:2525, Mode-11 spy on
  `precompute_psi_state`/`cell_contribution`, bulk `psi_view` zeroed) +
  `test_apply_transpose_is_the_in_sweep_seed_cells_bar` (:2591, spy on
  `angular_adjoint`).

So the "welded copies" are the **inline CALL SITES** in the walk
(`loss_representation/__init__.py`: `_seed_rows_forward:2853`, `_seed_rows_transpose:
2875`, `_apply_walk:3050`, `loss_action_transpose:3411`, `_run:3739`,
`_run_transpose:4387`) — NOT divergent kernels. 4e replaces each call site with a
call to the engine object (or a shared free-function kernel) and DELETES the inline
leg + the unified marching substrate.

**The reassociation source (HAZARD 1, the principled-equiv driver).** The shared
closure kernel `precompute_psi_state`/`cell_contribution` computes **A_AA + A_AB
SUMMED in one pass**. Extracting A_AB as a SEPARATE `RadialCharacteristicSeeding.apply`
means a second kernel pass (psi_view zeroed) + a field SUM `A_AA·ψ + A_AB·seed` — a
different reduction tree from the one-pass fused `(A_AA+A_AB)`. That is the ~5.5e-16
principled-equiv the carrying sphere pays. **The (A,A)-only path (no ray leg) is
untouched** — this is the honest claim boundary (§ Ask 2, H1).

---

## Ask 1 — THE SURVIVORS LEDGER

Legend: **STAYS** = green-unchanged, no edit; **RE-POINT** = accessor/import change,
value HOLDS; **RE-CAPTURE** = principled re-baseline (anchored, not old-vs-new ULP);
**RE-AIM** = the tooth's mutation target dies → re-point the mutation or go silently
toothless.

### A. The #284 / regression floor (psi_half `TestRegressionFloor`, :254-420)

Landed as O(1)-threshold structural pins + one principled-equiv oracle → **ROBUST to
4e's reduction-tree drift by design**. All six already ride `_joint_M`/`_coupled_
template` (B.2d re-point).

| gate | :line | claim | 4e fate | why |
|---|---|---|---|---|
| `test_loss_operator_is_block_triangular_in_the_ray` (#284 cert) | 265 | within (L+C): `A_sb=A_st<1e-12`; `A_bs>1`, `A_ss>1` | **STAYS** (LIVE 4e catcher) | block-triangularity is STRUCTURAL; O(1) thresholds absorb ULP drift. The step-4 memo's "re-expresses at 4e" = interpretation shift (augmented-(L+C)→CoupledOperator block), NOT a value re-baseline. Its tooth (a walk bulk→seed back-edge) is a genuine 4e catcher — an extraction that let A_BB read a bulk moment reds `A_sb` O(1). |
| `test_bulk_to_ray_coupling_lives_in_the_lagged_A_BA_gain` | 285 | `S_sb=0` structural; `A_BA_sb>1e-6`; N's (A,B)=∅ | **STAYS** | reads the gain grid N (not the walk M); 4e touches only M. |
| `test_outer_si_splitting_rate_is_bounded_by_scattering_ratio` | 334 | `ρ(M⁻¹N)<c+1e-6` | **STAYS** | ρ is a spectral bound of `M⁻¹N`; densifies `system.resolvent`(M)+`gains`(N). M's reduction tree drifts ~1e-15 ≪ the O(1) bound. |
| `test_folded_ray_seed_is_a_direct_solve_zero_spectral_radius` | 354 | `ρ(lag)<1e-10` (nilpotent) | **STAYS** (LIVE catcher) | `A_sb=0` ⟹ no back-edge ⟹ nilpotent. Tooth (a back-edge) is a 4e catcher. |
| `test_welded_sweep_is_exact_direct_inverse` | 377 | `‖solve∘apply − ψ‖<1e-12` | **STAYS** + doc reword | post-4e "welded"→"extracted"; solve∘apply is principled-equiv ~1e-15 ≪ 1e-12. Reword the docstring "welded sweep IS the exact inverse" → "extracted block solve" (tense update, no re-baseline). |
| `test_extract_to_dense_is_principled_equivalent_not_bit_identical` (**row 6, THE oracle**) | 397 | `‖extracted_solve − LU(_dense(M.apply))‖<1e-11` (~5.5e-16) | **STAYS — becomes literal** | the durable principled-equiv floor. Pre-4e "welded vs dense-LU"; post-4e "extracted vs dense-LU-of-extract" — it was DESIGNED as the anchor every EXTRACT step pins against. **The re-capture anchor for the sphere walk_matvec row (§B).** |

### B. The frozen walk_matvec baselines (`tests/sn/regression/test_walk_matvec_baselines.py` + `_generate_walk_baselines.py`)

`assert_regression(kind="direct", reduction_depth=1)` → `nulp=1` (tightest). Measures
`fwd_bulk/fwd_trace/adj_bulk/adj_trace` = the System-A OUTPUT (the ray FEEDS the bulk
via A_AB on carrying meshes).

| snapshot | 4e fate | bar | anchor |
|---|---|---|---|
| `walk_matvec_slab_2g` | **STAYS bit-exact** (H8 CONTROL) | `nulp=1` | Cartesian, no ray — the reduction tree is untouched. |
| `walk_matvec_cyl_2g` | **STAYS bit-exact** (CONTROL) ⚠ verify | `nulp=1` | cyl is NON-carrying by the #229 `τ_raw=0` clamp ⟹ `_draw_with_seed` sets `seed_leaf=None` ⟹ seedless path. **FINDING C1**: the generator's own docstring says "sphere/cylinder rows are RE-CAPTURED" — that prose predates the #229 non-carrying fact. Confirm `sn.radial_characteristic_space is None` for the generator's cyl (`_make_cyl`); if non-carrying (expected) → bit-exact CONTROL; if a product-cyl carries → RE-CAPTURE like the sphere. Do NOT trust the docstring — check the object. |
| `walk_matvec_sphere_2g` | **RE-CAPTURE** (the ONLY carrying member) | re-baseline the `.npz`, drift ~5.5e-16 documented | at 4e the A_AA+A_AB fold reassociates → `fwd_bulk`/`adj_bulk` drift > 1 nulp. **RE-CAPTURE anchored to `_dense(lc.apply)` probe** (structurally-independent — NOT old-vs-new ULP, the L18/snapshot-inheritance-needs-anchor rule): the re-captured sphere values MUST agree with `_dense(lc.apply)` at ≤ 5.5e-16 (the row-6 floor). The generator ALSO re-points its `radial_characteristic_flux=` leg to the split-native surface (the unified leaf retires). |

### C. The affine-carve bit-identity baselines (the (A,A)-untouched CONTROL)

- `tests/sn/sweep/core/test_affine_carve_baseline.py` + `tests/sn/solve/test_affine_
  carve_bit_identity.py` — `(L+C).apply(ψ)` / `transport_sweep` on het σ_t + non-flat
  ψ, **NO ray legs** (the ray is zeroed / (A,A)-only action). `kind="direct"`, the
  `DriftWarning` tripwire fires on ANY ULP move.
- **STAYS bit-exact on slab+sphere+cyl.** This is the load-bearing 4e CONTROL: it
  isolates "ray-leg extraction only" from "accidental (A,A) perturbation." The
  extraction touches only the ray legs; A_AA stays fused/unchanged in the walk, so the
  no-leg (A,A) output is byte-identical **on the sphere too**. **If 4e reds this on the
  sphere, the walk-slot rewrite folded arithmetic into the (A,A) path** — that is the
  bug, not an accepted drift. Run it `-W error::DriftWarning` as the e1 gate (§ Ask 5).

### D. The coherence swap-law home (`test_inverse_adjoint_coherence.py`)

| gate | :line | 4e fate | why |
|---|---|---|---|
| fused **G1/G2/G4** (slab, cyl_product) | 195-314 | **STAYS** (CONTROL) | on `_loss(sn)` for SEEDLESS meshes; 4e touches only the carrying-sphere walk. Bit-exact. |
| `test_carrying_fused_wrap_refuses_and_the_coupled_sibling_carries` | 479 | **STAYS** | the predicate contract (fused refuses / coupled carries); 4e does NOT flip `is_adjointable`. |
| **G1c** `test_coupled_gate1_forward_matvec_g_reciprocity` | 502 | **STAYS** | `⟨Mψ,x⟩_G=⟨ψ,b⟩_G < 1e-10` is a MATHEMATICAL identity (x solves M*x=b); holds ~1e-16 regardless of M's reduction tree. |
| **G2c** `test_coupled_gate2_swap_law_bit_identical` | 521 | **STAYS bit-identical** | `array_equal(M.H.inverse().apply(b), M.inverse().H.apply(b))` — BOTH spellings run the SAME `M.inverse().H` graph. Object identity (L19), bit-identical BY CONSTRUCTION even as M's internals change at 4e. |
| **G4c** `test_coupled_gate4_reverse_scan_executes` | 544 | **STAYS** (RE-POINT-safe) | Mode-11 wrap on `CoupledInvertibleOperator.solve_transpose` (the OUTER surface) — robust to the internal `to_unified` retirement. |
| **A2a fwd** `test_a2a_grid_forward_h_reciprocity` | 624 | **STAYS** (LIVE catcher) | `⟨A·ψ,x⟩_G=⟨ψ,A.H·x⟩_G < 1e-12` on the loss grid; identity holds ~1e-16. The M-ADJ-metric tooth (strip the CoupledSpace metric) still reds O(1). |
| **A2a inverse** `test_a2a_grid_swap_law_inverse_arm` | 665 | **STAYS xfail** | `xfail(strict=False)` until step-5 flips `grid.is_invertible`; 4e does NOT cross that boundary. |

### E. The driver / builder gates (psi_half `TestCoupledBuilder`, `TestWithinGroupSystem`, `TestWithinGroupSystemAnchors`)

| gate | :line | 4e fate | why |
|---|---|---|---|
| **G-c2.3** `test_grid_equals_m_minus_n` | 2998 | **STAYS — tightens** | grid(explicit blocks) ≡ M−N(welded); at 4e M's legs become extracted, so the two paths CONVERGE. Per-row bars (bulk rtol=1e-11) absorb the drift. |
| centrepiece teeth / M2 layout | 3047/3080 | **STAYS** | structural (off-diagonal swap unconstructable; layout coextensive). |
| **G-c2.5** `test_forward_block_adjoint_reciprocity` | 3134 | **STAYS** (LIVE catcher) | block-`.H` reciprocity < 1e-12 + M-ADJ-metric tooth; identity, holds ~1e-16. |
| **G-d1.4** `test_g_d1_4_m_leg_plumbing_consistency` | 3388 | **RE-POINT** | zero-leg ≡ no-leg + round-trip. ⚠ exercises the `to_unified`/`from_unified` round-trip IN M — when the bridge retires (§ Ask 4) the round-trip becomes the split-native pass-through; RE-POINT the accessor, value HOLDS. |
| **G-d1.5** `test_g_d1_5_gain_grid_matches_the_pieces_with_teeth` | 3455 | **STAYS** | N is the gain grid (S+B_a, Emission, B_b); 4e does not touch N. |
| **G-d3.1** `test_g_d3_1_flat_flux_equilibrium_on_the_carrying_sphere` (E4) | 3701 | **STAYS — the redistribution catcher** | φ=Q/Σ_t END-TO-END, rtol 1e-10, both drivers. A PHYSICS invariant a principled-equiv extraction preserves; a dropped/sign-flipped ray fold moves φ off Q/Σ_t O(1). **KEEP as the primary 4e correctness anchor** (flux-shape layer). |
| **G-d2.3** honest-DOF / ERR-053 restart | 3616 | **STAYS/RE-POINT** | the coupled `to_flat` count. If 4e shrinks the ray member's flat size (unified→split-native) the honest count could shift by the seed/corner layout delta — assert `restart == coupled.to_flat().size` recomputed from the NEW space (not a frozen integer). |

### F. The test_282 cold-residual acceptance suite (`tests/sn/sweep/curvilinear/test_282_direct_seed_fixed_point.py`) — THE L18 acceptance number

| gate | :line | 4e fate | why |
|---|---|---|---|
| **C(i)** `test_ci_cold_residual_is_machine_zero` | 170 | **STAYS — THE acceptance gate** | `‖A·solve(b)−b‖_∞ < 1e-11` slab/sphere/cyl. Certifies solve≡apply CONSISTENCY (not coefficients) — invariant to the extraction's reduction-tree drift. **The single headline number that certifies the un-weave didn't break the direct solve** (L18: 5.18e5 pre-fix → <1e-11). |
| **C(ii)** `test_cii_sphere_solve_is_seed_insensitive_bitwise` | 188 | **STAYS bit-identical** | `array_equal` on two solves with different seeds — compares two runs of the SAME post-4e code (structural seed-insensitivity), not a frozen value. |
| **C(ii)-probe6 / C(iv)** | 209/260 | **STAYS** | cold solve recovers preimage / pure-absorber exact — physics, principled-equiv. |
| **C(iii)** finite+positive | 233 | **STAYS** | physicality (NaN/neg → RED); coarse LS S4 end-to-end. |
| **Mode-11** `test_mode11_direct_solver_executes_on_sphere_not_slab` | 291 | **RE-POINT** — the CLOSEST precedent for the § Ask 2 sentinel | today it spies the direct-solver entry; post-4e RE-POINT the wrap onto the EXTRACTED engine (`RadialCharacteristicOperator.solve`). Its "sphere-not-slab" asymmetry IS the single-source proof pattern 4e generalizes. |
| **Mode-10** `test_mode10_seed_source_activation_q_half_moves_the_sphere_solve` | 432 | **STAYS — re-verify activation** | proves the seed block is CONSUMED (a zero q½ leaves the solve unmoved; a nonzero one moves it). After the extraction re-verify the injected q½ still moves the solve (an extraction that dropped the feed nulls this — silent coverage loss, §0.6). |
| **Mode-12** `test_mode12_g_reciprocity_catches_a_seed_row_flip` | 339 | **⚠⚠ RE-AIM — THE HEADLINE FINDING** | see below. |

**FINDING M1 (the headline) — the Mode-12 ERR-067 closure goes silently toothless.**
`test_mode12_g_reciprocity_catches_a_seed_row_flip` (:339) mutates
`_OneDimScanWalk._seed_rows_forward` (:413/:419 `monkeypatch.setattr(_OneDimScanWalk,
"_seed_rows_forward", _flipped)`) as its TOOTH — the seed-row forward flip that the
SPD `G_sd=V_cell` reciprocity catches (the ERR-067 Mode-12 closure, L18). **`_seed_
rows_forward` is EXACTLY a welded leg 4e DELETES** (retirement entry, the A_AB forward
seed rows). Post-4e the `monkeypatch.setattr` either raises `AttributeError` (symbol
gone) OR — worse — the walk no longer calls `_seed_rows_forward`, so the flip is a
NO-OP and the gate goes **GREEN-but-vacuous** (Mode-11: the tooth's target is off the
call graph). **RE-AIM the tooth onto the successor** — flip the seed-row inside the
extracted `RadialCharacteristicSeeding.apply` (or the shared kernel the walk now routes
through), and confirm the SAME `G_sd=V_cell` reciprocity reds `>1e-6` (the 0.107
measured defect). This is the "a tooth whose target symbol dies must be re-aimed" case;
without it the ERR-067 closure (the highest-value adjoint catcher in the campaign) is
lost. **The `test_282` file will not even IMPORT past 4e** (it patches a deleted
symbol) — this is a compile-time break, land the re-aim IN the 4e commit.

### G. STAYS UNCHANGED (do NOT touch)

- `tests/numerics/test_coupled_operator.py` (M1–M5, `TestBlockMatvec`/`TestAssembleProbe`/
  `TestBlockAdjoint`) — semantics-agnostic SYNTHETIC toy blocks, ZERO SN-walk reference.
  4e is SN-only. Grep-confirmed no ψ½/walk refs.
- `tests/diffusion`, `tests/cp` — 2-block base algebra, the eviction's bit-id oracle;
  4e touches only the ray. Assert ZERO re-points (a diffusion red = 4e broke the base).
- The Phase-C L0 gates (`tests/sn/sweep/core/test_phase_c_gates.py`:355/829 per-ordinate
  flat-flux residual `catches("ERR-026")`; :413 reciprocity) — L0 physics invariants
  (per-ordinate residual = 0), reduction-tree-independent; STAY green, but **re-run as
  the redistribution re-verification** (they are the ERR-026 curvilinear catchers).

---

## Ask 2 — THE SINGLE-SOURCE PROOF (the Mode-11 discipline)

The extraction's claim is "production routes through the block engines, the welded
copies are gone." A green value gate measures the UNCHANGED sibling (Mode 11) — the
claim needs a committed catcher on BOTH halves.

**S1 — the engine-execution sentinel (the newly-routed line MUST fire).** On a REAL
production `solve_sn_fixed_source`/`solve_sn` carrying-sphere run (BOTH
`inner_solver="source_iteration"` and `"krylov"`), an in-process autouse/`monkeypatch`
**WRAP counter** on the entry point the walk now routes to returns `count > 0`:

| routed line (4e implementer's choice) | wrap target | assert |
|---|---|---|
| (A) route through the engine OBJECTS | `RadialCharacteristicOperator.solve`/`.apply` + `RadialCharacteristicSeeding.apply`/`.apply_transpose` | each `> 0` on a real carrying solve |
| (B) route through a shared FREE-FUNCTION kernel | that kernel (as `radial_characteristic_forward_residual` is today) | `> 0`, AND both the walk AND the engine hit it (the existing `test_walk_forward_routes_through_the_shared_kernel` is the template) |

This is the L19 gold-standard (wrap the internal call in-process, not a file-write).
Tooth: an inlined walk that recomputes the seed leg leaves the counter at 0 → RED.
**Precedent:** `test_282::test_mode11_direct_solver_executes_on_sphere_not_slab`
(:291) is the exact pattern — RE-POINT + generalize it to all four legs.

**S2 — the welded-symbols-are-GONE tripwire (the reverse).** An AST/grep-tier absence
test (the `test_one_dim_loop_walk.py` spy+AST-tripwire precedent): assert the walk
module no longer DEFINES `_seed_rows_forward`/`_seed_rows_transpose` (the welded A_AB
legs) — `not hasattr(_OneDimScanWalk, "_seed_rows_forward")`, and grep-tier
`"_seed_rows_forward" not in loss_representation source`. This is what forces FINDING M1
to be discharged: the deleted symbol is the proof the extraction happened, and the
Mode-12 tooth re-aim is unavoidable once it's gone.

**H1 — the fused-A_AA+A_AB claim boundary (the honest scope).** `precompute_psi_state`/
`cell_contribution` compute A_AA + A_AB **summed in one pass**. Two honest 4e designs,
each with its pin:

- **H1-narrow (extract A_BB only; keep A_AB fused with A_AA).** The forward matvec
  (A,A)+(A,B) stays BIT-IDENTICAL (only A_BB solve/apply extract → drift lives in the
  ray self-block). More bit-identity survives; the sphere walk_matvec `fwd_bulk` may
  even STAY bit-exact if the A_AB feed is untouched (only `adj`/solve drift). Pin:
  `test_apply_matches_the_in_sweep_seed_injection` (:2525) stays the A_AB↔walk
  single-source proof; the sphere re-capture (§B) shrinks to the solve/transpose blocks.
- **H1-wide (extract A_AB too, as a separate `…Seeding.apply` + sum).** `A_AA·ψ +
  A_AB·seed` reassociates vs the one-pass fused sum → the sphere `fwd_bulk` drifts
  ~5.5e-16 → RE-CAPTURE. The claim boundary is "A_AB is isolable (the zeroed-psi_view
  proof) AND now separately routed"; the sentinel S1 must fire on `…Seeding.apply` in
  the forward walk.

**State which design 4e picks (AskUserQuestion at 4e start) and name the pin.** The
memo's default recommendation: **H1-narrow** (extract the genuinely-separable A_BB ray
self-block; leave A_AB fused with A_AA — it is already single-sourced via the shared
kernel, and forcing a second-pass sum buys reassociation drift for no architectural
gain). Under H1-narrow the sphere forward matvec can STAY bit-exact and only the
solve/transpose re-capture — maximizing surviving bit-identity.

---

## Ask 3 — THE EQUIVALENCE ARCHITECTURE FOR THE CARVE

**Control / experiment split.**

- **CONTROL (bit-exact, `array_equal`/`nulp=1`, both drivers):** non-carrying `slab` +
  non-carrying `cyl` — no ray, the walk is untouched. Gates: `walk_matvec_slab_2g`/
  `cyl_2g` (nulp=1), the fused G1/G2/G4 (slab/cyl), the affine-carve baseline (all
  geoms, incl. the sphere's (A,A)-only path — the ray-less control ON the carrying mesh).
  **A slab/cyl red = the extraction leaked into the shared kernel** (should be impossible;
  it is the Mode-9-style tripwire).
- **EXPERIMENT (carrying sphere, per-surface bars):**

| surface | deviation mechanism | bar | structurally-independent anchor |
|---|---|---|---|
| **matvec fwd** `M.apply` | A_AA+A_AB fold reassociation (H1-wide) / NONE (H1-narrow) | 5.5e-16 rtol / 0-ULP | `_dense(M.apply)` self-probe is TAUTOLOGICAL — anchor to **G-d3.1 φ=Q/Σ_t** + the reciprocity `⟨Mψ,x⟩=⟨ψ,b⟩` (G1c), NOT to old-vs-new ULP |
| **matvec transpose** `M.apply_transpose` | seed-row pullback reassociation | 5.5e-16 rtol | **A2a forward reciprocity** (`⟨Aψ,x⟩_G=⟨ψ,A.H·x⟩_G`) — the transpose's independent check |
| **solve fwd** `M.solve` | forward-substitution reassociation (2/Δr vs Δr·σ+2, L7) | 5.5e-16 rtol | **row-6 dense-LU oracle** `test_extract_to_dense…` (:397) — LAPACK LU ⟂ WDD sweep; THE keystone |
| **solve transpose** `M.solve_transpose` | reverse-scan reassociation | 5.5e-16 rtol | **G3 round-trip** `A*(A*)⁻¹b=b` + the dense `solve(M.T,b)` cross-check |

The row-6 oracle (`‖extracted − LU(_dense)‖ < 1e-11`, ~5.5e-16) is the **durable
principled-equiv floor** the plan names; it is `.solve`-based and LAPACK-independent →
the keystone for all four surfaces' solve legs. It STAYS-GREEN (§A) and becomes literal
at 4e.

**Cold-vs-converged coverage (ERR-053/#282 lesson: converged-only gates mask seed
defects).** The `test_282` C(i) COLD residual `‖A·solve(b)−b‖<1e-11` is the cold
catcher — it certifies the extracted solve is a genuine inverse from a COLD start
(no iterate-threading hides a seed defect). PAIR with the converged G-d3.1/G-d1.6
end-to-end. **Do NOT rely on the converged anchors alone** — the un-weave touches the
seed march, exactly the ERR-053 family's blind spot.

**Eigenvalue-level smoke.** The B.2d d2 precedent: sphere-Krylov `Δ1.8e-10` was the
anticipated honest-DOF arm; sphere-SI was BIT-IDENTICAL. For 4e assert `k(SI) ≡
k(Krylov)` cross-driver at `Δ ≤ 1e-10` on the carrying sphere (the honest mechanism =
the extracted-solve reassociation + Krylov-subspace change), NAMED — NOT "bit-identical
k" (the sphere solve genuinely reassociates). k_inf vs the derivations closed form
(`test_kinf_homogeneous`, the G-d3.2 cross-reference) stays the eigenvalue-layer anchor
at rtol 1e-10, **Mode-12-paired** (k never credited against a shape-class mutation — the
pack leg-swap / ray-relabel is caught OBJECT-level by G-d1.4 + the het FLUX FIELD by
G-d1.6/G-d3.1, never by keff alone).

---

## Ask 4 — THE UNIFIED-RETIREMENT COMPLETENESS GATES

4e retires: the `to_unified`/`from_unified` bridge, the unified `RadialCharacteristic
Field` leaf (`fields/_bases.py:1105`), the unified `RadialCharacteristicSpace`
(`numerics/spaces/radial_characteristic_space.py:313`), and the F2 space-identity
aliasing. This is the biggest retirement in the campaign (~40 files reference
`radial_characteristic_space`).

**R1 — the bridge retirement + the split-fidelity replacement.** `RadialCharacteristic
Composite.to_unified`/`from_unified` (`radial_characteristic_composite.py:258/:226`)
were the "B.1d licence" — the role-preserving exact re-label, pinned by a
`to_unified∘from_unified == id` gate (the DEMOTION licence). Post-4e the composite IS
the native representation → the round-trip has no referent. **The split-fidelity pin
that replaces it:** the B.1b arange-fingerprint / member-slice gate on the
`radial_characteristic_composite_space` (interior⊕corner offsets + `np.arange`
fingerprint of the member slices) — assert the composite's `to_flat`/`from_flat`
round-trips bit-exactly WITHOUT the unified detour. The `_dense_ray` probe helper
(psi_half :464, "probing the composite basis IN THE UNIFIED LAYOUT — bit-comparable
with the pre-B.2b `_dense_seed` matrices") RE-POINTS to probe the composite's NATIVE
flat basis (drop the `.to_unified().values` currency; use `.to_flat()`). Every helper
using `.to_unified()`/`from_unified` as the comparison currency (psi_half `_ray_
composite`:455, `_dense_ray`:464, `_seed_leaf`:431; the coherence `_ray_composite`:461;
the walk_matvec generator's `radial_characteristic_flux=` leg) RE-POINTS.

**R2 — the F2 aliasing dissolution (deletion makes it structural — NO new gate).** B.2c
minted the F2 runtime-apply proof precisely BECAUSE the unified and composite spaces
COLLIDE on `(name, shape)=("radial_characteristic",(ni+nb,))` ⟹ `FunctionSpace.__eq__`
is Mode-12-blind to a still-unified-typed block (`build_coupled_system` G-c2.1 runtime
catcher; the memo called it "transitional aliasing, dissolves at 4e"). **When the
unified space is DELETED, the collision is structurally impossible** — there is only
ONE `("radial_characteristic", …)` space (the composite). So: (a) the F2 runtime-apply
proof (G-c2.1's collision arm) RE-POINTS to assert the composite space is now the ONLY
one (`sn.radial_characteristic_composite_space` is the sole carrier; `sn.radial_
characteristic_space` attribute is GONE — `hasattr` False), and (b) any gate that used
IDENTITY-not-`==` to dodge the collision (e.g. `test_b2c_member_composite_block_
boundary` :1334 "== cannot see a unified-typed block") can RELAX to `==` post-deletion
— but KEEP the identity check as a stronger-form regression guard (cheap, no downside).
No NEW gate needed; the deletion is the proof — assert the ABSENCE.

**R3 — the 3-search retirement audit (the implementer's obligation; enumerate).** Per
the retirement-means-test-migration rule, run all THREE searches and RE-POINT vs RETIRE
each. Known consumers of the unified leaf/space that PIN the layout (must re-point vs
retire):
- `tests/sn/mesh/test_radial_characteristic_carrier.py`, `…_slot_coordination.py`,
  `…_split_spaces.py` — the carrier algebra on `sn.radial_characteristic_space`;
  largely the unified-layout probes → RE-POINT to the composite native layout or RETIRE
  the unified-only rows (the composite has its own carrier suite).
- `tests/sn/operators/test_ray_operator.py`, `…test_radial_characteristic_metric.py`,
  `…test_loss_transpose_solve.py`, `…test_invertible_operator.py` — `.radial_
  characteristic_space` probes → RE-POINT.
- `tests/sn/sweep/test_assembly_mode.py`, `…curvilinear/test_282…`, `…curvilinear/
  test_psi_half_angle_seed.py`, `tests/numerics/test_face_layout_typed_key.py`/
  `test_face_streaming_normal.py` — mixed; audit each.
- `tests/sn/regression/_generate_walk_baselines.py` — RE-POINT the `radial_
  characteristic_flux=` leg (§B).
- **The DOC blast radius (search 2, the silent one):** `:class:`/`:func:` cross-refs to
  `RadialCharacteristicSpace`/`RadialCharacteristicField`/`to_unified` in `docs/theory/`
  render as plain text with NO `-W` warning — grep `docs/` for the retired symbols and
  re-point (the B.2d NIT-1 precedent: 13 stale `_within_group_triple` doc refs rendered
  silently). Nexus `impact`/`callers` on `RadialCharacteristicSpace` for the graph-reached
  leaves the text-grep misses (property-reached, class-name-bypass).

**R4 — the `CoupledInvertibleOperator` bridge-body rewrite (the M surfaces).** M's four
surfaces (`coupled_system.py:370/390/418/436`) currently `to_unified` the System-B
member, call the walk, `from_unified` back. Post-4e they pass the composite member
DIRECTLY to the split-native walk. **G-d1.4's round-trip RE-POINTS** (§E); G4c's Mode-11
sentinel on `solve_transpose` (the outer surface) is robust (STAYS). Correctness gate:
the M surfaces' OUTPUT is principled-equiv (5.5e-16) vs the pre-4e captured baseline
(`--capture-baseline` before 4e) — the same sphere-re-capture bar as §B.

---

## Ask 5 — SEQUENCING (verification-first)

**Recommendation: e1 (walk-slot rewrite, bit-id-ABLE) FIRST, then e2 (extraction,
principled-equiv).** This keeps the frozen 0-ULP baselines valid LONGEST and localizes
the bit-identity break to e2's sphere solve/transpose rows only.

- **e1 — WALK-SLOT REWRITE (split-native marching + unified leaf/space/bridge
  retirement). BIT-ID-ABLE.** The rewrite changes only the STORAGE/marching container
  (unified buffer → composite native), preserving the recurrence order and the (A,A)
  reduction tree. The ray-slot DOFs are separate from the bulk, so removing them from
  the System-A march does NOT touch the (A,A) arithmetic; the ray recurrence is
  order-fixed by the sweep direction, not by storage. **Gate e1 (the tripwire that
  PROVES it was a pure representation change):** run EVERY frozen baseline under
  `-W error::DriftWarning` — `walk_matvec_{slab,sphere,cyl}_2g` at nulp=1, affine-carve
  bit-exact ALL geoms, `test_282` C(i)/C(ii) bit-exact, the coherence G2c
  `array_equal`. **If e1 drifts ANY of these (especially the sphere walk_matvec or the
  affine-carve sphere), the rewrite folded arithmetic — that is the bug, not an accepted
  drift.** e1 is likely ATOMIC (retiring the unified space breaks `cells_view`/`corner_
  view` which the walk marches — the same kwarg-coupling that forced the d2 eviction
  atomic). ⚠ **Caveat:** if the composite's internal member order forces a reassociation
  the recurrence cannot avoid (a sum-over-members regroups), e1 is NOT bit-id → the
  sphere re-captures AT e1 (anchored to row-6). The `-W error::DriftWarning` run tells
  you which — do NOT assume; measure.
- **e2 — EXTRACTION (route the split-native legs through the engines + delete the
  welded copies). PRINCIPLED-EQUIV.** Gate e2: the sphere solve/transpose RE-CAPTURE
  (anchored to the row-6 dense-LU oracle at 5.5e-16); slab/cyl STAY bit-exact (the
  control — an e2 red here = the extraction leaked into the shared kernel); the S1
  engine-execution sentinel fires; the S2 welded-symbols-GONE tripwire passes; the
  FINDING-M1 Mode-12 tooth is RE-AIMED onto `…Seeding.apply`; `test_282` C(i)
  cold-residual STAYS <1e-11.

**Why NOT e2-first (extract behind the bridge, then rewrite):** extracting while the
unified buffer still exists forces the engines to round-trip through `to_unified` — the
extraction reassociation entangles with the bridge, making the drift harder to localize
and requiring the sphere to re-capture TWICE (once at e2 behind the bridge, again at the
rewrite). e1-first localizes the reassociation to e2's clean split-native path.

**Why NOT one atomic landing:** e1 (rewrite) is genuinely bit-id-able and e2 (extract)
is genuinely principled-equiv — landing them atomically CONFLATES the two failure modes
(a representation bug and a reassociation drift both show as sphere drift), destroying
bug-localization. The d2 atomicity was FORCED by kwarg-coupling (carrier+space+protocol);
here e1↔e2 are separable (extract is additive-then-delete on top of the rewritten walk).
Keep them separate unless the implementer finds a hard coupling — then flag it.

---

## Mutation ledger (monkeypatch-only; NEVER `git checkout` an uncommitted file — L4)

| gate | tooth (in-process) | target `file:line` (⚠ best-estimate) | expected RED |
|---|---|---|---|
| **M1 RE-AIM** Mode-12 seed-flip | flip the seed-row inside `RadialCharacteristicSeeding.apply` (the extracted successor of `_seed_rows_forward`) | `radial_characteristic.py` `…Seeding.apply` (NOT the deleted `loss_representation:2853`) | `G_sd=V_cell` reciprocity `>1e-6` (~0.107) + control holds <1e-12 |
| S1 sentinel | inline the walk's ray leg (bypass the engine) | NEW 4e walk route | engine wrap counter = 0 |
| S2 tripwire | leave `_seed_rows_forward` defined | walk module | `hasattr(_OneDimScanWalk,"_seed_rows_forward")` True → RED |
| row-6 oracle | naive dense `M⁻¹` ignoring the inflow/seed rows | (the gate's own reference) | `‖extracted − LU‖ > 1e-11` |
| G-d3.1 (E4) | drop / sign-flip the extracted ray fold | NEW 4e route | φ moves off `Q/Σ_t` O(1) |
| #284 triangular | ray reads a bulk moment (back-edge) in the extracted engine | `RadialCharacteristicOperator` | `A_sb > 1e-12` O(1) |
| affine-carve (e1) | fold the (A,A) march into the ray rewrite | e1 walk-slot rewrite | sphere bit-exact RED under `-W error::DriftWarning` |
| Mode-10 activation | zero the injected q½ | (the gate's own input) | the solve MUST still move (else silent coverage loss) |

## Refutations (fire before the ink dries)

1. **The Mode-12 ERR-067 closure tooth targets a DELETED symbol** (FINDING M1) —
   `test_282::test_mode12…seed_row_flip` patches `_seed_rows_forward`; 4e deletes it →
   `AttributeError` (compile break) or silent no-op (Mode-11 vacuous-green). RE-AIM onto
   `RadialCharacteristicSeeding.apply` IN the 4e commit. Highest-signal survivors finding.
2. **The row-6 oracle is the durable floor, NOT a re-baseline victim** — it was designed
   as the anchor every EXTRACT step pins against; it STAYS-GREEN at 5.5e-16 and anchors
   the sphere walk_matvec re-capture (do NOT re-baseline old-vs-new ULP — L18).
3. **G2c bit-identity SURVIVES the extraction** — it compares two spellings of the SAME
   `M.inverse().H` graph (object identity, L19), bit-identical BY CONSTRUCTION even as M
   drifts. Do NOT loosen it to rtol.
4. **The affine-carve sphere row is the (A,A)-untouched CONTROL** — it must STAY bit-exact
   ON THE SPHERE (no ray leg); a sphere red there = the walk-slot rewrite folded (A,A)
   arithmetic, the e1 bug tripwire.
5. **cyl is the CONTROL, not a re-capture** (FINDING C1) — cyl is non-carrying (#229
   `τ_raw=0` clamp); the generator's "cylinder re-captured" prose predates that fact.
   Verify `radial_characteristic_space is None` for the generator's cyl; expect bit-exact.
6. **4e does NOT flip `grid.is_invertible`** — the A2a inverse arm STAYS `xfail`; the
   block-solve spelling is step 5 (#41). A 4e that greens `test_a2a_grid_swap_law_inverse_
   arm` is out of scope.
7. **cold-residual, not converged-only** — the un-weave touches the seed march (the
   ERR-053/#282 blind spot); `test_282` C(i) COLD `‖A·solve(b)−b‖<1e-11` is the acceptance
   number, PAIRED with the converged G-d3.1 (converged-only masks seed defects).
8. **the unified-space deletion makes F2 STRUCTURAL** — assert the ABSENCE
   (`hasattr(sn,"radial_characteristic_space")` False; the composite is the sole carrier),
   NOT a new runtime collision gate.

## Result contract

4e lands in TWO sub-commits: **e1** (walk-slot rewrite split-native + unified leaf/space/
bridge retirement, BIT-ID-ABLE — every frozen baseline STAYS at nulp=1 under
`-W error::DriftWarning`) then **e2** (extraction, principled-equiv — sphere solve/
transpose RE-CAPTURE anchored to the row-6 dense-LU oracle at 5.5e-16; slab/cyl bit-exact
CONTROL). NEW gates: **S1** engine-execution sentinel (both drivers, the L19 wrap counter
on the newly-routed line) + **S2** welded-symbols-GONE tripwire (AST/grep absence of
`_seed_rows_forward`/`_seed_rows_transpose`). RE-AIM: **FINDING-M1** the Mode-12 ERR-067
seed-flip tooth (`test_282:339`) onto `RadialCharacteristicSeeding.apply` (IN-commit —
the file won't import past 4e). RE-CAPTURE: `walk_matvec_sphere_2g.npz` (anchored to
`_dense(lc.apply)`, NOT old-vs-new ULP). RE-POINT: G-d1.4 M round-trip + every
`.to_unified()`-currency probe helper (psi_half/coherence/walk-generator) + the doc
`:class:`/`:func:` refs to the retired symbols. STAYS (do NOT touch): the numerics
synthetic M1–M5, diffusion/cp walls, the Phase-C L0 ERR-026 catchers, the fused-slab/cyl
coherence controls, the row-6 oracle, G-c2.3/G-c2.5/G-d1.5/G-d3.1, test_282 C(i)-C(iv).
**Load-bearing deliverables:** the **row-6 dense-LU oracle** (the 5.5e-16 keystone every
surface's solve leg anchors to) + **test_282 C(i) cold-residual <1e-11** (the L18
acceptance number) + **G-d3.1 φ=Q/Σ_t** (the redistribution catcher) + **the FINDING-M1
Mode-12 re-aim** (or the ERR-067 closure goes silently toothless) + **the affine-carve
sphere (A,A)-control** (the e1 bit-exact tripwire). End-to-end acceptance per sub-commit:
full `tests/sn -m "not slow"` + `tests/numerics` + `tests/diffusion` + `tests/cp` +
ratchet `transport:1` + sphinx -W; the wall STAYS 3733/0 · 37 xfailed through e1 (bit-id),
and after e2 only the sphere walk_matvec `.npz` moved (re-captured) — every count else holds.
