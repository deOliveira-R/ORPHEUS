---
name: radial-characteristic-forward-extract-verification
description: The 4b gate spec for the Coupled Block Operator campaign — completing A_BB (RadialCharacteristicOperator) forward apply/apply_transpose via a SHARED-KERNEL EXTRACTION (relocate _seed_residual_march + extract radial_characteristic_forward_residual to psi_half_angle_seed.py; the walk's _seed_rows_forward/_transpose become thin wrappers; A_BB.apply wraps the SAME fn) + flip is_invertible/is_adjointable True + inverse() involution. Covers: the walk-bit-id frozen baseline + Mode-11 anti-twin routing sentinel; the involution round-trip is PRINCIPLED-EQUIV ~3 ULP NOT 0-ULP (refutes "bit-exact by construction" — the dr·(2/dr)≠2 reassociation; ONLY the outflow-corner closes bit-exact) AND identity only on the CONSISTENT subspace; the transpose seam (walk method = shared_fn + seed_cells_bar A_AB term, A_BB.apply_transpose = shared_fn ONLY, feed empty seed_cells_bar); Euclidean-NOT-metric adjoint (L19/Mode-12 is the composite's job); the predicate-flip contract is EMPTY of assertions (net-new positive pins + doc-drift only). Sibling of coupled_operator_step4 + aba_schur_fold.
metadata:
  type: reference
---

# Verifying a shared-kernel forward EXTRACTION — A_BB.apply (Coupled Block Operator Step 4b)

**When this recipe applies.** An operator's forward `apply` is currently
WOVEN into a monolithic walk (`(L+C).apply`) as a per-leg residual march; a
carve COMPLETES the operator by (a) relocating the inner engine + extracting
the multi-leg orchestration into ONE shared free function beside the SOLVE
engine, (b) making BOTH the walk row-emitter AND the operator's `apply` thin
wrappers of that one function (zero twin — user ruled extract-now, not
tracked-twin), (c) flipping `is_invertible`/`is_adjointable` True + wiring the
`inverse()` involution. Sibling of `coupled_operator_step4_verification.md`
(the ASSEMBLE/LIFT/EXTRACT parent) and `aba_schur_fold_unweld_verification.md`
(the un-weld Mode-11 recipe). This is a DIFFERENT decomposition axis than the
step-4 memo's 4a/4b/4c (WRAP/LIFT/EXTRACT) — here "4b" = the A_BB.apply carve.

**The object (worked instance).** `A_BB = RadialCharacteristicOperator`
(`orpheus/sn/operators/radial_characteristic.py`). Its `solve`/`solve_transpose`
(the Carlson march resolvent + Euclidean adjoint) already landed (step 1c). 4b
adds the FORWARD: `A_BB.apply(ψ) = (μ∂_r + σ_t)ψ` (the DD residual `m = σ·c +
(2/Δr)(c−f)`), the exact algebraic inverse of the march. Extraction map:
- `_seed_residual_march` (`loss_representation:170`, ONE caller `_seed_rows_forward`)
  → relocate to `psi_half_angle_seed.py` (beside `carlson_inward_sweep_from_source`).
- NEW `radial_characteristic_forward_residual(values, space, sigma, dr)` (the
  two-leg orchestration + corner rows, lifted verbatim from `_seed_rows_forward:2761-2789`)
  → `psi_half_angle_seed.py`.
- NEW `radial_characteristic_forward_residual_transpose(values, space, sigma, dr)`
  (the PURE reverse straight-line program, lifted from `_seed_rows_transpose:2823-2843`).
- `_seed_rows_forward` → thin wrapper: `return ..._forward_residual(seed.values, seed.space, sigma, dr)`.
- `_seed_rows_transpose` → wrapper + coupling: `out = ..._transpose(chi_seed.values, space, sigma, dr)`
  THEN `for p in seed_cells_bar: cells_view(out,p,-1) += seed_cells_bar[p]` (2844-2846 — the A_AB term STAYS in the walk).
- `A_BB.apply` wraps `..._forward_residual`; `A_BB.apply_transpose` wraps `..._transpose` (PURE, no coupling).
- Flip `is_invertible=True`, `is_adjointable=True`; `inverse()` → involution
  (`inverse().apply ≡ solve`, `inverse().solve ≡ apply` — the `InverseOperator`/`InverseWrapMixin` leaf pattern, requires `is_invertible=True`).

Home: `tests/sn/operators/test_psi_half_coupling.py::TestA_BB_RadialBVP` (extend;
`foundation`). Carrying = sphere-GL S4 (`_sphere`/`_graded_sphere`, ≥2G); cyl/slab =
non-carrying CONTROL (constructor refuses — existing `:936`). Reuse helpers
`_ray_sigma`/`_ray_source`/`_ray_cotangent`/`_two_leg_reference`/`_install_engine_spy`/
`_euclid_adjoint_defect`. Canonical `-O`: sentinels/teeth = `np.testing.*`/`pytest.fail`/raise.

## The 4 load-bearing refutations (fire before the ink dries)

1. **"bit-exact by construction" is FALSE for the cell round-trip — it is
   principled-equiv ~3 ULP** (MEASURED 6.66e-16, nx=7). The forward residual
   uses `two_over_dr=2/Δr`; the inverse march divides by `denom=Δr·σ+2`;
   composing reassociates `Δr·(2/Δr)≠2` in FP. ONLY the outflow-corner defect
   closes bit-exact 0.0 (subtracting a same-ops-recomputed face from its stored
   copy). ⟹ gate 2 cells at `rtol=1e-11` (L7 SAFETY ≫ ~10 ULP), corner at
   `array_equal` 0. Try `array_equal` FIRST to CONFIRM it fails (documents the claim).
2. **`solve∘apply=id` is identity ONLY on the CONSISTENT subspace.** The +1
   outflow-corner is a FREE datum `apply` MEASURES (residual `f_outer−k_plus`)
   and `solve` OVERWRITES (`corner_out=f_outer`). Testing on an ARBITRARY ψ
   FALSELY REDS (corners differ, non-bug). Build ψ0=`solve(q0)` (consistent), OR
   test `apply∘solve` on a PHYSICAL source (corner+1=0, which `_ray_source`
   satisfies — the q½ fold never writes the +1 source corner). §0.6-adjacent:
   the convenient arbitrary-ψ config activates a term outside the identity claim.
3. **The transpose seam: the walk method is NOT a pure wrapper.**
   `_seed_rows_transpose = ..._transpose(...) + seed_cells_bar` (the A_AB M-M
   coupling term, `:2844-2846`). `A_BB.apply_transpose` wraps ONLY the pure
   `..._transpose`. Gate 3c/4 MUST feed `seed_cells_bar={}` so the walk method
   reduces to the pure A_BB transpose (aligns with the existing A_AB-isolation
   at `test:1556`, zero χ_seed). A gate that let the A_AB term leak would falsely
   diverge A_BB.apply_transpose from the walk.
4. **Gate 4 is EUCLIDEAN, NOT the V_cell metric.** `⟨apply(u),v⟩=⟨u,apply_transpose(v)⟩`
   is the plain flat dot. The `G_sd=V_cell` Hilbert reciprocity (`.H`) is realized
   ONCE at the COMPOSITE (L19/Mode-12/ERR-067) — the step-4-ASSEMBLE `.H` gate, NOT
   4b. Do NOT gate a metric reciprocity on the leaf A_BB here.

## The 6 gate classes

### Gate 1 — walk stays bit-identical + Mode-11 routing sentinel (regression/characterization)
- **Existing coverage — SUFFICIENT for VALUE, INSUFFICIENT for bit-identity.**
  `TestRegressionFloor::test_welded_sweep_is_exact_direct_inverse` (`:266`,
  `(L+C).solve(apply(ψ))≈ψ` @3.5e-16) + `test_extract_to_dense...` (`:286`,
  dense-LU of `(L+C).apply` @5.5e-16) + `tests/sn/sweep/**` end-to-end + the
  ERR-067/L18 Mode-12 catcher `test_282_direct_seed_fixed_point.py:367`
  (monkeypatches `_OneDimScanWalk._seed_rows_forward` — survives IFF the method
  stays as a wrapper). All are TOLERANCE gates → keep the walk GREEN at current
  tol (bit-identical relocation) but cannot PROVE 0-ULP unchange.
- **1a (NEW, RECOMMENDED) — frozen pre-carve bit-identity baseline** (L17/L19). BEFORE
  touching production, `assert_regression --capture-baseline` on `(L+C).apply(ψ).to_flat()`,
  `(L+C).apply_transpose(χ).to_flat()`, AND `_seed_rows_forward(σ_t, seed).values`
  (the isolated unit) on the fixed-seed sphere-GL S4 ≥2G het composite (non-flat, §0.6).
  Post-carve: `array_equal`. The honest byte-exact regression the tol-floor can't give.
- **1b (NEW, ESSENTIAL) — Mode-11 anti-twin routing sentinel.** The value floor
  stays GREEN if the walk still calls a LEFTOVER `_seed_residual_march` copy in
  `loss_representation` (L16 "only-new-reds→twin"/Mode-11 "gate-never-executes-the-
  rewired-path"). Spy the relocated `radial_characteristic_forward_residual`
  (patched in the WALK's namespace) → `(L+C).apply(seed-carrying ψ)` enters it
  count>0; same for `..._transpose` via `(L+C).apply_transpose`. Structural pair:
  grep-assert NO `_seed_residual_march` def survives in `loss_representation`
  (retire-as-you-go — it MOVED).

### Gate 2 — apply correctness via the involution (PRINCIPLED-EQUIV, refutations 1+2)
- **2a `apply∘solve≈id`:** q=`_ray_source` (physical, corner+1=0); `apply(solve(q)).values`
  vs `q.values` → cells+corner(−1) `assert_allclose(rtol=1e-11)`; corner(+1) `array_equal` 0.0.
  Graded sphere, ≥2G, ≥2 draws.
- **2b `solve∘apply≈id`:** ψ0=`solve(q0)` (consistent subspace, refutation 2);
  `solve(apply(ψ0)).values` vs `ψ0.values` `rtol=1e-11`. Graded sphere, ≥2G.
- **2c corner-closure bit-exact:** `apply(solve(q))` outflow +1 corner == 0.0 `array_equal`
  (the walk's precise "streamed−stored on a consistent state" claim; 0.0 IS achievable — L7).
- **Tooth** (monkeypatch the relocated fn, `-O`): face-chain sign flip `2·c−f → 2·c+f`
  OR magnitude `σ·c → 3σ·c` in `radial_characteristic_forward_residual` → round-trip
  defect O(1) ≫ 1e-11.

### Gate 3 — A_BB.apply ≡ walk forward (bit-id by construction) + Mode-11 anti-twin
- **3a bit-identity:** `A_BB.apply(seed).values == _loss(sn)._seed_rows_forward(σ_t, seed)`
  `array_equal` (SAME shared fn, σ=σ_t — feed the walk σ_t, NOT the σ-free 0 of `streaming_action`;
  A_BB is the FUSED L+C ray block). Sphere-GL S4 ≥2G.
- **3b Mode-11 CENTERPIECE:** spy `radial_characteristic_forward_residual`; BOTH
  `A_BB.apply(seed)` AND `_seed_rows_forward(σ_t, seed)` enter it (count>0 each, EXACT 1
  orchestration call/invocation — the fn loops levels internally). Green-but-count-0 = divergent copy.
- **3c transpose bit-id + anti-twin:** `A_BB.apply_transpose(v).values ==
  _loss(sn)._seed_rows_transpose(σ_t, v, seed_cells_bar={})` `array_equal` (refutation 3 —
  empty coupling ⟹ walk method == pure transpose); Mode-11 both enter `..._transpose`.

### Gate 4 — apply_transpose Euclidean-adjoint consistency (refutation 4)
- **4a:** `⟨A_BB.apply(u),v⟩ = ⟨u,A_BB.apply_transpose(v)⟩` < 1e-11 (plain flat dot; NEW
  helper `_euclid_apply_adjoint_defect` mirroring `_euclid_adjoint_defect:673` but on
  apply/apply_transpose). u,v random ray fields, ≥2 draws, ≥2G, graded (het σ). Reverse-mode
  of a linear SLP ⟹ exact transpose, defect ~1e-14.
- **Tooth:** monkeypatch relocated `..._transpose` — sign-flip incoming face-cotangent
  (`f_bar=-f_bar → +f_bar`, the `:2831/2838` sign; mirrors existing solve tooth `:744`) → defect O(1)>1e-3.

### Gate 5 — inverse() involution + predicate report
- **5a bit-id involution:** `A_BB.inverse().apply(q).values == A_BB.solve(q).values` `array_equal`;
  `A_BB.inverse().solve(ψ).values == A_BB.apply(ψ).values` `array_equal` (InverseOperator delegates
  apply→inner.solve, solve→inner.apply — same graph).
- **5b identity:** `A_BB.inverse().inverse() is A_BB` (InverseWrapMixin object identity, L15) or bit-id.
- **5c predicate:** `A_BB.is_invertible is True` AND `A_BB.is_adjointable is True` (net-new pins).
  ⚠ MUST land WITH `apply` — `InverseOperator.__init__` raises `NotInvertible` if `inner.is_invertible`
  False (the module's faithfulness contract: is_invertible=True ⟺ working solve+apply+inverse()).
- **5d Mode-11 delegation:** spy `A_BB.solve`; `A_BB.inverse().apply(q)` enters it count>0 (proves
  delegation, not a reciprocal twin). Tooth: an inverse() re-deriving the march → count 0.

### Gate 6 — the predicate-flip contract (test-architect flag-2 blast radius) — EMPTY of assertions
- **RESULT (live grep, 2026-07-08): the assertion blast radius is EMPTY.** NO test asserts
  `RadialCharacteristicOperator(...).is_invertible is False` / `is_adjointable is False`; NO
  `pytest.raises(NotImplementedError)` on `A_BB.apply` (the only such raise in the file is
  `RadialCharacteristicBoundaryOperator` white/albedo `:542`, unrelated); A_BB is NOT in
  `test_capability_survival.py`'s operator table (it enumerates SNBoundaryOperator/Incoming
  Source/Invertible/Streaming/Sweep/Fission/Scattering/N2N/Multiplication — NOT A_BB); A_BB is
  NOT composed in production (only docstring refs + the `loss_representation:4095` `.solve`
  comment) ⟹ no COMPOSITE `is_invertible/is_adjointable` AND-predicate flips. **The flip is
  additive-safe — no existing assertion breaks.**
- **(i) Doc-drift (Cardinal Rule 3 — prose, won't fail the run, MUST update):** module docstring
  `radial_characteristic.py:53-100` (the "Scope — resolvent action only" section + the
  is_invertible/is_adjointable-False explanation), `:175` (class docstring), `:260-281` (the
  `apply` NotImplementedError method + deferral prose); `test_psi_half_coupling.py:692` ("apply/
  inverse defer to step 4"). ⚠ DO NOT touch `:512-513`/`:584-589` — those are
  `RadialCharacteristicSeeding` (A_AB), whose `is_invertible` CORRECTLY stays False (rectangular).
- **(ii) RECOMMENDED net-new positive pins (L15 — the flip should propagate to a CONTRACT):**
  gate 5c is the minimal pin; the durable home is a positive row in `test_capability_survival.py`
  (else a future silent flip-back has no catcher).
- **Re-run at landing (the contract returns ZERO):**
  `grep -rn "is_invertible is False\|is_adjointable is False" tests/ | grep -i "radial\|ray\|A_BB"` → EMPTY;
  `grep -rn "RadialCharacteristicOperator" tests/ | grep "NotImplementedError\|raises"` → EMPTY;
  `grep -rn "RadialCharacteristicOperator" orpheus/ | grep -v radial_characteristic.py` → docstring-only.

## Claim-layer note (§1.5) — no NEW L1 apply-anchor needed
The forward `apply`'s VALUE is verified TRANSITIVELY: solve's L1 closed-form-ODE convergence
(`test_ray_operator.py::TestA_BB_RadialBVP`, already landed) ∘ the involution (gate 2, apply =
exact inverse of the convergence-verified solve) ∘ the walk equivalence (gate 3, apply ≡ the
production forward the whole sweep suite value-verifies). A standalone MMS/closed-form for the
residual would be redundant. The E4-style structural anchor (`Q/Σ_t` fixed-source, k_inf) lives
at the COMPOSITE (step-4-ASSEMBLE), not on the leaf residual.

## Mutation ledger (monkeypatch-only; NEVER `git checkout` an uncommitted file)
| gate | tooth | target `file:line` (post-relocation) | expected RED |
|---|---|---|---|
| 2 fwd | face sign `2c−f→2c+f` | `psi_half_angle_seed.py` NEW `..._forward_residual` (from `loss_representation:189`) | round-trip O(1) |
| 2 mag | `σ·c→3σ·c` | same | round-trip O(1) |
| 3b | (structural) walk inlines a copy | spy `..._forward_residual` in walk namespace | count 0 |
| 3c | (structural) transpose inlines a copy | spy `..._forward_residual_transpose` in walk namespace | count 0 |
| 4 | reverse face-cotangent `-f_bar→+f_bar` | `psi_half_angle_seed.py` NEW `..._transpose` (from `loss_representation:2831/2838`) | defect O(1)>1e-3 |
| 5d | inverse() re-derives, not delegates | spy `A_BB.solve` | count 0 during `inverse().apply` |
| 1b | leftover `_seed_residual_march` def survives | grep + spy | Mode-11 count 0 |

## Result contract
Extend `test_psi_half_coupling.py::TestA_BB_RadialBVP` (foundation) with gates 1b/2/3/4/5/6; add
the 1a frozen baseline (capture PRE-carve). Every tooth mutation-verified in-process under `-O`.
Acceptance: full `tests/sn -m "not slow"` 0 reds + `TestRegressionFloor` unchanged + the ERR-067
`test_282` Mode-12 catcher green + ratchet `transport:1` + sphinx -W. The `is_invertible`/
`is_adjointable=True` flip may ripple into the step-4-ASSEMBLE swap-law gates
(`test_inverse_adjoint_coherence.py`) — OUT of 4b scope, flag as a downstream note.
