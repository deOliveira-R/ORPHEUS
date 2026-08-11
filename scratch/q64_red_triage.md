# Q5.6.4 — red triage after the ω-partition carve

Slice: `tests/sn/sweep tests/sn/verification/mms/test_mms_ordering_blindness.py
tests/sn/mesh/test_cylindrical_quadrature_admission.py -m "not slow"`, serial,
`python -O`. **25 failed / 784 passed / 1 skipped / 31 xfailed** in 290.67 s.

## A. PRE-EXISTING — the #51 ledger, NOT caused by this carve (8 of the 25)

| test | ledger entry |
|---|---|
| `test_2d_octant_sweep_equivalence` × 5 (LS4/LS6) | the 5 octant rows, verified identical at `1f220c41` |
| `test_affine_carve_baseline::…_unmoved[CYL]` × 2 | sha256 affine-carve digests frozen at `63719a2` |
| `test_diamond::test_spherical_inward_bit_identical` | the spherical-inward diamond row |

⚠ verify each against the ledger before the merge gate; do not re-triage from
scratch.

## B. CAUSED BY THE CARVE — 17, and every one is a PREDICTED category

### B1. The absorber's TWIN in a test oracle (1) — blast-radius audit item (3)
* `test_cell_visit_c_stamp::test_multilevel_cylinder_production_stamp_matches_inline`
  ⟹ `tests/sn/sweep/core/_c_surrogate.py:137` still does
  `np.clip(tau_raw, 0.5, 1.0)` under a comment claiming the reference is
  UNCLAMPED. **Fix: delete the clip.** The "independent" oracle was replicating
  the absorber, so it is now a WRONG reference, not a broken SUT.

### B2. Gates whose THESIS was the clamp (3)
* `test_tau_producer_equivalence::test_cyl_tau_clamp_is_the_only_difference_from_reference[8,16]`
  — "the clamp is the only difference" is now vacuous. **Do NOT delete** (the
  audit: it is the cylinder's only vv-L11 producer gate). Re-pose to the
  sphere's 0-ULP shape: production τ ≡ an independently-written P2 expression.
* `test_tau_arc_wellposedness::test_the_folded_tau_is_bounded_with_the_reversal_identity[8,16,32,64]`
  (counted in B4 — the τ⊂[1/5,4/5] leg still passes; the BIT-EXACT reversal leg
  is what reds.)

### B3. ⭐ THE SCOPE QUESTION — theorems about rules production now REFUSES (9)
The new arc guard in `angular_cell_edges_per_level` refuses a level whose ω is
not monotone (the σ_y double cover). These rows deliberately drive the τ
producer with NON-ARC rules to establish theorems about it:

* `test_march_start_structure::test_the_tau_trichotomy_is_a_theorem_about_the_facts`
  × 7 — `product_node_aligned_{even_4,even_8,odd_5,odd_7}`,
  `product_staggered_full`, `level_symmetric_{4,8}`
* `test_tau_arc_wellposedness::test_shipped_families_pass_the_guard_including_the_closed_endpoints`
* `test_mms_ordering_blindness::test_alpha_and_tau_are_bit_identical_across_tie_breaks`

**The τ trichotomy (`on_edge_node ⟹ τ₀=0`, `degenerate ⟹ τ₀=1`) is a theorem
about the RETIRED η-midpoint partition.** Under the ω partition those rules have
no angular cells at all, so the trichotomy is not false — it is VACUOUS for
them. `[M]` R12a no longer depends on it: the seed-presence predicate is
`march_start_structure_per_level` (integer facts, Q5.4/T26); the trichotomy has
been a *gated consequence for documentation* since then.

⟹ Coherent with the 6.3 flip, which already made these rules inadmissible in
production. Retiring the theorems is arguably cleanup 6.3 owed. **But it is a
scope call, so it is the user's.**

⚠ **Is the guard even load-bearing?** `[M]` WITHOUT it, P3 still refuses
`product(2,8)` (`τ[2] = 22.43` outside `[0,1]`). So the guard changes the
MESSAGE, not accept/reject — *for that rule*. It is still the better design: a
non-arc rule could coincidentally land τ ∈ [0,1], and then the trichotomy rows
would assert WRONG VALUES rather than get a clear refusal. Keep it.

### B4. Predicted, already flagged in the plan (4)
* `test_the_folded_tau_is_bounded_with_the_reversal_identity[8,16,32,64]` —
  `[M]` τ⊂[1/5,4/5] STILL HOLDS (τ₀ = 0.2599 vs the old 0.2195); what reds is
  `assert_array_equal(τ + τ[::-1], 1)` — the reversal identity is now
  **1.11e-16, not bit-exact**. The ω partition trades bit-exactness for
  correctness (the chord partition's symmetry was exact *because* both end
  cells were stretched symmetrically). **Re-pose to 1-ULP**, and re-pose the
  `[1/5,4/5]` leg's numbers: the bound is now `[0.2599, 0.7401]` at n_φ=8.

### B5. Message-provenance (1)
* `test_a_mis_ordered_level_is_refused_as_outside_its_cell` — regex did not
  match: the mis-ordered level now trips the ARC guard before P3. Refusal still
  happens, one frame earlier, with a different (better-targeted) message.
  Re-point the regex, and keep a row that still exercises P3 itself.

## C. What is VERIFIED GREEN about the carve

`[M]` production τ on `folded_product(4, n_φ)` level 0:
* equals the P2 closed form `τ = ½ + ½cot(ω)tan(Δω/4)` to **1.67e-16**
* **ν-closure = 1.000000000000** exactly (the clamp gave 1.016389, τ≡½ 1.164784)
* P3 holds; reversal identity 1.11e-16
* **the SPHERE is untouched** — τ still equals the cumulative-weight reference
  (1.8e-15), and `[M]` 4/8 of its τ are below ½ at S₈ (8/16 at S₁₆), which is
  the literature's own point that `[½,1]` was never the admissible range.
