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

---

## D. ⚠ OPEN — an INDEPENDENT-REFERENCE cross-check went red after the re-baseline

`tests/sn/verification/analytical/test_phase_c_crosscheck.py::test_phase_e_trajectory_resolvent_flux_shape_crosscheck[cyl_2g_3reg_folded_4x8_dd_n40]`

`[M]` per-group max `|Δφ_norm|` = `[0.1268, 0.0134]`, overall `1.268e-01`
against `tol_per_cell = 1e-01`. The other 7 rows in that module PASS (the
k_eff literal update is correct — `1.2302082296342958 → 1.2310212585879858`).

**Why this is different from every other red in this file.** It is not a
snapshot, not a tolerance, not a re-pose. The test compares the SN cylinder
flux SHAPE against the **trajectory-resolvent** reference — a genuinely
different method (L1-class, structurally independent). So it is a claim about
correctness, not about frozen state.

**Mechanism:** the test reads the FROZEN SNAPSHOT's `scalar_flux`, so
re-baselining swapped the new flux in. `[M]` old→new flux change is
`4.38e-02` max-relative, which is the right order to push a ~10 % agreement
over a 10 % bar.

⚠ **DO NOT relax `tol_per_cell`.** If the carve moved the flux AWAY from an
independent reference, that is evidence to weigh against the partition change,
and the honest options are (a) accept with a stated, understood reason, (b)
re-open the partition question, or (c) show the resolvent reference itself
carries error at this level. Loosening the bar would erase the signal.

**IN FLIGHT:** probe `$CLAUDE_JOB_DIR/tmp/q64_phaseE_old_vs_new.py` computes
BOTH snapshots against the SAME resolvent reference, to establish direction and
magnitude. Until it reports, the sign of this finding is UNKNOWN — do not
assume the carve caused it (the old margin may already have been ~9.9 %) and do
not assume it did not.

⭐ Note the asymmetry that makes this worth chasing rather than filing: the
carve's justification is *structural* (P2/P3 satisfied, ν-closure exact), and
its accepted cost was an MMS-floor regression at n_φ≥16. An INDEPENDENT-METHOD
disagreement is a different and stronger class of evidence than an MMS floor,
because it cannot be dismissed as truncation-order bookkeeping.

### D-ANSWER `[M]` 2026-08-11 — the carve moved the flux AWAY from the reference, by 1.92×

Probe `$CLAUDE_JOB_DIR/tmp/q64_phaseE_old_vs_new.py`, BOTH snapshots against the
SAME resolvent reference (tolerance `1.000e-01`):

| τ convention | per-group max `|Δφ_norm|` | max |
|---|---|---|
| OLD chord + `[½,1]` absorber | `[0.06593, 0.01335]` | **`6.593e-02`** |
| NEW ω partition, no clamp | `[0.12676, 0.01341]` | **`1.268e-01`** |

⟹ the old margin was **6.6 %, comfortably inside** — not the ~9.9 % I
speculated. The carve **nearly doubled** the disagreement with an independent
method. Two independent instruments now agree the new τ is less accurate (MMS
floor ~1.8–2× worse at n_φ≥16; resolvent shape 1.92× worse), and the resolvent
one is NOT dismissible as truncation-order bookkeeping on a manufactured
fixture.

⛔⛔ **THE STRUCTURAL ARGUMENT FOR THE ω PARTITION RESTS ON A FALSE PREMISE, and
the literature memo said so in a sentence I under-weighted.**
`scratch/q64_tau_edge_convention_literature.md:943`:

> *"the **recursion-defined edge** — which is what α *is*, by definition, in
> both sources — and the geometric arc-half-angle edge, **which no source
> uses**."*

My case was: T3 gates α at the arc half-angle boundaries ⟹ τ must use the same
partition. **But T3's closed form is `α = −w_gl·κ·ξ(ω_{k−1/2})`, and the κ IS
the discrepancy between α's actual (recursion-defined) edge and the geometric
arc edge.** So α does NOT live at the geometric arc boundary; it lives at a
κ-scaled one. I unified τ onto the partition α *doesn't* use.

⭐ **This re-reads the whole session.** The chord partition may be *closer to
α's recursion-defined edges* than the arc one is — in which case it was
approximately RIGHT, and the `[½,1]` absorber was patching its endpoint
pathology (`τ→0` at an edge node) rather than compensating a wrong partition.
That is the opposite of this campaign's stated finding.

⚠ Caveat, stated so it is not over-read: the resolvent reference carries its own
discretisation error, and this is an L∞-normalised SHAPE comparison on a
heterogeneous 3-region closed cylinder at a loose 1e-1 bar. It is not a gold
standard. But it is INDEPENDENT, and 1.92× is not noise.

**NEXT MEASUREMENT (do this before any decision):** build τ on the
**recursion-defined** edges — the ν ladder from α's own recursion, BMC Eq. 43
read forward — and re-run both instruments. That is the partition that makes α
and τ genuinely share one object, which was the campaign's actual goal; the
geometric arc was a mis-identification of it. Three candidates now have to be
measured on the same two instruments: chord (+absorber), geometric arc,
recursion-defined.

⚠ **DO NOT relax `tol_per_cell`** and do not close 6.4 until this is settled.
The commits `3dda18ca`/`d5067c4d`/`39b46a31`/`c9bb61b4` are all on the branch,
so a revert is cheap and available — but the recursion-defined partition may
well be the answer both the structure AND the numbers want, which would make a
revert premature.
