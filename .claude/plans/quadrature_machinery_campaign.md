# Quadrature machinery campaign — the plan of record

> **⏱ COLD-PICKUP ANCHOR — read this block first, then §0. Trust `git`, never this
> file, for merge status** (`.claude/rules/process-discipline.md`).
>
> **What this is.** A deliberate **interlude**: a general improvement of the
> quadrature machinery, opened 2026-08-01 by user ruling. It is NOT a fix for one
> bug. It happens to dissolve #326 when the machinery is right — that is the
> *test* of the design, not its purpose.
>
> **Branch** `refactor/operator-strategy-layers`. **Q0, Q1, Q2, Q3, Q4 ALL
> LANDED** (2026-08-02). **#325 is no longer blocked**: it was stuck behind a
> sort-key ruling that Q4 dissolved (see T19).
>
> **IN FLIGHT: Q5 — THE FOLD.** Design of record written 2026-08-02, reviewed
> and RULED. Read §5's Q5 block, then **T22b / T23 / T24 / T25 / T27 / T28**.
> All three open questions are ruled: the fold reaches the **SOLVE** (every
> cylindrical snapshot moves, and that is *scheduling*, not a cost — "it is a
> matter of WHEN we rebuild the snapshot, not IF"); Q7's mirror-plane naming is
> **pulled forward** into Q5.0.2; and the **exactness space is pulled forward**
> as Q5.E because T28 showed it is the root.
>
> **⏸ CHECKPOINT — 2026-08-02, extended 2026-08-07 twice (compaction before
> Q5.3; compaction before Q5.5) and 2026-08-08 four times (compaction before
> #327; compaction before 6.3 — after the #325/#337 interleave; compaction
> MID-6.3 after legs 1/2a/2b; compaction POST-FLIP before leg 6 — legs
> 2c/2d/2e/3/4 AND leg 5 THE FLIP are LANDED, `2befb14d`..`1689faf4`; the
> "6.3 LANDED SO FAR" + "SECOND HALF" + "LEG 5 LANDED" ledger blocks in §5's
> Q5.6 block are the per-hash authority for the flip's own legs).
> Landed so far, in order:**
>
> | step | commit | what |
> |---|---|---|
> | **Q5.0.1** | `a7695148` | the reflection partner map is **certified** — caught a live defect (odd-`n_φ` products had an involutive-but-wrong `σ_x` map feeding the `r=0` pole continuation) |
> | **Q5.E / E1** | `9e74faa1` | `numerics/exactness.py` — a claim names the SPACE its degree indexes |
> | **Q5.E / E2** | `34a97f43` | `DiscreteMeasure` carries ONE claim; `tensor_with` vs the direct sum; the "6.2832 bug" re-measured and partly **inverted** |
> | **Q5.E / E3** | `647ece5a` | `rules_circle.py` — the periodic trapezoid on `S^1`, shift as a `Fraction`, nodes as roots of unity |
> | **Q5.E / E4** | `255a9f1c` | the sphere product theorem + the azimuthal substitution; #325's generator half closed |
> | **Q5.0.2** | `bfedc621` | `Z2` RETIRED → `Mirror(axis)`; the 1-D and 3-D arms unified *(hash filled at the 2026-08-07 reconciliation — the row shipped saying "this commit")* |
> | **Q5.1** | `681bc49b` | `DiscreteMeasure.quotient(group)` — the composite named, weights orbit-stabilizer-derived, refusal off the certificate; the "idempotent" gate realized as TWO arms (see the block's landing note) |
> | **Q5.2** | `e0554df9` | `spherical_product(polar, azimuthal)` — the combinator seam; group DERIVED at construction (3 `preserves` generator checks, refusal names the unrealized family); STAGGERED ⟹ `Σ = ∅` ⟹ the fold composes with NO flag |
> | **Q5.3** | `65ff5bb0` | `fiber()` RETIRED → `LevelStructure.quotient(parent=, onto=)` — the structure DESCENDS by selection (charts bit-exact, level order = the parent's η-order RESTRICTED to survivors); fiberwise certified by per-level mass; ⛔ the plan's `==` gate spelling REFUTED (the two orders are exact REVERSES — see the block's landing note) |
> | **Q5.4** | `679cba61` | the R12a predicate re-posed on TWO named integer facts — `MarchStart(on_edge_node, degenerate)`, carrying ⟺ neither; `_SEED_TAU_EPS` + the `(0,1)` interval RETIRED; the τ_raw trichotomy demoted to a bit-exact gated THEOREM; ⛔ T26's 5.6e-16 flip was PRE-E3 (see the block) — post-E3 the trichotomy is bit-exact at every parity; the theory page's "odd count would carry" falsehood corrected in place |
> | **Q5.5** | `0ec701d4` | the clamp SPLIT per T27 — the `[0,1]` membership PROMOTED to a raising guard on the raw producer (`_assert_tau_within_unit_interval`, both arms; T22's ω-mis-ordering stops at source); the `[½,1]` absorber STAYS until Q5.6 (docstrings re-posed on the two-objects finding); the T27 mechanism gates at n_φ ∈ {8..64} (Σ=∅ COMPUTED + τ_raw ⊂ [1/5,4/5] + reversal identity BIT-EXACT — the Q5.3 selection descent upgraded T27's ~1e-15 to 0.0); label `morel-montry-folded-arc`; ⭐ the guard immediately CAUGHT a live latent defect (7 sphere_2g fixtures fed LS(4) raw 3-D to the sphere arm, τ ∈ [−20.3, 1.13] consumed silently — re-propped to GL(4), #336 filed) |
> | **Q5.6a** | `8416253c` | `Quadrature.quotient(group)` (the thin lift) + `folded_product(n_mu, n_phi)` (offset DERIVED per T25, odd n_φ refused, every level a carrying arc; the σ_x pole map gated as surviving the fold) |
> | **Q5.6b** | `70b5b6d2` | ⭐⭐ the System-B seed machinery was SPHERE-HARDCODED in 3 places, caught by the wiring's own L0 battery (82 %/158 % flat-equilibrium errors, Mode-12-masked in every two-spelling test): the q̄½ Legendre fold at μ=±1 → the arc's GC1 endpoint synthesis; the \|μ\|=1 march → PATH widths Δr/sinθ_p; the ½ emission → 1/Σw. Sphere byte-identical (243 regression green); **folded cylinders solve at machine L0 incl. scattering** |
> | **Q5.6c** | `b99fcbc3` | the folded rule's harmonic machinery binds the σ-EVEN sub-basis (`MirrorEvenSphericalHarmonicBasis`, rectangular layout, odd columns zeroed, parity DERIVED by mirror-evaluation); `Quadrature.folded_by` provenance marker; the +6.49 ξ-moment garbage channel bit-closed; even-block Gram == continuum at 1e-15 |
>
> **▶ NEXT = Q5.6 RESUMES at 6.3 (the flip) — `SNMesh(CYLINDRICAL)`
> admits exactly the carrying rules.** Existence-checked at
> pointer-write (2026-08-08): the mesh already READS the march-start
> facts — `augmented_mesh.py:803` calls
> `march_start_structure_per_level` as the R12a seed-presence filter
> (which levels carry ψ½ blocks) — so 6.3's deliverable is NOT that
> call; it is the ADMISSION refusal (a non-carrying rule on a
> cylinder mesh REFUSES at construction, off the same two named
> facts) plus the call-site migration onto `Quadrature.folded_product`
> (both symbols exercised this session: the fold gate file imports
> both). The ladder + standing warnings are in the block below and
> §5's Q5.6 block.
>
> **#327 CLOSED @ `414f2cb6`** (the citation half; Fork A was
> `cf2e8a07`/`df33913d`/`c3475228`): the seed `μ₁² = 4/(N(N+2))` is a
> **PROJECT CONVENTION** — the source leaves μ₁ FREE (LA-3251-MS
> printed p. 32, scan-verified; the report is the 1968 Greenspan
> chapter's report form), the published tables are LA-3186's
> moment-matched family (S4: `0.3500212` vs ours `0.4082483`, exact
> contrast — single orbit forces the weights), and the
> positivity-frontier hypothesis is **CONFIRMED**: theirs positive
> through n = 22 vs our S12 ⟹ the seed, not the shape, limits the
> family. Corpus note resolved in place; memo
> `scratch/issue_327_mu1_provenance.md`; bib + `LathropCarlson1964`;
> two present-tense-false docstrings fixed with it
> (`level_symmetric_sn`'s equal-weight paragraph; both "Common
> values: 4, 8, 16" traps — S16 refuses). The family upgrade
> (moment-matched μ₁ per order, frontier S22) is **#337** — a NODE
> change moving every LS consumer; NOT in this campaign's scope
> unless the user pulls it in.
>
> **THE QUADRATURE-ISSUES INTERLEAVE (user-ruled 2026-08-08, before 6.3): #325 + #326-hygiene + #337 — ALL LANDED.**
> `bad5223e` #325 remainder (MoC rides the group action; last approximate
> mirror → index arithmetic; derivations twin repointed; eps-gate
> family-complete; #326's ords[0] second-tie-break collapsed) ·
> `e76fd656` #337 step-0 pre-carve captures · `59bb38a0` **the
> moment-matched seed** (μ₁ = smallest root of ∫μ_z^N exact; degrees
> N+1 at S4–S10/S14; frontier S12→S18, 50-digit-confirmed; G1–G8 gates
> + the derived per-order tolerance ladder) · `36f81d8f`/`be8f9dfb`/
> `027ae222`/`ec5b905e`/`be52e3b8` the five re-baseline families (the
> standing GL 1-ULP 7th red + the affine-CYL deliberate reds + the T4b
> pair ABSORBED BY NAME — task #51's known-red set shrinks to the
> spherical-inward diamond) · `ec076008` the T4 capture script was a
> STALE TWIN of its gate's arm (missing #282 route-(a)'s block sum) —
> caught because a re-capture greened every changed arm and reddened
> the unchanged one · `0e102421` the closing docs (Closes #325, #337).
> ⭐ The delegated S20/S22 re-evaluation closed BY MEASUREMENT: LA-3186's
> axis-weight decomposition is derivable-systematic but reaches full 3-D
> degree 11 at every order ≥14 (ours: 15/15/17); frontiers INTERLEAVE;
> S22 LP-certified infeasible. Extractions:
> `scratch/issue_337_la3186_la4058_extraction.md` + the verification
> plan + `scratch/issue_325_326_remainder_map.md` (the #326 remainder
> rides 6.4/6.5 as planned). Batteries: 247 gate rows + 4433-test
> medium slice + the 108-row artifact re-run all green; mutation
> battery M0/M1/M2 = 52/33/3 reds (the plan's predicted pattern).
> **▶ NEXT remains Q5.6's 6.3 flip** — unchanged by the interleave;
> cylinders still admit LS today and the CYL re-baselines are marked
> re-baseline-now-retire-soon for exactly that step. Pointer symbols
> re-checked at this compaction (2026-08-08): `augmented_mesh.py`
> reads `march_start_structure_per_level` (the R12a filter — 6.3's
> deliverable is the ADMISSION refusal, NOT that call);
> `Quadrature.folded_product` live and gated.
>
> ⚠ **NEW standing corrections for the 6.3 pickup (from the
> interleave):** (a) `_CYL_LS` in
> `tests/sn/mesh/test_radial_characteristic_slot_coordination.py:101`
> now builds cylinder+LS meshes at orders `(2, 4, 8, 18)` — S18 was
> added at #337; ALL its rows become refusal cases at the flip and
> must re-pose or migrate with the call-site sweep. (b) The T4 capture
> script (`tests/sn/_fixtures/wave_t_t4/_capture_pre_t4_snapshots.py`)
> and its T4c gates consume an LS4 CYLINDER — part of the migration
> list, and the script must be RE-RUN (not hand-edited) after the
> migration, per the stale-twin lesson at `ec076008`. (c) The
> LS family now serves S2..S18, but NOTHING about the flip's refusal
> logic changes: LS cylinder levels are `degenerate` (η 4-to-1) at
> every order — the level structure is seed-independent (`[M]` gated
> by the step-0 literals).
>
> **Batteries at this compaction (HEAD `34355c2a`):** 247 gate rows +
> the 4433-test medium slice (numerics/geometry/sn-operators/
> primitives/mesh, `-m "not slow"`, 7:15) + the 108-row artifact
> re-run + the 54-row streaming file ALL GREEN; sphinx `-W` clean
> every leg; xref 0 dead; ratchet green; mutation battery M0/M1/M2 =
> 52/33/3 (the plan's predicted teeth). ⚠ tests/sn's FULL not-slow
> slice was NOT re-run this leg (last full: 2890/6 at `fb38ab31`) —
> the known-red set SHOULD now be the spherical-inward diamond alone
> (T4b ×2 + affine-CYL ×2 + the GL 7th absorbed by the re-baselines),
> but that count is INFERRED from the absorptions, not re-measured;
> the 6.3 landing's own battery re-measures it. tests/numerics
> not-slow ran green inside the medium slice; the WITH-slow battery
> stays owed at the flip's close, as before.
>
> **Q5.6's remaining ladder (was parked at 6.3)** — its landed legs +
> remaining ladder (6.3 flip: march-start refusal + call-site
> migration + snapshot re-baselines + LS sites → `folded_product`;
> 6.4 acceptance: #229 falls + T3 α gate + the absorber retirement +
> τ ∉ {0,1} + the three xfails re-posed; 6.5 retirements + theory
> pass + #326 close) live in §5's Q5.6 block. ⚠ **§6b for the
> absorber retirement (RULED at Q5.5):** its retirement and the fold
> wiring are ONE call-site set — production NODE_ALIGNED cylinders
> have `τ_raw,0 = 0` bit-exact and divide by zero unclamped, so the
> absorber falls only when the wired fold makes every consumed
> cylinder rule an arc. ⚠ The 2026-08-07 existence-check corrections
> still hold: `_XFAIL_326` lives in
> `tests/sn/sweep/curvilinear/test_azimuthal_mirror_symmetry.py` (NOT
> `_test_helpers`); the `eta_edge` twin at
> `orpheus/derivations/discrete/sn/contamination.py:45+`;
> `product_level_ordering`/`PRODUCT_LEVEL_ORDERINGS` at
> `tests/sn/_test_helpers.py:983+`.
>
> `[M]` **Batteries at the 2026-08-08 pre-#327 compaction checkpoint**
> (HEAD `fb38ab31`, serial `-O`): `tests/sn -m "not slow"` **2890
> passed / 6 failed** (35:10) — the six are EXACTLY the known
> deliberate reds (T4b streaming snapshots ×2, affine-carve
> bit-identity ×3, spherical-inward diamond), zero new; the +19 over
> the Q5.5 close are the Q5.6 gates (fold lift/factory 10, folded
> harmonics 4, System-B cylinder arm 5). `tests/numerics` untouched by
> the Q5.6 legs' test additions (all under tests/sn; the
> numerics-side PRODUCTION edits — directional.py, the basis, the
> mirror_axis accessor — are exercised by the 81-item targeted
> regression at the 6.2c close: SH-space ERR-039 Gram gate, scattering
> kernel crosscheck, harmonic-moment-flux, scattering adjoint, all
> green; the full tests/numerics battery was last run 2005/0 at the
> Q5.3 close and is OWED at the Q5.6 flip's own close). `sphinx -W`
> clean at every leg; pyright ratchet green (49-71 s runs at 6.2c and
> Q5.5); xrefs 0 dead at the Q5.5 close (no docstring-xref-bearing
> edits since were left unbuilt — every leg rebuilt sphinx -W clean).
> The #33 Lambertian snapshot (7th red, slow set) stands unmeasured
> since checkpoint #7. Uncommitted-by-policy set unchanged (vv SKILL +
> catalog, qa lessons ×2).
>
> ---
> **⏱ RECONCILED AGAINST THE TREE 2026-08-07 — the G/P interlude landed
> BETWEEN Q5.0.2 and Q5.1, and this plan predates all of it.** Track G
> (geometric transformation consolidation, carved out AT `bfedc621` to land
> before Q5.1 — G1–G6.5 COMPLETE) and track P (the affine boundary source
> channel — COMPLETE) put ~170 commits on this branch that this anchor does
> not know. The premise audit (every §5 forward step checked against the
> tree) found the FOLD LADDER ITSELF SOUND — Q5.1–Q5.6 consume the
> measure-side machinery (orbit certificate, Σ, `pushforward`+`consolidate`,
> `periodic_trapezoid`, `AngularSymmetry`), none of which G/P touched
> destructively — with these updates a resuming session MUST carry:
>
> 1. ⛔ **The reflection table is GONE** (G §7d.3, 2026-08-07):
>    `reflection_index` / `reflection_partners` /
>    `_compute_sphere_reflection_partners` are DELETED. The single source for
>    "which permutation does a motion induce" is
>    **`Quadrature.ordinate_permutation(motion)`** — derivation-at-need
>    through the G2-verified `RigidMotion.preserves`, same embedding
>    (`_embedded_nodes`) and windows as `_orbit_closure`, same three
>    certifications (ERR-074/073/042), returning `None` for a non-symmetry.
>    Every T-item citing `reflection_index` describes the pre-G tree
>    (annotated in place at T8, T23, §7).
> 2. ⭐ **T23's must-precede cleanup is DISCHARGED — beyond its ask.** The
>    `r=0` pole map is now DERIVED per-quadrature at the mesh
>    (`_ensure_pole_mirror` → `ordinate_permutation(SelfPairedDeck.mirror("x").motion)`,
>    4 sites) and REFUSES when no pairing exists ("cannot seed the r = 0
>    pole", pinned at `product(4,5)`+cylinder). Consequence for the fold: if
>    folding ever produced a rule without σ_x-closure, the cylindrical sweep
>    refuses LOUDLY instead of shipping the silent wrong pole map T23 feared.
>    (Geometric expectation: the ξ-fold quotients by σ_y, and σ_x maps the
>    half-range onto itself, so closure should SURVIVE folding — but the net
>    exists if it does not. Gate it in Q5.6 rather than assuming.)
> 3. **`symmetry.py` was rewired under Q's feet — for the better** (G3/G4):
>    it speaks `RigidMotion` end-to-end; the 1-D/3-D arm split is DELETED
>    (`_check_invariance_1d`/`_3d` no longer exist — a new family needs NO
>    per-dimension arm, superseding the Q5.0.2 template's warning); the
>    checker's `C_n`/`σ_v` are re-based on `roots_of_unity`; `Mirror(axis)`
>    SURVIVED and is the 1-D invariance tag. `_NAMED_LATTICE`, `_GROUP_CACHE`
>    (repr-keyed) and `_contains` survive as the template describes.
> 4. **Three mirror spellings now coexist by design, one per tier**:
>    `Mirror(axis)` (the SYMMETRY-GROUP member, numerics), 
>    `SelfPairedDeck.mirror(axis)` (the boundary LAW's deck slot, geometry),
>    `RigidMotion.reflection(normal)` (the MOTION, geometry core). The first
>    two both RESOLVE to the third; the fold speaks the first
>    (`AngularSymmetry` carries `SubgroupOfO3` values). Axis-letter →
>    normal resolution is single-sourced per tier with ONE deliberately-local
>    test-side cross-check (`tests/_harness/references.py::_AXIS_INDEX`,
>    vv #22) — do not "deduplicate" it.
> 5. **The red baseline is unchanged but is now TRACKED elsewhere**: the six
>    SN reds + the geometry seventh below are task **#51**'s merge gate
>    (with #33 for the geometry snapshot and **#333** for the three affine
>    sha256 rows). The full wide slice has been externally stopped three
>    times from background sessions — it runs at #51's merge, `| tee`d.
> 6. **G-side contracts a Q session must not violate**: the SN deck arm
>    REFUSES space-less method spaces (the refusal is the contract, G6.5);
>    `TraceRestrictionOperator` cross-checks its gather against a bound
>    codomain's `ordinate_indices` elementwise (ERR-077) — a fold that
>    mints half-trace-like spaces must keep declared rows ≡ gather rows.
> 7. **Sibling-matcher note for Q7/Q6 work**: `ordinate_permutation`
>    (quadrature tier) and `_orbit_closure` (measure/group tier) are
>    LAYERED, not twins — both consume the same embedding and windows, and
>    the quadrature method documents the relationship in its docstring. A
>    one-grep verification that `_orbit_closure` also routes through
>    `RigidMotion.preserves` post-G3 is owed at Q7 pickup before touching
>    either.
> ---
>
> ⛔ **UNCOMMITTED AND UNCOMMITTABLE: ERR-074** was appended to
> `.claude/skills/vv-principles/error_catalog.md` (user-authorised, 2026-08-02).
> That file is on the **forbidden-to-commit** list and carries other uncommitted
> user state. It is the quadrature-layer sibling of ERR-042 — an uncertified
> partner map that is *involutive but wrong*. **Do not `git checkout` /
> `restore` / `stash` / `clean` it; the content is irrecoverable** (lesson L28).
> The same account is in `a7695148`'s message if the file is ever lost.
>
> ⛔ **THIS ANCHOR PRESCRIBED A FALSE NEXT STEP UNTIL 2026-08-02.** It said
> "**NEXT: #326** — swap `level_indices` to the fiber ordering Q4 built, and see
> whether the three `xfail(strict=True)` rows XPASS". **Measured: that swap
> makes the cylindrical solve produce NaN**, and the same anchor's own §7 said
> "do not fix #326 directly". See **T22** for the measurement, **T22b** for why
> the mechanism is *implementation, not inherent*, and **T23** for why those
> three rows cannot serve as the acceptance test at all.
>
> **The one-paragraph version of Q5.** ξ is the coefficient of `∂ψ/∂ω`, so
> `Σ = {ξ = 0}` (Q1's own object) is where angular redistribution vanishes and
> the closure `α = 0` is *physically* right — the endpoints of the arcs a level's
> circle is cut into. On one arc `η` is strictly monotone in `ω`, so the march
> order and the fiber order **coincide** and #326 has nothing left to decide.
> `Σ = ∅` makes the fold a **free** action, and that condition **derives** the
> azimuthal offset (`δ = π/n_φ`, i.e. Gauss–Chebyshev in `cos φ`) rather than
> leaving it a preference. Every piece already exists and is unconsumed.
>
> ⚠ **A SEVENTH RED lives OUTSIDE `tests/sn`, and the campaign had not been
> running it.** `tests/geometry/test_bc_equivalence_snapshot.py::
> TestWhiteXminPartial03GLSnapshot::test_matches_the_frozen_scaled_lambertian`
> — `[M]` 8 of 60 elements off by **2.2e-16 abs / 1.1e-15 rel** against an
> `rtol` of `8.9e-16`, i.e. over by 1.26×. Bisected 2026-08-02: it fails at
> `292a1ba5` too, so **neither E3/E4 nor Q5.0.2 caused it** — it is older, and
> the likely origin is Q3's GL re-baseline (`bc89b62e`), which re-captured the
> SN slab snapshots and evidently missed this geometry one. Same class as the
> other re-baseline items; needs its own call. **Every "the same six" report in
> this plan is scoped to `tests/sn`** — run `tests/geometry` too.
>
> ⚠ **SIX SN GATES ARE RED**, deliberately, all documented in `81689a58`; each
> needs its own call, none is a physics failure:
> two on a frozen historical `.npz` (`test_streaming_operator`), three on a
> hard-coded **sha256** frozen at pre-carve `63719a2`
> (`test_affine_carve_bit_identity` — re-hashing is a different act from
> re-capturing an array), and `test_diamond`, which is NOT a stored snapshot:
> production vs the test's own hand-written algebra, now `[M]` **exactly 1 ULP**
> apart (2.06e-16 rel). That one wants a small nulp gate, not a re-baseline.
>
> Still open by design, non-blocking: `align_to` (no consumer), the ε retirement
> (mis-scoped per T18b).
>
> | commit | what |
> |---|---|
> | `1f9d4818` | **Q0** — `roots_of_unity` + the three promoted gates + `_test_helpers` |
> | `49bd7314` | this plan + the 8 `scratch/` evidence reports + agent memory |
> | `93287e57` | the boundary campaign's 13 orphaned evidence reports |
> | `35a1fdb1` | `.gitignore` the 215 MB local literature corpus |
> | `d6bc55e4` | **ERR-072 + ERR-073** — two false certifications, with 4 proven catchers |
> | `63d0b234` | **the subgroup GRAPH + walk** (T17) and the `O2 → D_∞h` rename (T15d) |
> | `c9c57152` | **the orbit certificate + `Σ`** (T14, T18) + 2 perf fixes |
> | `4fd0c7b3` | the 1-D polar-marginal action; `Z2` is improper; dead code retired |
> | `e6f01d7e` | **`DiscreteMeasure.consolidate()`** — the second half of a quotient |
> | `9915f15b` | the two sub-agent evidence reports + agent memory |
> | `e7d44f3c` | **the Stage-1 re-pose** (T16/T16b) — `AngularSymmetry`, the computed gate, `product`'s false tag corrected |
> | `c630153e` | **Q3** — `GeneratingMeasure` (Golub–Welsch, families as VALUES) + the `min()` soundness fix |
> | `6a04fcdf` | the slab rule's reference measure; the drift gate given teeth; a doc falsity |
> | `bc89b62e` | **⭐ IMPOSE the symmetry** (T19) — from reading numpy/scipy source |
> | `d6f53afe` | `level_symmetric` advertised `O_h`, realized `D_2h` — the index fix (T20) |
> | `579d5eaf` | **ONE Gauss-Legendre construction** — the consolidation + re-baseline |
> | `81689a58` | the SN baselines re-captured through the project's own generators |
> | `3afb52c2` | **Q4** — a level is a FIBER of an invariant (T21) |
> | `fadb827a` | the Q0–Q4 compaction checkpoint |
> | `6a4a64af` | **T22/T22b/T23** — the α-march is a march in ω ARC BY ARC |
> | `32809d98` | **Q5 design of record** — T24 (`Σ = ∅`), T25 (the derived offset), T26 (the R12a float) |
> | `b54c1066` | **T27** — the clamp is TWO objects fused; both Q5 rulings recorded |
> | `a7695148` | **⭐ Q5.0.1** — CERTIFY the reflection partner map (+ **ERR-074**) |
> | `63049c84` | Q5.0.1's anchor + substep update |
> | `58fe2dc3` | **T28** — the product rule is a WELDED PAIR that re-derives its own properties |
> | `9e74faa1` | **⭐ Q5.E / E1** — an exactness claim names the SPACE its degree indexes |
> | `89af5d15` | the Q5.E ladder (E1 landed, E2–E4 specified) |
> | `34a97f43` | **⭐⭐ Q5.E / E2** — ONE claim on `DiscreteMeasure`; the "6.2832 bug" inverted |
>
> **`[M]` Gate at THIS checkpoint (before E3) — measured, not asserted:**
> ```
> tests/numerics                        1585 passed              (1:29)
> tests/sn -m "not slow"                2553 passed, 6 failed   (16:14)
>                                       the SAME six as 81689a58 — no new ones
> sphinx -b html docs -W --keep-going   build succeeded
> npx pyright orpheus/numerics/            0 errors
> ```
> `tests/sn` is the real blast-radius check and its counts are **unchanged
> across Q5.0.1 and E1–E2** — three carves through `measure.py` (163 external
> callers), `symmetry.py` (91) and `product_mu_phi` (17) that moved **no
> solver's behaviour**. `tests/numerics` grew 1546 → 1585 purely by added gates
> (+18 Q5.0.1, +21 E1).
>
> Earlier in the session, untouched since: `tests/geometry + tests/diffusion +
> tests/data` 909 passed, 4 skipped, 1 xfailed.
> ⚠ `tests/numerics/test_symmetry.py` costs ~**80 s** on its own. Budget it.
>
> `tests/sn` is the real blast-radius check: `symmetry.py` has 91 external
> callers, `measure.py` 163, `product_mu_phi` 17. Its counts are **identical to
> the pre-re-pose baseline** — which is the evidence that matters here: the change
> altered which rules a *geometry* admits and changed no solver's behaviour.
> ⚠ `tests/numerics/test_symmetry.py` now costs ~**80 s** (was ~5 s) — it verifies
> the lattice, the walk, the certificate, `Σ`, orbit-stabilizer and 2 ERR
> catchers. Budget it.
>
> **What Q0 landed:**
>
> | path | note |
> |---|---|
> | `orpheus/numerics/roots_of_unity.py` | #325's primitive. `npx pyright` 0. **Not exported, no consumer repointed** — the repoint is **Q5**'s business (it touches the RANGE/RULE factorization) |
> | `tests/numerics/test_roots_of_unity.py` | the primitive's defining laws |
> | `tests/sn/sweep/curvilinear/test_alpha_closed_form.py` | 20 `l1` + 15 `foundation`; drives **production** `cylindrical_streaming` |
> | `tests/sn/sweep/curvilinear/test_azimuthal_mirror_symmetry.py` | 5 `l1` (**3 `xfail(strict=True)`** = the live defect) + 7 `foundation` |
> | `tests/sn/verification/mms/test_mms_ordering_blindness.py` | 10 `foundation` |
> | `tests/sn/_test_helpers.py` | `product_level_ordering`, `PRODUCT_LEVEL_ORDERINGS` |
>
> `[M]` **Gate, re-verified at the Q0 commit — measured, not asserted:**
> ```
> 305 passed, 3 xfailed, 1 warning in 6.05s
>
> $ pytest --runxfail tests/sn/sweep/curvilinear/test_azimuthal_mirror_symmetry.py
> E  assert 0.11908767547973154 < 1e-11
> E  assert 0.05144789466670699 < 1e-11
> E  assert 0.3083565899072802 < 1e-11
> 3 failed, 9 passed
> ```
> Each xfail reds on **its own** documented assertion — not a swallowed setup
> error (`vv-principles` Mode 8 class 4).
>
> **Evidence reports**, now committed under `scratch/`:
> `quadrature_landscape_survey.md` (the 4-axis inventory + #327's discovery),
> `issue325_main_agent_measurements.md`, `issue325_blast_radius.md`,
> `issue325_verification_plan.md`, `issue326_alpha_theory.md`,
> `issue326_mms_adjudication.md`, `issue326_halfrange_map.md`,
> `issue326_quotient_frames.md`.
>
> **Issues opened by this campaign:** **#326** (the double-cover; restructured
> twice — read ALL comments, the last one is the current statement) and **#327**
> (`level_symmetric` is degree-3 at every order — a LIVE L1 defect).
>
> **ON HOLD, resumes after this campaign** — the boundary-machinery campaign
> (`.claude/plans/boundary_machinery_review.md`, B3.4c landed at `79b5affe`;
> next was **#325 → B4**, with B3.3 and B3.5 also open). #325 and #326 are both
> *downstream* of this campaign now: see §7.
>
> **Two `vv-principles` files are FORBIDDEN to commit** and carry uncommitted
> working-tree state git cannot recover (lesson L28). Path-scope every commit to
> exclude `.claude/skills/vv-principles/{SKILL.md,error_catalog.md}`.

---

## 0. Why this campaign exists — the chain that produced it

Recorded because the chain is the argument, and a fresh session that sees only
the conclusion will re-litigate it.

| step | what happened |
|---|---|
| B3.4a | An `AngularAverageOperator` classifier used `> 0.0` on a direction cosine. Removing it exposed that `product`'s tangential cosines are **not exactly zero**. |
| #325 | Root cause: a node set defined by a symmetry group was being **generated by evaluating a parameterization** instead of by the group action. Primitive built (`roots_of_unity`). |
| #326 | Making the nodes exact would make the per-level `argsort` ties exact — which exposed that the ordering was under-determined. Chasing *that* found the real defect. |
| **here** | The real defect is that the cylindrical level **double-covers the azimuthal range**, and the machinery has no way to say so. Every layer of the chain was a missing *abstraction* surfacing as a numerical symptom. |

**The user's ruling that opened this campaign (2026-08-01, verbatim):**

> *"This is not a fixing the solution problem. This is a quadrature machinery
> design problem, that happens to fix the solution when the machinery is fixed.
> Should half range be a restriction over full-range? Do we need some abstraction
> class for something like Levels? For symmetry? Don't map what is missing in THIS
> quadrature. Map what is missing **in machinery design**."*

And on the two forks that were posed as choices:

> *"Again, this should not be an arbitrary choice. Why do we need to choose this?
> There should simply be a mathematically formal correct machinery."*

**Standing acceptance criterion for every abstraction in this campaign:** a choice
earns first-class status when there are **≥2 genuine realizations** in the
codebase's real or documented-foreseeable needs — because ≥2 realizations let a
mathematical invariant be tested *across* them, which **proves** the implementation
rather than pinning it. (The user's sharpening of `coding-elegance` Pattern 6:
the second instance is a verification instrument, not just evidence of a pattern.)

---

## 1. ⭐ THEOREMS AND INVARIANTS ESTABLISHED — the durable math

These survive any implementation decision. Each is either proved, or measured and
marked `[M]`. **Nothing in this section should ever need re-deriving.**

### T1 — The azimuthal rule is the periodic trapezoid, and its exactness is TRIGONOMETRIC

`rules_product.py:115-116` is `np.linspace(0, 2π, n, endpoint=False)` with equal
weights `2π/n`. That is the **periodic trapezoidal rule**. Its equal weights are
not a compromise: for a periodic analytic integrand the periodic trapezoid is
spectrally accurate (Trefethen & Weideman, *The exponentially convergent
trapezoidal rule*, SIAM Rev. 56(3), 2014).

`[M]` Exactness is exactly `k < n`, with total aliasing at `k = n`:

```
  n=  4:  int cos( 1phi)=-2.89e-16  int cos( 3phi)=+6.44e-16  int cos( 4phi)=+6.28e+00
  n=  8:  int cos( 1phi)=-2.22e-16  int cos( 7phi)=+2.89e-15  int cos( 8phi)=+6.28e+00
  n= 16:  int cos( 1phi)=-3.53e-16  int cos(15phi)=+7.70e-16  int cos(16phi)=+6.28e+00
```

**The three fused choices.** `linspace(0, 2π, n, endpoint=False)` + `2π/n` fuses:
**(a) the RANGE** `[0,2π)`, **(b) the SPACING** equispaced, **(c) the RULE**
periodic-trapezoid. All three are unexamined and none is nameable in the current
types. This decomposition is the user's and it is the campaign's main axis.

### T2 — ⚠ The `degree_of_exactness` tag is CORRECT and TIGHT. The defect is different from what it looks like.

`rules_product.py:148` ships `degree_of_exactness=min(2*n_mu - 1, n_phi - 1)`, a
`min()` over a *polynomial* degree and a *trigonometric* degree — two different
units. That LOOKS like a category error. **It is not, and the plan must not claim
it is.** `[M]` Integrating every spherical monomial `x^a y^b z^c` against the
closed form:

```
  product(2,4) tag=3 | deg3:2.1e-16  deg4:1.1e+00     <- tight
  product(4,8) tag=7 | deg7:1.1e-16  deg8:7.3e-02     <- tight
  product(4,4) tag=3 | deg3:2.3e-16  deg4:8.4e-01     <- tight
```

The tag is exact and **tight**, not "conservative" as its own docstring claims
(that word is the falsehood — a small one, fix it in place). It works because a
degree-`d` spherical monomial's azimuthal factor is a trig polynomial of degree
`≤ d`, so the two bounds coincide.

**The real defect is representational:** a single integer cannot say *which
space* the claim is about. A symmetry quotient CHANGES the space (T5), and the
tag survives the change unaltered — which is precisely the `+2.94` failure mode
of T6. **What is missing is a subspace attached to the degree, not a corrected
number.**

### T3 — The α-recursion has an exact closed form, and it is the literature's DEFINITION

`α_{m+½} = α_{m−½} − w_m η_m`, `α_½ = 0` (production: `reduced_operator.py:781-786`
cylinder, `:688-695` sphere) is a **cumulative quadrature in the azimuthal angle**.
Via the Dirichlet identity `2 sin(t₀/2)·Σ_{m'=0}^{m} cos(m' t₀) = sin((m+½)t₀) + sin(t₀/2)`:

> **α_{m+½} = −sinθ · (π/n)/sin(π/n) · [ sin(ω_{m+½}) + sin(π/n) ]**,  `ω_{m+½} = ω_m + π/n`

So **α is `−ξ` (the tangential cosine) at the half-angle boundary**, up to
`κ = Δω/(2 sin(Δω/2)) → 1`. This is not an inference — Hébert (2009) §3.9.3
Eq. (3.399) *defines* it: *"Defining `α_{p,q±1/2} = 𝒲_p η_{p,q±1/2}`…"* with his
`η` = the tangential cosine = ORPHEUS `ξ = μ_y`; `α_½ = 0` is just "the starting
boundary is the plane `ξ = 0`".

`[M]` The recursion reproduces it **iff** ordinates are enumerated in increasing ω:

```
     n |  omega-order      eta-order (legacy)   lexsort (eta,xi)
     8 |   7.003e-16                3.044e+00          3.044e+00
    16 |   4.466e-16                2.970e+00          2.970e+00
    32 |   8.604e-16                2.905e+00          2.905e+00
```

⚠ **κ is 2.6 % off 1 at the production `n_phi = 8`.** Never gate a bare `α == −ξ`;
it is false. Pinned by `test_kappa_prefactor_is_not_unity_and_converges_to_one`.

### T4 — ψ is EVEN in ξ (a geometric theorem), and production violates it

A 1-D cylinder is mirror-symmetric about the r–z plane, so `(η, +ξ, μ)` and
`(η, −ξ, μ)` are physically equivalent. The product rule with equal azimuthal
weights preserves the symmetry, so the semi-discrete problem inherits it exactly.
**Needs no MMS, no benchmark, no reference solution.**

`[M]` Stock tree, converged, `inner_tol=1e-13`, `max|ψ_n − ψ_mirror(n)|/max|ψ|`:

```
          product(4,4)  hetero   9.797828e-03        product(4,8)  hetero  1.957716e-02
         product(4,16)  hetero   1.175783e-02       product(4,32)  hetero  7.318974e-03
   product(4,8) homogeneous (non-degenerate)        1.190877e-01
   level_symmetric(4)                               3.083566e-01
```

⚠ A homogeneous config with reflective/reflective + uniform source is
**H2-degenerate** — flat exact solution, trivially symmetric, reads `1e-15`. That
is `vv-principles` anti-pattern #4 and it fooled the first probe. Falling in
`n_phi`, **flat in `n_mu`** ⟹ an *azimuthal* defect.

### T5 — Half-range is a QUOTIENT, not a restriction; the weights are DERIVED

Not an ad-hoc fold: the quotient of a G-invariant measure by a symmetry subgroup,
**defined only when the measure is G-invariant**. Orbit-stabilizer **derives** the
weights — `W = w·|G|/|Stab|` gives the trapezoid `[w, 2w, …, 2w, w]`. It is the
**orbifold measure**, not a quadrature taste. Burnside's `(N+F)/2` reproduces the
measured orbit counts:

```
                quad  orbits  fixed  pairs   sum w preserved  xi-even integ
product(4,8)            20      8     12              True           True
product(4,16)           36      8     28              True           True
```

The in-tree discriminator between *restriction* and *quotient* is **total mass**:
`restrict` drops it, `partition_by` preserves it, the fold preserves it.

### T6 — The ξ-even/odd split IS the trivial/sign isotypic decomposition of ℤ₂

Not an analogy. The character count `dim V₊^{(ℓ)} = ℓ+1` reproduces the measured
`ℓ=1: 2 even, 1 odd`. `[M]` A ξ-odd integrand on the quotient:

```
  full sphere  int(xi) = +4.441e-16     quotient  int(xi) = +5.993e+00
```

That is **not an error** — it is the quotient reporting that ξ-odd functions are
**not in its space**. The reported `φ_1^ξ : −1.3e-16 → +2.94` is the signature of
quotienting the *measure* but not the *basis*: an inconsistent pair, i.e. a
**domain check firing**, not a cost of the design.

⚠ `dim V_{A₁}^{(ℓ=1)} = 1` for `C_{2v}` vs `2` for `⟨σ_y⟩` ⟹ **a 4-fold and a
2-fold rule on the same geometry have DIFFERENT trial spaces.** The restriction
MUST be derived from the group the measure carries. Never hardcode it.

### T7 — A level is a FIBER of an invariant, not an orbit — and the two producers fiber differently

Corrects an earlier main-agent claim. `product` fibers over **signed** `μ_z`;
`level_symmetric` fibers over **|μ_z|**. One type (`LevelStructure`), two
meanings — the direct cause of the two distinct `τ_raw` fingerprints.

Consequently **level-partition and symmetry-fold are NOT two rungs of one chain**:
one is a fibration by an invariant, the other a quotient by isotropy. **Fibered,
not towered.** Inside a fiber the folds DO commute (`C_{2v}` is abelian), so
`level_symmetric`'s 4-to-1 needs **no new object** — same object, bigger `G`.

### T8 — The angular fold and the spatial `r=0` pole map are ONE decision

`σ_x σ_y = C₂(ẑ)` exactly. The physical centreline continuation is `C₂(ẑ)`; the
code uses `σ_x` (`loss_representation/__init__.py:4189` via
`reflection_index("x")`). `[M]` They differ on **56 of 64 ordinates** and coincide
**only on the quotient**. The tree is right by a coincidence whose precondition it
has never taken.

⛔ **MECHANISM SUPERSEDED by track G (annotated 2026-08-07); the THEOREM
stands.** `reflection_index` is DELETED (§7d.3). The σ_x pole map is now
DERIVED per-quadrature at the mesh — `_ensure_pole_mirror` →
`Quadrature.ordinate_permutation(SelfPairedDeck.mirror("x").motion)`, four
sites in `loss_representation` — certified at derivation (ERR-042/073/074)
and REFUSING when the pairing does not exist ("cannot seed the r = 0 pole",
pinned). So the "coincidence whose precondition it has never taken" now has
its precondition TAKEN at every derivation: a rule not closed under σ_x
cannot seed the pole at all, loudly. The fold's Q5.6 acceptance should gate
that the FOLDED rule still derives its σ_x pairing (expected: yes — σ_x
maps the ξ ≥ 0 half-range onto itself — but gate it, don't assume it).

Corollary, and it corrects a plausible-sounding wrong claim: ERR-042 checks `σ_x`
(`φ→π−φ`, needs even `n_φ`) — the **centreline's** condition. The ξ-fold closes
under `σ_y` for **every** `n_φ`, so the fold exists unconditionally. ORPHEUS
checks the wrong group element and agrees with the right one only because
`σ_y`-closure is automatic.

### T9 — Two structural criteria that fall out for free

- **Angular redistribution (the α term) exists ⟺ the spatial-reduction group acts
  non-trivially on the fiber.** Translations act flat (slab: no α); rotations do
  not (cylinder: α). `SO(2)_ẑ`'s fiber action **IS** the α term.
- **`n_levels > 1` ⟺ the point isotropy group is finite.**

One criterion explains slab-has-no-α, sphere-has-one-level, cylinder-needs-levels.

### T9b — ⛔ `level_symmetric` is degree-3 at EVERY order — filed as #327

The campaign's first independent finding, and the most severe thing in this file.
`Quadrature.level_symmetric(sn_order=N)` advertises `degree_of_exactness = N − 1`;
`[M]` its measured degree is **3 for every N**, an over-claim of **12** at S₁₆.
Verified against the closed-form monomial integral, with controls:

```
         level_symmetric    4           3         3           0
         level_symmetric    8           7         3           4
         level_symmetric   16          15         3          12
   -- controls --
            product(4,8)      advertised=7    MEASURED=7
           product(6,12)      advertised=11   MEASURED=11
             lebedev(17)      advertised=17   MEASURED=17
```

Robust from `atol=1e-12` to `1e-2` — not a tolerance artifact. Mechanism: `[M]`
exactly **ONE distinct weight** in the whole rule at every order
(`w_octant = 4π/(8·n_octant)`), i.e. **not** the Carlson–Lathrop moment-matched
construction the docstring cites. The degree-3 it reaches is free from the `O_h`
orbit symmetry, not from the weights.

Production consequence — the discrete SH Gram vs the analytic `4π/(2ℓ+1)`, which
`quadrature.angular_frame(L).gram` assumes:

```
  N= 8 l=2: 9.687%    l=3: 10.726% (off-diag 8.769e-02)
  N=16 l=2: 18.000%   l=3: 24.970% (off-diag 1.118e-01)
```

⚠ **It gets WORSE with N.** `ℓ=0,1` are exact — which is why isotropic and P1
look healthy and anisotropic `ℓ ≥ 2` does not. No test catches it: the exactness
tests stop at degree 2 (inside the sector `O_h` makes free) and the rest pin the
**tag**, not the property.

**This is the campaign's thesis proved on independent ground.** The defect is not
arithmetic — it is that an integer with no attached function space (T2 / gap 3)
cannot be checked against anything, so a false claim shipped and a selector
consumes it.

### T12 — ⭐ EVERY Gaussian rule is ONE construction: Golub–Welsch on a weight function

User question, 2026-08-01: *"if Gauss-Legendre and Gauss-Chebyshev are either
polymorphism at construction of some rule, or thin sub-classes of some parent
class, that seems like something we should build."* The answer is **stronger than
subclasses — they are not even polymorphism. They are DATA.**

A weight function `w` on an interval fixes a family of orthogonal polynomials,
which obey a three-term recurrence `(α_k, β_k)`. Assemble the symmetric
tridiagonal **Jacobi matrix** `J = tridiag(√β_k, α_k, √β_k)`; then

> **nodes = eigenvalues of `J`;  weights = μ₀ · (first component of each unit eigenvector)²**,  `μ₀ = ∫w dx`

That is **Golub & Welsch (1969)**. The named families differ ONLY in `(α, β, μ₀)`:

| family | interval | weight `w(x)` |
|---|---|---|
| Legendre | `[-1,1]` | `1` |
| Chebyshev-1 | `[-1,1]` | `(1-x²)^(-1/2)` |
| Chebyshev-2 | `[-1,1]` | `(1-x²)^(+1/2)` |
| Jacobi | `[-1,1]` | `(1-x)^α (1+x)^β` |
| Laguerre | `[0,∞)` | `x^α e^{-x}` |
| Hermite | `ℝ` | `e^{-x²}` |

`[M]` One body reproduces every family the tree hand-rolls separately:

```
  Legendre  n= 16: vs orpheus.gauss_legendre   dnode=4.44e-16 dweight=1.26e-15
  Chebyshev n= 16: vs orpheus.gauss_chebyshev  dnode=2.22e-16 dweight=1.47e-15
  Laguerre  n=  8: vs numpy laggauss           dnode=1.78e-15 dweight=3.00e-15
```

**Gauss-Lobatto / -Radau are NOT a family.** They are the same construction with
**prescribed nodes**, obtained by modifying the last entries of `J`
(Golub 1973). So the taxonomy is exactly **two orthogonal axes**:

1. **The measure** (weight function) → the recurrence → nodes + weights.
2. **The node constraint** — free (Gauss) / one endpoint (Radau) / both (Lobatto).

### T12b — ⭐⭐ The generating measure IS the exactness space. One object, not two.

This is the campaign's central unification, and it collapses gap 3 into T12.

`[M]` `gauss_legendre` and `gauss_chebyshev` both ship
`degree_of_exactness = 2n − 1` **in the same field, meaning different things** —
GC's own docstring even says *"with respect to the weighted integral"* and *"the
weight function is in the integral, not in the quadrature"*. The type erases the
distinction the prose states (`coding-elegance` anti-pattern #20). Measured at
`n=4`, `q(x) = x⁶`:

```
    GL  sum w*q = 0.285714285714   int q dx            = 0.285714285714
    GC  sum w*q = 0.981747704247   int q (1-x^2)^-1/2  = 0.981747704247
    GC against the UNWEIGHTED integral: err = 6.960e-01   <- NOT exact
```

> **The weight function that GENERATES the rule is also the measure its exactness
> claim is ABOUT.** So "the missing Gauss abstraction" and "the missing exactness
> subspace" (gap 3) are the SAME object seen from two sides.

**Corollary — this is why #327 (T9b) was possible.** `level_symmetric` has **no
generating measure at all**: it assigns one weight to every ordinate by hand, so
its exactness claim is an unconstrained integer someone typed. A rule constructed
*from* its weight function **cannot over-claim**, because its exactness follows
from the construction. That is `coding-elegance` Pattern 4 applied to quadrature:
**a rule built from its measure cannot lie about its exactness space.**

### T12c — ⭐ The REALIZATION decision: families are VALUES, not subclasses

User question, 2026-08-02: *"polymorphism at construction of an object or
sub-classing? … decided based on benefits, such as elegance, performance."*
Decided by measurement, not taste.

**Performance — measured, and it removes the only argument for an override hook.**
The generic Golub–Welsch body is **2–6× FASTER** than the specialized routine it
would replace, and both sit at machine precision:

```
      n |  Golub-Welsch max err  numpy leggauss max err     GW time   leggauss time   ratio
      8 |             3.214e-16               5.553e-17     34.5 us       211.5 us   0.16x
     32 |             2.468e-16               6.473e-17    153.9 us       679.6 us   0.23x
    128 |             2.665e-16               6.045e-17   2444.9 us      4174.6 us   0.59x
```

`leggauss` is ~4–5× more accurate **in the last ULP** (0.3 ULP vs 1–2 ULP) because
it runs the same eigen-decomposition and then Newton-refines. Irrelevant here:
quadrature is built ONCE at setup and feeds a solver running at `1e-13`. `[M]`
The generic body preserves exactness to `2n−1` (max err `4e-16`–`8e-16` over
degrees `0..2n−1`).

⚠ Honest scope on the failure side: at `n=8` the first failing degree is
unmistakable (`4.7e-5`), but at `n=32/64` the *absolute* error at degree `2n`
reads `~1e-16` because `x^128` on `[-1,1]` collapses in magnitude. Do NOT gate
"fails at `2n`" on an absolute error at high `n` — use a relative measure or a
weighted integrand.

**Therefore there is NOTHING to override, and that decides the shape.** A subclass
exists to vary behaviour; here every family supplies the same three things
`(α_k, β_k, μ₀)` and ONE generic body consumes them. A subclass whose entire
content is returning three arrays is `coding-elegance` Pattern 6's *"the instance
with extra ceremony"*.

**But the weight IS a genuine type** by the project's own decidable criterion
(`.claude/rules/coding-standards.md`, "Type vs property"): (a) ≥2 non-isomorphic
realizations — yes; (b) a **non-identity morphism** is actually applied — yes,
three of them: the **affine remap** `[-1,1] → [a,b]` (changes `μ₀` by the
Jacobian, shifts `α`, scales `β`); the **Christoffel modification** (given the
recurrence for `w`, derive it for `(x−c)·w`); and the **endpoint constraint**
(a rank-modification of `J`).

**The realization, by axis:**

| concept | realization | why |
|---|---|---|
| **weight family** | **VALUES of one frozen type** — module constants `LEGENDRE`, `CHEBYSHEV_T`, `CHEBYSHEV_U`; factory functions `jacobi(a,b)`, `laguerre(a)`, `hermite()` | many cases, no varying behaviour ⟹ data. Parameterized families are *parameterized values*, which is how the math reads them |
| **construction** | **ONE generic body** (Golub–Welsch) | measured faster AND machine-accurate ⟹ no fast-path hook is earned |
| **node constraint** (Gauss / Radau / Lobatto) | **closed sum type + exhaustive match**, applied as a morphism on the rule | 3 closed cases, operations may grow — the expression problem's *dual*, so a match, not polymorphism |
| **domain topology** (interval / circle) | **closed sum type + exhaustive match** | few closed cases; selects WHICH construction (Golub–Welsch vs the circle's uniform/DFT rule — T1) |

**No subclassing anywhere.** The expression-problem discriminator
(`coding-elegance`, "a repeated conditional is a missing type"): weight families
are the axis that GROWS while operations stay fixed *and generic*, so they are
data; node-constraint and domain-topology are small CLOSED sets whose operations
may grow, so they are exhaustive matches.

**Retire, don't wrap.** The three hand-rolled bodies (`gauss_legendre`,
`gauss_chebyshev`, `gauss_laguerre`) become the ≥2 instances that PROVE the
generic body (measured agreement `1e-16`–`3e-15`), then retire. If anyone later
wants `leggauss`'s extra ULP it is an *optimization* gated by an equivalence test
against the generic body — do NOT build it now (Pattern 6).

### T12d — ⭐ Enumeration/naming: the registry EXISTS and is already the right shape

User question, 2026-08-02: a UI wants a two-level drop-down — pick the **family**,
then the **member** — and asks whether the codebase's *"abstract class that
subclasses register with during construction"* pattern fits.

**The mechanism already exists.** `orpheus/numerics/quadrature/registry.py` has
`QuadratureSpec` — a frozen dataclass carrying `name`, `factory`,
`parameters: dict[str, type]`, `invariance_group`, `degree_of_exactness_for`,
`expected_node_count`, four structural trait flags, and `structural_flags()` —
plus `quadrature_registry: list[QuadratureSpec]`, a module-level list of
**VALUES**. That is exactly the shape T12c derives. Do not mint a new mechanism
(clean-before-extend).

**The project uses BOTH patterns, and where each is used IS the discriminator:**

| pattern | where | what is registered |
|---|---|---|
| `__init_subclass__` | `numerics/registry.py` (generic), `transport/displacements/_displacement.py`, `transport/spatial/scheme.py` | classes with **behaviour to override** |
| explicit value registry | `quadrature/registry.py`, `CLOSURE_REGISTRY`, `BOUNDARY_OPERATOR_REGISTRY` | **data / specs** |

> **Discriminator: does the registered thing have behaviour to override?**
> Yes ⟹ class + `__init_subclass__`. No (it is data) ⟹ value registry.

Per T12c the Gauss families have **no** behaviour to override — one generic
Golub–Welsch body consumes their `(α, β, μ₀)`. So: **value registry**. Forcing
them into `__init_subclass__` would (a) require classes whose entire content is
three arrays, and (b) inherit **import-order fragility** — a subclass in an
unimported module is invisible, which is disqualifying for a UI that must
enumerate *everything*.

**But `__init_subclass__`'s one real benefit is worth keeping: you cannot
define-and-forget.** Obtainable for values, by making the registry the single
**declaration** site rather than a second one. NOT "define the family, then
register it" (two sites ⟹ Pattern-2 drift), but a registering constructor so the
constant and the registry entry come from ONE expression:

```python
LEGENDRE = _family("legendre", "Gauss–Legendre", ...)   # declares AND registers
```

### What is actually missing for the two-level UI

1. **No FAMILY grouping** — the registry is flat (4 entries). And the grouping is
   **not cosmetic**: it IS the construction axis (Golub–Welsch on a weight /
   the circle's uniform rule / tabulated orbit rules / product constructions),
   which T12c already establishes as a closed sum type. **The two-level UI
   mirrors the mathematics**, which is the master standard, not a concession to
   presentation.
2. **`parameters: dict[str, type]` is too weak, and it is ALREADY a twin source of
   truth.** It cannot distinguish "any `int ≥ 1`" (`n`) from "one of a **discrete
   tabulated set**" (Lebedev `order ∈ {3,5,…,47}`) from "any `float > −1`"
   (Jacobi `a, b`). Those are three different widgets AND three different
   validations. The admissible domain is currently encoded *implicitly* in
   `degree_of_exactness_for` (it returns `None` past Lebedev 47) and **not at all**
   in `parameters` — so a `ParameterSpec` is earned on **existing evidence**, not
   on UI taste, and it serves constructor validation and the degree inverter as
   well as the UI (≥2 consumers, T-criterion satisfied).
3. **Display label vs identifier** — `"GaussLegendre1D"` is an identifier, not a
   human label; there is no description field.
4. **The Gauss families are not registered at all.** Registering them gives
   `gauss_chebyshev` / `gauss_laguerre` / `equispaced` — currently shipping with
   **zero production consumers** — their consumer, and closes that waste.

### T10 — Closure is Mode-12 blind and CANNOT gate any of this

`α_{M+½} = 0` is a **telescoping sum**: `[M]` it holds to ~1e-16 under every
ordering including a **random shuffle**. `vv-principles` anti-pattern #8 / Mode 12
exactly. Worse than blind — the cylinder has **no closure assertion at all** (the
sphere has one at `reduced_operator.py:693`), `verify_alpha_closure` rebuilds η
from a fresh `linspace` and never touches `level_indices`, and the three α
identity tests at `test_pole_angular_closure.py:142-188` build from
`Quadrature.gauss_legendre` only — the product rule is never instantiated.

### T11 — The MMS is exactly Mode-7 blind to all of it

Both cylindrical MMS ansatzes are functions of `η` and `ξ²` only (`mms/sn.py:3778`
ψ, `:3820-3824` Q) — and the mirror pair shares both. `[M]` 22/22 tests PASS
identically under legacy / lexsort / stable, the two rate gates agreeing to
`3e-12` and `9e-15`, both reporting order 2.003. **The MMS cannot adjudicate a
ξ-odd-sector question.** A ξ-odd companion ansatz sees it but does not adjudicate
either (two orderings converge to the same angular floor from opposite sides).

### T13 — ⭐ The root of `Σ` is the GROUP. Not the space, not the quadrature.

User question, 2026-08-02: *"It's very unlikely that this is an orphan convention.
It certainly originates somewhere. […] Where is the root of this?"* Answered by
tracing the import graph, and it settles the siting for good.

**`orpheus/numerics/symmetry.py` is the root**, and both the quadrature layer and
the space layer descend from it:

- It is the **lowest** of the three. Its only import is `.measure`, and that edge
  is one-way by construction — `measure.py:82` carries the comment *"symmetry
  imports FROM this module, so a runtime import would cycle"*, so its own
  `SubgroupOfO3` reference sits under `TYPE_CHECKING`.
- It **already owns the mirrors**: `_reflections(axis)` (`:750`) returns exactly
  σ_x, σ_y, σ_z; `_vertical_mirrors(n)` (`:787`), `_octahedral_ops` (`:802`),
  `_icosahedral_ops` (`:827`) sit beside it. The reflections whose fixed-point
  sets ARE `Σ` are already there.
- The edge is **already walked**: `quadrature/rules_1d.py:30` does
  `from ..symmetry import SubgroupOfO3`. `directional.py` simply never took it.

> `Σ = Fix(σ_n̂) = {x : σ_n̂ x = x}`, and the orbifold singular set is
> `{x : Stab(x) ≠ {e}}` — a function of **(group, point set)** and nothing else.
> A face normal and a quadrature rule are both *inputs* to it, never its home.

**Corollary — the duplication was caused by SITING** (`.claude/rules/coding-standards.md`
L32). The shared constant lives in `numerics/spaces/angular_trace_space.py`, which
is ABOVE one of its consumers: `numerics/quadrature/directional.py` knows nothing
about faces and must not import a face-trace space. Unable to read it legally,
`directional.py` minted `_OCTANT_SIGN_EPS = 1e-15`. **The twin exists because the
primitive was sited too high**, which is L32's failure mode exactly.

**Naming (user ruling, 2026-08-02).** Collection = **`singular_set`**; individual
member = **`singular_point`**; the hyperplane itself = **`mirror`** (Thurston's
orbifold term for a reflection's fixed set; `wall` in Coxeter vocabulary).
`singular_point` beats `mirror_point` because it stays true for **cone points** —
a pure-axis ordinate `(0,0,±1)` is fixed by every rotation about ẑ, so it is
singular without lying on any mirror. `tangential` stays at the boundary layer as
the honest face-normal **specialization**, reading down to `singular_set`.

⚠ Do NOT re-inherit the earlier framing *"six ad-hoc ε-detectors for one object."*
It is wrong twice — see §3 gap 1 as rewritten.

### T14 — ⭐⭐ `Σ` is NOT a floating-point question. It is `Fix(π)`, an integer identity.

The campaign's sharpest single result, and it dissolves gap 1 rather than tidying it.

`_orbit_closure` (`symmetry.py:904`) already computes, for each node `i` and each
group element `M`, the matched index `j` — then discards it to return `bool`
(that is L-013). But observe what `j == i` **means**: `M·x_i = x_i`, i.e.
`x_i ∈ Fix(M)`.

> **`Σ` is the fixed-point set of the orbit permutation. Membership is
> `π_M(i) == i` — an integer identity. Exact. No tolerance.**

The only surviving tolerance is the one that matches nodes while *building* `π`,
which is `symmetry.py`'s existing single `atol` — the one place the question is
honestly numerical. The three scattered ε-detectors do not get consolidated onto
one ε; they **stop being ε-questions**.

`[M]` The measurement that independently confirms the ε was never doing real work
— across 29 production rules (`gauss_legendre` ×6, `product` ×7,
`level_symmetric` ×6, `lebedev` ×10):

```
  ordinates in the disputed band [8.88e-16, 1e-15) : 0
  largest cosine BOTH constants call zero          : 1.806e-16
  smallest cosine BOTH constants call nonzero      : 2.435e-02   (2.7e13x separation)
```

Any constant in `[2e-16, 2e-2]` classifies identically on every production rule.
So the ε disagreement (`TANGENTIAL_EPS = 4·eps ≈ 8.88e-16` vs the two private
`1e-15`) is **latent, not live** — a single-source-of-truth defect, never a
numerical one. Probe: `$CLAUDE_JOB_DIR/tmp/q1_tangential_band.py`.

**Corollary — `Σ` is only expressible where the certificate exists**, i.e. on a
G-invariant set. That is exactly T5's precondition (*half-range is defined only on
a G-invariant measure*) reappearing as `coding-elegance` Pattern 4: the
illegal state is unrepresentable because the type that names `Σ` cannot be
constructed without the closure proof.

### T15 — ⭐ Containment vs SUBCONJUGACY: two standard relations, and the module had neither

User ruling, 2026-08-02: *"An arbitrary choice is not the answer. We need a
concept that is mathematically grounded in group theory (how do group theorists
tackle this?) and that works for curvilinear and cartesian geometries."*

Group theory keeps **two** named relations, and conflating them was the defect:

1. **Containment `H ≤ K`** — literal, set-theoretic, **setting-dependent**. This is
   what "subgroup" means. Crystallography (International Tables Vol. A) always
   lists subgroups *with their setting* for exactly this reason.
2. **Subconjugacy `H ≼ K` ⟺ `∃g ∈ O(3): gHg⁻¹ ≤ K`** — **setting-independent**. The
   poset of *conjugacy classes* ordered by subconjugacy is what indexes orbit /
   isotropy types (Bredon, *Compact Transformation Groups*, Ch. I §5) — i.e. the
   same object as T14b's stratification.

They genuinely differ. `O ≅ S₄` has element orders 1,2,3,4 only:

| relation | literal (z-axis setting) | subconjugate |
|---|---|---|
| `C₁, C₂, C₄ ≤ O_h` | ✓ | ✓ |
| `C₃ ≤ O_h` | ✗ (the cube's 3-fold is the **body diagonal**) | ✓ |
| `C₆ ≤ O_h` | ✗ | ✗ (no order-6 proper rotation in `S₄`) |
| `D_2h, D_4h ≤ O_h` | ✓ | ✓ |
| `D_3h ≤ O_h` | ✗ | ✗ (`O_h` carries `D_3d`, not `D_3h`) |

**Which one the selection gate needs — and why it is not a preference.** A rule and
its geometry are used in **one frame**: the ordinates stream against a mesh whose
axes the problem fixes. So the gate asks **literal containment**, which is also the
only sense `is_invariant` can test (it applies z-oriented operators to actual
nodes). Making `contains` literal is what makes `contains` and `is_invariant`
answer the **same question** — they previously did not, which is how the lattice
and the checker could disagree with nothing noticing.

The re-orientation freedom is real (rotating every ordinate yields a valid rule)
and is **exactly** subconjugacy. It must be explicit, not folded into `contains`.
Per L-013/T14, its honest shape returns the **aligning rotation**, not a bool:

> `contains(H) -> bool` — literal, same setting. `align_to(K) -> Rotation | None`
> — the subconjugacy **certificate**, so a caller that re-orients a rule gets the
> rotation it must apply and selection stays honest about having done so.

`align_to` is NOT yet built — recorded here as the named next step.

### T15b — ⭐ The lattice is COMPUTED, not tabulated (and that closed a deliberate gap)

`_contains` hand-tabulated its relations, and the source admitted the hole:
*"Cn ⊂ Oh would need n ∈ {1,2,3,4,6}; we do not encode this until a consumer
needs it."* `[M]` That note is itself wrong under **either** convention — 6 is
never right. The gap made `O_h.contains(X)` return **False for every** `C_n` and
`D_nh`, including `C₁`, `C₂`, `C₄`, `D_2h`, `D_4h`.

Fixed by **computing** finite-vs-finite containment from the SAME operator sets
`_check_invariance_3d` applies (`_realized_ops` → `_close_group` → matrix-set
containment, memoised). A hand-maintained table is a claim with no construction
behind it — the T12b lesson applied to the symmetry lattice. `[M]` The computed
group orders are textbook: `|Z₂|=2, |C₃|=3, |C₄|=4, |D_2h|=8, |D_3h|=12,
|D_4h|=16, |O_h|=48, |I_h|=120`, and the lattice now satisfies reflexivity,
antisymmetry, transitivity and **Lagrange** (a subgroup's order divides the
group's) — laws a computed lattice can be checked against and a tabulated one
cannot.

Continuous groups keep the static lattice: they have no finite realization to
compare against.

### T15c — ⛔ `_vertical_mirrors` placed the PLANES at the normals' angles

`[M]` A fourth live defect, found while verifying the derivation and **missed by
the generator audit's clean bill**. `_vertical_mirrors(n)` put the k-th mirror
**normal** at `kπ/n`, hence every **plane** at `kπ/n + π/2`. The standard `D_nh`
setting (principal axis z, vertex on x — which is exactly what
`linspace(0, 2π, n, endpoint=False)` builds) requires the planes AT `kπ/n`.

For **even** `n`, `π/2` is a multiple of `π/n`, so the shift maps the normal set
onto itself and is invisible. For **odd** `n` it is a different set of planes:

```
  product(4,n) closed under the phi=0 mirror : True at n = 3,4,5,6,7,8
  Dnh(n).is_invariant(product(4,n))  BEFORE  : False at n = 3,5,7
                                      AFTER  : True  at every n
```

⚠ **Why the audit could not see it:** orthogonality, determinant, group closure
and group ORDER are all preserved by a *rotated* mirror set. Only comparing the
angles against the setting can catch it. Generalises: **a generator audit that
checks algebraic properties cannot see a CONVENTION error.**

### T16 — ⭐ The Stage-1 gate was fed the SPENT half of the symmetry group

> **✅ LANDED `e7d44f3c`** (2026-08-02). Read T16b immediately after this section:
> executing the derivation changed four of its conclusions, and the ⚠ below is the
> reason — every number here was hand-derived by an agent with no shell.

The derived answer to *"what replaces `G_geom ⊆ G_rule` once the checker is
truthful?"* ⚠ **Derived, not executed** — the deriving agent had no shell, so
every number in `scratch/q2_stage1_predicate_derivation.md` is hand-derived or
cited. The three prerequisites below WERE verified against the tree (T15b/T15c).

**The chain that made the re-pose necessary.** `GEOMETRY_GROUPS["cylinder"] = SO2`
is right *as geometry*; `registry.py:438` declares `product → invariance_group =
SO2`, which is **false** (`[M]` its node set is `D_{n_φ h}`: `C_{n_φ}` ✓, `σ_v` ✓,
`σ_h` ✓, `C_{2n_φ}` ✗); `registry.py:677` gates on
`geom_group.is_subgroup_of(spec.invariance_group)`, passing **only** because of
that false tag; and `is_invariant(SO2)` — the checker that should have caught it —
was broken in exactly the certifying direction (ERR-072).

**The structural root:** no finite point set on S² is SO(2)-closed, so Stage 1 as
posed is **unsatisfiable by any discrete azimuthal rule**. It could only ever pass
on a false tag.

**The derived resolution.** Split the geometry's group by how the action is used:
`G⁰` acting freely (it REDUCES the spatial dimension — and its non-trivial fiber
action **IS the α term**, T3/T9), continuous isotropy `G⁰_{x₀}` (the angular
domain IS the quotient `S²/G⁰_{x₀}`), and the **discrete residual** `Γ = G/G⁰`,
which must be realized as an ordinate permutation. `Γ` is always finite.

> **`GEOMETRY_GROUPS` recorded the half of the symmetry group that the dimensional
> reduction had already SPENT (paid as the α term) and omitted the half still
> owed.** The relation `G ⊆ Sym(Q)` was never wrong — it was fed the wrong group.

Proposed: `slab, sphere → Z₂`; `cylinder, cartesian2d → D_2h` (+ a Stage-0 domain
column). `O_h` for cartesian2d over-claims by 6× — it demands `x↔z`, never a
symmetry of a z-uniform problem. `invariance_group` becomes
`invariance_group_for(params)`, mirroring `degree_of_exactness_for`; three of four
rules are constant and `product` is the parameter-dependent one — **the one that
lied, because the field's type could not spell it** (the THIRD instance of
"the type erases a distinction the prose states": T2/T12b degree, T12d
`parameters`, and now this).

`[M]` **Decisive discriminator, verified:** odd `n_φ` for cylinder. Today ADMITted;
under `D_2h` **REJECTed** — and that rejection *is* ERR-042, today a hand-written
guard. (`product(4,3)` genuinely fails `σ_x`: `φ→180°−φ` sends the 0° node to
180°, which is not a node.) Even `n_φ` stays admitted.

⚠ **It does NOT catch #326, and must not** — the rule *is* σ_y-invariant, so
Stage 1 passes and should. ⛔ **The reason given here was wrong**: this said "#326
is a sweep/ordering defect". `[M]` T22b establishes it is a **double-cover**
defect, and the fold is a *quadrature-layer* operation, so the fix does live at
this layer — just not in this conjunct. The predicate that sees it is not
"is the rule σ_y-invariant?" (yes, trivially) but **"has the owed mirror been
REALIZED as a quotient, or only as a permutation index?"** — Q5.

**New for Q5:** the azimuthal **offset is exactness-invisible but `Sym`-visible** —
`φ_m = δ + 2πm/n` puts mirror planes at `δ + kπ/n`, so `D_2h ⊆ Sym` needs `n` even
**and** `δ ≡ 0 mod π/n`. Stage 1 is the only conjunct **in the selection gate**
that can see the offset. ⛔ **But it is NOT the only thing that can** — and the
sharpening matters, because this paragraph was used to conclude the offset is a
free choice. `[M]` **Σ sees it decisively**: `δ = 0` ⟹ `|Σ| = 2` per level,
`δ = π/n` ⟹ `|Σ| = 0`, which is exactly the free-vs-non-free fold. See **T25**:
the offset is **derived**, not chosen.

### T16b — ⭐⭐ What EXECUTING T16 changed (and the general lesson)

T16 was derived without a shell. Four of its conclusions moved when measured, and
the direction of the movement is the lesson: **a derivation predicts the defect it
was looking for; measurement finds the ones it was not.**

**1. T16 predicted ONE behaviour change. `[M]` There are FIVE, and the one it
predicted is not the load-bearing one.** Across target degrees 2–9:

| change | conjunct | why it is a correction |
|---|---|---|
| slab/sphere **reject** `Product` | domain | an S² rule handed to a 1-D solver |
| cylinder **rejects** `GaussLegendre1D` | domain | a μ-marginal cannot carry 2 angular DOF |
| cylinder **admits** `Lebedev`/`LS_N` | symmetry | the `SO(2) ⊄ O_h` impossibility |
| cylinder **rejects** `Product` at odd `n_φ` | symmetry | ERR-042 — *the predicted one* |
| cartesian2d **admits** even-`n_φ` `Product` | symmetry | the `O_h` 6× over-claim |

**2. The DOMAIN conjunct does most of the work, and a group-only gate has almost
none.** `[M]` At `d=5` a group-only gate admits **every rule for every geometry**
— 16/16 — because `Z₂` is satisfied by any symmetric rule and `D_2h` by every S²
rule *and* by the 1-D marginal. T16 flagged the Stage-0 domain column
parenthetically; it is in fact the conjunct that separates the geometries.

**3. The root cause is one field doing two jobs.** `invariance_group` was
carrying *domain* and *symmetry* at once — `SO2` meant "1-D-ish", `O_h` meant
"on S²" — so a wrong domain answer and a wrong symmetry answer **cancelled**, and
the selector looked healthy. That is why the false tag survived: it was load-
bearing for a job it was not named for. Same family as the twin-errors-cancel
pattern; the tell is a field whose *values* correlate with something other than
its name.

**4. `invariance_group_for(params)` was the wrong fix, and so was the lattice.**
T17 already superseded the first. Measurement killed a second candidate: gating
as `residual.is_subgroup_of(measure.invariance_group)` **rejects Gauss-Legendre
for a slab**, because GL's tag `SO(2)` is *true but not maximal* (it names the
group the domain was quotiented BY) and `Z₂` is a reflection, hence not in the
rotation group `SO(2)`. ⟹ **the gate must ask the NODES, not the tag.** A
declaration is allowed to be true-and-not-maximal; only a *computed* answer can
be gated on.

**5. A reasoning error worth keeping.** `rules_product.py` defended its `SO(2)`
tag as *"a conservative upper bound"*. For an **invariance** claim a larger group
is a **STRONGER** claim, so an upper bound is an over-claim — never conservative.
(Contrast `degree_of_exactness`, where `min` genuinely is conservative: a lower
degree is a weaker claim.) The same paragraph *also* undersold the truth, calling
the rule "`C_{n_φ}`-invariant strictly" — an index-4 subgroup of the real answer.
**Check which direction "conservative" runs for the specific claim type**; it is
not a property of the word.

**6. The reason Q2 stopped here was itself wrong.** The checkpoint warned this was
"a live behaviour change in rule selection" needing a ruling. `[M]` Nexus:
`select_quadrature` has **15 callers, all in `tests/numerics/test_registry.py`** —
**zero production**. Production builds quadratures directly
(`Quadrature.level_symmetric(order)`, `Quadrature.lebedev(...)`). The selector is
machinery-in-waiting, so the re-pose changed test expectations only. **Measure the
consumer set before escalating a change as user-facing.**

### T18b — ⛔ The ε retirement is MIS-SCOPED: two of the three sites do not want `Σ`

T18 concluded "pass `Dnh(2)` explicitly, not the full group". `[M]` **That half is
confirmed exactly** — `Σ(D_2h)` reproduces the octant-zero set node-for-node on
every rule tested (`product(4,4)` 16/16, `product(4,8)` 16/16,
`level_symmetric(8)` 0/0, `lebedev(11)` 18/18). But measuring the *consumers*
(not just the answer) shows the retirement as scoped would **lose information at
one site and widen the answer at another**. L30, one level deeper than T18 saw.

| site | what it actually asks | vs `Σ(D_2h)` |
|---|---|---|
| `_OCTANT_SIGN_EPS` | the **chamber LABEL** — a per-component sign vector driving `partition_by` | `Σ` is only its degenerate stratum |
| `_MU_DIRECTION_EPS` | `Σ` under the **single mirror** `⟨σ_x⟩` (sweep-leg assignment) | strictly SMALLER — 8 vs 16/18 |
| `TANGENTIAL_EPS` | the boundary's Γ₊/Γ₋/tangential split on a **face normal** | unrelated (T18 already carved it out) |

**Why `_OCTANT_SIGN_EPS` cannot become `singular_set`.** `[M]` The predicate
emits **8–26 distinct sign vectors**, not 2. It is the T14b **covector**, and the
octants are the full-dimensional chambers; `Σ` is the union of the lower-
dimensional cells. A boolean membership set cannot drive a partition into
chambers. The primitive to extract here is the **chamber label**, not `Σ` — and
that is Q5/Q7 business (the chamber/RANGE factorization), not a swap.

**Why `_MU_DIRECTION_EPS` must not be folded onto it.** `[M]` `{|μ_x| < ε}` is a
strict subset of `Σ(D_2h)` (8 vs 16 on `product(4,8)`, 8 vs 18 on `lebedev(11)`).
Folding it would mislabel 8–10 ordinates as radially degenerate and put them on
NO sweep leg. Its group is `⟨σ_x⟩`, rank one.

**Bonus finding.** `product(4,4)` has **ZERO full-dimensional octants** — all 16
nodes lie on walls, because `n_φ=4` puts every azimuth on a coordinate plane.
Any consumer assuming "most ordinates are interior to an octant" is wrong for
that rule.

⟹ **The item is NOT a mechanical swap and is deliberately left open.** It
re-enters as part of the chamber work, with `Dnh(2)` and `⟨σ_x⟩` as the two named
groups. The general rule this sharpens: *measuring what a detector SELECTS is not
enough — measure what its CONSUMER does with the answer.* A set and a labelling
can agree on membership and still not be interchangeable.

### T19 — ⭐⭐ IMPOSE a structural property; never inherit it to round-off

User instruction, 2026-08-02: *"go into the source of the numpy or scipy function
and check their implementation."* It paid, and the lesson is general.

Both `numpy.polynomial.legendre.leggauss` and scipy's
`_gen_roots_and_weights` end the SAME way, for a measure whose weight is even:

```
x = (x - x[::-1]) / 2 ;  w = (w + w[::-1]) / 2     # symmetry
w *= mu0 / w.sum()                                 # zeroth moment
```

They do not hope the rule comes out symmetric — they **make** it symmetric. `[M]`
Adopting the two lines in the generic body: symmetry defect `8.6e-16 → 0.0`
(exact, every family, every `n`); mass defect `4.4e-16 → 0.0`; and exactness
*improves* `7.9e-16 → 3.3e-16`, because the mirror average cancels the
antisymmetric part of the eigensolver's error. The banded solver and the Newton
step, also in both references, buy nothing measurable here (`6.1e-16 → 6.1e-16`)
and were NOT adopted — read the source, then measure which parts earn their place.

> **Three tiers of exactness, in preference order.** (1) **Generate by group
> action** when the operators are IEEE-exact — sign flips and coordinate swaps.
> Lebedev does this: `[M]` 48/48 operators exact. (2) **Use the closed-form index
> relation** when a constraint determines the value. (3) **Impose by averaging**
> when neither applies. Gauss rules are tier 3; `level_symmetric` was doing
> arithmetic where tier 2 applied (T20).

**And derive the condition, don't flag it.** scipy passes `symmetrize` as a
hand-set boolean. But `α_k ≡ 0 ⟺ w` is even is a *theorem* (α_k is the first
moment of an even function against an odd integrand), so `is_symmetric` reads it
off the recurrence. `[M]` It agrees with scipy's flag on every family including
the parameterised ones — `jacobi(a,b)` is symmetric exactly when `a == b`, which
the derivation gets right without being told and a declared flag could not.

### T20 — ⛔ Recomputing what a closed-form index determines (a RECURRING pattern)

`level_symmetric_sn` recovered its third direction cosine as
`sqrt(1 - mu_z² - eta²)`. But Carlson's construction fixes
`p + k + j = N/2 - 1`, so the third cosine is `mu_levels[j]` — an **index**.
`[M]` The square root landed within ~1e-16 of the level value but not ON it, so
at `N=16` the y-axis carried **22 distinct magnitudes where the level array has
8**, and only **8 of 48** `O_h` operators were exact — precisely the pure sign
flips (negation is exact in IEEE; a coordinate swap only if the values match).
**The rule advertised `O_h` and realized `D_2h`** — the same over-claim shape as
T16's geometry table, one layer down. Fixed: 48/48 at every order, and the
`xi_sq < -1e-14` guard retired because the quantity it guarded stopped existing.

The pattern recurs: MoC's `_reflected_azi_index` is an `argmin` computing exactly
`n_azi-1-a`; `_find_reflections` is an `argmin` with **no distance guard at
all**, load-bearing for every reflective/albedo BC. Inventory:
`scratch/q3_symmetry_exactness_inventory.md`.

### T21 — ⭐⭐ An ordering by a NON-INJECTIVE key is not an ordering

User correction, 2026-08-02, on being asked to rule on a sort-key tie-break:
*"Do we have proper level machinery? If we did, I feel like this problem would
disappear. Seems like an 'architecture incomplete' problem, rather than a
decision problem."* Correct, and measurement settled it.

`LevelStructure` stored `level_indices` **already sorted by**
`η = sinθ cosφ`. That key is EVEN in φ, hence **2-to-1 on a circle**. `[M]` On
one level of `product(2,8)`, 8 ordinates give only **5 distinct η**, and the
stored order does not determine `sign(ξ)` — that information is not in the key.

> **So the "ties" are not ties.** They are the fibers of a non-injective map, and
> no choice of tie-break makes a 2-to-1 projection into an ordering of the
> circle. There was nothing to rule on. **This IS #326** — "the cylindrical level
> double-covers the azimuthal range" — and the double cover is `cos φ`.

The type stored **the RESULT of an ordering decision, and a lossy one**, instead
of the structure that determines the ordering — build-the-product-not-the-
primitive. Q4 gives it the invariant (`PolarInvariant.SIGNED_MU_Z` vs
`ABS_MU_Z` — T7's two producers, finally distinguishable) and the fiber's own
coordinate.

**The coordinate takes BOTH `azimuth` and `hemisphere`, and writing it proved
why**: under `ABS_MU_Z` a level carries two circles, so φ alone repeats — the
type failed its own injectivity test on the first draft. `[M]` The pair is
injective on every level under either invariant.

⚠ The corrected order is NAMED (`fiber(level)`), not swapped in: the stored order
is what the cylindrical sweep consumes and changing it moves results. That swap
is #326's business, with its own gate.

**Meta-lesson, and the sharpest one of the session: round-off can DISGUISE a
structural degeneracy by manufacturing fake distinctions.** With inexact nodes,
1e-16 noise broke the ties and made a 2-to-1 map look 1-to-1. Making the nodes
exact collapsed the disguise. Exactness is not only an end — it is a *diagnostic*.

### T22 — ⛔ The ordering is CLOSED as an adjudication: every option is measured

T21 ended by deferring "the swap" to #326 "with its own gate", and the cold-pickup
anchor then promoted that to **the** next step. `[M]` **Measured 2026-08-02** on the
cylindrical fixed-source problem the three xfail rows use (`product(4,8)`, nx=20,
isotropic source, `inner_tol=1e-13`), varying ONLY the per-level ordering through
`tests/sn/_test_helpers.product_level_ordering` with node VALUES held identical:

```
ordering      exact-nodes   homogeneous       heterogeneous
PRODUCTION    (n/a)         1.190877e-01      5.144789e-02
lexsort       False         1.190877e-01      5.144789e-02
lexsort       True          1.190877e-01      5.144789e-02
stable        False         1.190877e-01      5.144789e-02
stable        True          1.190877e-01      5.144789e-02
azimuthal     False         nan               nan
azimuthal     True          nan               nan
```

**⚠ Do NOT read row 6/7 as "the fiber ordering is wrong".** That was this
theorem's first draft and it is a misreading — the fiber ordering is the
*physically* natural march order. User challenge, 2026-08-02: *"Don't simply
accept the overflow. Investigate the mechanism and reason if the mechanism is a
principled and inherent mechanism or a problem of implementation."* It is the
latter, and running that down is what produced T22b.

1. **η-ascending orderings are bit-identical to each other.** The tie-break is not
   a free variable that moves the answer — not even with exact nodes. (This also
   re-confirms independently that **#325 does not gate #326**.)
2. **The NaN localises to ONE expression** — the cylinder's edge convention at
   `pole_angular_closure.py:599-603`:
   ```python
   eta_edge[0] = -sin_theta                        # the march STARTS at −sinθ
   for m in range(M - 1):
       eta_edge[m + 1] = 0.5 * (eta[m] + eta[m + 1])
   eta_edge[M] = sin_theta                         # …and ENDS at +sinθ
   ```
   That hard-codes **"this level is ONE monotone arc from −sinθ to +sinθ"**. Under
   the ω-order `eta[0] = +sinθ`, so `τ_raw[0] = 2sinθ/(0.434+0.508) = 1.079` —
   the measured value, reproduced exactly from the convention alone.
3. **Nothing catches it.** `τ_raw`'s only guard is `abs(deta) > 1e-15` (exact
   collapse), and `morel_montry_tau_per_level` then applies the structural
   `max(0.5, min(1.0, ·))` clamp, which **absorbs the out-of-range value into a
   plausible one**. So a violated premise is laundered into a finite wrong answer;
   the NaN is downstream divergence, not a guard firing. `[M]` A mere *reversal*
   of a level yields `τ_raw ∈ {+1.079, −0.079}` with **no NaN and no raise** at
   all. The clamp is the thing that hides mis-ordering.

### T22b — ⭐⭐ The α-march is a march in ω ARC BY ARC; the arcs end on Σ = {ξ = 0}

The principled statement the overflow was pointing at, and it makes the fold
*necessary* rather than merely elegant.

ξ is the coefficient of `∂ψ/∂ω` in the cylindrical streaming operator, so
**ξ = 0 is exactly where the angular-redistribution flux vanishes** — the only
place the recursion's closure `α_{1/2} = α_{M+1/2} = 0` is physically right, not
merely conventional. Those points are `ω ∈ {0, π}`: each level's circle carries
**exactly two**, and they cut it into two arcs. This is the set the campaign
already named in **Q1** — `Σ = {ξ = 0}`, the orbifold singular set — arriving
independently from the sweep side.

On ONE arc, `η = sinθ·cos ω` is **strictly monotone in ω**. So `[M]` (measured,
every level of `product(4,8)`):

> **the η-ascending order and the ω-order are THE SAME ORDER.** `[4 3 2 1 0]`
> both ways. They are not competitors; they are one order seen through two
> coordinates. They can only *appear* to compete on a level carrying two arcs.

**Verdict on the mechanism: implementation, not inherent — and the convention
encodes the bug.** The edge rule is *correct code for one arc* applied to a
two-arc level. Sorting by η makes a two-arc level **impersonate** one arc (a
monotone sequence from −sinθ to +sinθ) by interleaving ordinates from *different*
arcs — which is precisely why mirror partners land adjacent and split `{1, ½}`.
So neither ordering is right, because the LEVEL is wrong: the η-order returns a
finite wrong answer, the ω-order a divergent one, and they are the same disease.

**The strongest result: the fold needs NO change to the closure at all.** `[M]`
The unmodified production producer, run on a half-range level:

```
                            tau_raw                          degenerate
production (2 arcs)  [0, 1, 0, 1, 0, 1, 0, 1]                 8/8
half-range (1 arc)   [0, 0.292893, 0.5, 0.707107, 1]          2/5
```
`0.292893 = 1 − 1/√2`, `0.707107 = 1/√2` — smooth, monotone, in range, identical
on all four levels. **The alternation is gone.** The two remaining degenerate entries are
the arc's own endpoints (`η = ∓sinθ` *is* the boundary, by definition, and they
are the Σ points), and only one of the five is touched by the clamp. `[M]` Mass
is exact: `12.566370614359172` before and after `= 4π` bit-for-bit, `32 → 20`
ordinates, `|Σ| = 8`, with T5's orbit-stabilizer trapezoid `[w, 2w, …, 2w, w]`
supplying the weights (Σ points have `|Stab| = 2`).

⟹ `morel_montry_tau_raw_per_level` was written for a half range all along. The
defect is that nobody ever gave it one.

⚠ The lesson for planning, not just for #326: **an anchor written at a compaction
checkpoint is a claim, and it can be falsified by the repo it summarises.** This
one contradicted its own §7 ("do not fix #326 directly") and the promoted test
file's own docstring ("WHY NO RE-ORDERING FIXES IT"), and survived a checkpoint
because prose was reconciled against prose. Cost: one probe. Reconcile a
next-step claim against a MEASUREMENT before it becomes the anchor.

### T23 — ⛔⭐ The declared acceptance test CANNOT adjudicate the fold — and the reason is a live Pattern-4 gap

§7 declared the three `xfail(strict=True)` rows the campaign's proof: "when the
machinery is right they XPASS and force their own deletion". They measure
`max|ψ − ψ[reflection_index("y")]| / max|ψ|`.

**The fold makes that expression meaningless**, because the quotient carries ONE
representative per mirror pair — the partner is *not in the node set*. `[M]`
Measured on a half-range restriction (ξ ≥ 0) of `product(4,8)`, 32 → 20 ordinates:

```
                          full sphere      half-range quotient
involution?               True             False
identity (slab case)?     -                False
max |ξ[partner] + ξ|      2.303e-16        9.404e-01     <- maximally wrong
raised?                   -                NO
partner map               -                [0 0 0 4 4 5 5 0 9 9 10 10 ...]
```

Every node is mapped onto some **ξ = 0** node — the nearest thing to a reflection
that exists in the set — silently, many-to-one, non-involutive.

**Root cause, and it is a bug independent of #326.**
`directional._find_reflections` (`directional.py:125`) is a bare `argmin` over
squared distance with **no distance threshold and no closure check**. It *asserts*
"these are the reflection partners"; it never verifies it. Today every production
rule happens to be closed under all three axes, so the map is correct — a latent
trap, not a live defect. **The fold breaks exactly that unstated precondition**,
so the guard is a MUST-PRECEDE cleanup (`coding-standards.md` "clean before
extending"), not follow-up polish: fold first and the cylinder ships a silently
wrong pole map.

The fix is already in the tree and is one call: `[M]`
`SubgroupOfO3.Dnh(2).is_invariant(measure)` returns **True** for the full rule and
**False** for the half range — the campaign's own computed predicate (T14/T17)
gets the right answer where the nearest-neighbour matcher gets a silent wrong one.
Same shape as ERR-072/073: *a certification claimed and never computed.*

⛔ **THE MUST-PRECEDE IS DISCHARGED — TWICE OVER (annotated 2026-08-07).**
Q5.0.1 certified the table (`a7695148`, catching the live odd-`n_φ` defect —
ERR-074); then track G RETIRED the table entirely (§7d.3) onto
`Quadrature.ordinate_permutation(motion)` — derivation-at-need through the
G2-verified `RigidMotion.preserves`, carrying the same three certifications
and returning `None` for a non-symmetry, with each caller owning its refusal
("no specular pairing" at certification tier, "cannot seed the r = 0 pole"
at sweep tier). `_find_reflections`, `reflection_index` and the
`reflection_partners` field NO LONGER EXIST. T23's conclusion is UNCHANGED:
the three xfail rows cannot adjudicate the fold, and the #229 floor is the
acceptance.

**The replacement acceptance test.** The non-vacuous prediction is already on the
issue: a half-range level should **remove the #229 azimuthal error floor** and
restore a convergence order to the anisotropic curvilinear MMS gate (today flat at
≈1.9e-2). That is a manufactured-solution measurement against a structurally
independent reference — not a self-comparison, and not something the fold can
satisfy by construction. Secondary, both substantive: `τ_raw` must leave `{0,1}`
(T22's mechanism), and the α closed form of **T3** must hold exactly.

⚠ **Do NOT** re-pose the three xfail rows as "ψ even in ξ on the quotient". That
IS satisfied by construction and would be a gate that cannot red — `vv-principles`
Mode 8, the tautological-guard class.

### T24 — ⭐⭐ `Σ = ∅` is the FOLD'S WELL-POSEDNESS CONDITION, and it is a property of the RULE

User direction, 2026-08-02, declining a menu of implementation sites: *"this is a
question of principles and design … grounded in mathematics, invariants and
theorems … we should consider the possibility of leveraging the singular set
(since the quadrature already knows it). When considering this possibility, we
should consider other quadratures that could be used in the same problem … and if
they can also define, in a mathematically clean way, their own singular set."*

They can, and the answer selects the design. Σ is already a *computed* object
(`symmetry.singular_set`, T14: an integer identity `π_M(i) == i`, unrepresentable
without the invariance proof). Asking it of each candidate cylindrical rule under
the owed mirror `σ_y` — `[M]`, all four levels of a GL(4) polar factor:

| azimuthal rule | \|Σ\| per level | why |
|---|---|---|
| equispaced **left**, `n_φ` even | **2** | nodes land on `ω = 0` AND `π` |
| equispaced **left**, `n_φ` odd | **1** | only `ω = 0` is a node |
| equispaced **midpoint** | **0** | nodes straddle both |
| `level_symmetric(4)` | **0** | no ordinate has `ξ = 0` |

**`Σ = ∅` ⟺ the mirror acts FREELY ⟺ every orbit has length exactly `|G|`.** That
is the condition under which the quotient has no fixed points to special-case: the
orbit-stabilizer weight `W = w·|G|/|Stab|` collapses to a uniform `2w`, so mass is
preserved **bit-exactly** rather than by a mixed `w`/`2w` sum. `[M]` folded:

```
                 ordinates    mass                      tau_raw           degen
equi-LEFT(8)     32 -> 20     12.566370614359169 (3ulp) [0,.293,.5,.707,1]  2/5
equi-MID(8)      32 -> 16     12.566370614359172 = 4pi  [.2195,.4142,       0/4
                                                         .5858,.7805]
```
`0.4142 = √2−1`, `0.5858 = 2−√2` — strictly interior, symmetric about ½.

**This unifies a decision §4 already made on other grounds.** Gauss–Lobatto was
characterised and *declined* there because "`τ_raw=0` breaks the M-M recurrence AND
the R12a predicate" — i.e. because Lobatto **puts nodes on the interval endpoints,
which are Σ**. Equispaced-left is doing the very thing Lobatto was rejected for,
on the circle. One criterion, stated once: **a rule is admissible for a fold iff it
places no node on Σ.**

⚠ **A hypothesis this REFUTED, recorded so it is not re-derived.** The natural
guess — "the `τ_raw ∈ {0,1}` degeneracy IS Σ-occupancy" — is **FALSE**. `[M]`
`equi-MID(8)` has `|Σ| = 0` and is *still* **8/8 degenerate** unfolded. The
degeneracy comes from the **double cover** (mirror partners share `η` whether or
not anything sits on Σ), so the FOLD is what cures it, in every rule; Σ decides
only the folded arc's **endpoint** behaviour. Two independent structural facts,
one float.

### T25 — ⭐ The azimuthal OFFSET is exactness-invisible, `Sym`-visible, and Σ-DECISIVE — so it is DERIVED, not chosen

`[M]` **Gauss–Chebyshev-1 in `x = cos φ` IS the midpoint rule on `(0, π)`**:
`x_k = cos((2k−1)π/2n)` ⟺ `φ_k = (2k−1)π/2n`, agreeing to **5.55e-17** at `n = 4`
and `n = 6`. So the azimuthal rule that natively lives on the half range is a
*Gauss* rule — `CHEBYSHEV_T`, which **Q3's `GeneratingMeasure` already
constructs**. The campaign built the object two steps before it was needed.

This closes the offset question the plan had left open in three places, and it
closes it by **derivation rather than by vote** (the campaign's own standing
ruling: *when a fork looks like a preference, the abstraction is missing*). The
three visibility levels, now complete:

* **exactness — BLIND.** The periodic trapezoid is exact for trig degree `< n` at
  any offset (T1). `§4` already recorded this.
* **`Sym` — VISIBLE.** `φ_m = δ + 2πm/n` puts the mirror planes at `δ + kπ/n`, so
  `D_2h ⊆ Sym` needs `n` even **and** `δ ≡ 0 mod π/n`. Recorded at T16.
* **Σ — DECISIVE, and new.** `δ = 0` ⟹ `|Σ| = 2` per level; `δ = π/n` ⟹
  `|Σ| = 0`. Nothing else in the machinery separates them, and it is the
  separation that decides whether the fold is free.

⟹ the offset is not a knob on the RANGE axis. **It is fixed by requiring the
quotient to be free**, which is the only reading under which "half-range" is a
quotient (T5) rather than a restriction.

### T26 — ⛔ The R12a carrying predicate decides a STRUCTURAL question with a 5.55e-16 float

`SNMesh.radial_characteristic_levels` — which decides whether a mesh carries
independent ψ½ state, and therefore its **field block structure and coupled
operator grid** — is `τ_raw[0] ∈ (0,1)` exclusive, documented as a "bit-exact
trichotomy". `[M]` at full precision:

```
equi-LEFT(4)      0.0                    1 - t0 = 1.000e+00   -> not carrying
equi-LEFT(5)      0.9999999999999994     1 - t0 = 5.551e-16   -> CARRYING
equi-LEFT(7)      0.9999999999999989     1 - t0 = 1.110e-15   -> CARRYING
equi-LEFT(8)      0.0                    1 - t0 = 1.000e+00   -> not carrying
level_symmetric   1.0                    1 - t0 = 0.000e+00   -> not carrying
```

For **odd** `n_φ` the exact value **is** 1 — the level's two most-negative `η`
coincide, so `eta_edge[1]` lands on them — and trig round-off alone puts it
`5.6e-16` below, flipping the classification. The trichotomy is bit-exact only for
even `n_φ`; the docstring's claim is true of the cases anyone ran.

**And the predicate conflates two different structural facts**, which is why one
float cannot carry it:

* `τ_raw[0] = 0` ⟺ a node sits at the march-start endpoint — **Σ occupancy**;
* `τ_raw[0] = 1` ⟺ the first two ordinates share `η` — **the double cover**;
* `τ_raw[0] ∈ (0,1)` ⟺ neither.

Both degenerate cases read "not carrying" today, for unrelated reasons. This is
exactly the disease **T14** already diagnosed and cured for Σ — *not a
floating-point question, an integer identity* — recurring one layer up, in the
sweep. The certificate answers both conjuncts exactly.

### T27 — ⭐⭐ The `[½,1]` clamp is TWO objects fused; retire one, promote the other

User instruction, 2026-08-02: *"check if twisting its perspective, function and
location adds some value. If we have reached a way to structurally avoiding the
clamp, then retiring it is the right thing. If we have reached a mechanism to
structurally make the clamp unnecessary, then a test must assert this condition in
the configuration that would be most likely to activate the clamp and show that it
passes because our structural mechanism is working."*

Twisting it does add value, and the measurement **corrects a hint this plan
carried**: I had suggested the fold puts `τ_raw` "in the sphere's unclamped
regime". `[M]` It does not, and the sphere is not in `[½,1]` either:

```
folded, n_phi =  4 .. 64    tau_raw  [0.292893, 0.707107] -> [0.200289, 0.799711]
sphere GL(8/16/32)          tau_raw  [0.390354, 0.610160]   -- UNCLAMPED, also
                                                               outside [1/2, 1]
```

So the clamp still *bites* after the fold. Three findings settle it anyway.

**1. The `[½,1]` box was never the invariant.** The sphere runs **unclamped** at
`[0.390, 0.610]` — outside the box, and correct. The clamp exists for exactly one
thing: `τ = 0` makes the recurrence `(ψ − (1−τ)ψ)/τ` divide by zero.

**2. The fold bounds `τ` away from the singularity, ANALYTICALLY.** With midpoint
nodes on the arc the smallest-`η` node sits at `φ = π − ε`, `ε = π/2n`, so
`η_0 = −1 + ε²/2`, `η_1 = −1 + 9ε²/2`, `η_{1/2} = −1 + 5ε²/2`, giving

> `τ_0 = (ε²/2)/(5ε²/2) = 1/5` exactly, and by reversal `τ_{M−1} → 4/5`.

`[M]` approached monotonically from inside (`0.2929 → 0.2003` over `n_φ = 4…64`).
`τ = 0` becomes **structurally unreachable**, so the clamp's reason is gone even
though its range is still violated. This is the user's first branch: *structurally
avoided ⟹ retire.*

**3. Keeping it would now be ACTIVELY HARMFUL.** On a reversal-symmetric arc the
node's fractional position in its cell flips under `ω → π−ω`, so
`τ_m + τ_{M−1−m} = 1` is an exact structural identity. `[M]`

```
n_phi= 8   RAW max|tau + tau_rev - 1| = 0.000e+00    CLAMPED = 2.805e-01
n_phi=16                              6.661e-16                2.953e-01
n_phi=32                              8.216e-15                2.988e-01
   raw     : [0.219545, 0.414214, 0.585786, 0.780455]
   clamped : [0.5,      0.5,      0.585786, 0.780455]
```

> `[M]` **Configuration note added at Q5.5 (2026-08-07):** the residuals
> above were measured on this T-item's HAND-BUILT fold (pre-Q5.3 — no
> production fold machinery existed). On the LANDED fold
> (`LevelStructure.quotient`, selection descent, bit-copied charts) the
> residual is **0.0 bit-exact at every n_φ ∈ {8..64}**, so the shipped
> gate asserts the identity with NO epsilon. Same conclusion, stronger
> — the design choice at Q5.3 is what upgraded it.

The clamp is **asymmetric** (`[½,1]`, not `[0,1]`), so it destroys an identity the
fold creates — **a symmetry-breaking defect of the same kind as #326, one level
down.** Retiring it is not merely permitted; keeping it would reintroduce the
disease the campaign is curing.

**The twist — the fused expression is TWO objects at TWO locations.**
`max(0.5, min(1.0, τ_raw))` welds an *absorber* to a *range check*:

| half | verdict | why |
|---|---|---|
| `[½,1]` **absorption** | **RETIRE** | reason structurally removed (2); actively breaks an exact identity (3) |
| `[0,1]` **membership** | **PROMOTE to a guard that RAISES** | `τ_raw ∉ [0,1]` means a node lies **outside its own angular cell** — impossible on a monotone arc, so it certifies arc well-posedness |

Today the absorber *hides* the range violation (T22), which is why a mis-ordering
laundered into a finite wrong answer instead of stopping. `[M]` the promoted
guard's honest teeth:

```
production: full circle, eta-sorted (THE #326 DEFECT)  [0.0000, 1.0000]  MISSES
full circle, omega-ordered (T22's NaN case)            [0.2929, 1.0790]  CATCHES
folded + midpoint (the Q5 target)                      [0.2195, 0.7805]  silent
```

⚠ **State the limit plainly: the guard does NOT catch the double cover** — the
`[0,1,0,1,…]` fingerprint is entirely *inside* `[0,1]`. It catches the ordering
violation (and would have stopped T22's NaN at source instead of 400 lines
downstream). The #326 detector is the `Σ`/fold criterion (T24), not this.

**The gate the user's instruction requires** — it must pass *because* the
mechanism works, at the most-activating configuration:

* **config**: folded + midpoint, swept to the largest supported `n_φ` (min
  `τ_raw` falls monotonically toward `1/5`, so large `n_φ` is the worst case;
  `[M]` `τ_raw` is `n_μ`-independent — `sinθ` cancels in the ratio);
* **assert the MECHANISM**: `Σ = ∅`, computed via `singular_set` (not declared);
* **assert the CONSEQUENCE**: `τ_raw ⊂ [1/5, 4/5]`, strictly bounded away from
  `{0, 1}`, and the reversal identity `τ_m + τ_{M−1−m} = 1` to round-off;
* **reddenable, on one mutation**: revert Q5.2's offset to `δ = 0` and *both*
  legs fail together — `Σ` becomes non-empty AND `τ_raw` hits `0`. That is what
  attributes the pass to the mechanism rather than to luck.

### T28 — ⭐⭐⭐ The product rule is a WELDED PAIR that re-derives its own properties; it should be a COMBINATOR over registered rules

User diagnosis, 2026-08-02, on being shown the odd-`n_φ` reflection defect:
*"the most likely reason why we have a problem in the product rule is … machinery.
(1) is the product rule hard-coded to a product pair? If so, we must become capable
of selecting the pair we want to instantiate. (2) does it derive its properties in a
principled way based on the properties of the rules we're composing + the invariants
and theorems of the product rule? If the rules we can make a product off are
registered, the product rule should be capable of using a registered rule to make a
product, and all intrinsic properties come downstream of this."*

Both halves confirmed by measurement. This is the root T24–T27 are symptoms of.

**(1) It is a welded pair.** `product_mu_phi(n_mu, n_phi)` hard-codes
`gauss_legendre_on_mu(n_mu)` for the polar factor and
`np.linspace(0, 2π, n_φ, endpoint=False)` for the azimuthal one. The azimuthal
factor **is not a measure object at all** — it is a raw array plus a scalar
weight, so it can carry no degree, no group, no generating measure, and cannot be
substituted.

⛔ **And the rule the fold needs is ALREADY IN THE TREE, registered and unwired.**
`measure.equispaced(a, b, n)` is the **midpoint** rule, and its own docstring says
why it is midpoints:

> *"the project's existing code uses left-endpoints, but this primitive offers
> midpoints because they integrate constants exactly **while preserving symmetry
> under reflection through the centre of the interval**"*

That is **T25's `Σ = ∅` criterion**, written down by whoever built the primitive
and then never reached, because the welded pair left no seam to reach it through.
`[M]` `equispaced(0, 2π, 8)` gives exactly the `π/n`-offset nodes T24 measured as
`|Σ| = 0`. §4 already recorded that `equispaced` / `gauss_chebyshev` /
`gauss_laguerre` ship with **zero production consumers** — "the abstraction is
being paid for and not used". This is what the non-use cost: the double cover.

**(2) It re-derives what its factors already know.**

| property | today | what it should be |
|---|---|---|
| nodes / weights | hand-built nested loop | `(polar ⊗ azimuthal).pushforward(embedding)` |
| `degree_of_exactness` | `min(2*n_mu - 1, n_phi - 1)` **hand-written** | composed from the factors' own degrees |
| `invariance_group` | `Dnh(n_phi)` **declared** | computed / derived from the factors |
| `Σ`, level structure | — | downstream of the nodes and the polar factor |

`[M]` `polar.degree_of_exactness` **already equals** `2*n_mu − 1` at every order —
so the product recomputes a number its own factor is holding. That is the **fourth
in-tree instance of T20** (recomputing what a closed-form relation determines).

⚠ **Be precise about the group tag: it is CORRECT, not false.** `[M]` declared
`== computed` for `n_φ = 4,5,6,7,8` — `D_{n_φ h}` every time, odd `n_φ` included.
The objection is not that it lies; it is that **a declaration cannot be falsified
by construction**, and this exact module has already shipped two declarations that
*were* false (ERR-072, ERR-073) plus a partner map that was (ERR-074). T17's rule
stands on its own: a computed property cannot lie about the object it was computed
from.

**⛔ The real prerequisite, and it is not cosmetic.** A naive composition
`min(p_polar, p_azimuthal)` would give **degree 1**, because
`equispaced.degree_of_exactness = 1` is the **algebraic** degree of the midpoint
rule *on an interval* — while on the **circle** the same nodes are the periodic
trapezoid, exact to **trigonometric** degree `n − 1`. Same nodes, same weights,
two different exactness claims, and the integer cannot tell them apart. This is
§4's "**EXACTNESS SPACE — EARNED MOST**" axis and T2/T12b's gap, arriving as a
hard blocker: **the product cannot derive its degree until each factor carries the
SPACE its degree is measured in.** T12b already named the fix — the *generating
measure* IS the exactness space.

**Consequence for the plan.** Q5.2 ("the derived offset") stops being a change to
`product_mu_phi` and becomes a **selection**: the azimuthal factor is a registered
rule VALUE, and the fold's `Σ = ∅` requirement *picks* `equispaced`(midpoint) —
equivalently `CHEBYSHEV_T` in `cos φ` (T25) — over the welded left-endpoint array.
No flag, no boolean, no `half_range=True`. The campaign's standing ruling holds:
**when a fork looks like a preference, the abstraction is missing** — here the
missing abstraction is the seam that lets a registered rule be the factor.

### T17 — ⭐⭐ WALK the subgroup graph; stop declaring the symmetry group

User ruling, 2026-08-02: *"when I was studying crystallography, there was a literal
graph that we could walk to determine [the] group. That seems like essentially the
machinery we need. Why don't we implement the graph, and walk the graph?"*

The remembered object is the **Hasse diagram of the subgroup lattice** — nodes are
groups, edges are **maximal-subgroup** relations (`H ⋖ G` iff `H < G` with no `K`
strictly between). Crystallography walks it downward from high symmetry to find
the symmetry a structure actually has; International Tables Vol. A1 renders it as
the **Bärnighausen tree**. It is exactly the machinery this campaign needs, and it
**supersedes** the `invariance_group_for(params)` proposal of T16.

**Definition.** For a measure `μ` and an expressible candidate set `C`:

> `Sym_C(μ) := the MAXIMAL elements of { G ∈ C : μ is G-invariant }`

**The theorem that makes the walk sound** — *invariance is DOWNWARD-CLOSED*: if `μ`
is `G`-invariant and `H ≤ G`, then `μ` is `H`-invariant. So the invariant set is an
**order ideal**, and two prunings follow: if `μ` is **not** `G`-invariant, no
supergroup of `G` can be invariant; if `μ` **is** `G`-invariant, every subgroup is
implied and need not be tested. This is precisely the law the audit measured
**68 violations** of (S-16) — so the walk both *requires* a correct lattice and
*is a test of it*.

**Edges are DERIVED, never tabulated.** `H ⋖ G` is computed from T15b's
matrix-set containment. A hand-drawn Bärnighausen tree would be the same
"claim with no construction behind it" this campaign keeps finding.

**Two realizations — and the second is the verification instrument** (the user's
standing ≥2 rule):
- **R1 brute force** — test every candidate, take the maximal elements.
- **R2 Hasse walk** — descend from the tops with both prunings.
The invariant tested ACROSS them: `R1 ≡ R2` on every rule × candidate set. R2 is
the fast path; R1 is what proves it.

**Why this beats a declared tag** (and dissolves T16's item 3):

1. `invariance_group` stops being a **declaration** and becomes a **computed
   property of the object**. The #327 disease — a tag with no construction, which
   nothing can check — is cured on the **G axis** exactly as T12b cures it on the
   **V axis**. A declaration, if any survives, is gated by *equality* against the
   computed value.
2. `product`'s parameter-dependence needs **no** `invariance_group_for(params)`
   function at all: `D_{n_φ h}` simply falls out of the walk. The measure tells
   you; the spec cannot lie.
3. Stage 1 becomes `G_ang(geom) ⊆ Sym(Q)` with **`Sym(Q)` computed** — so the
   false `product → SO2` tag becomes structurally unspellable rather than merely
   corrected.

**Bounding the infinite families.** `C_n` / `D_{nh}` are infinite families, so `C`
is bounded by the measure: an order-`n` rotation acting on `N` nodes needs `n ≤ N`,
giving ~`2N` candidates plus the named entries — cheap at production sizes.

⚠ The walk answers about the group's **realization in the standard setting**
(T15), not up to conjugation. A rule whose symmetry axis is not `z` reports a
smaller group — correctly, since the gate compares against a geometry in the same
frame. `align_to` (§6) is what makes re-orientation expressible.

### T18 — ⭐ `Σ` is NOT a tidier spelling of the ε-detectors. It is a LARGER, correct set.

`[M]` The measurement that reshaped the retirement, taken once the certificate
existed (landed `c9c57152`). The ε-detectors test **cosine magnitudes**, so they
see only **coordinate-plane** mirrors and are structurally blind to the
**diagonal** mirrors that `O_h` and `D_nh` (n>2) also carry:

```
  rule                 Σ|D_2h   ε-answer   identical   Σ|full group
  product(4,4)             16         16       True             16
  product(4,8)             16         16       True             32
  product(4,16)            16         16       True             64
  level_symmetric(8)        0          0       True             32
  level_symmetric(16)       0          0       True             96
  lebedev(11)              18         18       True             50
  lebedev(17)              30         30       True            110
```

`level_symmetric` is the sharpest case: **no node has a zero cosine**, so the
ε-answer is EMPTY while 32 of its ordinates lie on `O_h` mirrors.

> **The retirement warrant is not "swap the spelling" — it is "NAME the
> subgroup".** The ε-detectors were asking `Σ` under **`D_2h`**: the group
> generated by the three coordinate reflections, `≅ (ℤ₂)³`, whose chambers are
> exactly the octants — **T14b's reflection group, predicted before it was
> measured**. `Σ|D_2h` reproduces the ε-answer EXACTLY on every production rule,
> which is what licenses retiring them.

⚠ So a blind swap to "the singular set" would have **changed behaviour** — silently
widening every octant/degenerate classification. The consumer-side retirement must
pass `Dnh(2)` explicitly. Still OPEN (the consumers are hot paths).

### T15d — `O2` is realized as `C_∞h`, not `O(2)` — RESOLVED: rename to `D_∞h`

`[M]` `_check_invariance_3d` built `O2` as *rotations about z* + `σ_h`, i.e.
`C_∞h`. True `O(2)` embedded in 3-D is `C_∞v` (rotations + **vertical** mirrors).
Consequence for the audit's S-7 (`D_nh ⊆ O2`, asserted true, **pinned by a
committed test** at `test_symmetry.py:216`): `D_nh` contains `C₂` about *in-plane*
axes, which lies in **neither** `C_∞h` nor `C_∞v` — so the claim was false under
both readings. It is true only for `D_∞h`, the full cylindrical group.

**RESOLVED (user ruling, 2026-08-02): rename `O2 → D_∞h` and re-realize it as the
full cylindrical group** — `C_∞` rotations + `σ_h` + the vertical mirrors + the
in-plane `C₂` axes. Under `D_∞h` the relation `D_nh ≤ D_∞h` is **true for every
n**, so the pinned assertion becomes correct rather than deleted — the test was
right about the *relation* and wrong about the *name*. `D_∞h` is also what a
cylinder actually carries, which is why the geometry table wanted it (T16).

### T14b — The octant sign label is already a named mathematical object

`_octant_sign_predicate` (`directional.py:103`) returns `{+1, 0, −1}` per axis.
The eight octants **are** the chambers of the reflection group `(ℤ₂)³` acting on
ℝ³; the three coordinate planes are its **walls**; a node with `k` zero components
lies on a codimension-`k` face of a chamber, with isotropy `(ℤ₂)^k`. The sign
vector is the **covector** of that hyperplane arrangement (oriented-matroid
vocabulary), and the collection of them is its **face lattice**.

So the label is not an ad-hoc encoding — it is the **orbit-type stratification**,
and the zero components name *which* mirrors fix the point. Two zero components is
Thurston's **corner reflector**, which is precisely the
`|μ_x| < ε ∧ |μ_y| < ε` case `directional.py:116` describes as the 2-D wavefront
short-circuit.

---

## 2. Diagnostics and gates that now exist

| gate | what it pins | state |
|---|---|---|
| `tests/sn/sweep/curvilinear/test_alpha_closed_form.py` | T3 — α ≡ its closed form, with κ; drives **production** `cylindrical_streaming` | 20 `l1` + 15 `foundation`, GREEN |
| `tests/sn/sweep/curvilinear/test_azimuthal_mirror_symmetry.py` | T4 — ψ even in ξ | 5 `l1` (**3 xfail-strict = the live defect**) + 7 `foundation` |
| `tests/sn/verification/mms/test_mms_ordering_blindness.py` | T11 — the MMS's Mode-7 declaration, made executable | 10 `foundation`, GREEN |
| `tests/numerics/test_roots_of_unity.py` | #325's primitive's defining laws | GREEN |

**No `verifies(...)` on the mirror module** — the property is currently BROKEN, and
a `verifies` edge on a red gate is a phantom coverage edge in the audit. Add it
when the xfails flip.

**Two Mode-8 lessons the promotion produced, both `[M]`:**
1. The α diagnostic re-implemented the recursion locally and never executed
   production (**Mode 11**). Rewritten to drive `cylindrical_streaming`; a
   production sign-flip mutation now reds exactly the 20 `l1` rows and leaves the
   15 `foundation` rows green — which independently validates the marker split.
2. **`xfail` swallows fixture SETUP errors too** — measured, falsifying the
   received belief that a setup error surfaces as an ERROR the xfail cannot hide.
   The working structure is: the xfail body asserts only its documented
   inequality, **plus a reddenable un-xfailed sibling** consuming the same
   fixtures. Both directions re-verified.

**No `catches("ERR-NNN")` anywhere yet — #326 has no ERR entry.** Minting one is
an open close-out item (§6).

---

## 3. The machinery gaps — what is actually missing

Ordered by (reach ÷ cost). **Item 1 first** — smallest change, largest reach.

1. **⚠ REWRITTEN 2026-08-02 — the original claim was wrong twice.** It read *"the
   orbifold singular set `Σ = {ξ = 0}` has no name and six ad-hoc ε-detectors."*
   Reading the code refutes **both** halves. Do not re-inherit it.

   **(a) It HAS a name — `tangential`** — with an exported, empirically-calibrated
   constant `TANGENTIAL_EPS = 4·eps ≈ 8.88e-16`
   (`numerics/spaces/angular_trace_space.py:164`), read across `sn/boundary/`,
   `geometry/boundary/`, `numerics/operator.py`, `numerics/face_layout.py`,
   `numerics/spaces/full_field_space.py`. The **boundary campaign minted it** and
   did the calibration (its docstring cites N=2..64 + Lebedev 3..53).

   **(b) It is THREE sites, not six** — and the other three are different
   OPERATIONS, so folding all six would be the L30 mistake (shared *data* ≠ shared
   *operation*):

   | plan's original site | actual predicate | what it really is |
   |---|---|---|
   | `directional.py:120` | `\|c\| ≤ _OCTANT_SIGN_EPS` ⟹ label 0 | **Σ** ✓ per coordinate axis |
   | `loss_representation:2903` | `\|μ_x\| < _MU_DIRECTION_EPS` | **Σ** ✓ for σ_x only |
   | `geometry/boundary/_specular.py:223` (*not in the original list*) | `\|μ_axis\| > TANGENTIAL_EPS` | **Σ** ✓ canonical |
   | `augmented_mesh.py:815` | `eps < τ_raw[0] < 1−eps` | τ-**degeneracy** — a closure property, NOT membership |
   | `rules_sphere.py:213` | `\|\|μ_z\|−level_μ_p\| < 1e-12` | **fiber membership** ⟹ that is **Q4** |
   | `symmetry.py:952` | `\|w_j − w_i\| > atol` | **orbit-partner weight equality** ⟹ part of the Q2 certificate |

   **The real gap** is that ONE object has three private tolerances
   (`TANGENTIAL_EPS` ≈ 8.88e-16, `_MU_DIRECTION_EPS` = 1e-15, `_OCTANT_SIGN_EPS`
   = 1e-15) because it is **sited above one of its consumers** (T13). The
   disagreement is **latent, not live** — measured, `2.7e13×` separation, zero
   ordinates in the disputed band (T14).

   **⟹ DISSOLVED INTO Q2, not fixed here.** Per T14, membership is `π_M(i) == i`,
   an integer identity — so the ε-detectors do not get consolidated, they *cease
   to be ε-questions*. Naming and retirement ride on the Q2 certificate.
2. ✅ **`DiscreteMeasure.consolidate()` — one missing verb.** `pushforward` already
   exists and its own docstring documents the atom-merging case as
   *valid-but-unreduced*. The gap is one method, not a class hierarchy.
   *(LANDED `e6f01d7e` — Q2.6.)*
3. ✅ **`degree_of_exactness` carries no subspace** (T2). `restrict`/`pushforward`
   drop it, losing a provable claim. **That gap IS the `+2.94` failure mode.**
   *(DISSOLVED by Q5.E — E1/E2 minted `ExactnessClaim(reference, degree)` and
   wired it into `DiscreteMeasure`.)*
4. **`SubgroupOfO3` has containment but no residual-group operation**, and
   `invariance_group=None` conflates *trivial* with *unknown* — which makes the
   mathematically-genuine second fold **non-composable in code**.
5. **No operator-side equivariance predicate.** `is_invariant` is one conjunct of
   THREE; the operator-equivariance and data-invariance conjuncts have no home
   anywhere, and without them a descent silently returns the group average
   `(1/|G|)Σ gAg⁻¹`. Owes a negative test: a ξ-odd source must RAISE.
6. **"half-range" already names the RESTRICTION** in `restrict`'s own docstring.
   Two operations, one word.
7. ✅ **The two `LevelStructure` producers fiber over different invariants** (T7).
   Clean-before-extend: fix this BEFORE folding anything. *(LANDED `3afb52c2` —
   Q4.)*
8. ✅ **`SO2` is documented in-tree as a "conservative upper bound" that is false at
   finite `n_φ`**, and `_so2_representatives()` checks four rotations.
   *(RESOLVED — annotated 2026-08-07: `_so2_representatives` is RETIRED
   (`symmetry.py`'s own docstring records the replacement), and the
   "conservative upper bound" prose survives only as a PAST-TENSE correction
   note in `rules_product.py` ("Until 2026-08-02 this rule tagged SO(2) and
   this paragraph defended it…"). Gap 8's content is done; Q7 keeps gaps 4,
   5, 6.)*

### ⭐ The lesson that shrinks the work (new, L-013)

**`is_invariant` is not the checking face.** `_orbit_closure` **already computes
the permutation** and discards it to return `bool`. The machinery computes the
quotient and throws it away.

> **A `-> bool` predicate that internally builds the permutation / partition /
> certificate IS the "missing" primitive. Widen the return type before minting a
> class.**

This is why the fix is one verb rather than a class hierarchy — and it is a
general review heuristic worth carrying beyond this campaign.

### The `Frame` question — SPLIT, do not force it

Measure-side quotient is **not** a Frame; basis-side isotypic restriction **is**.
Discriminators: `conjugate` has no equivariance precondition it can honour; no
operator *owns* the group (it is in the commutant of every equivariant operator at
once, and the project's rule is "an operator owns its frame IFF the frame is its
eigenbasis"). ⚠ The earlier claim *"both composites are the identity"* is **FALSE**
on the full space (`R∘M = Π ≠ I`) and true only after domain restriction — which
is exactly what the type must carry. Strongest counter, recorded and not hidden:
the fold *is* expressible as an `IndicatorBasis` Galerkin frame (`W = M𝟙`).

---

## 4. The three axes — what earns abstraction

The campaign's decision table. **Fill the "≥2 realizations?" column from the
landscape survey** (`scratch/quadrature_landscape_survey.md`, in flight at time of
writing) — do NOT fill it from intuition.

**SETTLED against the landscape survey** (`scratch/quadrature_landscape_survey.md`).
The expression fuses **FOUR** choices, not three — generation method is the fourth
(that is #325's axis).

| axis | today | ≥2 realizations? | verdict |
|---|---|---|---|
| **RANGE** | `[0,2π)` hardcoded | **YES — 8 realizations.** MoC's `[0,π)` is the same physical angle under a symmetry quotient; `angular_range` is **already a polymorphic property** in the derivations tree | **EARNED.** And it is T5's quotient, not a free parameter |
| **SPACING** | equispaced hardcoded | **YES** | **EARNED** — but see T25: the one degree of freedom that remains (the *offset*) is **DERIVED by `Σ = ∅`**, not exposed as a knob |
| **RULE — on the CIRCLE** | periodic trapezoid hardcoded | **NO — there is exactly one, and mathematically there should be** | **NOT earned.** The periodic trapezoid **IS the circle's Gauss rule**. The real variation there is *offset* and *generation*, and `[M]` the offset provably does not change exactness — ⚠ **which is NOT the same as "the offset does not matter"**: it is exactness-blind, `Sym`-visible, and **Σ-decisive** (T25). Reading exactness-invisibility as freedom is what left `δ = 0` in place |
| **RULE — on an INTERVAL** | GL / tabulated, per call site | **YES** | **EARNED.** Gauss-Lobatto is fully characterised and *declined* — and the decline reason (`τ_raw=0` breaks the M-M recurrence AND the R12a predicate) is itself the proof the rule choice is not encapsulated today |
| **EXACTNESS SPACE** | one bare integer | **YES — the most realizations of any axis**; nine distinct function spaces already share the one field | **EARNED MOST.** The only axis with a **measured falsehood shipping today** (T9b / #327) |
| **GENERATION** | trig evaluation | **YES** (#325's group action) | **EARNED** |

⚠ **The axes are NOT independent, and that is the design's crux.** For a periodic
integrand, equispaced + equal-weight IS the trapezoid and IS optimal (T1) — on the
circle, spacing and rule are welded **by mathematics**. Three orthogonal knobs
would admit meaningless combinations: `coding-elegance` Pattern 4 run backwards.
**The right factorization is "a rule on a domain" as ONE object** — the domain
determines what optimal means, which is exactly why *rule-on-the-circle* has one
realization and *rule-on-an-interval* has several.

**Already-unwired second realizations** (they exist, nothing consumes them —
so the abstraction is being paid for and not used): `equispaced`,
`gauss_chebyshev`, `gauss_laguerre` all ship with **zero production consumers**.

**Documented foreseeable needs** (grounded in issues/docs, not inferred): #325
(generation), #326 (range), **#128 QMC — a Hardy–Krause variation space**, #109
(τ-coordinate + Laguerre + Jacobi), **#123 — the project's own ratification of the
≥2 criterion**, #265, #235. `D_6h` is a documented deferral with the group already
implemented and **no `D_6h`-invariant rule in tree**.

---

## 5. Phases

> Each phase ends at a commit boundary with the tree green. Per
> `feedback_compaction_points_in_campaign_plans`, ⏸ marks a compaction point:
> **commit first, then checkpoint this file, then compact.**

- ~~**Q0 — Land what is already staged.**~~ ✅ **DONE 2026-08-02, `1f9d4818`.** The
  three promoted gates + `roots_of_unity` + `_test_helpers`, at
  `305 passed, 3 xfailed`, with `--runxfail` confirming each xfail reds on its own
  documented assertion. The plan + evidence + agent memory landed at `49bd7314`.
- **Q1 ⇒ MERGED INTO Q2** (user ruling, 2026-08-02). They are ONE carve: per T14,
  widening `_orbit_closure` to return its permutation **is** what names `Σ`, because
  `Σ = {i : π_M(i) == i}`. Splitting them would build the certificate in Q2 and then
  not use it for the question it directly answers.
- **Q2 — The orbit certificate, and `Σ` with it** (L-013; gaps 1, 2, 5; T13, T14).
  Sited in **`numerics/symmetry.py`** — the root (T13), below every consumer.
  1. ✅ **Widened `_orbit_closure`** to return `OrbitCertificate` — **LANDED
     `c9c57152`**. `is_invariant` keeps its `-> bool` face by asking whether the
     certificate exists. (The audit REFUTED the fingerprint worry: `O_h`/`I_h`
     already run full orbit closure; only the *docstring* claimed a fingerprint.
     On any d≥2 measure every non-`Trivial` tag routes through `_orbit_closure`.)
  2. ✅ **Named `Σ`** — `singular_set`, plus `stabilizer_order` and `orbits()`.
     Membership is `π_M(i) == i`: integer, exact, `[M]` stable across
     `atol` 1e-15…1e-11. Expressible only on a G-invariant measure, so the
     quotient's precondition is enforced by construction.
  3. ⬜ **Retire** `_OCTANT_SIGN_EPS` and `_MU_DIRECTION_EPS` — **OPEN, and now
     KNOWN MIS-SCOPED; see T18b.** T18's "pass `Dnh(2)`" is `[M]` confirmed as an
     *answer* (Σ(D_2h) reproduces the octant-zero set exactly), but the two sites
     do not want that answer: `_OCTANT_SIGN_EPS` needs the **chamber label** (8–26
     distinct sign vectors — Σ is only its degenerate stratum), and
     `_MU_DIRECTION_EPS` is Σ under the **single mirror `⟨σ_x⟩`**, strictly
     smaller (8 vs 16–18). Deferred into the chamber work.
     `TANGENTIAL_EPS` survives as the boundary's **matching** tolerance.
  4. **Fix what the correctness audit finds** in `symmetry.py` (user ruling: fix it
     while we are in the module). ✅ ERR-072, ERR-073, the computed lattice
     (T15b), `_vertical_mirrors` (T15c). ⬜ the `O2 → D_∞h` rename (T15d), the 1-D
     path, the tolerance asymmetries, the dead generators.
  5. ✅ **The subgroup GRAPH and its walk** (T17) — **LANDED `63d0b234`**, with the
     `O2 → D_∞h` rename (T15d). `maximal_invariance_groups(measure)`; Hasse edges
     derived from the computed containment; two realizations that verify each
     other. `[M]` `product(4,n) → D_nh`, `level_symmetric/lebedev → O_h`, walk ≡
     bruteforce on every rule, and **downward-closure now 0 violations over 1088
     pairs** (was 68). **Superseded** `invariance_group_for(params)`.
     ⚠ Cost: `test_symmetry.py` 5s → 57s, because `_orbit_closure` is an O(N²)
     Python loop per operator. Vectorising it is an open follow-up (§6).
  6. ✅ **`DiscreteMeasure.consolidate()`** — **LANDED `e6f01d7e`**. `[M]` The
     quotient's weights come out as orbit-stabilizer predicts
     (`W/w_parent ∈ {1,2}` for `|G|=2`, 8 on-mirror + 12 free), orbits `32→20`
     matching Burnside `(N+F)/2`, mass preserved to `3.6e-15`. Preserves
     `invariance_group`/`degree_of_exactness` (it changes no integral), unlike
     `pushforward`. **One method, not a class hierarchy** — as predicted.
  7. ⬜ **`align_to(K) -> Rotation | None`** — the subconjugacy certificate (T15).
     STILL open, and now **genuinely unconsumed**: the Stage-1 re-pose landed
     without needing it, because the gate wants *literal* containment in ONE
     frame (T15). Re-orientation freedom is real but has no caller yet — this is
     the honest "defer until a consumer exists" case, not a blocked one.
  8. ✅ **The registry Stage-1 re-pose** (T16) — **LANDED `e7d44f3c`**. Executing
     it moved four of T16's conclusions; see **T16b**. `AngularSymmetry` (spent +
     owed, `support` derived as S²/G⁰), `GEOMETRY_ANGULAR_SYMMETRY`, Stage 0
     (domain) + Stage 1 (symmetry **computed from nodes**), V evaluated first,
     `QuadratureSpec.invariance_group` **retired**, and `product_mu_phi`'s false
     `SO(2)` measure tag corrected to `Dnh(n_phi)`. `[M]` `tests/numerics` 1250
     passed; `sphinx -W` clean; pyright 0.

  **⚠ THE STOP-REASON RECORDED HERE WAS WRONG, and that is worth keeping.** The
  previous checkpoint held Q2 open pending a user ruling because the re-pose was
  "a live behaviour change in rule selection". `[M]` It is not: `select_quadrature`
  has **15 callers, all tests** — production never calls it. The escalation was
  written from the *shape* of the change (a selection gate) instead of from its
  *consumer set*. **Measure who calls it before escalating.**
- ⏸ **COMPACTION POINT**
- ✅ **Q3 — The measure-parameterized rule** — **LANDED** `c630153e`+`bc89b62e`+`579d5eaf`. (T12 / T12b; absorbs gap 3 and T2).
  ONE Golub–Welsch body; the families become `(α, β, μ₀)` **data**, not
  subclasses. The weight function is carried on the rule, and
  `degree_of_exactness` becomes a claim *with respect to that measure* — which
  is what makes the quotient's domain check expressible instead of a `+2.94`
  surprise, and what makes #327's over-claim structurally unspellable.
  Node-constraint (Gauss / Radau / Lobatto) is the second axis, a modification of
  the Jacobi matrix — **not** a new family, so the declined Gauss-Lobatto study
  was declining a *constraint*, not a rule.
  **Verification is free and is the point:** three existing hand-rolled rules
  (`gauss_legendre`, `gauss_chebyshev`, `gauss_laguerre`) become the ≥2 instances
  that prove the generic body, at `1e-16`–`3e-15`. Retire the hand-rolled bodies
  as it lands (`.claude/rules/coding-standards.md`).
- ✅ **Q4 — Reconcile the two `LevelStructure` producers** (gap 7, T7, T21) —
  **LANDED** `3afb52c2`. It was upstream of #325 AND #326, as the plan said.
- ⏸ **COMPACTION POINT**
- **Q5 — ⭐ THE FOLD: RANGE realized, SPACING derived** (§4; T5, T22b, T23, T24,
  T25, T26). *Design of record below — written 2026-08-02 for review BEFORE any
  code moves. Nothing in it has been implemented.*

  **The thesis in one line.** RANGE is not a parameter — it is T5's quotient by
  the geometry's **owed** residual mirror; the fold's well-posedness condition is
  `Σ = ∅` (T24); and that condition **derives** the azimuthal offset (T25)
  instead of leaving it a preference. #326 then dissolves because a folded level
  is a single arc on which the march order and the fiber order coincide (T22b).

  **Why this is the campaign's shape and not a #326 patch.** Every piece already
  exists and is unconsumed: `AngularSymmetry` names the owed residual (T16),
  `orbit_certificate`/`singular_set` compute Σ exactly (T14), `pushforward` +
  `consolidate` **are** the two halves of the quotient (Q2.6), `CHEBYSHEV_T` is
  the half-range azimuthal rule (Q3). Q5 is the verb that composes them.

  ---
  **Q5.0 — Must-precede cleanup** (`coding-standards.md` "clean before
  extending"). The fold breaks two unstated preconditions; fix them FIRST or the
  fold ships a silent wrong answer.

  1. ✅ **Certify the reflection partner map** (T23) — **LANDED `a7695148`.**
     `_find_reflections` retired; `_compute_sphere_reflection_partners` routes
     through `symmetry._orbit_closure`, which already computes the permutation
     while proving closure and requires a **bijection** with matched positions
     AND equal weights (ERR-073) — so the certification is free and the
     permutation it returns IS the partner map (L-013). An uncertified axis is
     OMITTED, so `reflection_index` raises rather than returning garbage.

     ⛔ **It caught a LIVE defect, not a hypothetical one.** `[M]` A product
     rule's mirror planes sit at `kπ/n_φ`, so `σ_x` needs `k = n_φ/2` — integer
     only for **even** `n_φ`. At `product(4, 5/7/9)` the shipped axis-0 map was
     wrong by `0.58 / 0.42 / 0.33` in the direction cosines **and was still an
     involution**, so `test_q4_2` passed on it. That map feeds the `r=0` pole
     continuation. `σ_y` is the `k=0` plane and survives at every `n_φ` — which
     is why the **fold is unconditional while the centreline map is not**.

     `[M]` No behaviour change anywhere else: every `product(4, even)`,
     `level_symmetric(4/8/16)`, `lebedev(5/11/17)`, `gauss_legendre(8)` keeps
     all three axes. `tests/numerics` 1546 passed (== baseline); `tests/sn -m
     "not slow"` 2553 passed with **the same six** deliberately-red gates and no
     new ones; sphinx `-W` clean; pyright 0.

     Teeth verified by monkeypatch back to the legacy body (positive control
     first): **6 gates red**. ⚠ Note for future gate design: with only
     even-`n_φ` rows the three-check gate stayed **GREEN** under the legacy body
     — every map it inspected happened to be correct. The odd rows are what give
     it teeth, and they are in `_SHIPPED_RULES` deliberately.
  2. **Name the mirror's PLANE** — ✅ **PULLED FORWARD from Q7 by user ruling.**
     `Z2` is realized as `_reflections("z")` — σ_z only — so the cylinder's σ_y
     is **not nameable**, an axis convention hidden inside an "abstract order-2
     group". Same class as T15c.

     **Blast radius, mapped** (`scratch/q5_mirror_plane_blast_radius.md`): 29
     sites over 4 production files, 3 test files, 1 doc page, carrying **three
     different meanings** — σ_z-specific (3-D containment + invariance), a
     plane-free `x → −x` (the 1-D arm), and `μ → −μ` (the geometry table).

     **The plumbing already exists.** `_reflections(axis: str)` accepts exactly
     `x`/`y`/`z` and raises otherwise; σ_y is constructible today. Only the enum
     member is missing — all five callers pass the literal `"z"`.

     ⛔ **A shipped rule already gets a FALSE certification.** `[M]`
     `product(4, 3)`: `σ_x = False`, `σ_z = True`, and
     `SubgroupOfO3.Z2.is_invariant(...) = True`. It is harmless *today* only by
     **dispatch accident** — two independent guards (the 1-D branch at
     `symmetry.py:706-713`, and Stage 0 rejecting every S² rule on domain before
     `Z2` is asked) — not by any designed invariant. The slab embedding really
     is `(μ, 0, 0)` (written twice in `directional.py`), so the slab's mirror is
     **σ_x** while `Z2` realizes **σ_z**; `symmetry.py:816-818`'s *"Any single
     reflection works; the choice is convention"* is simply **false**. On a
     symmetric GL μ every embedding reports `True` for all three axes, so the
     natural fixture **structurally cannot see the difference** — only an
     asymmetric μ, or an odd-`n_φ` product, exposes it.

     ⭐ **A convergence worth keeping.** A prototype `Mirror(axis)` got every
     finite relation right for free, including `Dnh(3).contains(Mirror("x")) =
     False` — which **breaks `test_dnh_reflection_in_dnh`**, whose prose ("a
     single reflection sits inside every `D_nh`") is a false generalization.
     That is **the same parity fact** Q5.0.1 measured on the node set: `D_nh`'s
     vertical mirrors sit at `kπ/n`, so `σ_x` needs `k = n/2`. The group lattice
     and the point set agree, from opposite directions — the sharpest available
     evidence that the parity is structural and not an artifact of either.

     **Template + hazards for the new family** (from the `Cn`/`Dnh` arms):
     * `_contains` decides finite×finite by **computed matrix containment**; the
       static `_NAMED_LATTICE` is consulted only when a side is CONTINUOUS. So a
       new family needs **2** relations (`⊆ D_∞h`, `⊆ O(3)`), not 5 — `[M]`
       **3 of the 5 existing `Z2` edges are dead code.**
     * ⚠ `_contains` ends in a bare `return False`, so `O3.contains(Mirror("x"))`
       measures **`False`** — wrong, and silent. Fix this in the same pass or the
       new family inherits it. *(✅ fixed in this step's own landing — the
       `Oh/O3/Dinfh ⊇ σ_a` facts are right BY THE ARM; see the Q5.0.2 record.)*
     * `_NAMED_LATTICE` is typed enum→enum, so parameterized families live
       entirely in `isinstance` arms — that is the pattern to follow.
       *(✅ still true post-G, verified 2026-08-07.)*
     * `repr` is **load-bearing**: `_GROUP_CACHE` and the lattice walk's
       `visited` set both key on it. *(✅ still true post-G.)*
     * The 1-D arm (`_check_invariance_1d`) never touches `_realized_ops` /
       `_orbit_closure` (`[M]` 0 calls) and ends in `return False`, so a new tag
       **silently reports not-invariant** on every 1-D measure until its arm is
       written. *(⛔ STALE since track G's G3, 2026-08-03 — annotated
       2026-08-07: the 1-D/3-D arm split is DELETED;
       `_check_invariance_1d`/`_3d` no longer exist and matching routes
       uniformly through the `RigidMotion` core over `_embedded_nodes`. A new
       family needs NO per-dimension arm.)*

  ---
  **Q5.1 — `DiscreteMeasure.quotient(group)`.** ✅ **LANDED `681bc49b`
  (2026-08-07).** Name the composite
  `pushforward(orbit representative).consolidate()`. Precondition already
  enforced by construction (Σ is unrepresentable on a non-invariant measure).
  *Gates:* mass preserved exactly; weights equal orbit-stabilizer's
  `W = w·|G|/|Stab|`; idempotent; and **`Σ = ∅` ⟹ every orbit has length `|G|`**
  (the free-action certificate, T24).

  **Landing note.** The representative is the orbit's **first-appearing
  member** — the only group-generic section (`ξ → |ξ|` exists only for a
  mirror; the pre-existing manual-composite test keeps that geometric section
  as the independent reference realization, gated equal orbit-by-orbit).
  Refusal = `orbit_certificate(...) is None ⟹ ValueError`, message family
  matching `singular_set`'s. `[M]` product(4,8) folds 32 → 20 (Burnside
  `(N+|Σ|)/2`, `|Σ| = 8`), weights `{w, 2w}` == `w·|G|/|Stab|` read off the
  certificate; a hand-built midpoint ring (Σ = ∅, because the SHIPPED product
  rule is equispaced-LEFT until Q5.2 derives the offset) folds 8 → 4 at
  uniform `2w` per-atom BIT-exact. Result honestly drops BOTH claims:
  `invariance_group=None` (the section is not equivariant) and
  `exactness=None` (a claim would be against `φ_*λ`, no consumer — the E2
  direct-sum precedent).

  ⭐ **The one-word "idempotent" gate was underdetermined; realized as TWO
  arms, derived from the API.** Literal idempotence holds exactly where the
  action is TRIVIAL (every node fixed ⟹ the quotient is the identity on
  `(nodes, weights)`). On a free action the fold **consumes the symmetry** —
  the folded measure keeps one member per pair, so it is no longer
  G-invariant and a SECOND `quotient(group)` **refuses**: T5's "the quotient
  CHANGES the space", enforced at runtime rather than faked as a no-op. Both
  arms gated. Teeth `[M]` (in-process monkeypatch): skip-consolidate 4 red,
  skip-fold 6 red; the invariant-integral gate is blind to skip-consolidate
  BY THEOREM (consolidate changes no integral) — the count/weight gates are
  its designed catchers.

  ⚠ **A docstring `:label:` is invisible in an un-`automodule`'d module** —
  the first `-W` run failed with `equation not found` at both `:eq:` sites
  because `numerics/measure.py` is among the 218 unrendered modules (the
  corpus-integrity lesson, hit live). The labeled equation now lives on the
  theory page (`discrete_measures.rst`, with its vv-status pair); the
  docstring keeps the math unlabeled. The page pass also closed owed drift:
  consolidate's propagation-table row (owed since Q2.6), the rotting
  "four"/"fifth" ordinals, and stale `μ.space` cells → `μ.support`.

  ---
  **Q5.E — ⭐ THE EXACTNESS SPACE, pulled forward** (T2, T12b, T28). User
  ruling 2026-08-02: *"pull the exactness space forward, it's the root. Let's
  implement the principled product machinery."* It is a **prerequisite of
  Q5.2**, not a later phase — without it the combinator's derived degree is
  wrong by construction (`[M]` a naive `min()` gives **1**).

  The unification: `V` is the span of the first `m` functions of the reference
  measure's own **orthogonal system**, and which system that is follows from
  the measure — weight-on-an-interval → orthogonal polynomials (algebraic);
  uniform-on-the-circle → the Fourier basis (trigonometric); uniform-on-`S²` →
  spherical harmonics. T12b's "the generating measure IS the exactness space"
  holds in general once the system is **read off** the measure instead of
  assumed to be polynomials.

  1. ✅ **E1 — mint the claim** — **LANDED `9e74faa1`.**
     `numerics/exactness.py`: `OrthogonalSystem` (a VALUE, per T12c),
     `ReferenceMeasure` (a Protocol), `ExactnessClaim(reference, degree)` with
     `combined_with` as a **partial meet**. `GeneratingMeasure` satisfies the
     protocol structurally and reports `ALGEBRAIC` as a *theorem of the
     construction* (a three-term recurrence generates polynomials by
     definition). 21 law tests; `tests/numerics` 1585 passed; pyright 0.

     ⭐ **The fork, resolved and recorded so it is not re-opened.**
     `GeneratingMeasure` is a Golub–Welsch *generator* — more than a claim
     needs, and less than every rule has (the circle's Fourier system has no
     such recurrence; Lebedev's reference generates nothing). So the claim is
     typed against the **protocol**, and the generator is the sub-case that
     *also* constructs. Widening `GeneratingMeasure` would weld "what the claim
     is about" to "how the rule was built" — the very welding T28 indicts.

  2. ✅ **E2 — wired into `DiscreteMeasure`.** One `exactness` field;
     `degree_of_exactness` and `generating_measure` are now **derived
     read-only views** (no separate storage — the `Quadrature.mu_x`
     precedent), so the read sites survive while the storage that let a degree
     float free is gone. `with_metadata`'s `degree_of_exactness: int`
     parameter **retired, not renamed** — it attached a half-claim after the
     fact and `[M]` had zero callers.

     ⭐⭐ **The suite caught a real design error on the first run, and the
     correction is the most valuable thing in E2.** `_combined_degree` was ONE
     law serving **two different operations**:

     * `__add__` (direct sum) lands on the **same** space — summing two rules
       for `λ` gives a rule for **2λ** (`[M]` its weights sum to 4, not 2), so
       its reference is `λ₁ + λ₂`, not the shared `λ`;
     * `__mul__` (tensor) lands on the **product** space, so its reference is
       `λ₁ ⊗ λ₂`.

     Keeping a factor's reference in either case asserts exactness against a
     measure the composite is not exact against. Split: `tensor_with` mints a
     `ProductMeasure`; the direct sum carries **no claim** (nothing consumes
     it, and a `SumMeasure` built on speculation is worse than silence).
     `combined_with` then had no consumer and was **retired** — it was never
     the direct sum's law, only a plausible-looking one.

     ⛔⛔ **AND THE CAMPAIGN'S OWN "6.2832 BUG" WAS PARTLY A MISREADING.** The
     tree cited `gauss_legendre(4) * gauss_chebyshev(4)` "advertising degree 7
     while integrating the constant 1 to **6.2832** instead of **4**" as proof
     that a mixed-reference product must claim nothing. `[M]` Re-measured with
     the reference nameable:

     > `2π = 6.2832` **IS** the correct mass of `legendre ⊗ chebyshev_t`, and
     > the product **IS** exact to degree 7 per axis against
     > `dx ⊗ (1−y²)^(−1/2)dy` — worst error **4.16e-13** over every
     > `(a,b) ≤ 7`, with degree 8 missing by **1.5e-2**, so the bound is TIGHT.

     The expected `4` was the **Lebesgue** product — the wrong reference. So
     the real defect was always *the claim was unfalsifiable for want of a
     reference*, and the "drop the degree" guard was a **conservative
     workaround for a missing type**, not a mathematical law. Naming the
     reference makes a correct, tight claim representable that the tree could
     previously only suppress. Two tests were **inverted** accordingly.

     **Lesson worth carrying beyond quadrature:** when a claim is compared
     against the wrong reference, the *measurement* is real and the
     *conclusion* can still be backwards. "Suppress the claim" and "name what
     the claim is about" both make the red go away; only one of them is right,
     and the cheap one hides a correct result for years.

     Also fixed in passing: `_build_level_symmetric_arrays`' return annotation
     said 7-tuple while Q4 (`3afb52c2`) had made the body return 9 — a stale
     signature no test could fail on and `sphinx -W` could not see.

  3. ✅ **E3 — the periodic trapezoid on `S^1`, as a first-class rule.**
     `numerics/quadrature/rules_circle.py`: `periodic_trapezoid(n, *, shift)`
     carrying `ExactnessClaim(UNIFORM_ON_CIRCLE, n − 1)`, with the classification
     constants `NODE_ALIGNED` / `STAGGERED`. 30 tests × parameterisation = 150;
     `tests/numerics` 1735 passed; `sphinx -W` clean; pyright 0.

     ⭐ **Three design calls, each derived rather than chosen:**

     * **The nodes are POINTS, not angles** — `(cos φ, sin φ)` from
       `roots_of_unity`, not a stored angle chart. This is the call the rest of
       the step hangs off, and the reason is `Σ`: on-axis `sin` is **exactly
       `0.0`**, so `Σ = {ξ = 0}` is decided by an EQUALITY. Under the
       construction it replaces, `np.sin(np.pi) = 1.22e-16` and the membership
       test would need an `eps` — i.e. **the fold's well-posedness condition
       would be a tuning parameter**. `[M]` The mutation control is brutal: the
       `linspace+cos` nodes differ by **2.22e-16** and that alone destroys the
       bit-exact mirror AND one of the two exact zeros.
     * **The shift is a `Fraction`, and REQUIRED.** Rational because the exact
       generation rests on integer arithmetic (`Fraction(0.1)` would ask for a
       root of unity of order `n · 3.6e16`); required because the shift is
       **exactness-invisible and yet `Σ`-decisive** — a default is exactly how a
       parameter no gate can see goes unexamined (L19 route 2). It lives in
       `ℚ/ℤ` (a whole step is the identity relabelled), so the reduction
       canonicalises rather than restricts.
     * **`NODE_ALIGNED` / `STAGGERED` are theorems, not conventions.** The node
       set is mirror-symmetric iff `2s ∈ ℤ`, so within `[0,1)` there are
       **exactly two** such shifts — the classification T25 wanted. Which is why
       the shift is a *parameter* and the two values are *named constants*,
       never a boolean flag (anti-pattern #3 avoided by deriving, not by taste).
       Caveat pinned by its own test: staggered empties `Σ` only for **even**
       `n` — at odd `n` the node `m = (n−1)/2` lands on `φ = π`.

     ⭐ **The step found a live mis-count in the shipped product rule.** `[M]`
     `product_mu_phi(4, 8)`: the node set meets the axis **twice per level**, but
     `|Σ| = 4` by equality against `8` by any sane tolerance — because
     `np.sin(np.pi) ≠ 0`. So today `Σ`'s SIZE depends on the tolerance used to
     measure it, which is not a defensible state for a quantity that is a fold's
     well-posedness condition. Generating the azimuths as roots of unity makes
     both counts `8`. Recorded in `rules_product`'s docstring; the fix is Q5.2's
     substitution, not a patch here.

     Also fixed in passing: `product_mu_phi` called its azimuthal factor the
     "midpoint rule" — false, it is left-endpoints (shift `0`, not `½`). The
     *exactness* sentence beside it was right anyway, **because the shift cannot
     change the degree** — a small worked instance of the same lesson.

     Renamed `SPACE_CIRCLE` from `"[0,2π)"` to `"S^1"` (4 sites). It was the
     off-pattern member of a family whose sibling is `SPACE_SPHERE = "S^2"`, and
     it made the tag assert a *coordinate* where the rest assert a *space* —
     which is precisely the distinction E3 exists to draw.

  4. ✅ **E4 — the sphere product theorem + the azimuthal substitution.**
     `spherical_product_claim(polar, azimuthal)` in `rules_product.py` carries
     the theorem and its derivation; `product_mu_phi` now composes **two
     registered rules** and DERIVES its degree. Value unchanged and tight.
     `tests/numerics` 1753 passed; `tests/sn -m "not slow"` **6 failed / 2554
     passed — the SAME six**; `sphinx -W` clean; pyright 1 (the #288 residual).

     `[M]` **What the substitution moved** (isotropic cyl MMS, nx=40, S(4,8)):
     MMS L2 `5.389276744519e-04 → …525e-04` (1.2e-12); scalar flux **6.7e-16**;
     per-ordinate ψ **3.0e-05**; and the degree-7 monomial residual **improved**
     `1.1e-16 → 1.1e-17`. The movement is entirely in the ξ-odd sector, where
     `test_azimuthal_mirror_symmetry.py` already measured the defect magnitude
     as **ordering-INVARIANT** — so this is not correct→incorrect, it is
     **noise-decided → named**.

     ⭐ **Two committed tests went red, both CORRECTLY, and both were good
     tests doing exactly their job.** Worth carrying because the pattern
     generalises: *a fixture can be discriminating for a reason nobody pinned.*

     * `test_domain_is_the_tangential_band_gamma_plus` — its **ACTIVATION
       guard** fired. It discriminated only because `product(2,4)`'s `±π/2`
       azimuths came out of `cos(linspace(…))` as `±6.1e-17` instead of `0`;
       exact nodes made them exactly tangential and the two classifiers agreed.
       `[M]` Re-measured across EVERY shipped rule (`product` ×3,
       `level_symmetric` 4/8, `lebedev` 3/5, all three axes): **none
       discriminates any more** — every projection is exactly `0.0` or far
       above the 4-ULP band. So `TANGENTIAL_EPS` existed to absorb the product
       rule's round-off, and the product rule was its only source. Re-posed to a
       **constructed** sub-EPS measure, which cannot decay that way again.
     * `test_trig_nodes_hide_the_tie_break_behind_rounding_noise` — a
       characterization test of the pre-#325 state that **predicted this change
       in its own docstring**. Split into a pair: the historical claim stays as
       a control, its inversion (production ties are now exactly `0`) is the new
       row. Keeping both is what makes the second legible as a defect removed.

     ⚖ **On that file's "#326 blocks #325".** The blocking condition was that
     exact nodes hand the ordering decision to a tie-break, so one must be
     CHOSEN rather than inherited from noise. It was — `kind="stable"` (η
     ascending, ties by increasing φ) — so the condition is **discharged, not
     violated**. Choosing ≠ ranking: that file's headline result stands (no
     ordering is more correct), and Q5.3 may re-pose it as the azimuthal march.

     ---
     *Original E4 specification, kept for its derivation:*

     T2 already proved the shipped
     `min(2n_μ − 1, n_φ − 1)` is **tight**, and *why*: a degree-`d` spherical
     monomial's azimuthal factor is a trig polynomial of degree `≤ d`, so the
     two bounds coincide. So E4 does not change the formula — it gives it
     correctly-typed inputs and a name, turning a hand-written `min()` over two
     bare integers into the theorem it always was.

     **The derivation, worked out at E3 and worth keeping** — the mechanism is
     sharper than "both bounds coincide", and it is the reason the theorem is
     NOT `tensor_with`. On `S²` the product rule FACTORISES on a separated
     integrand, so for `Ω^a Ω^b Ω^c` of total degree `d`, split off
     `sin^k θ` (`k` = the combined transverse power):

     * **`k` odd** ⟹ the polar factor `(1−μ²)^{k/2}·μ^{d−k}` is **not a
       polynomial in μ** — but its azimuthal partner `cos^aφ sin^bφ` with
       `a+b = k` odd contains only ODD harmonics, so every `m ≠ 0` and the
       azimuthal rule reproduces the exact `0`. *The azimuthal factor
       annihilates exactly the terms whose polar factor is non-polynomial.*
     * **`k` even** ⟹ the polar factor IS a polynomial, of degree exactly `d`,
       so the polar rule needs `2n_μ − 1 ≥ d`; the azimuthal trig degree is
       `k ≤ d`, so the azimuthal rule needs `n_φ − 1 ≥ d`.

     Hence `min`, and hence why the two systems (ALGEBRAIC × TRIGONOMETRIC) do
     NOT compose by `ExactnessClaim.tensor_with` — that method correctly returns
     `None` here, because the product's claim is about **spherical harmonics on
     `S²`**, reached through the EMBEDDING `(μ, φ) ↦ Ω`, not about a mixed
     tensor system on the square `[−1,1] × S¹`. E4 is the theorem the
     `ProductMeasure` docstring says "belongs with the rule that applies it".

     ⚠ **E4's honest form needs Q5.2 step 1 first, and that is a scheduling
     call, not a technical one.** To hand the theorem its two typed inputs, the
     azimuthal factor must BE `periodic_trapezoid(n_φ, shift=NODE_ALIGNED)` —
     constructing that rule *only for its claim* while still building nodes from
     `np.linspace` would be a twin path asserting a claim about a measure it
     does not build. So E4 ⟹ substitute the azimuthal nodes. Precedent: Q3 took
     exactly this trade for the GL polar nodes (user ruling, 2026-08-02) and
     re-baselined the SN **slab** snapshots. **User ruled 2026-08-02:
     substitute now.**

     ⛔ **THE COST QUOTED FOR THIS RULING WAS INCOMPLETE — recorded because the
     understatement is the lesson.** It was presented as "a pure-ULP node move,
     `2.2e-16 … 6.1e-16`". The node move IS that. But the exact azimuths ALSO
     change the **level ordering**, which is discrete, not ULP:
     `η = sinθ·cos φ` and `cos φ_m = cos φ_{n_φ−m}` bit-exactly, so a level has
     only `⌊n_φ/2⌋+1` distinct η — `[M]` **5 of 8** at `n_φ = 8`. Under
     `linspace` round-off manufactured `n_φ` fake distinctions, so **the
     intra-pair order was previously decided by NOISE**; with real ties it is
     decided by the sort's tie-break. `[M]` The orderings differ at every `n_φ`
     tested (4, 8, 12). That is the exactness-as-diagnostic lesson landing on
     the campaign's own change, one step after the campaign wrote it down.

     `[M]` What does NOT block it: `test_product_bit_identical_to_legacy_adapter`
     compares `Quadrature.product(...)` against `product_mu_phi(...)`, and the
     wrapper CALLS the rule — so both sides move together. (It is the same
     compare-a-value-with-itself-through-a-wrapper shape `test_rules_1d`'s own
     docstring already flags; worth re-posing when touched.)

  ---
  **Q5.2 — ⭐ UNWELD the product: a COMBINATOR over registered rules** (T28).
  Re-posed 2026-08-02 by user diagnosis — this is no longer "add an offset", it
  is the root the fold's other symptoms hang off.

  1. **Make the azimuthal factor a MEASURE.** Today it is `np.linspace` plus a
     scalar — an object that can carry no degree, no group, no generating
     measure, and cannot be substituted.

     ⛔ **This step named the WRONG factor until E3.** It said "`measure.
     equispaced` (the **midpoint** rule) is the factor the fold wants, and its
     own docstring already gives T25's reason". Both halves are now false, and
     the way they are false is the whole point of E3: `equispaced` is a rule on
     an **interval**, carrying an ALGEBRAIC degree of **1**. Handing it to the
     product would be the exact confusion the exactness carve exists to prevent
     — same nodes, wrong space, wrong claim. The factor the fold wants is
     `periodic_trapezoid(n_φ, shift=…)` (E3, `647ece5a`), on `S^1`, carrying
     TRIGONOMETRIC degree `n_φ − 1`. Its `shift` is the parameter that used to
     be spelled "midpoints vs left-endpoints", and T25's reason now lives there
     as a **classification** (exactly two shifts are mirror-symmetric) rather
     than as a passing remark in an interval rule's docstring.
  2. **Take the pair as arguments, not as literals.** `product_mu_phi` welds
     `gauss_legendre_on_mu × linspace`. The seam is what lets `Σ = ∅` *select*
     the azimuthal rule instead of a `half_range=True` flag existing.
  3. **Derive every intrinsic property downstream of the factors** — nodes and
     weights through `⊗` + `pushforward` (both exist), the group **computed**
     (not declared — `[M]` today's `D_{n_φ h}` tag is *correct*, but a
     declaration is unfalsifiable and this module has shipped three false ones:
     ERR-072/073/074), `Σ` and the level structure downstream of the nodes.
     Retire the hand-written `2*n_mu - 1`: `[M]` `polar.degree_of_exactness`
     already holds it — the fourth in-tree instance of **T20**.

  ⛔ **BLOCKED on the exactness space — do NOT compose degrees before fixing
  it.** `[M]` a naive `min(p_polar, p_azimuthal)` yields **degree 1**, because
  `equispaced.degree_of_exactness = 1` is the *algebraic* degree on an interval
  while the same nodes on the *circle* are the periodic trapezoid, exact to
  *trigonometric* degree `n − 1`. One integer, two incompatible claims. So the
  **EXACTNESS SPACE** axis (§4's "earned most", T2/T12b) is a **prerequisite of
  Q5.2**, not a later phase — pull the minimum of it forward, or land Q5.2's
  nodes/weights/group and leave the degree hand-written with a `[M]`-justified
  comment until it lands. Decide explicitly; do not let it drift.

  ✅ **UNBLOCKED — Q5.E landed in FULL (E1–E4), and E4 already absorbed part
  of this step** (annotated 2026-08-07). `[M]` against the tree:
  `product_mu_phi` now composes the two registered rules internally and
  DERIVES its degree through `spherical_product_claim` (steps 1 and 3's
  degree-half are DONE), but its signature still takes the ORDERS
  (`def product_mu_phi(` at `rules_product.py:215`) — so step 2 (the pair as
  ARGUMENTS, the actual combinator seam) and step 3's group-computed /
  Σ-downstream halves REMAIN, and the §6 "conservative"-word fix at `:238` is
  still owed. Re-scope the residual at pickup against `rules_product.py`
  itself (§7.2 — the plan's step boundaries no longer match the tree's).

  ✅ **RESIDUAL LANDED `e0554df9` (2026-08-07) — Q5.2 COMPLETE.**
  `spherical_product(polar, azimuthal) -> (DiscreteMeasure, LevelStructure)`,
  the constructive twin of `spherical_product_claim`; `product_mu_phi` is the
  named family (order guards + the two registered default factors + delegate).
  The fold's spelling, flag-free:
  `spherical_product(gauss_legendre_on_mu(n), periodic_trapezoid(m,
  shift=STAGGERED))`.

  * **Step 3's group-half:** `_derived_product_group` — three
    `RigidMotion.preserves` generator checks ON THE FACTORS (C_n rotation +
    vertical mirror on the circle rule; μ-mirror on the polar rule) ⟹
    `Dnh(n_φ)`; an unrealized family (`C_nv`/`C_nh`/`C_n`) is REFUSED naming
    it (the field means MAXIMAL — an under-claim is as false as an
    over-claim). Inputs are MEASURED on the factors, never inherited from
    their tags (ERR-072/073/074). The Hasse walk over the assembled nodes is
    the committed independent realization; teeth `[M]`: hardcoding the
    derivation reddens the refusal gate while walk-agreement stays green —
    each gate catches its own failure direction.
  * **Step 3's Σ-half — ⭐ REFUTED-AS-OWED:** Σ was ALREADY downstream (a
    query via `singular_set`, nothing stored). Recorded in the combinator's
    docstring; no code was owed.
  * **The stated gates, all discharged `[M]`:** same-factors bit-identity —
    pre-carve capture on (4,8)/(2,4)/(4,12)/(3,5), every array `array_equal`;
    `Σ` computed as ∅ under σ_y at STAGGERED (integer equality) with the fold
    composing end-to-end through `quotient` (32→16, mass exactly 4π — T24's
    number) and `|Σ| = 2n_μ` at NODE_ALIGNED; derived group == the previous
    declaration on every shipped order (the pre-existing 5-config
    parameterized gate now exercises the derivation). **T16b respected: the
    DEFAULT stays NODE_ALIGNED** — the offset became *selectable*, not
    re-defaulted, so no consumer moved; flipping any consumer's default is
    the fold-wiring steps' business, under the ruled snapshot re-baseline.
  * E3's same-nodes/wrong-space confusion is now unrepresentable END TO END
    (an interval rule as the fiber is refused BY ITS SUPPORT); off-circle
    nodes refused; the "conservative" wording corrected to TIGHT; the module
    head's "not verified at construction" deferral contract updated (it now
    is — at construction, from the factors).

  *Gates:* `Σ` **computed** as ∅ under the geometry's owed mirror; the group
  computed and equal to the previous declaration on every shipped order; nodes
  bit-identical when the same factors are passed (the un-welding alone must move
  nothing).
  ⚠ **Measure the consumer set before calling the offset user-facing** (T16b):
  changing which azimuthal rule `Quadrature.product` defaults to reaches **every**
  consumer, 2-D Cartesian included, not just the cylinder.

  **Q5.3 — A level becomes an ARC.** On a folded measure `η` is strictly monotone
  in `ω`, so `level_indices` and `fiber()` **are the same order** (T22b, `[M]`
  `[4 3 2 1 0]` both ways). The two accessors MERGE — this is the retirement Q4
  deliberately deferred.
  *Gate:* `level_indices[p] == fiber(p)` for every level of every folded rule —
  T22b's measurement promoted to a permanent test.
  ⛔ **The gate's `==` spelling REFUTED 2026-08-07 at pickup** — `[M]` probe on
  every level of folded staggered (4,8)/(2,4)/(3,6): the two accessors were
  exact **REVERSES**, never equal. Both halves of the plan sentence were true
  and the composition wasn't: T22b's "`[4 3 2 1 0]` both ways" compared
  η-ascending against the **march** order (ω-DESCENDING, from Σ at ω=π toward
  Σ at ω=0), while `fiber()`'s lexsort was the ω-ASCENDING chart. Same total
  order, two charts, opposite orientation — `∂η/∂ω = −sinθ sinω < 0` is the
  theorem, and it fixes the orientation the sentence elided. The landed gate
  asserts the sharper form: strict η-injectivity AND strict ω-descent along
  the stored order, per level, per config.

  ✅ **LANDED `65ff5bb0` (2026-08-07) — Q5.3 COMPLETE.** `fiber()` deleted;
  minted `LevelStructure.quotient(parent=, onto=)` — the structure **descends
  by pure selection** along the quotient the measure already took (a quotient
  never moves a node ⟹ charts descend bit-identically; each level's order is
  the parent's η-order **restricted** to survivors, so the sort convention
  stays spelled once, in the producer — no re-sort, no twin). Fiberwise
  precondition certified per level by mass conservation; a level-merging fold
  (σ_z) refused loudly; keyword-only measures (a positional swap of two
  same-typed measures would be silent). Gates +7: the arc gate (staggered
  (4,8)/(2,4)/(3,5)-with-Σ-endpoint + T22b's own NODE_ALIGNED (4,8), carrying
  `verifies("folded-level-arc")` → the new labeled equation on the
  discrete-measures page), the selection contract against an independently
  re-derived node-bit map, the two refusal arms; the two `fiber()` tests
  REWIRED (they assert chart injectivity — order-free). Teeth `[M]` 4/4
  mutations caught by exactly their designed catchers (control 7 reds;
  ascending-ω re-sort = the retired orientation → arc×4+selection;
  guard-off → ONLY the merging refusal; arctan2 recomputation → ONLY the
  selection gate, moving 4/16 azimuth entries on staggered(4,8), 2/9 on
  (3,5) — measured, not assumed). ⚠ The teeth driver's first grep scanned
  for `FAILED` lines that `-q --tb=no` never emits — vv #17's catalogued
  instrument failure re-hit and re-fixed before any verdict was trusted.
  `azimuth`/`hemisphere` FIELDS deliberately survive the accessor merge:
  the ω-chart's imminent consumers are Q5.4–Q5.6 (arc endpoints, Δω, T3's
  α closed form).

  **Q5.4 — Re-pose the R12a predicate on the integer** (T26). Replace
  `0 < τ_raw[0] < 1`, which decides a structural question on a `5.6e-16` gap AND
  conflates Σ-occupancy with the double cover.
  *Gate:* classification invariant under trig round-off; odd `n_φ` no longer
  flips. ⚠ Requires first deciding **which** of the two facts R12a actually
  wants — do not port the conflation forward.

  ✅ **LANDED `679cba61` (2026-08-07) — Q5.4 COMPLETE.** The adjudication
  answered: R12a wants **BOTH facts, as separately-named conjuncts** (T26's
  own closing line) — `MarchStart(on_edge_node, degenerate)` minted in
  `pole_angular_closure`, carrying ⟺ neither; the consumer
  (`radial_characteristic_levels`) reads the facts; `_SEED_TAU_EPS` and the
  interval retire; the trichotomy (`on_edge_node ⟹ 0`, `degenerate ⟹ 1`,
  neither ⟹ strict interior) is a bit-exact gated THEOREM with no epsilon.
  * ⛔ **The step's "5.6e-16 gap" premise was STALE at pickup** — `[M]`
    post-E3 (roots-of-unity azimuths) the trichotomy is ALREADY bit-exact at
    every shipped `n_φ` including 5 and 7: T26's `0.9999999999999994` was
    measured on the pre-E3 `linspace`+cos rule. The fragility half of T26 was
    cured upstream by E3; what Q5.4 cured is the CONFLATION half + the eps
    guard + the encoding. The theorem gate now pins the exactness (a
    `linspace` regression reds it bit-exactly).
  * ⛔ **A theory-page falsehood found and fixed in passing**: the
    circle-vs-interval section claimed an odd azimuthal count "would then
    carry a seed" — that documented the pre-E3 BUG as structure. `[M]` odd
    `n_φ` is `degenerate` (the mirror pair straddling π shares η bit-exactly),
    not carrying. Sharp corrected form now on the page: **on a full circle
    σ_y closes the march at EVERY parity** (node ON the plane at even, tie
    ACROSS it at odd); **only the folded arc opens a genuine off-node start**.
  * `[M]` **The fold future pinned before wiring**: folded staggered reads
    `τ_raw,0 ≈ 0.22` strictly interior on every level — CARRYING (the
    sphere's structure); folded NODE_ALIGNED starts ON Σ — not carrying.
    Ten-config classification gate carries
    `verifies("sn-direct-seed-r12a-predicate")` — the label went from ZERO
    verifying tests to ten.
  * `[M]` Teeth: control/swap/drop-degenerate each caught by exactly their
    designed catchers (the fact-swap leaves the carrier suite GREEN — the
    consumer reads only the conjunction, which is why the fact-level gates
    exist); the OLD eps-interval pose restored is measured **BLIND (35/35
    green)** — extensionally equal on every shipped post-E3 realization, so
    the step is behavior-neutral by measurement (`tests/sn -m "not slow"`
    2861 passed / the 6 checkpoint reds, zero new).

  **Q5.5 — The `[½,1]` clamp: SPLIT it in two, retire one half, promote the
  other.** ✅ **ADJUDICATED — see T27.** The absorption retires; the `[0,1]`
  membership becomes a guard that RAISES. The accompanying gate is specified in
  T27 and is reddenable by reverting Q5.2's offset.

  ✅ **LANDED `0ec701d4` (2026-08-07) — the guard + gates half; the
  absorber half is fold-coupled and falls at Q5.6.** The guard is
  `_assert_tau_within_unit_interval` on `morel_montry_tau_raw_per_level`
  (BOTH arms; closed interval — the {0,1} endpoints are legal march
  starts per Q5.4 — NaN-catching; message fragment "outside its own
  angular cell"). Gates: `tests/sn/sweep/test_tau_arc_wellposedness.py`
  (10 items) — guard positive (the NODE_ALIGNED full circle ATTAINS
  both closed endpoints) + negative (ω-ordered members ⟹ τ_raw = 1.079
  ⟹ raises) + the T27 pair at n_φ ∈ {8,16,32,64}: mechanism (Σ = ∅
  computed via `singular_set`) and consequence (τ_raw ⊂ [1/5, 4/5] +
  the reversal identity `np.testing.assert_array_equal` — NO epsilon).
  New label `morel-montry-folded-arc` (`structured_geometry.rst`, the τ
  doctrine home; the SN page cross-references it) ×4 in the matrix.
  `[M]` teeth: δ=0 (patch `orpheus.numerics.quadrature.STAGGERED` →
  NODE_ALIGNED BEFORE collection) reds BOTH legs at every n_φ (8/10
  red); a no-opped guard reds exactly the negative (1/10).
  `seam_quad` promoted to `tests/sn/_test_helpers` (2 consumers).

  ⭐ **The landed fold UPGRADED T27's identity.** Reversal residual
  **0.0 bit-exact** at every n_φ — T27's `[M]` table (0 / 6.66e-16 /
  8.22e-15) was measured on its HAND-BUILT pre-Q5.3 fold; the Q5.3
  selection descent bit-copies the charts, so the landed machinery is
  exact and the gate asserts with no epsilon (the campaign's
  exact-by-construction discipline, same as the Q5.4 trichotomy).

  ⭐⭐ **The guard's first catch was a LIVE latent defect, found by the
  landing battery itself.** Seven operator-equivalence tests
  (`test_removal_form_matvec_sweep.py` + `test_inverse_operator_equivalence.py`,
  the `sphere_2g` rows) built `SNMesh(spherical 1-D, level_symmetric(4))`
  — a raw 3-D rule (24 UNSORTED `mu_x` with duplicates, Σw = 4π) fed to
  the sphere arm's weight-sum edge convention: `[M]` τ_raw ∈
  [−20.310, 1.130] with **23 of 24 outside [0,1]**, consumed SILENTLY
  by the *unclamped* sphere closure (a recurrence thread weight of
  +21.3). The tests stayed green pre-guard because both compared
  spellings SHARE the τ — the measured functional's invariance group
  contains the defect (vv Mode 12). Fixtures re-propped onto
  `gauss_legendre(4)` (their claims are quadrature-incidental; no
  snapshots involved); **#336** tracks the refuse-or-reduce design for
  `SNMesh(SPH)` + non-μ-line rules, incl. the march-start family
  table's missing LS-on-sphere row.

  ⚠ **Instrument (the vv #17 family's 4th appearance, NEW face):** the
  first δ=0 teeth run measured **10/10 BLIND** — the driver had
  imported `tests.sn.sweep.test_tau_arc_wellposedness` and patched ITS
  `_FOLD_SHIFT`, while pytest (`--import-mode=importlib`; `tests/sn/`
  has no `__init__.py`) built a FRESH twin module and ran that. Fix:
  patch the mutation's SOURCE (`orpheus.numerics.quadrature.STAGGERED`)
  before collection — in-process ≠ same-module; verify the patched
  object is the one pytest runs.

  ---
  **Q5.6 — Acceptance** (T23 — the three xfail rows CANNOT serve).

  ⏳ **IN FLIGHT (2026-08-08) — ladder + landed-so-far.** Rulings taken
  (user, this session): the wiring shape is **primitive + named factory +
  mesh refusal**; the factory is **`Quadrature.folded_product`** (the
  quotient understanding ratified over "half-range" — T4/T5/T6 grounds);
  the `SNMesh(CYLINDRICAL)` refusal predicate is **the Q5.4 march-start
  facts** (every level carrying: `degenerate` catches the double cover,
  `on_edge_node` the non-free fold).

  * ✅ **6.0/6.1 LANDED `8416253c`** — `Quadrature.quotient(group)` (the
    thin lift; delegates measure Q5.1 + structure Q5.3, owns no refusal)
    + `folded_product(n_mu, n_phi)` (offset DERIVED per T25; odd n_φ
    refused — `[M]` staggered odd puts |Σ| = 4, one node per level ON
    the mirror); 10 gates incl. the σ_x pole map surviving the fold as
    an involution (the reconciliation note's owed gate). `[M]` en route:
    **σ_y-folded level_symmetric is NOT carrying** (its |μ_z|-keyed
    levels keep η 2-to-1 from the ±μ_z pair) — the LS-on-cylinder
    ruling (fold by the full vertical pair vs retire) is OPEN, owed to
    the user before the migration.
  * ✅ **6.2b LANDED `70b5b6d2`** — ⭐⭐ **the System-B seed machinery was
    SPHERE-HARDCODED in three places**, caught by the wiring's own L0
    battery (82 % / 158 % flat-equilibrium errors), diagnosed by the
    explorer with a closed-form reproduction of the exact garbage seed
    values (+3.7217/−0.5801 = the two GL level families under the
    sphere's Legendre fold AT μ=±1, a point outside the arc;
    Mode-12-masked in every two-spelling test): (1) the q̄½ source fold
    → the arc's OWN Gauss family (GC1 in x = cosω, T25;
    `T_k(±1) = (±1)^k` mirrors `P_ℓ(±1) = (±1)^ℓ`); (2) the DD march's
    baked |μ|=1 → PATH widths `Δr/|η_start|` via
    `_march_start_cosines` (flat-INVISIBLE — its catcher is the
    linear-ψ½ term gate); (3) the Emission's hard ½ → the constant's
    reproducing weight `1/Σw`. `[M]` flat L0 82 % → 2.8e-13; c = 0.4
    → 158 % → 2.2e-13; sphere byte-identical (243 regression tests
    green). The false "(the sphere)" guard texts corrected across 5
    files + the theory ladder passage.
  * `[M]` **consumer-audit measurements banked**: T3's α closed form
    holds on the folded arc to ≤ 2.2e-16 per level with the dome
    closing exactly (the 6.4 gate's content); ⚠ **the ξ-carrying P1
    harmonic slot reads +6.49 GARBAGE on the folded rule for a FLAT
    flux** (T6: ξ-odd is not in the quotient's space — must be
    EXCLUDED structurally, not computed).
  * ✅ **6.2c LANDED `b99fcbc3`** — the folded rule's harmonic
    machinery binds the σ-EVEN sub-basis:
    `MirrorEvenSphericalHarmonicBasis` (rectangular layout, odd
    columns structurally zeroed — the same mechanism as the |m|>l
    padding, so every shape contract survives; parity DERIVED by
    mirror-evaluation, never hand-listed — ⚠ the basis measures
    azimuth FROM μ_y, so "mask the sin branch" would zero the WRONG
    slots); `Quadrature.folded_by` provenance marker (the measure
    stays provenance-free per Q5.1); one seam covers the kernel + its
    adjoint + future curvilinear DSA. `[M]` even-block Gram ==
    continuum to 1.1e-15 on the quotient; the garbage slot bit-zero.
    ⭐ **The explorer's stale-premise catch re-scoped the acceptance**:
    the #229 aniso MMS is the aniso-ANSATZ case with ISOTROPIC
    scattering — the ξ-odd defect never touched it; NO converged
    P1-scattering cylinder test exists (the only P1 cylinder gate is
    route-equivalent = folded-blind); the kernel-tier gates here are
    the fix's coverage, and a converged P1 cylinder solve gate is 6.4+
    material.
  * **Remaining ladder**: 6.3 the flip (mesh refusal on the
    march-start facts + call-site migration + snapshot re-baselines +
    LS sites → folded_product per the user's ruling, with the
    double-fold capability FILED not built; #327 is the NEXT Q-track
    step after Q5.6 per the same ruling) → 6.4 acceptance (#229 falls
    + recovers an order; T3 α gate; absorber retires — Q5.5's
    deferred half; τ ∉ {0,1}; the three xfails re-posed) → 6.5
    retirements + theory pass + #326 close.
  * **6.3 execution order (designed 2026-08-08 from the call-site map
    `scratch/q5_6_3_flip_call_sites.md` — ~48 fixture sites / ~40
    files / 13+ artifacts / 3 capture scripts / 2 production
    builders).** Migrations FIRST, the refusal LAST — §6b: the wiring
    commit must follow every call site's migration; the only green
    ordering.
    1. Production MMS builders → `folded_product(4,8)` + a σ_y-parity
       gate (the ansatz `A + Bη` and its `ξ²B/r` source are ξ-even —
       VERIFY by mirror evaluation on the parent rule, never assume:
       the quotient computes GARBAGE, not zero, for excluded
       functions) + the 5 consuming MMS gate files re-posed (#229
       re-pinned at measured on the folded fixture; its floor FALLS
       at 6.4, not here; the n_phi sweep keeps PARENT counts — the
       factory's own parameter semantics).
    2. Mechanical fixture sweeps by directory (operators, sweep/core,
       sweep/curvilinear, eigenvalue, solve, primitives, mesh) — the
       no-artifact families.
    3. Snapshot families, each WITH its in-commit re-capture:
       regression snapshots + walk baselines (`_make_cyl` is shared
       with `_generate_walk_baselines` — coordinate, don't fork) +
       affine-carve + bc_extraction + T4 (script migrated then
       RE-RUN, never hand-edited — the `ec076008` stale-twin lesson)
       + `test_phase_c_crosscheck` name-keyed dicts/literals. ⚠ All
       6.3 captures are stamped INTERIM in provenance: `[M]` the
       [½,1] absorber still clamps 4/4 folded levels (8 entries in
       [1/5,½) per folded_product(4,8) — map §5), so every folded
       cylinder artifact moves AGAIN at 6.4. Kept deliberately: the
       un-clamp is a physics change and lands WITH its acceptance
       instruments (plan-authoring §8), and 6.4's re-capture delta
       then measures the absorber ALONE.
    4. The three non-mechanical redesigns (map §9.1/9.2/9.3):
       slot-coordination (its subject INVERTS — folded multi-level
       is the admitted cylinder norm; its refusal rows land with
       step 5), `test_cyl_direct_seed_fold` (subject unreachable
       post-flip — re-pose or retire WITH marker migration),
       `test_azimuthal_mirror_symmetry` (#326 — its fixtures refuse
       at the flip, so its handling cannot ride to 6.4 unchanged;
       decide at the file, re-poses recorded).
    5. THE FLIP: a raising helper beside
       `march_start_structure_per_level` + the CYLINDRICAL
       `_init_core` arm wiring (after `cylindrical_streaming`, whose
       structure-less guard fires first with the more specific
       message) + the refusal gate module (positive folded admission
       AND negative per family — vv #11; refusal is total by
       structure: full NODE_ALIGNED = edge-node-ONLY (⛔ this line
       read "edge-node+degenerate" until 2026-08-08 — REFUTED by the
       test-architect's §1.1 measurement: `(on_edge_node=True,
       degenerate=False)` on every product level; a gate transcribing
       the old parenthetical would pin a FALSE reason on every
       product row), full STAGGERED = mirror-tie degenerate, LS =
       hemisphere-tie degenerate at every order, slab-GL = fake
       single level on Σ; AND `[M]` leg 2c: all-tangential
       folded(n,2) is refused EARLIER still, by the trace tier's
       `build_omega_dot_n` — an order-discrimination row like
       slab-GL, never a MarchStart-refusal row, gate_design §8.1)
       + `_CYL_LS`/`_CYL_PRODUCT` refusal rows + the staleness sweep
       (`cylindrical_streaming` guard text + docstring recommend
       exactly the refused families; `MarchStart` bracketed instance
       lists; `radial_characteristic_levels` docstring).
    6. Docs + close: theory pages (curvilinear_one_group R12a
       trichotomy rows, angular_quadrature, structured_geometry),
       the LS double-fold capability issue FILED, V&V matrix, plan +
       batteries (tests/sn full not-slow re-measures the known-red
       set; tests/numerics WITH-slow owed here).
    Ruled at design (2026-08-08): admission stays **SNMesh-only**
    (map §9.6 — the carrying predicate is M-M sweep semantics;
    `cylindrical_streaming` is a shared geometry primitive and keeps
    its family-agnostic contract; its 6 direct test consumers stay).
  * **6.3 LANDED SO FAR (2026-08-08):** `479c2a19` fix-on-sight
    (loss_representation's R12a float spelling + "every cylinder
    carries none") · `ce6607f5` **leg 1** (MMS builders →
    folded_product; the σ_y-parity gate two-legged — parity on the
    fold's PARENT + the restriction leg BIT-EQUAL; `[M]` #229 folded
    floor 3.538e-3 → 6.782e-4 at n_φ 8→16, ratio 5.22× — most of the
    full-circle floor falls at the fold, before 6.4's honest τ) ·
    `f8ecb4f6` ⭐⭐ **ERR-078 FOUND AND FIXED** — the ψ½ march's
    `solve` DROPPED the outflow-row rhs (ERR-071's mechanism verbatim
    one system down; found by leg 2a's first coupled random-composite
    round-trip, `back = 0.0` exact on one slot per level; production
    bit-unchanged since physical q_out ≡ 0; THREE gates had pinned
    the wrong contract, incl. an apply∘solve bit-zero closure that
    was green BECAUSE of the bug — catalog entry ERR-078) ·
    `b1539468` **leg 2a** (13 operator-tier files; ⭐ the degenerate
    pure-azimuthal class SURVIVES the flip at folded n_φ ≡ 2 (mod 4)
    — φ = π/2 exact by roots-of-unity, μ_r = 0.0 BIT-EXACT, stronger
    than the retired 6e-17 NODE_ALIGNED fixtures; sphere_gl coupled
    row added; loss_transpose G4 retired + G2's mask exclusion
    retired; b3's `_sn_2d` is 2-D-Cartesian — the map overcounted).
    `384d62e4` **leg 2b** (6 sweep-core/top files; ⭐ the #209
    resonance re-posed to its measured truth — the old resonant node
    was the OLD LS8 seed, ALREADY retired by #337 (a latent red the
    inferred known-red count missed), and `[M]` the pole-cell a = 0
    resonance is UNREACHABLE at physical Σ_t on the admitted folded
    family (Σ* ≤ 0 for every inward ordinate, every probed mesh) —
    now pinned as the 6.4 TRIPWIRE that reds if the absorber
    retirement makes it reachable again; phase_c_gates' caller-less
    `prod_2x4` arm retired; unified :424/:543 verified Cartesian).
    **Remaining legs:** 2c sweep/curvilinear (⚠ the P(4,5) odd row →
    (4,6); + the @slow `test_si_cyl_20cell_nan_regression` fixture —
    same #209 family as 2b's re-pose) · 2d eigenvalue+solve+analytical
    · 2e primitives (⚠ mesh-simple carrier/split_spaces/split_leaves
    are FLIP-COMMIT family, not 2e — their cylinder fixtures are the
    R12a NEGATIVE instances whose assertions invert; they ride leg 5
    with slot-coordination) → leg 3 snapshots (interim-stamped) →
    ⚠ 2c decide-at-file: `test_coupled_pole_mu_level_invariant`'s
    P(4,5) row is SUBJECT-coupled (its docstring: odd n_φ has no
    x-mirror closure — that property IS the row's point), and the
    @slow `test_si_cyl_20cell_nan_regression.py` is the #209 family
    2b re-posed (its k_inf pin's fixture migrates coherently with
    2b's unreachability finding). ⚠ Renames of shared test fixtures
    owe the L20 importer sweep — `[M]` the `_cyl_product` rename
    broke `test_inverse_adjoint_coherence`'s cross-file import,
    caught only by full COLLECTION (targeted per-file runs never
    collect the importer; the V&V-matrix regen inside sphinx -W is
    the gate that surfaced it). **Batteries at this checkpoint
    (mid-6.3, HEAD after the coherence fix):** targeted per-leg
    batteries green (leg 2a 461 rows / leg 2b 75 / MMS 40 / the
    ERR-078 production-regression slice 140); full COLLECTION 9206
    tests / 0 errors; pyright ratchet green; xref checker 0 dead;
    sphinx -W clean (re-run after the coherence fix). NOT run:
    tests/sn full not-slow (owed at 6.3 close — re-measures the
    known-red set, already known wrong by ≥1: the latent #209
    resonance red) and tests/numerics WITH-slow (owed at 6.3 close,
    as before). →
    leg 4 redesigns → leg 5 THE FLIP — with the test-architect's
    three corrections (`scratch/q5_6_3_gate_design.md`; read its §6
    SIX user questions before leg 5): (i) the `folded_by`-provenance
    pincer — a hand-built tag-less folded-equal rule MUST construct,
    a `quotient(NODE_ALIGNED product)` (tag present, non-carrying)
    MUST refuse; (ii) `[M]` NODE_ALIGNED product is edge-node-ONLY
    (degenerate=False) — the two facts split the negatives, so
    per-conjunct mutations each red exactly half; (iii) split
    `non_carrying_levels()` out of the raiser or the ∀-quantifier
    ships ungated; slab GL is the ORDER-discrimination row (the
    structure-less guard fires first), never a refusal row → leg 6.
  * **6.3 LANDED, SECOND HALF (2026-08-08 continuation — legs
    2c/2d/2e/3/4 ALL LANDED; the "Remaining legs" block above is now
    HISTORY):** `2befb14d` **leg 2c** (7 curvilinear files; ⭐ `[M]`
    `test_both_quadratures_agree` RETIRED — its homogeneous
    REFLECTIVE fixture is k_inf-blind, five rules all measured keff =
    1.500000000000 exactly, so the 1e-6 cross-family agreement was
    Mode-12-degenerate all along AND the subject dissolves post-flip;
    ⭐ `[M]` folded(n,2) is INADMISSIBLE one tier BELOW the
    march-start predicate — the all-tangential rule (every μ_x ==
    0.0 bit-exact) is refused by the trace tier's
    `build_omega_dot_n` while the MarchStart facts say carrying —
    gate_design §8.1: an order-discrimination row like slab-GL,
    never a MarchStart-refusal row, and "carrying ⟹ admitted" must
    never be claimed (admission = the tier guards' CONJUNCTION in
    firing order); the P(4,5) pre-march refusal pin is flip-commit
    family (§8.2); 282's cyl row PROMOTED to the coupled branch
    incl. the ERR-078 outflow row; the #206 battery stays honestly
    xfail 27/27 on folded; apply_matvec's n_phi grid {2,4}→{4,6}
    keeps the surviving bit-exact μ_r=0 degenerate class in the L0
    batteries) · `569a75f2` **leg 2d** (5 files; keff_curvilinear's
    two product-vs-LS param lists → two folded SPLITS 4x8/8x4 for
    level-bookkeeping diversity; `[M]` the touched @slow cylinder
    set 12 passed in 1:13:35 — sweep vs trajectory_resolvent in the
    3% budget, twin-path 1e-5, refinement 20/40/80 both paths,
    SI≡Krylov eigenpair, CP het) · `4afc8313` **leg 2e** (5
    primitives files; α-dome over five folded splits incl. the
    degenerate 4x6; contamination β machine-zero on folded;
    quadrature-LEVEL family tests keep their LS/product rows — they
    test those rules AS RULES, which stay constructible) ·
    `c39b7d44` **leg 3** (23 files — every artifact family
    re-captured on the fold, stamped INTERIM for the 6.4 absorber
    double-move: the 3 regression snapshots RENAMED to honest folded
    names + regenerated, `[M]` new keffs 1.5 exact / 1.5−2ulp /
    3-region 1.2302082296342958 (was 1.2284281074857448 on LS4,
    0.14%); walk_matvec_cyl_2g on folded(4,8) with `[M]`
    slab/sphere/cart2d re-captures BYTE-IDENTICAL; bc_extraction 6 +
    affine-carve 3 CYL baselines re-captured via --capture-baseline
    (⚠ that flag re-captures EXISTING baselines too — verified no
    non-CYL bytes moved); T4 script re-run, `pre_t4_walltime.json`
    restored to HEAD (it measures an UNCHANGED slab fixture; the
    fresh numbers were warm-machine noise); phase_c keys renamed +
    literals re-measured; `[M]` the 4 phase_c @slow cylinder
    cross-checks passed in 55:46) · `7c0db72f` **leg 4**
    (slot-coordination's subject INVERTED — its own "if a future
    quadrature ever carries >1 level" tripwire fired; folded rows
    assert the multi-level carrier as the admitted norm; NEW
    per-(p_idx,sign) disjoint-marker gate with LEVEL-ASYMMETRIC
    markers 10p±1 defeating the bit-palindrome Mode-12 hazard;
    product/LS rows keep () and INVERT into step 5's refusal
    negatives; the other two redesigns RULED in gate_design §8.4 —
    test_cyl_direct_seed_fold DEFERRED to leg 5 coupled to §6 Q3
    (its six gates ARE Q3's "six tests"; green today),
    azimuthal_mirror #326 rides the flip as retire-with-tombstone +
    quadrature-level survivors, ⛔ superseding the old "re-poses at
    6.4" schedule). **Batteries at this checkpoint (HEAD `7c0db72f`
    + the matrix regen):** all touched-file batteries green (2c
    92+27xf not-slow + 3 @slow 10:44 · 2d 44 + 12 @slow 1:13:35 ·
    2e 122 · leg-3 consumers 111 + 4 @slow 55:46 · leg-4 45); full
    COLLECTION 9221 tests / 0 errors; pyright ratchet green; xref 0
    dead; sphinx -W clean + V&V matrix regenerated. NOT run:
    tests/sn full not-slow + tests/numerics WITH-slow (owed at 6.3
    close, as before). ⛔ the "STOPPED AT THE USER BOUNDARY" state
    that stood here RESOLVED same-day — the user answered §6 and
    ruled; superseded by the leg-5 landing record below.
  * **LEG 5 LANDED `1689faf4` (2026-08-08) — THE FLIP IS IN.** The
    user's §6 answers: Q1 YES (below-admission sites migrated,
    `1f220c41`; the map missed the two kinf standoff files — found
    by the tree-wide straggler sweep); Q2 {S4, S18} adopted as
    principled (mechanism representative + the #337-pinned range
    boundary); ⭐ Q3 RETIRE the #280 2.5b fold — ruled after the
    disambiguated seed briefing + `[M]` the foreign-rule check
    (`scratch/probe_foreign_gl_arc_rule.py`: a hand-built GL-in-φ
    arc rule — non-uniform azimuthal nodes AND weights, no factory,
    no `folded_by` tag — classifies carrying, admits, and the
    route-(a) coupled COLD solve is a true inverse at 4.5e-16 with
    k_inf exact; the user's GL×Gauss-Chebyshev example IS the
    shipped staggered rule read in the angle variable). Landed:
    `non_carrying_levels` + `assert_carrying_quadrature` (public,
    Q4) wired AFTER `cylindrical_streaming`; the fold + transpose
    twin + every `is_seed_ord` branch DELETED (behavior-neutral on
    every admitted input — the fold never fired on carrying
    meshes); the 16-row admission module with the architect's three
    corrections + the §8.1 trace-tier row; `[M]` the 10-mutation
    battery 10/10 with MA2 first — the control caught its own
    broken parser TWICE (in-process `-q` emits no FAILED lines;
    forced ANSI broke the match — vv anti-#17 verbatim); all
    riders in the one commit (mesh-simple trio re-posed,
    slot-coordination's product/LS rows → the admission module,
    psi_half's four non-carrying-cylinder controls → slab-only +
    a folded 2×2 presence leg, the P(4,5) pin re-posed onto the
    admission fragment through the public funnel, azimuthal_mirror
    tombstoned per §8.4, direct_seed_fold deleted). `[M]` the wide
    slice's 20 reds ALL accounted: 14 = one straggler fixture
    (fixed, 95 green with both kinf importers), 5 = the 2-D octant
    family PRE-EXISTING (identical reds at `1f220c41` in a
    throwaway worktree — the #51 known-red ledger, "known wrong by
    ≥1", is now measured wrong by +5:
    `test_2d_octant_sweep_equivalence` rows 01–05), 1 = the known
    spherical-inward diamond. Collection 9197/0. The
    Gauss-Lobatto × admission interaction FILED as #338 (the
    sphere-side node-aligned case needs its own decision).
  * **LEG 6 LANDED `6837a429` (2026-08-08) — 6.3 IS CLOSED.** The
    theory corpus follows the flip, retirement-by-tense: the four #280
    2.5b fold references go past-tense; `curvilinear_one_group`'s "only
    the sphere carries the block; every cylinder inlines the 2-point
    extrapolation" was the **exact inverse** of the tree and is now
    "every ADMITTED level carries"; the R12a trichotomy's first two rows
    are marked refused-at-admission (they survive as quadrature-level
    classifications, not constructible meshes); the circle-vs-interval
    conclusion inverted with it — the fold **renounces** the circle's
    free edge-inclusion deliberately, because that freeness came bundled
    with Σ-on-the-mirror, the #326 η-ties, and the τ=0 division block, so
    the folded cylinder JOINS the sphere in paying one seed per level;
    `angular_quadrature` gained the missing **Folded Product** section
    (the admitted family) + an LS warning (`level_indices` is necessary,
    NOT sufficient); the Gauss-Lobatto study's pole node re-framed as the
    SAME non-carrying class, arising sphere-side (#338);
    `structured_geometry`'s "the absorption must survive until the fold
    wiring" ⛔ REFUTED in place — the wiring landed and dissolved that
    reason; what keeps the clip alive to 6.4 is **baseline discipline,
    not structure**. ⭐ `[M]` the CONCEPT grep caught what seven symbol
    greps missed: `sn/index.rst:516` opened *"For cylindrical 1-D radial
    sweeps with a product or level-symmetric quadrature"* — an
    instruction naming exactly the two refused families (the
    `coding-standards` "grep the CONCEPT, not only the symbol" clause,
    earning its keep again).
    `[M]` **both owed batteries measured**: tests/sn FULL not-slow
    **2882 passed / 10 failed** (30:59) and tests/numerics WITH-slow
    **2075 passed / 0 failed** (16:16). Nine of the ten reds are the
    known ledger (3 affine sha256 = #333 · 5 octant · 1 spherical
    diamond). ⚠ **the tenth is NEW and the full battery is what found
    it** — `test_d3_admission.py::test_d3_pure_absorber_per_ordinate_psi_exact`,
    `[M]` max rel **3.287e-10** vs its `rtol=1e-10`, and `[M]` **bit-
    identical at every `inner_tol` from 1e-9 to 1e-15** ⟹ NOT a
    convergence artifact; it lives on exactly the 8 |μ_x|-dominant
    ordinates (the |μ_x|=0.3500 eight are exact at 1e-15) and only in the
    low-Σt group. `git bisect run` ⟹ first-bad **`59bb38a0`** (#337's
    moment-matched LS seed: S4 |cos| moved 0.4082/0.8165 → 0.3500/0.8689).
    ⛔ REFUTED, with reasons: a broken/many-to-one mirror map (ERR-073
    class) — the cosine sets are bit-exactly negation-closed, the x-mirror
    is an injective involution, partner weights differ by 0.0; and
    source-normalisation drift — Σw − 4π = 0.0 exactly. Under
    ✅ **RESOLVED `11e78430` — and it was NEITHER hypothesis.** The
    all-reflective pure absorber needs `[M]` **1631** sweeps against a
    default `max_inner=1000`; the solve returns the 999th iterate flagged
    `history.converged = False`, and the gate **never read the flag**. The
    exact uniform field's residual is 1.06e-15 (the discretization is
    fine) and a converged solve delivers 2.79e-14 — so `rtol=1e-10` was
    the RIGHT number all along and is untouched; the fix grants the budget
    and asserts convergence before reading any value. ⭐⭐⭐ **the ruling
    banked three paragraphs above is ⛔ REFUTED by this**: "bit-identical
    across a 6-decade tolerance sweep ⟹ not a convergence artifact" is
    exactly BACKWARDS — the residual at k=999 (1.185e-09) sat above every
    tol tested, so all four runs hit the same cap and returned the same
    bytes. **Tolerance-insensitivity means the tolerance never BOUND; read
    `n_inner` against the cap before concluding "floor".** The same trap
    dressed the error MAP as a bias — the 8 |μ_x|-dominant ordinates are a
    decaying MODE (shape cosine ≥ 0.9928 across `max_inner` 600→1400 at
    constant ρ = 0.9853). ⭐⭐ **#337 EXONERATED**: the gate had ALWAYS
    ridden a truncated exit, green pre-#337 only because the error landed
    at 2.05e-11, inside `rtol=1e-10` by 5×; the ratified 6.4 % cosine rise
    lifted the budget 1369→1631 and the error 16×, across the line. A
    correct change stopped a test-design defect from hiding. Its docstring
    premise WAS the defect — "c=0 needs no iteration" kills the SCATTERING
    iteration, not the REFLECTIVE-BOUNDARY one. ⭐⭐ an audit plugin over
    four reflective modules (52 solves, exactly 2 unconverged) caught the
    SIBLING as a **latent false green** —
    `test_d3_scattering_infinite_medium_matches_multigroup_balance` rode
    the same truncation under a looser `rtol=1e-7`; both now carry the
    budget and the flag check. `[M]` budget laws pinned in
    `derivations/diagnostics/diag_d3_absorber_02`: ~an order per dimension
    (d=1 32 · d=2 258 · d=3 1631) and `Σ_t · n_inner` INVARIANT. Filed,
    not fixed (public-API contract, needs the wide slice): **#340** the
    best-effort exit is indistinguishable from a certified one, **#341**
    boundary G-S ~2× SLOWER than Jacobi at d=3 (its documented rate gain
    was measured at d=2 and INVERTS), **#342** `solve_sn` hardcodes
    `converged=True` while its adjoint twin is honest. Report:
    `scratch/d3_absorber_diagnosis.md`. ⟹ the #51 ledger returns to NINE. Gates at this stamp:
    sphinx -W exit 0 on a FULL `-E` re-read, xref 0 dead / 940 files,
    pyright ratchet green, V&V matrix re-derived identical (9197 tests,
    336 equations). **#339 FILED** — the LS double-fold capability
    (quotient the vertical pair too; the σ_y fold alone does NOT repair
    LS, since the `|μ_z|` grouping keeps the hemisphere pair's η-tie
    through the ξ-fold).
    **▶ NEXT: 6.4 acceptance.** ⚠ existence-checked 2026-08-08 against
    the tree, and **one item of the list below is already DONE**: the
    three `xfail(strict=True)` rows (item 4) were retired-with-tombstone
    AT THE FLIP — `_XFAIL_326` now grep-hits *only* inside
    `test_azimuthal_mirror_symmetry.py`'s tombstone docstring, past-tense.
    Genuinely pending and verified live: the #229 floor (item 1); T3's α
    closed form (item 2 — `test_alpha_closed_form.py` still calls it "the
    Q5.6.4 T3 α" gate); τ ∉ {0,1} (item 3); and the **absorber
    retirement**, confirmed live at
    `orpheus/sn/sweep/pole_angular_closure.py:929`
    (`tau[m] = max(0.5, min(1.0, ...))`) — it owns its own re-baseline
    window because `[M]` it still moves 4/4 folded levels' τ at n_μ = 4.

  ✅ **Q5.6.4's PARTITION IS DECIDED AND STAYS** (2026-08-11, verdict in
  `.claude/plans/q64_tau_partition_memo.md` **§9bis.10 — read that section
  first, not the memo's §1–§9**). The carve is **BMC Eq. (43) verbatim** and
  the first time the cylinder satisfies it exactly: `[M]` the closure
  reproduces μ at the faces to `2e-15`, where chord+absorber reads `2.4e-3`
  and `τ ≡ ½` `2.4e-3`.

  ⛔ **Two frames were tried and BOTH refuted** — do not re-open either:
  candidate (3), "τ on α's recursion-defined edges", is **not a definition**
  (α is the arc-edge dome scaled by κ — an amplitude, and on every even-M rule
  one edge has no real position); and re-posing τ into ω is **BMC's own
  "diamond scheme"**, which their cylinder section (Eq. 53 + line 172) proves
  is leading-order only. Hébert §3.9.3 ships the plain diamond; BMC 2010 exists
  to show it insufficient. Do not cite Hébert against BMC.

  ⚠ **THE ONE RED IS REAL BUT IS NOT THE CLOSURE.** `[M]` the accuracy ranking
  of the four candidates is in perfect rank correlation with the angular
  recurrence's **transient error amplification** and INVERTED against closure
  fidelity. The retired `[½,1]` absorber bought damping by silently forfeiting
  Eq. (43); the carve removed a mask and exposed a pre-existing problem.
  ⟹ **the remaining work is the amplification / the march SEED
  (`MorelMontryAngularSweep.psi_half_seed`), not τ.**

  ⚠ Also `[M]`: the resolvent cross-check that produced the red is
  **REFERENCE-limited** — its floor is `≈3e-2` and refining the SN side moves
  it the WRONG way, so it cannot grade a closure past `n_φ = 16`.

  Everything below this banner is the ORIGINAL item plus its in-place
  refutations — accurate as history, but design from §9bis.10.

  1. **The #229 azimuthal floor** — ⛔ **THIS ITEM'S TEXT WAS STALE IN TWO WAYS;
     both corrected 2026-08-11 by an existence-check, original kept per
     `plan-authoring` §3.** It read: *"today flat at `≈1.9e-2` on the
     anisotropic curvilinear MMS with no convergence order. It must fall and
     recover an order. … This is the acceptance test."*

     ⛔ **(a) `#229` is CLOSED** (2026-06-13), and it never was the work item —
     it is the **measurement record** that named the floor and attributed it to
     half-angle-thread interpolation. A fresh session running `gh issue view
     229` sees CLOSED and concludes this item is done. ⟹ read every bare `#229`
     here as *"the floor measured in #229"*, never as an open issue.

     ⛔ **(b) `1.9e-2` is the PRE-6.3 configuration** (`[M]` 2026-06-13, FULL
     `NODE_ALIGNED` product rule) and quoting it as "today" is the
     §4 numbers-carry-their-configuration failure. The tree's own gate
     docstring (`tests/sn/verification/mms/test_curvilinear_aniso_convergence.py`)
     carries the live ladder:

     | when | fixture | `n_phi` 8→16 | ratio |
     |---|---|---|---|
     | `[M]` 2026-06-13 | full `NODE_ALIGNED` product | `1.90e-2 → 7.37e-3` (→ `3.10e-3` at 32) | 2.58× |
     | `[M]` 2026-08-08, at the 6.3 flip | `folded_product`, **absorber still live** | **`3.538e-3 → 6.782e-4`** | **5.22×** |

     ⟹ **the fold ALREADY took the floor down 5.4× at `n_phi=8` and steepened
     the azimuthal scaling**, before 6.4 touches anything. The floor is
     azimuthal, not polar: `[M]` holding `n_phi` and sweeping `n_mu` 4→8→16
     leaves it flat (1.90e-2, 1.91e-2, 1.91e-2).

     **What 6.4 actually owes, restated as an outcome:** re-measure this ladder
     under the **RETIRED `[½,1]` absorber** and show the floor falls again AND
     an order is recovered. The gate is structurally independent
     (manufactured-solution, not satisfiable by construction) and it is still
     the acceptance test — but the baseline it moves from is `3.538e-3`, not
     `1.9e-2`.

     ⛔ **REFUTED 2026-08-11 — retiring the absorber makes the floor WORSE, and
     the acceptance test as stated cannot be met.** `[M]` the aniso cylinder
     MMS, `nx=80`, volume-weighted L2, probe
     `$CLAUDE_JOB_DIR/tmp/q64_absorber_ladder.py` (the ABSORBER-LIVE row
     reproduces this gate's own docstring to every printed digit, so the
     instrument is sound):

     | τ convention | `n_phi`=8 | 16 | 32 |
     |---|---|---|---|
     | chord-midpoint + `[½,1]` absorber (HEAD production) | `3.5384e-3` | `6.7824e-4` | `2.4837e-4` |
     | chord-midpoint, **unclamped** (the naive retirement) | `6.2244e-3` | `2.3020e-3` | `6.0065e-4` |
     | arc ω half-angle edges, unclamped | `3.1503e-3` | `1.1326e-3` | `3.2611e-4` |
     | cumulative-WEIGHT edges (the SPHERE's convention) | `1.9252e-2` | **diverges (nan)** | **diverges (nan)** |
     | ⚠ CONTROL `τ ≡ ½` (angular diamond) | `3.4055e-3` | `3.7443e-4` | `1.3279e-4` |

     ⟹ the naive retirement is **1.8–3.4× worse** at every rung. And the
     CONTROL discriminates the tempting rescue: `τ ≡ ½` beats every
     *principled* convention everywhere, so the absorber is **not**
     "accidentally approximating the right thing" — it is partial movement
     toward angular diamond, and diamond is simply more accurate on this
     fixture. (Divergence of the sphere's convention was announced by the
     #340 N4.7 warning, unprompted: `hit max_inner=500 … last residual nan`.)

     ⛔ **THAT LAST SENTENCE IS ITSELF REFUTED, same day, by two further
     controls.** `τ ≡ ½` is not "angular diamond beating the principled
     conventions" — **it IS the principled convention**, and the table above
     was read on a contaminated fixture. Both corrections below.

     ⛔ **(i) `nx=80` DOES NOT ISOLATE THE ANGULAR FLOOR** — the table above,
     and this gate's whole "fixed-fine `nx=80`" framing, mix a spatial error
     into a claim about angular closure. `[M]` probe
     `$CLAUDE_JOB_DIR/tmp/q64_is_the_floor_spatial.py`, HEAD convention:

     | `n_phi` | nx=40 | 80 | 160 | 320 | spatial orders |
     |---|---|---|---|---|---|
     | 32 | `6.1936e-4` | `2.4837e-4` | `1.6358e-4` | `1.4449e-4` | 1.32 · 0.60 · 0.18 |
     | 64 | `5.3570e-4` | `1.5503e-4` | `6.1967e-5` | `4.0491e-5` | 1.79 · 1.32 · 0.61 |
     | 128 | `5.1562e-4` | `1.3397e-4` | `3.8772e-5` | `1.5488e-5` | 1.94 · 1.79 · 1.32 |

     At `n_phi=128` refining the mesh still drops the error 8.6× from `nx=80`,
     so the `≈1.3e-4` that every convention "saturated" to at `nx=80` is the
     **MESH**, not the closure. The `n_phi` 8→16 leg the shipped gate asserts
     IS angular-dominated (`3.5e-3 ≫ 1.3e-4`), so the gate is not wrong — but
     anything read at `n_phi ≥ 32, nx = 80` is mixed, and **"the cylinder has
     no O(h²) window / the floor dominates" needs re-examination**: at
     `nx=320` the angular error converges at a clean ~O(n_φ⁻²) with no flat
     floor in range.

     ⭐⭐ **(ii) THE CLEAN COMPARISON, and it names the real defect: THE
     CLOSURE INTERPOLATES IN THE WRONG VARIABLE.** `[M]` probe
     `$CLAUDE_JOB_DIR/tmp/q64_angular_isolated.py`, `nx=320` (spatial
     contribution ≤ `1.5e-5`, so 1e-4…3e-3 is angular):

     | τ convention | 8 | 16 | 32 | 64 | orders in `n_phi` |
     |---|---|---|---|---|---|
     | chord-midpoint + absorber (HEAD) | `3.5111e-3` | `5.8890e-4` | `1.4449e-4` | `4.0491e-5` | 2.58 · 2.03 · 1.84 |
     | chord-midpoint, unclamped | `6.2063e-3` | `2.2824e-3` | `5.7024e-4` | `1.4175e-4` | 1.44 · 2.00 · 2.01 |
     | arc ω half-angle edges (in η) | `3.1281e-3` | `1.1078e-3` | `2.8285e-4` | `7.1658e-5` | 1.50 · 1.97 · 1.98 |
     | **fractional position in ω** | `3.4258e-3` | `3.4907e-4` | **`3.9321e-5`** | **`9.1485e-6`** | 3.29 · 3.15 · 2.10 |

     **4.4–15× better at `n_phi=64`, with a better order, across the whole
     clean ladder — not a pre-asymptotic accident.** And it is not a tuned
     constant: `[M]` probe `$CLAUDE_JOB_DIR/tmp/q64_omega_fractional_position.py`

     | `n_phi` | `max｜τ_ω − ½｜` | `max｜τ_η,arc − ½｜` | `max｜τ_η,chord − ½｜` |
     |---|---|---|---|
     | 8 | **`0.000e+00`** | 0.240108 | 0.280455 |
     | 16 | `5.551e-16` | 0.247575 | 0.295311 |
     | 32 | `3.442e-15` | 0.249397 | 0.298839 |
     | 64 | `2.276e-15` | 0.249849 | 0.299711 |

     ⟹ **τ is the fractional position of the ordinate inside its own angular
     cell, measured in the variable the march marches in.** T22b already
     established that the azimuthal march is *a march in ω, arc by arc*; Q5.2's
     STAGGERED offset then puts every node exactly at the midpoint of its own
     ω-cell (`ω_j = (j+½)Δω`, cell `[jΔω,(j+1)Δω]`), so the fractional
     position in ω is `½` **bit-exactly** — while the same quantity measured
     in `η` is wrong by a permanent O(1) `0.25`–`0.30`, because `η = sinθ·cosω`
     is nonlinear. The `[½,1]` absorber was dragging the η-value crudely back
     toward the ω-value; that is the whole of what it was doing, and why
     retiring it alone regresses.

     ✅ Bonus consistency: `τ ≡ ½` keeps `τ_m + τ_{M−1−m} = 1` **bit-exactly**
     (`½+½`), so the landed T27 gate
     (`test_tau_arc_wellposedness.py::test_the_folded_tau_is_bounded_with_the_reversal_identity`,
     which asserts with `assert_array_equal`, NO epsilon) stays green. ⚠ the
     arc-in-η candidate would RED it at `1.11e-16` — a concrete reason to
     prefer the ω pose over the η-arc repair beyond the accuracy margin. But
     ⚠ that gate ALSO asserts `τ_raw ⊂ [1/5, 4/5]`, which becomes vacuous at
     `τ ≡ ½`; re-pose it, do not just watch it stay green (`vv-principles` #22).

     ⚠ **The fix is NOT "hardcode ½".** `½` is the *value on this rule family*
     (staggered fold ⟹ midpoint nodes). The DELIVERABLE is that the closure
     computes the fractional position in the march variable, so a
     non-midpoint azimuthal rule gets the right answer too — otherwise the
     campaign trades a wrong constant for a lucky one.

     ⚠ **UNSETTLED — do NOT touch the sphere on this evidence.** The sphere
     marches in `μ` with cumulative-weight edges. Whether ITS march variable
     should be `θ` rather than `μ` is the same question one dimension over,
     and nothing here measures it. (`[M]` the reverse transplant — the
     sphere's convention onto the cylinder — DIVERGES, so the two arms are not
     interchangeable in either direction.)

     ---

     ⛔⛔ **(iii) THE LITERATURE REFUTES (ii)'s CONCLUSION, AND THE MMS WAS
     MODE-12 BLIND THE WHOLE TIME.** Deliverable:
     `scratch/q64_tau_edge_convention_literature.md` (Bailey–Morel–Chang 2010
     + Hébert 2009 Ch.3 + Lathrop 2000, all local, load-bearing equations
     scan-verified). Four corrections, in descending importance:

     1. ⛔ **`τ ≡ ½` is NOT "the principled convention in disguise" — it is
        HÉBERT'S DIAMOND SCHEME, and it has the WRONG DIFFUSION LIMIT.**
        Lathrop's `δ` is τ affinely (`τ = (1+δ)/2`), and Lathrop states both
        halves two pages apart, both scan-verified: p.245 *"only with
        `μ_m = μ̄` (`δ=0`) is the truncation order `O(Δμ²)`"*; p.249→250
        *"`δ = 0` … would not have `c = ⅔` and hence would not give the
        correct diffusion limit"*. ⟹ **`τ=½` optimises exactly the quantity my
        L2 ladder measures and breaks exactly the quantity τ exists to fix.**
        The `nx=320` table above is not evidence for `τ=½`; it is a textbook
        Mode-12 result — the measured functional's stabiliser contains the
        error class. `[[feedback-nsweep-discriminator-closure-repose]]` said
        so in advance and I ran the ladder anyway.
     2. ⛔ **"α and τ should read ONE partition object" is REFUTED as a
        MEANS** (the goal survives; the mechanism does not). They *share* the
        partition but impose **different conditions on it**: α accumulates the
        FIRST moment, the τ edges the ZEROTH (Lathrop (22)/(25)/(26):
        `α = 1−μ_edge²` **iff `δ=0`**, else `α = 1−μ²+β` with
        `β_{m+1/2}−β_{m−1/2} = −δΔμ²/2`). **Unifying them naively silently
        re-imposes `τ=½`** — i.e. the elegant-looking single-source carve is
        the wrong-diffusion-limit scheme wearing a Pattern-2 badge. ⚠ Also:
        Hébert's `η_{p,q±1/2}` is a *constant-flux conservation recursion*
        (3.398, after Alcouffe & O'Dell), **not** a trig evaluation at a
        bisected ω — so T3's closed form is a THEOREM about our equispaced-ω
        rule, not the literature's definition of the partition.
     3. ⛔ **The code's own citation is wrong.** It cites Hébert Eqs.
        3.437/3.439 for the M-M τ; those are **`τ=½` extrapolations** carrying
        zero information about angular cell edges (they differ only in the
        SPATIAL half, inward vs outward sweep). The τ actually implemented is
        BMC-shaped. Same disease as #327 — a scheme wearing another paper's
        name. Fix the citation regardless of which way the closure lands.
     4. ✅ **The absorber is refuted by the source directly**: no reference
        prescribes ANY limiter, the admissible range is `[0,1]`, and BMC's own
        S₂ example gives `τ₁ = μ₁+1 = 1−1/√3 ≈ 0.4226 < ½` (Eq. 47,
        scan-verified). `max(0.5,·)` contradicts the published scheme and
        re-introduces exactly the `β` contamination τ exists to remove. **The
        retirement is CORRECT; only my justification for it was wrong.**

     ⭐⭐⭐ **(iv) THE ACTUAL BLOCKER, and it is bigger than τ: OUR AZIMUTHAL
     RULE IS INCOMPATIBLE WITH THE LITERATURE'S CLOSURE.** The literature
     convention is **(C) cumulative-weight edges in the radial cosine, for
     BOTH arms** — sphere BMC Eq. (12) `[VERIFIED-ON-SCAN p.5]`, cylinder BMC
     Eq. (52) per ξ-level with weights renormalised to `Σw̄ = 2√(1−ξ²)`
     `[p.10]`, τ by Eq. (74) ≡ (42). So **ORPHEUS's sphere is verbatim
     correct**, ORPHEUS's cylinder chord-midpoint is **nobody's convention**,
     and **the sphere/cylinder asymmetry is ORPHEUS's own, not the
     literature's**.

     But (C) cannot simply be adopted. `[M]` probe
     `$CLAUDE_JOB_DIR/tmp/q64_convention_C_applicability.py`, `folded_product(n_mu=4)`,
     level 0 — Lathrop's admissibility condition is that the ordinate lie
     INSIDE its own cell:

     | `n_phi` | M | min `τ_C` | max `τ_C` | ordinates OUTSIDE their own cell |
     |---|---|---|---|---|
     | 8 | 4 | 0.1522 | 0.8478 | 0 — admissible |
     | 16 | 8 | **−0.3259** | **1.3259** | **4/8** |
     | 32 | 16 | **−1.1841** | **2.1841** | **12/16** |
     | 64 | 32 | **−2.8552** | **3.8552** | **28/32** |

     **Why, structurally:** (C) presupposes an azimuthal rule whose per-level
     weights ARE the η-cell widths. Ours is equispaced in **ω** with EQUAL
     weights ⟹ equal η-cells — while the nodes `η = sinθ·cos ω` are
     cos-CLUSTERED toward `±sinθ`. Equal cells + clustered nodes ⟹ from
     `n_phi=16` most ordinates sit outside their own cell. (This also explains
     the `nan` in (i)'s table: not a probe bug, and NOT the missing
     renormalisation the literature memo hypothesised — the probe already
     renormalised to `2sinθ`.) The incompatibility WORSENS with `n_phi`.

     ⟹ **the fork is between the CLOSURE and the RULE, and it is the user's
     call** — the equispaced-ω rule is what the entire fold / `Σ=∅` / orbit /
     T22b–T25 machinery is built on, so "just switch to an η-partitioning
     azimuthal rule" is not a local edit. Options, none yet chosen:
     (a) keep the ω rule, accept the closure is NOT BMC's, and **name +
     verify it on its own terms**; (b) change the azimuthal rule so weights
     partition η, and pay the fold machinery's re-derivation; (c) a
     per-level re-weighting that reconciles the two.

     ⭐ **AND THE INSTRUMENT IS NOW KNOWN — STOP USING THE MMS FOR THIS.**
     BMC Eq. (75)'s **`β`** is a SOLVE-FREE oracle: a pure function of the
     quadrature and the α/edge recursions, `≈` round-off for a correct
     convention (BMC Table II: `−4.12E−16 … −1.55E−15`) and `O(10⁻²–10⁰)` for
     a wrong one. It sees the diffusion limit the L2 ladder is blind to. **The
     next measurement in this item is β, not another solve.**

     ⭐⭐⭐ **(v) THE RESOLUTION — THE FOLD ALREADY DOES THE M-M τ'S JOB, SO τ
     IS FREE AND GOES TO TRUNCATION ORDER.** (User steer, 2026-08-11: *"fix our
     own before changing to any other implementation … what we need might not
     be their closure or method but the CONCEPT they used to make it accurate.
     That is probably transferable."* — and the user then pointed at the
     in-tree instrument: `orpheus/derivations/discrete/sn/contamination.py`,
     written when M-M was implemented.)

     The transferable concept is **β = 0 ⟺ diffusion-limit consistency**
     (BMC Eq. 41 sphere / Eq. 75 cylinder). That module's own docstring states
     the design intent: *"Implementing Morel–Montry angular weights (Bailey
     Eq. 74) forces β = 0 exactly"*. β is a pure function of
     (nodes, weights, EDGES) — **no solve** — so "which edges?" becomes
     "which edges give β = 0?".

     ⭐ **It also already CONTAINS the correct diagnosis of (iv)**, at
     `contamination.py:64-66`: *"the weight-sum approach is wrong for
     cylindrical because weights are uniform in φ-space, not η-space."* The
     original implementer saw the incompatibility I re-measured in (iv). The
     DIAGNOSIS was right; the chord-midpoint SUBSTITUTE was never checked
     against β. So check it.

     `[M]` probe `$CLAUDE_JOB_DIR/tmp/q64_beta_per_convention.py`,
     `folded_product(n_mu=4)`, max |β| over levels:

     | `n_phi` | chord (HEAD) | arc ω | weight-sum | node-centred (τ≡½) |
     |---|---|---|---|---|
     | 8 | `5.55e-17` | `6.94e-18` | `1.11e-16` | `5.55e-17` |
     | 16 | `6.94e-18` | `5.55e-17` | `3.47e-18` | `5.55e-17` |
     | 32 | `1.39e-16` | `7.49e-16` | `2.78e-17` | `5.55e-17` |
     | 64 | `6.94e-17` | `2.78e-16` | `4.86e-16` | `2.36e-16` |

     **All four at round-off — including the weight-sum convention that
     `[M]` DIVERGES the solve at `n_phi≥16` with τ outside [0,1].** An
     instrument that certifies a divergent scheme is not measuring, so per
     `vv-principles` #17 the instrument was validated on a known-positive
     before any negative was trusted. `[M]` probe
     `$CLAUDE_JOB_DIR/tmp/q64_beta_is_blind.py`, `n_phi=16` level 0:

     | edge set | β |
     |---|---|
     | chord (production) | `+6.94e-18` |
     | GARBAGE — edges scaled 0.5× | `+3.47e-18` |
     | GARBAGE — edges CUBED | `+1.73e-18` |
     | GARBAGE — random, antisymmetrised | `−3.47e-18` |
     | one edge nudged (breaks antisymmetry) | **`−3.53e-03`** |

     ⟹ **β is a SYMMETRY IDENTITY on the folded arc.** The proof is two lines:
     the fold makes nodes antisymmetric (`[M]` `max|η + η[::-1]| = 0.000e+00`)
     and the α dome symmetric (`[M]` `2.78e-17`), so for ANY antisymmetric
     edge set `e`,
     `term_{M−1−m} = (−η_m)(α_m(−e_m) − α_{m+1}(−e_{m+1})) = −term_m`
     and the sum cancels pairwise. β sees **only** whether the edges are
     antisymmetric; random garbage passes.

     ⭐⭐ **The physics reading, and it is the item's answer.** β = 0 is the
     entire reason the M-M τ exists (BMC choose τ precisely to zero it). On
     the σ_y-folded arc, **the fold delivers β = 0 structurally, for free,
     independent of τ** — so on OUR rule the M-M τ has no job left to do, and
     τ is a FREE parameter that must be fixed by the NEXT criterion. That
     criterion is truncation order, and Lathrop p.245 names the optimum:
     *"only with `μ_m = μ̄` (`δ=0`) is the truncation order `O(Δμ²)`"* —
     i.e. **`τ = ½`**.

     ⟹ the two independent lines CONVERGE: the structural argument says τ is
     free and truncation order picks `τ=½`; the `nx=320` angular-isolated
     ladder in (ii) measured `τ≡½` best by **4.4–15×** with a better order.
     And **Lathrop's objection to `δ=0` is DISARMED HERE, not ignored**: it
     says `δ=0` "would not give the correct diffusion limit" — true on the
     UNFOLDED sphere, where β=0 must be BOUGHT with `δ≠0`. On the folded arc
     the symmetry buys β=0 for free, so `δ=0` costs nothing. ⟹ (ii)'s
     measurement was not Mode-12 blind after all; it was measuring the only
     criterion still live. **(iii).1's refutation of (ii) is itself now
     narrowed to the unfolded case** — the τ=½ ⟹ wrong-diffusion-limit
     implication does NOT transfer to a rule whose fold already zeroes β.

     ⛔ **(v) IS REFUTED — its premise (b) was false, and the τ≡½ conclusion
     is DEAD.** The literature follow-up resolved the tension as **two objects
     sharing a letter**, `[SCAN-VERIFIED]`, and the agent **reproduced BMC
     Table I to every printed digit** (step `7.698004e-01…`; diamond
     `2.06E-01, −3.57E-03, −4.57E-05, −1.21E-05`; MM `~1e-16`) — so the
     mechanism is verified, not argued:

     | | BMC Eq. 41/75 | Lathrop Eq. 25 |
     |---|---|---|
     | is | **one scalar** — the `J⁽²⁾` contamination coefficient | **a sequence** — α's pointwise defect |
     | zero iff | `τ =` Morel–Montry | `δ ≡ 0 ⟺ τ ≡ ½` |

     And BMC's `η_{m±½}` inside β are **not** cumulative-weight edges — they
     are the edges the scheme's **τ IMPLIES** via Eq. (43). So my β probe fed
     it the wrong quantity. (The blindness finding SURVIVES — the agent
     re-measured it independently and got round-off for (A),(B),(C) *and*
     `τ≡½`, blind by Lathrop's own oddness argument on a symmetric
     equispaced-ω level. ⟹ β stays unusable as a gate here, for the right
     reason.) Also: Lathrop's `δ=0 ⟹ c≠⅔` is about **the NODES `δ=0`
     implies** (midpoints), and `c ≡ Σ w_m η_m² = ⅔` is **quadrature-only —
     no edges, no τ**. Our nodes are fixed and satisfy it regardless, so
     "the fold frees τ" was answering a question β never asked.

     ✅ **(vi) THE RESOLUTION, RATIFIED — use (B), and P2 DETERMINES τ.**
     Q1's rule-agnostic predicate ladder (Lathrop Eqs. 53/54, `[SCAN p.12]`;
     the P₁ pair is recovered iff `3c/2 = 1`):

     | | condition | order |
     |---|---|---|
     | **P0** | `Σw = |range|`, `Σ w η = 0` (α closes) | — |
     | **P1** | `c ≡ Σ w_m η_m² = ⅔` | LEADING; **quadrature only** |
     | **P2** | `η_m = τ_m e_{m+½} + (1−τ_m) e_{m−½}` pointwise (BMC 43 = Lathrop 23) | FIRST |
     | **P3** | `η_m ∈ [e_{m−½}, e_{m+½}] ⟺ τ ∈ [0,1]` | well-posedness |
     | **P4** | `α_{m+½} = α_{m−½} − k w_m η_m`, seeded geometric (Hébert 3.397-9) | — |

     ⟹ **cumulative-weight edges are ONE SOLUTION, not the system** — BMC
     Eq. (52) is merely the statement that *in their* quadrature the weight
     equals the cell's η-measure. `[M]` on OUR arc the cell's η-measure is
     `2sinθ·sin(ω_m)·sin(Δω/2) ∝ sin ω_m`, **not constant** (ratio 0.30→1.53),
     which is why (C) violates P3 and worsens with refinement — the same
     4/8→28/32 ladder measured in (iv), now explained. **This is exactly the
     user's "take the CONCEPT, not their closure" instruction discharged.**

     **The decision:** the **arc ω half-angle partition (B)** — chosen because
     it satisfies **P3** for our nodes, not on BMC's authority — and **P2 then
     DETERMINES τ**, closed form ⛔ (this line first read *"`[M]` verified to
     `1.7e-16`"` — a SINGLE-FIXTURE reading published as if it were a
     machine-epsilon identity. `[M]` corrected: the agreement DEGRADES two
     orders under refinement — `1.1e-16 / 2.2e-16 / 7.8e-16 / 7.4e-15 /
     2.3e-14` at `n_φ = 4/8/16/32/64`, max over all four levels. The
     shipped gate already knew, at `atol=1e-13`; the plan did not. Same
     defect class as §4's configuration rule, one level down: a number
     measured at one order, quoted without its order):

     > **`τ_m = ½ + ½ cot(ω_m) · tan(Δω/4)`**

     giving `[0.259892, 0.458804, 0.541196, 0.740108]` at `n_phi=8` — the arc
     values measured in (ii). `[M]` **(A) is (B) in disguise**: its interior
     edges are `cos(Δω/2) ×` (B)'s to `1e-16` with the END CELLS STRETCHED —
     precisely the ~17.5 % ω-width defect of (iii)'s table. Prefer (B).

     ⭐⭐ **THE GATE — ν-CLOSURE, and it is the instrument this item was
     missing.** The march implied BY τ (`ν_{½}=−sinθ`,
     `ν_{m+½} = (η_m − (1−τ_m)ν_{m−½})/τ_m`) must land on `+sinθ`. Solve-free,
     no MMS, and it separates a DERIVED τ from a fabricated one. `[M]` probe
     `$CLAUDE_JOB_DIR/tmp/q64_verify_verdict.py`, `ν/sinθ` at close:

     | `n_phi` | arc | chord | **clamped (HEAD)** | **τ ≡ ½** |
     |---|---|---|---|---|
     | 8 | `1.000000` | `1.000000` | **`1.016389`** | **`1.164784`** |
     | 16 | `1.000000` | `1.000000` | `1.001930` | `1.039182` |
     | 32 | `1.000000` | `1.000000` | `1.000238` | `1.009677` |
     | 64 | `1.000000` | `1.000000` | `1.000030` | `1.002412` |

     ⟹ **the `[½,1]` absorber and `τ≡½` correspond to NO partition of the
     level** — their implied march overshoots the level's own endpoint by
     1.6 % and 16.5 %. That condemns the absorber on principle with zero
     reference to any MMS, and it is the honest replacement for the refuted
     "the floor falls" acceptance test. **Ship this as 6.4's gate.**

     ⚠ **ACCEPT THE ACCURACY TRADE, DO NOT RE-LITIGATE IT.** `[M]` at
     `nx=320`, (B) is BETTER at `n_phi=8` (`3.128e-3` vs `3.511e-3`) and
     ~1.8–2× WORSE at 16/32/64 than the clamped baseline. Per
     `[[feedback-principled-over-bit-identical]]` and
     `[[feedback-nsweep-discriminator-closure-repose]]` (*principled ≠ more
     accurate*), the scheme that satisfies P2/P3 wins over the one with the
     smaller number on one manufactured fixture. **Re-baseline; document the
     ladder honestly in the gate's docstring; do NOT keep the clamp for its
     L2.**

     ⚠ **Retirement is now doubly-motivated, and `τ ∈ [0,1]` gets PROMOTED to
     P3 as a named predicate** — the existing raising guard IS P3; say so.
     `[M]` at S₈ Gauss–Legendre **four of eight** MM τ are below ½, so the
     `[½,1]` box was never the admissible range in either arm.

     ⚠ **Still open / not ours yet:** whether a *third* floated partition
     beats (B) is an **original-work gap, unsettled anywhere in the corpus**;
     **Morel & Montry (1984) TTSP 13(5):615** remains the one unacquired
     primary source. And a printed-constant discrepancy: Lathrop Eq. (26)
     gives `−δΔμ²/2` where its own Eqs. (19)+(23) yield `−δΔμ²` (structure
     unaffected).

     ⚠ **And β must NOT be shipped as the gate for this** — on the folded arc
     it is green for random garbage. If a gate cites β for the cylinder, it
     is a coverage claim an audit will trust and it cannot fail
     (`coding-standards`' tautological-gate rule). Either gate β on the
     UNFOLDED parent (where antisymmetry does not cancel it) or gate the
     antisymmetry itself and say that is what it tests.

     ⚠ `contamination.py`'s sphere arm has **API bit-rot** (#347 family):
     `contamination_beta(gauss_legendre_on_mu(N), "spherical")` raises
     `AttributeError: 'DiscreteMeasure' object has no attribute 'mu_x'`. The
     sphere control could not be run. Fix before citing the module anywhere.

     ⚠ **Two literature hazards recorded** (`scratch/q64_tau_edge_convention_literature.md`):
     BMC Eq. (52) line 1 is a **typo in the printed journal** (`μ_{m+1/2} =
     μ_{m+1/2} + w̄`; RHS must be `m−1/2`, proven by the correctly-printed
     sphere twin Eq. 12 + non-vacuity) — a scan match proves the OCR faithful,
     not the equation correct. And **`μ`/`η` are SWAPPED vs ORPHEUS**: BMC's
     and Hébert's `μ` is our `η`. Morel & Montry (1984) TTSP 13(5):615 and
     Alcouffe & O'Dell are **not local** — MM-1984 is the one remaining place
     a statement about τ's admissible range could live.

     ⚠ **And the MMS is the WRONG INSTRUMENT for this choice anyway** —
     `[[feedback-nsweep-discriminator-closure-repose]]`: *principled ≠ more
     accurate*, and a τ change is a CLOSURE re-pose, so an MMS-floor
     comparison measures error cancellation, not correctness. Do not let this
     table pick a convention; it is here to refute a claim, not to found one.

     ⭐⭐ **THE STRUCTURAL FINDING, which needs no accuracy argument.** α and τ
     both reference ONE physical object — the boundary between azimuthal cell
     `m` and `m+1` — and they DERIVE IT INDEPENDENTLY, in disagreement
     (Cardinal Rule 2 / elegance Pattern 2, two spellings of one concept):

     * **α** (gated by T3, `test_alpha_closed_form.py:113`) puts the boundary
       at the real half-angle `ω_{k−1/2} = (k−½)Δω`, so its radial cosine is
       `sinθ·cos(ω_{k−1/2})`.
     * **τ** (`pole_angular_closure.py:840-844`) puts it at the CHORD midpoint
       `(η_m + η_{m+1})/2`. `cos` is nonlinear, so these are different numbers.

     `[M]` probe `$CLAUDE_JOB_DIR/tmp/q64_cell_partition_disagreement.py`,
     `folded_product(n_mu=4)`, level 0:

     | `n_phi` | max rel η disagreement | implied ω-width spread of τ's partition | quadrature's own weight spread |
     |---|---|---|---|
     | 8 | `5.3825e-2` | **18.71 %** | `0.00e+00` |
     | 16 | `1.7752e-2` | **17.59 %** | `0.00e+00` |
     | 32 | `4.7227e-3` | **17.48 %** | `0.00e+00` |
     | 64 | `1.1987e-3` | **17.46 %** | `0.00e+00` |

     The quadrature's cells are bit-exactly equal (45.000° each at `n_phi=8`,
     weight spread `0.0`); τ's partition implies **49.211°, 40.789°, 40.789°,
     49.211°**. The η disagreement vanishes as `Δω→0` but **the ω-width spread
     does NOT** — it converges to ≈17.45 %, an O(1) inconsistency at every
     azimuthal order. The algebra: for interior edges the chord midpoint is
     exactly `cos(Δω/2) ×` the arc edge (the κ prefactor's sibling), while the
     ENDPOINTS are pinned at `∓sinθ` and left unscaled, so the outermost cells
     stretch to absorb the shrink. **That is what the `[½,1]` absorber was
     compensating for**: the chord partition drives `τ_0` to `0.2195` where the
     arc value is `0.2599`, and clamping to `½` crudely shoves it back.

     ⚠ **Fixing the partition alone does NOT recover the floor** (`[M]` the arc
     row above is still worse than production at `n_phi` 16/32), so the floor
     has a second, unattributed contributor. Do not close 6.4 on the partition
     repair and claim the floor.

     **⟹ 6.4's goal, RESTATED IN DOMAIN TERMS** (`plan-authoring` §5, §1): *a
     polar level's azimuthal cells are ONE partition, derived once from the
     quadrature that defines them, and every coefficient that references a cell
     boundary reads it from there.* The absorber retirement is a CONSEQUENCE of
     that, not a step of its own. **Open, and blocking the means (not the
     goal):** which partition the literature prescribes — the sphere uses
     cumulative-WEIGHT edges and the cylinder chord-midpoint, and it is not yet
     established whether that asymmetry is the literature's or ORPHEUS's own
     (dispatched; deliverable `scratch/q64_tau_edge_convention_literature.md`).
  2. **T3's α closed form** holds exactly on the arc
     (`α_k = −w_gl·κ·[ξ(ω_{k−1/2}) − ξ(ω_{−1/2})]`, `κ = Δω/(2 sin(Δω/2))`).
     ⚠ Do NOT gate a bare `α == −ξ`: `κ` is 2.6 % off 1 at `n_φ = 8`.
  3. `τ_raw` leaves `{0,1}` (T22b's mechanism).
  4. The three `xfail(strict=True)` rows are **re-posed or retired with reasons
     recorded** — never silently XPASSed, and never re-posed as "ψ even in ξ on
     the quotient" (tautological, `vv-principles` Mode 8).

  **Retirement list** (a numbered deliverable, per `coding-standards.md`):
  the `[½,1]` **absorber** (RULED at Q5.5 `0ec701d4`: the `[0,1]` membership
  half is ALREADY promoted to the raising guard; the absorption half retires
  HERE, call-site-coupled to the wiring — production NODE_ALIGNED cylinders
  divide by zero unclamped until the wired fold makes every consumed
  cylinder rule an arc) · ~~`LevelStructure.fiber()`~~ ✅ retired at
  Q5.3 `65ff5bb0` ·
  ~~the τ_raw trichotomy as the R12a predicate~~ ✅ retired at Q5.4 `679cba61`
  (the trichotomy survives as a gated theorem; the predicate is `MarchStart`) ·
  `derivations/discrete/sn/contamination.py:45-90` (the `eta_edge` twin, whose
  docstring writes the double cover down as a feature) · `_XFAIL_326` + its three
  rows · `_test_helpers.product_level_ordering` / `PRODUCT_LEVEL_ORDERINGS` (the
  adjudication they existed for is closed by T22).

  **✅ BOTH RULED, 2026-08-02 — do not re-litigate:**
  * **The fold reaches the SOLVE.** A cylinder gets 16 ordinates where it had 32,
    and every cylindrical snapshot moves. User: *"we had already settled on the
    quotient at the very start of the quadrature interlude … There is ZERO point
    in having a snapshot that is not reflecting our best understanding of
    correctness, when correctness is our first cardinal rule. Remaking a snapshot
    is expected when our understanding of correctness improved in a principled
    way. It's just a matter of WHEN we will rebuild the snapshot, not IF."*
    ⟹ **a snapshot is never a reason to decline a principled correction.** Plan
    the re-baseline as scheduling, not as a cost to weigh against the fix.
  * ✅ **Q5.0.2 — DONE.** `Z2` **RETIRED**, `Mirror(axis)` minted. Blast-radius
    map `scratch/q5_mirror_plane_blast_radius.md` (its §6 checklist held).
    `tests/numerics` 1753 + 9 new; `sphinx -W` clean; pyright 1 (#288).

    ⭐ **Why retire rather than add.** `Mirror('z') == Z2` was `False` while
    `contains` said each held the other — **two unequal spellings of one
    group**. And `Z2` meant two DIFFERENT things by node shape: the 3-D arm
    tested σ_z, the 1-D arm tested plane-free `x → −x`. The tree's canonical
    embedding of a polar marginal is `(μ, 0, 0)`, so the slab mirror `μ → −μ`
    **is σ_x** — a different matrix. `[M]` On an ASYMMETRIC μ-set the two arms
    disagreed, the 3-D one CERTIFYING a set that violates `μ → −μ` (ERR-072
    family, the dangerous direction). And the false docstring — *"any single
    reflection works; the choice is convention"* — is falsified by a SHIPPED
    rule: `product(4, 3)` is σ_z-closed and NOT σ_x-closed.

    ⭐ **The 1-D arm is now DERIVED, not declared.** Under `(μ, 0, 0)`, σ_y and
    σ_z fix every node POINTWISE (so they hold trivially) and σ_x carries the
    whole content. That is what makes the two arms answer one question —
    verified with an asymmetric discriminating leg AND a symmetric control (on
    a symmetric set every plane is True and the arms agree even while asking
    different questions, which is exactly how this hid).

    `[M]` The map's predictions all held: `Dnh(n) ⊇ σ_x` iff `n` even, `⊇ σ_z`
    always; `Oh/O3/Dinfh ⊇ σ_a`; `SO2/SO3 ⊉ σ_a` — and the last two are now
    right BY THE ARM rather than by the bare `return False` fallthrough. Three
    of the five old `Z2` lattice edges were already dead code; the family owes
    exactly two facts, both on `_contains` arms since `_NAMED_LATTICE` is typed
    enum-to-enum.

    Re-posed `test_dnh_reflection_in_dnh`, whose prose *"a single reflection
    sits inside every D_nh"* was a false generalisation the parameter-free tag
    made unfalsifiable — two of its five orders answer False once named.

    ⛔ **Method note, and it cost a full mutation round.** The first mutation
    harness monkeypatched the PARENT process and ran pytest in a SUBPROCESS,
    which re-imports a clean module — all four mutations read GREEN. The
    positive controls were real and proved the mutation bit **in the wrong
    process**. Redone with source-level mutation (copy to tmp, mutate, restore
    — never `git checkout`): M1/M3/M4 red immediately, and M2/M5 red only after
    removing a `-k` filter that had deselected the very tests written to catch
    them. **Five for five once the mutation actually reached the tests.**
- **Q6 — Enumeration & naming** (T12d). Extend the EXISTING `QuadratureSpec` /
  `quadrature_registry` — do not mint a mechanism. Adds the **family** grouping
  (= the construction axis, a closed sum type), a **`ParameterSpec`** replacing
  `parameters: dict[str, type]` (it cannot express Lebedev's discrete admissible
  set vs `n ≥ 1` vs Jacobi's `> −1`, and the domain is ALREADY duplicated inside
  `degree_of_exactness_for`), a display label distinct from the identifier, and
  **registers the Gauss families** — which finally gives `gauss_chebyshev` /
  `gauss_laguerre` / `equispaced` a consumer. Use a **registering constructor**
  so the constant and the entry are ONE expression.
- **Q7 — `SubgroupOfO3` residual groups + the `None` conflation** (gaps 4, 8, 6).
- **Q8 — #327**: `level_symmetric` either implements the Carlson–Lathrop
  moment-matched weights it claims, or is renamed to what it is. Under T12b this
  becomes structural: a rule built from its generating measure cannot over-claim.
- **Q9 — Retire the campaign's own scaffolding**; docs; the corpus has **zero**
  mention of T4 or #326 today. Plus the §6 close-out list.

---

## 6. Open close-out items

- **`align_to(K) -> Rotation | None`** — the subconjugacy certificate (T15). Named
  and derived, NOT built. Without it the re-orientation freedom is unspellable and
  a caller that rotates a rule has no way to say so.
- ~~`O2` means `C_∞h`~~ — **RESOLVED**, see T15d: rename to `D_∞h` and re-realize
  as the full cylindrical group. Closes the audit's S-7.
- **Vectorise `_orbit_closure`** — it is an O(N²) Python loop per operator, and
  the graph walk (T17) now calls it many times. Cost is real and measured:
  `test_symmetry.py` 5s → 57s, and a brute-force scan over every candidate of a
  110-node Lebedev rule took **150 s** before the gate's rule set was trimmed. The
  fix is one `(N,N,3)` broadcast + `argmin` per operator.
- **Remaining `symmetry.py` audit findings** (`scratch/q1q2_symmetry_audit.md`),
  not yet addressed: the 1-D path returns `True` for all ten tags
  (`SO3`/`Cn` unconditionally, though `R_x(π)` maps `μ→−μ`); the node window is
  `100·atol` while the weight window is `1·atol` (and a third constant `10·atol`
  in the 1-D helper), all absolute; `_rotation_x`/`_rotation_y` are dead code and
  `SPACE_SPHERE` is an unused import.
- **Mint the ERR entry for #326** and attach `catches(...)` to the mirror gates.
  Deliberately not minted yet — a parallel investigation was live on the same
  issue and an ID collision is worse than the gap.
- **`rules_product.py:40-44` cites the RETRACTED** Bailey-Adams-Yang-Zika 2009
  paper that `curvilinear_one_group.rst` explicitly corrects, attributing to
  "Eq. 50" a claim it does not make (Eq. 50 is the closure `α_{1/2} = α_{M+1/2} = 0`).
- **L33 doc falsification**: `curvilinear_one_group.rst`
  §`sn-direct-seed-circle-vs-interval` asserts *"a periodic redistribution axis
  gives edge-inclusion for free; an interval axis makes you pay."* Its *reason* is
  wrong either way (the α recursion is a cumulative integral on the **quotient
  interval**). The downstream Gauss-Lobatto study inherits the correction.
- **`rules_product.py` docstring calls the degree tag "conservative"** — `[M]` it
  is **tight** (T2). One-word fix.
- **MoC `_reflected_azi_index`** (`moc/geometry.py:222`) is an `argmin` search
  computing what is exactly `n−1−k`; `[M]` verified for `n_azi ∈ {2,3,4,5,6,8,16,32,64}`,
  and both of its wrap guards are dead (`π − φ_k ∈ (0,π)` always). ERR-042-class
  latent defect on a grid where `π−φ` is not a node.

---

## 7. What is on hold, and why it is downstream

- **#325** (symmetry-exact node generation) — the primitive is BUILT and **landed**
  (`1f9d4818`), but is **not exported and has no consumer**. Repointing its three
  consumers is Q5's business, because the repoint touches the RANGE/RULE
  factorization.
- **#326** (the double-cover) — **do not fix it directly.** Still true, and now
  measured rather than argued: every re-ordering either changes nothing
  (`lexsort`/`stable` are bit-identical to production, with exact AND inexact
  nodes) or destroys the march (the fiber ordering NaNs). See **T22**.
  ⛔ **What has changed is the acceptance test.** This section used to say the
  three `xfail(strict=True)` rows XPASS "when the machinery is right". `[M]`
  **They cannot adjudicate the fold at all** — on the quotient the ξ-mirror
  partner is not in the node set, and `reflection_index("y")` returns a silent
  many-to-one garbage map rather than raising. *(⛔ That mechanism sentence
  was true when written and is doubly superseded — annotated 2026-08-07:
  Q5.0.1 made an uncertified axis RAISE, then §7d.3 deleted the accessor;
  the ξ-mirror question is now `ordinate_permutation(σ_y motion)`, which
  honestly answers `None` on the quotient. The CONCLUSION stands.)* See
  **T23** for the measurement
  and for the non-vacuous replacement (the #229 azimuthal floor).
  ⟹ #326 is now **subsumed by Q5** (the fold), not a separate item: it dissolves
  when a level becomes a single arc. See §5's Q5 block for the design of record.
- **The boundary campaign** (B3.3, B3.5, B4–B7) — resumes after Q7.

### ⚠ Two claims that must NOT be inherited

- **"#325 gates #326"** — FALSE, and I relayed it before measuring. `[M]` The
  mirror pair's η differ by ~1 ULP **today** (`1.110e-16`), not bit-equally, but
  `eta_edge = 0.5*(η_m + η_{m+1})` collapses onto the node anyway, so
  `tau_raw` is **already** `[0,1,0,1,…]` on `main`. The defect does not wait for
  exact nodes. The two issues are independent.
- **"#326 is a 0.6–7.2 % error"** — that was the angular flux compared
  ELEMENTWISE, conflating a mirror **relabeling** with real error. `[M]` Against
  an ordinate-invariant functional (scalar flux) it is `1.4e-4 – 6.6e-4`, and at
  `n_phi = 8` it is a **pure relabeling** (`1.2e-16`; un-permuting collapses the
  angular discrepancy to `3.5e-16`). The corrected magnitude is on the issue.
