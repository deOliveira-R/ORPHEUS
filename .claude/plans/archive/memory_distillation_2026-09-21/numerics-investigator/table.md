# Audit trail — numerics-investigator memory distillation (2026-09-21)

Verdicts: **KEPT** · **CUT** (lossless) · **DISTILLED** into a meta-lesson · **DUPLICATE** of a
named rule/skill clause (adjudicated by READING the file, clause named) · **ARCHAEOLOGY** →
`archive_additions.md` · **STALE** (landing named).

Every DUPLICATE verdict below was made by opening `.claude/rules/<name>.md` or
`.claude/skills/<name>/SKILL.md` in this session and reading the clause, never from memory.

## A. `lessons.md` (1 228 lines → 715)

| lines | entry (first words) | verdict | where its correction lives now |
|---|---|---|---|
| 1–8 | header + "## L1 Run the diagnostic cascade in order" | DUPLICATE + KEPT (2 lines) | `AGENT.md` § Diagnostic Cascade ("Execute in order … Do NOT skip steps"); `probe-cascade` § Anti-patterns, "Don't skip probes". The six-wrong-hypotheses cost is kept as the founding teeth. |
| 10–17 | "## L2 Curvilinear redistribution is the prime suspect" | DUPLICATE (→ 4-line stub) | `numerical-bug-signatures` Signature 1 (symptom / mechanism / diagnostic probe / catching test / why it hides), and `AGENT.md` Step 5. |
| 18–38 | "## L3 Rank-N closure … ≥2 quadrature schemes" | KEPT (compressed 21 → 12) + STATUS verified | Digest L3. `#123` verified OPEN (`gh issue view 123`, 2026-09-21) and its title literally is the protocol. Gate `tests/cp/test_peierls_rank_n_protocol.py` verified present, `assert_rank_n_structural_win` at line 90. |
| 29–38 | the Direction-C / Direction-Q / Direction-N narrative, "#121 closed, #122 closed" | ARCHAEOLOGY (pointer only) | Already cold: `_archive/direction_c_pca_rich_adaptive.md`, `_archive/direction_q_lambert_marshak_derivation.md`; plus the live topic file `frame_5_qmc_quadrature.md`. No new archive file needed. |
| 39–71 | "## L4 Convergence-rate fingerprints" | KEPT (table intact, prose 33 → 19) | Digest L4. `numerical-bug-signatures` Signature 6 carries rows 1–2 and the fix family; rows 3–5 exist nowhere else, so the table stays whole. The behavioural payload (a wrong ATTRIBUTION picks the wrong fix) is kept. |
| 72–107 | "## L5 Read the paper's stated approximation level" | KEPT (36 → 19) | Digest L5; discriminator table intact. Not in any rule or skill (checked `vv-principles` #6 and `numerical-bug-signatures` Signature 7 — both are about instrument/quadrature, not a paper's stated order). |
| 108–130 | "## L6 A curvilinear matvec … NON-FLAT" | DUPLICATE (→ 6-line stub) | `AGENT.md` Step 5 carries the rule verbatim (the entry already said so); `vv-principles` Mode 7 carries the enumerate-the-exact-regime check, citing main-lessons L27. Forensics (two ERR-049-family bugs, the Pattern-2 cure) kept. |
| 131–161 | "## L7 A direction-dependent per-ordinate moment" | KEPT point 1+3, DUPLICATE point 2 (31 → 13) | Point 2 ("matvec self-consistency is never sufficient") is `vv-principles` § Bit-identity, the paragraph ending "**tell:** a matvec-sweep round-trip at 1e-16 offered as the floor (ERR-061)". |
| 162–198 | "## L8 The project's own theory page can be the contaminated reference" | KEPT headline + 1 + 3, DUPLICATE point 2 (37 → 15) | Point 2 (prove a no-re-baseline claim by a DIRECT old-vs-new value comparison) is `vv-principles` anti-pattern #12. The headline is a genuine sharpening of #6, whose tells are all solvers — so it is an UPLIFT candidate, see `uplift.md` U1. |
| 199–216 | "## L9 A reported hang … FIXTURE cost" | KEPT (18 → 8) | Digest L9. `numerical-bug-signatures` Signature 8's catalog entry CITES this lesson ("bound the solver cost before declaring a hang") but carries no procedure, so the how-to-apply stays here. |
| 217–246 | "## L10 Error grows with refinement + a discarded info-flag" | DUPLICATE of the headline (30 → 6) | `numerical-bug-signatures` Signature 8 carries symptom, mechanism, all three discriminators and the probe. Two sharpenings kept: derive `restart` from the domain size (`coding-elegance` Pattern 7); an L1 anchor can pass by numerical coincidence at one mesh. |
| 247–265 | "L10b — recorded-but-unread … TOLERANCE-INSENSITIVITY" | KEPT (19 → 14) | Digest L10b. Not in Signature 8 (checked): the vacuous-tolerance-sweep rule, the `max_iter−1` plateau, the error-MAP shape discriminator and the `sessionfinish` blast-radius plugin. UPLIFT candidate U2. |
| 266–281 | "L10c — an all-reflective absorber is SI-HARD" | KEPT budget law + inversion + trace twin, DISTILLED the sawtooth mechanism (16 → 9) | The sawtooth MECHANISM now lives once, in L19·2 / L23·6 / L24·1 (it was derived there). The budget law `Σ_t·n_inner` and the dimension-inversion warning have no other home. |
| 282–303 | "## L11 ρ-honest stopping … residual not increment" | DUPLICATE (→ 9-line stub) + STALE | `AGENT.md` Step 3 (verbatim standing rule) and `numerical-bug-signatures` Signature 9. STALE: `FluxDisplacement.contraction_ratio` no longer exists — `grep -rn FluxDisplacement orpheus/` returns nothing. Landing named: `IterationRecord.increment_norms` / `.contraction_ratios` / `.true_error_estimate()` in `orpheus/numerics/convergence.py:1012/1076/1104`, `AngularResidual` in `orpheus/transport/residuals/angular_residual.py:89`; the retirement is recorded in `coding-elegance` anti-pattern #18. |
| 304–330 | "## L12 An OFFLINE-isolated error is only THE floor…" | DUPLICATE points 1–2, KEPT point 3 (27 → 10) | Points 1–2 are `vv-principles` § Bit-identity, "NEVER call an OFFLINE-isolated error 'the floor' … (1) end-to-end swap; (2) term-silent control; (3) AMPLIFY" — verbatim, including the ERR-061 tell. Point 3 (measure the ORDER at the cell, not in the volume-weighted norm) is in no rule: UPLIFT candidate U3. |
| 331–353 | "## L13 greedy `(Ellipsis, *idx)`" | KEPT (23 → 13) | Digest L13. A code-level signature in no skill; proposed as a new `numerical-bug-signatures` Signature (U4). Commit hash kept, dead diag path CUT. |
| 354–393 | "## L14 curvilinear `(L+C).solve` seed-lag" | ARCHAEOLOGY (verbatim → archive) + KEPT imperatives (40 → 11) | Cold at `_archive/curvilinear_seed_metric_and_ordering.md` (lossless, `diff`-verified). The kept half is the behavioural correction: my own headline was over-generalised from ONE probe, and the mechanism I published was a mis-attribution. |
| 394–439 | "## L15 PRINCIPLED-vs-REGRESSION … N not h; Krylov restart" | ARCHAEOLOGY (verbatim → archive) + KEPT imperatives (46 → 13) | Cold at `_archive/krylov_composite_restart_and_stale_references.md`. Worktree-recipe commit archaeology (`5170f20`) CUT; `[LANDED a29ab2d]` kept. |
| 440–468 | "## L16 compare two ANGULAR QUADRATURES" | ARCHAEOLOGY (verbatim → archive) + KEPT the 5-step recipe (29 → 14) | Cold at `_archive/curvilinear_seed_metric_and_ordering.md`. Dead artefact paths CUT: `scratch/experimental/glob_sphere_study/` and `derivations/diagnostics/diag_glob_0{1..5}_*.py` do not exist (`ls`, 2026-09-21). |
| 469–519 | "## L17 A STATE DOF's Hilbert metric" | ARCHAEOLOGY (verbatim → archive) + KEPT 3 points (51 → 16) | Cold at `_archive/curvilinear_seed_metric_and_ordering.md`. Dead path `derivations/diagnostics/diag_gsd_0{1,2,3}_*.py` CUT (does not exist). |
| 520–579 | "## L18 LABELING/ORDERING degeneracy" | ARCHAEOLOGY (verbatim → archive) + KEPT 5 points + 2 corrections (60 → 18); STATUS corrected | Cold at `_archive/curvilinear_seed_metric_and_ordering.md`. STALE→corrected: #326 is CLOSED and the remediation LANDED (`dde93b64`, ancestor of HEAD). The method warning ("verify the knob is REACHABLE") DISTILLED into spine M5; the `xfail`-swallows-setup-errors point is kept as the sharpening of `vv-principles` Mode 8(4), and the re-wire-at-promotion point cites Mode 11 (both read this session). |
| 582–635 | "## L19 RATE question is a spectrum question" | ARCHAEOLOGY (verbatim → archive) + KEPT 5 points + verdict (54 → 31) | Cold at `_archive/rate_spectrum_certificate_and_angular_axis.md`. STATUS corrected: #341 CLOSED, docstring repairs LANDED `adc887d6`, and the "octant-order issue owed" was FILED — it is #343, OPEN (`gh issue view 343`). |
| 638–695 | "## L20 A RESIDUAL cannot gate an EIGENVALUE contract" | ARCHAEOLOGY (verbatim → archive) + KEPT 5 points (58 → 27) | Cold at `_archive/rate_spectrum_certificate_and_angular_axis.md`. Point 5 (gate the fixture on the mixture's consistency identity) is UPLIFT candidate U5; point 1 (TRANSFER GAIN) is U6. |
| 696–748 | "## L21 ANGULAR-consistency … `h→0`" | ARCHAEOLOGY (verbatim → archive); point 2 DUPLICATE; rest KEPT (53 → 29) | Cold at `_archive/rate_spectrum_certificate_and_angular_axis.md`. Point 2's regime check LANDED as `vv-principles` anti-pattern #24(e) (commit `916a93ef`, "the four durable lessons from #235 Phase 0"; the clause is present in the skill I hold) — only the fixture-LIVENESS column survives here. STATUS: the Phase-0 gates LANDED `a3121cfe`; #319 and #235 remain OPEN. |
| 749–800 | "## L22 A frozen reference is stale by MAGNITUDE CLASS" | ARCHAEOLOGY (verbatim → archive); points 2+5 DUPLICATE; 1,3,4,6,7 KEPT (52 → 29) | Cold at `_archive/krylov_composite_restart_and_stale_references.md`. Point 2 (the sibling-pass discriminator) is `numerical-bug-signatures` Signature 10's first discriminator; point 5 (a re-baseline's radius is the frozen REFERENCES by kind, not the `.npz` files) LANDED as `vv-principles` anti-pattern #25. Points 1, 4, 6 and 7 are in no rule: UPLIFT candidates U7, U8, U9. |
| 803–945 | "## L23 SINGULARITY is a two-question object" + both addenda | ARCHAEOLOGY (verbatim → archive); points 7, 8-first-half, 10-partial DUPLICATE; rest KEPT (143 → 45) | Cold at `_archive/issue_344_singularity_kernel_and_gauge.md`. Point 7 + point 10 LANDED as `vv-principles` Mode 9's "false-RED premise" (singular `A` ⟹ the fixed point is a MANIFOLD, with the three discriminators (a)/(b)/(c)) — read in the skill this session. Point 8's parity split is `vv-principles` #13's third disguise (break the congruence class); the `‖Ad‖/‖d‖`-beside-the-error method point is not there and stays. STATUS corrected: #344 CLOSED; the kernel + gauge LANDED `f934ff57`/`b51bc802` and the 24 gates LANDED `1a2be025` — they are NOT "awaiting promotion". |
| 948–1012 | "## L24 KERNEL is a CLOSED-FORM problem" | ARCHAEOLOGY (verbatim → archive) + KEPT 6 points (65 → 31) | Cold at `_archive/issue_344_singularity_kernel_and_gauge.md`. Nothing in any rule or skill covers the counting-law/sign-character route: UPLIFT candidate U10 (a new `numerical-bug-signatures` entry was considered and rejected — it is a DERIVATION method, not a bug signature). |
| 1015–1082 | "## L25 ARITY question is a THEOREM question" | ARCHAEOLOGY (verbatim → archive); point 2 DUPLICATE; rest KEPT (68 → 29) | Cold at `_archive/curvilinear_tau_ld_and_gram_ownership.md`. Point 2 ("a SIGNATURE is not an invariance proof") is `vv-principles` Mode 8(3), the signature-tautological class — the entry already cited it; it is now one clause inside spine M7. |
| 1083–1165 | "## L26 OWNERSHIP question" | ARCHAEOLOGY (verbatim → archive) + KEPT 9 points (83 → 34) | Cold at `_archive/curvilinear_tau_ld_and_gram_ownership.md`. Point 4 (reciprocity cannot adjudicate `G`) DISTILLED with L17·3, which is its instance: L17 keeps the instance and points at L26·4 for the theorem. |
| 1168–1228 | "## L27 A NOISE mode reading is a HYPOTHESIS" | ARCHAEOLOGY (verbatim → archive) + KEPT 6 points (61 → 25) | Cold at `_archive/rcond_threshold_rederivation.md`. Point 2 cites `vv-principles` #24(d) (zero-set) itself; the two-edges sharpening is not in #24(d) and stays. UPLIFT candidate U11. |
| — | the SPINE (new) | NEW, 62 lines | Eight meta-lessons M1–M8, each naming its instances by `L` number. This is the standard's "find the meta-lesson that drives multiple failure points and encode THAT". |

**Renumbering refused, with the reason.** `numerical-bug-signatures/SKILL.md:434,493` cite
"numerics-investigator lessons L10/L9" and "lesson L11" by NUMBER, and
`derivations/diagnostics/diag_282_sphere_repose_convergence.py:23` cites "L14/L15". Renumbering
would dangle three external citations, so retired entries keep their number as a stub.

## B. `MEMORY.md` (162 lines → 86)

| lines | entry | verdict | where its correction lives now |
|---|---|---|---|
| 12–62 | §1, the 52-line keyword lookup table for L1–L27 | CUT (→ 8 lines) | It is a second copy of the digest's headings (index discipline 5: a hook is ≤ 15 words). The new §1 names the eight SPINE meta-lessons and stops; the digest is the authority, as §1's own banner already said. |
| 66–70 | "**None.** Every SN / curvilinear / Peierls campaign … merged to main" | KEPT + re-verified | Re-verified 2026-09-21: `a1c90aac`, `1a2be025`, `dde93b64`, `adc887d6`, `a3121cfe`, `916a93ef`, `317944ec` are all ancestors of HEAD. |
| 72–74 | the merge-status-goes-stale banner | KEPT (3 → 2) | It is the reason four of the lines below were wrong; kept as the discipline. |
| 76–78 | "#326 (cylindrical ordering — adjudicated, remediation open)" | STALE → collapsed | #326 is **CLOSED** and the remediation LANDED (`dde93b64` "a level's order is the FIBER's", plus the 6.3 leg series). Record: the digest L18 + `_archive/curvilinear_seed_metric_and_ordering.md`. |
| 77–79 | "**#341** … docstring fixes + the octant-order issue owed" | STALE → collapsed, successor named | #341 **CLOSED**; the docstring repairs LANDED `adc887d6` ("three published claims that were false, and the CI gate that never ran"); the octant-order issue WAS filed — **#343, OPEN**. |
| 79–83 | "**#344** (ANSWERED …; 24 gates awaiting promotion)" | STALE → collapsed | #344 **CLOSED**; the gates LANDED `1a2be025` ("the #344 characterization runs in CI") and the fix LANDED `f934ff57`/`b51bc802`. "Awaiting promotion" was false. |
| 84–85 | "**#319/#235** … 3 gates awaiting promotion" | STALE → collapsed, remainder kept | The gates LANDED `a3121cfe` ("promote the Phase-0 angular diffusion-limit gates; retire the diagnostic") and Phase 0 CLOSED `317944ec`. #319 and #235 themselves are still OPEN for later phases, so they stay in the open list with the corrected reason. |
| 85 | "#123, #128, #132/#100, #129, #170" | KEPT, all verified OPEN | `gh issue view` 2026-09-21: all five OPEN. **#200** added (OPEN, and L14/L15 both point at it). |
| 89–162 | §3, 20 "durable reference" lines | KEPT, rewritten to ≤ 15-word hooks; six `_archive/` lines added | The six new archive files are the durable half of five records that are UNTRACKED in `scratch/`. |
| 92, 104, 111, 117 | "⭐ `scratch/…md` (in-repo, not agent memory)" ×4 + `scratch/rcond_rederivation.md` | KEPT with a correction | `git ls-files --error-unmatch` 2026-09-21: `issue_344_null_space_structure.md`, `adjoint_gram_ownership_audit.md`, `tau_under_ld_dip_analysis.md`, `issue_344_kernel_basis.md` and `rcond_rederivation.md` are **UNTRACKED** — "in-repo" over-claims. The index now says so, and the archive move makes their lessons durable. (`issue_341_gs_jacobi_mechanism.md`, `n5_outer_certificate_measurement.md` and `q68_flux_dip_discriminator.md` ARE tracked.) |

## C. Topic files judged ARCHAEOLOGY (listed only — nothing retired in this pass)

Inbound-reference counts measured over `docs/`, `.claude/plans/` and `.claude/agent-memory/`
(973 files, excluding `_build` and the file itself), 2026-09-21.

| file | inbound refs | reason it reads as archaeology |
|---|---|---|
| `phase5a_moment_consuming_scatter_derisk.md` | **0** | Phase-codename title ("phase5a"), no referrer anywhere in docs, plans or memory. |
| `r1_step_d_sphere_preconditioner_oscillation.md` | **0** | Step-codename title ("R1 step D"); the live successor is #200 and the seed taxonomy. |
| `sn_sig_t_layout_drift_indexerror.md` | **0** | A single fixed layout-drift IndexError; ERR-055 catalogues it and the catalogue is the durable home. |
| `phase_f_step2_mesh_refinement.md`, `phase_f_step3_diagnostic.md`, `phase_d_gate_1_1_sphere_mms_diagnosis.md` | 9 / 8 / 12 | NOT archaeology despite the phase codenames — they are cited from `docs/theory/`, so they are load-bearing graph nodes. Listed here only so the orchestrator does not mistake the codename for disuse. |

The orchestrator runs the full `retirement-audit` blast radius before any of the three
zero-reference files is touched: they are Nexus graph nodes, and `-W` warns on none of it.

## D. Counts (X2)

**Lines, per file:**

| file | before | after | note |
|---|---|---|---|
| `lessons.md` | 1 228 | **715** | −42 %. 865 of the removed lines are not deleted: they moved verbatim into six `_archive/` files (below), `diff`-verified lossless. |
| `MEMORY.md` | 162 | **86** | −47 %. §1's 52-line keyword table collapsed to 8 lines; §2's stale status block to 10. |
| `_archive/` (new files) | 0 | **+904** | Cold, never loaded: six files, each a verbatim `sed` extraction with a 5-line provenance header. |

**Hot-surface accounting:** 1 390 loaded lines (`lessons.md` + `MEMORY.md`) → **801**. The figure
the standard names is ≤ 400 for the digest; this proposal lands at 715 and says so rather than
cutting corrections. The remaining 315 lines are the nine dense entries L19–L27 (their numbered
point lists). **A further cut to ≈ 400 is available and is the orchestrator's call**: demote
those nine numbered lists to their archive files entirely, leaving a 6-line entry each. The cost
is that the discriminators (which measurement distinguishes which mechanism) become one `Read`
away instead of in-context, which is exactly the trade the standard's hot/cold split makes
everywhere else — I did not take it unilaterally because these nine are the newest and least
internalised material.

**Rows per verdict** (28 digest rows in table A, counting the new SPINE row):

| verdict | rows | which |
|---|---|---|
| KEPT (whole or compressed) | 12 | L3, L4, L5, L9, L13, L14, L15, L16, L17, L19, L24, L26 |
| DUPLICATE (whole entry → stub) | 4 | L1, L2, L6, L11 |
| DUPLICATE (part of the entry) | 8 | L7·2, L8·2, L10 headline, L12·1–2, L18 (two clauses), L21·2, L22·2+5, L23·7+8+10, L25·2 |
| DISTILLED into a meta-lesson | 3 | L10c's sawtooth → L19·2/L23/L24; L18's knob-reachability → M5; L26·4 ↔ L17·3 |
| ARCHAEOLOGY → `_archive/` | 9 | L14, L15, L16, L17, L18, L19+L20+L21, L22, L23+L24, L25+L26, L27 (six files) |
| CUT (dead path / commit archaeology, nothing lost) | 5 | the `~/.claude/jobs/…` probe paths; `diag_glob_*`; `diag_gsd_*`; `diag_344_*`; `diag_276_*`; `scratch/experimental/glob_sphere_study/` — each verified absent by `ls` |
| STALE (landing named) | 5 | L11's `FluxDisplacement`; and `MEMORY.md`'s #326, #341, #344, #319/#235 status lines |

**Rules and skills READ this session to adjudicate a DUPLICATE** (7): `vv-principles`,
`numerical-bug-signatures`, `probe-cascade`, `instrument-doctrine`, `coding-elegance` (skills,
held in-context); `.claude/rules/vv-testing.md`, `.claude/rules/coding-standards.md`,
`.claude/rules/code-search.md`, `.claude/rules/process-discipline.md` (read from disk).

## Orchestrator's notes at apply time (2026-09-21, the main agent with the session's context)

Applied as proposed; the six `_archive/` files copied after re-running the nine `diff`s (9 of 9 empty). **Ruling on 715 vs 400: the digest stays at 715.** L19–L27's numbered points are each a corrected mistake with its measurement, the newest and least-internalised material; the count is a reading and the per-line test is the criterion (the other passes landed at 1 013, 721, 405 and 394 by the same test). Applied from `agent_md.md`: P1 ("Before the cascade: what KIND of question is this?") ahead of the cascade; P2 not (X1 covers it; spine M5 keeps the instances). Applied from `uplift.md`: U1 (`vv-principles` #6 check and tell), U2 (`numerical-bug-signatures` Signature 8's caveat), U3 (`vv-principles` #5 check and tell), U5 (`vv-testing`, the mixture identity), U6 (the `instrument-doctrine` skill's X1 transfer gain), U7 and U8 (`vv-principles` bit-identity, one paragraph), U9 (`code-search`, two probe traps), U11 (`vv-principles` #24 clause (f)). HELD in the uplift queue: U4 (the new signature owes an ERR entry first, per the skill's own protocol) and U10. The digest's entries whose corrections moved gained pointers. The five untracked `scratch/` records are not tracked here (the archive copies carry the lessons; a scratch file is not a memory surface); the three zero-reference topic files are not retired (graph nodes; the blast-radius audit is a separate pass).

Counts as applied (`wc -l`): lessons.md 1228 → 721; MEMORY.md 162 → 86; `_archive/` 5412 → 6316 lines.
