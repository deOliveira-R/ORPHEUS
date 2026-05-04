# Plan — Peierls documentation cleanup (post-#138 close-out)

**Date filed:** 2026-04-30
**Predecessor:** #138 full-collapse (commits 8438df9, 1940b7f, 99a05ab, 742d3b0 — landed locally on `refactor/peierls-138-full-collapse`, **not yet pushed to origin/main**)
**Trigger:** User directive after #138 audit — "documentation should reflect code in detail, but should not necessarily reflect a bunch of failed experiments. The right place to document failed experiments would be in the GitHub issues that started them."
**Branch (new):** `docs/peierls-cleanup` (cut from `main` AFTER pushing #138 cascade)
**Estimated sessions:** 3 sessions (issue post-mortems → doc edits → verification)

---

## Pre-conditions (verify before starting)

1. **#138 cascade is on `origin/main`.** This plan assumes the 134-commit cascade landed via `git push` from the prior session. If not yet pushed, push first then start this branch.
2. **Sphinx -W clean.** Run `rm -rf docs/_build && sphinx-build -W -b html docs docs/_build/html` — must be EXIT 0 with 0 warnings (#139 fix in 8438df9 made this true).
3. **Issue #138 status.** Verify it's CLOSED with the full-collapse summary comment posted.

---

## What this plan does

Cleans up the two largest Peierls docs by relocating ~3500 LoC of failed-experiment narrative to GitHub issues (their natural home), keeping the docs as a "Sphinx-as-brain" reference for current code. Net: `peierls_nystrom.rst` shrinks 11479 → ~8344 LoC (−27%); `collision_probability.rst` shrinks 4342 → ~3880 LoC (−11%). Total: ~3587 LoC moved to issues.

## What this plan does NOT do

- Does not touch the slab native E₁ Nyström module (`peierls_slab.py`) — preserved as V-of-V per `§theory-peierls-slab-polar-retirement`.
- Does not modify any code outside docs.
- Does not relitigate Issue #138 (the API-posture rewrite in commit 8438df9 stays).

---

## Source-of-truth audits

Four Archivist agents produced detailed line-by-line audits in the prior session. Their classifications are the authoritative input for this plan. Re-read them before starting:

| Audit partition | Agent ID | Focus |
|---|---|---|
| `peierls_nystrom.rst` lines 1–3697 | a23cb57a2c3218c15 | Architecture, capabilities, glossary, slab-polar, API posture (preserve), MG, Moment-form ARCHIVED, Sections 1-7 |
| `peierls_nystrom.rst` lines 3698–7212 | a13809078c2f05e93 | §8 White-BC, F.5 #119 close-out, F.6 #132 falsification, Class B Hébert (1-P_ss)⁻¹, Specular BC + multi-bounce |
| `peierls_nystrom.rst` lines 7213–11479 | ad5f66a205d114f0f | Sections 9-17 (test-bed, hierarchy, FlatSourceCPGeometry, BickleyTables retirement), §22 coordinate transforms, §23 MC, Apps D+E, Part IV tensor structure, §29 deferred Marshak, §30 BoundaryClosureOperator |
| `collision_probability.rst` Peierls sections | a18306500e7781213 | Slab/cyl/sph Peierls reference subsections; Issue #100 retraction duplicated content |

If those agents have been garbage-collected (likely across session boundary), use this plan as the source of truth — line ranges and destination issues are reproduced below.

---

## Phase 0 — Cut branch + lock baseline

1. From `main` (post #138 cascade push):
   ```
   git checkout main && git pull --ff-only
   git checkout -b docs/peierls-cleanup
   ```
2. Verify baseline: `rm -rf docs/_build && sphinx-build -W -b html docs docs/_build/html` → EXIT 0, 0 warnings.
3. Snapshot `wc -l docs/theory/peierls_nystrom.rst docs/theory/collision_probability.rst` for delta tracking.
4. `gh issue list --state all --label module:peierls --json number,title,state` to confirm issue states.

---

## Phase 1 — Issue close-out comments (the load-bearing prerequisite)

**Cardinal Rule 1 (Sphinx is the brain) only holds if relocated knowledge lives somewhere durable.** This phase posts 12 close-out / investigation-log comments to GitHub issues BEFORE the docs are edited. Each comment carries the relocated narrative + math + tables + citations.

### Closed issues — need post-mortem comments

| Issue | Topic | Source lines | Content to relocate |
|---|---|---|---|
| **#94** | BickleyTables retirement | `peierls_nystrom.rst:8314–8413` | Replacement table (5 rows: legacy → canonical) + retirement sequence (Phase B.1/B.2/B.4 commit chain + 4-bullet measured impact: k_inf shifts, agreement, table-size test retirement, latency) |
| **#100** | Sphere white-BC retraction | `peierls_nystrom.rst:3824–3915` + `collision_probability.rst:3922–4051` | The retraction debate (pre-correction "ratio varies 40 % so rank-1 fails" argument + post-correction "that conflated u/v with J⁻"); preserve numerical tables at CP:3946–4011 |
| **#113** | Cylinder failed-scheme forensic | `collision_probability.rst:2324–2347` | Failed-scheme narrative (full forensic detail) |
| **#117** | Moment-form Nyström archived | `peierls_nystrom.rst:1355–2289` | ENTIRE 935-line section: slab specialisation, Vandermonde solve, moment recursion gates, GL convergence tables, performance numbers, **plus** the τ-Laguerre polar-form failure narrative within (lines 1745–1942 — falsified-experiment-of-falsified-experiment with the convergence table + α-scan analysis) |
| **#119** | F.5 rank-N per-face close-out | `peierls_nystrom.rst:4443–4538` | Villarino-Stamm'ler per-mode (falsified 2026-04-21) — ~95 lines of failed-extension narrative |
| **#122** | Lambert/Marshak gauge investigation | `peierls_nystrom.rst:4580–4793` | Lambert/Marshak gauge close-out — ~215 lines |
| **#131** | Slab MR E₂ closed-form | `peierls_nystrom.rst:803–859` | Probe A/B/D diagnostic cascade (1.5 % gap diagnosis) — ~57 lines |
| **#133** | Visibility-cone rollout (closed wontfix variant) | `peierls_nystrom.rst:7017–7196` | Phase 5 continuous-µ retreat — 6 rounds of investigation across Phase 5a/Round 1 (3 fronts)/Round 2 (PRIMARY+BACKUP)/Round 3 (PRIMARY+SECONDARY) — ~180 lines |
| **#133/#134/#135/#136** | Quadrature rollout completed | `peierls_nystrom.rst:9939–10169` | §22.9 per-primitive landing map (17-row table) + originating commits (Q1 → L3 in dependency order) + acceptance-criterion audit — ~245 lines |

### Open issues — need investigation-log comments (active research record)

| Issue | Topic | Source lines | Content to relocate |
|---|---|---|---|
| **#112** | Phase A/C canonical Marshak deferred | `peierls_nystrom.rst:11117–11205` | §29 empirical variant scan (sphere + cylinder k_eff error tables for V1/V2/V4/V6) + observations 1–3 — ~88 LoC |
| **#123** | Rank-N stability protocol L19 | `peierls_nystrom.rst:4794–4889` | Full L19 stability protocol — ~95 lines |
| **#132** | Class B MR×MG rank-N falsification | `peierls_nystrom.rst:5267–5409` (Probe G LEGACY/CANONICAL — ~145 LoC) + `peierls_nystrom.rst:5807–5970` (Davison + Aug Nyström + Cyl Hébert follow-up dirs — ~165 LoC) | Probe-cascade narrative (relocate; keep core empirical evidence at 5109–5217 in docs); follow-up direction investigation full forensic detail |

### Need to file (NEW issues)

| Topic | Source lines | Why a new issue is needed |
|---|---|---|
| **Planar-limit cross-check** | `peierls_nystrom.rst:885–907` | Currently sits as an OQ bullet inside §slab-polar with no issue number. Empirical probe at r₀=0.999R found 22 % disagreement; chord-distribution explanation is non-trivial. Per directive, this is failed-experiment-style content that needs an issue home before relocation. |
| **MG residual benchmark vs discrete CP** | `peierls_nystrom.rst:1313–1326` | "Multi-group residual benchmarks vs the discrete CP modules… have not been benchmarked as part of Issue #104 commit 2 — the parity gate is the planned follow-up". Needs its own tracking issue. |

### Verification artifact

| Item | Action |
|---|---|
| `tests/l0_error_catalog.md` ERR-NNN for `∫E₂` algebra bug | Verify catalog entry exists for the `peierls_nystrom.rst:3601–3624` "History of the algebra bug" forensic. If not, file the catalog entry per Cardinal Rule 3 BEFORE trimming the doc text. |

### Phase 1 commit cadence

One commit per issue close-out. Commits are docs-only (just the new issue comments aren't tracked in git, but the local plan + reference notes can be). If you want a git audit trail of "what was relocated", create a `docs/peierls-cleanup-log.md` scratch file inside this plan directory listing each relocation and its destination, then delete the scratch file in the final Phase 4 cleanup commit.

---

## Phase 2 — Doc edits, top-down

Edit in this order so cross-references resolve cleanly (later sections may reference earlier ones).

### Phase 2a — `peierls_nystrom.rst` lines 1–3697 (Audit a23cb57a)

Order matters: do the bigger relocations first, then the staleness sweeps.

1. **§Moment-form ARCHIVED → 25-LoC stub** (lines 1355–2289, −910 LoC).
   - Preserve anchors `theory-peierls-moment-form` and `theory-peierls-moment-form-failed-polar` (cross-referenced from line 871 and Key Facts).
   - Stub points to: (a) `derivations/archive/peierls_moments.py` + `peierls_slab_moments_assembly.py`, (b) Issue #117 closing comment, (c) git commit `investigate/peierls-solver-bugs`.
   - One sentence on τ-Laguerre failure ("plateaued at 9e-4; see #117 for the diagnostic table + α-scan").

2. **§G.5 routing → relocate Issue #131 probe cascade** (lines 803–859, −52 LoC).
   - Replace the "How the 1.5 % gap was diagnosed (Issue #131)" subsubsection with: "Issue #131 closed the original 1.5 % multi-region gap — see [issue link] for the diagnostic cascade (Probes A/B/D, the underconvergent GL branch, the closed-form fix). Both paths now agree bit-exactly on the shipped fixture."
   - Remove the dead cross-reference to §G.5-diagnosis in the MG section (lines 1144–1158).

3. **§Slab-polar related-OQs → trim planar-limit narrative** (lines 885–907, −20 LoC). Replace with: "**OQ — Planar-limit cross-check against hollow cylinder?** Probed empirically and found structurally non-trivial; tracked in Issue #NNN. The cylinder's Ki₁ has already integrated axial direction, giving a different chord-length distribution than slab — there is no simple thin-shell equivalence at matched optical parameters."

4. **§White-BC analytical slab → trim algebra-bug forensic** (lines 3601–3624, −18 LoC). After verifying the ERR-NNN entry exists in `tests/l0_error_catalog.md`, replace the 24-line forensic with ~6 lines: state the algebra mistake (`∫E₂ = ½ - E₃, not 1 - E₃`), state the lesson ("two independent derivations agreeing at 1e-39 is worthless if both share a factor-of-2 algebra bug"), reference the catching test name (`TestSlabKernelRowSum`). Forensic narrative lives in catalog + git history.

5. **§Key Facts milestone-tag trim** (lines 162–214, −37 LoC). Three dated milestones bleed into Key Facts:
   - Trim "Historical note. Standard references…" to keep technical content (dimensional-reduction observation), drop milestone framing.
   - Trim "Production status (2026-04-20)" to one sentence: "Slab K is production via the unified adaptive `mpmath.quad` path".
   - Trim "Vacuum-BC verification milestone (2026-04-20)" — keep the *fact* (link to §peierls-vacuum-bc-analytical-references), drop the milestone-narrative wrapper.

6. **§Capabilities — staleness sweep** (lines 359–375, ~−9 LoC). Verify status of Issues #103, #101, #112; trim closed-issue bullets, reframe as "tracked in #NNN" for any still active.

7. **§MG section — staleness sweep** (lines 1313–1326, ~−10 LoC). After filing the new "MG residual benchmark" issue (Phase 1), trim to "Multi-group residual benchmarks vs the discrete CP modules are tracked separately under Issue #NNN."

8. **§Motivation and scope — staleness fix** (lines 2306–2308, ~−5 LoC). Replace the stale "(Planned) `peierls_sphere`" bullet with present-tense: "`orpheus.derivations.peierls_sphere` — sphere `e^{-τ}` reference (registry-only façade per Issue #138; canonical entry is `solve_peierls_mg(_pg.SPHERE_1D, ...)`)".

9. **§API posture — DO NOT TOUCH** (lines 912–1117). This was just rewritten in commit 8438df9 as the deliberate Issue #138 posture statement. Preserve.

10. **§Retirement candidates — minor cleanup** (lines 1051–1056, ~−3 LoC). Replace empty-section placeholder with one sentence pointing at V-of-V: "None as of 2026-04-29; the slab native E₁ Nyström is classified as verification-of-verification (see below), not transitional."

**Sub-phase 2a verification**: `sphinx-build -W` → EXIT 0. Commit: `docs(peierls): #138 close-out — relocate moment-form + Probe cascade + milestone trims (lines 1-3697)`.

### Phase 2b — `peierls_nystrom.rst` lines 3698–7212 (Audit a13809078)

This is the heaviest partition (−1455 LoC, 42% reduction). Execute in section-natural order:

1. **§8 Sphere Issue #100 retraction → 25-LoC stub** (lines 3824–3915, −70 LoC). One sentence: "Sphere white-BC rank-1 was probed (Issue #100); the closure is structurally correct, the original 'rank-1 fails because ratio varies 40 %' argument was wrong (conflated u/v with J⁻). See Issue #100 for the historical debate." Preserve any cross-reference labels.

2. **§8 Rank-N skeleton + cyl/sph magnitude limitations → trim** (lines 3918–4006, −45 LoC). Keep skeleton derivation (formulas), relocate magnitude-limitation narrative to Issue #112.

3. **§8 (ρ_max/R)² Jacobian fix narrative → trim heavily** (lines 4008–4207, −155 LoC). Keep one paragraph: the canonical Jacobian formula + which test gates it (`peierls-rank-n-jacobian-derivation` label survives). Relocate pre-fix/post-fix tables and retraction debate to Issues #112/#132.

4. **§F.5 status + motivation trim** (lines 4211–4248, −25 LoC). Compress to one paragraph: "Issue #119 closed F.5 rank-N per-face as falsified beyond F.4. The Stamm'ler IV Eq. 34 rank-2 closure is the production decision."

5. **§F.5 five-reference synthesis trim** (lines 4248–4364, −70 LoC). Keep the F.4 residual table (4-3 row block — production data); relocate literature synthesis details to Issue #119.

6. **§F.5 structural obstruction (`c_in` remap) → KEEP** (lines 4365–4442). This is the durable lesson — what physically prevents rank-N from working. Preserve `:label:` `c-in-remapping`.

7. **§F.5 Villarino-Stamm'ler per-mode (falsified) → relocate to #119** (lines 4443–4538, −95 LoC). Preserve label `peierls-change-of-basis` if cross-referenced from elsewhere.

8. **§F.5 production decision + F.4 residual scan → KEEP** (lines 4540–4578).

9. **§F.5 Lambert/Marshak gauge → relocate to #122** (lines 4580–4793, −155 LoC). Keep ~60-LoC stub: state the gauge ambiguity, point at #122 for the empirical scan + the Stepanek 1981 calibration plan.

10. **§F.5 L19 stability protocol → relocate to #123** (lines 4794–4889, −95 LoC).

11. **§F.5 infrastructure retained → KEEP** (lines 4891–4931). Documents what dead-code lives in `peierls_geometry.py` (e.g. compute_P_esc_mode_marshak) for future Phase A/C session.

12. **§F.5 open research bullets → relocate to #120/#121** (lines 4933–4970, −35 LoC). Verify these issues exist (they may need filing).

13. **§F.5 session trail → DELETE** (lines 4972–5022, −50 LoC). Pure investigation log; git history + issue comments cover it.

14. **§F.6 status + hypothesis → trim** (lines 5025–5108, −60 LoC). One paragraph: Issue #132 is OPEN, rank-N on Class B MR×MG empirically falsified, the Hébert (1-P_ss)⁻¹ partial fix is the active resolution path.

15. **§F.6 MR×MG empirical evidence → KEEP** (lines 5109–5217). This IS the active evidence for Issue #132; preserve the 3 observations.

16. **§F.6 Probe-cascade + Probe G LEGACY/CANONICAL → relocate to #132** (lines 5267–5409, −95 LoC). Keep ~50 LoC: state the root cause (Probe G normalisation mismatch), reference issue for the full cascade.

17. **§F.6 production decision + Issue #132 fix paths → KEEP** (lines 5410–5476).

18. **§F.6 infrastructure + diagnostic-script list → trim** (lines 5478–5560, −50 LoC). Keep infrastructure-retained list; remove the diagnostic-script enumeration (lives in `derivations/diagnostics/`, discoverable by `ls`).

19. **§F.6 lesson + session trail → trim** (lines 5562–5605, −33 LoC). Keep ERR-030 reference; delete session-trail listing.

20. **§Class B Hébert (1-P_ss)⁻¹ resolution → KEEP** (lines 5608–5805). All current production code; preserve `peierls-class-b-Jn-canonical` label.

21. **§Class B follow-up directions → trim** (lines 5807–5970, −120 LoC). Keep the cylinder Hébert direction (active code in `compute_P_ss_cyl`); relocate Davison + augmented-Nyström forensic to Issue #132.

22. **§Class B synthesis → trim** (lines 5972–6033, −50 LoC). Compress 3 conclusions + 3 candidate paths to 6 bullets pointing at Issue #132.

23. **§Specular BC defn + R_specular derivation + ladder → KEEP** (lines 6088–6224). Production code (`compute_T_specular_*`).

24. **§Specular no-Jacobian P primitives → KEEP** (lines 6225–6263). Active gotcha; preserve the rationale for why specular branch uses no-Jacobian primitives.

25. **§Specular convergence sphere/cyl/slab → KEEP** (lines 6265–6525). Numerical evidence for current production code.

26. **§Specular verification posture + tests → trim** (lines 6526–6605, −40 LoC). Keep the "verification posture" framing; remove the per-test inventory (rediscoverable via `pytest --collect-only`).

27. **§Specular MG/MR convergence (Phase 2) → KEEP** (lines 6606–6727).

28. **§Multi-bounce specular (Phase 4) → KEEP** (lines 6729–7016). Production code (`compute_T_specular_multibounce`); the multi-bounce derivation + best-use envelope are load-bearing.

29. **§Phase 5 continuous-µ retreat → relocate to #133** (lines 7017–7196, −150 LoC). Keep ~30 LoC: state that `compute_K_bc_continuous_mu_*` raises NotImplementedError, kernel is hypersingular, see #133.

**Sub-phase 2b verification**: `sphinx-build -W` → EXIT 0. Commit: `docs(peierls): #138 close-out — relocate F.5/F.6/Phase 5 narrative to issues (lines 3698-7212)`.

### Phase 2c — `peierls_nystrom.rst` lines 7213–11479 (Audit ad5f66a)

1. **§9 Phase 4.2 test-bed → minor trim** (lines 7213–7298, −3 LoC). Delete commit-hash parentheticals (line 7274 `(Phase-4.3; commits 435c0b3, 9d03948, cad2f0b)` — `git log` is authoritative).

2. **§14 FlatSourceCPGeometry → staleness sweep** (lines 7982–8176, −10 LoC). Code IS shipped in `cp_geometry.py`. Replace ALL "Phase B.2 / not yet shipped / will deliver / target / intended skeleton" forward-looking phrases with present-tense. Caption `# orpheus/derivations/cp_geometry.py  (Phase B.2, not yet shipped)` → `# orpheus/derivations/cp_geometry.py`.

3. **§15 escape probabilities L2 → minor staleness fix** (lines 8177–8267, −5 LoC). Convert "available for Phase B.3" to "tracked in Issue #NNN" or just remove if no follow-up tracker exists.

4. **§12 derivation-source-of-record → trim** (lines 7833–7853, −15 LoC). Keep first paragraph (cites SymPy script in `collision_probability`). Decide on second paragraph: if `derive_second_difference()` was lifted into `cp_geometry.py`, document it; if not, delete the deferred-work paragraph or shrink to a one-line "Future: lift script into `cp_geometry.derive_second_difference()`. See Issue #NNN."

5. **§16 BickleyTables retirement → aggressive trim** (lines 8269–8445, −95 LoC).
   - Keep status note + intro (8269–8313, ~45 LoC) — the 20 000-point lru_cache + scipy.quad context is load-bearing.
   - Relocate replacement table (8314–8341, ~28 LoC) to Issue #94 close-out.
   - Trim retirement sequence + Phase B.4 postmortem (8343–8413, ~70 LoC) to ~10 LoC of measured-impact bullets; full sequence narrative to Issue #94.
   - Keep "Why Chebyshev not direct mpmath" (8414–8445, ~32 LoC) — explains active design.

6. **§22.3 τ-coordinate → trim MC-redundant text** (lines 8716–8975, −20 LoC). Keep math (Issue #109 OPEN, doc-of-record for deferred Phase H). Trim ~20 LoC of MC-tracking explanation that §23.1 covers verbatim.

7. **§22.6 Davison u=r·φ → minor trim** (lines 9139–9221, −15 LoC). Compress the Case-de Hoffmann-Placzek bullet.

8. **§22.7 visibility-cone provenance + planned-rollout table → relocate** (lines 9224–9688, −75 LoC). Keep math (production primitive) + numerical-evidence table. Delete the "Planned rollout sites" Phase 1A-1F table (9626–9663) — superseded by §22.9 (post-rollout audit). Relocate "Provenance" historical detail (9665–9687) to Issue #133 close-out.

9. **§22.9 quadrature rollout audit → relocate to #133/#134/#135/#136** (lines 9910–10216, −245 LoC). Keep ~60-LoC stub: intro paragraph + "intentional residual consumers" bullet list (slab `compute_T_specular_slab`, `_compute_K_bc_sanchez`, `build_volume_kernel` non-adaptive — these are load-bearing for future Phase work) + cross-link to issues.

10. **§29 deferred Marshak → trim heavily** (lines 11058–11298, −120 LoC).
    - Keep intro + continuous derivation + structural-entanglement lesson + Stepanek path-forward + deferred-work status block (~120 LoC total).
    - Relocate empirical variant scan tables (V1/V2/V4/V6 sphere + cylinder, lines 11117–11205, ~88 LoC) to Issue #112.
    - Trim observations 1–3 (lines 11210–11225, ~24 LoC → ~6 LoC) — keep observation 4 only as one sentence.

11. **§seealso closing trim** (lines 10577–10638, −3 LoC). Delete commit-hash parentheticals.

12. **Sections 9, 10, 11, 13, 17, 22.1, 22.2, 22.4, 22.5, 22.8, 23, App D, App E, Part IV §§24-28, §30, References → KEEP entirely**. These are foundational math, active gotchas, or current code documentation — no relocation needed.

**Sub-phase 2c verification**: `sphinx-build -W` → EXIT 0. Commit: `docs(peierls): #138 close-out — BickleyTables/§22.9/§29 relocations (lines 7213-11479)`.

### Phase 2d — `collision_probability.rst` Peierls sections (Audit a18306500e)

1. **Stale cross-ref fix immediate** (line 3046, +0 LoC, −0 LoC). Replace `solve_peierls_cylinder_1g(..., boundary="vacuum")` → `solve_peierls_1g(geometry=_pg.CYLINDER_1D, ..., boundary="vacuum")`. Already partially done in #138 commit 8438df9; verify no stragglers remain.

2. **Issue #100 retraction duplicate → cut to ~25-LoC stub** (lines 3922–4051, −105 LoC). The same retraction debate exists in `peierls_nystrom.rst:3824–3915` (also being cut to a stub in Phase 2b step 1). The CP-page version should be a one-paragraph cross-reference to the unified page's stub. **Preserve the numerical tables at lines 3946–4011** — they're the production data, not narrative.

3. **Issue #113 failed-scheme forensic → relocate** (lines 2324–2347, −25 LoC). Move full forensic detail to Issue #113 close-out comment.

4. **Phase-4.2 tracked-work narrative → trim** (lines 3078–3091, ~−10 LoC; lines 3115–3123, ~−7 LoC). Compress to "tracked in Issue #NNN" pointers.

5. **V&V harness label retrofit → file as separate issue** (line 2738). The audit flagged this as needing its own tracking issue; file it before the cleanup PR merges.

6. **CaseZweifel R_c table programmatic ingestion → file as separate issue**. Audit flagged — file it.

7. **Slab leg of #138 → file as separate issue or roll into Green's-function-trigger plan**. The slab native solver is the carved-out V-of-V exception per the user's directive in this session; track its eventual retirement under the Green's function plan rather than as a #138 followup.

**Sub-phase 2d verification**: `sphinx-build -W` → EXIT 0. Commit: `docs(cp): #138 close-out — collapse Issue #100 retraction duplicate, trim #113 forensic`.

---

## Phase 3 — Final verification + audit

1. `rm -rf docs/_build && sphinx-build -W -b html docs docs/_build/html` → EXIT 0, 0 warnings.
2. Spot-check critical cross-doc citation links (`Sanchez1982`, `Stamm1983`, etc. from `collision_probability.html` resolve to `peierls_nystrom.html#stamm1983`).
3. `python -m tests._harness.audit` → ERR coverage 30/31 (baseline match).
4. Spot-check Nexus: `mcp__nexus__provenance_chain` from `peierls-cylinder-equation` / `peierls-sphere-equation` to a code node — should still land on `_pg.solve_peierls_mg`.
5. **Full peierls test suite (chunked)** — same chunked approach as #138 verification:
   - batch 1: `tests/derivations/test_peierls_{cylinder,sphere}*.py + tests/cp/test_peierls_*_flux.py` (~80 tests, ~3 min)
   - batch 2: `tests/derivations/test_peierls_{assembly_drivers,closure_operator,geometry,unified_verification,rank_n_*}.py + test_peierls_multigroup.py::TestMG{Ng1BitMatch1G,2GHollowRegistration,SlabPolarMatchesNativeSlabMG} + test_peierls_convergence.py` (~130 tests, ~10 min)
   - **No code changes — tests should pass identically to baseline.** Any failure here would indicate a doc edit broke a doctest or example block.

6. LoC delta verification:
   - `peierls_nystrom.rst`: should be 8000–8500 LoC (down from 11479)
   - `collision_probability.rst`: should be 3800–3900 LoC (down from 4342)

---

## Phase 4 — Merge + close-out

1. FF-merge `docs/peierls-cleanup` → `main`:
   ```
   git checkout main && git pull --ff-only
   git merge --ff-only docs/peierls-cleanup
   git push origin main
   ```
2. **Issue close-out comments**: ensure all 12 Phase-1 comments are posted before pushing. Each comment should link to the corresponding doc commit for a back-pointer.
3. **2 new issues filed**: planar-limit cross-check, MG residual benchmark. Plus any side issues found during audit (V&V harness label retrofit, CaseZweifel R_c programmatic ingestion).
4. **Delete** `.claude/plans/peierls-docs-cleanup.md` (this file) — done. Optionally archive to `.claude/plans/archive/` if the team values plan-history retention.
5. **Update `.claude/lessons.md`** with any meta-lesson learned (e.g. "Failed-experiment narrative belongs in GitHub issues, not in evergreen theory docs. Sphinx-as-brain ≠ Sphinx-as-history. The L21-style 'we tried X and it failed because Y' content is load-bearing for the issue, not for the next agent reading the theory page.")
6. Delete branch local + remote.

---

## Risk register

### R1 — Issue close-out comments grow unbounded

12 comments, some carrying ~250 LoC of math/tables. GitHub issue comments handle large markdown but are cumbersome to navigate. **Mitigation**: each comment leads with a 3-line summary + collapsible `<details>` block for the long content. Use the `.claude/scratch/` directory to draft each comment locally first, then paste-and-post.

### R2 — Anchor preservation across cuts

Many `:label:` and `:ref:` anchors are referenced from sections that survive. Killing surrounding prose can orphan an anchor. **Mitigation**: before each Phase 2 sub-phase, `grep -n ":label: <anchor-name>" docs/theory/*.rst` to map every anchor in the cut-target range, and confirm survivors. Specific known-load-bearing anchors:

- `theory-peierls-moment-form` (referenced from line 871) — preserve in §Moment-form stub.
- `theory-peierls-moment-form-failed-polar` (referenced from Key Facts) — preserve in §Moment-form stub.
- `c-in-remapping` (in F.5 structural obstruction) — preserve.
- `peierls-rank-n-bc-closure`, `peierls-rank-n-jacobian-derivation`, `peierls-rank-n-P-esc-moment` — survive the §8 trim.
- `hebert-3-323` — survives in Class B Hébert section.
- `peierls-change-of-basis` (Villarino-Stamm'ler) — verify if cross-referenced; if so, preserve in stub.
- `peierls-M-rank-1`, `peierls-M-rank-2` — survive in §8 white-BC.
- `peierls-class-b-Jn-canonical` — survives in Class B Hébert.
- `peierls-vacuum-bc-analytical-references` — survives at line 3390.
- `peierls-scattering-convention` — survives in MG section.

### R3 — Numerical tables disappearing as collateral

Two tables in CP page (3946–4011) and several in F.5/F.6 narratives are PRODUCTION DATA, not failed-experiment evidence. **Mitigation**: Phase 2 step descriptions explicitly call out "preserve numerical tables at X-Y" — the cuts are for narrative, not data. When in doubt, KEEP the table.

### R4 — Cardinal Rule 1 violation (Sphinx as brain)

If a future agent picks up Issue #112 and the empirical Marshak scan tables are gone from docs AND the Issue #112 comment was lost, the agent re-runs the same scan. **Mitigation**: Phase 1 (issue comments) is gating for Phase 2 (doc edits) — explicit dependency. Each Phase 2 sub-step references the corresponding Phase 1 commit hash in the relocation comment. If Phase 1 is incomplete, Phase 2 stops.

### R5 — Sphinx -W gate stays clean

The cleanup may surface broken cross-refs not visible today (e.g. cuts that orphan an `:eq:` or `:ref:` target). **Mitigation**: clean-rebuild Sphinx after every sub-phase commit. The Phase 0 baseline EXIT 0 / 0 warnings is the contract.

### R6 — Issue #138 confusion

Issue #138 was just closed (per prior session) with the API-posture flip + cyl/sph wrapper retirement. This cleanup is a SEPARATE, parallel close-out — it does NOT reopen #138. Each Phase 1 issue comment / Phase 2 commit message references "post-#138 doc cleanup" not "#138 phase B" to keep the scope distinction clear.

---

## Estimated session count

- Phase 0 + Phase 1 (close-out comments): **session 1** (~2-3 hours, ~12 comments to draft and post; mostly mechanical extraction)
- Phase 2a + 2b: **session 2** (~3 hours; the heaviest narrative cuts)
- Phase 2c + 2d + Phase 3 + Phase 4: **session 3** (~2 hours; smaller cuts + verification + merge)

Total: ~7-8 hours of focused work spread over 3 sessions. Compresses to 2 sessions if Phase 1 + 2a are bundled.

---

## What this plan does NOT cover (deferred to future sessions)

1. **Slab native solver retirement.** Per user directive, `peierls_slab.py` survives until the parallel Green's function approach lands. That work is tracked in `.claude/plans/peierls-greens-function-approach.md` (existing scratch plan) and is out of scope here.
2. **Issue #112 actual fix.** Phase 1 relocates the empirical scan tables to #112 as the active research record. Solving the canonical-Marshak problem (the next Phase A/C session) is a separate effort.
3. **Phase 5 continuous-µ K_bc.** Issue #133 holds the retreat narrative; the function still raises NotImplementedError. If a future session resurrects it, a new tracking issue is needed.
4. **`derivations/diagnostics/` directory cleanup.** This plan only touches `docs/`. The diagnostic scripts (~150+ files) are an entirely separate cleanup target — file as a new tracking issue.

---

## Hand-off checklist for fresh next-session agent

1. Read this plan top-to-bottom.
2. `gh issue view 138` for the predecessor close-out context.
3. `git log --oneline 742d3b0..HEAD` (after #138 cascade is on main) to see if anything has landed since.
4. Run Phase 0 baseline (sphinx -W clean, ERR 30/31).
5. Start Phase 1 — open all 12 issues in browser tabs, draft each close-out comment in `.claude/scratch/`, post + verify.
6. Then Phase 2 in sub-phase order. Sphinx -W after each sub-phase commit.
7. Phase 3 + Phase 4: ship.

If at any point the audit line ranges have drifted (because someone landed unrelated doc commits between sessions), re-derive them via `grep -n "section-title-text" docs/theory/peierls_nystrom.rst` — the section TITLES are stable even when line numbers shift.
