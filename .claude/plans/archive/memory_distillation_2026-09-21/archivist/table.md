# Audit trail — archivist memory distillation (PROPOSE phase, 2026-09-21)

Row granularity is the **contiguous bullet block** of the current `lessons.md` (blocks are the
blank-line-separated clusters the file already uses); `n=` is the block's top-level bullet count,
so every row carries its own population. All 483 bullets over lines 16–3250 are covered, plus the
header (1–15) and the coda (3251–3256).

**Verdict vocabulary** — KEPT · CUT (lossless) · DISTILLED into `<meta-rule>` · DUPLICATE of
`<rule/skill § or item>` · ARCHAEOLOGY → archive `§L-0NN` · STALE (landing named) · REFUTED
(measurement named).

**One standing fact for every DISTILLED row:** the war story, the `[M]` numbers and the `file:line`
detail of each collapsed bullet already live in the archive section its own `→ L-0NN` pointer
names, and the proposed digest carries every one of those pointers — `[M]` all **112** archive
sections are cited by the proposal, **0** lost (checked by set difference on `L-\d{3}`). Nothing is
retired into a gap.

---

## A. `lessons.md`

| lines | entry (first words) | n | verdict | where its correction lives now |
|---|---|---|---|---|
| 1–15 | header + THE SPINE | — | KEPT, sharpened | Proposal header. The SPINE is unchanged in content; one clause added ("a line that restates a rule is retired on sight") and the pointer to AGENT.md's Quality Checklist as the itemised form of the same bar. |
| **§1 Ground truth is the LIVE tree — 226 bullets, lines 16–1568** | | | | |
| 18–34 | "RENAMED EVERYWHERE" IS A CLAIM ABOUT A *PREDICATE* | 3 | DISTILLED into §1a "A relayed claim is a PREDICATE, not a fact" | §1a bullet 1; the dead-xref-during-a-rename half and the mine-the-commit-messages half → archive L-112 |
| 36–62 | A CARVE CAN CORRUPT A HISTORY CLAIM IN THE **CODE** | 4 | DISTILLED into §1a + §1c ("two honest `[M]`s", "census the DIFF") | §1a, §1c; the `28435e11:solver.py:594` detail and the `n_inner`/`n_outer` census → archive L-111 |
| 64–85 | A "the artefacts were RE-BASELINED" claim is a claim about the TREE | 3 | DISTILLED into §1c ("when a quoted figure will not reproduce, publish the CURVE") + §2b (the private-name gate blindness) | §1c, §2b → archive L-110 |
| 87–110 | RE-DERIVE EVERY NUMERIC LITERAL … IN ONE SCRIPT | 4 | KEPT as §1c bullet 1 (and PROPOSED for AGENT.md — see `agent_md.md`) | §1c; the four instances → archive L-109 |
| 112–166 | WHEN A RELAYED `[M]` WILL NOT REPRODUCE, COMPUTE THE NORM LADDER | 8 | DISTILLED into §1c (statistic/predicate reconciliation; draw-dependence), §1f (the refuted-prediction rule), §3 (the call-signature label) | §1c, §1f, §3 → archive L-108 |
| 168–255 | WHERE a new cross-method TYPE is documented is a PACKAGE question | 13 | DISTILLED into §1a (brief-target-measures-zero), §1c (two honest counts), §1f (the three-part deferral row, EVENT-phrases), §4a (changelog hygiene) | §1a, §1c, §1f, §4a → archive L-106, L-107 |
| 257–285 | A BRIEF'S "DOCS FALSIFIED" LIST IS SCOPED TO THE CURRENT CARVE | 4 | DISTILLED into §1a + §4 ("sort by tense"; the retired-kwarg-in-a-code-block case) + §4a | §1a, §4, §4a → archive L-105 |
| 287–374 | WHEN A CARVE AXIS-IFIES A FACTOR, RE-DERIVE WHICH *ARM* | 13 | DISTILLED into §1d (re-derive F's MEMBERS; the STAND-IN), §2c (the markup regexes), §4 (retitling keeps the label), §5 (the mid-paragraph block replacement) | §1d, §2c, §4, §5 → archive L-100, L-101 |
| 376–445 | A PAGE'S OWN ⚠ CAVEAT CAN *BE* THE DIAGNOSIS A CARVE LANDS | 9 | DISTILLED into §1e (the caveat rule, kept nearly verbatim), §1a (a review's verdict is a recommendation), §1c (publish the scale-free statistic), §4a | §1a, §1c, §1e, §4a → archive L-103, L-104 |
| 451–613 | ON A SAME-COMMIT DOCS+CODE TASK, RE-RUN `git status` AFTER EVERY FILE | 22 | DISTILLED into §1b (the tree moves under you — 5 bullets), §1e (a table's CAPTION owns its columns), §1a (a forecast), §1d (census by AST), §2c, §6 | §1a, §1b, §1d, §1e, §2c, §6 → archive L-092…L-099 |
| 615–1097 | A CORPUS PASS FOR A CARVE THAT LANDED WITHOUT ONE IS A BUILD REPAIR FIRST | 74 | DISTILLED into §1c (numbers: 8 bullets), §1d (read the object: 4), §1e (surfaces nothing gates: 6), §1f (landings: 6), §1g (instruments: 4), §9 (the `-E` baseline is the carve's own red) | §1c–§1g, §9 → archive L-043, L-048…L-051, L-055, L-057, L-060, L-062…L-071, L-077…L-094 |
| 1101–1178 | A gap YOU REPORTED upward has the shortest shelf life | 10 | DISTILLED into §1b ("a gap YOU report", kept), §1a (run the brief's own census), §1f (a landing gives a phrase a real name; a refusal half-falsified), §1d (the finite-roster miss) | §1a, §1b, §1d, §1f → archive L-072, L-074, L-075 |
| 1180–1184 | Verify the memo's `file:line` ATTRIBUTION | 1 | DISTILLED into §1a ("the brief is the FLOOR" meta-rule, which is the general form) | §1 meta-rule → archive L-072 |
| 1186–1567 | A GATE'S NAME IS A UNIVERSAL — run its predicate over the whole shipped family | 58 | DISTILLED into §1d (2 bullets), §1c (4), §1e (2), §1g (2), §4 (the class split, the SYMBOLS block), §6 | §1c–§1g, §4, §6 → archive L-076…L-091 |
| **§2 The build is BLIND — 71 bullets, lines 1569–2079** | | | | |
| 1574–1586 | `:noindex:` ON AN `automodule` MINTS NO CROSS-REFERENCE TARGET | 2 | KEPT as §2a, with its count **re-measured this session** (`[M]` 48 automodules / 24 noindex, 2026-09-21); the length-changing-underline half moved to §9 | §2a, §9 → archive L-112 |
| 1588–1596 | THE PROJECT XREF GATE IS STILL BLIND TO A DEAD `:class:`… | 1 | KEPT as §2b bullet 1, verbatim in substance — it is the load-bearing ⛔⛔ and MEMORY.md's copy of it retires in its favour | §2b → archive L-111 |
| 1598–2052 | EVERY CORPUS GREP OWES `\| grep -v _build` | 65 | DISTILLED into §2a (the silent/warn split), §2b (4 bullets on the gate, nexus, page convention), §2c (6 bullets on render-vs-source), §3 (declarations), §9 | §2a–§2c, §3, §9 → archive L-002, L-044…L-047, L-052…L-054, L-058…L-061, L-068…L-071, L-078…L-082, L-088…L-090, L-104 |
| 2056–2069 | KEEP THE PRE-EDIT `-E` BUILD AND DIFF IT | 2 | KEPT, merged into §2c bullet 3 (with the self-contradicting-index half moved to §4a) | §2c, §4a → archive L-072 |
| 2071–2078 | A SOURCE regex CAN gate nested markup | 1 | KEPT, merged into §2c bullet 2 (the working pattern set) | §2c → archive L-074, L-076 |
| — | *(the "cross-doc dangling `:ref:` renders plain text" premise carried inside §2a, §3, §4 and §2b's neighbourhood)* | — | **REFUTED** `[M]` 2026-09-21 | New §2a bullet 2: a 2-page throwaway project with live-beside-dead controls shows `WARNING: undefined label [ref.ref]` at DEFAULT severity, bare and with explicit text, plus `[ref.doc]` for a dangling `:doc:`; and all five dead py-domain roles on an `.rst` page warn under `-n` while none warn by default. The label-rename caution argued from silence is void; see `uplift.md` and `agent_md.md`. |
| **§3 A `:label:` is a V&V edge — 16 bullets, lines 2080–2162** | | | | |
| 2082–2101 | A DOCUMENTED-sentinel label adds NO test | 3 | DISTILLED into §3 ("classify every label", which now carries the sentinel arithmetic) + §3's declaration bullet | §3 → archive L-077 |
| 2103–2159 | NEVER rename or delete a label a `verifies(...)` targets | 13 | DISTILLED into §3's seven bullets; the KEEP-fate rationale gains the ⚠ refutation note | §3 → archive L-003, L-004, L-024, L-027, L-030, L-032, L-035…L-037, L-039, L-049, L-064, L-065, L-070 |
| **§4 Retirement & staleness — 65 bullets, lines 2163–2582** | | | | |
| 2165–2174 | A CHANGELOG PAGE'S OWN PREAMBLE CAN LICENSE LEAVING A STALE SPELLING | 2 | DISTILLED into §4a ("a dated row keeps the spelling current on its date"; the merge-hash contract) | §4a → archive L-111 |
| 2176–2181 | A MEMBER-LIST CLAIM HAS TWO SHAPES | 1 | DISTILLED into §1c's QUANTIFIER line + §4's "charter BOTH predicates" reading inside the census bullet | §1c, §4 → archive L-105 |
| 2183–2191 | When a residue census is LARGE and the residue is a SIMPLIFICATION | 1 | KEPT as §4's DECLARE-it bullet | §4 → archive L-077 |
| 2193–2234 | A MUTATION MAGNITUDE IS A DRAW | 6 | DISTILLED into §1c ("publish the MECHANISM, not the value"), §4 (the step-label collision), §4a | §1c, §4, §4a → archive L-075, L-104 |
| 2239–2255 | A FIELD SPLIT is not a rename | 2 | KEPT as §4's split bullet (merged with the class-split instance and the remint/homonym sibling) | §4 → archive L-017, L-077, L-091 |
| 2257–2287 | An ONTOLOGY OVERTURN is not a retirement sweep | 4 | KEPT as §4 bullet 1 (merged with the design-unification-corollary case) and §3's four fates | §3, §4 → archive L-013, L-020, L-063 |
| 2289–2347 | A retirement's stale REASON outlives its stale NAME | 8 | DISTILLED into §4's stale-REASON bullet, the hidden-PARAMETER bullet, the symbol-collision bullet, and §6's ERR-chapter rule | §4, §6 → archive L-051, L-065, L-067, L-069 |
| 2349–2457 | A brief's SITE CENSUS is a sample; run the windowed CONCEPT grep | 16 | DISTILLED into §4 (census-before-repair, HISTORY as a category error, fix-half-a-claim), §1a, §9 (the flag-deleted-not-updated case) | §1a, §4, §9 → archive L-054, L-056, L-058…L-062 |
| 2459–2562 | A DELETION (unlike a MOVE) leaves a stale PARAGRAPH | 23 | DISTILLED into §4's five-register tense bullet, the tombstone bullet, the DEMOTE citation, the stale-status bullet, and §7 | §4, §7 → archive L-007, L-013…L-015, L-017…L-020, L-037, L-040…L-042, L-044…L-047, L-049, L-056 |
| 2566–2574 | A VOCABULARY retirement has FOUR classes | 1 | KEPT, folded into §4's five-register bullet (the ADDRESS and genuine-referent registers are its content) | §4 → archive L-072 |
| 2576–2581 | A section RENAME is cheap when you count citers FIRST | 1 | KEPT in §4's tombstone bullet, with the silence premise corrected per the §2a refutation | §4 → archive L-076 |
| **§5 Page surgery — 20 bullets, lines 2583–2697** | | | | |
| 2585–2679 | A mid-task scope REVOCATION on a file you already edited | 18 | DISTILLED into §5's ten bullets; the `git checkout` prohibition is a DUPLICATE of `process-discipline` § "Mutation-testing an uncommitted file" and is now a citation inside §5 bullet 1 | §5 → archive L-011, L-012, L-022, L-023, L-026…L-030, L-034, L-054, L-056, L-058, L-060…L-062, L-071 |
| 2683–2696 | A REFUSAL BECOMING A CAPABILITY is its own arc | 2 | KEPT, relocated to §6 (it is an event-class shape, not page surgery); the per-page changelog-direction half stays in §4a | §4a, §6 → archive L-076 |
| **§6 Match the doc SHAPE to the event class — 43 bullets, lines 2698–3005** | | | | |
| 2700–2715 | A NAMING PASS ADJUDICATES SENTENCES, NEVER WORDS | 2 | KEPT as §6 bullets 1–2 | §6 → archive L-112 |
| 2717–2737 | AN ERR FOR A DEFECT OF *SILENCE* EARNS ITS ENTRY | 4 | KEPT as §6 (the ERR rule, the deferral table, the arm-asymmetry dissolution); the ledger-row property clause folded into §1f | §1f, §6 → archive L-110 |
| 2739–2758 | THE ANSWER TIER IS NOT A NEW LAYER | 4 | KEPT as §6's ARITY bullet; the "shipped witnesses" tense clause → §1e, the retraction `r = σ∘π` and the generic-primitive `implements::` clauses → §3 and §6 | §1e, §3, §6 → archive L-109 |
| 2760–2766 | A LEDGER GAINING A FIELD splits across two pages by REGISTER | 1 | KEPT as §6 | §6 → archive L-091 |
| 2768–2990 | A KERNEL CHANGING HOUSE splits into TWO registers | 31 | KEPT as §6's remaining event-class bullets, one line each (kernel-changing-house, branch-becomes-one-formula, new layer, zero-consumer mint, route re-point, truncation correction, new theorem, un-weld, structural-IFF, dialectical seed page, retrodiction table, SSOT register, success-resolution chapter, ontology-overturn changelog, second-subject, stub→narrative, plus the eight one-pointer classes) | §6 → archive L-005, L-012…L-014, L-016, L-018, L-025, L-036, L-038, L-039, L-057, L-063, L-064, L-068, L-070…L-072, L-075, L-078, L-079, L-089, L-090 |
| 2994–3004 | The POSING-CONTRACT section has six parts | 1 | KEPT as §6 | §6 → archive L-072 |
| **§7 V&V vocabulary — 9 bullets, lines 3006–3049** | | | | |
| 3011–3021 | Never "MMS verifies the eigenvalue" … NAME the pillar | 1 | **DUPLICATE** of `vv-principles` § "The three pillars of verification" (the pillar table and its "MMS does not prove eigenvalues" clause) and anti-patterns #1 and #3 | §7's opening sentence cites the skill and orders a verbatim match; the never-list itself is the skill's → archive L-010 |
| 3022–3030 | Never upgrade a `foundation` gate to an L-level in prose; a doc sentence "gates X, Y pin claim C" IS a coverage claim | 2 | KEPT (the level definition is `vv-principles` § "V&V level taxonomy", ORTHOGONAL; the *don't-upgrade-in-prose* imperative and the coverage-claim-in-PROSE analogue are in no skill) | §7 → archive L-040, L-047 |
| 3031–3035 | The SAME gate cited for TWO claims | 1 | KEPT | §7 → archive L-047 |
| 3036–3039 | Distinguish the EUCLIDEAN transpose from the metric HILBERT adjoint | 1 | KEPT | §7 → archive L-010, L-034 |
| 3040–3043 | A Mode-10 sub-floor term is closed by STRUCTURAL teeth; Mode-12 blindness boundary; the FIRST iterative member | 3 | **DUPLICATE** of `vv-principles` test-design modes 10 and 12 (mode 10 carries the "no value-improvement leg — do not manufacture one" clause verbatim; mode 12 carries the stabiliser boundary and the leaf-transpose exception) | §7's opening pointer → archive L-010, L-015, L-036 |
| 3044–3046 | Skill-uplift duty | 1 | KEPT (it is AGENT.md Directive 5's operational half) | §7 → archive L-010 |
| **§8 Code-prose rebalance — 6 bullets, lines 3050–3080** | | | | |
| 3052–3077 | Expect ZERO MOVED … prove the edit is doc-only | 6 | KEPT, all five rules, tightened | §8 → archive L-033, L-034, L-041, L-045, L-073 |
| **§9 Gates, generated artefacts, tooling — 27 bullets, lines 3081–3250** | | | | |
| 3083–3090 | `:label:` ON `.. math::` IS SPHINX-ONLY | 1 | DISTILLED into §2c's differential-docutils bullet (where the subtraction rule belongs) + §1d's probe ladder | §1d, §2c → archive L-076, L-096, L-105 |
| 3092–3209 | Generated artefacts are NEVER hand-edited | 19 | KEPT as §9's fourteen bullets; three collapse as DUPLICATE: the xref-gate `head_role` mechanics (→ §2b, one home) and the `-F`-not-`-m`/backtick family (→ `process-discipline` § "A commit message with backticks", cited from §9's heredoc bullet) | §2b, §9 → archive L-006, L-008, L-025, L-026, L-027, L-029, L-030, L-031, L-035, L-041, L-051, L-054, L-059, L-061, L-062, L-063, L-068, L-069, L-081, L-082, L-089, L-091, L-102, L-103 |
| 3211–3247 | A CAPABILITY FLIP stales DEFERRAL CONTRACTS | 7 | KEPT, relocated to §4 (a staleness class, merged with the stale-status blast radius) and §9 (the AST doc-only proof, the `raise`-string rule, the `git archive` baseline) | §4, §9 → archive L-051, L-070, L-073 |
| 3251–3256 | Quality self-assessment (Directive 3) | — | KEPT, 3 lines; the rubric list is AGENT.md Directive 3's, cited | Coda → AGENT.md Directive 3 |

### Lines retired with no surviving digest line

**None.** Every block above resolves either to a named proposal section, to a named rule/skill
clause, or to the archive section its own pointer already named. The three DUPLICATE rows (§7's
pillar list, §7's mode 10/12 boundaries, §9's `-F` backtick family) are the only content whose
home is now exclusively a rule or skill, and each row names the heading — read, not recalled:
`vv-principles` § "The three pillars of verification" and its anti-patterns #1/#3; `vv-principles`
§ "The 6 test-design failure modes", items 10 and 12; `process-discipline` § "A commit message with
backticks is written with `-F`, never `-m`".

---

## B. `MEMORY.md`

| lines | entry (first words) | verdict | where its correction lives now |
|---|---|---|---|
| 1–11 | header + index disciplines | KEPT, re-worded | Proposal header; discipline (4) restated as "a landed campaign is ONE line" and given its reason |
| 13–30 | §1 Lessons — hot digest over cold archive | KEPT | Proposal §1, plus the next free archive number (**L-113**), which the brief's "after L86" guessed from another agent's tree |
| 32–45 | §2 header: "⚠ This list is a SNAPSHOT and it has frozen on landed work EIGHT times" | **STALE — it froze a ninth time.** | Proposal §2 replaces the list with a single measured statement. `[M]` 2026-09-21: `git status --porcelain -- docs/` = **0**; `1ce64371`, `628997b1`, `deacd897`, `2c1667b0`, `7b4d2b78`, `6379e9ab`, `0c6978a4`, `3b8e0591`, `bcd9c83c`, `207d2b07`, `93225c65`, `0221e19e` are all ancestors of `main`; `git branch -a` lists **no** `refactor/consumers-*` branch, local or remote |
| 46–67 | rows #412 · step-3 U6 · U2 · U1 · step-2 C3b-2 · C3b-1 · C3a — each "uncommitted on branch X" | **STALE (landing named)** | All seven merged; the #412 prose half is `0c6978a4` and the typed residue `3b8e0591`, both on `main`. Their lessons are digest §1–§9 over archive L-106…L-112; the branches are gone. Collapsed into proposal §2's one-line MERGED statement |
| 52–55, 59–61 | "ERR-086 … its catcher is owed by the main agent, so `nexus errors` reads 1 uncaught"; "ERR-085 minted, catcher marker owed" | **STALE** | Both landed: `tests/sn/solve/test_subcritical_multiplying_source.py:169,189` carry `catches("ERR-086")` and `tests/sn/sweep/core/test_cache.py:1000` carries `catches("ERR-085")`; the generated catalogue index reads **86 entries · 320 catchers · 0 uncaught**. Recorded as a ⛔ in proposal §2 so the pattern is not repeated |
| 68–84 | rows C2 · C1 · step 1 · CS4c ×7 · #425 · #448 | ARCHAEOLOGY → archive L-095…L-105 (each row already names its own section) | Proposal §2's pointer to `git log --oneline -- docs/` and to the archive span. The per-row diffstats (`9 .rst, +642/−73`) are `git show --stat`, not memory |
| 85–88 | "Everything older — #428, #434, #432, #429 …" with `[[lessons-L39]]`…`[[lessons-L94]]` wikilinks | CUT (lossless) | Those wikilinks resolve to no file; they were archive-section pointers in wikilink clothing. Proposal §2 names the archive directly (`L-095…L-112` for 2026-09, and the archive itself for older) |
| 89–90 | ⚠ ERR-026 history block, branch gone | KEPT as a ⛔ void-claim note | Proposal §2 bullet 3 — the branch is confirmed absent from `git branch -a`, so the "still OPEN" claim asserts nothing and must not be acted on |
| 91–97 | ⛔⛔ the dotted-target blindness note + its acceptance rider | KEPT, RELOCATED | Proposal §2 bullet 4 points at digest §2b, which now holds it once, with its measurement and its two controls. An index is not the home of a measured lesson |
| 98–104 | §3 header + landed-milestone record | KEPT | Proposal §3 intro and §2's "where the record lives", including the ⚠ that the SN history page is the post-split `methods/sn/history.rst` and that a stale-ref grep must `test -f` the source |
| 105–153 | §3's 23 reusable doc-architecture pointers | KEPT, hooks tightened to ≤ 15 words | Proposal §3. `[M]` all **25** markdown links in the proposal resolve to existing files; **0** wikilinks remain |
| 154–159 | "Doc-architecture redesign (#231, OPEN) … Phase 2 ACTIVE" | **STALE** | `[M]` 2026-09-21 `gh issue view 231 --json state` → **CLOSED**. Proposal §3's last bullet states CLOSED, keeps the conventions as the standing target, keeps the `phase2_code_prose` map pointer and the five file-classes pointer (digest §8), and says to file a new issue rather than reopen |

---

## C. Topic files judged archaeology (listed only; not retired here)

`[M]` 2026-09-21: **80** files in `.claude/agent-memory/archivist/` besides `MEMORY.md`,
`lessons.md` and `lessons_archive.md`. **23** are cited by the proposed index §3 (the reusable
recipes); **57** are uncited by it — and were uncited by the current index too, so their
per-dispatch cost is already zero.

⚠ **The predicate is by NAME, not by reading:** a file whose name encodes a phase, task, wave, PR
or issue-leg codename and no doc-shape noun. I did not open these 57 files this session, so this is
a candidate list for the orchestrator's blast-radius audit, not a verdict. **28 of 57** match:

`425_outside_chapter` · `425_sn_chapter` · `feedback_bc_trace_law_wave_12` ·
`feedback_d5b_d6_campaign_closeout` · `feedback_err058_success_closeout_supersedes_phase_chain` ·
`feedback_issue_196_eigenvalue_verification_closeout` · `feedback_issue_247_legA_mode10_closeout` ·
`feedback_issue_251_legB_boundary_trace_closeout` · `feedback_issue_257_consolidated_taxonomy_pass` ·
`feedback_issue_257_s9_coherent_promise_criterion_seam` · `feedback_misc_cluster_triage` ·
`feedback_phase2a_large_relocation` · `feedback_phase2b_label_dense_partition` ·
`feedback_phase2c_staleness_sweeps` · `feedback_phase2d_table_preservation_constraint` ·
`feedback_phase_c_stub_expansion` · `feedback_phase_d_carlson_seed_narrative` ·
`feedback_plan_triage` · `feedback_plan_triage_quadrature_cluster` ·
`feedback_post_wave_cleanup_docs` · `feedback_s10a_emission_spectrum_value_object` ·
`feedback_wave_o_operator_algebra_docs` · `issue_196_pr_cleanup_docs_closeout` ·
`issue_196_pr_index_6_closeout` · `peierls_greens_phase123_rich_narrative` ·
`task55_transport_sweep_docs_pass` · `task57_psi_bc_matvec_docs_pass` ·
`task58_staleness_sweep_fix`.

The other **29** carry a doc-shape or subject noun (`feedback_label_reconciliation_sweep`,
`feedback_unified_trace_rewrite`, `feedback_type_retirement_concept_survives`, …) and may hold a
reusable recipe that merely lost its index line; they are neither promoted nor retired here.

**Blast radius, for whoever runs the audit:** `[M]` 2026-09-21
`grep -rn "agent-memory/archivist" docs/ .claude/plans/ | grep -v _build | wc -l` → **0**, which
reproduces the brief's measurement. The memory-file retirement law still applies
(`feedback_memory_distillation_standard.md`: agent-memory files are Nexus nodes, so a delete owes
the full `retirement-audit` blast radius, not a wikilink check) — and a 2026-06 `git rm` of 333
closeouts broke 98+ in-repo references and had to be reverted.

---

## D. Counts

**Lines before / after**

| file | before | after | ratio |
|---|---:|---:|---:|
| `lessons.md` | 3 256 lines / 272 956 bytes | 720 lines / 56 574 bytes | 4.5× lines, 4.8× bytes |
| `MEMORY.md` | 159 lines / 12 206 bytes | 113 lines / 8 037 bytes | 1.4× lines, 1.5× bytes |
| `lessons_archive.md` | 13 714 lines (cold) | +1 section (**L-113**, 60 lines) — see `archive_additions.md`; nothing the distillation retired needs one | — |

**Rows per verdict** — `[M]` counted by parsing this file's own tables, after the parse caught two
errors in the numbers I had first written by hand (a `n=2` that is `n=1`, and a five-STALE tally
that is four). The mis-parse of the `` `\| grep -v _build` `` row was the PARSER's fault, not the
table's: the pipe is escaped and renders correctly.

Section A: **53** rows over **483** bullets (the `n` column sums to 483). Section B: **12** rows.

| verdict | A rows | B rows |
|---|---:|---:|
| KEPT (verbatim in substance, tightened, or relocated) | 25 | 6 |
| DISTILLED into a named meta-rule | 24 | 0 |
| DUPLICATE of a named rule/skill clause | 3 | 0 |
| ARCHAEOLOGY → a named archive section | 0 | 1 |
| STALE (landing named) | 0 | 4 |
| CUT (lossless) | 0 | 1 |
| REFUTED (measurement named) | 1 | 0 |
| retired with no surviving home | **0** | **0** |

**Coverage checks**

- `[M]` archive sections cited by the proposed digest: **112 of 112**; lost: **0**; bogus (a
  pointer naming no section): **0** — set difference on `L-\d{3}` between the proposal, the current
  digest and `^## L-\d{3}` in the archive.
- `[M]` markdown links in the proposed index: **25**, all resolving to existing files; wikilinks
  remaining: **0**.
- `[M]` lines over 100 columns in either proposal: **0** (the corpus's own width).
- Bullets covered by section A's rows: **483 of 483**, plus lines 1–15 and 3251–3256.

## Orchestrator's notes at apply time (2026-09-21, the main agent with the session's context)

Applied as proposed. Additions of the orchestrator: pointer notes in the digest for what moved (checklist items 0 and 6 into `AGENT.md`; the `-n` scope into `retirement-audit` A.2; the prose coverage claim into `vv-principles`; the input-count clause into the `instrument-doctrine` skill; the unquoted-heredoc check into `code-search`, whose digest bullet shrank to its founding case); `L-113` added to the refutation bullet's pointers; the index's next free number set to L-114 after the append. The two `AGENT.md` corrections were accepted on the transcript in §L-113 and on Sphinx's documented behaviour (an undefined label warns; nitpicky mode reports unresolved domain roles on rendered pages). Part C's 28 archaeology candidates are not retired here (a topic file is a graph node; its retirement is the retirement-audit's blast radius, a separate pass).

Counts as applied (`wc -l`): lessons.md 3256 → 721; lessons_archive.md 13714 → 13793; MEMORY.md 159 → 113.
