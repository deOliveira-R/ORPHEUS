# Audit table — explorer memory distillation, 2026-09-21

Every verdict below was made by reading the source named (`AGENT.md`,
`.claude/skills/nexus-exploring/SKILL.md`, `.claude/skills/nexus-guide/SKILL.md`,
the brief's five rule lines and the rule files behind them). The war story of
EVERY lesson moves verbatim to `lessons_archive.md` §L-NNN (`archive_additions.md`,
a lossless move: one hunk of `diff`, the `L-034b` heading). The verdict column
says what the DIGEST does with the entry; "lines" are `lessons.md` at HEAD.

## `lessons.md` (1 641 lines → `lessons.proposed.md` 406 lines)

| lines | entry (first words) | verdict | where its correction lives now |
|---|---|---|---|
| 1–21 | preamble: "Behavioral corrections only … spine … codified as OP 4–7" | CUT (lossless) | digest preamble (rewritten: archive pointer, stable IDs, the `[uplift]` convention) |
| 24–67 | L-001 A retirement/rename blast radius = graph callers AND grep AND constructors AND doc nodes | DUPLICATE of `AGENT.md` Operating Principles §4 (the four searches, verbatim there) | digest L-001 stub → OP4; the two graph-blinding constructs (`@singledispatchmethod` alias, Protocol receiver) distilled into M-4 and kept in the stub; war story archive §L-001 |
| 69–87 | L-002 The issue text is a stale premise | DUPLICATE of `AGENT.md` OP §5 | digest L-002 stub → OP5; M-1 generalises; archive §L-002 |
| 89–109 | L-003 Separate the DURABLE subsystem-shape from line numbers | DUPLICATE of `AGENT.md` OP §7 (+ `MEMORY.md` preamble for home placement) | digest L-003 stub (home-placement sentence kept); archive §L-003 |
| 111–128 | L-004 A clean carve verdict names BOTH the retire case and the keep-as-anchor case | KEPT (cut 18 → 5 lines) | digest L-004 → M-6 |
| 130–145 | L-005 Git is authoritative for merge-status | DUPLICATE of `AGENT.md` OP §6 and `process-discipline` § "Trust git for merge status" | digest L-005 stub; archive §L-005 |
| 147–162 | L-007 On a branch under ACTIVE edit, re-run the census | DISTILLED into M-2 | digest L-007 (three-step imperative) + M-2; campaign narrative (F2 recon, `refactor/pyright-burndown`) archive §L-007 |
| 164–172 | L-008 zsh: an unquoted separator starting with `=` | KEPT (cut) — not in `code-search` (read § "`grep` here is ugrep": quote-eating, unquoted heredoc, `$var` only) | digest L-008, marked `[uplift]` → `uplift.md` U-1 |
| 174–207 | L-009 A dataclass-FIELD rename audit is a grep problem … substring of an English word | DISTILLED into M-3 (anchored vs unanchored) and M-4 (fields mint no node) | digest L-009 (anchored greps + the `\b` delta hazard); archive §L-009 |
| 209–229 | L-010 "Complement" ≠ "the named sibling": a THIRD bucket | KEPT (cut) | digest L-010 → M-5 |
| 231–249 | L-006 A "shape probe" is not always a missing predicate | KEPT (cut); reordered numerically | digest L-006 → M-6 |
| 251–303 | L-012 On a "blast radius ahead of a carve" brief, `git diff --stat` FIRST (+ two sharpenings) | DISTILLED into M-2 | digest L-012 (open/close protocol, `git ls-files --error-unmatch`, runtime-probe reconciliation, emptiness greps at close); B3.4c / #326 map / `roots_of_unity` narrative archive §L-012 |
| 305–338 | L-011 A docstring that DELEGATES is the highest-yield falsity shape | KEPT (cut) | digest L-011 → M-1 |
| 340–384 | L-013 SWAP IT AND RUN — a grep-classification of exact assertions is guesswork | KEPT (cut 44 → 12 lines) | digest L-013 → M-5 |
| 386–429 | L-014 Read the source's DERIVATION and its INDEX-DOMAIN sentence | KEPT (cut) | digest L-014 |
| 431–470 | L-015 Enumerate the FUNCTIONALS … fold of the ALGORITHM or of the STATE | KEPT (cut) | digest L-015 → M-6 |
| 472–514 | L-016 A stored NUMERIC PROPERTY is a claim: sweep it … what the SYMMETRY gives for free | KEPT (cut) | digest L-016 → M-1, M-5 |
| 516–559 | L-017 A NON-target sibling shares the target's name; a test's self-description | DISTILLED into M-3 (homonym) and M-4 (captured function object) | digest L-017 (two-number probe, RHS one hop, `.npz` find); archive §L-017 |
| 561–589 | L-018 Which rows of a STATIC TABLE are consulted; the discriminating fixture FAILS | KEPT (cut) | digest L-018 → M-5 |
| 591–637 | L-020 The BRIEF's own timeline is a claim; "cannot express X" EXPIRES | DISTILLED into M-1 and M-2 | digest L-020 (the `git log -1` + mtime command, read-the-new-module, domain/codomain tiers); archive §L-020 |
| 639–681 | L-019 Hunting a hidden TRANSFORMATION: the chart-defining ASSIGNMENT; COUNT the parts | KEPT (cut) | digest L-019 |
| 683–726 | L-023 "N spellings" is a SYMPTOM: the primitive that DISCARDED the information | KEPT (cut) | digest L-023 → M-1 |
| 690–701 | L-023 finding "`power_iteration` returns a 3-tuple with no flag" | STALE — #340 CLOSED (gh, 2026-09-21); `IterationRecord` landed | the lesson (read the callee's `return`) is kept without the finding; archive §L-023 |
| 728–757 | L-022 The campaign's own MID-FLIGHT PROSE is the staleness class | KEPT (cut) | digest L-022 → M-1 |
| 747–752 | L-022 "(a) … eps-gap gate … still unfixed 7 days later" | ARCHAEOLOGY — #325 CLOSED; a dated gap status | archive §L-022 |
| 759–813 | L-021 A brief's TYPE table is a claim about MATERIALIZED objects; INERT, not SAFE | KEPT (cut) | digest L-021 → M-1, M-5, M-6 |
| 815–858 | L-024 A solver's NESTING SHAPE is per-ENTRY-POINT; trace it by RUNNING | KEPT (cut) | digest L-024 → M-4, M-5; the #340 narrative archive §L-024 |
| 860–899 | L-026 "DELIBERATE?" is decided by a COUNTERFACTUAL; the unit is the SOLVE | KEPT (cut) | digest L-026 → M-5 |
| 901–963 | L-025 FRAME-LOCALS probe; a `[M]`-marked NEGATIVE claim is the most perishable | KEPT (cut); ⭐1 restates `plan-authoring` §2 [M]-SCOPE, which NO brief line carries → not retired | digest L-025, marked `[uplift]` → `uplift.md` U-3 |
| 965–1005 | L-027 Diff the PRIMITIVE the survey is about; a count that reproduces under no convention | KEPT (cut) | digest L-027 → M-1, M-2, M-3 |
| 1007–1052 | L-028 Sort consumers by WHICH GUARD they sit behind | KEPT (cut; the table dropped, its three classes kept in prose) | digest L-028 → M-5 |
| 1054–1097 | L-029 MEASURE the capability against the cheaper alternative; per TIER | KEPT (cut) | digest L-029 → M-5, M-6 |
| 1099–1134 | L-030 Count the seam's `return`s first — "the ONE construction site" is a claim | KEPT (cut) | digest L-030 → M-1, M-3 |
| 1103–1110 | L-030 finding "three `return Solution(...)` sites bypass `_package_solution`" | STALE — `[M]` 2026-09-21 `grep -n 'return Solution(\|return _package_solution' orpheus/sn/solver.py`: 5 `return _package_solution(`, 0 `return Solution(` | digest L-030 carries `[REMEDIED]`; the topic file `sn_solve_exit_and_reflective_default.md` still states it present-tense (owner audit) |
| 1136–1162 | L-031 Cite an anchor by GREPPING its string, never from the `sed` range | KEPT (cut) | digest L-031 |
| 1164–1211 | L-032 The authored rationale sits on a NEIGHBOURING label; NOTHING fails flatteringly | KEPT (cut); absorbs the method corrections the `MEMORY.md` §3 hooks carried (CP `pytestmark` comments 15 of 15; MC/operator-algebra KIND prior) | digest L-032 |
| 1213–1282 | L-033 The CODE may declare the label; the `verifies()` CLAIMANT; the principled NOTHING | KEPT (cut 70 → 10 lines) | digest L-033 |
| 1284–1318 | L-034 "Doctrine X is overturned" needs a TWO-SIDED inventory | KEPT (cut) | digest L-034 → M-6 |
| 1320–1363 | L-035 OPENS with `git log --since`; a cited site vs its CONSUMERS | DISTILLED into M-1, M-2 | digest L-035 ("CONFIRMED (substance) / REFUTED (site)"); archive §L-035 |
| 1365–1409 | L-036 A RETAINED section keeps its ORIGINAL vintage | DISTILLED into M-1, M-2, M-6 | digest L-036 (per-SECTION `--since`, "conflict to adjudicate"); archive §L-036 |
| 1411–1442 | L-034 (second) Reconciling a PLAN includes the OTHER PLANS; `git log -S` | KEPT (cut); RENUMBERED L-034b (no explorer file cites L-034/035/036: `grep -n 'L-03[4-6]' explorer/*.md` → 0 outside `lessons.md`) | digest L-034b → M-1; archive §L-034b |
| 1444–1470 | L-037 The AST route has its own viewport; a same-named field/method pair | KEPT (cut); half (a) is an instance of X1 as the brief line carries it — cited, not retired | digest L-037 → M-3 |
| 1472–1497 | L-038 Split EXECUTABLE sites from PROSE citations in the FIRST pass | KEPT (cut); restates `instrument-doctrine` X2 check "a membership question is parsed with an AST, never grepped" which the brief line does NOT carry → not retired | digest L-038, marked `[uplift]` → `uplift.md` U-2 |
| 1490–1493 | L-038 (b) "the CLI `nexus impact` prints its depth-1 list then dies with a traceback" | ARCHAEOLOGY — a tool-state observation of 2026-09-02 with no behavioural rule beyond "use the MCP tool" | archive §L-038 |
| 1499–1511 | L-039 Nexus `callers` on a METHOD returns empty + `unresolved`; read the COUNT | KEPT (cut); the base fact is the skill's (`nexus-exploring/SKILL.md:78` "`callers` empty ≠ dead: read the `unresolved` block") — cited; the addition (count = census, pair with AST) is the lesson | digest L-039 → M-3, M-4 |
| 1513–1540 | L-040 A retirement's PROSE census has three populations | KEPT (cut) | digest L-040 → M-3 |
| 1515 | L-040 `scratch/_step6/explorer_prose_census.md` pointer | ARCHAEOLOGY — a scratch path (exists today, untracked; will not survive) | archive §L-040 |
| 1542–1562 | L-041 An `is`-identity on an OPTIONAL constituent is vacuous on the `None` arm | KEPT (cut) | digest L-041 → M-5 |
| 1546–1553 | L-041 finding "`is_same_phase_space` … `True` for different meshes" | STALE — #459 CLOSED; `[M]` 2026-09-21 `grep -rn same_phase_space orpheus/` → `SNProblem.same_phase_space`, comment "Until 2026-09-12 `is_same_phase_space` compared CONSTITUENT …" | digest L-041 carries `[REMEDIED 2026-09-12]`; two topic files and two §3 hooks still name the old predicate present-tense (owner audit) |
| 1564–1582 | L-042 Dataclass introspection LIES about hashability | KEPT (cut) | digest L-042 → M-5 |
| 1584–1605 | L-043 A kwarg's census is a RECEIVER split; the clamp is per-ENTRY | KEPT (cut) | digest L-043 → M-3 |
| 1607–1633 | L-044 A "bit-identical" verdict = ULP probe + the WALLS' PREDICATE; spy on the CALLER's binding | KEPT (cut) | digest L-044 → M-4, M-5 |
| 1635–1641 | L-045 A `Class.attr (file:line)` citation names a LINE; a `\.method\(` grep counts docstrings | KEPT (cut) | digest L-045 → M-1, M-3 |

Rows: 53. KEPT (cut) 35 · DISTILLED 7 · DUPLICATE 4 · CUT 1 · STALE 3 · ARCHAEOLOGY 3.
(46 lesson headings + the preamble = 47 entry rows; 6 sub-rows for a finding inside a kept entry.)

## `MEMORY.md` (152 lines → `MEMORY.proposed.md` 55 lines)

| lines | entry (first words) | verdict | where its correction lives now |
|---|---|---|---|
| 1–9 | preamble "One line per entry … 2026-08-11: that block had grown to 47 lines" | CUT (lossless) | proposed preamble (adds the ≤ 15-word hook cap and the graph-node retirement caveat); the 2026-08-11 anecdote is archaeology |
| 13–19 | §1 `lessons.md` line carrying one-sentence restatements of L-039…L-045 | CUT — violates the preamble's own "never restate a lesson's content here" | digest headings; §1 becomes two pointer lines (digest, archive) |
| 22–26 | §2 preamble "Merge-status in memory goes STALE (L5) …" | DUPLICATE of `AGENT.md` OP §6 | one sentence pointing at OP6 |
| 28–34 | §2 "Every SN campaign … MERGED … Only surviving local SN branch: `feature/sn-adjoint-transport` (the paused #276 campaign)" | STALE — `[M]` 2026-09-21: `git branch -a \| grep -i adjoint` → none; `gh issue view 276` → CLOSED; `git merge-base --is-ancestor` yes for `b23e972e` and `588f2429`; the "#34 ray-leg retirement" number is a diffusion issue (OPEN), i.e. a task id, dropped per `process-discipline` "task numbers COLLIDE" | ONE merged line naming the two hashes, #276 CLOSED and its record `.claude/plans/archive/adjoint_completion_campaign.md` (`9ff1ee3f`) |
| 35–39 | §2 "Durable post-#280 facts: A_BB = `RadialCharacteristicOperator` … walk executors …" | KEPT one line, flagged — verified true (`radial_characteristic.py:565/:570` call `carlson_inward_sweep_from_source`; `_OneDimScanWalk` at `loss_representation/__init__.py:3032`) but it is subsystem SHAPE, not active state | `agent_md.md` item 5 (durable-shape section); the index line deletes when that lands |
| 41–44 | §3 preamble | CUT | proposed §3 preamble (vintage caveat) |
| 45–51 | §3 CP hook (7 lines, ⭐ method correction "pytestmark COMMENTS answered 15 of 15") | CUT to ≤ 15 words; the method correction is a lesson | digest L-032 |
| 52–57 | §3 operator-algebra hook (KIND split) | CUT to ≤ 15 words; the KIND split is a lesson | digest L-032 |
| 58–63 | §3 spatial transform hook | CUT to ≤ 15 words | topic file |
| 64–67 | §3 angular layer hook | CUT to ≤ 15 words | topic file |
| 68–72 | §3 quadrature landscape hook "`level_symmetric_sn` … (tag says N−1)" | STALE — `[M]` 2026-09-21 `orpheus/numerics/quadrature/registry.py:665-669`: "The wrapper has recorded a build-measured degree since #337" | hook tensed ("LS degree tag since build-measured (#337)"); topic file to be tensed (owner audit) |
| 73–77 | §3 convergence-knob hook | CUT to ≤ 15 words; #364 verified OPEN | topic file |
| 78–83 | §3 SN solve exit hook "the exit is THREE `Solution`-construction sites … BYPASS `_package_solution`" | STALE (same measurement as the L-030 row) | hook tensed ("exit sites since folded"); topic file to be tensed |
| 84–88 | §3 MoC hook | CUT to ≤ 15 words | topic file |
| 89–93 | §3 loss_representation hook | CUT to ≤ 15 words | topic file |
| 94–99 | §3 Monte Carlo hook (KIND prior ⭐) | CUT to ≤ 15 words; the KIND prior is a lesson | digest L-032 |
| 100–103 | §3 multigroup hook | CUT to ≤ 15 words | topic file |
| 104–107 | §3 GENDF hook | CUT to ≤ 15 words | topic file |
| 108–113 | §3 HarmonicMomentField + harmonic-frame hooks | CUT to ≤ 15 words | topic files |
| 114–119 | §3 flux torsor hook "⚠ snapshot 2026-08-19 pre-ruling — if the cone campaign landed since, this maps the RETIRED side" | STALE — `[M]` 2026-09-21 `grep -rln 'FluxDisplacement\|class FluxRole' orpheus/` → 0; `AGENT.md` durable-shape: "campaign 1 CS3, 2026-08-19, retired the affine `FluxRole` gate and the `FluxDisplacement` mint" | hook reads "maps the RETIRED torsor side (CS3 landed 2026-08-19)" |
| 120–129 | §3 affine-operator + SN α ends hooks | CUT to ≤ 15 words | topic files; the α-ends "bare `assert`" claim is UNVERIFIED (NEEDS) |
| 130–142 | §3 cylindrical + non-SN census hooks | CUT to ≤ 15 words; #326 verified CLOSED | topic files |
| 143–147 | §3 phase 5 + pyright hooks | CUT to ≤ 15 words | topic files |
| 148 | §3 Problem → Solution hook "`is_same_phase_space` is not a content identity (L-041)" | STALE (renamed `same_phase_space`, #459 CLOSED) | hook cut to the durable half; topic file to be tensed |
| 149 | §3 identity-step hook "every generating datum … RAISES on `==`/`hash` today; exactly ONE `is_same_phase_space` pin, and it FLIPS" | STALE (same) | hook "identities as of 2026-09-12 (#459 CLOSED)" |
| 150–151 | §3 k-solve walls + census-predicates hooks | CUT to ≤ 15 words | topic files |

Rows: 26. CUT 18 · STALE 6 · DUPLICATE 1 · KEPT (flagged) 1.

## Topic files judged archaeology (listed only; the orchestrator's blast-radius audit decides)

`[M]` 2026-09-21, `grep -rl <file stem> docs/ .claude/plans/` → 0 for each; the
single memory reference each carries is `explorer/MEMORY.md` itself.

| file | reason |
|---|---|
| `flux_torsor_vs_cone_inventory.md` | maps the RETIRED torsor side; CS3 landed 2026-08-19; its "⚠ pre-ruling" caveat has resolved |
| `identity_step_census_durable.md` | #459 CLOSED; the predicate it censuses (`is_same_phase_space`) is renamed and repaired |
| `problem_solution_split_census.md` | a campaign OPENER census; the split has landed (the `SNProblem` hub exists) |
| `sn_solve_exit_and_reflective_default.md` | its headline (three exit sites) is remedied; the unset-BC default half may still hold |
| `phase5_mu_resolved_primitive_inventory.md` | a "Phase 5" codename inventory of one file; no lesson, no convention |

## Counts

| file | before | after |
|---|---|---|
| `lessons.md` → `lessons.proposed.md` | 1 641 | 406 |
| `lessons_archive.md` (new) → `archive_additions.md` | 0 | 1 636 (16-line header + 1 620 body lines moved verbatim; one heading renamed) |
| `MEMORY.md` → `MEMORY.proposed.md` | 152 | 55 |

Verdict totals over both tables (79 rows): KEPT 36 · CUT 19 · DISTILLED 7 ·
STALE 9 · DUPLICATE 5 · ARCHAEOLOGY 3. Every DUPLICATE names its source
heading; every STALE names the command that measured it; every retired
war-story line is in `archive_additions.md` under its own L-number.

## Orchestrator's notes at apply time (2026-09-21, the main agent with the session's context)

Applied as proposed; the lossless split verified by `diff` (one hunk, the `L-034b` heading). Accepted from `agent_md.md`: OP5 broadened (item 1), OP8 added (item 2), the `:446` line number dropped (item 4), the curvilinear block added to the durable-shape section (item 5) and the index's curvilinear line deleted with it. OP9 (item 3) NOT added: U-4's open/close protocol went to the workflows page's brief template instead, one home. Accepted from `uplift.md`: U-1 as a `code-search` clause; U-2 into the `instrument-doctrine` brief line; U-3 into the `articulation` brief line; U-4 into the brief template (L38's sentence); U-5 into the brief template as the fixed sentence before the generated list. The digest's `[uplift]` marks became pointers to those homes. The five topic files judged archaeology are not retired here (a topic file is a graph node; its retirement is the retirement-audit's blast radius, a separate pass); the `sn_alpha_and_psi_half_ends.md` hook's bare-`assert` claim stays UNVERIFIED as the owner reported.

Counts as applied (`wc -l`): lessons.md 1641 → 405; lessons_archive.md 0 → 1636 (new); MEMORY.md 152 → 54.
