# Audit trail — test-architect memory distillation

Two-part audit: Part A traces the ORIGINAL `.claude/agent-memory/test-architect/lessons.md`
(3 514 lines, the file this whole distillation started from) to `lessons.proposed.md`
(1 007 lines). Part B traces `MEMORY.md` §2 (49 entries) to `MEMORY.proposed.md` §2.
Verdicts follow the brief's vocabulary: KEPT · CUT (lossless) · DISTILLED into
`<meta-lesson>` · DUPLICATE of `<rule/skill § or item>` · ARCHAEOLOGY → archive
§LNN · STALE (landing named).

Part A is at BLOCK granularity (a contiguous original heading's line range),
not per-bullet: the original file already carried this two-tier shape —
families §1–§9 (the running digest) plus ~19 chronologically-appended
"additions" blocks, each one campaign's raw output not yet folded in. A prior
dispatch of this agent (interrupted by four server-side API failures, not by
its own error) did the fold that is this table's Part-A subject; this
session's own incremental action is the two rows marked **(this session)**.
Every original line's content is traceable: the 19 addition blocks were
verbatim appended to `lessons_archive.md` as sections L58–L86 by the
concluding campaigns themselves (verified: `grep -n "^## L58\|^## L67\|^## L78\|^## L86"
lessons_archive.md` finds all four at the expected lines), and a citation
sweep (`grep -oE "\`L<N>[a-z]?\`" lessons.proposed.md`, run for every L1–L86)
found at least one live pointer for every section — nothing orphaned.

## Part A — `lessons.md` (3 514 lines) → `lessons.proposed.md` (1 007 lines)

| lines | entry (first words) | verdict | where its correction lives now |
|---|---|---|---|
| 18–363 (346) | `## 1. Gates that cannot red` | KEPT — reorganized into family 1, cite-preamble ("`vv` Mode 8 / anti-#17") added, duplicates of rule clauses already excised by the prior pass | `lessons.proposed.md` §1 |
| 364–745 (382) | `## 2. Harness discipline` | KEPT — reorganized into family 2, cite-preamble ("`vv` anti-#17, `instrument-doctrine` X1") added | `lessons.proposed.md` §2 |
| 746–1013 (268) | `## 3. Config blindness` | KEPT — reorganized into family 3, cite-preamble ("`AGENT.md` §0.6, `vv` §H2/anti-#3/anti-#4/Mode 7") added | `lessons.proposed.md` §3 |
| 1014–1230 (217) | `## 4. Reference, claim layer` | KEPT — reorganized into family 4 | `lessons.proposed.md` §4 |
| 1231–1312 (82) | `## 5. Tolerance is a claim` | KEPT — reorganized into family 5 | `lessons.proposed.md` §5 |
| 1313–1667 (355) | `## 6. Carve archetypes` | KEPT by prior pass into family 6 verbatim (355 lines); **ARCHAEOLOGY → `lessons_archive.md` §L87 (this session, 2026-09-21)** — the 60+-shape lookup table is reference material, not a per-dispatch rule; digest keeps only the meta-rule (8 lines) + a pointer | `lessons_archive.md` §L87; digest meta-rule at `lessons.proposed.md` §6 |
| 1668–1707 (40) | `## 7. Snapshots, generators, and exactness` | KEPT — reorganized into family 7 | `lessons.proposed.md` §7 |
| 1708–1755 (48) | `## 8. Verifying a pure-math PRIMITIVE` | KEPT by prior pass into family 8 verbatim (78 lines incl. its obligation preamble); **ARCHAEOLOGY → `lessons_archive.md` §L88 (this session, 2026-09-21)** — the fourteen further gate shapes are reference material; digest keeps only the pillars-differ fact + a pointer | `lessons_archive.md` §L88; digest fact at `lessons.proposed.md` §8 |
| 1756–1778 (23) | `## 9. Pointers` | KEPT — unchanged as family 9 | `lessons.proposed.md` §9 |
| 1779–1883 (105) | `## CS4b additions (2026-08-21)` | DISTILLED into families 1–5/7 per its own "grouped by the families above" heading; war story independently preserved | `lessons_archive.md` §L62 (CS4a) / §L59–L66 range as cited by pointer in each folded bullet |
| 1884–1956 (73) | `## CS5 additions (2026-08-29)` | DISTILLED into families per its own grouping heading | `lessons_archive.md` §L65 |
| 1957–2030 (74) | `## P4-remainder additions (2026-08-29)` | DISTILLED into families per its own grouping heading | `lessons_archive.md` §L66 |
| 2031–2113 (83) | `## CS4c binding ladder — pre-carve (2026-08-30)` | DISTILLED into families 1–5 (cited 5×: `L67a`–`L67h` found in §§1–5) | `lessons_archive.md` §L67 |
| 2114–2176 (63) | `## FUSED step A+B (#429 / ERR-080) — pre-carve (2026-09-02)` | DISTILLED (cited 2×) | `lessons_archive.md` §L69 |
| 2177–2275 (99) | `## O(2)_a stabiliser additions (2026-09-02, #429 tracker 1.9 / #432)` | DISTILLED (cited 5×) | `lessons_archive.md` §L70 |
| 2276–2360 (85) | `## #434 R1 additions (2026-09-03, the realization carve)` | DISTILLED (cited 4×) | `lessons_archive.md` §L72 |
| 2361–2399 (39) | `## #434 R4 additions (2026-09-03, the lift-as-derivation-output carve)` | DISTILLED (cited 10×, incl. §0 meta-lesson M3) | `lessons_archive.md` §L73 |
| 2400–2503 (104) | `## #434 R3 additions (2026-09-03, the three-field registry ledger)` | DISTILLED (cited 9×) | `lessons_archive.md` §L75 |
| 2504–2568 (65) | `## #426 additions (2026-09-03, the (n,2n) anisotropy carve)` | DISTILLED (cited 6×) | `lessons_archive.md` §L76 |
| 2569–2655 (87) | `## CS4c step-5 additions (2026-09-04, ends-select-the-body carve)` | DISTILLED (cited 8×) | `lessons_archive.md` §L77 |
| 2656–2788 (133) | `## #448 additions (2026-09-05, the eigenvalue-finalize reconstruction)` | DISTILLED (cited 9×, incl. §4's own free-oracle rule) | `lessons_archive.md` §L78 |
| 2789–2864 (76) | `## CS4c step-6 additions (2026-09-07, the CS2 residue — pre-carve)` | DISTILLED (cited 5×) | `lessons_archive.md` §L79 |
| 2865–2935 (71) | `## CS4c step-6 item 6.2c additions (2026-09-07/08, axis-built moment head)` | DISTILLED (cited 7×) | `lessons_archive.md` §L80 |
| 2936–3012 (77) | `## CS4c coda additions (2026-09-08, homogeneous path re-points last)` | DISTILLED (cited 7×) | `lessons_archive.md` §L81 |
| 3013–3114 (102) | `## Consumers campaign step 1 additions (2026-09-12, Problem-identity — pre-carve)` | DISTILLED (cited 5×) | `lessons_archive.md` §L82 |
| 3115–3210 (96) | `## Consumers campaign step 2 additions (2026-09-12, the terminal object — pre-carve)` | DISTILLED (cited 11×) | `lessons_archive.md` §L83 |
| 3211–3417 (207) | `## Consumers campaign step 2 DELTA additions (2026-09-13, after the design was ruled)` | DISTILLED (cited 13×, the largest single source of §0's meta-lessons) | `lessons_archive.md` §L84 |
| 3418–3514 (97) | `## Consumers campaign step 3 additions (2026-09-17, the Solution carries its posing)` | DISTILLED — this block's own sub-headers ("→ family 1/2/3/6/7") name the fold explicitly (cited 11×) | `lessons_archive.md` §L86 |
| n/a | §0 "Eight meta-lessons" | NEW this campaign (prior pass) — cross-campaign synthesis of the above, not present in the original file at all | `lessons.proposed.md` §0; two of the eight (M3, M6) are proposed rule-uplift candidates, see `uplift.md` |
| n/a | Cite-never-restate preamble + THE SPINE + Maintenance note | NEW this campaign (prior pass) — makes the "duplicate retires, rule cited" law explicit for this file | `lessons.proposed.md` lines 1–27 |

Duplicate check performed on the still-hot families (0–5, 7, 9, ~930 lines):
read every clause of `vv-principles`, `instrument-doctrine`, `retirement-audit`,
`coding-elegance`, `algebra-of-record`, `numerical-bug-signatures`,
`nexus-verification`, `nexus-impact`, `cardinal`, `articulation`,
`process-discipline`, `workflows`, `code-search`, `nexus-tools`, `vv-testing`,
`coding-standards`, `plan-authoring` in full this session (all pasted into
context by the skill/rule loaders) and scanned the digest for entries that
restate a clause without ORPHEUS-specific mechanism. **No additional
high-confidence DUPLICATE was found beyond what the prior pass had already
excised** — every remaining entry either names a concrete ORPHEUS
file/measurement/mechanism the cited rule does not carry, or extends a cited
clause with a project-specific instance (e.g. digest §2's NEXUS-resolver-blind
entry cites no rule because none of the loaded rules states the 94.8 %
`OperatorSum` figure). This is a targeted read, not an exhaustive
per-clause cross-product; see `NEEDS`.

## Part B — `MEMORY.md` §2 (49 entries, 9 279 chars) → `MEMORY.proposed.md` §2 (49 entries, 5 492 chars)

One row per campaign name (all 49 kept — none retired, since §2's own
discipline is "active/in-flight", not archaeology). Verdict is uniform:
**KEPT, slimmed** — every entry's inline finding/refutation sentence and
artefact-path/row-count detail was cut (lossless: the finding is a lesson
already cited by the SAME `LNN` pointer in `lessons.proposed.md`, or an
artefact path already held by `active_campaigns.md`, which this section's own
intro already points at). No campaign's terminal status changed: a
`git log --oneline -30 main` spot check found the last ~30 commits on `main`
are all harness/docs/CI work (none touching `orpheus/sn/`), consistent with
every "PRE-carve" SN-campaign line still being un-landed — no stale "unmerged"
claim found. Examples (full 49-row list is the bullet-for-bullet diff between
the two files' §2, trivially re-derivable with `diff`):

| entry (name) | verdict | where its finding lives now |
|---|---|---|
| Consumers step 3 | KEPT, slimmed | finding ("no reported number changes is FALSE for 1/5") → `lessons.proposed.md` §4 `L10`, `L86a`; artefact paths/row counts → `active_campaigns.md` |
| Consumers step 2 DELTA | KEPT, slimmed | findings (0-D law not bit-identical; every fissile 0-D mixture supercritical; `A−F` never positive-stable) → `lessons.proposed.md` §3 `L7`, `L84d`, `L84e` |
| #434 R4 | KEPT, slimmed | finding (kernel Mode-12 BLIND 0/9925) → `lessons.proposed.md` §6 meta-rule pointer `L73`; full table → `lessons_archive.md` §L73 |
| #448 | KEPT, slimmed | finding (`φ = ∫ψ dΩ` is the definition) → `lessons.proposed.md` §4 `L78a` |
| *(remaining 45 entries)* | KEPT, slimmed | same pattern — each entry's own `→ LNN` pointer is unchanged and is where both the finding and the war story live |

## Counts

| file | lines before | lines after | chars before | chars after |
|---|---:|---:|---:|---:|
| `lessons.md` → `lessons.proposed.md` | 3 514 | 1 014 | 254 702 | 71 746 |
| `MEMORY.md` → `MEMORY.proposed.md` | 109 | 116 | 16 335 | 13 069 |
| `lessons_archive.md` (+ `archive_additions.md`, to be appended) | 10 906 | 11 404 | 690 000 | 725 921 |

**Rows per verdict, Part A** (28 original blocks + 2 preamble/meta additions = 30 rows):
KEPT 7 (§1,2,3,4,5,7,9) · KEPT-then-ARCHAEOLOGY-this-session 2 (§6 → L87, §8 → L88)
· DISTILLED 19 (the addition blocks) · NEW-this-campaign 2 (§0 meta-lessons,
cite-preamble). CUT (lossless) 0, DUPLICATE 0, STALE 0 — none found in the
duplicate-check pass described above (see `NEEDS`).

**Rows per verdict, Part B** (49 entries): KEPT-slimmed 49. CUT 0 (no campaign
line was archaeology — §2 is explicitly "active/in-flight", and none was
found merged/closed against git that the index still lists as open).

`MEMORY.md` line count rose (109→116) because §1/§2's intro prose grew by a
few sentences explaining where content moved (the L87/L88 split, the "findings
live in the digest" rule); its CHARACTER count fell 20 % because §2's 49
bullets each lost their inline finding/path/count. Line count is the wrong
instrument for this file's cost since every §2 bullet was already
one physical line; character/token count is what the standard's "index tax"
argument is about.

## Orchestrator's corrections at apply time (2026-09-21, the main agent with the session's context)

1. **M5 trimmed (law 2).** Its checkpoint protocol (`git status --porcelain` and `git log` at start and end; `cp -a` before the first mutation, `diff -q` after) restated `process-discipline` "Trust git for merge status" and "Mutation-testing an uncommitted file" verbatim; the entry now cites them and keeps only what this tree adds.
2. **Part B redone by git evidence (the NEEDS item).** The owner's instrument, `git log --oneline -30 main` reading no `orpheus/sn/` commits, read the harness campaign's tail, not the carves. `gh issue view` on 2026-09-21: #459, #448, #426, #434, #429, #325, #337, #2, #280, #340, #344, #290 CLOSED; #432, #235, #358, #21, #41 OPEN, and "#41" in the index was a plan-internal task number (issue #41 is a fuel issue). The CS-ladder campaign is COMPLETE 2026-09-18 @ `c27373b9` (main memory). Consumers step 3 is the one PRE-carve line that survives: its anchor file still carries 7 strict-xfail rows. §2 collapses from 49 lines to 5 (four open items and one merged line); the DSA line ("PRE-carve, runs after #280") was stale by eight weeks.
3. **§1 of the index no longer inlines the gap-landing lesson** (`L40a`): the digest's Maintenance paragraph carries it, and the index's own rule is that no lesson content is inlined.
4. The §0 "pending orchestrator review" notes became plain pointers: M1 → `AGENT.md` §1.5 (applied), M3 → `vv-principles` Mode 12 (applied), M6 → `plan-authoring` §2 AXIS-CHOICE (applied).

Counts as applied (`wc -l`): lessons.md 3514 → 1013; lessons_archive.md 10906 → 11389; MEMORY.md 109 → 66.
