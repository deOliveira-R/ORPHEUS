# QA Agent Memory — index

## 1. Lessons — a HOT digest over a COLD archive

Two files, read at different times. Do **not** summarize lesson content here.

- **[lessons.md](lessons.md)** — the HOT digest (**1076 lines** at 2026-09-21 — ⛔ OVER the distillation trigger by 2.7x; distil BEFORE the next append). Behavioral rules
  only: one imperative rule + its failure→correction core + a
  `→ lessons_archive.md L-0NN` pointer. Nine sections: **A** mutation mechanics,
  **B** structural blindness, **C** structural independence, **D**
  re-baseline/bit-identity, **E** markers & audit surface, **F** claim-scope,
  **G** doc-correctness, **H** mechanics/environment, **I** the map of what is
  already in `vv-principles`/`numerical-bug-signatures` (point, don't restate).
  **Read this before every review.**
- **[lessons_archive.md](lessons_archive.md)** — the COLD archive (~5940 lines,
  L-001..L-089, append-ordered). War stories, evidence, `file:line`, measured
  tables, verdicts. **Open only the specific `L-0NN` a digest rule points at** —
  never read it whole (it is ~60K tokens).

Maintenance: a new lesson appends `L-0NN` to the archive AND lands a 2–5 line
rule in the digest. Sharpen the digest in place; **it is now past the ~400-line
distillation trigger — distill before the next append** (encode a shared
meta-lesson once, list its instances); never truncate.

## 2. Active / in-flight state

⏹ Nothing in flight. Every entry below is COMPLETE: one line each — the verdict
lives in its report, the behaviour in its lesson, the tense in git.
⛔ **Every branch/hash/line below is a SNAPSHOT — reconcile with git and a grep
before acting on any of it** (`process-discipline`; measured to have lied twice).

### Harness RE-EVALUATION, 2026-09-21 (branch `docs/harness-reeval`)

- **W4-P3 R5 second-reader RECOUNT** (`f5f2cf2d` vs `main` `17b31b11`) ->
  `scratch/_harness_eval/reeval/qa_recount.md`. The claim "no clause was lost" across R2
  consolidation, R3 re-tiering and R4 prose: **CONFIRMED, 323 of 323, 0 LOST** (25 of 25
  audit items in the new `retirement-audit` skill, 21 byte-identical; §6b's census
  byte-identical). 4 WEAKENED, 5 findings, 10 rejected candidates — all REPAIRED in-tree by
  the parent. L-089 / **E19**,**E20**. ⭐ It also LANDED the gloss debt: "a gloss beside a
  link asserts the target's scope" is now `instrument-doctrine` rule X3 (tell) AND skill X3
  (remedy), so that row left the debt table.

### Harness substrate campaign, 2026-09-20 (branch `docs/development-substrate`)

Five dispatches, all READ-ONLY, all reports under `scratch/_harness_eval/review/`.

- **W1-P3 `tools/harness/` package refactor** (`refactor/harness-generator`, uncommitted)
  -> `scratch/_harness_eval/k3/qa_refactor.md`. All 8 claims CONFIRMED except the
  always-on LABEL (REFUTED: 20341 omits 7479 of hand-maintained always-on rules).
  14 findings; the two to act on are an unguarded second `yaml.safe_load` that crashes
  `sphinx-build`, and a skill front matter allowed to be empty. L-087 / **B16**.
- **W1-P3 generator MUTATION** -> `qa_generator_mutation.md`. All 3 claimed arms
  have teeth; 7 findings (the DRIFT message's own repair duplicates a role block
  then blesses it; `budget_tokens` ignored on `[[agent]]`; nothing runs `--check`).
  L-082 / **A21**,**A22**.
- **W4-P3 rule-core DISTILLATION fidelity** -> `qa_fidelity_rules.md`.
  `plan-authoring` 61 of 67 full, 4 LOST; 1 unsafe merge; appendix 110 not
  reproducible. **[REMEDIED @`0d376e4f`]** — all repaired, verified below.
  L-081 / **A20**,**E9**,**E10**.
- **W4-P3 skill-core fidelity** -> `qa_fidelity_skills.md`. Mode 8 at 1 of 11 was
  the finding. **[REMEDIED @`0d376e4f`]**.
- **W4-P3 LESSONS RETIREMENT claim** -> `qa_lessons_retirement.md`. Retired set
  EXACTLY 20, 18/1/1 exact; 12 of 20 carried, 8 PARTIAL, 0 losses; 44/44 anchors
  resolve; 9 stale pointers / 7 targets. L-083 / **A23**,**A24**.
- **W4-P3 T5 RESIDUE re-review** (`df3e0f31` @ HEAD `ae381c4b`) -> `scratch/_harness_eval/review2/qa_t5_closures.md`.
  R1 (both arms MUTATION-proven), R2, R3, R4, R6 **CLOSED**; R5 closed on the pointer,
  **R5′ open** (`explorer.md:3`'s gloss names 2 of the clause's 4 requirements, staled by
  its own commit). 2 new: an evidence paraphrase dropped "in a code fence"; `rules/
  workflows.md:70`'s `[M]` "≈29.5K" (from `a5545ced`) matches no recorded instrument
  while its "≈71K" half is keep−omit exactly. Census 85/37, unchanged. L-086 / **E16**,**E17**.
- **W4-P3 K3b RULE-GENERATION fidelity** (`docs/k3b-rules`, uncommitted)
  -> `scratch/_harness_eval/k3/qa_k3b.md`. The three always-on rules brought under
  generation and `delegation` retired into the workflows rule + W7. **82 of 82 clauses
  KEPT, 0 LOST**; 23/23 links; the 4 generated rules body-identical to source. The finding
  was not textual: the Zotero-liveness clause lost its always-on AUDIENCE and was absent
  from the AGENT.md W7 names as its carrier, which held the opposite instruction.
  **[REMEDIED @ the carrier]**, with the four dropped support details restored.
  L-088 / **E18**.
- **T5 SUPPORT-BRIEF paragraphs** (W4, uncommitted) -> `scratch/_harness_eval/t5/qa_support_briefs.md`.
  12 findings; the headline REFUTED by a zero-tool explorer probe: a Support agent
  keeps its OWN memory index and skills, only the PROJECT index is dropped. After
  repair 11 closed, 1 half (the new "how a budget is set" comment is false for 4 of
  8 rows). L-085 / **B-DISPATCH**.
- **W4-P3(b) RESTORATION fidelity** (`0d376e4f`) -> `qa_restorations.md`.
  **23 of 30 faithful**, 1 lost a discriminating step (Mode 11's surrogate
  qualifier + its two-sided mutation check), **0 invented numbers** — all four
  defects MODAL. ⭐ the rewritten appendix convention IS reproducible (13 of 13
  rows, Total 110) though the actual-before was 108, so the published total is
  inert across a commit that grew the corpus. 179/179 links; `--check` clean.
  L-084 / **E11**–**E15**.

### Earlier reviews — one line each

- **#428 four-solver (n,2n)** 2026-09-03 -> `scratch/_428_four_solver_check.md`;
  all SIX families handle MT=16, ERR-023 FIXED. L-079 / **E7**,**A17**.
- **#448 SN eigenvalue finalize** 2026-09-06 -> algebra correct on all 5 arms, 6
  findings, no numerical defect. L-080 / **A18**,**A19**.
- **#429 symmetry/quotient carve** 2026-09-02 -> confirmed by brute force; 8
  findings. L-077 / **A15**,**A16**.
- **#426 (n,2n) anisotropy reproduction** 2026-09-03 -> `scratch/_426_repro.md`;
  REPRODUCED to 9 dp; the pcm gap was a UNIT and it inverts the study's
  conclusion. L-078 / **F21**–**F23**.
- **CS4c step 0 feeding census** 2026-08-30 -> L-076 / **A14**.
- **SN specialization audit** 2026-08-26 -> `scratch/specialization_audit.md`;
  findings 1-2 REMEDIED. Live residue: DD-hardcoded `2.0`s in
  `psi_half_angle_seed.py:180-185`; `cylindrical_streaming`'s docstring
  recommends two quadratures its guard refuses. L-075 / **A13**.
- **CS4a-R Phase-1 gate review** 2026-08-21 -> L-074 / **A12**.
- **CS4a round-1 design assemblies** 2026-08-20 ->
  `scratch/cs4a_attack_{algebra,physics,parsimony}.md`; all three miss
  `La13511Case` (`sood_registry/la13511.py:171`). L-072 / L-073.
- **Task 51 — 7 CYL snapshot reds** -> `scratch/task51_cyl_snapshot_audit.md`;
  RE-BASELINE all 7. Two doc repairs still stand (the FALSE `[M]` in
  `tests/sn/regression/_generate_snapshots.py`; `test_streaming_operator.py:869`).
  L-069.
- **Q5.6.4 SN cylindrical τ** -> `scratch/q64_attempt2_qa_review.md`; five
  tree-carried findings. L-068.
- L-001..L-062 SN campaigns: all merged to `main`.

⛔ 2026-09-20 ruling: a HARNESS lesson is NOT a vv-principles item (it lives on
`docs/development/harness.md`). Re-check each row below against that before landing.

### ⚠ Standing debt — `vv-principles` §Anti-patterns drop-ins NOT landed

Every brief so far forbade tracked-file edits. Land these before the next review
closes; each drop-in text is at the cited digest rule or lesson.

| owed | source |
|---|---|
| N1 inferred-relation-under-a-DECLARED-name; N3 a recall counter DOWNSTREAM of a filter | L-070 §12 / **A11** |
| Γ-reuse as a licence to fold; the α-normalised-AST control check | L-077 |
| the overloaded-unit-name rationale | L-078 |
| a `catches` marker on a `slow` test; a two-stage census needs a control per STAGE | L-079 / **E7**,**A17** |
| two mechanisms under ONE check are two clauses | **E10** |
| a retirement note is a CARRIER claim | **A23** / L-083 |
| #17 (i) *prescribed repair* + (a)'s *per CONSUMER kind* | ✅ LANDED @`0d376e4f` — but with **no evidence entry**, and (i)'s `[M]` carries no configuration (L-084) |
| a before/after `[M]` pair states ONE instrument for both halves (`plan-authoring` §4, the temporal twin of RATIO-NEEDS-ITS-POPULATIONS) | L-086 / D2 |

## 3. Durable reference (topic files)

- [field_role_typing_apply_sourcesink_contract.md](field_role_typing_apply_sourcesink_contract.md)
  — re-checkable SN role contract (`.apply`=AngularSourceSink,
  `.solve`=AngularFlux) + the A2D-1 source-hash-pin update procedure +
  affine-gate test-migration playbook. Cited by `qa/AGENT.md` enforcement
  rule #10 — **durable**.
- [phase1_moment_space_review.md](phase1_moment_space_review.md) — the
  ERR-039 moment-space P1.1–P1.7 verification-of-record; cited as the
  verification artifact by 3 plan/agent files — **durable**.
- [issue_247_legA_review.md](issue_247_legA_review.md) — full #247 Leg A
  (slope-source) review; distilled into L-037. **Retire candidate** (merged
  campaign; the reusable behavior is in the lesson).
- [issue_251_legB_review.md](issue_251_legB_review.md) — full #251 Leg B
  (boundary face-slope) review; distilled into L-038. **Retire candidate**
  (merged campaign; the reusable behavior is in the lesson).
