# QA Agent Memory — index

## 1. Lessons — a HOT digest over a COLD archive

Two files, read at different times. Do **not** summarize lesson content here.

- **[lessons.md](lessons.md)** — the HOT digest (~665 lines). Behavioral rules
  only: one imperative rule + its failure→correction core + a
  `→ lessons_archive.md L-0NN` pointer. Nine sections: **A** mutation mechanics,
  **B** structural blindness, **C** structural independence, **D**
  re-baseline/bit-identity, **E** markers & audit surface, **F** claim-scope,
  **G** doc-correctness, **H** mechanics/environment, **I** the map of what is
  already in `vv-principles`/`numerical-bug-signatures` (point, don't restate).
  **Read this before every review.**
- **[lessons_archive.md](lessons_archive.md)** — the COLD archive (~4000 lines,
  L-001..L-071, append-ordered). War stories, evidence, `file:line`, measured
  tables, verdicts. **Open only the specific `L-0NN` a digest rule points at** —
  never read it whole (it is ~60K tokens).

Maintenance: a new lesson appends `L-0NN` to the archive AND lands a 2–5 line
rule in the digest. Sharpen the digest in place; **it is now past the ~400-line
distillation trigger — distill before the next append** (encode a shared
meta-lesson once, list its instances); never truncate.

## 2. Active / in-flight state

**Nexus V&V surface — the demand memo (L-070), 2026-08-15.** ⏹ **DESIGN MEMO
DELIVERED** (scratchpad only; the parent agent holds the text). Headline, all
`[M]` at `a1c90aac` / nexus 0.16.1: the graph has **no test→code edge**;
`implements` is 100 % inferred (81 % on ONE token); `verified` has **no
evidence floor** (351/692 have no declared test); static `calls` had **0 %
recall** vs runtime's 7/21 on the one ground-truthable equation. Ranked demand
D1–D5 (D1 = a per-test `exercises` runtime family — I prototyped it, ~15 lines,
23/23 contexts joined). ⚠ **Three novel rationales are OWED to `vv-principles`
and were NOT landed** (write fence + the skill files were already dirty):
N1 inferred-relation-under-a-declared-name, N2 the Mode-8 DECODER dual,
N3 recall-counter-downstream-of-a-filter — drop-in text is in the memo §14 and
summarized in L-070 §12. **Land them before the next review closes.**

**Task 51 — the 7 CYL snapshot reds (L-069).** ⏹ **AUDIT COMPLETE**;
`scratch/task51_cyl_snapshot_audit.md`. Verdict **RE-BASELINE all 7** — they
are the REMAINDER of `39b46a31`'s re-baseline (already in the tree, scoped to
one directory). Sole cause bisected to `3dda18ca` (ω-partition + absorber
retirement); `c33178ef` and the fold are both REFUTED as causes. Two blocking
doc repairs stand regardless: the `[M]`-marked *"at M=2 τ did not change at
all"* in `tests/sn/regression/_generate_snapshots.py` is FALSE (`Δτ =
2.071e-01`), and `test_streaming_operator.py:869` claims the CYL arms are
still bit-identical. Re-verify against git before acting.

**Q5.6.4 — the SN cylindrical τ partition (L-068).** ⚠ Live on branch
`refactor/operator-strategy-layers` (NOT main; the harness's session-start git
snapshot said `main` and was wrong). My review
`scratch/q64_attempt2_qa_review.md` **refutes link C6** of the proposal and
fills the memo's own PENDING decisive row (`τ≡½ = 1.0181e-01`, 1.54× worse than
the retired chord+absorber). Findings the tree still carries regardless of how
the design lands: the production Hébert 3.437/3.439 citation is false three
ways; the single-source partition producer has **no value gate**; the cylinder
τ has **exactly one** catcher; `alpha_defect_beta` / `nu_closure_residual` have
**no test consumer at all**; no `ψ̂`-positivity gate on either arm. Re-verify
against git before acting.

All SN review campaigns behind L-001..L-062 are merged to `main` (verified via
`git merge-base --is-ancestor`). Git is authoritative for merge status —
re-verify before acting on any "in-flight" claim.

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
