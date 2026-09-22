# QA Agent Memory — index

## 1. Lessons — a HOT digest over a COLD archive

Two files, read at different times. Do **not** summarize lesson content here.

- **[lessons.md](lessons.md)** — the HOT digest (**726 lines** at 2026-09-21, from 1076).
  Behavioral rules only: one imperative, the check that makes it decidable, and a
  `→ L-0NN` pointer. Ten sections: **A** mutation mechanics, **B** structural
  blindness, **C** structural independence, **D** re-baseline/bit-identity,
  **E** markers & the audit surface, **E′** auditing a rewrite or distillation,
  **F** claim-scope, **G** doc-correctness, **H** mechanics/environment, **I**
  what retired INTO the skills and rules. **Read before every review.**
- **[lessons_archive.md](lessons_archive.md)** — the COLD archive (L-001..L-090,
  append-ordered): war stories, evidence, `file:line`, measured tables. **Open
  only the `L-0NN` a digest rule points at** — never read it whole.

Maintenance: a new lesson appends `L-0NN` to the archive AND lands a 2–5 line
rule in the digest. Before appending, apply the per-line test to a neighbouring
rule: a rule that restates a clause of a preloaded skill or an always-on rule
retires into §I, citing the clause. Sharpen in place; never truncate.

## 2. Review record

Every campaign below is **merged** (`git merge-base --is-ancestor <hash> HEAD`,
2026-09-21). The verdict lives in its report, the behaviour in its digest rule,
the tense in git. Reports under `scratch/` are untracked and may be gone; the
lesson is the durable half. ⛔ Re-verify any path or symbol here before acting on
it (`process-discipline`).

| review | lesson → digest rule |
|---|---|
| harness re-evaluation, R5 second-reader recount (2026-09-21) | L-089 → E19, B-DISPATCH |
| harness substrate: `tools/harness/` package refactor | L-087 → B16 |
| harness substrate: generator mutation battery | L-082 → A21 (retired to #17), A22 |
| harness substrate: rule-core and skill-core distillation fidelity | L-081 → A20, E9, E10 |
| harness substrate: lessons-retirement claim | L-083 → A23, A24 |
| harness substrate: T5 residue re-review | L-086 → E16, A3 |
| harness substrate: K3b rule-generation fidelity | L-088 → E18 |
| harness substrate: T5 Support-brief paragraphs | L-085 → B-DISPATCH |
| harness substrate: restoration fidelity | L-084 → E11+E12, E13, E14, E9+E15 |
| #428 four-solver (n,2n) | L-079 → §I (#36, #17g) |
| #448 SN eigenvalue finalize | L-080 → A18, §I (#13) |
| #429 symmetry/quotient carve | L-077 → A15, §I (#34) |
| #426 (n,2n) anisotropy reproduction | L-078 → F23, §I (#35) |
| CS4c step-0 feeding census | L-076 → A14-r |
| SN specialization audit | L-075 → A13-r |
| CS4a-R Phase-1 gate review; CS4a design assemblies | L-072, L-073, L-074 → F16, F17, F20, §I |
| Task 51 CYL snapshot reds | L-069 → A10, D14-r, H13 |
| Q5.6.4 SN cylindrical τ | L-068 → C1, F14, H12 |
| L-001..L-067, the SN campaigns | the digest's §A–§H |

## 3. Standing debt — drop-ins NOT yet landed in a rule

None. The five rows this table carried on 2026-09-21 landed the same day: the
inferred-relation and downstream-counter clauses in the `instrument-doctrine`
skill X1; two-mechanisms-one-check in its X3; the retirement-note carrier claim
in `retirement-audit` item 19; BEFORE-AND-AFTER-ONE-INSTRUMENT in `plan-authoring`
§4; the `.claude/` sweep in `retirement-audit` item 4. A new row is a candidate
only after reading the current rule file (they change), and a HARNESS lesson's
home is `docs/development/harness.md`, never `vv-principles` (ruled 2026-09-20).

## 4. Durable reference (topic files)

- [field_role_typing_apply_sourcesink_contract.md](field_role_typing_apply_sourcesink_contract.md)
  — the SN role contract (`.apply` = AngularSourceSink, `.solve` = AngularFlux),
  the A2D-1 source-hash-pin update procedure, the affine-gate migration
  playbook. Cited by `qa/AGENT.md` #10 — **durable**.
- [phase1_moment_space_review.md](phase1_moment_space_review.md) — the ERR-039
  moment-space verification-of-record; cited from three files outside this
  memory — **durable**.
