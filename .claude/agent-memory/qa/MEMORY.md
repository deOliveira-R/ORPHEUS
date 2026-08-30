# QA Agent Memory — index

## 1. Lessons — a HOT digest over a COLD archive

Two files, read at different times. Do **not** summarize lesson content here.

- **[lessons.md](lessons.md)** — the HOT digest (~760 lines). Behavioral rules
  only: one imperative rule + its failure→correction core + a
  `→ lessons_archive.md L-0NN` pointer. Nine sections: **A** mutation mechanics,
  **B** structural blindness, **C** structural independence, **D**
  re-baseline/bit-identity, **E** markers & audit surface, **F** claim-scope,
  **G** doc-correctness, **H** mechanics/environment, **I** the map of what is
  already in `vv-principles`/`numerical-bug-signatures` (point, don't restate).
  **Read this before every review.**
- **[lessons_archive.md](lessons_archive.md)** — the COLD archive (~4800 lines,
  L-001..L-076, append-ordered). War stories, evidence, `file:line`, measured
  tables, verdicts. **Open only the specific `L-0NN` a digest rule points at** —
  never read it whole (it is ~60K tokens).

Maintenance: a new lesson appends `L-0NN` to the archive AND lands a 2–5 line
rule in the digest. Sharpen the digest in place; **it is now past the ~400-line
distillation trigger — distill before the next append** (encode a shared
meta-lesson once, list its instances); never truncate.

## 2. Active / in-flight state

**CS4c step 0 — the RUNTIME FEEDING CENSUS, 2026-08-30.** ⏹ **DELIVERED** —
`scratch/cs4c_feeding_census.md` (HEAD `2f44ed4e`, tree clean; read-only, no
tracked file edited; instruments `scratch/cs4c_step0_{spy,drive}.py` + 2 probes).
11 production entries × the 13-site roster; 23-verb denominator; every arm fired
by an activation control; 11/11 headline numbers bit-identical instrumented vs
control. Load-bearing, all `[M]`: (1) **SN's C binding
(`sn/coupled_system.py:446`) is minted 20–25× per k-solve and applied ZERO
times** — `StreamingCollisionOperator` overrides `OperatorSum.apply` and reads
`self.diagonal.coefficient.values`; (2) the S `FullField` arm **re-dispatches**
on `psi.interior`, so "4 arms" over-counts bodies; (3) on **2-D Cartesian**, S's
arriving composite is a *different* `FullFieldSpace` (moment interior) from the
bound one (angular interior) — a shipped non-endomorphism corroborating the §1
two-space ruling; (4) `S.apply[ScalarFlux]` = **0 prod traffic** (#205 arm);
(5) `C.apply[ndarray]` = **2 calls, one site** (`homogeneous/solver.py:194`);
(6) the iso pair returns **bare ndarray from typed input**, 4816 applies, bound
`space=None`. Refuted: the windowed carrier is routed by **dimensionality**, not
by `inner_schedule` (12-row probe); and SN C is silent under Krylov too.
✅ `vv-principles` #29 sharpening (the (d) NO-arm way + bodies-not-arms) **LANDED**
— nothing owed. Lesson **L-076**, digest rule **A14**.

### ⚠ The ONE still-open debt

**Nexus V&V demand memo (L-070, 2026-08-15) — two novel rationales still OWED to
`vv-principles`, `[M]` re-verified 2026-08-30 (0 hits in the skill):** N1
*inferred-relation-under-a-DECLARED-name*, N3 *a recall counter placed DOWNSTREAM
of a filter cannot count what the filter dropped*. Drop-in text: L-070 §12 (they
live in digest **A11** only). N2 (the Mode-8 DECODER dual) ✅ landed.
**Land N1/N3 before the next review closes.**

### ⏹ Complete — one line each; the evidence is the lesson, the tense is git

- **CS4c step 0 feeding census** (2026-08-30) → above. L-076 / **A14**.
- **SN specialization audit** (2026-08-26) → `scratch/specialization_audit.md`;
  findings (1)+(2) REMEDIED by P1 (`ebe5d22f`, `37d6d1af`); the class is
  `StreamingCoefficientCache` now. Still-live residue: the DD-hardcoded `2.0`s in
  `psi_half_angle_seed.py:180-185`; `cylindrical_streaming`'s docstring
  recommending two quadratures its own guard refuses; the `Bailey 2009` citation
  its module header condemns. L-075 / **A13**.
- **CS4a-R Phase-1 gate review** (2026-08-21) → 12 findings; the `#17` sharpening
  ✅ landed. L-074 / **A12**.
- **CS4a round-1 design-assembly reviews** (2026-08-20) →
  `scratch/cs4a_attack_{algebra,physics,parsimony}.md`; all four OWED skill items
  ✅ landed (`vv#29`, `#28`'s 7-of-13 correction, `#30`, `plan-authoring` §10
  third shape). Still-live: all three assemblies miss `La13511Case`
  (`sood_registry/la13511.py:171`). L-072 / L-073.
- **Task 51 — the 7 CYL snapshot reds** (L-069) →
  `scratch/task51_cyl_snapshot_audit.md`; verdict RE-BASELINE all 7. Two doc
  repairs still stand: the FALSE `[M]` *"at M=2 τ did not change at all"*
  (`tests/sn/regression/_generate_snapshots.py`; `Δτ = 2.071e-01`) and
  `test_streaming_operator.py:869`'s CYL bit-identity claim.
- **Q5.6.4 SN cylindrical τ** (L-068) → `scratch/q64_attempt2_qa_review.md`;
  refutes link C6. Tree-carried findings: the false Hébert 3.437/3.439 citation;
  no value gate on the partition producer; **one** catcher on the cylinder τ;
  `alpha_defect_beta`/`nu_closure_residual` have no test consumer; no
  `ψ̂`-positivity gate.
- L-001..L-062 SN campaigns: all merged to `main`.

⛔ **Every branch/hash/line above is a SNAPSHOT — reconcile with git and a grep
before acting on any of it** (`process-discipline`; measured to have lied twice).

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
