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

**CS4a round-1 adversarial review — the ALGEBRA assembly, 2026-08-20.**
⏹ **DELIVERED**: `scratch/cs4a_attack_algebra.md` (`feature/cs1-energy-space` @
`71515847`). (R1) **the bound space does not determine the carrier** — `[M]`
per-INSTANCE runtime census over one solve per family: **6 of 12** production
instances refute construction-time arm selection; sharpest, `SNSolver.fission_op`
is bound to `full_field_space` and fed **`ndarray` ×17 and nothing else** on the
k-eigenvalue path. This kills the mechanism **all three** assemblies share, and
the algebra assembly had listed the census as its own strongest self-attack and
deferred it. (R2) the flip's true §6b unit is **10 production + 165 space-less
TEST constructions in ~50 files**, and the one new *numerical* gate it ships (the
counting-measure adjoint theorem) reads `0.0` under the defect, `0.0` under the
fix and `≤4.4e-16` at `V_cell` spread 3358× — Mode-12 dud. Per-objective: O2/O7/O9
BREAK, O1/O6/O8 HOLD, O3/O4/O5 UNDERSPEC. Duty-B adjudications: **parsimony's a1
breaks the tree** (deletes `MaterialXSField.mesh`, its 3 readers are scheduled
for a5 ⟹ its own "D5 8/8" gate is un-runnable); its `scattering.py:713` blocker
is REAL as a fact, **self-inflicted** by that same ordering, and its stated C8
cause is wrong; **physics has the safest ordering** (K2 *adds* a constructor arm,
so K1–K3 cannot break a call site). ⚠ **Owed, NOT landed** (charter forbade
edits): `vv-principles` **#29** (drop-in text in memo §5.2) + #28's number.
Lesson: **L-073**.

**CS4a round-1 adversarial review — the PHYSICS assembly, 2026-08-20.**
⏹ **DELIVERED**: `scratch/cs4a_attack_physics.md` (same branch/HEAD). (R1) its
CENTRAL claim — *spatial locality ⟹ "a mesh is never data of the interaction"* —
has a **true premise and a false corollary**: `[M]` all four in-scope operators
ARE spatially diagonal and all still read `mesh` off the OPERAND to build the
codomain carrier (≥11 sites; `FullFieldSpace` carries no mesh on any block), and
`streaming.py:589` asserts the opposite in a production guard. (R2) O7's tell —
*"grep `SigS` finds the datum owned once, viewed once"* — is **DESIGNED-RED**:
`[M]` 70 hits, direct readers in cp/moc/mc/sn+derivations, all inside its own
Untouched list. ⛔ **Corrects the sibling-landed `vv-principles` #28: its `[M]`
"8 of 13" is 7 of 13** (`from_mesh`'s chain NAMES `full_field_space` but
`MaterialMesh` has none — runtime-probed). Duty-B ledger, rule applied per
channel: scattering **mint**, **n2n mint** (parsimony REFUTED — `N2NMomentOperator`
is a second realization, and its fold silently overturns a twice-stated physics
ruling), fission **property** (physics/algebra wrong), collision **property** —
each assembly 3-of-4, on different rows. Funk–Hecke/`Σ_s0ᵀ`: **no divergence**
(the ᵀ is the storage convention, shared). ⚠ **Owed to the skills, NOT landed**
(brief forbade edits): #28's number, a codomain-constructor discriminator, and a
`plan-authoring` §10 third shape (designed-red tell) — text in the memo's §D.

**CS4a round-1 adversarial review — the PARSIMONY assembly, 2026-08-20.**
⏹ **DELIVERED**: `scratch/cs4a_attack_parsimony.md` (branch
`feature/cs1-energy-space` @ `71515847`). Two refutations, both `[M]` mine:
(R1) **`BulkField.mesh` is MANDATORY on all 10 concrete field leaves** — probed:
`CrossSectionField(values, space)` → `TypeError: missing … 'mesh'`, `mesh=None`
→ `AttributeError: 'NoneType' … 'ng'`, no `from_space` constructor — so
retiring the degenerate carrier while listing `CrossSectionField` *Untouched*
leaves the homogeneous `MultiplicationOperator` **no constructible
coefficient**; (R2) the "**529 files**" rename blocker is **26** (503 in
gitignored `docs/_build/`). Cross-cutting: **no** assembly discharges O4 (`[M]`
**28** mesh-bound output mints inside apply bodies — the operand's mesh is a
*constructor argument of the result*, not dispatch), O5 (tell unfalsifiable
today), or O6 (`[M]` production reach 0). ⭐ **All three miss `La13511Case`**
(`sood_registry/la13511.py:171`) — a shipped `geometry_kind`-tagged
(spatial-description | infinite) × materials pairing with `to_geometry()`, i.e.
the medium's existing third spelling AND its only real first consumer
(`builders.py:34-53`). Landed `vv-principles` **anti-pattern #28** +
lessons **H15**, and ⛔ **corrected H14** (the `grep` wrapper does NOT honour
`.gitignore`). Re-verify against git before acting.

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
