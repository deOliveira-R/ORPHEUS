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

**SN streaming/redistribution SPECIALIZATION audit (report-only), 2026-08-26.**
⏹ **DELIVERED** — `scratch/specialization_audit.md` (938 lines), branch
`refactor/unweld-phase-b` @ `226cc6ca` with 3 corpus files DIRTY (the un-weld
Phase-B carve is uncommitted; the memo describes the WORKING TREE). C1–C4 table
over ~45 symbols + D1/D2/D3. Load-bearing, all `[M]`: (1) `mu_start` is a
**three-link dead chain** — `GeometryCoefficients.mu_start` has zero readers of
any kind, so `StreamingTerms.mu_start`'s only production consumer is the write
into it, while its docstring names a closure that reads the OWNER instead;
(2) `requires_upstream_angular_state`/`angular_marching_axis` — 0 production
readers, and a **997-operator flip over 2591 rows reddens exactly the 6
assertions that name them** (the founding measurement for lesson **A13**);
(3) the DD-hardcoded `2.0`s in `psi_half_angle_seed.py:180-185` sit in the
ANGULAR factor of a product `pairing.py` declares orthogonal — latent, guarded
only by `_require_slab` in another package on the other factor; (4) ⛔
`cylindrical_streaming`'s docstring recommends `Quadrature.level_symmetric`
AND `.product`, **both refused** by the guard its only caller runs 12 lines
later (`folded_product` appears only in an error string); (5) the module header
condemns the `Bailey 2009` citation its own class docstring still makes
(`:12` vs `:499`). **D1 verdict**: α is ANGULAR (`alpha_dome` takes no geometry
arg; `angular_redistribution(quad, coord)` needs no mesh and is CALLED that way
at `augmented_mesh.py:417`) — ownership already fixed by Phase B, the HOME is
not. ⭐ Reusable: **a boundary Protocol's SIZE measures a misplacement** — trace
which members survive moving the suspect object; 4 of `AngularMeasure`'s 6 go
with α. ✅ `vv-principles` #17 ⭐ sharpening **LANDED** (nothing owed).
Lesson **L-075**, digest rule **A13**.

**CS4a-R Phase-1 GATE review (report-only), 2026-08-21.** ⏹ **DELIVERED** as the
final message (no file; `feature/cs1-energy-space` @ `a9a2d55a`). Re-run: all
three gate files + D5 byte gate **green** (163 passed / 14 xfailed; D5 8/8), no
tracked file touched. Twelve findings, three load-bearing, all `[M]` by
in-process plugin mutation: (F1) the shared `assert_energy_extent_conforms` has
**4 call sites and 1 witness** — disabling C → 0 reds, IsoS+IsoN2N → 0 reds, F →
1; `"energy extent"` occurs in exactly ONE assertion tree-wide, and the
un-witnessed C site is the one passing a *different* ng expression. (F2) G1.3's
"independently assembled" licence is false — kernel and carrier run the SAME
`np.asarray(s.todense())`; transposing BOTH leaves **51/51 green** (the tree
catches it only in `tests/homogeneous/`, 17 reds). (F3) G2.7's named falsifier is
unreachable — `apply_metric(x) **is** x`, so `.H` literally calls
`apply_transpose`; a dense AFFINE transpose keeps it green (its one live
falsifier is the wrapper ceasing to delegate). Plus: the strict-xfail gate's
denominator is 3 hand-named lists (new list ⟹ silent evasion, measured); the R1
annotation gate reads `domain` only, never `codomain`; G2.1's "all 8 D5 cases" is
**1 of 8** discriminating. Two of my own attacks WITHDRAWN by probe (G2.4's
monkeypatch is live; G2.3's frozen literals reproduce independently). ⚠ **OWED,
NOT landed** (report-only charter): one `vv-principles` **#17 sharpening** — *a
hoisted guard's arms are its CALL SITES* — drop-in text in the memo §OWED.
Lesson **L-074**, digest rule **A12**.

**CS4a round-1 design-assembly reviews (algebra / physics / parsimony),
2026-08-20.** ⏹ ALL DELIVERED — memos `scratch/cs4a_attack_{algebra,physics,
parsimony}.md` @ `71515847`; full evidence in **L-072/L-073**. The three
still-live carry-forwards: (i) `[M]` per-INSTANCE runtime census — **6 of 12**
production instances refute construction-time arm selection from the bound space
(`SNSolver.fission_op` is bound to `full_field_space` and fed `ndarray` ×17
only), which kills the mechanism all three assemblies shared; (ii) ⛔
`vv-principles` **#28's `[M]` "8 of 13" is 7 of 13** (runtime-probed —
`MaterialMesh` has no `full_field_space`); (iii) all three miss `La13511Case`
(`sood_registry/la13511.py:171`), the medium's existing third spelling and its
only real first consumer. ⚠ **OWED to the skills, NOT landed** (charters forbade
edits): `vv-principles` **#29** (text in algebra memo §5.2), #28's corrected
number, a codomain-constructor discriminator, and a `plan-authoring` §10
third-shape (designed-red tell) — text in the physics memo §D. Landed at the
time: `vv-principles` #28, lessons **H15**, and the ⛔ correction of **H14**.

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
