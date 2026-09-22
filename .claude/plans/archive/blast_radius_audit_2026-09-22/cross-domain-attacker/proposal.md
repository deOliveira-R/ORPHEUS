# Blast-radius proposal — cross-domain-attacker, 1 candidate

Dispatch 2026-09-21. Read-only on tracked files; this is the only file written.

## Verdict

| file | verdict | reason (one line) | referrer edits |
| --- | --- | --- | --- |
| `.claude/agent-memory/cross-domain-attacker/elegance_smell_rank_non_monotone.md` | **SALVAGE, then RETIRE** | Every behavioural clause (smell, why, how-to-apply, trigger, ORPHEUS precedent) is in `cross-domain-frames` reference.md Part C Smell #15 lines 252-266 — except the file's "Diagnostic first test", which is **not** in the skill and was **falsified** as RH13 on 2026-04-22, so the file states a wrong instruction in the present tense. | 5 (R1-R5 below) |

The file's own survival condition — "remains as evidence / precedent until the
live rank-N work threads close" (line 53-54) — is met: its research log sits in
`.claude/plans/archive/`, and the only non-archive plan matching `rank-?n` is
`.claude/plans/dsa_saddle_point_frame.md:1238`, which is the MPO refutation and
does not concern this file. `[M]` `Grep(pattern="rank-N|rank_n|rank-n", path=".claude/plans", -i)`
→ 19 files, 18 under `archive/`.

### Why SALVAGE and not plain RETIRE

The file's lines 34-41 prescribe: compute `R[ψ_b]` of the converged rank-1
ansatz; if it matches the production eigenvalue to 1e-5 the method IS
variational, then build rank-2 as a Gram-Schmidt-nested subspace. That test was
run and **rejected**:

> `.claude/plans/archive/rank-n-closure-research-log.md:1867-1898` — "RH13 —
> rejected … `R[1] ≈ 2×10⁻⁴` at the anchor — four orders of magnitude below
> `k_eff^{F.4} = 1.4963` … **the volume term is load-bearing and absent from the
> stated Rayleigh quotient** … F.4's integral CP is a Schur-reduced eigenproblem
> … No self-adjoint Rayleigh quotient on the **boundary trace alone**
> reproduces `k_eff`. Issue #126 step 2 (nested rank-2 Ritz) NOT dispatched."

Deleting the file without salvage loses the correction and leaves the skill
entry with no first test at all. The SMELL itself survives the falsification and
was strengthened, not weakened: same log, lines 1900-1906 and 1922-1927 —
"rank-N non-monotonicity … remains valid … doubly confirmed by the
anisotropic-BC scan — rank-N non-monotonicity survives on BOTH kernels."

## Salvaged

### S1 → `lessons.md` Part 4 (Refuted-frame ledger), as a new bullet

Destination chosen because Part 4 exists to carry a frame's NAME plus the
trigger it actually needs, so an UNEXPLORED line can be written without
re-deriving the refutation; `[M]` `Read(lessons.md, 367-396)` — no
Rayleigh-Ritz entry there today.

```markdown
- **Rayleigh–Ritz on a REDUCED variable** — the frame is sound and its smell is
  promoted (Part C #15), but the OBVIOUS reformulation target is the wrong one.
  Refuted FOR "is the F.4 white-BC closure secretly rank-1 Ritz on the boundary
  trace?" (RH13, 2026-04-22): `R[1] ≈ 2e-4` against `k_eff = 1.4963`, because
  the CP eigenproblem is SCHUR-REDUCED and the eliminated volume block is
  load-bearing. The FACT it establishes, and the standing rule: **a Schur
  complement is not self-adjoint in the reduced variable's inner product, so no
  Rayleigh quotient exists on the reduced variable alone** — before claiming a
  variational principle, ask which block was eliminated. The frame needs the
  FULL (volume + trace) eigenproblem with a nested ladder
  `V_n = V_{n-1} ⊕ µ·V_{n-1}`, Galerkin on the whole operator; the smell itself
  survived and was doubly confirmed on both kernels by the anisotropic-BC scan.
  Literature: Courant & Hilbert 1953 Vol. I §VI (min-max); Case & Zweifel 1967
  §6 and Wendroff 1961 (Boltzmann variational theory).
```

Nothing else in the file is unique. Its two "why"/"how to apply" paragraphs are
reference.md 252-266 verbatim in substance; `Case-Zweifel 1967 §6, Wendroff 1961`
already survive at research-log line 1918; the PCA-sector M = 2/3/5 numbers are
a spent VALUE claim whose sighting the skill already names ("Direction-C/Q
failures").

### S2 → OPTIONAL uplift, `cross-domain-frames` reference.md Part C Smell #15

Not required for retirement; offered because Smell #15 is the only Part C entry
with a promoted status and no discriminator, and my role contract requires a
fail-able first test. One clause, appended after line 266:

> First test (corrected; the boundary-trace-only version was falsified as RH13,
> 2026-04-22): form the Rayleigh quotient on the FULL eigenproblem, never on a
> Schur-reduced variable, and build rank-(N+1) as a NESTED subspace
> (`V_n = V_{n-1} ⊕ µ·V_{n-1}`); a Ritz method then cannot be non-monotone, so a
> surviving non-monotonicity refutes the variational reading rather than the
> truncation order.

If the orchestrator takes S2, the Version History block (reference.md 307-315)
gains its dated line.

## Referrer edits

`referrers.md` lists 4. `[M]` re-run: `Grep(pattern="elegance_smell_rank_non_monotone", path=<repo>)`
→ 11 hits over 8 files; excluding this audit's own directory and
`scratch/_memory_distillation/` (the pre-apply staging copies of R4/R5), the
consumer set is **5**, not 4 — see the finding on R5.

| # | referrer | class | exact edit |
| --- | --- | --- | --- |
| R1 | `.claude/plans/archive/rank-n-closure-research-log.md:1592` | repo (archived session artifact list) | Keep as history, append the forwarding pointer. Replace the line with:<br>``- `.claude/agent-memory/cross-domain-attacker/elegance_smell_rank_non_monotone.md` (retired 2026-09-21; content in the `cross-domain-frames` skill, reference.md Part C Smell #15)`` |
| R2 | `.claude/plans/archive/rank-n-closure-research-log.md:1757` | repo (present-tense DEFINITION pointer) | Re-point to the skill. Replace line 1757:<br>``(`.claude/agent-memory/cross-domain-attacker/elegance_smell_rank_non_monotone.md`):``<br>with<br>``(`cross-domain-frames` skill, reference.md Part C Smell #15; promoted 2026-04-30 from a memo retired 2026-09-21):`` |
| R3 | `.claude/plans/archive/memory_distillation_2026-09-21/cross-domain-attacker/table.md:77-79` and `:123` | repo (the distillation's own candidate list — not a consumer, per `referrers.md`) | Both are dated records and stay. Per `plan-authoring` REMEDIED-FACT, close the trail: append to line 79 (after "which this dispatch cannot run.")<br>`[REMEDIED 2026-09-21] the audit ran; SALVAGE-then-RETIRE, salvage in lessons.md Part 4.`<br>and append to the clause at line 123 ("…is not retired (a graph node; the blast-radius audit is a separate pass)")<br>`[REMEDIED 2026-09-21] retired; and the "graph node" premise was wrong — `.claude/` is outside the Nexus graph.` |
| R4 | `.claude/agent-memory/cross-domain-attacker/issue_168_phase_c_sweep_frame.md:77` | own-memory (live cross-link) | Re-point. Replace<br>`- elegance_smell_rank_non_monotone.md (Smell 15) — Smell 16 is a sibling.`<br>with<br>``- `cross-domain-frames` Part C Smell #15 (rank-N non-monotone) — Smell 16 is a sibling.`` |
| R5 | `.claude/scratch/open_fronts_audit.md:189` | repo — **NOT in `referrers.md`** | Close the open box. Replace<br>`- [ ] **J9. cross-domain-attacker `elegance_smell_rank_non_monotone.md`** — dangling memo; subject not actioned this session.`<br>with<br>`- [x] **J9. rank-N non-monotone smell** — CLOSED 2026-09-21: promoted 2026-04-30 to `cross-domain-frames` Part C Smell #15 and shipped as a gate (`tests/cp/test_peierls_rank_n_protocol.py::test_protocol_catches_non_monotone`); the memo retired, its falsified Ritz test salvaged to `lessons.md` Part 4.`<br>Apply only if the file is tracked — I could not check (no Bash). If untracked, no edit. |
| R6 | `.claude/agent-memory/cross-domain-attacker/MEMORY.md:93-95` | own-index | Delete the clause. Replace the three lines with the single line:<br>`Files here not listed above = archive: early Peierls / Variant-α / phase attacks.` |

`scratch/_memory_distillation/cross-domain-attacker/{MEMORY.proposed.md:94, table.md:77,123}`
are untracked pre-apply staging copies of R6/R3 and need no edit.

## Findings (discrepancies against the brief's claims)

1. **The census missed a referrer inside `.claude/`.** `referrers.md`'s predicate
   excludes `scratch/`, which also excludes `.claude/scratch/` — where R5 lives.
   R5 is the one referrer that is an OPEN tracker item (`- [ ]`) naming this file
   by stem, i.e. the most consequential of the five. `[M]`
   `Grep("elegance_smell_rank_non_monotone", path=<repo>)` → `.claude/scratch/open_fronts_audit.md:189`.
   Recommend the other owners' censuses be re-read with `.claude/scratch/`
   included.
2. **The stated retirement blocker does not exist.** The apply-time note at
   `table.md:123` withheld retirement because the file is "a graph node". It is
   not: `[M]` `mcp__nexus__file_brief(".claude/agent-memory/cross-domain-attacker/elegance_smell_rank_non_monotone.md")`
   → `"error": "… is not in the graph"`; `query("elegance_smell_rank_non_monotone")`
   → `[]`; `query("agent-memory")` → `[]`, i.e. `.claude/` is an excluded tree.
   No graph edge is broken by this retirement, for any owner.
3. **The smell is load-bearing in the tree, the memo is not.** Smell #15 ships as
   an executable acceptance protocol: `assert_rank_n_structural_win` in
   `tests/cp/test_peierls_rank_n_protocol.py`, with three L1 gates
   (`test_protocol_catches_non_monotone:322`, `…_sign_flip:295`,
   `…_non_beat_any_quad:348`) all `@pytest.mark.verifies("peierls-rank-n-stability")`.
   None references the memo by stem, so this is a consumer of the CONTENT, not a
   referrer of the FILE — it argues for keeping the skill entry, not the file.
4. **The file carries a present-tense falsehood today.** Lines 34-41 prescribe a
   test refuted 17 months of campaign-time ago, with no banner
   (`plan-authoring` §3, READ-THE-END-FIRST). Retiring it resolves this; if the
   orchestrator overrides to KEEP, the banner is mandatory.

## NEEDS:

- **No Bash in this dispatch** (the tool set was Read/Grep/Glob/Write/Nexus).
  So: (a) the brief's requirement to verify status lines against `git`/`gh` could
  not be met — `#121`, `#122` (cited closed in reference.md:315), `#125`, `#126`,
  `#127` are unverified here; the retirement does not rest on any of them, since
  the survival condition was checked against the archive location of the plans
  instead. (b) I could not determine whether `.claude/scratch/open_fronts_audit.md`
  is tracked, which gates R5. (c) Line counts are unmeasured.
- Verification owed by the applier: `gh issue view 126 --json state` (expect
  CLOSED — the log recommends it for closure at line 1897-1898) before the R2
  wording ships, and `git ls-files --error-unmatch .claude/scratch/open_fronts_audit.md`
  for R5.
