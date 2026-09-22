# Topic-file blast-radius audit — owner `qa`, proposal

**Dispatch:** the harness campaign's close-out, owed item (1). READ-ONLY on
tracked files; this file is the only deliverable. Both candidates read in
FULL (114 and 171 lines). Every status line below was re-verified against
git or GitHub on 2026-09-21, never against a memory's claim.

## Verdicts

| file | verdict | reason (one line) | referrer edits |
|---|---|---|---|
| `.claude/agent-memory/qa/issue_247_legA_review.md` | **RETIRE** | Campaign narrative of a MERGED review (#247 CLOSED, `d9396a20` is an ancestor of HEAD, branch `refactor/sn-foundation-cleanup` deleted): every behavioural lesson is already `lessons_archive.md` **L-037** (six-step Mode-10 closeout recipe + the WATCH doc-nits), and both of its non-blocking follow-ups are REMEDIED in the tree. | (1) delete `MEMORY.md:72–76` (the joint entry, both files); (2) `issue_251_legB_review.md:11` `[[issue-247-legA-review]]` — its referrer is itself a RETIRE, so no edit. |
| `.claude/agent-memory/qa/issue_251_legB_review.md` | **RETIRE** | Same shape (#251 CLOSED): behaviour is `lessons_archive.md` **L-038** (five numbered lessons + the scope-note verbatim); its one novel finding LANDED in `vv-principles` Mode 10; its one standing follow-up is OPEN issue **#252**, whose body carries the storage-verified / sign-unverified disposition in more detail than the memo. | (1) delete `MEMORY.md:72–76` (same entry — one edit retires both); (2) no other inbound referrer. |

Nothing to salvage: no line of either file is behaviour that is not already in
the archive, a rule, a skill, an open issue, or the production tree (evidence
below). No KEEP.

## The exact `MEMORY.md` edit

Delete these five lines (72–76) from `.claude/agent-memory/qa/MEMORY.md`, under
`## 4. Durable reference (topic files)`, leaving the two `**durable**` entries
above them untouched:

```
- [issue_247_legA_review.md](issue_247_legA_review.md) and
  [issue_251_legB_review.md](issue_251_legB_review.md) — full slope-source and
  face-slope reviews of merged campaigns, distilled into L-037/L-038; `git grep`
  finds **0** references outside this memory. **Archaeology — retire candidates
  for the orchestrator's blast-radius audit.**
```

No other line in `MEMORY.md` mentions either file (`[M]` `grep -n` on both
stems: 2 hits, both inside this block). The digest's two pointers to the
LESSONS (`lessons.md:63` → L-038, `lessons.md:479` → L-037) point at the
ARCHIVE, which stays; they need no edit.

## Evidence that nothing is brought forward

Status lines, each re-verified (`[M]` 2026-09-21):

| memory claim | instrument | reading |
|---|---|---|
| "#247 … UNCOMMITTED, branch `refactor/sn-foundation-cleanup`" | `git merge-base --is-ancestor d9396a2 HEAD`; `git branch -a --list` | merged (`d9396a20`, 2026-06-18); branch deleted. The memo's front-matter `UNCOMMITTED` is a frozen-forward lie (`process-discipline`). |
| "#247 / #251" | `gh issue view --json state` | both **CLOSED**. |
| "#252 is filed OPEN with correct labels" | `gh issue view 252` | **OPEN**, `module:sn, module:tests, level:L1, type:improvement`. Body carries the physics expectation, the three closing requirements and the Mode-1 rationale. |
| Leg A nit 1, "lift docstring overstates the eigenvalue sharing" | read `orpheus/sn/solver.py:3190–3218` | **REMEDIED** — the docstring now reads "for the fixed-source path (#240 D5b-S3 / #247)… kept as a single-source helper so a *future* eigenvalue external-source hook reuses the same lift/thread policy." |
| Leg A nit 2, "`linear_discontinuous.py:304` slope-label transposition" | `orpheus/transport/spatial/linear_discontinuous.py:276` | **REMEDIED** — reads `{ψ̄, ψ̂_y, ψ̂_x, ψ̂_xy}`, the canonical Kronecker order (file also re-homed `sn/` → `transport/spatial/`). |
| Leg B scope-note, "DD identity at `n==1` ≠ spec D6" | `orpheus/sn/loss_representation/__init__.py:755–790` | still true AND now self-documented at the arm ("single-moment closure (DD/Step, `n == 1`): identity"); also verbatim in L-038's closing note. |
| Leg B self-improvement, "recommend the one-line Mode-10 row addition" | `docs/development/skills/vv-principles.md:233` | **LANDED** — "When no isolating regime exists the verification is a STRUCTURAL pair… There is then NO value-improvement leg — do not manufacture one." |
| Leg A's last distinctive datum, the per-slot consumption signals | `tests/sn/verification/mms/test_mms_ld_2d.py:867,881` | in the tree at the gate itself: `_CONSUMPTION_TOL = 1e-8`, "probed ~3e-3 x-slope, ~1e-2 y-slope, ~6e-5 xy". |
| Leg B's L-007 tagging observation | digest `lessons.md:371` (**E1**) | covered as a general rule, citing L-007. |
| Leg B item "confirm the follow-up issue exists before crediting 'filed as #NNN'" | `process-discipline` (always-on) + digest **G1** | covered by the rule floor. |

What remains in each memo after this accounting is a per-claim narration of a
merged review plus ephemeral gate readings (`460 passed / 1 skipped / 4 xfailed`,
slot shapes `(24,2,6,2)`, mutation transcripts) — the definition of archaeology.

## Findings outside the verdicts (for the orchestrator, not for me to fix)

1. **The census predicate is blind to the link spelling the memory format
   prescribes.** `referrers.md` keyed on the file STEM (`issue_247_legA_review`,
   underscores); memory files cross-link by their front-matter `name:` SLUG in
   `[[...]]` (hyphens). Re-running the census on the slug spelling found one
   referrer the stem census missed: `issue_251_legB_review.md:11`,
   `[[issue-247-legA-review]]`. Harmless here (both candidates retire), but a
   stem-only census under-reports every intra-memory `[[slug]]` link across the
   other five owners' audits. `[M]` 2026-09-21,
   `grep -rn -- "issue-247-legA-review" . --exclude-dir=.git/_build/.venv`: 2
   hits, one the definition, one the link. The rest of `referrers.md` reproduced
   exactly (own-index, 1 line each) and its positive control re-ran clean
   (`phase_f_step2_mesh_refinement`, 6 non-scratch files).
2. **OPEN issue #252 carries a dead path reference.** Its body cites
   `orpheus/sn/boundary_operator.py:219-227`; that file does not exist —
   `_reflect_trace` now lives in `orpheus/sn/operators/boundary.py` (and
   `orpheus/sn/boundary/realizer.py`). Cardinal Rule 4: an open issue is the
   plan, and this one sends its executor to a missing file.
3. **#252 also cites three OTHER owners' memory files by path**
   (`method-implementer/issue_251_legB_boundary_trace_closeout.md`,
   `explorer/issue_251_trace_space_widening.md`,
   `test-architect/issue_251_legB_boundary_gate_spec.md` — all three exist
   today). Those owners' blast radii are therefore **not** clean the way mine
   is: retiring any of them breaks a live issue. Worth passing to those audits.
4. **A live X3 nit found while verifying nit 1** (new, not the memo's): the
   repaired `_lift_external_source_to_moments` docstring asserts "One production
   caller (`_build_fixed_source_rhs`)", but the symbol now has **two**
   production call sites — `orpheus/sn/solver.py:3157` in `_build_fixed_source_rhs`
   and `orpheus/sn/solver.py:2939` in `solve_sn_adjoint_fixed_source`. Prose
   asserting a count the tree contradicts; a one-word fix.

NEEDS:
- Nothing blocking. Nexus MCP tools were available to this dispatch but were
  not needed: the question is a text-reference census over `.claude/` and
  `docs/`, which `grep` answers with a stated predicate, and agent-memory
  inbound edges add nothing over the two-spelling census run above.
- Items 2 and 4 are tracked nowhere; they need the orchestrator's decision
  (in-session fix vs issue) — I am read-only on tracked files and on GitHub.
