---
name: issue-236-phase2-stepc-tau-retirement
description: #236 Ph2 Step C geometry-τ retirement — PASS-WITH-NITS; the migrate-then-delete carve's stale-doc blast radius lands in untouched-but-adjacent hot files (present-tense claims about deleted data contracts).
metadata:
  type: project
---

# #236 Ph2 STEP C — geometry-τ retirement (the terminal τ-carve step)

**Verdict: PASS-WITH-NITS** (uncommitted, branch `feature/sn-spatial-angular-product`, MAIN checkout). ONE production file changed (`reduced_operator.py`, +56/−113); the other 16 diff files are 11 tests + 5 memory/skill. All nits are stale-DOC / dangling-reference (anti-#11), zero code rework.

## What landed (the surgical excision — exemplary)
Deletes the DEAD geometry-side τ producer + 5 fields (`StreamingTerms.tau_mm`/`alpha_in`/`alpha_out`; `ReducedStreamingOperator.tau_mm`/`tau_mm_per_level`) after B3 made the LIVE sweep/scan/matvec read closure-owned τ/c off `CellVisit`. Migrate-then-delete: test oracles re-pointed onto the structurally-independent `contamination.morel_montry_weights` (a DIFFERENT BMC-2010-Eq.43 path, imports nothing from `pole_angular_closure.py` — L11 VERIFIED) BEFORE deletion. DriftWarning-escalated gates 0 failures + count delta reconciles EXACTLY (762→758 = 4 retired Leg-1 legs; 554→552 = 2 retired bit-id tests) = no live numeric changed, no silent test loss.

## The 5 ASSESS questions — all PASS
- **Q1 surgical**: only τ removed; `mu_edge`/`eta`/`M`/`eta_edge` (τ-only locals) correctly dropped with it; α-dome recursion + GL-antisymmetry assert + `redist_dAw`/`face_areas`/`delta_A` survive, read like the domain. Smoke-verified fields gone + geometry intact.
- **Q2 `mu_start=-1.0`**: NOT magic — Hébert §3.9.4 inward pole μ₁/₂=−1, NAMED + field-docstringed + site-commented; old `float(mu_edge[0])` was *itself* −1.0 (array seeded `mu_edge[0]=-1.0`), so the swap removed a τ-serving indirection, buried no derivation. Cyl analog `-sin_theta` is the genuine per-level derived expr, kept computed.
- **Q4 `_c_surrogate.py`**: clean Pattern-5 split — `c_from_constants(τ,α_in,α_out)` pure formula primitive (reads like math, matches `balance.py:265` algebra-of-record) + `mm_constants_for_ordinate(op,...)` per-ordinate resolver. L11-independent, cyl-clamp `np.clip(...,0.5,1.0)` in right place (mirrors closure `morel_montry_tau_per_level` `:883`), `del cell_idx` honest dead-weight. Mutation-verified non-tautological (τ×1.1 RED, c-swap RED).
- **Q5 no twin**: all 7 migrated test files route through the ONE surrogate API; L12 grep found 5 readers beyond the audit list, all migrated consistently.
- **Q6 retirement complete + T5**: τ genuinely GONE (no half-retired field-then-drop). T5 (legacy `__call__`-arg `tau_mm` at `pole_angular_closure.py:320/577/1468/1777`) is a DELIBERATE documented separate-surface deferral (unbound `MorelMontryAngularSweep(sn_mesh=None)` runtime path), tracked + #248 (orphaned Protocol). Not an accidental inconsistency.

## ⭐ THE STANDING LESSON (the recurring tell this carve instances)
A migrate-then-delete of a field DEAD in production but described by comments across the codebase leaves a STALE-DOC BLAST RADIUS in **untouched-but-adjacent hot files** — comments the carve makes false but that live OUTSIDE the diff. Discriminate by TENSE/CITE:
- **PRESENT-TENSE claim about a deleted data contract = MUST-FIX-before-commit** (false-invariant bug habitat; a maintainer re-adds the field "to match the docstring," re-opening the twin). Here: `scheme.py:84-85` ("`StreamingTerms` carries `alpha_in=alpha_out=0.0, tau_mm=1.0`"), `scheme.py:273` ("Step C **will** retire" — future tense for what THIS commit does), `reduced_operator.py:245-246` ("the `alpha_in is None` test discriminates slab from curvilinear" — live discriminator is actually `requires_upstream_angular_state`, `scheme.py:123`).
- **DELETED LINE-NUMBER cross-references = MUST-FIX** (anti-#11; my standing "line-cites→symbol-cites" tell). Here: `pole_angular_closure.py:845/:863/:793` `# Replicates spherical_streaming (reduced_operator.py:681-688) EXACTLY` — cites deleted ranges, "Replicates EXACTLY" asserts a live geometry-side twin that no longer exists (it WAS the bit-id de-risk bridge in Steps A/B; now a phantom). Prefer symbol-ref "formerly replicated ..., retired in #236 Step C".
- **HISTORICALLY-framed contrast against a deleted thing = FOLLOW-UP** (still true as history, only the contrast target is gone). Here: the "uses THIS τ **instead of** `StreamingTerms.tau_mm`" comments at `pole_angular_closure.py:310/485/560/1065` + `geometry.py:1317/1551` + `diamond.py:204/240/315` + `cell_balance.py:46/286`. Fold into the #236 Step-6 archivist doc sweep.

Mechanism: the brief scopes the review to the changed files, but Cardinal Rule 2 + anti-#11 require grepping the WHOLE tree for refs to the deleted symbols/line-ranges — the false claims hide where the diff does not reach. (Mirrors the method-implementer's own LESSON that oracle-readers hide outside the audit's named files; the SAME blast-radius logic applies to DOC refs, not just test refs.)

## Approval conditions issued
Must-fix: `scheme.py:84-85`+`:273`; `pole_angular_closure.py:793/845/863`. Follow-up: `reduced_operator.py:199/206-211/245-246` + the historical-contrast cluster. Production `reduced_operator.py` deletion APPROVED as-is (conditions are doc-correctness only).
