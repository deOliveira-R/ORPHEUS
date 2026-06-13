---
name: err058-success-closeout-supersedes-phase-chain
description: ERR-058 #195 — archiving a FIX-FOUND close-out that supersedes a long multi-phase OPEN-loop (Phase A-F / "ERR-026 PARTIAL CLOSURE"); preserve-history + retraction-tombstone discipline, vv-status on representational eq-labels, audit-strict-pre-existing report-not-fix
metadata:
  type: feedback
---

ERR-058 / Issue #195 in `docs/theory/discrete_ordinates.rst`: a curvilinear
closure-seed fix that CLOSES a family the page had tracked across ~3000 lines
of Phase A-F + Wave-E narrative under the name "ERR-026 PARTIAL CLOSURE". The
distinctive archival shape: a **fix-found close-out that supersedes a long
multi-phase OPEN loop**, not a falsification close-out and not a fresh feature.
The 9-step CLOSED arc in AGENT.md is for FALSIFICATION; this is its sibling
("the chase was real, the resolution arrived, the chase's terminal claims are
now wrong but its REASONING is preserved").

**Rule: a fix that closes a multi-phase chase = ADD one close-out section +
RETRACTION-TOMBSTONE the chase's terminal claims; do NOT rewrite the chase.**
**Why:** the Phase A-F reasoning (why anyone tried a Carlson inward sweep, an
apply-vs-sweep twin audit, a Krylov-default flip) is pedagogically load-bearing
— a future session asking "why was this so hard?" must find the path. Only the
TERMINAL claims each phase reached (PARTIAL CLOSURE, default→krylov,
"pre-asymptotic transient") were falsified. **How to apply:** new H1/H2
close-out section at the END of the phase block (here after Phase F, before the
next chapter); inside it, "Motivation PRESERVED" flips tenses ("is expected to"
→ "was expected to") on the phase claims without deleting them; then a
`.. note:: **Retraction/Superseded (date, #N).**` tombstone immediately above
EACH stale terminal claim site, with (a) what the claim was, (b) one-line why
it's wrong, (c) forward-pointer to the close-out. Numerical tables STAY (bug-era
evidence); only the INTERPRETATION gets the tombstone.

**Tombstone placement for a deeply-nested phase chain — cover the scan path
AND the deep site.** A reader scanning a phase block's top-of-section step
summary (item "default flipped to krylov") would hit the stale claim ~700 lines
before the deep tombstone. So: one forward-pointer `.. note::` at the TOP of the
phase block (the "Phase D scope" entry point) listing which terminal decisions
were reverted, PLUS inline tombstones at each deep site. Don't rely on a single
close-out section to catch a scanning reader.

**Two close-out flavors in one fix:** the wrong-fixed-point family CLOSED (#195),
but the SAME fix surfaced a residual ANISOTROPIC angular floor that stays OPEN
under a NEW issue (#229). Pattern: a closed close-out can re-pin existing xfail
markers from the closed issue to the new one. The isotropic MMS tests' xfail
came OFF (verify on disk: markers removed, `catches("ERR-058")` added); the
ANISO + prescribed-inflow-sphere tests STAY xfail but their reason-string +
issue-pointer moved #195→#229 (verify on disk too — the test files carry the
authoritative re-scope). Doc must mirror BOTH transitions: "markers OFF" for
isotropic, "re-pinned to #229" for aniso, with the floor-vs-quadrature evidence
table (floor drops under quadrature refinement = the structural proof the
residual gap is angular, not a wrong fixed point).

**vv-status on the new eq-labels is MANDATORY or `audit --strict` regresses.**
New `.. math:: :label:` blocks that are structural/representational identities
(a boundary-continuity condition; a FALSIFIED proxy-source definition kept as
the diagnosed defect; an operator-consistent linear-map seed) are NOT solver
claims → each needs `.. vv-status: <label> documented` (a comment directive
placed AFTER the math block, blank-line-separated) + a `.. (vv-status rationale)`
comment naming the verifiable content (the per-ordinate operator-admission gate
+ the strategy-owned adjoint bit-identity). Without it the label counts as a
testable orphan and trips `--strict`. A prose `.. note::` EXPLAINING the
vv-status intent is NOT enough — the harness reads the DIRECTIVE, not prose.
Confirmed: documented-count 216→219 after adding the 3 directives; my labels
became 0-of-N orphans.

**`audit --strict` exit code: rebuild `-E` BEFORE baselining it.** The audit
reads `docs/_build/html/_nexus/graph.db`. A first audit run against a STALE
graph (built before the branch's test markers landed) gave EXIT=0; the TRUE
baseline after a `-E` rebuild was EXIT=1 (55 orphan-equation labels + a
MISSING-ERR list, all pre-existing, none ERR-058). Lesson: to know whether YOUR
docs work moved the strict gate, stash your edits, `-E` rebuild, run audit
(true baseline), pop, `-E` rebuild, run audit (post). Compare orphan/documented
COUNTS, not the raw exit code. Report a pre-existing exit-1 as pre-existing
(task said "report rather than fix") — confirmed my edits added 0 orphans.

**Read order that worked (authoritative-first):** error_catalog.md ERR-058 entry
(the canonical two-manifestation mechanism + post-fix evidence verbatim) →
`gh issue view 195 --comments` (the re-scope premise-refutation) → the new
code's RICH docstrings (psi_half_angle_seed.py AngularEdgeExtrapolation class
docstring carried the full math + the "why this replaces Carlson" narrative;
operator.py `_compute_LpC` coupled-pole comment block; solver.py docstring
History block) → the new test module docstring (the per-ordinate-vs-scalar
blindness story) → `gh issue view 229`. The docstrings + error_catalog were so
complete the close-out was 90% synthesis, ~10% new derivation (the
recurrence-threads-exactly linear-in-μ step I expanded from the one-line
docstring claim). Cardinal Rule 3 working as designed: rich docstrings ARE the
prose seed.

**Stale test-PATH sweep alongside the claim sweep.** The page cited
`tests.sn.test_mms_curvilinear...` but the live path is
`tests.sn.verification.mms.test_mms_curvilinear...`. When a verification section
is already being edited for the claim flip, repoint the test `:func:` paths in
the same edit (grep the cited dotted path against the live filesystem). These
render plain-text if wrong (L1 — no warning), so the build won't catch them.

Quality self-score: derivation depth 5 (full linear-in-μ threading proof from
the docstring seed), cross-refs 5 (every code symbol + test + eq-label verified
to exist before citing), numerical evidence 5 (post-fix ladders + per-ordinate
residual + #229 floor-vs-quadrature, all verbatim from error_catalog/#229),
failed approaches 5 (3 refuted hypotheses recorded), code traceability 5,
derivation source 5 (from psi_half_angle_seed.py + error_catalog, not
hand-written). Related: [[feedback_retirement_docs]] (flip will→did),
[[feedback_phase_d_carlson_seed_narrative]] (the ORIGINAL Phase-D archival this
supersedes), AGENT.md "Close-Out Narrative Arc" (the falsification sibling).
