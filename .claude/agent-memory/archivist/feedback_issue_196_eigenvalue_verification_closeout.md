---
name: issue-196-eigenvalue-verification-closeout
description: Issue #196 — archiving the VERIFICATION + regression-gate step that CLOSES the last ERR-026 manifestation (#7), one session after the ERR-058 (#195) fix close-out; the bit-identical(fixed-source)-vs-floor-equivalent(eigenvalue) distinction, the two-layer L0/L1 close history (ERR-048 vs ERR-058), retitle-keep-anchor, stale-canary-path repoint
metadata:
  type: feedback
---

Issue #196 in `docs/theory/discrete_ordinates.rst`: the **verification +
regression-gate** close-out that VERIFIES the prior-session ERR-058 (#195) fix
actually closed ERR-026 manifestation #7 (the curvilinear SI-vs-Krylov O(h)
eigenvalue asymmetry) and PINS it. Distinct archival shape from its predecessor
[[err058-success-closeout-supersedes-phase-chain]]: that was the FIX close-out
(replaced wrong seeds, superseded the Phase A-F chase); THIS is the
**verification layer stacked on top** — measure the equivalence, add the gate,
formally retire the family. A fix close-out and a verification close-out are
two SEPARATE archival acts on the SAME ERR family across two sessions; the
verification one ADDS a focused subsection inside the existing fix-close-out
section + retracts the bug-era O(h)-gap tables, it does NOT re-write the fix
close-out.

**The load-bearing nuance: bit-identical (fixed-source) ≠ floor-equivalent
(eigenvalue).** The prior #195 close-out blanket-said "SI ≡ Krylov bit-identical
on the curvilinear MMS ladders." TRUE only for fixed-source
(`solve_sn_fixed_source` — one operator, one quadrature, within-group `L.solve`
vs Krylov-on-`apply` = same L⁻¹ arithmetic = same bits). FALSE for the
EIGENVALUE path (`solve_sn`): there the inner is wrapped in POWER ITERATION, and
SI vs Krylov are different iteration schemes converging to the SAME correct
fixed point only to ~inner_tol (~1.9e-11 keff / ~2.4e-10 flux-shape) — same
PHYSICS, not same ARITHMETIC. **Rule:** when archiving any "two inner solvers
agree" claim, ALWAYS ask whether the entry point is fixed-source (→ bit-
identical possible) or eigenvalue (→ floor-equivalent, NEVER call it bit-
identical). Gave it a dedicated `^^^^` H3 + a `.. _<anchor>-bit-identical-vs-floor:`
target so other sites can `:ref:` the distinction. The existing blanket bit-id
claims elsewhere were ALL in fixed-source/MMS-ladder contexts (so correct) — did
NOT flip them; only added the eigenvalue caveat where the eigenvalue path is
named (status banner + the new subsection).

**Two-layer L0/L1 close history is the WHY a fix didn't suffice the first time.**
ERR-048 (Phase G Step 2) closed only the L0 FLAT-FIELD twin-agreement (SI sweep
≡ apply-matvec on the homogeneous streaming gauntlet); the L1 HETEROGENEOUS
EIGENVALUE O(h) asymmetry PERSISTED (→ #196 stayed open) because the shared
seeds were still exact-on-flat / O(1)-wrong-on-non-flat (the ERR-058 defect).
ERR-058 was the TERMINAL fix. The doc MUST spell this ladder out — a reader
asking "didn't ERR-048 already fix the sweep?" needs the L0-vs-L1 distinction or
they'll think #196 was redundant. Use vv-principles level terms VERBATIM (L0
term / L1 equation), never paraphrase.

**Retraction-tombstone the bug-era O(h)-gap TABLES, keep the numbers.** The page
carried a "Post-Phase-F SI-vs-Krylov convergence" table reading the residual
O(h) gap as a benign discretisation artefact of "two methods now solving the
same equation." That INTERPRETATION is wrong (they did NOT yet solve the same
equation — shared seed defect). Added a `.. note:: **Retraction (date, #196).**`
above the table + retitled the table caption "(bug-era — gap closed by
ERR-058/#196)" + flagged the interpreting prose "(Bug-era reading, retracted)".
Numbers STAY as historical evidence; only the interpretation gets the tombstone
(same discipline as the PARTIAL/OPEN variant, here applied to a CLOSED follow-up).

**Retitle the section, KEEP the anchor.** The section `.. _sn-phase-f-residual-o-h-open:`
"What stays open after Phase F (ERR-026 manifestation #7)" was retitled to
"ERR-026 manifestation #7 — CLOSED by ERR-058 (#195), verified + pinned by #196"
WITHOUT renaming the anchor. Two existing intra-doc `:ref:`s point at the anchor;
keeping it means they auto-pick-up the new title text (verified in rendered HTML:
the `std-ref` link text became the new title). Renaming the anchor would have
forced a tree-wide sweep for zero benefit. Close-out label-rename hygiene from
AGENT.md, INVERTED: when the CONCEPT survives and only the STATUS flips, keep the
anchor, change the heading.

**Stale canary test-PATH sweep (L1 — renders plain-text, NO warning).** The page
cited `tests.sn.test_phase_c_crosscheck.<fn>` but the live module is
`tests.sn.verification.analytical.test_phase_c_crosscheck.<fn>` (3 sites). Same
trap as the [[err058-success-closeout-supersedes-phase-chain]] sweep — grep the
cited dotted path against the live FS (`find tests -name <file>.py`) and repoint
in the same edit. The canary (`test_phase_e_..._crosscheck`) is now a PLAIN L1
test (xfail removed on disk — verify: no `@pytest.mark.xfail`, just `l1`+`verifies`)
and is the structurally-INDEPENDENT anchor (Variant-α semi-analytical reference
at 8%/12%): SI≡Krylov alone is twin-path agreement (vv L11 necessary-not-sufficient);
the k_inf homogeneous legs (closed-form) + the Variant-α cross-check (semi-
analytical) supply the independent ground. ALWAYS pair a twin-path-equivalence
close-out with its structural-independence leg in the prose.

**No new eq-labels this task → no vv-status directives needed.** Unlike the #195
fix (3 new representational eq-labels needing `.. vv-status: documented`), #196
added only SECTION anchors (`.. _sn-issue-196-*:`), resolved by `:ref:` not
`:eq:`. Section anchors are NOT in the eq-label namespace, carry no vv-status,
and don't trip `audit --strict`. The new test's `verifies(...)` set reused
PRE-EXISTING labels (transport-spherical, matrix-eigenvalue, …) so no orphan
verifies-targets created.

**Build gate: MAIN-checkout baseline = 1 warning** (the mesh.py `from_geometry`
`paramref` ERROR, out of scope). Forced `-E` pre+post; warning SET unchanged.
This was on a feature branch in the MAIN checkout (not a worktree) so Nexus WAS
available — `staleness()` correctly flagged `theory/discrete_ordinates` stale
(code landed, doc hadn't recorded #196); this edit clears it. Confirmed the new
`:ref:` anchors resolve as live hyperlinks in rendered HTML (grep `href="#<anchor>"`),
not the silent plain-text render.

Quality self-score: derivation depth 4 (verification close-out, not a new
derivation — the measured-evidence + gate-threshold tables are the depth),
cross-refs 5 (every code symbol/test/anchor verified on disk; `SNSolver.inner_solver`
checked against solver.py:757), numerical evidence 5 (eigenvalue + homogeneous +
gate-threshold-vs-bug-era tables, all verbatim from task background + test
docstrings), failed approaches 5 (Option a/b/c kept as clearly-marked bug-era
history with WHY-each-was-wrong), code traceability 5, derivation source N/A
(verification task, no derivation). Related: [[err058-success-closeout-supersedes-phase-chain]]
(the sibling FIX close-out this verifies), [[feedback_retirement_docs]] (flip
will→did + keep-tests-pinned-labels).
