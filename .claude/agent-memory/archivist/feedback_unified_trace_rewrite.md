---
name: unified-trace worked-example rewrite (Issue #223)
description: Closing the deferred holistic rewrite that C4/C5 left as a stale-narration .. note:: — modernizing boundary_conditions.rst's BC-resolution worked example from the pre-#188 split-trace (InflowTraceSpace/OutflowTraceSpace, _inflow_trace/_outflow_trace, for_face(inflow_trace=...)) onto the one-space-two-selectors unified TraceSpace. The DIRECT sequel/payoff of the #223-pointer in [[feedback_axis_primary_3d_admission_docs]].
type: feedback
---

# Unified-trace worked-example rewrite (Issue #223)

This task EXECUTED the deferral that the C5 memo
([[feedback_axis_primary_3d_admission_docs]]) explicitly predicted and
pinned: it said "the worked-example CODE BLOCK (~886-909) is #223's,
NOT trivial to half-fix — the signature itself changed, so renaming
JUST the method leaves a self-inconsistent block; the C5-scoped move
was a `.. note:: Stale narration ... (tracked by #223)` ABOVE the
block". #223 is the holistic rewrite that note promised. The triage
the brief carried had ALREADY located all 7 edit sites + the
leave-alone set — the C5 memo's stale-ref-sweep section was the
ground-truth map.

## The decisive doc moves (reusable for a "modernize a stale worked
example onto a unified type" rewrite)

1. **The stale `.. note::` becomes the rewrite's seed, not its
   leftover.** C5 had parked a `.. note:: Stale narration ... tracked
   by #223` above the dead code block. The rewrite REPLACES the whole
   step body (header→last code-comment) with the CURRENT narration,
   and the stale note is REWRITTEN into a compact `.. note::
   Historical — the pre-#188 split trace` APPENDED to the end of the
   now-current step. Don't keep both the "this is stale" warning AND
   the stale code — collapse: current narration up top, one historical
   tombstone at the bottom naming the old API as code literals.

2. **Narrate the CURRENT path step-by-step from the code, with the
   WHY.** The Sphinx-as-brain bar for a worked example: (a) the actual
   factory call (`TraceSpace.from_quadrature_and_layout(self.quad,
   self.boundary_face_layout)`), (b) the data shape it produces (the
   `(n_faces, N)` SIGNED-projection table — emphasize *one signed
   table*, NOT two boolean masks), (c) the predicate arithmetic worked
   for one named face (`xmin`: Ω·n = −μ_x, inflow ⇔ −μ_x < −ε ⇔ μ_x >
   ε), (d) the single `trace=` arg threaded to `for_face` + the
   on-demand `inflow_indices_for_face` extraction, (e) the WHY ONE
   SPACE rationale (predicate-not-identity; removes the inflow/outflow
   drift bug class; gives #208 a single |Ω·n| boundary inner product).
   Read `_resolve_bcs`/`_resolve_one`/`for_face` FIRST — the prose
   mirrors the code line-for-line.

3. **face= flip is `"left"`→`"xmin"` ONLY where presented as current;
   INPUT-mesh `bc_left=BC(...)` STAYS.** Same two-surface discipline
   as the C4 memo ([[feedback_face_name_carve_docs]]): Step 1's
   `Mesh1D(... bc_left=BC("vacuum"), bc_right=...)` is the USER-FACING
   geometry DECLARATION (input-mesh dataclass field) → KEEP. Only the
   SOLVER-resolved `face="left"` kwarg in the `for_face` call is the
   current-API spelling → flip to `"xmin"`. Add a parenthetical at the
   attr-doc site noting `"left"`/`"right"` were pre-C4 aliases of
   `"xmin"`/`"xmax"`.

4. **Two-typed-subclasses → one-type-two-selectors is a `.. note::`,
   not a silent swap.** The Layer-1 trace structure presented two LIVE
   `:class:` typed FunctionSpace subclasses. Rewrite to ONE
   `:class:`TraceSpace`` whose `inflow_indices_for_face` /
   `outflow_indices_for_face` are SELECTORS, + a `.. note:: One space,
   two selectors (#205/#201)` carrying the unification observation
   (incoming-vs-outgoing is a PREDICATE on shared data, not a property
   of the space's identity). The module docstring of trace_space.py is
   the canonical source for this prose — it has the exact framing
   ("inflow and outflow are operations on a single space, not two
   spaces"); mirror it, don't reinvent.

5. **De-role dead `:class:` refs in scope; LEAVE the ones the triage
   walled off.** A dead `:class:`InflowTraceSpace`` (class never
   existed in code post-unification) renders PLAIN TEXT (no -W
   warning) but is still a Cardinal-Rule-1 staleness bug AT A LIVE
   SITE. Fix the in-scope ones to literals (``InflowTraceSpace``
   "then-named") or to the live `:class:`TraceSpace``. But the triage
   explicitly walled off the deep historical sections
   (curvilinear-realizer-unification ~2237-2300, close-out ~2748,
   gate-only ~3036) — those keep their dead `:class:` xrefs because
   they sit in unambiguously historical prose and the acceptance
   criterion is "no LIVE (non-historical) reference remains", not
   "zero lexical hits". Respect the triage's leave-alone boundary;
   list each remaining hit + its justification in the L12 closeout.

## Verification grep (L12) — the acceptance criterion is non-live,
not zero-hit

`grep -n "_inflow_trace\|_outflow_trace\|InflowTraceSpace\|inflow_trace="`
WILL still return hits after a correct rewrite — every one must be a
code literal (double-backtick) inside explicitly-historical prose OR a
dead `:class:` in the triage's walled-off historical sections.
`probe_inflow_trace` (a variable name in the PrescribedInflow
realization-map row) is NOT a split-trace type — it's a false-positive
of the grep, leave it. Enumerate each surviving hit with its one-line
justification in the closeout.

## Build gate (MAIN checkout, this session)

Worked in the MAIN checkout (NOT a worktree — a prior worktree attempt
hit Edit permission denials). Baseline = **1 warning** (the
pre-existing `orpheus/geometry/mesh.py` `Mesh1D.from_geometry`
`:paramref:` ERROR — needs sphinx-paramlinks, out of scope) — matches
the AGENT.md standing baseline, NOT the C4/C5 worktree's 11 (those 11
were `_generated`/`.h5` env-artifacts of a fresh worktree; the MAIN
checkout's env is materialized → 1). Forced `-E` cold build: EXIT 0,
"build succeeded, 1 warning". My 7 edits added ZERO new
WARNING/ERROR/CRITICAL; grep of the `-E` log for
`boundary_conditions|InflowTraceSpace|inflow_trace|Title underline|
Inconsistent title|trace_space` returned only `reading/writing`
progress lines. `trace_space` module is NOT automodule'd → my
`:class:`TraceSpace``/`:meth:` refs render plain text (matches how the
original `InflowTraceSpace` refs rendered — same render, correct
name). Referenced anchors (`bc-trace-structure`, `bc-face-name-carve`,
`sn-c5-geometry-blind-trace`) all exist intra-doc — grep
`^\.\. _<anchor>:` BEFORE building.

## Quality scores (this task)

| Dimension | Score | Note |
|---|---|---|
| Derivation depth | 5 | current-path step-by-step: factory call + signed-projection table shape + xmin predicate arithmetic (Ω·n=−μ_x → μ_x>ε) + WHY-one-space rationale (predicate-not-identity, drift-bug-class removed, #208 inner product) |
| Cross-references | 4 | every live SNMesh/SNMethodSpace/TraceSpace symbol repointed; trace_space refs plain-text (unsurfaced module — matches pre-edit render); 3 intra-doc :ref: anchors verified-exist |
| Numerical evidence | n/a | doc-staleness rewrite, no numbers change (structural narration) |
| Failed approaches | 5 | the pre-#188 two-space split preserved as an explicit historical tombstone naming the old API as literals + WHY it was retired (inflow/outflow drift) |
| Code traceability | 5 | every claim mirrors _resolve_bcs/_resolve_one/for_face read line-for-line; signature change (inflow_trace=/outflow_trace= → trace=) verified against current dataclass fields |
| Derivation source | n/a | architecture/staleness, no SymPy derivation needed |
