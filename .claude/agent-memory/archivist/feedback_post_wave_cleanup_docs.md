---
name: Post-wave architectural-cleanup docs sweep (Issue #188+#176 / C176.5)
description: 5-element pattern for a small follow-up cleanup that closes a deferral introduced by a larger refactor wave; uses CLOSE-OUT narrative arc with motivation preserved + algebra of where-they-agree-vs-diverge + Option-X/Y rationale + new equation labels for diverging-semantics decomposition
type: feedback
---

# Post-wave architectural-cleanup docs sweep

Branch `feature/bc-curvilinear-realizer-cleanup` 2026-05-11 commit
`188bf9a`. The branch sits on top of `refactor/sn-operator-algebra`
(Wave 12 close of the 12-wave BC trace-law refactor) and closes a
specific deferral that Wave 12 documented as an open follow-up (the
curvilinear `InflowTraceSpace` and the 2-arg shim simplification it
unblocks).

## Rule

A two-issue cleanup that closes a deferral introduced by an earlier
wave (where Wave-N documented the deferral, and the cleanup commits
under "Issue #188 unblocks Issue #176" finally land) gets a
**post-cleanup C-N close-out commit** to the docs, NOT a full
re-narration of the earlier wave. The docs lift load is much smaller
than a wave-summary doc: 5 elements, ~500 LoC.

**Why:** The earlier wave's docs already carry the architecture
narrative. The cleanup commit's job is to flip stale text
("deferred to Issue #176" → "closed by Issue #176"), to add the
algebraic decomposition that ONLY becomes load-bearing after the
cleanup (because before the cleanup, the divergence was masked by
the bypass), and to install a closure note pointing the cleanup
plan. The "novel-extension falsified" step of the close-out arc is
absent because nothing was falsified — Issue #176 was always
intended; Issue #188 was the prerequisite.

**How to apply:**

1. **Key-Facts vacuum-bullet (or equivalent semantic-correction
   bullet) flip.** The post-cleanup state is "uniform realizer
   routing"; the pre-cleanup state was "Cartesian-only realizer +
   curvilinear bypass". The Key-Facts entry has been the
   load-bearing summary since Wave 12; flip it with explicit
   "Post Issue #X+#Y (date)" framing so future sessions can
   timestamp the change.

2. **New section: the algebraic decomposition.** When the cleanup
   makes two previously-equivalent code paths converge into one,
   the precise algebraic statement of where the two semantics
   AGREE vs DIVERGE becomes load-bearing. Add new
   `:label:`s for the decomposition equations:
   - `ordinate-partition-inflow-outflow` (the trace partition
     :math:`\{1..N\} = I_f \sqcup O_f \sqcup T_f`)
   - `vacuum-legacy-vs-trace-correct` (the two functions side by
     side, with annotated agreement on :math:`I_f` and divergence
     on :math:`O_f`)
   These appear in the V&V audit as "documented but no
   implementing code" — correct, because they're decomposition
   statements, not solver claims.

3. **Replace the deferred-work section with a close-out section
   (using the CLOSE-OUT narrative arc).** Keep the original
   motivation as past-tense ("Why the split existed... The Wave 2
   factory deferred curvilinear support because..."), then add
   "Why the split has been removed" + "The architectural sequence
   (Issue #X unblocks Issue #Y)" + "Closure" (commits + plan
   pointers). The forward-pointer (formerly "see Issue #N for
   future work") becomes a backward-pointer (formerly "...will
   ship when..." → "...shipped 2026-05-11 in commits ABC+DEF").

4. **New section: the Option-X-vs-Option-Y rationale.** When the
   cleanup adopts a non-trivial design (Option A keyword-optional
   parameter vs Option B strict, in this case), document the
   rationale at a dedicated `:label:` (`bc-option-a-signatures`)
   with a per-BC behaviour table (3-column: BC name × `quadrature`
   semantics × body). The table is the discoverable summary;
   visitors who want only "what does X.apply look like" land here
   and find the answer without re-reading the close-out narrative.

5. **Anti-pattern catalog extensions.** The earlier wave's
   anti-pattern catalog gets 2 new entries (one for the bypass
   that's been retired, one for the rejected Option B). Numbered
   continuation (#6, #7 etc.) of the existing catalog, not a
   parallel cleanup-specific catalog. Anti-patterns are durable
   because the seductive trap to "preserve the dual mode for
   flexibility" survives in future sessions until explicitly
   flagged.

## Module-docstring touch-up

The cleanup commits themselves already updated per-BC docstrings
(via the implementation work). The archivist's commit only
touches:

- The package `__init__.py` top docstring (one bullet flip in the
  three-layer §16A.3 decomposition list).
- The base ABC (`_base.py`) top docstring + class Notes section
  (rewrites the "transitional `(psi_out, *args, **kwargs)`" prose).
- The trace-space module docstring References section (adds
  pointer to the cleanup plan that lifted the original
  NotImplementedError).

Cross-verify the per-BC files (vacuum / reflective / white / albedo
/ periodic / prescribed_inflow) but expect them to be correct
already — the implementation commits do this work as part of the
Option-A signature migration. If any reference the old
curvilinear-deferral language, fix it; otherwise leave them alone.

## Forward-pointer hygiene

Searches to run after the close-out edits:

```
grep -rn "curvilinear.*deferred\|deferred.*curvilinear\|curvilinear realizer ships\|Issue #176\|2-arg" docs/theory/ orpheus/
```

The expected output AFTER cleanup:

- Every `Issue #176` mention is now historical-context ("Issue #176
  / C176.X" with past-tense framing) or in the close-out section
  itself. NONE should say "until Issue #176 ships" or "deferred to
  Issue #176".
- Every `2-arg` mention is historical ("the pre-#176 2-arg form
  was...") or in the Wave-7 Option-a-vs-b-c retrospective (which
  predates this cleanup and is genuinely about a different
  decision).
- `curvilinear.*deferred` should appear ONLY in references to the
  2-D cylindrical `Mesh2D` deferral (which IS still deferred —
  no 2-D cylindrical sweep exists today).

## Sphinx warning gate

This kind of cleanup MUST hold the warning count at the pre-edit
baseline. The Wave-12 baseline was 9 warnings; the 2026-05-11
baseline observed at C176.5 was 8 (one warning got fixed somewhere
in the Wave 12 ship-state work after the prior session). The
acceptance gate is "count is unchanged from immediately-before-this-
commit baseline", NOT "absolute 9". Document the pre/post in the
commit body so reviewers can verify.

## Files this commit touched (C176.5 reference inventory)

- `docs/theory/boundary_conditions.rst` — primary, ~+360 LoC.
- `docs/theory/discrete_ordinates.rst` — secondary, the BC
  resolution table + abstract-apply description.
- `docs/theory/operator_algebra.rst` — small, the trace-spaces
  consumers list.
- `orpheus/geometry/boundary/_base.py` — top docstring + class
  Notes.
- `orpheus/geometry/boundary/__init__.py` — Layer 1 bullet flip.
- `orpheus/numerics/trace_space.py` — References section in top
  docstring (small).

Per-BC files (`vacuum.py`, `reflective.py`, `white.py`,
`albedo.py`, `periodic.py`, `prescribed_inflow.py`) were already
correct from the implementation commits — Archivist verified but
did not touch.

## When to use this pattern

- A small follow-up commit closes a deferral introduced by a
  larger refactor wave.
- The deferral was load-bearing for some piece of the larger
  refactor's transitional code (here: the dual-mode bound-quadrature
  shim).
- The cleanup commits adopt a non-trivial design (Option X over
  alternatives) — needs a rationale block.
- The cleanup makes two previously-equivalent code paths converge
  into one, and the algebraic decomposition of where-they-agree-
  vs-diverge becomes load-bearing.

Inapplicable if:

- The cleanup is bit-identical (no new algebraic statement) — just
  flip stale text.
- The cleanup is large enough to be a full wave (then use the
  Wave-12-style synthesis pattern instead).
- The cleanup is a falsification close-out (a path that was tried
  and didn't work) — use the 9-step CLOSE-OUT narrative arc
  instead.
