---
name: type-confinement docstring sync (subtype confined to one role; supertype admitted at the boundary)
description: Syncing docstrings/comments/error-strings to a surgical type change that CONFINES a subtype to one role (the driver iterate) while the boundary now admits the SUPERTYPE — the timed iterate passes via inheritance (`TimedFullField(FullField)`), a bare supertype is the degenerate (history_depth=0), and a long-retired sibling arm (`AngularFlux`) is dropped from the signatures. The keep-vs-flip rubric: flip INPUT-boundary prose to the supertype; PRESERVE OUTPUT prose describing the subtype (the driver carrier stays). Plus the live-code-corrects-the-stale-docstring instance (`rhs.boundary` not `initial_guess.boundary`) and the two dead-ref kill.
type: feedback
---

A code-final docstring/comment/error-string sync task (NOT a theory-page
write). Branch `refactor/operator-inverse-algebra`, campaign P4.5 "W-C —
TimedFullField→FullField boundary confinement". The shape recurs whenever
a surgical carve CONFINES a subtype to one role and admits the supertype
at a boundary. The durable discipline:

**The semantic to internalise (so prose is accurate):** a subtype
(`TimedFullField`, the history-bearing comonad carrier) was confined to
its genuine role — the DRIVER ITERATE OUTPUT of `solve`/`solve_moments`.
The INPUT boundaries (solve/residual `rhs`/`initial_guess`/`q_ext`) now
speak the timeless SUPERTYPE (`FullField`). The timed iterate still flows
in via inheritance (`class TimedFullField(FullField)`), and a bare
supertype is admitted as the `history_depth=0` degenerate. A long-retired
sibling (`AngularFlux`, killed earlier by a runtime guard) is now DROPPED
from the signatures entirely. The residual OUTPUT also becomes the bare
supertype (a one-shot balance defect carries no history — it was already
built with `history_depth=0`).

**The keep-vs-flip rubric (the load-bearing decision):**
- **FLIP to the supertype:** every INPUT-boundary param prose
  (`rhs : AngularFlux or TimedFullField` → `rhs : FullField`), the
  dual-dispatch narrative ("Both branches dispatch via runtime
  isinstance"; "Output type matches input type" — DELETE wholesale, the
  arm is gone), the residual OUTPUT/return-type prose, and any
  user-facing error-string naming the boundary type
  (`...solve(TimedFullField): ...` → `...solve(FullField): ...`).
- **PRESERVE (do NOT flip):** prose describing the SUBTYPE in its
  surviving role — the `-> "TimedFullField"` iterate output, "the
  history-bearing TimedFullField comonad", "A TimedFullField IS a
  FullField" inheritance notes, the driver-owns-history threading note,
  and `psi : TimedFullField` (the reconstructed full-angular flux
  argument, which IS the driver-side iterate). The brief will say
  "PRESERVE prose describing the X iterate/comonad" — take it literally.
- A param that reads the subtype's data CONTAINER-AGNOSTICALLY (a sweep
  reading only `.bulk.values` via an extractor) flips to the supertype
  AND gains a "history-blind" note — the extractor sees only the bulk
  ndarray, so the timed history is irrelevant to it.

**Live-code-corrects-the-stale-docstring (L-001, the sharpest instance):**
the dead docstring claimed the boundary seed came from
`initial_guess.boundary` via two now-nonexistent helpers. The LIVE method
body seeded from `rhs.boundary` (a W-C inline comment even said "REPLACES
the pre-extraction seed-from-`initial_guess.boundary`"). I rewrote to the
TRUE source + the TRUE mechanism (`boundary_buf.face_view(name)[:] =
seed_boundary.face_view(name)`, found by grepping the method body for
`face_view`). NEVER transcribe the stale docstring's claim — read the
body, name the real source variable.

**The two dead refs were the headline deliverable (L-002).** The stale
paragraph named `:meth:`TimedFullField.to_legacy_angular_flux`` and
`:func:`_copy_boundary_face_state`` — both grep to ZERO defs in
`orpheus/`, i.e. dangling refs rendering plain-text with NO `-W` warning.
The build gate is BLIND to them; the grep gate (`grep -rn "<symbol>"
orpheus/` → only the two docstring sites) is the real check. Confirm
"GONE ✓" by grep post-edit, not by the clean build.

**Scope-adjacent dead module-path repointing (in-scope only).** The
class-level "The .solve API" section carried the dead path
`:class:`~orpheus.sn.angular_flux.AngularFlux`` (the module is
`orpheus.transport.fields.angular_flux`; `orpheus/sn/angular_flux.py`
does not exist). Where a line I was ALREADY editing under brief scope
carried that dead path, I repointed it (Cardinal-Rule-1, within the
edit). 4 MORE instances of the same dead path live in OTHER solver.py
methods (eigenvalue/SI driver docstrings, not the W-C solve/residual
boundary) — I FLAGGED them as a pre-existing follow-up, did NOT fix
(scope creep into untouched methods). Flag-don't-fix the out-of-scope
sibling staleness; fix only what your edit already touches.

**Out-of-scope stale Notes flagged, not fixed.** A matvec `apply` Notes
line ("Depth B D-I.3d — `apply` accepts ONLY `TimedFullField`") is stale
vs the live `apply(psi: "FullField")` + the #257 S8a guard — but it's the
MATVEC contract, not the W-C solve/residual boundary. Flagged, left.

**Gate (this was a CODE-FINAL task → import sanity matters):** baseline
`-E -W` was EXIT=0 / zero WARNING|ERROR|CRITICAL (this branch has NO
mesh.py `:paramref:` ERROR — the AGENT.md "baseline=1" note is checkout-
specific; re-measure, don't assume). Post-edit identical. ALSO ran
`.venv/bin/python -c "import <the 3 modules>"` — a malformed `r"""`
docstring breaks import but can still build Sphinx from a stale pickle;
the import check is the cheap proof the prose edits didn't corrupt a
raw-string. The combined `git diff` mixes the pre-session W-C code with
my prose; prove "docstrings/strings only" by confirming every flagged
signature/guard/return line is W-C-origin (it was final before the
session per the brief), not something an Edit of mine touched.

Quality self-assessment (Directive 3): this is a sync, not a derivation
— Derivation depth / Numerical evidence are structurally N/A (no math
authored, no flux moves). Code traceability + Cross-references = 5 (every
flipped type names the canonical `:class:` path; the new
`_initial_guess_values` cross-ref resolves). The discipline that carried
it: L-001 (the `rhs.boundary` correction) and L-002 (the dead-ref grep
gate). Weakest applicable dimension: none below 4 — the task was narrow
and the live-code-read caught the one trap (stale seed source).
