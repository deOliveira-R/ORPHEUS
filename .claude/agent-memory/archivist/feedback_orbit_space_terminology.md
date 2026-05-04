---
name: Orbit-space M/G terminology sweep
description: Pattern for replacing loose "topological similarity"/"topology family" with the precise "orbit-space M/G classification" term across DocStrings + Sphinx + memory; preserve `topology` attribute as a code-name with upgraded definition; add ONE definitional aside in landing page that downstream pages :ref: to
type: feedback
---

When the user introduces a precise mathematical term to replace a
loose informal one already entrenched in the codebase, follow the
**add-aside-then-bridge-then-sweep** pattern, not blind find/replace:

**Why:** Codebases evolve a working vocabulary that gets used
internally as a code-name (e.g. `CurvilinearGeometry.topology`
returning string `"two_surface"`). Renaming the attribute is a
method-implementer concern (potentially breaking change). Renaming
the prose ABOUT the attribute is an Archivist concern (safe). The
two must be sequenced — prose-precision first establishes the
canonical definition; attribute rename can follow without breaking
the conceptual chain.

**How to apply** when a user says "X is actually Y, fix the
docstrings + Sphinx":

1. **Triage the existing usage**: search for the loose term across
   `orpheus/`, `docs/theory/`, `.claude/agent-memory/`. Inspect each
   hit's context — some are CORRECT usages of the loose term in its
   precise sense (e.g. "boundary topology" of a multi-region
   problem); others are loose usages where the precise term would
   be clearer (e.g. "6-geometry × 2-topology family"). Case-by-case
   judgment, not blanket find/replace.

2. **Pick the landing page** (the canonical place for one
   definitional aside) — usually a `theory/reference_solvers.rst`
   or `theory/index.rst` style master page. Add a `.. _<label>:`
   anchor + a section with the formal definition + a table mapping
   each instance (geometry/case/etc.) to the formal structure +
   prose explaining why this signature drives the existing
   architectural decisions.

3. **Bridge in adjacent pages**: in pages that consumed the loose
   term, replace ONE early occurrence with the precise term + cross-
   reference (`:ref:` to the new label) so a reader landing there
   knows where to find the canonical definition. Don't replace EVERY
   instance — just the introductory ones.

4. **Upgrade the code-name docstring**: where the loose term was
   used as a Python attribute name (e.g. `topology`), keep the
   attribute but rewrite its docstring to explicitly cite the
   precise term: "X is the human-readable label for Y" or "X is the
   M/G endpoint count made human-readable". This frames the
   attribute as the API surface of the precise concept, not the
   precise concept itself.

5. **Light-touch memory updates**: hit only the most-prominent
   memory descriptions/index entries, not every body line. Memory
   files are durable; a single-line precision uplift in the
   description is enough for the next agent reading it to know the
   correct vocabulary. Body content (numerical evidence, ranked
   candidate lists) stays unchanged.

6. **Verify Sphinx -W clean** + key tests pass + the new label
   appears in `docs/_build/html/searchindex.js`. The label IS the
   contract — downstream pages and future agents will :ref: to it.

**Concrete from this session (2026-05-03)**:
- Loose: "6-geometry × 2-topology family" / "topological analog of
  the asymmetric slab" / "Topologically rank-1/rank-2"
- Precise: "6-geometry × 2-orbit-space-class family" / "orbit-space
  M/G analog of the asymmetric slab" / "Closure rank-1/rank-2
  (orbit-space class: one-surface-compact / two-surface)"
- Definitional aside: `docs/theory/reference_solvers.rst`
  §`orbit-space-m-g-classification` with a (G, M/G dimension, M/G
  shape, endpoints) table per geometry — 102 LoC table + prose.
- `CurvilinearGeometry.topology` property docstring upgraded to
  cite M/G endpoint count + new label; attribute stays.
- Two clean commits, Sphinx -W clean, 1376 tests collect, 22+3
  spot-check tests pass.

**Anti-pattern caught here**: don't replace `Topologically rank-N`
naively with `Orbit-space-class rank-N` — that loses the
"closure-rank" semantic. Better: split the concept — keep "closure
rank-N" as the rank vocabulary, add "orbit-space class:
one-surface-compact / two-surface" as a parenthetical that names
the M/G structure that DETERMINES the rank. Preserves both threads.

**Self-assessment for this task** (Directive 3 rubric, 1-5):
- Derivation depth: 4 — the M/G aside derives endpoint-count → rank
  for each geometry, with the symmetry group named explicitly. Could
  be 5 with worked geometric examples (e.g., a literal R²-translation
  trace).
- Cross-references: 5 — every relevant theory page now :ref:s
  to the new label; docstrings cite Sphinx section name.
- Numerical evidence: n/a — pure terminology fix, no new numerical
  claims.
- Failed approaches: n/a — no failed approaches in this session.
- Code traceability: 5 — every code attribute that uses the term
  is now linked to the definitional aside.
- Derivation source: n/a — no equation derivation; the M/G
  framework is standard differential geometry, not project-derived.

Weakest dimension: numerical evidence — but this is structurally
absent for a terminology fix, not a deficit.
