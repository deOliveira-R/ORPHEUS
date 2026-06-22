---
name: Stub-expansion sweeps for family-of-geometries pages
description: Patterns for end-of-session rich-narrative expansion of method-implementer stubs across a multi-geometry family page (peierls_greens.rst Phase 1/2/3 sweep, 2026-05-02)
type: feedback
---

End-of-session "expand all stubs" sweeps on a multi-geometry family page (sphere/cylinder/slab/slab-asym/hollow-sphere/annulus on the same page) work best when the per-geometry sections share a SHAPE rather than just a topic. Five subsections per geometry, in a fixed order:

1. **Lead paragraph** referencing the closeout memo by file path and stating the verdict in one sentence (e.g., "shipped clean on first try, V_α1 at machine precision in 1 iter, zero ERR-NNN").
2. **Phase-space and conserved invariants** (the foundation-tagged Branch-1 SymPy module's algebraic content).
3. **Resolvent / closure** (the load-bearing rank-1 vs rank-2 algebra; cross-link to `.variant_alpha_core` shared primitives via `:func:`).
4. **Operator-level identities V_α1/V_α2/V_α3** (one paragraph per `derive_*` SymPy function, naming the function as `:func:` and the proof structure).
5. **Verification status** (k_inf-exactness, MG, vacuum, plus geometry-specific gates like method-of-images / R_in→0 reduction).
6. **Source code, tests, and provenance** (Branch-1 SymPy module, Branch-2 production module, foundation tests, L1 tests, closeout memo, cross-domain frame memo if applicable, plan).

**Why:** When all six per-geometry sections share this shape, a future reader can navigate any single section, then spot the *structural* analogues across geometries trivially. Cross-cutting sections (Phase 2 unification, V_α2 strengthening) sit between the shape and the data, naming the unifying objects.

**How to apply:**
- BEFORE the sweep, decide which cross-cutting sections will be added (e.g., the operator-theoretic unification + V_α2 strengthening across geometries) and where they go in the page (between Phase 1 and Phase 3 sections is natural — they reference Phase 1 evidence and forward-prepare for Phase 3 algebra).
- WHEN expanding each per-geometry stub, use the closeout memo as the SOURCE; the memo records numerical evidence, closeout verdicts, and bug-history. The Sphinx page presents this evidence with full algebra and cross-references; do NOT re-derive math that's already in the closeout's "what was shipped" section.
- PRESERVE every `:label:` from the original stub. Method-implementer stubs allocate equation labels for cross-reference; expansion replaces the surrounding prose without touching the labels.
- AFTER the sweep, ALSO update the page's Key Facts / "What this prototype does NOT cover" bullets — they are typically stale because the original page was written when only one geometry existed. Family-overview tables at the top of the page replace single-geometry "Scope" headers.

**The duplicate-equation-label trap.** When expanding a stub by inserting a NEW subsection with the SAME `:label:` as an existing equation, Sphinx errors: `duplicate label of equation X`. The fix: pick the canonical placement (typically inside the new richer subsection), DELETE the original bare equation. Verified 2026-05-02 in the annulus expansion: I added a "Through-ray rank-2 resolvent" subsection with the rank-2 resolvent equation, but the original stub already had the same equation with the same `:label:` 30 lines earlier. Sphinx caught it; the fix was to delete the original.

**Title-underline length re-trap when expanding subsection titles.** Adding a subsection title like "Operator-level identities (V_α1_hollow_sph / V_α2_hollow_sph / V_α3_hollow_sph)" — 79 code points (with Greek α counting as 1 code point each). The underline must be ≥ 79 dashes. Don't trust `awk '{print length}'` (counts bytes, with α taking 2 bytes); use `python3 -c "print(len('...'))"`.

**Per-geometry verdict subsections preserve the closeout memo's "what's load-bearing" voice.** When the closeout memo says "the load-bearing structural-independence anchor is V_α1_cyl numerical k_inf-exactness", the Sphinx section should reproduce that framing verbatim ("This is the load-bearing L1 structural-independence anchor"). Don't dilute the memo's clear hierarchy of evidence by demoting the load-bearing test to "one of several sanity checks". The closeout memo's hierarchy is the V&V map; Sphinx is its presentation layer.

Quality self-assessment for the Phase 1/2/3 rich-narrative expansion (peierls_greens.rst, 5 stub sections + 2 cross-cutting sections, +1418 lines):
- Derivation depth: 5/5 (cylinder phase-space foundation, slab 2-bounce-per-period structural distinction, ERR-035 (1−αe^{-τ})(1+αe^{-τ}) algebraic identity, hollow-sphere 4-case partition + composability proof, annulus byte-equal-shared discussion all derived in full)
- Cross-references: 5/5 (all `:func:`, `:mod:`, `:class:`, `:eq:`, `:ref:` resolve; 6 new equation labels added; family-overview table cross-links every per-geometry section)
- Numerical evidence: 5/5 (cylinder k_inf-exactness 4-row table; cylinder vacuum convergence 5-row table; Phase-2 bit-equality preservation 5-row table; method-of-images convergence 3-row table; slab vacuum 56× improvement 4-row table; annulus convergence 4-row table)
- Failed approaches: 5/5 (full ERR-035 forensic with mechanism + corner-agreement explanation + bug fingerprint + delegation fix architecture + cascade of side-effects; ERR-034 trajectory bug fingerprint with uniform-source-blindness explanation; analogical-generalisation anti-pattern called out as the algebra-of-record teaching point)
- Code traceability: 5/5 (every per-geometry section ends with Source/tests/provenance subsection naming Branch-1 SymPy module + every `derive_*` function + Branch-2 production module + every public solver function + foundation tests + L1 tests + closeout memo + cross-domain frame memo + plan)
- Derivation source: 5/5 (all algebra cross-references the SymPy `derive_*` functions; no hand-derivation introduced; the algebra-of-record discipline is respected throughout — Sphinx prose presents what the SymPy modules already prove)

The session's overall quality score is 5/5 across all six dimensions, achieved by:
- Treating the closeout memos as the algebra-of-record source.
- Preserving the existing sphere V_α1/V_α2/V_α3 sections (already at depth 5).
- Adding cross-cutting Phase-2 unification + V_α2 strengthening sections that name the family-level mathematical objects without overwriting per-geometry detail.
- Catching all 6 stale Key-Facts claims and updating them to reflect the post-Phase-3 family scope.
