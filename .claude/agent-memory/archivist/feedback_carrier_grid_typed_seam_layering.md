---
name: carrier-grid-typed-seam-layering
description: Doc recipe for documenting a completed NxM typed-carrier grid + a typed-seam layering rationale (which layer owns the cast); plus the sibling-section staleness trap when completing one path splits a previously-fused one.
metadata:
  type: feedback
---

Recipe for documenting a **completed typed-carrier grid + the typed-seam
layering rationale** (the "which layer owns the cast" decision). Instance:
Frame-campaign P4 — the `(angular ⊗ moment) × (flux ⊗ source)` scattering
carrier grid + `HarmonicFrame` typed seam, in `operator_algebra.rst` under
the §5.6 kernel material (new H3 `scattering-carrier-grid`).

**Why:** these phases close a 2×2 (or N×M) type grid: N leaves + the edges
between them, where some edges are role-PRESERVING (a change of basis/axis —
the frame faces M/R) and exactly one is role-CHANGING (physics — here Λ,
flux→source). The doc must make the type signatures carry what was buried in
an ndarray chain. Three load-bearing pieces a fresh reader needs:

1. **The grid diagram as a labelled `.. math::`** (array env + `\xrightarrow`/
   `\xleftarrow`/`\downarrow`): N leaf TYPES at the nodes, the edge verbs on
   the arrows. Follow with a `.. list-table::` (one row per leaf, one per
   edge) giving Type + role/axis semantics. Tag the grid eq-label
   `vv-status: <label> documented` — it is a named-field-typing identity, NOT
   a solver claim (lesson L-004); it auto-files under Documented-only.
2. **The layering rationale** — WHICH layer owns the cast between the generic
   primitive faces and the typed carriers, stated as a FORCED constraint not a
   preference: the carriers share their deepest primitive (`Field`) in the LOW
   layer (numerics), but their *castability* (the `mesh` binding +
   `from_mesh*` factories) lives in the HIGH layer (transport `BulkField`), and
   the low layer cannot import the high one → the typed seam (the `IS-A
   GalerkinFrame` wrapper with the typed verbs) MUST live one layer up. Spell
   the three bullets: shared-primitive-low / castability-high / can't-invert-
   layer-order. This is the reusable shape for any "typed wrapper over a
   generic numerics primitive" seam.
3. **Explicit-typed vs fused-composed design choice** — when the SAME math has
   two realizations (a hot fused `np.ndarray` composed-operator path kept as
   the 0-ULP canary, AND an explicit typed edge-by-edge path for a consumer
   that holds a typed iterate), document that the choice is **legibility at
   the call site, not math** — both route the SAME kernel + SAME face, so they
   agree numerically; the fused path keeps the reduction tree single (canary
   meaningful), the typed path makes the role flow read off the signatures.

**THE TRAP (most valuable):** completing ONE path of a previously-FUSED pair
SILENTLY STALES a sibling section. Here P4 gave the windowed arm its own
explicit typed path, which invalidated an existing angular-windowing section's
claim "*both arms call the single `_aniso_source_from_moment_values`
primitive*". **How to apply:** when a phase completes/splits one path of an
operator chain, GREP the sibling sections that describe the SAME chain for now-
stale "shared primitive / single path / both arms call" claims, and update them
with a `.. note::` (what changed, that the old primitive is RETAINED as the
crosscheck oracle not deleted, forward-ref to the new section). The fused-path
helper usually survives as the **0-ULP crosscheck oracle** — say so; do not
imply it was removed. (L-001 verify-vs-live + L-007 update-the-wrong-claim,
specialized to the path-split case.)

**Mechanics that recurred:** API `automodule` for `orpheus.transport.*` is
NEVER added (page convention — surfaced via `:class:` xrefs only; 119 such
cites already on the page); matching that is correct, NOT half-surfacing
(L-002). New module docstrings being `:label:`-free is what *would* make them
automodule-safe IF a transport home existed — it doesn't, so no API change.
The matrix auto-regen also picks up P-phase TEST additions (count drift on
`test_scattering_operator`) — report the real number, it is not a defect
(L-008). Quality self-assessment: derivation-depth 4 (layering argument is
the "derivation"), cross-refs 5, numerical-evidence n/a (structural typing
change — no flux moves; the evidence is the 0-ULP crosscheck, cited not
re-run), failed-approaches 4 (the "Λ-is-not-a-frame-verb" + explicit-vs-fused
rejected-alternative notes), code-traceability 5, derivation-source n/a.
