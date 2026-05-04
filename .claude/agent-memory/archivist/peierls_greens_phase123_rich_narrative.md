---
name: Peierls Variant α — Phase 1/2/3 rich-narrative expansion (end-of-session)
description: 2026-05-02 archivist sweep replacing the 5 stub Sphinx sections (Phase 1 cylinder, Phase 2 unification, Phase 3A slab, Phase 3B slab-asym, Phase 3C-1 hollow sphere, Phase 3C-2 annulus) with full theory-page-depth rich narrative; also added family-overview + V_α2-strengthening cross-cutting sections. Sphinx -W clean; 178/178 tests still pass.
type: project
---

# Peierls Variant α theory page — Phase 1/2/3 rich-narrative expansion

**Branch**: `feature/peierls-greens-cylinder` (commit `87a3fed` baseline).
**Date**: 2026-05-02.
**Verdict**: shipped clean. 178/178 tests still pass; Sphinx `-W` build clean.

## Sections expanded

The expansion replaced 5 `.. todo:: Archivist expansion needed.`
stub blocks plus inserted 2 new cross-cutting sections. Page grew
from 2807 → 4225 lines (+1418 net).

| Section | Pre (lines) | Post (lines) | Status |
|---|---:|---:|---|
| Phase 1 cylinder (`peierls-greens-cylinder`) | ~210 | ~270 | full algebra + tables + V_α identity proofs + cross-check disclaimer + Source/tests |
| Phase 2 unification (`peierls-greens-unification`) | NEW | ~170 | new section: cross-domain frame + bit-equality preservation table + why theory-page space |
| V_α2 strengthening (`peierls-greens-valpha2-strengthening`) | NEW | ~85 | new section: 1A/1B/1C states across sphere/cylinder/slab + L11 lesson |
| Phase 3A slab (`peierls-greens-slab`) | ~175 | ~310 | post-ERR-035 delegation forensic + Source/tests |
| Phase 3B slab asymmetric (`peierls-greens-slab-asym`) | ~200 | ~390 | rank-2 monodromy + method-of-images + ERR-035 forensic + Source/tests |
| Phase 3C-1 hollow sphere (`peierls-greens-hollow-sph`) | ~165 | ~330 | impact-parameter partition + 4-case analysis + composability + R_in→0 reduction + V_α identities + Source/tests |
| Phase 3C-2 annulus (`peierls-greens-annulus`) | ~190 | ~365 | partition + 3D chord scaling + byte-equal-shared discussion + family-closeout roadmap + Source/tests |

Plus the page title change (`sphere, α∈[0,1]` → `family`) and the
family-overview note section near the top, and updates to the Key
Facts to remove stale "sphere geometry only / no cylinder" claims.

## Commits made (all signed by Archivist co-author)

To be made — see "Recommended commit message" below.

## Sphinx build status

`sphinx-build -W --keep-going` exits clean. No new warnings.

The 4 pre-existing SN MMS soft-warnings (`sn-mms-spherical-psi`,
`sn-mms-spherical-qext`, `sn-mms-cylindrical-psi`,
`sn-mms-cylindrical-qext`) appear as informational `... — skipping`
messages (not warnings), unrelated to peierls.

## Test status

Full peierls suite: **178 passed in 332s** (5.5 min). Identical
to pre-edit baseline. No source files touched; tests confirm
no accidental regression.

## Honest verdict on the page's overall coherence

**Does the 6-geometry × 2-orbit-space-class family read as a unified
architecture?** Yes, with three structural threads tying it
together that future readers can follow:

1. **The family-overview table at the top** lists all 6 geometries
   × 2 orbit-space classes × closure rank in one place, names the unifying
   axis (`variant_alpha_core`), and forward-references each
   per-geometry section.
2. **The Phase 2 unification section** names the load-bearing
   mathematical object (the boundary-to-boundary scattering
   resolvent `T = (I − S)^{-1}`) and explicitly states that the
   per-geometry V_α1 proofs are specialisations of the same
   chord-cancellation algebra to specific chord algebras.
3. **Each per-geometry section ends with a Source/tests/provenance
   subsection** in identical format, making it trivial to navigate
   between Branch-1 SymPy / Branch-2 production / closeout memo /
   tests for any given prototype.

**Are there still gaps that future readers would stumble on?**
Three minor ones, deliberately not fixed in this sweep because
they require structural decisions outside Archivist scope:

1. **Sphere V_α1/V_α2/V_α3 sections were already deeply developed
   pre-session** (the existing post-archivist sphere section was
   the depth target). I did NOT modify these sections, only
   forward-referenced them. A future archivist sweep could harmonise
   the structural framing — e.g., add a small forward-pointer
   inside the sphere V_α1 section to "the same chord-cancellation
   pattern is the load-bearing structural identity for sphere/
   cylinder/slab/slab-asym/hollow-sphere/annulus per
   :ref:`peierls-greens-unification`". This would tighten the
   narrative thread but is incremental polish, not a gap.

2. **The "Cross-verification matrix (Plan 2 B5)" section is
   sphere-only** in scope, by historical accident. Cross-verification
   matrices for cylinder/slab/slab-asym/hollow-sphere/annulus
   would be valuable additions but require new analytical or
   structurally-independent reference values to compare against.
   The sphere matrix uses Phase 4 N=1 ≡ white_hebert (V_α2 algebraic
   anchor) and PS-1982 — neither has direct analogues for the
   sister geometries. Future cross-verification work needs new L1
   reference solvers, not just archival expansion.

3. **The Key Facts header still has sphere-centric V_α2 / V_α3
   bullets** (lines 185-197). I deliberately preserved them
   because they're load-bearing for the sphere prototype's history
   and the V_α2/V_α3 sphere sections immediately below the Key
   Facts block; rewriting them as "family-level claims" would
   obscure the sphere-specific evidence. The new family-overview
   note + the cross-cutting V_α2 strengthening section together
   carry the family-level framing without overwriting the
   sphere-specific Key Facts.

The page is functionally coherent as a reference for any given
geometry. A reader landing on the cylinder section gets the cylinder
algebra plus pointers to the unification section + family overview;
a reader landing on annulus gets the same structural pattern with
the cylinder axial-correction identity called out explicitly.

## Documentation drift / discrepancies flagged

**No source-code drift.** Every closeout-memo function name and
file path I cross-referenced exists at the path claimed. The
ERR-NNN entries in `error_catalog.md` exist with the IDs claimed
(ERR-034, ERR-035). The cross-domain frame memo exists at the
path claimed.

One minor stale claim was found and FIXED in the Key Facts:
- "The prototype now covers **sphere geometry only** (no cylinder)..."
  → updated to family-overview pointer.
- "What this prototype does NOT cover: ... multi-region cylinder,
  multi-region cylinder fixed-source (Phase 1b). ... Phase 2
  unification with sphere is deferred until both standalone
  reference solvers have captured accuracy floors."
  → updated to current 6-geometry coverage list with explicit
  "Phase 2 unification SHIPPED at commits efbae9c/92d4f10/166b9ae;
  rank-2 framework SHIPPED at commit 7cb5bc6 (Phase 3B)" closure.
- "**Scope as of 2026-05-01**" header → updated to `2026-05-02`
  with reference to the new family-overview section.

## Recommended commit message

```
docs(peierls): expand Variant α theory page — Phase 1/2/3 rich narrative

Replace 5 stub `.. todo:: Archivist expansion needed.` blocks across
the cylinder (Phase 1), slab symmetric (Phase 3A → ERR-035
delegation), slab asymmetric (Phase 3B), hollow sphere (Phase 3C-1),
and annulus (Phase 3C-2) sections with full rich-narrative
expansions. Add cross-cutting Phase-2 unification section
(operator-theoretic frame, bit-equality preservation) and V_α2
strengthening section (sphere/cylinder/slab three-state pattern)
that together unify the 6-geometry × 2-orbit-space-class family on the
shared `variant_alpha_core` resolvent + closure primitives.

Per Cardinal Rule 3 (Sphinx is the LLM's brain): every per-geometry
section now carries full algebra (chord conservation arguments,
phase-space reasoning, V_α1/V_α2/V_α3 identity proofs), numerical
evidence tables (k_inf-exactness, vacuum convergence floors,
method-of-images agreement), design rationale and what-was-tried
narratives (e.g. ERR-035 forensic), and a Source/tests/provenance
subsection cross-referencing the Branch-1 SymPy module, Branch-2
production module, foundation tests, L1 tests, closeout memo, and
cross-domain frame memo.

Updates the family-overview at the top of the page (table of 6
geometries × 2 orbit-space classes × closure rank), the page title (now
"family" rather than "sphere"), and the Key Facts to remove stale
"sphere only / Phase 2 deferred" claims now that Phase 3 is closed.

Sphinx -W build clean. 178/178 peierls tests still pass at
unchanged tolerances.
```

Single commit recommended (one cohesive thematic edit) per the
brief; could be split into 2 (per-geometry sections + cross-cutting
sections + Key-Facts updates) if the user prefers atomic-by-section
commits, but the cross-cutting sections include forward references
into the per-geometry ones and vice versa, so a single commit is
semantically cleanest.

## Self-assessment quality scores (per Archivist Directive 3)

| Dimension | Score | Notes |
|-----------|------:|-------|
| Derivation depth | 5 | Full chord-conservation arguments, 4-case partition analysis, V_α1 chord-cancellation pattern derived for each geometry, ERR-035 forensic with algebraic identity (1−αe^{-τ})(1+αe^{-τ}) explicit |
| Cross-references | 5 | Every `:func:`, `:mod:`, `:class:`, `:eq:`, `:ref:` resolves; new sections cross-link to existing sphere V_α sections and to closeout memos by file path |
| Numerical evidence | 5 | 5 new convergence/agreement tables (cylinder k_inf-exactness, cylinder vacuum floor, slab vacuum 56× improvement, method-of-images convergence, annulus convergence floor, Phase-2 bit-equality preservation) |
| Failed approaches | 5 | ERR-035 heuristic origin + corner-agreement mechanism + bug fingerprint + delegation fix architecture + cascade of side-effects all narrated with the (1−αe^{-τ})(1+αe^{-τ}) algebraic identity making the failure visible; ERR-034 trajectory bug fingerprint also documented |
| Code traceability | 5 | Every per-geometry section has a Source/tests/provenance subsection naming Branch-1 SymPy module, Branch-2 production module, every public solver function, every internal helper, every test file, every closeout memo, with absolute paths |
| Derivation source | 5 | All math derives from existing SymPy `derive_*` functions (which are explicitly named and cross-referenced); no hand-written derivations introduced; the algebra-of-record discipline is respected |

**Overall**: full-depth rich-narrative expansion as Cardinal Rule 3
demands. The weakest dimension would be the "Cross-verification
matrix" coverage gap noted above (sphere-only matrix; sister
geometries don't have analogous L1 reference solvers yet) — but
that's a verification-coverage gap, not an archival gap.

## Lessons learned for `vv-principles` skill

The session repeatedly highlighted three patterns that the
`vv-principles` skill could capture more sharply:

1. **The "analogical generalisation of a closure formula" anti-pattern**
   (ERR-035's mechanism). The skill's anti-pattern list could gain
   an explicit entry: "NEVER analogise a closure formula across
   geometries by parameter substitution — instead derive the new
   geometry's closure from first principles via Branch-1 SymPy on
   a non-uniform source, with structurally-independent cross-check
   against the rank-N first-principles form." This is the discipline
   that prevented Phase 3C-1 / 3C-2 from re-introducing the ERR-035
   class of bugs.

2. **The "uniform-source-test blindness" pattern** (ERR-034's
   mechanism). When a bug only manifests on non-uniform spatial
   profiles, the V_α1-style closed-α-with-constant-source tests are
   structurally blind to it. Method-of-images on asymmetric BC is
   the canonical structurally-independent gate for slab; R_in → 0
   reduction is the canonical gate for hollow geometries. The skill
   could call this out in the "what V_α1 does NOT verify" section.

3. **The "convergence-rate-as-bug-fingerprint" pattern** (Phase-3A
   slab pre-fix). When apparent convergence is much slower than
   theoretical expectations, that's a fingerprint for a sub-leading
   linear bias from a latent bug (ERR-034 in Phase-3A) — NOT a
   quadrature limitation, until proven otherwise. The skill's
   "diagnostic fingerprints" section could pin this with the
   Phase-3A 5e-4 → 9e-6 (56× improvement after fix) as the
   canonical example.

Will dispatch a `feedback_*.md` memory entry summarising these
lessons after the commit.

## Status

Branch state: clean (only `docs/theory/peierls_greens.rst`,
`docs/verification/matrix.rst` (auto-regenerated), and
`.claude/agent-memory/cross-domain-attacker/MEMORY.md` (sub-agent
side-effect from earlier session) modified). 178/178 tests pass.
Sphinx -W clean. Ready to commit.
