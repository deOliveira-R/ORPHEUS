---
name: Phase 2c staleness sweeps and small-cut relocations
description: Pattern for distributed staleness fixes + audit-table relocations across many small sections; orphan-citation trap when cutting MC sub-subsections; new-issue-before-edit rule for "deferred-work" pointers.
type: feedback
---

Phase 2c (peierls_nystrom.rst, ~−536 LoC across 5 commits) was the small-cut counterpart to Phase 2a/2b. The pattern differs from large relocations:

**Why:** Earlier phases moved 100s-of-LoC contiguous failed-experiment narrative; Phase 2c was 12 small edits scattered across 9000 LoC, each <100 LoC. Repeating the 2a/2b commit-per-step cadence would produce noise. Group adjacent same-theme cuts (e.g. all "Phase B.X forward-looking framing" into one commit) — 5 commits for 12 sub-steps.

**How to apply:**

1. **Staleness-sweep grouping.** Cluster sub-steps by theme + by section neighbourhood. §§9/14/15/12 (all "Phase B.x → present-tense" or "deferred → tracker pointer") share a single commit. §§22.3/22.6/22.7 (all "trim X-LoC of redundant prose / move historical detail elsewhere") share a single commit. §22.9 alone (~245 LoC stub) is its own commit. §29 alone is its own commit.

2. **Forward-looking-framing flips.** When code is shipped (`grep -rn` confirms function/class lives in production), the doc should say so present-tense:
   - "Phase B.2 will deliver" → "ships" or just "."
   - "not yet shipped" → drop entirely
   - "intended skeleton" → "actual code (excerpt)"
   - Code-block caption "Phase B.2, not yet shipped" → ""
   - "available for Phase B.3" → drop the deferred-work framing entirely if there's no tracker (or replace with `Tracked in Issue #NNN`).

3. **Deferred-work paragraphs need a tracker BEFORE the doc edit.** If a §X paragraph reads "future Phase Y will lift Z into module W", and the lift hasn't happened, **file a new issue first**, then edit the doc to point at it. Don't delete the deferred-work narrative — it's load-bearing intent for a future investigator. Example: §12 derive_second_difference deferred lift → filed Issue #141 with full deliverable spec, then trimmed §12 paragraph to a 4-line pointer.

4. **Orphan-citation trap when cutting MC / discipline-specific sub-subsections.** Cutting a sub-subsection sometimes orphans citations that were referenced ONLY from that sub-subsection. After the cut, run `sphinx-build -W` — orphaned `[CitationName]` will surface as `WARNING: Citation [X] is not referenced`. Resolution paths:
   - **Move citation to another section** that references the same topic (e.g. §22.3 MC sub-subsection cut → move [MartinBrown2003], [Leppanen2010] to §23 if §23 covers MC).
   - **Delete the citation entirely** if it's only relevant to the cut sub-subsection AND there's a sibling section that covers the topic with its own citations (§23 already had [Woodcock1965] etc., so deleting the orphans was clean).
   - DO NOT introduce a new `[Sanchez1986]_` citation that doesn't exist in References — `grep -n "^\.\. \[Sanchez" docs/theory/peierls_nystrom.rst` to verify before adding any citation marker.

5. **Audit-table cut to issue-cross-link block.** §22.9 had a 17-row per-primitive landing table + commit-dependency narrative + (a)-(h) acceptance-criterion audit (~245 LoC). All four destination issues (#133/#134/#135/#136) already had verbatim copies in close-out comments. The Sphinx stub that replaces this should:
   - Name the recipes (`chord_quadrature`, `observer_angular_quadrature`, `surface_centred_angular_quadrature`) without enumerating per-primitive landing.
   - Cross-link to ALL relevant close-out issues with explicit comment URLs.
   - Preserve the "intentional residual consumers" subsection — these are LOAD-BEARING for future Phase work (slab geometric immunity, Sanchez verification benchmark, build_volume_kernel non-adaptive). They are not legacy leftovers.
   - Preserve the section anchor (`section-22-9-rollout-outcome`) since §22.7 stub references it.

6. **Empirical-table relocation pattern (§29 V1/V2/V4/V6 scan).** When the table is already in the issue close-out comment, replace the inline table with a **headline-finding paragraph** that names the variants + the headline number(s) the reader needs to know without clicking through. Cross-link to the issue comment URL. Drop observations 1-3 (covered in issue) but keep one teaching observation (Observation 4 — the conceptual entanglement) inline because it's the "what this teaches" load-bearing lesson.

7. **§seealso commit-hash trims.** Same pattern as §9 Phase 4.2: drop "(commits aaa → bbb → ccc)" parentheticals — git log is authoritative. Keep narrative pointers like "Phase B.4" (a phase name, not a hash); drop the actual hash trailer.

8. **Anchor preservation guards (Phase 2c specific).** After every step, `grep -n "^\.\. _" docs/theory/peierls_nystrom.rst | grep <suspect>` to confirm anchor-bearing labels and section anchors all survive. The Phase 2c critical anchors were `theory-peierls-multigroup`, `peierls-rank-n-bc-closure-section`, `peierls-part-iii`, `section-22-7-visibility-cone`, `section-22-9-rollout-outcome`, plus all `:label:` equation labels in §§22.3/22.7/§29 (`peierls-tau-coordinate-transform`, `gauss-legendre-visibility-cone`, `peierls-davison-urho`).

**Quality scores (Phase 2c, 1-5 scale):**
- Derivation depth: 4 — kept all teaching math (§22.3 Step 1/Step 2 derivation, §22.7 Bernstein-ellipse analysis, §29 continuous Marshak derivation).
- Cross-references: 5 — every cut content has explicit comment-URL cross-link.
- Numerical evidence: 4 — kept §22.7 spectral-vs-plain-GL evidence table (production data); relocated §29 V1/V2/V4/V6 tables (failed-experiment).
- Failed approaches: 5 — §29 "What this teaches" + V2-cylinder-rank-3 = 0.45 % proof preserved inline.
- Code traceability: 5 — every relocated section retains `:func:` / `:class:` / `:mod:` cross-refs to current production code.
- Derivation source: 4 — §22.7 derivation tied to `derivations/diagnostics/diag_phase5_round3_visibility_cone_quad.py`; §22.3 has no derivations/ script (the τ-substitution is too elementary to warrant one).

**Final delta:** −536 LoC across 5 commits, exceeding the −415 target. The over-delta came from §22.9 (−233 vs planned −245) plus §16 (−78 vs planned −95) + §29 (−97 vs planned −123) + §§22.3/22.6/22.7 (−123 vs planned −110) + §§9/14/15/12 (−5 vs planned −33). Group A came in well below target because §14's architectural prose ("Option (a) vs Option (b)" decision narrative) was load-bearing teaching content that should NOT be cut despite the staleness-sweep premise.
