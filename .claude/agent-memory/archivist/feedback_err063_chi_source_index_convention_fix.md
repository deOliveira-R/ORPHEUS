---
name: ERR-063 chi source-index convention fix
description: Correcting a DOCUMENTED-WRONG physics convention (Peierls MG fission χ sink→source swap) across theory page + code docstring, with a "Sphinx is the brain" rationale note recording the bug so it can't be re-introduced
type: feedback
---

A WRONG-physics-convention correction (Mode-6 convention drift / Mode-2
variable swap): a numerics fix (ERR-063 / #264 / #257 S10a) re-indexed
the multigroup Peierls fission spectrum χ from the SINK node r_i to the
SOURCE node r_j (`χ shares the fission-point argument with νΣf and φ`;
Hébert 2009 *Applied Reactor Physics* Ch.3 Eq.3.57/3.58). My job: flip
3 doc sites that still showed the sink form + add a rationale note so a
future session can't re-introduce it. Do-NOT-commit (main agent
integrates with the S10a bundle).

**Why / How to apply — the load-bearing moves:**

1. **VERIFY THE CODE, NOT THE PRE-FIX DOCSTRING.** The brief warned the
   pre-fix docstrings carried the bug. Read the canonical fission
   assembly FIRST (`geometry.py:6416/6426` `chi_je = chi_n[j, ge]`;
   `slab.py:313/368` `chi_at_node[j][ge]`) to confirm the CORRECT index
   before editing. The code already had ERR-063 comments naming the
   convention — those are the ground truth, the math docstring was the
   straggler. [[feedback-drift-verification]] generalises: diff against
   the code, not the prose.

2. **error_catalog.md was the RICHEST prose seed** — the ERR-063 entry
   carried the full physics ground (Hébert Eq.3.57/3.58), the
   independent structural witness (sibling Variant-α trajectory_resolvent
   already source-indexes χ via local `chi_nodes·F_r`), the "how it hid"
   story (every fixture region had IDENTICAL χ → χ_i≡χ_j → Δ=0.000
   byte-identical; surfaced only when S10a zeroed a non-fissile region's
   χ → O(1) k_eff move), AND the verification gate name. ALWAYS read the
   ERR entry before writing a bug-history note — it is the algebra-of-record
   for the bug.

3. **The "Sphinx is the brain" rationale note = 4 beats:** (a) state the
   convention (χ at r_j, kernel is sole spatial coupler); (b) WHY (birth
   spectrum fixed at fission point, can't depend on observation point —
   causal argument, not just "Hébert says so"); (c) `.. warning::`
   showing the WRONG sink-indexed form explicitly + how it hid + "do NOT
   re-introduce"; (d) verification gate (the promoted intrinsic-property
   test "non-fissile region's χ must not affect k_eff", mutation-verified
   4.7% move). Show the wrong form INSIDE the warning so r_i survives in
   the page ONLY there — grep `chi.*(r_i)` afterward to confirm it's
   walled to the admonition + the "causally meaningless" prose.

4. **STALE :doc: TARGET TRAP.** Cited `:doc:`...<peierls_greens>`` from
   the error_catalog's phrasing — but `peierls_greens.rst` was RENAMED
   to `trajectory_resolvent.rst`. Cross-doc `:doc:` to a missing target
   renders PLAIN TEXT with NO `-W` warning (Cardinal-Rule-1 staleness).
   `ls docs/theory/ | grep -i green` + grep the page for the live target
   BEFORE trusting a memo's doc name. Verified the repointed href
   actually resolved in the built HTML (`grep href=...trajectory_resolvent`).

5. **[Hebert2009] vs [Hebert2020] both exist** — 2009 defined in
   collision_probability.rst, 2020 (3rd ed.) in peierls_nystrom.rst. The
   brief/catalog specify "Hébert 2009 Eq.3.57/3.58" → cite `[Hebert2009]_`
   (resolves cross-doc, NO duplicate-bib warning since I didn't redefine
   it). Match the citation key to the edition the source names.

6. **New eq-label vv-status.** Minted `peierls-mg-fission-source-local`
   (the local birth-source product) → tagged `.. vv-status:
   peierls-mg-fission-source-local documented` (a representational
   identity, not a solver claim; the verifiable content is χ's
   source-indexing pinned by the intrinsic-property gate). Lands in the
   matrix's Documented-only bucket, off the orphan gate.

7. **`verifies("peierls-mg-operator")` edge was EMPTY in the rebuilt
   graph** — because the promoted test is in the not-yet-committed S10a
   bundle. Expected & correct: the edge draws once the main agent lands
   the test + rebuilds, and it then points at the NOW-corrected equation
   node (degree 68). The whole point of doing the doc fix INSIDE the
   bundle. Report this explicitly so the orchestrator knows the edge is
   pending-on-commit, not broken.

8. **Title underline = len() code points.** "...(χ shares the fission
   node)" with the χ glyph = 66 code points; set the `~` run to exactly
   66 via a python rewrite (NOT eyeballing). An over-long underline is
   tolerated by docutils but I matched it for cleanliness.

**Build gate:** forced `-E -W` gave EXIT=1 with the SOLE warning being
the standing-baseline `mesh.py` `:paramref:` ERROR (out of scope);
non-`-W` `-E` build gave EXIT=0 "build succeeded, 1 warning". Warning
SET unchanged from baseline = acceptance. Staleness still lists
peierls_nystrom (git-timestamp compares last-committed state; clears on
S10a commit — not actionable in an uncommitted working-tree edit).

**Quality self-assessment:**
- Derivation depth: 5/5 (causal birth-spectrum argument + local-source
  eq + sink-form-is-meaningless, not just "Hébert says r_j")
- Cross-references: 5/5 (`:eq:`, `:doc:` repointed+resolved, `:func:`,
  `[Hebert2009]_` cross-doc, new label registered)
- Numerical evidence: 4/5 (Δ=0.000 byte-id, 4.7% mutation move, both
  cited from catalog; no fresh table — the fix IS byte-identical so
  there's nothing to tabulate)
- Failed approaches: 5/5 (the sink-indexed bug IS the failed approach,
  fully documented in the warning with how-it-hid)
- Code traceability: 5/5 (geometry.py/slab.py sites + test gate named)
- Derivation source: 4/5 (physics from Hébert + error_catalog; no SymPy
  derivation script for this — the convention is a one-line indexing
  fact, a script would be ceremony)
