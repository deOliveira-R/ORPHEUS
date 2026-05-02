---
name: Label-dense partition relocation (Phase 2b pattern)
description: How to execute 7-commit ~−1455 LoC partition cleanup when many labels are externally referenced and load-bearing for tests; 4-class anchor triage + stale-text-flip during KEEP edits
type: feedback
---

**For label-dense partition cuts (Phase 2b style — 25 sub-steps,
mixed KEEP/TRIM/RELOCATE/DELETE across one H1 section, ~1500 LoC
total), do an anchor-triage pass FIRST and group the cuts so each
commit's stub preservation pattern is uniform.**

**Why:** Task 2026-04-30 Phase 2b peierls_nystrom.rst cleanup
relocated §F.5 / §F.6 / §Specular-multibounce / §Phase-5 to
Issues #100/#112/#119/#122/#123/#132/#133. The partition contained
~24 `:label:` blocks. Naïve "cut each step in its own commit and
hunt warnings" approach produced 25 micro-commits and invited
broken `:ref:`s. Better: 7 grouped commits, each preceded by an
anchor-triage classifying the labels in the to-be-modified region
into 4 classes:

1. **External-cited from outside this partition** — preserve
   in stub OR retain in KEEP-content. Examples:
   `peierls-rank-n-stability` (10+ test verifies markers),
   `peierls-phase5-retreat` (cited from peierls_geometry.py code),
   `peierls-class-b-Jn-canonical`, `peierls-M-rank-2`,
   `peierls-f4-rank-1-gauge-why` (all cited from §F.6).
2. **Cross-referenced within this partition** — preserve in the
   KEEP-content; the citing site is also being kept.
   Examples: `peierls-rank-n-jacobian-derivation` cited from §F.6's
   `peierls-class-b-Jn-canonical` block.
3. **Auto-listed in `verification/matrix.rst` only** — safe to
   drop with the cut content. Examples: `hebert-3-350`,
   `c-in-jacobian`, `peierls-phase5-sanchez-A6`.
4. **Internal references within the cut block** — drop with the
   cut.

The triage tells you up-front which labels are negotiable. Trying
to figure this out commit-by-commit is what produces dangling
:refs.

**How to apply:**

1. **Triage script** before any edit:
   ```bash
   for L in label1 label2 ...; do
     count=$(grep -rn "$L" docs/ tests/ orpheus/ derivations/ 2>/dev/null \
             | grep -v _build \
             | grep -v "verification/matrix.rst" \
             | grep -v "the_doc_being_edited.rst" \
             | wc -l)
     echo "$L external-citation count: $count"
   done
   ```
   Anything > 0 must be preserved.

2. **Stale-text-flip during KEEP edits.** When a section is
   marked KEEP intact, scan it for stale text from prior phases.
   The Class B Hébert "Code reference" said
   `cylinder-1d raises NotImplementedError` — but in the
   immediately preceding §"Cylinder Hébert + Issue #112 Phase C"
   subsection, ``compute_P_ss_cylinder`` and
   ``compute_G_bc_cylinder_3d`` are documented as shipped. KEEP
   with stale-text-flip is more honest than KEEP-as-is:
   ```rst
   - The closure is currently **sphere-only**. ``cylinder-1d``
     raises :class:`NotImplementedError`...
   ```
   should become
   ```rst
   - :func:`compute_P_ss_cylinder` — cylinder analog
     (Issue #112 Phase C resolution).
   ```

3. **Commit-grouping rule.** Group adjacent steps where the cuts
   share an issue destination + similar stub pattern. Phase 2b's
   29 sub-steps grouped naturally into 7 commits:
   - A: §8 sphere/rank-N → #100/#112/#132 (3 sub-steps, 3 stubs)
   - B: §F.5 V-S falsification → #119 (4 sub-steps)
   - C: §F.5 Lambert/Marshak + L19 → #122/#123 (3 sub-steps)
   - D: §F.5 housekeeping → #120/#121 + DELETE (3 sub-steps)
   - E: §F.6 Probe G cascade → #132 (6 sub-steps)
   - F: Class B Hébert + Specular → #132 (6 sub-steps)
   - G: Specular cleanup + Phase 5 → #133 (4 sub-steps)
   ~5-7 commits per phase keeps the diff reviewable while leaving
   the per-commit blast radius small.

4. **Sphinx -W after every group, not every step.** A grouped
   commit's diff fits in a single Sphinx warning scan; the cost
   of re-running between every micro-step is wasted (and triggers
   noise from in-flight half-edits).

5. **Avoid the "section anchor renamed = ALL old refs broken"
   trap.** Section titles change naturally during trims (e.g.
   "Why F.4 works at rank-1 but does not generalise (Direction Q,
   Issue #122)" → "(Issue #122 close-out)"). Section *anchors*
   (`.. _peierls-f4-rank-1-gauge-why:`) MUST remain stable across
   the title change. Sphinx auto-generated anchors from titles
   would silently break — explicit `.. _label:` directives are
   non-negotiable for cross-referenced sections.

**Counter-pattern (avoid):** "I'll fix the dangling :ref:s when
Sphinx -W complains" — only catches anchors that were used in
:ref: form. Anchors used in :eq: form (LaTeX equation references)
or that propagate via `verifies(...)` test markers won't be
caught by Sphinx -W; you discover the breakage at next test run
or next Nexus query. Triage-first finds them before the cut.

**Net stats Phase 2b:** 10478 → 9517 LoC (−961 LoC, target was
−1455 but the plan's line ranges were pre-Phase-2a so under-trim
is acceptable when KEEP-content is correctly identified). 7 new
commits, 0 Sphinx warnings, 0 broken :refs.
