---
name: issue-257-consolidated-taxonomy-pass
description: #257 CONSOLIDATED 5-stub archivist pass on operator_algebra.rst (S3b/S5/S6 todos + S7 plumbing note + Pattern-M idiom) — the §5.5-5.7 Operator/Kernel/Functional taxonomy; deferred-for-one-rich-pass pattern; vv-status on structural eq-labels; intra-doc :ref: created same-edit
metadata:
  type: feedback
---

# #257 consolidated taxonomy pass — the deferred-multi-stub rich pass

The campaign-wide doc debt: a method-implementer ran S3b→S8c shipping clean
Sphinx STUBS (`:label:`+`:vv-status:`+`.. todo:: Archivist expansion needed`)
at each stage, DELIBERATELY deferring the rich narrative for ONE consolidated
pass once the taxonomy stabilized. This task WAS that pass. The pattern is the
algebra-of-record skill's stub-vs-narrative separation taken to a campaign
scale: 3 `.. todo::` markers (S3b/S5/S6) + 2 smaller owed items (S7 plumbing
note, the Pattern-M idiom) on ONE page.

## What worked (repeat)

- **Read closeout memos as the prose seed, READ CODE to confirm the spelling.**
  Each stage had a `issue_257_s<N>_*_closeout.md`. They are dense and accurate
  (multiplier-algebra laws, the variance hinges, the kernel-shape rejections),
  but they use LOOSE NAMES in prose. The S3b memo said "make_within_group_operators"
  (a `:func:` I cited) — it DOES NOT EXIST; the real builder is the private
  `_within_group_triple`. `-W` does NOT catch a dead `:func:` (renders plain-text)
  → grep `def <name>` for EVERY symbol a closeout names before citing. The
  page ALREADY cited `_within_group_triple` as `:func:...<orpheus.sn.solver._within_group_triple>`
  (private-func-full-path convention) → matched the established convention.
- **Source-of-truth for an idiom = the CODE + the deciding issue comment.**
  Pattern M (`if TYPE_CHECKING: @overload def apply... else: apply = _apply_impl`)
  is grounded in `fission.py:466-482`/`scattering.py:1338-1358` (read the actual
  blocks) + `gh issue view 261 -c` (the dispatch-spelling A-vs-B decision +
  WHY parked on #261). The plan's S8c entry gave the master-standard rule-2>rule-4
  rationale verbatim. Three sources triangulated.
- **vv-status DIRECTIVE on the 2 NEW derivation eq-labels I minted**
  (`multiplication-operator-embedding`, `fission-as-composition`). These are
  STRUCTURAL/DEFINITIONAL identities (the L∞→B(L²) embedding; the §5.6
  composition reading), NOT solver claims → `.. vv-status: <label> documented`
  with a rationale comment naming the verifiable content (the law-suite / the
  0-ULP bit-identity). They're NOT verifies-targets (grep tests/ = 0) → land in
  the matrix orphan-list, which auto-regenerates (main agent rebuilds matrix).
  This keeps the page's convention (EVERY :label: math block in this campaign
  carries a vv-status) consistent.
- **Intra-doc :ref: created in the SAME edit.** The S6 section's prose flows
  into "see :ref:`heteromorphic-apply-typing`" (the Pattern-M anchor). `-W`
  CATCHES a dangling intra-doc :ref: → I placed the Pattern-M section (with its
  `.. _heteromorphic-apply-typing:` anchor) BEFORE running any build. Confirmed
  resolved (anchor at :281, ref at :1654).
- **Placement by page structure, not by stub order.** Pattern M is about
  TYPING the `apply` primitive when heteromorphic → placed as an H2 (`---`)
  under the Definitions H1 block, right after the apply/solve/apply_transpose
  primitives. The brief said "near the apply primitive definition OR a focused
  subsection" — the Definitions block was the natural home (reads in flow with
  the endomorphism contract it qualifies).

## Page-specific facts (operator_algebra.rst)

- Marker ladder: `===`(H1) / `---`(H2) / `~~~~`(H3) / `^^^^`(H4 — NONE used in
  file yet, but `~~~~` present). My new H3s use `~~~~`; sub-points use plain
  bold-run prose, not H4.
- Bib keys `[LewisMiller1984]`/`[AdamsLarsen2002]`/`[Hebert2009]` are defined
  in OTHER theory pages (discrete_ordinates / collision_probability), NOT here.
  They ARE already cited on THIS page (pre-existing, lines 816/819/849/4529/4571)
  — rendering plain-text by cross-doc convention. I added NO new citations
  (the taxonomy is internal-algebra prose; literature lives in the cited code
  closeouts/issues). Did NOT re-define any `.. [Key]` (would be a dup-cite warning).
- `source_sinks.ScalarSourceSink`/`.AngularSourceSink` are PACKAGE re-exports
  (`source_sinks/__init__.py __all__`) → cite via the package path (page convention).
- `transport.IntegralKernelOperator` / `transport.production_rate_functional.
  ProductionRateFunctional` re-exported from `transport/__init__`.

## Gate

baseline = MAIN-checkout 1 warning (the `mesh.py Mesh1D.from_geometry :paramref:`
ERROR, needs sphinx-paramlinks, out of scope). Forced `-E -b html` to /tmp:
EXIT=0, "build succeeded, 1 warning" — that EXACT paramref ERROR, **zero new**.
All 4 anchors render (`id="..."` grep =1 each); zero `admonition-todo` /
"Archivist expansion needed" remain in the HTML. Label count 66→68 (the 2 new
derivation eq-labels, untagged-but-vv-status'd). Did NOT regenerate the V&V
matrix (main agent batches it).

## Quality self-assessment (1-5)

- Derivation depth: 4 (multiplier-algebra laws as a full list-table; the
  variance argument derived not asserted; the two kernel-shape REJECTIONS with
  the rank-changing reason — but these are TAXONOMY/type-system narratives, not
  PDE derivations, so "full derivation" is bounded by the subject).
- Cross-references: 5 (every operator/field/functional/test linked; private
  funcs via the page's full-path convention).
- Numerical evidence: 3 (the 0-ULP bit-identities + the CAP_SOLVE-NaN
  before/after are cited from closeouts; this is a type-surface campaign, no
  convergence tables to add — the evidence is the law-suite + bit-identity gates).
- Failed approaches: 5 (the literal 3-factor fission composition rejected for
  rank-change; SumOfTensorProductsOperator kernel rejected as a re-derivation;
  Pattern C vs Pattern M; the raising-base NoReturn poison — all narrated WHY).
- Code traceability: 5 (every eq linked to its operator/test).
- Derivation source: 4 (closeout memos + the deciding #261 comment + the code
  blocks; no SymPy module — correct, it's a type-surface carve per the
  closeouts' "NO algebra-of-record Branch-1 manifest owed" posture).

Extends [[feedback-affine-in-sigma-stub-expansion]] / [[feedback-scanmarch-coefficient-model-stub]]
(same #257-family loss_representations.rst stub-expansions; THIS one is the
operator_algebra.rst TAXONOMY consolidation — Operator/Kernel/Functional +
the §5.7 promotion + the first @overload idiom).
