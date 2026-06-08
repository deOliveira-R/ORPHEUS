# Archivist Memory Index

Universal build-gating, cross-ref reality, venv/worktree, and close-out
arc discipline now live in `.claude/agents/archivist/AGENT.md` — memory
notes carry only task-specific doc-craft. Read `lessons.md` at start.

## Core lessons & rubrics
- [lessons.md](lessons.md) — Sphinx build traps (L1-L9) + content quality rules + per-page quality baseline
- [Verify drift before fixing](feedback_drift_verification.md) — quoted "stale text" may pre-date a fix on disk; diff the committed file before editing
- [Retirement doc sweeps](feedback_retirement_docs.md) — flip "will retire"→"retired in commit X" postmortem; preserve historical WHY, delete orphan-test eq labels, preserve tests-pinned labels

## Wave-O / SN role-typing docs (Issue #208/#201)
- [Wave-O operator-algebra docs](feedback_wave_o_operator_algebra_docs.md) — CONSOLIDATED playbook for the bc-extraction family: sibling-extraction, output role-typing, WavefrontFlux cochain, fold→variadic, composite-source API, affine algebra, eigenvalue+G-adjoint, landed-carve, alias-collision; shared section template + breadcrumb-resolution + per-step kernels
- [SI rate-recovery docs](feedback_si_rate_recovery_docs.md) — Phase 3 boundary-G-S in discrete_ordinates: σ_r-fold trap (#215), ERR-056 shared-face rule, Case-Zweifel-vs-within-group c nuance, honest-scope c-mode=Krylov/DSA
- [SN-internal architectural rewrite](feedback_sn_internal_architectural_rewrite.md) — within-module rewrites closing smoking-gun loops: L7-trap fix subsection, bit-identity-vs-principled-equivalence taxonomy, mesh-time-precompute, §15A invariant-set call-out

## Stub-to-rich-narrative expansion
- [Non-vacuum MMS stub expansion](feedback_nonvacuum_mms_stub_expansion.md) — BC-trace-property novelty (a0>0): source-reading order (closeout→attacker-frame→SymPy→tests), vv-status -psi=documented/-qext=verified, Mode-7 activates/nulls
- [Phase C stub expansion](feedback_phase_c_stub_expansion.md) — 12-part shape for a method-implementer's multi-phase-rewrite stub; deviation-from-plan subsection; verify refs EXIST before citing; :label:(eq) vs :ref:(section)
- [Phase D Carlson-seed narrative](feedback_phase_d_carlson_seed_narrative.md) — 4-memo fix overriding plan architecture (proposed site falsified, diagnostic site confirmed); PARTIAL 3-sub-claim breakdown; duplicate-citation grep gate
- [Stub→rich narrative expansion](feedback_stub_to_rich_narrative_expansion.md) — 3-page sweep doubling LoC from TODO stubs; SymPy derive_* docstrings + ERR closeout memos as prose seeds; cross-doc citation duplicates → plain-text inline
- [Stub-expansion family pages](feedback_stub_expansion_family.md) — 6-subsection per-geometry shape; cross-cutting sections BETWEEN per-geometry; preserve method-implementer :labels; duplicate-eq-label trap
- [Stub-expansion cylinder MR](feedback_stub_expansion_cylinder_mr.md) — 8-todo-stub expansion via closeout-memo + verification-plan + Branch-1 SymPy; per-label math+:vv-status+ancestor/impl :func: links; skills are .md not :ref:; :func: resolves by definition-site

## Peierls / continuous-reference & primitive docs
- [Peierls docs patterns](feedback_peierls_docs.md) — RST pitfalls (csv-table LaTeX commas, title-underline-as-codepoints, citations resolve cross-doc but not inside .. math::), vv-status tagging, quality scores
- [Peierls Greens rich narrative](peierls_greens_phase123_rich_narrative.md) — Phase 1-3 Green's-function rich-narrative archival record
- [Numerics primitive docs (Wave 0)](feedback_numerics_primitive_docs.md) — generic primitives lifted from per-solver code: NEW theory pages + existing-page extensions; DO NOT automodule modules with rich :label: docstrings (duplicate-label) — cross-link instead
- [Orbit-space M/G terminology](feedback_orbit_space_terminology.md) — replace loose "topology family" with precise orbit-space M/G classification across docstrings + Sphinx; ONE definitional aside in landing page

## Issue #196 (storage-layout convention flip)
- [Canonical-convention page](feedback_canonical_convention_page.md) — 13-section index_convention.rst; 6-class sweep-audit rubric; "PR-N era docstring" present-tense flip idiom; :doc: uses actual filename
- [Label reconciliation sweep](feedback_label_reconciliation_sweep.md) — resolve :verifies: missing-label warnings: 3 resolution paths; section-vs-eq label via *-section suffix; in-scope vs out-of-scope by test-file module label
- [Memo correction + Sphinx promotion](feedback_memo_correction_and_promotion.md) — correct stale planning memo + promote its catalog to a new Sphinx H1; banner-at-top; replace_all gotcha; bidirectional memo↔Sphinx cross-link
- [#196 PR-INDEX-6 closeout](issue_196_pr_index_6_closeout.md) — SN principled-storage-flip closeout record
- [#196 PR-CLEANUP-DOCS closeout](issue_196_pr_cleanup_docs_closeout.md) — docs-side cleanup closeout record

## BC trace-law / descriptor cleanup
- [BC trace-law Wave 12 synthesis](feedback_bc_trace_law_wave_12.md) — master page absorbs per-wave closeout memos; anti-pattern catalog + Option-a-vs-b-c retrospective; §16A.5 vacuum-semantic correction; baseline-9-warning gate
- [BC descriptor-cleanup docs](feedback_descriptor_cleanup_docs.md) — retire a same-page interim (Option A → pure-descriptor): replace canonical-section anchor wholesale; cross-ref-rewrite-on-retire; type-system-enforcement framing; :pydata:/:pyattr: not valid roles

(NOTE: a sibling durable note `feedback_post_wave_cleanup_docs.md` — the
post-wave-cleanup close-out pattern, KEPT — lives in the MAIN-checkout
`archivist/` dir, not this worktree; indexed there. Cross-ref'd from
[Wave-O operator-algebra docs](feedback_wave_o_operator_algebra_docs.md).)

## Auto-generated tables & landing pages
- [Auto-generated Sphinx tables](feedback_autogen_tables.md) — registry as single source of truth: metadata-only function + generator + builder-inited hook; label-must-precede-title trap
- [Landing-page split](feedback_landing_page_split.md) — split a flat theory toctree without git mv via index.rst + sister landing pages; existing :doc: refs keep working; landing page IS the taxonomy document

## Standalone build-hygiene
- [Build-hygiene staleness sweep](feedback_build_hygiene_sweep.md) — standalone "make Sphinx clean": module-move autodoc repoint, generated-artifact materialization (.h5 converter / generate_rst = ENV not staleness), malformed-grid-table column-drift, matrix.rst autogen drift, sphinx.ext.todo for stubs

## Audit / relocation / triage (read-only & cross-session)
- [Audit-then-edit partitions](feedback_audit_partition.md) — KEEP/RELOCATE/TRIM/REMOVE table schema for read-only doc-cleanup audits; preserve numerical-evidence tables even when surrounding narrative is a wrong attempt
- [Issue close-out relocation comments](feedback_relocation_comments.md) — 5-block template relocating failed-experiment narrative Sphinx→GitHub; RST→Markdown math/anchor/citation conversions; OPEN-log vs CLOSED-post-mortem framing
- [Large relocation cuts (>900 LoC)](feedback_phase2a_large_relocation.md) — preserve every cross-referenced :label: as a stub anchor; 25-LoC stub structure; rewrite citing sites since tables moved to the issue
- [Label-dense partition (2b)](feedback_phase2b_label_dense_partition.md) — 4-class anchor triage BEFORE the cut; group ~25 sub-steps into 5-7 commits sharing issue destination; flip stale text during KEEP edits
- [Phase 2c staleness sweeps](feedback_phase2c_staleness_sweeps.md) — group ~12 scattered edits into 5 theme-clustered commits; "will deliver"→present-tense flips when grep confirms shipped; orphan-citation trap; file new-issue BEFORE editing deferred pointer
- [Phase 2d table preservation](feedback_phase2d_table_preservation_constraint.md) — preserve numerical tables inside a relocation cut: real LoC delta = section−tables−stub; retitle from debate-narrative to production-evidence
- [Plan-file triage rubric](feedback_plan_triage.md) — DELETE/UPDATE/POST/CONSOLIDATE/SUPERSEDED 5-action framework; archive-not-rm; Cardinal-Rule-2 supersession; cluster-aware comment routing
- [Triage with built-in close-outs](feedback_plan_triage_quadrature_cluster.md) — when umbrella issue carries a §22.9 audit with per-criterion landing rows, default collapses to DELETE for row-mapped plans; fresh issue for genuinely-unshipped tooling
- [Misc-cluster plan triage](feedback_misc_cluster_triage.md) — 4 sub-patterns: RETRACTED banner→DELETE, ERR-cataloged→DELETE, shipped-milestone→DELETE, multi-stage phase→POST close-out; gh-api-PATCH-existing-comment; issue-OPEN-despite-Closes-#NNN recommend-not-close
