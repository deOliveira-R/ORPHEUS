# Archivist Memory Index

Universal build-gating, cross-ref reality, venv/worktree, and close-out
arc discipline now live in `.claude/agents/archivist/AGENT.md` — memory
notes carry only task-specific doc-craft. Read `lessons.md` at start.

## Core lessons & rubrics
- [lessons.md](lessons.md) — Sphinx build traps (L1-L9) + content quality rules + per-page quality baseline
- [Verify drift before fixing](feedback_drift_verification.md) — quoted "stale text" may pre-date a fix on disk; diff the committed file before editing
- [Retirement doc sweeps](feedback_retirement_docs.md) — flip "will retire"→"retired in commit X" postmortem; preserve historical WHY, delete orphan-test eq labels, preserve tests-pinned labels
- [ERR-058 success-close-out supersedes phase chain](feedback_err058_success_closeout_supersedes_phase_chain.md) — #195: a FIX-FOUND close-out closing a long Phase A-F OPEN loop ("ERR-026 PARTIAL"); ADD close-out section + RETRACTION-tombstone the chase's terminal claims (preserve REASONING, flip terminal claims); top-of-block forward-pointer note covers scan path + deep tombstones; isotropic xfail OFF / aniso re-pinned #195→#229; vv-status DIRECTIVE (not prose) on representational eq-labels or --strict regresses; `-E` rebuild BEFORE baselining audit-strict (stale graph gave false EXIT=0); pre-existing exit-1 = report-not-fix; stale test-PATH sweep alongside claim sweep; docstrings ARE the prose seed

## Wave-O / SN role-typing docs (Issue #208/#201)
- [Wave-O operator-algebra docs](feedback_wave_o_operator_algebra_docs.md) — CONSOLIDATED playbook for the bc-extraction family: sibling-extraction, output role-typing, WavefrontFlux cochain, fold→variadic, composite-source API, affine algebra, eigenvalue+G-adjoint, landed-carve, alias-collision; shared section template + breadcrumb-resolution + per-step kernels
- [SI rate-recovery docs](feedback_si_rate_recovery_docs.md) — Phase 3 boundary-G-S in discrete_ordinates: σ_r-fold trap (#215), ERR-056 shared-face rule, Case-Zweifel-vs-within-group c nuance, honest-scope c-mode=Krylov/DSA
- [SN-internal architectural rewrite](feedback_sn_internal_architectural_rewrite.md) — within-module rewrites closing smoking-gun loops: L7-trap fix subsection, bit-identity-vs-principled-equivalence taxonomy, mesh-time-precompute, §15A invariant-set call-out
- [Sweep-walk-collapse re-layering docs](feedback_sweep_walk_collapse_relayering.md) — S6.4(e) #222 sequel to Wave-2: 4-method direction×storage product → 2 walks × level-op objects × storage-free kernel pair; "Architecture history" note preserving retired-symbol literals; de-role via grep (NOT -W); build-gate proven by filename/symbol grep returning only progress-line matches
- [Type-retirement-but-concept-survives docs](feedback_type_retirement_concept_survives.md) — S6.4(f) #222 sequel: a TYPED FIELD (WavefrontFlux+InteriorFaceSpace) git-rm'd but concept C¹_int survives in 2 realizations (_MovingFrontier front + _octant_face_cochain history); retitle to concept not type, KEEP anchor (cross-doc), prominent 4-para succession note, de-role dead module-path roles to literals (grep gate NOT -W), past-tense type/repoint present-claims to realizations, reconcile now-landed future tense; KEEP derivation/grids/tables as history
- [Axis-primary inversion + 3-D admission docs](feedback_axis_primary_3d_admission_docs.md) — C5 #225 sequel to C4 face-name-carve: constructor data-flow INVERSION (axes→mesh→axes round-trip retired, axes tuple PRIMARY) + axis-native 3-D admission NO Mesh3D; NEW C5 section in boundary_conditions (6-substep arc); lossy-round-trip ASCII before/after; Mode-9 gate-retarget 2-row table (reduced-is-None coincidence-proxy → genuine is_cartesian∧ndim); d=3 value-gate table mapped 1:1 to V&V pillars (k_inf L1-closed-form NEVER MMS + structural-independence note); method-rename stale-ref sweep (from_mesh_and_quadrature→from_quadrature_and_layout) LIVE-vs-HISTORICAL triage + #223-pointer note on out-of-scope worked-example; index_convention rank-generic note (NOT rewrite); incremental-build IS zero-warning (reconciles AGENT.md vs C4's 11)
- [Unified-trace worked-example rewrite](feedback_unified_trace_rewrite.md) — #223 EXECUTES C5's deferred holistic rewrite: modernize boundary_conditions.rst worked example off the pre-#188 split trace (InflowTraceSpace/OutflowTraceSpace, _inflow_trace/_outflow_trace, for_face(inflow_trace=...)) onto one-space-two-selectors TraceSpace; stale-note→rewrite-seed (current up top + ONE historical tombstone bottom); narrate CURRENT path from code w/ WHY-one-space; face=left→xmin ONLY in solver call (INPUT bc_left STAYS); two-subclass→one-type+selectors .. note: (#205/#201); de-role dead :class: in-scope but LEAVE triage-walled historical; L12 acceptance is NON-LIVE not zero-hit (probe_inflow_trace=false-positive); MAIN-checkout baseline=1 NOT worktree's 11

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
- [Face-name carve docs](feedback_face_name_carve_docs.md) — C4 #220 STORAGE-KEYING carve (bc_xmin/bc_left attrs → face-name-keyed bc dict): two-surface bc_* triage (INPUT-mesh KEEP vs SNMesh-resolved MIGRATE); curvilinear bc_right→bc["xmax"]; structural-absence-vs-null 3-col table; latent-d3-bug-closed-by-construction = section + "NO ERR entry (type unconstructible)" note; absorb retired-feature section into carve record + repoint its removed anchor; historical-code-block keep-verbatim+retirement-note; insert-H1-mid-arc re-parents trailing Closure

- [Post-wave architectural-cleanup docs sweep](feedback_post_wave_cleanup_docs.md) — 5-element pattern for a small follow-up that closes a deferral left by a larger refactor wave: CLOSE-OUT arc + motivation preserved + where-they-agree-vs-diverge algebra + Option-X/Y rationale + new eq labels for diverging-semantics decomposition. Sibling of [Wave-O operator-algebra docs](feedback_wave_o_operator_algebra_docs.md).
- [Capstone architecture page](feedback_capstone_architecture_page.md) — S5.5 #222: NEW page documenting the LAYER (loss-representation = N schedules over ONE lower-triangular operator) not the methods; 9-H1 shape; cross-ref-not-duplicate (cell math lives in discrete_ordinates/operator_algebra); measurement+VERBATIM-user-decision admonition = the WHY of the default; bit-id-vs-principled-equiv verification section maps 1:1 to vv; production-docstring staleness OUT of scope; grep-gate :ref:/:eq:; worktree baseline=11 NOT AGENT.md's 1.

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
