# Archivist — Memory Index

Slim index. Behavioral lessons live in `lessons.md` (read FIRST each
dispatch). The mechanical build-gating / cross-ref / venv-worktree /
close-out-arc procedure lives in `AGENT.md` ("Build-Gating & Cross-Ref
Reality", "Close-Out Narrative Arc"). The V&V vocabulary lives in the
`vv-principles` / `algebra-of-record` skills. This index holds only
(1) the lessons pointer, (2) git-true active/doc-debt state, (3) durable
doc-architecture reference. Campaign play-by-play is retired — its
behavioral lesson is in `lessons.md`; its landed milestones are in the
SN theory page's "Development history" section.

## 1. Lessons (read first)

- [lessons.md](lessons.md) — 27 documentation/Sphinx/knowledge-architecture
  lessons (L-001…L-027), read FIRST each dispatch. The spine: a page is done
  when every cross-ref resolves against the LIVE tree, every claim's V&V level
  matches the skill verbatim, every retired symbol leaves no dangling ref, and
  the build's WARNING/ERROR/CRITICAL set is unchanged from the `-E` baseline.
  The per-lesson catalog + war-stories live in lessons.md — do NOT re-summarize
  them here (this index stays slim). Newest: L-027 — a "RELOCATE section to page
  X" brief that the mandated CLOSE READ reveals is ALREADY-fully-on-X →
  Cardinal-Rule-2 DE-DUPLICATE (replace with a `:doc:` pointer preserving the
  conceptual bridge, merge NOTHING, carry `.. _` aliases onto the canonical
  content), NOT relocate+merge; FLAG the inversion (the brief's "additive parts"
  list was the scoping estimate). NEW `-W`-caught class: `ref.ref` "A title or
  caption not found" fires on a BARE `:ref:` to a `.. _label:` before a
  paragraph (no title/caption) — fix by anchoring before a titled/CAPTIONED
  element (section title / `.. list-table:: Caption`) OR explicit-text
  `:ref:`text <label>``. Reframe-consistency applies to a FOLD too (a fold is a
  move — spell the sub-composite `(L+C)`, result is full `A=L+C−S−B`). Stale
  deferred-issue tags distill by `git merge-base --is-ancestor`, keeping the
  design rationale. L-026 (prior): the SPLIT-to-new-page pattern — locate by
  STABLE title + prove contiguity, grep-inventory traveling labels, slice
  PROGRAMMATICALLY; the f-string-over-authored-LaTeX brace-mangle trap (`-W`
  blind) — don't f-string math-bearing prose; discriminate orphaned `-E` HTML
  by "does the source `.rst` exist?".

## 2. Active / doc-debt state — git-true

**No owed Sphinx pass on merged work.** Every SN documentation campaign
this index once tracked — field-typed operator algebra, LD-on-the-DAG,
the foundation-cleanup cluster, Wave-O role-typing, the matvec carve,
the eigenvalue-verification close-outs, the χ source-index fix — is
**MERGED to main** (git-verified 2026-06-21). Their landed milestones
live in the SN theory page's **"Development history"** section
(`docs/theory/discrete_ordinates.rst:15338`); their behavioral lessons
are in `lessons.md`.

The only OPEN SN branch is **#236** (`feature/sn-spatial-angular-product`,
NOT an ancestor of main). I carry no pending #236 doc work; if a #236
docs slice arrives, follow the stub-expansion / capstone patterns in §3.

> Merge-status in memory goes STALE — a note frozen mid-flight merges in
> a later session. ALWAYS reconcile any "resume / in-flight / NOT pushed"
> claim against `git merge-base --is-ancestor <hash> HEAD` before acting;
> never trust a frozen "NOT merged". (This index was rebuilt on that rule.)

## 3. Durable reference (reusable doc-architecture)

Each entry is a ONE-LINE pointer; the full recipe lives in the linked topic
file (this index stays slim — detail belongs in the `feedback_*.md`).

- **Landed-milestone record:** the SN theory page's "Development history" section
  (`docs/theory/discrete_ordinates.rst`) — architectural-changelog of every landed
  SN refactor. POINT here instead of re-listing merged campaigns.
- [canonical-convention-page](feedback_canonical_convention_page.md) — 13-section
  anatomy for a multi-PR migration's canonical theory page + sweep-audit keep/flip
  rubric. Instance: `index_convention.rst`.
- [canonical-axis-convention SSOT section](feedback_canonical_axis_convention_ssot_section.md)
  — SSOT convention SECTION for a just-landed axis-flip enforced once at a data-ingest
  boundary + the overloaded-"group" terminology note. Instance: energy group-0-fast
  (`cross_section_data.rst`).
- [double-category architecture insight](feedback_double_category_architecture_insight.md)
  — 6-move recipe documenting a deep categorical framing of a SHIPPED type system
  (0-ULP crosscheck = interchange-coherence witness; impossibility = numbered-obstruction
  table). Instance: (Rep×Role) carrier grid, `operator_algebra.rst` (#268/#261).
- [orientation-axis two-frames doc](feedback_orientation_axis_two_frames_doc.md) —
  8-move recipe for a 2×2-face operator unification where ORIENTATION (fwd↔adjoint) is
  the coherence axis, execution a non-free 3rd axis; Euclidean-transpose `.. warning::`;
  swap law `A.H.inverse()≡A.inverse().H` as object identity. Instance: #280 Phase 2.5e
  (`loss_representations.rst`).
- [carrier-grid-typed-seam-layering](feedback_carrier_grid_typed_seam_layering.md) —
  completed NxM typed-carrier grid + typed-seam layering (seam one layer up as a forced
  constraint); completing one fused path silently stales a sibling "shared primitive"
  claim. Instance: Frame P4 scattering grid (`operator_algebra.rst`).
- [capstone-architecture-page](feedback_capstone_architecture_page.md) — shape for a
  NEW capstone page documenting the LAYER above per-method pages (cross-ref-not-duplicate).
  Instance: `loss_representations.rst`.
- [capstone-completion-status-reaudit](feedback_capstone_completion_status_reaudit.md) —
  the COMPLETION phase: re-audit every ship-state STATUS claim vs what merged since phase
  N; document a designed-but-unbuilt sibling as a SEAM (literal not `:class:`); the
  documented-future-seam `.. important::` block. Instance: frame-projection P7 (`frame.rst`).
- [capstone-root-cause-ruling](feedback_capstone_root_cause_ruling.md) — 6-move recipe
  retrofitting a capstone ruling that supplies the structural WHY (a theorem) behind a
  split the docs only ASSERTED; literature as a negative-space table. Instance:
  Funk-Hecke/Schur eigenbasis-ownership (#268, `frame.rst`).
- [operator-classes→frame-faces re-homing](feedback_operator_classes_to_frame_faces_rehoming.md)
  — doc-sweep when a refactor retires standalone M/R operator classes into two FACES of one
  frame; records the #268 reversal (discipline is a TYPE; homog/cond are Petrov-Galerkin).
  Instance: P1 discipline-type carve (#268).
- [operator-reification/retype doc pattern](feedback_operator_reification_retype_doc_pattern.md)
  — 6-move recipe for a #226 carve reifying a duck-typed operator (matrix splitting; block
  coisometry `analysis∘reconstruction=4π·I`, never `=I`/ERR-051). Instance: #226 step 2.
- [named-family-member theory section](feedback_named_family_member_theory_section.md) —
  6-move recipe for a NEW § documenting a named member of an invariant-keyed operator family
  ("name = a promise backed by a test"); iterative-member V&V framing (L-010). Instance:
  #226 step 4 GreenOperator (`operator_algebra.rst`).
- [step-5b first-consumer close-the-loop](feedback_issue_138_step5b_first_consumer_closeloop.md)
  — recipe when a follow-on wires the FIRST production consumer of a verified-but-unwired
  operator type; the spectral-invisibility→object-gate V&V unit (→ vv Mode 12). Instance:
  #226 step 5b (`homogeneous.rst`).
- [consumption-mode + capability-axis](feedback_consumption_mode_and_capability_axis.md) —
  two-page recipe for a NEW consumption mode/capability axis on an operator algebra
  (solve/apply/ASSEMBLE). Instance: stencil-assembly 2b (#272/#284/#282).
- [algebra-of-record stub→narrative](feedback_stub_to_rich_narrative_expansion.md) —
  SymPy-module-as-canonical-source; per-geometry 6-subsection shape; stub/expand separation.
  (Also `lessons.md` L5 + the `algebra-of-record` skill.)
- [solver-replacement campaign close-out](feedback_solver_replacement_campaign_closeout.md)
  — P8 recipe when a legacy "island" solver is replaced by the operator-algebra family
  (theory overhaul + mis-named-law re-attribution + LIVE/MOOT history split). Instance: #290.
- [type-confinement docstring sync](feedback_type_confinement_docstring_sync.md) —
  code-final docstring sync to a carve confining a subtype to one role; keep-vs-flip rubric.
  Instance: P4.5 W-C TimedFullField→FullField.
- [Petrov-Galerkin homogenization reframe](feedback_petrov_galerkin_homogenization_reframe.md)
  — THE LIVE recipe (#268 P3) reframing "flux-weighted homogenization" from forward-only
  "Galerkin in L²(φV)" to honest PETROV-GALERKIN (test=φ·1_R ≠ trial=1_R; L²(φV)-fold is the
  Galerkin DEGENERATE of the φ*≠φ bilinear ⟨φ*,Σφ⟩ case); 9-piece order. Instance:
  `Solution.homogenize` (#268).
- [Galerkin-natural-metric reframe](feedback_galerkin_natural_metric_reframe.md) —
  **SUPERSEDED 2026-06-24** by the entry above (its "Galerkin in L²(φV)" framing was reversed
  to PG). Kept ONLY for why-it-was-tried (the metric-fold IS a real identity for φ*=φ). Do
  NOT follow its "Galerkin" framing.
- [domain-op + L2-promotion + asymmetry-law](feedback_domain_op_l2_promotion_asymmetry_law.md)
  — 6-part section shape for a domain OPERATION (transform on a `Solution`, not a solver step)
  born from an L2 promotion + data/behavior split; lead with THE preservation identity as the
  verifies-target. Instance: `Solution.homogenize`+`MaterialMesh` (#267).
- [orbit-space terminology sweep](feedback_orbit_space_terminology.md) —
  add-aside-then-bridge-then-sweep for introducing a precise math term replacing a loose
  code-name.
- [auto-generated tables](feedback_autogen_tables.md) — registry-as-SSOT: `capability_rows()`
  metadata fn + generator tool + `builder-inited` hook. (Also `lessons.md` L8.)
- [audit-then-edit partitions](feedback_audit_partition.md) — the KEEP/RELOCATE/TRIM/REMOVE
  partition-table shape for a read-only doc-cleanup audit.
- [cross-solver unified-law doc architecture](feedback_cross_solver_unified_law_doc_architecture.md)
  — ONE correctness LAW spanning N solver families: canonical derivation on the prototype page
  + SHORT sibling spellings cross-ref'd; retired seam as "dead-by-DESIGN" backed by a
  consistency theorem. Instance: k-estimator #259/#291.
- **Doc-architecture redesign (#231, OPEN):** theory-page template, machine header, prose
  rebalancing, V&V slices, bibtex. Spec lives in the issue; coupled to Nexus semantic search.
  The standing target for any "modernize a theory page" task. **Phase 1a DONE** (frame.rst
  rename + PG-first reorg + SN homog/cond theory eviction — L-022). **Phase 1b+1c DONE**
  (SN page additive skeleton: nexus-meta machine-header dropdown + Key-Facts/Overview→Synopsis
  fold + §8 Gotchas consolidation + §5/§6 automation-pending stubs + evicted the two
  Investigation-History essays; clean `-W` build; 19788→19408 lines — L-023). NEXT = Phase 1d
  (decompose the ~10k-line Sweep mega-section — CHECKPOINT WITH USER before large deletions),
  then Phase 2 code-prose rebalancing. Not committed by me — main agent commits.
