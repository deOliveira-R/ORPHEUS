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

- [lessons.md](lessons.md) — 38 lessons (L-001…L-038), read FIRST each dispatch.
  The spine: a page is done when every cross-ref resolves against the LIVE tree,
  every claim's V&V level matches the skill verbatim, every retired symbol leaves
  no dangling ref, and the build's WARNING/ERROR/CRITICAL set is unchanged from
  the `-E` baseline. Per-lesson detail is in lessons.md — do NOT re-summarize here.
  Newest by number (one-line hooks; full detail in lessons.md):
  - **L-038** — AUDITING a "is the terminal docs phase done?" charter (frame-projection P7): the answer is usually EFFECTIVELY-DONE — each earlier phase's doc-pass (P3/P5/P6) landed its slice INTO the eventual capstone, so the last phase was executed piecewise before it was reached. Verify by the page's OWN self-identification (grep the intro/notes for the phase tag — here a note literally titled "What shipped since (P3/P5/P7)" said "This page (P7) is the capstone"), NOT an open plan line (stale tracking artifact; trust git). The plan's task-`#N` ≠ GitHub issue-`#N` (plan's "#46–52" collided with unrelated TH/Kinetics issues; real trackers #268/#226/#281). Per-item method: locate label → READ (judge vs articulation std: rejected alts? structural WHY? honest seam?) → `-E -W` clean → grep-gate cross-doc `:ref:`/`:eq:` (the -W-BLIND class). A documented SEAM (#275 LeastSquaresFrame; anisotropic-order Σ_{s,ℓ}) is the OPPOSITE of a gap; a charter's "the condensation PAGE" correctly delivered as a SECTION of a shared page is DRY, not missing.
  - **L-037** — FLIPPING a "documented-future seam" to LANDED across an existing rich page (#281 P6 frame.rst/sn.rst): the stale-status blast radius is the WHOLE page — grep `blocked|not built|pending|future seam|lands with P6` (brief named 3, grep found 7); a "one remaining not-built discipline is X" bullet must RE-POINT to the still-unbuilt sibling, not just drop X; the wrong (φ→φ*) rule hides in a LITERATURE-TABLE cell (tell = test φ*·1_R against an INDICATOR trial ⟹ bare-φ* `∫φ*Σ/∫φ*`, not bilinear — the fix is the PRODUCT weight φ*⊙φ); wired⟹no-sentinel VERIFIED against the LIVE test's stacked verifies() (not the brief) then fast-theory-scan (label_exists=True/documented=False); grow via `.. include::` the generated fragment + supporting math UNLABELED (0 net `:label:` change) + the ONE verifies-target label byte-identical.
  - **L-036** — GROWING a thin honest-stub chapter to full at campaign close (A6/ch15 adjoint): flip the stub's stale "in flight" status (verify merge, not prose); PRESERVE the landed-earlier section verbatim (its `:label:`s are live wiring-backlog), grow AROUND it; RECONCILE a NEW canonical taxonomy that subsumes ≥2 sibling framings (walk-orientation⊂Euclidean, μ-reversal=continuous's signature) — reconcile, don't contradict; deferred-wire verifies-target labels UN-sentineled (tests wait for the label; audit is the main agent's separate gate) while definitional siblings get `vv-status documented`; correctness catches (no 1/k-fission+q* fusion; KEigenvalue's 1st arg is the RESOLVENT (L+C), not A_loss); Mode-12 EXACTNESS (k blind to factor-order/vector/**G-metric**; DROPS not blind — F†=F 1.488→0.153); xref: only non-`:noindex:` automodule links, solution.py has `:label:` docstrings so plain-text-by-convention is CORRECT.
  - **L-035** — V7 orphan-slice adjudication (WIRE/SENTINEL/GAP): WIRE iff a test's PRIMARY assertion IS the equation vs a structurally-independent ref (sign-flip reds it); 3 SENTINEL shapes (schema/continuous-def/literature tested-under-a-DIFFERENT-label · native-vs-legacy bit-id · code-not-built); ROOT narrative page = all-SENTINEL (verified downstream under method labels); foundation-file resolve per-test (computed→WIRE, invariant→SENTINEL, module-foundation+class-verifies coexist); doc's named catcher beats a stale line-range (spectral test is Mode-12-blind — pin the OBJECT); fast `_scan_theory_equations` self-check, no pytest collection.
  - **L-034** — #231 P2 rebalance CONTRACT-DENSE file classes (machinery/driver/ABC/mesh/contract-operator/ψ½): honest cut ≪ the teaching-operator pilot & that's CORRECT; cut SURFACE differs by class (comments dominate driver/mesh); automodule even `:noindex:` makes the `-E -W` gate live; a rebalance READ is a free staleness-audit.
  - **L-033** — #231 P2 rebalance PILOT: teaching ALREADY TWIN → expect ZERO MOVED (Haiku MOVED-column is ~noise); CONTRACT = "would a file-local modifier err without this?"; docstring-only proven by token-invariance; not-`automodule`'d ⟹ no gate.
  - **L-032** — #304 P10 `:label:` re-namespacing: label follows its heading's ruling; self-description oracle; delimiter-anchored COUNTED replace.
  - **L-026…L-031** — #231 corpus split→de-dup→metadata→surface-taxonomy→`:label:`-backfill→bibtex-migration.

## 2. Active / doc-debt state — git-true

**No owed Sphinx pass on merged work.** Every SN doc campaign this index
once tracked is **MERGED to main** (git-verified); landed milestones live in
the SN theory page's **"Development history"** section
(`docs/theory/discrete_ordinates.rst`), behavioral lessons in `lessons.md`.
The active track is the **#231 doc-architecture corpus** (see §3) — the
`operator_algebra.rst` reframe + Phase-3 splits + Phase-4/5a cleanup; main
agent commits, I stage + gate. Only OPEN SN branch: **#236**
(`feature/sn-spatial-angular-product`, not a main ancestor); no pending #236
doc work.

> Merge-status in memory goes STALE. Reconcile any "resume / in-flight / NOT
> pushed" claim against `git merge-base --is-ancestor <hash> HEAD` before acting.

## 3. Durable reference (reusable doc-architecture)

Each entry is a ONE-LINE pointer; the full recipe lives in the linked `feedback_*.md`.

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
  **SUPERSEDED** by the PG-homogenization entry above (kept only for why-it-was-tried).
- [domain-op + L2-promotion + asymmetry-law](feedback_domain_op_l2_promotion_asymmetry_law.md)
  — 6-part section shape for a domain OPERATION (transform on a `Solution`, not a solver step)
  born from an L2 promotion + data/behavior split; lead with THE preservation identity as the
  verifies-target. Instance: `Solution.homogenize`+`MaterialMesh` (#267).
- [orbit-space terminology sweep](feedback_orbit_space_terminology.md) —
  add-aside-then-bridge-then-sweep for a precise math term replacing a loose code-name.
- [auto-generated tables](feedback_autogen_tables.md) — registry-as-SSOT: `capability_rows()`
  metadata fn + generator tool + `builder-inited` hook. (Also `lessons.md` L8.)
- [audit-then-edit partitions](feedback_audit_partition.md) — the KEEP/RELOCATE/TRIM/REMOVE
  partition-table shape for a read-only doc-cleanup audit.
- [cross-solver unified-law doc architecture](feedback_cross_solver_unified_law_doc_architecture.md)
  — ONE correctness LAW spanning N solver families: canonical derivation on the prototype page
  + SHORT sibling spellings cross-ref'd; retired seam as "dead-by-DESIGN" backed by a
  consistency theorem. Instance: k-estimator #259/#291.
- **Doc-architecture redesign (#231, OPEN):** the standing target for any "modernize a
  theory page" task — template, machine header, prose rebalancing, V&V slices, bibtex
  (spec in the issue; coupled to Nexus semantic search). Phase 1a–1c DONE (L-022/L-023);
  **Phase 2 code-prose rebalancing ACTIVE** — done: P2-A `scattering.py` PILOT
  (teaching-operator, L-033), P2-B `sn/loss_representation/` (machinery, L-034), P2-C
  `sn/solver.py`+`numerics/iteration.py` (driver, L-034), P2-D `numerics/operator.py`
  (ABC, L-034), P2-G `sn/operators/{streaming,boundary}`+`sn/mesh/augmented_mesh`
  (contract-heavy operator + MESH, L-034). Maps in `.claude/plans/phase2_code_prose/`.
  Five file-classes calibrated (teaching-operator=aggressive TWIN-cut; machinery/driver/
  mesh=small, COMMENTS dominate; ABC=leanest; contract-heavy-operator=small). Main agent
  commits; I stage + gate.
