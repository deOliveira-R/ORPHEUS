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

- [lessons.md](lessons.md) — 42 lessons (L-001…L-042), read FIRST each dispatch.
  The spine: a page is done when every cross-ref resolves against the LIVE tree,
  every claim's V&V level matches the skill verbatim, every retired symbol leaves
  no dangling ref, and the build's WARNING/ERROR/CRITICAL set is unchanged from
  the `-E` baseline. **Per-lesson detail lives ONLY in lessons.md — never
  re-summarize a lesson here.** Topic index of the newest:
  - L-042 — audit the corpus against a just-landed multi-commit REFACTOR
  - L-041 — DOC-ONLY "retire the false promises" B0 pass; prove doc-only by AST
  - L-040 — retiring a per-X flag from the docs (grep the CONCEPT, not the symbol)
  - L-039 — authoring a campaign-CAPSTONE theory page from an algebra-of-record
  - L-038 — auditing an "is the terminal docs phase done?" charter
  - L-037 — flipping a documented-future SEAM to LANDED across a rich page
  - L-036 — growing a thin honest-stub chapter to full at campaign close
  - L-035 — V7 orphan-slice adjudication (WIRE / SENTINEL / GAP)
  - L-033/L-034 — #231 P2 code-prose rebalance: pilot, then contract-dense classes
  - L-032 — `:label:` re-namespacing onto design families
  - L-026…L-031 — #231 corpus split→de-dup→metadata→taxonomy→labels→bibtex

## 2. Active / doc-debt state — git-true

**No owed Sphinx pass on merged work.** Every SN doc campaign this index
once tracked is **MERGED to main** (git-verified); landed milestones live in
the SN theory page's **"Development history"** section
(`docs/theory/discrete_ordinates.rst`), behavioral lessons in `lessons.md`.
The active track is the **#231 doc-architecture corpus** (see §3); main
agent commits, I stage + gate.

**Uncommitted doc work awaiting main-agent review (verify with git before
acting — these claims go stale):**

- **Boundary B3.0–B3.2 doc repair** (branch `refactor/operator-strategy-layers`,
  2026-07-31). 46 claims adjudicated across 6 `docs/theory/**` files; findings
  table at `scratch/b3_doc_repair.md`. Added `bc-domain-narrowing` +
  `bc-method-realizability` / `bc-equivariance` / `bc-refusal-axes` and three
  `documented` eq-labels (matrix auto-regen 521→524). `-E -W` EXIT 0 / 0
  warnings, unchanged from baseline. See [[lessons-L42]].
- **DSA #2 docs close-out** (branch `feature/sn-dsa`, 2026-07-27). New capstone
  `docs/theory/methods/sn/acceleration.rst` + verification/sn.rst gate table +
  13 `sn-dsa-*` labels + 12 `verifies()` markers + refs.bib ×5 (each
  `% AWAITING ZOTERO BACK-PORT`). See [[lessons-L39]].

> Merge-status in memory goes STALE. Reconcile any "resume / in-flight / NOT
> pushed" claim against `git merge-base --is-ancestor <hash> HEAD` before acting.

## 3. Durable reference (reusable doc-architecture)

Each entry is a ONE-LINE pointer; the full recipe lives in the linked `feedback_*.md`.

- **Landed-milestone record:** the SN theory page's "Development history" section
  (`docs/theory/discrete_ordinates.rst`). POINT here instead of re-listing campaigns.
- [canonical-convention-page](feedback_canonical_convention_page.md) — 13-section anatomy
  for a multi-PR migration's canonical theory page + keep/flip rubric (`index_convention.rst`).
- [canonical-axis-convention SSOT section](feedback_canonical_axis_convention_ssot_section.md)
  — SSOT section for an axis-flip enforced at a data-ingest boundary (`cross_section_data.rst`).
- [double-category architecture insight](feedback_double_category_architecture_insight.md) —
  documenting a categorical framing of a SHIPPED type system; impossibility as an
  obstruction table. Instance: (Rep×Role) carrier grid (#268/#261).
- [orientation-axis two-frames doc](feedback_orientation_axis_two_frames_doc.md) — 2×2-face
  operator unification with ORIENTATION as the coherence axis (#280 P2.5e).
- [carrier-grid-typed-seam-layering](feedback_carrier_grid_typed_seam_layering.md) — NxM typed
  grid + seam one layer up; completing one path silently stales a sibling claim (Frame P4).
- [capstone-architecture-page](feedback_capstone_architecture_page.md) — a NEW page for the
  LAYER above per-method pages (cross-ref, don't duplicate). `loss_representations.rst`.
- [capstone-completion-status-reaudit](feedback_capstone_completion_status_reaudit.md) — the
  COMPLETION phase: re-audit ship-state claims; document an unbuilt sibling as a SEAM (P7).
- [capstone-root-cause-ruling](feedback_capstone_root_cause_ruling.md) — retrofitting the
  structural WHY (a theorem) behind a split the docs only ASSERTED (#268 `frame.rst`).
- [operator-classes→frame-faces re-homing](feedback_operator_classes_to_frame_faces_rehoming.md)
  — sweep when standalone operator classes retire into two FACES of one frame (#268 P1).
- [operator-reification/retype doc pattern](feedback_operator_reification_retype_doc_pattern.md)
  — reifying a duck-typed operator; block coisometry `= 4π·I`, never `= I` (#226 step 2).
- [named-family-member theory section](feedback_named_family_member_theory_section.md) — a NEW
  § for a named member of an invariant-keyed operator family (#226 step 4 GreenOperator).
- [step-5b first-consumer close-the-loop](feedback_issue_138_step5b_first_consumer_closeloop.md)
  — wiring the FIRST consumer of a verified-but-unwired type; → vv Mode 12 (#226 step 5b).
- [consumption-mode + capability-axis](feedback_consumption_mode_and_capability_axis.md) — a
  NEW consumption mode on an operator algebra (solve/apply/ASSEMBLE) (#272/#284/#282).
- [algebra-of-record stub→narrative](feedback_stub_to_rich_narrative_expansion.md) —
  SymPy-module-as-canonical-source; stub/expand separation (also lessons L5).
- [solver-replacement campaign close-out](feedback_solver_replacement_campaign_closeout.md) —
  a legacy island solver replaced by the operator-algebra family; LIVE/MOOT split (#290 P8).
- [type-confinement docstring sync](feedback_type_confinement_docstring_sync.md) — code-final
  sync when a carve confines a subtype to one role (P4.5 W-C).
- [Petrov-Galerkin homogenization reframe](feedback_petrov_galerkin_homogenization_reframe.md)
  — THE LIVE recipe: flux-weighting is a TEST weight, not a measure (#268 P3). Supersedes
  [Galerkin-natural-metric](feedback_galerkin_natural_metric_reframe.md) (why-it-was-tried only).
- [domain-op + L2-promotion + asymmetry-law](feedback_domain_op_l2_promotion_asymmetry_law.md)
  — section shape for a domain OPERATION born from an L2 promotion (#267).
- [orbit-space terminology sweep](feedback_orbit_space_terminology.md) —
  add-aside-then-bridge-then-sweep for a precise math term replacing a loose code-name.
- [auto-generated tables](feedback_autogen_tables.md) — registry-as-SSOT: metadata fn +
  generator + `builder-inited` hook (also lessons L8).
- [audit-then-edit partitions](feedback_audit_partition.md) — the KEEP/RELOCATE/TRIM/REMOVE
  partition table for a read-only doc-cleanup audit.
- [cross-solver unified-law doc architecture](feedback_cross_solver_unified_law_doc_architecture.md)
  — ONE law spanning N solver families: canonical derivation + short sibling spellings (#259/#291).
- **Doc-architecture redesign (#231, OPEN):** the standing target for any "modernize a theory
  page" task — template, machine header, prose rebalancing, V&V slices, bibtex (spec in the
  issue). Phase 1a–1c DONE; **Phase 2 code-prose rebalancing ACTIVE** — P2-A/B/C/D/G done
  (maps in `.claude/plans/phase2_code_prose/`). Five file-classes calibrated:
  teaching-operator = aggressive TWIN-cut; machinery/driver/mesh = small, COMMENTS dominate;
  ABC = leanest; contract-heavy-operator = small. Main agent commits; I stage + gate.
