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

## 1. Lessons — a HOT digest over a COLD archive

Same hot/cold split as this index: read the digest always, page the
archive on demand. **Never re-summarize a lesson here or in the digest —
each layer points down, it does not copy up.**

- [lessons.md](lessons.md) — **HOT digest, read FIRST every dispatch**
  (one `Read` fits it). Every lesson as one imperative rule + its
  failure→correction core, clustered into 9 themes: (1) verify against the
  LIVE tree · (2) the build is blind, grep is the gate · (3) a `:label:` is
  a V&V edge · (4) retirement & staleness · (5) page surgery · (6) doc
  SHAPE per event class · (7) V&V vocabulary curation · (8) code-prose
  rebalance · (9) gates & tooling. Every entry carries a `→ L-0NN` pointer
  into the archive.
- [lessons_archive.md](lessons_archive.md) — **COLD, load on demand**
  (~250 KB — never open whole). One `## L-0NN` section per lesson: war
  stories, evidence tables, `file:line` detail, per-instance calibrations.
  Open ONLY the section a `→ L-0NN` pointer sent you to. New lessons are
  appended HERE first, then distilled into the digest.
  (Counts are deliberately NOT quoted here — a frozen number rots, and
  this line claimed "44 lessons / L-001…L-044" while the archive was at
  L-047. `grep -c '^## L-0'` answers it in one command.)

## 2. Active / doc-debt state — git-true

**No owed Sphinx pass on merged work.** Landed SN milestones live in the
**"Development history"** changelog at `docs/theory/methods/sn/history.rst` (the
`docs/theory/discrete_ordinates.rst` path this index once named was SPLIT away;
orphaned July HTML survives in `_build`, so a stale-ref grep must `test -f` the
SOURCE). Active track = the **#231 doc-architecture corpus** (§3); main agent
commits, I stage + gate.

**In flight (ONE line each; the evidence is in the lesson, the state is in git):**

- **The campaign-1 history rows** — branch `docs/sn-history-campaign1-landing`, 2026-08-24,
  COMMITTED `68d265ef`, not pushed. Discharged the held merge-hash contract: 5 rows on
  `methods/sn/history.rst` (+287) for the `feature/cs1-energy-space` ff-merge at **`55bb47b9`**
  — one per THESIS (fields-are-space-elements · the frame's Parseval metric + the mint ·
  the operator base demands its two spaces · kernels/born-bound · the axis layer), each
  `(merged @ 55bb47b9)`. Plus the merge's own blast radius: 3 stale `*(in development)*` cells
  on `spaces`/`field_algebra` → `merged @ 55bb47b9`, and `frame.rst`'s F-1 heading re-dated
  2026-08-23 (git). Gates: `-E -W` EXIT=0, W/E/C/Syntax 0 ↔ 0; `DEAD TARGETS 0`; nexus
  `dead_references` 0; vv violations 0, sentinels 549 unchanged; my own probe over the 36 added
  roles = 0 dead; the new rows render 0 visible backticks / 0 dead roles. ⛔ **open, REPORTED not
  fixed:** RST-cannot-nest-inline-markup is corpus-wide — `[M]` 125 nested runs / 25 pages and
  **104 roles that do not render at all** / 28 pages (`:math:`mu`` ×10, `:math:`tau`` ×6);
  widened onto **#379** with a retitle proposal. `history.rst` itself carries 18 pre-existing
  nested runs I deliberately left. → [[lessons-L68]]

- ⏹ **MERGED `55bb47b9` (Campaign 1: CS1 · CS4a+R · CS4b S1–S7 incl. F-0/F-1) — verified
  2026-08-24.** The doc halves are IN the tree and their milestones are in `history.rst`;
  `spaces.rst` / `field_algebra.rst` / `frame.rst` carry the per-page narratives.
  ⛔ **open follow-ups I reported and did NOT fix — all CODE/SKILLS, not docs:**
  `orpheus/numerics/face_layout.py:363`'s dead role (`FaceField.from_face_arrays` →
  `BoundaryField`'s) · `transport/fields/scalar_flux.py:45`'s unqualified dead `from_mesh` /
  `from_ndarray` roles in a correctly-past-tensed line · `_require_typed_composite`'s
  "mesh identity" docstring over a space-CONTENT body · `AxisSectionOperator`'s docstring citing
  a probe that no longer prints its number · `IntegratedReactionRate`'s class docstring naming
  `volume_measure`/`.mesh` over a body reading the SPACE's axes · `spherical_harmonic_space.py`'s
  class docstring still claiming `inner_product_weights` holds `4π/(2ℓ+1)` ·
  `_coefficient_role.py`'s 8 present-tense-false lines + 1 dead `:class:` ·
  `cross_section_field.py:33` · the retired flux ontology in the `cross-domain-frames` A.1
  worked example + `numerical-bug-signatures` · `infinite_medium.rst:1153`'s retired
  `basis_shape=(ng, 1)` teaching · `orpheus.numerics.space` / `…mesh.material_mesh` are
  `automodule`'d nowhere (#302). → [[lessons-L68]] [[lessons-L67]] [[lessons-L66]]
  [[lessons-L65]] [[lessons-L64]] [[lessons-L63]]

- **ERR-026 history block: 29 roles → 13, 15 dead → 0** — branch
  `docs/err026-history-is-not-a-crossref`, 2026-08-18. ⛔ open, re-confirmed unlanded
  2026-08-24: the `head_role` one-liner (blindness is ROLE-scoped, not `.rst`-scoped — see
  L-067) and **40 of 100** stale raw file paths in the catalogue. → [[lessons-L62]]

- ⏹ **MERGED — re-verified against git 2026-08-24 (FOURTH time this list had frozen "awaiting
  review" on landed work).** Durable record = the lesson; doc changes are in the tree:
  MD→corpus catalogue port `a79f57aa` [[lessons-L61]] · nexus #82 implementer declarations
  `144cdf51` [[lessons-L60]] [[lessons-L59]] · nexus graph-path retirement [[lessons-L58]] ·
  #344 loss-kernel-gauge [[lessons-L57]] · Q5.6.4 τ/partition carve [[lessons-L54]] ·
  Boundary B3.0–B3.2 [[lessons-L42]] · DSA #2 close-out [[lessons-L39]].
  ⚠ `git status --porcelain -- docs` before believing ANY entry above it.

> Merge-status in memory goes STALE. Reconcile any "resume / in-flight / NOT
> pushed" claim against `git merge-base --is-ancestor <hash> HEAD` before acting.

## 3. Durable reference (reusable doc-architecture)

Each entry is a ONE-LINE pointer; the full recipe lives in the linked `feedback_*.md`.

- **Landed-milestone record:** `docs/theory/methods/sn/history.rst`. POINT here instead of
  re-listing campaigns. (This line named the pre-split `discrete_ordinates.rst` until 2026-08-18,
  contradicting §2 four lines up — an index can go stale against ITSELF.)
- **Ontology-overturn rewrite** (a page whose THESIS was refuted): the recipe is
  [[lessons-L63]] — argument-unit not symbol-unit, the 4-way eq-label fate rubric, the
  unlabelled-history-equation trick, the two-sided illegal-states rule. Instance:
  `field_algebra.rst` affine → cone (CS3).
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
- [orbit-space terminology sweep](feedback_orbit_space_terminology.md) — add-aside-then-bridge-then-sweep for a precise math term.
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
