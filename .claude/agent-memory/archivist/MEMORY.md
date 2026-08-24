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

- **CS4b S5 — the field-factory sugar tier retires, doc half** — branch `feature/cs1-energy-space`,
  2026-08-24, NOT committed (main agent commits). 5 hand-edited doc files, +360/−106 (a 6th,
  `verification/matrix.rst`, is the `builder-inited` regeneration of the LANDED test tree —
  `10063 → 10067` collected — not a hand edit; report it, never revert it). Rewrote the
  `indexing_and_layout.rst` allocator surface to the SPACE-primary story (carrier mints table ·
  measured shapes · `angular_trial_space` construct-general/select-narrow · presence-gated ψ½
  composites · a SECOND dated correction note beside the #346 one), re-worded live guidance in
  `foundations/boundary_conditions.rst` + `verification/sn.rst` (incl. a runnable code block that
  was ALSO carrying a stale `TimedFullField(bulk=…)` — the field is `interior=`), re-tensed the
  `cartesian_multid.rst` S3 select-narrow subsection (landed-as-future) and repaired the
  `operator_algebra.rst` layering clause. Every surviving `zeros_on`/`from_mesh`/`from_ndarray`
  mention is a history-``literal`` or a survivor-API role. Gates: `-E -W` EXIT=0, W/E/C/Syntax
  set unchanged (0 ↔ 0); `DEAD TARGETS 0`; my own import probe over the 37 qualified roles on
  added lines = 0 dead; vv violations 0, sentinels 545 unchanged. ⛔ **open follow-up I reported
  and did NOT fix (code is out of scope):** `orpheus/transport/fields/scalar_flux.py:45` is the
  ONE dead-role survivor tree-wide — its module-docstring history line spells
  ``:meth:`from_mesh` / :meth:`from_ndarray``` as ROLES (the sentence is correctly past-tensed and
  even names the S5 retirement; only the markup is wrong, and being UNQUALIFIED it is skipped by
  `check_docstring_xrefs.py` and renders plain text). ⚠ I first wrote this follow-up from memory as
  "two files, present-tense capability lines" and the live read refuted BOTH halves — `_bases.py:7`
  is already correct plain literals in a "Before B.1, …" sentence.
  → [[lessons-L66]]

- **F-0 — the metric truth, doc half (the three-metric contradiction)** — branch
  `feature/cs1-energy-space`, 2026-08-23, NOT committed (main agent commits). 6 doc files,
  +1246/−81: `foundations/frame.rst` (+656-line `frame-parseval-metric` chapter: the
  theorem φ = Gc ⟹ metric = G⁻¹, the induced-not-a-constant principle, declared-vs-measured
  Gram, the slab DENSE refusal, the frame square, a 7-frame residual table, a 3-shield
  "why nothing caught it", + Key Facts bullet + history entry) · `foundations/
  spherical_harmonics.rst` (the metric×adjoint table; `hilbert-adjoint-equals-metric-times-S0`
  corrected IN PLACE with the pre-F-0 form preserved unlabelled; `implements::` re-derived
  7 → 8; `sh-space-metric` reframed + declared 3) · `conventions/normalization.rst` (the
  ledger's adjoint row was the THIRD contradiction — and its `(2ℓ+1)/W` "unification the
  canon misses" IS the Parseval metric) · `verification/error_catalog.rst` (ERR-039 F-0
  chapter — **no new ERR number**, the landed gates already carry `catches("ERR-039")`) ·
  `foundations/operator_algebra.rst` (1 aside). Gates: `-E -W` EXIT=0, warning set unchanged
  (0 ↔ 0); vv violations 0, sentinels 541 → 545; `DEAD TARGETS 0`; my own import probe over
  the 16 roles on added lines = 0 dead; 17 declared `implements::` + 3 `no-implementation`
  verified in the rebuilt graph. ⛔ **open follow-ups I reported and did NOT fix (code/tests
  off limits):** `frame.py:116-119` says the slab off-diagonals are "~0.5 of the
  Cauchy–Schwarz scale" — `[M]` **0.9347** relative to C–S, 0.5774 relative to the max
  diagonal (same wording copied into 2 test docstrings); `spherical_harmonic_space.py`'s
  CLASS docstring still says `inner_product_weights` holds `4π/(2ℓ+1)` (`from_L`'s was
  updated, the class's was not, and the frame-dressed instance IS that class);
  `scratch/probe_f1_parseval*.py` no longer reproduce their own headline post-repair.
  → [[lessons-L65]]

- **CS1 step 5 — the `spaces.rst` seed (campaign 1, "operators born bound")** — branch
  `feature/cs1-energy-space`, 2026-08-20, NOT committed (main agent commits). NEW 1158-line
  `docs/theory/foundations/spaces.rst` (axis taxonomy · counting-measure theorem in the METRIC
  register · the collapse doctrine written DIALECTICALLY with both refuted versions · fences ·
  dev history) + toctree/list-table registration + `orpheus.numerics.axis` automodule + 3
  `cone_violations` sites in `field_algebra.rst`. Gates: `-E -W` EXIT=0, warning set unchanged
  (0 ↔ 0); 141 role targets 0 dead; `DEAD TARGETS: 0`; vv violations 0, sentinels 540 → 541
  (`spaces-axis-product`, the ONE new label). ⛔ **open follow-ups I reported and did NOT fix
  (out of scope, both measured):** `infinite_medium.rst:1153` + its `:1243-1244` code block still
  teach the retired `basis_shape=(ng, 1)` spelling as current (`[M]` the block still RUNS —
  stale description, not broken code); and `orpheus.numerics.space` / `…mesh.material_mesh` are
  `automodule`'d nowhere while `api/homogeneous.rst` is `:noindex:`, so `FunctionSpace.of_axes`,
  `has_coordinate_cone`, `MaterialMesh.bulk_space` and `solve_homogeneous_infinite` render plain
  text corpus-wide. → [[lessons-L64]]

- **CS3 step 5 — the flux-ontology overturn, doc half** — branch `refactor/cone-field-algebra`,
  2026-08-19. `field_algebra.rst` rewritten as the cone chapter (602 → ~1540 lines) + 7 citer
  pages + `coding-elegance` #18 reversed. 32 dead refs → 0; labels 4 → 5; sentinels 539 → 540;
  `-E -W` EXIT=0, warning set unchanged. ⛔ **open follow-ups I reported and did NOT fix
  (all out of my scope):** `_coefficient_role.py` has 8 present-tense-false lines + 1 DEAD
  `:class:` role (`_flux_role.FluxRole`) and `cross_section_field.py:33` says "unlike the
  affine flux" — the step-3 commit fixed that file's TABLE and REFERENCES and left its
  comparative essay; the cone-witness gate's own docstring mis-states its fixture
  (both legs are `Δx·Σ_t = 100`); `cross-domain-frames/reference.md` (A.1 row, §192, §201)
  and `numerical-bug-signatures/SKILL.md` (§479, §488) still teach the retired ontology;
  `boundary_conditions.rst`'s `SNMesh.axis_widths` is dead (pre-existing, unrelated).
  → [[lessons-L63]]
- **ERR-026 history block: 29 roles → 13, 15 dead → 0** — branch
  `docs/err026-history-is-not-a-crossref`, 2026-08-18, `error_catalog.rst` +63/−23.
  ⛔ **open follow-ups I reported and did NOT fix:** the `check_docstring_xrefs.py` `.rst`
  blind spot (one-line `head_role` fix, `docs/` 49/71 → 207/255) and **40 of 100** stale raw
  file paths in the catalogue. → [[lessons-L62]] ⚠ **re-confirmed unlanded 2026-08-19** — the
  gate saw 1 of 32 dead refs on this task for exactly that reason.

- ⏹ **MERGED — verified against git 2026-08-18, THIRD time this list had frozen "awaiting review"
  on landed work.** Their durable record is the lesson; their doc changes are in the tree:
  MD→corpus catalogue port `a79f57aa` [[lessons-L61]] · nexus #82 implementer declarations
  (both pages) `144cdf51` [[lessons-L60]] [[lessons-L59]] · nexus graph-path retirement
  [[lessons-L58]] · #344 loss-kernel-gauge [[lessons-L57]] · Q5.6.4 τ/partition carve
  [[lessons-L54]] · Boundary B3.0–B3.2 [[lessons-L42]] · DSA #2 close-out [[lessons-L39]].
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
