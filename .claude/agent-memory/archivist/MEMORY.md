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

Same hot/cold split as this index: read the digest always, page the archive on
demand. **Never re-summarize a lesson here or in the digest — each layer points
down, it does not copy up.** Counts are deliberately NOT quoted (a frozen number
rots; `grep -c '^## L-0'` answers it).

- [lessons.md](lessons.md) — **HOT digest, read FIRST every dispatch** (one `Read`
  fits it). Every lesson as one imperative rule + its failure→correction core, in
  9 themes: (1) verify against the LIVE tree · (2) the build is blind, grep is the
  gate · (3) a `:label:` is a V&V edge · (4) retirement & staleness · (5) page
  surgery · (6) doc SHAPE per event class · (7) V&V vocabulary curation · (8)
  code-prose rebalance · (9) gates & tooling. Each entry carries a `→ L-0NN`
  pointer into the archive.
- [lessons_archive.md](lessons_archive.md) — **COLD, load on demand** (~250 KB —
  never open whole). One `## L-0NN` section per lesson: war stories, evidence
  tables, `file:line` detail. Open ONLY the section a pointer sent you to; new
  lessons are appended HERE first, then distilled into the digest.

## 2. Active / doc-debt state — git-true

**No owed Sphinx pass on merged work.** Landed SN milestones live in the "Development
history" changelog at `docs/theory/methods/sn/history.rst` (⚠ NOT the pre-split
`discrete_ordinates.rst`, and orphaned July HTML survives in `_build`, so a stale-ref
grep must `test -f` the SOURCE). Active track = **#231** (§3); main agent commits, I
stage + gate.

**In flight — ONE line each. The evidence is in the lesson; the STATE is in git
(`git status --porcelain -- docs/`), never in this list. All are branch
`fix/angular-phantom-support`, docs UNCOMMITTED (I stage + gate; the main agent commits).**

- **#434 R1 — every question about a group is COMPUTED from its realization** (2026-09-03;
  3 pages +832/−85; sentinels 587→**591**; all gates green; a concurrent rename + a
  thrice-moving test count) → [[lessons-L88]]
- **#429 2.2b — the Γ-slot: a symmetry is asked ON the orbit space** (2026-09-02; 4 pages
  +1318/−107; sentinels 584→587) → [[lessons-L87]]
- **#432 — an orbit space is named by its STABILISER** (2026-09-02; 9 pages +1144/−355;
  sentinels 582→584; one dead xref BY INSTRUCTION) → [[lessons-L86]]
- **#429 fused commit — a 1-D rule's frame binds the basis its ORBIT SPACE admits; ERR-080
  CLOSED** (2026-09-02; 12 pages +1756/−240; sentinels 579→582) → [[lessons-L85]]
- **#429 2.5 — the angular moment space is READ off the frame** (2026-09-02; sentinels
  578→579) → [[lessons-L84]]
- **#429 3.1 — a catalogue entry carries its own ARROW and the measure it pushes forward**
  (2026-09-02; sentinels 577→578) → [[lessons-L83]]
- **#429 2.3 — a map carries its two point sets, so a codomain cannot be forged**
  (2026-09-02; sentinels 576→577) → [[lessons-L82]]
- **#429 2.1b — a basis declares its symmetry by naming what it EATS** (2026-09-01;
  sentinels 576→576) → [[lessons-L81]]
- **#429 2.0a — a quotient carries TWO coordinate systems, the chart and the section**
  (2026-08-31; code landed `b55bba56` mid-session) → [[lessons-L80]]
- ⛔ **ERR-026 history block** — branch `docs/err026-history-is-not-a-crossref`, 2026-08-18,
  still OPEN and re-confirmed unlanded 2026-08-24 (the `head_role` one-liner; 40 of 100 stale
  raw file paths in the catalogue) → [[lessons-L62]]

⚠ **Every entry above is a snapshot.** Reconcile with git FIRST — this list has frozen
"in flight" on landed work SIX times. `git merge-base --is-ancestor <hash> HEAD`, and
`git branch --list <branch>` (a vanished branch means merged).

**⏹ MERGED — the durable record is the lesson + the tree; this is only a pointer.**
Every docs pass from 2026-08 back to the Boundary/DSA work is archived as
[[lessons-L39]]…[[lessons-L79]] (one `## L-0NN` section each, with its commits).
Nothing here needs re-listing: `git log --oneline -- docs/` is the index, and each
lesson names its own hashes. Their CODE-side reports are GitHub's to track; the
corpus-wide RST-nested-markup finding lives on **#379**.

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
