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

**In flight (ONE line each; the evidence is in the lesson, the state is in git):**

- **P4-remainder — the producer binds the axis; the courier dies** — 2026-08-29, branch
  `feature/p4rem-producer-binds-axis`, **UNCOMMITTED** (mine; main agent commits). 5 doc
  files + `matrix.rst` (10215 → **10236 = +21**, exactly predicted); CS5's third seam
  DISCHARGED on `spaces.rst` + the new `spaces-generator-route-gate` doctrine section.
  `-E` 0 → 0; xref gate 0 dead (live positive control); `dead_references` 0/52; sentinel 0.
  ⭐⭐ **My own CS5 page shipped a claim already false when it landed** (my 3 reported gaps
  were fixed in `cb3cd15b`, same commit-second as my `f8c69117`). ⭐⭐ The briefed courier
  sweep was `[M]` **1 site**; every real find was un-swept **P4.5–P4b** staleness.
  ⚠ **REPORTED, code-only:** the keystone helper's decoy-refusal attribution is wrong (the
  α-dome refuses only the ROLL; negate/reverse die at the closure's P3 τ-guard).
  → [[lessons-L75]]

- **CS4c step 4 — the fission channel becomes TWO bindings of one datum** — 2026-08-31, branch
  `refactor/cs4c-step4-fission-binding` (HEAD `fadad026`), **UNCOMMITTED** (mine; main agent
  commits). 19 authored `.rst` + the regenerated `matrix.rst`. `-E -W` **0/0/0 both sides**,
  EXIT=0; directive edges 412 → **415** (predicted exactly); documented-sentinel labels
  571 → **574** (3 NEW eq-labels: `sn-gain-channels-one-shape`,
  `sn-gain-transposes-one-shape`, `energy-condensation-fission-dyad`); xref gate 0 dead
  (patched, live positive control); `dead_references` 0/52; harness audit 16 passed/5 xfailed.
  NEW H2 `sn-fission-binding-adjoint` in `adjoint.rst`. ⭐⭐ **A published CODE BLOCK in
  `infinite_medium.rst` was FALSE** (`FissionOperator.from_solver_data` vs the live
  `IsotropicFission.from_material_xs`). ⭐⭐ The corpus's **"same Python classes" thesis** was
  refuted at 3 sites + a machine header; the AST construction-site census (`FissionOperator`
  **1 site, sn only** vs `IsotropicFission` **4, all three packages**) is the repair's evidence.
  ⭐⭐ Two sibling harmonizations described alike differ in KIND: N2N **bit-identical 1000/1000**,
  F **≤5 ULP, 0/200 bit-equal** — a theorem of ℓ=0, not fixture luck.
  ⚠ **REPORTED, code-only:** `n2n.py`'s module docstring understates its own result
  ("a pure IEEE-754 order change") — measured, there is none. ⚠ **DEFERRED:** 37 SN-chapter
  sites still spell the pre-N₂ₙ algebra (step-3 residue) — DECLARED at the chapter root
  rather than swept; a scoped follow-up. → [[lessons-L77]]

- **P7 — the metric becomes an OBJECT; a refusal becomes a capability** — 2026-08-30, branch
  `feature/p7-nondiagonal-metric`, **UNCOMMITTED** (mine; main agent commits). 8 doc files
  (brief named 4 items over 3; my census found 5 more) + the regenerated `matrix.rst`.
  `-E -W` **0/0/0 both sides**, EXIT=0; generated artefacts moved exactly as predicted
  (documented labels 567 → **568**, no-implementation 17 → **18**, tests 10236 → 10266 from the
  code side); `dead_references` 0/52; xref gate 0 dead/988 files; corpus ref/eq/doc 0 dangling.
  New H1 `spaces-metric-object` (9 subsections, 1 new eq-label `spaces-pseudo-inverse-parseval`
  + `no-implementation :kind: identity`); `frame-parseval-dense-refusal` **renamed** →
  `frame-parseval-dense-arm`; ERR-039 gains **chapter 3**. ⭐⭐ Chasing one outlier row turned a
  gate's NAME ("a sphere-family property", `[M]` false twice) into a theorem —
  `M* = R/W ⟺ Y(G⁺ − diag(d)/W) = 0`. ⭐⭐ Three published point-values were one-draw readings
  and were replaced by measured bands / draw-free operator figures. → [[lessons-L76]]

- **(n,2n) isotropy: a physics claim vs a model claim** — 2026-08-31, branch
  `docs/n2n-isotropy-claim`, commit **`6906f2a2`** (mine; not merged/pushed). 11 files
  (4 `.rst` + 5 production docstrings + 2 test-prose). `-E -W` EXIT=0 / 0 warnings both
  sides; xref gate 0 dead (stock AND head-role-patched, live positive control);
  `dead_references` 0/52; `matrix.rst` untouched (predicted: no new eq-label); 53 tests
  pass. New anchor `sn-n2n-p0-truncation`. ⭐⭐ The sort is (a) PHYSICS vs (b) MODEL —
  opposite treatments, and one docstring carried both across a dash. ⭐⭐ Two files
  contradicted THEMSELVES and the **hedge was the true half**. ⭐⭐ A relayed CONTRAST
  (`μ̄ +0.278 vs +0.094 elastic, "~3×"`) summed two DIFFERENT energy windows — over the
  same 50 groups elastic is **+0.4264** and the claim inverts; every other number
  reproduced exactly. ⚠ **REPORTED, not fixed:** `cross_section_data.rst:582-700` says
  (n,2n) is not extracted and the balance is 1-in-1-out — `[M]` false at 3 `file:line`s,
  a DIFFERENT claim class, deliberately not swept. → [[lessons-L78]]

- **ERR-026 history block: 29 roles → 13, 15 dead → 0** — branch
  `docs/err026-history-is-not-a-crossref`, 2026-08-18. ⛔ still open, re-confirmed unlanded
  2026-08-24: the `head_role` one-liner (blindness is ROLE-scoped, not `.rst`-scoped) and
  **40 of 100** stale raw file paths in the catalogue. → [[lessons-L62]]

**⏹ MERGED — collapse to one line each; the durable record is the lesson + the tree.**
⚠ This list has frozen "in flight" on landed work FIVE times; reconcile with `git` first.

- 2026-08-29: **CS5** axis-generator doctrine (`4e7b8977`/`b0bfc06c`/`cb3cd15b` + docs
  `f8c69117`; all 3 reported code gaps repaired at `cb3cd15b`) [[lessons-L74]] ·
  **§5b P0** carrying-prose sweep (`8d093334`) [[lessons-L73]] · **P4.9b** operator poses
  (`9c3eb60a`) [[lessons-L72]] · **P4.9a** closure owns its march (`ca852c44`/`7a0f434c`)
  [[lessons-L71]].
- 2026-08-28 and earlier: α-dome citation retraction + `alpha-dome-recursion` rename
  [[lessons-L70]] · τ-arity + LD-curvilinear-Padé [[lessons-L69]] · campaign-1 history rows
  `68d265ef` [[lessons-L68]] · Campaign 1 `55bb47b9` [[lessons-L68]]…[[lessons-L63]] ·
  MD→corpus catalogue port `a79f57aa` [[lessons-L61]] · nexus #82 declarations `144cdf51`
  [[lessons-L60]] [[lessons-L59]] · nexus graph-path retirement [[lessons-L58]] · #344
  loss-kernel-gauge [[lessons-L57]] · Q5.6.4 τ/partition carve [[lessons-L54]] · Boundary
  B3.0–B3.2 [[lessons-L42]] · DSA #2 close-out [[lessons-L39]].
- Their CODE-side reports are GitHub's to track, not this file's; the corpus-wide
  RST-nested-markup finding lives on **#379**.

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
