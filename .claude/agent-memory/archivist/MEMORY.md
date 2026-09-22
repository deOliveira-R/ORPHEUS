# Archivist — Memory Index

An index, not a memory: one line per entry, each a pointer. Behavioural lessons live in
`lessons.md` (read FIRST each dispatch) over `lessons_archive.md` (cold). Mechanical procedure —
build-gating, cross-ref reality, venv/worktree facts, the 9-step close-out arc — is `AGENT.md`.
V&V vocabulary is the `vv-principles` / `algebra-of-record` skills. Rules are cited by ID, never
copied here.

**Index disciplines**, each written after this file bloated by violating it: (1) a campaign's
rulings live in its TOPIC FILE, never here; (2) no `NEXT = <step>` pointer (`plan-authoring` §6);
(3) merge status comes from git, never from a frozen claim here (`process-discipline`); (4) a
landed campaign is ONE line — name, terminal status, pointer — because its commits are in git,
its lessons in the digest and its open items in GitHub; (5) a hook is ≤ 15 words: enough to decide
relevance, never the content.

## 1. Lessons — a HOT digest over a COLD archive

Same hot/cold split as this index: read the digest always, page the archive on demand. **Never
re-summarise a lesson upward** — each layer points down. Counts are deliberately not quoted
(`grep -c '^- \*\*' lessons.md` answers it).

- [lessons.md](lessons.md) — **HOT digest, read FIRST every dispatch.** Every lesson as one
  imperative plus its failure→correction core, in 9 themes: (1) the LIVE tree is ground truth ·
  (2) the build is blind, grep is the gate · (3) a `:label:` is a V&V edge · (4) retirement &
  staleness · (5) page surgery · (6) doc SHAPE per event class · (7) V&V vocabulary curation ·
  (8) code-prose rebalance · (9) gates & tooling. Each entry carries a `→ L-0NN` pointer.
- [lessons_archive.md](lessons_archive.md) — **COLD, load on demand** (~900 KB — never open
  whole). One `## L-0NN` section per lesson with the war story, the `[M]` numbers and the
  `file:line` detail. Open ONLY the section a pointer names; new lessons append HERE first, then
  distil into the digest. Next free number: **L-115**.

## 2. Landed work — nothing owed

Every archivist pass through 2026-09-18 is **MERGED**. `[M]` 2026-09-21:
`git status --porcelain -- docs/` = **0**; every hash this file used to call "uncommitted" is an
ancestor of `main`; every `refactor/consumers-*` branch is gone; `gh issue view` reads **CLOSED**
for #412 and #231; `nexus errors` reads **0 uncaught** with ERR-085 and ERR-086 both carrying
catchers. There is no owed Sphinx pass and no owed marker.

- **Where the record lives:** landed SN milestones are the "Development history" changelog at
  `docs/theory/methods/sn/history.rst` (⚠ NOT the pre-split `discrete_ordinates.rst`; orphaned
  July HTML survives in `_build`, so a stale-ref grep must `test -f` the SOURCE). Commits are
  `git log --oneline -- docs/`. Per-pass lessons are the digest, §1–§9, over archive sections
  L-095…L-112 for the 2026-09 passes.
- ⛔ **Two claims this section carried until 2026-09-21 were frozen and false, and both are the
  index lying forward** (`process-discipline` "Trust git for merge status"): a list of nine
  "uncommitted on branch X" passes, all merged; and "ERR-085/086 catchers owed by the main
  agent", both landed. Do not re-introduce a per-pass status list here — `git status` answers it
  in one command and cannot go stale.
- ⛔ **The ERR-026 history-block entry is void**: its branch
  `docs/err026-history-is-not-a-crossref` is gone locally and remotely, so the 2026-08-24 "still
  OPEN, unlanded" claim says nothing.
- ⛔ **The xref-gate blindness note moved into the digest** (§2b) — it is a lesson with a
  measurement and two controls, not index state. Its rider is there too: acceptance evidence for
  a page is your OWN import probe with a live AND a retired control.

## 3. Durable reference (reusable doc-architecture)

One-line pointers; the recipe lives in the linked file. A reusable recipe earns a line; a campaign
pass does not.

- **Ontology-overturn rewrite** — archive **L-063**: argument-unit, the 4-way eq-label fate rubric,
  the unlabelled-history-equation trick, the two-sided illegal-states rule.
- [canonical-convention page](feedback_canonical_convention_page.md) — 13-section anatomy + the
  keep/flip rubric for a multi-PR migration.
- [axis-convention SSOT section](feedback_canonical_axis_convention_ssot_section.md) — an axis flip
  enforced at a data-ingest boundary.
- [double-category insight](feedback_double_category_architecture_insight.md) — a categorical
  framing of a shipped type system; impossibility as an obstruction table.
- [orientation-axis two frames](feedback_orientation_axis_two_frames_doc.md) — 2×2-face operator
  unification with ORIENTATION as the coherence axis.
- [carrier-grid typed seam](feedback_carrier_grid_typed_seam_layering.md) — an N×M typed grid, seam
  one layer up; one path completed stales its sibling.
- [capstone architecture page](feedback_capstone_architecture_page.md) — a NEW page for the LAYER
  above per-method pages: cross-ref, never duplicate.
- [capstone completion re-audit](feedback_capstone_completion_status_reaudit.md) — re-audit
  ship-state claims; document an unbuilt sibling as a SEAM.
- [capstone root-cause ruling](feedback_capstone_root_cause_ruling.md) — retrofitting the structural
  WHY behind a split the docs only ASSERTED.
- [operator classes to frame faces](feedback_operator_classes_to_frame_faces_rehoming.md) — the
  sweep when operator classes retire into two FACES of one frame.
- [operator reification / retype](feedback_operator_reification_retype_doc_pattern.md) — reifying a
  duck-typed operator; block coisometry `= 4π·I`, never `= I`.
- [named family member](feedback_named_family_member_theory_section.md) — a NEW section for a named
  member of an invariant-keyed operator family.
- [first-consumer close-the-loop](feedback_issue_138_step5b_first_consumer_closeloop.md) — wiring
  the FIRST consumer of a verified-but-unwired type.
- [consumption mode + capability axis](feedback_consumption_mode_and_capability_axis.md) — a NEW
  consumption mode on an operator algebra (solve / apply / ASSEMBLE).
- [stub to rich narrative](feedback_stub_to_rich_narrative_expansion.md) — the SymPy module as
  canonical source; the stub/expand separation.
- [solver-replacement close-out](feedback_solver_replacement_campaign_closeout.md) — a legacy island
  solver replaced by the operator algebra; the LIVE/MOOT split.
- [type-confinement docstring sync](feedback_type_confinement_docstring_sync.md) — the code-final
  sync when a carve confines a subtype to one role.
- [Petrov-Galerkin reframe](feedback_petrov_galerkin_homogenization_reframe.md) — THE LIVE recipe:
  flux-weighting is a TEST weight, not a measure. Supersedes
  [Galerkin natural metric](feedback_galerkin_natural_metric_reframe.md) (why-it-was-tried only).
- [domain op + L2 promotion](feedback_domain_op_l2_promotion_asymmetry_law.md) — the section shape
  for a domain OPERATION born from an L2 promotion.
- [orbit-space terminology](feedback_orbit_space_terminology.md) — add-aside, bridge, then sweep for
  a precise mathematical term.
- [auto-generated tables](feedback_autogen_tables.md) — registry-as-SSOT: metadata function,
  generator, `builder-inited` hook.
- [audit-then-edit partitions](feedback_audit_partition.md) — the KEEP/RELOCATE/TRIM/REMOVE table
  for a read-only doc-cleanup audit.
- [cross-solver unified law](feedback_cross_solver_unified_law_doc_architecture.md) — ONE law over N
  solver families: canonical derivation plus short sibling spellings.
- **Doc-architecture redesign — #231 is CLOSED** (`[M]` 2026-09-21, `gh issue view 231`). Its
  template, machine header, prose-rebalancing, V&V-slice and bibtex conventions remain the standing
  target for any "modernize a theory page" task; the Phase-2 maps are
  `.claude/plans/phase2_code_prose/` and the five calibrated file-classes are digest §8. File a new
  issue rather than reopening this one.
