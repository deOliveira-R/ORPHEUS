# SN Documentation Architecture — issue #231 prototype (campaign plan)

**Branch:** `docs/sn-doc-architecture` (off `main`)
**Issue:** #231 — docs(theory): documentation architecture. This plan is the
**execution** of #231's settled design comment, scoped to the SN flagship.
**Started:** 2026-07-13.

This file is the durable, cross-session plan (Cardinal Rule 4). The authoritative
DESIGN is issue #231's "Settled design" comment; this file is the SEQUENCED BUILD +
the three scoping decisions the user made 2026-07-13.

---

## Scoping decisions (user, 2026-07-13)

1. **Entry point = vertical slice on the flagship SN page.** Restructure
   `docs/theory/discrete_ordinates.rst` into the 9-section template AND migrate the
   top SN-file prose + kill task-ID archaeology in one coherent (phased) pass. The
   destination sections and the relocated prose land together.
2. **Extensions = bibtex + sphinx-design NOW.** The two Nexus-independent wins
   (central `refs.bib`; dropdowns/cards for machine header + gotchas). Defer
   `graphviz` + `sphinx-proof` — they are coupled to unshipped Nexus features
   (flow-graph nexus#20, typed-node nexus#21).
3. **Nexus-blocked sections (§5 Implementation map, §6 Verification slice) = STUB +
   TRACK.** A placeholder note in each, citing the blocking Nexus issue; no
   auto-generated content until the feature ships. `tests/_harness/` exists, so a
   later upgrade to a real slice is a follow-up, not a rewrite.

## The 9-section template (target section order)

1. Machine header (`nexus-meta`, rendered collapsed) — AUTHOR the fields Nexus
   can't derive (conventions, invariants, operator glossary, aliases); leave
   equation/entry-point lists for Nexus to derive later. Ingestion waits on nexus#1
   Phase 2 — render now, machine-consume later.
2. Synopsis — dense, named, retrieval-targeted (built from current Key Facts + Overview).
3. Formulation — continuous equations, all `:label:`ed.
4. Discretisation — discrete forms + closures, all labeled.
5. Implementation map — **STUB+track** (auto flow-graph needs nexus#20). `plot_directive`
   figures may be added by hand later; not this pass.
6. Verification — **STUB+track** (auto slice needs Nexus label↔test linking).
7. References — `sphinxcontrib-bibtex`, auto-filtered.
8. Gotchas — consolidated (currently scattered as `~`-subsections); inline boxes
   preferred where a gotcha scopes cleanly to one part (micro-decision default).
9. History — ONE collapsed changelog. `date — decision — commit — issue`. Narrative
   essays relocate to issues (L10).

---

## Audit baseline (explorer, 2026-07-13) — where the work is

**`discrete_ordinates.rst`: 21,881 lines.** Has Key Facts (L12) + Development history
(L20765). Monolith: the `Sweep Algorithm` section (L2534→12572) is ~10k lines = 46% of
the page, itself full of `Wave 2`/`Phase D`/`Phase F`/`S2`/`S3`/`#282` subsection essays.
- Gotchas NOT top-level — scattered (`~` L11835 under Sweep, L15445 under SNSolver).
- THREE history surfaces: two "Investigation History" narrative essays (L18075, L18211)
  + the "Development history" changelog (L20765). Template wants ONE (changelog);
  essays → issues.
- TWO extra sections that belong on OTHER pages: `Spatial homogenization` (L18599) →
  homogenization/galerkin page; `Energy condensation` (L19657) → condensation page.

**Code prose (post-reorg tokenize measurement):** `orpheus/sn` 59.9%, `orpheus/transport`
69.1%. Dominant ABSOLUTE-prose sinks (this is the ROI):
| prose lines | ratio | file |
|---|---|---|
| 2332 | 53.4% | `orpheus/sn/loss_representation/__init__.py` |
| 1787 | 56.9% | `orpheus/sn/solver.py` |
| 1320 | 73.2% | `orpheus/transport/operators/scattering.py` |
| 1125 | 70.1% | `orpheus/sn/sweep/pole_angular_closure.py` |
| 1117 | 82.4% | `orpheus/transport/spatial/scheme.py` |
| 964  | 60.1% | `orpheus/sn/mesh/augmented_mesh.py` |
| 896  | 61.5% | `orpheus/sn/operators/radial_characteristic.py` |
| 834  | 71.8% | `orpheus/sn/operators/streaming.py` |

**Archaeology sweep** (`grep -rnE`, code only): ~527 high-confidence campaign-ID lines
(`Phase [A-Z]` 203, `Wave [A-Z]` 128, `Step [0-9]` 78, `D[0-9][a-z]` 118). Top files:
`sn/solver.py` 93, `sn/loss_representation/__init__.py` 65, `sn/operators/streaming.py` 30,
`sn/sweep/pole_angular_closure.py` 26, `sn/mesh/augmented_mesh.py` 25,
`transport/operators/scattering.py` 25. **`solver.py` + `loss_representation/__init__.py`
are top-2 on BOTH axes — the two dominant targets.** `derivations/` hits are a distinct,
lower-priority cluster (dated/benchmark IDs are more defensible as provenance there).

**Concrete relocation exemplars** (file:line, from the audit):
- MOVES-TO-PAGE: `streaming.py:61-79` (WDD symmetric-closure derivation), `:81-119`
  (BC block-matrix essay — already `:ref:`s bc-extraction, so collapse to pointer),
  `radial_characteristic.py:27-67` (ψ½ characteristic derivation, Hébert §3.9.4),
  `loss_representation/__init__.py:1-96` (loss_action architecture essay).
- MOVES-TO-ISSUE/DELETE: `streaming.py:42-59` (dated D-I/D-J retirement log),
  `solver.py:3-24` ("Wave E Round 2… Wave A-D…" module-docstring narration),
  `pole_angular_closure.py:97-115` (PR-TYPED-6c/Issue-#248 retirement diary),
  `augmented_mesh.py:274-302` (Issue #168 Phase C/D inline comments),
  `radial_characteristic.py:339-343` (step-6 sentence polluting a contract docstring).
- STAYS (contract-grade, healthy): `radial_characteristic.py:352-362` + `:115-135`
  (typed carrier, shapes, Hébert cite), `streaming.py:36-40` (operator forms + typed I/O).
- `:eq:` role is THIN in production: 42 total, only `transport/operators/scattering.py`
  uses it outside `derivations/`. Most cites are free-text "Hébert §3.9.4" — weakly
  machine-linkable. Backfilling `:eq:` where cheap is a Phase-2 sub-goal.

---

## Phased execution (compaction points marked ⏸)

### Phase 0 — Infra foundation (small, low-risk, mergeable on its own)
- [ ] pyproject `docs` deps: add `sphinxcontrib-bibtex>=2.6`, `sphinx-design>=0.6` (DONE in-plan).
- [ ] `pip install -e .[docs]` into `.venv`.
- [ ] `conf.py`: add both extensions; `bibtex_bibfiles = ["refs.bib"]`.
- [ ] `docs/refs.bib`: seed with SN's own citations (Hébert 2009, Lewis & Miller,
      Bailey-Morel-Chang 2010, Morel-Montry, Carlson, Sood registry keys). Source of
      truth is Zotero (manual Better-BibTeX export) — seed by hand now, mechanize later.
- [ ] `docs/theory/glossary.rst`: seed canonical terms (albedo, white/reflective/vacuum
      BC, diamond difference, lethargy, scattering ratio, optical thickness, ordinate,
      sweep, WDD, starting-direction). Wire into a toctree.
- [ ] Sphinx clean `-W` build; Nexus graph rebuilds.
- **⏸ COMPACTION POINT** — commit `chore(docs): adopt bibtex + sphinx-design; seed refs.bib + glossary`.

### Phase 1 — `discrete_ordinates.rst` template restructure (BIG; internally phased)
Order chosen for reversibility + early net-line-reduction wins.
- [ ] **1a — Evict the two EXTRA sections.** Move `Spatial homogenization` (L18599) and
      `Energy condensation` (L19657) to their proper pages (galerkin_projection /
      condensation). Verify no dangling `:ref:` (grep docs/ per retirement-audit search #2).
      ⏸ commit.
- [ ] **1b — Collapse three history surfaces → one.** Relocate the two "Investigation
      History" narrative essays (L18075, L18211) into their originating GitHub issues
      (ERR-025/curvilinear); keep only the one-line "Development history" changelog (L10).
      ⏸ commit.
- [ ] **1c — Impose the template skeleton + Gotchas consolidation.** Introduce the 9
      section headings in order; reorganize existing content under them; author the
      `nexus-meta` machine header; build the Synopsis from Key Facts + Overview;
      consolidate the scattered Gotchas; STUB+track §5 and §6. ⏸ commit.
- [ ] **1d — Decompose the 10k-line Sweep Algorithm mega-section.** Strip the page's OWN
      campaign-ID subsection essays (Wave 2/Phase D/F/S2/S3/#282) to production math +
      relocate narrative → issues. Largest single task; sub-checkpoint internally.
      **CHECKPOINT WITH USER before large deletions.** ⏸ commit(s).

### Phase 2 — Code-prose rebalancing (the "improve code quality" deliverable)
File-by-file in ROI order; each file = one commit; docstring-only, but run the SN test
subset + Sphinx clean after each (guard doctest/`:ref:` breakage). The Phase-1 page now
has the destination sections for relocated derivations.
- [ ] **2a — `sn/solver.py`** (1787 prose / 93 archaeology). Kill campaign narration;
      relocate design-rationale → Implementation-map/Formulation sections; trim to
      contract-grade synopsis. ⏸ commit.
- [ ] **2b — `sn/loss_representation/__init__.py`** (2332 / 65). Relocate the loss_action
      architecture essay → theory; trim. ⏸ commit.
- [ ] **2c — operator-leaf trio**: `streaming.py`, `radial_characteristic.py`,
      `pole_angular_closure.py`. Relocate the exemplar derivations above; kill retirement
      diaries; add `:eq:` roles where a free-text Hébert cite is cheap to upgrade. ⏸ commit(s).
- [ ] **2d — `sn/mesh/augmented_mesh.py`** + spillover. ⏸ commit.
- Rubric (issue #231): STAYS = args/returns/shapes/units/conventions/invariants +
  `:eq:`/ERR one-liner. MOVES-TO-PAGE = derivations/alternatives/teaching. MOVES-TO-ISSUE
  = phase-history narration. COMMENTS = constraint-stating only. Guardrail: contract tier
  stays self-sufficient — no anemic docstrings forcing a query round-trip for basics.

### Phase 3 — SN glossary + equation-label backfill (partial, incremental)
- [ ] Backfill unlabeled `.. math::` blocks on the SN page with module-prefixed labels
      (`sn-…`). Full-corpus backfill (389 blocks) is a separate incremental effort.
- [ ] Expand glossary terms surfaced during Phases 1–2.
- ⏸ commit. Then PR: `Closes #231` is premature (corpus-wide); this PR is the SN prototype
  — reference #231, leave it open for the remaining ~36 pages' mechanical pass.

---

## Guardrails / hazards
- **Correctness (Cardinal Rule 1):** trimming a docstring must NOT delete a load-bearing
  convention/gotcha/invariant. When unsure whether prose is contract-grade, KEEP it or
  relocate — never delete outright. Review every archivist output with session context.
- **Retirement-audit for doc moves:** an unresolved `:ref:`/`:func:`/`:eq:` renders as
  plain text with NO `-W` warning (coding-standards search #2). grep `docs/` for every
  moved anchor/label.
- **`.claude/*` + uncommitted-state hazard (L28):** never `git checkout`/`restore` a path
  with uncommitted work to revert; re-edit instead.
- **Nexus worktree hazard (L22):** this branch is NOT a worktree, so the session's Nexus
  graph tracks it after a `sphinx-build` — but rebuild before trusting structural queries.
- **Delegation:** archivist is the executor for page restructure + code-prose (vv-principles
  + algebra-of-record + nexus preloaded); main agent reviews every output. Delicate
  "is-this-contract-grade" judgment may stay main-agent.

## Status log
- 2026-07-13: plan created; branch cut; pyproject deps added. NEXT = Phase 0 install + conf.py + seeds.
