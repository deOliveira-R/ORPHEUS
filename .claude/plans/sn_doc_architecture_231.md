# SN Documentation Architecture — issue #231 prototype (campaign plan)

**Branch:** `docs/sn-doc-architecture` (off `main`)
**Issue:** #231 — docs(theory): documentation architecture. This plan is the
**execution** of #231's settled design comment, scoped to the SN flagship.
**Started:** 2026-07-13.

> ## ⚠ SEQUENCING SUPERSEDED — read `documentation_corpus_architecture.md` FIRST
>
> **2026-07-14.** A whole-corpus architecture proposal now exists:
> **`.claude/plans/documentation_corpus_architecture.md`**. It is the layer ABOVE this plan
> (#231 designed the *page*; that proposal designs the *corpus*), and **it re-scopes this
> plan's remaining phases**. It is the authority on WHAT TO DO NEXT.
>
> - **Phase 1d below is SUPERSEDED** → it becomes corpus **Phase C** ("split the monoliths into
>   chapters"), which strictly contains it: *splitting SN into chapters* ⊃ *decomposing the
>   mega-section in place*. The Haiku fan-out catalog (recorded in 1d) is still the instrument —
>   it now classifies chunks into **target chapters** as well as KEEP/ARCHAEOLOGY/DISTILL/RELOCATE.
> - **USER RULING (2026-07-14): corpus Phase A → B → C.** Phase A (the mechanical skeleton:
>   directories, `git mv`, ref fixes, raise `:maxdepth:`, fix `theory/index`'s one-section defect,
>   exclude `**/*.inc.rst`) and Phase B (de-duplicate by the `:label:` oracle) both come
>   **BEFORE** any monolith split. **Do NOT start 1d/C first.**
> - **Phase 2 (code-prose rebalancing) is unchanged and still valid**, but is downstream of the
>   corpus skeleton for the same reason (relocation needs its destinations to exist).
>
> Phases 0 / 1a / 1b+1c below are **DONE and merged into this branch** — that record stands.

This file is the durable, cross-session plan (Cardinal Rule 4). The authoritative
DESIGN is issue #231's "Settled design" comment; this file is the SEQUENCED BUILD +
the three scoping decisions the user made 2026-07-13.

---

## Scoping decisions (user, 2026-07-13)

1. **Entry point = vertical slice on the flagship SN page.** Restructure
   `docs/theory/methods/sn/index.rst` into the 9-section template AND migrate the
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
- [ ] **1a — Frame page restructure + evict the two SN sections** (user-directed 2026-07-13).
      RENAME `galerkin_projection.rst` → `frame.rst`, retitle "The Discrete Frame — projection
      machinery", reorganize into the pedagogical tree: (1) frame theory general → (2) Petrov-
      Galerkin frame [2a general, 2b **Posing the adjoint** — promote the L482 "not Galerkin-in-
      weighted-metric" note to a subsection: forward reaction-rate folds into L²(φV) but the
      eigenvalue-consistent bilinear ⟨φ*,Σφ⟩ has test=φ*≠integrand=φ so no measure-metric
      reproduces it → discipline on the test side, adjoint falls naturally; + short design-
      rationale, 2c applied-to-homogenization, 2d applied-to-condensation] → (3) Galerkin frame
      [3a general, 3b applied-to-spherical-harmonics, ALL terms named: M/R/S₀/Λ/(2ℓ+1)/4π/Funk–
      Hecke] → (4) advanced retained (ownership/unifying/consumer-table/seam/discipline-as-type/
      evidence/impl-map/history/refs). Code hierarchy `GalerkinFrame(PetrovGalerkinFrame)`
      validates PG-first. MOVE the SN page's `Spatial homogenization` (L18599→19657) + `Energy
      condensation` (L19657→20765) GENERAL theory into 2c/2d; leave a short SN "Consuming the
      frame in SN" subsection + `:doc:`frame`` link. Rename audit: index.rst toctree L75,
      spherical_harmonics.rst:213, operator_algebra.rst:2247 (internal anchors already `frame-*`).
      ⏸ commit.
- [ ] **1b — Collapse three history surfaces → one.** Relocate the two "Investigation
      History" narrative essays (L18075, L18211) into their originating GitHub issues
      (ERR-025/curvilinear); keep only the one-line "Development history" changelog (L10).
      ⏸ commit.
- [ ] **1c — Impose the template skeleton + Gotchas consolidation.** Introduce the 9
      section headings in order; reorganize existing content under them; author the
      `nexus-meta` machine header; build the Synopsis from Key Facts + Overview;
      consolidate the scattered Gotchas; STUB+track §5 and §6. ⏸ commit.
- [ ] **1d — Decompose the 10k-line Sweep Algorithm mega-section + the deferred physical
      re-level.** Largest, riskiest task. Two coupled jobs: (i) strip the Sweep section's OWN
      campaign-ID subsection essays (Wave 2/Phase D/F/S2/S3/#282) to production math + relocate
      narrative → issues; (ii) physically re-level the 23 top sections into template order
      (Scattering→§3, Architecture→§5, etc. — the reorder deferred from 1c), incl. folding the
      in-place `sn-282-gotchas` block into §8 and reconciling the References/Gotchas/History order.

      **EXECUTION APPROACH (user-directed 2026-07-13/14):**
      1. **Compact context FIRST** — before executing 1d, /compact so the main agent has budget to
         genuinely review the archivist's large diff (not rubber-stamp it). This plan file + git log
         are the re-anchor (per compaction-points discipline).
      2. **Haiku fan-out pre-catalog** — before the archivist touches the Sweep section, run a
         cheap-model fan-out (Workflow or parallel `model: haiku` Agents — needs user opt-in) that
         CHUNKS the ~9,360-line Sweep section and CLASSIFIES each chunk/paragraph into buckets:
         KEEP (production math/invariant) · ARCHAEOLOGY→ISSUE (Wave/Phase/Step/#NN narrative) ·
         DISTILL→GOTCHA (current trap) · RELOCATE (→ frame.rst / another page) · with file:line
         spans. Aggregate → a catalog the archivist EXECUTES against (deterministic, reviewable).
         Makes both the archivist's job tractable AND the main-agent review a diff-vs-catalog check.
      3. **Archivist executes** from the catalog; **main agent reviews** the diff against the catalog
         with fresh (post-compact) context; then commit.
      **CHECKPOINT WITH USER before large deletions** (confirm the catalog's KEEP/ISSUE/DISTILL/
      RELOCATE calls before the archivist deletes). ⏸ commit(s).

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
- 2026-07-13: plan created; branch cut; pyproject deps added.
- 2026-07-13: **Phase 0 DONE** @ `8dc3dfca` — bibtex + sphinx-design wired; refs.bib (12 SN
  cites) + glossary (15 terms) seeded; clean -W build.
- 2026-07-13: **Phase 1a DONE** (archivist, reviewed) — `galerkin_projection.rst` → `frame.rst`
  restructured into the frame pedagogy (PG-first; §2b "Posing the adjoint" promoted; §3b all
  terms named); SN homog/condensation general theory (~2166 lines) evicted → §2c/2d; SN reduced
  to the "Consuming the frame in SN" stub. Rename audit clean; independent -W build rc0.
  frame.rst 1521→3706; SN page 21881→19788. NEXT = Phase 1b (collapse the 3 SN history surfaces).
  NOTE: archivist also compacted its OWN agent-memory (MEMORY.md 395→140, +L-022) as a side
  effect — user-approved KEEP, committed @ `c7f74c32`.
- 2026-07-13: **Phase 1b+1c DONE (skeleton scope)** (archivist, reviewed) — machine-header
  nexus-meta dropdown (§1, unregistered-directive-safe code-block form) + Synopsis (§2, replaces
  Key Facts, folds in Overview) + consolidated Gotchas (§8, incl. the distilled ERR-025
  homogeneous-rescale-degeneracy gotcha with verified catcher tests) + §5/§6 automation-pending
  notes + removed the 2 Investigation-History essays (§9, L10; full text preserved in closeout for
  optional issue relocation). NO big-block reorder (deferred to 1d per user). One dangling intra-
  page `:ref:`investigation-err-025`` repointed to the distilled gotcha. Page 19788→19408.
  Independent -W build rc0.
- 2026-07-14: **CORRECTNESS FIXES** (found while grounding the corpus proposal — the corpus was
  asserting a RETIRED operator algebra as current fact in four places):
  `275a753a` the SN machine header + Synopsis (shipped in `afb5571d`) encoded the pre-B-extraction
  fold `L: streaming + boundary` / `(L+C−S−F/k)`; `018ecb7b` `numerics/eigenvalue.py`'s docstring
  claimed CP/diffusion/homogeneous have "no (L,S,F) factorization" (**all three do** — it
  manufactured a wrong taxonomy in-session before the code was read) + made the docstring raw;
  `0ca0d378` retired the "four-operator algebra" spelling corpus-wide (6 sites / 3 pages; B added
  as the fifth leaf). Honest algebra: **A = L + C − S − B**, posed `Aψ = (1/k)Fψ` / `Aψ = q`.
- 2026-07-14: **CORPUS PROPOSAL written** → `.claude/plans/documentation_corpus_architecture.md`
  (`bb56a33a`). Grounded by a 71-page inventory + a verified literature survey + an adversarial
  review of the path-integral root. **Issues filed: #298 / #299 / #300.**

### ⏭ NEXT ACTION (post-compaction, authoritative)

**Read `.claude/plans/documentation_corpus_architecture.md` — §7.1 (what Phase A landed) +
§8 (decisions).** It is the authority on sequencing.

**Corpus Phase A is ✅ DONE @ `08e58ee6`** (2026-07-15): `docs/theory/` is no longer flat —
it is `conventions/ · foundations/ · methods/{sn/} · references/`, every `:doc:` under it is
absolute, and `-W` is clean. **This page's SN content now lives at
`docs/theory/methods/sn/index.rst`.**

**⏭ NEXT = corpus Phase B** (de-duplicate by the `:label:` oracle — move the four mis-filed
blocks to their labelled homes; re-namespace labels as they land). **Then** Phase C, which is
this plan's old Phase 1d re-scoped: split the SN monolith into §5's chapters, driven by the
Haiku fan-out catalog. Because `methods/sn/` already exists and its index is the final URL,
Phase C is **pure addition** — it adds chapters and shrinks `index.rst` into a router; it
moves nothing and breaks no refs. Do NOT start C before B. Phase 2 (code-prose) remains
valid but is downstream.

**Parked, no decision needed to proceed:** the two removed Investigation-History essays (full text
in the Phase-1b/1c archivist closeout; git history preserves them; ERR-025/ERR-026 catalog entries
carry the outcomes). Optional L10-canonical relocation: essay-2 → #95 (clean origin fit); essay-1
has no clean issue home.
