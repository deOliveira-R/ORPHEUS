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

- [lessons.md](lessons.md) — 12 documentation/Sphinx/knowledge-architecture
  lessons. The spine: a page is done when every cross-ref resolves
  against the LIVE tree, every claim's V&V level matches the skill
  verbatim, every retired symbol leaves no dangling ref, and the
  build's WARNING/ERROR/CRITICAL set is unchanged from the `-E`
  baseline. Highlights: verify-against-live-code-not-quoted-prose
  (L1); code-xrefs render plain-text with no warning (L2);
  grep-the-matrix-before-touching-an-eq-label (L3); vv-status DIRECTIVE
  on representational labels (L4); stub-expansion source-reading order
  (L5); cross-doc citations resolve-not-redefine (L6);
  retirement-tombstones-evidence-preserves-WHY (L7); generated
  artifacts never hand-edited (L8); file-local marker ladder +
  code-point underlines (L9); V&V vocabulary curation (L10);
  overloaded-symbol convention sweep — inventory every meaning of the
  letter, verify the target spelling vs LIVE code, replace_all only
  unambiguous strings (L11); re-staged-branch doc-merge into a diverged
  tree — programmatic verbatim splice + pre-reorg path translation +
  place-by-anchor + same-merge forward-ref (deferred→landed) tense-flip
  + literal-vs-dead-xref retirement audit (L12).

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

These hold a reusable doc-design RECIPE, not a campaign log:

- **Landed-milestone record:** the SN theory page's "Development
  history" section (`docs/theory/discrete_ordinates.rst`) — the
  architectural-changelog of every landed SN refactor. POINT here
  instead of re-listing merged campaigns.
- [canonical-convention-page](feedback_canonical_convention_page.md) —
  13-section anatomy for a multi-PR migration's single canonical theory
  page (the documentation-of-record a refactor leaves behind); the
  sweep-audit keep-vs-flip rubric. Instance: `index_convention.rst`.
- [canonical-axis-convention SSOT section](feedback_canonical_axis_convention_ssot_section.md) —
  recipe for an SSOT convention SECTION (not a whole page) documenting a
  just-landed axis-ordering flip enforced once at a foreign-data ingest
  boundary (banner+label → precise statement → ORTHOGONALITY `.. important::`
  to the sibling-axis page → enforcement (`:func:` the reverser, which
  arrays flip/don't) → weighted rationale with the contraction `:label:`
  +vv-status:documented backed by the real L2 permutation gate), PLUS the
  orthogonal-axis TERMINOLOGY note that must ride with it (overloaded word
  "group" = energy-bin vs octant-group, grounded on the live
  `_GaussSeidelResolvent` docstring). Traps: extraction-internals constant
  indexes NATIVE order (flag-not-renumber, verify call order in live code);
  flag-don't-fix the 1-based-vs-0-based satisfied claims; verify the L2
  gate's docstring matches the rationale verbatim (NOT over-claimed as
  solver correctness). Instance: energy group-0-fast (#`refactor/group-0-fastest-convention`,
  `cross_section_data.rst` §`canonical-group-convention`).
- [double-category architecture insight](feedback_double_category_architecture_insight.md) —
  6-move recipe for documenting a deep CONCEPTUAL/categorical framing of an
  ALREADY-SHIPPED type system: grow the seeded section; categorical-part→code→
  consequence list-table; the 0-ULP crosscheck RE-FRAMED as the interchange
  coherence witness; census + PRINCIPLED HOLES; the impossibility argument as a
  numbered-obstruction table (success-story analogue of the close-out's
  structural-obstruction step); realization status via a Development-history
  changelog line (grep the LIVE branch — plan tense lags code). Instance: the
  (Rep × Role) carrier grid as a DOUBLE CATEGORY, `operator_algebra.rst` (#268/#261).
- [carrier-grid-typed-seam-layering](feedback_carrier_grid_typed_seam_layering.md) —
  recipe for a completed NxM typed-carrier grid (N leaves + role-preserving
  vs role-CHANGING edges, grid-as-labelled-`.. math::` + leaf/edge list-table)
  + the typed-seam LAYERING rationale (shared-primitive-low / castability-high
  → seam lives one layer up, as a FORCED constraint) + explicit-typed-vs-fused
  design-choice framing (legibility-not-math). THE TRAP: completing one path
  of a fused pair SILENTLY STALES a sibling "both arms share one primitive"
  claim — grep sibling sections, `.. note::`-update, the fused helper survives
  as the 0-ULP crosscheck oracle (not deleted). Instance: Frame-campaign P4
  scattering carrier grid + `HarmonicFrame` (#refactor/operator-inverse-algebra,
  `operator_algebra.rst` §`scattering-carrier-grid`).
- [capstone-architecture-page](feedback_capstone_architecture_page.md) —
  shape for a NEW capstone page documenting the LAYER above existing
  per-method pages (cross-ref-not-duplicate; measurement+decision
  section; verbatim user-decision admonition). Instance:
  `loss_representations.rst`.
- [capstone-completion-status-reaudit](feedback_capstone_completion_status_reaudit.md) —
  the COMPLETION (final) phase of a multi-phase campaign's capstone: FIRST
  pass = re-audit every ship-state STATUS claim ("X is the only frame
  shipping today", "no concrete Y yet", consumer-table Live/Pending) vs
  what MERGED since the page was written at phase N (grep the consumer
  verbs; walk the table row-by-row) — stale-status correction is the
  load-bearing half, new content the easy half. PLUS: document a
  designed-but-UNBUILT sibling type as a SEAM (literal not `:class:`, tied
  to the real guard — here `GramStructure.DENSE`→`MissingCapability`+#275,
  NOT a `LeastSquaresFrame` class; kill the dangling `:class:` ref, verify
  0 py-class spans in built HTML), AND the success-capstone
  "documented-future seam" `.. important::` status block (what is NOT
  wired / blocked on what / "don't read forward green gates as evidence
  for the lift"; generalises a shipped φ†=1 degenerate; bidirectional
  capstone↔derivation wiring; sharpen "is built later"→"deferred, blocked
  on X"). Instance: frame-projection P7 (`galerkin_projection.rst` —
  `frame-composed-verbs`/`frame-least-squares-discipline`/`frame-adjoint-weighted-seam`).
- [capstone-root-cause-ruling](feedback_capstone_root_cause_ruling.md) —
  6-move recipe for retrofitting a CAPSTONE ruling that supplies the
  structural WHY (a theorem) behind a split the docs only ASSERTED
  axis-by-axis (distinct from the capstone-PAGE entry + the close-out
  arc): bold biconditional → structural leg as the theorem-written-out
  (RE-FRAME existing identities as the spectral resolution, don't
  re-derive) → ownership-fixing asymmetry `.. list-table::` → literature
  corroboration as a NEGATIVE-SPACE table ("N refs, ZERO
  cross-validation") → unifying-principle subsection (axis × symmetry ×
  eigenbasis? × discipline) with the negative consequences → relocation
  tripwire as a forward `.. _label:` section; vv-status:documented on
  every theorem-transcription label; sibling echoes
  cross-ref-not-duplicate. Instance: Funk-Hecke/Schur eigenbasis-
  ownership unifying Galerkin-vs-PG (#268, `galerkin_projection.rst`
  + sh/operator_algebra echoes).
- [operator-classes→frame-faces re-homing](feedback_operator_classes_to_frame_faces_rehoming.md)
  — doc-sweep recipe when a refactor retires standalone projection M +
  reconstruction R operator classes into the two FACES of one discrete
  frame: re-home onto the abstraction + add the frame-theory framing
  (T/T*/S/tight-frame), don't find-replace; KEEP documented-only eq-labels
  by name. NOW ALSO records the #268 REVERSAL — discipline is a TYPE
  (`FrameBase→PetrovGalerkinFrame→GalerkinFrame`), NOT a property/role-marker;
  homogenization/condensation are PETROV-GALERKIN (test=φ/φ*-weighted basis),
  NOT "Galerkin in L²(φV)"; measure carries axis+fixed-L²-metric, never
  discipline. Plus: sibling-page-staleness-FLAG-not-silent-flip (repoint to
  discipline-neutral `FrameBase`, flag the substantive prose reversal as a
  separate task — the `discrete_ordinates.rst` homogenization §); M-not-Π
  total relabel. Instance: P1 discipline-type carve (#268,
  `refactor/operator-inverse-algebra`).
- [operator-reification/retype doc pattern](feedback_operator_reification_retype_doc_pattern.md)
  — 6-move recipe for a #226-taxonomy carve that REIFIES a duck-typed operator
  (confused apply/solve pairing → regular matrix splitting `M=(L+C)−B_lower`;
  output-mode arg that changes codomain → typed composition `P∘A⁻¹`, P a block
  COISOMETRY `analysis∘reconstruction=4π·I` NEVER `=I`/ERR-051, no CAP_SOLVE):
  category-error framing FIRST → reifying math UNLABELLED (foundation gates,
  no verifies-target) → source-subspace HONESTY note (M⁻¹ exact on the
  production {outflow-rows=0} subspace; round-trip STRONGER than falsifier) →
  tombstone the "named methods" design in-place → repoint dead refs +
  rejected-alternatives (whole-face-overwrite dropped y_row; moment-proxy
  RESIDUAL gate CATEGORY-CONFUSED on a non-invertible composition) → numbers +
  mutation gates by name (M-SPLIT-DIR/M-SPLIT-PART) + changelog row. Trap:
  verify the gate predicate vs LIVE `_maybe_window`/`_select_si_resolvent`
  (`reduced is None` was the pre-C5.4 proxy; live = `is_cartesian and ndim==2`).
  Instance: #226 step 2 (`refactor/inverse-as-operator`,
  `si-gauss-seidel-reification`/`windowing-retyped`).
- [named-family-member theory section](feedback_named_family_member_theory_section.md)
  — 6-move recipe for a NEW § documenting a NAMED member of an invariant-keyed
  operator family ("name = a promise backed by a test", §13): labeled
  identity+vv-status → "name is earned" invariant `.. list-table::` (asserts /
  why the generic parent FAILS it) → wraps-the-driver → ordering-ruling edge
  table → promise/"what-failed" (the §18.A ρ/(1−ρ) delta + 2 divergence FP
  shapes = Rule-3 gotchas) → verification. Carries the ITERATIVE-MEMBER V&V
  framing (first iterative member of an all-exact family → NO bit-id twin,
  structural-independence anchors ONLY, NO eigenvalue/MMS, solve≡inverse.apply
  a TAUTOLOGY to exclude; L-010), Euclidean-transpose-vs-`.H` reciprocity, and
  the OPEN→RESOLVED issue tense-correction (L-007 on a DECISION). Instance:
  #226 step 4 GreenOperator (`operator_algebra.rst` §`green-operator` + SN
  Development-history changelog).
- [step-5b first-consumer close-the-loop](feedback_issue_138_step5b_first_consumer_closeloop.md)
  — recipe when a landed follow-on wires the FIRST production consumer of a
  verified-but-unwired operator type: 3-page tense-flip (theory=main rewrite /
  operator consumer-ruling="follow-on #138"→"landed" / API=fix stale caller),
  DON'T over-flip "not wired"→"wired" (spelling YES / factory-routing STILL
  waits + deliberately bypassed); strategy-choice-as-TYPE (explicit
  `MatrixInverseOperator(loss)` vs `loss.inverse()`→iterative GreenOperator);
  shared-extraction-home reframe (old convenience engine `direct_eigenvalue`→
  `(A,F)`-posed SIBLING delegating to `dominant_eigenpair`, 0 consumers, NOT
  fuller-view-oracle); the reusable **spectral-invisibility→object-gate** V&V
  unit (factor-swap `F·A⁻¹`~`A⁻¹F` similarity + transpose `eig(Mᵀ)=eig(M)` →
  |Δk|=0, every k-gate blind → pin the OBJECT `K.as_matrix≡solve(A,F)`;
  vv-status:documented on the identity). Instance: #226 step 5b / task #138
  (`refactor/inverse-as-operator`, `homogeneous.rst` §`spectral-invisibility`).
- [algebra-of-record stub→narrative](feedback_stub_to_rich_narrative_expansion.md)
  — the SymPy-module-as-canonical-source discipline; per-geometry
  6-subsection shape; the method-implementer-stubs / archivist-expands
  separation. (Generalised in `lessons.md` L5 + the `algebra-of-record`
  skill.)
- [solver-replacement campaign close-out](feedback_solver_replacement_campaign_closeout.md)
  — P8 recipe when a legacy "island" solver is REPLACED by the operator-algebra family:
  theory overhaul (+ rejected-alternative `.. warning::`, refs untouched) · mis-named-law
  re-attribution (rename the LAW word, KEEP math+labels) · investigation-history LIVE/MOOT
  split · sibling forward-ref "expected→now real, consumer still unbuilt" · brief-vs-live
  registry catch · per-phase Dev-history table. Instance #290. Traps: verifies-targets grep
  TREE-WIDE; `:noindex:`-whole-package plain-text (L-002).
- [type-confinement docstring sync](feedback_type_confinement_docstring_sync.md)
  — code-FINAL docstring/comment/error-string sync to a surgical carve
  that CONFINES a subtype to one role (the driver iterate OUTPUT) while
  the INPUT boundary admits the supertype via inheritance + the
  history_depth=0 degenerate, dropping a retired sibling arm. The
  keep-vs-flip rubric (FLIP input-boundary/dual-dispatch/error-string
  prose to the supertype; PRESERVE output/comonad/inheritance prose);
  live-code-corrects-the-stale-docstring (`rhs.boundary` not
  `initial_guess.boundary`, found by grepping the body for `face_view`);
  the two dead-ref kill via grep gate (L-002); flag-don't-fix
  out-of-scope sibling dead-paths; import-sanity check on a code-final
  prose edit. Instance: P4.5 W-C TimedFullField→FullField
  (`refactor/operator-inverse-algebra`, `operator.py`/`solver.py`/`loss_representation.py`).
- [Petrov-Galerkin homogenization reframe](feedback_petrov_galerkin_homogenization_reframe.md)
  — **THE LIVE recipe** (#268 P3) for re-framing a "flux-weighted
  homogenization" theory § from the FORWARD-ONLY "Galerkin in L²(φV)"
  reading to the honest PETROV-GALERKIN framing: test=φ·1_R ≠ trial=1_R;
  the L²(φV)-fold is the Galerkin DEGENERATE of the eigenvalue-consistent
  adjoint-weighted (φ*≠φ) bilinear ⟨φ*,Σφ⟩ case (k stationary by 1st-order
  perturbation theory); measure carries axis+fixed-L²-metric, NEVER the
  discipline. 9-piece order (headline G⁻¹M-of-PG-frame → retraction
  `.. warning::` → trial/test/cross-Gram → why-PG-not-Galerkin → measure-
  never-discipline → mesh-yields-basis+n-D → G⁻¹M verb+MP-zero → two-frame
  `.. list-table::` (Σ=reaction rate / χ=emission rate, via project_through)
  → Mode-11 TEST-side-reader sentinel). Traps: cited test NAMES reframe
  CONCURRENTLY (re-inventory `grep "def test_"`); anchor rename
  `...-galerkin-frame`→`...-petrov-galerkin-frame` + the 1 cross-doc `:ref:`
  in `galerkin_projection.rst`; retire `-galerkin-equals-petrov`, mint
  `-test-functions`/`-metric-fold`/`-bilinear`. Instance: `Solution.homogenize`
  (#268, `refactor/operator-inverse-algebra`).
- [Galerkin-natural-metric reframe](feedback_galerkin_natural_metric_reframe.md)
  — **SUPERSEDED 2026-06-24 by the entry above** (the "Galerkin in L²(φV)"
  framing it documented was the forward-only reading, now reversed to PG).
  Kept ONLY for the WHY-it-was-tried (the metric-fold IS a real identity for
  φ*=φ). Do NOT follow its "Galerkin" framing.
- [domain-op + L2-promotion + asymmetry-law](feedback_domain_op_l2_promotion_asymmetry_law.md)
  — 6-part section shape for a domain OPERATION (transform on a
  `Solution`, not a solver step) born from an L2 module promotion +
  data/behavior type split: lead with THE preservation identity as the
  verifies-target (the rest are vv-status:documented decomposition
  steps), derive each weight forced-not-chosen, `.. warning::` the
  source-group variable-swap trap (vv Mode 2), special-case the simplex
  χ channel, one-line balance argument, asymmetry-law `.. list-table::`
  contracting the deferred sibling (mesh-COUPLED vs DECOUPLED → different
  return types). Instance: `Solution.homogenize`+`MaterialMesh` (#267).
- [orbit-space terminology sweep](feedback_orbit_space_terminology.md) —
  add-aside-then-bridge-then-sweep pattern for introducing a precise
  math term to replace a loose one entrenched as a code-name.
- [auto-generated tables](feedback_autogen_tables.md) — registry-as-SSOT:
  `capability_rows()` metadata fn + generator tool + `builder-inited`
  hook. (Generalised in `lessons.md` L8 + `algebra-of-record`
  "Capability matrices".)
- [audit-then-edit partitions](feedback_audit_partition.md) — the
  KEEP/RELOCATE/TRIM/REMOVE partition-table shape for a read-only
  doc-cleanup audit, with the three-categories-of-narrative rubric.
- **Doc-architecture redesign (#231, OPEN):** theory-page template,
  machine header, prose rebalancing, V&V slices, bibtex. Spec lives in
  the issue; coupled to Nexus semantic search. The standing target for
  any "modernize a theory page" task.
