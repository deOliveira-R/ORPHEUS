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

- [lessons.md](lessons.md) — 10 documentation/Sphinx/knowledge-architecture
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
  code-point underlines (L9); V&V vocabulary curation (L10).

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
- [capstone-architecture-page](feedback_capstone_architecture_page.md) —
  shape for a NEW capstone page documenting the LAYER above existing
  per-method pages (cross-ref-not-duplicate; measurement+decision
  section; verbatim user-decision admonition). Instance:
  `loss_representations.rst`.
- [operator-classes→frame-faces re-homing](feedback_operator_classes_to_frame_faces_rehoming.md)
  — doc-sweep recipe when a refactor retires standalone projection M +
  reconstruction R operator classes into the two FACES of one discrete
  `Frame` (`frame.analysis`/`frame.reconstruction`): re-home onto the
  abstraction + add the frame-theory framing (T/T*/S/tight-frame), don't
  find-replace; KEEP documented-only eq-labels by name; the
  concurrent-code-carve-vindicates-re-verify-live lesson; the
  co-located-separate-retirement scoping discipline. Instance:
  Frame/Basis carve (#263, `refactor/operator-inverse-algebra`).
- [algebra-of-record stub→narrative](feedback_stub_to_rich_narrative_expansion.md)
  — the SymPy-module-as-canonical-source discipline; per-geometry
  6-subsection shape; the method-implementer-stubs / archivist-expands
  separation. (Generalised in `lessons.md` L5 + the `algebra-of-record`
  skill.)
- [Galerkin-natural-metric reframe](feedback_galerkin_natural_metric_reframe.md)
  — doc-recipe when a carve proves a "flux-weighted average" is the
  L²(weighted-metric) ORTHOGONAL (Galerkin) projection, NOT
  Petrov-Galerkin (the PG reading is the unweighted-dV-metric artifact;
  flux multiplier lives in test-fn OR measure = same map; Frame reads it
  off the measure → Galerkin). 8-piece derivation order (P0 trial space →
  normal eqns → measure-is-DERIVED-not-chosen → two-readings-same-map →
  Radon–Nikodym μ_φV=φ·μ_V & why DiscreteMeasure stays 1-D → mesh-YIELDS-
  IndicatorBasis → R∘G⁻¹∘M with 1/Φ_R in the space metric → Mode-11
  sentinel); + the sibling-page staleness TOMBSTONE (L-007); + asymmetry-
  law gains a frame reading (geometric-K=mesh-coupled vs spectral-K=
  decoupled). φV-vs-dV discriminator is the load-bearing guard. Instance:
  `Solution.homogenize` (#268, `refactor/operator-inverse-algebra`).
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
