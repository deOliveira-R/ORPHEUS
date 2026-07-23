# Elegance Enforcer — Memory Index

Slim index. The pattern/anti-pattern CATALOG is the preloaded `coding-elegance`
skill; the cross-cutting institutional SMELLS (twin-delivery plumbing, role-grid,
fuller-view-oracle exception, tells-to-grep) are in AGENT.md §"Institutional
knowledge". This index holds only (1) my review-PROCESS lessons, (2) git-true
active state, (3) durable design-review reference. Campaign play-by-play is
retired — the behavioral lesson is lifted to `lessons.md`.

## 1. Lessons (read FIRST each review)

- [lessons.md](lessons.md) — 12 review-process lessons. The spine: a VIOLATION is
  earned by a bug-habitat argument + proof the two things are coextensive TODAY (so
  it's a NIT) + verification against the LIVE tree, not the diff's docstring.
  Highlights: docstring-names-a-primitive-the-code-doesn't-call is my #1 recurring
  catch (L1); single-source debt collapses on the THIRD divergent consumer not the
  second (L2); two parallel predicates with different return types are NOT a unify
  trigger (L3); stale-doc blast radius lands OUTSIDE the diff (L4); prove the teeth
  bite before crediting a re-baseline (L5); Mode-10 sub-floor terms need structural
  teeth NOT a value band (L6); Pattern-6 TRIM a predicate whose own docstring says
  prod won't use it (L7); emergent invariant-gate over hand-written dunder (L8);
  stash-baseline pyright delta for typing-honesty verdicts (L9); sibling-repo latent
  crash = the guard a sibling method already has (L10); the briefing's DIFF SCOPE is
  a claim to verify against live `git status` — a stale snapshot manufactured a false
  "ungated / promise-not-backed-by-test" finding (L11); PROVE a dead `# type: ignore`
  via a throwaway `reportUnnecessaryTypeIgnoreComment` config — the error-count ratchet
  is blind to it, so a carve can leave dead ignores while claiming it retired them (L12);
  a doc gather/split (monolith → chapter) is a RETIREMENT review — verify twin-content by
  label-uniqueness+prose-grep, math-fidelity by git-HEAD extract+char-diff, account every
  deleted hunk, and flag stale `A=L+C`/`A^{-1}`-for-sweep spelling in NEW text only (the
  honest algebra is `A=L+C−S−B`, sweep=`(L+C)^{-1}`) (L13).

## 2. Active / in-flight state

**None of my own — every reviewed campaign is MERGED to origin/main** (git-verified
2026-07-22). #257/#247/#251/#245/#246/#249/#240/#158/#208/#20/#206 AND the whole
**#226 inverse-as-operator carve** (steps 1–6, incl. the step-6 frozenset/CAP_*
retirement `f4919b1` — the 3-layer predicate/operator-method/realization-verb surface
with PEP-647 `invertible()`/`adjointable()` bridges, confirmed live in operator.py).
Verdicts went to the main agent; durable rulings in the two #226 topic files below;
behavioral lessons in `lessons.md`. (Step-6 delivered 2 POLISH + 1 doc-gate FILE-AS-ISSUE;
on a merged carve those are the archivist's / a GitHub issue's, not my active state.)

One genuinely OPEN branch (reconcile against git before trusting):
- **#236** (`feature/sn-spatial-angular-product`, tip `6409328` — NOT in origin/main).
  My Phase 1b/2/3 reviews COMPLETE and delivered; no pending #236 work. New slice ⟹
  lead with the SN-carve institutional smells in AGENT.md.

> Merge-status in memory goes STALE — a note frozen mid-flight merges in a later
> session. ALWAYS reconcile a "resume/pending X" against
> `git merge-base --is-ancestor <hash> origin/main` before acting; never trust a
> frozen "NOT pushed". (This whole index was rebuilt on that rule.)

## 3. Durable reference (reusable design-review pointers)

- [doc_prose_rebalance_certification.md](doc_prose_rebalance_certification.md) — method for
  certifying a "docstring/comment-only" doc-prose-rebalance batch (#231 Phase 2 P2-*; certified
  P2-D/E/F/G 2026-07-22, all CERTIFIED). The 4 checks: dual token+AST invariance (strip ALL bare
  string-exprs, not just leading docstrings — attribute docstrings are bare Expr too); the
  dropped-contract net (grep `-` diff lines for shape/raise/mutation tokens, confirm each
  survives elsewhere); in-pass fixes vs LIVE code (Hilbert→Euclidean, stale-trait, repointed
  `:ref:`, `:label:`-cut orphan-safety); pointer honesty (resolve ALL + content-check ≥3/batch).

- [pyright_carrier_generic_carve.md](pyright_carrier_generic_carve.md) — #226 pyright
  carrier-generic carve rulings, C2→C4 (reviewed through 2026-07-03). C2: generic-public-surface
  + `Any`-typed-PRIVATE-realization split; `default=Any` REQUIRED under INVARIANCE. C3: F2
  read-mesh-off-the-narrow-boundary-leaf; stale-`# type: ignore`-after-Optional→required smell.
  **C4 (composition-leg generics): the COVARIANT-leg keystone is PRINCIPLED — the forced chain
  retire-casts→pin-legs→pinned-subclass-upcasts-to-defaulted-base→covariance→read-only-`Final`;
  `default=covariant-leg-type` CORRECT (contrast C2's invariant `Any`); it GENERALIZES (#289
  cites it). Operator-leg `__init__` guard vs static pinning = NOT a twin (parse-once/propagate,
  runtime-guard survives erasure). Typed-field ROLE parses = honest #289-tracked family at a
  type-erasure seam. `_CellSolve` Optionals+XOR-guard→ABC+2 kw_only subclasses = Pattern 4,
  tagged-union worse. `other:T→other:"FullField"` = anti-#20 honest direction.** And THE recurring
  approval condition (verified C2+C4) — the ratchet baseline (`pyright_ratchet --update`) must land
  WITH the carve or `test_pyright_ratchet.py` reds `main` on IMPROVEMENT (grep `git status` on it).

- [issue_226_inverse_as_operator_rulings.md](issue_226_inverse_as_operator_rulings.md)
  — the inverse-family design rulings: delegation-not-reciprocal for value-bearing
  leaves (involution-by-identity, no ULP twin, no `1/Σ` units lie); the keystone
  FORCES `solve=inner.apply` on inverse objects; the SweepOperator↔InverseOperator
  wrap-delegate back-half TWIN (collapses at the 3rd sibling Green/Matrix); base
  `LinearOperator` has no `solve`/`inverse` so `_SolveBackedLeaf` is needed + rightly
  private; the shim forwards `CAP_SOLVE` with no `solve` method (keystone-armed,
  latent). Reusable when reviewing the Green/Matrix inverse siblings.

- [coupled_block_boundary_unweld_rulings.md](coupled_block_boundary_unweld_rulings.md) —
  reusable rulings for reviewing a per-system boundary block `B = Σ_x ι_x∘B_x∘π_x`
  (the SN 1b un-weld, PRECEDENT for DSA/multiphysics). What HELD: present-zero-not-None
  is FORCED (presence law) + IS the ι∘B∘π embedding; the `isinstance` narrowing before
  `.split()` is honest (grading lives on B_a not the composite); role grid + twin-delivery
  single-sourced at the `_reflect_*` cores. Recurring FINDINGS to check next time:
  role-parse-guard symmetry (L-010), `_RULED_CORNER_KINDS` single-source, dead moved
  `type: ignore` (L-012), in-diff stale delivery-helper docstring, B_a naming asymmetry.

- [frame_projection_coarsening_shape.md](frame_projection_coarsening_shape.md) —
  homogenize/condense coarsening machinery shape review (P3 merged / P5 draft).
  Six rulings reusable for ANY new coarsening axis: axis-yields-its-views law
  (the keystone), collapse-verb-on-the-XS-container-not-the-frame, the
  HomogenizationFrame/CondensationFrame FALSE symmetry, `project()` is the
  diagonal/PoU special case (type the Gram-structure seam, don't build the dense
  solve), `Mixture.eg` bare-numpy, defer `OverlapBasis`-for-space.

- [doc_label_naming_certification.md](doc_label_naming_certification.md) — method +
  rulings for certifying descriptive equation `:label:` NAMING on the #231 theory
  corpus (greppability / one-concept-one-spelling). Durable ruling: "one domain term,
  two legitimate concepts" (`mg-multiplication-operator` K=A⁻¹F vs
  `multiplication-operator-*` M[f]) = SHOULD-CONSIDER not MUST-FIX — the fix trades
  token-overload for label-vs-prose divergence on entrenched vocab. Includes the
  collision-scan-by-concept-token method, the section-anchor `-eq` rule, and
  non-collisions that look like collisions (P_ii vs r_ii; areas vs fractions).

Sibling-repo (`sphinxcontrib-nexus`, the project's Nexus engine) reviews — kept
because the invariants are reusable and the repo shares the elegance discipline:

- [nexus_runtime_overlay.md](nexus_runtime_overlay.md) — #26 execution-flow overlay
  review. The "MCP tools never raise out" invariant; the missing-stale-node-guard
  latent crash; the orthogonal-capabilities-as-flat-fields (not tagged-union) call.
- [nexus_workspace_resolution.md](nexus_workspace_resolution.md) — workspace/worktree
  wiring; `_switch_workspace` single-source; the `q.knowledge_graph` SSOT accessor
  that retired a private-attr poke; "one server = one workspace" process-local state.
- [nexus_elegance_diagnostics.md](nexus_elegance_diagnostics.md) — a CLI/MCP
  diagnostic family (native_place / twin_paths / discriminations / dead_functions /
  protocol_conformers) that mechanizes my structural-smell first pass, mapped to my
  review axes. NB (2026-06-21): NOT confirmed merged/installed — `nexus` CLI absent
  and these tools are not in the current MCP surface. Verify `nexus <cmd> --help`
  before relying on them; treat as candidate tooling, not available today.
