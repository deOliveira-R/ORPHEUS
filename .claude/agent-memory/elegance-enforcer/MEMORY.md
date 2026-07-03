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
  is blind to it, so a carve can leave dead ignores while claiming it retired them (L12).

## 2. Active / in-flight state

**None of my own.** Every campaign I reviewed — #257 (field-typed operator algebra),
#247/#251/#245/#246/#249 (foundation cleanup), #240/#158 (LD-on-the-DAG), #208/#20
(role-typing), #206 (matvec carve) — is **MERGED to origin/main** (git-verified
2026-06-21: their HEAD commits are ancestors of origin/main). Those reviews are done;
their verdicts went to the main agent at the time; their lessons are in `lessons.md`.

Two genuinely OPEN branches (reconcile against git before trusting):
- **#236** (`feature/sn-spatial-angular-product`, tip `6409328`). My Phase 1b/2/3
  reviews are COMPLETE and delivered — no pending #236 work. New slice ⟹ lead with
  the SN-carve institutional smells in AGENT.md.
- **#226** (`refactor/inverse-as-operator`, step-1 `e7115a2`, step-2 `cc293ef`, step-3
  COMMITTED `1ab7429`, step-4 UNCOMMITTED working tree). Reviewed step-1 (2026-07-01) +
  step-2 (2026-07-01, BLOCK-on-doc-retirement) + step-3 (2026-07-02, SHIP-WITH-FIXES) +
  **step-4 "GreenOperator + OperatorSum.inverse() + InverseWrapMixin extraction"
  (2026-07-02) — verdict SHIP-WITH-FIXES: code+tests SHIP-QUALITY, zero code defects; the
  mixin fired EXACTLY as my step-1 ruling predicted (back-half extracted at the 3rd
  sibling; leaves keep guard/apply/repr). All 4 attack points PASS: refinement-loop is
  bit-identical-to-hand-SI (0-ULP) so it wraps not re-implements + the §18.A near-critical
  mutation gate (budget BETWEEN increment-stop ~460 & refinement-close ~920) proves it's
  driven-not-checked; `is_invertible`=leading-term honest (constructibility, not
  convergence-guarantee → loud ConvergenceFailure); NO prod consumer of `OperatorSum.solve`
  so fresh-Green-per-call is off-hot-path; `_SolveBackedLeaf→_InvertibleForward` rename +
  0 residue. 77 green/0.31s, pyright 0-new (ratchet 148 holds).** The prior doc-retirement
  gate is now CLOSED (no broken xrefs, historical framing). Remaining = 3 doc-sync fixes
  (NOT code): sweep_operator.py:12 stale-"Green deferred/Krylov" (inline before commit);
  operator_algebra.rst:8285 frames #285 as open when carve resolved it structural
  (archivist); Green theory-page section additive (Rule 3). Step-4 rulings + step-1/2/3 in
  [issue_226_inverse_as_operator_rulings.md]. Step 5 (`MatrixInverseOperator`) reviewed
  2026-07-02 (SHIP-WITH-FIXES, 1 docstring drift). **Step 6 (carve P4 — the frozenset
  retirement) reviewed 2026-07-02 @ `f4919b1` (W1+W2 landed; MY review IS the W3 gate):
  verdict PASS — ZERO MUST-FIX code defects, ZERO VIOLATIONS.** The `capabilities:
  frozenset[str]` + CAP_*/MissingCapability RETIRED tree-wide; each axis now 3-layer
  (predicate / operator-method / realization-verb) with PEP-647 TypeGuard checked bridges
  (`invertible()`/`adjointable()`). Verified against LIVE tree: retirement complete (exactly
  6 docstring-history mentions, 0 functional caps reads); 0 `# type: ignore` ADDED (9
  removed); ratchet 148→145 honored; all 12 bridge sites §44.E-clean (LinearOperator-typed,
  narrow to a downstream `.inverse()`/`.apply_transpose()`); message family axis-uniform;
  keystone v2 single-sourced + `-O`-safe raises + bridge-consistency-leg gives the TypeGuards
  teeth (95 gates green live); RankOne-bug-class did NOT recur (every caps-removal got a real
  `is_adjointable` override / role-base inherit). Findings = 2 POLISH + 1 FILE-AS-ISSUE only.
  POLISH-1 (brief-requested judgment): `OperatorSum.is_invertible` (operator.py:1151)
  `getattr(self.a,"is_invertible",False)` is a SMELL not honest duck-typing — its OWN sibling
  `is_adjointable`:1124 + EVERY other composite predicate read operands DIRECTLY, so the
  `False` fallback is provably DEAD (coextensive today → NIT); fix `return self.a.is_invertible`.
  POLISH-2 (softer): the 2 ctor-guard getattr (KEigenvalue iteration.py:1090, GreenOperator
  green_operator.py:256) are a 3rd spelling of "is X invertible?" alongside the canonical
  `invertible()` bridge (used at `_seeded_inverse`:263 same file) — defensible as boundary-parse
  but prefer direct read. FILE-AS-ISSUE (the recurring doc gate, route to W3 archivist Sphinx
  pass): 13 BROKEN `:class:MissingCapability` xrefs (operator_algebra.rst 9 / galerkin 2 /
  discrete_ordinates 2) + stale CAP_/OperatorSum.solve prose across 6 theory pages;
  operator_algebra.rst CAP_SOLVE-gate section needs a REWRITE not find-replace.

> Merge-status in memory goes STALE — a note frozen mid-flight merges in a later
> session. ALWAYS reconcile a "resume/pending X" against
> `git merge-base --is-ancestor <hash> origin/main` before acting; never trust a
> frozen "NOT pushed". (This whole index was rebuilt on that rule.)

## 3. Durable reference (reusable design-review pointers)

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

- [frame_projection_coarsening_shape.md](frame_projection_coarsening_shape.md) —
  homogenize/condense coarsening machinery shape review (P3 merged / P5 draft).
  Six rulings reusable for ANY new coarsening axis: axis-yields-its-views law
  (the keystone), collapse-verb-on-the-XS-container-not-the-frame, the
  HomogenizationFrame/CondensationFrame FALSE symmetry, `project()` is the
  diagonal/PoU special case (type the Gram-structure seam, don't build the dense
  solve), `Mixture.eg` bare-numpy, defer `OverlapBasis`-for-space.

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
