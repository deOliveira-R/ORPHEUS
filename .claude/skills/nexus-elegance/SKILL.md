---
name: nexus-elegance
description: "Use when reviewing recently-changed code for architectural decay (twin paths, duplicated concepts, broken math-to-code provenance, incomplete retirement) and you want the knowledge graph to scope and corroborate the review. The graph map for the coding-elegance discipline. Examples: \"Review this diff for elegance\", \"Is this a twin path or one source?\", \"Did the refactor actually retire the old pattern?\", \"What's the blast radius of this structural smell?\", \"Does this code still link to the equation it implements?\""
---

# Elegance Review with Nexus

This skill is the **graph map for the `coding-elegance` discipline**. It does
NOT restate the patterns — load `coding-elegance` (`.claude/skills/coding-elegance/SKILL.md`)
for those. It does NOT re-document the tools — see
[../nexus-exploring/reference.md](../nexus-exploring/reference.md) for full
schemas. It owns exactly one thing: **the map from each elegance review axis
to the graph query that corroborates it, plus the discriminators that stop a
graph signal from becoming a false finding.**

IMPORTANT: The graph is a WITNESS, never the accuser. Every finding originates
from reading the changed code against the math/domain. Nexus *scopes* the read
(what changed, what it touches) and *corroborates* the hypothesis (is this
really a twin? does this concept really live in two places?). A finding whose
entire basis is a graph metric — with no code-level bug-habitat argument — is
not a finding.

## Worktree precondition (READ FIRST)

The MCP server may answer from the MAIN checkout's graph. Before any query:
`session_briefing` → confirm the graph's branch matches your checkout. If you
are in a worktree, build Sphinx there, then `use_workspace(<worktree root>)`.
A stale graph produces phantom "severed provenance" and "orphan" findings.

## Workflow

```
0. detect_changes(scope="branch")        → changed lines → graph node worklist
1. impact(<node>, "upstream")  per node  → blast radius = severity multiplier
2. per-axis corroboration (table below)  → confirm/refute the code-read hypothesis
3. provenance_chain + verification_coverage → math↔code↔test chain intact?
4. callers(<predecessor>) + retest        → retirement complete + tests rewired?
5. assemble: severity = (code-read finding) × (blast radius)
```

## Tool → elegance axis map

Each row: the axis from the Elegance Enforcer methodology, the tool, the
**signal** (what the graph shows), the **smell it indicates**, and the
`coding-elegance` pattern it serves.

### Axis 1 — Data structures (Pattern 4: illegal states unrepresentable; anti-pattern 13: bare numpy across boundaries)

- **MCP:** `context(node_id="py:function:...")`,
  `neighbors(node_id, direction="out", edge_types="type_uses")`
- **CLI:** `nexus context <node_id> --db <path>` ·
  `nexus neighbors <node_id> --db <path> --direction out --edge-types type_uses`
- **When it fires, this is the smell:** a function whose `type_uses` edges are
  all `np.ndarray` / `dict` / `tuple` and that crosses a module boundary —
  distinct physical quantities sharing one representation, no static check
  against an argument swap. (anti-pattern 13, the densest bug cluster.)
- **Usage example:** reviewing a new `solve_*` that takes
  `(materials, mesh, quad)` positionally — `context` shows three `type_uses`
  to bare containers and `callers` shows 4 call sites → positional-swap bug
  habitat. Demand a typed `PhaseSpace` / keyword-only signature.
- **Serves:** Pattern 4 (make illegal states unrepresentable).

### Axis 2 — Path multiplicity (Pattern 2: single source of truth; anti-pattern 1: twin paths)

- **MCP:** `twin_paths(min_similarity=0.7)` for the whole-graph sweep;
  `callees(node_id=<suspect A>)` + `callees(node_id=<suspect B>)` and
  `shortest_path(source=<A>, target=<B>)` to adjudicate a specific pair.
- **CLI:** `nexus twin-paths --db <path> --exclude derivations,scratch` ·
  `nexus callees <A> --db <path>` · `nexus shortest-path <A> <B> --db <path>`
- **When it fires, this is the smell:** two functions whose AST bodies share a
  high shingle fraction (`twin_paths` similarity ≥ 0.8) and that do not call
  each other → two implementations of one computation. The fingerprint catches
  the array math (`@`, einsum, slicing) the call graph cannot see. (anti-pattern
  1, the load-bearing failure mode — Phase F ERR-026.)
- **Usage example:** `nexus twin-paths` flags
  `CPMesh._resolve_bc ≈ MOCMesh._resolve_bc` (0.90, cross-module) → confirm by
  reading: same registry-lookup + error string, only the default differs →
  VIOLATION, extract a shared `resolve_bc_factory`. Contrast
  `StreamingOperator.apply / apply_transpose` (0.92): symmetric-by-design, NOT
  a twin (see false-positive table).
- **Serves:** Pattern 2 (composition over duplication).

### Axis 3 — Procedural transcription (anti-pattern 2: for-loops over a domain with algebra)

- **MCP:** `callees(node_id=<fn>, transitive=True, max_depth=2)`
- **CLI:** `nexus callees <fn> --db <path> --transitive --max-depth 2`
- **When it fires, this is the smell:** the transitive callee set is entirely
  numpy/stdlib primitives with no domain operator (`.apply`, `.solve`, `+`/`@`
  on operator types) → the code is *about* the algebra, not *the* algebra.
- **Usage example:** a function whose callees are only `np.sum`, `np.zeros`,
  `range` and that lives beside an operator algebra → it transcribes the math
  procedurally instead of invoking `L.solve(q)`. Flag for operator-expression
  rewrite.
- **Serves:** Pattern 1 (match the algebra via dunder methods).

### Axis 4 — Single source of truth / missing types (Pattern 7: normalise at the definition site; "a repeated conditional is a missing type")

- **MCP:** `discriminations(min_sites=2)` for tag-dispatch fan-in;
  `native_place(min_callers=1)` for free functions coupled to one class;
  `graph_query("* -implements-> equation")` for one concept in two code nodes.
- **CLI:** `nexus discriminations --db <path> --exclude derivations` ·
  `nexus native-place --db <path>` ·
  `nexus graph-query "* -implements-> equation" --db <path>`
- **When it fires, this is the smell:** the SAME tag branched on at ≥2 sites
  (`discriminations` — `if geometry == "spherical"` re-asked everywhere = a
  missing type / absent single dispatch); a module-level function whose every
  caller is a method of one class (`native_place` — Feature-Envy, belongs
  inside it); one equation with ≥2 `implements` edges (a concept in two
  places).
- **Usage example:** `nexus discriminations` shows `coord` discriminated at 13
  sites → the geometry partition is a missing polymorphic dispatch (ties to the
  dimension-agnostic carve). `nexus native-place` flags
  `_subdivide_zone → Mesh1D` (cross-module, private, 0 excluded) → move it in;
  but `gauss_legendre_on_mu → Quadrature` carries `likely_free_primitive`
  (public, tested) → leave it free.
- **Serves:** Pattern 7, Pattern 5, and "a repeated conditional is a missing
  type — discriminate once, at the boundary" (anti-patterns 3/4/7).

### Axis 5 — Math/domain alignment (the master standard: "reads like the math")

- **MCP:** `provenance_chain(node_id=<solver fn or equation>)`,
  `verification_coverage(status_filter="")`
- **CLI:** `nexus trace <test_node_id> --db <path>` (for the failing-test
  direction) · `nexus audit --db <path>`
- **When it fires, this is the smell:** `provenance_chain` returns empty or
  broken for a function that visibly implements a cited equation → the
  code↔equation↔citation chain was severed by the change. A refactor that
  breaks provenance is a Cardinal-Rule-1 correctness regression even when
  tests pass.
- **Usage example:** after a carve, `provenance_chain` on the relocated kernel
  no longer reaches its `:label:` equation → either the docstring label was
  dropped (fix it) OR the graph is stale (rebuild + `use_workspace`, then
  re-check — see Worktree precondition).
- **Serves:** the Master Standard + Cardinal Rule 1.

### Axis 6 — Architectural forwardness / retirement (coding-standards: aggressive retirement)

- **MCP:** `dead_functions()` for orphaned predecessors;
  `callers(node_id=<predecessor>)`,
  `impact(target=<predecessor>, direction="upstream")`,
  `retest(scope="branch")`
- **CLI:** `nexus dead-functions --db <path> --exclude derivations` ·
  `nexus callers <predecessor> --db <path>` ·
  `nexus retest --db <path> --scope branch`
- **When it fires, this is the smell:** the "retired" symbol still has live
  (non-test) callers → retirement incomplete, the predecessor pattern lingers
  and invites accidental extension of the wrong path. `dead_functions` lists
  the orphans directly (private + undecorated ranked first); `retest` not
  reaching the successor's tests → test migration skipped.
- **Usage example:** a refactor claims to replace `_compute_LpC` with a walk
  method. `callers(_compute_LpC)` returns 2 live operator nodes → not retired.
  Demand deletion + caller migration before approval. (Beware the
  retained-oracle exception in the false-positive table — a fuller view kept
  for an equivalence test is NOT dead weight.)
- **Serves:** coding-standards "retire as you go"; Architectural forwardness.

### Axis 7 — Conformance / unused weight (Pattern 4: make conformance explicit; anti-pattern 11: no "for future use")

- **MCP:** `protocol_conformers(min_methods=2)` for undeclared structural
  conformers; `dead_functions()` / `callers(node_id=<symbol>)` for dead weight.
- **CLI:** `nexus protocol-conformers --db <path> --exclude derivations` ·
  `nexus dead-functions --db <path>`
- **When it fires, this is the smell:** classes that satisfy a `Protocol`'s
  method-set without declaring it (`protocol_conformers`) — the `inherits` edge
  only sees explicit subclassing, so a structural conformer is invisible;
  declaring it (or confirming the Protocol is load-bearing) makes the contract
  checkable. A symbol/argument with zero consumers = dead weight kept "for
  future use".
- **Usage example:** `nexus protocol-conformers` shows `LinearOperator` ←
  StreamingOperator, OperatorSum, … (they inherit a mixin, not the Protocol) →
  candidates to declare conformance. NOTE: this is a method-NAME heuristic; the
  authoritative check is the type checker (pyright / LSP `goToImplementation`).
- **Serves:** Pattern 4 (illegal states unrepresentable); anti-pattern 11.

## Severity from blast radius

`impact(<node>, "upstream")` does not *create* findings — it *scales* them.

```
severity = (code-read finding strength) × (upstream blast radius)
```

A confirmed twin path with 14 upstream dependents is a VIOLATION (14× the
divergence surface). The same twin with 1 caller is a CONCERN. NEVER open a
finding whose only basis is a blast-radius number.

## Corroboration across tools

Two tools flagging the SAME symbol is a stronger escalation than either alone.
`native_place` and `discriminations` independently flag `axes_from_legacy_mesh`
and `_subdivide_zone` (a free function coupled to `SNMesh` that also branches
on geometry) — the agreement says "this symbol is the seam of a missing
abstraction," not two unrelated nits. Cross-reference the tool outputs before
ranking.

## False-positive table (graph says "smell"; READ before flagging)

| Graph signal | Looks like | Discriminator (what the graph can't see) | Correct verdict |
|---|---|---|---|
| `twin_paths` high similarity, or 2 nodes `implements` one equation | Twin path (Pattern 2) | One is apply (`Aψ`, source/sink), the other is solve (flux). matvec≡sweep is a code FACT; two genuine walks sharing one denominator leaf. | Shared leaf → PASS. Inlined-each → VIOLATION |
| 2 `calls` edges reach the SAME operator node | Twin path | Twin DELIVERY (driver-fold vs Krylov-inline), single-sourced at the operator. | Both routes consume one operator + byte-identical overlap → CONCERN, demand cross-ref + removal trigger |
| `twin_paths` `apply`/`apply_transpose`, `domain`/`codomain` | Twin path | Symmetric-by-design forward/adjoint or dual accessors. | PASS unless they inline divergent arithmetic |
| `native_place` with `likely_free_primitive=true` | Feature-Envy | Public + independently tested = a primitive that is *correctly* free. | PASS (leave free) |
| High `impact` upstream on a `god_nodes`/`bridges` hit | Risky change | The node is the single source doing its job. | Amplifies a confirmed finding; NOT a finding alone |
| `provenance_chain` empty after a refactor | Severed math chain | **Usually nothing is wrong** — the symbol implements no equation, which is the COMMON case (`[M]` 639 on one corpus). Or a stale graph (Sphinx not rebuilt in the worktree). | Check `also_on_these_pages` first, then compare with the PRE-change reply; rebuild + `use_workspace` if in doubt. **Flag only on a DELTA**, never on emptiness alone |
| `dead_references` hit on a dynamic attribute | Dead doc reference | `__getattr__` / metaclass magic creates it at runtime | Static analysis cannot see it — verify the attribute really is gone before flagging |
| `dead_functions` / `callers(<fuller view>)` = test nodes only | Orphan / dead weight | Deliberate retained oracle (aggressive-retirement EXCEPTION). | Wired to `optimized ≡ oracle` test (`tests` edges) → PASS. No test edges → VIOLATION |
| `dead_functions` with `decorated=true` or `public=true` | Dead code | Registry/route/property decorator invokes it indirectly; public = entry point. | Read for dynamic dispatch / external callers before flagging |
| `protocol_conformers` match | Undeclared conformer | Method-NAME match only (signatures ignored). | Confirm with the type checker (pyright / LSP goToImplementation) |
| `discriminations` single dispatcher reads one tag | Missing type | A guard + its dispatcher SHARING one predicate cannot drift. | Flag only a SECOND independent spelling on one side |

## Closing checklist (before returning the review)

1. Every VIOLATION has BOTH a `coding-elegance` pattern citation AND a
   bug-habitat argument. Graph metric alone is never the basis.
2. Every twin-path finding was checked against the apply-vs-solve and
   twin-DELIVERY discriminators before being called a VIOLATION.
3. Provenance/orphan findings were re-checked after confirming the graph
   matches the checkout (`session_briefing`), not against a stale MAIN graph.
4. Retained-oracle nodes were checked for `tests` edges before being called
   dead weight; `protocol_conformers` hits were confirmed with the type checker.
5. Blast radius scaled severity; it did not manufacture a finding.
6. Confirmed smells hand off to `nexus-refactoring` for the safe rename/extract;
   required changes name the destination (the typed structure, the shared
   primitive, the predecessor to delete) — not the path.
