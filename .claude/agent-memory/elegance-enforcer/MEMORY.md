# Elegance Enforcer — Memory Index

Slim index. The pattern/anti-pattern CATALOG is the preloaded `coding-elegance`
skill; the cross-cutting institutional SMELLS (twin-delivery plumbing, role-grid,
fuller-view-oracle exception, tells-to-grep) are in AGENT.md §"Institutional
knowledge". This index holds only (1) my review-PROCESS lessons, (2) git-true
active state, (3) durable design-review reference. Campaign play-by-play is
retired — the behavioral lesson is lifted to `lessons.md`.

## 1. Lessons (read FIRST each review)

- [lessons.md](lessons.md) — 10 review-process lessons. The spine: a VIOLATION is
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
  crash = the guard a sibling method already has (L10).

## 2. Active / in-flight state

**None of my own.** Every campaign I reviewed — #257 (field-typed operator algebra),
#247/#251/#245/#246/#249 (foundation cleanup), #240/#158 (LD-on-the-DAG), #208/#20
(role-typing), #206 (matvec carve) — is **MERGED to origin/main** (git-verified
2026-06-21: their HEAD commits are ancestors of origin/main). Those reviews are done;
their verdicts went to the main agent at the time; their lessons are in `lessons.md`.

The only genuinely OPEN branch is **#236** (`feature/sn-spatial-angular-product`,
tip `6409328`, NOT an ancestor of origin/main). My Phase 1b/2/3 reviews of it are
COMPLETE and were delivered — I carry no pending #236 review work. If a new #236
slice arrives, lead with the SN-carve institutional smells in AGENT.md.

> Merge-status in memory goes STALE — a note frozen mid-flight merges in a later
> session. ALWAYS reconcile a "resume/pending X" against
> `git merge-base --is-ancestor <hash> origin/main` before acting; never trust a
> frozen "NOT pushed". (This whole index was rebuilt on that rule.)

## 3. Durable reference (reusable design-review pointers)

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
