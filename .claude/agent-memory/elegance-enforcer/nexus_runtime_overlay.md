---
name: nexus-runtime-overlay
description: sphinxcontrib-nexus (sibling repo) #26 dynamic execution-flow overlay review — runtime.py pure/wrapper split, 3 GraphQuery overlay methods, sidecar-survives-rebuild stale-node hazard
metadata:
  type: project
---

Reviewed `feat/runtime-overlay` in `/Users/rodrigo/git/sphinxcontrib-nexus`
(#26 dynamic execution-flow overlay; cProfile + coverage.py trace joined onto
static node-IDs via a `_nexus/traces/<run>.json` sidecar, joined at QUERY time —
never written into `graph.db`). 516+27 tests green, behavior pinned. Substrate =
ORPHEUS SN solver traced under cProfile. Design doc =
`docs/notes/issue26_trace_overlay_spike.md`.

**Verdict: PASS-with-fixes (one should-fix latent crash + two nits).** The
core design is sound and TRUE to intent.

What good looked like (reinforce, do NOT re-flag):
- **The sidecar invariant holds.** `runtime.py` has no Sphinx import; the store
  writes `_workspace.db_path.parent / "traces"` (sibling of graph.db, not inside
  it); nothing in the module touches graph.db. Re-bind by node-ID at query time.
- **Pure/wrapper split is clean (Pattern 5).** `overlay_cprofile(stats_dict,
  index, ...)` / `overlay_coverage(json_dict, index, ...)` are pure over parsed
  structures; `ingest_*` are thin path-loaders (`pstats.Stats`/`json.loads` then
  delegate). Testable without the artifact producers installed. `build_node_index`
  + `resolve_node` are separately reusable primitives.
- **`resolve_node` decorator-window logic is CORRECT.** `bisect_right(starts,
  firstlineno + DECORATOR_WINDOW)` bounds candidates; exact `ln==firstlineno`
  short-circuits; enclosing-body fallback keeps innermost (latest start) via
  ascending-sort overwrite of `best`. The window asymmetry (cProfile
  `co_firstlineno`=first-decorator vs AST `node.lineno`=def-line) is documented
  with the measured 68%→97% join lift. `None` for lambdas/closures is by design.
- **`runtime_edges(mode=...)` is the RIGHT factoring, not a latent twin.**
  dynamic_only/fired collapse to one body via `want_static=(mode=="fired")` +
  `((u,v) in static_calls)==want_static`; `dead` is a genuinely different set op
  (iterate static_calls ∩ reachable) — co-located in one method because all three
  are "edge-overlay views". Correct parameterized primitive over 3 twin methods.
- **3-layer convention mirrored faithfully.** GraphQuery method → `@nexus_tool`
  returning `to_json(to_dict(results))` → CLI subparser+`_run_*` handler, exactly
  like `dead_functions`/`protocol_conformers`. `runtime_runs`/`runtime_ingest`
  correctly return raw `to_json(...)` (summary dicts, not dataclass lists — no
  to_dict needed). `--all`→`partial_only=not args.all` CLI inversion is clean.

THE SHOULD-FIX (latent crash; the sidecar's whole reason to exist triggers it):
- **`runtime_edges` omits the `if node_id in self._g` guard that BOTH sibling
  overlay methods have** (`runtime_hotspots` line 2794, `runtime_branches` line
  2884 both guard; `runtime_edges` does not). Mechanism: a trace-resolved node can
  be RENAMED/REMOVED by a `sphinx-build` rebuild between ingest and query — the
  exact stale-binding the sidecar design embraces. `_node_result` calls
  `self._g.degree(node_id)`; for a MISSING MultiDiGraph node that returns a
  `DiMultiDegreeView` (NOT an int, does NOT raise) → carried through `asdict` →
  `to_json` (which has NO `default=`) → **`TypeError: ... not JSON serializable`**
  at the server/CLI boundary. Violates this repo's "MCP tools never raise out"
  invariant. dynamic_only/fired endpoints come from `run.edges` (ingest-time
  index) so are the vulnerable ones; `dead` mode is safe (endpoints from current
  static_calls). Fix: guard both endpoints `if u in self._g and v in self._g`
  before emitting, in the dyn-iterating branch (the `dead` branch already only
  emits current static endpoints).

Nits (cosmetic, non-blocking):
1. `overlay_coverage` redefines a `def in_range(line)` closure INSIDE the
   `for lineno,end,node_id in spans` loop (query.py runtime.py ~254). Called
   same-iteration so late-binding is NOT a bug, but it's per-iteration redefinition
   + a procedural smell. A module-level `_in_range(line, lo, hi)` or inlined
   `lo<=ln<=hi` reads cleaner. (Functional, low priority.)
2. Stale-node-after-rebuild has NO test (the design's headline scenario). Add a
   test that ingests, then drops a node from the graph, then runs all 3 overlay
   queries — pins the guard above and documents the re-bind contract.

Cross-ref: [[nexus-workspace-resolution]] — same repo; `q.knowledge_graph` SSOT
accessor (used by ingest here) was introduced there retiring a private-attr poke;
"MCP tools never raise out / return error payloads" is that review's documented
invariant, which the should-fix above would violate.

NB: the Phase-1+2 should-fix above (missing `runtime_edges` stale-node guard) WAS
FIXED — `query.py:2877` now guards `if u not in self._g or v not in self._g: continue`
in the dyn branch, with the documented comment. Confirmed during the Phase-3 review.

---

**PHASE 3 review (uncommitted working tree vs HEAD `f7fb5d8`, 2026-06-17). Verdict:
CONCERNS-RAISED, no blockers.** 528 tests green (behavior pinned; review = STRUCTURE).
Three additions: (1) `merge_runs` multi-run union, (2) `_is_accessor` edge classifier +
`substantive_only`, (3) viztracer backend + `runtime_timeline`.

HIGH-STAKES STRUCTURAL CALL (got it RIGHT): `RuntimeRun` now carries FOUR metric families
(calls/edges/coverage/timeline) as flat per-family empty-dict-default fields — NOT a tagged
union/per-kind subtype. CORRECT: the families are ORTHOGONAL CAPABILITIES a run may have, not
mutually-exclusive variants. `merge_runs` unions WHATEVER families inputs carry (a heterogeneous
cprofile+coverage merge legitimately holds both) → a union would make the merged run
unrepresentable + force per-variant branching. `kind` is NOT dispatch in the query hot path
(each query reads its family unconditionally, empty-safe) — it's provenance/label. Same shape
as ORPHEUS per-leaf field-role-vs-forced-base.

What good looked like (reinforce):
- **`runtime_ingest` kind→ingester DICT dispatch** (`server.py:759`): `{KIND_*: ingest_*}` +
  membership-check + single uniform call site replaced if/elif/else string-branch. Textbook
  anti-pattern-#4 avoidance (3 ingesters share one signature). Nit: error string could
  `"|".join(ingesters)` to SSOT the valid-kinds list.
- **`runtime_timeline` carries the Phase-1+2 stale-node guard** (`query.py:2989`) — same re-bind
  guard as the siblings, not forgotten; test pins it.
- **`substantive_only` filter BEFORE the `limit` slice** (`query.py:2882`) — filter-then-truncate
  yields N substantive edges, not "N raw minus accessors"; opposite ordering is the common bug.
- **either-endpoint accessor classification** correct for the #16 goal (goal is docstring-stated,
  so checkable against intent).
- **`merge_runs` merge LAWS all correct**: ncalls/tottime=SUM, cumtime=MAX (cross-run cumtime not
  additive; CONSISTENT with overlay_cprofile's per-run max — Pattern 7), edges count=SUM, coverage
  missing_arcs = set.INTERSECTION ("still-missing iff missing in EVERY run" = exact dead-code dual
  to edges-union). Arc semantics exact + test-pinned.
- **`@property`/`@cached_property` decorator signal is correct-by-construction**: verified
  ast_analyzer.py:874 renders decorators via `ast.unparse` (stored as TUPLE), so all three
  spellings contain substring `"property"`; `_is_accessor`'s `(... or [])` iter works for tuple
  AND the test's list.

APPROVAL CONDITIONS (4; defensible-contract / one-twin; no blocker):
1. **SHOULD-FIX twin path**: `server._load_runs` (`server.py:165`) ≈ `cli._runtime_load_many`
   (`cli.py:1355`) — comma-split + `merge_runs(..., name=",".join(wanted))` duplicated verbatim;
   only the per-name loader differs. PLUMBING-TWIN (single-sourced at `merge_runs`, two parse-and-
   fold delivery routes — institutional note #1). Habitat: change separator/label → edit one,
   other diverges silently → server/CLI accept multi-run differently, tests pass (each pinned
   apart). FIX = extract `runtime.load_and_merge(names, load)`; wrappers carry only the loader
   binding (store vs db_path = legit layer-symmetry). The SINGULAR `_load_run`/`_runtime_load` pair
   differs substantively (no shared convention) → leave.
2. **SHOULD-FIX viztracer depth stack** (`runtime.py:397-404`): `min_depth` from an end-times-only
   `open_ends` stack + `while open_ends[-1] <= ts: pop`. Wrong at corners: zero-dur child at parent
   boundary (`parent_end == child_ts`) pops parent BEFORE measuring depth → child at grandparent
   depth. Root: stack discards the containment it claims (anti-pattern #16 — tree from flatten+sort
   +end-time). Habitat: `min_depth` IS the `runtime_timeline(max_depth=1)` filter input (the stated
   viztracer value-add) → misdepthed helper masquerades as top-level stage, defeats the filter.
   Test probes only clean separated intervals. FIX = push `(start,end)` + never pop same-`ts`
   parent + add corner tests; OR prove strict containment empirically + NARROW the docstring claim.
3. **CONCERN `_is_accessor` ≤2-line fallback** (`query.py:2908`): `(end-ln)<=2` magic threshold
   (unnamed; span ≠ body-line-count). Self-defeating false-negative: a 2-line GENUINE dispatcher
   (`return self._dispatch[type(n)](n)` — the exact #26 signal) → classified accessor → DROPPED by
   `substantive_only` → filter hides what it reveals. Decorator path immune; only line-count net
   does this. FIX = name `_ACCESSOR_MAX_SPAN=2` + derivation comment; document the false-negative
   (current "without risking real small methods" overclaims); consider exact-set decorator membership.
4. **CONCERN self-disproving docstring** (`runtime.py:120-122`): "Exactly one metric family
   populated, per kind" REFUTED by the same commit's `merge_runs` (heterogeneous union) +
   `timeline`/`KIND_MERGED`. Institutional note #7 (false-invariant comment = bug habitat); most
   insidious because "just a docstring" = the contract the next agent reads. Cardinal Rule 1. FIX =
   rewrite to the bag-of-orthogonal-overlays contract.

Nits: SSOT ingest valid-kinds error from the dict; `merge_runs` lines_hit=MAX is an approximation
not a set-union (arc intersection IS exact; one-line docstring note). RETRACTED an early "eager
accessor compute wasteful" objection — eager keeps `RuntimeEdgeResult` a clean value object (flag
is part of identity, available to any JSON consumer) = the more elegant choice.

Arch opportunity (spike-doc, not issue): depth reconstruction is the ONE place the overlay rebuilds
a TREE (call-stack nesting) from a FLATTENED trace; a future "viztracer call-tree edges" Phase-4
needs the interval-containment tree, not the end-time stack.
