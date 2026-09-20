---
name: nexus-workspace-resolution
description: sphinxcontrib-nexus (sibling repo) four reviews — workspace/worktree wiring, the ambient brief layer, the token-budget caps, and the LSP-parity oracle (set-equality-with-reducer-encoded-tolerance)
metadata:
  type: project
---

Reviewed `feature/workspace-worktrees` (8 commits `d6f9674..a7483ae`) in
`/Users/rodrigo/git/sphinxcontrib-nexus` — a sibling project to ORPHEUS that
re-uses the coding-elegance discipline. Its own CLAUDE.md states conventions:
plain pytest (no `-O`/markers), pyright via `pyrightconfig.json`, producer-side
normalization invariant, "one server process = one workspace" process-local state.

**Verdict: PASS-with-nits (cosmetic).** 384 tests green, pyright clean.

What good looked like here (reinforce, don't re-flag):
- **`_switch_workspace` is a genuine single source** (Pattern 2). Both the
  `use_workspace` tool and `_auto_align_workspace` delegate to it; the switch
  core has ZERO duplication. `use_workspace` thinned to resolve+delegate;
  `_briefing_payload` shared by the async tool and the `nexus://briefing`
  resource (the prior twin between them collapsed). The `ingest` tool's
  hand-rebuilt `KnowledgeGraph(); kg._graph = q._g` (private-attr poke) replaced
  by `q.knowledge_graph` property — that retired a real SSOT violation.
- **`KnowledgeGraph(graph=None)` wrap constructor is the right boundary** (not a
  classmethod): the edge-key continuation (`max(int_keys, default=-1)+1`) is an
  INVARIANT of holding a graph, so it belongs in `__init__` where every entry
  (empty or wrapped) passes through it. Fixed a latent key-collision bug
  (`_edge_key=0` over a wrapped graph silently updates existing parallel edges).
- **`resolve_checkout_root` 3-form order is principled**: absolute path → worktree
  name/branch match → relative dir (last resort). Ambiguity (name matches >1
  checkout) RAISES `WorkspaceResolutionError` listing candidates — illegal-states
  surfaced, not papered. Worktree-match wins over relative-dir by design
  (documented). Parse-at-boundary (Pattern 4): tool catches the typed exception,
  returns an error payload, never raises out of an MCP tool (their invariant).
- **`direction` validation** at `impact`/`neighbors` tool boundary + `Literal`
  type on `assemble_neighbors` — illegal-states-unrepresentable internally,
  validated at the stringly-typed MCP edge. Correct split.
- **README↔FastMCP registry drift guard** (`test_server_registry.py`) — pins a
  Pattern-7 convention (the tool count lived in 3 places, all disagreeing) to the
  registry as SSOT. Same shape as the version single-source (`dynamic=["version"]`
  from `__init__.__version__`).

Nits (cosmetic, non-blocking, did NOT require rework):
1. `_auto_align_workspace` re-derives an `info` block from `_switch_workspace`'s
   return dict by string-key probing (`outcome.get("switched")`/`"error"`/`"hint"`).
   This couples to the switch payload's dict SHAPE — a stringly-typed seam. Low
   bug-habitat (one consumer, one producer, both in one file) so CONCERN-not-
   VIOLATION; a typed `SwitchOutcome` result would remove the re-probe but is
   premature at one consumer (Pattern 6). Watch for a 2nd consumer.
2. `checkout_containing` "deepest match" via `max(key=len(root.parts))` is correct
   for the nested-worktree case (`.claude/worktrees/<name>` under main) — tie on
   identical depth is impossible (distinct paths). Sound.

The `default_branch` (origin/HEAD → main/master) retired the old `for base in
("main","master"): ... if files: break` fallback in `query.py` that conflated
"ref absent" with "no .py changed" — a genuine convention bug, not just style.

---

**AMBIENT-LAYER review (later 4 commits `700df12→72b36b9→cff9113→2d7d13f`, vs
parent `9b1fb68`, 2026-06-12). Verdict: PASS-with-nits.** 423 tests green,
pyright 0/0. The "ambient push channel": `brief.py` (NEW, edit-time `file_brief`
reads SQLite DIRECTLY — no networkx, ~100ms warm budget), `GitProvenance.from_stamp`
+ `STAMP_*` vocab, `changed_files` None/`frozenset()` tri-state, server
position-staleness warning on `node_at`, `NodeResult.file_path/lineno` reverse bridge.

What good looked like (reinforce):
- **`STAMP_*` + `from_stamp` is a real SSOT win** (Pattern 2+7). Collapsed 4 reader
  sites that hand-indexed raw `prov["git_branch"]`/`.get("git_commit")` into one
  vocabulary owned by the writer's module. Deferred to the 2nd/3rd reader (Pattern 6,
  per commit msg) — correct timing.
- **Tri-state staleness** (`None`=unknown ≠ `frozenset()`=verified-unchanged ≠
  non-empty) carried in `bool | None`, each state independently tested. Docstring
  forbids the collapse callers would make ("a commit this clone lacks is MORE
  suspect"). Pattern 4 in the small.
- **Performance comments state CONSTRAINTS not narration** — "per-type SQL predicate
  tempts the planner onto the huge type index, measured 15× slower"; "must not pay
  ~150ms graph-stack import". They explain why the obvious refactor is wrong.
- **Lazy-networkx split + additive `idx_node_attrs_key_value`** verified genuinely
  additive (no prior node_attrs index to retire; SCHEMA_VERSION stays 1, old DBs
  readable). No retirement obligation.

The flagged `brief._norm` vs `query.node_at._norm` "twin" = **JUSTIFIED BIFURCATION,
NOT a twin path.** Byte-identical bodies but opposite substrates: node_at normalizes
in NetworkX-space (graph already in memory), brief normalizes in SQL-space (builds
`json.dumps` literal spellings to push the match into the indexed `WHERE value IN`,
because never-load-the-graph IS its reason to exist). Unifying forces brief to load
the graph (defeats budget) — same "one contract, two substrates" shape as the SN
SI-vs-Krylov / apply-vs-residual fold. CONCERN-not-VIOLATION: bound it (reciprocal
cross-ref comment + a corner-probing equivalence test: symlink/mixed-sep through
BOTH), do NOT collapse (Pattern 6 / FaceField-deferral precedent).

Open nits (cheap, in-session, non-blocking): (1) the two `_norm` share only a
one-directional prose pointer + no equivalence pin → drift habitat when one hardens
(symlink/casefold) and the other doesn't → silent: ambient brief reports wrong file's
blast radius. (2) `brief.py`'s 2nd-tier basename `LIKE...ESCAPE` fallback (the
subtlest code — hand-rolled `\\`/`\%`/`\_` escaping + load-bearing trailing-`"`
anchor) has ZERO test coverage; it only runs on the corner inputs no exact-match test
hits. Both silent-failure shaped → defend with a test.

---

**TOKEN-BUDGET review (`feature/token-budgets`, single commit `53dda60`, vs parent,
2026-06-12). Verdict: PASS clean.** 438 tests green, pyright 0/0 on `sphinxcontrib/`.
The "context-bomb budgeting" carve: `assemble_context` gains `per_type_limit=25`
(per-edge-type buckets, most-connected-first, honest `omitted`+`hint` iff dropped);
NEW `assemble_impact` (`per_depth_limit=50`) single-sources the MCP+CLI impact path;
`_compact_node` strips falsy NodeResult fields. Motivation MEASURED (2.7 MB degree-3429
`numpy.array` context → 12 KB).

Durable RULINGS (a future session WILL re-litigate these — they are the four questions
the user posed):
- **`_compact_node` falsy-strip (`{k:v ... if v}`) is PRINCIPLED, not a hack** — and
  the discriminator is the TYPE. `NodeResult` (query.py:22) has `id:str` required +
  EVERY other field defaulting to a sentinel non-value (`degree=0`, `lineno=0` docstring
  "0 when unknown", rest `""`). The falsy values literally ARE the absence encoding, so
  `if v` drops exactly the absents and never a live value. SAFE ONLY because the type was
  designed sentinel-as-absence; would be a lossy hack against a type where 0/"" are live.
  Two guards make it sound: context sort reads `e.get("degree",0)` (degree may be
  stripped); `id` (secondary sort key + every `e["id"]` index) is required→never stripped.
- **The `X if X > 0 else None` "0=uncapped" mapping at 4 call sites is NOT a Pattern-7
  smell** — it follows the ESTABLISHED package boundary convention (pre-existing
  `processes`/`verification_coverage` use identical spelling at server.py:626,903 +
  cli.py:1150). Assemblers are `None`-native (honest semantic type); integer-zero is an
  MCP/argparse protocol affordance (int params can't carry clean `None`) converted to
  `None` once per boundary. Producer-side `None`, edge-side sentinel = the RIGHT split.
- **capped-by-default context/impact INVERTING the uncapped-by-default
  processes/verification contract is CORRECT** — "select narrow on MEASURED cost"; only
  these two are the measured bombs, docstrings+commit cite the bytes. Not arbitrary.
- **`omitted`/`hint` present IFF something dropped + `total_affected` read PRE-cap** =
  illegal-states-unrepresentable for a budget layer (truncation un-silenceable; small
  node carries no confusing empty `omitted:{}`; blast-radius can't lie). The single most
  important budget-layer property, built in not asserted.
- **cross-test-module import `from test_serialize import _build_hub_graph` in test_cli
  is the ANTI-twin choice** — the degree-spread builder is correctness-load-bearing
  (it's what makes "most-connected survives cap" meaningful); one builder + 2 consumers
  beats a drift-prone copy. Mild coupling < duplicated fixture risk.

Two micro-nits (record-and-move-on, NOT rework): (1) 2 near-identical hint strings —
leave at 2 instances (Pattern 6); single-source `_budget_hint(unit,limit)` IF a 3rd
budgeted tool appears. (2) the two sort keys differ textually (`e.get("degree",0)` over
compacted dict vs `n.degree` over raw NodeResult) for a REAL reason (pre- vs
post-compaction); a half-line comment at `_serialize.py:84` would inoculate against a
future "tidy" to `e["degree"]` that KeyErrors on zero-degree neighbors.

---

**LSP-PARITY ORACLE review (`feature/lsp-parity-oracle`, single commit `bf04373`, vs
parent, 2026-06-12). Verdict: PASS-with-nits.** 430 tests green, pyright 0/0 on
`sphinxcontrib/`. New `tests/test_lsp_parity.py` (393 lines) + a CI step. An
LSP↔graph structural-independence drift guard: pyright re-derives symbols + call edges
for a 2-file fixture project; the test asserts the AST analyzer's graph agrees. Would
auto-catch the v0.12.0 51%-worktree-contamination symbol-drift class. Two probes:
documentSymbol-vs-graph-def-nodes (set EQUALITY per file, granularity in the REDUCERS
not a weakened assertion) + incomingCalls-vs-graph-callers (static `{"top","unused"}`
exact; dynamic dispatch two-strength: same-class `self.run()` capability-pinned, annot-
mediated `s.run()` only pyright resolves → graph ⊆ pyright, the gap is Phase-F4 input).

What good looked like (reinforce):
- **Granularity tolerance encoded in the REDUCERS, not by weakening equality→subset**
  is the headline elegance win. `_lsp_def_symbols` strips closures (container ∈
  `function_names`) so Probe-1 stays set-EQUALITY. A subset assertion would silently
  tolerate a DROPPED top-level symbol — exactly the v0.12.0 bug class the oracle exists
  to catch. Keeping the assertion at equality + moving the known divergence into a named
  reducer is Pattern-4-in-tests: the tolerated difference is REPRESENTED (a line you can
  read), not ASSUMED-AWAY (a weaker operator).
- **The two dynamic-dispatch strengths are correctly typed as different CONTRACTS** —
  `self.run()` = capability pin (`assert "_step" in graph_callers`), `s.run()` =
  directional `graph ⊆ pyright`. NOT collapsed to one loose assertion. The gap is named
  as a deliverable (Phase-F4 pyright-enriched edges), not hidden.
- **Closures pinned FROM BOTH SIDES** (`test_graph_excludes_closure_pyright_sees`).
- **`LspClient` is correctly SIZED for test infra, and the no-dep call is RIGHT.**
  ~6 verbs (initialize/didOpen/documentSymbol/prepareCallHierarchy/incomingCalls/
  shutdown); single-outstanding-request synchronous pump = simplest correct shape (Hoare
  "obviously no deficiencies"). Three robustness details earn their place, NOT
  over-building: timeout-guarded reads (wedged server FAILS not HANGS), answering
  `workspace/configuration` with `[None]*len(items)` (spec-correct per-item count, pump
  can't deadlock), dropping notifications. Minimum viable correct client; pytest-lsp/pygls
  rejection sound (more surface to own than the thing under test).
- **Module-scoped fixtures + module-level skipif + CI-runs-only-in-pyright-job** = right
  cost structure; matrix jobs skip BY DESIGN (documented); oracle step after type-check
  warms the npm cache (stated constraint, not narration).

Nits (cheap, non-blocking, did NOT require rework):
1. **`q._g` private-attr access in two reducers (lines 283, 336) breaks this repo's
   test convention.** Every other test reaches the graph via PUBLIC
   `q.knowledge_graph.nxgraph` (test_query.py:666-673 pins the property; test_merge.py
   uses `kg.nxgraph` ~10×). `query.py:425` is `self._g = self._kg.nxgraph`, public
   `knowledge_graph` prop at :428 — `graph.knowledge_graph.nxgraph` is byte-equivalent,
   in-convention. CONCERN-not-VIOLATION (private access in tests is soft; bug habitat low
   — a `_g`→`_kg.nxgraph` rename breaks this module while the suite rides the property).
2. **Reducer's def-like triple `("function","method","class")` is a hardcoded
   string-literal DUPLICATE of a production SSOT.** `ast_analyzer.py:962`
   `_CANONICAL_TYPES = frozenset({CLASS,FUNCTION,METHOD,MODULE,EXCEPTION,TYPE})` from the
   `NodeType` enum (graph.py:12). CURRENTLY CORRECT for the AST-only fixture (verified:
   analyzer assigns `NodeType.CLASS` to EVERY `class` regardless of base; `EXCEPTION`/
   `TYPE` only from the Sphinx `:py:exc:` role / ID-prefix path a pure `analyze_directory`
   never hits). Latent Pattern-7: the day the analyzer base-inspects to assign
   `EXCEPTION`/`TYPE` to AST nodes, the reducer SILENTLY drops them — equality fails
   spuriously or the drift guard itself drifts. SAME `_FACE_NORMALS` hand-list-vs-enum
   shape as C5.3. Cheapest defense: `NodeType.{FUNCTION,METHOD,CLASS}.value` (rename
   caught by checker) OR assert the set ⊆ `_CANONICAL_TYPES` (vocab extension fails LOUD).

Granularity ruling (durable, cross-repo): set-EQUALITY-with-reducer-encoded-tolerance is
STRICTLY STRONGER than subset-and-hope, and is the correct default for any
structural-independence oracle. Reach for subset ONLY where the two views genuinely
disagree by CAPABILITY (the `s.run()` type-inference gap) — and there, NAME the gap as a
deliverable. A subset assertion papering over a granularity choice is anti-pattern #17
(loosened contract) in test clothing.
