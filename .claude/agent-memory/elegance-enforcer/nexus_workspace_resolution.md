---
name: nexus-workspace-resolution
description: sphinxcontrib-nexus (sibling repo) workspace/worktree wiring review — _switch_workspace single-source, resolve_checkout_root 3-form order, KnowledgeGraph(graph) wrap, briefing payload share
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
