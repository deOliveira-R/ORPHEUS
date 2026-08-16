# A graph query returns what an agent needs to ACT, not everything the graph knows

**Status: PARTLY APPROVED.** Sections marked *proposed* are still hypotheses;
sections marked `[M]` carry the command or query that produced them.

## ⏸ COMPACTION POINT #3 — 2026-08-16, Track 0 COMPLETE (0.1–0.6)

**Track 0 has landed**; Track 1 has not started. Work is on branch
`feat/config-and-ontology` in `~/git/sphinxcontrib-nexus` (ours — a folder
change, not a hand-off) and `chore/nexus-project-config` in ORPHEUS. `main` is
untouched in both. Terminal commits: nexus `c51672c`, ORPHEUS `4c3c46b3`.
⚠ Trust `git merge-base --is-ancestor`, not this line.

⚠ **Every graph number in this plan is an instance count from one build and
moves on every rebuild.** The before/after pairs in the resume table are the
exception — each was measured across a real rebuild and is a *delta*, which is
what makes it re-checkable.

⭐ **The graph store MOVED on 2026-08-16 (Track 0.6).** It is now
`<project root>/.nexus/graph.db`, derived — never declared. Any older text in
this plan naming `docs/_build/html/{_nexus,graph}/graph.db` is describing a
path that no longer exists; ask `nexus config db`. `.nexus/traces/` holds the
runtime overlays and is DURABLE — it is not in the build tree precisely so a
clean build cannot delete it.

**Environment, verified after the 0.6 landing:** the MCP server was restarted
and now reports `db = <root>/.nexus/graph.db`, `has_graph: true`, provenance
stamped at `f81e3258`, and **23013 nodes / 206868 edges** — matching disk. If a
future session sees the server's counts disagree with `nexus status`, it is
serving a frozen snapshot; reconnect it, `use_workspace` cannot self-heal that.

<details><summary>superseded — compaction points #1 and #2</summary>

**#1, "before any code is written":** Nothing in this plan has been
implemented. What exists is measurement, four filed issues, and two user
decisions.

**#2, "Track 0 COMPLETE":** said 0.1–0.5. Track 0.6 landed the same day, from
a question asked at the Track 1.1 design checkpoint — see the 0.6 row.

</details>

### User decisions, 2026-08-16

| decision | status |
|---|---|
| ⭐ **The nexus RESHAPE goes FIRST** — before #68 and before `layer` (§6ter) | ✅ **APPROVED 2026-08-16**, supersedes the ordering below |
| Fix the prose-symbol leak (**#68**) | ✅ approved — **folded into Track 0.4**, no longer standalone |
| `layer` as a node property (§5.4 / §6.1) | ✅ **APPROVED** — now behind the reshape |
| §5.1 retype test→code `calls` as `exercises` | ✅ **ACCEPTED** |
| §5.3 observed state is a **DATABASE**, not a JSON blob | ✅ **CHOSEN** — only the CI story is deferred |
| §5.2 halo demotion | *proposed*, not yet ruled on |

⚠ I recorded the first two of these wrongly on 2026-08-16 and the user
corrected both: `exercises` had been accepted ("a good suggestion"), and §5.3's
deferral was only ever about **CI**, not about DB-vs-JSON. Kept visible per
`plan-authoring` §3 — a plan that silently repairs its own misreadings teaches
the next reader nothing about how they happen.

### ✅ Track 1.1 LANDED 2026-08-16 — a query knows which working tree it answers about

Two commits on `feat/config-and-ontology`, both mutation-verified per arm:

| | outcome | `[M]` |
|---|---|---|
| **1.1a** `35d5d31` | the binding — `GraphQuery` carries its `Workspace` | `project_root` off all 7 public signatures; server globals **2 → 1**; `_changed_cache` retired (the query IS the key); `graph_db_in` single-sources the 0.6 convention; suite **777 → 786** |
| **1.1b** `8de24ec` | every returned position says when it may be wrong | staleness applied by **1 of 40 tools → the `@nexus_tool` boundary, once**; fresh-graph payload byte-identical (identity, not equality); **0.00 ms** fresh vs **35.2 ms** on a 2059 KiB payload with 30 files dirty; suite **786 → 790** |

`[M]` ORPHEUS rebuild after the carve: **23013 / 206868** — unchanged from
the 0.4/0.6 baseline, so it is content-neutral by measurement, not by
argument. Sphinx clean under `-W`.

✅ **VERIFIED ON THE LIVE MCP SERVER 2026-08-16, both directions.** A green
reading could not have discriminated a working pass from a dead one
(vv #19), so a positive control was forced. Same server, same edit, the
only variable being the age of the query object:

| step | result |
|---|---|
| `session_briefing`, clean tree | no flags — and this call CACHES the changed-set |
| edit `orpheus/numerics/measure.py`, then `query("DiscreteMeasure")` | ⛔ **no `stale` key** on any of 3 |
| `touch .nexus/graph.db` (content unchanged → new `GraphQuery`), same edit | ✅ all 3 flagged, naming commit `7593db5e` |
| restore the file, `touch` again | ✅ flags gone |

⟹ the pass is alive and per-file; what is frozen is the **changed-set**.
`files_changed_since_build` is cached per query object, and a query object
outlives working-tree edits — so in a live session the set freezes at the
FIRST tool call. **Filed as nexus #70**, with the measured cost of the
obvious fix (`changed_files` is `[M]` 25.1 ms on ORPHEUS, 1.00 s across a
40-call session if run unconditionally).

⚠ **NOT a regression** — the retired `_changed_cache` key was
`(root, db_mtime, commit)`, which does not include working-tree state
either. What changed is reach: 1 tool → 40. ⚠ And note the existing gates
**cannot** see it: every one builds a fresh query after the edit, which is
exactly the configuration that hides it. Same shape as the byte-identity
blindness above, one layer out — [[lessons-L58]] earned its keep twice in
one step.

⚠ **A reconnect is needed after ANY nexus commit** — the server serves the
code it imported at startup, and `use_workspace` does not re-import.

⭐ **Two things this step got wrong first, both caught by instruments and
both worth carrying:**

1. **The subdirectory control passed against the configuration it exists
   to condemn.** `rich_graph` spells `file_path` *relatively*
   (`"alpha.py"`), which happens to equal git's repository-relative
   output, so its nodes match whatever root they are asked about. `[M]`
   **16527 of 16527** `file_path` values in ORPHEUS's real graph are
   ABSOLUTE — 0 relative. The fixture was blind to the very axis the
   test varies (`vv-principles` #24(d)/(e); the configuration includes
   *what the fixture cannot see*).
2. **The byte-identity done-when had NO WITNESS**, and the first mutation
   battery said so — GREEN on two arms. Identity cannot distinguish "did
   not parse" from "parsed and returned the original", and the
   re-serialisation branch is unreachable from the fresh case. Promoted
   to `vv-principles` **#26**: gate a claim about the PATH by
   instrumenting the path, never by asserting the output.

**Track 1.2 is LANDED** (see its block below). Next is **`NodeId`**, whose
evidence is ✅ **verified, mine** — but it was measured *before* Track 0.3
landed, and 0.3's own outcome line is "one symbol, one id, on both producers".
⟹ **reconcile the NodeId evidence against the tree before designing**
(`plan-authoring` §3's REMEDIED-fact clause: a fact can die by being FIXED, and
this campaign is the thing that fixed it). `Evidence` and `Diagnostic` remain
⬜ **relayed, not reproduced**.

> ⚠ **This paragraph said "`NodeId` … remain ⬜ relayed" until 2026-08-16 and
> that was wrong** — the evidence table at `:1306` marks R2 ✅ verified with my
> own measurement. Written in the same session that landed 1.2, i.e. by an
> author who had just read the table. Kept per §3: the failure mode is a
> *summary line* drifting from the table it summarises, which is what every
> resume pointer is.

> ⛔ **This heading read *"Track 1.1 — `ProjectView`, graph + working tree as
> one object"* until 2026-08-16, and the mechanism is REFUTED.** Do not mint
> `ProjectView`. `workspace.Workspace` is already documented, verbatim, as
> *"One checkout paired with its graph database"* — the identical concept, so a
> new wrapper would be the **fourth** spelling in this file family
> ([[lessons-L57]]) and would leave `project_root` still travelling per-call.
> By the concept-count test it ADDS a layer. The goal in the title is the
> durable half and it stands; the means is to **bind the type that exists**.

**The design, decided 2026-08-16 and measured, not proposed.** ✅ **BUILT AS
WRITTEN** at `35d5d31`/`8de24ec` — every FACT below describes the tree
*before* that landing and is now past tense. Kept per `plan-authoring` §3:
the rows are what the design was reasoned from, and the ones marked
✅ REMEDIED are the ones the campaign repealed rather than refuted.

`GraphQuery.__init__(self, graph)` takes no workspace (`query.py:711`), so the
server holds `_query` and `_workspace` as **two module globals for one
object** — and they have already disagreed in production: the server answered
from a snapshot while `_workspace.db_path` pointed at a deleted file, so
`_reload_if_stale()` failed its `stat`, logged at `debug`, returned, and the
server could never reload while reporting `has_graph: false`. That is this
step's motivating defect, observed live rather than argued.

⟹ `GraphQuery` carries the `Workspace` it was loaded from; `project_root`
leaves the public signatures; `_active_root()` becomes a read off the query.

**The precedent, with every adjective checked at HEAD** (§1's PRECEDENT
clause — a cited shape is what drifts):

| claim about `Workspace` | verified |
|---|---|
| `@dataclass(frozen=True)`, `workspace.py:318` | ✅ `w.root = …` raises `FrozenInstanceError` |
| two fields: `db_path: Path` required, `root: Path \| None = None` | ✅ `[f.name for f in fields(w)] == ['db_path', 'root']` |
| carries `db_relpath` + `sibling()` (the worktree transplant) | ✅ both present and working |
| `root` is OPTIONAL by design — a bare `--db` server has a graph but no tree | ✅ its own docstring says why |

⚠ That last row is the constraint on the binding: `workspace=None` must stay
legal, because **176 of the 178 `GraphQuery(...)` construction sites are in
tests** and most pass a bare graph. The parameter is additive or the change is
a 176-site migration.

⭐ A fresh precedent for the *shape* of the repair landed in 0.6:
`RuntimeStore.beside(db_path)` — the CLI and the server each had their own
`db_path.parent / "traces"`, and the fix was one home both read.

**Evidence re-derived at HEAD 2026-08-16, with the FILTER beside each
fraction** (a fraction without its predicate is not re-runnable — the
`10 of 56` below read as refuted last compaction for exactly that reason):

| claim | at HEAD | how to re-run |
|---|---|---|
| `ProjectView` does not exist | **0 hits** — ✅ and it never will; `Workspace` was bound instead | `grep -rn ProjectView sphinxcontrib/ tests/` |
| `GraphQuery` methods taking `project_root` | **10 of 56** — 7 public (`detect_changes`, `node_at`, `rename`, `retest`, `session_briefing`, `staleness`, `verification_audit`) + 3 private (`_apply_renames`, `_git_changed_files`, `_git_file_timestamps`) — ✅ **REMEDIED 2026-08-16 @ `35d5d31`**: all 7 public dropped it; the 3 private stay, they are pure functions of a root and take `self.project_root` | AST-walk `query.py`, arg name exactly `project_root` |
| MCP tools | **40** carry `@nexus_tool` | AST-walk `server.py` decorators |
| …that consult staleness | **1** — `_position_staleness_warning`, defined `server.py:217`, called at `:428` only, inside `node_at` — ✅ **REMEDIED 2026-08-16 @ `8de24ec`**: that helper is RETIRED; `_mark_stale_positions` runs in the `@nexus_tool` wrapper, so the count is the boundary itself | `grep -n _mark_stale_positions server.py` |
| …that take a position as INPUT | **3** — `node_at`, `ingest`, `runtime_ingest` | AST-walk for an arg in `{file, file_path, artifact}` |
| producers EMITTING `(file_path, lineno)` | **51** `_node_result` call sites | `grep -c _node_result query.py` |
| `GraphQuery(...)` construction sites | **178** — 176 tests, 3 `server.py`, 1 `cli.py`. This is why the parameter had to be ADDITIVE, and it held: ✅ 2026-08-16, **9** sites migrated, the rest untouched | AST-walk the repo for `GraphQuery(` calls |

⛔ **The step's own done-when was REFUTED before design began.** §6ter said
*"staleness is applied by construction, and the count is 40 of 40."* Only
**3 of 40** tools take a position as input, so 37 have nothing for a
position-staleness warning to apply to; the number silently changed what it
counted, from tools that APPLY it to tools that SHOULD. The real exposure is
the **51** emitters, whose `(file_path, lineno)` `NodeResult`'s own docstring
invites feeding straight to `Read`. Logged as a `plan-authoring` surprise.

✅ **USER RULING 2026-08-16 on the replacement scope:** bind, then attach
staleness **once at the `@nexus_tool` boundary** — the wrapper scans the
outgoing payload and appends a line only when a returned file is genuinely
stale. One site, not 51; absent entirely when the graph is fresh, so it costs
nothing against **#67**. New done-when: *every tool that returns a position
says so when the position is stale, and the fresh-graph payload is byte-identical
to today's.*

⚠ **A control this step needs, and it is not a mutation.** `[M]` at the time
of writing **11564 of 11564** stored positions still land on their own
`def`/`class`, because the graph is fresh — so a green staleness gate proves
nothing until a file is *made* stale. The positive leg is worthless here
(vv #19); the gate must be built around an edited file. Note also the one
existing consumer had been **inert** on this project: `_files_changed_since_build`
returns `None` without a provenance stamp, and the server was reporting none.

✅ **HEEDED, and it was not enough** (2026-08-16 @ `8de24ec`). Every gate is
built around an edited file, as instructed — and the mutation battery still
came back GREEN on two arms, because the *byte-identity* half of the
done-when is a claim about the FRESH path, where an edited file cannot
reach. ⟹ the caution above covers the flag's presence; the promise that it
costs nothing needed its own instrument (a parse counter) and its own
configuration (dirty tree, unaffected file). `vv-principles` **#26**.

### ▶ Track 1.2 `PositionIndex` — the gate is SATISFIED, and it moved the step

`[M]` 2026-08-16, `scratchpad/probe_position_index.py` + `…_theft_exact.py`,
against ORPHEUS's real graph (`23013 / 206868`). ⚠ The relayed measurement
never recorded WHICH 4 positions, so this is a **re-measurement, not a
replication** — positions were derived from the AST of
`orpheus/numerics/measure.py`, not copied.

| probe | line | `query.node_at` (A) | `resolve_node` (B) | |
|---|---|---|---|---|
| inside a method body | 267 | `DiscreteMeasure.__post_init__` | `DiscreteMeasure.__post_init__` | ✅ |
| class body, no method | 144 | `DiscreteMeasure` (class) | `None` | ≠ |
| module scope (import) | 75 | `orpheus.numerics.measure` | `None` | ≠ |
| decorator line | 288 | `DiscreteMeasure` (class) | `DiscreteMeasure.dim` | ≠ |

⟹ **3 of 4 — the relayed count reproduces on independently-chosen positions.**

⛔ **But the relayed FRAMING is wrong, and it would have misdirected the
design.** It is not "three implementations of `(file, line) → node`". It is
**two concepts, each duplicated**, and the third cited site answers a
different question entirely:

| concept | spellings | agree today? |
|---|---|---|
| **path → comparable key** | 3 — `node_at._norm`, `_in_file_node_ids._norm` (its own comment says *"keep the two in lockstep"*), `NodeBinder._abs` — ✅ **REMEDIED 2026-08-16 @ `7db466d`**: one `workspace.canonical_path` | ✅ yes, by hand |
| **span search** | 2 — A (all types, smallest span, module fallback) vs B (function/method only, latest start, decorator window, no fallback) — ✅ **REMEDIED 2026-08-16 @ `c497ddb`**: one `PositionIndex._innermost`, two verbs over it | ⛔ no |
| `brief._in_file_node_ids` | is `(file) → nodes`, **not** a position lookup | `[M]` **56 = 56**, set-equal to A's `_norm` match |

`[M]` B returns `None` where A and C succeed when the query is spelled
**relatively** — B's index keys are the raw stored paths, so its
normalization lives in its *caller*. Invisible on this graph
(`[M]` **16527 of 16527** stored `file_path` are absolute), which is why it
needs saying rather than probing.

⛔ **The step's own done-when is REFUTED, the same way Track 1.1's was.**
*"one implementation, and the 4 probed positions agree"* — two of the three
disagreements are B's **deliberate design**: classes are excluded because they
would shadow their own methods, and a trace record has no module to bind to.
Forcing agreement would break the runtime join. What must be single-sourced is
the **key space and the span data**, not the resolution policy.

⭐⭐ **And the third disagreement is a LIVE DEFECT, in the tool this campaign
runs on.** `resolve_node` writes the window and the body test as one
condition — `ln - WINDOW <= ask <= end` — so a *later sibling* can claim a
line above an earlier def, and "latest start wins" hands it the answer:

```
ask = 288 (the decorator of n_points)
   289-291  n_points   <- the decorator's OWN def
   294-301  dim        <- 294-8 <= 288 <= 301, scanned later, WINS
```

`[M]` cProfile reports `co_firstlineno`, which for a decorated def is the
**first decorator line** — so this is the production path, not a corner.
Measured over the whole graph, with the predicate stated (`plan-authoring` §2):

| population | count |
|---|---|
| indexed `function`/`method` nodes | 10207 |
| …of which **decorated** (by AST, own def line indexed) | 3530 |
| …… bound to the **right** node | 2969 |
| …… bound to the **wrong** node | **456** — 12.9 % of decorated |
| …… bound to nothing | 105 |

291 of the 456 are `@property` — the exact node class `runtime_hotspots`
advertises (*"a property called 10k×/run = a caching opportunity"*). The theft
is always to the **next sibling**, so a run's counts are shifted one property
down the file: plausible numbers on the wrong node, i.e. [[lessons-L57]]'s
invisible-wrong-answer again.

⚠ **Its gate is green and cannot fail — the fixture is blind to the only axis
that matters.** `test_resolve_decorator_line_above_def` spaces its defs at
**10, 30, 52** (gaps of 20 and 22) against `DECORATOR_WINDOW = 8`, so no
sibling can ever reach another's decorator line. Real `@property` blocks sit
5 lines apart (`mixture.py`: 108/113/118/123). Third blind fixture in this
campaign, after `rich_graph`'s relative paths and the byte-identity witness —
`vv-principles` #24(d)/(e).

**Replacement done-when** (the refuted one is struck above):
1. **one** path→key implementation — `grep -c "def _norm" sphinxcontrib/` is 0;
2. **two** named span verbs, and no third — navigation vs trace-binding, each
   saying which question it answers;
3. `[M]` the 456 misbindings go to **0**, and the decorator fixture is re-spaced
   so the existing gate *can* fail (`plan-authoring` §6c — the gate lands with
   the case it catches).

### ✅ Track 1.2 LANDED 2026-08-16 — all three met

| | outcome | `[M]` |
|---|---|---|
| **1.2a** `7db466d` | one contract decides when two paths are the same file | 3 private realizations → `workspace.canonical_path`; `grep "def _norm"` **0**; 8 hand-written contract gates replace the demoted copy-comparison; suite **790 → 798** |
| **1.2b** `c497ddb` | a position resolves one way, and decorators land on their own def | misbindings **456 → 0**, unbound **105 → 6**; `decorator_lineno` recorded at the producer; graph **23013 / 206868** unchanged; suite **798 → 803** |

`[M]` both joins run against the **same rebuilt graph**, so the comparison
isolates the join, not the graph:

| join | RIGHT | WRONG | NONE |
|---|---|---|---|
| the retired `resolve_node` | 2969 | **456** | 105 |
| `PositionIndex.defined_at` | **3524** | **0** | 6 |

⭐ **The repair was at the PRODUCER, not in a better heuristic.** The analyzer
already walked `decorator_list` and kept only the rendered *names* — it knew
where the definition started and threw the number away, so the consumer had to
guess with an 8-line window. This is [[feedback-lossy-return-type-is-the-root-cause]]
exactly: triage one hop UP, and the guess becomes an exact match.

⭐⭐ **Two duplications were found by MEASURING rather than arguing, and both
were mine, written the same hour:**
1. the exact-match-first pass I put in `defined_at` — `[M]` **0 of 1 830 000**
   realizable positions changed answer, so it was a copy of the search beneath
   it. Deleted, with the measurement in the docstring so it is not "restored";
2. the two verbs' innermost rules (smallest extent vs latest start) — `[M]`
   agree on **all 1 830 000**; they can differ only for extents that overlap
   without nesting, which no source file produces. One `_innermost`.

⚠ **The blind-arm count was 2 of 6, and both were the FIXTURE again** — the
third and fourth instances in this campaign. Every graph fixture stored
ABSOLUTE paths, so mutating the index's key canonicalisation was a *no-op*;
and no fixture had a decorator stack longer than the window, so removing
`decorator_lineno` was *rescued by the fallback*. Witnesses built for both. The
pattern is now unmistakable enough to state as a rule: **a fixture inherits the
one spelling its author had in front of them, and that is the axis the gate
cannot see.**

⚠ **The old decorator gate could not fail.** Its fixture spaced defs at
10/30/52 — gaps of 20 and 22 against `DECORATOR_WINDOW = 8` — so no sibling
could reach another's decorator lines. Re-spaced to real `@property` spacing.

⭐ **Filed, not fixed here: nexus #71.** The residual 6 unbound are all
`source = both` nodes, and the merge drops their whole decorator block —
costing `[M]` **36 of 420** `@pytest.mark.verifies` declarations their `tests`
edge, i.e. **8.6 % of ORPHEUS's declared V&V links are missing from the
coverage matrix**. Error direction is the safer one (false MISSING, not false
COVERED), but `verification_audit` is exactly the tool an auditor trusts.
Separate defect, separate blast radius, own issue.

**Landed 2026-08-16, each mutation-verified, all on `feat/config-and-ontology`:**

| | outcome | `[M]` on ORPHEUS's rebuilt graph |
|---|---|---|
| **0.2** `3e137ff` | the graph never claims a call the source does not make | `calls` pairs **53117 → 51376** (−1741, **0 added**); 291 of the removed pointed at a real symbol; fabricated self-loops **57 → 5**; 94 real symbols had no un-fabricated incoming call; `dead_functions` candidates **2919 → 2967** |
| **0.3** `2ddf61b` | one symbol, one id, on both producers | short-prefix ids **316 → 0**; symbols carrying both spellings **263 → 0**; nodes 24214 → 23971 as the halves merged |
| **0.4** `090a793` | the equation namespace holds declared labels, not LaTeX | `math:equation:*` **1860 → 903**, all 903 declared (was 51 % LaTeX); newline ids **13 → 0**; `dead_references` **0 → 0**; nodes added **0** |
| **0.5** `7882dba` | a project may WIDEN a base edge, never narrow it | `[extend.edge.X]` unions `domain`/`range`/`sources`/`attributes`; scalars and `forbid_source_attr` refused; monotonicity gated over every node type |
| **0.6** `c51672c` | the graph lives with the project, not with the build output | store → `<root>/.nexus/`; `[graph].db` RETIRED (derivable); `traces/` out of the directory a clean build empties; `find_project_root` bounded at `$HOME` + the checkout; `DEFAULT_OUTPUT` single-sourced. `[M]` rebuild **23013 / 206868**, identical to the 0.4 baseline, 5 runtime runs still binding; suite **777 passed / 1 skipped** |

Net over 0.2–0.4: nodes **24307 → 23013**, edges **215226 → 206868**, and every
removed item was something the graph had asserted and could not support. 0.6
moved the store without touching content — the rebuild lands on **exactly**
that 23013 / 206868, which is what makes it a re-runnable check rather than a
claim. Sphinx clean throughout; nexus suite **777 passed / 1 skipped** (was 773
before 0.6 added four boundary gates).

⚠ 0.5 ships the **mechanism only** — the extension tier still has zero
production consumers. Its first are §5.1 `exercises` and §5.4 `layer`.

⚠ Those 94 symbols are **not** thereby dead — most are protocol-typed methods
reached by dynamic dispatch, which the static graph cannot see either way (the
co-execution ruling). What changed is that it stopped presenting invention as
evidence: a candidate list may be wrong, a call edge may not.

**Read the ⛔⛔ block under the Track 0 table before sizing anything.** All three
of 0.2/0.3/0.4 refuted their stated size the same way, and 0.4 refuted this
plan's *description* of its own defect as well.

> ⛔ This line read *"Next is 0.5 (#69), then Track 1"* until 2026-08-16. 0.5
> landed the same day, at `7882dba`, two rows above in this very table — the
> `plan-authoring` §1 pointer rot, inside a compaction point whose own table
> already contradicted it. Kept because the shape is the point: a pointer and
> the table it sits beside keep different clocks, and the pointer is the half
> that is read.

<details><summary>superseded pointer — Track 0.2 (kept per <code>plan-authoring</code> §3)</summary>

⚠ **The sequencing CHANGED 2026-08-16 and this pointer supersedes the earlier
one.** #68 was first; the **nexus reshape is now first** (user: *"better to do
this first so we can build on top of better foundations"*). #68 has been folded
into **Track 0.4**, because it and the newline-bearing ids are the same defect at
the same producer. **Read §6ter for the full order and the critique of the
order the review proposed** — I do not accept that order and the reasons are
measured.

**Track 0.1 is ✅ LANDED** (nexus `37bad88`, bit-identical, 195 refusals).

**Next is Track 0.2.** *Goal (outcome):* the graph never claims a call that the
source does not make. `ast_analyzer._unparse_attribute` fabricates `calls` edges
— `get_thing().method()` mints an edge to a bare `method` — and it fails in the
false-**alive** direction, so it inflates `impact`/`retest` and *rescues*
genuinely dead functions from `dead_functions`. `[M]` 8 of 13 incoming `calls`
on one function in nexus's own graph are fabricated. Reported as a one-line fix;
**verify that before believing it**, and check what the suite pins — I have not
re-measured this one myself.

Then 0.3 (apply `py_type_map` on the xref path, 2 lines), 0.4 (id
normalisation = #68 + the 13 newline ids), then Track 1 in the stated order.

⛔ **REFUTED 2026-08-16, and the caution above is what caught it.** The
mechanism was right; *"a one-line fix"* was wrong, and so was "2 lines" for
0.3 — see the ⛔⛔ block under the Track 0 table. The relayed *"8 of 13"* was
never reproduced as stated: it was measured on **nexus's own** graph, and my
re-measurement was on **ORPHEUS's**, where the worst single case is 6 of 8
(`_CompositeLeaf.values`). Not a contradiction — a different denominator — but
it is why the number is not quoted forward.

</details>

Work is in `~/git/sphinxcontrib-nexus` (ours; a folder change, not a hand-off),
branch `feat/config-and-ontology`.

⚠ `[M]` **run the nexus suite with the NEXUS venv**, not ORPHEUS's: the latter
lacks the optional `sphinx-proof` and `test_proof_relations.py` importorskips at
module level, dropping **36 tests** from collection — a green 691 that looks
identical to a green 727.

### Where the work is tracked

⚠ Every `#N` below is a **nexus** issue (`deOliveira-R/sphinxcontrib-nexus`),
never an ORPHEUS one — the two number spaces collide.

**CLOSED by this campaign:** **#56** (`5796c6d`), **#64** (the config file —
`.nexus/{config,ontology}.toml`, plus the `nexus config` verb), **#68** (the
prose leak — `090a793`), **#69** (the ontology `extend` verb — `7882dba`).

**Still open, and the ones this plan will reach:** **#67** (payload is 5× its
information content — §5.5), **#55** (call-resolution fragmentation — Track 1's
`NodeId`), **#59** (an empty result is indistinguishable from a broken one),
**#57** (per-test attribution), **#63** (`catches` is a string, not an edge —
wants the ontology extension tier 0.5 just built).

**Open and adjacent, not scheduled here:** **#65** (the tool reference documents
29 of 40 tools), **#66** (`--db` declared 42× per-subparser), **#16**
(annotation-mediated dispatch), **#58** (`bridges`/`communities` do not
complete), **#60** (a cone needs an antisymmetric relation), **#61** (markers
from a collection manifest), **#62** (`retest` truncates at depth 3), **#20**
(filtered subgraph → graphviz render).

**Started 2026-08-16** from the question *"is a query a context bomb because of
the tests, or naturally?"* — and the measured answer refuted the premise of the
first design that followed from it. Both are kept, per `plan-authoring` §3.

---

## 1. Goal (outcome, not mechanism)

An agent asks the graph a structural question and gets back **the context it
needs to act**, at a cost proportional to the question. The graph's job is to
answer what grep cannot — *relationships* — and it must do so in a form that is
cheaper to read than the code it describes. Today a hub query returns megabytes
and the agent pays to learn mostly nothing.

**This is an ERGONOMICS goal, not a size goal.** Slimming is a means. The test is
not "how many bytes" but *"did the answer change what the agent did next?"* A
100 kB payload that surfaces the one caller that matters beats a 10 kB payload
that surfaces twenty-five arbitrary ones.

**Done when:** the routine structural questions (*who calls this, what does this
depend on, what tests it, what does it implement*) are answerable at a cost an
agent will pay **by default**, without a human first knowing which flag to pass.

⚠ The evaluation process is part of the work (user, 2026-08-16). Design changes
here get judged by whether an agent's actual workflow improves, not by the
byte counts alone.

---

## 2. Measured baseline — `[M]` 2026-08-16, pre-Track-0

Graph: **24,307 nodes / 215,226 edges**. Queries below are against that build;
re-measure after any rebuild, since these are instance counts, not schema facts.

⚠ Measured at `docs/_build/html/graph/graph.db`, **a path that no longer
exists** — Track 0.6 moved the store to `<root>/.nexus/graph.db`. Kept as the
provenance of these numbers, not as a location: ask `nexus config db`. The
current graph is **23013 / 206868**; every number in this section predates
Track 0's removals and is a *starting* count, not a description of the tree.

### 2.1 Node and edge population `[M]`

| node type | n | | edge type | n |
|---|---:|---|---|---:|
| function | 5,632 | | `calls` | 121,788 |
| **unresolved** | **5,169** | | `contains` | 26,706 |
| method | 4,575 | | `references` | 17,959 |
| attribute | 2,387 | | `type_uses` | 17,901 |
| data | 1,676 | | `implements` | 13,968 |
| class | 1,359 | | `imports` | 9,011 |
| module | 919 | | `documents` | 4,031 |
| equation | 903 | | `tests` | 2,748 |
| external | 746 | | `equation_ref` | 543 |
| section | 680 | | `inherits` | 326 |
| file, tag, term, exception | 261 | | `discriminates_on` | 245 |

`is_test` is stored **only when true** — 7,305 nodes (30.1%). Absence of the
attribute means "not a test", not "unknown".

### 2.2 The layer matrix — `[M]` this REFUTED the 3-way-split proposal

Layers assigned as: `ext` = type ∈ {external, unresolved}; `docs` = type ∈
{file, section, equation, term}; `test` = `is_test`; `code` = the rest.

| | → code | → test | → docs | → ext | row total |
|---|---:|---:|---:|---:|---:|
| **code →** | 40,037 | 9,387 | 14,134 | 48,815 | 112,373 |
| **test →** | 29,440 | 6,483 | 2,799 | 56,093 | 94,815 |
| **docs →** | 3,714 | 7 | 3,436 | 881 | 8,038 |
| **ext →** | 0 | 0 | 0 | 0 | **0** |

- **Intra-layer 23.2% / inter-layer 76.8%.**
- ⛔ **REFUTED 2026-08-16 — "split the graph into docs / code / tests".** Three of
  every four edges cross a boundary, so nearly every query becomes a cross-part
  join. The seam is real as a *concept* and wrong as a *partition*.
- ⭐ **The `ext` row is all zeros** — `[M]` external/unresolved nodes have
  **0 outgoing edges**. The halo is strictly TERMINAL: nothing traverses
  *through* it; it only inflates degree.

```
full : 24,307 nodes  215,226 edges
core : 18,392 nodes  109,437 edges     -24% nodes, -49% edges
```
Dropping the halo is **provably lossless for reachability**, because a terminal
node cannot lie on any path except as its endpoint.

### 2.3 Layers are ALREADY encoded in the edge type — with exactly one exception `[M]`

Excluding the halo, per edge type, the fraction crossing a layer boundary:

| edge type | core edges | % cross | dominant pair |
|---|---:|---:|---|
| `implements` | 13,968 | **100%** | code→docs |
| `documents` | 3,378 | **100%** | docs→code |
| `tests` | 2,748 | **100%** | test→docs |
| `equation_ref` | 543 | 0% | docs→docs |
| `imports` | 5,731 | 0.0% | code→code |
| `type_uses` | 4,134 | 5.6% | code→code |
| `references` | 12,320 | 14.8% | code→code |
| `inherits` | 230 | 14.8% | code→code |
| `discriminates_on` | 245 | 22.0% | code→code |
| `contains` | 26,706 | 35.4% | code→code |
| **`calls`** | **39,434** | **70.5%** | **test→code (27,647)** |

⭐⭐ **Every inter-layer relation already has its own edge type, except test→code
`calls`.** The layer structure is ~90% readable off the edge type today. This is
the finding the whole design rests on.

### 2.4 Payload cost `[M]`

`neighbors()` emits, per neighbour, a 9-field node dict **plus** a 4-field edge
dict: **488–572 B/pair**, stable across nodes (`py:module:orpheus` 488, `BC` 570,
`solve_sn` 572, `LinearOperator` 511).

`solve_sn`, 374 neighbours:

| | bytes | vs full |
|---|---:|---:|
| as shipped | 213,962 | — |
| test neighbours removed | 83,734 | −61% |
| **slim node repr, tests INCLUDED** | **46,402** | **−78%** |

Redundant *in this context*: `source`/`target` (one is the query node, the other
is the neighbour's own `id`), `key` (internal multigraph key), `display_name`
(usually `== name`), and `domain`/`docname`/`file_path`/`lineno` (empty for every
external node — costing bytes exactly where they carry nothing).

Across the **172** ORPHEUS production nodes with degree ≥ 100 (38,719 neighbour
slots, ≈18 MB at measured density): test neighbours **27.1%**, `calls`
(production) 20.7%, `implements` 13.9%, `imports` 13.8%, rest 24.5%.

⟹ **Tests are ~27% of the bomb.** Removing every one leaves 73%. Filed as **#67**.

---

## 3. What the 5,169 phantoms are — `[M]` 2026-08-16

26,582 edges point at `type='unresolved'` nodes (12.4% of the graph). They are
**not one defect**. Split by dominant inbound edge type:

| origin | nodes | what they are |
|---|---:|---|
| `calls` | 3,579 | call-resolution failures — see below |
| `references` | 1,282 | **prose, not code** |
| `documents` | 307 | doc-side reference failures |
| `imports` | 1 | — |

### 3.1 The `calls` phantoms — four distinct causes

The user's premise — *"in principle the AST should resolve packages ORPHEUS uses,
such as numpy"* — is **correct, and nexus already does it**: `numpy.array`,
`numpy.zeros`, `float`, `pytest.raises` all resolve to real `external` nodes
(746 nodes absorbing 79,207 edges). Resolution fails only where the callee is
reached **through a value whose type the AST does not know**:

| cause | example (inbound calls) | size |
|---|---|---|
| **local-variable receiver** | `rng.standard_normal` (640), `rng.uniform` (577), `rng.random` (215), `op.apply` (265), `monkeypatch.setattr` (278) | dominant |
| **bare method name**, receiver unknown | `apply` (312), `ravel` (205), `to_flat` (184), `realize` (178), `reshape` (170) | large |
| **inherited member — MRO not walked** | `AngularBoundaryFlux.zeros_on` (252) resolves to `_bases.FaceField.zeros_on`; `AngularFlux.from_mesh` (178) to `_bases.AngularField.from_mesh` | `[M]` 74 nodes / 942 edges |
| **type-prefix mismatch** (`py:function:` vs `py:method:` vs `py:meth:`) | `zeros_on` exists under 4 spellings, all unresolved | `[M]` 218 nodes / 359 edges |

`rng = np.random.default_rng()` then `rng.uniform(...)` needs to know the **type
of `rng`** — that is type inference, not AST walking. Nexus already names it:
**#16** (Pyright-enriched call edges). **#55** is the fragmentation half (one
method scattered across receiver spellings).

⚠ **The MRO and prefix causes are cheap and were NOT what I first guessed.**
`[M]` I hypothesised a constructor-call/type-prefix collision
(`py:function:X` vs `py:class:X`) and it is **refuted** —
`py:function:…AngularBoundaryFlux` does not exist at all, and the prefix class
totals only 359 edges. The real cause was inheritance.

### 3.2 The `references` phantoms — FILED as **#68**, and first in the queue

⛔ **My first statement of the cause was WRONG, and the correction matters.** I
wrote: *"single-backtick text in RST/docstrings is being bound as a Python
cross-reference and minted as a code node when it fails"*, and inferred a
`default_role` misconfiguration. `[M]` **`default_role` is not set in
`docs/conf.py` at all**, and the phantom ids are not py-domain: they are
`math:equation:R`, `math:equation:[0, 1]`, `math:equation:c > 1`,
`math:equation:SO(2)`.

`[M]` the true breakdown of the 5,276 `references`-origin phantom edges:

| sub-defect | edges | share | what it is |
|---|---:|---:|---|
| **inline math minted as an equation LABEL** | 4,383 | 83.1% | ``:math:`R` `` → `math:equation:R` |
| **unqualified py-domain role** | ~659 | 12.5% | ``:meth:`apply` `` → `py:method:apply` (33) |
| **citations typed `unresolved`** | 156 | 3.0% | `cite:p:Hebert2009` (12) |

`math:equation:c > 1` is the tell — that is not a label anyone could declare, it
is the **body** of an inline expression. A real labelled equation is
`.. math:: :label: sn-within-group-system`, cited with `:eq:`. Inline `:math:`
has no label and should mint nothing.

⚠ **Why it costs more than its node count:** `math:equation:*` is the identifier
space that `implements`, `equation_ref`, `provenance_chain` and
`verification_coverage` all traverse. ORPHEUS has **903 real equation nodes**;
these phantoms add ~4.4k more, so any "how many equations / how many verified"
denominator computed from node counts is wrong.

Citations already carry a distinct `cite:p:` prefix, so that third of the fix is
a node-**typing** change, not a resolution change — and it would make `cites` a
traversable edge (currently 0 instantiated despite the project having a
bibliography). `ontology.toml` already flags citation-as-`unresolved` as a smell.

⟹ Distinct from #55/#16, which are `calls`-origin and need type inference. This
is entirely `references`-origin and needs none.

⚠ **Consequence for §5.2:** the halo is NOT homogeneous, so "drop terminal
non-ORPHEUS nodes" is too blunt. `external` is safe to demote; `unresolved`
contains real ORPHEUS symbols that failed to bind (dropping them LOSES edges) and
prose that should never have been there.

---

## 4. The framing: three axes, not one split

The 3-way split conflates three independent questions. Separating them is the
core insight of this document.

| axis | separates | edges crossing | what it delivers | right mechanism |
|---|---|---:|---|---|
| **layer** (docs/code/test) | who authored it | 76.8% | query scoping | edge type + node property |
| **reachability** (core/halo) | terminal sinks | n/a (terminal) | −49% edges, free | query default |
| **lifetime** (derived/observed) | survives a rebuild? | ~0 | **state ownership** | separate store |

⭐⭐ **The axis that forces a physical boundary is LIFETIME, not layer.**

Test *structure* is **derived** — rebuilt from AST on every `sphinx-build`,
exactly like code. Splitting it into its own file does not change that; the file
is still overwritten. Test *state* (green / red + the error / untested) is
**observed**. Only the observed half needs to survive a rebuild, and nexus
already has that pattern — `runtime.py`'s own docstring:

> *never written into `graph.db` (which is rebuilt on every `sphinx-build`). It
> lives in a sidecar … and re-binds to the live graph at query time because node
> IDs are stable across rebuilds.*

---

## 5. The proposals

### 5.1 Retype test→code `calls` as `exercises` — ✅ **ACCEPTED 2026-08-16**

> **User:** *"We do need a type that differentiates code→code from code→test."*

**Goal.** `calls` means one thing (code→code architecture), and "what exercises
this?" is a different question with a different answer.

⚠ **Direction, because getting it backwards is a real bug and the phrasing is
easy to flip.** `[M]` the edge runs **test → code**: source is the test function,
target is the production symbol it calls (`calls` dominant pair `test->code`,
27,647 edges). The existing `tests` edge type already runs the same way
(`test->docs`, 2,745), so test→X is the established convention. The relation
reads *"this test exercises that symbol"*, and its **antisymmetry is the whole
point** — that is what `calls` reachability lacks (#60) and what a dependence
cone needs.

**Why.** `[M]` §2.3 — `calls` is the only edge type that is majority
cross-layer, and its cross half is exactly test→code (27,647 edges). Fixing it
makes the layer boundary readable from the edge type, so *inter-layer query
becomes `edge_types=`* — an API that already exists. No split needed.

**Also unlocks:** #358's subject relation (#60 — a cone needs an antisymmetric
"rests on"; `exercises` is antisymmetric where `calls` is not), and #57 would
populate a per-test version of it from coverage contexts.

**Open:** is the static `exercises` (test function → production symbol it calls)
the same edge as the runtime one (#57, from coverage dynamic contexts), or two
edges with different confidence? *Leaning: same type, different `confidence` and
`source`, exactly as `implements` distinguishes declared (1.0) from inferred (0.7).*

⭐ **The general principle the user drew from this, worth applying beyond `calls`:
look for nodes and edges doing DOUBLE DUTY and separate the duties.** Candidates
already visible: `unresolved` (four unrelated defects, §3), the node `id`
(carries identity *and* type), `contains` (35.4% cross — file→symbol and
class→method are arguably different relations).

### 5.2 Demote the halo — *proposed, needs §3's refinement*

**Goal.** A traversal query is not charged for nodes it can never traverse
through.

**Means.** Traversal queries (`impact`, `communities`, `bridges`,
`shortest_path`, transitive `callers`) exclude terminal nodes unless asked.
`[M]` −49% edges, provably lossless for reachability.

⚠ **Refined by §3:** apply to `external` (safe — genuinely foreign, 79,207
edges). Do **not** blanket-drop `unresolved`: it holds real ORPHEUS symbols that
failed to bind. Those want **repair** (§5.1's double-duty principle), not
demotion. Likely also unblocks **#58** (`bridges`/`communities` not completing at
215k edges) — that is a *compute* blow-up, and this is the input-size lever.

### 5.3 Observed state is a DATABASE — ✅ **CHOSEN 2026-08-16**; only CI is deferred

> **User:** *"If we chose DB, we can use ATTACH when we're developing (like now).
> With JSON we don't seem to gain any benefit. So even if … the CI has partial
> functionality compared to when we're developing and have both DBs, it at least
> provides more benefits when we're developing. So we don't seem to have lost
> anything by choosing DB."*

**The argument, and it is a dominance argument rather than a trade-off.** During
development both files sit on one filesystem, so `ATTACH` is available and the
DB strictly beats the blob (joins, indices, partial reads). In CI the join may
not be available — but a JSON blob would not have offered it *either*. There is
no configuration in which JSON wins, so the choice costs nothing to make now and
the CI story can be solved later with more information.

⏸ **DEFERRED, specifically:** how CI gets at the store when the two files are not
co-located (artifact fetch, export view, a merge step, or simply doing without
the join). Do not design that now.

**What is decided:**
- observed state (runtime overlays **and** test status) lives in a SQLite DB;
- `ATTACH` is the development-time join mechanism `[M]` (verified this session:
  cross-file join works, default limit 10 `SQLITE_MAX_ATTACHED`, sqlite 3.51.3);
- `RuntimeRun` becomes a **query interface over the store**, not the store
  itself (user, 2026-08-16) — today it is a Python object serialised wholesale,
  so a CI pane wanting 30 rows loads a ~900 kB blob.

**Obligations this takes on**, stated so they are not rediscovered as surprises:
a schema and its migrations (the blob had none), and `ATTACH`'s same-filesystem
constraint (the deferred CI question).

### 5.3a ✅ REMEDIED 2026-08-16 by Track 0.6 (`c51672c`) — the store is in the WRONG PLACE today, and choosing a DB is when to fix it

> ✅ **The relocation half of this section has LANDED, ahead of the DB choice.**
> `traces/` now lives at `<root>/.nexus/traces/`, outside the build output, and
> the 5 runs were copied and byte-compared across the move. The section's
> diagnosis was right and its *timing* was wrong — it said the fix should ride
> along with choosing a database; it did not need to, and waiting would have
> left the near-miss below live for another campaign. What remains open here is
> only the DB-vs-blob question (§5.3), not the location.
>
> Read the rest of this section as the argument that produced the move, not as
> a description of the tree.

`[M]` 2026-08-16 — the trace sidecars live at
`docs/_build/html/graph/traces/*.json`. That path is **inside `docs/_build/`**,
i.e. the Sphinx **output** directory, for data whose entire design rationale is
*"survives a rebuild"*. `.gitignore:7` ignores `docs/_build/`, which is correct
for derived artefacts and exactly wrong for observed ones.

⚠ **This is not hypothetical — it nearly happened in the session that found it.**
When `[graph].output` moved from `_nexus/` to `graph/`, the four existing traces
(a 251-solve coverage run and three cProfile runs) sat under the old directory
and would have been destroyed by the `rm -rf docs/_build/html/_nexus` that
cleaned it up. They survived only because they were copied out by hand first.

⟹ The observed store belongs **outside the build output**. `.nexus/` is the
natural home — it is the project-owned, non-derived directory the config and
ontology already live in — with the DB itself gitignored (`.nexus/*.db`), since
it is machine-local and potentially large. A tracked directory holding an
untracked file is fine and is the normal shape for this.

**Open (first design question when this starts):** one DB for all observed
families with a `kind` column, or one per species (runtime / test-state)? The
`ATTACH` limit of 10 argues for few; the two consumers differ (runtime → Sphinx
docs, test state → CI visualization) which argues for two.

⛔ **My earlier reasoning here was bad and the user caught it.** I wrote
*"leaning one: `RuntimeRun` is already documented as a bag of orthogonal
overlays, and test status is a fifth family in that bag."* That argues from a
docstring which is **itself the smell** — see §5.3b. A bag of optional fields is
a missing type, not a justification for adding a sixth optional field.

### 5.3b Look in the bag first — `[M]` 2026-08-16

> **User:** *"If `RuntimeRun` is a bag of orthogonal overlays, first you need to
> look into the bag and check if something is misplaced. `RuntimeRun` has a name
> strongly associated with runtime, so we don't want to mix a bunch of concepts
> in a single place."*

`RuntimeRun` (`runtime.py:299`) holds `name`, `kind`, `meta`, and five payloads.
Three measured defects, before any new family is added:

**(a) `edges` is the misplaced thing.** `calls`, `coverage` and `timeline` are
all `dict[node_id, metrics]` — **node** overlays. `edges` is
`list[(src, tgt, count)]` — an **edge** overlay. So the "orthogonal overlays"
claim is already false: the bag holds *two species*, and the difference is not
cosmetic — `merge_runs` special-cases it (`runtime.py:393-398` builds an
`edge_counts` dict then rebuilds the list) and `from_dict` must re-`tuple()` it.

**(b) `kind` is a discriminant's NAME on a field nothing discriminates by.**
`[M]` `grep -rn "run\.kind|r\.kind|\.kind =="` over the package returns **one**
hit — `server.py:804`, which merely echoes it back in the response. Zero
behavioural reads. The `KIND_*` constants are used at `server.py:784-786`, but as
keys of a table dispatching on the **CLI argument**, not on `run.kind`. And
`merge_runs` overwrites it with `KIND_MERGED` (`runtime.py:392-393`), destroying
the provenance of what was merged — so the field fails at its only remaining job.

**(c) The metric records are a missing type, already duplicated.** `[M]`
`{"ncalls": 0, "tottime": 0.0, "cumtime": 0.0}` is written literally at
`runtime.py:399` and `:481`, with a third statement of the same shape in the
field docstring at `:315`. Three copies of one record. `run.calls[nid]["ncals"]`
is a runtime `KeyError`, and `{"ncalls": "hello"}` is constructible.

**(d) The deeper one — `merge_runs` knows every family's algebra.** `[M]` 64
lines (`runtime.py:386-450`) containing the sum/max rule for `calls`, the sum
rule for `edges`, the intersection rule for `coverage`, and the *refusal* rule
for `timeline`. **Adding a fifth family means editing it** — an open-closed
violation, and precisely the edit I was about to make for test status.

⭐⭐ **Why this converges with the DB decision instead of competing with it: in a
DB, `merge_runs` DISSOLVES.** A run becomes a row (`name`, tool, command,
captured_at, ledger); an overlay becomes a table keyed by `(run, node_id)` or
`(run, src, tgt)`; and merging becomes a `GROUP BY` across runs. Every current
rule is expressible as an aggregate — `SUM(ncalls)`, `MAX(cumtime)`,
`SUM(count)`, and coverage's intersection as *"no run in the set hit this arc"*,
which reads more clearly in SQL than in the current Python. Test status ("red in
any run") is one more aggregate. **A new family then adds a table, not a branch.**

⟹ So the answer to *"one DB or two?"* is not yet the right question. The right
first question is **what the species are** — run-provenance, node-overlay,
edge-overlay — and the type that is missing is `Overlay` (a payload plus its own
merge algebra), not a container. Settle that, and the file-count question
answers itself.

⚠ **Scope check before acting:** `RuntimeRun` is 745 lines of `runtime.py`, but
the store, the CLI's five `runtime-*` verbs, six MCP tools and four `GraphQuery`
methods consume it. A restructure is a real retirement audit, not a rename —
and four ingested sidecars exist on disk in the old JSON shape.

**Goal.** Test status (green / red + the error / untested) and runtime overlays
survive a `sphinx-build`, and are queryable *with* the graph.

**The user's sharpening (2026-08-16):** now that `ATTACH DATABASE` is known to
work, the sidecar could be a **database** rather than a JSON blob — and
`RuntimeRun` might **query** the store rather than *be* the store.

`[M]` verified this session: SQLite `ATTACH DATABASE` joins across separate
files in one query; default limit 10 attached (`SQLITE_MAX_ATTACHED`);
sqlite 3.51.3.

**Argument for a DB over the JSON sidecar:**
- joins against `graph.db` in one query instead of load-all-then-filter in Python
- indices on node_id — the overlay queries are all node_id lookups
- partial reads: a CI status pane wants 30 rows, not a 900 kB blob
- concurrent writers (a test run streaming results) are a real use case for
  status in a way they never were for a profile dump

**Argument against / open questions:**
- the current sidecars are 400 kB–1.5 MB JSON and load fine; is the join actually
  the bottleneck, or is this a solution looking for one?
- `ATTACH` requires both files on one filesystem — fine locally, a constraint for
  a CI artifact fetched separately
- schema migration becomes a real obligation the JSON blob does not have

⭐ **Design constraint the user flagged, and it changes the answer: both of these
become VISUAL.** The runtime graph is destined for Sphinx docs; the test-state
graph for a CI visualization. So the store's consumer is not only an agent — it
is a renderer that wants *filtered subgraphs* with attributes attached. That
favours a queryable store over a blob, and it connects to **nexus #20**
(flow-graph directive: filtered subgraph → graphviz render for theory pages).

**Open:** one store for all observed data (runtime + test state + coverage) with
a `kind` column, or one per species? *Leaning: one — `RuntimeRun` is already
documented as "a bag of orthogonal overlays, not a tagged union", and test status
is a fifth family in that bag.*

### 5.4 `layer` as a node property + `scope=` on queries — ✅ **APPROVED 2026-08-16**

Design detail in §6.1. Approved as a direction, not yet as an implementation —
the open sub-questions are still open:

- **What are the layer values?** `{code, docs, test, external}` is what the
  measurement used. Is `external` a layer, or is it the *reachability* axis
  (§5.2) wearing a layer label? They coincide today; they are different concepts.
- **Where is it computed?** Derived at build time from `is_test` + node type +
  source directory. It must stay derived — a hand-maintained layer is a second
  source of truth for something the extractor already knows.
- **`scope=` on which queries?** Everything that traverses. Note `--db`'s history
  (#66): declaring a parameter per-subparser 42× is how the last one drifted.
- ⚠ **Does `scope=` compose with `edge_types=` or duplicate it?** §5.1 makes
  layers readable from the edge type; if that lands, `scope="test"` and
  `edge_types=["exercises"]` may answer the same question two ways. They are
  claimed to be complementary (node-kind vs relation-kind) — verify that on a
  real query before shipping both, or this becomes a double-duty seam of its own.

### 5.5 Payload ergonomics — **FILED as #67**

`[M]` −78% on `solve_sn` with a slim representation, tests included. But per §1
the goal is *relevance*, not size — see §6.2 for what "ergonomic" must mean here.

---

## 6. Detail requested by the user

### 6.1 What "`layer` as a node property" actually means

**The claim:** everything the 3-way file split was for — *query control* — is
available from **one column and one parameter**, with no split, no join, no
rebuild transaction.

**The mechanism.** Every node gains `layer ∈ {code, docs, test, external}`,
computed at build time from what nexus already knows (`is_test`, node type,
source directory). It is *derived*, so it costs nothing to maintain and cannot
drift. Then every query takes `scope=`:

```
callers("solve_sn")                        # today: everything, tests included
callers("solve_sn", scope="code")          # architecture: who really depends on this
callers("solve_sn", scope="test")          # coverage: what exercises this
impact("solve_sn", scope=["code","docs"])  # blast radius excluding the suite
```

**Why a property beats a partition, concretely:**

- **A node belongs to one layer; an edge does not.** 76.8% of edges cross. A
  property attaches the label where it is unambiguous (the node) instead of
  forcing every edge to choose a home.
- **Scoping becomes a filter on a traversal, not a choice of which file to open.**
  `scope=["code","docs"]` is one predicate; with three files it is a decision
  about which two to `ATTACH` before you know what the answer touches.
- **It composes with §5.2.** `external` is just another layer value, so
  "demote the halo" and "scope to code" are the same mechanism, not two.
- **It is reversible.** A wrong layer assignment is a rebuild. A wrong file split
  is a migration.

**What it does NOT give you:** state ownership (§4 — that is the lifetime axis),
and it does not by itself reduce payload (§5.5).

⚠ **Relationship to §5.1:** these are complementary, not alternatives. The edge
type says *what kind of relation this is*; the node layer says *what kind of
thing this is*. You want both — `edge_types=["exercises"]` and `scope="test"`
answer different questions, and today neither is expressible.

### 6.2 What "ergonomic" has to mean (user, 2026-08-16)

> *The objective is not just slimming, but to make it highly ergonomic for you
> and other agents to use it. It should deliver highly relevant context in a way
> that is helpful for your actions and helpful in a way that grep is not.*

Consequences for the design, to be held against every proposal here:

1. **Relevance beats completeness.** A capped list of 25 arbitrary neighbours is
   worse than 10 chosen ones. Today's `limit_per_type=25` truncates by insertion
   order — it should rank. *Open: rank by what? Degree is available and wrong
   (it favours hubs). Candidates: same-module first, then `type_uses`/`inherits`
   over `calls`, then runtime-confirmed (§5.3) over static.*
2. **Answer the question, not the schema.** "Who calls this" should not require
   knowing that test callers arrive as `calls` while doc mentions arrive as
   `references`.
3. **The differentiator vs grep is RELATIONSHIPS AND ABSENCE.** Grep finds text.
   The graph should answer *"what would break"*, *"what is untested"*, *"what
   implements this equation"* — and, per **L56**, must distinguish "nothing"
   from "I could not tell".
4. **Cheap by default, complete on request.** The default should be affordable
   enough that an agent reaches for it without being told to.

---

## 6bis. Nexus's own architecture — adversarial review, 2026-08-16

Full report: `scratch/nexus_architecture_review.md` (untracked working space —
the verified essentials are distilled here so they survive it). A second copy
of the survey is at
`.claude/agent-memory/elegance-enforcer/nexus_architecture_survey.md`.

⛔⛔ **DO NOT re-derive the ordering from that report.** It was written under
the constraint described in the method note below, so it still carries
*"sound layering, loose vocabulary, do NOT restructure"* as its verdict —
the framing that was **corrected mid-flight and superseded**. A cold reader
landing in it will find an authoritative-sounding conclusion that this plan
contradicts, and nothing inside the report says so. **§6ter's order is
MINE, not the review's**, and its three objections to the review's order are
measured. Read §6ter; use the report only for the raw findings it enumerates.

**Method note, because it changed the answer.** The first pass was briefed with
*"say what is well-factored"* in the same breath as the findings, and returned
*"sound layering, loose vocabulary, do NOT restructure"*. That verdict was then
(by me) fed back as a **constraint** on the adversarial pass — a Phase-2
conclusion bounding Phase 1. Both errors corrected mid-flight; with the
constraint lifted, the answer changed. See
[[feedback-adversarial-phase-before-balance]].

### Verified independently (I re-measured these; they are not relayed on trust)

**R1 — `GraphQuery` has TWO states, and one is an optional keyword argument.**
`[M]` by AST over `query.py`: **56 methods, 10 take a working-tree argument**
(`node_at`, `detect_changes`, `rename`, `verification_audit`, `staleness`,
`retest`, `session_briefing`, + 3 private). So it is not "a namespace over
`self._g`" — it is a procedural module over *two* states, one held and one
threaded. That is the user's *"data and functions that share a context,
un-consolidated"* case exactly.

⭐ **The causal payoff:** a staleness check needs *both* states. With the tree as
an optional kwarg, applying it is something every tool author must remember — so
it landed at **1 of 40** tools. With graph+tree as one object it is paid at
construction and becomes unforgettable. `Workspace`, `GitProvenance` and
`stamp_provenance` already exist, unassembled.

**R2 — the node id is a Sphinx *artefact* key, not a symbol identity.** `[M]` on
ORPHEUS's graph: ids carrying a Sphinx **role** spelling rather than the
canonical **objtype** — `py:func` 206, `std:doc` 94, `py:property` 68, `py:meth`
23, `py:attr` 22, `py:obj` 8. `[M]` **527 names are carried by more than one node
id**; `ValueError` exists under four (`class` / `exc` / `exception` /
`function`).

`[M]` **Cause is exactly one-sided, and it is two lines.** The role→objtype map
lives at `ast_analyzer.py:920` (`{"func": "function", "meth": "method", …}`) and
is applied on the **AST** minting path at `:947`. The **Sphinx xref** path at
`extractors.py:402-404` puts the raw `reftype` straight into the id:
`_node_id(refdomain or "std", reftype or "any", reftarget)`. One normalisation,
two doors, applied at one.

⟹ This is the same root as **§3.1's** type-prefix fragmentation and **§5.1's**
double-duty principle: the id carries identity AND type, and the type has two
spellings.

**R3 — ids are built from raw text with no normalisation at all.** `[M]` 13 node
ids contain a **literal newline**, e.g.
`math:equation:(P_1 P_2) x =\nP_1(P_2 x)`. All are `math:equation:` — i.e. the
same defect as **#68**, wearing a different symptom. An id is a key; keys are
normalised at the producer (Pattern 7) or every consumer pays.

### The revised verdict

*File* rearrangement stays rejected — `[M]` the layering has no inversions, and a
seven-file split of `query.py` adds seven concepts and removes none. **Object
extraction is not the same question**, and it survives: `NodeId`, `ProjectView`,
`PositionIndex`, `Evidence`, `Diagnostic` — five objects latent in
`GraphQuery` / `workspace.py` / `server.py` / `runtime.py`, all already
half-built, with the concept count going **down**.

⚠ **Two live correctness defects, both invisible to the suite** (`[M]` 739
passed / 1 skipped at review time):
- `ast_analyzer._unparse_attribute` fabricates `calls` edges —
  `get_thing().method()` mints an edge to a bare `method`. `[M]` 8 of 13
  incoming `calls` on one function in nexus's own graph are fabricated. Fails in
  the false-**alive** direction, so it inflates `impact`/`retest` and *rescues*
  dead functions from `dead_functions`. One-line fix.
- Three implementations of "(file, line) → node" disagreeing on 3 of 4 probed
  positions (`query.node_at`, `runtime.build_node_index`+`resolve_node`,
  `brief._in_file_node_ids`).

⛔ **And one thing I shipped this week that must not stay as-is:** `ontology.py`
has `[M]` **zero production consumers** (only its own test imports it), while
`merge.py:248` still ships the hardcoded `not tgt_attrs.get("in_test_file")` rule
its docstring claims was dissolved into a declared constraint. I took the
concept count from 2 to **3** and wrote a present-tense-false docstring saying
otherwise. Either wire `check_edge` into `_finalize_graph` **and** retire
`merge.py:248` in one commit, or hold the module.

⚠ Likewise `NodeBinder`'s docstring (mine, this week) says *"**Every** backend
performs this join"* — measured over the three **runtime** backends;
`query.py:757 node_at` is a fourth. `plan-authoring` §2's quantifier defect, in a
docstring, by its own scribe.

## 6ter. ▶ THE RESHAPE IS THE FIRST TASK — ordering, critiqued

**User ruling 2026-08-16:** the nexus reshape goes **before** #68 and before the
`layer` property — *"better to do this first so we can build on top of better
foundations."* Everything in §5 now sits behind it.

The review proposed: `NodeId` → `ProjectView` → `PositionIndex` → `Evidence` →
`Diagnostic` → wire `Ontology`. I do **not** accept that order. Three objections,
each measured.

**(a) It bundles a 2-line correctness fix into a 68-site refactor.** R2's
*defect* is one-sided normalisation — the `py_type_map` at
`ast_analyzer.py:920` is applied on the AST path (`:947`) and not on the Sphinx
xref path (`extractors.py:402`). Applying it on both is **two lines** and fixes
the fragmentation today. Minting a `NodeId` *type* is the reshape that prevents
recurrence. They are separable, and the fix must not wait for the refactor.

**(b) `NodeId` first maximises time on unstable ground.** `[M]` 68 id-construction
sites across 10 modules — the largest blast radius of the five — plus heavy
test-literal pinning (`test_dead_references.py` 111, `test_ast_analyzer.py` 67,
… and 17 in ORPHEUS itself). It also changes the graph's *identity space*, which
every stored sidecar and every committed expectation rests on. Doing it first
means the other four are built on a foundation that is still moving.
⟹ when it does happen, scope it **producers-first** (the ~45 minting sites in
`extractors.py` + `ast_analyzer.py`), because normalising at the producer is what
makes the defect unspellable; consumer adoption can follow incrementally.

**(c) `Ontology` last is backwards — it is a foundation, not a capstone.** It is
the *admission authority*: it decides which edges may exist. Landing it after the
extractions means each is built without it and may then be rejected by it. It is
also the natural home of the vocabulary the review says is missing (`NodeId`'s
valid type set, `implements`'s domain), so `NodeId` and `Ontology` are the **same
layer**, not opposite ends. ✅ Partly done already — see below.

### The order I would run instead

**Track 0 — correctness now, prejudging no reshape.** Each is small, local, and
independently valuable:

| # | item | size | done when |
|---|---|---|---|
| 0.1 | ✅ **Ontology wired** — `_infer_implements` consults it | LANDED `37bad88` | `merge.py` holds no copy of the rule |
| 0.2 | ✅ **A call edge is never fabricated** — the truncating twin RETIRED | LANDED `3e137ff` | `get_thing().method()` mints no edge |
| 0.3 | ✅ **A Python id is spelled with the objtype on both producers** | LANDED `2ddf61b` | `[M]` short-prefix ids **316 → 0**, split symbols **263 → 0** |
| 0.4 | ✅ **The equation namespace holds declared labels, not LaTeX** (**#68** CLOSED) | LANDED `090a793` | `[M]` `math:equation:*` **1860 → 903**, all declared; newline ids **13 → 0** |
| 0.5 | ✅ **A project may WIDEN a base edge, never narrow it** (**#69** CLOSED) | LANDED `7882dba` | `[extend.edge.implements] range = [...]` widens; narrowing raises; monotonicity asserted over every node type |
| 0.6 | ✅ **The graph lives with the project, not with the build output** | LANDED `c51672c` | `nexus config db` → `<root>/.nexus/graph.db`; `[graph].db` is a REFUSED key; a clean build no longer destroys `traces/` |

**0.6 was not in this plan — it came from the user asking, at the Track 1.1
design checkpoint, "why is the graph in `docs/`? Shouldn't it be in
`.nexus/`?"** The answer was "because nexus is a Sphinx extension and
`output` resolves against `app.outdir`", and checking it turned up a
data-loss bug: `RuntimeStore` sat at `db_path.parent / "traces"`, i.e. inside
the directory `rm -rf docs/_build` empties, while `runtime.py`'s own docstring
says the sidecar exists precisely to survive the rebuild. Three lifetimes in
one directory — derived / served / **durable** — and a directory's lifetime is
its most-derived member's.

⭐ **`project.py`'s module docstring had already said the graph belongs there**
(*"the directory also holds `ontology.toml` and (once built) the graph and its
overlays"*). The convention was documented and never implemented, so this was
honouring a stated design, not inventing one. ⟹ when a layout question comes
up, read the module docstring of the thing that owns the layout *before*
designing — the intent may already be written down, and then the work is to
make the code match rather than to choose.

⛔⛔ **Two pre-existing defects fell out, and one of them was mine.**

| | what | how it was caught |
|---|---|---|
| **pre-existing** | `find_project_root` walked past every boundary into `$HOME`, because `~/.nexus/` exists on any machine that has run the MCP server (it holds `usage.jsonl`). `[M]` `find_project_root("<nexus repo>/tests/roots")` → `/Users/rodrigo`. Every unconfigured tree under `$HOME` had been READING the wrong project's config, silently. | only when the store move turned that read into a **write** and `graph.db` appeared in `~`. Promoted to [[lessons-L57]] |
| **pre-existing** | the `output` default was spelled **twice** — the `pick()` fallback and the `nexus_output` config value. | changing one silently relocated every unanchored project's graph; 26 tests reddened. Now one `DEFAULT_OUTPUT` |
| **mine, caught by the tests** | my first design anchored the store to the root `find_project_root` discovers. For a docs *fixture* that walks out of the fixture entirely, so a test build would write into the repository source tree — or, given the bug above, into `~`. | the suite. Corrected to: a project is **ANCHORED** only when it declares itself with a `.nexus/`; unanchored, `root` is a guess and the store stays with the build |

⛔⛔ **BOTH LANDED ITEMS REFUTED THEIR OWN SIZE ESTIMATE, THE SAME WAY.** 0.2
was priced at "1 line", 0.3 at "2 lines". Neither was a line count — each was a
**duplication**, and the reported defect was the copy that happened to be wrong:

| item | the reported fix | what it actually was |
|---|---|---|
| 0.2 | make `_unparse_attribute` return `None` | `_dotted_name`, **360 lines above in the same file**, already did — retire the twin (−28 lines) |
| 0.3 | apply `py_type_map` on the xref path | the map was a dict **local to one producer** — hoist to `_mappings.REFTYPE_OBJTYPE_MAP`, both read it |

⟹ **Do not inherit 0.4's "small".** It is filed as *"an id built from raw text
with no normalisation"* — which is already the language of a producer with more
than one spelling. Enumerate the id-minting sites FIRST (0.4's own note says
~45 across `extractors.py` + `ast_analyzer.py`), and expect the repair to be
"one home" rather than "one edit". Promoted to [[lessons-L57]].

⛔ **0.4 CONFIRMED IT, and refuted this plan's own description of the defect.**
It is **not** a normalisation defect — it is a *semantic* one, plus a
duplication. `:eq:` REFERENCES a label; `:math:` TYPESETS an expression and
references nothing. The scanner routes both into the equation namespace on
purpose (forgiving a common authoring slip); its guard was a **blocklist**
(reject `\`, `{`, `}`) sitting three lines above a Python branch asking the
opposite, stronger question (`_is_dotted_identifier`). Inline math clears a
blocklist trivially. The 13 newline ids were **12 + 1**, not 13 of a kind: 12
inline math, and **1** the doctree walker missing the wrap-normalisation the
docstring scanner has always had — the *same* one-of-two-producers split as
0.3, now the fourth instance in this file family.

⚠⚠ **And the repair introduced a regression that only the arithmetic caught.**
Sharing `_normalize_wrapped_target` with a second producer exposed that it was
already wrong: it collapsed ALL whitespace whenever the result matched the
dotted shape, and a sentence of letters plus a full stop matches once its
spaces are gone. `_classify_unresolved` refuses a non-identifier as napoleon
noise — so prose was already being dropped correctly, and collapsing first
*disguised prose as an identifier* and walked it through that gate. `[M]` 48
junk classes minted (`py:class:allkeyvariables.`). Found by noticing 957 nodes
dropped but net −910, i.e. **48 appeared**. ⟹ **when a pass is meant to be
purely subtractive, check that its node delta is purely subtractive** — the
totals close or they do not, and it costs one query.

⭐ **0.4 absorbs #68.** The prose leak (inline math minted as an equation
*label*) and the 13 newline-bearing ids are the **same defect**: an id built from
raw text with no normalisation. They were filed separately because I found them
separately; they are one fix at one producer.

> ⛔ **REFUTED 2026-08-16 by doing it** — the conclusion (fold #68 into 0.4)
> was right, the reasoning was wrong in both of its clauses. They are **not**
> the same defect: 12 of the 13 newline ids were inline math (a *semantic*
> error — `:math:` is not a reference), and **1** was the doctree walker
> missing a wrap-normalisation (a *normalisation* error). Nor is it **one
> producer**: it is two, and the fix is two, in two different modules. Kept
> because the shape recurs — "these were filed separately, so they are
> probably one thing" is a hypothesis, and here it was 12/13 right and
> 0/2 right about the mechanism.

**Track 1 — the reshape, reordered by leverage and stability:**

1. **`ProjectView`** (R1) — graph + working tree as one object. Highest
   *correctness* leverage of the five: staleness needs both states, and with the
   tree as an optional kwarg it is something every tool author must remember, so
   it landed at **1 of 40**. Self-contained; touches no ids.

   ✅ **DONE 2026-08-16** as Track 1.1 — `35d5d31` (bind `Workspace`) +
   `8de24ec` (staleness at the tool boundary). The diagnosis in this row is
   exactly right and is why the step paid: "something every tool author must
   remember" is the whole defect, and the repair was to make it something no
   author has to. Read the ✅ block at the top of this plan, not this section.

   > ⛔ **MECHANISM REFUTED 2026-08-16** — the *outcome* stands, the named type
   > must not be built. `workspace.Workspace` already IS "one checkout paired
   > with its graph database" (its own docstring), so `ProjectView` would be a
   > fourth spelling and would leave `project_root` travelling per-call. Bind
   > `Workspace` into `GraphQuery` instead. Full reasoning, the verified
   > precedent table, and the live defect that motivates it are in the
   > **RESUMES AT** block at the top of this plan; do not design from this
   > section alone.

   ✅ **RE-DERIVED AT HEAD 2026-08-16, at the compaction point** — both halves
   reproduce, and the units are now written out because the first re-derivation
   *appeared* to refute the number:

   | claim | at HEAD | how to re-run |
   |---|---|---|
   | `ProjectView` does not exist | **0 hits** in `sphinxcontrib/` + `tests/` | `grep -rn "ProjectView"` |
   | 10 of 56 `GraphQuery` methods take a working-tree arg | **10 of 56** — 7 public (`detect_changes`, `node_at`, `rename`, `retest`, `session_briefing`, `staleness`, `verification_audit`) + 3 private (`_apply_renames`, `_git_changed_files`, `_git_file_timestamps`) | AST-walk `query.py`, arg name `project_root` |
   | staleness is applied by 1 of 40 MCP tools | **1 of 40** — `_position_staleness_warning` is defined at `server.py:217` and called at `:428` only, inside `node_at`; 40 functions carry `@nexus_tool` | `grep -c "@nexus_tool"`, `grep -n "_position_staleness_warning"` |

   ⚠ **The unit was ambiguous and it cost a reconciliation cycle.** "56 methods
   / 10" does not say whether *methods* means public or all. Filtering to
   **public** gives **7**, which reads as a refutation of the plan; the 10 is
   right and counts private helpers too. A second trap sits beside it: an
   arg-name heuristic wide enough to catch `start` also catches
   `_dominant_call_chain`, whose `start` is a **graph node id**, not a path.
   ⟹ state the FILTER beside the fraction, not just the two numbers.

   > ⛔ *Done when:* **"staleness is applied by construction, and the count is
   > 40 of 40"** — REFUTED 2026-08-16, before design. `[M]` only **3 of 40**
   > tools take a position as INPUT, so 37 have nothing for a position-staleness
   > warning to apply to; the target silently changed the predicate its own
   > `1 of 40` measurement had counted. The real exposure is the **51**
   > `_node_result` producers that EMIT `(file_path, lineno)`. Replacement
   > done-when (user-ruled) is in the **RESUMES AT** block.
2. **`PositionIndex`** (F2) — ✅ **LANDED 2026-08-16** @ `7db466d` + `c497ddb`.
   **A prerequisite for the test-state work**, which needs exactly this binder;
   adding it as a fourth consumer of a three-way disagreement would have made it
   four-way. ⛔ *Done when: the 4 probed positions agree* — **REFUTED before
   design**: two of the three disagreements are the trace join's deliberate
   policy (no module to bind to, classes shadow their methods), and forcing
   agreement would have broken it. See the Track 1.2 block at the top for what
   replaced it and what landed.
3. **`NodeId`** — producers first, per (b). ▶ **NEXT**, and **re-chartered
   2026-08-16 after reconciling against the tree** — see the block below; the
   original evidence was measured before Track 0.3 repaired half of it.
   *Done when:* an id cannot be constructed except through the type, **and the
   143 same-domain duplications are 0**.

#### `NodeId` re-chartered — `[M]` 2026-08-16, ORPHEUS graph 23013 / 206868

Track 0.3 (`2ddf61b`) fixed the **spelling** disagreement in the `py` domain
(role `func` vs objtype `function`), and said so in its own commit: *"applied
only when `refdomain == 'py'`"*. `[M]` role-spelled `py:` ids are now **0**.
What survives is not a spelling problem at all:

| shape | names | both halves carry edges | what is split |
|---|---|---|---|
| `py:method` + `py:property` | **68** | **68 of 68** | AST node holds `calls`/`type_uses`; Sphinx node holds `documents` |
| `py:class` + `py:function` | 43 | 43 | builtins — `int` as annotation (3271 `type_uses`) vs `int()` called (809 `calls`) |
| `py:function` + `py:method` | **16** | **16 of 16** | same split as the properties, e.g. `AngularBoundaryFlux.zeros_on`: 252 `calls` on one, `documents` on the other |
| `py:function` + `py:module` | 16 | 16 | a test module documented as a function |
| `doc` + `std:doc` | 94 | — | the `std:doc:` half is **inert** (`out={}`), i.e. 94 orphan duplicates |
| `math:equation` + `std:label` | 10 | — | |

⚠ **CORRECTION 2026-08-16 — "143" was the top FOUR shapes, reported as the
total.** `[M]` re-derived independently: **174** `py` names carry more than one
objtype; the four shapes tabulated above sum to **144** (`function + method` is
17, not 16). The tail is 30 more names across 11 further shapes. ⭐ This is
`plan-authoring` §2's quantifier defect — *mine*, in a design document, hours
after writing [[lessons-L59]] about exactly this class. The tell was available
and I did not look: I summed a `most_common(4)` and wrote the sum as though it
were `len()`.

⭐⭐ **The consequence, and it is the one that matters for this project: for
85 symbols the CALL graph and the DOCUMENTATION graph hang off different
nodes.** `orpheus…Mixture.ng` has `calls`/`type_uses` on `py:method:…ng` and
`documents` on `py:property:…ng`. `[M]` aggregated over all 68 property pairs:
the method-spelled half holds `references 109 / type_uses 84 / calls 55`, the
property-spelled half holds `documents 55`. **Even the `contains` parents
differ** — the class contains the method node, the doc page contains the
property node, so a class's member list and its own doc page's member list are
disjoint views of one symbol.

⛔ **This means minting a `NodeId` type is NOT sufficient, and the original
charter mis-states the fix.** The producers do not disagree about how to
*spell* an objtype — 0.3 settled that. They disagree about **what the objtype
IS**: the AST sees a `method`, the Sphinx domain sees a `property`. A type that
validates spelling would happily mint both. What is needed is a **single
authority on a symbol's objtype**, and that is the ontology — which is
objection (c)'s point (*"`NodeId` and `Ontology` are the same layer"*) arriving
with a measurement behind it.

⟹ *Proposed means (2026-08-16, NOT yet verified):* the ontology declares the
objtype lattice and which producer wins when two disagree (`property` refines
`method`; `module` is never a `function`); `NodeId` is the type that cannot be
constructed without asking it. Producers-first, ~45 minting sites.

⚠ Related and already filed: **nexus #71** — ⛔ **NOT the same defect**, and I
guessed that it was. `[M]` **1 of 185** lost-decorator names is a split name,
against a 2 % base rate among the 3345 that kept theirs. Independent; neither
fix helps the other, and #71 is **not** a prerequisite. (The issue's own stated
mechanism was also too strong and has been corrected in a comment: `source=both`
is necessary, not sufficient — 185 of 2812, 6.6 %.)

#### The design — `[M]` 2026-08-16

**The 143 are FOUR defects, not one.** Decomposed by what the halves actually
*are* (the `type` attribute), rather than by id prefix:

| n | the two halves | family |
|---|---|---|
| 94 | `file` + `file` | **A — two id NAMESPACES for one concept** (`doc:X` vs `std:doc:X`) |
| 73 | `attribute` + `method` | **B — two producers TYPE one symbol differently** (`@property`) |
| 44 | `external` + `external` | D — one builtin used two ways (`int` annotated vs `int()` called) |
| 17 + 7 + 6 | placeholder + real | C — a phantom that never folded onto its definition |
| 16 | `unresolved` ×2 | D |
| 10 | `equation` + `section` | out of scope — a labelled equation is also a `:ref:` target |
| 2 | `class` + `exception` | B |

⟹ **A and B are this step. C and D are not** (C is `_canonicalize_phantoms`'s
job and D may be legitimate), and saying so is the point — the original charter
would have had me fix all 143 with one instrument.

**The invariant that makes A checkable, and it is one number.** A node's `type`
is EITHER the type segment of its own id, OR a declared placeholder state.
`unresolved` is a *resolution state*, not a kind — `py:function:foo` with
`type=unresolved` is correct by design, and the ontology already marks those
`placeholder = true`. `[M]` at HEAD:

| | count |
|---|---|
| id type-segment **==** node type | 17456 |
| differs, but type is a **placeholder** | 4595 — correct by design |
| **VIOLATIONS** | **962**, in exactly 5 shapes |

| n | shape |
|---|---|
| 680 | id says `label`, type is `section` |
| 94 | id has **no type segment** (`doc:api/data`), type is `file` |
| 94 | id says `doc`, type is `file` |
| 68 | id says `property`, type is `attribute` |
| 26 | id has no type segment (`citation:<key>`), type is `unresolved` |

⭐ `property`, `label` and `doc` are **not declared node types** — the ontology
declares `section`, `file`, `attribute`, `method`. So the producers are minting
ids in a vocabulary the ontology does not have, and nothing refuses them. That
is the whole defect in one sentence, and it is Pattern 4: the illegal state is
representable because an id is a **string built at ~45 sites**.

**Two parts, separable, in this order** (the plan's own objection-(a) precedent,
vindicated by 0.3 landing the fix ahead of the refactor):

**1 — `NodeId`, the type.** A frozen `(domain, type, name)` whose constructor
refuses a `type` the ontology does not declare, plus `parse()`/`__str__` so the
string form round-trips. Producers translate their own vocabulary through the
map Track 0.3 already established as its one home (`REFTYPE_OBJTYPE_MAP` in
`_mappings.py`) *before* minting; an unmapped producer word then fails **loudly
at mint time** instead of silently opening a new id namespace. Fixes family A
by construction — `doc:X` and `std:doc:X` cannot both exist when one type builds
both.
*Done when:* the invariant gate reads **0** violations, and `grep` finds no
`f"py:{...}:"`-style construction outside `NodeId`.

**2 — the arbiter.** Fixes family B, which part 1 does **not**: once `property`
maps to a declared type, both halves become individually *valid* and still two
nodes. Something must decide which type a symbol IS when producers disagree.
The ruling, from `coding-standards`' type-minting criterion (mint a type iff
≥2 non-isomorphic realizations AND a non-identity morphism): a Python property
has one realization and no morphism, so **it is a `method` with a flag, not a
type**. It has a body, a file position and callers; `attribute` nodes have none
of those. So `property` → `method`, `is_property` in metadata, and the two
halves merge — uniting 68 symbols' `calls`/`type_uses` with their `documents`.
*Done when:* `[M]` the 75 family-B duplications are 0, and
`orpheus…Mixture.ng` has both its call edges and its doc edges on one node.

**The invariant is not a new idea — it is a PUBLISHED PROMISE the tree breaks.**
`server.py:1606` advertises to every MCP client:
`"node_id_format": "<domain>:<type>:<qualified_name>"`. That says the middle
segment IS the type. 962 nodes say otherwise. So part 1 is not "impose a new
rule", it is "make the tree honour the contract it already ships" — which is
also why no consumer can be blamed for trusting it.

**`[M]` the id↔type correspondence is currently encoded in NINE places**, which
is the argument for the type in one line: three producer-private helpers
(`extractors._node_id` `:32`, `extractors._doc_node_id` `:36`,
`PyASTAnalyzer._node_id` `:857` — none imported by the other producer), **27**
inline construction sites across 8 modules, **38** helper call sites, plus
lookup tables that must stay in step (`_ID_PREFIX_TO_TYPE`
`ast_analyzer.py:1621`, the candidate-objtype loop `_mappings.py:361`, and a
consumer at `directives.py:133` that *guesses* by trying
`py:function:`/`py:method:`/`py:class:` in turn — which is why it can never
find a `@property`).

⭐ **The one site that IS the whole family**, `extractors.py:71-76`: the type is
the MAPPED value and the id is built from the RAW one.
```python
node_type = DOMAIN_TYPE_MAP.get((domain_name, obj_type), obj_type)
node_id   = _node_id(domain_name, obj_type, name)   # <- raw, unmapped
```
`py:property:`, `std:doc:` and `cite:p:` all come from this line pair.

⚠ **Blast radius — SMALLER than I warned, and measured.** I wrote that 962
nodes changing id would cost committed test literals. `[M]` `tests/` carries
1171 id literals over 423 distinct strings, but the two spellings this step
actually retires appear **0 times**: `py:property:` is 0 in `tests/`,
`sphinxcontrib/` and `docs/`; `std:doc:` is 0 in `tests/` and `sphinxcontrib/`.
Retiring them pins against nothing. The costly ones are the spellings we do NOT
have to touch — `py:method:` (69 literals), `py:attribute:` (24), `doc:` (~68).

⟹ **that reopens a design choice I had closed.** Family A can be fixed two ways:
- **strict** — every id gets a type segment, so `doc:api/data` becomes
  `std:file:api/data`. Honours the published grammar exactly; costs ~68 test
  literals and touches the 120-node `doc:` namespace that currently carries all
  the traffic.
- **minimal** — retire the inert `std:doc:` twin (94 nodes, `out={}`) onto
  `doc:`, and declare `doc:<name>` a documented two-segment exception in the
  grammar.
`[M]` the minimal option costs 0 test literals and removes 94 duplicate nodes;
the strict option additionally fixes 680 `std:label:`→`section` and the 26
`citation:` ids, i.e. it is the only one that reaches **0 violations**.
**Proposed (2026-08-16, NOT ruled on): strict, because a published grammar with
a documented exception is a grammar nothing can be validated against** — and an
exception is exactly how `py:property:` got in.

⚠ Also found, informational: **72 `cite:p:<key>` nodes** exist because
`extractors.py:363` special-cases `refdomain == "citation"` while
`sphinxcontrib-bibtex` emits `refdomain="cite", reftype="p"` — so the branch
never fires and the generic path mints a fourth namespace. Same defect class,
own fix.

#### The vocabulary itself — evaluated `[M]` 2026-08-16 (user ruling: strict, and audit the types)

⭐⭐ **The ontology already encodes the answer, and nothing checks it.** Every
`[node.*]` declares an `origin` — `sphinx`, `ast`, or `derived`. **Two are
FALSE**, and both falsehoods are the defect this step is chasing:

| type | declared `origin` | who can actually assign it | `[M]` result |
|---|---|---|---|
| `exception` | **`ast`** | **Sphinx only** — `("py","exception")` in `DOMAIN_TYPE_MAP`. `NodeType.EXCEPTION` appears in `ast_analyzer.py` only inside `_CANONICAL_TYPES` and `_ID_PREFIX_TO_TYPE`; **no assignment site** | **2** nodes — and 24 more exception classes in `orpheus/` are typed `class`. Both of the 2 are `source=None` + a `docname`, i.e. Sphinx-only, and both are halves of a `class`+`exception` split |
| `type` | **`ast`** | Sphinx only — `("py","type")`, plus C/C++ typedefs | **0** nodes |
| `property` | **not declared at all** | Sphinx only | 68 nodes, all split |

⟹ **the missing rule, and it is checkable statically: the producer that OWNS
the fact must be able to assign the type.** Whether a class is an exception is
a fact about its *bases*; whether a method is a property is a fact about its
*decorators* — and the AST analyzer already records `decorators`
(`ast_analyzer.py:1290`) and throws the conclusion away. The ontology says
`origin = "ast"` for both; the AST delivers neither.

**Rulings, each by `coding-standards`' type-minting criterion** (mint a type
iff ≥2 non-isomorphic realizations AND a non-identity morphism is applied):

- **`exception` → RETIRE as a type; it is a property of `class`.** An exception
  *is* a class: one realization, no morphism. It is currently assigned by
  whether autodoc happened to document it, which is not a fact about the
  object. `[M]` 2 of 26.
- **`property` → never mint it; it is a property of `method`.** Already ruled.
- **`type` → KEEP, and fix its `origin` to `sphinx`.** A type alias is used in
  annotations and never called or read as a value — that *is* a different
  morphism, and the C/C++ typedef entries make it multi-language. Its only
  defect is a mislabelled origin, not existence. ⚠ This is where the user's
  point bites: a declared-but-unused type is not waste, it tells an author
  `.. py:type::` will be understood.
- **`attribute` vs `data`** — NOT duplicates: class/instance scope vs module
  scope, both populated (2387 / 1676).
- **`external` vs `unresolved`** — NOT duplicates: "defined elsewhere" vs
  "defined nowhere". But `unresolved` is CONFLATED (below).

**Two types are genuinely MISSING, and one of them is what this project runs
on:**

1. ⭐⭐ **`error` / defect.** `[M]` 224 nodes carry a `catches` marker naming
   **78 distinct `ERR-NNN` entries**, and **0 of those 78 exist as a node**.
   Compare the sibling marker: `verifies` (856 nodes) resolves to **2639
   `tests` edges** into real `equation` nodes. So the V&V triangle is
   **equation ✅ / test ✅ / defect ❌** — the error catalog is a first-class
   concept in this project (`vv-principles/error_catalog.md`) and the graph
   cannot express it. `catches` should be an EDGE into an `error` node, not a
   string in an attribute. This is the gap that makes "which tests catch
   ERR-051?" unanswerable by the graph today.
2. **`citation`.** 72 `cite:p:<key>` nodes typed `unresolved`. The ontology's
   own comment under `[node.unresolved]` already says this is wrong — *"a cited
   work is a real bibliographic entity, not a dangling reference, so it inflates
   the `unresolved` count and falls inside `dead_references`' phantom set…
   Modelling it as its own node type is a candidate change, not a described
   one."* The comment has been waiting for this step.

**What "test needs types properly" resolves to** — `[M]` the test layer today
is `is_test` (7305) / `in_test_file` (9530) / `vv_level` (1490) / `verifies`
(856) / `catches` (224), all **node attributes**. A test function IS a function
(same realization, has callers and a body), so `test` fails the minting
criterion and stays a property. What is missing is not a *test* type — it is
the `error` node its `catches` marker points at, plus the `exercises` edge
already accepted in §5.1. ⟹ the test layer needs **one node type and one edge
retype**, not a parallel vocabulary.

**The enforcement the user asked for — three gates, cheap:**
1. **mint-time** — `NodeId` refuses an undeclared type (part 1). Turns a silent
   new namespace into an error at the producer.
2. **build-time warning** — an id whose type segment is not a declared type
   warns during `sphinx-build`. This is the half that serves doc authoring: an
   author writing a directive nexus does not know gets told so.
3. **origin gate** — for every `[node.X]`, some producer of its declared
   `origin` must be able to emit it. `[M]` this fails today on `exception` and
   `type`, which is how a type gets declared, documented, and never produced.
4. **`Evidence` / `Diagnostic`** — the diagnostics family, once the vocabulary
   exists. `[M]` 285 lines of adapter across two surfaces, with each diagnostic's
   contract written twice in different words (`native_place`'s ranking rule is in
   `query.py` and **absent** from the MCP docstring an agent actually reads).
   *Done when:* a new diagnostic is one object, registered once.

⚠⚠ **GATE ON TRACK 1 — evidence status per object.** These five come from ONE
review. Two I re-measured myself; three I have only relayed, and a relayed
measurement travels with the reviewer's authority and none of its fixture
(`plan-authoring` §4). **Do not design an object marked ⬜ until its own
measurement is reproduced.**

| object | evidence | status |
|---|---|---|
| ~~`ProjectView`~~ → **bind `Workspace`** (R1) | `[M]` mine, by AST: **10 of 56** `GraphQuery` methods (7 public + 3 private) take a `project_root` arg; **178** construction sites, 176 in tests | ✅ verified, re-derived at HEAD 2026-08-16; ⛔ the TYPE is refuted, the outcome is not |
| `NodeId` (R2) | ~~`[M]` mine: `py:func` 206 / `std:doc` 94 / 527 duplicated names; one-sided map at `ast_analyzer.py:920`+`:947` vs `extractors.py:402`~~ — ✅ **PARTLY REMEDIED 2026-08-16 @ `2ddf61b` (Track 0.3)**, which landed exactly the "two lines" objection (a) said must not wait. `[M]` re-derived at HEAD: role-spelled `py:` ids **206 → 0**. The defect CLASS survives one domain over and is re-chartered below | ✅ verified, then reconciled |
| `PositionIndex` (F2) | ~~relayed: 3 implementations disagreeing on 3 of 4 probed positions~~ → `[M]` mine 2026-08-16: **3 of 4** reproduces, but the framing was wrong — **2 concepts each duplicated** (path-key ×3, span search ×2), the third cited site answers `(file) → nodes` and agrees exactly (56 = 56). Live defect found: **456 of 3530** decorated defs mis-bind | ✅ **reproduced** — see the Track 1.2 block at the top |
| `Evidence` | relayed only | ⬜ **re-derive** |
| `Diagnostic` | relayed: 285 lines of adapter, contract written twice per diagnostic (docstring similarity 0.05–0.67) | ⬜ **re-derive** |

The same applies to Track 0.2 (`_unparse_attribute`, "8 of 13 fabricated",
reported as a one-line fix) — **relayed, not reproduced.** It is the next work
item, so reproducing it is the first thing that happens, not an optional check.

### ✅ 0.1 LANDED — the ontology has a production consumer

`_infer_implements` (`merge.py`) now asks the ontology what it is allowed to
infer, instead of restating it. **Three duplicates removed**, all of which were
copies of `[edge.implements]`:

| was hardcoded in `merge.py` | now read from the declaration |
|---|---|
| `code_types = {"function", "method", "class"}` | `spec.domain` |
| `and not tgt_attrs.get("in_test_file")` | `forbid_source_attr` via `check_edge` |
| `confidence=0.7` | `spec.default_confidence` |

The check sits **at the producer**, not as a post-hoc validation pass: an edge
that should not exist must never be created, and the ontology's `implements`
enforcement is `warn`, so a validate-after design would have logged the violation
and shipped the edge anyway — strictly worse than the hardcoded filter it
replaces. No declaration ⟹ no inference. Refusals are counted and logged, because
a filter that drops silently is indistinguishable from one that never fires
([[lessons-L56]]).

`[M]` **Verified two ways, because a bit-identical diff is compatible with "the
check works" AND with "the check is inert"** (`vv-principles` #19):
- differential rebuild of the whole ORPHEUS graph — **24307 nodes symdiff 0,
  215226 edges equal, 13968 `implements` equal**, confidence values `{0.7}`
  unchanged, 0 Sphinx warnings;
- the guard is **load-bearing** — re-running the inference on the stripped
  production graph reports **195 candidates REFUSED** and re-creates exactly
  13968 edges. A direct probe refuses a test-file source *and* an out-of-domain
  source while admitting a production class.

Landed nexus `37bad88`. `[M]` 739 passed / 1 skipped.

#### ⚠ …and the same defect one level down: the EXTENSION is still inert

`[M]` 2026-08-16, prompted by the user asking why `ORPHEUS/.nexus/ontology.toml`
does not exist.

**It does not exist by design** — the ontology has two tiers: a **base**
shipping inside the package (`sphinxcontrib/nexus/ontology.toml`, the 16 node /
16 edge vocabulary) and an **optional project extension** at
`<root>/.nexus/ontology.toml` (`PROJECT_ONTOLOGY`, `ontology.py:46`). ORPHEUS
has declared nothing of its own yet, so it has only `config.toml`. Nothing is
missing.

But the question exposed that **the extension mechanism has no reachable
production effect**, which is the defect I had just "fixed" one level up:

1. `[M]` the only production call was `merge.py:219 Ontology.load()` — **no
   argument**, which the method documents as *"loads the base alone"*. Every
   call passing a root was in the test. Fixed at `c6c51dc` by threading
   `settings.project.root`.
2. ⚠ **That fix is behaviourally inert today, and I am not claiming otherwise.**
   `[M]` `_guard_redefinition` (`ontology.py:281`) enforces that an extension
   *"may add, never silently redefine a base entry"*, and `implements` is
   base-declared — so both loads yield an identical spec. The threading is
   correct (a consumer should hold the full vocabulary, and it would have been a
   live silent bug the moment this pass consulted an addable type), but the
   extension tier still has **zero** production consumers: it can only ADD
   types, and no production path consults an added type.

⟹ **The extension gets its first real consumer from the plan itself**, which is
the answer to *"is it part of what we execute?"* — **yes**, at these points:

| plan item | what ORPHEUS's `.nexus/ontology.toml` would declare |
|---|---|
| §5.1 `exercises` | the new test→code edge — `domain`, `range`, antisymmetry, `sources` (static AST vs runtime coverage). ⭐ `[M]` `tests/test_ontology.py:145` already uses `exercises` as its worked extension example |
| §5.4 `layer` | the layer vocabulary and which node types carry which value |
| §6ter Track 1 `NodeId` | the valid type set an id may be built from — the review's "missing vocabulary layer" is ontology data |
| nexus #63 | ORPHEUS's `err:ERR-NNN` nodes and a real `catches` edge — project-specific, so extension-tier by definition |

⚠ Which raises a design question to settle **before** writing one: is
`exercises` a *base* edge (nexus-wide — every project has tests) or a *project*
extension? Base, most likely — and then the `_guard_redefinition` rule means
ORPHEUS must not declare it. The extension tier is for genuinely
project-specific vocabulary (`err:`, V&V levels), not for anything nexus should
know about all projects.

### 5.3d The two tiers, derived from nexus-as-a-generic-project — *proposed*

> **User, 2026-08-16:** *"Nexus needs a relatively basic ontology.toml, and
> ORPHEUS expands it. At the very least, there should be 1 intra-connection for
> each [doc, code, test] and 1 inter-connection between [doc, code, test]…
> Nexus is dogfed to itself, so it's a good generic base to derive the base
> ontology from."*

**The method works, with one caveat stated below.** `[M]` nexus's self-graph
(3,309 nodes, 826 of them tests) is a Python project with docs and tests and no
scientific domain, so what it populates is what any project populates.

#### The criterion, measured against nexus

| cell | edges | via | verdict |
|---|---:|---|---|
| **INTRA** code→code | 2520 | `contains`, `calls`, `type_uses` | ✅ |
| **INTRA** docs→docs | 19 | `contains`, `references` | ✅ |
| **INTRA** test→test | 50 | `contains`, `calls` (helpers/fixtures) | ✅ |
| inter docs↔code | 43 | `documents` | ✅ |
| inter test→code | 1588 | `calls` | ⚠ present but **mistyped** — this is §5.1's `exercises` |
| **inter test↔docs** | **0** | — | ❌ **the missing pairwise link** |
| inter docs→test | 0 | — | not needed; the pair is covered by test→docs |

⟹ **All three intra-connections already exist generically.** Of the three
pairwise inter-connections, doc↔code holds, test↔code exists but wears the wrong
type, and **test↔docs is genuinely absent** in a generic project. That third one
is `tests` (ORPHEUS 2,748 / nexus 4) — *"this test verifies this documented
claim"* — and it is the edge the base tier most needs to name properly.

⚠ **Artefact to fix, not a finding:** `[M]` code→test shows 799 edges, of which
796 are `contains`. That is a *test module* containing test functions — the
module node does not carry `is_test` while its functions do. **The layer
assignment is inconsistent at module level**, which matters directly for §5.4's
`layer` property: it must be derived from the file, not from a per-symbol flag.

#### What ORPHEUS adds, measured

`[M]` types with **zero** instances in nexus and many in ORPHEUS:
`implements` (0 / 13,968), `equation_ref` (0 / 543), `exception` (0 / 2),
`term` (0 / 15). Everything else is exercised by both.

⚠ **The caveat on the method: "nexus does not instantiate it" ≠ "not
universal".** `term` is a Sphinx **glossary** entry — universal to any documented
project; nexus simply has no glossary. `exception` likewise. The base/extension
cut is *conceptual universality*, with nexus's graph as evidence that a type is
**reachable** generically, not as the definition.

By that test the genuinely domain-specific vocabulary is the **mathematical**
one: `equation`, `proof_object`, `equation_ref`, and `implements` in its
"code *is* this equation" reading.

#### ⛔ Two things block this being a config change

**(a) `equation` is not a vocabulary entry — it is woven through the package.**
`[M]` 220 mentions across **14 of nexus's 22 modules** (`query.py`, `merge.py`,
`extractors.py`, `directives.py`, `registry.py`, `visualize.py`, `server.py`,
`brief.py`, `ingest.py`, `_mappings.py`, …), plus `implements` at 62 across 14.
`provenance_chain` and `verification_coverage` are *named for* it. So "move the
math vocabulary to the extension tier" is a **reshape**, not a TOML edit — and
it is the deepest instance of the review's *"still shaped by what it was"*: nexus
was built for ORPHEUS and its base vocabulary is ORPHEUS's domain.

**(b) The extension tier can ADD a type but cannot EXTEND a base one.**
`[M]` `_guard_redefinition` (`ontology.py:281`) refuses any project redefinition
of a base entry. So ORPHEUS can declare an `equation` **node**, and then cannot
say *"`equation` is a valid target of `implements`"* — that would redefine the
base edge's `range`. The only escape today is the `ANY = "*"` wildcard
(`ontology.py:49`, already used by 7 base entries), which buys openness by giving
up all range checking.

⟹ **The mechanism needs a third verb.** Today: `add` (new name) and `redefine`
(refused). It needs **extend** — a project may *widen* a base edge's
`domain`/`range` with its own node types, never narrow it or change its meaning.
Widening is **monotone**: anything the base admitted it still admits, so no
consumer's expectation can be invalidated — which is why it is safe where
redefinition is not, and is the property to assert in a test
(`base.admits_target(t) ⟹ extended.admits_target(t)`, ∀ t).

✅ **FILED as nexus #69** and approved for implementation (user, 2026-08-16).
Scope is the **mechanism only** — the vocabulary split itself is larger and
blocked on (a).

#### Proposed base (nexus) — *not yet ratified*

- **nodes** — code: `module`, `class`, `function`, `method`, `attribute`,
  `data`, `exception`; docs: `document`, `section`, `term`; resolution:
  `external`, `unresolved`; meta: `tag`
- **edges** — intra-code: `contains`, `calls`, `imports`, `inherits`,
  `type_uses`; intra-docs: `contains`, `references`; intra-test: `contains`,
  `calls`; **inter: `documents`** (docs→code), **`exercises`** (test→code, new),
  **`verifies`** (test→docs, generalising today's `tests`), and `implements`
  (code→docs) in the weaker universal reading *"this code realises this
  documented statement"*

#### Proposed ORPHEUS extension — *not yet ratified*

- **nodes**: `equation`, `proof_object`, `err` (the `ERR-NNN` catalogue, #63)
- **edges**: `equation_ref`, `catches` (#63)
- **extends**: `implements`.range += `equation`; `verifies`.range += `equation`
  — ⛔ **blocked on (b) above**

⚠ Sequencing: this sits **after** Track 1's `NodeId`, because the valid type set
an id may be built from is the same vocabulary question. Attempting the tier
split before then means doing it twice.

## 7. Open questions

**Answered 2026-08-16:** §5.3's shape is DEFERRED; §5.4 is APPROVED; §3.2 is
filed as **#68** and goes FIRST.

**Answered 2026-08-16 (second pass):** §5.1 is ACCEPTED; §5.3's store is a
**DB** — only the CI story is deferred.

Still open:

1. **§5.1 scope** — retype `calls`→`exercises` only for test→code, or take the
   double-duty principle further in the same pass (`contains` at 35.4% cross,
   the `id`'s type prefix)? The prefix one is implicated in #55's fragmentation,
   so they may not be separable.
2. **§5.1 identity** — is static `exercises` (test calls production) the same
   edge as runtime `exercises` (#57, from coverage contexts), distinguished by
   `confidence`/`source` as `implements` already distinguishes declared 1.0 from
   inferred 0.7? *Leaning yes.*
3. **§5.2 boundary** — `external` demotes safely, but `unresolved` holds real
   ORPHEUS symbols. Does #68 (and then #55/#16) shrink `unresolved` enough that
   the question dissolves, or is a `terminal` flag needed regardless?
4. **§6.2 ranking** — `limit_per_type` truncates by insertion order. What should
   it rank by? Degree is available and **wrong** (it favours hubs, which is the
   opposite of relevance). Candidates: same-module first, `type_uses`/`inherits`
   over `calls`, runtime-confirmed over static. **Unmeasured — this needs an
   eval, not an opinion**, and the eval is the §1 acceptance test.

---

## 8. Insight log — what this conversation produced

| date | insight | where it landed |
|---|---|---|
| 2026-08-16 | Tests are 27.1% of hub payload, not the cause | §2.4 |
| 2026-08-16 | The halo is TERMINAL — 49% of edges, zero reachability value | §2.2 |
| 2026-08-16 | Layers are already in the edge type, except test→code `calls` | §2.3 |
| 2026-08-16 | The 3-way split is the wrong seam (76.8% crossing) | §2.2 ⛔ |
| 2026-08-16 | Lifetime, not layer, is what forces a file boundary | §4 |
| 2026-08-16 | `unresolved` is FOUR unrelated defects wearing one type | §3 |
| 2026-08-16 | 1,282 phantoms are prose math symbols and citations, not code | §3.2 |
| 2026-08-16 | Payload is 5× its information content | §2.4, #67 |
| 2026-08-16 | Separate the duties, wherever a node or edge does two jobs | §5.1 ⭐ |
| 2026-08-16 | 83% of the prose leak is inline math minted as an equation LABEL — **not** the py-domain cause I first claimed | §3.2 ⛔, #68 |
| 2026-08-16 | The phantoms pollute the namespace `verification_coverage` counts over | §3.2 ⚠ |
| 2026-08-16 | DB-vs-JSON is a **dominance** argument, not a trade-off — JSON wins in no configuration | §5.3 |
| 2026-08-16 | The observed store sits INSIDE the build output it was designed to outlive | §5.3a ⛔ |

### Surprises this conversation produced (`plan-authoring` §-loop)

Two of my own claims were refuted by the next measurement, both stated with
more confidence than they had earned:

1. **The constructor/type-prefix hypothesis** for the `AngularBoundaryFlux`
   phantom. `[M]` `py:function:…AngularBoundaryFlux` does not exist; the prefix
   class totals 359 edges. The real cause was **inheritance** (MRO not walked).
2. **The `default_role` hypothesis** for the prose leak. `[M]` `default_role` is
   not set in `docs/conf.py`; the ids are `math:equation:*`, not py-domain.

⭐ Both share one shape, and it is worth carrying: **I named a plausible
mechanism from the SHAPE of the symptom before querying the actual node ids.**
The id was one query away in both cases, and in both cases it named a different
subsystem than the one I had reached for. ⟹ *when a defect is identified by a
node's NAME, read its full ID before proposing the cause* — the id carries the
namespace, and the namespace is the subsystem.
