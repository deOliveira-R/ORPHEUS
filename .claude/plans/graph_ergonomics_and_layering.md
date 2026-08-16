# A graph query returns what an agent needs to ACT, not everything the graph knows

**Status: PARTLY APPROVED.** Sections marked *proposed* are still hypotheses;
sections marked `[M]` carry the command or query that produced them.

## ⏸ COMPACTION POINT — 2026-08-16, before any code is written

**Nothing in this plan has been implemented.** What exists is measurement, four
filed issues, and two user decisions. Re-verify against the tree before acting:
the numbers are instance counts from one build of
`docs/_build/html/graph/graph.db` and move on every rebuild.

### User decisions, 2026-08-16

| decision | status |
|---|---|
| Fix the prose-symbol leak (**#68**) **FIRST**, before the rest | ✅ APPROVED — next work item |
| `layer` as a node property (§5.4 / §6.1) | ✅ **APPROVED** |
| §5.1 retype test→code `calls` as `exercises` | ✅ **ACCEPTED** |
| §5.3 observed state is a **DATABASE**, not a JSON blob | ✅ **CHOSEN** — only the CI story is deferred |
| §5.2 halo demotion | *proposed*, not yet ruled on |

⚠ I recorded the first two of these wrongly on 2026-08-16 and the user
corrected both: `exercises` had been accepted ("a good suggestion"), and §5.3's
deferral was only ever about **CI**, not about DB-vs-JSON. Kept visible per
`plan-authoring` §3 — a plan that silently repairs its own misreadings teaches
the next reader nothing about how they happen.

### ▶ RESUMES AT: the equation namespace holds only real equations (nexus #68)

**Goal (outcome, not mechanism).** `math:equation:*` contains the 903 declared
equations and nothing else, so any count, traversal or coverage denominator over
that namespace is trustworthy.

**`[M]` unstarted** — verify before designing. Work is in
`~/git/sphinxcontrib-nexus` (ours; a folder change, not a hand-off), on branch
`feat/config-and-ontology`, which already carries the config/ontology layer and
the #56 runtime fix.

⚠ **Read §3.2 before designing, including its ⛔.** My first stated cause was
wrong (I blamed `default_role`, which is not set) and the real one — inline
`:math:` minted as an equation *label* — implies a different fix in a different
extractor. The issue body carries the corrected taxonomy.

⚠ `[M]` **run the nexus suite with the NEXUS venv**, not ORPHEUS's: the latter
lacks the optional `sphinx-proof` and `test_proof_relations.py` importorskips at
module level, dropping **36 tests** from collection — a green 691 that looks
identical to a green 727.

### Where the work is tracked

Filed from this conversation: **#67** (payload is 5× its information content),
**#68** (the prose leak — first). Adjacent, filed earlier the same day: **#65**
(the tool reference documents 29 of 40 tools), **#66** (`--db` declared 42×
per-subparser). Pre-existing and relevant: **#55** (call-resolution
fragmentation), **#16** (annotation-mediated dispatch), **#57** (per-test
attribution), **#58** (`bridges`/`communities` do not complete), **#59** (an
empty result is indistinguishable from a broken one), **#60** (a cone needs an
antisymmetric relation), **#20** (filtered subgraph → graphviz render).

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

## 2. Measured baseline — `[M]` 2026-08-16, `docs/_build/html/graph/graph.db`

Graph: **24,307 nodes / 215,226 edges**. Queries below are against that file;
re-measure after any rebuild, since these are instance counts, not schema facts.

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

### 5.3a ⛔ The store is in the WRONG PLACE today, and choosing a DB is when to fix it

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
