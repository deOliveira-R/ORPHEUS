# A red partitions the suite, and a coverage claim can be falsified

---

## ⏸ COMPACTION POINT — 2026-08-15

**Read this section, then `git log`, then §6. Do not reconstruct state from a
conversation summary.**

### What LANDED (verify with `git merge-base --is-ancestor <hash> HEAD`)

**`~/git/sphinxcontrib-nexus`, branch `feat/config-and-ontology`** (4 commits,
`[M]` 720 passed / 1 skipped):

| hash | outcome |
|---|---|
| `47693d8` | a project's graph settings and vocabulary get a home — `project.py`, `ontology.py`, `ontology.toml`; Python floor 3.10 → 3.11 so `tomllib` is stdlib |
| `dd36241` | the extension, the CLI and the server all read the same file |
| `0a209dc` | two settings were reached by `getattr` and my census missed them |
| `206f2d6` | six `--db` defaults survived the sweep in a different spelling |

**ORPHEUS, branch `chore/nexus-project-config`** (3 commits, `main` untouched):

| hash | outcome |
|---|---|
| `046f0a5e` | the graph's settings live in `.nexus/`, declared once |
| `e2b25e3c` | retire `conf.py`'s nexus options, and measure the one knob left |
| `f59d107e` | this compaction point + `lessons` L55 |

⚠ `[M]` every hash above verified with `git merge-base --is-ancestor <h> HEAD`
in its own repo. An earlier draft of this table carried `89f4e0c3` for the
retirement — the pre-amend hash, dead the moment trailers were added. A
compaction point with a wrong hash is worse than none, so re-verify rather
than copy forward.

### Measured state — the transition is CLOSED and neutral

`[M]` Three full ORPHEUS rebuilds, each compared on five dimensions
(node id set / node rows / edge multiset / node attrs / edge attrs):

| comparison | result |
|---|---|
| conf.py-driven **vs** config-driven | **0** difference ×5 |
| conf.py options present **vs** RETIRED | **0** difference ×5 |
| inference ON **vs** OFF | ONLY `implements` changes (13968); 0 node difference |

Only `build_time`/`provenance` ever differ — which is the differential's own
positive control, proving it is not comparing a file with itself.

`[M]` Current graph: `docs/_build/html/graph/graph.db`, **24307 nodes /
215226 edges**, resolved with no `--db` anywhere (CLI from any subdirectory,
and the MCP server via `--project-root`). Four runtime traces re-bound at
`graph/traces/` after a copy AND two rebuilds.

⚠ `docs/_build/html/_nexus/` (frozen pre-config baseline) and
`graph_noinfer/` (the inference experiment) still exist — deletion was
declined at the permission prompt. Both gitignored; traces already copied
out, so `_nexus/` is safe to remove.

### Issues filed (nexus repo unless noted)

**#55** call resolution mints a phantom per receiver SPELLING · **#56** ingest
can join nothing and exit 0 · **#57** coverage ingest discards dynamic
contexts · **#58** `bridges`/`communities` time out at 24k nodes · **#59** an
empty result is indistinguishable from "not applicable" · **#60** the
antisymmetric `subject` relation (the blocker) · **#61** lift markers from a
pytest COLLECTION manifest · **#62** `retest` truncates at `max_depth=3` ·
**#63** `catches` is a string while `verifies` is an edge · **#64** the config
surface (this work). ORPHEUS **#358** carries the synthesis comment; **#377**
is the `RigidMotion.determinant` perf finding.

### ✅ #56 CLOSED 2026-08-16 @ nexus `5796c6d` — an ingest that joins nothing says so

`NodeBinder` + `JoinLedger` now own the key space, the scope test and the
drop accounting for all three backends. `[M]` the same ORPHEUS artifact goes
**0 → 1691** bound nodes (158 of 158 coverage keys were relative); a zero-join
exits 1, names its root, and is **not stored**. `[M]` 739 passed / 1 skipped
(was 726); both structural changes mutation-verified, one red each.

⛔ **Two of the issue's own premises were REFUTED** — recorded on the issue
rather than routed around. (a) The suggested fix said to read the rundir from
`coverage json`'s `meta`; `[M]` meta carries only
`format/version/timestamp/branch_coverage/show_contexts` and records it
**nowhere**, so it must be supplied — hence `--root`. (b) The issue said the
multi-prefix trap was "filed separately"; `[M]` no such issue exists, so it
was folded in. Also: `coverage json --branch` is rejected — `--branch`
belongs to `coverage run`.

Two further instances of the same defect fell out: `--kind viztracer` reached
the **coverage** ingester (a binary dispatch over three registered choices),
and `len(calls) or len(coverage)` read `nodes: 0` for a *successful*
viztracer run. Filed while here: **#65** (the tool reference documents 29 of
40 tools, missing both worktree-hazard tools) and **#66** (`--db` declared 42×
per-subparser; `status` rejects `--project-root`).

### ▶ RESUMES AT: one method is one node, however its receiver is spelled (nexus #55)

**Goal (outcome, not mechanism).** A call to `X.foo()` binds to the same node
regardless of how the receiver was spelled at the call site, so `callers`
answers a question instead of returning `[]`.

**`[M]` unstarted** as of 2026-08-16 — verify before designing; the fix for
#56 touched `runtime.py`, `cli.py` and `server.py`, and **not** the AST call
resolver.

**Then** #57 (`exercises` from coverage contexts). The order is deliberate:
instrument, then resolver, then capability — measuring with a broken
instrument is worse than not measuring, because the result gets recorded.

⚠ Read §2.2 here first. The diagnosis is *a phantom per receiver spelling* —
one method fragmenting into five nodes, 20 948 of 121 788 calls landing on
`unresolved`. Nexus already solved this class for **doc** references in #37
and never applied it to the call resolver; that prior art is the starting
point, not a general "name resolution is hard" problem.

⚠ And #55 is where [[lessons-L56]] bites hardest: `callers` returning `[]`
is the flagship instance of "an empty answer and a broken one look alike"
(nexus #59). Fixing the resolver without also making `[]` self-describing
leaves the next reader with the same unanswerable output.

---

**Campaign home for ORPHEUS #358.** Work happens in **two repos**:
`~/git/sphinxcontrib-nexus` (the graph capability) and `ORPHEUS` (the authoring
that uses it). Nexus is ours, so this is a folder change, not a hand-off.

**Status 2026-08-15.** Design complete, nothing built. Nexus upgraded
`0.16.1 → 0.17.0`, graph rebuilt clean (`-E`, 0 warnings). All numbers below are
`[M]` on that rebuilt graph unless marked otherwise.

---

## 1. Goals — stated in the domain, separately from any means

**G1.** A red test partitions the suite into *valid / flip point / invalidated /
not-run*, so "fix one thing and resume" is a defined operation rather than a
re-run of the tree.
**Done when:** `invalidated_by(<the root test>) == {}` — a forward-only
assertion passes on the symmetric relation that already exists, so only a
genuine partial order can satisfy this.

**G2.** A coverage claim can be *falsified*. Today "test T verifies equation E"
is unfalsifiable: nothing records whether T ever executed E's implementation.
**Done when:** the claimed set and the exercised set are separately queryable
and their difference is computable.

**G3.** A frozen reference can be adjudicated *before* it is re-baselined —
"which non-`RECORD` test also constrains this subject?" is a query, not a hunt.
**Done when:** the query returns a candidate list, and its docstring states the
honest limit (it finds subjects with *no* pin, never a *blind* one).

⚠ **Explicitly NOT a goal.** Finding gates that cannot fail. Only mutation does
that. #358 says so, `vv-principles` #17–#24 back it, and no deliverable here may
claim that credit.

---

## 2. The finding that reorders everything

Two sub-agents worked opposite halves with no shared context and hit the same
wall:

> **`qa`:** it measures **co-execution, not co-constraint**.
> **`test-architect`:** `calls` gives "A and B both touch P" — **symmetric**,
> hence not a partial order. #358 needs antisymmetric "B rests on A".

`[M]` A naive cone invalidates **31–36 %** of the suite on one red and is *not
filterable* — blocking `external` traversal moves it 0.6 pp.

⟹ **#358's premise that the dependence DAG can be derived from the call graph is
structurally wrong, not merely incomplete.** `exercises` (what ran) is necessary
for G2 and **insufficient** for G1. G1 additionally needs a **subject** relation.

### 2.1 Static blindness, measured twice, independently

| probe | static | runtime |
|---|---|---|
| tests reaching `Quadrature.ordinate_permutation` | **0 / 21** (same at depth 4/8/20/50) | 7 / 21 |
| selection for `OperatorSum` | **13** | **229** (94.8 % miss) |

`[M]` All four of `OperatorSum.{apply,inverse,apply_transpose,assemble}` have
**zero** static in-edges while all fire. The static graph sees the operator
algebra **at construction, never at use**.

⭐ **And that is a consequence of Cardinal Rule 2, not an accident.**
`coding-elegance` Pattern 1 spells every operation as a dunder — exactly what a
syntactic resolver cannot follow. **The better the architecture gets, the
blinder the call graph becomes.**

⛔ `[M]` **0.17.0 does NOT fix this** — `callers` on that method is still
`total: 0` post-upgrade. 0.17.0 repaired *documentation* reference resolution,
not *call-edge* resolution.

### 2.2 Root cause, measured — it is not a dropped edge

The edge is **minted onto a phantom keyed by the caller's local variable name**:

```
quadrature.ordinate_permutation  unresolved  calls-in= 4
quad.ordinate_permutation        unresolved  calls-in=16
q.ordinate_permutation           unresolved  calls-in=16
good./broken.ordinate_permutation unresolved calls-in= 2 each
…Quadrature.ordinate_permutation method      calls-in= 0
```

One method, five phantoms, 40 edges. Corpus: **20 948 of 121 788 `calls`
(17.2 %)** land on `unresolved`, over 3 580 phantoms.

This is the **same defect nexus fixed for doc references in #37** ("it retargets
the *edge*, not the node") — applied to the reference resolver and never to
`_resolve_call_target` (`ast_analyzer.py:537`).

| unresolved call shape | count | share |
|---|---|---|
| simple `recv.method` | 13 830 | 66.0 % — **897 recoverable** by a one-hop param-annotation join |
| bare name | 3 851 | 18.4 % |
| chained `a.b.method` | 3 267 | 15.6 % |
| `self.x.method` | 0 | **dropped upstream** — returns `None`, never even unresolved |

⚠ **Do not overstate the annotation join.** 897 edges is **6.5 %** of the simple
class. It *does* cover the motivating case (all 4 callers carry
`type_uses[param="quadrature"] → Quadrature`), and it is build-time and always
fresh — but it is a **complement to** the runtime overlay, never a substitute.

---

## 3. The label vocabulary — coarse

Three layers, and **four** arrow classes (test→code is two questions, not one).

| arrow | asserts | answers | state `[M]` 2026-08-15 |
|---|---|---|---|
| CODE→MATH `implements` | this code **is** this statement | provenance | ships — **0 declared**; 13 968 inferred @0.7, **87.5 % single-token** |
| TEST→MATH `tests` | this test **verifies** this statement | *claimed* coverage | ships — 2 748, all claim @1.0; **50.5 % from multiplicative file-level markers** |
| **TEST→CODE `exercises`** | this code **ran** under this test | *measured* coverage (G2) | **MISSING** |
| **TEST→CODE `subject`** | this code is what the test is **about** | the cone (G1) | **MISSING — the blocker** |
| MATH→MATH `discretizes`/`derives_from`/`approximates` | the statement spine | derivation order | ships 0.17.0 — **0 authored** |
| CODE→CODE `calls`/`imports`/… | reachability | — | ships — reachability ≠ dependence (§2) |

## 4. The label vocabulary — fine

### 4.1 `claim_kind` — provenance of the EXPECTED VALUE

L51's CONSTRAINT/RANKER/DIAGNOSTIC is **refuted for tests**: every collected
pytest test is a CONSTRAINT by construction (pass/fail is the only outcome), so
the partition has one non-empty cell.

| coarse | the expected value comes from | if it reds | fine |
|---|---|---|---|
| **THEOREM** | a law holding for every admissible input | **voids every other claim on that subject** | `symbolic`, `invariant` |
| **REFERENCE** | a structurally-independent external route | the implementation disagrees with truth | `analytical`, `convergence` |
| **RECORD** | what this code produced on a chosen day | something *changed*; no information about which side is right | `empirical`, `regression` |

`[M]` It must be **declared, not derived**: `assert_allclose` appears in **218
files** serving all three.

**Payoff (G3), mechanical:** *every RECORD subject must also carry a THEOREM or
REFERENCE test.*

### 4.2 The other fine labels

| label | kind | goes on | why |
|---|---|---|---|
| `scope` | edge attr on `tests` | module- vs function-level marker | `[M]` `vv_level` resolves on **1524 of 5273** because **254 files** use module-level `pytestmark`; a 9-equation blanket marker scores `confidence=1.0` identically to a single-purpose L0 gate |
| `role` | node attr on code | `production`/`reference`/`oracle`/`diagnostic`/`harness` | enables the structural-independence check (`vv-principles` #1) |
| `layer` | node attr on module | the existing L0/L1/input/L2/L3 contract | a **move** of `tests/test_layer_imports.py:34`, not a mint; extends the gate to `tests/` for free |
| `error` node + `catches` edge | node + edge | ERR catalog | `[M]` `gaps` reports `error_catalog_size: None`; `catches` is a filterable attr on 239 nodes, **not traversable** |
| `prf_type` | node attr | statements | ✅ 0.17.0 — *what kind of statement*, **NOT** claim kind |

⛔ **REJECTED — a `quantity` node.** G3 wants "which test constrains the same
quantity". Minting it means hand-labelling — the thing that killed `cap()`. Once
`exercises` exists, "same quantity" ≈ "exercises the same code node", refined by
the equation label where one exists. Derive it.

⛔ **REJECTED — `body_shingles` as a marker-staleness detector.** `[M]`
bit-identical under tolerance, budget, re-baselined-expectation and fixture-arg
edits (it normalizes `Constant→"C"` by design) — blind to **every** decay cause,
in the flattering direction.

---

## 5. The honest ceiling

`[M]` For the worked equation: **CLAIMED 21 → EXERCISED 7 → ASSERTED ≤2 →
MUTATION-VERIFIED 0.**

⛔ **All 7 exercisers are identically connected in any exercised-coverage graph,
yet 5 cannot fail for a permutation error.** No edge quality separates rungs 2
and 3. **This design reaches rung 2.** Nothing in it may claim rung 3.

(#334's "~2 exercise it" is refuted in magnitude — 7, so 3.0× not 10× — but is a
good estimate of the *asserting* rung.)

---

## 6. Work — nexus side (`~/git/sphinxcontrib-nexus`)

Ordered. **P0 fixes the instrument before funding any measurement** — a
measurement taken with a broken instrument is worse than none, because it is
recorded.

| id | outcome | `[M]` evidence | touches |
|---|---|---|---|
| **P0-1** | an ingest that joins nothing says so | printed `nodes: 0`, **exit 0, no warning**; `coverage json` emits relative keys against an absolute index, dropped *upstream* of the `unresolved` counter. Absolute keys → 2892 | `runtime.py` |
| **P0-2** | scope is a set of prefixes, not one | root 21.9 % join, `orpheus/` 67.2 %, `tests/` 88.6 % — and `tests→orpheus` edges exist **only** in the root run | `runtime.py`, `cli.py` |
| **P1-1** | a call binds to the symbol, not to the caller's spelling | §2.2 — 5 phantoms/1 method; 20 948 edges (17.2 %) on `unresolved`; 897 recoverable by the annotation join already in `type_uses[param=…]` | `ast_analyzer.py:537` |
| **P1-2** | per-test **`exercises`** | coverage `dynamic_context` emits it; `grep -ci context runtime.py` → **0**. Prototyped on nexus's own `build_node_index`: 23/23 contexts, 1353 pairs, **~15 lines** | `runtime.py:overlay_coverage` |
| **P1-3** | markers arrive **resolved**, not re-parsed | collection-manifest sidecar `runtime-ingest --kind pytest`. Fixes 254 module-level `pytestmark` files, conftest precedence, granularity, and gives `claim_kind` a home — one mechanism. ⭐ **supersedes** a build-time `nexus_marker_attributes` config: AST cannot see resolved precedence | `runtime.py`, `cli.py` |
| **P1-4** | **`subject`** — an antisymmetric "rests on" | **the G1 blocker**. Gate: `invalidated_by(test_root) == {}` | `ast_analyzer.py`, `query.py` |
| **P2-1** | an incomplete answer says it is incomplete | `callers` → `total: 0` while 40 edges sit on phantoms; `dead_functions` then flags a live method | `query.py` |
| **P2-2** | `retest` walks to a fixed point | `[M]` hard-coded `max_depth=3`; 18 tests (5.1 %) reach `warn_if_unconverged` at depth 4 and are reported `safe_to_skip`. Symbol-dependent — `solve_sn` converges at depth 2, so the obvious spot-check shows nothing | `query.py:retest` |
| **P2-3** | `error` nodes + `catches` edge | `gaps` → `error_catalog_size: None` | `query.py`, `registry.py` |
| **P3** | mutation verdict as staleable data | ⛔ **mechanism refuted** (§4.2) — needs a new one before it is scheduled | — |

## 7. Work — ORPHEUS side (authoring, not capability)

Several "missing labels" are **not missing from nexus; they are unused by us.**

- statement relations `discretizes`/`derives_from`/`approximates` — **0 authored**
- `.. implements::` — **0 declared**, against 13 968 inferred @0.7 / 87.5 % single-token
- `claim_kind` declaration per test
- `sphinx-proof` — not installed; `proof_object`/`prf_type` unavailable until it is
- extend the layer gate to `tests/` (`test_layer_imports.py:34` scans `orpheus/` only) — #358's own suggested first step, still undone

---

## 8. Method notes for whoever picks this up

- ⚠ **`is_test` ≠ collectable.** `[M]` `is_test` = 7305 nodes, but restricted to
  `function`+`method` it is **5273 — exactly pytest's collected set**, 100 % join
  both directions (0 in either gap), 7.0 s to collect. Use `in_test_file` (9530)
  for "lives in the test tree". Both are 0.17.0.
- ⚠ **SQLite double quotes bind to COLUMNS.** `a.key="source"` joined against
  `edges` resolves to `edges.source` and silently returns empty. Single-quote
  string literals. (Cost me one false "no provenance" reading.)
- ⚠ **Python buffers stdout off a tty** — an empty log is not a dead run.
- ⚠ **zsh `nomatch` aborts the whole command line**, so a failed glob can
  silently skip the `ls` that follows it. (Cost me a false "no nexus checkout".)
- ⚠ **cProfile's `co_filename` is absolute**, so it does *not* hit P0-1's
  relative-key bug. **A clean cProfile ingest is not evidence that coverage
  ingest works.**

## 9. Refuted, with the structural reason

| candidate | why it fails |
|---|---|
| a finer `l0`–`l3` | ordinal quality scale ≠ dependence relation; two `l1` tests in unrelated packages are incomparable. Category mismatch (#358's own ruling) |
| backfilling `cap()` | 4 of 456 files, 0 outside `tests/sn/`, plans archived, being replaced |
| a hand-maintained marker DAG | a missing hand-edge *is* the silent-promotion failure |
| topological **sort** | discards the branching the partition is made of |
| static-`calls` cone | §2 — symmetric, and 0 % recall on the measured relation |
| `implements` as the equation→code relation | 100 % inferred, 87.5 % single shared token |
| `provenance` as the coverage source | answers "what else is on this page" |
| `verified` status | `query.py:1479` sets it iff `len(tests) > 0`, **no floor** — `[M]` 351 of 692 "verified" equations have no declared test |
| deriving `claim_kind` from the graph | nothing observes what an assertion *computes*; `assert_allclose` in 218 files serves all three |
| cProfile as the exercised-coverage source | process-global union; safe for cone invalidation, **unsafe** for a claimed-vs-exercised audit |
| `body_shingles` as staleness detector | §4.2 |
| a `quantity` node | §4 — derivable once `exercises` exists |
| "80.3 % of production nodes have zero tests" | measures the resolver, not coverage. Nearly shipped; caught only by implausibility |
