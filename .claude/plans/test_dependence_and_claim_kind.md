# A red partitions the suite, and a coverage claim can be falsified

---

## ⏸ COMPACTION POINT #5 — 2026-08-18 · the ledger SHIPS, the declaring path is guarded, and CI can be believed again

⚠ **Everything below this section is HISTORY unless a hash says otherwise.**
Both repos on `main`, **pushed**, 0 dirty / 0 unpushed:
nexus **`253581f`**, ORPHEUS **`c3349828`**.
`git merge-base --is-ancestor <hash> main` is the authority.

⚠ **`/mcp` needs reconnecting** — the server serves the code it imported at
startup, and five commits since CP#4 changed the reply shape.

### ▶▶ WHERE THE RULED SEQUENCE STANDS

The user's 1 → 2 → 3 ruling (CP#4) still holds, and the licensing framing
this session added a **1b** to it. Titled as outcomes (`plan-authoring` §1).

| # | outcome | state |
|---|---|---|
| 1 | `retest` answers from EVIDENCE, not popularity | ✅ **LANDED** — see below |
| **1b** | **the code side of every claimed equation is DECLARED, not guessed** | ⏳ in progress: guard landed, inventory taken, `:kind:` SHIPPED + first 11 declared; **151 remain — ORPHEUS #382** |
| 2 | a capture wide enough to matter exists | ⏳ not started |
| 3 | a test declares what it is about (`.claude/plans/test_architecture_redesign.md`) | ⏳ not started |

⭐ **1b exists because of the user's standing ruling**, verbatim: *"to
demonstrate correctness we must be able to demonstrate the ledger of tests,
which tests exercise which code, and why we can trust them by dependency…
bullet-proof V&V at a licensing code standard"*. That is what turned step 1
from "improve `retest`" into "build the ledger and make its claims
falsifiable", and what makes documentation-via-ontology campaign work rather
than a filing.

### ✅ nexus#85 — DONE, 2026-08-18. nexus `8d87abd` · ORPHEUS `bb075c93`

⛔ **This section SUPERSEDES the "▶ NEXT: nexus#85" pointer that stood here.**
Original text, kept per §3 because its reasoning is still the reasoning:
*"`[M]` existence-checked at CP#5 time … 0 hits for `nothing` in
`directives.py` or `ontology.toml`. It does not exist; #85 is still a sketch,
still OPEN. ⛔ Why it blocks: 150 of the inventory's 332 rows are NOTHING
verdicts with nowhere to land."*

**What shipped.** `.. no-implementation:: <label>` `:kind:` — its OWN
directive, not a `:nothing:` option on `.. implements::` (user ruling,
2026-08-18: it writes no edge, so `:by:` + `:nothing:` together becomes
unspellable rather than a runtime refusal). It sets `no_implementation_kind`
on the statement; `_infer_implements` reads that alongside declared edges, so
there is ONE stand-down, not two suppression paths. Closed vocabulary —
`identity` / `law` / `canonical-form` / `definition` — extendable per project
via `[extend.attribute.…] values` (user ruling, same turn).

**`[M]` the measured effect**, ORPHEUS graph rebuilt at `bb075c93`:

| | before | after |
|---|---:|---:|
| corpus `implements` edges | 13450 | **13206** (−244) |
| labelled set answered | 45 of 56 | **56 of 56** |
| inferred edges left on it | 244 | **0** |
| coverage `documented` / `no_implementation` | 188 / 0 | **177 / 11** |
| ⭐ RULE half `guesses_on_unimplementable_equations` | 239 | **239** |

The last row is the anti-cheat invariant and it is the one that matters: the
scorer's RULE half REIMPLEMENTS the inference rather than calling it, so a
declaration cannot flatter it. ⚠ Its `false_positives` DID move, 1305 → 1314,
and that is not a regression — the candidate pool is the symbols a page
DOCUMENTS, and these blocks cross-reference symbols in order to say they are
NOT the implementer. The scorer's own comment predicts that mechanism.

⚠ **§10, third instance this campaign, and the first inside the instrument
built to be the honest one.** `implements_ground_truth.py`'s corpus half
derived "declared" from having a non-inferred `implements` edge — a PROXY
this work removes. It would have read `45 of 56` forever however thoroughly
the eleven were declared, i.e. failed in the direction that reads as
unfinished work. Fixed in the same change; it now reports the two kinds of
answer separately.

### ✅ 1b's MAIN PUSH LANDED — 2026-08-19, ORPHEUS `596093fa`

⛔ **This SUPERSEDES the "▶ NEXT: the 151 remaining" pointer that stood here**,
which named the wrong population. `[M]` only **2 of those 140** carry any
coverage claim (6 claims); they are corpus hygiene. The population the
licensing argument turns on was a different 94 — the equations that carry
claims but had no implementer — and it is now done.

| | before | after |
|---|---:|---:|
| claims blocked by "equation has no implementer" | 976 over 94 eq | **3 over 3** |
| declared `implements` edges | 85 | **388** |
| ⭐ anti-cheat `guesses_on_unimplementable_equations` | 239 | **239** |

⭐⭐ **THE CONVENTION, ruled by the user 2026-08-19**: `implements` means
**every site that executes the equation's arithmetic**, not only the canonical
or assembling one. The measurement that forced it: `mg-balance`'s 70 claims
split 26 SN / 26 CP / 12 homogeneous / 5 MC / 1 MoC over **four independent
transcriptions**; a canonical-only declaration refutes **58 of 70** — tests
that DO exercise an in-scatter sum, just not that one. ⟹ narrow is not the
conservative choice, it is the one that manufactures false refutations.

⭐⭐ **METHOD, and it is the transferable half**: six `explorer` agents, one
page cluster each, deriving from the tree — **92 of 94 declarable**, 302 of 303
targets passing a mechanical resolve+ontology check. ⛔ The Haiku inventory
called 67 of these same 94 "nothing implements it"; it is refuted from six
independent directions and its header says so. **A fanout is an inventory
engine, not a judge** — and the judging is cheap when the agent reads the page.

⚠ **The authored-but-inert surfaces**, in the order they actually pay:
1. the **claiming test module's `pytestmark` comment block** — names the gate
   and often the symbol; found independently by three agents, and it is what
   works when the page has no rationale;
2. **labels declared IN THE CODE** — `orpheus/derivations/discrete/sn/balance.py`
   names labels term-by-term (`(Eq. dd-recurrence)`), and
   `sn/sweep/psi_half_angle_seed`'s docstring is called *the canonical
   algebra-of-record* by the page itself. Two greps resolved four labels.
3. `.. (vv-status rationale)` — ⛔ `[M]` only **2 of 94** equations have one
   within 12 lines of their label. It is dense on pages an earlier V&V
   campaign visited and absent elsewhere; do NOT brief it as universal
   (`plan-authoring` §2, logged).
4. `equation_labels=` on a `VerificationCase` — 58 refs / 16 files / 44 labels.
   ⛔ NOT an `implements` channel: its own docstring says *"IDs this reference
   **exercises**"*, i.e. a `tests` edge, and those already exist via
   decorators. Its docstring also claims *"Nexus builds test ↔ equation edges
   from these"* — `[M]` **0 hits** for `equation_labels` in the whole nexus
   package.

### ▶ NEXT

1. **Step 2 — the wide capture** is RUNNING as of 2026-08-19 (17 slices,
   `-m "not slow"`, contexts). ⛔ Two bugs cost v1 a 27-min slice: the CLI flag
   is `--note` not `--command`, and the script deleted the report before
   checking rc. Size is NOT a blocker — a 2633 MB report ingests in ~10 s.
   Keep the `.coverage` SQLite (3.5 MB) until ingest returns 0.
2. **#385** — `[M]` 0 of 4 MC solver kernels is imported by any test; four L0
   gates replicate the logic inline. Those claims WILL adjudicate REFUTED once
   the capture covers `tests/mc`, and that verdict is correct.
3. **#383** — `phase-f-q-bar-twin-forms` asserts an equivalence whose second
   twin was retired. The only one of the 94 left undeclared, deliberately.
4. **#384** (doc drift, 8 items + a duplicate label), **#386** (MoC: twin
   denominator for `t_s^eff`; trig on the angle chart), **#382** (the 140
   hygiene rows).


---

## ⏸ COMPACTION POINT #4 — 2026-08-18 · the cone's edges are honest, the graph was 33 % duplicate, and the SEQUENCE is ruled

⚠ **Everything below this section is HISTORY unless a hash says otherwise.**
Both repos on `main`, **pushed**, 0 dirty / 0 unpushed:
nexus **`6e469e9`**, ORPHEUS **`3f3cc9ca`**.
`git merge-base --is-ancestor <hash> main` is the authority.

⚠ **`/mcp` needs reconnecting** — the server serves the code it imported at
startup, and `6e469e9` changed the reply shape.

### ✅ STEP 1 LANDED — 2026-08-18, nexus `a0f038e` (main, pushed)

Re-scoped mid-flight by the user's **licensing framing**: *"to demonstrate
correctness we must be able to demonstrate the ledger of tests, which tests
exercise which code, and why we can trust them by dependency."* So step 1's
deliverable became the **ledger + the claim-falsification audit**, not a
better `retest`. Plan: `~/.claude/plans/peppy-conjuring-thimble.md`.

| commit | outcome |
|---|---|
| `f0a6e4c` | `ExecutionLedger` — a capture joined to the graph. Class DESCENT: `[M]` 0 of 438 production classes bind (coverage attributes lines), so 269 classes gained evidence they were structurally unable to carry |
| `1ca01f8` | `verification_coverage/audit(run=)` — 4 verdicts + `CaptureScope`; `orphan_code` stopped meaning "the call graph could not see it" (**1590** nodes moved to `tested`) |
| `2104ce4` | `retest(run=)` — evidence where a capture can speak; rows are runnable selectors (`[M]` 225 071 → 78 518 chars) and `limit` finally exists |
| `a0f038e` | the `run` arg reached the MCP tools (it had reached NEITHER), and the tally stopped counting evidence as its own corroboration |

⛔⛔ **The diagnosis that reframed everything — both instruments are wrong.**
Static `("calls","type_uses","inherits")` has **12–15 % recall** against
execution, and `[M]` **0 of 300** proven test↔symbol pairs have ANY path over
it (property 27 %, upstream break 25 %, dunder 21 %, no-source-caller 17 %,
phantom 9 %). It also over-claims 944 in-capture pairs — and ⛔ **that is NOT
a builtin-hub bug**: excluding non-project traversal leaves **250 of 250**
still reachable project-internally. My "ban `int`/`float` hops" hypothesis was
REFUTED by its own counterfactual. Coverage evidence is execution fact; its
limit is **scope**, not accuracy (def-line artifact refuted: 915 of 931 nodes
have `lines_hit ≥ 2`, within-module Jaccard 0.247 not ≈1.0).

⛔ **CORRECTION to this campaign's own earlier number.** An interim probe said
"**2720** of 2748 claims OUT-OF-CAPTURE". That conflated two causes with
opposite repairs: `[M]` **1751 out_of_capture** (needs a wider capture) +
**976 no_implementation** (needs a declared `implements` link — no capture can
fix it). 94 distinct equations, concentrated in `monte_carlo` (22),
`collision_probability` (15), `infinite_medium` (12), `method_of_characteristics` (11).

⛔ **And the finding that gates the licensing argument: 0 refutations are
trustworthy today.** `[M]` **all 10** refuted rows carry
`code_evidence=inferred` — the `implements` link was name-matched, so the
refutation lands on the guess, not the test. **nexus#82 is a prerequisite**,
alongside step 2's capture. Filed as ORPHEUS **#381**; #334 commented (its own
witness is 0/21 in capture).

⭐⭐ **Rulings that transfer.**
1. **A partial capture must not certify a test it never RAN.** I had "the
   capture covers the symbol ⟹ trust it" and it reported `safe_to_skip = 5161`
   for a geometry change on a 1499-of-5278 capture — 3779 tests blessed for
   never having been looked at. A capture speaks for the tests it *collected*,
   nothing more.
2. **A row minted from evidence cannot corroborate anything** — its verdict is
   a tautology. Counting them took `claims_corroborated` 5994 → 36466, a
   number tracking capture SIZE rather than suite quality (`plan-authoring`
   §10, in a metric I had just introduced). Now `executed_unclaimed`: 30 472.
3. **Gate an arm where it is FALSIFIABLE, not where it lives.** Two arms were
   blind through their feature and load-bearing one level down; two others
   were genuinely no-ops and the docstrings now say so rather than implying a
   gate exists.
4. ⚠ **Six for six: the MCP layer is where this campaign's defects hide.**
   `run` reached `GraphQuery` and neither tool; `assemble_*` payloads are
   hand-built, so a new field does NOT arrive for free.

▶ **Still open from step 1**: the nexus fixture has no contexts-carrying
coverage run, so `exercised_by` has no end-to-end witness in nexus's own suite
(every gate hand-builds a graph). And `verification_audit(run=)` returns
149 k chars — it has no `limit`.

### ▶▶ THE AGREED SEQUENCE — user ruling, 2026-08-18

Titled as outcomes (`plan-authoring` §1). Do them **in this order**; the
reasoning is that (3) must not land before (1) exists, or it repeats `cap`.

1. **`retest` answers from EVIDENCE, not from popularity.** Point the cone at
   coverage attribution (`exercised_by`, landed `ca6ccb0`) instead of — or
   ahead of — static `calls`. Needs no claims, no reorganisation, no new
   capture, and degrades to static where no capture exists.
2. **A capture wide enough to matter exists.** Decide the strategy; (1) tells
   you which directories matter first.
3. **A test declares what it is about** — the redesign, `.claude/plans/test_architecture_redesign.md`.

⛔ **Why (1) is first, measured 2026-08-18 on the deduped graph.** `retest`
today walks `("calls", "type_uses", "inherits")`, and `calls` reachability
measures **popularity, not dependence**:

| root symbol | cone reaches |
|---|---:|
| `Quadrature.gauss_legendre` (a leaf utility) | 1307 tests, **24.8 %** |
| `solve_sn` (the central solver) | 117 tests, **2.2 %** |
| `Mesh1D.__post_init__` | **0** |

A leaf helper pulls **12×** more of the suite than the solver. This confirms
nexus#60 (OPEN) survives the dedupe — duplicate edges share `(u,v)`, so
reachability was never affected; now measured rather than reasoned.

⛔ **A correction to CP#3, preserved from the interim block this point replaces.** CP#3's ▶ NEXT said a cone "walked over `imports` today over-invalidates". `[M]` `_RETEST_DEPENDENCE_EDGES = ("calls",
"type_uses", "inherits")` — `retest` does not follow `imports` at all. The
stamp is a prerequisite for the module-precise cone #358 still has to
build, not a repair to one that ships. And the issue's headline "14 of 365"
counts CROSS-LAYER edges; the cone's exposure is **199 of 4045**
intra-project. Both right, different predicates.

### Landed since CP#3

| what | outcome | re-measure with |
|---|---|---|
| **nexus#88 CLOSED** `9ac3d4b` | an `imports` edge minted under `if TYPE_CHECKING:` says so; declared on `[edge.imports]`; both spellings | `stats()` |
| nexus **`924efbf`** | ⛔ an `extra_source_dirs` entry INSIDE a scanned root was analysed TWICE — **207 643 → 139 761 edges (−33 %)**, `calls` 118 126 → 69 217. `god_nodes` RANKING changed (SNMesh now #1) | `stats()`, `god_nodes()` |
| nexus **`6e469e9`** | the stamp reaches the REPLY; a runtime + type-only pair of one import stops folding into a false `times: 2` | `neighbors(<mod>, edge_types="imports")` |
| ORPHEUS `1cb0c76f` | the `.nexus/config.toml` comment claimed a causation it did not have | — |
| ORPHEUS `95199f0e` | `pyproject.toml`'s `cap`/`sentinel` marker help cited 2 archived plans — 10 dead paths | — |
| ORPHEUS `a17243f0` / `54cd765a` / `3f3cc9ca` | the redesign proposal + two experiments + the reframe | the plan |
| filed | nexus **#89** (`extra_source_dirs` resolves against `srcdir.parent`; its warning hides in `-q`) | — |
| new run | **`num_ctx`** — `tests/numerics` coverage w/ contexts, 2344 tests | `runtime_runs()` |

⭐ Do NOT copy figures out of those tools into here (`plan-authoring` §9).
The fidelity probes RE-MEASURE — run `evals/fidelity_probes.py --project <root>`,
never quote a number from a doc.

### `[M]` Numbers a pick-up needs that nothing else re-measures

- **Capture cost**: `tests/numerics`, 2344 tests → **567 s** → **467 MB**
  coverage report → **3.37 MB** sidecar (**139×**), 1303 attributed nodes, 0
  unresolved contexts. ⟹ the STORED artifact is a non-problem; the
  **intermediate** is the blocker for (2).
- **Subject derivability** (the experiment that reshaped the redesign): "narrowest
  executed production node" recovers the module-level subject **78.2 %**
  (`tests/numerics`, 531 tests / 25 files) and **73.6 %** (clean geometry);
  ceiling 94.2 %; truth ranks #1 in 364/500, p90 rank 3.
- **Granularity instability**: at SYMBOL granularity **17 %** of derived
  subjects are a private helper, vs **1 %** at module granularity.
- **Family structure**: `protocol_conformers` reports **39 classes** on one
  `apply`/`apply_transpose` protocol across 7 packages. An `inherits`-based
  family detector finds this in **1 of 20** scattered files — conformance is
  STRUCTURAL.

### ⭐⭐ Rulings that transfer

1. ⭐⭐ **A stamp nothing can READ is stored, not shipped — and the reply can be
   worse than silent.** #88 put `type_checking` in the DB; no tool exposed it,
   and `neighbors` folded a runtime + a type-only import of one target into
   `times: 2`, asserting a sameness that does not hold. **Fifth** reply-shape
   defect of this campaign invisible in-process. ⟹ for every new fact, name the
   tool that returns it, and call that tool.
2. ⭐⭐ **Adding an attribute to the graph buys a FREE consistency check against
   the source — take it, and chase the residue.** The stamp let an independent
   AST census be compared against the graph; `orpheus/` and `examples/`
   reconciled exactly and `tests/` was **exactly 2×**, which is the only reason
   a 33 %-duplicate graph was ever noticed. Nothing else reported it — a
   doubled graph looks like a big graph.
3. ⭐⭐ **A feature whose own fixture silently no-ops has no witness.** nexus's
   fixture declared two `extra_source_dirs` that resolved to nothing and warned
   into a `-q` build; that is why the doubling survived. ⟹ when a gate is
   suspiciously green, check the fixture REACHED the feature.
4. ⚠ **"No duplicated pair" is the wrong shape for a duplicate-scan gate** —
   two `from X import a` / `from X import b` statements legitimately mint two
   edges, so the gate reddened on a correct tree one commit after I wrote it.
   Count against an INDEPENDENT census of the source, and mirror the analyzer's
   own traversal (function bodies are never visited statement-wise, so
   `ast.walk` over-counts).
5. ⚠ **"Not a runtime attribute" is not "not a resolvable name" — ask what the
   map is FOR.** #88's re-export half was reverted on exactly this.
6. ⚠ **MINE, §2 twice**: a blast radius priced from an `__init__.py`-only probe
   (15) when the code path was tree-wide (**263 of 6917**); and a revert
   justified by a +5 delta that **survives the revert**. Both are logged in
   `plan-authoring`.
7. ⚠ **A mutation battery must be CRASH-safe, not exception-safe** — a harness
   SIGTERM skips `finally` and leaves production code mutated on disk.
   Promoted to `.claude/rules/process-discipline.md`.
8. ⭐⭐ **Do not use today's names as ground truth when the names are what you
   intend to change** — that scores one heuristic against another and counts
   every disagreement against the wrong half. Ask what the tree SHOULD be, and
   characterise rather than score.

---

## ⏸ COMPACTION POINT #3 — 2026-08-18 · a claim is falsifiable now; the cone is not yet honest

⚠ **Everything below this section is HISTORY unless a hash says otherwise.**
Both repos are on `main` and **pushed**, 0 dirty / 0 unpushed:
nexus **`4252989`**, ORPHEUS **`530eab29`**.
`git merge-base --is-ancestor <hash> main` is the authority.

⚠ **`/mcp` needs reconnecting AGAIN.** It was reconnected mid-session and
served `ca6ccb0`; two fixes landed after (`1ab2ebe`, `4252989`) and the
server holds the code it imported at startup.

### Landed this session

| what | outcome | re-measure with |
|---|---|---|
| **nexus#57 CLOSED** `ca6ccb0` | `RuntimeRun.exercised_by` + `runtime_exercisers` (query/CLI/MCP). REACH is evidence now, not inference | `nexus runtime-exercisers --db .nexus/graph.db --run geom_ctx` |
| `1ab2ebe` | the exerciser list is **runnable pytest selectors**, not NodeResults | — |
| `4252989` | an `.. error-entry::` node carries its declared line | `errors()` |
| **ORPHEUS `538b5f5c`** | ERR-026's archaeology block: 29 roles → 13, **15 dead → 0** | the probe in §"roles" below |
| **ORPHEUS `c0d79933`** | the xref gate could not import `tests`; **3 failed → 46 passed**, 0 dead tree-wide | `python -O -m pytest tests/test_docstring_xrefs.py` |
| filed | ORPHEUS **#380** (30 stale raw PATHS — no gate checks a path) | — |

⭐ Do NOT copy figures out of those tools into here (`plan-authoring` §9).

### ⭐⭐ Rulings that transfer

1. ⭐⭐ **Scoring and dependence are different facts.** `# pragma: no cover`
   removes a line from coverage's numerator AND denominator while coverage
   goes on stamping CONTEXTS on it — it ran. `[M]` gating attribution on the
   coverage guard drops 4 ORPHEUS nodes, all pragma'd guards;
   `DiscreteMeasure.__post_init__` is executed by **131** tests and would have
   reported none.
2. ⭐⭐ **Absence of a run is not absence of exercise** — hit on the feature's
   FIRST use, by its author. A geometry-only capture made **53 of 53**
   equations look like every claiming test executed nothing of their
   implementation. `[M]` **0** of `alpha-recursion`'s 39 claimants were in the
   capture; they are SN tests. The whole signal was the slice. ⟹ intersect the
   claimants with the tests the run CONTAINS before reading a zero as a
   refutation.
3. ⭐⭐ **A "hidden findings" number must be checked for FALSE positives before
   it is believed** — the §10 lesson in a new costume. A well-argued one-line
   repair to the xref gate was measured to reveal **158 hidden dead targets**;
   `[M]` all 155 in `docs/` were `tests.*` and **alive**. The diagnosis was
   right and the remedy was wrong: the real defect was `sys.path`, and with it
   corrected the same patch reveals **nothing** (5 tree-wide with and without).
   ⟹ sample the revealed items by hand before adopting the instrument.
4. **An unimportable target reads as a missing one.** The gate ran as a script,
   so `sys.path[0]` was `tools/` — `orpheus` resolved (pip-installed editable)
   and `tests` did not. **49 of 49** dead targets were `tests.*` and existed.
   A gate that cries wolf gets ignored, which costs more than no gate.
5. **A reply shape must be checked at the SERVER.** `runtime_exercisers` blew
   the 20 000-char budget on its first real MCP call (38 of 130 kept) while the
   in-process gate stayed green — it asserted a test NAME was present, which
   remained true as the list was truncated around it.

### ▶ NEXT — nexus#88, and it is a CORRECTNESS prerequisite, not polish

A `TYPE_CHECKING` import creates **no runtime dependence**, and the `imports`
edge does not say so — so a cone walked over `imports` today
**over-invalidates**, which is the exact failure the DAG exists to remove.
`[M]` 14 of 365 cross-layer edges, concentrated on the order-inverting L2→L3 /
L1→L3 pairs that pull the largest downstream sets. ORPHEUS's own
`tests/test_layer_imports.py` explicitly TOLERATES these, so the graph is
currently **stricter than the contract and cannot express why**.

The issue proposes stamping `type_checking: True` on the edge rather than
dropping it — "what does this module reference for typing" is how `#76`'s
attribute-mediated dispatch would be recovered.

**Then ORPHEUS#358.** Its three tiers: (a) coarse layer DAG — ruled out by the
claim-KIND ruling, (b) module-precise — `0d6bfdf`, (c) symbol-precise per-test
— `ca6ccb0`. What #358 still needs beyond #88 is a **capture wide enough to
matter**: `[M]` one directory produced a 265 MB report (1.44 MB stored), so a
whole-suite capture is not a scaled-up geometry capture.

**Also open:** ORPHEUS **#302** (its gate half is REMEDIED; its SUBJECT —
roles that resolve but have no autodoc target, hence no link — is untouched),
**#301**, **#379**, **#380**; nexus **#85/#86/#87**, **#76/#16**, **#55**.
Owed on the fidelity side: **F1-honesty**, **F4-recall**, and a real
multi-agent field trial — round 6 was probes plus targeted stress only.

---

## ⏸ COMPACTION POINT #2 — 2026-08-17 · the import DAG is visible; REACH is next

⚠ **Everything below is HISTORY unless a hash says otherwise.** Both repos are
on `main` and **pushed**: nexus `c68fe4d`, ORPHEUS `5ae4cd6d`, both 0 dirty /
0 unpushed. `git merge-base --is-ancestor <hash> main` is the authority.

⚠ **Reconnect `/mcp`** — the server serves the code it imported at startup, and
`errors` is new. ⚠ **The ORPHEUS graph was REBUILT** after nexus `0d6bfdf`; a
graph built before it still has the collapsed imports.

### Landed since the 2026-08-15 point

| what | outcome | re-measure with |
|---|---|---|
| **nexus#63 CLOSED** | `errors` query + CLI + MCP tool; `catches` 0 → **258 edges** | `nexus errors --db … --format text` |
| **ORPHEUS#308 first instance** | the 79-entry catalogue is 79 corpus nodes; skill index 374 KB → 9 KB | `tools/verification/generate_error_index.py --check` |
| ⛔ **`.. error-entry::` crashed every real build** | `apply_pending_edges` read `entry["label"]` on a declaration payload | `tests/test_directives.py` |
| **F9 — the granularity collapse** | `imports` kept only `split(".")[0]`; god node **5299 → 1**, cross-layer **1 → 365** | `evals/FIDELITY.md` round 6 |
| filed | nexus **#88** (TYPE_CHECKING edges), ORPHEUS **#379**; #301 has the derivability criterion | — |

⭐ Do NOT copy figures out of those tools into here (`plan-authoring` §9).

### ⭐⭐ Rulings that transfer

1. ⭐⭐ **The layer order is an ordering over CLAIM KINDS, not a scheduler.**
   *You cannot validate before you verify; you cannot verify numerics whose
   primitives are untested.* Correct for V&V, too coarse for test selection.
   It gives DIRECTION, never REACH. See the ⛔ in §2.
2. ⭐⭐ **The DAG the work wants is per-symbol.** Change `sn_solver` → the base
   node of *its* test cluster flips untested, and that cluster is itself a
   mini-DAG **ordered by increasing complexity**, re-executed in that order.
   The point is fail-fast: the cost being paid is *waiting for tests to clear*,
   so the cheapest discriminating test must red first.
3. **A collapse that preserves the EDGE COUNT is invisible to every aggregate**
   — and degree centrality *promotes* it to the top of the hub list, where it
   reads as architecture. Two sub-agents concluded "no antisymmetric relation
   exists" from a corpus that machine-enforces one.
4. **A metric adopted to CHOOSE between implementations needs the same check as
   one adopted as a target** — `plan-authoring` §10. I picked a regex by warning
   count (21 vs 27); a census found **830** mangled delimiters, ~2 % of which
   ever warned.

### ✅ LANDED 2026-08-18 — REACH is measurable (nexus #57, `ca6ccb0`)

`RuntimeRun.exercised_by` + `runtime_exercisers` (query / CLI / MCP). The
graph can now answer *"which tests rest on this symbol"* from EVIDENCE —
the half ruling 1 says the layer order can never supply.

`[M]` ORPHEUS `tests/geometry`: **426 contexts, 426 resolved, 0 unknown;
931 code nodes, 11 182 (code, test) pairs**. Stored as run `geom_ctx`.
Re-measure with `nexus runtime-exercisers --db .nexus/graph.db --run
geom_ctx` — do NOT copy figures out of it (`plan-authoring` §9).

⚠ **The capture, not the read, was the work.** The compaction point's own
`[M]` was right: `.nexus/traces/*.json` are sidecars with the contexts
already discarded, so this needed a re-capture (`dynamic_context =
test_function` + `coverage json --show-contexts`, config-only).

⭐⭐ **Rulings this produced, both transferable:**

1. **Scoring and dependence are different facts.** `# pragma: no cover`
   removes a line from coverage's numerator AND denominator while
   coverage goes on stamping contexts on it — it ran. `[M]` gating
   attribution on the coverage guard drops 4 nodes, every one a pragma'd
   guard; `DiscreteMeasure.__post_init__` is run by **131** tests and
   would have reported none.
2. ⛔ **Absence of a run is not absence of exercise** — hit on the
   feature's FIRST use, by its author. A geometry-only capture made
   **53 of 53** equations look like every claiming test executed nothing
   of their implementation. `[M]` **0** of `alpha-recursion`'s 39
   claimants were in the capture: they are SN tests. The whole signal was
   the slice. ⟹ before reading a zero as a refuted claim, intersect the
   claimants with the tests the run actually contains.

⚠ **A whole-suite capture is not a scaled-up geometry capture.** `[M]`
one directory produced a **265 MB** report (reducing to 1.44 MB stored).
Slice, or budget accordingly.

**Then, and only then, ORPHEUS#358.** Its premise is corrected in §2; the three
tiers are (a) coarse layer DAG — available today, **and ruled out by ruling 1**,
(b) module-precise — LANDED at `0d6bfdf`, (c) symbol-precise per-test — LANDED
at `ca6ccb0`. ⟹ **#358 is no longer blocked on evidence; it is blocked on
`nexus#88` and on a capture wide enough to matter.**

⚠ **nexus#88 is a correctness prerequisite, not a nicety.** A `TYPE_CHECKING`
import creates NO runtime dependence, and the edge does not say so, so a cone
walked over `imports` today **over-invalidates** — the exact failure the DAG
exists to remove. `[M]` 14 of 365 cross-layer edges, concentrated on the
order-inverting L2→L3 / L1→L3 pairs that pull the largest downstream sets.

**Also open:** ~~ORPHEUS **#302** (`test_docstring_xrefs` is RED on `main`, 71
dead sites, none mine)~~ ✅ **REMEDIED 2026-08-18 @ `c0d79933`** — the gate is
GREEN (46 passed, 0 dead tree-wide). ⛔ The 71 were **not dead references**: the
tool ran as a script, so `sys.path[0]` was `tools/` and it could not import
`tests` — 49 of 49 dead targets in `docs/` were `tests.*` and alive. #302's own
SUBJECT (roles that resolve but have no autodoc target, hence no link) is
untouched and still open. **#301** (algebra-of-record derivability), **#379**,
**#380** (30 stale raw PATHS in the catalogue — no gate checks a path); nexus **#85/#86/#87**, **#76/#16**, **#55**. Owed on the fidelity
side: **F1-honesty**, **F4-recall**, and a real multi-agent field trial —
round 6 was probes plus targeted stress only.

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

⚠ **PARTLY REFUTED 2026-08-17 — the conclusion about `calls` STANDS; the
universal about the corpus does NOT.** Both agents searched for an
antisymmetric relation, found `calls` symmetric, and concluded none existed.
The reasoning was sound and the corpus was lying: ORPHEUS **declares and
machine-enforces** a layered import contract (`tests/test_layer_imports.py`
— L0 `derivations` → L1 `numerics` → input → L2 `transport` → L3 `sn`/…),
and nexus had flattened it to a star before either agent could see it.
`[M]` **5298 of 5299** project `imports` edges pointed at bare
`py:module:orpheus` — the graph's #1 `god_nodes` hit, which was the collapse
rather than a hub. Fixed at nexus `0d6bfdf`: degree **5299 → 1**, cross-layer
edges **1 → 365**, of which 351 honour the contract and all 14 exceptions are
`TYPE_CHECKING` (which the contract tolerates). Recorded as fidelity class
**F9** (`evals/FIDELITY.md` round 6).

⛔ **But do NOT build the DAG on the layer order** (user ruling, 2026-08-17).
It is an ordering over **claim kinds** — *you cannot validate before you
verify, and you cannot verify numerics whose primitives are untested* — which
is epistemic and correct for V&V, and **too coarse for scheduling**. It gives
DIRECTION, never REACH: it says an L3 test cannot invalidate an L1 test; it
cannot say *which* L3 tests a given L1 red invalidates. Reach is still §2's
open problem, and is `nexus#57`.

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
