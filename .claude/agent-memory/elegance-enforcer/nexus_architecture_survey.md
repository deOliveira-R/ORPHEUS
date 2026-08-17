---
name: nexus-architecture-survey
description: Cardinal-Rule-2 survey of sphinxcontrib-nexus (2026-08-16) — the verdict, the three missing types, what is well-factored and must not be touched, and the rejected splits.
metadata:
  type: project
---

Whole-package architecture review of `/Users/rodrigo/git/sphinxcontrib-nexus`
(branch `feat/config-and-ontology`, HEAD `5796c6d`, 15,682 ln / 22 modules).
Report: `/Users/rodrigo/git/nuclear/ORPHEUS/scratch/nexus_architecture_review.md`.
Trigger was "where does *test state* (green/red/untested) live?" — the user was
about to add it as a fifth family on `RuntimeRun`.

**Why:** the user's standing instruction — *"if Nexus is loose and not well
structured, tackle it while it's still relatively small."* This is the baseline
survey a later restructure is planned against.

**How to apply:** re-verify against git before acting (nexus moves fast). Reuse
the *shape* of the verdict, not the line numbers.

## Verdict — reshape by OBJECT, not by file

The pipeline (extract → merge → store → query → serve) is real and the file
structure respects it, so a *file* split is still wrong (adds ~7 concepts,
removes none). But `[M]` `GraphQuery` has **two** states — `self._g` AND a
working tree passed as a kwarg to **10 of 56 methods** (2 shell out to `git`) —
i.e. a procedural module wearing a class, which is exactly the *"consolidate data
and functions that share a context into objects"* case. Five objects are latent
and already half-built across the package: **`NodeId`**, **`ProjectView`**
(graph+tree+provenance), **`PositionIndex`**, **`Evidence`**, **`Diagnostic`**.

⭐⭐ **The load-bearing move: B2 (only 1 of 40 MCP tools checks graph staleness)
is not an oversight — it is the mechanical CONSEQUENCE of the missing
`ProjectView`.** A check needing two states cannot be enforced when one is an
optional kwarg, so it becomes something 40 authors must remember, and 39 did not.
**Generalise this probe: when a guard is applied at 1 of N sites, ask what object
would make it unforgettable** — that object is usually the finding.

⚠ **My own two retractions, same error both times: a balanced-phase verdict was
allowed to BOUND a later adversarial pass.** (a) I called `GraphQuery` "a
namespace, one shared field" — incomplete measurement, and I had B2 and its
explanation in the same review without connecting them because "the layering is
sound" had closed the question. (b) I graded the missing `NodeId` "compensation
machinery, not wrong answers" — true on nexus's own 3,309-node graph, **false on
ORPHEUS**: `[M]` 171 Sphinx-role shadow nodes (`py:meth:`/`py:func:`/`std:doc:`),
221 edges, **126 provable duplicates**, `std:doc:` 100 %, and one node id
containing a **literal newline**. Cause: the role→objtype map
(`ast_analyzer.py:920`) is applied on the AST minting path and **not** on the
Sphinx-xref path (`extractors.py:402`). ⟹ **choose the denominator that can show
the defect** — the small self-graph hid it.

## ⭐ METHOD LESSON (user, via coordinator) — the adversarial phase comes FIRST

*"Code review for elegance needs adversarial intent and a clear 'how would I make
this 100x better / how would I break this' frame … Then later you can review your
own criticism and re-evaluate … But there must be a clear adversarial phase."*

⟹ **Phase 1 ADVERSARIAL** (no balance): *how would I break this* — silent
wrongness, weighted to defects that fail in the REASSURING direction — and *how
would I make this 100× better* (reframe the JOB, not the spelling). **Phase 2
RE-EVALUATION**: withdraw what does not survive, saying why you expected a
defect. "Well-factored / do not touch" belongs in Phase 2 **as a withdrawn
attack** — it carries far more information than a compliment written up front.

**Measured yield difference on this very review.** The balanced first pass asked
"what is misplaced?" and produced 12 findings, 2 of them live defects. The
adversarial pass, run afterwards on the same loaded context, produced the
flagship (B1 below) plus 3 more break-it findings **and inverted one Phase-1
verdict**: `workspace.py`, rated "coherent, do not touch" by the survey, is
actually *correct and under-wired* — its provenance stamp has `[M]` 2 readers
against 40 MCP tools. A balanced pass rates a module by internal quality and
cannot see that its product goes nowhere.

## The flagship break-it finding (B1) — the pattern to look for elsewhere

**An inference layer whose output is presented identically to a fact.**
`merge._infer_implements` (`merge.py:188`) links code→equation on **one shared
name token ≥3 chars**. `[M]` on ORPHEUS: **16,068 of 16,068 `implements` edges
are `source="inferred"`, confidence 0.7** — zero declared; median 12 / max 99
claimed implementers per equation; `Billiard` "implements" a cylinder Green's
trajectory equation on `shared=["trajectory"]` taken from its *module path*.
`verification_coverage` consumes them **with no confidence filter**
(`query.py:1528-1533`) while reading confidence for `tests` edges 8 lines later
(`query.py:1544`), then BFS-3-hops for an `is_test` ancestor ⟹ `[M]` **345 of 893
ORPHEUS equations (38.6%) read `verified` with ZERO declared evidence**.
`TestReference` carries `source`+`confidence`; `CoverageEntry.implementing_code`
is bare `NodeResult`, and **`EdgeResult` has no confidence field at all** — so no
query surface can express it.

⟹ **Generalisable review probe:** when a codebase stores `confidence`/`source` on
its data, grep whether the RESULT type carries them. The gap between "the model
knows" and "the answer says" is where unmarked guesses escape.

## The 100× reframe (for reuse)

Nexus today = a graph API with `[M]` 36 frozen analyses × 40 MCP tools × 44 CLI
verbs × 48 result dataclasses. **The 100× is: an evidence engine** — (1) every
answer is a claim with provenance/confidence/graph-commit; (2) one query language
+ a versioned query library replaces the 36 frozen questions (concept count goes
DOWN — the Pattern-6 test the module split failed); (3) stamped graphs gain
memory so nexus answers about CHANGE, which is what an agent actually asks.
⭐ The tell that the ARTEFACT is wrong (not the fields): every Part-A fix was
"add a field to one more bespoke dataclass", 36 times, and the 37th forgets.

## The two live defects (measured, not hypothetical)

1. **Twin AST dotted-name reconstructors** — `ast_analyzer._dotted_name` returns
   `None` on a non-`Name` chain root; `_unparse_attribute` silently returns the
   partial leaf, and it is the one wired into `_resolve_call_target`. `[M]` 8 of
   13 incoming `calls` edges on `project.resolve` are fabricated, every one from
   a `pathlib.Path.resolve()`. False-ALIVE direction: inflates `impact`, and
   *rescues* dead functions from `dead_functions`.
2. **Three implementations of "(file,line) → node"** — `query.node_at`,
   `runtime.build_node_index`+`resolve_node`, `brief._in_file_node_ids`. `[M]` 3
   of 4 probed positions disagree; `DECORATOR_WINDOW = 8` lets a `def` BELOW the
   position claim it (L707 class-docstring → `__init__` at L711). This is the
   shared mechanism the test-state overlay needs ⟹ **prerequisite, not cleanup**.

## The branch's one must-not-ship

`ontology.py` (292 ln) landed with **zero production consumers** and a docstring
claiming the hardcoded `merge.py:248 in_test_file` rule was dissolved into a
declared constraint. `[M]` both spellings ship; the concept count went UP
(Pattern 6). Wire `check_edge` into `_finalize_graph` AND retire `merge.py:248`
in the same commit, or hold the module.

## What is well-factored — the standards to hold new work to

- `_serialize.py` — the shared assembly layer; `[M]` 7/7 used by server, 6/7 by CLI.
- `brief._in_file_node_ids` — **the model for a NECESSARY twin**: structural
  reason + reciprocal comments in both files + a real both-realizations pin
  (`test_brief.py:140`). Cite this when grading any other duplicate in this repo.
- `EffectiveSettings`/`_effective` (Pattern 7, argued in its own docstring),
  `resolve_db` (one answer for CLI+server), `_mappings.candidate_rank`
  (6-tuple ranking key with a per-level rationale + `_IDENTITY_RANK_WIDTH`),
  `workspace.py`, and the enum↔`ontology.toml` join pinned at `test_ontology.py:36`.

## Rejected splits — do not re-derive these

- **Split `query.py` / `GraphQuery` by section**: fails the concept-count test.
  `[M]` 56 methods, only shared state is `self._g`, 33 % of intra-class calls go
  to one 15-line adapter, 28/56 methods call no sibling — it is a namespace, and
  splitting adds ~7 concepts and removes none.
- **`ast_analyzer.py` is a god module**: it is not — 4 coherent layers. Its only
  misplaced content is the ~500 ln of graph-*repair* post-passes, whose real
  problem is the triple `external`-vs-`unresolved` predicate, not the file.
- `callers`/`callees` twin (0.90): symmetric-by-design; unifying needs a
  `direction=` flag = anti-#3.
- `dead_functions`' 39 self-FPs: all the `cli.py:934` dispatch dict (a registry)
  + 2 Sphinx hooks. Tool limitation, documented, not a defect.
- `RuntimeStore` JSON vs `export.py` SQLite: **two lifetimes ⟹ two stores**
  (the graph is rebuilt every sphinx-build; sidecars must survive). Do not unify.

## The reusable measurements

`[M]` bare vocabulary strings vs enum refs: **340 : 187** package-wide;
`query.py` alone **140 : 20**. `NodeType(str, Enum)` gives ZERO static
protection — `== "function"` and `== NodeType.FUNCTION` type-check identically.
The repo already `Literal`-types `direction` (2–3 values) and leaves
`node_types: list[str]` (32 values) untyped — that contrast is the sharpest way
to make the point. `[M]` the vocabulary is enumerated 4× exhaustively and
`visualize.py` has already drifted (15/16 nodes, 13/16 edges) while 24 `tag`
nodes and 42 `discriminates_on`/`tests` edges render in the fallback colour.
See [[nexus-runtime-overlay]], [[nexus-workspace-resolution]].
