# Explorer — Lessons

Behavioral corrections only: "what mistake did I make exploring, and what
did I learn that improved my behaviour?" The HOW of each Nexus tool lives
in the preloaded skills (`nexus-exploring`, `nexus-guide`, and the wider
`nexus-impact` / `nexus-debugging` / `nexus-verification` / `nexus-refactoring`
family) — point there, never duplicate the workflow here. Per-campaign
`file:line` maps are archaeology; they go stale in days. A lesson earns its
place only if it changes how the NEXT exploration is run.

The cross-cutting spine: **an exploration answer is not "I found the symbol"
— it is "I found EVERY consumer the next action will touch, I verified the
premise against the current tree (not the issue text, not a frozen memory),
and I separated the durable subsystem-shape from the line numbers that will
drift."** This spine is now codified as standing directives in `AGENT.md`
Operating Principles 4–7 (blast-radius, premise-verification, git-merge-status,
durable-vs-line). L-001/L-002/L-003/L-005 below are RETAINED for forensic
value (the war-story behind each directive); the directive itself, not the
incident, governs behaviour. L-004 and L-006 stay lesson-only — they fire on
narrower question shapes (carve-verdict / probe-collapse), not every task.

---

## L-001 -- A retirement/rename blast radius = graph callers AND text grep AND direct constructors AND doc nodes

→ **Now AGENT.md Operating Principle 4** (the four-search discipline). War-story kept below for the specific misses each search catches.

`mcp__nexus__callers` / `impact` find the *graph* consumers, but a retirement
audit that stops there under-scopes — and under-scoping a retirement forces a
mid-session re-plan (the documented ~2× cost behind the proactive-explorer
trigger). The consumers a single `callers` query misses:

- **Property-reached leaves.** A method only reachable by reading a
  `cached_property` and calling `.apply` shows `callers() == 0` while still
  being live through the property. Audit the property's readers, not just the
  method's callers.
- **Bypass-trick / class-name consumers.** A test that uses an orphan via its
  CLASS NAME for a side purpose (a validation-bypass) is invisible to a
  method-level `callers`. A repo-wide grep for the `_ClassName` surfaces it.
- **Direct constructors of a guarded type.** A guard-at-the-data-source has
  blast radius = EVERY direct `Foo(...)` call, not just the factory path.
- **Doc nodes that will dangle.** A retired symbol referenced from a theory
  page leaves a broken `:ref:`. `graph_query` for the doc→symbol edge (or the
  symbol name in `docs/`) catches it so the archivist hand-off is complete.

How to apply: for ANY retirement/rename audit, run BOTH graph (`callers`/
`impact`) AND a text grep of the symbol/class name AND a constructor audit (if
a guarded type) AND a doc-node scan. Four searches, not one. (Reinforces the
proactive-explorer-before-retirement trigger; sibling to method-implementer
L-004.)

**Sharpening (two graph-blinding patterns endemic to the operator algebra).**
Two constructs make `callers()` systematically lie in this codebase, and BOTH
appeared in the W-F scope audit — when an `apply`-dispatch arm's liveness is the
question, grep is not a cross-check, it is the *primary* evidence:
(a) **runtime-aliased dispatch** — `apply = _apply_impl` with `@singledispatchmethod`
arms means the graph attributes calls to the alias, not to the per-type arm; a
`@_apply_impl.register` leaf shows near-zero callers while a whole solver feeds it.
(b) **Protocol-typed receivers** — a call `solver.compute_fission_source(...)` where
`solver: EigenvalueSolver` (a Protocol) is unresolvable to the concrete `SNSolver`
method, so the concrete method reads `callers==0` even though `power_iteration`
drives it every outer step. The liveness of a dispatch arm is decided by the
ACTUAL input TYPE at the production call site — trace `power_iteration →
solver.compute_X → op.apply(<what type?>)` by READING the chain, and let grep, not
the graph, enumerate the `op.apply` sites.

---

## L-002 -- The issue text is a stale premise; verify it against the current tree FIRST

→ **Now AGENT.md Operating Principle 5.** War-story kept below for the worked examples (diamond-coefficients / 2-D-matvec premises already landed).

Repeatedly, an audit's first deliverable was "the premise the issue describes
is STALE — that work already landed." An issue body is written at one moment
and the natural trigger for its work (a related carve) often lands it early
under a *different* campaign. Examples of the same shape: a "lift the inline
diamond coefficients" issue was already folded onto the scheme; a "2-D matvec
recomputes inline" concern was resolved by an earlier phase.

How to apply: before mapping HOW to do an issue, spend one query confirming it
still NEEDS doing. Grep the named symbol / read the current body of the named
function. If the premise is stale, the deliverable flips from "implement" to
"CLOSE-VERIFY (regression-pin + issue hygiene)" — say so up front. This pairs
with the git-authoritative discipline (L-005): code state, not the issue's
prose, is ground truth.

---

## L-003 -- Separate the DURABLE subsystem-shape from the line numbers that will drift

→ **Now AGENT.md Operating Principle 7.** War-story kept below for the home-placement detail (durable → AGENT.md durable-shape section; transient → topic file flagged with the HEAD it was current at).

Every audit I wrote mixed two things with opposite shelf-lives: the durable
STRUCTURE (what couples to what, which seam is polymorphic, which path is
canonical) and the perishable `file:line` map. The structure survives years of
churn; the line map is wrong within a sprint. A memory that fronts the line map
reads as authoritative long after it has rotted, and a future session trusts a
dead address.

How to apply: lead every finding with the durable claim ("the within-group
operator is the variadic `(L+C, S, B)`; the sweep reads `ψ.boundary.inflow` and
does NOT re-apply `R·G` internally"), then mark line numbers as
re-derive-via-Nexus, never as the headline. The durable subsystem-shape belongs
in `AGENT.md` (its "SN operator-algebra subsystem — durable shape" section is
the canonical home); transient maps belong in a topic file flagged
"line numbers current at HEAD X, re-derive if drifted" — and are deletable once
the campaign merges.

---

## L-004 -- A clean carve verdict names BOTH the retire case and the keep-as-anchor case, with the discriminator

The strongest audit verdicts were not "retire it" or "keep it" but "RETIRE-eligible
BY <discriminator>, AND here is the defensible documented-KEEP, AND the call is
the user's because it turns on a judgment the explorer can surface but not make."
The discriminator that decides it is the `coding-standards` aggressive-retirement
rule's own test: same-math-available-via-the-surviving-helper ⟹ retire (genuine
redundancy); genuine-independent-consumer-need (even a future one a named typed
leaf would serve) ⟹ keep-as-anchor is defensible. The honest counter-weight is
Cardinal Rule 2 cutting both ways: a clean typed surface a future preconditioner/
DSA would consume is an architectural asset, not just dead weight.

How to apply: when asked "does X earn its keep?", deliver the dependency surface
+ the retire-with-rewire map + the keep-as-anchor counter-weight + the explicit
discriminator, and hand the value judgment to the user. Do not pre-decide a
retirement that turns on "will a future open issue consume this."

---

## L-005 -- Git is authoritative for merge-status; a memory's "in-flight / NOT pushed" freezes mid-flight

→ **Now AGENT.md Operating Principle 6** (and the always-on `process-discipline.md` rule). War-story kept below for the SN-campaign pattern that motivated it.

Memory notes captured a campaign as "uncommitted on branch X / NOT pushed," but
nearly every SN campaign merged in a later session — the note froze the moment it
was written. A future dispatch that trusts the frozen "in-flight" wastes effort
re-deriving landed work or, worse, treats merged code as still-pending.

How to apply: NEVER trust a memory's merge-status. Reconcile every "resume X /
in-flight X" against `git merge-base --is-ancestor <hash> HEAD` (or
`... <branch> main`) before acting. Active-state in MEMORY.md should say only
what git confirms; when in doubt, the answer is "check git." (Now an always-on
rule: `.claude/rules/process-discipline.md` §"Trust git for merge-status".)

---

## L-007 -- On a branch under ACTIVE edit, re-run the census immediately before reporting

During the F2 cast-family recon (2026-07-02, `refactor/pyright-burndown`), a cast
site moved 1532 → 1552 BETWEEN two of my greps: the main session was editing the
same files concurrently (uncommitted C3 carve in flight), and part of my brief
(the scattering `apply_transpose` item) was being FIXED while I explored. This is
intra-session drift — a different failure shape from L-002's stale-issue drift.

How to apply: when `git status` shows uncommitted edits in the subsystem being
audited, (1) re-run the position census as the LAST step before writing the
report, (2) diff the uncommitted hunks against the brief's items — an item may
already be mid-fix, flipping that deliverable to "confirm the in-flight fix +
report the alternative," and (3) timestamp reported line numbers as "at final
read; tree moving."

---

## L-008 -- zsh: an unquoted separator starting with `=` aborts the WHOLE compound command

`echo ===` (or any unquoted word starting with `=`) triggers zsh's `=cmd`
expansion; the lookup fails ("== not found") and the ENTIRE command line is
aborted — including greps sequenced after the echo, silently costing a
round-trip. Quote separators (`printf 'NAME\n'` or `echo "---"`), never bare
`===`, when batching multiple searches into one Bash call.

---

## L-009 -- A dataclass-FIELD rename audit is a grep problem, not a graph problem — and the field name may be a substring of an English word

Two independent findings from the `WithinGroupSystem.resolvent`/`.gains` rename
audit, both of which change how the NEXT field-rename audit is run:

**(a) Nexus does not model dataclass fields as nodes.** `context` on
`py:class:…WithinGroupSystem` returns class-level edges (doc pages, implemented
equations, referencing functions) but the only `py:attribute:` node was
`…WithinGroupSystem.loss` (degree 2) — `.resolvent` and `.gains` had NO nodes at
all, despite ~75 consumer lines. The graph surface contributed **0 of the 75**
sites. This is L-001's "graph alone under-scopes" at its extreme: for a FIELD
(as opposed to a function/class/method) rename, **text-grep is the primary
evidence and the graph is at best a way to find the owning class's doc pages.**
Don't spend a round-trip on `impact`/`callers` for a field; spend it on
`grep -rn "\.<field>\b|<field>=|<field>:"` plus a `replace(obj, <field>=` /
`getattr` / `asdict` sweep for dynamic access.

**(b) Check whether the OLD token is a substring of a common word before
proposing any replace strategy.** `gains` is a substring of **`against`**
(a-**gains**-t) — 679 occurrences in `orpheus/`+`tests/` `.py` alone. A bare
`sed s/gains/…/g` or an unanchored `replace_all` corrupts every one. The same
class of trap: `loss` ⊂ `lossless`, `space` ⊂ `namespace`, `role` ⊂ `payroll`,
`gain` ⊂ `bargain`. It also poisons the CENSUS, not just the edit — my first
grep of the owning file reported 21 `gains` hits where 2 were `against` in
prose.

How to apply: for any rename audit, (1) skip the graph for fields and go
straight to anchored greps; (2) run one `grep -c "<newword-containing-old>"`
sanity probe — or simply grep the old token with `-w`/`\b` anchors AND without,
and report the delta as a hazard line in the deliverable. Report the anchoring
requirement as an explicit instruction to the implementer, since a rename is
usually executed by a different agent than the one that audited it.

---

## L-006 -- A "shape probe" is not always a missing predicate — split boolean-presence from integer-width before proposing a typed swap

Asked to collapse N value-based `arr.shape[-1] > 1`-style probes into one typed
predicate, the load-bearing finding was that the probes split into two KINDS with
opposite fates: Kind-A pure-presence ("does this axis exist?", boolean → swap to
the typed predicate) vs Kind-B width/count ("the actual `2^d`, needed for buffer
ALLOCATION → these are honest counts, KEEP them). Proposing to delete the width
derivations would have broken allocation. A second constraint that governs such
work: a typed factor may live on the FIELD, but the inner-walk sites that do the
probing often see only a bare ndarray + `mesh.scheme` — so the "clean predicate
swap" is really a small-plumbing change, not a one-line rename.

How to apply: before recommending "collapse these probes into one predicate,"
classify each probe boolean-presence vs integer-width, and check whether the
probe site even HAS the typed object in scope. Report the verdict as
"(B) small plumbing," not "(A) clean swap," when the factor isn't reachable at
the site.
