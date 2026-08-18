# A test declares what it is ABOUT, and blindness becomes computable

⚠⚠ **STATUS: PROPOSAL, 2026-08-18. Nothing here has landed.** Every design
statement below is a *hypothesis* unless it carries `[M]`. The measurements
are real and mostly come from issues filed earlier by other sessions; the
architecture built on them is not yet verified anywhere.

## Relationship to the running campaign

The live campaign is `.claude/plans/test_dependence_and_claim_kind.md`
("A red partitions the suite, and a coverage claim can be falsified"). That
plan **instruments** the current architecture — per-test attribution
(nexus#57, CLOSED), module-precise imports, type-only import edges
(nexus#88, CLOSED).

This document asks a different question, at the user's request: *if we
redeveloped the test infrastructure with no constraints, what would differ?*
It does **not** supersede the campaign. Everything the campaign has landed is
an input here — the instrumentation is what makes the redesign measurable
rather than aesthetic.

---

## 1. The goal, in the domain's terms

**A test is an argument that some proposition about the code is true.** The
infrastructure should therefore be able to answer, mechanically:

- *What is this test about?* (its **subject** — one thing, not the hundreds it touches)
- *What proposition does it support?* (its **claim**)
- *What makes the expected value true?* (its **oracle**)
- *What would have to be false for this test to be worth running?* (its **premises**)
- *What can this test NOT see?* (its **blindness**)

Today the first, fourth and fifth are unrepresentable, and the second and
third live in prose. That is the whole of the gap.

**Done when** (the outcome, not a mechanism): for any test in the suite, each
of the five questions above is answered by a query rather than by reading the
test.

---

## 2. What is wrong today — measured, by other sessions

⭐ These are frozen issue measurements. Cite the issue; do not re-copy the
numbers into new documents (`plan-authoring` §9).

### 2.1 The relation a cone needs does not exist

`[M]` **nexus#60 (OPEN)** — `calls` reachability is **symmetric**
("A and B both touch P"). A cone requires an **antisymmetric** "B rests on A".
Consequences measured on ORPHEUS: a naive cone invalidates **31–36 % of the
suite** on a single red and is **not filterable** (excluding `external`
traversal moves it 0.6 pp). Separately: the static call graph has **0 %
recall** on the relation that matters for protocol-heavy code — 0 of 21 tests
reach the method under test via `calls` closure at depths 4/8/20/50, while 7
reach it at runtime.

⟹ #60's own diagnosis, reached independently of this document: *"what seems to
be missing is a **subject** relation — which code node a test is about, as
distinct from the many nodes it merely touches."* **This proposal's `Claim`
object is that subject relation**, and it answers #60's stated open question
("derivable or declared?") with: **declared**.

⛔ And #60 supplies the reason not to try deriving it: `[M]` the analogous
*claim-kind* axis is **not** derivable — `assert_allclose` appears in 218
files serving three different kinds of claim.

### 2.2 The marker layer cannot carry the semantics

`[M]` **nexus#61 (OPEN)** — nexus lifts markers by AST-parsing decorators and
recognises exactly **four** facts (`vv_level`, `verifies`, `catches`, `slow`).
Everything else the project registers is invisible: `regression`,
`foundation` (1515 usages across 308 files), `sentinel`, `cap` → **0 nodes
each**. And AST parsing cannot see *resolved* marker state: `vv_level`
resolves on only **1524 of 5273** collected tests, because **254 files** use
module-level `pytestmark` and the conftest applies a 5-rule precedence at
collection.

⚠ **This refutes a claim made in conversation on 2026-08-18** ("the AST layer
already reads decorator arguments, so claims-in-code works today"). It reads
*four hardcoded markers*. The mechanism is narrower than it appeared.

### 2.3 The current axis credits the wrong file, in both directions

`[M]` **ORPHEUS#334 (OPEN)** — a file-level `verifies(...)` list naming nine
labels writes an edge from **every** test in the file to **every** label. One
label collected ~21 `tests` edges of which ~2 exercise it. `[M]` nexus#61
measures the general case: blanket file-level markers generate **50.5 %** of
all `tests` edges.

⭐ The sharper half of #334, and the strongest single argument for this
redesign: the **strongest** gates for that same label deliberately carry
`foundation` and no `verifies` (the verifies⊥level doctrine), so the audit is
simultaneously **inflated by incidental rows and blind to the real battery**.
No amount of marker discipline fixes that, because "this file mentions X" and
"this test is about X" are different relations and the marker layer has only
one.

`[M]` **ORPHEUS#378 (OPEN)** — **139 of 456** test files carry no `verifies`
or `vv_level` at all, and the audit cannot distinguish *"no marker because
nothing is claimed"* from *"no marker because someone forgot"*. The failure
direction is the safe-looking one: coverage is **understated**, so the
response is "write a test that already exists".

### 2.4 The authored path is less controlled than the guessing path

`[M]` **nexus#86 (OPEN)** — `_infer_implements` consults the ontology before
minting an edge; `apply_pending_edges` (the *directive* path) does not. A
directive minted `py:data:` targets that `[edge.implements]` forbids, and the
build was `-W` green with no warning at any severity.

⛔ **This is a standing obligation on the present proposal, not a curiosity.**
A claim layer is an authored path. If it ships without admission control it
inherits this exact bug class, at `confidence=1.0`.

### 2.5 The negative space is not representable

`[M]` **nexus#85 (OPEN)** — an equation that *cannot* have an implementer is
indistinguishable from one whose implementer nobody declared. On a
hand-verified ground-truth set, **11 of 56** equations have no implementer at
all and attract **213 guesses**, **17 %** of all guesses on the set. #85 also
supplies a ready-made vocabulary: the four kinds are *identity*, *law enforced
by absence*, *canonical form not realized*, and *definition with no
declaration site*.

---

## 3. The design

### 3.0 The one structural change everything follows from

**Today a pytest function is four things welded into one object**: the
proposition, the argument for it, the coverage claim (`catches(...)`), and its
own quality certificate (green ⇒ believed).

Every anti-pattern in `vv-principles` (26 and growing) is a failure of that
weld — a blind gate, an unverified `catches` marker, a tautological guard, a
Mode-12 invariant functional. Each reduces to *"the evidence does not bear on
the claim"*, and **that sentence is unrepresentable today**, which is why it
lives in a list agents must remember instead of in a check.

⟹ **Separate the CLAIM (authored) from the EVIDENCE (derived).** Most of the
list then becomes a type error or an unmet obligation.

### 3.1 Ordering — three distinct orders, currently one flat namespace

| order | question | source | status today |
|---|---|---|---|
| **reach** | which tests are implicated at all | coverage attribution | ✅ landed (nexus#57) |
| **macro-DAG** | which cluster before which | **premises** between claims | ⛔ absent |
| **mini-DAG** | which test first within a cluster | measured discriminating power | ⛔ absent |

**The macro-DAG comes from premises, not from categories.** *"Theorems before
observables"* is right, but not because a theorem outranks an observable — it
is because **a failed premise makes the observable's result uninformative**.
Order by declared entailment and the category ordering falls out as a
consequence. This also keeps clear of the standing ruling (user, 2026-08-17)
that the L0–L3 layer order gives DIRECTION and never REACH: premises give
reach honestly; claim-kind never could.

**The mini-DAG is measured, not authored.** *Proposed (2026-08-18, NOT
verified):* **kill-set inclusion is the partial order** — `A ⊑ B iff kills(A)
⊆ kills(B)`. Four consequences:

- the order is **derived**, so it cannot drift from the tests;
- run order = the DAG, **cheapest-first within each antichain**;
- a test whose kill set is a *subset* of another's **and** costs more is
  provably redundant — the architecture says what to delete, which nothing
  does today;
- **an empty kill set is a test that cannot fail** — the tautological-gate
  class stops being something reviewers must spot and becomes a computed
  property.

`[M]` The measurement is affordable: `tests/_mutation/README.md` records
cosmic-ray 8.4.6 at **~1.0 s/mutant**, **374 mutants for `diamond.py` ≈ 6 min**
single-process, **373/374 killed**, and survivors already map to their
enclosing function via `definition_name`. Scoping is per-module, so the matrix
is incremental.

⛔⛔ **THE MISTAKE TO AVOID — do not linearize.** `[M]` **ORPHEUS#358 (OPEN)**
states it directly: a topological *sort* discards the branching, and with it
exactly the information that makes the partition possible — you can no longer
tell which later reds are in the failing node's cone and which are in an
unrelated branch. **The DAG must survive to REPORT time, not just schedule
time. The deliverable is `items + edges`, not sorted items.**

⚠ This corrects the version of this design first stated in conversation
(2026-08-18), which emphasised *ordering*. Ordering is a means; the
**partition** — valid / flip point / invalidated / not-run — is the outcome.
`[M]` #358's stake: a **1:08:43** whole-tree gate (9325 passed, `-m "not
slow"`, 2026-08-12).

**The explicit objective:** order by **expected bits per second** —
`P(fail | this symbol changed) × discrimination ÷ cost`. `sentinel` is already
a hand-curated estimate of exactly that quantity; measuring it makes `sentinel`
*derived*.

### 3.2 Category — typed, with obligations

Do not tag categories. The discriminating axis is **what makes the expected
value true** — where the oracle comes from — and each category carries
*required fields*, so an incomplete claim is unconstructible:

| category | oracle | obligation it must discharge |
|---|---|---|
| definitional | the type's own construction | names the invariant; may not claim more than the type promises |
| **theorem** | derivation from named premises | **must name premises, each itself a claim** |
| law / invariant | a quantified property | names the quantifier domain and one witness outside the trivial set |
| **observable** | external mathematics | names its pillar and its independence chain |
| convergence | asymptotics | names the limit **and** the rate |
| characterization | the past | states what it does **not** claim |
| agreement (L4) | another implementation | names its L0–L2 backing |
| **declared-empty** | — | one of nexus#85's four kinds |

`[M]` Category is **not** derivable (nexus#60: `assert_allclose` in 218 files,
three kinds), so it is authored — and it is the *only* axis of the claim that
must be.

### 3.3 Epistemology — the reason to do this at all

Make the epistemic content first-class and machine-checked:

- **Strength is measured, not asserted** (kill set).
- **Blindness is computable**: evidence whose kill set does not intersect the
  claim's threat model.
- **Mode 12 becomes a query**: does the measured functional's stabiliser
  contain the threat class?
- **Belief propagates**: a theorem's confidence is bounded by its *weakest
  premise*. Today a theorem test can be green while its premises are untested
  and nothing notices.
- **The negative space is storable**: *"nothing tests this"* and *"this is
  untestable, and here is why"* become different states (nexus#85's four
  kinds), instead of both looking like a gap.

⟹ The suite answers **"what do I now believe, and why"**, and reports what it
*cannot* see.

---

## 4. Where claims live — DECIDED, and why the fork dissolved

An earlier framing (2026-08-18, same session) posed *corpus vs code* as a
fork. ⛔ **REFUTED — it is a false dichotomy**, and the user caught it.

`[M]` Sphinx **renders 48 of 921 modules (~5 %)** and **zero** of `tests/`,
yet the graph carries **2748 `tests` edges**, every one minted from a
decorator in a file Sphinx never renders. The AST channel and the directive
channel are independent, and **the graph does not care which side authored a
node**.

**Decision: claims are authored IN CODE.** The argument that settles it —
and which neither side of the original fork stated:

> **Premises must be references, not strings.** The macro-DAG is built from
> "claim A rests on claim B". In the corpus that edge is an ID string that
> resolves if someone spelled it right — the exact failure the error-catalogue
> campaign spent a session fixing for `catches`. In code it is an object
> reference: pyright checks it, a rename moves it, a deletion breaks the
> build. No corpus design can provide that.

**Consequences of the decision:**

- **A claim does not live on a test.** Because the relation is one-to-many
  (one subject, many tests of increasing strength), a claim is a standalone
  module-level object near the code it constrains. Tests *reference* it. This
  makes the one-to-many explicit and countable instead of emergent from
  directory layout.
- **Equations stay in the corpus.** They are prose-and-math artifacts with
  `:eq:`-resolvable labels (903 nodes). A claim *references* an equation; it
  never restates it. Clean seam: the corpus owns the mathematics, the code
  owns assertions about the implementation, the graph joins them.
- **Publication needs one new directive.** A claim in code reaches the graph
  but not the rendered artifact. An autodoc-style `.. autoclaim::
  <module>` imports a module and emits its claims. Arguably better than
  today's default, where a page must *remember* to mention a claim.
- ⛔ **Rejected middle option:** claim declared in code, prose in the corpus,
  joined by ID. Two homes for one concept — the twin-source smell. The prose
  belongs in the object's docstring.

### 4.1 `[M]` What the code channel captures today, and what it does not

Measured 2026-08-18 by running `CodeVisitor` on a representative claim
declaration:

| | result |
|---|---|
| module-level annotated assignment → node | ✅ `py:data:<module>.<NAME>` |
| the **category**, via the annotation | ✅ `metadata["annotation"] == "Theorem"` |
| `file_path`, `lineno` | ✅ |
| **constructor kwargs** (`subject=`, `premises=`, `equations=`) | ❌ **not captured** |
| docstring on the object | ❌ not captured |

⟹ **Identity and category are free today; the fields are the work.** The
precedent for that work already exists and is close: `_parse_pytest_markers`
(`ast_analyzer.py:130-245`) already extracts literal positional **and keyword**
arguments from decorator `ast.Call` nodes. The same literal-extraction applied
to an assignment's RHS `ast.Call` is the mechanism.

⚠ But per nexus#61 the **test↔claim link must NOT be lifted from decorator
AST** — it is marker-like, and 254 files apply markers via module-level
`pytestmark` plus a conftest precedence that only collection resolves. `[M]`
The collection manifest join is proven: **9861/9861 items, 5273/5273 keys,
100 %**, `--collect-only` costs **7.0 s**.

⟹ **The split:** claim objects → AST. Test↔claim links → collection manifest.

---

## 5. What the issue review changed

Recorded because the design before the review is preserved in the session, and
the deltas are the value of having read them.

| from | change |
|---|---|
| nexus#60 | ⭐ **Validated the whole design from an independent direction** — it names the missing *subject* relation and asks whether it is derivable or declared. This proposal answers: declared, and the Claim is it. Also supplied the measured reason not to derive category. |
| ORPHEUS#358 | ⭐ **Corrected the deliverable**: partition, not sort. The DAG must survive to report time; `items + edges`. |
| nexus#61 | ⛔ **Refuted** "AST already reads decorators": four hardcoded facts only. Forced the AST/manifest split in §4.1. |
| ORPHEUS#334 / #378 | Supplied the measured case that the marker layer *cannot* express subject — inflated and blind simultaneously. |
| nexus#86 | ⛔ Added a hard obligation: the authored path needs admission control at least as strict as the derived path. |
| nexus#85 | Supplied the declared-empty category and its four kinds. |
| nexus#87 | Confirmed the authored/derived boundary is the right variable, and named what authored edges need that derived edges get free: both-end queries, provenance, a dangling report. |

---

## 6. Phases

Each names an outcome. Means are proposed and unverified unless marked.

### P1 — a test's SUBJECT is a fact, not an inference
**Goal.** For a test, "what is this about?" is answered by one code node.
**Proposed means:** a `Claim` object + a link from tests to it.
**Done when:** nexus#60's mechanism criterion holds on a fixture —
`invalidated_by(<the root test of the tree>) == {}`. ⚠ #60 explains why the
obvious forward-only assertion is the wrong gate: it *passes on the broken
symmetric relation*, so only the antisymmetry check discriminates.

### P2 — every registered marker is queryable, at its true scope
**Goal.** A project can add a marker without a nexus release, and a
file-level marker is distinguishable from a per-test one.
**Proposed means:** nexus#61's collection manifest.
**Done when:** `regression` / `foundation` / `sentinel` / `cap` return
non-zero node counts (all **0** today), and module- vs function-scope is
distinguishable in the reply.

### P3 — a claim cannot be malformed
**Goal.** An incomplete claim does not exist.
**Done when:** a `theorem` with no premises raises at construction (with a
test proving it), the four declared-empty kinds are expressible, and — per
nexus#86 — a claim-authored edge that the ontology forbids is **refused**,
with a witness that is refused. (§6c: the gate lands with the case it catches.)

### P4 — strength is measured
**Goal.** "How strong is this test" stops being a human estimate.
**Done when:** the suite reports (a) every test with an empty kill set, and
(b) every test whose kill set is a subset of another's at higher cost.
⚠ Sizing: a full matrix is hours, not minutes; it must be incremental and
content-hash keyed from the start, not retrofitted.

### ⏸ COMPACTION POINT — after P4, before P5

### P5 — a red partitions the suite
**Goal.** ORPHEUS#358's outcome.
**Done when:** on a deliberate red, the run reports valid / flip point /
invalidated / not-run as **items + edges**, and re-running `descendants(X)` is
a defined operation.

---

## 7. Risks, and what would refute this

- ⛔ **The authored surface may not get authored.** `[M]` ORPHEUS#378 measures
  139 of 456 test files with no marker today. A claim layer is *more* work per
  test than a marker. If the authoring does not happen, this is strictly worse
  than the status quo — a second empty axis. **Mitigation to design in from
  P1, not later:** the claim layer must be introducible incrementally and must
  degrade to today's behaviour where absent.
- ⛔ **Kill-set measurement may not scale to the whole tree.** Priced per
  module (§3.1) and never measured suite-wide. If the matrix cannot be kept
  fresh, the mini-DAG order silently goes stale — and a stale derived order is
  worse than an honest authored one, because nothing signals it.
- ⚠ **This proposal has no consumer yet.** Same shape as `cap`, which is
  measured write-only today. **P1 should not land without the P5 consumer at
  least sketched**, or it becomes another well-designed taxonomy nothing reads.
- **What would refute the core bet:** if `subject` turns out to be derivable
  after all — e.g. the narrowest node a test exercises is, in practice, the
  node it is about — then the authored claim layer is unnecessary for the DAG
  and only the *category* axis needs authoring. `[M]` nexus#57's per-test
  attribution is exactly the data needed to test that, and it is landed.
  **This is the cheapest experiment in the document and it should run first.**

---

## 8. ✅ THE EXPERIMENT RAN — 2026-08-18. The strong bet is REFUTED; the
## design survives in a cheaper form

**Question.** Is a test's SUBJECT derivable from what it executes?

**Method.** Ground truth = the file-naming convention `tests/.../test_X.py`
↔ a production module named `X.py`, taken only where the basename is
UNIQUE. That truth is independent of coverage, which is what the candidate
rules read. Candidates restricted to PRODUCTION nodes — a first pass let
test-tree nodes compete and a helper defined in one test file is trivially
the "narrowest", which refuted a strawman (R1 read 0.0 %; repaired, 73.6 %).

**Slices.** `geom_ctx` (existing) and `num_ctx` (captured for this
experiment: `tests/numerics`, 2344 passed, 567 s, **467 MB** report →
**3.37 MB** sidecar, 1303 attributed nodes, 0 unresolved contexts).

| rule | geometry (146 tests / 5 files) | numerics (531 tests / 25 files) |
|---|---|---|
| ceiling — truth is executed at all | 79.5 %¹ | **94.2 %** |
| **R1 narrowest executed node** | 63.0 %¹ | **78.2 %** |
| R2 most-exclusive module | 61.0 % | 68.4 % |
| R3 modal (most nodes) | 56.2 % | 68.5 % |
| truth ranks #1 by node count | 82/116 | 364/500, median 1, **p90 3** |

¹ depressed by a known bad truth label — `tests/geometry/test_geometry.py`
maps by basename to `orpheus.moc.geometry`. Excluding it, R1 = **73.6 %**,
ceiling 92.8 %.

### What it decides

⛔ **The STRONG bet is refuted.** "Subject is not derivable, so it must be
authored from scratch" is false: a naive rule recovers it ~74–78 % of the
time, and the true subject ranks #1 in ~73 % of cases and within the top 3
at p90.

⚠ **The weak form is not licensed either.** 22 % wrong is unusable as an
*authority* for an invalidation cone — ORPHEUS#358's outcome is a partition
you can TRUST, and a cone rooted on the wrong subject silently invalidates
the wrong branch. Being right 3 times in 4 is the worst kind of instrument:
right often enough to be believed.

⟹ **The synthesis, and it changes the economics of the whole proposal.**
Derivation is a strong PRIOR, not an authority. Seed the authored subject
from R1, show the top-3, and make DISAGREEMENT the review queue. The plan's
own #1 risk was *"the authored surface may not get authored"* (`[M]` 139 of
456 files carry no marker today). That risk drops from *author 5273 subjects*
to *review the ~22 % where the rule is unsure* — an order of magnitude.

⭐ **And the real value of the authored layer is not the 78 %, it is the
cases that are not a node at all.** `[M]` `test_bc_universal_invariants`:
52 tests, **12 distinct derived subjects, top answer 19 %** — the file
asserts universal invariants over *every* boundary law, so the derived
subject tracks the parametrisation rather than the claim. That subject is a
statement about a FAMILY. No executed node is it, at any accuracy. This is
what a `Claim` expresses and a derived node cannot.
⚠ Honest counter-observation: `test_reduced_operator` scatters just as badly
(8 answers, 24 % top share) while being module-named — so "concept-named vs
module-named" is NOT the discriminator. Across files: concept-named average
4.1 distinct answers / 55 % top share, module-named 3.0 / 72 %. Real but
weak; do not build on it.

### ⚠ Limits of this measurement — read before quoting the 78 %

1. **The scored population is definitionally the EASY case.** Only tests in
   files whose stem uniquely names a production module are scoreable — 531 of
   ~2344 numerics tests. The unscoreable remainder is exactly where
   derivation is expected to be worst, so **78 % is an optimistic estimate
   for the suite.**
2. **Ground truth is itself a proxy** (a naming convention), not a hand
   label.
3. **The hardest regime was NOT tested.** nexus#60 measured 0 % call-graph
   recall on the protocol-heavy SN sweep spine; neither slice here is that.
   Notably the "adversarial" numerics slice scored BETTER than geometry, so
   my prior about which code is hard was wrong — which is itself a reason to
   distrust the extrapolation.
4. `[M]` R1 (narrowest) beats R3 (modal) on both slices. The subject is the
   SPECIFIC thing a test reaches, not the most-executed thing — the opposite
   of what the rule list in §3.1 implied by ordering.
