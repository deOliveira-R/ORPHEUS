# Knowledge-structure tiering — the agent harness as the artefact (2026-09-20)

**New problem class.** Every prior entry in this library attacks MATHEMATICS (an operator, a
carrier, a quadrature). This one attacks a DOCUMENT CORPUS: the ORPHEUS agent harness
(`docs/development/` → `.claude/`, W5 re-evaluation). The class recurs whenever the brief is
"does this rule/doc set earn its cost, is each concept defined once, is each clause in the right
tier" — a rule corpus, a style guide, a test-policy document, a prompt library. The frames below
are the ones that FIRED, with the trigger that selected each; they are not transport frames and
they transfer to any tiered knowledge corpus.

Memo: `scratch/_harness_eval/reeval/attacker_memo.md` (W5, 2026-09-20).

## The structural features that select these frames

1. A document/section loads under a PREDICATE on session state, and the predicate ranges over
   several INDEPENDENT base spaces (here: file PATH via globs, workflow PHASE, ACTOR/agent role).
2. Two or more CITATION systems coexist, only one of which the build checks.
3. A generic statement (`X1`–`X4`) with instances scattered across other documents.
4. A cost (tokens × loads) against a logged value signal (surprises caught).
5. One source tree, N generated views, a target Protocol with no stated law.

## The four frames, with first tests

### Sheaf on a site — base = path × phase × actor

**Trigger** feature 1. **Reformulation** a clause is a SECTION over its SUPPORT (the set of
session states where it is non-vacuous). "Always-on" is not a tier: it is Γ(X), the section over
the terminal open. Tier order IS support inclusion, so per-session cost is MONOTONE in support and
the tier is DECIDED by support rather than negotiated against a budget.

⭐ **The durable criterion this frame produces — GLUING, i.e. "write the restriction and diff".**
For a fact stated on several opens (the `python -O` fact: `[M]` 23 occurrences over 10 files of
`docs/development/`), restrict the would-be global statement to each open and diff it against the
local text. **Zero delta ⟹ it is one global section: one home plus a citation. Non-zero delta ⟹
the delta is the whole of what the local text should say.** A presheaf that will not glue (the
copies disagree) IS the drift, named. This settles "restatement vs legitimate local instance"
without a taste argument, and it generalizes to any corpus.

**What it made visible** a base space with NO mechanism: PHASE. Clauses whose support is a phase
(a retirement audit, a resume-a-plan procedure, a compaction rule, a look-at-CI-after-push rule)
ride Γ(X) because nothing loads at a phase transition — even though the orchestrator OWNS phase
transitions, so the hook point exists. Second: a clause supported on `docs/**`/`.rst` sitting in
Γ(X) while the path mechanism already ships (`harness.paths`, `[M]` 1 of 9 rules used it).
**First test** add `paths:` to a second rule and run the generator's `--check`: the always-on sum
MUST drop by that page's budget with the drift test green. A page carrying its tier anywhere but
front matter FAILS — its generated body still asserts the old scope in prose.

### Relational schema, third normal form

**Trigger** features 2–3: IDs as primary keys, citations as foreign keys. **Reformulation** each
ID registry is a table; a gloss beside a link is a COPIED NON-KEY ATTRIBUTE, and its update
anomaly was already measured in-tree (a role block naming two of a census clause's four
requirements, staled by the commit that added the fourth).

⭐ **Criterion for one home vs a local instance:** a fact has ONE HOME iff it is functionally
determined by that table's key ALONE; it is LOCAL iff determined by the COMPOSITE key
(statement, site). So an instance is written as the JOIN ROW — "X2 applied to a census's
exclusion list" — and only the site-dependent half is spelled out.

**What it made visible** an INTEGRITY ASYMMETRY: link+anchor citations are build-checked (the
generator reds on a missing file or a missing MyST heading) while `[M]` 89 plain-text ID citations
(`X1`–`X4`, `#NN`, `Pattern N`, `ERR-NNN`, `Lnn`) over 12 files carry no constraint at all. And
the generic→instance relation is stored ONE-WAY, in ONE place, at §-granularity, over a corpus
keyed at CLAUSE granularity (one cited section held 29 clauses by the page's own appendix census).
**First test** a `--check` resolving every ID against its registry; seed a non-existent `X5` AND a
renamed registry entry and require BOTH reds — a today-green pass proves only that nothing happens
to dangle. **Attack** a clause not writable as `Xk[binding]` is NOT an instance: it carries a
mechanism the generic statement does not entail, so it is a fifth generic statement or genuinely
local. That is a decidable test where the corpus had judgement.

### Screening with a base rate (PPV) — the tier criterion that is not a token budget

**Trigger** feature 4. **Reformulation** a clause is a DETECTOR over a population. Its yield
`catches / loads` FACTORS as (catches per load INSIDE its support) × (fraction of loads inside its
support), and **only the second factor is a tier question**. Criterion: place a clause where that
fraction ≈ 1; a clause with high sensitivity and support-fraction ≈ 0 has PPV → 0 whatever its
token cost. **What it made visible** a target of the form "surprises per campaign trending to
zero" has NO DENOMINATOR, so it cannot separate a clause that has stopped catching (quality) from
one loaded in the wrong tier (allocation). **First test** recount the catch log per clause ID and
per tier; ≥3 catches all inside one phase ⟹ a phase section, not a global one. Discriminates
against "it caught something once, so it earns always-on".

Note the self-application: the corpus's own metric rule (report rule quality over a FIXED
population separately from how much corpus is affected) IS this criterion — the move was to apply
the corpus's rule to the corpus. Look for that move first on any self-evaluation brief.

### Functor with envelope invariance — what a SECOND generated view may vary

**Trigger** feature 5. **Reformulation** a harness is `F : Src → Rt` over a category of pages
(objects) and links (morphisms); the relinker is the action on morphisms and already carries its
law (target must exist, anchor must be one the parser mints). Two laws are missing and both are
checkable: **(i) BODY INVARIANCE** — every harness-specific token is an ENVELOPE (a stamp, front
matter, block markers, an injected `!cat` line), so `strip_envelope(F(page)) == page.body` up to
relinking; **(ii) SUPPORT IS A SOURCE FACT** — `always_on` is computed from `kind` + `paths`, so a
second harness may vary the loading MECHANISM and never WHICH clauses are global. Together they
state exactly what a second runtime may change, as a law rather than a prose convention.
**First test** add a second target with NO path-scoped mechanism: it must fold the path-scoped
rule into another tier, after which that rule's BODY asserts its old scope ("Applies when working
under `tests/`") under a harness where the claim is false. `[M]` 3 of 13 sources failed
body-invariance at the time of the attack.

## UNEXPLORED — frames checked on this class, with the structural reason

Carry these forward; they are the class's dead ends, not this corpus's.

- **Homology / chain complex** — no `∂²=0`; citation composes to a non-zero relation. (Third
  independent refutation of this frame in the library; it is baited by the word "boundary" and by
  "section", and it has never fired.)
- **Category theory past the one functor** — the whole win (one source, many views) is captured by
  functor + envelope invariance; no further law produces a test. (L-001's standing verdict holds
  on a non-transport artefact too.)
- **Tensor network / MPO** — bond dimension 1: one source, one harness. Fires at N ≥ 3 targets.
- **Rate–distortion / compression** — tempting because the brief is "cost vs value", but the corpus
  defines no RECONSTRUCTION TARGET, so there is no distortion measure. The value signal is a
  DETECTION statistic, which is why the screening frame fired instead. ⭐ Generalizable: when a
  brief is phrased as compression, check for a reconstruction target before reaching for
  information theory; absent one, the right frame is detection with a base rate.
- **Dependent types indexing a clause by its tier** — the tier is a RUNTIME predicate on session
  state, not a static index; nothing composes over it. The sheaf frame carries the payoff.
- **Group theory / quotient** — no group acts on the clause set.
- **CRDT / merge semantics** for the generated copies — the one-way flow makes them a pure function
  of the source and a drift test already gates it; no concurrent-edit problem exists.
- **Knapsack under a token budget** — this is the framing already IN USE; screening dominates it
  because the screening weights are measured (the catch log) while the knapsack's values are not.
- **Queueing / Markov session model** — per-session cost is deterministic given the tier.

## Cross-discipline pollination that fired (the analogue of Part B for this class)

Document-engineering disciplines have solved this exact problem; reach for them before inventing.

- **ISO/IEC Directives Part 2** — normative clause vs informative annex = core rule vs evidence
  page (already matched in-tree); the verbal-forms table (shall/should/may, one meaning each) is
  ALREADY validated here as the epistemic-marker rule, which is evidence the pattern transfers.
  The missing half is a TERMS CLAUSE: each defined term one number, never repeated later.
- **Legislative drafting** — the general-definitions rule is a COUNT, not taste: a term used in
  more than one section goes to the definitions section; used in one, defined there. That is the
  one-home-vs-local criterion, decidable by grep. The CFR's "reasonably available" condition on
  incorporation by reference maps onto the on-demand tier's link-checked anchors.
- **Conformance clauses (W3C / IETF, RFC 2119)** — a spec names CLASSES OF PRODUCT and states
  which requirements bind each. The harness's classes exist (orchestrator / key / support); a
  hand-pasted "rules that apply to you" line in every brief IS the requirement × class table,
  denormalized. Borrowing the TABLE is the fix, and it is the ACTOR axis of the sheaf made explicit.
- **DITA `conref` + `ditaval` profiling** — build-resolved transclusion by ID (a dangling conref
  reds exactly as a dangling anchor does) and conditional-processing attributes, which are
  `paths:` generalized to arbitrary conditions with the same generator shape. This is the concrete
  mechanism that would put plain-text ID citations under the integrity link citations already have.

## How to open the next attack of this class

1. Enumerate the BASE SPACES of the loading predicate and, for each, ask whether a MECHANISM
   exists. A base space with no mechanism forces its clauses into the global section — that is
   the highest-value finding and it is found by inspection, not by reading prose.
2. Sort citations into build-CHECKED and unchecked, and count both. The asymmetry is always there.
3. Run the corpus's own metric/instrument rule against the corpus.
4. Grep the SOURCE pages for TARGET vocabulary (here: `.claude/`, and any body sentence restating
   its own loading predicate). Those are Smell #16 shape 2 — one fact, two representations — and
   they are exactly what a second generated view falsifies.
