---
name: instrument-doctrine
description: PROACTIVELY load when building, judging or citing evidence — a gate, a metric, a census, a canary, a mutation battery, a ratio, a timing. Procedures for the four always-on statements of the instrument-doctrine rule (X1 an instrument must be able to fail; X2 every claim carries its population and instrument; X3 prose is not enforcement; X4 one definition per quantity). Preloaded by qa, test-architect, numerics-investigator and archivist; vv-principles cites it.
harness:
  kind: skill
  budget_tokens: 3200
---

# Instrument doctrine — the procedures

The rule (`instrument-doctrine`, always-on) states four things every claim is
held to. This skill is how to meet them. Each procedure ends in a number with
its predicate, or in a red, never in an adjective.

## X1 — prove the instrument can fail

**Positive control, before any battery.** Include one mutation that MUST redden
many gates; read an all-blind verdict as *the harness is broken* until the
control reddens. Mutate INSIDE the object's algebraic class (a mutation that
also breaks linearity, symmetry, positivity or a shape contract reddens for the
wrong reason). The verdict is a per-arm TABLE — a multi-arm guard is N claims,
and an arm that reddens nothing is a guard with no witness.

**Positive control, for a filter.** A census grep, a name net or a regex over a
tree is an instrument, and its zero obeys the same law: before believing it,
show one known member of every shape it must find (uppercase labels, relative
imports, CamelCase, a string-form `getattr("name")`); a staged filter (a name
net, then a literal scan) needs a control per STAGE, named with the spelling
you are least sure the net catches; a member the filter would find for the
wrong reason is a null control. `grep` here is ugrep, and an anchor inside an
alternation group matches nothing, silently (`code-search`). The nine
battery-specific checks are `vv-principles` #17; the nine classes of a gate
that fires and cannot fail, `vv-principles` mode 8; the design-time question a
plan asks of a new gate, `plan-authoring` §6c. A gate that ranges over a LIST also reports its INPUT COUNT beside its finding count, because a correct filter over an empty list and a correct filter over a clean tree print the same zero and the first is a broken harness: print `len(inputs)`, assert it non-zero, and assert one input resolves to an existing path (`[M]` 2026-09-21, twice in one session: an unsplit shell `$FILES` made an xref probe and a markup scan each read one nonexistent path and print a clean 0).

**Stabiliser enumeration, at gate design time.** Write down the functional the
gate measures and the group of errors it is invariant under: spectra are blind
to similarity and transpose; balance and telescoping sums to any per-term error
that cancels; normalised shapes to global scaling; trace and determinant to
similarity. Intersect with the threat model. A threat inside the stabiliser is
designed-green: gate the object, not the functional.

**Activation count, for a canary or an acceptance artefact.** Instrument the
carved path (a counter or a file write — not a bare `assert`, not a print) and
run the artefact: zero executions means the artefact passes under success and
failure alike, whatever its name. Ask the inverse question too: *if the
campaign fully failed, would this artefact move?*

**Full-success reading, for a metric adopted as a target.** Read the metric's
implementation and answer *"if this campaign fully succeeds, what does it
print?"* — for the named target AND every other instrument that will read the
result. Three shapes fail it: a proxy the work removes; a population the work
empties (the score over "the ones still broken" degrades as you fix the easy
ones); a predicate ranging over a set the design is forbidden to touch
(designed-red — intersect the tell's hits with the UNTOUCHED set).

**A repair with no gate.** A fix that makes a defect unspellable feels to need
no test, and its pinning gate is green either way. Mutate every repair you
believe is structural; re-measure the repaired module's own suite separately,
to tell *no witness* from *not measured*.

**Execution evidence.** A gate is evidence only if it RUNS under the canonical
invocation: check the asserting statement survives `python -O` (the scope is
`coding-standards` § "A bare `assert`"), check the marker set is not deselected by `-m "not slow"`, and
for a re-routed path check with a counting spy that the gate still reaches the
changed line.

## X2 — state the population and the instrument

**The census protocol.** (1) Write the predicate before counting: what the
numerator counts, over which tree (`orpheus/`, `tests/`, `docs/` — `tests/` is
usually the majority and usually unranged). (2) A membership question is
parsed — `ast` for code, the doctree for docs — never a line window, a regex
that can land inside a nested literal, or `| head`. (3) Validate the filter
with a positive control per shape and per stage (X1). (4) State every
exclusion in the claim: a file, a package, a directory, the definition site.
(5) For a completeness claim, run a SECOND, independently vocabularied filter;
one validates the other. (6) Report `k of N <predicate>`; list several counts
only with their separate predicates and trees. (7) Verify the DECODER as well as
the filter: a production predicate reused as a detector inherits its other
meanings; enumerate every state it maps to True and give each state its own
control.

**The shapes a denominator hides in**, each a `plan-authoring` clause cited
there by its tag. QUANTIFIER: a universal is `k of N <predicate>`; the
unmeasured members are where the defect survives. QUANTIFIER-AT-WRITE-UP: run
the check when publishing, not when measuring. CONSTANT-DENOMINATOR: the
denominator least likely to be written is the one held fixed across every row;
ask what is the same in every row and whether the conclusion is about that.
PREDICATE: a fraction names what its numerator counts, or it is not
re-runnable. DENOMINATOR-THAT-IS-A-GATE: a guard named as the denominator
bounds the arm it is called on, not the question; `grep -c` its call sites,
check the prescribed measurement is not tautological, and check that a scope
inherited from a code comment did not drop the sentence that limited it.
STRUCTURAL-DENOMINATOR: for "the rebuild loses X" or "these are duplicates",
enumerate against the TYPE (`dataclasses.fields(T)`), never against the
concept chased. EXCLUSION-IS-A-PREDICATE: an excluded file hides in-module
consumers, an excluded package reads as "0 production consumers", an excluded
directory measures the issue's folder and not the goal's corpus; state the
exclusion in the claim or filter by line, and measure one level up.
DEFINITION-FILE-ONLY: a self-check whose whole population is the defining
module certifies "the def plus one call"; a self-check is a census and owes
its population. RESUME-BLOCK-M: an `[M]` count in a
compaction point is a claim with a predicate; a one-file census lies beside
`[R]` neighbours that were checked for free. A-LIST-IS-N-CENSUSES: counts
listed in one breath imply a shared scope that never existed; each owes its
predicate and its tree. CARVE-FORKS-THE-DENOMINATOR: when N classes share one
body, "N surfaces" is two counts, bodies and role × surface rows; state both
and the arithmetic between them. RENAME-SIZE: a rename counts the identifier
and, separately, the concept spelled without it; word-bounded, then triaged by
meaning. STATE-THE-CENSUS-PREDICATE: "complete" states its method ("complete
for literal-name calls"). DUPLICATES-IS-A-UNIVERSAL: a bit-identity claim owes
its denominator and is routinely asserted without being run.

**Configuration, for a number.** Fixture, settings, what the fixture is blind
to (its kernel and regime), the `-m` filter of a pytest count, the exclusion
list of a census. A relayed number travels with all of it or does not travel.

**Draw-stable statistics.** A seed, a wall-clock timing and a hand-typed
prediction are all draws. Sweep seeds before writing "bit-exact"; time with
the min of ≥15 interleaved repeats and say so; compare a transcribed
prediction with a scale-free statistic (alignment, relative residual, rank),
never component-wise to the digits you typed. Pin gates on the absolute
statistic that is stable across draws.

**Ratios.** State both legs' populations side by side; if they differ, there is
no ratio. A positive control validates the instrument, not the comparison.

**Predicted-then-measured.** Explain every unit of a gap before publishing the
explanation; a registry-driven gate gains rows from corpus surfaces the
arithmetic never modelled.

## X3 — replace prose with an assertion

For every claim in a docstring, marker, label, table header or plan row, find
the line that fails if the claim is false. A `verifies`/`catches` marker is a
claim until a coverage capture or a mutation adjudicates it. A docstring
naming a shape, a bijection, a compatibility law or a primitive that the body
does not assert or call is a single-source divergence: assert it, return the
structure (a permutation makes its own bijectivity assertable; a `bool` does
not), or delete the claim. A plan row in a rule's vocabulary carries the
rule's `[M]`, not its phrase. A gloss beside a link ("whose census clause is X
and Y") asserts the target's scope, reads as navigation so nobody audits it,
and goes stale when the target grows: point, and let the definition speak.

## X4 — find the shared upstream

Ask independence per axis: the DERIVATION axis (no shared identity, integrand
or closed form) and the INPUT axis (no shared constructed object handed to both
sides). A verified primitive both sides call (one dense eigensolver, one
integrator) is neither axis: the SUT and its oracle stay independent when they
ASSEMBLE its inputs by different routes, and the primitive is trusted once it is
pinned against a domain-free closed form (a matrix with chosen eigenpairs), never
against an in-domain solver. Compare two "independent" implementations by α-normalised AST; an
α-equivalent pair is one implementation under two names, and its agreement is
a tautology. Before retiring a duplicate, name the mechanism that kept the
copies equal and grep the shortest distinctive fragment of its message for a
witness; if none exists, write it in the same commit — the retirement created
the exposure.

## Pointers

- Rule: `instrument-doctrine` (always-on). Instances, each citing its
  statement by ID: `vv-principles` #7, #11–#14, #17–#20, #22–#24, #26, #31,
  #34, #36 and test-design mode 8; `plan-authoring` §2, §6c, §8, §10;
  `retirement-audit` D.16, D.18 and its three searches; `coding-elegance`
  Pattern 2, Pattern 7, anti-patterns #1 and #20; Cardinal Rule 2.
- Evidence: [V&V anti-patterns](../evidence/vv-anti-patterns.md),
  [test-design modes](../evidence/test-design-modes.md),
  [plan-authoring evidence](../evidence/plan-authoring.md).
