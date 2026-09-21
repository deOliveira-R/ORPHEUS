---
harness:
  kind: rule
  budget_tokens: 1500
  brief: >-
    a zero from a filter is evidence only after a positive control of each shape it must find (X1); every count states its predicate, its tree and its exclusions, a universal is `k of N`, and a completeness claim is re-run in Python (`re` + `pathlib.rglob`) so its denominator is stated (X2).
---

# Instrument doctrine — four statements every claim is held to

An instrument is anything whose reading you cite as evidence: a test, a metric,
a census, a canary, a mutation battery, a docstring, a plan row. The four
statements below apply to every artefact an agent writes — an issue, a commit
message, a brief, a plan, a report. The
procedures for meeting them are the `instrument-doctrine` skill; the
domain-specific instances live in `plan-authoring`, `coding-standards`,
`vv-principles` and `coding-elegance`; the skill's Pointers section maps each
instance to its statement.

## X1. An instrument is evidence only if some realizable state changes its reading

Before citing a gate, metric, canary or fix as evidence, name the input that
would make it read differently and show that the input exists in the tree today.
The same law read from the other side: a negative reading (a clean tree, a
zero, an all-green battery) is evidence only after a POSITIVE CONTROL has made
the same instrument read positive, because a broken instrument and a clean
tree print the same thing, and the broken one reads as "nothing to do".

- check: ask *"if the thing this guards were fully broken, what would it
  print?"* — answered by reading the implementation and counting activations,
  never by its name. A gate lands with the input, existing in the tree the
  moment it lands, that it rejects (its first red); for a metric or an
  acceptance artefact ask the inverse, *if the work fully failed, would this
  reading move?* For a repair you believe is structural, mutate it and
  require a red. Before believing a zero, name the known member the filter
  found or the mutation the battery reddened on.
- tell: a gate green before and after the change; a canary whose carved path
  executes zero times; a metric that moves the wrong way while the work
  succeeds; a docstring's claim standing in for a witness; a confident, empty
  answer with no control beside it.

## X2. Every claim carries its population and its instrument

A universal (*every, all, none, only*) is written as `k of N <predicate>`. A
number travels with the command, the fixture, the seed or repeat protocol and
the precision that produced it. A ratio states both legs' populations side by
side. A filter is validated against a known member of every shape before its
zero is believed.

- check: before publishing, ask *"how many did I actually look at, and is that
  number in the sentence?"* and *"what is held fixed across every row?"* A
  membership question is parsed with an AST, never grepped through a line
  window or a `head`.
- tell: a bare number; "every" with no denominator; several counts in one
  breath; two filters that disagree; an exclusion (a file, a package, a
  directory, the definition site) that is not stated in the claim.

## X3. Prose is not enforcement

A docstring, a marker, a label, a table header, a column's vocabulary, a rule's
phrase, a plan sentence or a gloss beside a link asserts nothing the code or
the gate does not. Assert the structure the prose names, return it so that it
asserts itself (a returned permutation makes its own bijectivity assertable; a
`bool` does not), or weaken the prose.

- check: for every claim in prose, find the line that would fail if the claim
  were false. A plan row written in a rule's own vocabulary is a summary of a
  check unless the row carries that check's `[M]`.
- tell: "by construction", "tautological", "unaffected because X", "no new
  edge either way" with no measurement beside it; a `verifies` or `catches`
  marker on a test the defect does not redden; a gloss that names part of its
  link's target and goes stale when the target grows.

## X4. One definition per quantity

Two spellings of one thing agree tautologically. Before calling agreement
evidence, find the shared upstream — an identity, an integrand, an input
object, an α-equivalent body, a single constant. Before retiring a duplicate,
name the mechanism that made it redundant and that mechanism's witness.

- check: independence is asked per axis, derivation and input; two bodies are
  one implementation iff they agree as α-normalised ASTs (every local renamed
  to a placeholder), whatever their names. The mechanism that keeps two copies
  equal has a witness, a test that greps the shortest distinctive fragment of
  its message; none means write one in the same commit as the retirement,
  which created the exposure.
- tell: `allclose(solver_a, solver_b)`; a gate whose two sides derive from one
  rule; a convention re-applied at N consumers.
