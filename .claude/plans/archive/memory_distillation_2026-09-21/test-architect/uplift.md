# Rule-uplift candidates — test-architect memory distillation, 2026-09-21

Two candidates, both drawn from the digest's §0 meta-lessons (the
cross-campaign syntheses the prior pass wrote). Both generalise past this
agent: Mode 12 and X2/§2 are consulted by qa (mutation-battery review),
method-implementer (gate design at build time) and any agent designing a
census, not only test-architect. Both are evidenced by ≥7 independent
campaign pointers in the digest, which is the bar this file uses for
"generalises, not a one-off."

## Candidate 1 — the stabiliser has TWO sides, not one

**Joins:** `vv-principles` skill, Mode 12 ("Invariant-functional gate").

**Current clause (reproduced for context, not restated in the digest):**
"enumerate the functional's stabiliser (spectra: similarity + transpose;
balance sums: cancelling per-term errors; normalised shapes: scaling;
trace/det: similarity), intersect it with the threat model... for the WHOLE
committed gate set... gate the OBJECT where a mutation sits inside the
stabiliser (DESIGNED-GREEN)."

**What it does not carry:** every worked instance in the current clause is
about the OPERATOR/functional side of a gate (`G` acting on `A`). The
digest's M3 meta-lesson (`lessons.proposed.md` §0, 7 archive pointers:
`L43e`, `L47a`, `L47e`, `L58b`, `L59a`, `L61a`, `L84i`) found a second,
independent stabiliser that Mode 12's own worked examples never enumerate:
the SPACE side. A rank-1 point axis makes every weight a one-element array —
a SCALAR — and a scalar commutes with every operator, so `[G, A] = 0`
regardless of what `G` or `A` are; a ratio functional is blind to a uniform
scale of the space; a spectrum, to a similarity transform of the space's
basis; a palindromic generating rule, to its own reversal; a symmetric
generating rule, to half the permutation group acting on the space. None of
these five involve the operator's algebraic class at all — they are facts
about the SPACE (or the group acting on it) the gate is evaluated on, found
by degenerating the fixture, not the functional.

**Proposed clause addition (Mode 12, as a new sentence in the existing
paragraph, not a new numbered mode — the mechanism is the same stabiliser
concept, one more place to look):**

> check: enumerate the stabiliser on BOTH sides — the OPERATOR side
> (`[G, Aᵀ] = 0`, the existing four shapes) AND the SPACE side: does the
> fixture's own axis collapse a weight to a one-element array (a scalar,
> which commutes with everything)? Is the gate's functional a ratio (blind to
> uniform scale), a spectrum (blind to similarity), a palindrome (blind to
> reversal) or built from a symmetric generating rule (blind to half the
> permutation group)? A space-side stabiliser is as designed-green as an
> operator-side one, and is invisible if only the operator's algebraic class
> is checked.
>
> tell: a rank-1 or otherwise degenerate fixture credited for a gate whose
> functional is a ratio, spectrum, palindrome or symmetric-rule construction.

**Founding case:** `lessons_archive.md` §L84i (Consumers campaign step 2
DELTA) is the clause's own worked example — a 1×1 adjoint lift on a rank-1
point axis read as bit-identical under a ratio-valued gate while the
underlying object had genuinely moved; the archive section carries the full
measurement.

## Candidate 2 — a universal's AXIS is a design choice, not a formality

**Joins:** `plan-authoring` rule, §2 (the census-clause list: QUANTIFIER,
CONSTANT-DENOMINATOR, PREDICATE, STRUCTURAL-DENOMINATOR, ...). This is the
better home than `instrument-doctrine` X2 itself: X2 states the law once,
§2 is where ORPHEUS keeps the enumerated, named failure shapes the law takes
here, and this candidate is exactly that shape of entry.

**What the existing clauses do not carry:** every existing §2 clause is
about what a denominator MUST state (its predicate, its exclusions, whether
it is held fixed). None of them asks the prior question: *which axis does
this universal range over at all* — and a corpus that is uniform along the
axis actually chosen leaves the claim witness-less even when every existing
§2 check (predicate stated, denominator explicit, exclusions named) passes
cleanly. The digest's M6 meta-lesson (`lessons.proposed.md` §0, 11 archive
pointers: `L59b`, `L60e`, `L64f`, `L65d`, `L66d`, `L70e`, `L71d`, `L74b`,
`L76f`, `L77b`, `L86e`) is the densest of the eight meta-lessons and found
the SAME failure — a technically well-formed `k of N` claim that is
witness-less because `N` was counted along the wrong axis — recurring across
eight distinct axes: per ARM, per MEMBER of a union, per BRANCH of a
dispatch, per CONSUMER, per CALL SITE, per ROW of a parametrize, per FAMILY
of the corpus, per SUB-FAMILY of a tolerance.

**Proposed clause (new §2 item, named in the file's own convention):**

> **AXIS-CHOICE** A universal's axis (per ARM / per MEMBER of a union / per
> BRANCH of a dispatch / per CONSUMER / per CALL SITE / per ROW of a
> parametrize / per FAMILY of the corpus / per SUB-FAMILY of a tolerance) is
> itself a design choice, prior to and independent of stating the
> denominator on the chosen axis. check: before publishing a `k of N`, ask
> whether the corpus is uniform IN THE DISCRIMINATING FIELD along the axis
> chosen — a corpus uniform there leaves an arm (or member, or branch, ...)
> witness-less no matter how precisely `N` is counted; a single non-uniform
> member can be what makes the invariant spellable at all. tell: a
> well-formed `k of N` whose `N` was the easiest axis to enumerate, not the
> one the claim is actually about.

**Founding case:** `lessons_archive.md` §L76f — a census correctly counted
"14 of 14 isotopes" while the discriminating axis was actually per-CHANNEL
within an isotope (2 of 13 isotopes carry an empty (n,2n) channel), so the
per-isotope `N` was witness-less for the per-channel claim being made.

## Considered and not proposed

- **M1** (construct-and-measure before designing the gate, 13 pointers — the
  most-evidenced of the eight): proposed instead as an **AGENT.md**
  promotion (see `agent_md.md`), not a rule uplift — it is a standing
  test-architect WORKFLOW step (what to do before drafting a matrix), not a
  general evidentiary law a rule states; `instrument-doctrine` X1 already
  states the general law ("an instrument is evidence only if some realizable
  state changes its reading") and M1 is this agent's operational form of it,
  which is what AGENT.md is for.
- **M2, M4, M5, M7, M8**: each generalises less cleanly — M2 (four
  null-hypotheses) and M7 (identity/route/count instruments) are close
  restatements of `instrument-doctrine` X1's existing "name the mutation"
  check with ORPHEUS-specific instrument names, not a new clause; M4 (a
  brief's claim is a hypothesis) and M5 (the tree moves under you) are
  general process discipline already substantially covered by
  `process-discipline`'s "Trust git for merge status" section and by this
  agent's own high pointer-density (9 and 7 respectively) reading more as
  "recurs often for test-architect" than "missing from the rule corpus." Not
  proposed; conservative per the brief.
