# Plan authoring — what a plan owes its own future author

A campaign plan in `.claude/plans/` is not a note to the reader. It is a note to
**you, after your context is gone**, and it will be read by someone with your
authority and none of your memory. That reader cannot tell a considered decision
from a placeholder, an outcome from a guess at the means, or a measurement from
a plausible number — unless the plan marks the difference.

## ⭐⭐ This file is a LIVING document, and it has a measurable target

A plan is the hardest case of the **articulation standard**: articulate = take a
concept apart so a fresh reader reassembles it *losslessly*, and **the LOSS is
the measure**. For prose the loss is hard to observe. For a plan it is not —
it shows up as a **SURPRISE**: a moment where a session, reading the plan in
good faith, believed something the tree contradicts. Every surprise is a
quantum of information that failed to survive the context boundary.

**So the loop is closed, and running it is part of the work:**

1. **Log the surprise** the moment it happens — in the plan (in place, per §3)
   and, if the mechanism generalises past that campaign, here.
2. **Ask what a plan could have carried that would have prevented it.** Not
   "what was wrong" — what *absent structure* let it through. A clause earns its
   place only if it would have caught its own founding surprise.
3. **Add the clause, with the surprise as its worked example**, dated and
   measured. A rule without its founding failure is advice; with it, it is
   evidence.

⟹ **the target is falsifiable: surprises per campaign trending to zero.** When a
session picks up a cold plan and finds nothing that contradicts the tree, the
transfer was lossless and this file is done growing. Until then it is not, and a
campaign that produced a surprise but no clause has thrown away its own
measurement.

⚠ **Two failure modes of the loop itself.** (a) Adding a clause for a
*one-off* — if the same absent structure cannot plausibly bite a different
campaign, it belongs in that plan, not here; this file is the floor, not an
archive. (b) Growing without distilling — clauses that have stopped catching
anything should merge or retire, exactly as the memory-distillation standard
requires, or the file becomes a hot surface that costs more to read than the
surprises cost to hit.

### Surprise log — each clause traces to what produced it

| date | the surprise | clause |
|---|---|---|
| 2026-08-06 | A phase title named a MECHANISM; its own author designed to it three compactions later and lost a cycle. | §1 |
| 2026-08-06 | An inherited `[M]` number (`\|B(x)\| = 1.320`) was a *different fixture's*; the resuming session measured `1.824`. Both correct, silently incompatible. | §2, §4 |
| 2026-08-06 | A plan's §5 blocker ("there is NO public entry point") was true when written and void two commits later, inside the same campaign. | §3, §7.3 |
| 2026-08-06 | A "forbidden to commit" note was a point-in-time snapshot; the work had been committed at `34af8474`. | §7.1 |
| 2026-08-06 | A plan split one signature change into three steps. Step 2 was unlandable alone: the retired probe was one of the two call sites of the signature being changed. | §6b |

Companion to CLAUDE.md **Cardinal Rule 4** (issues are the cross-session log) and
to the compaction-point discipline. Those say *where* state lives; this says what
a plan must CONTAIN so it does not lie to the next session.

> **The failure this rule exists to prevent, measured 2026-08-06.** Campaign
> phase **P6** was titled *"promote what P4 hand-rolled into production"*. That
> title named a **mechanism** (move the class), when the goal was an
> **outcome** (make the MMS expressible as an ordinary user-written source).
> The author of the title and the reader who designed against it were the same
> agent, three compactions apart. The reader spent a full design cycle on
> "where should the class live in `orpheus/`" before the user pointed out that
> moving the class was never the point — a verification artefact in the
> production tree would have been the *wrong* outcome. **A plan's phase title
> is the highest-leverage text in the document**, because it is what survives
> summarisation intact while the body is compressed away.

---

## 1. Every phase states its GOAL, separately from any proposed MEANS

The two must be visually distinguishable, because the goal is durable and the
means is a guess that later evidence may overturn.

```markdown
### P6 — a directional inflow is expressible without smuggling
**Goal.** A user can declare `q(Ω)` on a face without the source having to
carry trace knowledge through its constructor.
**Proposed means** (as of 2026-08-05, NOT verified): promote the MMS's
hand-rolled spec into production.
**Done when:** the MMS source is an ordinary user-written source, and
`tests/.../test_mms_declared_inflow.py` builds it with no per-face
`mu_inflow=` argument.
```

- **Title the OUTCOME, never the mechanism.** "Promote X into production" is a
  mechanism; "X is expressible without smuggling" is an outcome. If the title
  contains a verb that names an edit (*move*, *promote*, *rename*, *extract*),
  ask what the edit is FOR and title that instead.
- A means proposed before the investigation is a **hypothesis**. Label it, date
  it, and say it is unverified. An unlabelled means reads as a decision.
- **Done-when is a checkable predicate**, not a feeling. Prefer one that a
  grep or a test run can answer.

## 2. Mark the epistemic status of every claim

The reader cannot re-derive which is which, and the cost of guessing wrong is
asymmetric — acting on a guess labelled as a measurement is how a campaign
inherits a false premise.

| marker | means |
|---|---|
| `[M]` | **measured** — with the command or file:line that produced it |
| ⛔ / ⚠ | a hazard or a refutation, with what refuted it |
| *"proposed"* / *"hypothesis"* | not yet checked |
| ✅ LANDED `<hash>` | done, verifiable with `git merge-base --is-ancestor` |

A bare number with no marker will be read as measured. If you did not measure
it, say so or delete it.

## 3. A refuted premise is EDITED IN PLACE, never silently dropped

When investigation refutes something the plan asserted, leave the original text
and put the refutation beside it (`⛔ REFUTED <date> — …`). Deleting it loses
the reason and invites the next session to re-derive the same dead end;
rewriting it as if it were always right destroys the record of how the campaign
learned. The falsified text plus its refutation is worth more than either alone.

This is the plan-side twin of the `coding-standards` retirement rule: past-tense
history stays, present-tense falsehood is a MUST-FIX.

## 4. Numbers carry their CONFIGURATION

A measurement without its fixture is not reusable — and worse, it is *usable
wrongly*. State what was measured, on what, with what settings.

> ⛔ Measured 2026-08-06: a compaction point recorded `|B(x)| = 1.320` as P5's
> activation value. The next session's fixture measured `1.8239310798528774`.
> Both were correct; they were different fixtures. Had the plan pinned the
> inherited number in a gate, it would have been a false red with an
> authoritative-looking provenance.

Same rule for a sub-agent's number: **read its configuration before adopting its
value.**

## 5. State the goal in the DOMAIN's terms, not the tree's

"Make `InflowSourceSpec.evaluate` take a space" is a statement about the current
code, and it stops being meaningful the moment the code moves. "A boundary
source is a function of direction, so it must be told the directions" stays true
across refactors and lets a later reader re-derive the means. Write the second;
put the first under *proposed means*.

## 6. Sizing — a plan that is too big to read is not a plan

- **Compaction points** every ≥4 phases (see the compaction-point discipline),
  each carrying: the phase→commit table, corrections that supersede older text,
  the measured red baseline and gate costs, and the durable lessons.
- **Never carry a `NEXT = <step>` pointer.** Git and the task list hold it and
  update themselves; a hand-written pointer is the field guaranteed to rot.
- A plan's internal task numbers **collide with real GitHub issue numbers**.
  Never write a bare `#N` for an internal step.

### 6b. A step boundary must not cut across a signature's CALL SITES

Steps are usually drawn along conceptual lines — *change the signature*, then
*retire the guard*, then *write the contract*. That is how a plan reads best and
it is the wrong unit when the steps touch one interface: **the unit of work is
the call-site set, not the tidiness of the description.**

Before committing to a step decomposition that changes a signature, enumerate
every call site and ask which step each one lands in. A step that leaves any
call site speaking the old signature is not a step; it is half of one, and the
tree does not compile (or worse, compiles and fails at runtime) between them.

> `[M]` 2026-08-06, campaign P6. "Step 2 — change `evaluate(shape)` to
> `evaluate(space)`" and "Step 3 — retire the ERR-047 probe" read as two clean
> steps. They are one: the probe **is** one of the two call sites of
> `evaluate`, so landing step 2 alone leaves it calling `evaluate((N,))` on a
> source that expects a space. Discovered only because the pick-up honoured §7
> and re-read the call sites before designing; a session that trusted the
> decomposition would have found it by breaking the tree.

Corollary for the plan text: when steps ARE fused this way, say so where the
step is defined, not only in the commit. The next reader is planning against
the plan, not against your commit history.

## 7. Before resuming a plan, reconcile it against the tree

The plan is a snapshot; the tree is the fact. On pick-up:

1. `git log` + `git merge-base --is-ancestor <hash> HEAD` for every claimed hash;
2. read the **implementing class's first line** before designing to a phase's
   prose — a scope read from prose gets refuted by the realization;
3. re-check any "blocked by / not possible" claim, because the blocker may have
   been dissolved by a later phase of the same campaign. (P4 §5's "there is NO
   public entry point" was true when written and void two commits later; the
   banner saying so had to be added retroactively.)
