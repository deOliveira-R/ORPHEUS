# Plan authoring — what a plan owes its own future author

A campaign plan in `.claude/plans/` is not a note to the reader. It is a note to
**you, after your context is gone**, and it will be read by someone with your
authority and none of your memory. That reader cannot tell a considered decision
from a placeholder, an outcome from a guess at the means, or a measurement from
a plausible number — unless the plan marks the difference.

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

## 7. Before resuming a plan, reconcile it against the tree

The plan is a snapshot; the tree is the fact. On pick-up:

1. `git log` + `git merge-base --is-ancestor <hash> HEAD` for every claimed hash;
2. read the **implementing class's first line** before designing to a phase's
   prose — a scope read from prose gets refuted by the realization;
3. re-check any "blocked by / not possible" claim, because the blocker may have
   been dissolved by a later phase of the same campaign. (P4 §5's "there is NO
   public entry point" was true when written and void two commits later; the
   banner saying so had to be added retroactively.)
