# AGENT.md promotion candidates

**One candidate, plus one correction to an existing directive.** Read the
current `.claude/agents/qa/AGENT.md` first: its Enforcement list already carries
make-it-RED (#11), the field-role contract (#10), the multi-group/heterogeneous
demand (#4, #5) and the V&V-level classification (#1, #2), which is most of what
a promotion would otherwise propose. Be conservative — the rest of the digest is
failure→correction material and stays a lesson.

## P1 — promote: the instrument that reads clean is the one to distrust first

**Where:** a new Enforcement item #12, after #11 (which it generalises from
gates to every instrument).

> 12. **A clean reading is a claim about the INSTRUMENT before it is a claim
>     about the tree.** Whenever a census, a filter, a fingerprint, a counter or
>     a fidelity diff returns zero, name what it could not have seen and show one
>     known member it DID find. Three shapes recur here and each reads as good
>     news: a filter that dropped its whole input (an unsplit `$VAR`, a path
>     filter, a two-stage net missing a spelling); a detector normalised for
>     robustness, hence blind to the change class it is now asked about; and an
>     identity-keyed diff over text, which cannot see a check that moved away
>     from the imperative it belongs to. State the population and the instrument
>     with every zero you publish (`instrument-doctrine` X1, X2).

**Why identity-level rather than a lesson.** It is applied on essentially every
task, not in one situation: it is what I do with the SECOND half of every
review — the half where a measurement came back clean. It governs digest rules
in five different sections (A10, A11, A24, E19, H9, H13, H14, H15), which is the
signature of an operating principle rather than an instance. Its cost in
AGENT.md is six lines, loaded fresh per dispatch, and it sharpens #11 from
"make the gate red" to "make every instrument prove it can move".

**Leave in the digest:** each instance, with `→ now in AGENT.md #12` appended to
the §A header line, so the mechanics stay findable and the principle stays
unrepeated.

## P2 — correct an existing directive (not a promotion)

Enforcement **#9** reads *"Check `verification_coverage` — every equation should
have status `verified`"*. That is refuted by my own measurement: `verified` is
set iff `len(tests) > 0` with no confidence floor, all `tests` edges are
test→equation (there is no test→code edge), and `[M]` 351 of 692 "verified"
equations have no declared test (L-070; digest E7). As written, the directive
instructs me to chase a number whose weakest admissible evidence is a shared
name token.

**Proposed replacement:**

> 9. **Read `verification_coverage` for its PREDICATE, not its verdict.**
>    `verified` means only "some test edge exists", and most `implements` edges
>    are name-token guesses; ask of any status what predicate sets it and what
>    its weakest admissible evidence is, and adjudicate with a coverage capture
>    (`run=`) or a mutation (`nexus-verification`, its three ⛔ blocks).

This is a correction of a present-tense-false sentence in an always-loaded file
(Cardinal Rule 3), not a distillation product; flagged here because this pass is
the first time the two texts were read side by side.
