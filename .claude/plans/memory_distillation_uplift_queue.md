# Uplift queue — proposals from the 2026-08-03 fleet memory distillation

Six agents distilled their own memory (commits `89e81322`, `4ea700b8`, `db9585fa`,
`d97348b7`, `ce75650c`, `eec8a238`). Each returned promotion proposals. This file holds
the ones NOT applied on 2026-08-03, so they survive compaction.

**Applied that day** (so they are NOT in this queue): the `-n` correction to
`coding-standards.md` (`c5d7b7c2`); six additions to `coding-standards.md` /
`process-discipline.md` / `delegation.md`; the archivist AGENT.md baseline + `-n` fix.

---

## A. `vv-principles` — LANDED 2026-09-22 (the header below it is history)

The skill is generated from `docs/development/skills/vv-principles.md` and committed; the
"uncommitted state" this section was written against no longer exists. A1 landed as
anti-pattern #36 (2026-09-05); A5 was resolved in `coding-standards`; A2, A3, A4 and A6
landed 2026-09-22 as Mode 8 class (10), the "Log every caught bug" REASON half,
anti-pattern #37 and anti-pattern #17 rider (j). The original header, kept as the record
of why the section was blocked:

> BLOCKED, do not apply without checking the working tree.
> `.claude/skills/vv-principles/{SKILL.md,error_catalog.md}` carry irrecoverable
> uncommitted state and are forbidden to commit. Every item below is queued for a session
> where that is resolved. **Read the live file first — some may already be there.**

### ⚠ NUMBERING — #18 IS NOW TAKEN (updated 2026-08-03, later the same day)

**`18` is spoken for.** A THIRD anti-pattern #18 was written and LANDED that
evening: *never credit a mutation's reds as coverage when the mutation also
breaks a structural law the object obeys* — the over-powered-mutation rule, the
exact dual of #17. Both queued proposals below must therefore be renumbered
(**19** and **20**, in whatever order they land). Do not paste either as "#18".

That the same number was claimed three times by three agents in one day is
itself the argument for #308's generated-index direction: a hand-maintained
ordinal in a file many agents append to is a collision waiting to happen. Until
then, **re-read the live file's last number before writing a new entry** —
`grep -oE "^[0-9]+\. \*\*NEVER" SKILL.md | tail -1`.

The original two collisions, both still queued:

1. **`qa` and `test-architect` each proposed a DIFFERENT "anti-pattern #18".**
2. **Three agents proposed a DIFFERENT "eighth class" for Mode 8** (`qa`'s deselected
   catcher, `test-architect`'s no-op marker flip). Mode 8 already documents seven.

Renumber on landing; do not paste either verbatim.

### A1 — Mode 8, new class: the DESELECTED catcher (`qa`)
A gate marked `@pytest.mark.slow` (or any marker the run config excludes) cannot fire in
the gate that actually runs, so crediting it is the same false green as a stripped assert.
Review: check the catcher's markers against the invocation actually used, and confirm by
simulating the regression under that invocation.

> ⛔ **Fix the citation before landing.** `qa` justified this as "`-m "not slow"` is the
> fleet-wide canonical invocation per `.claude/rules/vv-testing.md`". **That file says the
> canonical invocation is `python -O -m pytest`**; `not slow` appears only inside an
> example. The 52-minute practical gate does use `not slow`, but that is recorded in the
> main agent's `reference_test_execution_env` memory, not in the rule. Cite it correctly
> or the rule ships a false reference — the exact defect class this anti-pattern is about.

### A2 — Mode 8, new class: the MARKER WHOSE FLIP IS A NO-OP (`test-architect`)
The existing class-4 defence ("prove the XPASS flip") assumes the flip-edit is sensitive
to the landing. It need not be: a strict xfail can be written so its own prescribed
flip-edit turns it into a character-for-character duplicate of the live flip-proof beside
it — green before AND after the phase, whose marker-deletion signals "phase complete"
while asserting nothing new. Discriminator: diff the xfail body against its own
flip-proof; if the documented edit makes them textually equal, the flip is ceremony.

### A3 — Mode 8 class 7 extended to the REASON, not just the fixture (`test-architect`)
Class 7 covers a `catches` marker whose fixture drifted. The same half-life hits a gate's
*justification*, and there the mutation run cannot see it because the gate stays correctly
green. When a phase falsifies a STRUCTURAL claim, rows asserting it can stay green on a
now-special-case fixture while the argument that made them meaningful is false. Review:
grep the claim's WORDS in `tests/`, not only its symbols.

### A4 — new anti-pattern: the reciprocity gate needs a one-sided partner (`qa`)
`⟨A.solve q, p⟩ = ⟨q, A.solve_transpose p⟩` pins the transpose RELATIONSHIP, not
correctness — satisfied by any genuine `(S, Sᵀ)` pair, so Mode-12 blind to a SYMMETRIC
regression dropping the same completion from both halves. Mutation battery is two-sided:
undo only the transpose half (must red) AND drop both halves (stays GREEN — that is the
finding). Only the one-sided `A ∘ A⁻¹ ≡ I` identity gate catches the symmetric half; the
two are non-redundant partners. (ERR-071.)

### A5 — the `legacy`-pin demotion (`test-architect` B3) — **RESOLVED ELSEWHERE**
Substance landed 2026-08-03 in `coding-standards.md` under "Retirement means test
migration", which is the right owner (it is a retirement consequence). If a
`vv-principles` cross-reference is wanted, make it a pointer, **not a second copy** — that
duplication is the exact disease this whole distillation was cleaning up.

### A6 — bit-identity teeth for a value-correct-by-coincidence twin (`elegance-enforcer`)
For a twin that is value-correct by coincidence, the teeth must be `array_equal` (0 ULP):
only bit-identity separates a leak from a genuine override. Currently parked in the
elegance-enforcer digest §A.

---

## B. AGENT.md promotions — held for review

Each is the owning agent's own identity file. All are proposals, none applied.

- **`elegance-enforcer` A1** — replace the §"Scope of Review" sentence that treats the
  dispatch brief as co-equal authority with "enumerate from a FRESH `git status` at review
  time; the brief's scope is a claim to verify". (This produced a retracted finding.)
- **`elegance-enforcer` A2** — leg-3 rider: a surprising ABSENCE is a tooling failure until
  re-verified with a differently-spelled grep. (Fired 3×, each one keystroke from a false
  MUST-FIX.)
- **`test-architect` A1** — new §0.7 "measure the premise before you gate it": run the
  proposed acceptance criterion as a probe first; trace the RUNTIME object not the type
  hint; never trust a count; `ls` the target path before delivering a pre-carve plan.
- **`test-architect` A2** — validate the mutation harness before believing any negative it
  reports; an all-blind verdict is a broken instrument until a positive control says
  otherwise.
- **`archivist` A2** — the retirement radius as a standing discipline. ⚠ **Now largely
  redundant** — the three-grep rule including the concept-paraphrase half landed in
  `coding-standards.md` on 2026-08-03. Land only the archivist-specific residue, as a
  pointer to the rule.
- **`archivist` A3** — never rename/delete an equation `:label:` that a
  `@pytest.mark.verifies(...)` targets; grep `orpheus/` and `tests/` first, and for a stale
  equation that IS a verifies-target, keep the label and rewrite only the body.
- **`cross-domain-attacker` A1** — new task-classification row: a naming/vocabulary
  adjudication is frame-detection work, not taste work (3 sightings).
- **`cross-domain-attacker` A2** — the backbone also says WHERE a foreign frame fires: a
  frame keyed to an operator's ALGEBRAIC SHAPE fires only on members whose shape matches.
- **`literature-researcher`** — proposed ZERO promotions, correctly: it checked and found
  its candidates already in AGENT.md/`delegation.md`, and declined to create a third copy.
  Optional low-priority sharpening: AGENT.md §6 does not state the phantom-citation
  MECHANISM (each forward reference makes the phantom look more established; the error
  compounds silently).

---

## C. `cross-domain-frames` — held

**Smell 17 — a `-> bool` predicate whose BODY builds the object.** A function named
`is_*` / `check_*` / `*_closure` whose body constructs an index map, permutation, matching,
partition, or certificate and returns `bool`. The capability is not missing; its WITNESS
is. Fix: widen the return type first; hand-rolled downstream re-implementations then delete
themselves. Distinct from Smell #16 shape 1 (that says "collapse two paths"; this says "one
path exists and is discarding its output").

Held because the agent honestly reported the evidence bar: **one** in-tree sighting (#326
`_orbit_closure`, with two downstream re-implementations) against Smell #16's seven. User
call whether Part C's stated bar is met.

---

## The finding that motivated all of this

The largest cut in every heavy agent was NOT compression — it was doctrine that had been
**uplifted into a shared skill and never retired from the private copy** (`qa` ~1165 lines,
`test-architect` ~460, `elegance-enforcer` 6 lessons, `cross-domain-attacker` 22 lines).

**So: when landing anything from this queue, delete the source lesson from the agent's
memory in the SAME change.** Retiring the original is part of the promotion. Otherwise
this file is just a recipe for regrowing what was cleaned up.

---

## Reconciled 2026-09-21 (the fleet-wide memory distillation, `.claude/plans/harness_context_budget.md`)

- **`test-architect` A1** — LANDED as `AGENT.md` §1.5 item 1 ("construct and measure first", the digest's meta-lesson M1), and its general form as `process-discipline` "Measure a brief's premise before arguing its scope".
- **`test-architect` A2** — RESOLVED by the rule: `instrument-doctrine` X1 (a positive control before any battery; an all-blind verdict is a broken instrument), on the always-on floor since 2026-08.
- **`archivist` A2, A3** — A2 stays redundant with `coding-standards`; A3's imperative is the archivist digest's §3 first bullet, not promoted (conditional on a label edit, not identity-level).
- **`cross-domain-attacker` A1, A2** — re-endorsed by the owner at four sightings (A1) and one more corollary (A2); still HELD for the user's review. Its new P1 (record the question a frame was refuted FOR) LANDED in `AGENT.md` and, as a check, in `process-discipline` "A refuted candidate".
- **Smell 17** — still HELD; the owner resubmits it as a two-shape smell (shape (a) the `-> bool` body that builds the object; shape (b) a correct predicate wired only to the advisory path). `[M]` 2026-09-21 on shape (b)'s sighting: an AST pass over `orpheus/` finds `admits_domain` at its definition only, zero callers of either kind, so the sighting is "unconsumed", stronger than "advisory only". The user's call on Part C's bar stands.
- **`elegance-enforcer` A1, A2** — untouched; that agent's digest (408 lines) was not in this pass.
- **New, surfaced and HELD** (from the main-memory census): a "Naming" section of `coding-standards` (six high-signal checks and the greppability law, `feedback_high_signal_names.md`, `feedback_naming_consistency_greppable.md`); the four lossy-return-type checks as a `coding-elegance` Pattern 4 corollary (`feedback_lossy_return_type_is_the_root_cause.md`).
- **`numerical-bug-signatures`, a new signature (HELD, owes an ERR entry first)** — the greedy `(Ellipsis, *idx)` spectator-axis index: bit-identical for the scalar-moment case, `IndexError` on a rectangular grid or a silent wrong value on a square grid with an asymmetric material map once a trailing axis is present; fix `cells = (slice(None), slice(None), *idx)`; blind test classes: every scalar-moment test. Founding case #276 A2, commit `0b3275d`, all four `MaterialXSField` moment-scatter verbs. The skill's add-a-signature protocol asks for the ERR entry first; none exists (`vv-principles` "Log every caught bug"), so the archivist owes it before the signature lands. Proposal: `.claude/plans/archive/memory_distillation_2026-09-21/numerics-investigator/uplift.md` U4.
- **`instrument-doctrine` X2, parameter-independence as a derivation hint (HELD)** — a fitted law measured independent of a parameter the object contains is combinatorial; derive it. One sighting (#344's kernel basis, 0.05 s closed form against a 23 s SVD). Proposal: the same file, U10.

## Landed 2026-09-22 (the user's ruling: "I accept all your recommendations on this implementation. Go ahead.")

Branch `chore/uplift-queue`, four commits, each item's home and founding case recorded in
`.claude/plans/harness_context_budget.md` at ⏸ COMPACTION POINT #12 and its "Uplift queue
landed" record.

- **§A** — LANDED: A2 as Mode 8 class (10); A3 as the REASON half of the "Log every caught
  bug" decay paragraph; A4 as anti-pattern #37 (evidence `## AP37 reciprocity partner`);
  A6 as #17 rider (j). A1 had landed as #36 (2026-09-05); A5 resolved in `coding-standards`.
- **§B** — LANDED: elegance-enforcer A1 (enumerate from a FRESH `git status`/`git diff`; the
  brief's scope is a claim); cross-domain-attacker A1 (the task-type row) and A2 (the
  backbone's WHERE clause), with "→ now in AGENT.md" pointers in its digest; the
  literature-researcher's optional §6 mechanism sentence. RESOLVED without an edit:
  elegance-enforcer A2 (the `code-search` positive-control check, which the Key agent
  loads). Already reconciled on 2026-09-21: test-architect A1/A2, archivist A2/A3.
- **§C** — LANDED: Smell #17, two shapes, in `cross-domain-frames` Part C, version history
  2026-09-22; the attacker's digest M1 item 6 points at it.
- **Naming section / lossy-return corollary** — LANDED in `coding-standards` "Naming" and
  `coding-elegance` Pattern 4; the three main-memory notes retired outside git.
- **Greedy-`Ellipsis` signature** — LANDED as Signature 11 after its ERR-087 entry and the
  verified catcher (`tests/transport/test_material_field.py::TestIndependentReference::
  test_moment_source`, `catches("ERR-087")`; the mutation reds exactly the four
  trailing-axis rows under `-O`); the numerics-investigator's L13 retired to a pointer.
- **Still HELD, on purpose:** the `instrument-doctrine` X2 parameter-independence heuristic
  (one sighting); archivist A3 (a lesson, not identity-level).

This file is now a record. Nothing in it is pending.
