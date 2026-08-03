# Uplift queue — proposals from the 2026-08-03 fleet memory distillation

Six agents distilled their own memory (commits `89e81322`, `4ea700b8`, `db9585fa`,
`d97348b7`, `ce75650c`, `eec8a238`). Each returned promotion proposals. This file holds
the ones NOT applied on 2026-08-03, so they survive compaction.

**Applied that day** (so they are NOT in this queue): the `-n` correction to
`coding-standards.md` (`c5d7b7c2`); six additions to `coding-standards.md` /
`process-discipline.md` / `delegation.md`; the archivist AGENT.md baseline + `-n` fix.

---

## A. `vv-principles` — BLOCKED, do not apply without checking the working tree

`.claude/skills/vv-principles/{SKILL.md,error_catalog.md}` carry irrecoverable
uncommitted state and are forbidden to commit. Every item below is queued for a session
where that is resolved. **Read the live file first — some may already be there.**

### ⚠ TWO NUMBERING COLLISIONS — reconcile before writing anything

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
