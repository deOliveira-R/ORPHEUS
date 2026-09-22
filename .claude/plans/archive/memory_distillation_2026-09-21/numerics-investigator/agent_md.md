# AGENT.md promotion candidates — numerics-investigator, 2026-09-21

`AGENT.md` (316 lines) was read whole before proposing. It already carries, as standing
directives: the cascade and its step order; the reference-pillar trace (Step 1); the residual-
not-increment rule (Step 3, the whole L11 paragraph); the token-adjacency sweep (Step 4.5); the
NON-FLAT per-ordinate rule (Step 5, the whole L6 paragraph); the scaling table (Step 6); the
structural-independence bar before promoting agreement evidence (Step 7); the ERR-entry and
skill-maintenance close-out; and Rules 1–7. **Two of the three things I would have promoted are
already there** — which is why the digest's L6 and L11 are stubs.

Being conservative, as the standard asks: **two candidates, one strong and one marginal.** No
war story, no codename, no numbers in either.

---

## P1 (recommended) — a new standing directive, "Ask what KIND of question this is"

**Where.** A short section immediately before "## Diagnostic Cascade", since it governs whether
the cascade is the right instrument at all.

**Why identity-level.** It is applied on essentially every task and it changes the *first*
action, not a later check. The cascade answers "which component is wrong"; a large and growing
share of my dispatches are not that question — they are "is this rate claim real", "must this
accessor widen", "who owns this quantity", "can this statistic gate that contract", "what is
this operator's kernel". For each, the cheap correct instrument is a theorem, a spectrum or a
closed form, and re-running the solver is the expensive wrong one. Six separate investigations
lost time to starting with a measurement, and the correction was the same each time.

**Proposed text (verbatim, ~120 words).**

> ## Before the cascade: what KIND of question is this?
>
> The cascade answers *"which component is wrong?"*. Not every dispatch asks that. Name the
> question first, because several kinds have a cheaper and stronger instrument than a solver run:
>
> - a **RATE** question is a spectrum question — build the iteration matrix and eigen-solve it;
>   never re-time the solver.
> - a **CONTRACT / arity** question ("must this widen?") is a theorem question — ask what the
>   defining conditions COMMUTE with.
> - a **KERNEL or counting** question is usually a closed form — a law independent of a
>   parameter the operator contains is combinatorial; derive it, do not fit it.
> - an **OWNERSHIP** question is answered by measuring the increment's structure, not by what
>   the quantity is called.
> - a **"can statistic X gate contract Y"** question is one number, the transfer gain `|Δy|/X`,
>   measured before any threshold.
>
> Only when the question really is "which component is wrong" does the cascade start at Step 1.
> The instances, with their measurements, are agent memory `lessons.md` (spine M1).

**Leave behind in the digest.** Spine M1 keeps the instance list and gains a
`→ now in AGENT.md` pointer on application.

---

## P2 (marginal — the orchestrator's call) — a clause added to Rule 5

**Where.** `AGENT.md` § Rules, appended to Rule 5 ("Write runnable evidence").

**Why marginal.** It is identity-level by the standard's test (a standing diagnostic discipline,
applied on every probe I write), and `instrument-doctrine` X1 already states the general law, so
this would be a domain restatement — exactly what law 2 of this pass retires. I propose it only
because X1's examples are gates and censuses, while the thing that repeatedly fails here is an
**analytic** instrument (a symmetry-annihilated diagnostic, a knob decided upstream by rounding
noise, a linearity check on a matrix that is linear by construction, a monkeypatch onto a
`property`). If the orchestrator judges X1 sufficient, drop this and keep spine M5.

**Proposed text (one sentence appended to Rule 5).**

> Every probe carries its control in the same run: for an ANALYTIC instrument (a symmetry
> functional, a closed-form diagnostic, a knob sweep), feed it deliberate garbage in the varied
> slot and require it to move, because a symmetry the design introduced can annihilate the very
> functional that would judge it.

**Leave behind in the digest.** Spine M5 keeps the five instances.

---

## Considered and NOT promoted (with the reason)

- **"A curvilinear matvec needs a NON-FLAT per-ordinate reference"** — already in `AGENT.md`
  Step 5, verbatim, with a pointer already recorded in the old digest.
- **"Measure the residual, not the increment"** — already in `AGENT.md` Step 3, verbatim.
- **"Never skip a cascade step"** — already in `AGENT.md` § Diagnostic Cascade and in
  `probe-cascade` § Anti-patterns.
- **"Trust git for merge status, never a frozen memory claim"** — a real and recurring failure
  here (four stale status lines found in this pass), but it is `process-discipline` § "Trust git
  for merge status", an always-on rule I already load. Duplicating it in `AGENT.md` is exactly
  the failure law 2 of this pass exists to remove; the digest carries it as a two-line banner in
  `MEMORY.md` §2, which is where the stale claims lived.
- **"A degenerate fixture annihilates a bug class"** — `AGENT.md` Rules 2, 3 and 4 already carry
  three of its instances (refinement, homogeneous-exact, 1-group). A fourth abstraction layer on
  top of three concrete rules would add words and no behaviour; spine M4 is the right home.
