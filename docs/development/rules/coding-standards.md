# Coding standards — the minimum-quality floor

Minimum standards every contributor (main agent and sub-agents) follows by default: the **floor**; `coding-elegance` is the ceiling (Cardinal Rules 1, 2). Each clause carries `check:` (the mechanical check), `tell:` (how the failure looks) and, where one exists, a `[case]` link into [the evidence page](../evidence/coding-standards.md). `[M]` measured, `[R]` reasoned.

## Clean before extending

Before adding a capability to a class/module, run a cleanup pass on that layer first: collapse double paths, move concepts to their native place, delete dead shims, fix twin sources of truth. The capability then lands as a **no-op extension through the one generic body**, not a third arm grafted onto debt.

- check: a plan proposing a capability extension inserts a **cleanup phase before** it; order findings into must-precede / independent-polish / explicit-WAIT; gate each cleanup substep bit-identical where possible.
- tell: the new arm needs a matching arm in the converter AND the constructor AND the gate. [case](../evidence/coding-standards.md#2026-06-11-from-axes-roundtrip)

## Type vs property — before minting a type

Mint a **type** iff (a) the concept has **two or more non-isomorphic realizations** AND (b) a **non-identity morphism** is actually applied to it. Otherwise it is a **property**: a field or flag on an existing type. One realization plus an identity change-of-basis makes a "type" theatrics: a conversion seam and no illegal state made unrepresentable (a single-basis spatial moment is a `property`, not a `SpatialOrder` type).

- check: count the realizations and name the morphism; identity-only change of basis means property.
- **An axis that changes the ARITHMETIC INTERFACE cannot be a phantom type parameter.** `Generic[Tag]` is erased at runtime and does not specialize dunders, so every instantiation shares ONE `__add__`; a torsor `A×V→A` that must forbid `A×A` and a vector `V×V→V` cannot share a body. Arithmetic or shape changes: a class; neither: a phantom parameter is allowed.
- tell: an implementation that "passes" only by branching on a stored tag at runtime is stringly-typed dispatch; `replace(obj, tag=Other)` type-checks and walks through the gate the type was minted to be.

## A bare `assert` in `orpheus/` is not a contract — the canonical runner strips it

`python -O -m pytest` is canonical; `-O` sets `__debug__ = False` and removes every `assert` at compile time. A contract written as a bare `assert` **does not run in the suite that matters**, and the code ships accepting the input the assert refuses.

- check: `grep -n "^\s*assert " orpheus/` and sort the hits. **Type-narrowing** (`assert x is not None` for pyright) may stay: it was never the guard. A **numerical / domain / admission contract** (a tolerance, an invariant, a shape law) **MUST be a real `raise`**, modelled on the nearest admission guard (`_assert_alpha_dome_closes`) so the vocabulary stays greppable.
- **Prove it, don't argue it.** Run the guard's own arithmetic on a deliberately-bad input under `python` and under `python -O`. tell: it returns instead of raising, so the contract is inert. [case](../evidence/coding-standards.md#2026-08-12-alpha-dome-assert)
- **Converting one is a retirement** (the audit below); tests pin the **shortest distinctive fragment of the OLD assert's message**: grep that, not your new wording.
- Cardinal Rule 2 first, then the guard: the founding recursion had **three copies**, which is why the contract could live on one arm only.

## A guard is elegance debt — tag it, and name what retires it

A runtime guard (`require_member`, `admit_composite`, a typed refusal on an alien carrier) is a **signal that the architecture failed** to make the mistake unspellable (`coding-elegance` Pattern 4). A legitimate protection *today*, not the target state: **the ultimate state is to not need the guard** ([R] user ruling, 2026-09-07: [case](../evidence/coding-standards.md#2026-09-07-guard-is-debt)).

- **Every guard that lands carries a greppable marker in its docstring:** the token **`ELEGANCE-DEBT[guard]`**, the issue number, and ONE sentence naming the structural change that makes the guarded mistake unspellable (e.g. *"retires when B is bound on its own trace end"*). check: `grep -rn "ELEGANCE-DEBT" orpheus/` is the debt ledger.
- The issue is filed **with the carve that lands the guard**, never before (a guard without its retirement plan is an unpriced debt; a plan without its guard is a promise). The step landing the structural change deletes guard AND tag in the same commit, and the mutation battery must show the mistake is now unspellable, not merely refused.
- tell: the guard's docstring justifies the *check* rather than naming the *shape that would make the check unnecessary*.

## Retire as you go — the audit as a numbered checklist

Superseded code is noise that invites extending the wrong path. **Retirement is a first-class deliverable.** Every refactor introducing a better pattern MUST retire its predecessor; shims live **one merge cycle only**; never keep backward-compat unless the user explicitly authorizes it. The audit is its own numbered substep, with a `file:line` retirement list.

**A. The three searches — run all three, then retire.**

1. **Graph callers.** check: `nexus impact` / `callers`; necessary, NOT sufficient. tell: `callers()==0` but live via a `cached_property`; class-name *bypass* consumers; direct constructors of a guarded type.
2. **Text-grep the symbol across code, tests, AND `docs/`.** An unresolved Python-domain xref (`:func:`/`:class:`/`:meth:`/`:mod:`) renders as plain text with **no `-W` warning, and `-n` does NOT save you**: Sphinx nitpicks only what it RENDERS, and only ~45 modules are `automodule`'d; nothing under `tests/` renders. check: is the module rendered at all? If not, **grep is the only gate**. tell: an unchanged warning count (it proves nothing). [M] 2026-07-15: `:doc:`/`:ref:` DO warn; raw path strings never do (grep a segment-assembled path by its **last segment**). [case](../evidence/coding-standards.md#2026-07-15-sphinx-severity)
3. **Direct constructors** of any guarded type. check: grep `T(`; a guard-at-source change reaches every `T(...)` caller, not just the factory path. tell: the factory path is guarded while `T(...)` callers are not.

**B. Surfaces a symbol grep cannot reach.**

4. **`dead_references` is the only instrument that reads the DOCSTRING surface.** check: run it before calling any retirement or re-home done, and after the fix. tell: a validated AST import audit returned **0** while 5 dead targets sat in docstrings; the filter pointed at the wrong SURFACE. [case](../evidence/coding-standards.md#2026-08-28-p44-dead-references)
5. **A math symbol has THREE spellings, and the NUMBER is a fourth**: ASCII identifier (`tau_raw`), Unicode prose (`τ_raw`), LaTeX role body (`\tau_{\rm raw}`). check: grep all three, then the number the claim carries (`tfrac15|tfrac45`). tell: a two-spelling sweep reports clean while a page still asserts the old figure. [case](../evidence/coding-standards.md#2026-08-11-tau-raw-spellings)
6. **A name inside a STRING.** `\.name\b` and `name\s*[:=]` miss `getattr(x,"name",None)`, `hasattr`, `setattr`, `__getattr__` keys. check: `grep -rnE "['\"]<symbol>['\"]"` on every retired name, then read what each default MEANS. tell: nothing raises; the call returns the DEFAULT and every branch keyed on it flips. [case](../evidence/coding-standards.md#2026-08-26-curvature-getattr)
7. **Grep the CONCEPT, not only the symbol.** A field is documented by NAME (greppable) and by PARAPHRASE (not): a column headed "Sweep-cycle flag" carries per-law values with no symbol in any cell. check: grep the hyphen/space variants too, then **triage every hit by MEANING**. tell: 7 exact hits beside 17 missed cells; 11 of 17 flagged pages false positives. [case](../evidence/coding-standards.md#2026-08-11-tau-raw-spellings)
8. **The paragraph that JUSTIFIED the retired thing inherits its wrongness.** check: re-verify that prose against the replacement. tell: the dead name deleted, its argument left standing.
9. **A LABELLED EQUATION is an API.** Correcting the prose around it does not correct it; Sphinx is silent because the `:eq:` reference RESOLVES. check: grep ``:eq:`<label>` ``, read every citer, fix the equation itself; confirm the label's `verifies()` marker still means what it says. tell: the page condemns the claim in prose while its Key Facts card still carries it. [case](../evidence/coding-standards.md#2026-08-13-pole-mm-recurrence)
10. **A CROSS-REFERENCE is a load-bearing dependency.** check: for every `:ref:`/`:doc:`/`:eq:` you ADD or LEAN ON while correcting, read the target and ask whether it still says the old thing. tell: the repair reads as *more* rigorous for its citation, and the cited section holds the half-claim it refutes. [case](../evidence/coding-standards.md#2026-08-14-q6c-cross-reference)
11. **A MESSAGE STRING is an API the moment a test pins it, and tests pin SUBSTRINGS.** check: grep the **shortest distinctive fragment**, never your own longer wording. tell: the audit reports only the definition site, then gates go red in the wide run. Prefer KEEPING the established vocabulary: it is load-bearing provenance. [case](../evidence/coding-standards.md#2026-08-06-message-fragment)

**C. Migration — a retirement that only deletes loses coverage.**

12. **Test migration.** Behavioral test (correctness contract): rewire to the successor. API-smoke test: delete. Characterization test: keep under `tests/<module>/characterization/`. check: `grep -rn "<symbol>" tests/`. tell: a delete-only diff.
13. **Marker migration.** A retired test takes its `catches(...)` / `verifies(...)` markers with it; re-tag the successor asserting the same invariant. check: grep the catalog and the `tests/_harness` registry with the code grep. tell: the catalog names a dead test class; the audit reads "MISSING".

**D. What a retirement silently does to the surviving gates.**

14. **DEMOTION.** Re-pointing a comparison target at the successor turns a two-implementation gate into a value compared with itself through a wrapper: green forever, invisible in review because BOTH sides are real production calls. check: *are the two sides still INDEPENDENTLY produced?* One-line mutation: replace the SUT's body with garbage; still green means never a pin. tell: a local still called `legacy` beside a brand-new API. Do not delete the gate; re-scope every doc crediting it and name the surviving pin, **then check that pin is real**: the cited snapshots may have been re-captured BY the same carve. [case](../evidence/coding-standards.md#2026-08-03-recaptured-snapshots)
15. **PROMOTION, the unhunted mirror.** A retirement can make a surviving gate STRONGER; nothing fails, so its docstring keeps advertising the weaker claim, which is how a real gate gets deleted as redundant. check: for every gate that survived UNTOUCHED, re-derive what its assertion now compares. tell: a docstring saying "tautological" or "by construction" about an implementation that may no longer exist. [case](../evidence/coding-standards.md#2026-08-17-loss-action-promotion)
16. **SINGLE-SOURCING a duplicate demotes every gate that compared its copies; the demotion is CORRECT, never back it out.** Neither tell above fires; both sides now derive from one constant. check: in order, (a) what input could still make the two sides differ? "none" means tautological, a *design-time* question, not a mutation one; (b) hunt an EXTERNAL hand-written pin; none means you owe a replacement; (c) keep the gate for what it still tests and say so **in its own docstring**. tell: an authoritative name on a comparison that cannot fail. [case](../evidence/coding-standards.md#2026-08-09-capability-rows)
17. **Check that the FIRST pin you reach for MOVES under the OLD value; it is the one most likely BLIND.** Applies to every re-baseline: the nearby candidate usually rests on the fixture that hides the change, and green is compatible with *loaded* and with *blind* (`vv-principles` #19). check: one in-process mutation at the old value. tell: an all-reflective flat medium reads the quadrature only through `sum(w)`, bit-exact before and after, while `μ₁` moved 14 %. [case](../evidence/coding-standards.md#2026-08-12-blind-pin-twice)
18. **ENFORCEMENT: retiring a duplicate makes whatever KEPT THE COPIES EQUAL load-bearing, usually a production `raise` with no test.** Its code does not change, so nothing prompts a re-look. check: in order, *what makes the copy provably redundant?* (the mechanism) and *does that mechanism have a witness?* (grep the shortest distinctive fragment of its message). None: write one **in the same commit** (the retirement created the exposure). tell: the grep returns only the production `raise` lines. The tests that asserted the retired copy (a stored literal read back) are the right place for the witness. [case](../evidence/coding-standards.md#2026-08-27-p41a-guard-witness)

**E. Retiring the retirement's own residue.**

19. **A retirement's past-tense NOTE is a confidence trap; sort surviving references BY TENSE.** check: grep the deleted name tree-wide. Past tense ("was retired") is history and STAYS; a **present-tense claim** ("provides a shortcut") or an **imperative** ("Apply `@verify.lN(...)`") is a MUST-FIX: a follower hits `ImportError`; a maintainer re-adds the symbol. tell: the offenders are PRE-EXISTING lines outside the diff, ~50 lines from the batch's own note.
20. **Routing a call site through the ALGEBRA raises "has the verb" to "IS the type".** `f.apply(g.apply(x))` to `(f @ g).apply(x)` is arithmetic-neutral, not a pure re-spelling: `@` needs `__matmul__`. Production objects come from a factory, so the breakage surfaces in a *test* double. check: make the surrogate honour the contract; do NOT add a runtime guard for a case the type system now covers. tell: `TypeError: unsupported operand type(s) for @` in a duck-typed double. [case](../evidence/coding-standards.md#2026-08-06-matmul-operand)
21. **A retirement onto a SHARED helper moves the raise's provenance one frame out.** check: re-point any frame-asserting gate at the helper AND assert the CALLER frame (or an owner tag in the message). tell: "this composite refused" has widened to "something refused".
22. **Mass-deletes are retirements too.** check: (a) grep the module/script **NAME**, not only its symbols (a worker import inside a `textwrap.dedent` string is invisible to the call graph); (b) a diagnostic consumed by a tracked test is **production infrastructure**, never debris; (c) a consumer still green after the delete: check WHERE its import resolves. tell: an untracked scratch shadow or stale `__pycache__` serves the import until it evaporates. [case](../evidence/coding-standards.md#2026-07-13-mass-delete-shadow)

**F. A rename is a retirement of the old spelling.**

23. **A LENGTH-changing rename owes the `.rst` corpus a section-underline scan.** A heading carrying the renamed word (`SNMesh` to `SNProblem`, +3) outgrows its underline; `sphinx -W` errors. check: one pass over `docs/**/*.rst` comparing each heading's code-point length to its underline's after ANY length-changing rename. tell: the mechanical pass ships green; the build finds the red an hour later. [case](../evidence/coding-standards.md#2026-09-18-underline-scan)

**G. The retirement's own scope and order.**

24. **Retire leaves first, in one commit.** Symbols that exist only to support the retirement target are part of the retirement: before the plan, write a dependency-audit table (rows the symbols, columns the three surfaces above) to `.claude/plans/<plan>_dependency_audit.md`; retire in its order — leaves first (zero production callers), then the helpers only the retired set called, then the top-level symbol — in ONE commit whose body references the table. tell: an orphaned helper surviving the delete; a retirement spread over commits with green in between. [case](../evidence/lessons.md#l20-retirement-dependency-audit)
25. **A move's residues.** Read every tool that WRITES into the moved tree (a generator, a hook) before declaring the move complete; a path constant gets ONE home that every consumer, tests included, imports; archaeology on a moved tree is a per-FILE judgement (what does this text DO in this file), never a per-directory one. [case](../evidence/lessons.md#l34-path-segments-grep)

## The mirror — landing a deferred capability stales its DEFERRAL CONTRACT

Flipping a case from *deferred* to *implemented* makes every docstring saying "raises / deferred / not yet supported" for it present-tense-FALSE; same three-search blast radius. check: grep the case's NAME and prose forms across the package, then discriminate **by arm** (a matvec-transpose does not un-defer the transpose-SOLVE). tell: the prose ledger is rewritten, the machine-facing contracts are not: the Protocol stub, the BASE class's default docstring, the sibling class, operators in untouched files. A campaign-closing `Closes #NN` contradicts a status the tree still tags deferred.

## Exception — keep a relinquished *fuller view* as a verification oracle

Retirement targets a *superseded* predecessor, never a **fuller view an optimization relinquished** (a full field behind a rolling window; a dense operator behind a factored form): that is the verification pathway pinning the optimized path's reference. The decision is EXPLICIT (a wired, exercised oracle, or retired); the oracle is not production-reachable, shares the optimized path's kernel, and is pinned by a permanent end-to-end `optimized == oracle` gate. check: keep only a genuine STRUCTURAL reference; the same math procedurally rearranged AND already verified by an independent oracle (MMS, closed form) is redundancy: retire it (`vv-principles` #7; [L11](../evidence/lessons.md#l11-structural-independence)). tell: a fuller view orphaned but undeleted. Worked: `_sweep_2d_full_field` kept as the rolling-window oracle; the per-ℓ scattering kernel retired.
