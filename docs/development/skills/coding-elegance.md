---
name: coding-elegance
description: PROACTIVELY load when writing or reviewing any production code, designing an API surface, choosing between abstractions, refactoring, or evaluating whether an implementation reads like the math/domain it claims to encode. This skill codifies what "elegance in coding" means — patterns to invoke, anti-patterns to flag, the prevention-by-construction argument for why elegant code has fewer bugs, and the elegance checklist used at code-write time. Preloaded by all sub-agents that produce code (method-implementer, numerics-investigator, qa, test-architect) and by the main agent when orchestrating implementation.
harness:
  kind: skill
  budget_tokens: 9000
---

# Coding Elegance core

The argument, the worked ORPHEUS cases
and the long code contrasts live on the
[evidence page](../evidence/coding-elegance.md), with its
[examples 1-5](../evidence/coding-elegance.md#example-1-four-operator-unification);
open a `[case]` link only when you need *why* or *how it was found*.

> "Code elegance literally prevents bugs by construction." Not aesthetics:
> elegance makes domain-illegal states **unrepresentable**. Bugs are
> notation/domain mismatches; eliminate the mismatch and the bug class becomes
> unspellable.
> [why](../evidence/coding-elegance.md#why-this-skill-exists);
> [the argument](../evidence/coding-elegance.md#bug-prevention-argument).

---

## One-line summary

**Code should be a notation for thought; the notation must match the domain.**
Every pattern here is a corollary.

**The acceptance test:** read the code aloud. Does it sound like the math, the
workflow, the protocol it implements? If you hear procedural for-loops,
special-case branches, stringly-typed dispatch, or "and now we multiply by
4π/3", the notation lost its grip on the domain.

---

## The two questions to ask BEFORE writing code

1. **What is the algebra of the domain?** Its primitives, compositions,
   identities, inverses. Write them on paper. SN transport is
   `(L + C − S − F/k) ψ = q`; HTTP is `Request → Response`; a parser is
   `String → AST | ParseError`.
2. **What syntax in the host language captures that algebra?** Operator
   overloading, sum types, Protocols, composition combinators, method chaining.

**The answers to these two questions ARE the architecture.** The implementation
is whatever syntax mechanics the language requires to spell them. Skip them and
you write code that is *about* the algebra rather than code that *is* the
algebra — the procedural-transcription failure mode, of which every twin-path
bug is a downstream consequence.

Beck's four rules of simple design, in priority order — **pass the higher rule
before optimising the lower one**: passes all tests; reveals intention; no
duplication; fewest elements. Clarify first, then compress. Here the principle
has a sharper form: **the code is the math, character for character**.
[lineage](../evidence/coding-elegance.md#first-principle-code-as-notation)

---

## The Master Standard

The user's verbatim acceptance criterion:

> "The implementation should read like math under the dunder methods
> implemented with operator expressiveness."

| Domain statement | Elegant code |
|---|---|
| Fixed-source `(L + C − S − F/k) ψ = q` | `(L + C - S - F/k) @ psi = q` |
| Adjoint `(L + C − S)^† ψ^† = R` | `A_loss.H.solve(R.as_source())` |
| k-eigenvalue `K = (L + C − S)^{-1} F` | `K = A_loss.inverse() @ F` |
| Reaction rate `r = ⟨Σ_a, ψ⟩` | `r = ReactionRateFunctional(sigma_a) @ psi` |
| Sweep `ψ = (L + C)^{-1} q` | `psi = (L + C).solve(q)` |

When the right column reads aloud the same as the left, the standard is met.
When it doesn't, the implementation is one level of indirection away from the
math — and that gap is the bug habitat this skill exists to close.

---
## Seven elegance patterns

Each makes a class of bugs unspellable, and each is backed by ≥3 historical
ORPHEUS bugs (the per-bug audit is the evidence memo under
[pointers](../evidence/coding-elegance.md#pointers)).

**1 — Match the algebra of the domain via dunder methods.** Implement the
domain's operators as dunders on the domain types (`__add__`, `__matmul__`,
`inverse()`, `.H`). **Trigger:** the moment you would write a *function name* for
an operation the domain writes as a *symbol*. **Domain example:**
`A_loss = L + C - S`; `psi_dagger = A_loss.H.solve(response)` — `.H` propagates
through `OperatorSum` to leaves, so the adjoint solver is six lines and **zero
new code beyond leaf `.H`**. **Counter:** `solve_adjoint(L, C, S, q)` needs a new
function the moment F appears; an unbounded explosion.
[case](../evidence/coding-elegance.md#pattern-1-dunder-algebra)

**2 — Single source of truth (composition over duplication).** Every concept
appears in exactly one place; two pieces of code computing the same mathematical
quantity are a bug in waiting (X4 decides "the same": α-normalised ASTs, never names). **Trigger:** the moment you reach for copy-paste,
OR write a parallel implementation "because the layout differs". **Signature:**
"I just had to apply this fix in two places." **Domain example:** the sweep and
the matvec both compose ONE `SNCellOperator` over `iter_cells_by_direction(±1)`,
so ERR-026 manifestation #7 (O(h) drift between them) dissolves by construction.
**Allowed exception:** a *same-role* broadcast injection across storages
(`A ⊂ B`, same physical role) may stay implicit — the embedding carries no
convention to drift. Cross-*role* combinations stay explicitly named.
[case](../evidence/coding-elegance.md#pattern-2-single-source)

**3 — Named intermediates with domain semantics.** Every value crossing a
function boundary, living across iterations, or appearing in a return type gets
a domain name. **The load-bearing physics statement: an unnamed quantity is
evidence that the physics is wrong** — every unit-bearing combination has a
literature name, and the only unnamed values are (a) an identity lost at a
boundary, (b) units the author never tracked, or (c) a wrong formula. Dimensionless
intermediates are the exception that proves the rule: they are constructed
DELIBERATELY and are themselves named (`k_eff`, `τ`, `c = Σ_s/Σ_t`, `Re`).
**Trigger:** you are about to write `np.sum(...)`, `.reduce(...)`, `for ... +=`,
or any inline multi-term arithmetic. **Domain example:**
`np.sum(nu_sigma_f * phi * V[:, None])` becomes
`compute_group_production_rate(phi).sum()` — the intermediate IS the per-group
fission production rate `[1/s]`, and its ratio to the absorption rate IS
`k_eff`. **Review rule:** a line that
cannot answer "what is this quantity, in physical units?" is a code-review
blocker: the physics is opaque, not merely the code ugly.
[case](../evidence/coding-elegance.md#pattern-3-named-intermediates)

**4 — Make illegal states unrepresentable.** Encode invariants in the type
system so violating values cannot be constructed ("parse, don't validate").
**Trigger:** you are about to write a runtime check whose rationale is "this
can't happen in valid input". **Domain example:** `L.solve(q)` raises
`MissingCapability` — inverting streaming without collision is not a thing that
can be typed; `(L + C).solve(q)` is the only spelling of the sweep. Compare
`streaming_op.solve(q, include_collision=False)`: the invalid combination is
expressible and the type system cannot help. **Other
instances:** capability `frozenset`s on `LinearOperator`; `BoundaryRealizer` as a
constructor parameter of L (so "sweep without BC" is not a value);
`DiscreteOrdinatesPhaseSpace` (one object, not three swappable args).
**Corollary (Pattern 4 ∩ 2):** an operation producing a NEW instance of the same
frozen type routes through `dataclasses.replace(...)`, which re-runs
`__post_init__` so the invariant re-fires for free — a hand-written dunder that
*restates* the invariant is a Pattern-2 duplicate of the law, and goes stale.
[case](../evidence/coding-elegance.md#pattern-4-illegal-states)

**5 — Build the right primitive, not the right product.** Decompose a complex
behaviour into small composable primitives; the product is their composition.
**Why:** N primitives express N! products; a product must be extended once per
requirement. **Trigger:**
you are about to write `solve_<specific_problem>` and can already see its future
siblings. **Domain example:** `SourceIteration(L, C, S, F, q_ext).solve()` serves
fixed-source (`F=None`), k-eigenvalue (Arnoldi outer), α-eigenvalue (substitute
the pencil), and adjoint (`.H` at each leaf) — one primitive, many products.
**Counter:** the `_solve_<problem>_<method>` family, an unbounded explosion.
[case](../evidence/coding-elegance.md#pattern-5-primitive-not-product)

**6 — Defer abstraction until you have evidence.** Do not unify until ≥2 working
concrete instances exist (the project's amendment of the rule of three, for its
higher correctness stakes).
**Trigger:** you reach for `Protocol` / `ABC` / a generic and the second concrete
consumer does not yet exist. **Domain example:** `LinearOperator` was extracted
from `SNStreamingOperator` + `SourceIteration` + the Peierls Nyström kernel,
three instances. Four checks, each catching a different failure:

- **Concept-count test (forward).** Count the *concepts* the codebase carries
  before and after. The right abstraction **shrinks** the count; the wrong one
  adds a generic layer on top of the instances it failed to unify. If the count
  does not go down, keep the instances.
- **The don't-over-defer guardrail.** "≥2 instances" blocks *speculative*
  abstraction, not a **primitive** whose benefit is already established (an
  illegal state made unrepresentable; methods the alternative type cannot carry;
  independent expert frames agreeing). Build those at the first consumer.
- **The check that guardrail owes** [REFUTED 2026-08-19: its own worked
  example, `FluxDisplacement`, was retired ten weeks after being cited here].
  Spell out WHICH illegal state, then grep for a production path that
  legitimately produces it. If one exists the benefit is inverted, not
  established: the type will refuse correct output. The premise was never
  checked against the *producers*.
  [history](../evidence/coding-elegance.md#ap18-reversal-history)
- **Self-conceding-docstring trim (retirement-side mirror, grep-decidable).**
  Trim on sight when both hold: **zero production consumers** (callers/importers
  return only tests or nothing) AND a docstring conceding its own irrelevance
  ("not used in production", "diagnostic only"). A
  *fuller-view oracle* with a permanent equivalence test consuming it is not
  zero-consumer ([coding-standards](../rules/coding-standards.md), the
  fuller-view-oracle exception).

[case](../evidence/coding-elegance.md#pattern-6-defer-abstraction)

**7 — Normalise at the definition site, not at every consumer.** Convention-
dependent values — signs, normalisations, weight sums, energy-grid orderings,
axis conventions — are fixed at the ONE place the value is defined. The
convention is a property of the PRODUCER; N consumers means N chances to drift,
and a *future* consumer gets the bug for free. **Trigger:** you are about to
depend on a sign / normalisation / ordering / axis convention, or you see the
same convention applied at multiple call sites (X4's tell). **Domain example:**
`lethargy_bin_widths(eg)` returning `np.abs(np.log(eg[1:]/eg[:-1]))` at the
definition site — ERR-022 was a signed `Δu` re-derived at three consumers, one of
which met a descending grid. Eliminates ERR-004/008/014/018/022/025/031 by
construction.
[case](../evidence/coding-elegance.md#pattern-7-definition-site)

**Convention crosswalk, the mechanical form of Pattern 7.** Any carve crossing
subsystem boundaries writes this table to `.claude/plans/<carve>_crosswalk.md`
BEFORE any code; the crosswalk IS the architecture. One row per subsystem,
`| Subsystem | Input convention | Internal | Output |`, plus a **Bridge** row
naming which way and which transform. The Bridge row is where the pattern
demands action: move the bridge to the producer, or record the load-bearing
reason it stays at the consumer. Seven axes: per-ordinate vs iso scalar
(`/sum_w`); `/W` normalisation; μ sign; packed vs typed layout; normal vs
adjoint (`.H` propagation, `apply_transpose`); signed vs unsigned lethargy;
group ordering. [M] R-1 Step 4 session 1: ~3x debug time (three convention
bugs, ≥1 h each) against a 15-minute table ([L17](../evidence/lessons.md#l17-convention-crosswalk-first), [L18](../evidence/lessons.md#l18-pattern-seven-producer)).
[case](../evidence/coding-elegance.md#convention-crosswalk-case)

---
## A repeated conditional is a missing type

Discriminate once, at the boundary. This is the unifying lens behind
anti-patterns #3, #4, #7 and Patterns 1 and 4.

**The rule.** The smell is not the *number* of branches — it is **repeated,
interior, tag-based discrimination**: the same distinction (`if geometry ==
"spherical"`, `inner == "krylov"`) re-asked at many sites deep
in the call tree. Each is a **type that was never made**. Resolve it ONCE, at
the boundary, into a value/type; the core then runs decision-free because the
object it holds already encodes the choice. The domain does not say "if
spherical then...", it says "the spherical operator".

**The check — essential vs accidental.** *Does adding the next case force me to
edit existing branches?* **YES means accidental** — a missing type, an
Open–Closed violation (a new geometry must not require editing the sweep AND the
matvec AND the BC AND the source); resolve to a type, dispatch once. **NO means
essential** — a single, local, one-off branch (guard clause, early return,
exhaustive `match` over a closed sum type) is the honest expression of a genuine
split; keep it. **The goal is not branchless code:** each distinction decided
*once*, as *early* as possible, *recoverable as a type*.

**Why.** Repeated discrimination is a divergence habitat (fix one site, its
twins survive). One typed dispatch also buys **exhaustiveness**: a sum-type
`match` flags the *missing* case; a scattered `if/elif` with no `else` silently
does nothing.

**Know the axis before reaching for polymorphism (the expression problem).**
Cases grow, operations stable: **polymorphism / Protocol family**. Operations
grow, cases stable: **closed sum type + exhaustive `match`** (yes, branches). In
ORPHEUS geometries are stable while operations (sweep, matvec, adjoint, DSA,
...) keep growing, so polymorphism is usually right *here*, for that reason,
not by dogma.

**Positive forms.** (a) *Resolve at the boundary*: decide at construction
(`StreamingOperator` already IS spherical), pass the typed object down, the core
never re-asks. (b) *Data-driven
dispatch*: a `dict`/registry mapping tag → behaviour, so the branch becomes one
extensible data structure.

**Trigger.** You are about to write `if <tag> == ...` / `elif kind == ...`, or
add a `mode=` / `is_x=` parameter, OR you notice the same conditional already
living elsewhere. STOP: is an existing type already carrying this distinction
(extend it), or is this a genuine one-off local split (keep the branch)?

**Tell.** McCabe cyclomatic complexity is a *tripwire*, never a target: it is
Goodhart-prone (relocating a branch shrinks the count) and essential-blind.
Use it to FIND candidates; use the essential/accidental check to JUDGE them.

---
## Anti-patterns to flag

All twenty, each a redirect: **NEVER** X — **instead** Y. Flag on sight in
review. 1–12 are construction-level discipline; 13–17 the densest empirical bug
clusters from the 47-entry catalog; #3, #4, #7 are instances of *a repeated
conditional is a missing type*. Each `[case]` is the item's full statement with
its ERR cross-cites.

1. **NEVER** two implementations of one mathematical quantity — **instead**
   factor the common math into a primitive both consumers call. The "twin path"
   is the load-bearing failure mode (X4).
   [case](../evidence/coding-elegance.md#ap1-twin-paths)
2. **NEVER** procedural for-loops over a domain that has its own algebra —
   **instead** invoke the dunders. `psi = L.solve(q)` IS the math.
   [case](../evidence/coding-elegance.md#ap2-procedural-loops)
3. **NEVER** thread a boolean flag through a call to switch behaviour —
   **instead** polymorphism / strategy / sum types. Four booleans = 16
   combinations, 12 of them bugs waiting.
   [case](../evidence/coding-elegance.md#ap3-boolean-flags)
4. **NEVER** a stringly-typed dispatch parameter — **instead** a singleton, enum
   or class: `inner_solver="krylov"` becomes `inner_solver=Krylov()`.
   [case](../evidence/coding-elegance.md#ap4-stringly-typed-dispatch)
5. **NEVER** reduce over an anonymous intermediate, period — **instead** name it.
   If you cannot say what the quantity is, in domain language, in five seconds,
   the physics is wrong, not just the code.
   [case](../evidence/coding-elegance.md#ap5-anonymous-intermediate)
6. **NEVER** transcribe MATLAB/Fortran line-by-line — **instead** read the math
   the legacy claims to encode and write from the math.
   [case](../evidence/coding-elegance.md#ap6-transcription)
7. **NEVER** special-case "the boundary cell" / "the first iteration" / "the
   singular point" — **instead** find the abstraction that makes the special case
   a value, not a branch. If the pole is special, you lack the primitive.
   [case](../evidence/coding-elegance.md#ap7-special-cases)
8. **NEVER** validate inputs deep inside a function — **instead** parse at the
   boundary and trust the type. `def f(x: NonNegativeFloat)` is structural; an
   inner `assert x >= 0` is procedural — and `python -O` strips it.
   [case](../evidence/coding-elegance.md#ap8-deep-validation)
9. **NEVER** write code you would have to comment to explain to a domain expert —
   **instead** name in the domain's vocabulary and let the code be the comment.
   Exception: cite the literature (`# WDD per Lewis-Miller Eq. (4.27)`).
   [case](../evidence/coding-elegance.md#ap9-explanatory-comments)
10. **NEVER** more abstraction than you have concrete instances to justify —
    **instead** wait for the second, then extract.
    [case](../evidence/coding-elegance.md#ap10-premature-abstraction)
11. **NEVER** mark code "temporary" without a removal trigger — **instead**
    remove it now, or write the trigger into a tracked artefact (issue, plan,
    test). `# TODO: fix later` with no date/condition/consumer is a memory leak.
    [case](../evidence/coding-elegance.md#ap11-untracked-temporary)
12. **NEVER** ship code that "works but is ugly" alongside intent to clean it up
    later — **instead** clean before merging.
    [case](../evidence/coding-elegance.md#ap12-ugly-shipped)
13. **NEVER** pass `numpy.ndarray` across a module boundary with shape and
    convention encoded only in a docstring or test fixture — **instead** wrap
    distinct physical quantities in `NewType` / frozen dataclasses / Protocols.
    **Densest single cluster: 13 evidenced bugs** (ERR-002, 009, 011, 022, 031,
    034, 040–047). Minimum: keyword-only args + constructor validation;
    preferred: typed wrappers.
    [case](../evidence/coding-elegance.md#ap13-bare-numpy)
14. **NEVER** hardcode a numerical constant derivable from a typed object in
    scope — **instead** compute it: `4π` becomes `quadrature.weights.sum()`
    (ERR-004: the level-symmetric weight sum is not 4π).
    [case](../evidence/coding-elegance.md#ap14-hardcoded-constants)
15. **NEVER** ship a specialised method before the general one when the general
    is reachable — **instead** derive general-first, then specialise by
    collapsing modes. Specialised-first calibrates to accident and
    then constrains the general form to reproduce it (ERR-030, ERR-035).
    [case](../evidence/coding-elegance.md#ap15-specialised-first)
16. **NEVER** rely on sequential ordering to encode a data dependency —
    **instead** encode it as a DAG. ERR-003 and ERR-044/045 both shipped "the
    loop happens to visit cells in a good order".
    [case](../evidence/coding-elegance.md#ap16-sequential-ordering)
17. **NEVER** loosen a test tolerance to paper over a known approximation gap —
    **instead** document the gap in docstring AND test, pin the residual with a
    structurally-independent reference proving it is approximation not bug, and
    file the closure issue. The tolerance is a CONTRACT (ERR-036/038 both shipped
    `tol=5e-2`).
    [case](../evidence/coding-elegance.md#ap17-loosened-tolerance)
18. **NEVER** strand an iterative method's convergence data (contraction ratio
    `ρ`, a-posteriori `‖Δx‖/(1−ρ)`, Aitken Δ²) — **instead** home it on the
    object that knows *"previous"*: the **iteration**, not the state and not a
    second field type. **Tell:** a loop that subtracts two iterates and keeps only
    a float; the whole `c→1` false-convergence diagnostic is gone.
    [REFUTED 2026-08-19] This item read the opposite (mint a displacement
    type) until its precedent was retired. **Check before minting a difference
    type, two grep-able questions:** (a) is there a
    canonical zero, distinguished by the domain rather than chosen? (vacuum yes;
    a frame origin no); (b) is superposition physical? (linear operator, two
    solutions in one medium add). Two YES means one type, signed differences,
    diagnostics on the iteration record. Two NO means the affine/torsor frame
    genuinely applies. **The two-sided test for "make illegal states
    unrepresentable":** mint the invariant iff (i) every value the type admits
    is legal AND (ii) every legal value is admitted. (ii) is the half that gets
    skipped, being a claim about the **producers**; when it fails the type
    **refuses correct output** ([M] a converged production solve ships
    `min ψ = −6.4e−1`).
    [history](../evidence/coding-elegance.md#ap18-reversal-history)
19. **NEVER** silence a type-checker red on a NEW symbol with `# type: ignore`
    before proving no principled spelling exists — **instead** declare the
    variance the type actually has, use `if TYPE_CHECKING:` where a subclass
    override obscured the signature, or spell a reduction as the contraction it
    IS (`weights @ stacked`). Only a proven upstream-stub gap
    earns `# type: ignore[specific-code]` + a one-line reason. A bare ignore on
    new code is the typed analogue of #17: a contract silently relaxed.
    [case](../evidence/coding-elegance.md#ap19-type-ignore)
20. **NEVER** let a docstring NAME a primitive (a shape, a named quantity, a
    closed form) that the body then **open-codes by hand** — **instead** call the
    primitive the docstring names, or delete the claim. Two spellings of one
    quantity, the prose looking like a verified contract while unenforced (X3). **Special case:** a
    `.shape` / axis-order claim on a bare-ndarray return is #13 wearing a
    docstring — verify the ACTUAL shape at the seam or wrap the return in a typed
    carrier.
    [case](../evidence/coding-elegance.md#ap20-docstring-open-codes)

> **Floor vs ceiling.** #11 and #12 are the *diagnostic* framing of rules stated
> prescriptively in [coding-standards](../rules/coding-standards.md); kept here
> as recognition signals, not a second copy.

---
## Decision flow

Apply at each decision point, in order.

| At this moment | STOP and apply |
|---|---|
| I do not understand the algebra of the domain | Write it on paper first: primitives, compositions, identities, inverses, invariants |
| The language has syntax for the algebra's operations | Pattern 1 (dunders). If it doesn't: methods on the domain type, never free functions outside it |
| About to write the same math in two places | Pattern 2 — factor into a primitive |
| About to reduce over an intermediate that HAS a domain name | Pattern 3 — name it (truly anonymous FP-noise passes) |
| About to assert an invariant that could be a type | Pattern 4 (parsing untrusted external input at the boundary passes) |
| About to write `solve_<specific>` with siblings coming | Pattern 5 — build the primitive |
| About to add a Protocol/ABC/generic with one user | Pattern 6 — defer, wait for ≥2 |
| About to add a boolean flag to a signature | Anti-pattern #3 — polymorphism / sum type |
| About to special-case a boundary cell / first iteration / singular point | Anti-pattern #7 — make the special case a value |

Then write the code and run the checklist **before committing**.

---

## The elegance checklist

Before committing code, walk this checklist. If any answer is NO, return to the editor.

1. **Reads like the domain?** Read the new code aloud. Does it sound like the math / workflow / protocol it implements? If you hear procedural for-loops where the domain has algebra, NO.

2. **No twin paths?** Is there a second piece of code in the codebase that computes the same mathematical quantity by a different procedure? If yes, NO — factor into a primitive.

3. **All intermediates named, with units?** Every value — every product, every reduction, every function return — should answer "what is this quantity, in physical units?" in under five seconds. In physics code, unnamed-and-unit-less is a category that does not exist; every combination of physical quantities has units, and every unit-bearing combination has a name. If any line of arithmetic does not pass the "name + units" challenge, NO. Dimensionless intermediates are permitted IFF they are deliberately constructed (`τ`, `k_eff`, `c`, `Re`, ...) and themselves named — the dimensionless reduction is the result you sought.

4. **Illegal states unrepresentable?** Every runtime assertion of an invariant — could the type system have expressed it instead? If yes, NO.

5. **Primitive vs product?** Is the code a primitive (small, composable, one concern) or a product (specific, large, multi-concern)? If a product, can it be reduced to a composition of primitives? If yes, NO — refactor.

6. **No unjustified abstraction?** Every Protocol / ABC / generic — does it have ≥2 concrete instances justifying it? If not, NO — inline.

7. **No boolean flag parameters?** No `def f(..., is_x=True, is_y=False)`. If yes, NO — use polymorphism / sum types.

8. **No stringly-typed dispatch?** No `func(..., kind="legacy")`. If yes, NO — use types.

9. **No transcription style?** No MATLAB/Fortran/legacy-language idioms imported wholesale. If yes, NO — re-derive from the math.

10. **No "TODO: cleanup later"?** Either clean now, or write a tracked artefact (issue, plan, test). If the comment exists without the tracked artefact, NO.

11. **No bare numpy across module boundaries?** Every `np.ndarray` that crosses a function or module boundary — does the shape + convention + units + physical-quantity-identity live in the type, not just the docstring? If a swapped argument or transposed matrix would type-check but be wrong, NO. Wrap in `NewType`, frozen dataclass, or Protocol. (This is the densest single bug cluster in the project — 13 evidenced bugs.)

12. **No hardcoded numerical constants where a typed object is in scope?** No `4 * np.pi` when a quadrature object's `.weights.sum()` is reachable. No `0.5236` when `(4/3) * np.pi * r_pole**3 / 8` reads as the geometric quantity it is. If a magic constant exists, NO — derive it from the typed source.

13. **No sequential-ordering-as-silent-contract?** Iteration order in a loop should not encode a data dependency. If "the loop happens to visit cells in the right order" is the dependency mechanism, NO — encode the dependency as a DAG and traverse it explicitly.

14. **No specialised method shipped before its general form?** If you're writing a rank-1 / 1-group / homogeneous / single-region case AND the general form is reachable, derive the general first and specialise. If the specialised case ships first, NO — refactor to general-first development order.

The checklist is fast — under a minute per file. The cost of running it is small; the cost of not running it is Phase F.

---

## Implicit principles the bug catalog teaches

- **Symmetry in math implies symmetry in code.** A symmetry in the math
  (apply/sweep, forward/adjoint, rank-N/rank-1) must be visible in the code;
  breaking it is a drift bug waiting. **Corollary:** if the sweep dispatches
  geometry through a protocol, the matvec must too; no `..._matvec_spherical`
  survives. [case](../evidence/coding-elegance.md#symmetry-in-math-and-in-code)
- **Bugs cluster at corners, not the interior.** The bug-rich region is the
  boundary between regimes (boundary, pole, group coupling, angular coupling);
  probe corners deliberately: the convenient MMS is isotropic-flat, the
  bug-rich one angularly-varying with discontinuous σ_t.
  [case](../evidence/coding-elegance.md#bugs-cluster-at-corners)
- **Cross-cutting concerns belong in the type system, not in coding
  discipline.** "Remember to apply the BC at every entry" is fragile;
  `L = StreamingOperator(V, boundary=B)` makes it unforgettable (ERR-040..047).
  [case](../evidence/coding-elegance.md#cross-cutting-concerns-in-types)
- **Dimensional analysis is a free verification tool.** Units must agree on
  both sides of every line, so an inconsistency IS the bug signal, by
  inspection, with no test run (why Pattern 3 is load-bearing). No `pint` here,
  so: a units comment on every physics declaration
  (`phi  # (nx, ng) [1/(cm²·s)]`), units in every `Returns:`, canonical σ names.
  [case](../evidence/coding-elegance.md#dimensional-analysis-free-verification)
- **Refactors extend forward, not to fit legacy.** An N+1 change reaches for the
  N-dimensional generalisation and lets existing cases fall out as
  specialisations, retiring the legacy shape as it extends; the cost is paid
  once. [case](../evidence/coding-elegance.md#refactors-extend-forward)
- **An unused argument is dead weight iff the math layers it OUTWARD for every
  family** (SN, CP, MoC, diffusion, MC): then drop it from the inner signatures;
  keep it inner only if some family needs it there.
  [case](../evidence/coding-elegance.md#unused-argument-dead-weight)
