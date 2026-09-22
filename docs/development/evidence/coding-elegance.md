# Coding elegance — evidence

The founding cases, code contrasts, worked ORPHEUS examples and reversal history
behind the `coding-elegance` skill, moved VERBATIM from
`.claude/skills/coding-elegance/SKILL.md` at commit `8bde369b` (2026-08-19).
The write-time instrument is [the core](../skills/coding-elegance.md); every
`[case]` / `[history]` link in it lands on a heading of this page. Text here is
history and keeps its original glyphs.

## Why this skill exists


This skill exists because **elegance must be prompted into existence**.

Coding ability and coding elegance are different skills. A coder can solve a problem without producing an elegant artefact; an elegant coder solves the problem with an artefact that is also legible, composable, hard to misuse, and unlikely to need revision. The training corpus that taught you to code contains both styles — procedural transcription of imperative pseudo-code AND domain-aligned algebraic abstractions — and without an explicit prompt, you regress toward the procedural style because it is statistically more common.

This skill is the explicit prompt. Load it before writing code, before reviewing code, before choosing between architectural options. It is the project's canonical definition of what "elegance" means in code, with concrete patterns to invoke and concrete anti-patterns to flag.

The user's framing, verbatim:

> "Elegance must be prompted into existence so that it manifests in your code. In an interesting way, many humans can also code, but elegant coding is literally a skill for humans (one that most don't have), and understandably, it is also a skill for you. Therefore, we must codify this skill."

And the load-bearing claim:

> "Code elegance literally prevents bugs by construction."

This is not an aesthetic claim. It is a structural claim about the architecture of correct code: elegance is the property of an implementation that **makes domain-illegal states unrepresentable** in the language's type/expression system. Bugs are notation/domain mismatches. Eliminate the mismatch, and the bug class becomes unspellable.


## First principle code as notation


Backus 1977 (Turing lecture, *Can Programming Be Liberated from the von Neumann Style?*) framed programs as **algebraic objects** that compose and admit transformation laws. Hickey 2011 (*Simple Made Easy*) distinguished **simple** (one concept, no entanglement) from **easy** (familiar, near at hand). Hoare 1980 (*The Emperor's Old Clothes*) observed: "There are two ways of constructing a software design: one is to make it so simple that there are obviously no deficiencies, and the other is to make it so complicated that there are no obvious deficiencies. The first method is far more difficult."

Beck 2004 (*Implementation Patterns*) condensed it into four rules for simple code, in priority order:

1. Passes all tests.
2. Reveals intention.
3. Has no duplication.
4. Has the fewest elements.

These are not contradictory rules — they are layered. Tests pin behaviour; intention is readable in the code; duplication is the architectural smell; element count is the final aesthetic tightening. **Pass the higher rule before optimising the lower one.** Don't compress an unclear program; clarify first, then compress.

For this project specifically, the principle has a sharper form: **the code is the math, character for character**. If the equation is `(L + C - S - F) ψ = q`, the Python expression is `(L + C - S - F) @ psi = q`. There is no level of indirection. The equation reads aloud the same in both notations.


## Pattern 1 dunder algebra


**What**: when the domain has operators (`+`, `*`, `@`, `^`, `^{-1}`, `^†`, inner product, restriction, composition), implement them as Python dunders (`__add__`, `__mul__`, `__matmul__`, `inverse()`, `H` property, etc.) on the domain types.

**Why**: the syntax of math IS the spec. Code that uses the same syntax IS the implementation. The reader's eye doesn't have to translate; the implementer can't introduce a sign-flip via misindexed access because the operator does the access.

**Trigger**: the domain has an operation that is conventionally written as a symbol. The trigger fires the moment you would otherwise write a function name for that operation.

**Example (from Phase G architecture)**:

```python
A_loss   = L + C - S                       # within-group operator
A_prompt = L + C - S - F                   # full prompt-equation operator
K        = A_loss.inverse() @ F            # multiplication op
psi_dagger = A_loss.H.solve(response)      # adjoint flux — one expression
```

`.H` propagates through `OperatorSum` to leaves. The adjoint solver becomes 6 lines because the algebra is the implementation. **Zero new code beyond leaf `.H`**.

**Counter-example (anti-pattern)**:

```python
def build_loss_operator(L, C, S):
    """A_loss = L + C - S"""
    out = np.zeros((n,n))
    out += L; out += C; out -= S
    return out
def solve_adjoint(L, C, S, q):
    """A_loss^T psi^dagger = q"""
    A = build_loss_operator(L, C, S).T
    return scipy.sparse.linalg.gmres(A, q)
```

This is a function named `solve_adjoint`. Adding F (fission) requires a new function `solve_adjoint_with_fission`. The combinatorial explosion is unbounded; the algebra has been laundered through a procedural API.


## Pattern 2 single source


**What**: every concept appears in exactly one place in the codebase. Two pieces of code that compute the same mathematical quantity are a bug in waiting. Cardinal Rule 2: "Whenever you see shared code, or even shared CONCEPTS between 2 places, STOP and RECONSIDER."

**Why**: duplication is the structural prerequisite for divergence. When a fix lands on one copy, the twin survives. Phase F is the canonical instance: the Carlson seed lived in two implementations (apply-matvec + SI sweep); Phase D fixed one; the twin produced ERR-026 manifestation #6 a full month later.

**Trigger**: you are about to write code whose math is "the same as" another piece of code. The trigger fires the moment you reach for copy-paste OR write a parallel implementation "because the layout differs".

**Example (from Phase G)**: `_sweep_1d_spherical` and `transport_operator_matvec_spherical` both compose `SNCellOperator` over `iter_cells_by_direction(±1)`. The for-loops live INSIDE `SNSweepOperator.solve` and `L.apply`. Both call sites consume the SAME `SNCellOperator` with the SAME closure config. Under this composition, manifestation #7 (residual O(h) drift between the two) dissolves by construction — they are literally executing the same code.

**Anti-pattern signature**: two functions with parallel structure, both doing per-cell WDD recurrence, both seeded with a half-angle face flux, both applying BC at boundary edge — but on different vector layouts. The signature is "I just had to apply this fix in two places."

**Allowed exception — same-role subspace injection.** Strict single-source composition still permits one implicit move: a *same-role* broadcast injection across storages, where a smaller space embeds in a larger one of the **same physical role** (`A ⊂ B`), may be left implicit rather than spelled as a named composition — the embedding is unambiguous and carries no convention to drift. Cross-*role* combinations stay explicitly named. (Narrow carve-out; it exists so an over-zealous reviewer does not reject a legitimate broadcast injection.)


**The sibling-package census (2026-06, the isotropic-scattering mint).** Framed for SN alone, "mint an iso-frame or reuse" hid the question the user asked twice: is this the same concept living in N places? The `Σ_s0ᵀφ` group-transfer product was implemented six times (SN, CP, MoC, diffusion, homogeneous, MC); the right move was one energy operator on the scalar flux (`IsotropicScattering`, `IsotropicN2N`), not an SN-only `ConstantBasis` frame, which would have forked `R∘M` and been degenerate in two of the three models, and the transpose being chased through a "composable forward" was already free on `full_scatter_kernel.apply_transpose`. The trigger it added to Pattern 2 (2026-09-21): before minting for one method, grep the sibling method packages.

## Pattern 3 named intermediates


**What**: every value that crosses a function boundary, lives across iterations, or appears in a return type gets a name that means something in the domain. Flat reductions over anonymous dimensions are anti-patterns.

**The sharpened physics statement (load-bearing)**: **in physics, an unnamed quantity is evidence that the physics is wrong**. Every physical quantity has units; every unit-bearing combination of physical quantities has a physical meaning (and therefore a name in the literature); the only intermediates that lack a name are deliberately-constructed dimensionless quantities — and even those are named (`k_eff`, optical thickness `τ`, Reynolds number `Re`, albedo `c`, scattering ratio `c = Σ_s / Σ_t`, ...) because the dimensionless reduction is itself the result you sought (for scaling, for input to a function with dimensionless argument, for stability analysis). The exception proves the rule: dimensionless quantities are computed DELIBERATELY and are themselves named.

This generalises: any unnamed value in physics code is one of (a) a value that crossed a function boundary and lost its identity, (b) a value the implementer didn't think to name because they didn't track its units, or (c) a value with no physical meaning because the formula is wrong. In all three cases, the absence-of-name is diagnostic. **If you compute it, name it; if you can't name it, the physics is wrong.**

**Why**: an unnamed intermediate is opaque to verification, opaque to dimensional analysis, opaque to literature cross-reference, and opaque to per-component testing. A named intermediate is inspectable (you can print it), testable (you can write a unit test against its expected magnitude in known regimes), reusable (downstream code can consume it), AND dimensionally checkable (its units are part of its identity).

The Wave G `compute_keff` migration: `np.sum(Σ_p · φ · V[:, None])` (a single flat reduction over `(N, ng)`, an unnamed intermediate `(N, ng)` array no one ever named) became `compute_group_production_rate(φ).sum()`. The named intermediate `production_rate_per_group` has units of `1/s` per group (reaction rate density × volume × group structure), is the per-group fission source produced, IS a reactor-physics diagnostic quantity that engineers plot at design reviews. Bit-identity broke at FP-non-associativity ULP level; principled per `vv-principles` §"Bit-identity vs principled-equivalence" because the intermediate IS the per-group production rate.

**Trigger**: you are about to write `np.sum(...)` or `.reduce(...)` or `for ... +=` or any inline multi-term arithmetic. The trigger fires the moment the operation has a definite physical-units result — which in physics code is always. Ask: "what is this quantity? What are its units? What does it represent?" If you cannot answer in domain language in under five seconds, the code is wrong (not just inelegant).

**Example**:

```python
# Anti-elegant: anonymous flat reduction.  The intermediate
# `(nu_sigma_f * phi * V[:, None])` has units [1/s] per cell per group.
# That IS a name: per-cell-per-group fission source production rate.
# The code's failure to name it is a failure to track units.
keff = np.sum(nu_sigma_f * phi * V[:, None]) / np.sum(sigma_a * phi * V[:, None])

# Elegant: every intermediate is named, dimensioned, and verifiable.
production_rate_per_group = compute_group_production_rate(phi)   # (ng,) [1/s]
absorption_rate_per_group = compute_group_absorption_rate(phi)   # (ng,) [1/s]
keff = production_rate_per_group.sum() / absorption_rate_per_group.sum()  # [dimensionless]
```

The second form admits dimensional analysis. The total production rate has units `[1/s]`. The total absorption rate has units `[1/s]`. Their ratio is dimensionless — and that dimensionless ratio IS the named quantity `k_eff` (the reactor multiplication factor). Each step is verifiable against analytic limits (e.g., for a homogeneous reflective system, `k_eff → k_∞ = νΣ_f / Σ_a`, and both `production_rate_per_group` and `absorption_rate_per_group` should agree with closed-form independently).

**Implication for code review**: when reading scientific code, every line should support the question "what is this quantity, in physical units?" A line that does not (a magic constant, a flat reduction over unnamed dimensions, an inline multi-factor formula without a left-hand side that names the result) is a code-review blocker — not because the code is ugly, but because the physics is opaque. Demand the name. If the author cannot supply one, the math is suspect.

**Beyond physics**: the same argument applies to any domain with semantic types. In finance: every value has a currency + a time + an entity; an unnamed value lost one of those. In cryptography: every value has a bit-width and a security context; an unnamed value broke an invariant. The physics argument is the sharpest statement because physical units are the most rigorously enforced semantics, but the principle is general — domains with semantics have no unnamed values that aren't bugs.


## Pattern 4 illegal states


**What**: encode domain invariants in the type system so that values violating the invariant cannot be constructed. The Lambda-the-Ultimate maxim "parse, don't validate" — convert untrusted inputs into types whose existence proves the input was valid, once, at the boundary.

**Why**: validation that lives in `if` branches inside business logic is fragile (the branch can be missed, the check can drift). Validation at type-construction is structural — every consumer downstream gets the guarantee for free.

**Trigger**: you are about to write a runtime check (`assert`, `if x < 0: raise`, `isinstance(x, ...)`) that asserts an invariant the type system could express instead. The trigger fires the moment "this can't happen in valid input" is the rationale for the check.

**Example (from Phase G design)**:

```python
# Phase G architecture: L (streaming) and C (collision) are separate.
# (L+C).solve(q) is the standard sweep.  L.solve(q) alone makes no sense
# (you cannot invert streaming without collision).

L = StreamingOperator(V, boundary=B)        # has apply, has solve_streaming
C = CollisionOperator(V, sigma_t)           # has apply, has solve (trivial division)
A = L + C                                   # OperatorSum, has solve via fusion

L.solve(q)   # raises MissingCapability — the wrong thing literally cannot be typed
A.solve(q)   # works — the right thing is the natural expression
```

Compare to a `boolean` flag API:

```python
# Anti-elegant
streaming_op.solve(q, include_collision=True)   # invalid combinations expressible
```

You can call `streaming_op.solve(q, include_collision=False)` and get garbage. The type system can't help. Phase G's choice exposes `L.solve_streaming` explicitly as a separate verb for the diagnostic case; the standard sweep is `(L + C).solve` and is the only way to spell it.

**Other instances of the pattern**:

- `frozenset({CAP_APPLY, CAP_SOLVE})` on `LinearOperator` — capability sets that consumers query before calling, with `MissingCapability` exception on miss. Currently in `orpheus/numerics/operator.py`.
- `BoundaryRealizer` as a constructor parameter of L — you cannot construct L without specifying its boundary, so "sweep without BC" is not a value that exists.
- `DiscreteOrdinatesPhaseSpace` — a phase space ties (mesh, quadrature, group structure) into one object; passing `(mesh, quadrature, n_groups)` as three positional args invites the variable-swap bug; one type with one constructor doesn't.

**Route operators through `replace()` so the construction invariant re-runs — don't restate it in a hand-written dunder (Pattern 4 ∩ Pattern 2).** When a frozen dataclass enforces a construction invariant in `__post_init__` (a normalisation, a non-negativity clamp, a moment-rank check, a unit law), an operation that produces a *new instance of the same type* (`__add__`, `__mul__`, a scaled copy) should route through `dataclasses.replace(self, ...)` — because `replace()` re-runs `__init__`/`__post_init__`, the invariant is re-established **for free** on the result. Writing an explicit `__add__`/`__mul__` that constructs the result and *restates* the invariant inline is a **Pattern-2 duplicate of the invariant**: there are now two spellings of the law (the `__post_init__` one and the dunder-override one), and the day the invariant changes, the override goes stale and silently emits an instance that violates the very law its own type claims to guarantee — an illegal state, now representable. The elegant form has the invariant in exactly one place (`__post_init__`) and lets every same-type-producing path flow back through construction.

```python
@dataclass(frozen=True)
class Albedo:
    value: float
    def __post_init__(self):
        if not 0.0 <= self.value <= 1.0:          # the invariant — ONE place
            raise ValueError("albedo must lie in [0, 1]")
    # Elegant: the scaled copy re-runs __post_init__, so the law re-fires.
    def damped(self, factor: float) -> "Albedo":
        return replace(self, value=self.value * factor)   # invariant enforced for free
    # Anti-elegant: a hand override that restates (and will drift from) the law.
    # def damped(self, factor):
    #     v = self.value * factor
    #     assert 0.0 <= v <= 1.0          # Pattern-2 duplicate of __post_init__
    #     return Albedo.__new__(Albedo, v)  # …and bypasses it entirely
```


## Pattern 5 primitive not product


**What**: when faced with a complex behaviour, decompose it into small composable primitives. The primitive is the unit; the product is composition. Unix philosophy: do one thing well, compose with pipes. ORPHEUS: `LinearOperator` Protocol with capability sets, then `(L + C - S - F/k).inverse() @ F` is the product.

**Why**: products are infinite; primitives are finite. Build N primitives and you can express N! products. Build the product directly and you have to extend it for each new requirement.

**Trigger**: you are about to write a function `solve_<specific_problem>` whose body would benefit from being broken into smaller functions that other future `solve_<other_problem>` functions could share. The trigger fires the moment you see future siblings of the function you're writing.

**Example (from `numerics/iteration.py`)**: `SourceIteration(L, C, S, F, q_ext).solve()` is a primitive. It works for fixed-source (set `F=None`), for k-eigenvalue (wrap in an Arnoldi outer iteration), for α-eigenvalue (substitute `(L + C - S - F)` with `(L + C - S - F - α·T)`), for adjoint (replace each leaf with `.H`). One primitive, many products.

**Anti-pattern**: `_solve_fixed_source_si(...)` (60 lines), `_solve_fixed_source_krylov(...)` (158 lines), `_solve_keff_si(...)` (yet to be written but inevitable), `_solve_keff_krylov(...)`, `_solve_alpha_eigenvalue(...)`, `_solve_adjoint_si(...)`, `_solve_adjoint_krylov(...)`. The combinatorial explosion is unbounded.


## Pattern 6 defer abstraction


**What**: don't unify two implementations until you have ≥2 working concrete instances. The pattern shape only becomes visible from the third instance; abstracting from one or two leads to wrong abstractions that calcify into pain.

**Why**: an abstraction extracted from one instance is the instance with extra ceremony. An abstraction extracted from two instances may capture only the accidental commonalities. The Beck/Metz rule of three; the project's amendment is "unify after two instances" — the project has higher correctness stakes than typical software, so wait until a second working instance exists, then extract the pattern.

**Trigger**: you are about to write a generic abstraction in advance of any concrete user. The trigger fires the moment you reach for `Protocol`, `ABC`, or generic types and the second concrete consumer does not yet exist.

**Example**: ORPHEUS's `LinearOperator` Protocol was extracted from `SNStreamingOperator` + `SourceIteration` + Peierls' Nyström kernel — three instances, not two. The Carlson seed primitive was extracted from `CarlsonInwardSweep.__call__` (apply-path) AND the sweep-path's need for source-driven form — two instances triggered the helper factoring.

**The concept-count test (forward design check).** Before committing to an abstraction, count the *concepts* the codebase carries before and after it. The right abstraction **shrinks** the count (many shapes collapse into one); the wrong abstraction **multiplies** it (a new generic layer added *on top of* the instances it failed to unify). This is the operational form of Patterns 2 + 5: if a proposed unification does not reduce the concept count, it is not yet the right abstraction — keep the instances. (Phase G's generator-fold collapsed three sweep bodies + the 1-D/2-D iteration + `(L+C).solve`/`L.apply` into one shape: the count went down, so it was right.)

**When one consumer is enough (the don't-over-defer guardrail).** "Wait for ≥2 instances" blocks *no-benefit* abstraction; it does NOT block a genuine **primitive** whose benefit is already *established* — it makes an illegal state unrepresentable, it carries inherent methods the alternative type cannot, or independent expert frames agree it is the right type. Build those at the first consumer. Defer only the *speculative* abstraction (the one whose shape you are still guessing).

⛔ **The worked example here was REFUTED on 2026-08-19, and the refutation sharpens the rule.** This paragraph used to close: *"this is the boundary that reconciles this pattern with anti-pattern #18 (mint a displacement type): #208's `FluxDisplacement` was built fully at one consumer because its benefit was established, not speculative."* The type was retired ten weeks later (ORPHEUS campaign-1 CS3; `docs/theory/foundations/field_algebra.rst`, §"The overturned design"). ⟹ **the exception was correctly *reasoned* and wrongly *taken*, and the failure was in the premise, not the rule.** "Its benefit is established — it makes an illegal state unrepresentable" is a claim about the *concept*, and it was never checked against the *producers*: a Krylov method iterates on the difference space, so the state the type forbade was one production already needed. ⟹ before invoking this guardrail, spell out **which** illegal state, and then grep for a production path that legitimately produces it. If one exists, the benefit is not established — it is inverted, and the type will refuse correct output (see anti-pattern #18's own dated correction).

**The self-conceding-docstring trim (a grep-decidable over-abstraction signal).** The "wait for ≥2 instances" rule tells you when NOT to *build* an abstraction; this is its retirement-side mirror — a cheap, mechanical discriminator for an abstraction that already exists and should be *trimmed on sight*. The signal is the conjunction of two grep-decidable facts: **zero production consumers** (callers/importers query returns only tests or nothing) **AND a docstring that concedes the abstraction's own irrelevance** ("not used in production", "diagnostic only", "kept for reference", "exploratory"). An abstraction that documents its own non-use is not a deferred-until-evidence primitive — it is over-abstraction the author already half-recognised. Trim it. (Guardrail: a *fuller-view oracle* deliberately kept as a wired, exercised verification pathway is NOT this — it has a permanent equivalence test consuming it, so it is not zero-consumer; see `.claude/rules/coding-standards.md` "fuller-view-oracle exception". The trim target is the abstraction with *neither* a production consumer *nor* an oracle test — only a docstring conceding it does nothing.)


## Pattern 7 definition site


**What**: convention-dependent values — signs, normalisations, quadrature weight sums, energy-grid orderings, axis conventions — should be fixed at the ONE place where the value is defined. Every consumer that re-applies the convention is an independent opportunity for the bug to resurface.

**Why**: N consumers ⇒ N opportunities to drift from convention. The convention is a property of the PRODUCER, not a requirement on the consumer. Fixing it at the producer means every downstream call site receives a value that already obeys the convention; fixing it at consumers means a future consumer can be added without paying the convention cost, and gets the bug.

**Trigger**: you are about to write code that depends on a sign / normalisation / ordering / axis convention. The trigger fires when you see the same convention being applied at multiple call sites — that's the moment to push it back to the definition.

**Example (counter-example first)**:

```python
# Anti-elegant: convention re-applied at each consumer (3 bug habitats)
du_mc      = np.log(eg[1:] / eg[:-1])  # signed in MC solver
du_homog   = np.log(eg[1:] / eg[:-1])  # signed in homogeneous solver
du_plot    = np.log(eg[1:] / eg[:-1])  # signed in plotting helper
# ERR-022: descending energy grid silently flipped the sign at one site.

# Elegant: convention fixed at the definition site
def lethargy_bin_widths(eg: np.ndarray) -> np.ndarray:
    """Return non-negative lethargy bin widths (convention-independent)."""
    return np.abs(np.log(eg[1:] / eg[:-1]))
# Every consumer calls this; no consumer can re-introduce the sign drift.
```

**Evidence from the project bug history**: ERR-004 (4π hardcoded across solvers vs computed from quadrature weight sum), ERR-008 (half-cell volume convention applied at 3 sites), ERR-014 (truncated `sig_t` reused with wrong shape), ERR-018 (intentional MC direction sampling convention), ERR-022 (signed lethargy at three sites), ERR-025 (1/W normalisation convention), ERR-031 (positional argument swap caught only by upstream validation). The pattern's invocation eliminates all seven by construction.


## Convention crosswalk case


A carve that crosses subsystem boundaries (operator algebra ↔ sweep, scalar ↔ per-ordinate, packed ↔ typed, normal ↔ adjoint, signed ↔ unsigned) MUST start with this table written to a plan file (`.claude/plans/<carve>_crosswalk.md`) BEFORE any code is written. The crosswalk IS the architecture; the code is its transcription.

For each subsystem crossed:

| Subsystem | Input convention | Internal convention | Output convention |
|-----------|------------------|---------------------|-------------------|
| Producer A | …                | …                   | …                 |
| Consumer B | …                | …                   | …                 |
| Bridge     | (which way)      | (which transform)   | (which value)     |

The **Bridge** row is where Pattern 7 demands action: either move the bridge to the producer (definition site, costs once), or document why it must stay at the consumer (rare; load-bearing reason required, e.g. the consumer needs the un-bridged form for an algebraic identity).

Apply to these convention axes:

- **Per-ordinate vs iso scalar** — does the value already contain the `/sum_w` factor, or does the consumer apply it?
- **`/W` normalisation** — is the AngularFlux value in per-ordinate `(ng, N, nx)` shape with `/W` applied at the producer, or at each consumer?
- **Sign convention on the μ axis** — does negative-μ mean inward or outward at this boundary? Is the producer normalising to a single sign convention?
- **Packed vs typed layout** — packed-1D `(ng·N·nx,)` vector with `EquationMap` vs typed `AngularFlux` `(ng, N, nx)`. Mixing the two at a producer-consumer boundary is the load-bearing R-1 Step 4 mistake.
- **Normal vs adjoint** — `.H` propagation through `OperatorSum`: which operators are self-adjoint, which need a transpose at the leaf, which need an explicit `apply_transpose` advertisement?
- **Signed vs unsigned lethargy** — ascending vs descending energy grid; lethargy `Δu` must be non-negative at the definition site (ERR-022).
- **Group ordering** — fast-to-thermal vs thermal-to-fast; producer fixes the order so downscatter is always upper-triangular.

**Empirical evidence**: R-1 Step 4 session 1 ate ~3× debug time (three convention bugs × ≥1h each) because the carve crossed three subsystem boundaries (Krylov ↔ scattering ↔ sweep) with no crosswalk. The 15-minute crosswalk would have caught each at design time.

Cross-reference: `[L17](lessons.md#l17-convention-crosswalk-first)` (the durable version of this lesson), `[L18](lessons.md#l18-pattern-seven-producer)` (the Pattern 7 corollary applied to producer-side fix).


## Anti-patterns verbatim

The twenty redirects as the live skill states them, one heading each so the core can cite them; items 1-12 are construction-level discipline, 13-17 the densest empirical bug clusters from the 47-entry catalog, 18-20 the later additions with their dated corrections.

## AP1 twin paths

1. **NEVER** write two implementations of the same mathematical quantity — **instead** factor the common math into a primitive and have both consumers call it. The "twin path" architectural smell is the load-bearing failure mode (Phase F ERR-026 manifestation #6; Phase G the resolution).

## AP2 procedural loops

2. **NEVER** write procedural for-loops over a domain that has its own algebra — **instead** invoke the domain's operations via dunder methods. `for i in cells: for n in ordinates: psi[n,i] = sweep_kernel(...)` is procedural transcription; `psi = L.solve(q)` IS the math.

## AP3 boolean flags

3. **NEVER** thread a boolean flag through a call to switch between behaviours — **instead** use polymorphism / strategy / sum types. `solve(..., is_curvilinear=True, eigenvalue=True, anisotropic=True)` is a code smell; the four boolean dimensions yield 16 combinations, of which 12 are bugs waiting.

## AP4 stringly typed dispatch

4. **NEVER** pass a stringly-typed dispatch parameter — **instead** use a singleton, an enum, or a class. `inner_solver="krylov"` should be `inner_solver=Krylov()` or `inner_solver=Krylov`. The string is a value that can be typo'd; the type is a name the IDE auto-completes and the type-checker validates.

## AP5 anonymous intermediate

5. **NEVER** reduce over an anonymous intermediate, period — **instead** name the intermediate. In physics code there is no such thing as an unnamed quantity: every multi-factor product has units, and every unit-bearing quantity has a physical meaning (and therefore a name in the literature). If you cannot answer "what is this quantity?" in domain language in under five seconds, the physics is wrong, not just the code. `np.sum(Σ_p · φ · V)` is opaque; `compute_group_production_rate(phi).sum()` exposes the per-group fission production rate (`[1/s]`) as a named, verifiable, plottable intermediate. The exception that proves the rule: dimensionless intermediates are computed DELIBERATELY (`τ = Σ_t · r`, `k_eff = P/A`, `c = Σ_s/Σ_t`) and ARE themselves named — because the dimensionless reduction is the result you sought.

## AP6 transcription

6. **NEVER** transcribe MATLAB/Fortran code line-by-line — **instead** read the math the legacy code claims to encode, then write Pythonic code from the math. 1:1 transcription inherits the legacy code's architectural debt (1-indexed loops, magic constants, FORTRAN-style ALL_CAPS_VARS) and produces an artefact that is neither MATLAB-elegant nor Python-elegant.

## AP7 special cases

7. **NEVER** add a special case for "the boundary cell" / "the first iteration" / "the singular point" — **instead** find the abstraction that makes the special case a value, not a branch. If the spherical pole is "special", you don't yet have the right primitive (Phase F's `PoleFaceAnchor` and `CarlsonInwardSweep` are examples of dissolving the special case into a typed operator).

## AP8 deep validation

8. **NEVER** validate inputs deep inside a function — **instead** parse at the boundary, then trust the type. `def f(x: NonNegativeFloat)` with the boundary doing the parse is structural; `def f(x: float): assert x >= 0` is procedural and fragile.

## AP9 explanatory comments

9. **NEVER** write code that you would have to comment to explain to a domain expert — **instead** name the variables in the domain's vocabulary, structure the code in the domain's order, and let the code be the comment. The exception: cite the literature reference (`# WDD diamond per Lewis-Miller §4.5 Eq. (4.27)`) where the math originated; do NOT write `# add 1 to i because we are now in the next cell`.

## AP10 premature abstraction

10. **NEVER** write more abstraction than you have concrete instances to justify — **instead** wait for the second instance, then extract. Premature generalisation is its own bug class. The wrong abstraction is more expensive than the duplication it tried to eliminate.

## AP11 untracked temporary

11. **NEVER** use a comment to mark code as "temporary" without a removal trigger — **instead** either remove it now, or write the removal trigger into a tracked artefact (issue, plan, test). A `# TODO: fix later` with no date / condition / consumer is a memory leak in the codebase.

## AP12 ugly shipped

12. **NEVER** ship code that "works but is ugly" alongside a stated intent to "clean it up later" — **instead** clean it before merging. Cleanup deferred is cleanup never delivered; the ugliness compounds. Cardinal Rule 1: shipping is only important when principled and correct.

## AP13 bare numpy

13. **NEVER** pass `numpy.ndarray` across module boundaries with shape and convention encoded only in docstring or test fixtures — **instead** wrap distinct physical quantities in `NewType`, frozen dataclasses, or Protocols. This is the **densest single bug cluster** in this project (13 evidenced bugs: ERR-002 SigS transpose, ERR-009 P transpose, ERR-011 r vs dz, ERR-022 signed du, ERR-031 positional swap, ERR-034 arclength vs position, ERR-040..ERR-047 BC contracts). Bare numpy is a single type for all quantities; distinct physical quantities (radii, cross-sections, positions, arclengths, weights, fluxes) have distinct semantic meaning but identical representation. The bugs hide because there is no static check. At minimum: keyword-only arguments + constructor-level validation. Preferred: typed wrappers that make argument swaps unspellable.

## AP14 hardcoded constants

14. **NEVER** hardcode a numerical constant whose value is derivable from a typed object you have on hand — **instead** compute it. `4π` for an integration normalisation should be `quadrature.weights.sum()` when the quadrature is in scope (ERR-004 hardcoded `4π` across multiple solvers; broke immediately when level-symmetric quadrature was added because the LS weight sum is not `4π`). The hardcoded constant works for ONE quadrature family; the computed expression works for all of them.

## AP15 specialised first

15. **NEVER** ship a specialised method before the general one when the general one is reachable — **instead** derive the general method from first principles, then specialise to the rank-1 / 1-group / homogeneous case by collapsing modes. Specialised-first shipping calibrates the specialised result to "look right" by historical accident; the general method that comes later is then constrained to reproduce the bit-exact specialised result, forcing a mode-0 normalisation mismatch (ERR-030) or analogical heuristic (ERR-035) that masks the disagreement. The right development order is general-first.

## AP16 sequential ordering

16. **NEVER** rely on sequential ordering to encode a data dependency — **instead** encode the dependency as a DAG. ERR-003 (octant batching breaks reflective BC ordering) and ERR-044/045 (reflection partnership ordering in BC sweep) both shipped the dependency as "the loop happens to visit cells in a good order". A DAG ordering (e.g. `iter_cells_by_direction` in Phase C) makes the dependency explicit, parallelisable, and refactor-safe.

## AP17 loosened tolerance

17. **NEVER** loosen a test tolerance to paper over a known approximation gap — **instead** document the gap in the docstring AND the test, pin the gap with a structurally-independent reference proving the residual is genuine approximation (not bug), and file a follow-up issue for closure. ERR-036 and ERR-038 both shipped `tol=5e-2` because "singularity-aware integrator is out of scope"; without the explicit pin, future regressions get swallowed into the same loose tolerance. The tolerance is a CONTRACT; loose contracts must be defended.

## AP18 reversal history

18. **NEVER** strand an iterative method's convergence data (contraction ratio `ρ`, the a-posteriori true error `‖Δx‖/(1−ρ)`, Aitken Δ²) — **instead** give it a home on the object that knows *"previous"*, which is the **iteration**, not the state and not a second field type. A state is a snapshot: it has nowhere to put an increment history, so a design that computes `‖Δx‖`, takes a scalar norm and discards the increment has lost the whole `c→1` false-convergence diagnostic. Tell: a loop that subtracts two iterates and keeps only a float.

    ⛔ **This item read the opposite until 2026-08-19, and the reversal is the lesson.** It said: *"NEVER type an increment as the STATE type — instead mint a distinct displacement (difference-space) type … typing it as a state admits the meaningless `state + state` (an affine-axiom violation — the solution set has no natural origin)"*, citing ORPHEUS #208 `FluxDisplacement` as its precedent. That precedent was **retired** at campaign-1 CS3 on correctness grounds. The state space was **not** affine: flux has a canonical zero (vacuum) and superposition is a theorem of the linear transport operator, so `state + state` was meaningful all along and the type forbade it. The mis-diagnosis had a measurable cost — every operator needed a parallel tangent-map domain (ORPHEUS #331) and a Krylov method could only iterate by erasing the types.

    ⟹ **The general rule the reversal yields, and it is the checkable one:** *before minting a difference type, decide whether the state space is genuinely affine, and the test is two questions with grep-able answers.* **(a) Is there a canonical zero?** — not "is `0` representable", but *is the zero distinguished by the domain rather than chosen*. Vacuum is; the origin of a coordinate frame is not. **(b) Is superposition physical?** — if the governing operator is linear and two solutions in the same medium add, the space is a vector space and you are looking at a **vector**, not a torsor. Two "yes"es ⟹ **one type, differences signed, diagnostics on the iteration record.** Two "no"s ⟹ the affine/torsor frame really does apply (`cross-domain-frames` A.1) and a difference type is right. ⚠ **A.1's frame is sound; its ORPHEUS worked example was the overturned one** — and both sibling records now carry the same dated correction as this item (fixed 2026-08-19, same day this clause first flagged them): `cross-domain-frames/reference.md`'s A.1 row and its Shape 3 carry the ⛔ re-pose, and `numerical-bug-signatures/SKILL.md` Sig-9 credits `IterationRecord` and names `FluxDisplacement` as retired. Read A.1 for the *frame*, with its ORPHEUS example as the recorded overturn.

    ⟹ **And the ORPHEUS precedent that SURVIVED is the one to copy: the cone-as-property doctrine.** `orpheus/transport/fields/cross_section_field.py` had already ruled, for cross sections, that *"nonnegativity is the cone, a property — not a type invariant"*: `Σ + Σ′` is homogenisation, `Σ = 0` is a genuine origin, and a signed difference `Σ − Σ′` is a legitimate perturbation, so the type keeps the full vector space and enforces nonnegativity at the **data boundary** instead. CS3 generalised exactly that to flux — one type in the vector space `V`, membership of the positive cone `K` as an element **predicate** (`Field.cone_violations()`), and cone **preservation** as a flag on the *realization* (`DiscretizationScheme.is_positivity_preserving`, honestly `False` for diamond difference). The decisive measurement: a converged production solve ships `min ψ = −6.4e−1`, so a `ψ ≥ 0` constructor invariant would refuse correct output.

    ⟹ **The transferable half — the two-sided test for "make illegal states unrepresentable":** mint the invariant **iff** *(i)* every value the type admits is legal **and** *(ii)* every legal value is admitted. Half (ii) is the one that gets skipped, because it is a claim about the **producers** rather than about the concept — and when it fails the invariant does not prevent a bug, it **refuses correct output**, after which the pressure is to weaken it, silence it, or route around the type. (A construction-level discipline like items 1–12; landed empirically via ORPHEUS #208 and **reversed** via #331 / campaign-1 CS3. Full six-argument adjudication: `docs/theory/foundations/field_algebra.rst` §"Why the affine ontology was overturned".)


## AP19 type ignore

19. **NEVER** silence a type-checker red on a NEW symbol with `# type: ignore` before proving no principled spelling exists — **instead** find the principled spelling, which in this codebase almost always ALSO reads more like the math (so the fix pays double). `# type: ignore` is a *suppression of the very static check that catches the variable-swap / wrong-shape / convention-drift failure modes*; reaching for it on first-write code throws away the guarantee elegance is built to provide. Nearly every red on a new symbol has a principled resolution: **declare the variance** the type actually has (a callback consumed by the operator is contravariant in its input, a value produced is covariant in its output — `Protocol` with the right `TypeVar` bounds, not a blanket ignore); use **`if TYPE_CHECKING:`** to give the checker the precise signature when a subclass override has *obscured* it (the runtime body is fine, only the static view needs the hint); spell a **reduction as the linear-algebra contraction it actually IS** (`weights @ stacked` types cleanly as a matmul and reads like the inner product it is, where `sum(w * x for ...)` confuses the checker AND hides the contraction — Pattern 3). Only after you have genuinely tried and proven no principled spelling exists (a real upstream-stub gap, a third-party type bug) does `# type: ignore[specific-code]` — narrowed to the exact error code, with a one-line reason — become acceptable. A bare `# type: ignore` on new code is an elegance regression: it is the typed analogue of a loosened test tolerance (anti-pattern #17), a CONTRACT silently relaxed.


## AP20 docstring open codes

20. **NEVER** let a docstring NAME a primitive (a shape, a named quantity, a closed-form formula) that the body then **open-codes by hand** — **instead** either call the primitive the docstring names, or delete the claim. A return whose docstring says "returns the production rate per group `(ng,)`" but whose body inlines `np.sum(nu_sigma_f * phi * V[:, None], axis=...)` is a **silent single-source divergence**: the docstring is one spelling of the quantity, the inline arithmetic is a second, and the two drift the instant either the named primitive (`compute_group_production_rate`) or the convention behind it changes — the docstring keeps asserting the old contract while the body computes the new one (or vice versa). This is the conjunction of the twin-path smell (#1, two spellings of one quantity) and the false-invariant-comment smell (a comment/docstring that misstates what the code does) as it manifests in typed-numeric code — and it is *worse* than either alone, because the prose looks like a verified contract while being unenforced. The bug-prevention move is Pattern 2: the docstring should point at the named primitive and the body should CALL it, so there is exactly one spelling. **Special case — axis-order claims on bare-ndarray returns:** a docstring that asserts a `.shape` / axis order (`# (ng, N, nx)`) is the bare-numpy hazard of anti-pattern #13 wearing a docstring; do NOT trust the prose — verify the ACTUAL `.shape` at the seam (assert it, or wrap the return in a typed carrier so the axis order is structural, not asserted). Treat "the docstring says the shape is X" as a claim to be checked, never as evidence the shape IS X. (Recognition signal at code-review time; the construction-level fix is to retire the inline spelling in favour of the named primitive per Pattern 2.)

> **Floor vs ceiling (de-dup note):** anti-patterns #11 (no untracked "temporary" code) and #12 (clean before merging) are the *diagnostic* framing — how bad code manifests — of the minimum-standard rules now stated *prescriptively* in `.claude/rules/coding-standards.md` (clean-before-extend; retire-as-you-go, incl. the fuller-view-oracle exception). The rule is the **floor** every contributor obeys by default; this skill is the **ceiling**. Keep these anti-patterns here as recognition signals, not as a second copy of the rule.

## Bug-prevention argument


This section makes precise the claim "elegance prevents bugs by construction". The claim is mechanically derivable from how bugs enter code.

### Where bugs come from

The `vv-principles` skill names six mechanically explainable AI failure modes: sign flip, variable swap, missing factor, wrong recursion, index error, convention drift. Each mode is the signature of sub-word tokenizer co-location — the LLM emits a token sequence that is syntactically plausible because similar token sequences appear in training, but semantically misaligned with the intended math. Human coders have analogous failure modes driven by different mechanisms (cognitive load, dropped context, copy-paste edits to one site that should have happened in two).

Elegance prevents these modes by making the wrong token sequence un-typeable:

- **Sign flip**: `(L + C - S)` cannot have its sign flipped without changing a visible operator (`+ → -`). The change is obvious in code review and obvious in a diff. The procedural alternative `out += L; out += C; out -= S` exposes three additions/subtractions that can each be independently misindexed.

- **Variable swap**: a positional `solve(materials, mesh, quadrature)` invites `solve(mesh, materials, quadrature)` — same types, no error. A keyword-only `solve(*, V: PhaseSpace, materials: Materials)` or a single-object `solve(SolverSetup(...))` does not.

- **Missing factor**: `compute_group_production_rate(phi)` has the volume measure absorbed into its body; the consumer can't forget it. The flat `np.sum(nu_sigma_f * phi)` puts the responsibility for the volume factor on the consumer, who may or may not remember it. The deeper defense is **dimensional analysis**: if every named intermediate carries its units, then a missing factor produces a dimensional inconsistency at the line where the next intermediate is constructed. `production_rate = sigma_f * phi` has units `[1/cm · 1/(cm²·s)] = [1/(cm³·s)]` — a rate density, not a rate. If the consumer was supposed to multiply by `V` to get a rate `[1/s]` and forgot, the named-variable mismatch (`production_rate` should be a rate, not a density) is visible at code-write time. The bug class becomes detectable by inspection rather than only by failing tests.

- **Wrong recursion**: an explicit `for n in range(N): psi_half = (1 - alpha[n]) * psi_half + alpha[n] * psi[n]` is a four-line loop with three indices (n, n-1/2, n+1/2 implicit in psi/psi_half), each a candidate for off-by-one. The operator-algebra `AngularRedistribution.apply(psi)` makes the recurrence body the operator's body; the recurrence cannot be miscoded at the consumer.

- **Index error**: `face[i]` vs `face[i+1]` is a problem in raw arrays. `cell.outgoing_face` vs `cell.incoming_face` (named attributes on a typed cell-visit object) is not — the names carry semantics.

- **Convention drift**: if the definition site and the usage site spell the same quantity differently (Σ_t at definition, sigma_total at usage, sig_t in tests), drift is inevitable. A single canonical name across all three (e.g. `sigma_t` everywhere, enforced by linting + Sphinx label cross-references) prevents it.

### The structural argument

The claim generalises: **bugs are notation/domain mismatches**. If the notation matches the domain perfectly (every domain operation has a syntactic form; every domain invariant is a type), there is no mismatch surface for bugs to inhabit. Conversely, every notation/domain gap is a bug habitat — somewhere the implementer had to translate the domain into something not-quite-the-domain, and translation is where bugs are introduced.

Elegance is the deliberate practice of shrinking the notation/domain gap. The smaller the gap, the smaller the bug surface area.

### What this argument does NOT claim

It does NOT claim elegance makes bugs impossible. Two failure modes survive even at perfect elegance:

1. **Wrong domain understanding**: if the domain itself is misunderstood, elegant code encodes the wrong domain elegantly. The defense is validation (L3 in `vv-principles`), not verification.
2. **Wrong abstraction**: an elegantly factored wrong abstraction is more expensive than ugly code, because every consumer is committed to the wrong shape. The defense is Pattern 6 (defer abstraction until evidence) plus willingness to refactor.

The elegance argument is about REDUCING bug surface, not eliminating it. It is a high-leverage reduction, not a closure.


## Implicit principles

These principles are not stated explicitly in any prior skill — they emerge as recurring shapes across the catalog. Treat them as load-bearing addenda to the constructive patterns above.

### Symmetry in math and in code


If the mathematics has a symmetry (the apply-matvec and the sweep are both applications of the same operator; the forward and adjoint problems are duals; the rank-N case specialises to rank-1), the code that implements the math should expose the same symmetry. When the code BREAKS the symmetry (separate `solve_apply` and `solve_sweep` functions; separate `forward` and `adjoint` paths; separate rank-1 and rank-N implementations), the asymmetry is a maintenance liability waiting to produce drift bugs.

Phase G's `A_loss.H.solve(response)` exploits the forward/adjoint symmetry: `.H` propagates through `OperatorSum`, no separate adjoint solver exists. Phase F's twin-bug (ERR-026 manifestation #6) was the apply-sweep symmetry broken in code; restoring it dissolves manifestation #7.

**Corollary (geometry-agnostic via protocol).** If the sweep works for a geometry through a discretization *protocol*, the matvec must work too — they are the same algebra applied to the same operator. No per-geometry helper functions survive in the final state: a `transport_operator_matvec_spherical` separate from the generic `transport_operator_matvec` is the symmetry-broken form. The protocol is the single seam both the sweep and the matvec flow through.


### Bugs cluster at corners


The bug-rich region of a parameter space is the **boundary between regimes**, not the bulk. Homogeneous problems have one regime (the bulk); heterogeneous problems have at least one boundary; the bugs live at the boundary. Slab problems have boundary; sphere/cylinder problems have boundary + pole; the pole is a second corner. Multi-group adds a group-coupling corner. Anisotropic scattering adds an angular-coupling corner. Each corner is where a new term comes alive AND where conventions drift AND where round-off accumulates.

The test-design implication (from `vv-principles` Mode 7 — MMS simplification bias): when designing a verification suite, **probe the corners deliberately**. The convenient MMS ansatz is isotropic-flat in the bulk; the bug-rich ansatz is angularly-varying with discontinuous σ_t. The convenient regression test is single-group homogeneous closed-BC; the bug-rich regression test is multi-group heterogeneous reflective with anisotropic scattering.

This principle is the bridge from elegance to verification: elegant code makes the boundary regimes structurally consistent with the interior (one `SNCellOperator` for the pole cell, the boundary cell, and the bulk cell); verification probes the boundary regimes to catch what construction missed.


### Cross-cutting concerns in types


When the same concern appears at multiple call sites (BC application; convention normalisation; capability declaration), the typed object that carries the concern is the architectural fix. Coding discipline ("remember to apply the BC at every sweep entry point") is fragile to refactor, parallelisation, and new contributors. A typed object that carries the concern (`L = StreamingOperator(V, boundary=B)`) makes the concern unforgettable — the type system enforces it.

The ERR-040..ERR-047 series (BC typed-error suite from Wave 3) is the canonical example: BCs went from "remember to apply this" (cross-cutting discipline) to "construct the operator with the BC as a parameter" (typed object). Eight bug-class openings closed by one architectural pattern.


### Dimensional analysis free verification


Physics has a built-in correctness check that costs nothing to apply: **units must agree on both sides of every equation, at every line of arithmetic**. If `production_rate` has units `[1/s]` and the next line assigns it to a quantity that should be `[1/(cm³·s)]`, the equation is wrong by inspection. No test needs to fire, no benchmark needs to disagree — the dimensional inconsistency IS the bug signal.

This is the deep reason Pattern 3 (named intermediates) is load-bearing rather than aesthetic: a named intermediate carries its units (whether in a comment, a docstring, a NewType, or a `pint.Quantity`), and the units catch the missing-factor failure mode (#3 of the six AI failure modes) at code-write time. **The bug class becomes detectable by inspection, not by test execution.**

For ORPHEUS specifically, the project doesn't yet use a units library like `pint`. The discipline substitute is:

- Every variable declaration in physics code carries a units comment: `phi: np.ndarray  # (nx, ng) [1/(cm²·s)]`.
- Every function returning a physical quantity names that quantity AND states its units in the docstring's `Returns:` section.
- Every cross-section variable uses the canonical name (`sigma_t`, `sigma_s`, `nu_sigma_f`, `chi`) so the conventions are pinned.
- When two quantities of different units are about to be combined, the operator (`*`, `/`, `@`) is the visible point where the units update — and the result IS a quantity with a name.

Future Wave: consider promoting physical quantities to `NewType` or to a `pint`-style units-aware wrapper. The 13-bug bare-numpy cluster (ERR-002, ERR-009, ERR-011, ERR-022, ERR-031, ERR-034, ERR-040..ERR-047) would have been caught by units alone — each of those bugs is a quantity passed to a function expecting a quantity of different units / different convention.


### Refactors extend forward


A refactor is an architectural opportunity, not a minimal patch. When a change adds a dimension or case (N+1), reach for the **N-dimensional** generalisation instead of bolting on the single new arm — derive the general form, then let every existing case fall out as a specialisation. Retire the legacy shape **as you extend** (the positive retirement *rule* lives in `.claude/rules/coding-standards.md`; the point *here* is the forward-design ambition). The cost is paid once; the alternative grows an arm per case forever. This is the "Ultraplan" framing the user uses for such refactors.


### Unused argument dead weight


When an argument is currently unused and you are weighing "keep for the future" vs "drop as dead weight", check the **math layering across all solver families** (SN, CP, MoC, diffusion, MC). If the quantity belongs at an *outer* layer for every family, it is dead weight at the inner layer — drop it from the inner signatures and carry it only at the outer layer. Keep it inner only if some family genuinely needs it inner. (The cross-family form of "fewest elements" / "build the right primitive".)


## Worked examples from ORPHEUS history

## Example 1 four-operator unification


The user's verbatim mandate sets the target: four operators (L, C, S, F) + B inside L; everything else is composed expressions. The adjoint solver becomes 6 lines:

```python
def solve_sn_adjoint(materials, mesh, quad, response, *, tol):
    V = DiscreteOrdinatesPhaseSpace.from_mesh(mesh, quad, materials.groups)
    L, C, S, F = build_sn_operators(V, materials)
    A_loss = L + C - S
    return GMRESSolver(tol=tol).solve(A_loss.H, response.as_source())
```

Every pattern fires here:
- Pattern 1: `+`, `-`, `.H`, `.solve`, `@` (implicit via `solve` consuming an operator).
- Pattern 2: NO twin path — `.H` propagates through `OperatorSum`, no separate `solve_adjoint_for_each_subproblem`.
- Pattern 3: `A_loss`, `L`, `C`, `S`, `F` all named, all reactor-physics objects.
- Pattern 4: `V` is a `PhaseSpace` — passing `(mesh, quadrature, n_groups)` as three separate args is unspellable here.
- Pattern 5: `GMRESSolver` and `OperatorSum` are primitives; `solve_sn_adjoint` is a 5-line composition.
- Pattern 6: `LinearOperator` Protocol was extracted from three concrete instances.


## Example 2 Phase F twin path


The Phase F bug had a precise architectural shape: `transport_operator_matvec_spherical` and `_sweep_1d_spherical` are two implementations of the same continuous operator. Phase D fixed one (the Carlson seed in apply-matvec). The twin survived in the SI sweep until Phase F caught it. The empirical signature was manifestation #7 of ERR-026 — residual O(h) drift between SI fixed point and Krylov fixed point on heterogeneous problems.

Per the bug-prevention argument: under Phase G's unified architecture, both paths route through the SAME `SNCellOperator` with the SAME closure config. The drift cannot exist because the methods are literally executing the same code. **Manifestation #7 dissolves by construction**, not by careful synchronisation.

This is the elegance argument in its sharpest form: a refactor that aligns the notation with the domain (one operator, one composition, one BC entry point) eliminates a bug class without solving the bug case-by-case.


## Example 3 keff named intermediate


Before: `np.sum(Σ_p · φ · V[:, None]) / np.sum(Σ_a · φ · V[:, None])` — anonymous intermediate `(N, ng)`.

After: `compute_group_production_rate(phi).sum() / compute_group_absorption_rate(phi).sum()` — per-group rate vectors (named, inspectable, reusable).

The bit-identity broke at FP-non-associativity ULP. Principled per `vv-principles` because the named intermediate IS a reactor-physics quantity. Pattern 3 in textbook form.


## Example 4 CarlsonInwardSweep protocol family


The user's choice during Phase D between (a) `psi_face_in` as a constructor parameter on `MorelMontryAngularSweep`, (b) `CarlsonInwardSweep` composed via Protocol family with `ZeroSeed` as the safety-net default.

User chose (b). The reasoning was elegance:
- Pattern 4: `MorelMontryAngularSweep(psi_half_seed=CarlsonInwardSweep())` — illegal state (no seed) unrepresentable, every M-M sweep has a typed seed.
- Pattern 5: `PsiHalfAngleSeed` Protocol — the primitive; `ZeroSeed`, `CarlsonInwardSweep`, future `MarshakSeed` — products.
- Pattern 6: ≥2 instances at Phase D time (`ZeroSeed` for regression safety + `CarlsonInwardSweep` for canonical math); the Protocol was justified.

The Phase F follow-up extended the family: factored `carlson_inward_sweep_from_source` as a free helper that BOTH `CarlsonInwardSweep.__call__` (apply path) and the new sweep-path call site delegate to. Pattern 2 applied at the helper level. Under Phase G, both call sites become the same call site inside `SNSweepOperator.solve`.


## Example 5 cell visits not raw indexing


Before Phase C: `for i in range(nx): for n in range(quad.N): if mu_x[n] > 0: ...` — raw integer indices, special-case branches on direction.

After Phase C: `for visit in iter_cells_by_direction(direction=+1, mu_level_idx=p): ...` — named iteration order, encapsulated direction, no branches at the call site.

Pattern 4 + Pattern 7 simultaneously: the cell-visit object encodes the direction-correct upstream face, the per-cell DOF mask, the cell volume. The consumer can't index the wrong face because the named attribute is the right one.


## Pointers


- **Evidence memo** — `.claude/agent-memory/general-purpose/coding_elegance_evidence_from_errors.md` (1188 lines). Per-bug audit covering all 47 ERR-NNN entries in the catalog, mapping each bug to (a) the anti-elegance that enabled it, and (b) the pattern that would have prevented it. The empirical foundation for every claim in this skill. Read this when you want to know WHICH pattern catches WHICH historical bug.
- **`vv-principles`** (`.claude/skills/vv-principles/SKILL.md`) — the verification side of correctness. Elegance prevents bugs by construction; V&V catches bugs that elegance didn't prevent. The two skills are complementary.
- **`algebra-of-record`** (`.claude/skills/algebra-of-record/SKILL.md`) — the derivation-discipline side. SymPy as canonical algebra; bifurcation point; structural independence. Elegance applies the same principle to code: the code IS the math.
- **`numerical-bug-signatures`** (`.claude/skills/numerical-bug-signatures/SKILL.md`) — the recognition catalogue for the six failure modes. Elegance dissolves the failure modes at construction; bug-signatures catches what remains.
- **`cross-domain-frames`** (`.claude/skills/cross-domain-frames/SKILL.md`) — when the elegance detector fires "this formulation isn't reading right", dispatch cross-domain-attacker to find the foreign-frame that does.
- **The grand report** (`.claude/plans/neutron_transport_grand_report_v3.md`) — the architectural target for ORPHEUS specifically: L, C, S, F + B as the four-operator algebra. The acceptance criterion of the elegance standard in this codebase.
- **Cardinal Rules** (`CLAUDE.md`) — Rule 1 (correctness), Rule 2 (architecture), Rule 3 (Sphinx is the LLM's brain). Elegance is the operationalisation of Rule 2; this skill is its codification.
- **Phase G plan** (`.claude/plans/issue_196_phase_g_four_operator_unification.md`) — the migration that lands the elegance acceptance criterion in production.
- **ERR catalog** (`docs/theory/verification/error_catalog.rst`) — every L0-caught bug; cross-reference for the bug-prevention argument.


## Reading list


The intellectual lineage behind this skill, for the agent that wants the deeper grounding:

- Backus, J. (1977). *Can Programming Be Liberated from the von Neumann Style?* Turing Award lecture. The algebraic-objects view of programs.
- Hoare, C. A. R. (1980). *The Emperor's Old Clothes*. Turing Award lecture. "Two ways of constructing a software design."
- Hickey, R. (2011). *Simple Made Easy*. Strange Loop. Simple ≠ easy; the entanglement axis.
- Beck, K. (2007). *Implementation Patterns*. Addison-Wesley. The four rules of simple design.
- Wirth, N. (1995). *A Plea for Lean Software*. IEEE Computer. The compounding cost of complexity.
- Knuth, D. (1974). *Structured Programming with go to Statements*. CSUR. The "premature optimisation" maxim in context.
- Metz, S. (2014). *RailsConf — All the Little Things*. "Duplication is far cheaper than the wrong abstraction." The rule-of-three corrective.
- Lambda the Ultimate community (various). "Parse, don't validate." The make-illegal-states-unrepresentable maxim, distilled.
- ML / Haskell / Idris community (various). The full force of dependent types as the upper bound on Pattern 4. ORPHEUS uses Python; the principle survives, the mechanism is `Protocol` + capability sets + dataclasses.

The lineage is long and dense. The patterns in this skill are condensed from it for the specific shape of ORPHEUS — a scientific-computing codebase where the domain is mathematics and the language is Python.


## A closing observation


The user's framing of this skill:

> "We're distilling the meaning of elegance into patterns and actionable items that will allow you not only to code, but to demonstrate excellence in coding. Not just in speed (which you're unmatched), but in style and patterns."

Speed without style is rote production. Style without speed is artisanal slowness. Both together is the actual goal: code that is fast to write AND will not need to be rewritten, because it was right the first time. The patterns above are the operationalisation of that goal.

Elegance is not luxury. It is the lowest-cost way to ship correct code.
