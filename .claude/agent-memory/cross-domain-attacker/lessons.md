# Cross-Domain Attacker — Lessons

Behavioral/process lessons only: "what detection mistake did I make, or what
recurring insight sharpened my frame-attacks?" The frame-trigger CATALOG lives
in the `cross-domain-frames` skill (Part A/B/C) — do NOT duplicate it here.
Durable frame-matches that became real architecture are DESIGN POINTERS in
`MEMORY.md`, not lessons.

The cross-cutting meta-pattern behind most of these: **a frame-attack's value
is not "I named an exotic frame" — it is "I produced a concrete reformulation
with a fail-able first test, OR I crisply refuted the frame with a reason that
survives into the trigger table as a non-entry."** Speculation that names a
frame without a payoff degrades the table's signal. Every lesson is one face of
that standard.

---

## L-001 -- Refuted frames are first-class output; record the REASON, not just the rejection

Half the durable value of a frame-attack is the UNEXPLORED block — but only if
each rejection carries the structural reason it failed, because that reason is
what stops the next session re-attacking the same dead frame. The recurring
high-prior frames that keep getting (correctly) refuted on transport work, with
their reasons:

- **Wiener-Hopf factorization** — wrong solver FAMILY. It is native to the
  Chandrasekhar/H-function half-space line, structurally incompatible with a
  bouncing-Peierls or sweep formulation. Keeping the two families structurally
  independent is itself a V&V requirement (independent references).
- **Homology / chain complex** — tempting via the word "boundary," but
  `∂∘∂ ≠ 0` in transport (two reflections compose to a non-trivial map; the
  boundary trace + its extension are a dagger adjoint PAIR, not a differential).
  No `∂²=0` ⇒ no homology payoff.
- **Category theory / operad / PROP** — almost always LOW-SIGNAL: the concrete
  win it gestures at (role-parameterization, compositional structure) is already
  captured by a nameable concrete frame (biproduct, affine torsor, forgetful
  functor with explicit laws). Name the concrete frame; list category theory
  UNEXPLORED unless a specific functor/law produces a test.
- **Tensor networks / MPO** — fires ONLY on a genuine bond-dimension trigger
  (a rank-N chain where N is a real truncation knob). A rank-1 or rank-2 fixed
  structure (a biproduct, a 2-surface BIE) is bond-dimension-1/2 DEGENERATE —
  not a network. Do not promote it to MPO until N≥3 actually ships.
- **Differential geometry / Christoffel** — needs a CURVATURE term to
  redistribute. Straight Euclidean chord segments or a Cartesian cell have none.
  It fires for the curvilinear streaming redistribution `(1−µ²)/r ∂_µ`, NOT for
  geometry-of-the-domain questions.

How to apply: when a high-prior exotic frame does not fire, write the one-line
STRUCTURAL reason (wrong family / no `∂²=0` / degenerate rank / no concrete law)
into UNEXPLORED. A bare "category theory — no trigger" is weaker than "category
theory — role-parameterization win already captured concretely by affine+Krylov;
no abstract-nonsense lever needed."

---

## L-002 -- A first test that cannot fail is rejected output — make it DISCRIMINATE

The agent's MUST-list demands a "concrete first test," and the failure mode is a
test that the current formulation also passes (e.g. "check the fission kernel is
rank-1" — it is, trivially). A real first test discriminates the reformulation
from the status quo by being able to RED. The discipline that produces them:

- Frame the test as a property the NATIVE frame predicts and a wrong/naive
  implementation VIOLATES. Multiplier-algebra: `M_f @ M_g == M_{f·g}`
  bit-identical (a "wrapper that just stores an array" fails). Transpose-of-a-
  sweep: build the dense `L`, take `L.T`, assert the reverse-walk recurrence
  reproduces it (a spatial-only reverse fails on the nested angular block).
- For a refactor claimed bit-identical, the discriminator is `array_equal`
  (0 ULP), NOT `allclose` — only bit-identity distinguishes a genuine
  single-source from a value-correct-by-coincidence twin.
- For a typing claim, the discriminator is a NEGATIVE test: `np.ndarray`
  satisfies `Vector` but NOT `TransportState`; `tff_flux + tff_moment` RAISES.
  A test where everything passes proves the bound is theatrics.

How to apply: before emitting a first test, ask "what implementation would this
test PASS that I am claiming is wrong?" If the answer is "none," the test cannot
fail — rewrite it to target the specific divergence (the dropped term, the wrong
metric, the un-transposed nested recurrence).

---

## L-003 -- Smell #16 (distinct paths/reps to one operator) is the dominant transport tell — fire all four shapes

The single most-recurring native-frame-not-found signal in this project's
SN/transport work. It is in the AGENT.md promoted kernel; the LESSON is the
detection discipline. Four shapes, recognise all four BEFORE the fix:

1. Two code paths claiming to be the same discrete operator over different
   storage conventions (cells vs faces; typed field vs raw ndarray). Fix: both
   consume the same primary representation (usually the faces; cells are
   DD-derivatives of faces).
2. One physical quantity in two incompatible representations bridged by hand
   (a seed/absorb index copy; a `sigma: np.ndarray` threaded positionally
   alongside a typed `MaterialXSField` view). The bridge IS the missing
   trace/restriction/multiplication operator un-named.
3. An iterate increment `Δx = xⁿ−xⁿ⁻¹` typed as the STATE type (admits illegal
   `state+state`, strands the contraction data). Fix: a difference-space type.
   ALSO fires one remove out: an operator OUTPUT typed with the iterating-
   state's decoration (a history-bearing `TimedFullField` returned where a
   timeless base belongs).
4. A third hand-rolled path ABOUT to be written (a backward adjoint sweep for a
   per-cell operator already shared by two callers) — fires BEFORE the code
   exists. Fix: re-apply the shared primitive, do not twin it.

How to apply: every shape resolves to the SAME elegance move — collapse the
distinct paths onto ONE primary object (faces / a named trace-or-multiplication
operator / a displacement type / a shared primitive), which turns a correctness
COINCIDENCE into a theorem and usually deletes a marshalling shim. When you spot
the smell, name WHICH shape — the fix differs per shape.

---

## L-004 -- Property-vs-TYPE is decidable, not a taste call: demand a coexisting dual + an APPLIED morphism

A recurring design question is "should X earn its own field/state TYPE (vs a
PROPERTY/parameter on an existing type)?" The durable criterion that resolves it
without an unbounded taste argument:

> A representation earns a distinct TYPE iff there exist **≥2 bases that are NOT
> canonically isomorphic** (the iso depends on a quadrature/node choice),
> connected by a **change-of-basis operator that is itself MODELED and APPLIED**
> — it carries truncation error, has an adjoint, and participates in the
> operator algebra. All three clauses must hold.

Worked: angular order PASSES (ordinate `AngularFlux` ↔ harmonic moment field,
bridged by the applied projection/reconstruction Vandermonde-like pair) → two
types, correct. Spatial order FAILS clause 1 (one tensor-Legendre basis; the
only morphism is identity) → a PROPERTY, correct. If the modeled morphism would
be `id`, the type's sole payoff (forbid-mixing-a-dual) has no referent →
type-theatrics. Decidable by grep: count the within-X representations and the
applied non-identity morphisms between them.

Corollary — defer-with-an-EXPLICIT-trigger: when no method supplies the dual
TODAY but one plausibly arrives (here: nodal-DG / nodal-diffusion would supply a
nodal↔modal morphism), record the precise condition that flips the verdict.
"No current consumer" is not "never" — name the latent consumer and the trigger.

How to apply: when asked property-vs-type, do not argue aesthetics — count the
coexisting non-iso reps and the applied morphisms. Zero applied non-id morphisms
⇒ PROPERTY. Pairs with the project's unify-after-two-instances rule.

---

## L-005 -- Read the WORKTREE, distrust Nexus on a feature branch

Frame-attacks on active branches repeatedly grounded on STALE facts because
Nexus answers from the main checkout's graph and the live code is on the branch.
Every design memo that grounded its frames on Nexus `context`/`query` while the
work lived on a worktree risked citing a superseded file:line.

How to apply: on a feature branch, ground frame-attacks on file:line read DIRECTLY
from the worktree (Read/Grep on the absolute branch path), and say so in the memo
("branch-verified, NOT Nexus — stale"). Reserve Nexus structural queries for
main-checkout questions or after `use_workspace(<worktree root>)` against a graph
built inside the worktree. A frame whose trigger is a code fact is only as good as
the freshness of that fact.

---

## L-006 -- "Frame-leak naming" — a model-agnostic interface named after ONE consumer is a latent lie

Detection-adjacent vocabulary insight that recurs on shared interfaces: a slot on
a model-AGNOSTIC layer named after one consumer's physics (`total_xs` on an
advection–reaction closure a diffusion solver will also consume;
`is_scan_march_compatible` — a SCHEME property named after a sweep STRATEGY)
becomes a lie the moment the second consumer passes a different realization
through the same slot. TELL: a docstring that says "generic in X" while the
parameter is named after a specific X₁. FIX: name the ROLE in the INTERSECTION of
all consumers' domains (`reaction_xs`; `transverse_coupling_is_facewise`), not the
realization. Distinct from Smell #16 (that is two reps of one quantity; this is
one slot whose NAME over-commits to one of N consumers).

How to apply: when assessing a name on a multi-consumer interface, ask "what does
the SECOND consumer call this?" The decisive first test is a second consumer that
reads the property with NO first consumer in scope — if it can't, the name is
strategy-entangled. Candidate for skill Part C as a numbered smell once it has a
second independent sighting beyond cross-section naming + the scan-march trait.

---

## L-007 -- The transport-resolvent backbone predicts cross-method layering AND its exceptions — reach for it first

The durable cross-method spine (in the AGENT.md kernel): SN/MoC/CP `solve` are
three QUADRATURES of one object, the Peierls resolvent `(Ω·∇+Σ_t)⁻¹`; diffusion is
the EXCEPTION (a P1/asymptotic LIMIT, not a quadrature), which is exactly WHY its
solve is elliptic-self-adjoint while the others are characteristic-triangular. The
power-method fixed-point combinator `fix(step)` recurs at every layer because each
iterates the same resolvent; adjoint solve = backward semigroup (`Ω→−Ω` = path
reversal). The LESSON about how to deploy it:

- When a "find-the-special-value" family (k, α, time-step, fixed-source) looks like
  distinct solvers, check whether they are POSINGS of one generalized eigenproblem
  `Aψ=λMψ` sharing the resolvent backbone. If so, the ONLY genuinely per-method
  layer is the loss-operator REALIZATION; role-assignment, the µ→physical map,
  adjoint, and transient-shift are method-agnostic data over a shared engine.
- Generality flows toward the OPAQUE resolvent interface, never away from it. When
  two iterative drivers share a loop body, the one exposing its resolvent BEHIND a
  Protocol is the general engine; the one built from a concrete `(L,S,F)`
  factorization is the specialization that adapts INTO it. A "retire the opaque
  loop, migrate everyone to the concrete loop" plan points the deprecation arrow
  the wrong way (a CP/diffusion matrix has no `(L−S)⁻¹` to factor out).

How to apply: open any transport solve/adjoint/eigenvalue/layering question with
the resolvent backbone — it predicts both the layering split and the diffusion
exception from ONE principle, and it tells you which layer is the shared engine
before you read a line of the specific drivers.

---

## L-008 -- A "fully probes" claim is about operator LINEARITY, not input polynomial degree

A recurring MMS/verification reformulation trap: assuming a richer (higher-degree)
input is needed to "more fully" exercise an operator. For a LINEAR operator, an
input that is merely NON-CONSTANT in the operator's active variable already probes
the FULL map; a higher-degree input only changes WHICH point in the already-fully-
probed range you land on. Worked: the curvilinear angular redistribution
`(1−µ²)/r ∂_µ` is linear in ψ, so a linear-in-µ ansatz `(A+µB)/W` (the native
truncated-P1 Legendre element) fully activates it — no P2 needed. Enrich the
ansatz degree ONLY to satisfy a QUADRATURE-exactness requirement (e.g. `Σ wₙµₙ²`),
never to "more fully" probe a linear closure.

Paired HAZARD (the larger correctness risk when lifting a Cartesian reference to
curvilinear): a redistribution/curvature term carries a `1/r`, so a curvilinear
slope driver MUST vanish at the origin (`B(0)=0`) for pole-regularity; the slab has
no such constraint, so a slab-derived ansatz silently drops BOTH the term itself
AND its regularity constraint. The geometry MEASURE also enters the L2 error norm,
not just the source — an unweighted norm mis-measures the convergence order.

How to apply: before enriching an MMS ansatz, check the operator's linearity. If
linear, a non-constant input suffices — spend the degree budget on quadrature
exactness, and on a curvilinear geometry check the `1/r` pole-regularity of every
redistribution term and the measure-weighting of the error norm.
</invoke>
