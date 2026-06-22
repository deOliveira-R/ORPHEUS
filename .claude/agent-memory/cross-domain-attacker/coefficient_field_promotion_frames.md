---
name: coefficient-field-promotion-frames
description: Frame attack on the CoefficientField → MultiplicationOperator promotion (§5.5–5.7 grand plan) — the multiplier-algebra embedding M: L^∞→B(L²), the coefficient cone/simplex algebra, the Kernel-vs-Coefficient locality criterion, and the honest TransportState generic that retires bare Vector[V].
metadata:
  type: project
---

# CoefficientField promotion — structural frame attack (branch `feature/field-typed-operator-algebra`)

DESIGN MEMO, not a frame-trigger enumeration. Feeds a design memo + a new
GitHub issue + an implementation plan. Four frames pinned, each grounded in
branch-correct code (read on this worktree, Nexus graph NOT trusted — stale).
The user's core complaint — *"if type hinting doesn't help, it's compliance
theatrics"* — is the through-line: every verdict below is judged by whether it
makes a reader's eye STOP hunting (`C = sigma_t.as_multiplication_operator()`
reads as the math; `op.apply(x: V) -> V` does not).

## Verified code facts (the ground truth these frames stand on)

- `orpheus/numerics/vector.py`: `Vector` = structural Protocol with
  `__add__/__sub__/__rmul__/__truediv__`; `V = TypeVar("V", bound=Vector)`.
  `apply(self, x: V, /) -> V` is the #256 signature. The docstring ALREADY
  admits the carrier is really `TimedFullField` and `np.ndarray` only appears
  at the scipy serialization boundary.
- `orpheus/numerics/field.py`: base `Field = (values, space) + UNITS`; dunders
  are a PLAIN vector space; `_check_partner` = same-class + same-space.
- `orpheus/transport/fields/_flux_role.py`: `FluxRole` SPECIALIZES the base to
  an AFFINE TORSOR (`flux+flux → TypeError`, `flux−flux → FluxDisplacement`,
  `affine_combination` with Σλ=1). Scalar scaling is UNTOUCHED (inherited).
- `orpheus/numerics/units.py`: FOUR field unit signatures (a 2×2 areal/volumetric
  × angular/scalar); cross sections are `1/cm` — a FIFTH signature that lives on
  NO field leaf today. The operator unit-gain (#208) is the `cm⁻¹` map
  `ANGULAR_FLUX_UNITS (cm⁻²s⁻¹sr⁻¹) → ANGULAR_RATE_UNITS (cm⁻³s⁻¹sr⁻¹)`.
- `orpheus/sn/material_xs_field.py`: `MaterialXSField` is a bare `@dataclass`
  (NOT a `Field`) holding per-material `Mixture` + lazy per-cell `(ng,*spatial)`
  views (`total_cross_section`, `absorption_cross_section`, `fission_production`,
  `emission_spectrum`) + per-material `(ng,ng)` dense `sig_s_legendre`/`n2n_matrix`.
  Cross sections are raw `np.ndarray`. The class docstring §"Composability framing"
  ALREADY names "Mixing" (homogenisation = weighted volume-average → single field),
  "Restriction" (per-region subset), "Action" (`mat_xs.fission_production * scalar_flux`
  "eventually reads as the math") as the THREE intended operations — unimplemented.
- `orpheus/sn/loss_representation.py`: `loss_action(self, sigma: np.ndarray, psi)`
  — Σ_t is ALREADY passed as an explicit positional `np.ndarray` coefficient
  (#240 Phase 2 Step B). The docstring states `C = σ⊙` (collision = pointwise
  multiply by Σ_t) and `L = (L+C) − C` (Resolution A). The collision operator is
  the diagonal-multiply FOLDED INTO the streaming apply, NOT a standalone leaf.
- `orpheus/numerics/operator.py`: `DiagonalOperator(weights: 1-D ndarray, axis)`
  = broadcast-multiply along ONE tagged axis (self-adjoint, `solve` iff nonzero).
  `RankOneOperator(left, right, axis)` = `|ℓ⟩⟨r|` per-position rank-1 along an axis.
  `MultiplicationOperator` DOES NOT EXIST. The fission kernel is
  `RankOneOperator(χ, νΣf, axis=0) & IdentityOperator()` over bare ndarrays.
- `orpheus/sn/fission.py`: `FissionOperator.apply` is `singledispatchmethod`
  over {legacy AngularFlux, TimedFullField, ScalarFlux, ndarray}; the rank-1
  math is delegated to `self.kernel` (the TP above). Output of the ScalarFlux
  arm is `ScalarSourceSink` (a SOURCE, not a flux). `block_role = BULK`.
- Grand plan §5.5–5.7 names the target: `CoefficientField` (sibling of flux),
  `MultiplicationOperator → {CollisionOperator, TimeMassOperator}`,
  `IntegralKernelOperator → {ScatteringOperator, FissionOperator}`, suffix law
  (Operator=field→field, Kernel=integrated-against-measure, Functional=field→scalar,
  Projection/Reconstruction=to/from coefficients, adjoint-related §4678).

PLAN IMPRECISION FLAGGED UP FRONT (§5.7 vs §21.1): the plan uses the NAME
`multiplication_operator` for TWO different things. §5.7 line 643 = the
`MultiplicationOperator` CLASS (field→operator promotion, this memo). §4158
`CriticalityEigenproblem.multiplication_operator()` = the eigenvalue ITERATION
operator `A_loss⁻¹ @ F` (the resolvent-fixed-point combinator, the promoted
library kernel). These are unrelated. The implementation MUST NOT reuse the
name across both — rename the §4158 verb to `iteration_operator()` or
`fixed_point_map()`. (Cross-domain-attacker library: this IS the resolvent
backbone `(Ω·∇+Σ_t)⁻¹` — §4158's `A_loss⁻¹F` is its power-method quadrature.)

---

## STRUCTURAL FEATURES

- A local (pointwise, diagonal) coefficient σ_t/νΣf/χ stored as raw ndarray,
  multiplied into a field — a multiplication by a function, the symbol of a
  zeroth-order operator.
- A 2×2 algebra split already half-built: `DiagonalOperator` (multiply by a
  coefficient that lives on a SUBSPACE — one tagged axis — broadcasting over the
  complement) vs the per-cell `(ng,*spatial)` full coefficient (multiply by a
  coefficient that lives on the WHOLE phase space minus angle).
- A commutative pointwise (+, ·) structure on the coefficients (homogenisation
  = `Σ_t = Σ_m f_m Σ_t^m`, a convex/volume-weighted blend — an ALGEBRA operation,
  not the affine-flux blend).
- A locality dichotomy: σ_t·ψ, νΣf·φ, χ·(·) are DIAGONAL (local symbol); Σ_s·φ
  integrates over the angular measure (a Schwartz kernel — nonlocal).
- A rank-1 separable object χ⊗νΣf sitting AMBIGUOUSLY between "integral kernel"
  and "multiplication ∘ functional ∘ multiplication".
- A single remaining genuine polymorphism: the cross-method carrier
  `TransportState` (SN `TimedFullField[AngularFlux]`, CP `[ScalarFlux]`, MoC
  `[RayField]`) typed today by the structural-but-unexplanatory `Vector`.
- A measure (the angular cubature weights `DiagonalOperator.from_measure`) that
  separates "diagonal multiplication" from "integrated against a measure" — the
  exact line between Operator and Kernel in §5.6.

## ELEGANCE DETECTOR HITS

- **Smell #16 shape 2 (one quantity in two incompatible representations).**
  Σ_t lives in TWO places bridged by hand: the per-cell `(ng,*spatial)` view ON
  `MaterialXSField` AND the positional `sigma: np.ndarray` THREADED THROUGH
  `loss_action`/`sweep`. The thread IS the un-named restriction/multiplication
  operator. The fix is to name `C = σ_t.as_multiplication_operator()` once and
  let the composite carry it, not re-pass it positionally at every leaf.
- **Smell #16 shape 1 (two paths claiming one operator).** "C = σ⊙" is realised
  TWICE: folded into `StreamingOperator.apply` (`L = (L+C) − C`) AND conceptually
  as the standalone `CollisionOperator` the plan wants. The double-subtract guard
  in the `loss_action` docstring is the tell — a correctness coincidence that a
  named `MultiplicationOperator` leaf makes a theorem.
- **Metric/role-blindness (the #208 candidate smell).** `MaterialXSField` is a
  bare `object` record while every flux/source/residual is a typed `Field` with
  `UNITS`. The coefficient is the ONLY operand in `(L+C−S−F−B)ψ=q` with no
  dimensional type — so the operator unit-gain check (#208) has nothing to read
  on the C/S/F side. The `1/cm` gain is asserted in prose, not in the type.
- **Stringly/positionally-typed coefficient.** `loss_action(sigma, psi)` takes a
  bare ndarray whose meaning (which σ? on which axis? what units?) is carried by
  a docstring, not a type. This is the user's "compliance theatrics" complaint in
  miniature: `sigma: np.ndarray` helps the reader NOTHING.

---

## FRAME CANDIDATES

### Frame 1 — Multiplier-algebra embedding M: L^∞ → B(L²) (field→operator as a *-homomorphism)

**Trigger.** A local coefficient (σ_t, νΣf, χ) multiplied pointwise into a
field — the symbol of a zeroth-order operator; `DiagonalOperator` is already
this restricted to one axis (skill A.3 functional-analysis: "multiplication
operator / spectral theorem for a bounded self-adjoint operator").

**Structural verdict.** The promotion `CoefficientField f ↦ MultiplicationOperator
M_f` (M_f ψ = f·ψ pointwise) IS the **multiplier-algebra embedding**
`M: L^∞(phase space) → B(L²(phase space))`, a unital injective **\*-homomorphism**
of commutative C\*-algebras onto the *maximal abelian subalgebra* of diagonal
operators. This is not a metaphor — it is the precise structure, and it dictates
the laws that BECOME the verification suite. The embedding is a functor from the
(commutative, pointwise) coefficient algebra to the (noncommutative, composition)
operator algebra; it is faithful and preserves the *-structure.

**Laws that must hold (and become L0/L1 tests):**

| Law | Statement | Test |
|---|---|---|
| Homomorphism (product) | `M_f @ M_g == M_{f·g}` | L0 — two coefficients, compose vs multiply-then-promote |
| Unit | `M_1 == IdentityOperator()` | L0 — `CoefficientField.ones(space).as_multiplication_operator()` |
| Zero | `M_0 == ZeroOperator(codomain_zero=…)` | L0 — note codomain ≠ domain (flux→source), so the zero is the CODOMAIN zero (the existing `codomain_zero` hook in `ZeroOperator` is exactly this) |
| Linearity | `M_{af+bg} == a·M_f + b·M_g` | L0 — links coefficient `+`/scalar-`·` to operator `+`/`·` |
| Commutativity of the image | `M_f @ M_g == M_g @ M_f` | L0 — the image is abelian (it is NOT for general operators; S∘C ≠ C∘S) |
| Adjoint (real ⇒ self-adjoint) | `M_f.H == M_{f̄} == M_f` for real f | L1 — under the space metric; ties to #208 G-adjoint. Real cross sections ⇒ C, the time-mass T, are SELF-ADJOINT — a structural fact the streaming-folded form HIDES |
| Spectrum | `spec(M_f) == ess-range(f)` | L1 — invertible iff f bounded away from 0; `solve` capability iff `min|f| > 0` (EXACTLY `DiagonalOperator`'s invertibility gate today) |

**Is `DiagonalOperator` exactly M_f restricted to a tagged axis?** YES, with one
caveat. `DiagonalOperator(weights, axis)` is `M_f` where `f` is a coefficient
whose support is the 1-D `axis` subspace, broadcast over the complement (the
weights are constant along every other axis). `MultiplicationOperator(CoefficientField)`
generalises it: the coefficient lives on the FULL `(ng,*spatial)` phase space
(σ_t varies per group AND per cell), broadcasting only over the ANGULAR axis. So:

- `DiagonalOperator` = M_f, f constant except on one axis (the angular weight W,
  the per-group-only removal in a homogeneous problem).
- `MultiplicationOperator` = M_f, f varying on a *sub-product* of axes
  (`(ng,*spatial)`), broadcasting over the complement (angle).
- Both are the SAME mathematical object (diagonal in the position/group/angle
  basis); they differ only in WHICH axes the coefficient is constant on. The
  clean design makes `DiagonalOperator` a *special case constructor* of
  `MultiplicationOperator` (coefficient = a measure on one axis), NOT a sibling.
  `DiagonalOperator.from_measure` becomes `MultiplicationOperator.from_measure`.

**Elegance payoff.**
- *Structure-exposing*: the `1/cm` unit-gain and the self-adjointness of C/T
  become readable facts on the type (`M_f.H == M_f`), not prose in a docstring.
  The collision double-subtract guard (`L = (L+C)−C`) is exposed as "C is a named
  leaf you add once", killing the coincidence.
- *Expressive*: `C = sigma_t.as_multiplication_operator()` /
  `C = MultiplicationOperator(sigma_t)` reads as "collision is multiplication by
  Σ_t" — the §5.7 promise verbatim. `loss = L + C − S − F − B` composes with C as
  a first-class term instead of a positional ndarray smuggled into `loss_action`.
- *Structurally-simpler*: collapses Smell #16 (one Σ_t in two representations) —
  the coefficient is named once, promoted once, and the `sigma: np.ndarray`
  positional thread RETIRES from the `LossRepresentation` Protocol.
- *Algorithmic-advantage*: none new — M_f is still an O(N) broadcast-multiply;
  the win is correctness-by-construction, not speed.

**First test.** `M_{σ_a} @ M_{σ_s0} == M_{σ_a · σ_s0}` on a 2-group 3-cell mesh,
bit-identical to the pointwise product promoted — discriminates a real homomorphism
from a "wrapper that just stores an array". A `MultiplicationOperator` that fails
this is not the embedding. Second discriminator: `M_f.solve(M_f.apply(x)) == x`
iff `min|f|>0`, and raises `MissingCapability` otherwise (proves the spectrum law).

**Structural attack on current.** The collision operator does not exist as a
leaf — it is a sign convention (`L = (L+C)−C`) folded into the streaming apply
and re-passed as `sigma: np.ndarray` at every `loss_action` call. The algebra
`(L+C−S−F−B)` is therefore a FICTION on the C term: there is no `C` object, only
a coefficient threaded by hand and a guard comment warning you not to subtract it
twice. The frame exposes that C, the simplest operator in transport (multiply by
a number per cell), is the ONLY term in the named algebra with no class.

---

### Frame 2 — The coefficient space is a commutative algebra / cone / simplex, NOT the flux affine torsor

**Trigger.** A commutative pointwise (+, ·) structure on σ with a positivity
constraint (σ ≥ 0) and, for χ, a partition-of-unity constraint (Σ_g χ_g = 1) —
distinct algebraic structures sharing a storage shape (skill A.2 algebra: "an
algebraic structure — group/ring/module/algebra — is present but the code treats
elements as bare arrays").

**Structural verdict.** `CoefficientField` MUST NOT inherit `FluxRole`. The two
spaces have OPPOSITE algebraic structure:

- **Flux** = affine torsor over a displacement space: NO origin, `flux+flux`
  forbidden, only Σλ=1 blends legal. (`_flux_role.py`.)
- **Coefficient** = a commutative **ring/algebra with multiplication**, AND a
  **module over ℝ**, AND a **cone** (σ ≥ 0): `σ + σ` is legitimate
  (homogenisation is a NON-convex sum when not volume-normalised: `Σ_t^mix =
  Σ_m N_m σ_t^m`, number-density-weighted — number densities do NOT sum to 1),
  `σ · σ` is legitimate (the M_f@M_g=M_{fg} product), `a·σ` legitimate. The
  coefficient KEEPS the base `Field` plain-vector-space dunders AND ADDS a
  pointwise `__mul__(other: CoefficientField)` (the algebra product the flux
  never has — flux `__mul__` is scalar-only).

So `CoefficientRole` is the COMPLEMENT of `FluxRole`: where FluxRole REMOVES
`flux+flux` and RETYPES `flux−flux`, CoefficientRole RESTORES nothing (keeps the
base `+`/`−`/scalar-`·`) and ADDS a field-field `*`. Its `_check_partner` is the
BASE one (same-class + same-space + same-mesh), no torsor gate.

**Is "CoefficientField" one role or several?** SEVERAL — this is the sharpest
finding and an amendment the §5.5 plan needs. The coefficients do NOT share one
sub-object:

| Quantity | Algebraic structure | Constraint | Role |
|---|---|---|---|
| σ_t, σ_a, σ_s0(diag) | nonneg cone element of the algebra | σ ≥ 0 | `CrossSectionField` (the multiplier; `1/cm`) |
| νΣ_f | nonneg cone element | νΣ_f ≥ 0 | `CrossSectionField` (same role — a production rate density `1/cm`) |
| χ (emission spectrum) | **probability simplex** element | χ_g ≥ 0 AND Σ_g χ_g = 1 | DISTINCT — `SpectrumField` (dimensionless, simplex-constrained) |
| 1/v (time-mass) | nonneg cone, units `s/cm` | 1/v ≥ 0 | `CrossSectionField`-like but DIFFERENT units (TimeMassOperator's coefficient) |

χ is a normalized emission spectrum — it lives on the **(ng−1)-simplex**, NOT the
cone. Adding two χ's leaves the simplex (Σ=2); the legal blend is a CONVEX
combination (Σλ=1) — which is EXACTLY the flux affine constraint, but on a
different object. So χ is affine-simplex-constrained, νΣf is cone, and they are
multiplied together in fission: a SIMPLEX element times a CONE element. This is
why fission's `χ⊗νΣf` is a genuine rank-1 (the two factors are algebraically
different kinds), and why a single `CoefficientField` role would mistype χ.

RECOMMENDATION: `CoefficientField` is an ABSTRACT base (the commutative-algebra
machinery: pointwise `*`, scalar `·`, `+`, `_check_partner`). Two concrete role
leaves: `CrossSectionField` (cone, `1/cm` or `s/cm`) and `SpectrumField`
(simplex, dimensionless, with a `_check_normalized` invariant + convex-blend
`mix` verb mirroring `FluxRole.affine_combination`). The simplex/cone distinction
is the coefficient-side mirror of the flux/displacement distinction — both are
"illegal-states-unrepresentable" (Pattern 4) applied to the coefficient algebra.

**Elegance payoff.**
- *Structure-exposing*: makes χ's simplex constraint (the ERR-039 normalization
  class — `Σχ=1`, the `sr`/`/4π` bug family the `units.py` doc obsesses over)
  a TYPE invariant. A χ that fails Σ=1 is unconstructable.
- *Expressive*: homogenisation reads as `Σ_t_mix = sum(N_m * sigma_t[m] for m)`
  — a literal algebra sum (the `MaterialXSField` docstring's "Mixing" operation
  finally has a home and reads as math). χ blends read as `chi.mix([(f, chi_m)])`.
- *Structurally-simpler*: removes the temptation to reuse `FluxRole` for
  coefficients (which would forbid `σ+σ` — breaking homogenisation by
  construction). One base algebra, two constrained leaves.

**First test.** `CrossSectionField(σ) + CrossSectionField(σ')` succeeds and is
the pointwise sum (cone closed under +); `SpectrumField(χ) + SpectrumField(χ)`
must EITHER raise (simplex not closed under +) OR be permitted only via
`.mix(Σλ=1)`. Discriminates "coefficient = flux torsor" (wrong — would forbid the
σ sum) from "coefficient = algebra with a simplex-constrained sub-leaf" (right).
Pick the convention with the user; the test pins whichever is chosen.

**Structural attack on current.** `MaterialXSField` stores σ_t, νΣf, χ in
identically-shaped `(ng,*spatial)` ndarrays with NO algebraic distinction — χ's
simplex constraint and σ's cone constraint are invisible. A depletion/feedback
update that perturbs χ off the simplex (the classic normalization bug) is
undetectable until a downstream k_eff comes out wrong. The frame exposes that the
coefficient record conflates three different algebraic objects into one untyped
blob.

---

### Frame 3 — Kernel-vs-Coefficient: the locality criterion, and where fission's rank-1 sits

**Trigger.** A locality dichotomy: σ_t·ψ is diagonal (the operator's Schwartz
kernel is `δ(r−r')δ(Ω−Ω')δ(g−g')·σ_t`), Σ_s·φ integrates over angle (a genuine
nonlocal kernel `Σ_s(Ω·Ω')`). Skill A.7 (integral equations): "an integral
kernel with compactness / separability structure is present."

**Structural verdict — the criterion.** The Operator/Kernel split in §5.6 is the
**locality / diagonality** of the Schwartz kernel:

- **Coefficient → MultiplicationOperator** iff the operator's kernel is a
  multiple of the identity-measure (diagonal in EVERY phase-space axis): a
  zeroth-order PDO whose symbol is the coefficient. σ_t, σ_a, νΣf, 1/v, χ
  (as a per-cell broadcast). Test: `(A x)(p)` depends only on `x(p)` — point p
  maps to point p.
- **Kernel → IntegralKernelOperator** iff the operator integrates the field
  against a measure on at least one axis (nonlocal in that axis): Σ_s integrates
  over the OUTGOING angle (`∫ Σ_s(Ω·Ω') ψ(Ω') dΩ'`), energy-group transfer
  integrates over g'. Test: `(A x)(p)` depends on `x(p')` for p'≠p — there is a
  contraction.

The line is precisely "diagonal (symbol) vs off-diagonal (Schwartz kernel with
support off the diagonal)". `DiagonalOperator`/`MultiplicationOperator` are the
diagonal side; `ScatteringOperator` (the `SumOfTensorProductsOperator`
`Σ_ℓ Pℓ⊗Σ_{s,ℓ}` form) is the kernel side. The `DiagonalOperator.from_measure`
factory is the HINGE: a measure on one axis that is INTEGRATED is a kernel; a
measure that MULTIPLIES (the angular weight W in `Mφ`, NOT summed) is a
multiplication. Same weights, different verb — the measure is integrated-against
(Kernel) vs multiplied-by (Operator).

**Where does fission's χ⊗νΣf sit? — the decisive question.**

It is **(b): a composition `Multiplication(χ) ∘ Σ_g ∘ Multiplication(νΣf)`**, NOT
a primitive integral kernel. Decompose the action
`(Fφ)_g = χ_g · Σ_{g'} νΣf_{g'} φ_{g'}`:

1. `Multiplication(νΣf)`: pointwise `νΣf_{g'}·φ_{g'}` — a `MultiplicationOperator`
   (diagonal, coefficient = νΣf, a cone element). Flux → rate-density.
2. `Σ_g` (group-sum, = integrate against the COUNTING measure on the group axis):
   `Σ_{g'} (·)_{g'}` — this is a **Functional** in the §5.6 sense, field → scalar
   PER CELL (the production rate `p(r)`, a scalar field on space). It collapses
   the group axis. This is the `compute_production_rate` / `ProductionRateSolver`
   the #226 pyright pass just split out (memory: that split was a structural fix
   the user demanded). Suffix law: it ends in a contraction to a per-cell scalar,
   so it is `ProductionRateFunctional` — NOT an Operator.
3. `Multiplication(χ)`: broadcast the per-cell scalar `p(r)` across groups
   weighted by the simplex χ — a `MultiplicationOperator` whose coefficient is
   the SpectrumField χ. Scalar → rate-density (re-fills the group axis).

So `F = M_χ ∘ R_groupsum ∘ M_{νΣf}` where the middle is the production-rate
functional (or its reconstruction). The `RankOneOperator(χ, νΣf, axis=0)` IS this
composition collapsed into one matvec — it is a *separable LowRankKernel*
(§5.6 `SeparableKernel`/`LowRankKernel`), the rank-1 factorisation
`F = |χ⟩⟨νΣf|` on the group axis. Both readings are correct and they are
adjoint-dual: the rank-1 operator IS `M_χ ∘ (sum) ∘ M_{νΣf}`.

DECISION for the §5.7 placement: `FissionOperator` is an `IntegralKernelOperator`
(it integrates over g' — nonlocal in energy, exactly the §5.7 placement), whose
kernel is a `SeparableKernel`/`LowRankKernel` (rank-1 = `χ⊗νΣf`). `RankOneOperator`
becomes the L1 numerics primitive that BACKS the separable fission kernel — it
stays where it is (a tagged-axis rank-1), and `FissionOperator.kernel` keeps
returning it. The composition reading (b) is the SEMANTIC decomposition that the
production-rate Functional makes explicit; the rank-1 is the EFFICIENT realisation
(one einsum, no intermediate per-cell scalar array). They must be cross-verified
bit-identically.

Σ_s, by contrast, is NOT separable in general (it is `Σ_ℓ Pℓ⊗Σ_{s,ℓ}` — a
SUM of rank-(2ℓ+1) separable pieces, i.e. a genuine `LegendreMomentScatteringKernel`),
which is why it is a kernel and not a multiplication, and why it has no `solve`.

**Elegance payoff.**
- *Structure-exposing*: names the production rate as a `Functional` (field→scalar
  per cell), exposing WHY fission is rank-1-in-energy and power iteration
  converges geometrically (the `fission.py` docstring already argues this in
  prose — the type would carry it). Connects to the #226 `ProductionRateSolver`
  split (already validated by the user).
- *Expressive*: `F = M_χ @ production_rate @ M_νΣf` reads as the physics
  ("multiply by νΣf, integrate over groups to get the production rate, redistribute
  by the spectrum χ"). The current `singledispatchmethod` over four input types
  hides this behind dispatch arms.
- *Structurally-simpler*: unifies fission and scattering under
  `IntegralKernelOperator` with the SEPARABILITY (rank) as the discriminator —
  rank-1 (fission) vs sum-of-low-rank (scattering) — one axis of variation, not
  two unrelated classes.

**First test.** Build `F_composed = M_χ @ ProductionRateFunctional @ M_νΣf` and
assert `F_composed.apply(φ) == FissionOperator.kernel.apply(φ)` bit-identically on
a 3-group fissile cell. Discriminates the genuine composition reading from the
opaque rank-1 — if they disagree, the production-rate Functional is mis-specified
(e.g. wrong measure on the group axis). An experiment that cannot fail would be
"check F is rank-1" (it is, trivially) — the discriminating test is the
composition equivalence.

**Structural attack on current.** Fission lives as a `singledispatchmethod` over
{AngularFlux, TimedFullField, ScalarFlux, ndarray} with the rank-1 math buried in
a `kernel` property and the `1/W` projection scattered across arms. The production
rate — the physically central scalar — has no type; it is an intermediate
`np.einsum("gxy,gxy->xy")` inside `RankOneOperator.apply`. The frame exposes that
the most important diagnostic quantity in criticality (the per-cell fission rate)
is an anonymous einsum reduction, not a named `Functional`.

---

### Frame 4 — The honest generic: `V` is `TransportState`, not `Vector`

**Trigger.** A single genuine polymorphism remains (the cross-method state
carrier) typed by a structural Protocol whose NAME (`Vector`) describes its
algebraic capability, not its domain meaning — the user's "type hint that doesn't
help" complaint. Skill A.2 / elegance Pattern 1 ("code should read like the
domain").

**Structural verdict.** After coefficients are typed and operators promoted, the
ONLY thing that varies across SN/CP/MoC at the `apply(x: V) -> V` signature is the
STATE the method carries. The right structure is a **bounded `TypeVar` over a
`TransportState` Protocol**, NOT genuinely-parametric-per-method and NOT a bare
`Vector`. Reasoning:

- The algebra is an ENDOMORPHISM on ONE carrier (the `vector.py` docstring proves
  this: everything added in the SI/Krylov loop is the same type). So a single `V`
  per solve is correct — NOT per-operator generics (memory `issue_226`: the user
  already rejected per-op `Generic[F]`; the same logic forbids per-method here).
- But `Vector` names the WRONG thing. `Vector` says "I support +, −, scalar·, /".
  That is the CAPABILITY needed by the iteration drivers, true — but it sends a
  reader hunting ("a vector of WHAT?"). `TransportState` says "I am the discrete
  state of a transport solve" — which is what the carrier IS (`TimedFullField`
  for SN: bulk⊕boundary⊕history; a `ScalarFlux`-carrier for CP; a `RayField`
  carrier for MoC). The domain meaning, not the algebraic capability.

RECOMMENDATION: introduce `TransportState(Protocol)` in `orpheus/transport/`
(NOT numerics — it is a transport concept) that REFINES `Vector` (it still needs
the four vector ops for the drivers) AND names the transport-carrier contract
(the bulk/boundary structure the algebra's block-roles dispatch on — `TimedFullField`
already satisfies this). Then `V = TypeVar("V", bound=TransportState)` at the
TRANSPORT operator layer. `numerics.Vector` STAYS (it is the correct numerics-layer
abstraction — `np.ndarray` and bare `Field` leaves satisfy it, and numerics cannot
import transport per the layering note). The honest generic is a LAYERED pair:
`Vector` (numerics floor: "supports vector ops") refined by `TransportState`
(transport ceiling: "is a transport solve's state"). The SN/CP/MoC `solve`
signatures read `solve(self, ...) -> TransportState` / `apply(x: V) -> V` with
`V bound=TransportState` — and NOW the reader knows what flows through.

This is NOT a rename of `Vector` to `TransportState` (the brief's either/or) — it
is a REFINEMENT. `Vector` is genuinely needed at the numerics layer where
`np.ndarray` (the scipy wire format) must satisfy it and `np.ndarray` is NOT a
transport state. The two-level Protocol is the honest answer to "make the generic
read as the domain": the numerics layer keeps the capability name, the transport
layer adds the domain name.

**Elegance payoff.**
- *Structure-exposing*: `apply(x: TransportState) -> TransportState` at the
  transport layer exposes that the carrier is a transport state (bulk⊕boundary⊕
  history), not just "something addable". The #208 block-role dispatch
  (BULK/FULL/BOUNDARY) only makes sense on a `TransportState` (it reads `.bulk`/
  `.boundary`) — so the bound is load-bearing, not cosmetic.
- *Expressive*: the reader stops hunting. `V bound=TransportState` answers "a
  vector of what?" at the type site.
- *Structurally-simpler*: one carrier type per solve, layered cleanly across the
  numerics/transport boundary; no per-method genericity, no per-operator genericity.

**First test.** Type-check (pyright) that `np.ndarray` satisfies `Vector` but does
NOT satisfy `TransportState` (it has no `.bulk`/`.boundary`), and that
`TimedFullField` satisfies BOTH. Then assert at runtime (`@runtime_checkable`)
that `isinstance(timed_full_field, TransportState)` and
`not isinstance(np.ndarray(...), TransportState)`. Discriminates a genuine
refinement (transport layer is stricter) from a useless alias (`TransportState =
Vector`, which would let `np.ndarray` through and re-create the "vector of what"
problem). If `np.ndarray` satisfies `TransportState`, the bound is theatrics.

**Structural attack on current.** `apply(self, x: V, /) -> V` with `V bound=Vector`
is honest about the endomorphism but mute about the domain — exactly the user's
complaint. `Vector` is the floor; at the transport layer it is the wrong altitude.
The frame exposes that the #256 generic is correct in ARITY (one carrier) but
under-specified in MEANING (capability, not domain) at the layer where the carrier
is always a transport state.

---

## CROSS-METHOD POLLINATION

Current method class: SN (the operator-algebra retype is SN-led; the carriers
generalise to CP/MoC).

**From diffusion / Pn (the multiplication-operator SIBLINGS).** §5.7 puts
`TimeMassOperator` (= multiply by 1/v) as the OTHER `MultiplicationOperator` leaf
besides `CollisionOperator`. Trigger: 1/v is a cone coefficient with units `s/cm`
— a FIFTH coefficient unit signature, same algebraic structure as σ_t (cone,
pointwise, self-adjoint), DIFFERENT units. Borrowing: the `MultiplicationOperator`
+ `CrossSectionField` design MUST be built so `TimeMassOperator =
MultiplicationOperator(inverse_velocity_field)` falls out for free (alpha/transient
eigenvalue, §5.6 `TimeMassOperator T = TimeMassOperator(xs.velocity).on(V_sn)`).
First test: `T = MultiplicationOperator(one_over_v)` and `T.H == T` (self-adjoint),
`spec(T) == ess-range(1/v)` — same law suite as C, different coefficient. Confirms
the embedding is method/leaf-agnostic.

**From CP/MoC (the carrier generalisation).** Trigger: CP and MoC solve the same
`(L+C−S−F)ψ=q` algebra on a DIFFERENT carrier (CP: region-averaged `ScalarFlux`
+ collision-probability matrix as the resolvent; MoC: `RayField` along tracks).
The transport-resolvent backbone (promoted library: SN/MoC/CP `solve` are three
quadratures of `(Ω·∇+Σ_t)⁻¹`) says C = `MultiplicationOperator(σ_t)` is IDENTICAL
across all three — only the resolvent (the `solve`) differs. Borrowing: build
`CrossSectionField`/`MultiplicationOperator` in `orpheus/transport/` (NOT
`orpheus/sn/`) so CP and MoC consume the SAME C leaf. The `MaterialXSField`
docstring's "Restriction" operation (per-region subset) is the CP/MoC entry point
— a CP region-average is a `Projection` (§5.6) of the coefficient field onto the
region-indicator basis, adjoint to a `Reconstruction`. First test: a
`CrossSectionField.restrict(region_map)` that returns a coarser `CrossSectionField`,
and `Projection.H == Reconstruction` (the §4678 adjoint law) — verified on a
2-region homogenisation where the volume-weighted average is the known answer.

**From the eigenvalue layer (naming-collision fix).** Trigger: §4158
`CriticalityEigenproblem.multiplication_operator()` returns `A_loss⁻¹@F` (the
iteration operator) — colliding with the §5.7 `MultiplicationOperator` class.
Borrowing: rename §4158 to `iteration_operator()` / `resolvent_map()`. This is the
resolvent backbone again — `A_loss⁻¹F` is the power-method quadrature of the
transport resolvent; the `fix(step)` combinator (`power_iteration` /
`power_iterate_variant_alpha`) consumes it. First test: none (a rename); the
discriminator is that grep for `MultiplicationOperator` must NOT return the
eigenvalue verb after the rename.

---

## UNEXPLORED

- **Tensor networks / MPO** — no compositional-bond-dimension trigger on the
  coefficient algebra (M_f is rank-full diagonal, not a low-rank network). The
  rank-1 fission kernel IS a degenerate (bond-dimension-1) separable, covered
  under Frame 3's `SeparableKernel`, not a network.
- **Group/representation theory** — no symmetry trigger on the coefficient
  promotion itself (σ_t carries no group action; the SO(3) angular symmetry lives
  on the quadrature/scattering kernel, not the multiplication operator).
- **Differential geometry / connection coefficients** — no curvilinear trigger in
  this scope (the coefficient is a 0-form scalar field; the curvilinear
  redistribution lives in the streaming operator L, already framed elsewhere —
  `phase4_o2b_4_6_mms_ansatz_frame`).
- **Category theory (functor/embedding)** — Frame 1 IS the categorical reading
  (M is a faithful unital *-homomorphism, a functor from the commutative
  coefficient algebra to B(L²)); pinned there with concrete laws, so not listed
  separately as speculative.
- **Probability / Feynman-Kac** — the simplex structure of χ (Frame 2) touches it
  (χ is a probability vector, the emission law of a branching kernel), but the
  promotion target is a multiplication operator, not a Markov kernel; the MC
  `BranchingKernel` reading is a §5.6 MC concern, out of scope for the
  field→operator promotion.
- **Optimization / variational** — no variational trigger; the promotion is
  algebraic, not a minimisation. (Distinct from Smell #15 rank-non-monotone,
  which fires on Galerkin-without-a-principle — not present here.)
