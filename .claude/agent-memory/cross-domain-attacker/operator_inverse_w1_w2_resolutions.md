---
name: operator-inverse-w1-w2-resolutions
description: Resolutions of the #226 operator-taxonomy self-audit weak points W1 (moment-emit codomain crack) and W2 (fold-config identity crack), plus three grain-sharpenings of the §13 name-earning bar. Companion to operator_inverse_taxonomy_frames. Branch-verified refactor/inverse-as-operator.
metadata:
  type: project
---

# #226 taxonomy — W1/W2 self-audit resolutions + §13 grain-sharpenings

Attack of `.claude/plans/operator_machinery_taxonomy.md` §16 weak points, on
branch `refactor/inverse-as-operator` (branch-read, L5 — Nexus graph was stale on
`feature/sn-adjoint-transport`). Companion to [[operator-inverse-taxonomy-frames]]
(the layer below: substrate/operator split, four inverse strategies). These are
the NEXT layer: the two cracks the authors themselves flagged, resolved.

## W1 — moment-emit is `P ∘ A⁻¹`, NOT "a SweepOperator with moment emit"

Native frame: **deforestation/loop-fusion of a composition through a coisometry.**
The moments-emitting object maps FullField→MomentField = `P ∘ (L+C)⁻¹`. Typing:
- IDENTITY = `OperatorProduct(P, A.inverse())`. `P = scattering_op.frame.analysis ⊕
  Id_∂` is a BLOCK COISOMETRY (`AnalysisOperator[AngularFlux,HarmonicMomentFlux]`
  on the bulk summand of `C¹=C¹_int⊕C¹_∂`, identity on the trace — the reflective
  B needs the full per-ordinate trace). Rank-deficient (N→(L+1)², `M∘R=I`,
  `R∘M`=projector) ⇒ NOT invertible, NO round-trip promise. Domain=FullField(q),
  codomain=`TimedFullField(bulk=HarmonicMomentFlux, boundary=full trace)`.
- FUSION (never materialize ψ) = SUBSTRATE application-context (§38): the base
  `SweepOperator` sweeps in `_SweepEmit` MOMENT mode (`_SweepEmit` ALREADY EXISTS,
  loss_representation/__init__.py:659). `A.inverse().apply(q, emit=MomentEmit(P))`
  — `emit` is a substrate kwarg like `initial_guess`. The SweepOperator stays
  ENDOMORPHIC + S-direct; the codomain change is owned by the P FACTOR.
- FACTORY = the SOLVER (`_maybe_window`), building `P @ A.inverse()`. NOT an
  `A.inverse_into(frame)` on the operator (repeats the cross-reach).

Honest invariant (the plan's "residual of WHAT, without materializing ψ?"): the
composite CANNOT certify `(L+C)ψ=q` (P not injective). Demanding a moment-proxy
residual `M[(L+C)Rφ]−Mq` is CATEGORY-CONFUSED — `Rφ≠ψ`, so that is the residual
of the DIFFERENT P_N-reduced system and would RED a perfect fused sweep. The
contract FACTORS: (fusion-correctness = the existing windowed≡post-projection
oracle, ≤4 ULP, RENAMED as deforestation) ∧ (S-direct on the base A⁻¹, its own
full-field round-trip) ∧ (coisometry `M∘R=I` on the frame). Three pieces, each on
its owner.

Structural attack: plan §7 ("SAME SweepOperator, emit-policy set to moments")
VIOLATES the plan's own §1 two-layer law — a substrate emit-config that changes
the OPERATOR codomain. The substrate emit is honest; the crack is letting a
SweepOperator silently change codomain while calling itself `A.inverse()`. The
`solve_moments(rhs, frame)` signature (streaming.py:862) is the confession — an
inverse needing a FOREIGN operator's (S's) projection is a product with P as an
un-named left factor. Evidence: `solve_moments` lives 3× (InvertibleOperator,
_GaussSeidelResolvent, _MomentWindowedResolvent-delegated) = Smell #16 shape-1;
all 3 + the `_MomentWindowedResolvent` CLASS dissolve into ONE product.

NEW SMELL candidate (Part C, 1st sighting) **"codomain-changing emit/config"**: an
OUTPUT-mode config (emit policy, "windowing") presented as a config of ONE
morphism while it silently changes the CODOMAIN ⇒ it is a DIFFERENT morphism (a
composition), not a config. TELL: a `solve`/`apply` variant taking a foreign
object's projection as a parameter. FIX: lift to an explicit product factor at the
operator layer; keep the emit purely substrate. Distinct from vestigial-forward
(that is a fake forward; this is a real morphism mistyped as a config). Held for a
2nd sighting.

## W2 — the G-S resolvent IS `M.inverse()`, `M = (L+C−B_lower)` a REAL forward op

Native frame: **regular matrix splitting `A=M−N`** (Varga). `(L+C−B) = M − N`,
`M = L+C−B_lower` (strictly-lower part of the reflective coupling in octant-
schedule order — STAYS triangular ⇒ direct sweep), `N = B_upper` (lagged). B has
NO octant-diagonal (μ→−μ maps each octant to a DIFFERENT octant) ⇒ `B=B_lower+
B_upper` exact + diagonal-free. Object = `M.inverse()` = a genuine SweepOperator on
M, S-direct, NO vestigial apply. Iteration `ψ_{k+1}=M⁻¹(q+Nψ_k)` lives in the
DRIVER.

Route verdict: REIFY `M.apply = (L+C).apply(ψ) − B_lower.apply(ψ.trace)` (route a),
NOT "type it as a preconditioner" (route b). Route b repeats the original sin
(naming by driver-role): the object is not an APPROXIMATE inverse of `(L+C−B)`, it
is the EXACT inverse of `M`. "Preconditioner" = 3rd mislabel after "resolvent". The
reification reuses EXISTING pieces (loss_action + the `SweepSchedule.gauss_seidel`
`_solve_scheduled` already builds + `_reflect_outflow_into_inflow`) — "two views,
one substrate" done right; B_lower = B masked to the schedule's strictly-lower
octant pairs. It is the fuller-view-oracle exception (a forward that pins the
inverse's keystone), not speculative.

Replacement gate (ranked): (1) round-trip `M.inverse().apply(M.apply(x))==x` to
MACHINE precision on a NON-EMPTY-B_lower config (≥2 reflective faces) — the
runnable keystone; current object REDs against `(L+C).apply` (off by B_lower). (2)
fixed-point-equivalence vs Jacobi SI on a DIAGONAL (non-axis-aligned) cubature
with shared faces (Mode-9, breaks degenerate coincidences) — DRIVER-level. (3)
splitting-law `ρ(M⁻¹N)<1` — convergence certificate, not correctness.

## Three grain-sharpenings of the §13 "name = promise backed by test" bar

1. **§11.3 (seed-as-context) ⟺ §13 S-direct — SAME claim, SAME scope.** "seed
   doesn't change the result" IS "solve is direct/seed-independent". So §11.3 is
   NOT universal — it inherits S-direct's geometry-split (slab/cyl TRUE, sphere
   FALSE until #282). Pre-#282 the sphere `.apply(q, initial_guess=s)` is a 2-arg
   RELAXATION step `Φ(s;q)` (one Picard iterate → a `GreenOperator`, L-007), not a
   1-arg inverse. The plan xfails the GATE but leaves the RULING TEXT universal — a
   pass-looking ruling over a failing gate = the plan promising what it cannot do.
   Test: `sphere_inv.apply(q,seed_A)!=…(q,seed_B)` off the FP (diag probe 4.6e-2)
   REDs the universal text; slab/cyl pass.
2. **`as_matrix()` partiality is RESOURCE (total-map-with-effect), NOT structural
   (restriction idempotent).** `.inverse()` partial = a singular op HAS no inverse
   (`is_invertible` = restriction idempotent, structural sub-object). `as_matrix()`
   = a TOTAL functor Op→Mat (every op HAS a matrix) with a resource EFFECT
   (`MatrixTooLarge`, like OOM); the matrix EXISTS, just not materializable. Do NOT
   add an `is_materializable` idempotent (reifies a resource bound as structure) —
   "will it fit" is a SIZE precheck (`dom.dim*cod.dim>thresh`), no structure read.
   Restriction categories (Cockett-Lack) is the ONE place CT fires concretely here.
   Test: identical op type at 2 sizes — `is_invertible` SAME (structural), `as_matrix`
   availability FLIPS with SIZE ⇒ the idempotent framing REDs.
3. **`ResolventOperator(z)` is a FACTORY-level name, not an instance type.** R-identity
   `R(z)−R(w)=(z−w)R(z)R(w)` is BINARY (the z-map); Sweep/Green/Matrix invariants are
   UNARY. A single instance = just `(A−zI).inverse()` (a Sweep/Green/Matrix on the
   shift), earning only "inverse of the shift". The name belongs to the FACTORY
   `A.resolvent(z)=(A−zI).inverse()` (which CAN run the family test). DROP the
   standalone class; "resolvent" is a parametrized family of the existing 3, not a
   4th strategy. §13 conflated instance-grain and family-grain invariants.

Cross-method (durable): the `P∘A⁻¹` deforestation recurs in MoC (track-sweep +
moment projection) ⇒ `emit` belongs on the SHARED inverse contract, 2nd witness.
Diffusion is the negative control (no angular P ⇒ genuinely endomorphic inverse)
⇒ `OperatorProduct(P,inverse())` belongs at `transport/`, not `numerics/`.
