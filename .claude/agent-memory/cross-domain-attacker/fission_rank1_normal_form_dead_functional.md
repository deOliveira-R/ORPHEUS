---
name: fission-rank1-normal-form-dead-functional
description: Fission F=|χ⟩⟨νΣf| is the single-mode frame-conjugate NORMAL FORM (S6 "unfold" is structurally empty); the Functional category is correctly modeled but its ONLY typed instance is production-dead while 5 live functionals stay untyped — a new "abstract-category-built-top-down, instance-population-bottom-up-and-disconnected" smell.
metadata:
  type: project
---

# Fission rank-1 normal form + the dead-instance-of-a-live-category smell

Frame-attack on the `ProductionRateFunctional` fate decision (branch
`refactor/operator-inverse-algebra`, #257 S5/S6). Two durable structural
findings, one a confirmed normal-form theorem, one a NEW detection smell.

**Why:** stress-test of "are FissionOperator (fused rank-1) and the
IntegralKernelOperator-with-ProductionRateFunctional (decomposed) competing
implementations?" — the answer reshapes how the Functional category should be
populated.

**How to apply:** fire BOTH findings on any `R∘Λ∘M` / co-vector / suffix-law
fate decision. Reconcile file:lines vs the live worktree (L-005) before acting.

## Finding 1 — the rank-1 frame-conjugate NORMAL FORM (a theorem, branch-verified)

`frame.conjugate(identity)` for a SINGLE-MODE frame IS the `RankOneOperator`.
Proof (byte-for-byte, not gesture): `RankOneOperator.apply` (operator.py:1811)
is `inner = (right·x).sum(axis,keepdims); return left·inner`. With left=χ,
right=νΣf, axis=0 the `inner` line is CHARACTER-IDENTICAL to
`ProductionRateFunctional.evaluate` (production_rate_functional.py:151
`(νΣf·φ).sum(axis=0,keepdims=True)`) and `left·inner` is the M_χ broadcast. So
`R∘I∘M` (single mode, Λ=id) collapses to the outer product with NO residual
structure → the fused rank-1 ALREADY IS `M_χ ∘ ProductionRate ∘ M_νΣf` in
normal form.

CONSEQUENCE: the S6 "rewire F's kernel through ProductionRate" step is
STRUCTURALLY EMPTY — there is no third object between "fused" and "decomposed";
the decomposition is the permanent READING of the normal form. A docstring that
promises a pending S6 unfold (fission.py:52-77) mis-states a no-op as future
work. DISCRIMINATING first test: `np.array_equal(F.kernel.apply(φ),
hand_unfold.apply(φ))` (0 ULP, NOT allclose) — bit-identity ⟹ S6 buys nothing.

This is the general move for ANY "fused operator vs its named-factor
decomposition" competition: the fused primitive is the frame-conjugate normal
form; the factors are the reading. "Which realizes the operator" is settled by
the normal form, NOT by benchmark. Fission is the rank-1 (bond-dim-1) DEGENERATE
of scattering's R∘Λ∘M — same theorem, N=1.

## Finding 2 — NEW smell: "dead typed instance of a live untyped category"

The §5.6 Functional category (functional.py:54-68) is CORRECTLY modeled — a
disjoint Protocol (`evaluate:V→R`, shares no member with LinearOperator), with
an intrinsic-property gate. But it has EXACTLY ONE typed instance
(`ProductionRateFunctional`) and that instance is production-DEAD (grep: the
`Functional` Protocol is imported nowhere in production `orpheus/`; the instance's
only consumer is `FissionOperator.production_rate`, read only by tests).
MEANWHILE five genuine functionals ship UNTYPED:
- `Solution.reaction_rate_density` (solution.py:292, ⟨σ,φ⟩ fiberwise) — LIVE
- `SNSolver.compute_group_production_rate` (solver.py:1070, ∫νΣf·φ dV) — LIVE
- `SNSolver.compute_group_absorption_rate` (solver.py:1120) — LIVE
- `_default_production_estimator` (iteration.py:272, Σ(Fψ)) — LIVE
- `_default_keff_estimator` (iteration.py:293, Rayleigh quotient) — LIVE
  (bare `Callable` aliases `ProductionEstimator`/`KeffEstimator`).

THE SMELL (distinct from frame-leak naming L-006, which is one slot
over-committed to one of N consumers): a category built TOP-DOWN from one
sub-structure (the fission decomposition's middle factor) and never connected to
the instance population that already exists BOTTOM-UP. The category is right; the
ONE typed inhabitant is the one nobody calls; the real inhabitants are untyped.
TELL: an abstract Protocol whose `isinstance` discriminator is exercised only by
its own category-partition gate (test_functional_category), never by a consumer
that branches on the Protocol. Promotion-status: held for a 2nd sighting before a
Part C slot (count = 1). Fire inline meanwhile.

## The fate verdict (option 2)

RETIRE the dead `ProductionRateFunctional` + `FissionOperator.production_rate`;
the 0-ULP `test_fission_kernel_crosscheck` is a PROCEDURAL TWIN (vv L11 — both
sides the same numpy primitive, νΣf identical on both, asserts no physics) so it
does NOT migrate as a correctness contract (delete/downgrade). RE-SEAT the
Functional Protocol on the live population: `ReactionRateFunctional(σ)` of which
production (σ=νΣf·V) and absorption (σ=Σa·V) are instances (the existing
`compute_group_*_rate`). GENERALIZE toward a `BilinearFunctional(A)`
evaluate(φ†,φ)=⟨φ†,Aφ⟩ — `_default_keff_estimator` IS its φ†=1 degenerate, and
the homogenization/condensation Petrov-Galerkin campaign's untyped
adjoint/flux-weighted inner products are the LIVE/IMMINENT consumer that makes it
load-bearing (the type the dead instance should have been).

THE CRUX (genuine vs mirage structural independence): the per-term analytic gain
(check production & absorption each vs k∞=νΣf/Σa SEPARATELY, catching a νΣf-only
shared-factor error the k-ratio masks) IS genuine structural independence — but
it does NOT need `ProductionRateFunctional`; it needs the LIVE
`compute_group_production_rate`/`absorption_rate`, checked against their OWN
closed forms. The named-co-vector IDEA is right; the specific dead instance
contributes nothing to it. The carve crosses ≥3 subsystems (touches compute_keff
live math) ⟹ test-architect verification plan BEFORE; ship the per-term
analytic-vs-k∞ test as the structural oracle replacing the retired twin.
