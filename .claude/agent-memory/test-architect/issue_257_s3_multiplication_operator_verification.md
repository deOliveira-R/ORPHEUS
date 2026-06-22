---
name: issue-257-s3-multiplication-operator-verification
description: #257 Stage S3 verification spec — CollisionOperator → MultiplicationOperator(CrossSectionField) promotion + generalize DiagonalOperator to N-D broadcast engine; the multiplier-algebra law-suite, broadcast oracle, resolvent-assembly check
metadata:
  type: project
---

# #257 S3 — MultiplicationOperator (C = M[σ_t]) verification spec (PRE-IMPL)

Branch `feature/field-typed-operator-algebra`, HEAD `1ce727a`. Pre-carve
test-architect dispatch. **Fact: this is an operator-algebra carve crossing
numerics↔transport — gate philosophy is PRINCIPLED-EQUIVALENCE not forced
bit-id (user 2026-06-19).**

## What S3 builds
- numerics: generalize DEAD `DiagonalOperator` (`operator.py:1436`, ZERO prod
  callers) from "1-D weights on ONE axis" → N-D coeff on a SUB-PRODUCT of axes,
  broadcast over complement (`expand_dims(coeff,bcast_axes)*carrier`). Shared
  broadcast ENGINE; 1-D is the rank-1 special case.
- transport: `MultiplicationOperator(CrossSectionField)` — FIRST operator in
  `orpheus/transport/`; the §5.7 promotion `M:L^∞→B(L²)`; DELEGATES raw
  broadcast to numerics engine on `ψ.bulk.values`, owns the typed codomain
  (AngularSourceSink, ANGULAR_RATE_UNITS). Folds dead DiagonalOperator in.

## Seam (verified, do NOT re-explore)
- `CollisionOperator` `sn/operator.py:511-673`: apply `σ[None]*ψ.bulk.values`→
  AngularSourceSink+`BoundarySourceSink.zeros_on`; solve `q/σ[None]`→AngularFlux
  +BoundaryFlux.zeros_on (NO σ=0 gate, IEEE NaN); apply_transpose=apply (self-
  adjoint, `.H` is inherited metric-blind `_AdjointOperator` domain=None →
  EUCLIDEAN); `__add__(StreamingOperator)→InvertibleOperator(L,C)` resolvent.
- 3 prod sites ALL `CollisionOperator(sn_mesh, solver.mat_xs.total_cross_section)`:
  `solver.py:217/925/1000`. Rewire σ source → S2 `mat_xs.total_cross_section_field`
  (CrossSectionField, `.values IS` the same `(ng,*spatial)` array).
- `DiagonalOperator(weights:1-D,axis=0)`: `_reshape(ndim)*x`, self-adjoint,
  CAP_SOLVE iff all w≠0, `from_measure`. 2 TEST callers ONLY:
  `tests/numerics/test_diagonal_operator.py` (170 lines, the direct pin),
  `tests/numerics/test_tensor_product_operator.py` (DiagonalOperator&Diagonal
  Kronecker via `&`). NO prod callers (grep clean).

## Structurally-independent references (the L14 four-legged correctness, NOT L4)
- **kinf homogeneous** `tests/sn/verification/analytical/test_kinf_homogeneous.py`
  (L1, ≥2G — 1eg/2eg/4eg): dominant eigenvector+eigenvalue of `A⁻¹F` by PURE
  linalg on XS matrices, NO transport discretisation → mathematically orthogonal.
  THE eigenvalue reference (closed-form pillar). Must STAY green — it routes σ_t
  through C in the loss operator.
- **SI resolvent → analytical kinf** `test_invertible_operator.py::TestInvertibleSolveBridge
  Regression::test_si_carve_recovers_analytical_kinf` (homogeneous reflective,
  ng≥2 → kinf rtol<1e-9). THE resolvent-assembly structural ref (L+C built via
  promoted C). Must stay green.
- **L0 streaming-equilibrium** φ=Q/(Σ_t(1−c)) — flat Q+flat Σ_t. Fixed-source.

## Gate verdicts (the 5 spec items)
1. BROADCAST ORACLE: generalized engine on σ_t (bcast_axes=ordinate) ≡ legacy
   `sigma[None]*ψ` at VALUES level. Expect 0 ULP (SAME op — both `expand_dims`
   on the leading axis), but FRAME via 3-criteria not blanket array_equal: if
   it diverges, accept iff (a) named principled intermediate, (b) struct-indep
   ref, (c) drift ≤ reduction_depth×ULP (reduction_depth=1 — pure broadcast-
   multiply, NO sum). RECOMMEND `assert_array_equal` as the PRIMARY (it IS the
   same op) with a documented fallback to `nulp=1` only if a future
   re-association appears.
2. LAW-SUITE (intrinsic-properties directive): M_1=I, M_0=ZeroOperator
   (CODOMAIN-aware — collision's zero is a SOURCE not flux, precedent
   `ZeroOperator(codomain_zero=…)` `:918`), linearity M_{af+bg}=aM_f+bM_g (≥2G
   asym het — NOT 1G/homog, anti-pattern #3/#4 BLIND), self-adjoint M.H=M (real
   coeff, Euclidean since domain=None), spectrum spec(M_f)=ess-range(f) →
   CAP_SOLVE iff min|f|>0 (f=0 → solve revoked, apply NaN per IEEE), homomorphism
   M_f∘M_g=M_{f·g} VALUES-level (σ·σ' is cm⁻² — build M_{f·g} from the
   numerics engine on the raw product array, NOT a CrossSectionField since units
   differ; OR build a synthetic CoefficientRole product. The homomorphism is a
   numerics-engine property, test it THERE).
3. RESOLVENT-ASSEMBLY: `L + M[σ]` still → InvertibleOperator; sweep/matvec
   through promoted C matches legacy resolvent AND hits kinf/Q-equilibrium ref.
4. INVENTORY (L20): STAY-GREEN `test_collision_operator.py` (the 457-line legacy
   pin — RE-POINT onto MultiplicationOperator OR keep CollisionOperator as a thin
   alias/subclass; decide at impl), `test_invertible_operator.py`,
   `test_streaming_operator.py::TestSumCapabilities/TestOperatorAlgebraComposition`,
   `test_kinf_homogeneous.py`. MIGRATE `test_diagonal_operator.py` +
   `test_tensor_product_operator.py` onto the generalized engine (1-D = rank-1
   case must STILL pass — backward-compat). NEW: law-suite + broadcast oracle +
   resolvent.
5. CROSSWALK (L17 Pattern 7): legacy `σ[None]·ψ` over leading ordinate axis,
   codomain AngularSourceSink ↔ engine `expand_dims(coeff,bcast_axes)·carrier`
   on the value-tensor, transport retyping. The BRIDGE row = transport
   MultiplicationOperator delegates broadcast to numerics, applies codomain type.

## Failure-mode exposure
- 1G/homogeneous BLIND (L2, anti-pattern #3 — k flux-shape-indep, flat ψ nulls
  nothing here but multi-AXIS broadcast bug needs ≥2 spatial axes + ng≥2 to see
  axis-ordering). LAW-SUITE linearity/homomorphism MUST be ≥2G asymmetric het.
- The N-D generalization's NEW risk = wrong bcast_axes / axis-ordering (variable-
  swap mode #2) → ONLY visible when the carrier has ≥2 complement axes of
  DIFFERENT sizes (so a transposed broadcast raises or mis-scales). 2-D mesh
  `(N_ord, ng, nx, ny)` with nx≠ny + ng=2 is the discriminating regime.
- σ=0 spectrum boundary (missing-factor-adjacent): CAP_SOLVE revocation is the
  catcher; legacy CollisionOperator had NO gate (IEEE NaN) — promotion is a
  CHANCE to add the honest capability gate (Pattern 4). Flag to method-implementer.

## NO new vv failure mode (matrix unchanged). NO ERR until a real bug (next free
ERR-063 per #251 note). Extends [[feedback_principled_over_bit_identical]] +
[[feedback_test_intrinsic_properties]] + [[feedback_vv_tagging]] +
[[feedback_regression_tolerance_design]].
