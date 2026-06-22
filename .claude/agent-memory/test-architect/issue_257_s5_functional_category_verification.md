---
name: issue-257-s5-functional-category-verification
description: #257 Stage S5 verification spec — Functional(Protocol[V]) category (§5.6 suffix law) + ProductionRateFunctional + estimators-as-Functionals; the category-distinctness gate, hand-derived production-density correctness ref, estimator bit-identity
metadata:
  type: project
---

# #257 S5 — Functional category verification spec (PRE-IMPL, DELIVERED)

Branch `feature/field-typed-operator-algebra`, HEAD `b404ae1` (S4.5).
BIT-IDENTICAL additive stage; gate by 3-criteria, **0 ULP expected for
the reductions** (user). NO production code written (method-implementer
follows). Spec files LANDED + dry-run-validated (correct stub → 21p/1skip;
3 mutations all RED the right gates; stub reverted).

## SUT (PRE-IMPL — none exist yet, grepped clean)
- `Functional(Protocol[V])` — numerics-floor (L1) `orpheus/numerics/functional.py`,
  mirrors `vector.py`. ONE method `evaluate(self,x:V,/)->float|V` (§5.6
  suffix law). `@runtime_checkable`, `Generic[V]` V-bound-to-Vector.
  STRUCTURALLY NOT a LinearOperator (which needs apply+capabilities+
  domain+codomain — all 4, runtime-checked).
- `ProductionRateFunctional` (L3, lives `sn/`) — carries νΣf as
  CrossSectionField (S2 `mat_xs.fission_production_field` accessor
  `material_xs_field.py:408`). `evaluate(phi)->scalar-field` =
  `p(r)=Σ_g' νΣf_g'·φ_g'` group contraction = EXACTLY the `inner` line of
  `RankOneOperator.apply` (`operator.py:1776`):
  `(self.right*x).sum(axis=0,keepdims=True)`, right=νΣf. NOT rewired in S5
  (that's S6: `F=M_χ∘ProductionRateFunctional∘M_νΣf`).
- estimators `iteration.py:251-300` retyped: `_default_production_estimator`
  =`float(np.sum(F.apply(ψ)))`, `_default_keff_estimator`=Rayleigh
  `Σ(Fψ)/(Σ(Lψ)−Σ(Sψ))`. ARITHMETIC UNCHANGED. They are bare module-level
  `(L,S,F,ψ)->float` callables (KeffEstimator/ProductionEstimator aliases),
  injected into KEigenvalue (`:945`).

## Deliverable files (4) — ALL collect under `-O` (canonical), guarded by
`pytest.importorskip`/duck-probe so PRE-IMPL = skip-not-error.
- `tests/transport/_functional_helpers.py` — shared. `require()` (-O-firing
  pytest.fail), `require_functional()`/`require_production_rate_functional()`
  (probes sn.production_rate_functional / sn.fission / sn.functional),
  meshes (slab + nx≠ny 2-D discriminator), ASYMMETRIC ≥2G het νΣf+φ
  (distinct group-trend + distinct per-cell ramp → swap/wrong-axis NOT
  null), `build_production_rate_functional(nu_sigma_f=…)` (THE single
  construction assumption — change ONE func if layout differs),
  `hand_derived_production_density()` (explicit Python double-loop, NO
  numpy reduction = structurally-indep ref), `squeeze_density()`.
- `tests/transport/test_functional_category.py` — Spec A (foundation).
  POSITIVE (SUT IS Functional + has evaluate) + NEGATIVE both directions
  (Functional NOT LinearOperator: isinstance + lacks apply + lacks
  capabilities; RankOne/Mult/Identity NOT Functional) + DISCRIMINATOR
  (Frame-4: not-just-Vector via bare ndarray-is-Vector-not-Functional;
  not-just-LinearOperator; disjoint surfaces).
- `tests/transport/test_production_rate_functional.py` — Spec B (foundation).
  B.1 CORRECTNESS = hand-loop ref (2G-2D-nx≠ny + 4G-slab) + wrong-axis
  mutation discriminator; B.2 EQUIVALENCE (clearly demarcated NOT
  correctness) = 0-ULP `assert_array_equal` vs live RankOne `inner`
  (de-risks S6); B.3 Mode-3 = no-volume-measure (shape group-collapsed-only
  + value==unweighted hand-loop on non-uniform-V mesh).
- `tests/numerics/test_estimators_as_functionals.py` — Spec C (foundation).
  C.1 bit-identity arithmetic (synthetic L0 triple, hand floats 6.0 /
  6.0/39.6 — PASSES TODAY, 3 green) + C.2 category-honesty.

## ⭐ THE RECOMMENDATION (Spec C, brief asked): KEEP bare Callable
estimators, do NOT wrap. Rationale: estimators are `(L,S,F,ψ)->float`
(4-arg, consume the operator triple) — NOT `Functional[V]` whose surface
is the 1-arg field→scalar `evaluate(x)`. The honest Functional is the
field→scalar CORE `ψ→Σ(Fψ)` once F is bound. C.2 asserts bare callables
are NOT misclassified as Functionals (lighter-touch = category stays
honest because they never CLAIMED to be Functionals) + an OPTIONAL
`ProductionFunctional(F)` wrapper leg that SKIPS if not shipped (probes
`functional_mod.ProductionFunctional`). Method-implementer free to ship
the wrapper or not; spec green either way.

## ⭐⭐ TEETH PROVEN (dry-run mutation, the load-bearing validation)
- Mode-2 wrong axis (`sum(axis=1)`): 6 RED (both hand-loop legs +
  equivalence + both shape gates). ✓
- Category leak (add apply+capabilities): lacks_apply + lacks_capabilities
  + disjoint-surfaces RED. ⚠ the headline `isinstance NOT LinearOperator`
  gate does NOT fire on a PARTIAL leak (runtime LinearOperator needs ALL 4:
  apply+capabilities+domain+codomain) — fires ONLY on a FULL leak
  (re-verified: +domain+codomain → isinstance flips True → gate RED). The
  3 lacks/disjoint gates are the DEFENSE-IN-DEPTH that catches partial
  leaks (anti-pattern #11 multi-angle). DESIGN NOTE for method-implementer:
  a category leak via partial operator-surface is caught by the surface
  gates, not the isinstance gate.
- Mode-3 measure fold (`dens*2.0`): density_unweighted RED. ✓

## vv classification
- Claim layer (1.5 gate): A = CATEGORY/structural (Protocol membership,
  no number → ref IS the Protocol defs, non-self-referential via
  pos+neg+discriminator). B = FLUX-SHAPE/value (struct-indep ref = the
  hand double-loop, NOT the RankOne it'll replace — that's the demarcated
  EQUIVALENCE leg, L11 honoured). C.1 = regression/value (bit-id inherited
  from unchanged arithmetic). ZERO eigenvalue/MMS/convergence claims (those
  ride the EXISTING test_iteration KEigenvalue gates).
- Mode-8: ALL via `require`/`np.testing.*`, zero bare assert (canonical
  `python -O` — note `-O` is the INTERPRETER flag `python -O -m pytest`,
  NOT a pytest arg).
- NO new vv failure mode. NO ERR (next free ERR-063).

## Gate list (regression subset — all VERIFIED green @ baseline)
STAY-GREEN inheritance pins (RankOne untouched in S5 → fission byte-id):
`tests/sn/operators/test_fission_operator.py` (18p), `tests/numerics/
test_operator.py -k rank/RankOne` (subset), `tests/transport/
test_multiplication_operator.py` + `test_kinf_homogeneous.py` (≥2G NOT 1G)
+ `test_invertible_operator.py` (74p/2xf together). NEW: the 3 S5 files.
Baseline reds = 7 (#250 SPHERE ×5 + #232 mu_y ×2) — route around, never
all tests/sn (#212). Recommended select:
`pytest -O tests/transport/test_functional_category.py tests/transport/
test_production_rate_functional.py tests/numerics/test_estimators_as_functionals.py
tests/sn/operators/test_fission_operator.py tests/transport/
test_multiplication_operator.py`.

Extends [[issue-257-s3-multiplication-operator-verification]] (the sibling
promotion; same nx≠ny discriminator + 3-criteria framing) +
[[feedback_test_intrinsic_properties]] (the category-distinctness intrinsic
gate) + [[feedback_principled_over_bit_identical]] + [[feedback_vv_tagging]]
(foundation, no verifies).
