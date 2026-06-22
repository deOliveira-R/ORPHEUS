---
name: issue-257-s6-integral-kernel-category-closeout
description: #257 S6 — ADD the §5.6 IntegralKernelOperator category (transport L2, a runtime_checkable REFINEMENT of LinearOperator exposing `kernel`) + reframe FissionOperator (new `production_rate` Functional) + ScatteringOperator (new `kernel` = R∘Λ∘M OperatorProduct); ADDITIVE + BIT-IDENTICAL (matvecs UNTOUCHED); 0 net-new pyright / 0 net-new type:ignore
metadata:
  type: project
---

# #257 S6 — the §5.6 IntegralKernelOperator category (FissionOperator + ScatteringOperator)

Branch `feature/field-typed-operator-algebra` 2026-06-20 (NOT committed,
on S5 HEAD `7b463bb`). ADDITIVE + BIT-IDENTICAL stage: ADD a base
Protocol + two properties + satisfy it; the operators' matvec arms stay
byte-identical. Extends [[issue-257-s5-functional-category-closeout]]
(S5 named the field→scalar Functional; S6 names the field→field NONLOCAL
Kernel — the §5.6 suffix-law completion).

## What landed (3 deliverables, all additive)

### 1. NEW `orpheus/transport/integral_kernel_operator.py` (L2, 0 pyright)
`@runtime_checkable class IntegralKernelOperator(LinearOperator[V], Protocol[V])`
— the §5.6 "Kernel" category: a `LinearOperator` whose action is NONLOCAL
(integrates a field against a measure on ≥1 axis), exposing its integral
structure via a `kernel` PROPERTY typed `-> LinearOperator`. ⭐ The
**inheritance form** `(LinearOperator[V], Protocol[V])` was chosen (brief's
preferred form) over a standalone Protocol — it STRUCTURALLY ENCODES
"is-a LinearOperator + has a kernel" and `@runtime_checkable` isinstance
then checks all 5 members (apply/capabilities/domain/codomain/kernel).
VERIFIED the partition is correct: `IdentityOperator`/`MultiplicationOperator`
(LinearOperator, NO kernel) → NOT IKO (strict refinement, not an alias);
`Functional` (no apply/kernel) → NOT IKO. Reused the INVARIANT `V`
(bound=Vector) verbatim from LinearOperator — a Kernel's `apply(x:V)->V`
both consumes+produces, so V is invariant-by-dual-use (NO variance
warning, UNLIKE S5's co-vector Functional which needed contravariant
`V_contra`+covariant `R_co`). Exported from `transport/__init__`.

### 2. `FissionOperator` reframe (`sn/fission.py`) — matvec UNTOUCHED
ADDED a `production_rate` property → `ProductionRateFunctional(self.mat_xs.fission_production_field)`
(the S5 functional + S2 typed accessor) = the §5.6 MIDDLE factor of
`F = M_χ ∘ ProductionRate ∘ M_νΣf` (Frame 3). It already had `kernel`
(`RankOneOperator(χ,νΣf,axis=0) & IdentityOperator()`, Wave T) → adding
`production_rate` completes the IKO Protocol. `kernel` KEPT as-is (the
efficient realization; user Q1). Docstrings (module + property) document
the semantic decomposition + the 0-ULP bit-identity (`production_rate.evaluate`
IS the RankOne `inner` line byte-for-byte; `χ·inner` = `left*inner`). The
4 dispatch arms + `kernel` UNCHANGED. **0 net-new pyright (8=8 baseline).**

### 3. `ScatteringOperator` reframe (`sn/scattering.py`) — 5 arms UNTOUCHED
ADDED a `kernel` property → `OperatorProduct(R, OperatorProduct(Λ, M))`
reproducing the EXISTING `_aniso_source_from_moment_values(M·ψ)` chain
`R(Λ(M·ψ))` byte-for-byte (R=`HarmonicMomentReconstruction.from_Y(self.Y)`,
Λ=`LegendreMomentScattering(mat_xs,L,skip_l0=True)`, M=`MomentProjection(weights,Y,L)`
— the SAME factors `_aniso_source_from_moment_values`+`build_aniso_source`
build). **DEFAULT (lower-risk) chosen: `kernel` is a SEMANTIC property +
cross-check; the 5 `@singledispatchmethod` arms NOT routed through it**
(guaranteed byte-identical; the brief authorized routing ONLY with provable
0-ULP, kept untouched per "any doubt → keep arms untouched"). Added a
`ValueError` guard for `scattering_order==0` (Pattern 4 — an isotropic
operator has NO anisotropic kernel; Y is None). Module + property
docstrings document the §5.6 Kernel reading + the `/W` lives OUTSIDE the
kernel (L18). **0 net-new pyright (22=22 baseline).**

## THE PYRIGHT HINGE (the load-bearing lesson, 0 type:ignore)

The scattering `kernel` `OperatorProduct(R, OperatorProduct(Λ, M))` is the
ONLY non-trivial typing battle. The numerics.projection primitives
(`MomentProjection`/`HarmonicMomentReconstruction`/`LegendreMomentScattering`)
subclass the **UNPARAMETRISED** `LinearOperatorMixin` with a CONCRETE
`apply(x: ndarray) -> ndarray` (NOT the generic `apply(x:V)->V`), so
pyright CANNOT unify the generic `V` of `OperatorProduct[V].__init__(a,b:
LinearOperator[V])` from them → 3 `reportArgumentType` "X cannot be
assigned to parameter a/b of type LinearOperator[V]". ⭐ **The `&`
(`__and__`) in fission's `RankOneOperator(...) & IdentityOperator()` types
cleanly at 0 errors DESPITE RankOneOperator ALSO being unparametrised — the
difference is the SECOND operand: `IdentityOperator` IS `LinearOperatorMixin[V]`
(PARAMETRISED), so V unifies; a `RankOneOperator @ RankOneOperator` between
TWO unparametrised mixins ALSO reddens (proven).** Tried + REJECTED:
(a) the `@` operator form `R @ Λ @ M` (cleaner Pattern-1 spelling) — the
`__matmul__` is on the MIXIN not the `LinearOperator` Protocol, so
annotating operands as the bare Protocol LOSES `__matmul__`; and the
unparametrised-mixin operands still don't resolve; (b) `LinearOperator`
bare annotation = `LinearOperator[Unknown]` → "capabilities not read-only"
+ "apply variance" assignment errors; (c) `LinearOperator[Vector]`
annotation → the projection ops' `cached_property` domain/codomain don't
match the Protocol's `property`. **FIX = `cast(LinearOperator, op)` on each
of the 3 factors** — the PEP-484 bridge for a known-correct runtime
interface the static analyser cannot infer (the runtime objects ARE
LinearOperators — `isinstance(kernel, LinearOperator)` passes in the test).
NOT a `# type: ignore` suppression (cast asserts the type, doesn't silence
an in-place error); it is the ESTABLISHED codebase convention
(`orpheus.data.micro_xs` uses `cast` in 5 files). The projection-operator
generic gap is #226 scope. Also narrowed `Y = self.Y; if Y is None: raise`
(Pattern 4) so the `from_Y(Y)`/`MomentProjection(Y=Y)` calls avoid the
`ndarray|None` arg errors my code would otherwise have ADDED (the existing
`_aniso_source_from_moment_values` does NOT narrow and inherits those as
baseline noise).

## GATES (all measured live, this session)

- **pyright** `npx pyright` (CLI oracle, NOT the IDE stream): **2311
  errors / 19 warnings = EXACT baseline, 0 net-new, 0 net-new `# type:
  ignore`.** ⚠ THE BASELINE IS 2311, NOT the brief's 2307 — the 4 S6 spec
  test files are untracked-but-on-disk and pyright analyzes them (the 2
  cross-check files carry 2 errors each = 4 pre-existing in the SPEC, which
  the method-implementer must NOT touch). Measured the real `7b463bb`-with-
  spec-files count BEFORE writing (same staleness pattern S5 flagged). New
  file `integral_kernel_operator.py` **0/0 standalone**; `fission.py` 8=8;
  `scattering.py` 22=22; `transport/__init__` 0=0. Per-signature
  line-independent diff confirmed NET-NEW=NONE in both touched files.
- **S6 suite** `python -O` 3 spec files: **20 passed / 0 skipped** (was
  7p/13skip PRE-IMPL — all skips flipped to passes).
- **Bit-identity regression** (matvecs untouched) + S6 suite + layer
  imports, ONE selection: **409 passed** (`test_scattering_operator.py`
  TestAnisoMomentSourcePath/TestProtocolCompliance/TestProducerSideNormalisation/
  TestP0AlgebraicIdentities + `test_fission_operator.py` incl.
  TestRankOneTensorProductKernel + `test_legendre_moment_scattering.py`
  TestBitIdenticalToLegacyInlinedMath + aniso MMS
  `test_curvilinear_aniso_scattering_p1.py` both l0 + S5
  `test_functional_category.py` + `test_layer_imports.py`). ZERO baseline
  reds introduced; routed around #250 SPHERE×5 + #232 mu_y×2 + #212 hang
  (never ran all of tests/sn).
- **`import orpheus`** clean; layer-imports green (`integral_kernel_operator.py`
  imports ONLY `orpheus.numerics` L1 + stdlib — clean L2, NO L3).
- **Sphinx** `sphinx-build` SUCCEEDED, NO warnings/errors (only pre-existing
  test-file SyntaxWarnings via the graph builder, not rST). NEW
  `.. _integral-kernel-category:` anchor + `:label: integral-kernel-category`
  (the `(Aψ)(x)=∫K(x,x')ψ(x')dμ` math) + `vv-status documented` + archivist
  `.. todo::` (#257 S6) in `operator_algebra.rst`, present in HTML
  (`id="integral-kernel-category"` rendered, V&V matrix picked up the
  documented-status label). Nexus graph picked up the `IntegralKernelOperator`
  class node + `kernel` method node.

## MUTATION-verified (the cross-checks have teeth, Mode-11)

- Scattering kernel DROP-R (`OperatorProduct(Λ,M)` only) → Spec C value
  gate `test_kernel_apply_equals_existing_R_Lambda_M_chain` + shape gate
  `test_kernel_apply_shape_matches_per_ordinate_source` BOTH RED.
- Fission `production_rate` wired to TOTAL XS (in-process monkeypatch) →
  B.2 `χ·production_rate==F.apply` (composed≠fused) + B.3
  `production_rate carries νΣf` (density≠Σνσf·φ) BOTH RED.
- ⚠ LESSON-FROM-A-NEAR-MISS: ran the FIRST scattering mutation via
  `git checkout -- orpheus/sn/scattering.py` to revert — which WIPED my
  edits (scattering.py is TRACKED at `7b463bb`, so checkout = revert to
  HEAD, NOT to my on-disk WIP). Had to re-apply all 3 scattering.py edits.
  **For mutation-testing an UNCOMMITTED tracked file, monkeypatch
  in-process (the fission MUT) or copy to a tmp, NEVER `git checkout`.**

## DEVIATIONS / decisions

- **Scattering arms NOT routed through `kernel`** (DEFAULT lower-risk;
  brief authorized either). The mild Cardinal-Rule-2 duplication (kernel
  parallel to `_aniso_source_from_moment_values`) is the same
  semantic-vs-realization pattern as fission, cross-checked at 0 ULP. Both
  share the SAME R/Λ factor types; the routing-to-single-source is a clean
  follow-up if desired but risked the 5-arm byte-identity.
- **NO algebra-of-record Branch-1 SymPy / L1-derivation manifest** — this
  is a type-surface carve naming an existing category + numerics-primitive
  delegation, NOT a new reference solver (same posture as S3a/S3b/S5; the
  fission cross-check's struct-indep ref is a TEST-SIDE Python double-loop
  `hand_derived_fission_emission`, the scattering physics L1 ref is the
  EXISTING aniso MMS gate). NO Branch-1/L1-xverif-in-derivations owed.
- **`kernel`/`production_rate` are bit-identity equivalence claims (0 ULP),
  NOT new correctness claims** (1.5 gate, L11 demarcated): the fission
  CORRECTNESS ref is B.1's hand-loop; the scattering physics is the aniso
  MMS gate.
- NO new ERR (the `scattering_order==0` guard is PROACTIVE Pattern-4
  illegal-state, not a found bug).
- NOT committed (main agent commits after elegance + qa review). archivist
  DISPATCH emitted (followup:false). Deferred follow-ups (file as issues):
  (a) un-orphan `SumOfTensorProductsOperator` for the full `Σ_ℓ Pℓ⊗Σ_{s,ℓ}`
  scattering form (principled-not-bit-identical re-derivation); (b)
  relocate carrier-agnostic C/F/S cores to `transport/` (needs CP/MoC
  carrier unification first, Pattern 6).

## LESSON (the durable one)

A §5.6 Kernel is a runtime_checkable REFINEMENT of LinearOperator
(`(LinearOperator[V], Protocol[V])` + a `kernel` property) — ASYMMETRIC to
the S5 Functional which is a DISJOINT sibling (shares no member); reuse the
INVARIANT `V` verbatim (a Kernel's `apply(x:V)->V` is invariant-by-dual-use,
NO variance warning, UNLIKE the co-vector Functional). When composing
existing unparametrised-`LinearOperatorMixin` primitives into an
`OperatorProduct[V]`, pyright cannot unify `V` unless at least ONE operand
is a PARAMETRISED `Mixin[V]` (this is why fission's `RankOne & Identity`
types clean but `R @ Λ @ M` between three unparametrised projection ops
does not, and the `@`-operator spelling does NOT help — `__matmul__` is on
the mixin, lost when annotating as the Protocol); the principled 0-ignore
fix is `cast(LinearOperator, op)` (the PEP-484 known-correct-runtime bridge,
NOT a suppression — the established `orpheus.data.micro_xs` convention) on
each factor, plus a `Y is None` Pattern-4 narrow so the `from_Y`/`Y=Y`
calls avoid the `ndarray|None` arg errors. NEVER `git checkout` an
uncommitted TRACKED file to revert a mutation — monkeypatch in-process.
Extends [[issue-257-s5-functional-category-closeout]] (S5=Functional/disjoint
sibling; S6=Kernel/refinement; both additive-bit-identical + 0-net-new).
