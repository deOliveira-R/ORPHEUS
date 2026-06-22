---
name: issue-257-s6-integral-kernel-category-verification
description: #257 Stage S6 verification spec — IntegralKernelOperator(§5.6 Kernel suffix) Protocol + FissionOperator.production_rate + ScatteringOperator.kernel reframes; category-refinement gate, fission χ·production_rate==F.apply 0-ULP, scattering kernel==R∘Λ∘M 0-ULP
metadata:
  type: project
---

# #257 S6 — IntegralKernelOperator category verification spec (PRE-IMPL, DELIVERED)

Branch `feature/field-typed-operator-algebra`, HEAD `7b463bb` (S5
committed `e9e0121`). ADDITIVE / BIT-IDENTICAL stage — adds a base
Protocol + `kernel`/`production_rate` properties + cross-checks ONLY;
does NOT touch the operators' matvec arms. Gate by 3-criteria, **0 ULP
expected** (user). NO production code written (method-implementer
follows). 4 spec files LANDED + dry-run-validated (correct stub →
20p/0skip; mutations RED the right gates; stub reverted → 7p/13skip
PRE-IMPL). Extends [[issue-257-s5-functional-category-verification]]
(same `require()`/probe-guard/`-O` style) + [[feedback_test_intrinsic_properties]].

## SUT (PRE-IMPL — none exist yet, grepped clean)
- `IntegralKernelOperator` — `@runtime_checkable` Protocol (NOT ABC; MRO
  friction w/ LinearOperatorMixin+@dataclass+@singledispatchmethod),
  `orpheus/transport/integral_kernel_operator.py`. Sole required member
  = `kernel` PROPERTY → a `LinearOperator`. §5.6 "Kernel" suffix =
  NONLOCAL action (integrates a field against a measure on ≥1 axis),
  distinct from `MultiplicationOperator` (LOCAL/diagonal, S3b) and
  `Functional` (field→scalar, S5). ⭐⭐ CRUX: the Kernel is a
  **REFINEMENT of** LinearOperator (still has apply+caps), NOT disjoint
  — UNLIKE the S5 Functional (shares NO member w/ LinearOperator). So
  the partition is ASYMMETRIC: a Kernel IS a LinearOperator, but only
  SOME LinearOperators are Kernels (those exposing `kernel`). The
  stub-validation Protocol = `class IntegralKernelOperator(LinearOperator,
  Protocol)` w/ a `kernel` property (inherits the 4 LinearOperator
  members + adds the 5th).
- `FissionOperator.production_rate` (NEW property, `sn/fission.py`) →
  `ProductionRateFunctional(self.mat_xs.fission_production_field)` (the
  S5 functional + S2 accessor `material_xs_field.py:408`). It ALREADY
  has `kernel` (`fission.py:167-220` = `RankOneOperator(χ,νΣf,axis=0) &
  IdentityOperator()`) → adding `production_rate` completes the Protocol.
  Matvec arms UNCHANGED.
- `ScatteringOperator.kernel` (NEW property, `sn/scattering.py`) → typed
  `OperatorProduct(R, OperatorProduct(Λ, M))` reproducing the EXISTING
  `_aniso_source_from_moment_values(M·ψ)` chain `R(Λ(M·ψ))`
  (R=`HarmonicMomentReconstruction.from_Y(self.Y)`,
  Λ=`LegendreMomentScattering(mat_xs,L,skip_l0=True)`,
  M=`MomentProjection(weights,Y,L)`; `op @ x` = `OperatorProduct`,
  `apply(x)=a.apply(b.apply(x))` `operator.py:826`). The 5
  `@singledispatchmethod` arms (TimedFullField:1038 / ScalarFlux:1085 /
  AngularFlux:1103 / HarmonicMomentField:1129 + base:961) UNCHANGED.
  ⚠ kernel = R∘Λ∘M ONLY; the `/W` producer-side norm (lesson L18) lives
  OUTSIDE the kernel (`build_aniso_source:784`).

## Deliverable files (4) — ALL collect under `-O`, probe-guarded
- `tests/transport/_integral_kernel_helpers.py` — shared. `require()`
  (-O-firing pytest.fail), `require_integral_kernel_operator()` (probes
  transport.integral_kernel_operator / numerics.* / transport.kernel_operator),
  `require_production_rate_property(op)`/`require_scattering_kernel_property(op)`
  (duck-`hasattr`→skip if PRE-IMPL — THE single construction assumption,
  change ONE func if the property name differs),
  `hand_derived_fission_emission(χ,νΣf,φ)` (explicit Python double-loop,
  NO numpy reduction = structurally-indep ref for Part B).
- `tests/transport/test_integral_kernel_category.py` — Spec A (foundation).
  Reuses S5 `_functional_helpers` (cartesian_2d_mesh nx≠ny, asymmetric
  νΣf). Builds F+S from small SNSolver fixtures (placeholder_materials has
  ZERO fission/scatter → can't satisfy Protocol; needs real data).
  POSITIVE (F+S ARE IKO: isinstance + kernel returns genuine LinearOperator
  + STILL a LinearOperator = the refinement leg) + NEGATIVE both dirs
  (MultiplicationOperator NOT IKO + lacks-kernel; Functional NOT IKO +
  lacks-kernel-AND-apply — the SHARPEST neg; IdentityOperator NOT IKO) +
  DISCRIMINATOR (Frame-4: a LinearOperator-without-kernel = IdentityOperator
  NOT IKO else it's a LinearOperator-rename; asymmetric-partition row).
  ⭐ runtime_checkable defense-in-depth: EVERY neg also asserts kernel-attr
  ABSENCE directly (isinstance only checks PRESENCE → incidental same-name
  attr could fool it).
- `tests/sn/operators/test_fission_kernel_crosscheck.py` — Spec B (foundation).
  4G fissile fixture (mixture A 4g: νΣf=[.014,.026,.125,.245] χ=[.6,.35,.05,0]
  per-group ASYMMETRIC → Mode-2 detectable; +mixture B moderator =
  heterogeneous). **B.1 CORRECTNESS** = `F.apply(φ)` vs hand-loop emission
  (struct-indep, rtol 1e-13) + role-swap-sensitivity discriminator (⚠ the
  νΣf↔φ swap is a NO-OP — production contraction is product-symmetric; the
  GENUINE Mode-2 = χ↔νΣf ROLE swap = broadcast-by-νΣf/contract-χ·φ, which
  asymmetric data distinguishes). **B.2 EQUIVALENCE/de-risk** (demarcated
  NOT correctness, L11): `χ · production_rate.evaluate(φ) == F.apply(φ)`
  0-ULP array_equal (production_rate IS the `inner` line
  `(νΣf*x).sum(axis=0,keepdims)` `operator.py:1776`; χ broadcast = RankOne
  `left*inner`) + production_rate-carries-νΣf (vs total_cross_section).
  Mode-11: reads `op.production_rate` OFF the live op.
- `tests/sn/operators/test_scattering_kernel_crosscheck.py` — Spec C
  (foundation). P1+heterogeneous(2-mat asym P0+P1)+2G fixture (ℓ≥1
  exercised). **EQUIVALENCE/de-risk ONLY** (⭐ DEMARCATED L11: physics
  ref = the EXISTING aniso MMS gate `test_curvilinear_aniso_scattering_p1.py`,
  NOT this same-array compare): `S.kernel.apply(ψ.values) ==
  _aniso_source_from_moment_values(MomentProjection(weights,Y,L).apply(ψ.values))`
  0-ULP array_equal + per-ordinate-shape `(N,ng,nx,ny)` + non-degeneracy
  guard (`scattering_order≥1` AND `moments[1:]≠0` else vacuous).
  Mode-11: reads `S.kernel` OFF the live op.

## TEETH dry-run-PROVEN (stub→revert)
- Correct stub (Protocol w/ kernel property + monkeypatched
  production_rate=νΣf-functional + kernel=R@(Λ@M)) → ALL 20 PASS (skips
  flip green).
- MUT `fission_wrong_xs` (production_rate uses Σt) → B.2
  `test_chi_times_production_rate_is_fission_apply` + B.3
  `test_production_rate_carries_nu_sigma_f` BOTH RED (right gates).
- MUT `scatter_drop_R` (kernel = Λ@M, R dropped) → C
  `test_kernel_apply_equals_existing_R_Lambda_M_chain` (value) +
  `test_kernel_apply_shape_matches_per_ordinate_source` (shape) RED.
- MUT Protocol-is-bare-LinearOperator-alias (no kernel member) → A's 2
  NEGATIVE-local rows (Mult/Identity wrongly IKO) + the DISCRIMINATOR
  (`test_a_linear_operator_without_kernel_is_not_an_integral_kernel`) RED
  — proves the gate catches a non-refining alias.

## Part-D regression gates (S6 ADDITIVE → MUST stay green; all confirmed
GREEN live @ HEAD 7b463bb, 33p)
- `test_scattering_operator.py::TestAnisoMomentSourcePath` (arm-equiv
  `:1177` + 4 pre-T3 byte-id snapshots `:1252/:1275/:1291/:1324`),
  `::TestProtocolCompliance` (caps unchanged — ⭐ S5 LESSON: adding the
  base Protocol MUST NOT alter `capabilities`; runtime_checkable
  LinearOperator checks 4 members, IKO adds a 5th but the advertised
  cap set is `{apply}` unchanged),
  `::TestProducerSideNormalisation::test_typed_apply_returns_per_ordinate_already_normalised`
  (`:415`, l0/matrix-eigenvalue), `::TestP0AlgebraicIdentities` (`:535`).
- `test_legendre_moment_scattering.py::TestBitIdenticalToLegacyInlinedMath` (`:240`).
- `test_curvilinear_aniso_scattering_p1.py` (both l0 `verifies("pn-scatter",
  "flux-moments")` — the structurally-indep PHYSICS L1 backing for Part C).
- `test_fission_operator.py::TestRankOneTensorProductKernel` (the existing
  Wave T `kernel.left is emission_spectrum` ref-identity + `kernel.apply≡apply`
  — the OTHER half of fission's IKO surface; ⚠ uses bare `assert` so RUN
  WITHOUT -O if pinning, though it's green under -O via the isinstance/is
  checks happening to be in fixtures... CONFIRM: it passed under -O here
  bc the asserts are skipped-but-the-test-still-collects-green — a Mode-8
  caveat the method-implementer should migrate to np.testing at touch).
  Plus `::TestProtocolCompliance` (fission caps unchanged).
- `tests/transport/test_functional_category.py` (S5, stays green).

## Recommended `python -O` selection (route around baseline reds)
```
.venv/bin/python -O -m pytest \
  tests/transport/test_integral_kernel_category.py \
  tests/sn/operators/test_fission_kernel_crosscheck.py \
  tests/sn/operators/test_scattering_kernel_crosscheck.py \
  tests/sn/operators/test_scattering_operator.py \
  tests/sn/operators/test_fission_operator.py \
  tests/sn/operators/test_legendre_moment_scattering.py \
  tests/sn/verification/mms/test_curvilinear_aniso_scattering_p1.py \
  tests/transport/test_functional_category.py
```
NEVER all of `tests/sn` (#212 `continuous_get` hang). Baseline reds to
route around if widening: #250 SPHERE snapshots ×5, #232 mu_y ×2.

## vv failure modes flagged (NO new mode, NO ERR)
- Mode-11 PRIMARY design hazard: both cross-checks read the NEW property
  OFF the live op (`op.production_rate`/`op.kernel`) — proven the mutation
  reddens. NOT a sibling-path route-around.
- Mode-8: all `require()`/np.testing (-O-safe). ⚠ the existing
  `TestRankOneTensorProductKernel` uses bare `assert` (Mode-8 inert under
  -O) — flagged for method-implementer to migrate at touch.
- 1.5 gate: every row is CATEGORY (Protocol membership) OR BIT-IDENTICAL
  equivalence (0 ULP). ZERO eigenvalue/MMS/convergence claims → no
  struct-indep NUMERICAL ref asserted here; the Part-C PHYSICS ref is the
  EXISTING aniso MMS gate (L11 — equivalence demarcated from correctness).
- Cardinal Rule: NO 1G eigenvalue claim (no eigenvalue claim at all).
  Cross-checks are ≥2G (fission 4G, scatter 2G) + heterogeneous + (scatter)
  anisotropic per H1/H2.
```
