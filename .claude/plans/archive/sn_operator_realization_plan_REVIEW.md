> **SUPERSEDED by operator_strategy_realization_campaign.md §10 — archived 2026-08-19 (plans triage).**

# Review — `sn_operator_realization_and_posing_plan.md`, grounded against `main` @ `8654d348`

**Reviewer:** main agent, 2026-07-28. **Method:** full read of the plan (1504 L), then three
parallel `explorer` grounding passes (maps retained at `scratch/review_map_{reaction_operators,
sn_assembly,agnostic_diffusion}.md`), plus main-agent independent verification of every claim
this review calls TRUE (L12 discipline). Sphinx rebuilt on `main` first so the Nexus graph
answers from this branch (the session-start briefing warned it was built on `feature/sn-dsa`).

**One-line verdict.** The plan's *mathematical taxonomy is sound and its instincts are the
project's own*, but it was written against `neutron_transport_grand_report_v3.md` rather than
against the live tree, the ~13k-line `docs/theory/foundations/` corpus, or the 60-issue
backlog — so roughly two thirds of it re-derives architecture that already exists (often
gate-pinned), three items would regress deliberate rulings, and its genuinely new content is a
small, sharp, high-value subset. **Keep ~7 items, retarget ~10 onto existing issues, drop the
rest.**

---

## 0. The process finding (most important, and cheapest to fix)

- The plan cites **zero** GitHub issues. Verified: `grep -o '#[0-9]\{1,3\}'` → 0 matches.
- It cites **one** theory page, `docs/theory/galerkin_projection.rst`, which **has not existed
  since ~2026-07-13** (renamed → `docs/theory/foundations/frame.rst`; the `_galerkin-projection:`
  label survives *inside* it, so the intended content is right and the pointer is dead).
- Its subject matter is covered by `docs/theory/foundations/`: `operator_algebra.rst` (4708 L),
  `frame.rst` (3975), `operator_inverse_family.rst` (1422), `operator_tensor_network.rst` (1064),
  `coupled_block_operator.rst` (710), `operator_adjoint.rst` (679), `field_algebra.rst` (591).
- And by ≥10 open issues (#296, #300, #261, #256, #205, #260, #213, #289, #219, #277) and ≥5
  existing plan files (`p4_5_two_param_operator_split.md`, `operator_inverse_algebra_carve.md`,
  `coupled_block_operator_campaign.md`, `typed_carrier_grid_carve.md`, `operator_algebra_reframe.md`).

Per **Cardinal Rule 4** (issues *are* the plan/log) and **Rule 2** (single source of truth), a
1500-line parallel roadmap is a **twin source of truth for the roadmap** — the very Pattern-2
defect the plan argues against, applied to planning artefacts. Not a reason to discard it; a
reason to *reconcile* it.

---

## 1. ALREADY DONE — the plan proposes existing architecture as new

| Plan proposal | Reality | Evidence |
|---|---|---|
| "`Operator` = one declared domain + one codomain" | **IS the implemented type system**, with a PEP-696 default so endomorphic sites need no change | `numerics/operator.py:504` `class LinearOperator(Protocol[Domain, Codomain])`; TypeVars `:117-118` (`Codomain` defaults to `Domain`) |
| Mint `ScalarSourceSink`, `HarmonicMomentSourceSink`, `AngularSourceSink` | **All exist**; `transport/source_sinks/` is a 7-module package | `source_sinks/scalar_source_sink.py:89`, `harmonic_moment_source_sink.py:104`, `angular_source_sink.py` |
| `S_SN = R·Λ·M` as named typed factors | **Documented as the design-of-record AND implemented** | `operator_algebra.rst:2941-2953` types `M`/`Λ`/`S` explicitly; realized as `frame.reconstruct(Λ.apply(φ))/W`, `scattering.py:1231-1248` |
| `L` σ-free, `C = M[σ_t]` separate | **Done**, with a foundation catcher + a Mode-11 σ-leak mutation | `operator_algebra.rst:530-574`; `tests/sn/operators/test_pure_L_sigma_free.py` (216 L) |
| Four-tier posing (leaves → posing → resolvent → algorithm) | **Documented + realized**; it *is* lesson L23(b) | `operator_algebra.rst:3938-4112` (`_eigenvalue-posing:`), commits `650032e`/`7603c8e` |
| `(n,2n)` as a separately-named leaf | **Two** already exist, and the code calls it "DISTINCT" | `scattering.py:299` `N2NMomentOperator`; `isotropic_scattering.py:345` `IsotropicN2N` |
| `Vector` protocol / typed algebra | Exists (#256 step 1 landed) | `numerics/vector.py:94` |
| `WithinGroupPosing` record | `WithinGroupSystem` docstring literally reads *"The POSED within-group system: the loss `A` and its named splitting `A = M − N`, constructed together"* | `sn/coupled_system.py:283-328` |
| "`.kernel` should become a first-class fixed-domain object" | It already **is** a typed composable `LinearOperator` | `fission.py:255-309` → `TensorProductOperator` |
| Extract method-neutral reduced operators | **Already extracted and consumed**: diffusion + homogeneous bypass the SN facade entirely and use the frozen, quadrature-free leaves | `diffusion/solver.py:236-237`, `homogeneous/solver.py:146` use `IsotropicScattering`+`IsotropicN2N` |
| Layering `numerics ← transport ← methods` | **ZERO violations today**, gated | `tests/test_layer_imports.py`; explorer audit found none |

---

## 2. WOULD REGRESS a deliberate ruling — do not execute as written

**2a. Phase 7 "remove runtime carrier polymorphism."** The multi-arm `apply` is **Pattern M**,
chosen deliberately over two named alternatives with a written rationale (Beck's rule ordering:
"reveals intention" > "fewest elements"), implemented as `singledispatchmethod` + a
`TYPE_CHECKING` `@overload` surface aliased to the same runtime object (so runtime is
byte-identical), and pinned by a **630-line** gate with static `assert_type` pins,
mutation-verified teeth, and a runtime dispatch-parity check.
`operator_algebra.rst:410-527`; `tests/sn/operators/test_operators_apply_typed.py`.

Decisive: `operator_algebra.rst:513-524` **parks this exact question on #261** —
> "The deeper *spelling* question — Pattern M versus a thin `@overload` + `match` router over
> shared primitives — is deliberately **parked on #261**, to be settled *together* with the
> C/F/S core relocation. … Deciding the spelling before the cores move would be premature; the
> sharing should dictate the form."

The plan decides it now, in a *third* direction, without engaging the ruling. Its Phase-1/2
sequencing (cores first) actually **agrees** with the ruling — it just doesn't know it is
executing #261.

**2b. Retiring the fused-kernel arm** (§5.3 "raw ndarray … should not be a second public
`apply` signature"; Phase 7.5). The corpus deliberately keeps **two** realizations of `RΛM`:
the fused `frame.conjugate(Λ)` ndarray kernel (the **0-ULP canary**, pinned by
`test_scattering_kernel_crosscheck.py`) and the explicit typed edges — *both routing through
the same `Λ` and the same `R`*, so they are **not** a twin path. This is precisely the
`coding-standards.md` **fuller-view-oracle exception**. `operator_algebra.rst:3072-3111`.

**2c. Moving the harmonic frame into `sn` (§6.2).** `HarmonicFrame`'s consumers are **all inside
`transport/`** — `scattering.py` (×7), `frames/__init__.py` (×3), `reaction_rate_functional.py`
(×2), `source_sinks/harmonic_moment_source_sink.py` (×1). Moving it to `sn` strands
`transport/reaction_rate_functional.py`, and `tests/test_layer_imports.py` forbids
`transport → sn`. As written this is a **layering violation, not a cleanup**. (The documented
reason it lives in transport: the generic ndarray faces are shared 0-ULP with the P3
indicator-homogenisation frame — `operator_algebra.rst:2999-3060`.) See §4c for the *real*
frame problem, which the plan sensed but mis-located.

---

## 3. NAME COLLISIONS — right idea, wrong word

**3a. `Kernel` — the big one.** The grand report's **own §5.6 suffix law** (line 606 ff.) says:
> "`Kernel` means **it is integrated against a measure**."

and §5.7 places `IntegralKernelOperator → {ScatteringOperator, FissionOperator}`. The code
implements exactly that: `transport/operators/integral_kernel_operator.py:164`
`class IntegralKernelOperator(LinearOperator[V], Protocol[V])`, whose module docstring states
> "A Kernel **REFINES** LinearOperator (unlike the disjoint Functional) … a Kernel **is** a
> LinearOperator (it has the full operator surface) **AND adds** the `kernel` member"

with a **locality** discriminator (local/diagonal ⇒ Operator; nonlocal/integrates ⇒ Kernel)
and a 13-test intrinsic-property gate carrying **negative** tests
(`tests/transport/test_integral_kernel_category.py` — Multiplication/Identity are *not*
Kernels; Functional is disjoint), satisfying `vv-principles` anti-pattern #11.

The plan's §1 asserts the **opposite**: kernels are "non-callable; do not inherit
`LinearOperator`; do not expose `apply`" — and reuses the node name `FissionKernel`, which in
§5.6's hierarchy already denotes the *integrable* kernel. Same word, inverted meaning.

⇒ The plan's concept (method-neutral, non-callable reaction **data + invariants**) is
legitimate and needs a **different suffix**. Project precedent is `BoundaryTraceLaw`, which the
plan itself cites approvingly as the model seam. The plan's §14.1 offers `Law`, then picks
`Kernel`; **the evidence inverts that choice — pick `Law`** (`FissionEmissionLaw`,
`LegendreScatteringLaw`). Also check first whether `CoefficientField` (#257) is already the
answer for "cross sections as fields".

**3b. `condensation`.** The plan uses it for the SCC quotient graph. In ORPHEUS the word is
overwhelmingly **energy-group condensation** (120 hits in `frame.rst` alone; `Mixture.condense`,
`Solution.condense`, `derive_energy_condensation_exactness`). Say "SCC quotient DAG".

**3c. `resolvent` — three live meanings.** (i) the corpus uses it for `K_pm ≡ A_loss⁻¹M`
(`operator_algebra.rst:3960`); (ii) `WithinGroupSystem.resolvent` **holds `M`**; (iii) the plan
uses it for `A⁻¹F`. Whatever the rename, it must pick one convention corpus-wide.

---

## 4. GENUINELY VALUABLE + NOVEL — the keep list

**4a. `SNSolver.self.L` holds `L + C`. TRUE, and worse than the plan states.**
`sn/solver.py:1040` and `:1102`: `self.L = build_streaming_collision(...)`, which returns
`StreamingOperator(...) + MultiplicationOperator(...)` (`coupled_system.py:374-377`).
Three sharpenings from the grounding pass:
- **`L` has three meanings inside one `__init__`**: the local `L` at `solver.py:990` is the
  *Legendre truncation order*; `self.L` at `:1040` is the `L+C` composite; and the honest σ-free
  streaming leaf `L` is the composite's left summand.
- `solver.py:1035` writes the misnomer out ("*the full transport operator L = Ω·∇ + Σ_t*") while
  `:1100` states the honest convention one line above its own rebinding.
- **`self.L` / `self.S` / `self.F` are PRODUCTION-DEAD.** `self.S`/`self.F` are read nowhere;
  `self.L` is read by exactly two tests, *both of which rename it at the call site*
  (`test_typed_residual_evaluation.py:233` labels it `"L+C"`;
  `test_curvilinear_operator_admits_mms.py:84` binds `LC, S = solver.L, ...`).

⇒ **The right remedy is aggressive retirement, not the plan's `SNDiscreteOperators` + aliases.**
Delete the vestigial `(L, S, F)` framing surface; cost is 2 test edits. Cheaper *and* more
correct than what the plan proposes.

**4b. `WithinGroupSystem.resolvent` holds `M`. TRUE.** `sn/coupled_system.py:327` — the field is
documented as "`M` — the sweepable part". Rename to `implicit_operator` (and `gains` →
`explicit_gains`). ⚠ Blast radius: **`resolvent` appears in 74 test files** — this needs the L20
three-search audit, and it should fix the §3c three-way overload at the same time.

**4c. `InvertibleOperator` names a *capability*, not the object.** `sn/operators/streaming.py:445`
— it wraps exactly a `StreamingOperator` + a `MultiplicationOperator`, i.e. `L+C`, "whose `solve`
is the WDD sweep". Naming a class after its capability is what **#213** argues against.
`FreeTransportOperator` / `AttenuatedStreamingOperator` reads like the physics. 36 test files
mention it; not in the module's `__all__`.

*Bonus defect found here:* `InvertibleOperator.solve(rhs, *, initial_guess=None)` **accepts and
drops** `initial_guess` (`del initial_guess`, `streaming.py:831`) — but lesson **L21** added that
kwarg deliberately so the sweep and matvec could share ONE M-M seed strategy. Either the seam
regressed or the parameter is now dead; it needs a ruling, not a silent `del`.

**4d. The held pencil `(A, M)` — the plan's single best idea, and now *justified* by the
project's own defer-until-2 rule.** There is no `Pencil`/`PosedEigenproblem` type (0 hits for
`GeneralizedEigenPencil`). What exists is the **same concept spelled three ways**:
- `KEigenvalue(A, S, F)` holds the triple **un-fused** as attributes (`iteration.py:1171-1173`) —
  but its *only* production consumer is the **adjoint** path (`sn/solver.py:2472-2474`);
- `DiffusionSolver` holds `self.loss` and `self.fission` **separately, un-materialized**
  (`diffusion/solver.py:240-243`) — then eagerly LUs `A⁻¹` alone
  (`MatrixInverseOperator(FlattenedOperator(self.loss, …))`, `:250-252`);
- `HomogeneousSolver` eagerly **fuses and densifies** `K = A⁻¹F`:
  `K = MatrixInverseOperator(loss, basis_shape=(ng,1)) @ production` (`homogeneous/solver.py:194`)
  then `dominant_eigenpair(K.as_matrix(...))` (`:202`);
- `direct_eigenvalue(A, F)` / `rayleigh_quotient_iteration(A, F)` take the pair as a **calling
  convention on bare ndarrays**; `direct_eigenvalue` then eagerly fuses at `eigenvalue.py:468`.

So **three solver families pose the SAME `A⁻¹F` eigenproblem three different ways through the
SAME shared leaves** (SN: no pair at all; diffusion: held pair, iterate; homogeneous: fused +
densified). That is the twin.
- The **SN forward** eigenvalue path holds **no pair at all**: `power_iteration` consumes the
  `EigenvalueSolver` Protocol, and the posing (`A_loss = L+C−S−B`, `M = F`, `k = μ`) exists only
  as **prose in a module docstring** plus arithmetic scattered across three `SNSolver` methods
  with a `str`-dispatch on `self.inner_solver`.
- The three engines "solve the same posed problem" **only by an assertion in a docstring**
  (`eigenvalue.py:392-394`).

That is a three-way twin of a method-neutral concept, and two genuine un-fused witnesses already
exist (diffusion's pair, `KEigenvalue`'s triple) — so **Pattern 6's "unify after two instances"
is satisfied now**, not speculatively. It also un-blocks #277 (RQI already keeps the pair
un-fused), shift-invert `(A − σM)⁻¹` (which `operator_algebra.rst:4049-4052` names as a
"documented future seam"), generalized Arnoldi, and the α-pencil.

> **Hard constraint the plan does not state.** Lesson **L23(a)** established that
> `power_iteration(solver: Protocol)` binds the resolvent **late**, and that late binding is
> *strictly more general* — it admits CP BiCGSTAB, diffusion's direct FD inverse, and the
> homogeneous dense solve, none of which have an `(L,S,F)` factorization. A pencil must therefore
> be an **additional Layer-2a posing object**, never a replacement for the Protocol. Getting this
> backwards would repeat exactly the mistake L23 was written to prevent.

**4e. The leaf-assembly record — motivated by measurement, not aesthetics.**
- `build_within_group_system` is called **fresh on every inner solve** (`solver.py:1568` SI,
  `:1756` Krylov, `:2140` finalize), and each call re-runs `build_streaming_collision`, re-builds
  `SNBoundaryOperator`, and on a carrying mesh re-builds `Seeding`/`Emission`/`B_b`/`march` plus
  two `CoupledSpace`s with fresh `zeros` lambdas. Nothing is memoized. (`self.L` — the one cached
  composite — is the dead attribute of §4a.)
- `scattering_op` is threaded by an **optional kwarg at ≥5 call sites** purely to preserve leaf
  identity (`solver.py:1569,1757,2141,3172,3403`).

⇒ A record that builds the leaves **once** fixes both. This is the plan's §5.5, but the
justification should be the measured re-construction, not the naming.

**4f. The "regular splitting" claim — a Cardinal-Rule-1 issue, sharper than the plan says.**
The phrase appears at 11 sites (5 code, 6 docs), cited to **Hackbusch 2016 §11**
(`sn/coupled_system.py:78,290`; `sn/solver.py:172,724`;
`loss_representation/sweep_schedule.py:173`; + `coupled_block_operator.rst`,
`boundary_conditions.rst`, `cartesian_multid.rst`, `history.rst`). Varga's classical *regular
splitting* requires `M⁻¹ ≥ 0` **and** `N ≥ 0`. `N = S + B ≥ 0` holds (non-negative cross
sections). But `M⁻¹ ≥ 0` is exactly what **diamond differencing does not guarantee** — DD
famously admits negative fluxes. So the claim is *plausibly false for the production scheme*
(and would hold for a strictly positive scheme such as step characteristics).
⇒ Not cosmetics. Dispatch `literature-researcher` to pin Hackbusch §11's exact definition
(**check `scratch/literature/` + the OCR sidecars in `scratch/literature_ocr/` FIRST** per
`.claude/rules/delegation.md`), then: prove the cone conditions, restrict the claim, or rename to
`stationary splitting`. Note `ρ(M⁻¹N) = 0.371` is already measured, so convergence is
characterized regardless.

**4g. `foldable_part` / `residual_part` — the plan aims at DEAD code; retire instead of rename.**
Both survive verbatim (`scattering.py:978`, `:1016`) with **zero production callers** (the only
non-test mention outside the file is prose in `sn/operators/streaming.py:180`). The DSA campaign
that closed #215 did **not** rename or re-pose them — it went *around* them to
`mat_xs.foldable_sigma()` / `mat_xs.residual_sig_s()` (`sn/acceleration/dsa.py:233,236`) and
**planted a sentinel forbidding the σ_r-fold wiring** (`tests/sn/acceleration/test_dsa_rate.py:1063-1134`,
`_FOLD_ACCESSORS`), backing the in-file trap note at `scattering.py:968-977` (*"do NOT wire the
σ_r-SWEEP … ships 46–56 % silent flux errors on anisotropic problems"* — vv Mode 9).
⇒ So the plan's §7.5 rename would polish a test-only API. The honest options are **retire** (with
test migration per `coding-standards.md`) or **keep as a documented anchor**; that value call is
the user's. Either way the plan's premise ("rename before exposing an implicit-fold strategy") is
moot — the fold is *fenced off by a gate*, which is stronger than any name.

---

## 5. ALREADY TRACKED — fold into the issue, don't duplicate

| Plan section | Issue | Note |
|---|---|---|
| §7.2 `PartitionedOperatorView`, §7.3 SCC/schedule | **#296** | Far more developed: already maps sweep = block-**triangular** over space, source iteration = block-**Jacobi** over angle, multigroup = block-**Gauss-Seidel** over energy, and is explicitly framed as the *forward quest* after ψ½ proved the machinery. **Real fork — see §6.** |
| §7.6 boundary cycles | **#300** | Woodbury/Sherman-Morrison on the rank-1 `B`; the SCC criterion, the CP precedent (`cp/solver.py:389-397` is S-M written longhand), and a discriminating test with named failure modes are all already designed. |
| §5–§6.2 core extraction | **#261** | Has an explorer-verified per-operator core/adapter split **and names the lone L3 blocker** (`MaterialXSField`'s `SNMesh` type vs `tests/test_layer_imports.py`) — which the plan omits entirely. |
| flux/source role typing | **#205** (4×3 storage×role matrix) + **#256** | #256's **open fork #3** is verbatim the plan's central question: *"S/F `singledispatchmethod` apply: keep the multi-arm union or split into typed siblings."* |
| ndarray arms → `apply_values` | **#256** | #256's answer is **better**: make the container the universal vector type + single-source the serialization into `as_scipy_linop_typed(op, template)`, so operators never carry ndarray. That *removes* the arm; the plan only *renames* it. |
| `Λ = ⊕_ℓ Σ_{s,ℓ}` | **#260** | `SumOfTensorProductsOperator`. |
| capability sets → one domain/codomain | **#213** | Capabilities as a morphism class (Iso ⊂ PartialIso ⊂ General). |
| `π_bulk` / `j_bulk` named injections | **#289** | FullField generic over leaf types; retire the role-erasure `isinstance` parses. |
| MethodSpace / realizer | **#219** | ⭐ **The plan improves the backlog here.** #219's F-C proposes a `MethodSpaceBuilder` **registry**; the realizer registry was **DISSOLVED** in #290 P7b (`44d583e`) in favour of `TransportMethod`. The plan is **right** to forbid a registry — #219 F-C is stale and should be updated to say so. |

Also stale in the plan: **§7.1 / §14** claim `SNMethodSpace` is a plausible reaction-realizer
input. It exists, but it is a **per-face boundary-realizer payload** carrying no σ, no scheme,
and no function space — not a method space in the §7-of-the-report sense.

---

## 6. The forks to DECIDE (not to build in parallel)

- **F1 — block structure.** Plan: a *new lazy* `PartitionedOperatorView` deriving `A_ij = R_i A J_j`
  from a monolithic `A`, leaving `CoupledOperator` untouched. #296: *reify the existing SN axes
  onto the existing explicit `CoupledOperator` grid*, which already does block matvec, block
  back-substitution, `inverse()`, and `assemble()` at offsets. Both are defensible; building both
  is the one outcome to avoid.
- **F2 — the `apply` spelling.** Parked on #261 **by ruling**, with the sequencing stated: move
  the cores, let the sharing dictate the form. Recommend honouring that order.
- **F3 — the data-vs-operator suffix.** `Law` (per §3a) vs checking whether `CoefficientField`
  (#257) already covers it.

---

## 7. Defects found during grounding that the plan does NOT contain

1. **`.H` silently degrades to a bare Euclidean transpose — via TWO independent routes.**
   The *mechanism* is metric-correct: `_AdjointOperator.apply` is `A* = G_V⁺ ∘ Aᵀ ∘ G_W` with the
   metric delegated to `FunctionSpace.apply_metric` / `apply_inverse_metric`
   (`numerics/operator.py:1204-1229`; `space.py:248,268` — Moore–Penrose, so zero-weight trace
   slots zero out instead of dividing by zero). But it degrades to the identity when **(a) the
   inner spaces are `None`** (`operator.py:1221-1226`) — and both reaction facades hardcode
   `is_adjointable = True` while `homogeneous/solver.py:193` constructs `FissionOperator`
   space-anonymously — **or (b) the space carries `inner_product_weights=None`**
   (`space.py:262-264`). So "is this `.H` genuinely weighted?" is a per-instance fact, not a
   hierarchy guarantee, and it is documented on neither facade. Latent correctness trap; it
   validates the plan's §11.3 insistence on *metric-aware* adjoint tests.
2. **Fission and scattering are NOT symmetric — the plan treats them as twins.**
   `FissionOperator` holds **zero** quadrature state and is genuinely cross-method (SN +
   diffusion + homogeneous all construct it). `ScatteringOperator` takes `quadrature` as a
   **required positional field with no default** (`scattering.py:386`, no `None` branch anywhere
   in 1326 lines) and has **zero** non-SN production consumers. A single symmetric carve is the
   wrong shape.
3. **The angular frame leaks out of the scattering facade.** `sn/solver.py:552` does
   `BulkAnalysisOperator(scattering_op.frame, sn_mesh) @ sweep` — the 2-D angular-windowing path
   reaches into `S` for its frame, making `ScatteringOperator` the *de facto owner and
   distributor* of the SN angular frame. Note `frame.rst:3283` already rules "an operator owns
   its frame IFF the frame is its eigenbasis" (Funk–Hecke), and the SH basis **is** scattering's
   eigenbasis — so *scattering owning it is correct*; the smell is a **third party sourcing it
   from `S`**. The fix is a frame home on the method/angular space, not the plan's §6.2 move.
4. **`apply` is `singledispatchmethod` but `apply_transpose` is a hand-rolled `isinstance` chain
   — on both facades**, and `isotropic_scattering.py` uses a third spelling (top-of-method
   `isinstance`). Three mechanisms for one job. This is a *sharper* version of the plan's §2
   complaint: the risk is not the multi-arm forward (ruled Pattern M) but that the forward and
   adjoint surfaces can drift because they dispatch differently.
5. **`gains: tuple[LinearOperator, ...]` has presence-dependent arity and positional meaning** —
   one `CoupledOperator` grid on a carrying mesh, the `(S, B_a)` tuple seedless, with "`B_a`
   LAST — the boundary-gain convention the G-S schedule arm parses"
   (`sn/coupled_system.py:331-339`). A tuple whose arity and positional semantics depend on mesh
   presence-structure, parsed positionally downstream, is illegal-states-representable — sitting
   right next to the rename of §4b.
6. **`_Y_cached` is mutable state on a mutable dataclass** (`scattering.py:398-400`) while all
   four sibling leaves are `frozen=True`; `ScatteringOperator` and `FissionOperator` are the only
   non-frozen members of the family (L15 immutability-strata smell).
7. **`assemble()` exists on the leaves but not the facades.** `IsotropicScattering`,
   `IsotropicN2N`, `MultiplicationOperator` are assemblable; `ScatteringOperator` and
   `FissionOperator` are not — so the DSA stencil-assembly `Op → Mat` functor cannot reach S or F.
8. **The `EigenvalueSolver.converged` body is spelled FIVE times** (`KEigenvalue`,
   `DiffusionSolver`, `SNSolver`, `CPSolver`, MOC) for a concept with **zero** method content.
   `KEigenvalue`'s version is strictly best (carrier-generic; handles `norm == 0`, which the four
   copies do not — they compute `Δ/1e-30`), and its docstring says the method solvers
   "never routed through this class" *by design* (`iteration.py:1121-1125`). Differing tolerance
   *defaults* are legitimate; the re-spelled *body* is not. An L23-shaped finding: the general
   layer exists and the specific ones don't delegate.
9. **Doc drifts (#302 class).** `ScatteringOperator.kernel_summands` is `:attr:`-referenced 3× in
   `docs/api/numerics.rst` (`:225,:252,:273`) and **does not exist**. `orpheus/sn/angular_flux.py`
   **does not exist** yet 6 production xrefs to `orpheus.sn.angular_flux` dangle. `fission.py:8-17`
   still writes the **retired** C-folded algebra `(L − S − F)ψ = q`, violating the ratified
   page-wide `A = L + C − S − B` internal-consistency directive. All silent under `-W` (Python-domain
   roles need `-n`).
10. **`str`-flag posing family on one entry point**: `inner_solver`, `inner_schedule` (validated
    twice — `solver.py:747` and `:937`), `acceleration ∈ {None,"dsa"}`, `boundary_condition`,
    `eigenvalue_method`. The DSA merge **added** a third string flag where a type was the elegant
    move — a data point *for* the plan's thesis and *against* its premise that the posing
    vocabulary is already in place.

---

## 8. Recommended path

**Tier 1 — naming honesty / retirement. No math changes; independent of every fork.**
1. Delete the vestigial `SNSolver.(L, S, F)` triple (§4a). 2 test edits.
2. `InvertibleOperator` → a name for the object (§4c); rule on the dropped `initial_guess`.
3. `WithinGroupSystem.resolvent → implicit_operator`, `gains → explicit_gains` (§4b), fixing the
   three-way `resolvent` overload; L20 three-search audit first (74 test files).
4. Adjudicate `foldable_part`/`residual_part`: retire-with-test-migration, or keep as a
   documented anchor (§4g).
5. Fix the §7.9 doc drifts; consider folding into #302.

**Tier 2 — the one genuinely new piece of machinery.**
6. Reify the pencil `(A_loss, M)` as an *additional* Layer-2a posing object (§4d), honouring the
   L23(a) late-binding constraint. Companion: `operator_inverse_algebra_carve.md` Phases 2–5
   (`.inverse()` returns an operator).
7. The leaf-assembly record, justified by the measured per-solve reconstruction (§4e).

**Tier 3 — correctness claim audit.**
8. Hackbusch §11 / regular-splitting positivity (§4f) — `literature-researcher`, local library first.
9. Document (and decide) the `.H`-without-space Euclidean degradation (§7.1).

**Tier 4 — decide, don't build.** F1/F2/F3 (§6).

**Do NOT:** run Phase 7 now (parked by ruling; would need to preserve the 0-ULP oracle); move
`HarmonicFrame` to `sn` (layering violation); mint `ScalarSourceSink` et al. (they exist); mint
`ScalarSourceSpace` — the source role differs from the flux space **by name only**, with no extra
structure and no non-identity morphism, so it **fails the project's own type-vs-property test**
(`coding-standards.md`); the mechanism already in use is the `_SPACE_NAME` string.

---

## 9. What to do with the plan file itself

Rewrite it as a **short reconciliation memo** that (a) keeps §§4a–4g as its scope, (b) points
each remaining item at its issue (§5), (c) records the three forks (§6), and (d) states the
three rulings it must not regress (§2). Then delete the re-derived sections — they are a second
spelling of `docs/theory/foundations/` and of the backlog, and the plan's own Pattern-2 argument
applies to them.
