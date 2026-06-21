# CoefficientField + operator-as-promotion (#257) — implementation plan

> **POST-COMPACTION REGROUNDING (read in this order before any code):**
> 1. THIS plan (the roadmap + the embedded census facts below).
> 2. GH **issue #257** (`gh issue view 257`) — the durable thesis.
> 3. Frame memo: `.claude/agent-memory/cross-domain-attacker/coefficient_field_promotion_frames.md`
>    (the 4 structural frames + the multiplier-algebra law-suite).
> 4. The campaign anchor `.claude/plans/field_typed_operator_algebra_campaign.md` (the #256→#257 re-scope + the baseline reds).
> 5. **Carrier design + layering verdict** (behind S4.5):
>    `.claude/agent-memory/cross-domain-attacker/issue_257_carrier_typing_layering_frames.md`
>    (the `Vector`-Protocol-is-irreducible forgetful-functor finding + the cofree-comonad
>    `TimedFullField = Cofree(FullField, d)` split). The DESIGN REVISION block below summarizes it.
> 6. Then the code, in the order the stages touch it.
>
> **Branch** `feature/field-typed-operator-algebra` (off `main`, NOT pushed; local `main` @ `05fa1ef`).
> The plan is committed (`99f108f`); implementation is UNDERWAY through **S8c (`35f5612`) — S8 COMPLETE** —
> see the ⭐ CURRENT STATUS block immediately below (S9 DONE; campaign behavioral work S1–S9 COMPLETE).

## ⭐ CURRENT STATUS (2026-06-21 — S9 COMPLETE; CAMPAIGN BEHAVIORAL WORK DONE (S1–S9); only deferred S10/S11 tails remain)

> **✅ S9 DONE + committed** (`7180f72` test(sn) + `bc1ebec` docs(sn) + `aee60d4` chore matrix), branch
> `feature/field-typed-operator-algebra`, **NOT pushed**. The original premise was CORRECTED TWICE this
> session: (a) there is **no `(moment_buf, None)` boundary drop** — that sentinel is a BULK scalar; the
> moment path already returns a complete `TimedFullField(bulk=HarmonicMomentField, boundary=BoundaryFlux)`;
> (b) the user named **fork #1 = TRANSVERSE-SPATIAL** (the LD within-cell slope), motivation = make the
> **MMS boundary source** moment-valued so the LD slope is properly verifiable.
> **Verification verdict (double-confirmed: test-architect + numerics-investigator, across optical depth +
> transverse frequency + amplitude):** the boundary transverse-slope moment is **SUB-FLOOR** for the
> converged flux in every regime → **NO value/order gate keyed on the slope** (vv **Mode-10**,
> companion-isolating gate unavailable — the **4th** such instance: #240 D5b-S4 → #247 → #251 → S9). The
> coherent promise ("LD 2nd-order EVERYWHERE incl. boundary") is **already TRUE** — delivered by the
> AVERAGE moment; S9 LOCKS it, does not fix it.
> **Shipped:** production `SN2DCartesianLDStressMMSCase.prescribed_inflow` now emits the moment slot for
> LD (DD/Step BYTE-IDENTICAL — gated on `face_moment_count`) via an L11 `leggauss` `_project_inflow_to_face_moments`
> (closes the #251 producer-blindness); the coherent-promise gate `test_first_cell_row_already_second_order`
> + the sub-floor verdict pins + the Mode-11 sentinel (promoted to `tests/sn/verification/mms/test_ld_2d_boundary_promise.py`);
> GATE-B/C re-targeted onto the production producer; `_solve_with_boundary_slope` re-baselined (flat leg
> decoupled from the now-honest producer). Sphinx narrative landed (`discrete_ordinates.rst` coherent-promise
> subsection + `operator_algebra.rst` #263 criterion). Gates: 35 + 590 + GATE-D 520 baseline; pyright 2282
> = baseline, 0 net-new. elegance PASS-WITH-NITS (applied), qa SUPPORTED. **NO `BoundaryMomentField`/new type**
> (the boundary moment is a PROPERTY), **NO new ERR** (closed a producer-blindness — the slope was UNVERIFIED,
> not wrong).
> **⭐ #263 FILED** — the type-vs-property criterion: a moment earns a first-class TYPE only when a
> NON-CANONICAL dual coexists (the integrating QUADRATURE is what makes the angular ordinate↔harmonic pair
> a type → `HarmonicMomentField`); spatial order is modal-only → a PROPERTY (`AngularFlux.spatial_moments`);
> the first-class `SpatialMomentField` is DEFERRED to the **collocation trigger** (nodal-DG / Lagrange-FEM —
> NOT NEM/ANM, which are modal). **#256 fork #1 RESOLVED + the `BoundaryMomentField`-as-subclass decision
> REVISED** (comment posted). Carry-forwards unchanged: **#261** (C/F/S relocation + the dispatch-spelling
> decision), **#262** (bulk-accessor under-typing).
> **NEXT:** only **S10** (χ-simplex at `Mixture.chi`, DEFERRED end-of-plan) + **S11** (`Field`→L2, DEFERRED /
> UNDER REVIEW) remain — both deferred tails; the campaign's behavioral work is COMPLETE. The owed
> CONSOLIDATED archivist pass ✅ **DONE** (`64d98fc` docs + `c9f5af5` matrix, 2026-06-21): the 3 `.. todo::`
> stubs (S3b `multiplication-operator-promotion`, S5 `functional-category`, S6 `integral-kernel-category`)
> expanded to full §5.5–5.7 narratives + the S8c **Pattern-M idiom** (`heteromorphic-apply-typing`) + a brief
> S7 scipy note; 2 new definitional eq-labels; build clean (paramref baseline only). **Campaign doc debt CLOSED.**
> ⭐ **COMPACTION LANDING POINT:** behavioral work S1–S9 DONE + doc debt closed, all committed on
> `feature/field-typed-operator-algebra` (HEAD `c9f5af5`), **NOT pushed** (whole campaign branch unpushed).
> **NEXT = DISCUSS S10 + S11 with the user** (per user, deferred to post-compaction) — do NOT implement them
> unprompted; both are end-of-plan, decide-with-hindsight (S10 = χ-simplex at `Mixture.chi`, needs an L20
> upstream-fixture audit; S11 = `Field`→L2, the decision-rule + for/against are in the S11 entry below).

> ✅ **DOC NUANCE LANDED (2026-06-20, `8d98714` + matrix `23b8eb0`).** The user-requested
> **apply-linear / solve-nonlinear asymmetry** subsection (`:label: apply-solve-asymmetry`,
> `operator_algebra.rst`) is committed: `(L+C)ψ = Lψ+Cψ` but `(L+C)⁻¹ ≠ L⁻¹+C⁻¹`; "apply & solve are
> two views of the same operator" holds for the bundled `InvertibleOperator` (= L+C, advertises
> `CAP_SOLVE` = the sweep), NOT the leaves (`StreamingOperator` = `{CAP_APPLY, CAP_APPLY_TRANSPOSE}`,
> verified); the WDD cell-denominator proof `1/(k_stream+σ)`; the operator-splitting **Neumann series**
> `(L+C)⁻¹ = C⁻¹(I+LC⁻¹)⁻¹ = …` + the parallel/resistor identity + the transport-native source-iteration
> series; + the EXPANDED S8b `streaming-action-pure-l` stub. archivist-verified against the code;
> incremental `sphinx -W` clean. **Working tree FULLY CLEAN (S1–S8c + the doc, all committed). RESUME
> DIRECTLY AT S9.** (NOTE: the deferred S3b/S5/S6 + S8b/S8c stubs remain owed for the consolidated
> archivist pass — now incl. the **"Pattern M" idiom** [the first `@overload` in `orpheus/`, the
> elegance do-now nit]; that is separate, not S9-blocking.)

**DONE + committed** on `feature/field-typed-operator-algebra` (NOT pushed):
- Foundation: `cfb651b` (`Vector` Protocol) + `41a92cb` (`apply(x:V)->V`); `e3f90d5`/`99f108f` (re-scope + this plan).
- **S1** `505e1b7` — `CrossSectionField` (1/cm cone) + `CoefficientRole` + `CROSS_SECTION_UNITS`. ⚠ `SpectrumField` DROPPED (see S1 entry); χ-simplex → S10.
- **S2** `1ce727a` — `MaterialXSField.{total_cross_section,absorption_cross_section,fission_production}_field` typed accessors (zero-copy: `.values IS` the raw cached array).
- **S3a** `c1da42d` — `DiagonalOperator` generalized → N-D broadcast ENGINE: `DiagonalOperator(coefficient, broadcast_axes=None, *, axis=0)`; 1-D mode byte-identical; the σ_t case `broadcast_axes=(0,)` on a `(N,ng,*spatial)` carrier ≡ `sigma[None]*carrier`. `weights` now a property (raises in broadcast mode). `from_measure`/self-adjoint/`CAP_SOLVE`-iff-nonzero preserved.
- **S3b** `4ebade6` — the §5.7 promotion `C = M[σ_t]`. `orpheus/transport/multiplication_operator.py` (FIRST operator in `transport/`): `MultiplicationOperator(CrossSectionField)`, the multiplier-algebra embedding `M:L^∞→B(L²)`; delegates the multiply to a STORED S3a engine (built once over the immutable coefficient — elegance store-once nit), typed codomain (apply→AngularSourceSink, solve→AngularFlux, `M.H=M`). `CollisionOperator` → thin SUBCLASS (keeps name + `(sn_mesh, sigma:ndarray|CrossSectionField)` back-compat ctor + `.sigma`/`.sn_mesh` + `L+C→InvertibleOperator` `+`-dispatch; duplicated apply/solve bodies DELETED — inherited). `solver.py:217/925/1000` source σ from `total_cross_section_field`. ⭐ **Behavioral strengthening**: `CAP_SOLVE iff min|σ|>0` (Pattern 4, inherited from engine; legacy was always-on→silent NaN; AUDITED safe — no prod `.solve` on C, σ_r path is #200/#215 not live). Gate = principled-equivalence: broadcast oracle 0 ULP (2-D nx≠ny ng=2), law-suite ≥2G-asym-het, k_∞ + streaming-equilibrium refs green (Mode-11 sentinel-confirmed route σ through promoted C). elegance PASS (2 nits fixed: dead CAP imports + store-once engine), qa SUPPORTED. Sphinx §5.7 STUB added (`:label:multiplication-operator-promotion`, `.. todo::` for the archivist).
- **S4** `93aa016` — `TransportState(Vector, Protocol)` (`orpheus/transport/state.py`), the §5.5 honest generic dissolving the opaque `Generic[V]`: refines numerics `Vector` with read-only `@property` members `bulk: BulkField`/`boundary: BoundaryField`/`history_depth: int` (the #208 composite accessors; read-only because a writable Protocol attr is invariant → rejects the frozen `TimedFullField`). `runtime_checkable`; discriminating test (`test_transport_state.py`, foundation): `TimedFullField` IS, `np.ndarray`+bare `AngularFlux` are `Vector` but NOT `TransportState` (strictly stronger). **PURE ADDITION** — numerics stays `Generic[V: Vector]`; NO operator re-pointed (see ⚠ OPERATOR RE-POINT below). pyright 2295=baseline NO offset; qa SUPPORTED (mutation-verified teeth). ⚠ **SUPERSEDED by S4.5** (the `TransportState` Protocol was retired).
- **S4.5** `b404ae1` (**closes #217**, SUPERSEDES S4) — extracted the timeless concrete `FullField` base (`orpheus/transport/full_field.py`) out of `TimedFullField` (cofree split `TimedFullField=Cofree(FullField,d)`): the 6 vector dunders + `from_flat`/`copy` live ONCE on `FullField`, routed through a polymorphic `_recombine(self:T,*,bulk,boundary)->T` hook; `TimedFullField(FullField)` overrides `_recombine` (empty history + preserved depth) + adds `_history`/`advance`/`at_lag`; `zeros` delegates to `FullField.zeros`; `from_flat` generic-over-template routed through `_recombine`; `_check_partner` widened to `isinstance(other,FullField)` (the `at_lag(0)-at_lag(1)` timed−timeless stencil). **S4's `TransportState` Protocol RETIRED** (`state.py` deleted; discriminating test now NOMINAL in `test_full_field.py`). Operators stay `["TimedFullField"]` (IS-A FullField); the `apply(x:FullField)` re-point folds into S8. BIT-IDENTICAL (dedicated TimedFullField suite byte-unchanged + green; `_recombine` two-way mutation teeth); pyright 2295=baseline 0 net-new + **0 net-new `type: ignore`** (trial `from_flat` override-ignore removed structurally; `zeros` ignores de-duped to original 2). elegance PASS (3 doc nits applied), qa SUPPORTED.

- **S5** `e9e0121` (+ matrix chore `7b463bb`) — the §5.6 `Functional` category. `orpheus/numerics/functional.py` (L1): `Functional(Protocol[V_contra, R_co])`, `@runtime_checkable`, one method `evaluate(x)->R`; the co-vector companion of `Vector`, a SIBLING of `LinearOperator` (NOT a subclass — carries NONE of `apply`/`capabilities`; the disjoint surface IS the category's defining property). Contravariant input / covariant UNBOUNDED result (a `float | V` union would mistype the per-cell scalar-field). `orpheus/transport/production_rate_functional.py` (L2 — shared SN/CP/MoC, NOT sn/): `ProductionRateFunctional(nu_sigma_f: CrossSectionField)`, `evaluate(phi)->(1,*spatial)` density `p(r)=Σ_g νΣf_g φ_g` (group-collapsed, keepdims=True, NO volume measure) — byte-identical to the anonymous `inner` in `RankOneOperator.apply`, so S6's `F=M_χ∘ProductionRate∘M_νΣf` inherits bit-identity. **Additive + bit-identical**: `RankOneOperator`/`FissionOperator` NOT rewired (S6). Estimators stay bare `(L,S,F,ψ)` callables (they consume the operator TRIPLE, not a lone field → not `Functional[V,R]`; the category just NAMES their field→scalar core). Gate: intrinsic-property (Functional≠LinearOperator both directions + discriminator foils, Mode-8 clean, teeth mutation-verified — axis-flip reds 6, ×1.5 measure-fold reds 5); `evaluate` vs a STRUCTURALLY-INDEPENDENT hand-derived double-loop (L11) + pinned 0-ULP LITERAL-RANK byte-identity vs the legacy `inner` (B.2 tightened post-qa — was squeeze-agnostic). pyright 0 net-new + 0 net-new `type: ignore`. elegance PASS-WITH-NITS (NIT-1 rationale-only: the `isinstance(phi,Field)` idiom's real precedent is `scattering.py:618/665`, NOT `MultiplicationOperator`; NIT-3 declined), qa SUPPORTED. Sphinx §5.6 STUB + `:label:production-rate-functional`/`functional-category` (full narrative DEFERRED to post-S6 archivist). ⚠ **Filed #259**: the production-rate `Σ_g νΣf_g φ_g` is coded 3× (SN `solver.py:1086` vol-weighted group-last+n2n / CP `solver.py:674` / numerics `iteration.py:269`) + the `KEigenvalue` estimator-injection seam is DEAD in production (SN/CP bypass via `power_iteration(self)`) — unify after the campaign reaches S6.

- **S6** `f509e74` (+ matrix chore `9ccfec7`) — the §5.6 **`IntegralKernelOperator`** category (completes Operator/Kernel/Functional). `orpheus/transport/integral_kernel_operator.py` (L2): `@runtime_checkable IntegralKernelOperator(LinearOperator[V], Protocol[V])` requiring a `kernel` property; a REFINEMENT of `LinearOperator` (still apply/caps + adds kernel — strict: a local `MultiplicationOperator`/`IdentityOperator` is NOT a Kernel), contrasted with the disjoint S5 `Functional`. **User scope (Q1/Q2 2026-06-20): keep the bit-identical matvecs (composition=semantic reading + cross-check); create base + reframe BOTH fission & scattering IN PLACE in sn/.** `FissionOperator` (additive): adds `production_rate`→S5 `ProductionRateFunctional` (the §5.6 middle factor of `F=M_χ∘ProductionRate∘M_νΣf`); `kernel`=fused RankOne χ⊗νΣf kept (Pattern 5). `ScatteringOperator` (additive): adds `kernel`=`OperatorProduct(R,Λ,M)` reproducing `_aniso_source_from_moment_values(M·ψ)` byte-for-byte; the kernel is the nonlocal ℓ≥1 aniso redistribution (P0 iso/n2n/1-W are the local components — a STRICT sub-component, pinned). 5 dispatch arms UNTOUCHED. ⚠ findings: `M_χ` is rank-CHANGING (not S3b's `MultiplicationOperator`) → literal composition rejected; scattering's R/M are rank-changing einsums (forbidden as TP factors) → `SumOfTensorProductsOperator` kernel rejected as a re-derivation. **3 `cast(LinearOperator,...)` bridge the unparametrised-`LinearOperatorMixin` generic gap (NOT `type:ignore`; #226 scope).** Gate: intrinsic-property (both ARE Kernels, both-directions + discriminator + direct kernel-attr-absence for the runtime_checkable loophole — qa-confirmed); two 0-ULP cross-checks (fission vs hand-loop+RankOne; scattering vs existing aniso path); the aniso-only-sub-component gate (qa-flagged misread guard); matvec byte-identical (pre-T3 snapshots + aniso MMS). pyright **2307=S5 baseline, 0 net-new** (helper-typing cleaned the 4 skeleton errors too), 0 net-new `type:ignore`. elegance PASS-WITH-NITS (N1+N2 docstring fixes applied), qa SUPPORTED. Sphinx §5.6 Kernel STUB + `:label:integral-kernel-category`. ⚠ **Filed #260** (un-orphan `SumOfTensorProductsOperator`) + **#261** (relocate shared C/F/S cores to transport/ + CP/MoC carrier unification — the DEFERRED relocation). ⚠ qa crossed L28 (git stash to measure baseline) but restored cleanly — tree verified intact.

- **S7** `1991d46` (+ matrix chore `f00d482`) — **scipy single-source** (bit-identical plumbing). The live ORPHEUS↔scipy Krylov boundary was TWO inline `spla.LinearOperator` closures in `KrylovAcceleration.solve` (`A_scipy` system matvec + `M_scipy` preconditioner), each duplicating the ravel-wrap (`_ravel`/`_unravel_like`); plus the public zero-prod-caller `as_scipy_linop` (`operator.py`, flat-only twin). THREE construction sites for one concept. Consolidated to ONE module-private `_as_scipy_linop(carrier_matvec, template, n)` in `iteration.py` (forced home: `iteration.py→operator.py` is one-way, so the ravel-aware adapter can't live in `operator.py`; the flat L0 case is the degenerate case of the ravel-aware one → one subsumes both). System matvec extracted as a named `loss_minus_gains(psi)->V` (reads like the within-group operator, not buried under plumbing). **RETIRED** `as_scipy_linop` (def + docstring ¶ + orphaned `spla` import + `__all__` ×2 across operator.py+`numerics/__init__`) + its 5 bare-ndarray tests (the one unique guard — `MissingCapability` on missing `CAP_APPLY` — is covered strictly-stronger at composition time by `KrylovAcceleration.__init__`); 3 doc xrefs repointed to the internal adapter. ⚠ explorer audit FIRST (public-API-retirement proactive trigger) confirmed zero prod callers (Nexus `callers`+`impact`+grep) + the layering verdict + the gate-test set. Gate = bit-identical: matvec computation character-identical (same `L.apply(psi)` then `for g: out-g.apply(psi)`, same order); A↔`solution_template` (flux space), M↔`q_ext` (source space) — qa **sentinel-confirmed BOTH branches fire** on the typed `TimedFullField` carrier (A built2/fired160, M built2/fired161), NO template swap (the Mode-2 hazard). pyright **2307→2297** (−10 from deleted pre-existing-noise lines, **0 net-new**; **−1 `# type: ignore`**); Krylov+round-trip gates 139 green; broad regression at the exact 7 baseline reds; Sphinx -W clean. elegance **PASS** ("cleanest carve in the S-series"; 2 optional doc/comment nits applied), qa **SUPPORTED**. ⚠ **archivist rich-narrative DEFERRED** into the consolidated S3b/S5/S6 pass (method-implementer's mechanical doc edits already keep -W green + factually correct; S7 is plumbing, not taxonomy — no point fragmenting the consolidated pass).

**Health**: baseline reds = **7** (pre-existing #250 SPHERE ×5 + #232 mu_y ×2 — NOT ours; every stage must keep exactly these). ⚠ **pyright baseline corrected to 2307 errors / 19 warnings** (the host-tree `b404ae1` count; the "2295" recorded through S4.5 was a stale S4.5-WORKTREE figure — re-measured this session, S5 files individually clean ⇒ no masked offset). Sphinx clean. Filed **#258** (units.py pint `Unit`/`PlainUnit` stub debt) + **#259** (production-rate fragmentation). Regression subset:
`.venv/bin/python -O -m pytest -q tests/sn/operators tests/sn/spatial tests/sn/sweep/core tests/sn/solve tests/numerics tests/transport --deselect tests/sn/solve/test_keff_slab.py::test_heterogeneous_absolute_keff`

## ⭐⭐ DESIGN REVISION (2026-06-19 — `TransportState` → `FullField`; explorer + cross-domain-attacker)

User flagged `TransportState` as mis-named ("not everything is a state") + asked to fold in **#217** (a
timeless `FullField` base, no history) + asked whether the `numerics ↛ transport` barrier is a
folder-hierarchy artifact. Dispatched **explorer** (layering/import map) + **cross-domain-attacker**
(structural frames). Both CONVERGED:

- **The layering is sound; the `Vector` Protocol is IRREDUCIBLE — not a hierarchy artifact.** `Vector`
  is the object-image of the forgetful functor `U: C_carrier → C_vec`; `numerics` lives in `cod(U)` so
  it structurally cannot name the carrier. The barrier is the *shadow of the layer DAG*; the only lever
  that removes it is moving the whole algebra+drivers to L2 (abandon method-agnosticism — rejected).
  The user's "bad-hierarchy" hypothesis is **disconfirmed**.
- **#217 split is structurally FORCED (cofree comonad).** `TimedFullField = Cofree(FullField, depth=d)`;
  operators are base arrows `FullField → FullField`; only the iteration drivers see the comonad. Typing
  an operator OUTPUT (`Cψ`=source, `b−Ax`=residual) as `TimedFullField` hands it an unused history tail
  ("a type error of altitude"). So the timeless **`FullField`** (bulk⊕boundary + Vector algebra, NO
  `history_depth`) is the correct operator-algebra carrier — confirming the user's "drop history_depth".
- **`Field`-in-numerics is genuinely a smell** (module-over-ring: `Field` is overloaded as flux-module
  base AND coefficient-ring set; zero numerics consumers; `numerics/__init__` re-export tell) — BUT
  relocating it is ORTHOGONAL to `Vector` (won't dissolve the barrier). → **deferred to end-of-plan**
  (user: "why should Field go to L2? push to the end with detail, review at the end"). See S11 below.
- **Fibration finding (act on in S6/S8):** `ScatteringOperator` (`scattering.py:286/618/665`) +
  `FissionOperator` (`@singledispatchmethod` over TFF/ScalarFlux/ndarray) branch on carrier type → the
  single-`V` typing is a partial lie for those two; resolve when they're re-typed.

**LOCKED (user, 2026-06-19):** `FullField` = a CONCRETE base class (close #217), `TimedFullField(FullField)`
inherits. Retire the S4 `TransportState` Protocol. Keep the `TimedFullField` name (strong grep-signal;
document as `Cofree(FullField, d)`).

**S4.5 ✅ DONE (`b404ae1`, closes #217, SUPERSEDES S4) — the record of what landed (was the realized
design):** Behavioral-NEUTRAL for `TimedFullField` (internal extraction):
- `orpheus/transport/full_field.py` (NEW): `@dataclass(frozen=True, kw_only=True) FullField` — the
  timeless carrier. Lift OUT of `TimedFullField` (DRY, the algebra lives ONCE): `bulk: BulkField`,
  `boundary: BoundaryField`, the vector-space dunders (`__add__/__sub__/__neg__/__mul__/__rmul__/
  __truediv__`), `to_flat`/`from_flat`, the `_check_partner` mesh/class guards, `zeros`, `copy`, the
  shared `__post_init__` validation. NO history.
- `TimedFullField(FullField)`: adds `_history`/`history_depth`/`advance`/`at_lag`/`history_length` + the
  history-aware `__post_init__` extension. ⚠ the lifted dunders must still return `TimedFullField`
  (empty history, preserved `history_depth`) for a `TimedFullField` operand — use a polymorphic
  recombine hook (base builds via `type(self)`-aware constructor, or TimedFullField overrides the 6
  dunders). #217: "algebra results carry empty history."
- Retire `transport/state.py`; `transport/__init__` exports `FullField` + `TimedFullField` (drop
  `TransportState`). Migrate `test_transport_state.py` → `test_full_field.py`: nominal discriminating
  checks (np.ndarray NOT a `FullField`; bare `AngularFlux` NOT a `FullField`; `TimedFullField` IS a
  `FullField`; `FullField` IS a `Vector`) + the timeless-vs-timed distinction.
- Operators STAY `LinearOperatorMixin["TimedFullField"]` for now (TimedFullField IS-A FullField →
  consistent). The operator re-point + codomain-timeless folds into S8 (below).
- Gate: BIT-IDENTICAL `TimedFullField` behavior + algebra results (the extraction is internal); the
  full TimedFullField suite stays green; pyright 0 net-new; the FullField discriminating test. Standard
  cycle (method-implementer → elegance + qa → L12 → commit).

**⭐ S8 IN PROGRESS — 3 sub-stages** (test-architect gate spec = `.claude/agent-memory/test-architect/issue_257_s8_streaming_pure_L_verification.md`, AUTHORITATIVE for S8b/S8c). Sequence: **S8a re-point → S8b drop (L+C)−C → S8c fibration.** Equivalence verdict for the carve = **PRINCIPLED-equivalent, NOT bit-identical** (the recomposition re-associates the FP tree: probe-grounded bulk drift slab 0 / sphere ≤4 / cyl ≤2 ULP; boundary trace strict 0-ULP). ⚠ Mode-11 landmine: `StreamingOperator.apply` has ZERO graph callers — the matvec leaf is reached ONLY by direct-`L.apply` gates (C1/C2/snapshots); the sweep `(L+C).solve` routes around it, so a solve-only gate is VACUOUS for the matvec change.

- **S8a ✅ DONE** (`f7d7135` + matrix `9316321`) — `apply(FullField)` re-point + helper widening. VALUE-NEUTRAL/bit-identical: operator matvec leaves (Streaming/Invertible/Multiplication/Boundary + the S/F TimedFullField-input arms) now consume+emit the timeless `FullField` (base arrows); the driver reattaches the timed type via `TimedFullField.__add__`/`_recombine` (NO `advance` — steady-state never advances, zero prod `advance` callers; the SI fold is timed-first so `_recombine` wins). `_require_typed_composite` widened `TimedFullField→FullField` (IS-A). ⭐ `InvertibleOperator.solve` STAYS `TimedFullField` (the iterate-building resolvent, covariant/Liskov-safe); `MultiplicationOperator.solve` flips timeless (zero callers, base-arrow inverse) — the discriminator is **iterate-builder vs base-arrow**, not solve-vs-apply. `full_field_space._rebuild`→`_recombine` (the one file beyond the leaf set; the metric analogue, resolved 7 g_adjoint reds). `(L+C)−C` arithmetic + the S/F `@apply.register` INPUT dispatch UNTOUCHED (S8b/S8c). Gate C5 (`tests/sn/operators/test_apply_full_field_codomain.py`): output history-free FullField across geom×depth (C5a, calls `.apply` DIRECTLY — Mode-11; mutation-verified teeth: reverting one leaf reds 4 C5a), k_inf=νΣf/Σa driver-reattach (C5b), advance type-guard (C5c). pyright 2297=baseline 0 net-new + 0 net-new ignore. elegance **PASS** (no nits — "textbook-clean"), qa **SUPPORTED** (mutation-grade). ~41 operator-suite tests re-pointed (type-surface only). ⚠ cosmetic nit deferred: a few re-pointed test FUNCTION names still read `..._timed_full_field` (arguably still correct — they test the timed INPUT arm).

- **S8b ✅ DONE** (`08cb9cf` + matrix `3e4fccb`) — dropped the `(L+C)−C` fold → pure σ-free `StreamingOperator`. ⭐ STEP-0 verified the WDD matvec is AFFINE in σ (two probes: `loss_action(0)==loss_action(σ)−σ·ψ` AND `streaming(σ_a)==streaming(σ_b)` for wildly different σ, both ≤72 ULP rel~1e-16 — pure streaming is GENUINELY σ-free; ERR-058/#195's σ-independent Carlson seed is what licenses it; the decomposition file's "3-13% wrong at σ=0" TOP docstring was STALE). CHOSE (a) done correctly: NEW named σ-free `LossRepresentation.streaming_action(psi) = loss_action(_zero_sigma_for(psi), psi)` + `streaming_action_transpose` on the Protocol + `_LossRepresentation` base (the base's `loss_action`/`_transpose` abstract sigs under `if TYPE_CHECKING:` to dodge `reportRedeclaration`); `StreamingOperator.apply`→`streaming_action`, `apply_transpose`→`streaming_action_transpose`; **`sigma_t` FIELD REMOVED** from `StreamingOperator` (Pattern 4); `loss_action` UNTOUCHED (so `InvertibleOperator.apply`=`loss_action(σ)` is BYTE-IDENTICAL — the sweep/solve path untouched; σ single-sourced from C's diagonal). MEASURED DRIFT: pure-L matvec leaf CART 16 / SPH 10 / CYL 8 ULP bulk, **boundary STRICT 0**; the **(L+C) composite matvec 0 ULP bulk+bdry** (byte-id). Migrated ~106 ctor sites across ~20 test/diag files (drop σ arg). NEW catcher C1 `test_pure_L_sigma_free.py` (`@foundation`, 9P: σ-freedom DIRECT `L.apply` Mode-11 + no-σ-surface + the Mode-11 σ-leak MUTATION TEETH); re-baselined B's `TestSubtractiveDefinition`(→nULP 256)/`TestT4b`(snapshots→`_PURE_L_NULP=256`)/`TestResolutionADifferentFromPriorWrong`(→`TestPureLIsLossActionAtZeroSigma`, `L.apply==loss_action(0)` byte-exact)/`test_loss_action_convention`(→`(L+C).apply==loss_action` byte + `L.apply≈loss_action−C·ψ` FP)/`test_invertible_operator`(nULP 256 bulk, bdry strict); RE-CAPTURED the 3 `bc_extraction_2d_baseline/*.npy` foundation snapshots to pure-L (verified `==(L+C)−σ·ψ` ≤64 ULP, strict gate RESTORED). pyright **2297/19 = EXACT baseline 0 net-new 0 type:ignore**; broad regression `-O` **EXACTLY the 7 baseline reds (#250 SPH×5 + #232 mu_y×2, confirmed identical on the `9316321` worktree), 0 non-baseline**; C2/C3/C4 + the §D MMS/k_inf/SI-rate backstop GREEN; G-adjoint reciprocity 12P (the `apply_transpose`→`streaming_action_transpose` path); Sphinx -W `build succeeded`. NEW Sphinx `:label: streaming-action-pure-l` (Eq `M(σ)ψ=streaming_action(ψ)+σ⊙ψ ⟺ streaming_action=loss_action(0)`) + `vv-status documented` + archivist `.. todo::` in `operator_algebra.rst`. NO new ERR (value-preserving; the stale-docstring fix is a doc correction). NO algebra-of-record Branch-1 manifest (type-surface + named-primitive carve over an existing verified discretization; σ-freedom pinned by C1, affine relation by the decomposition gate — S4.5/S8a posture). Closeout = `issue_257_s8b_pure_L_closeout.md`. archivist DISPATCH owed (consolidated S3b/S5/S6/S8b pass). ⭐ **POST-REVIEW (committed):** elegance **PASS-WITH-NITS** (both applied: function-local `import numpy as _np`→module `np`; the `.npy` capture-branch now ENFORCES the `(L+C)−σ·ψ` structural ground at capture-time — no self-referential freeze, verified the re-save is byte-stable); qa **SUPPORTED** (mutation-grade, all 5 claims independently verified — esp. CLAIM-1 `(L+C)` composite byte-identity via the `9316321` worktree + CLAIM-4 the 3 `.npy` == `composite−σ·ψ` ≤64 ULP, a load-bearing not-cosmetic re-baseline). Main-agent L12: C1/SI-rate/ERR-026-curvilinear-aniso re-run GREEN. ⚠ qa flag: leaf drift CYL **117 ULP** on 2G-het random-ψ (the implementer's "≤16" understated it; still FP-reassoc on the zero-caller leaf, under the 256-ULP test bound).

- **S8c ✅ DONE** (`35f5612` + matrix `08d0fd0`) — honest `@overload` typing of the `Fission`/`Scattering` `apply` fibration; **S8 COMPLETE**. The two operators are HETEROMORPHIC multi-carrier (`@singledispatchmethod`: ScalarFlux→ScalarSourceSink, AngularFlux→AngularSourceSink, TimedFullField→FullField, ndarray→ndarray, HarmonicMomentField→AngularSourceSink), but inherited the mixin's nominal `apply(x:V)->V` endomorphism AND the dispatcher base raises → `singledispatchmethod[NoReturn]`: callers saw `NoReturn` (poison) + the override + every `@apply.register` arm errored (fission 8, scattering 9 baseline dispatch errors). **Pattern M (the fix, RUNTIME BYTE-IDENTICAL):** rename dispatcher `apply`→`_apply_impl` (`-> Any` base so `.register` accepts real-typed arms), keep arms+bodies at natural indentation, append `if TYPE_CHECKING: @overload def apply(...)-><Out> ...` (one/carrier) `else: apply = _apply_impl` (public `apply` IS the SAME singledispatchmethod object, aliased — `__dict__['apply'] is ['_apply_impl']`→True). ⭐ **Pattern M chosen over the TYPE_CHECKING/else-split (Pattern C)**: master standard "reveals intention" (Beck rule 2) > "fewest elements" (rule 4) — C buries ~150/215 lines of source-assembly math in `else:`; M extends the mixin's EXISTING `if TYPE_CHECKING: def apply(self,x:V)->V` idiom (internally-consistent precedent; FIRST `@overload` in `orpheus/` → document the idiom in the archivist pass). Corrected stale dispatcher docstrings (fission's nonexistent AngularFlux arm; scattering's nonexistent ndarray arm + missing HarmonicMomentField; both →TimedFullField outputs now →FullField per S8a). Gate **C6**: `_c6_static_typing_pins` (pyright-only `assert_type`/carrier, MUTATION-VERIFIED teeth — breaking an overload reds the pin) + `test_c6_apply_dispatch_parity` (runtime, Mode-11 sentinel-confirmed on the aliased public `apply`, Mode-8 `pytest.fail`). 3 honest `cast()` (1 prod scattering `psi.bulk`→`AngularFlux|HarmonicMomentField`, 2 test) over the **#262** root gap; **0 net-new `type:ignore`**. pyright **2297→2282 (NET −15)**: removed the 16 NoReturn-poison errors; removing the mask UNMASKED **3 pre-existing latent errors in a NON-COLLECTED one-shot capture script** (`_fixtures/wave_t_t3/_capture_pre_t3_snapshots.py:191/204` — the bulk-accessor under-typing root, NOT regressions, NO runtime defect) → **filed #262** (`FullField.bulk:BulkField` too broad + `AngularFlux.integrate_angular()→object` + `build_aniso_source` return union; retires the 3 casts on fix). Runtime BIT-IDENTICAL: operator+C6 111P, §D MMS/k_inf/SI-rate 77P+2xfail(#195/#252). elegance **PASS-WITH-NITS** (decisively endorses M; nits: document "Pattern M" [owed→archivist pass] + #262 [filed]; optional C6-runtime HarmonicMomentField case skipped — windowing tests cover it), qa **SUPPORTED** (runtime bit-id + the −15 reconciliation + C6 mutation teeth, all independently re-verified; reconstructed baseline via `git show HEAD:` NOT stash). archivist rich-narrative + "Pattern M" doc DEFERRED into the consolidated S3b/S5/S6/S8b/S8c pass.

- **S9 ✅ DONE** (`7180f72`/`bc1ebec`/`aee60d4`, NOT pushed; premise corrected — boundary slope sub-floor, no value gate; LD 2nd-order-at-boundary LOCKED; #263 filed; #256 fork #1 resolved) — see the ⭐ CURRENT STATUS block. **NEXT = only deferred tails: S10** (χ-simplex at `Mixture.chi`, DEFERRED end-of-plan); **S11** (`Field`→L2, UNDER REVIEW). See the S10/S11 entries below.

(historical S8 detail follows.) Drop the `(L+C)−C` fold in
`StreamingOperator.apply` (`operator.py:417-418`): StreamingOperator computes pure streaming (σ-free);
`C` is the shared `MultiplicationOperator` (S3b); the model-specific loss representation is `L+C` (the
sweep, which keeps σ in its cell discretization). ⭐ S8 ABSORBS (from S4/S4.5): (a) re-point ALL operator
leaves to `apply(x: FullField) -> FullField` (timeless codomain — drop `history_depth=psi.history_depth`;
the driver `advance`s the timeless result into the timed state, the cofree decouple); (b) widen the
threaded helpers — `_require_typed_composite` (1 helper/4 sites) + the **14** `loss_action`/`_transpose`
defs + **15** `psi:"TimedFullField"` param sites in `loss_representation.py` (the `LinearOperator[V]`
invariance seam at `operator.py:641/794` RESOLVES once both ends read `FullField`); (c) resolve the
Scattering/Fission fibration finding (`@singledispatchmethod` over carrier type). Gate: SI rhs / Krylov
matvec value-equal (the composition graph changes — verify the field-level `L+C−S−B`); `InvertibleOperator`
σ single-source consistent; zero net-new pyright (the seam errors RESOLVE, not offset). **THEN S9
(BoundaryMomentField, BEHAVIORAL, test-architect FIRST), S10 (χ at source, DEFERRED end-of-plan),
S11 (Field→L2, under review).**

⚠ **S7 closeout note:** the scipy boundary is now internal to `KrylovAcceleration`; there is NO public
operator→scipy adapter anymore (was `as_scipy_linop`, retired). If a future caller needs to expose a
standalone `LinearOperator` to scipy (e.g. a DSA preconditioner built as an operator, #2), the move is to
generalize `_as_scipy_linop` (accept an `op.apply` callable + optional `rmatvec`) — NOT to resurrect the
flat twin. The ravel helpers `_is_ravellable`/`_ravel`/`_unravel_like`/`_zeros_like`/`_l2_norm` are
white-box-imported by `tests/transport/test_timed_full_field.py:515-524` — keep those names stable.

⚠ **S6 archivist DEFERRED (consolidated pass):** three clean Sphinx STUBS now await the archivist — S3b
`:label:multiplication-operator-promotion`, S5 `:label:functional-category`/`production-rate-functional`,
S6 `:label:integral-kernel-category` (each + a `.. todo::`). The §5.5–5.7 taxonomy (Operator/Kernel/Functional)
is now COMPLETE + stable (S7 is plumbing; S8 reshapes streaming, not these), so the consolidated
multiplier-algebra + suffix-law narrative can land post-S6 (per the campaign schedule) OR fold into the
post-S9 pass. Not an orphan — tracked here + on #257.

⚠ **S8 ABSORBS the operator re-point — now to `FullField` (timeless codomain).** S8 restructures
`StreamingOperator`→pure L + the loss representation anyway, so it also: (a) binds ALL operator leaves
`LinearOperatorMixin["FullField"]` (incl. the S3b `MultiplicationOperator`); (b) makes operator OUTPUTS
timeless `FullField` (drop `history_depth=psi.history_depth`; the driver `advance`s the timeless result
into the timed state — the cofree decouple); (c) widens the threaded helpers — `_require_typed_composite`
(1 helper, 4 sites) + **14** `loss_action`/`_transpose` defs + **15** `psi:"TimedFullField"` param sites
in `loss_representation.py` (the `LinearOperator[V]`-invariance seam at `operator.py:641/794` RESOLVES
once both ends read `FullField`); (d) resolves the fibration finding for `ScatteringOperator`/
`FissionOperator`. The S4 re-point trial was reverted (masked +2/−2 pyright offset) — S8 does it
coherently against the concrete `FullField`.

⚠ **Archivist DEFERRED (not an orphan):** S3b shipped a clean Sphinx STUB
(`:label:multiplication-operator-promotion` + `.. todo::`); the FULL multiplier-algebra narrative is
scheduled for the post-S6/S9 archivist pass (per the Verification section below + the campaign plan),
to avoid doc churn while S6 (FissionOperator composition) + S8 (StreamingOperator → pure L) reshape the
surrounding algebra. Tracked here + on #257.

## Context

`apply(self, x: V) -> V` over a bare `Generic[V]` was diagnosed (by the user) as unexplanatory —
"if type hinting doesn't help, it's compliance theatrics." That thread surfaced that grand report
§5.5–5.7 already specifies the target the code never reached: **cross-sections are `CoefficientField`s**
(a `Field` sibling), and **operators are those fields promoted** (`C = MultiplicationOperator(σ_t)`,
the multiplier-algebra embedding `M: L^∞ → B(L²)`). The opaque generic becomes the honest
**`TransportState`** bound. This re-scopes #256 (steps 1+3 are the foundation that stands); the rest
(scipy single-source, `Functional`, `BoundaryMomentField`) re-homes into the right taxonomy.

## The operator model (user-confirmed) — what is MODEL-SPECIFIC vs SHARED

- **Model-specific** (live in `sn/`, and each future model — CP/MoC — has its own):
  - `StreamingOperator` = **L** (pure streaming, `Ω·∇`; σ-free).
  - The **loss representation** = **L + C** = the `InvertibleOperator` / the sweep `(L+C)⁻¹`
    (`orpheus/sn/loss_representation.py`). This is where σ legitimately lives (the cell optical depth).
- **Shared across models** (live in `transport/`):
  - `CoefficientField` (+ `CrossSectionField`, `SpectrumField`), `MultiplicationOperator`,
    `C = CollisionOperator = M[σ_t]`, `ScatteringOperator`, `FissionOperator`, the `Functional` category,
    `TransportState`.

```
Field (numerics base = plain vector space; FluxRole adds the affine torsor)
├── CoefficientField (transport/; CoefficientRole = base algebra, NO affine gate)
│   ├── CrossSectionField  (cone σ≥0, units 1/cm)
│   └── SpectrumField      (probability simplex Σχ=1, dimensionless)
LinearOperator
├── MultiplicationOperator (transport/; = M[f], promotion of a CoefficientField)
│   └── CollisionOperator   (C = M[σ_t]; today ALREADY `σ[None]*ψ`)
├── StreamingOperator (sn/; = L, model-specific) ; loss rep = L+C (sn/, model-specific)
└── IntegralKernelOperator (transport/) → ScatteringOperator / FissionOperator
```

## Embedded census facts (so this plan regrounds without re-exploring)

- `CollisionOperator.apply` (`orpheus/sn/operator.py:592`) **is already** `self.sigma[None]*psi.bulk.values`;
  `.solve` (`:620`) `q/σ`; `apply_transpose`=`apply` (self-adjoint, currently untyped). The cleanest M[σ] target.
- The **C-fold** is in `StreamingOperator.apply` (`operator.py:417-418`): `loss_action(σ,ψ).bulk − σ[None]*ψ`
  — i.e. it computes (L+C) then subtracts C to *present* L. Target: StreamingOperator computes pure L
  (σ-free); the (L+C) action stays in the loss representation. `InvertibleOperator.apply` (`:884`)
  single-sources σ from `self.diagonal.sigma`.
- The **σ `np.ndarray` thread**: `LossRepresentation.loss_action(sigma, psi)` — 8 sites in
  `loss_representation.py` (`:252,272,429,442,773,933/947,1126/1147,1397/1418,1667`) + 4 callers in
  `operator.py` (`:417,461,884,903`) + `StreamingOperator.sigma_t` field (`:306`). σ in the loss-rep
  SWEEP is irreducible (cell discretization); σ in StreamingOperator (pure L) is removable.
- Primitives (`numerics/operator.py`): **only** `RankOneOperator` (`:1554`, fission, 1 prod site
  `fission.py:218`) and the σ-multiply are promotion targets. `DiagonalOperator` (`:1436`) has **NO prod
  caller** (`from_measure` `:1501` is the angular-W factory, no live caller). The 4 BC primitives
  (`Permutation`/`IncomingOrdinateMask`/`PeriodicWrap` + their `TensorProduct` folds, all in
  `boundary_realizer.py`) and `SumOfTensorProductsOperator` are basis-maps/projections → **stay flat**.
- `ScatteringOperator` already delegates to typed `mat_xs.apply_legendre_scattering_moments`/
  `apply_p0_in_scatter`/`apply_n2n` (`scattering.py`) — mostly a re-presentation, not a rewrite.
- `FissionOperator.kernel` (`fission.py:168`) = `RankOneOperator(χ,νΣf,axis=0) & IdentityOperator()`;
  `apply` is `@singledispatchmethod` over {TFF, ScalarFlux, ndarray}.
- `MaterialXSField` (`sn/material_xs_field.py`) = bare `@dataclass` (NOT a Field); one prod builder
  (`solver.py:858`). Views: `total_cross_section`/`absorption`/`fission_production`/`emission_spectrum`
  (ng,*spatial) + `sig_s_legendre`/`n2n_matrix` (group→group matrices). σ_s = a Kernel, not a Coefficient.
- `Functional` candidates: `iteration.py:251/272` estimators; the per-cell production rate is the
  anonymous `inner` einsum in `RankOneOperator.apply` (`operator.py:1671`).
- `units.py`: 4 signatures (`:130,136,140,144`); add `CROSS_SECTION_UNITS = UREG.cm**-1` (5th).
- Layering: `numerics ↛ transport` holds (zero `from orpheus.transport` in numerics) — `Vector`/`V` STAY
  in numerics; `TransportState` + coefficient leaves + `MultiplicationOperator` go in `transport/`.
- §4158 `.multiplication_operator()` is grand-plan VAPORWARE (`CriticalityEigenproblem` doesn't exist)
  — **no code rename needed**; reserve the name (the future eigen-iteration verb is `iteration_operator()`).
- `as_scipy_linop` (`operator.py:1680`): ZERO prod callers; live scipy boundary is the inline closure
  `iteration.py:744-766`.

## Staged implementation (each stage = one reviewed [elegance + qa] + gated commit)

**S1 — `CrossSectionField` + `1/cm` unit + `CoefficientRole`. ✅ DONE (`0565a24`, amended).**
Pure addition. `units.py` 5th signature `CROSS_SECTION_UNITS`; `transport/fields/` `CoefficientRole`
(base plain-vector-space algebra, NO affine gate), `CrossSectionField(CoefficientRole, ScalarField)`
(cone). Intrinsic-property tests: the cone closure (σ+σ ✓, λσ≥0 ✓, σ=0 origin), the role algebra
differs from `FluxRole` (no `flux+flux→TypeError`-style gate).

⚠ **`SpectrumField` DROPPED (user decision, 2026-06-19) — do NOT reintroduce in S1.** The per-cell
`emission_spectrum` view stores χ=0 in non-fissile cells (`mixture.py:192`), so a strict per-cell
simplex field is wrong; and a `SpectrumField` would have NO native downstream behaviour beyond the
simplex invariant. **The simplex (Σ_g χ_g=1, χ≥0) is a property of the SOURCE** (`Mixture.chi`, the
per-fissile-material spectrum that generates the broadcast) and is **enforced there** — deferred to
S10 (end-of-plan, decide-with-hindsight). `DIMENSIONLESS_UNITS` also dropped (no consumer). χ stays
the raw per-cell broadcast through S6 (its simplex guaranteed upstream).

**S2 — `MaterialXSField` typed accessors.** Bit-identical. Typed field accessors alongside the raw
ndarray views; `.values` bit-equal; existing consumers untouched this stage.

**S3 — `MultiplicationOperator` (transport/) + promote `CollisionOperator`. FOLD DiagonalOperator
(user decision 2026-06-19).** Two layers (the numerics↛transport boundary forbids a single class):
- **numerics**: generalize the (dead, zero-prod-caller) `DiagonalOperator` from "1-D weights on one
  axis" → an N-D coefficient on a *sub-product* of axes, broadcast over the complement
  (`np.expand_dims(coeff, broadcast_axes) * carrier`) — the shared broadcast ENGINE. Subsumes the 1-D
  case. Migrate its TEST callers (L20; production callers = 0). NO speculative axis modes — exactly
  the sub-product σ_t needs.
- **transport** (`transport/`, the first operator there): `MultiplicationOperator(CrossSectionField)`,
  the §5.7 promotion, DELEGATES the raw broadcast to the numerics engine on `ψ.bulk.values` and does
  the typed codomain (flux→`AngularSourceSink`; `solve`=`/`→`AngularFlux`; `apply_transpose=apply`,
  self-adjoint). `CollisionOperator` becomes one (or a thin subclass keeping its `+`-dispatch).

Rewire `solver.py:217/925/1000` to `MultiplicationOperator(mat_xs.total_cross_section_field)` (the S2
accessor). `L + M[σ]` MUST still assemble `InvertibleOperator` (the `+`-dispatch). **Law-suite tests
(the math concept's intrinsic properties):** `M_1=I`, `M_0=ZeroOperator`, `M_{af+bg}=aM_f+bM_g`,
`M.H=M`, `spec=ess-range` (invertible iff min|f|>0), and the homomorphism `M_f∘M_g=M_{f·g}` tested at
the VALUES level (σ·σ has units cm⁻² — the units-grading that deferred field·field `*`).

⚠ **GATE = principled-equivalence, NOT forced bit-identity (user, 2026-06-19).** The generalized
`expand_dims(coeff, axes)·x` reduces to `sigma[None]·ψ` for σ_t (likely bit-identical), but the gate
is the vv 3-criteria (named principled intermediate + STRUCTURALLY-INDEPENDENT reference [the
multiplier laws, `k_∞=νΣf/Σa`, the streaming-equilibrium analytic] + FP-bounded drift) — ACCEPT a
non-bit-identical result that arises from the more principled construction; do not force-fit 0 ULP.
**Dispatch test-architect FIRST** (operator-algebra carve crossing numerics↔transport — the MUST
proactive trigger): the values-level broadcast oracle, the law-suite, the `L+M[σ]≡InvertibleOperator`
assembly check, and the legacy-vs-new convention crosswalk.

**S4 — `TransportState(Protocol)` (transport/, refines `Vector`).** Bit-identical annotations. Names
`.bulk`/`.boundary`; re-points the transport/SN operator-leaf annotations to read like the domain;
`numerics` keeps `Generic[V: Vector]` (L0-ndarray + scipy path). Intrinsic-property test: `np.ndarray`
is NOT a `TransportState`; `TimedFullField` IS.

**S5 — `Functional` category.** Bit-identical reductions. keff/production estimators + the per-cell
`ProductionRateFunctional` as a `Functional` `evaluate(x)->scalar/field` (§5.6 suffix; #256-step-5).
Test: a `Functional` is NOT a `LinearOperator`.

**S6 — `FissionOperator` as `IntegralKernelOperator` = `M_χ ∘ ProductionRateFunctional ∘ M_νΣf`.**
⚠ **May be principled-not-bit-identical** (the composition changes the reduction order vs the fused
einsum) — judge by the 3-criteria, keep `RankOneOperator` as the structurally-independent cross-check.

**S7 — scipy single-source.** Bit-identical. Route the inline closure (`iteration.py:744-766`) through
one adapter; retire/fold the zero-caller `as_scipy_linop`. #256-step-4. Gate: keff unchanged to solver
tol; `from_flat(to_flat(x))==x`.

**S8 — restore `StreamingOperator` to pure L; loss rep = L+C. BEHAVIORAL.** Drop the `(L+C)−C` fold:
StreamingOperator computes pure streaming (σ-free); `C` is the shared `MultiplicationOperator`; the
loss representation (model-specific) is `L+C` (the sweep, which keeps σ in its cell discretization).
Retire σ from `StreamingOperator`'s surface; σ stays in the loss-rep sweep. **Dispatch test-architect
FIRST** (operator-composition seam). Gate: SI rhs / Krylov matvec value-equal (the composition graph
changes — verify the field-level `L + C − S − B`); `InvertibleOperator` σ single-source consistent.

⭐ **FOLD IN here (from S4): the `apply(x: TransportState)` operator re-point.** Since S8 already
restructures `StreamingOperator` + the loss representation, re-point ALL operator leaves to read the S4
`TransportState` Protocol atomically as part of this stage: `LinearOperatorMixin["TransportState"]` on
`StreamingOperator`/`MultiplicationOperator`/`InvertibleOperator` + override params/returns
`"TimedFullField"`→`"TransportState"`, and WIDEN the threaded helpers — `_require_typed_composite(field)`
(1 helper, 4 call sites) + the **14** `loss_action`/`_transpose` defs + **15** `psi:"TimedFullField"`
param sites in `loss_representation.py` (the invariance seam at `operator.py:641/794` closes once both
ends read `TransportState`). This is why S4 left the operators on `TimedFullField` (a `TransportState`-
bound `C` against a `TimedFullField`-bound `OperatorSum` is an invariant-mismatch pyright error). Gate
includes: zero net-new pyright (the seam errors RESOLVE, not offset); `apply(x: TransportState)` reads
like the domain.

**S9 — `BoundaryMomentField` + close the moment-state boundary drop. BEHAVIORAL.** #256-step-6, fork
#1 = general moment-tail (the #251 `boundary_face_layout` lever). The `(moment_buf, None)` drop is in
`apply_windowed`, not public `solve_moments`. **Dispatch test-architect FIRST** (scalar↔angular↔moment).
Gate: moment tensor byte-identical; new boundary block provably == old `None`.

**S10 — χ probability-simplex enforcement at the SOURCE (`Mixture.chi`). DEFERRED to end-of-plan
(user decision 2026-06-19); decide with hindsight when we reach it.** The simplex (Σ_g χ_g=1, χ≥0)
is a per-FISSILE-MATERIAL invariant on `Mixture.chi` (the `(ng,)` source that generates the per-cell
broadcast), NOT a per-cell field property — non-fissile materials store χ=0 (`mixture.py:192`).
Enforce it where χ is born (`data/macro_xs/mixture.py`), with the intrinsic-property test living
there. ⚠ This adds an invariant to EXISTING data → needs an upstream dependency audit (L20): every
existing fissile `Mixture` / test fixture must already normalise or it surfaces (correctly) as a
violation. At pickup, re-decide whether to enforce in `Mixture`, in `assemble_cell_xs`, or as a
standalone validator — with the benefit of having seen S3/S6 (where χ is actually consumed). This is
the home the dropped S1 `SpectrumField` was wrongly trying to be.

**S11 — relocate `Field` (data storage) out of `numerics/` (L1) to `transport/` (L2) or a new `fields`
layer. DEFERRED to end-of-plan + UNDER REVIEW (user, 2026-06-19: "why should Field go to L2? push to
the end with detail, review at the end").** Surfaced by the S4 layering investigation; decide with
hindsight. This is HYGIENE, ORTHOGONAL to everything else (it does NOT change `Vector`, the operator
promotion, or any solver) — so it is genuinely optional. The case, both sides, for the end-of-plan
review:
- **FOR:** (a) `Field` has ZERO numerics consumers — nothing in `numerics/` constructs or annotates it;
  it is ONLY ever subclassed from above (the L2 `transport` leaves). The `numerics/__init__.py` re-export
  of `Field` is the visible tell that it faces upward. (b) Module-over-ring frame: `Field` is OVERLOADED
  — simultaneously the additive base of the flux-MODULE and (via `CrossSectionField`) the underlying set
  of the coefficient-RING; co-locating it with its users + the #257 coefficient algebra clarifies both.
  (c) Co-locates a type with its only subclasses.
- **AGAINST:** (a) It does NOT dissolve the `Vector` barrier (the explorer + the forgetful-functor frame
  proved the barrier is forced by the carrier-generic ALGEBRA+DRIVERS staying in L1, NOT by `Field`'s
  location) — so the headline motivation (the user's original "is the hierarchy the root cause?") does
  NOT apply. (b) The documented L1/L2 principle ("L1 = knows-no-neutrons pure math; L2 = method-agnostic
  transport") is DEFENSIBLE for `Field`: it is "values + space + algebra", depends only on the
  unambiguously-L1 `FunctionSpace`, and "knows no neutrons". (c) `Field` is foundational — the move
  touches `tests/test_layer_imports.py` + every leaf's import path (mechanical but wide blast radius).
- **DECISION RULE at review:** relocate ONLY if, with hindsight from S1–S10, the module-over-ring
  overload (FOR-b) is judged a real clarity win that outweighs the churn; otherwise CLOSE as
  "principled-as-is" and keep `Field` in numerics. Check (L20) for any NEW numerics→Field coupling
  introduced by S5–S10 before moving. No cycle hazard found today (numerics drivers never construct a
  `Field`; `Field`→`FunctionSpace` is a valid downward L1 dep that survives the move).

## Verification (Cardinal Rule 1)

- **Per stage, judge the gate principledly — NOT a blanket 0 ULP.** For each touched output ask the
  vv-principles 3 criteria: (a) is every intermediate a NAMED principled quantity? (b) is it verified
  against a STRUCTURALLY-INDEPENDENT reference (the multiplier-algebra law, k_∞=νΣf/Σa, the
  RankOne cross-check, mpmath), not just old-vs-new? (c) is any drift FP-non-associativity bounded
  (`reduction-depth × ULP`)? 0 ULP is the *expectation* for pure re-typing/wrapping (S1–S5,S7); a
  principled reduction-order change (S6) may legitimately differ — accept ONLY if all 3 hold, then
  narrow the regression contract (`assert_array_almost_equal_nulp`) for that output with the documented
  justification. Be prepared for surprises and classify them, don't force-fit 0 ULP.
- **Every mathematical concept gets a test of its intrinsic properties** (user directive): the cone /
  simplex (S1), the multiplier-algebra laws (S3), the Functional-≠-Operator (S5), the
  `TransportState`-≠-ndarray (S4). These are the L0/L1 gates, not afterthoughts.
- Regression subset (route around the 7 baseline reds #250/#232 + #212): `.venv/bin/python -O -m pytest
  -q tests/sn/operators tests/sn/spatial tests/sn/sweep/core tests/sn/solve tests/numerics
  --deselect tests/sn/solve/test_keff_slab.py::test_heterogeneous_absolute_keff`.
- CLI `npx pyright` per stage: the opaque-generic surface shrinks; **zero new `# type: ignore`**.
- Sphinx clean; archivist documents the CoefficientField/promotion algebra after S6/S9.

## Execution discipline

Commit per verified stage (no push without an ask); `Co-Authored-By: Claude Opus 4.8 (1M context)
<noreply@anthropic.com>`; L28 (never `git checkout/restore/stash` on tracked paths); explicit paths;
never `.claude/skills/*`/`CLAUDE.md`/`.claude/rules`/`.claude/hooks`. First implementation action:
copy this plan into the repo + update issue #257 with the staged roadmap (durability before compaction).
