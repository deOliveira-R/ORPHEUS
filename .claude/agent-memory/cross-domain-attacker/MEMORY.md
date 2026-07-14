# Cross-Domain Attacker — Memory Index

Slim index. Behavioral/process lessons live in `lessons.md` (read FIRST each
dispatch). The frame-trigger CATALOG lives in the `cross-domain-frames` skill
(Part A/B/C) and the promoted library kernel (Smell #16 four-shape detector +
transport-resolvent backbone) lives in AGENT.md — fire it from there, do NOT
re-derive. This index holds only (1) the lessons pointer, (2) git-true active
state, (3) durable design pointers (frame-matches that became real architecture).

## 1. Lessons (read first)

- [lessons.md](lessons.md) — 11 process lessons. Spine: a frame-attack's value is
  a concrete reformulation with a FAIL-ABLE first test OR a crisply-reasoned
  refutation, never a named-but-payoff-free frame. Refuted-frame reasons are
  first-class (L1); first tests must DISCRIMINATE (L2); Smell #16 four shapes
  (L3); property-vs-type is decidable by counting applied morphisms (L4); read
  the worktree not Nexus on a branch (L5); frame-leak naming (L6); the
  transport-resolvent backbone predicts cross-method layering (L7); "fully
  probes" is about linearity not degree (L8); a change-of-basis frame's
  OWNER+Galerkin/PG discipline are set by the operator's SYMMETRY, not the first
  caller (L9); a conserved-COLLAPSE splits by what's conserved — rate⟹average
  (G⁻¹M) vs probability/mass⟹marginalize (M) — fixing the morphism not a weight
  (L10); a coupled/nested block system is a FREE re-association of an existing
  biproduct (Mat∘Mat≅Mat) not a new object, and defer-until-2 counts coupling
  KINDS not instances-of-the-word (L11).

## 2. Active / in-flight state

**None.** The campaigns these attacks fed — SN operator-algebra / role-typing /
Wave-O, the eigenvalue-posing layering, the LD-on-the-DAG spatial closure, the
Variant-α Green's-function family, the Δψ affine-typing — are landed work; their
frame DESIGN VERDICTS are the durable pointers in §3 and the skill, not in-flight
code state. The field-typed-operator-algebra DESIGN MEMOS (§3, carrier-typing /
coefficient-promotion / property-vs-type) are GUIDANCE for the operator-algebra
line, decoupled from any single branch's merge status.

> Merge-status in memory goes STALE — a "branch X, NOT merged" note freezes
> mid-flight and the work lands in a later session. A frame VERDICT does not
> expire, but a code FACT it cites can; re-ground against the live worktree
> (L5) before acting on a file:line. Only #236 (spatial⊗angular product) is a
> known-open branch as of this curation.

## 3. Durable design pointers (frame-matches that became, or should become, architecture)

These are frame-attacks whose reformulation produced a real refactor or a
load-bearing design verdict — keep them as pointers; the detection HABIT is in
`lessons.md`, the structural CONTENT is here.

### Operator algebra (the SN/transport algebra spine)
- [assembly_mode_functor_dag_frames.md](assembly_mode_functor_dag_frames.md) —
  Phase-2b ASSEMBLY verdict (Cartesian; curvilinear OUT). Q1 `assembled≡probed` IS FUNCTORIALITY of the
  additive-monoidal `F:Op→Mat` (`as_matrix`=F-dense, `assemble`=F-sparse; the flat C-ravel order IS the
  Kronecker-factor/monoidal structure, per-carrier via `to_flat`; `SumOfTensorProductsOperator L=Σ_a D_a⊗Ω_a⊗I_g`
  already names it). Q2 per-octant `(L+C)` = the weighted INCIDENCE matrix of the sweep DAG: triangular ⟺
  acyclic (`levels`=P, `Pᵀ(L+C)P` strictly-lower-tri); #282 pole-seed = a BACK-EDGE (above-diagonal on the
  permuted sparsity); assemble per-ordinate blocks + lift (block-diag ⊕ `I_space⊗Λ_moment` scattering), NOT
  monolith; DSA `R·A·P` = clean triple product on the ANGLE factor. Q3 emission earns NO type (scipy COO/CSR
  carrier; sum-concat/scale/compose laws INHERITED from the operator algebra via the sparse-Mat functor — a
  COO-builder-with-laws twins OperatorSum); local-to-global map = the EXISTING carrier ravel
  (`FullField.to_flat`/`FlattenedOperator` template), SHARE it (Smell #16 shape-2), reify a type only for an
  UNSTRUCTURED (MoC/DG) consumer. Diffusion `LeakageOperator` FD (probed today) = the 2nd consumer / unify
  trigger; `assemble_ubld` already emits FEM local element matrices. Smell #16 shape-4 (assemble = 3rd path,
  fires before code). REFUTED: DEC-homology (`∂²≠0`, belongs to the trace algebra), circulant/FFT (sweep is
  triangular not periodic), MPO (fixed rank-d, no bond knob).
- [operator_inverse_taxonomy_frames.md](operator_inverse_taxonomy_frames.md) —
  #226 inverse-as-operator TAXONOMY verdict. TWO LAYERS: SUBSTRATE `LossRepresentation`
  (ONE object, MANY action-views — loss_action/transpose/sweep; "3 actions of 1 object" CORRECT here,
  it's not a morphism) vs OPERATOR (each ONE morphism ONE promise; A/A.H/A⁻¹ are DISTINCT objects
  sharing the substrate). BREAKS the Q1 slogan "fwd-sub & matvec are VIEWS of ONE operator" → "two
  operators, one substrate" (the slogan is what licenses the resolvent confusion). BREAKS Q5
  "as_matrix = co-equal 4th view" → it's a functor OUT (Op→Mat, dim-gated); the co-equal 4th is
  DenseInverseOperator (already = `direct_eigenvalue`'s dense A⁻¹F, un-promoted). FOUR inverse
  strategies (Sweep/Green-Richardson(=SI)/Green-Krylov/Dense), structure-keyed factory. NEW SMELL
  (Part C candidate, 1st sighting) "layer-confusion/vestigial-forward": `_GaussSeidelResolvent`
  (solver.py:338) apply=(L+C) "NEVER called" + solve=(L+C−B)⁻¹ ⇒ round-trip `inv.apply(A.apply(x))≠x`
  REDs. Frames: dagger-partial-inverse category / LU-triangular / resolvent-Neumann / matrix-free↔explicit.
- [operator_inverse_w1_w2_resolutions.md](operator_inverse_w1_w2_resolutions.md) —
  #226 self-audit W1/W2 RESOLUTIONS (next layer of ↑). W1 moment-emit = `OperatorProduct(P, A.inverse())`,
  P=`scattering.frame.analysis⊕Id_∂` a BLOCK COISOMETRY (NOT invertible, no round-trip); fusion=substrate
  `_SweepEmit` moment-mode (`apply(q,emit=…)`, §38); SweepOperator stays ENDOMORPHIC. §7 "SweepOp with
  moment emit" VIOLATES the plan's own §1 (a codomain-changing config is a different morphism). Honest
  invariant FACTORS (deforestation ∧ base-S-direct ∧ coisometry `M∘R=I`); a moment-residual proxy is
  category-confused (`Rφ≠ψ`). W2 G-S resolvent = `M.inverse()`, `M=(L+C−B_lower)` a REAL forward op
  (matrix-splitting `A=M−N`; B diagonal-free ∵ μ→−μ off-octant); REIFY `M.apply` (route a, reuses the
  schedule), NOT "preconditioner" (3rd mislabel). Gate: round-trip vs reified M.apply > FP-equiv on
  DIAGONAL cubature (Mode-9) > ρ(M⁻¹N). 3 §13 grain-fixes: §11.3⟺S-direct (same geo-scope; sphere seed
  is a relaxation-step arg pre-#282); `as_matrix` is TOTAL-with-resource-EFFECT not restriction-idempotent
  (Cockett-Lack CT fires here); `ResolventOperator(z)` is a FACTORY-level name (family R-identity binary),
  DROP the class = `(A−zI).inverse()`. NEW SMELL cand "codomain-changing emit/config" (1st sighting).
- [iso_source_frame_conjugation_unification.md](iso_source_frame_conjugation_unification.md) —
  VERDICT on "every iso source (P0/n2n/fission) = iso_frame.conjugate(K)=R₀∘K∘M₀ in EVERY model". 3-PART:
  (1) conjugation HOME correct + ALREADY REALIZED in SN — `full_scatter_kernel` (scattering.py:850) IS
  `frame.conjugate(Λ_{ℓ≥0}+N2n)` thru the HARMONIC frame, iso "made nice like aniso" w/ free transpose;
  the rank reduction is captured (Λ in moment space, M run once) — naive P₀⊗Σ_s0 TP loses it. (2) the
  separate rank-1 `iso_frame` is WRONG-OBJECT: it is NOT a new ConstantBasis, it is `angular_frame(0)` =
  GalerkinFrame(SphericalHarmonicBasis(L=0),S²) = the ℓ=0/V₀-trivial-SO(3)-irrep SUB-block of the EXISTING
  harmonic frame (M₀=∫=ℓ=0 analysis row, R₀=broadcast=P₀=1); minting ConstantBasis forks R∘M (Smell #16
  shape-2, ULP-drifts vs the 0-ULP fused kernel); diffusion/CP/MoC have M₀=R₀=Id ⟹ DEGENERATE (L-004,
  0 applied non-id morphism). (3) dyad-as-K=Id-degenerate-of-frame.conjugate(Id): R∘Id∘M vacuous middle,
  REJECTED — `outer(χ,⟨νΣf|)` (operator.py:1857) IS the single-mode normal form w/ its OWN dual-dyad
  transpose (operator.py:1831); the dyad↔frame relation is DEGENERATE-CASE (frame=STACK of dyads), carried
  by the shared IntegralKernelOperator Protocol, NOT a frame wrap. THE REAL WIN (model-indep claim TRUE +
  pays): collapse the 3 hand-rolled cross-model fission transcriptions (diffusion solver.py:260, CP
  solver.py:496 byte-identical `chi*rate[:,None]/keff`, SN fission.py:308 typed) onto the ONE #261-relocated
  FissionOperator (diffusion/CP use the bare-ndarray arm fission.py:530, M₀=R₀=Id) — Smell #16 shape-1
  ACROSS MODELS; discriminating test = diffusion sums axis=1 vs dyad axis=0 (functional.py:252) ⟹ array_equal
  REDs unless the relocation maps the axis convention. DISANALOGY (L-009/L-010): spatial-homog/energy-condense
  are NOT frame.conjugate — they are PG `frame.project`=G⁻¹M (analysis-only EXTRACTION, test≠trial, owned by
  no operator/eigenbasis); iso-source=`conjugate` (Galerkin, SH=scattering eigenbasis) vs homog=`project` (PG,
  solution-weighted): DIFFERENT frame VERBS, do NOT unify. SECOND iso/angular-frame memo (1st = harmonic_frame_ownership).
- [fission_rank1_normal_form_dead_functional.md](fission_rank1_normal_form_dead_functional.md) —
  fission `F=|χ⟩⟨νΣf|` IS the single-mode frame-conjugate NORMAL FORM (`frame.conjugate(I)` on a 1-mode frame
  = `RankOneOperator`; byte-proven operator.py:1811 `inner` line ≡ production_rate_functional.py:151) ⟹ S6
  "unfold F" is STRUCTURALLY EMPTY (`array_equal` discriminator, NOT allclose). NEW smell: a correctly-modeled
  §5.6 `Functional` category whose ONLY typed instance (`ProductionRateFunctional`) is production-DEAD while 5
  LIVE functionals stay untyped (`reaction_rate_density`, `compute_group_production/absorption_rate`, the 2
  estimator Callable aliases) — "abstract-category-top-down, instance-population-bottom-up-disconnected".
  VERDICT (option 2) = retire dead instance + its procedural-twin crosscheck (vv L11), re-seat Functional on
  the live pop as `ReactionRateFunctional(σ)` (prod/abs=instances), generalize to `BilinearFunctional` ⟨φ†,Aφ⟩
  (`_default_keff_estimator`=its φ†=1 degenerate; homogenization-PG adjoint-weights=the LIVE/IMMINENT consumer).
  CRUX: per-term-vs-k∞ analytic test IS genuine structural independence but needs the LIVE `compute_group_*_rate`,
  NOT the dead instance (which only carries a procedural twin).
- [issue_261_cross_method_operator_relocation.md](issue_261_cross_method_operator_relocation.md) —
  #261 relocate C/F/S out of `sn/` into `transport/`. H1 CONFIRMED (C's `sn_mesh` is session-born — the
  `MultiplicationOperator` base is ALREADY mesh-free; collapse to a held `full_field_space` like F/S, Smell
  #16 shape-1). H2 REFUTED (the W-D guard checks `(name,shape)` ONLY ⊉ object-identity; the essential
  invariant is the MIDDLE strength `object-id ⊋ geometric-consistency ⊋ shape-eq` — σ from `diagonal`
  fuses with geometry from `streaming` in the sweep; replace `is` with `C.coefficient.mesh is L.sn_mesh`,
  not `==`). H3 CONFIRMED (C = general multiplier `M[f]`, f=σ_t OR derived σ_r — `CrossSectionField`'s signed
  algebra is the referent; σ_r path is LATENT not live, #200 unwired). H4 = DO NOT mint `ReactionOperator`
  (zero shared apply body ⇒ Protocol not base; `IntegralKernelOperator` + `LinearOperator[Flux,SourceSink]`
  already ARE the shared contract — L-004 type-theatrics). `__add__`: one-sided dispatch on `StreamingOperator`
  FORCED (incl. a NEW `__radd__` override, ∵ base `__radd__` returns plain `OperatorSum` not `InvertibleOperator`);
  C/F/S become dispatch-free. Minimal symmetric shape = hold data + optional `full_field_space`, read mesh off
  carrier — the pattern the `transport/` base + F/S already follow; the carve makes C conform.
- [rep_role_grid_double_category_frames.md](rep_role_grid_double_category_frames.md) —
  the (Representation × Role) carrier grid IS a DOUBLE CATEGORY (objects=cells, horizontal=Rep-morphisms
  M/R Role-fixing, vertical=Role-morphisms C/Λ/F Rep-fixing); scattering `S=(1/W)·R∘Λ∘M` = the 2-cell
  `frame.conjugate(Λ)` (`scattering.py:696`), "natural in the untouched axis" = the interchange law, the
  `:587` windowed-vs-full bit-identical crosscheck = its COHERENCE WITNESS. VERDICT on the brief's 4 carrier
  type-machinery options: (d) the current MI leaves `Leaf(RoleMixin,RepBase)` is the UNIQUE permitted form
  (NOT a baseline) — Role changes the `__add__` SIGNATURE (torsor‖vector) ⇒ must be a CLASS (HARD-refutes
  phantom-param a,c via Frame-2 affine-torsor); Rep changes SHAPE ⇒ must be a CLASS (refutes b); MI is the
  only both-classes form, role-arith once + rep-data once, ALREADY in normal form (~16 thin leaves, no
  duplication). A `Field[Rep,Role]` phantom carrier is IMPOSSIBLE (params erased ⇏ specialize dunders). The
  (Rep,Role) parametrization the brief wants BELONGS on `LinearOperator[Din,Cout]` (Din/Cout ARE the
  carrier's Rep,Role) — already HALF-built; THE CARVE = finish it (M Role-generic, Λ Rep-generic), RETIRE
  the 7 `@overload` stubs. Distilled the class-vs-phantom-param decider into lessons L4 corollary.
- [operator_protocol_mixin_collapse_frames.md](operator_protocol_mixin_collapse_frames.md) —
  P4.5 verdict on the `LinearOperator(Protocol[Din,Cout])` / `LinearOperatorMixin(Generic[V,W])` split:
  the stateless dunder-installer Mixin is Smell #16 shape-1 (TWO co-required contracts for ONE concept,
  bridged by the `TYPE_CHECKING` apply stub at `operator.py:446`). Two-param ENABLES the collapse. Frame-1
  = dunders-as-default-methods ON the Protocol (operators ARE morphisms; hom-set carries the algebra) GATED
  on a pyright spike (contravariant `Din` in `__add__->OperatorSum[Din,Cout]`?); Frame-2 fallback = single
  `Operator(Generic[Din,Cout],ABC)`, drop variance (composers ALREADY invariant ⇒ declared variance unused;
  operators are NEVER ndarray ⇒ the Protocol justification is MISCOPIED from `Vector`). Frame-3 (fibered
  category) RULES leaf-vs-composite: type the LEAF `[Flux,SourceSink]`, composite `[FullField,FullField]` is
  DERIVED from source-fiber closure. RETIRE ranked: 7 `@overload` stubs > Mixin apply stub > Mixin class >
  variant-TypeVar machinery + `_AdjointOperator [W,V]` workaround. Run the variance spike FIRST.
- [issue_208_operator_algebra_frames.md](issue_208_operator_algebra_frames.md) —
  the algebra IS a dagger inverse biproduct category with a CO-EQUAL non-trivial
  boundary metric `G` (`†=G⁻¹AᵀG`): block matrices + adjoint-for-free are
  THEOREMS, `solve` is a partial-inverse morphism (not a string tag),
  triangularity is an SN/MoC-only refinement mixin. Smell candidate:
  metric-blindness.
- [bordered_vs_biproduct_facefield_frames.md](bordered_vs_biproduct_facefield_frames.md) —
  VERDICT on the `BorderedOperator`(Hyp A)+`FaceField`-codim1-dual(Hyp B) two-primitive proposal for
  `V=V_bulk⊕V_trace⊕V_sd`. A is an OVER-ABSTRACTION: "append rows" fuses 3 structures w/ different
  invariants (discriminator=corner block+metric): (1) saddle-point/KKT/bordered (RQI corner=0 EXACT
  `eigenvalue.py:545`, indefinite, inf-sup; RQI is the ONLY instance→defer), (2) biproduct/dagger-trace
  (the composite; `A_ss=−B`≠0, PSD block-diag; "border⊣un-border"=the biproduct inject/project dagger
  pair ALREADY named), (3) Galerkin coarsening (condense/homog=lossy PG `G⁻¹M`, NOT exact Schur
  `D−CᵀA⁻¹B`). B is RIGHT as mimetic/staggered/Hodge (cell scalar↔face current) NOT cohomology (`ι_*ι*=id`,
  no d²=0 — L-001 holds); ONE face measure = |streaming-flux normal thru face|, vanishing tangent-to-face
  (`|Ω·n|` & `(1−µ²)|_{µ=±1}` are 2 instances — built twice today). But ψ½ STRAINS B: it is an ANGULAR
  trace (full-in-space, keyed by µ-level = `ι*_angle` on the `V_space⊗V_angle` product, already
  `SumOfTensorProducts`), NOT a spatial-incidence face → BAD first consumer. NOW: name the one measure +
  generalize `FaceLayout` to a typed key (kills the `StartingDirectionSpace` 2nd flat-layout, Smell #16
  shape-2). DEFER: structured spatial `FaceSpace` (1st consumer=RT0 mixed-form/CP interface currents);
  `SaddlePointOperator` (trigger=mixed-form 2nd instance, then retire RQI inline border onto it). SEQUEL =
  [[coupled-system-field-bc-frames]] (the sharpened coupled-2×2 refinement of this composite).
- [coupled_system_field_bc_frames.md](coupled_system_field_bc_frames.md) —
  VERDICT on the SHARPENED "augmented SN = coupled 2×2 `[[A_AA,A_AB],[A_BA,A_BB]]`, system=field+BC" (ψ½ ray
  = System B, transport = System A; sequel to [[bordered-vs-biproduct-facefield-frames]], C2 FaceField landed
  @ `4081c0d`). CRUX: the 2×2 is NOT a new object — the EXISTING 3-block composite biproduct RE-PARTITIONED
  (`Mat₂∘Mat₂≅Mat₄`, ⊕ coherently associative); the user CONVERGED onto the native biproduct frame, "different
  object" is true only vs the block-DIAGONAL strawman. G-adjoint composes block-wise FREE (G=diag(G_A,G_B),
  `A†` = `_AdjointOperator`'s A_bs↔A_sbᵀ at 2-subsystem granularity). DO NOT mint `CoupledOperator` (VIEW or
  Smell-#16 twin); minimal object = lift the `BlockRole` 3-enum FREEZE to the mesh-enumerated block SET
  (`_join_block_roles` already set-union) + NAME A_BA (the welded fold, Smell #16 shape-4). Q1 system=field+BC
  = biproduct-carrier + boundary block-ROW (NOT comma/slice — no morphism beyond); RECURSES, System B is a 1-D
  TWO-POINT BVP (r=R Dirichlet inward RIGOROUS + r=0 pole-reflection outward — "BC=r=R corner" INCOMPLETE). Q3
  defer JUSTIFIED (3 instances = 3 coupling KINDS: ψ½ linear/triangular/metric-adjoint, DSA linear/two-way/R⊣P,
  multiphysics NONLINEAR — no 2 match; DROP multiphysics, it under-reaches a Mat(𝒞) block matrix). Q5 over-reach:
  PSD metric EXCLUDES RQI=KKT (defer, mixed-form trigger), triangularity SN-only, linearity under-reaches
  multiphysics; DSA's R⊣P = OUTSIDE but the DEFINING trigger for the coupled-ITERATIVE machinery. Q4 guards
  retire via the fibration/mesh-block-list = FaceField §5.3 (NAME not new power). Pollination: the A_BA weld IS
  an un-named Schur elimination `A_AA−A_AB A_BB⁻¹ A_BA`. Watch: `carlson_inward_sweep_transpose` is EUCLIDEAN —
  the coupled seed G-adjoint needs V_cell (ERR-067 family). Distilled into lessons L-011.
- [issue_208_delta_psi_affine_frames.md](issue_208_delta_psi_affine_frames.md) —
  the affine/torsor + Banach-contraction + Krylov-dual triad for typing a
  solver's two "→0" quantities (increment Δx vs residual r). The affine frame +
  state-typed-increment anti-pattern are PROMOTED TO THE SKILL.
- [field_role_typing_faceflux_frames.md](field_role_typing_faceflux_frames.md) —
  the SN face flux is a DEC 1-cochain; the boundary trace is its restriction `ι*`
  ("absorption = identity" IS `ι_*∘ι*=id`); `C¹=C¹_int⊕C¹_∂` is the biproduct.
  G-S/Jacobi = back-edge split of the reflective coupling on the (octant×face)
  graph.
- [streaming_apply_transpose_frame.md](streaming_apply_transpose_frame.md) — `Lᵀ`
  = triangular-transpose (reverse DAG walk + transposed local coupling + face
  roles swapped) + Lewis-Miller adjoint-ordinate. THREE gotchas a spatial-only
  reverse silently violates (nested angular recurrence, pole seed→outer face,
  trace inflow↔outflow swap).
- [flux_to_sourcesink_operator_contract_frames.md](flux_to_sourcesink_operator_contract_frames.md)
  — the UNIVERSAL `Operator: flux-STATE→rate-density SOURCE/SINK` contract across
  S/F/L: native = LINEAR PART OF AN AFFINE MAP `A_flux→W_source` (bundle morphism
  `E_flux→E_source`, `1/cm` gain = fiber transition); codomain is a DIFFERENT typed
  space (vector role w/ origin) from domain (affine flux, no origin). `apply(x:V)->V`
  is the ENDOMORPHIC LIE (the in-code `@overload` block CONFESSES "NOT an
  endomorphism") → genuine contract = `Operator[Flux,SourceSink]`, retiring the
  per-leaf overload stacks. Fission = rank-1 DEGENERATE of scattering's
  analyse→act→reconstruct frame (both already satisfy `IntegralKernelOperator`).
  CORRECTS the "fused-vs-explicit split" premise: ALREADY collapsed
  (`kernel=frame.conjugate(Λ)` typed `OperatorProduct(R,Λ,M)`; 0-ULP crosscheck
  is DEFINITIONAL not a twin-guard).

### Curvilinear / phase-space boundary
- [psi_half_seed_angular_trace_frames.md](psi_half_seed_angular_trace_frames.md) —
  #282/#280 R10 verdict: ψ½(μ=−1) starting-direction seed IS the ANGULAR-INFLOW
  TRACE of the (r,μ) phase-space rectangle (∂angular, dual to the existing ∂spatial;
  the M-M `(1−µ²)/r ∂_µ` is a first-order advection in µ with inflow at µ=−1). Third
  `FullField` block; ZERO-metric GHOST (angular flux (1−µ²)=0 at µ=−1 = the α_{1/2}=0
  fact; F2's α dual-weight self-refutes; reuse tangential-trace pseudo-inverse); NEVER
  a fabricated volume metric (over-weights a zero-measure µ-point). Euclidean coupling
  via `apply_transpose` (G only at outer wrap) ⇒ 2.5b reverse-scan = DAG reversal, 2.5c
  = existing G-wrap, no new metric math. SPHERE-only (µ=−1∉GL nodes); CYLINDER seed=ψ₀
  (η₀=−sinθ, τ_raw=0) = DEAD DOF ⇒ curvature+quadrature-keyed presence. F3 REJECTED
  (trace-as-interior-ordinate, shape earthquake); F4 CONFIRMS index-0 differential ⇒
  forward-substitutable. New smell candidate "metric-invisible-yet-active DOF" (1st
  sighting, distinct from zero-in-zero-out tangential slot).

### Eigenvalue / iteration layering
- [eigenvalue_posing_layering_frames.md](eigenvalue_posing_layering_frames.md) —
  K/α/source/transient = ONE generalized eigenproblem `Aψ=λMψ`, backbone =
  resolvent `A_loss⁻¹M`, layered leaves→posing→resolvent→algorithm; posing
  bifurcates into method-agnostic role-assignment (2a) + method-specific
  realization (2b).
- [power_iteration_vs_keigenvalue_morphism.md](power_iteration_vs_keigenvalue_morphism.md)
  — a power-iteration loop and an operator-triple loop are the SAME `fix(step)`
  combinator at two layers; the Protocol-opaque-resolvent layer is STRICTLY more
  general (admits monolithic-matrix resolvents) → it is the engine, the triple
  loop adapts in. (Deprecation arrows point toward the opaque interface.)

### Field-typed operator algebra DESIGN MEMOS (operator-algebra line guidance)
- [issue_226_container_algebra_design.md](issue_226_container_algebra_design.md) —
  uniform container-typed algebra: structural `Vector` Protocol +
  `apply(x:V)->V` (NOT per-op `Generic[F]`); the `(L+C−S−F−B)` algebra is an
  endomorphism on one carrier; 7 numerics primitives stay flat ndarray axis-ops
  (composite-state⇒container, single-axis⇒ndarray); scipy = serialization
  boundary.
- [issue_257_carrier_typing_layering_frames.md](issue_257_carrier_typing_layering_frames.md)
  — `Vector` is the forgetful-functor image (irreducible by the layer DAG, not a
  bad-hierarchy artifact); `apply:V→V` is a fibration cartesian morphism;
  `TimedFullField = Cofree(FullField, d)` forces the timeless-base split.
  Strong promotion candidates: forgetful-functor + Cofree-comonad (skill A.2).
- [coefficient_field_promotion_frames.md](coefficient_field_promotion_frames.md)
  — `f↦M_f` IS the multiplier-algebra embedding `M:L^∞→B(L²)` (laws→tests);
  `CoefficientField` is a commutative algebra+cone (NOT a flux torsor) splitting
  into `CrossSectionField` (cone, 1/cm) vs `SpectrumField` (simplex, χ); locality
  criterion = diagonal-symbol⇒MultiplicationOperator vs integrated-against-
  measure⇒IntegralKernelOperator. ⚠ flags the `multiplication_operator` name
  collision (rename the eigenvalue verb to `iteration_operator`).
- [spatial_order_type_vs_property_criterion.md](spatial_order_type_vs_property_criterion.md)
  — the DECIDABLE property-vs-type criterion (≥2 non-iso bases + an applied
  change-of-basis morphism); spatial order = PROPERTY (one basis, morphism=id),
  angular order = TYPE (the applied projection/reconstruction pair). The criterion
  itself is the durable kernel (now distilled into lessons L4).
- [harmonic_frame_ownership_funk_hecke.md](harmonic_frame_ownership_funk_hecke.md)
  — CAPSTONE (the OWNERSHIP question, one level up from the rest of this cluster):
  WHO owns the angular SH frame? VERDICT the claim "M is intrinsically scattering's"
  is CONFIRMED (not non-falsifiable): M is LITERALLY the change-of-basis into Σ_s's
  EIGENSPACE — Funk–Hecke (zonal Σ_s(Ω·Ω') ⟹ {Y_ℓ^m} diagonalize it, eigenvalues =
  Legendre moments = Λ's diagonal) + Schur (Σ_s·∈SO(3)-commutant ⟹ scalar-per-ℓ,
  (2ℓ+1) degeneracy); streaming = ℓ=1 ladder (PN ℓ↔ℓ±1 by Clebsch–Gordan), does NOT
  diagonalize → the asymmetry OWNS the basis to scattering. 4 falsifiers adjudicated:
  (A) ψ∈L²(S²) licenses INFINITE basis not the TRUNCATED ℓ≤L M/R (L=Σ_{s,ℓ} support)
  →derived; (B) aniso external source →not live (iso ℓ=0-embedded); (C) angular
  detector-response order L_d>L_scatter →THE genuine falsifier-in-principle but ABSENT
  (only ℓ=0 scalar via integrate_angular, off-frame); (D) PN closure →scattering-rooted.
  Architecture ALREADY CORRECT (HarmonicFrame in transport/frames/, S HOLDS it
  `scattering.py:510`, Λ-on-S/M,R-on-frame spectral-theorem split A=R·Λ·M, L=scatter_order
  `:545`). RELOCATION TRIPWIRE = 2nd consumer with L independent of scatter_order (C OR
  PN L_flux>L_scatter) ⟹ constructor home moves off S to neutral `Quadrature.angular_frame`
  (exists, anticipates it). DISANALOGY: energy condensation = NOT eigenbasis (G×G transfer
  has no symmetry/Funk–Hecke) → solution-WEIGHTED PG; ONE cause for the campaign's
  Galerkin-vs-PG split = angular Galerkin ∵ SO(3)-eigenbasis-orthogonal, spatial/energy PG
  ∵ no symmetry.
- [projection_reconstruction_frame_pair.md](projection_reconstruction_frame_pair.md)
  — (PRE-CAPSTONE, frame LANDED as `GalerkinFrame`) M/R asymmetry was a HALF-APPLIED
  P1.3 refactor; the three weight families are ONE convention datum, `(2ℓ+1)=4π·g_C⁻¹`
  DERIVED. Smell #16 shapes 1+2 (metric ‖ its inverse independently; `S0` realised 3×).
  POST-LANDING naming verdicts (durable): SH `L` keep param / fix prose `order`→`degree`;
  `Quadrature.angular_frame` names the AXIS not the basis (greppable family + 2nd-basis is
  a signature change not a rename); phase-axis ON the measure; Galerkin/PG discipline →
  the FRAME TYPE (superseded by the GalerkinFrame/PG hierarchy that landed).
- [homogenization_measure_derivation_frames.md](homogenization_measure_derivation_frames.md)
  — DERIVES (normal equations) that spatial XS homogenization is the WEIGHTED-L²(φV)
  GALERKIN projection onto the coarse P0/indicator basis: the measure φV is FORCED
  (not chosen), the `/Φ_R` denominator IS the inverse coarse Gram = the reconstruction
  dual factor (exact analogue of SH `(2ℓ+1)=4π·g_C⁻¹`). CORRECTS the prior "first
  oblique consumer" slip: in L²(φV) it is ORTHOGONAL (L=K), obliqueness was a dV-metric
  artifact → no oblique adjoint path. Mesh is NOT a Basis — it YIELDS an indicator
  basis-view + the existing measure-view (symmetric); membership = per-axis searchsorted
  on coarse faces (n-D-clean) = `partition_by`/`pushforward`. Group axis rides the
  TRAILING tensor axis through ONE frame (φ as `M_φ` multiplier, NOT a measure) →
  DiscreteMeasure 1-D + Basis ABC blast radius ZERO. Energy condensation = SAME shape
  (K on EnergyGrid); asymmetry explained by K's mesh-coupling; numerator = region-
  resolved ReactionRateFunctional (share it). ⚠ The "Galerkin-in-L²(φV) is the
  NATIVE reading" verdict is SUPERSEDED by [[condensation-nonnested-fractional-overlap-frames]]:
  the campaign REVERSED it to PG-with-flux-as-TEST (the falsifier = adjoint-weighted
  homogenization's bilinear ⟨φ†,Σφ⟩, φ†≠φ, no metric fold). The DERIVATION (φV
  measure forced, /Φ_R = inverse Gram, mesh-yields-views, n-D membership) STANDS.
- [condensation_nonnested_fractional_overlap_frames.md](condensation_nonnested_fractional_overlap_frames.md)
  — NON-NESTED condensation (421→WIMS69/172, straddling groups). RESOLVES the
  "flux=test-weight-not-measure" ruling vs the user's "spectrum-as-measure
  least-squares" floated alternative: condensation is PETROV-GALERKIN with a
  FRACTIONAL-OVERLAP trial table `T[g,G]∈[0,1]` (partition of unity), NOT a new
  LeastSquaresFrame. User's least-squares REFUTED 3 ways: (1) folds φ into the
  metric → breaks under adjoint-weighting (the bilinear, φ†≠φ); (2) for P0/indicator
  basis least-squares COINCIDES with the flux-weighted average (same PG map, not a
  different one); (3) LSQ's trigger (test=A·trial, dense SPD Gram, a SOLVE) absent
  (diagonal Gram → reciprocal). CONFIRMS the brief's partition-of-unity gram claim
  (`gram=analysis(reconstruction(ones))=Φ_G=Σ_g T[g,G]φ_g` = the right denominator,
  `frame.project` returns the rate-preserving average unchanged). Within-group flux
  model `f_{g,G}` = the TRIAL TABLE column-split (basis/partition geometry), NOT a
  measure, NOT the test weight — the three-way trial/test/measure separation holds.
  MOST ELEGANT = generalize `IndicatorBasis.evaluate` one-hot→fractional; nested =
  the {0,1} special case; ZERO new discipline/frame. Cross-method: ONE fractional
  `evaluate` serves non-aligned SPATIAL homogenization too; the overlap op = the
  discrete SOFT pushforward (one atom fans out). Refuted: LeastSquaresFrame,
  Galerkin-in-L²(φ), conservative-remapping/supermesh (P0-collapses to fractional
  overlap), optimal transport (no cost metric).
- [xs_coarsening_collapse_marginalize_vs_average.md](xs_coarsening_collapse_marginalize_vs_average.md)
  — CROSS-VERB synthesis (homogenize P3 merged ∥ condense P5 draft, "what is the
  SHAPE"). The two are ONE op `Collapse(axis,table,test_weight,normalize?)` =
  `PetrovGalerkinFrame.project`. LOAD-BEARING (L-010): the "1-vs-2-frame" and "χ
  production-weighted vs bare `@T`" asymmetries are NOT "same slot ± a weight" —
  (a) χ collapses DIFFERENT AXES (homog=spatial-AVERAGE of χ; condense=energy-BIRTH-
  MARGINALIZE of χ) and (b) two MORPHISMS: average=`G⁻¹·M` (project, conserves a
  RATE) vs marginalize=`M` alone (`chi@table`=analysis-without-G⁻¹, conserves
  `Σχ=1`; a weight=1 project would ÷ bin-COUNT and BREAK Σχ=1). ONE machinery
  `(weight,normalize?)`; exposing marginalize(=analysis)/average(=project)
  DISSOLVES the frame-count asymmetry. Q1: 3 layers — `project`✓ + axis-adapter✓ +
  channel-taxonomy DUPLICATED (`material_xs_field.py:398` ∥ `mixture.py:288`, the
  ONE thing to unify). Q4: declare Gram structure ON THE BASIS (Indicator→DIAGONAL,
  Overlap→POU, GEC→DENSE); `FrameBase.gram`'s 30-line prose-precondition on a 3-line
  body wants to be a TYPE (negative test: a DENSE stub makes `.project` RAISE, not
  silently return the wrong row-sum — same landmine as WeightedIndicatorBasis raising
  on its unbuilt synthesis side). Q2 HAZARD: n-D `OverlapBasis` is a 2^d-DENSE
  outer-product table (corner overlaps), INVISIBLE on the ordered ≤2-straddle energy
  axis → MUST 2-D-test before homogenization adopts it. Q5: measure-view & basis-view
  from OPPOSITE ENDS (`fine.as_measure`+`coarse.as_basis`); the overlap table is a
  BINARY (fine,coarse) frame-construction step (`GroupCondensation` IS the energy
  instance), NOT a unary view. FLAGS the brief's row (d) as an axis category error.
  Discriminator: `project(Σ@T) ≠ (project Σ)@T` (normalization keyed on one axis).
- [projection_discipline_hierarchy_frames.md](projection_discipline_hierarchy_frames.md)
  — the projection-DISCIPLINE hierarchy (one level UP from M/R): (a) "projector"
  is RESERVED for the idempotent `R∘Π`; Π=analysis operator (Christensen)/
  projection (FEM), NOT a projector; (b) Galerkin/PG are SIBLINGS not
  `Galerkin(PG)` (Liskov: Galerkin STRENGTHENS `Π*=R`), but NEITHER earns a
  class; (c) ONE generalizing object = Saad's `(K,L)` trial/test pair = ORPHEUS
  `(basis,measure)`; collocation=Dirac measure, LSQ=`A*K` (the frame's
  expressivity boundary — needs `A`); (d) discipline = a PROPERTY decidable as
  `measure==basis.canonical_measure` → RETIRE both ABCs (type-theatrics, morphism
  =id, zero subclasses), matches `frame.py:27-35` CAP_SOLVE iso/non-iso ruling;
  (e) condensation/homogenization PG-ness IS a non-canonical (spectrum/flux-
  weighted) MEASURE, not a method — and the §17+§18 second oblique instance that
  justifies the path. Extends the M/R pair's verdict 3.2.

### Unified Frame API (the projection/reconstruction machinery design)
- [unified_frame_api_design.md](unified_frame_api_design.md) — LIFTS the two derivation
  memos (projection-discipline + homogenization-measure) into 3 API verbs. THE DEEP OP =
  Christensen canonical-dual coefficient map. (1) `project_weighted(f,w)=G_w⁻¹·M(w⊙f)`
  (homogenize/condense; the `/Φ_R` IS `G_w⁻¹`; diagonal Gram ⇒ reciprocal not solve);
  (2) `conjugate(A)=R∘A∘M` typed face (+`reconstruct_after` sub-op for the windowed arm;
  needs NO 2-param Din/Cout split — just type Λ's spaces=basis_space, delete the
  `scattering.py:663` cast); (3) `analyze`=bare M keep. Q2 = Option-B per-basis Gram split
  (keep SH's measure-free reconstruct, 0-ULP). Q4 = `Frame(basis,measure,*,test:Basis|None)`,
  `is_galerkin=(test is basis)` PROPERTY, RETIRE the 2 PG ABCs (#268), build the SEAM not
  the PG cross-Gram path, `canonical_measure` NOT needed. Q5 LS = flag-and-defer (`test=A·K`
  non-diagonal Gram ⇒ CAP_SOLVE coeff space). flux-as-MULTIPLIER is durable (retires
  `flux_volume_measure` consumption seam, DiscreteMeasure stays 1-D); condensation = same
  verb + `EnergyGrid.indicator_basis` greenfield leaf.

### Variant-α Green's-function family + spatial closure
- [variant_alpha_family_hindsight.md](variant_alpha_family_hindsight.md) +
  [trajectory_resolvent_foreign_frames.md](trajectory_resolvent_foreign_frames.md)
  — the 6-geometry × 2-orbit-class family: top frame = fiber bundle
  (BaseAtlas/AngularFiber/ChordOracle); rank-N IS the bond dimension of an open
  MPO (user-confirmed); ready-now refactors = shared ChordOracle (with Nyström),
  unified `power_iterate` driver, single GreensResult; MPO/bundle wait for the
  N≥3 / instance-N+1 tripwire.
- [d5_trait_and_mms_frames.md](d5_trait_and_mms_frames.md) — the DD/Step-vs-LD
  scan-march discriminator is TRANSVERSE-COUPLING ORDER (0th-order face trace ⇒
  separable d-D closure vs 1st-order slope moment ⇒ irreducible (1+d)-block); name
  the trait for the SCHEME property (`transverse_coupling_is_facewise`), not the
  strategy. The multi-D LD MMS ansatz `[A+µx·B+µy·C]/W` activates the bilinear
  slope rows IF the slope drivers vary along their own axis (add x↔y-broken cross
  terms to defeat a same-sign slope-row bug).

### Diffusion integration (#290 — pre-implementation frame verdicts)
- [diffusion_integration_frames.md](diffusion_integration_frames.md) — the diffusion
  solver IS the transport algebra hand-inlined (`_matvec` = `(L_diff+C−S)φ`,
  `compute_fission_source` = rank-1 `F/k`); the carve COLLAPSES it onto the shared
  leaves, minting only 2 new objects. Q1 trace = **partial currents J±** (= ℓ=0
  half-range moment of the SN trace under the SAME `G_s=|Ω·n|⊙w` metric the
  `TraceSpace` already carries; inflow=BC-input/outflow=solution-output mirrors SN;
  albedo BC = A_ss `BoundaryOperator`, retires the `BC_REGISTRY` string dispatch;
  kills Smell #9 + #16-shape2 — `_boundary_gradient` IS the un-named DtN map). Q2 =
  option (c) `A_diff=L_diff+C` (L_diff is the FULL primitive, elliptic sibling of
  streaming L — the L-007 diffusion exception NAMED), `S`=K_iso, `F`=shared rank-1
  dyad; C=σ_t so in-group CANCELS via S (removal as theorem); `A_diff.H` FREE ⟹
  adjoint-φ* for #281; DSA doc independently commits to this (`operator_algebra.rst
  :4210,4296`). Q3 = FD-as-is honest (face-D `1/(3·AM(σ))` = harmonic-mean-D = RT0
  stiffness, bit-identical; H(div)/RT0 is native, mint only the boundary J_face
  now). Q4 = 2nd ScalarFlux consumer ⟹ promotes shared TransportMethod; borrows
  K_iso/rank-1 F/eigenvalue-engines/FullField-trace-triple. First tests DISCRIMINATE
  (3-group breaks `σ_s[::-1]`; α=0.5 breaks string BC; `direct_eigenvalue` needs the
  in-algebra `A_diff.inverse()`).
