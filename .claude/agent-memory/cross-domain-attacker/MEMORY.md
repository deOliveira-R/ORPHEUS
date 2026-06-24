# Cross-Domain Attacker — Memory Index

Slim index. Behavioral/process lessons live in `lessons.md` (read FIRST each
dispatch). The frame-trigger CATALOG lives in the `cross-domain-frames` skill
(Part A/B/C) and the promoted library kernel (Smell #16 four-shape detector +
transport-resolvent backbone) lives in AGENT.md — fire it from there, do NOT
re-derive. This index holds only (1) the lessons pointer, (2) git-true active
state, (3) durable design pointers (frame-matches that became real architecture).

## 1. Lessons (read first)

- [lessons.md](lessons.md) — 8 process lessons. Spine: a frame-attack's value is
  a concrete reformulation with a FAIL-ABLE first test OR a crisply-reasoned
  refutation, never a named-but-payoff-free frame. Refuted-frame reasons are
  first-class (L1); first tests must DISCRIMINATE (L2); Smell #16 four shapes
  (L3); property-vs-type is decidable by counting applied morphisms (L4); read
  the worktree not Nexus on a branch (L5); frame-leak naming (L6); the
  transport-resolvent backbone predicts cross-method layering (L7); "fully
  probes" is about linearity not degree (L8).

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
- [issue_208_operator_algebra_frames.md](issue_208_operator_algebra_frames.md) —
  the algebra IS a dagger inverse biproduct category with a CO-EQUAL non-trivial
  boundary metric `G` (`†=G⁻¹AᵀG`): block matrices + adjoint-for-free are
  THEOREMS, `solve` is a partial-inverse morphism (not a string tag),
  triangularity is an SN/MoC-only refinement mixin. Smell candidate:
  metric-blindness.
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
- [projection_reconstruction_frame_pair.md](projection_reconstruction_frame_pair.md)
  — the M/R asymmetry (M generic-by-codomain, R specific-by-name+carrier) is a
  HALF-APPLIED P1.3 refactor; cleanest frame = analysis/synthesis pair of a
  `4π`-TIGHT FRAME on S² (one `SphericalHarmonicFrame` owning Y,W,L + space-pair);
  the three weight families are ONE convention datum with `(2ℓ+1)=4π·g_C⁻¹`
  DERIVED. Symmetric fix: rename R→`MomentReconstruction`, give it `.domain`,
  read `(2ℓ+1)` from it. Smell #16 shapes 1+2 (metric ‖ its inverse stored
  independently; `S0` realised 3×). Unify-NOW gate = the §10 PN second instance.
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
  resolved ReactionRateFunctional (share it).
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
