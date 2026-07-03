---
name: flux-to-sourcesink-operator-contract-frames
description: Frame attack on the universal "Operator: flux-STATE → rate-density SOURCE/SINK" contract across S/F/L (the balance (L+C−S−F/k)ψ=q). VERDICT the native structure is the LINEAR PART OF AN AFFINE MAP carrying the flux affine-space to the source vector-space (equivalently a bundle morphism E_flux→E_source with the 1/cm unit-gain as fiber transition); the codomain is a DIFFERENT typed space (vector role, has origin) from the domain (affine flux role, no origin). The shared contract should be Operator[Flux,SourceSink], NOT the endomorphic LinearOperator[V] apply(x:V)->V (whose TypeVar is in-code CONFESSED a lie via per-leaf @overload "NOT an endomorphism"). Fission = rank-1 DEGENERATE of scattering's analyse→act→reconstruct frame (both already satisfy IntegralKernelOperator; discriminator = middle rank). CORRECTS the brief's "fused-ndarray vs explicit-typed split" premise: ALREADY COLLAPSED in live code (kernel=frame.conjugate(Λ) typed OperatorProduct(R,Λ,M); both prod paths consume it; 0-ULP crosscheck is DEFINITIONAL not a twin-guard). Branch-verified on main worktree 2026-06-24.
metadata:
  type: project
---

# Flux→SourceSink operator contract — structural frame attack

Read DIRECTLY on the main worktree (L-005: file:line is ground truth, not Nexus).
This is the universal lift of `coefficient_field_promotion_frames.md` Frame-3
(fission `F = M_χ∘ProductionRate∘M_νΣf`) and `issue_208_delta_psi_affine_frames`
(affine torsor) to the WHOLE operator contract `Operator(flux)→source` across S/F/L.

## Verified code facts (the ground truth)
- `_flux_role.py`: flux = affine torsor (`flux+flux`→TypeError `:94`, `flux−flux`
  →Displacement, `affine_combination` Σλ=1). Domain has NO origin.
- `transport/source_sinks/*`: SourceSink = PLAIN vector space, `source+source`
  CLOSED ("a rate density, NOT a flux state"). Codomain HAS an origin (zero source).
- `numerics/operator.py:252`: shared base `LinearOperator(Protocol[V])`,
  `apply(x:V)->V` (ENDOMORPHIC), + `domain`/`codomain`/`capabilities`/`block_role`.
- `fission.py:466` + `scattering.py:1388`: `TYPE_CHECKING` `@overload` blocks whose
  comment CONFESSES *"NOT an endomorphism V→V (the mixin's nominal contract): it
  maps each input carrier to a DISTINCT output carrier."* The single TypeVar `V`
  spanning both role-spaces is the category error, patched per-leaf by overloads.
- `loss_representation.py`: `L.loss_action`→`FullField(bulk=AngularSourceSink,
  boundary=BoundarySourceSink)` — L also emits SOURCES, uniform with S/F.
- `integral_kernel_operator.py`: `IntegralKernelOperator(LinearOperator[V])`
  Protocol, REFINES (not disjoint), adds `kernel`. Fission `kernel`=`RankOneOperator
  (χ,νΣf,axis=0)&Identity` (rank-1 separable); scattering `kernel`=`frame.conjugate(Λ)`
  =`OperatorProduct(R,OperatorProduct(Λ,M))` (composite). Local ops (collision,
  Identity) have NO `kernel` → strict refinement, pinned by
  `test_integral_kernel_category.py`.
- `fission.py:267` `production_rate`→`ProductionRateFunctional` (the §5.6 middle
  factor, group-contraction to per-cell `p(r)`); rank-1 realization 0-ULP-pinned.
- `scattering.py:611-704`: `kernel = frame.conjugate(Λ)`; Λ (`LegendreMomentScattering`)
  carries REAL spaces (`domain==codomain==SphericalHarmonicSpace.from_L(L)`
  `:349-368`), the `None`-spaces `cast` RETIRED, `OperatorProduct` validates R∘Λ∘M
  natively. `:984-990` `build_aniso_source = (1/W)·self.kernel.apply(ψ.values)`.
  Windowed arm `:1378-1381` = `frame.reconstruct_after(Λ)` (M already done). SAME
  R,Λ,M. Crosscheck described IN-CODE as *"its DEFINITIONAL identity (it WAS a
  twin-guard between two parallel chains)"* — PAST TENSE.

## The three native frames (each with a FAIL-ABLE first test)

**Frame 1 — Affine map, linear-part-on-the-tangent (A.1).** Transport operators
are NOT linear `V→V`; they are affine arrows `T: A_flux → W_source` whose codomain
is the source VECTOR space, with linear part `dT: V_flux → W_source`. For the
genuinely-linear terms `T(0_A)`=zero source so `T=dT`, but the SIGNATURE crosses two
distinct typed spaces. The single `TypeVar V` identifies `V_flux` with `W_source`
ONLY because they share storage shape `(N,ng,*spatial)`. DISCRIMINATOR (L-002,
negative): `S.apply(ψ)` returns a source; `source+source` SUCCEEDS while `ψ+ψ`
RAISES — a `V→V` impl with one bound CANNOT express this. The killer:
`S.apply(S.apply(ψ))` is a TYPE ERROR in the affine frame (source is not a flux)
but type-checks under `V→V`. Native fix = `Operator[Domain,Codomain]` retiring the
per-leaf `@overload` stacks (4 scattering + 3 fission).

**Frame 2 — Bundle morphism E_flux→E_source (A.1).** Flux = section of state bundle,
source = section of PRODUCTION bundle (different fiber: the `1/cm` cross-section
gain = `ANGULAR_FLUX_UNITS→ANGULAR_RATE_UNITS`, the #208 operator unit-gain). `L`
= first-order differential (connection-like); C/S/F = zeroth-order bundle morphisms.
Balance = a vanishing-section equation in E_source. DISCRIMINATOR: `S.apply(ψ).UNITS
==ANGULAR_RATE_UNITS` while `ψ.UNITS==ANGULAR_FLUX_UNITS`, gain `codomain.UNITS/
domain.UNITS==1/cm` readable off `S` WITHOUT applying it — a `V→V` endomorphism has
no place for a domain≠codomain gain (same V, same units).

**Frame 3 — analyse→act→reconstruct, fission = rank-1 degenerate (A.2/A.7).** S =
`R∘Λ_{ℓ≥1}∘M`; F = `M_χ∘ProductionRate∘M_νΣf`. Fission IS the rank-1 degenerate of
scattering's frame: analysis basis = single production-rate measure, middle = 1×1
block, synthesis = broadcast-by-χ. Both satisfy `IntegralKernelOperator`; the
discriminator is MIDDLE RANK (rank-1 χ⊗νΣf vs sum-of-rank-(2ℓ+1) Σ_ℓ Pℓ⊗Σ_{s,ℓ}),
NOT a class difference. ALREADY BUILT (both expose `kernel`; fission also
`production_rate`). DISCRIMINATOR: `isinstance(F,IntegralKernelOperator)` and
`isinstance(S,...)` True; `isinstance(C,...)` (local multiply) False — the
strict-refinement test the existing category test pins.

## The brief-premise CORRECTION (the highest-value finding)
The brief asked "is the fused-ndarray-vs-explicit-typed split a defect?" — premise
is STALE. In live code BOTH aniso paths consume ONE typed composed operator
(`frame.conjugate(Λ)` / `frame.reconstruct_after(Λ)`), same R,Λ,M, 0-ULP crosscheck
DEFINITIONAL. "The fused composed operator IS the typed operator" is ALREADY TRUE.
The only residual is a CALL-SURFACE difference (full-angular applies to `ψ.values`
ndarray because the moment intermediate is internal; windowed materializes the typed
`HarmonicMomentSourceSink` because its iterate IS the moment tensor) — a legitimate
INPUT-SHAPE difference, not a divergent realization. The dissolution the brief sought
has shipped; the remaining move is one level UP: type the OPERATOR contract.

## Net verdict delivered to main agent
1. Native structure = linear part of an affine map (Frame 1) / bundle morphism
   (Frame 2); codomain = different typed space (vector role w/ origin) from domain
   (affine flux role, no origin). NOT cotangent-vs-tangent (no metric pairing) —
   affine-domain-vs-vector-codomain.
2. Balance IS a typed operator algebra = a SECTION equation; sources sum in the
   CLOSED source vector algebra (why the codomain must be the vector role) = q.
   S/F/L/C SHOULD share `Operator[Flux,SourceSink]` — the base/refinement machinery
   already exists; the `apply(x:V)->V` TypeVar is the endomorphic lie, patched by
   `@overload`. Genuine uniform abstraction = `Operator[Domain,Codomain]`.
3. Fused-vs-explicit = NOT a defect, ALREADY COLLAPSED; remaining elegance =
   type the Flux→SourceSink contract so the carrier-grid trip is legible for all
   S/F, fission as rank-1 degenerate of the SAME frame.

## Refuted (durable UNEXPLORED for this problem class)
- Homology — balance is a vanishing-section, not `∂²=0`; `B∘B≠0` (L-001).
- MPO/tensor-networks — fission rank-1, scattering rank-(2ℓ+1) FIXED, no bond knob.
- Category-theory adjunction (M⊣R) — captured by frame analysis/synthesis + `.H`;
  the affine-map/bundle-morphism frames ARE the concrete categorical content (functor
  flux-affine→source-vector cat), pinned with tests. No abstract-nonsense lever.
- Group/rep theory — SO(3) lives on Λ/quadrature, not the role-typing contract.
- Diff-geom/Christoffel — no curvature/1/r in flat S/F/L→source; curvilinear L only.

Cross-refs: [[coefficient_field_promotion_frames]] (Frame-3 fission composition, the
parent), [[issue_208_delta_psi_affine_frames]] (affine torsor, the domain-side mint),
[[issue_208_operator_algebra_frames]] (biproduct/dagger, the block algebra),
[[field_role_typing_faceflux_frames]] (the role grid as cochain decomposition),
[[projection_reconstruction_frame_pair]] (the analyse/synthesis frame faces M/R).
