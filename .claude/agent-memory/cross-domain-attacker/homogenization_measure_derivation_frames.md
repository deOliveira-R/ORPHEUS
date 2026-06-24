---
name: homogenization-measure-derivation-frames
description: DERIVED (not asserted) function-space structure of spatial XS homogenization / energy condensation on Frame(basis,measure). VERDICT the homogenization map H is the WEIGHTED-L²(φV) ORTHOGONAL (Galerkin) projection onto the coarse P0/indicator basis — the normal equations FORCE the measure μ=φV (flux·volume), NOT dV; the "/Φ_R" denominator IS the inverse coarse Gram G⁻¹=(R's dual factor), so H = G⁻¹∘M with M the un-normalized weighted moment the Frame's analysis face already computes. Galerkin-in-L²(φV) ≡ oblique-PG-in-L²(dV): SAME map, two readings; canonical one for the Frame = Galerkin (flux IS a measure factor, μ=φ·dV). Mesh is NOT a Basis: it YIELDS a coarse-indicator basis-view (P0) + a measure-view (volume_measure), symmetric with how it already yields volume_measure; cell membership = a restriction/pushforward π:fine→coarse built from coarse FACES (np.searchsorted is the 1-D shadow of a general cell-locate, n-D-clean via per-axis searchsorted on the coarse axis edges). Group axis rides the trailing tensor axis through ONE frame (analysis already einsum 'i,i...->...'); φ enters as integrand multiplier M_φ, NOT a per-group measure → DiscreteMeasure 1-D-weights constraint UNTOUCHED. Energy condensation = SAME shape, K=coarse-group indicator on EnergyGrid, mesh-decoupled (asymmetry explained: space's K lives on a mesh, energy's K is material-keyed). Numerator = region-resolved ReactionRateFunctional → share the primitive. Read on main, clean tree — file:line ground truth (L-005 N/A).
metadata:
  type: project
---

# Homogenization as a weighted-L²(φV) Galerkin projection on the coarse P0 basis

DESIGN VERDICT feeding the `.homogenize` slice (#36/#37, branch
`refactor/operator-inverse-algebra`). Read on `main`, clean tree —
`frame.py`/`basis/base.py`/`measure.py`/`space.py`/`functional.py`/
`material_mesh.py`/`production_rate_functional.py` file:line is GROUND TRUTH
(no Nexus/worktree staleness; L-005 N/A). Extends and OPERATIONALISES
`projection_discipline_hierarchy_frames.md` verdict-5 ("homogenization's
PG-ness IS a non-canonical measure"): that memo NAMED the measure as
flux-weighted; THIS memo DERIVES it from the normal equations and locates the
`/Φ_R` denominator as the inverse coarse Gram.

## The single load-bearing derivation (Q1+Q2 answered together)

H collapses fine XS `Σ_i` (i ∈ fine cells) to coarse XS `Σ_R` (R ∈ coarse
regions) preserving reaction rate. Treat `Σ` as a function on the spatial
domain. The coarse representation is **piecewise-constant on the coarse
regions** — the trial space is `span{1_R}` (the P0 / indicator basis). So H is
an APPROXIMATION: project `Σ ∈ L²(fine)` onto `span{1_R}`. The question is "in
WHICH inner product," and the reaction-rate constraint FORCES the answer.

**Normal-equations derivation.** Seek coarse coefficients `c_R` minimising the
projection residual in a weighted inner product `⟨f,g⟩_w = Σ_i w_i f_i g_i`.
The best `span{1_R}`-approximant satisfies the normal equations: residual ⟂ each
trial function `1_R`:

    ⟨Σ − Σ_c, 1_R⟩_w = 0   for every R,   where Σ_c = Σ_R c_R 1_R.

Because the `1_R` have disjoint support, the Gram `⟨1_R, 1_{R'}⟩_w = δ_{RR'}·(Σ_{i∈R} w_i)` is DIAGONAL, and the normal equations solve in closed form:

    c_R = ⟨Σ, 1_R⟩_w / ⟨1_R, 1_R⟩_w = (Σ_{i∈R} w_i Σ_i) / (Σ_{i∈R} w_i).

Match against the target `Σ_{R,g} = (Σ_i V_i φ_{i,g} Σ_i) / (Σ_i V_i φ_{i,g})`:
the weight is `w_i = V_i φ_{i,g}`. **The measure is φV (flux·volume), NOT the
geometric volume dV — DERIVED, not chosen.** dV would give the wrong
(volume-only) average and would NOT preserve reaction rate. The reaction-rate
preservation property `Σ_{R,g}·Φ_{R,g} = Σ_{i∈R} V_i φ_{i,g} Σ_{i,g}` is exactly
the statement `⟨Σ_c, 1_R⟩_{φV} = ⟨Σ, 1_R⟩_{φV}` — the projection preserves the
functional `⟨·, 1_R⟩_{φV}`, which is the region-resolved reaction rate.

**So H = G⁻¹ ∘ M:**
- `M = ⟨·, 1_R⟩_{φV}` = the un-normalized weighted moment `(MΣ)_R = Σ_{i∈R} V_i φ_i Σ_i` — EXACTLY what `Basis.analyze` computes (`base.py:124`: `(Mf)_k = Σ_n w_n φ_k(x_n) f(x_n)`, no Gram division). The numerator.
- `G⁻¹ = diag(1 / Σ_{i∈R} V_i φ_i)` = the inverse coarse Gram. The `/Φ_R` denominator. This is the RECONSTRUCTION dual factor: for the indicator basis the canonical-dual factor `d_R = 1/⟨1_R,1_R⟩_μ = G⁻¹` is the EXACT analogue of the SH basis's `(2ℓ+1) = 4π·g_C⁻¹` (`spherical_harmonic_basis.py:166`, `reconstruct` = `Σ (2ℓ+1) Y φ`). The SH dual factor is the inverse Gram up to 4π; the indicator dual factor is the inverse Gram, full stop.

H is therefore `frame.reconstruction ∘ frame.analysis` restricted to the
COEFFICIENTS (`G⁻¹·M`), NOT the full round-trip back to nodes — homogenization
WANTS the coarse coefficients `Σ_R` (one value per region), which IS the
`basis_space` element. The coarse XS field is the coefficient vector.

## Q2 verdict — Galerkin, not Petrov-Galerkin. The flux is a MEASURE factor.

The two readings are genuinely the SAME map:
- **Galerkin in L²(φV)**: trial = test = `{1_R}`, orthogonal projection in the
  flux·volume inner product. `measure = φ·dV`.
- **Petrov-Galerkin in L²(dV)**: trial `{1_R}`, test `{φ·1_R} ≠ {1_R}`, oblique
  projection in the plain-volume inner product.

They coincide because moving the φ from the test functional into the measure is
the identity `⟨Σ, φ1_R⟩_{dV} = ⟨Σ, 1_R⟩_{φdV}`. **The NATIVE reading for the
Frame is Galerkin-in-L²(φV)**, because `Frame.analysis` (`frame.py:142`) reads
the WEIGHTS off the measure (`self.frame.measure.weights`) — the machinery is
built to put the weighting in the measure, and the measure already owns the
inner product (`measure.py:242` `space.inner_product_weights = weights`). The
test space `{φ·1_R}` has no home in the current types (the basis is the test
space and it is `{1_R}`); the measure `φ·dV` has a home (it IS a
`DiscreteMeasure`). Encode the flux-weighting as the measure, and the
projection is Galerkin BY CONSTRUCTION (`measure == basis.canonical_measure` for
the indicator basis under that measure — see below).

**CORRECTION to the prior memo's "first genuinely-oblique consumer" claim.**
The prior `projection_discipline_hierarchy_frames.md` verdict-5 called
condensation/homogenization "the FIRST genuinely-oblique (L≠K) production
instance." THAT IS THE WRONG FRAME. In L²(φV) the projection is ORTHOGONAL
(L=K): the obliqueness was an artifact of insisting on the dV inner product. The
right statement: homogenization is Galerkin in a NON-CANONICAL (flux-weighted)
measure. The discipline criterion `measure == basis.canonical_measure` still
fires "non-canonical" (the indicator basis's canonical/self-dual measure is the
COUNTING or volume measure, not φV), so the prior verdict's CONCLUSION (PG-ness
= a non-canonical measure, retire the two ABCs, no class) STANDS — but the
"oblique" word was a category slip. It is Galerkin-in-a-reweighted-metric, which
is why it composes cleanly with the existing all-Galerkin Frame machinery
instead of needing a separate oblique adjoint path. This is the cleaner result:
NOTHING about the Frame's analysis/reconstruction machinery needs to change for
the oblique case because there IS no oblique case — only a different measure.

## Q3 verdict — the mesh is NOT a Basis. It YIELDS a basis-view and a measure-view.

The mesh carries BOTH roles, and they are DIFFERENT mathematical objects:
- **Basis role (trial)** = the coarse cells' indicator functions `{1_R}` + the
  FACES that define membership (which fine cell sits in which coarse region).
  This is the codomain coefficient structure (one coeff per region).
- **Measure role (test/integration)** = the cell VOLUMES → `dV`, already
  realised as `mesh.volume_measure` (`material_mesh.py:308`, `mesh.py:323/601`).

These are NOT the same object, so the mesh must NOT inherit `Basis` (that would
conflate trial and test — and `Basis` is the choice-free, measure-FREE half
`base.py:6`; a mesh that also carried the volume measure would violate the
basis's measure-freeness). The mesh is the COARSE geometry; it presents to
`Frame` via TWO accessors, symmetric with the existing one:
- `coarse_mesh.indicator_basis(fine_measure) -> Basis` — the P0 basis whose
  `evaluate(fine_centers)` returns the membership table `1_R(x_i)` (an
  `(n_fine, n_coarse)` 0/1 matrix). This IS a concrete `Basis`: `evaluate` =
  membership tabulation, `synthesize` = scatter coeffs back to fine cells,
  `analyze` = the region-sum moment, `reconstruct` = `G⁻¹·` region-sum, `space`
  = the coarse `RegionSpace` (`space.py:43` future-direction, now forced),
  `mass_matrix` = `diag(region volume-or-φV sums)`.
- `coarse_mesh.volume_measure -> DiscreteMeasure` — ALREADY EXISTS.

This satisfies the user's "the mesh is the basis, do not mint a free-floating
`IndicatorBasis`": the `Basis`-conforming object is NOT free-floating — it is a
VIEW the mesh yields (the mesh owns it, names it, builds its table from its own
faces), exactly as `volume_measure` is a view the mesh yields. The
`Frame(basis: Basis, ...)` signature is satisfied by `coarse_mesh.indicator_basis(...)`
without a new top-level type the user would have to construct by hand. The
indicator basis is mesh-OWNED, parameterised by the coarse faces, born from the
mesh — the antithesis of free-floating.

WHY the basis needs the fine measure to build (`indicator_basis(fine_measure)`):
`evaluate` tabulates `1_R` at the fine NODES, so it needs the fine cell centres
(the fine measure's nodes). This is the `Frame.table = basis.evaluate(measure.nodes)`
pattern (`frame.py:89`) — the table IS the membership matrix. Clean.

## Q5 verdict — n-D cell membership is per-axis searchsorted on coarse edges.

The prototype's 1-D `np.searchsorted(coarse_edges, fine_centers)` is the SHADOW
of the general operation: locate each fine cell centre in the coarse mesh. For a
STRUCTURED coarse mesh (the `Axis1D`-tuple `MaterialMesh`, `material_mesh.py:167`),
membership factorises per axis:

    region_index_along_axis_d = np.searchsorted(coarse.axes[d].edges, fine_centers[:, d]) - 1
    flat_region = ravel_multi_index(per-axis indices, coarse.spatial_shape)

This is n-D-clean with NO 1-D special-casing: loop over `d in range(ndim)`, one
`searchsorted` per axis against `coarse.axes[d].edges`, combine with
`np.ravel_multi_index` (the SAME `indexing='ij'` order `volume_measure` already
uses, `mesh.py:617`, `material_mesh.py:331`). The membership table the indicator
basis's `evaluate` returns is `(n_fine, n_coarse)` one-hot from `flat_region`.
The faces ARE `coarse.axes[d].edges` — the mesh already owns them. n-D for free.

This membership map is precisely a measure-theoretic operation already in the
vocabulary: it is the `pushforward` of the fine measure under the cell-locate map
π: fine_centre ↦ coarse_region (`measure.py:538`), OR equivalently the
`partition_by` of the fine measure keyed by `flat_region` (`measure.py:665`,
whose canonical use is the angular octant partition — here it is the
spatial-region partition). The region-sum numerator `M` is `partition_by(π)`
followed by per-partition `integrate`. Either primitive expresses it; both are
n-D-clean and already shipped.

## Q4 verdict — group axis rides the trailing tensor axis through ONE frame. DiscreteMeasure UNTOUCHED.

The flux `φ_{i,g}` is per-group, so the φV weight is group-dependent. DO NOT
loop ng frames, and DO NOT batch the measure weights (which would break the
hard-1-D `DiscreteMeasure.weights` constraint `measure.py:207`). The principled
factorisation: **the measure stays the GEOMETRIC volume measure `dV` (1-D
weights, group-independent); the flux enters as an INTEGRAND MULTIPLIER `M_φ`,
not as a measure.** Then:

    numerator_R,g = M[ M_φ_g(Σ_g) ]  with M = region-sum against dV
    denominator_R,g = M[ φ_g ]       (same M, integrand = φ itself)
    Σ_R,g = numerator / denominator

The group axis is a TRAILING tensor axis on the integrand. The Frame's analysis
already supports this: `analyze`'s einsum (`spherical_harmonic_basis.py:313`
"trailing axes broadcast"; `DiscreteMeasure.integrate` `measure.py:362`
`einsum("i,i...->...")`) contracts ONLY the node axis and broadcasts trailing
axes. So `M` applied to an `(n_fine, ng)` integrand returns `(n_coarse, ng)` in
ONE call — fully vectorised, ONE frame, group axis trailing. NO per-group frame,
NO batched measure.

This means the flux-weighting is `M_φ` (a `MultiplicationOperator` / coefficient
field, the `coefficient_field_promotion_frames.md` multiplier-algebra
`M:L^∞→B(L²)`), composed BEFORE the region-sum, not folded INTO the measure
weights. The "Galerkin-in-L²(φV)" reading and the "dV-measure + M_φ-multiplier"
reading are the SAME (φ·dV = M_φ then dV); the IMPLEMENTATION puts φ in the
multiplier (keeps the measure 1-D and group-independent), the MATH names it the
φV measure. The matrix-channel case (scattering `Σ_s[g_from,g_to]` weighted by
SOURCE group `g_from`) is the same with `M_φ` keyed on `g_from` — still a
trailing-axis multiplier, still one frame.

**BLAST RADIUS on the two foundational types: ZERO.**
- `DiscreteMeasure.weights` stays hard-1-D. The flux is NOT a measure weight; it
  is an integrand multiplier. Constraint UNTOUCHED. (This is the decisive payoff
  of the multiplier reading over the batched-measure reading: the foundational
  type does not churn for a single consumer.)
- `Basis` ABC: the indicator basis CONFORMS to the existing ABC unchanged
  (`evaluate`/`synthesize`/`analyze`/`reconstruct`/`mass_matrix`/`space`). The
  ABC was promoted (`base.py:13`) precisely "the Frame needs the interface, not
  a second basis" — the indicator basis IS that anticipated second basis,
  validating the un-deferral. No ABC method changes.

The only NEW code is the indicator `Basis` subclass (mesh-yielded) + the
cell-locate π (per-axis searchsorted). Both are leaves; neither touches a
foundational type.

## Augmentation answers

(a) **Energy condensation falls out as the SAME Frame shape — CONFIRMED.**
`Σ_G = Σ_{g∈G} φ_g Σ_g / Σ_{g∈G} φ_g` per material is the IDENTICAL
normal-equations projection with: trial basis K = coarse-energy-group indicator
`{1_G}` on the `EnergyGrid` (the energy analogue of `{1_R}` on the mesh); test =
trial (Galerkin in L²(φ-spectrum)); measure = the energy base-measure × spectrum
φ; the `/Σφ_g` denominator = the inverse coarse-group Gram, same `G⁻¹`. The
(K,L)-is-layered-per-axis finding holds: K is a domain object per axis (coarse
`Mesh` for space, coarse `EnergyGrid` for energy, `SphericalHarmonicBasis` for
angle), L is the per-axis base measure × the weight field. SAME `Basis` ABC,
SAME analysis=region-sum / reconstruction=G⁻¹·region-sum, SAME trailing-axis
multiplier for the weight (spectrum φ_g instead of flux φ_{i,g}). The energy
indicator basis's `evaluate(fine_group_centers)` = group-membership table
`1_G(g)`; cell-locate π becomes group-locate (which fine group ∈ which coarse
group) = `searchsorted` on coarse group boundaries. Identical structure, energy
axis. The frame architecture generalises cleanly; flag = NONE.

(b) **The asymmetry law is EXPLAINED by where K lives.** Homogenization's trial
basis `{1_R}` is defined BY the coarse mesh's faces — K is mesh-COUPLED, so its
output (coarse XS bound to coarse regions) is born on a mesh → `MaterialMesh`
(geometry+materials together). Condensation's trial basis `{1_G}` is defined by
ENERGY-group boundaries, which are MATERIAL-keyed and mesh-agnostic — K is
mesh-DECOUPLED, so its output (condensed XS per material) is a portable
`dict[int, Mixture]`. The asymmetry is precisely "is the coarse trial basis K
geometric (faces, on a mesh) or spectral (group boundaries, per material)." The
function-space framing RESPECTS and EXPLAINS the asymmetry: same Frame shape,
different K-carrier, and the K-carrier's mesh-coupling determines the output
container. This is the `feedback_domain_concrete_over_abstract` ruling in action
(K = a domain object per axis: coarse Mesh vs EnergyGrid).

(c) **The numerator IS the region-resolved ReactionRateFunctional — SHARE the
primitive.** The existing `⟨Σ,ψ⟩ = ∫ Σφ dV` (whole-domain scalar,
`functional.py:19`; `ProductionRateFunctional` `production_rate_functional.py`
contracts the GROUP axis to a per-cell density) is the WHOLE-DOMAIN reduction of
the homogenization numerator. The homogenization numerator `Σ_{i∈R} V_i φ_i Σ_i`
is that SAME `∫ Σφ dV` integrand REGION-RESOLVED instead of domain-collapsed —
i.e. the reaction-rate functional COMPOSED WITH the coarse-indicator ANALYSIS
(`M = ⟨·, 1_R⟩` resolves the single domain scalar into one-per-region). Smell
#16 shape-1 (two paths to one operator): the whole-domain reaction rate and the
region-resolved homogenization numerator are TWO reductions of ONE integrand
`Σφ·(volume measure)`; they should share the `∫ Σφ dV` primitive, differing only
in the FINAL contraction axis (all nodes → scalar, vs per-region partition →
vector). The clean factoring: `reaction_rate_density = M_Σ(φ)` (a per-cell
density, the existing functional), then the region-sum `M` against `dV` gives
the numerator, and the whole-domain rate is the same density summed over ALL
cells (the trivial one-region partition). RECOMMENDATION: the homogenization
numerator reuses `ReactionRateFunctional`'s `Σφ` density primitive, then applies
the coarse-indicator analysis `M` — do not write a second `Σφ` contraction. (The
denominator is the SAME `M` with integrand `φ` alone — `ReactionRateFunctional`
with `Σ ≡ 1`.)

(d) **Blast radius on the two foundational types: ZERO change required** (see Q4
verdict above). The user's wariness of foundational churn for a single consumer
is SATISFIED: neither `DiscreteMeasure` (1-D weights) nor `Basis` (ABC methods)
changes. The flux-weighting is a multiplier (`M_φ`), not a measure; the indicator
basis conforms to the existing ABC. The ONLY foundational-adjacent observation:
`space.py:43`'s anticipated `RegionSpace` is now FORCED (it is the indicator
basis's coefficient `space`), and `measure.py:291`'s reserved "spatial
homogenisation first Frame-binds a spatial measure" is now CONSUMED (the coarse
volume measure binds). Both were anticipated; neither is churn.

## Refuted frames (durable UNEXPLORED for this problem class)

- **Petrov-Galerkin / oblique projection (the prior memo's frame)** — REFUTED as
  the native frame: in L²(φV) the projection is ORTHOGONAL (L=K). Obliqueness is
  an artifact of the dV metric. The native frame is Galerkin-in-reweighted-metric
  → no separate oblique adjoint path needed. (This is the CORRECTION above, not a
  dead frame — PG is the SAME map in the wrong inner product.)
- **General finite-element / mixed FEM** — the coarse basis is P0 (cell
  indicators); the Gram is DIAGONAL (disjoint supports), so the normal equations
  solve in closed form (no linear solve, no mass-matrix inversion beyond a
  per-region reciprocal). A general FEM frame (continuous shape functions,
  non-diagonal mass matrix, assembly) is the WRONG weight class — it would
  introduce an assembly + solve the disjoint-support structure makes
  unnecessary. P0/finite-volume is the native basis; fire general FEM only if a
  higher-order coarse representation (P1 homogenisation) ever arrives.
- **Generalized Equivalence Theory (GET) / discontinuity factors** — REAL in the
  reactor-physics homogenization literature (Smith 1986: flux-volume weighting
  PLUS discontinuity factors to preserve surface currents). LOW-SIGNAL for THIS
  task: the brief's defining property is reaction-rate preservation ONLY
  (flux·volume average), not current/leakage preservation. GET adds a SECOND
  functional to preserve (surface current) → a SECOND projection constraint. It
  is the trigger for a FUTURE extension (a second test functional `{Ω·n on ∂R}`
  bound to the same coarse basis), not the current map. Name it UNEXPLORED with
  the precise flip condition: when leakage/current preservation joins
  reaction-rate preservation, the Frame gains a second (boundary-trace) measure.
- **Radon-Nikodym / measure absolute-continuity** — tempting via "the φV measure
  is dV reweighted by φ" (`dμ_φV/dμ_V = φ` is literally a Radon-Nikodym
  derivative). It is TRUE but adds no lever: the codebase's `pushforward`
  explicitly does NOT apply a Radon-Nikodym Jacobian (`measure.py:558`), and the
  φ-reweight is already cleanly the multiplier `M_φ`. The R-N derivative IS `M_φ`;
  naming it R-N adds vocabulary, not a test. (L-001 pattern: concrete frame —
  multiplier algebra — already captures it.)
- **Homology / chain complex** — no `∂²=0`; the fine→coarse restriction and its
  adjoint are a frame analysis/reconstruction pair, not a differential. (L-001.)
- **Category theory / Kan extension** — the fine→coarse collapse smells like a
  left Kan extension (restriction along the cell-locate functor π). Role
  fully captured by `measure.pushforward`/`partition_by` + the frame pair; no
  discriminating test the concrete frame lacks. (L-001 — name the concrete
  pushforward, list Kan extension UNEXPLORED.)

## First tests (each DISCRIMINATES — L-002)

- **DERIVED-MEASURE (Q1, the measure is φV not dV):** build the homogenization
  frame, homogenize a fine field with a SPATIALLY-VARYING flux φ_i that is NOT
  constant within a coarse region, assert `Σ_R == (Σ V_i φ_i Σ_i)/(Σ V_i φ_i)`
  bit-identical (`array_equal`). A dV-measure impl computes `(Σ V_i Σ_i)/(Σ V_i)`
  — DIFFERENT whenever φ varies within the region → REDs. (A constant-φ region
  makes the two agree → the test MUST use a varying φ to discriminate; a
  constant-φ test cannot fail and is rejected.)
- **REACTION-RATE PRESERVATION (Q2, the defining functional is preserved):**
  assert `Σ_R · Φ_R == Σ_{i∈R} V_i φ_i Σ_i` (the coarse rate equals the fine
  region rate) bit-identical, where `Φ_R = Σ_{i∈R} V_i φ_i`. This is
  `⟨Σ_c,1_R⟩_{φV} == ⟨Σ,1_R⟩_{φV}` — the projection-preserves-the-functional law.
  An impl that divides by the wrong denominator (e.g. volume `Σ V_i` instead of
  `Σ V_i φ_i`) violates the rate-preservation → REDs.
- **NUMERATOR = REACTION-RATE-FUNCTIONAL REUSE (augmentation c):** assert the
  homogenization numerator for the TRIVIAL one-region coarse mesh (the whole
  domain as a single region) equals the existing whole-domain `ReactionRateFunctional`
  value bit-identical. A separate-`Σφ`-contraction impl drifts at ULP (FP
  non-associativity, the `base.py:57` warning) → `array_equal` REDs unless the
  SAME primitive is reused. (This is the discriminator that proves the share, not
  just value-agreement.)
- **GROUP-VECTORISATION (Q4, one frame not ng frames):** homogenize an
  `(n_fine, ng)` field in ONE call, assert the result equals the per-group loop
  result bit-identical AND assert the measure's weights stayed `(n_fine,)` 1-D
  (a NEGATIVE test: constructing a `DiscreteMeasure` with `(n_fine, ng)` weights
  RAISES at `measure.py:207`). An impl that batched the measure weights would
  need the 2-D-weights measure that the constraint forbids → the negative test
  proves the multiplier reading, not the batched-measure reading.
- **n-D MEMBERSHIP (Q5, no 1-D special-casing):** build a 2-D coarse mesh,
  assert the membership table from per-axis `searchsorted` + `ravel_multi_index`
  matches a brute-force per-fine-cell point-in-region check, AND that the SAME
  code path runs for 1-D and 2-D (no `if ndim==1` branch). A 1-D-hardcoded impl
  has no 2-D path → fails to construct the 2-D frame.

## Elegance assessment

- **Structure-exposing (strong):** the `/Φ_R` denominator is EXPOSED as the
  inverse coarse Gram `G⁻¹` = the reconstruction dual factor, the exact analogue
  of the SH `(2ℓ+1) = 4π·g_C⁻¹`. Homogenization-normalisation and
  SH-reconstruction are the SAME frame operation (analysis then dual-factor),
  which the "flux·volume-weighted average" prose hides.
- **Structurally-simpler (strong):** the flux-as-multiplier reading keeps
  `DiscreteMeasure` 1-D and group-independent (ZERO foundational churn), and the
  P0 diagonal Gram makes the normal equations a per-region reciprocal (no solve).
- **Expressive (strong):** ONE Frame shape covers spatial homogenization AND
  energy condensation AND the existing SH angular projection — the (K,L)-per-axis
  layering. The numerator reuses the reaction-rate functional.
- **Algorithmic-advantage:** group axis vectorises as a trailing tensor axis
  through one frame (no ng-loop); n-D membership is per-axis searchsorted (no
  dimension special-casing).
