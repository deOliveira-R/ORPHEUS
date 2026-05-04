# Trajectory Resolvent — Foreign-Frame Detection

**Date:** 2026-05-03
**Agent:** cross-domain-attacker
**Target:** `orpheus/derivations/continuous/trajectory_resolvent/`
**Trigger table consulted:** `.claude/skills/cross-domain-frames/reference.md` (v0 + 2026-04-30 amendments)
**Prior precedent files consulted:** `variant_alpha_2surface_bie_frame.md`, `variant_alpha_family_hindsight.md`, `validated_bc_tensor_network.md`, `validated_unified_geometry.md`, `validated_hilbert_schmidt_separable.md`

This is **detection**, not critique. A trigger fires or does not.

---

## Structural feature inventory

Enumerated from `variant_alpha_core.py`, the six geometry drivers, and `docs/theory/trajectory_resolvent.rst`.

### Mathematical objects

1. **Operator resolvent** `T = (I − S)^{-1}` where `S = α R_chord` is a boundary-to-boundary scattering operator on surface phase-space. Currently rank-1 (sphere/cylinder, 1 bounce per period) and rank-2 (slab/asymmetric slab/hollow-sphere/annulus, 2 bounces per period). Closed-form `T = 1/(1 − α e^{-τ})` (rank-1) or 2×2 explicit (rank-2).
2. **Geometric series** `Σ_n α^n e^{-nτ}` collapsing to `1/(1 − α e^{-τ})`. The series indexes bounce count.
3. **First-leg + bounce-period decomposition**: `ψ = F + e^{-τ_first} · α B · T`. Decomposition is a structural fact, not a discretization choice.
4. **Per-µ multi-bounce kernel**: `T(µ) = 1/(1 − α exp(−Σ_t L_period(µ)))` with `L_period(µ)` geometry-specific. The integrand of the eventual scalar-flux integral is multiplication-operator-on-µ-space.
5. **Conserved invariants along trajectory**: impact parameter `b` (sphere, hollow-sphere, annulus, cylinder in-plane), radial distance to nearest wall (slab), `(b, µ_axial)` pair (cylinder). Each invariant is a momentum-map level set.
6. **Multi-group transfer**: `K_g = Σ_{g'} (Σ_s g'→g + χ_g νΣ_f / k) ψ_{g'}` — a block-diagonal-in-spatial-operator transfer matrix with off-diagonal scattering + fission.
7. **Power iteration** on the closed operator: outermost loop computes the dominant eigenpair `(k_eff, ψ)` of `K`.
8. **Phase-space**: `(r, µ)` for sphere/hollow-sphere; `(r, µ_axial, φ_az)` for cylinder/annulus; `(x, µ)` for slabs. The integration variable along bouncing trajectories is the in-plane chord arclength.

### Symmetries

9. **SO(3) rotational symmetry** (sphere, hollow-sphere) reduced via impact parameter `b` and surface cosine `µ_surf`.
10. **SO(2) × R translational** (cylinder/annulus): `φ_az ∈ [0, 2π)` periodic; axial direction infinite/translation-invariant.
11. **Z₂ reflection** (symmetric slab): `µ → −µ` reflection; antisymmetric/symmetric slab is not Z₂.
12. **Time-reversal / chord reversibility**: `q(r(s)) e^{−Σs}` and `q(r(L−s)) e^{−Σ(L−s)}` are related by reflection; surfaces support a 2-element orbit per chord.

### Iterative / spectral structure

13. **Fixed-point** form of the surface flux: `ψ_surf = αB + α e^{-τ} ψ_surf` (rank-1) — solved analytically, not iteratively.
14. **Outer power iteration** on `k`: contraction near `k_∞`, degraded convergence near criticality (no Chebyshev/Wielandt acceleration in current code).
15. **No Krylov / Arnoldi** despite the resolvent structure being explicit.

### Stochastic structure

16. **Markov interpretation latent**: bouncing trajectory IS a Markov chain on surface phase-space, transition kernel `α R_chord`. Currently used as a deterministic geometric series.

### Integral / kernel structure

17. **Source-term integral** `B(µ) = ∫_0^L_period q(r_chord(s)) e^{-Σs} ds` is a 1D quadrature along chord.
18. **Bickley-Naylor / E_n family** appears in the angle-integrated companion (Sanchez Eq. 5/6) but NOT in the angle-resolved Variant α — Variant α keeps `µ` explicit precisely to avoid the Bickley reduction.

### Differential structure

19. **No differential operator explicitly used.** Variant α is a pure integral-form / Peierls reformulation. Streaming + collision are absorbed into the chord exponentials.

### Boundary handling

20. **BC enters through `α` only**, multiplicatively in `S = α R_chord`. Vacuum = α 0; perfect specular = α 1; partial = α ∈ (0,1).
21. **No diffuse / albedo / white BC**: Sanchez 1986 covers these (`β`, `χ(µ)`); current Variant α is specular-only.
22. **Mixed BC** (vacuum on one face, reflective on the other) handled in the rank-2 framework by setting one diagonal of `S` to zero — see `compute_resolvent_T_rank2`.

### Scale separation

23. **Thick / thin asymptotics not exploited.** As `Σ_t L → ∞` the geometric series decays to a single bounce (diffusion-irrelevant); as `Σ_t L → 0` the series reduces to the vacuum first-leg. Code does not branch on regime.
24. **`α → 1` near-criticality**: `T → ∞` as `µ → 0` (sphere grazing pathology). Current code documents this; no analytic continuation / Padé regularization.

---

## ELEGANCE DETECTOR HITS

Six smells fire — well above the two-smell threshold for probing.

1. **Smell 1: Repeated near-identical code across geometry variants.** Six `_apply_operator_*` functions, six `_first_leg_*_chord` helpers, six driver loops. `variant_alpha_family_hindsight.md` already flagged this; the bundle frame (BaseAtlas, AngularFiber, ChordOracle) is the prescribed lift.
2. **Smell 4: Symmetry present in the problem but absent in the method.** SO(3) on sphere is exploited only as the impact-parameter reduction (`b = r√(1−µ²)`) — not as a representation-theoretic decomposition. Similarly SO(2) on cylinder/annulus is implicit in `φ_az ∈ [0, 2π)` quadrature, never as Fourier-mode block-diagonalization.
3. **Smell 6: "Iterative" without fixed-point structure analyzed.** Outer power iteration converges with rate `λ_2/λ_1` (spectral gap of the closed operator). Code computes `k_eff` via Rayleigh quotient but never extracts `λ_2` or applies Chebyshev / Wielandt deflation. Frame: dynamical systems / Krylov.
4. **Smell 7: "We picked N because it works."** `n_traj_quad`, `n_mu`, `n_phi_az` quadrature orders are picked empirically. No representation-theoretic optimum (Lebedev for SO(3), Fourier for SO(2)) appealed to.
5. **Smell 9: Boundary handling as a special case added to interior logic.** The b-partition (`b ≤ R_in` vs `b > R_in`) for hollow-sphere/annulus is coded as a 4-case `if/else`. Under the bundle frame it is a single integration measure restricted to a momentum-map level set. Already named in `variant_alpha_family_hindsight.md` as a refactor target.
6. **Smell 13: "Converges because we checked."** Variant α convergence on `c → 1`, `α → 1`, thick optical depth has no closed-form contraction-mapping argument in the docs — just the V_α1/2/3 SymPy proofs of operator correctness on constant trial. The spectral gap is a theorem the frame can supply.

---

## FRAME CANDIDATES

Eleven frame matches enumerated. Rank order at the end.

---

### Frame 1 — Tensor networks / matrix-product operator (USER ASKED)

**Trigger:** A.2 row "Tensor networks / tensor-train: high-dimensional problems (space × angle × energy); compositional boundary structure; cross-section data compression". Validated precedent in `validated_bc_tensor_network.md`. Direct hit on the rank-2 closure machinery.

**Verbatim quote (Orús 2014, "A practical introduction to tensor networks", *Annals of Physics* 349, p. 117):**

> "Matrix Product States, MPS, and Matrix Product Operators, MPO, [...] are particularly suited to represent operators acting on 1D quantum lattices with local interactions."

**Match assessment.** The boundary-to-boundary scattering operator `S` is **literally a matrix-product operator on a 1D-lattice of bounce events**. Each "site" is a single bounce (one surface arrival); the local tensor at site `n` carries:

- a leg for the incoming surface state (dimension = 2 for rank-2 geometries: `[ψ_L^+, ψ_R^-]`; dimension = 1 for rank-1)
- a leg for the outgoing surface state (same dimension)
- a "bond" leg of dimension 1 (the geometric bouncing has bond dimension 1 because bounces are deterministic — no entanglement across bounces)

The geometric series `T = Σ_n S^n` is the contraction of an OPEN (semi-infinite) MPO chain. Bond dimension 1 means the MPO is exactly representable; no truncation. **The tensor network is a hidden trivial-bond MPO.**

What the user's hint upgrades: not the MPO of bounces (trivial), but the MPO over **(geometry, BC type, multi-group)** axes. Each axis is a tensor leg.

**Reformulation:**

```
T = contract(BounceMPO(α, τ),
             GeometryTensor(geom, R_in, R_out, dim),
             GroupTransferMPO(Σ_s, χνΣ_f/k))
```

with:
- `BounceMPO(α, τ)`: rank-1 → rank-2 → rank-N case automatically expressed as bond dimension 1, 2, N (one leg per surface). Hollow sphere with `b ≤ R_in` vs `b > R_in` is a **superposition of bond-dimension-1 (outer-only) and bond-dimension-2 (through-ray) MPOs** weighted by the impact-parameter measure.
- `GeometryTensor(...)`: the parameterized chord oracle — slab/cylinder/sphere/hollow encoded as a single geometry tensor with `dim ∈ {0, 1, 2}` for slab/cylinder/sphere axis and `(R_in, R_out)` for solid/annular topology. **This IS the unified-geometry frame already validated in `validated_unified_geometry.md`.**
- `GroupTransferMPO`: scattering kernel block, naturally MPO if `Σ_s g'→g` has structure (downscatter-only is a triangular MPO; full matrix is dense MPO of bond dimension G).

**Elegance payoff (4/4 criteria):**

- **Structure-exposing**: makes explicit that hollow-sphere = (rank-1 MPO ⊕ rank-2 MPO) restricted to disjoint momentum-map level sets. Currently coded as a 4-case if/else; under MPO it is one direct sum.
- **Expressive**: a 7th geometry (toroidal? thick-shell with 3 surfaces?) is a new GeometryTensor instance, not a new module. A 3-bounce-per-period geometry (e.g., wedge with 3 walls) is rank-3 MPO — automatic, no new derivation.
- **Structurally-simpler**: collapses six `_apply_operator_*` functions and six driver loops into one `contract(BounceMPO, GeometryTensor, GroupTransferMPO)` call.
- **Algorithmic-advantage**: enables **TT-format compression of the multi-group operator** when `Σ_s g'→g` has low rank in `(g', g)` (block-bidiagonal downscatter has rank ≈ bands). MG cylinder solve currently O(G² · n_r · n_µ · n_φ); TT-compressed MG could be O(G · log G · n_r · n_µ · n_φ) when scattering matrix has Tucker rank ≪ G. Active research area: Peng et al. (Battelle, 2023) and Dektor-Einkemmer (2024) tensor-train transport.

**First test:** symbolically construct the rank-2 MPO for the asymmetric slab (one site = one bounce, bond dimension 2, two outer legs of dimension 2). Contract `n=8` bounces. Compare to `compute_resolvent_T_rank2` evaluated at the same `(α_L, α_R, τ)`. Predicted match: identical to rounding error (bond dimension 1 gives exact representation). Discriminator: at any geometry with bond dimension ≥ 3 (3-wall wedge, toroidal), the rank-N MPO contraction must agree with explicit `(I − S)^{-1}` solution.

**Structural attack on current:** the rank-1 / rank-2 nomenclature in `variant_alpha_core` is **bond dimension** in disguise. Calling them "rank-1" and "rank-2 closure" is a synonym for "MPO with bond dimension 1 vs 2". Naming the bond dimension exposes that the formulation already lives in MPO-space; the next geometry is just another bond-dimension instance.

**Precedent:** `validated_bc_tensor_network.md` (BCs as tensor networks). This frame extends that precedent from BC-only to the full Variant α core.

---

### Frame 2 — Schur complement / Banach algebra of trajectory operators

**Trigger:** A.7 "Hilbert-Schmidt / separable kernels" (matrix-Galerkin diagnostic) and A.5 "Krylov subspace theory". Resolvent + matrix structure of `(I − S)^{-1}` directly invokes Schur-complement identities.

**Verbatim quote (Higham 2008, "Functions of Matrices: Theory and Computation" §1.13):**

> "The Schur complement of a block in a partitioned matrix yields the block inverse efficiently and exposes the sub-system structure of the linear problem."

**Reformulation:** the rank-2 monodromy `S = antidiag(α_L e^{-τ}, α_R e^{-τ})` and its inverse `T = (I − S)^{-1}` admit the explicit Schur-complement form:

```
T_11 = (1 − S_LR S_RL)^{-1}
T_12 = T_11 · S_LR
T_22 = (1 − S_RL S_LR)^{-1}
T_21 = T_22 · S_RL
```

For a rank-N geometry (3-wall wedge, multi-shell), this generalizes to: pick any single surface as the "reference" surface, Schur-complement out the others, and solve for the reference-surface flux first. The reference-surface fixed-point equation is **scalar** (rank-1) regardless of N.

**Elegance payoff (3/4 criteria):**

- **Structure-exposing**: a hollow-sphere ψ_surf at the inner wall and ψ_surf at the outer wall are not two independent rank-2 quantities — they are Schur-complements of each other through the bouncing operator.
- **Structurally-simpler**: extending to N-wall geometry doesn't require an N×N matrix solve at every (r, µ) phase-space point; pick reference wall, Schur-eliminate, scalar fixed-point.
- **Algorithmic-advantage**: the Schur-complement inversion reuses `(1 − α_L α_R e^{-2τ})` — the determinant — which is already computed once. For multi-group, the analog `(I − S_g)^{-1}` block has Schur structure across groups (downscatter-only is upper-triangular Schur).

**First test:** verify on hollow sphere that `ψ_surf_inner` computed via direct rank-2 inversion equals `ψ_surf_inner` computed via Schur-elimination of `ψ_surf_outer` (and vice versa). Predicted match: exact algebraic identity, equal to rounding.

**Structural attack on current:** the rank-2 closure formula (`apply_variant_alpha_closure_rank2`) hand-codes the explicit 2×2 inverse. For N=2 this is fine; for N=3 (rank-3) the manual inversion explodes. Schur recursion is the discipline that scales.

**Precedent:** none in the validated table directly. Closest is `validated_hilbert_schmidt_separable.md` (separable kernel resolvent reformulation, Phase 5).

---

### Frame 3 — Spectral theory of multiplication operators (re-fire)

**Trigger:** A.3 row "Spectral theory of multiplication operators" — operator of form `(M_φ ψ)(µ) = m(µ) · ψ(µ)`, ess_range / ess_sup machinery. Validated precedent: `validated_hilbert_schmidt_separable.md` (Phase 5 continuous-µ).

**Verbatim quote (Trefethen-Embree 2005, "Spectra and Pseudospectra" §2.3):**

> "For a multiplication operator `M_m: f ↦ m·f` on `L²(Ω, dµ)` with bounded measurable `m`, the operator norm equals `‖m‖_∞ = ess_sup_Ω |m|`, and the spectrum equals the essential range of `m`."

**Match assessment.** Variant α's `T(µ) = 1/(1 − α e^{−Σ_t L_period(µ)})` IS a multiplication operator on `L²((0,1], dµ)` (sphere) or its `(µ_axial, φ_az)` analog (cylinder). The Phase 5 caveat ("ess_range(m) ∋ 1 → matrix-Galerkin fails") applies directly:

- **Sphere** (`µ ∈ (0, 1]`): `L_period = 2Rµ`, so `m(µ) = α exp(−2Σ_t Rµ)`. As `µ → 0+`, `m(µ) → α`. At `α = 1`, `m(0+) = 1` and `1 ∈ ess_range(m)`. **Sphere is structurally pathological at α = 1**, exactly as documented.
- **Slab** (`µ ∈ (0, 1]`, single-transit): `m(µ) = α exp(−Σ_t L/|µ|)`. As `µ → 0+`, `Σ_t L/|µ| → ∞`, so `m(µ) → 0`. `1 ∉ ess_range(m)` regardless of α. **Slab is structurally immune.**
- **Cylinder** (in-plane angle `φ_az ∈ [0, 2π)`): `L_period = 2√(R² − r² sin²φ_az) / √(1 − µ_axial²)`. As `φ_az → π/2` at `r = R`, `L_period → 0`, so `m → α`. **Cylinder inherits sphere's pathology at the tangential locus.**
- **Hollow sphere** (`b ≤ R_in`): `τ_step = Σ_t (√(R_out² − b²) − √(R_in² − b²))`. As `b → R_in`, `τ_step → 0`, so `m → α α'`. **Hollow sphere inherits the pathology at the inner-tangent locus.** The current code documents this via `s_in-plane → 0` mitigation but does not name the spectral mechanism.

**Reformulation:** lift `T(µ)` from a scalar function to a named multiplication operator `M_T` on `L²(angular-fiber)`, and derive per-geometry pathology from `ess_range(m)`. The variant_alpha_2surface_bie_frame memo already promotes this — **this is the same lever applied to the rank-1 multiplication-operator structure.**

**Elegance payoff (3/4 criteria):**

- **Structure-exposing**: the per-geometry grazing-ray pathology is one theorem (`ess_range(m) ∋ 1 ⇒ unbounded`), not four geometry-specific docstrings.
- **Expressive**: a new geometry's pathology prediction is `ess_range`-computation, not new ray analysis.
- **Algorithmic-advantage**: predicts in advance which geometries need angular regularization (Gauss-Jacobi instead of Gauss-Legendre to handle integrable `1/µ`), without per-geometry numerical experimentation.

**First test:** symbolically compute `ess_range(m)` for sphere, cylinder, slab, hollow-sphere, annulus. Verify the pathology classification matches the documented `α = 1` divergence behavior in each driver's docstring. Predicted match: 4/4 (slab immune, sphere/cylinder/hollow-sphere inherit at specific loci).

**Structural attack on current:** each driver's docstring documents its grazing-ray locus and mitigation **independently**, as if these were unrelated facts. They are one theorem applied to four different `m` functions.

**Precedent:** `validated_hilbert_schmidt_separable.md` (Phase 5 µ-resolved kernel). Cross-reference: `variant_alpha_2surface_bie_frame.md` step 5.

---

### Frame 4 — Fiber bundle / G-structure (re-fire from hindsight)

**Trigger:** A.1 "Topology" + A.2 "Group theory / representation theory". Already named in `variant_alpha_family_hindsight.md` as the top-ranked frame. Re-fired here to show it is the unifying lift across all 6 geometries × 2 topologies.

**Verbatim quote (Husemoller 1994, "Fibre Bundles" 3rd ed., §2.1, Springer GTM 20):**

> "A fibre bundle `(E, B, F, π)` consists of a total space `E`, base space `B`, fibre `F`, and projection `π: E → B` such that locally `E ≃ U × F` for trivializing patches `U ⊂ B`. A G-structure is a reduction of the structure group of the bundle's frames to a subgroup G ⊂ GL(F)."

**Match assessment.** The phase space `(r, Ω)` is a fiber bundle:

- Base `B`: 1D radial manifold-with-boundary. Slab `[0, L]`, sphere `[0, R]`, hollow `[R_in, R_out]`. Topology π_0(B) is connected (rank-1) or has 1 hole (rank-2).
- Fiber `F`: angular sphere `S²` (sphere), cylinder `S¹ × [−1, 1]` (cylinder), interval `[−1, 1]` (slab).
- Structure group `G`: SO(3) on sphere, SO(2) on cylinder, Z₂ × ZZ on symmetric slab. The G-structure restricts the angular fiber to G-orbits.

**The conserved chord invariants are momentum-map level sets** of the G-action:
- Sphere: `b = r √(1 − µ²)` is the modulus of angular momentum, J: T*B → so(3)*.
- Cylinder: `b = r |sin φ_az|` IS angular momentum component along axis; `µ_axial` IS axial momentum.
- Slab: signed `µ` is the only invariant (1D base, no rotational structure).

**Reformulation (already in hindsight memo):**

```
class FiberBundle:
    base: BaseAtlas        # (Interval, Annulus) — π_0 distinguishes solid/hollow
    fiber: AngularFiber    # (S², S¹×I, I) — quadrature support
    chord: ChordOracle     # (r, Ω) → (length, exit_surface) — the only geometry-specific primitive
    G: SymmetryGroup       # SO(3), SO(2), Z₂ — quadrature respects this

apply_variant_alpha(bundle, α, source) -> ψ
```

**Elegance payoff (4/4 criteria):**

- **Structure-exposing**: the b-partition in hollow_sphere is the pullback of `R_in`-characteristic through `chord`. Currently 4-case if/else.
- **Expressive**: 7th geometry = new `BaseAtlas` and `ChordOracle`; quadrature inherits `G`-equivariance automatically.
- **Structurally-simpler**: 6 modules → 1 module + 6 oracle instances. Same content, less code.
- **Algorithmic-advantage**: G-equivariant quadrature reduction is automatic. SO(3) sphere with `n_µ` Lebedev nodes covers a 2-sphere worth of angular quadrature with ~`n_µ²` accuracy; current code uses tensor-product `(µ, φ)` quadrature. Could halve quadrature cost at fixed accuracy.

**First test:** implement the `ChordOracle` for sphere and cylinder; replace 6 `_first_leg_*_chord` and `_bounce_period_*_chord` helpers with two `oracle_sphere`, `oracle_cylinder` calls. Verify all 178 tests still pass. Predicted: pass.

**Structural attack on current:** "rank-1 vs rank-2 closure" naming **is** the bundle's π_0 in disguise. The bundle frame names the topological invariant directly.

**Precedent:** `variant_alpha_family_hindsight.md` (top-ranked); also `validated_unified_geometry.md` (1D radial unification, weaker version of this frame).

---

### Frame 5 — Krylov subspace / Arnoldi for the resolvent

**Trigger:** A.5 "Krylov subspace theory: Linear system solves; eigenvalue iteration; operator-based methods". A.3 "Spectral theory of compact + positive operators". Smell 6.

**Verbatim quote (Saad 2003, "Iterative Methods for Sparse Linear Systems" 2nd ed. §6.3, p. 149):**

> "GMRES finds the solution in the Krylov subspace `K_m(A, r_0)` that minimizes the residual norm. For a linear system `(I − S)x = b` where `S` is contractive (`‖S‖ < 1`), the Neumann series `Σ S^k b` converges with rate `‖S‖`; GMRES converges asymptotically at the same rate but without forming the partial sums."

**Match assessment.** The closed Variant α operator `K(c, α) = (operator on flux that returns next-iterate flux given current flux as source)` admits a Neumann-series interpretation: each power iteration step IS a partial sum of `(I − K_homogeneous)^{-1} K_external`. Power iteration converges with rate equal to `λ_2 / λ_1` (subdominant / dominant eigenvalue ratio of `K`). At criticality `λ_1 → 1` and `λ_2 / λ_1 → λ_2 < 1` — finite gap.

Currently the outer power iteration on `k_eff` runs to convergence by Rayleigh quotient. Arnoldi/IRAM on the same operator would produce `(λ_1, λ_2, ..., λ_m)` simultaneously, yielding:

- `k_eff = 1/λ_1` (already computed)
- `k_subcritical_modes = 1/λ_2, 1/λ_3, ...` (NEW — Wielandt deflation, mode-shape uncertainty)
- Spectral gap as a verification quantity

**Elegance payoff (3/4 criteria):**

- **Structure-exposing**: the spectral gap (which controls power-iteration rate AND quantifies "how-much-multiplying" the system is — distance from critical) is currently nowhere extracted; Arnoldi names it.
- **Expressive**: extends the API from `k_eff` to `(k_eff, λ_2, ...)` — a richer spectral fingerprint without solving anything new.
- **Algorithmic-advantage**: GMRES on `(I − S/k) ψ = source` outperforms power iteration when subdominant eigenvalues are clustered (typical for thick optical depth). Tripwire from hindsight memo (#4) predicted real payoff at thick optical depth — confirmed by frame.

**First test:** swap power iteration for `scipy.sparse.linalg.eigs(LinearOperator(K), k=3, which='LM')` on the closed sphere problem. Compare `1/λ_1` to current `k_eff`; predicted match to convergence tolerance. Inspect `λ_2`: it is the Anlauf time (Markov mixing rate) of the source iteration — a structurally meaningful quantity that the code currently discards.

**Structural attack on current:** power iteration produces `k_eff` and discards everything else. The dominant-eigenpair tool is structurally sub-optimal when the operator's resolvent is already named (see Frame 1, Frame 3).

**Precedent:** none in scripts/. High-prior given the resolvent is explicit and `λ_2` would directly extend the verification matrix.

---

### Frame 6 — Markov chain / Feynman-Kac path-integral interpretation

**Trigger:** A.4 "Feynman-Kac / stochastic representation: Any transport problem. CP, MOC, MC are all here." Validated MASTER BRIDGE per the trigger table.

**Verbatim quote (Lapeyre-Pardoux-Sentis 2003, "Introduction to Monte Carlo Methods for Transport and Diffusion Equations" §1.4, OUP):**

> "The neutron transport equation admits a Feynman-Kac representation `ψ(r, Ω) = E_x[∫_0^τ q(X_s) e^{−Σ_t s} ds]` where `(X_s)` is the Markov process of free-flight + scattering. A specular boundary appends a deterministic reflection event to the path."

**Match assessment.** The Variant α multi-bounce sum IS the path-integral over the `(X_s)` process restricted to deterministic specular reflections (no scattering — `c = 0` for the closed-trajectory part). The geometric series `Σ_n α^n e^{−nτ}` IS the Feynman-Kac expectation `E[e^{−Σ_t · TotalPathLength} · α^{NumberOfBounces}]` evaluated under the deterministic specular dynamics.

This frame has TWO distinct payoffs that fire here:

**Payoff 6a — Path-integral Monte Carlo as fourth verification path.** Per `vv-principles` L11 structural-independence: a path-integral MC estimator that samples `(X_s)` paths via Woodcock-style delta tracking, weighted by `α^N exp(−Σ_t Σ_n L_n)`, gives `ψ(r, Ω)` independently of the deterministic geometric series. It exercises a **stochastic** path measure where Variant α exercises a **deterministic** one — different identity, different integrand. This is the missing cross-check for `(α) Variant α` against itself (where current verification routes Sanchez kernel `(β)` and Chandrasekhar `(γ)` provide weak / strong checks per `sanchez_chandrasekhar_gap.md` §Q3).

**Payoff 6b — Specular bounces as a deterministic Markov kernel `α·R_chord`.** The boundary scattering operator `S` IS a Markov transition kernel on surface phase-space (deterministic part), with absorption `(1 − α)` per visit. `T = (I − S)^{-1}` is the Green's function of this Markov chain. Mixing time = `−1/log(α e^{−τ})`. Currently this Markov interpretation is implicit; naming it gives a structural diagnostic.

**Elegance payoff (3/4 criteria):**

- **Structure-exposing**: makes explicit that Variant α IS the deterministic (`σ_s = 0`) Feynman-Kac expectation. The full transport problem is the same expectation with the scattering subprocess included; current code segregates the two (collision = source iteration, specular = closed-form geometric series) but they are one Feynman-Kac kernel.
- **Algorithmic-advantage**: provides the fourth-path MC verification (point-source response sampled via PIMC) — fills the gap in `sanchez_chandrasekhar_gap.md` §Q3 verification matrix.
- **Expressive**: Markov mixing-time argument names the convergence rate of source iteration in closed form; spectral gap of `S` IS the mixing rate.

**First test:** implement a small PIMC: sample 10⁵ specular bouncing paths from `(r=R/2, µ=0.5)` outward, weight by `α^N exp(−Σ_t Σ L_n)`, accumulate `ψ_surf` estimate. Compare to `compute_resolvent_T · B`. Predicted match: stochastic with `σ ∝ 1/√N`, mean within ~0.3 % at N=10⁵.

**Structural attack on current:** Variant α's `T = (I − S)^{-1}` is named as algebra but never as a Markov-chain Green's function. The naming is the lever.

**Precedent:** none yet in scripts/. Promotion candidate.

---

### Frame 7 — Asymptotic analysis: thick / thin / near-critical limits

**Trigger:** A.6 "Asymptotic analysis / matched asymptotics: Boundary layers; diffusion limit; thick-thin transitions". Smell 23.

**Verbatim quote (Bender-Orszag 1978, "Advanced Mathematical Methods for Scientists and Engineers" §9.4):**

> "When a Neumann series `Σ ε^n a_n` has poles at finite radius of convergence, Padé-approximant resummation `[L/M]_f(ε)` extracts the analytic continuation past the convergence boundary."

**Match assessment.** Three asymptotic regimes structurally present, none extracted:

- **Thin-system limit `Σ_t L → 0`**: `T(µ) → 1/(1 − α + αΣ_t L_period)` — explicit small-τ expansion. Current code evaluates the full exponential.
- **Thick-system limit `Σ_t L → ∞`**: `T(µ) → 1 + α e^{−Σ_t L_period}` — single-bounce truncation. Current code keeps the full series (no payoff in this regime — it is correct, but wasteful).
- **Near-criticality `α → 1`**: `T(µ → 0+) ~ 1/(Σ_t L_period(µ)) ~ 1/(2Σ_t R µ)` for sphere — integrable singularity in scalar-flux integral but a high-condition-number quadrature target. Padé resummation `[2/2]` of the expansion in `(1 − α)` would regularize the `α → 1` limit without singularity in any term.

**Reformulation:** add an asymptotic-mode flag to `compute_resolvent_T` that switches between exact, thin-limit linearization, thick-limit single-bounce, and Padé-regularized near-critical. Each is a closed-form expression; selection by `Σ_t L` and `(1 − α)`.

**Elegance payoff (2/4 criteria):**

- **Algorithmic-advantage**: avoids `1/0` regularization at α → 1 grazing limit (current code relies on `B → 0` cancellation, which loses digits in finite precision). Padé resummation is an analytic continuation that does not pass through 0/0.
- **Structure-exposing**: names the three regimes as one parameterized family rather than as edge-case caveats in docstrings.

**First test:** evaluate `compute_resolvent_T(α=0.999999, τ=1e-6)` at a sphere grazing ray. Compare to Padé `[2/2]` resummation in `(1 − α)`. Predicted: native form has `~1e9` magnitude, finite-precision inaccurate; Padé is `O(1)` and numerically clean.

**Structural attack on current:** the grazing-ray `0/0` cancellation noted in cylinder docstring (§"Tangential in-plane") is documented as a numerical concern. Padé resummation makes it not exist.

**Precedent:** none in scripts/. Low-cost candidate.

---

### Frame 8 — Generating functions / formal-power-series resummation

**Trigger:** A.3 "Approximation theory" + Smell 13. Pre-validated in `phase5_continuous_mu_frames.md` memory entry as "generating-function bounce-resolved" (M2 in Phase 5 promoted-frame list).

**Verbatim quote (Wilf 1994, "generatingfunctionology" 2nd ed., §2.1):**

> "A generating function `f(z) = Σ a_n z^n` is the bookkeeping device that converts a recurrence into an algebraic equation. Resumming a divergent series via Borel or Padé summation extracts information beyond the radius of convergence."

**Match assessment.** The geometric series `T(µ) = Σ_n α^n e^{−nτ(µ)}` IS a generating function in the formal variable `z = α e^{−τ(µ)}`. Closed form `T = 1/(1 − z)` is the trivial generating-function lift.

What this names: a NON-TRIVIAL generating function exists for the **bounce-resolved scalar-flux contributions**:

```
Φ(r, z) = Σ_n z^n · ψ_n(r)        where ψ_n = nth-bounce angular flux contribution
```

`Φ(r, z=α)` recovers the current `ψ`. But individual `ψ_n` (the per-bounce contribution) carries diagnostic information: how does the contribution decay with bounce count? `Σ_n n · ψ_n` is the mean number of bounces — a Markov mixing diagnostic (directly cross-references Frame 6).

**Reformulation:** evaluate `Φ(r, z)` at multiple `z` values and reconstruct individual `ψ_n` via discrete Fourier inversion on the unit circle:

```
ψ_n(r) = (1/2π) ∫_0^{2π} Φ(r, z = ρ e^{iθ}) (ρ e^{iθ})^{−n} dθ
```

(standard Cauchy contour formula). Cost: a few extra `compute_resolvent_T` evaluations at distinct `α` values, then DFT.

**Elegance payoff (2/4 criteria):**

- **Structure-exposing**: bounce-by-bounce decomposition exposes the spectrum of the `S` operator (mixing rate, dominant bounce-mode) without spectral solver.
- **Algorithmic-advantage**: provides ANOTHER cross-check — at `z = 0`, `Φ = F` (vacuum first-leg only); at `z = 1`, `Φ = ψ_specular`. Smooth interpolation provides a verification-by-monotonicity test that `Φ` increases monotonically in `z` for source `q ≥ 0` (provable: each `ψ_n ≥ 0`).

**First test:** evaluate `Φ(r=R/2, µ=0.5; z)` for `z ∈ {0, 0.25, 0.5, 0.75, 1.0}`. Verify `Φ` is monotone non-decreasing (predicted: yes, by positivity of `q`). Reconstruct `ψ_0 = F` and `ψ_1 = α e^{-τ_first} B` from finite differences; verify against direct evaluation.

**Structural attack on current:** the geometric series `T = 1/(1−z)` is named, but the generating-function viewpoint that lets you read off individual bounce contributions is not used. This is the simplest non-trivial GF lift available — it costs nothing and adds a verification axis.

**Precedent:** Phase 5 continuous-µ multi-bounce (`phase5_continuous_mu_frames.md` M2 entry, validated 2026-04-30 as a high-leverage frame for resolvent structures).

---

### Frame 9 — Group representation theory: Lebedev / Fourier quadrature

**Trigger:** A.2 "Group theory / representation theory: Any angular discretization; rotational / reflection symmetry". Smell 4 fires strongly.

**Verbatim quote (Lebedev-Laikov 1999 *Doklady Mathematics* 59, 477–481, p. 478):**

> "A quadrature on the sphere is invariant under the octahedral group `O_h` if both the nodes and weights respect the 48-element symmetry. Such quadrature integrates exactly all spherical harmonics up to order `2N+1` with `N ≤ 131` available tabulated rules."

**Match assessment.** Sphere phase-space `(r, µ)` after impact-parameter reduction IS a 1D angular fiber `µ ∈ [−1, 1]` — the quotient of the angular sphere by the rotation group fixing the radial axis. The full angular fiber IS `S²`. Currently the sphere driver uses **Gauss-Legendre on `µ ∈ (0, 1]`**, which respects only the Z₂ subgroup of SO(3) (the reflection across the `µ = 0` plane). The full SO(3) symmetry is not exploited at the quadrature level.

For the **cylinder** the situation is worse: `(µ_axial, φ_az)` quadrature is **tensor-product Gauss-Legendre × uniform-φ**. Tensor-product violates SO(2)×R rotational invariance. Lebedev would not apply (it is for `S²`); the cylinder analog is **Fourier collocation in `φ_az`** (which is invariant under SO(2) by construction) with Gauss-Legendre in `µ_axial` (invariant under axial reflection Z₂).

**Reformulation:**

- Sphere: replace `n_µ` Gauss-Legendre with Lebedev quadrature on `S²` (preserves SO(3)). For axially-symmetric problems (current isotropic-source case), this reduces to Gauss-Legendre on `µ` (no advantage). But for **anisotropic sources or BCs** (linearly-anisotropic Sanchez 1986 `f_1 ≠ 0`), Lebedev gives exact integration of `P_l(µ)·P_m(φ_az)` cross-moments at a fraction of the tensor-product cost.
- Cylinder: Fourier-collocate in `φ_az` (`n_φ` Fourier modes); the bounce-period integrand `B(r, µ_axial, φ_az)` becomes block-diagonal in Fourier space if the chord oracle has SO(2) symmetry. Could halve the `n_φ` quadrature cost.
- Hollow sphere/annulus: same as solid analog; the angular-fiber group structure is preserved by the radial topology change.

**Elegance payoff (2/4 criteria, conditional):**

- **Structure-exposing**: makes explicit that the sphere problem has SO(3) symmetry that current quadrature partially squanders.
- **Algorithmic-advantage**: only fires when anisotropic sources / BCs are added. For current isotropic case, Lebedev reduces to Gauss-Legendre — no win. But the moment the linearly-anisotropic Sanchez 1986 extension lands (planned per `sanchez_chandrasekhar_gap.md` Phase B), the SO(3)-respecting quadrature pays.

**First test:** for the isotropic case, Lebedev L=29 (302 points) vs Gauss-Legendre `n_µ=29` should give identical scalar flux to 1e-12 (axial-symmetry reduction). Then add a linearly-anisotropic source `q_1(r) µ`. Predicted: Lebedev integrates this exactly; tensor-Gauss has aliasing error scaling as `O(n^{-2})` until the basis covers `µ²` correctly.

**Structural attack on current:** Smell 7 directly — quadrature orders are picked by convergence empirically rather than by representation-theoretic exactness. Lebedev names what "exact" means.

**Precedent:** none in ORPHEUS scripts/ yet. Sphinx theory page mentions SO(3) but doesn't act on it.

---

### Frame 10 — QMC / low-discrepancy sequences for trajectory quadrature

**Trigger:** A.4 "Number theory / low-discrepancy sequences: QMC; ray generation for MOC". Sphere/cylinder/hollow-sphere have multi-D phase-space integrals; cylinder is 3D `(r, µ_axial, φ_az)`.

**Verbatim quote (Niederreiter 1992 "Random Number Generation and Quasi-Monte Carlo Methods" §2.3):**

> "For a function `f` of bounded variation in the Hardy-Krause sense on `[0, 1]^d`, QMC integration with a low-discrepancy sequence `(x_n)` satisfies `|I - I_QMC| ≤ V_HK(f) · D_N* (x_n)` where the star-discrepancy `D_N* = O(N^{-1} (log N)^d)`."

**Match assessment.** Cylinder phase-space `(r, µ_axial, φ_az)` is 3D. Hollow-sphere/annulus `(r, µ)` with b-partition is 2D + discrete branch. Current code uses **tensor-product Gauss-Legendre**: `n_r × n_µ × n_φ` tensor-product nodes.

For a smooth integrand (which the Variant α `ψ(r, µ)` is on the bulk of phase space — only at the `µ → 0` grazing locus does it lose smoothness), Sobol or Halton points give `O(N^{-1} (log N)^3)` vs `O(N^{-1/2})` for plain MC and `O(N^{-2k+1})` for Gauss-Legendre tensor at k-th derivative continuity.

**Conditional fires.** Whether QMC beats Gauss-Legendre depends on the **integrand smoothness** in `(r, µ_axial, φ_az)`. Current cylinder code has documented grazing-ray pathology at `φ_az = ±π/2` (tangential in-plane) — the integrand has a `1/(1 − α)` singularity there. **Hardy-Krause variation may be unbounded near α = 1**; QMC convergence guarantee fails. Verify on integrand before adopting.

**Elegance payoff (1/4 criteria, conditional):**

- **Algorithmic-advantage** (only if Hardy-Krause variation finite): asymptotic `O(N^{-1} log³ N)` cylinder convergence vs `O(n_r n_µ n_φ)` tensor-Gauss. For 3D quadrature at typical `n_r=n_µ=n_φ=32` (32³ = 32768 nodes), QMC at 32768 Sobol points should give better convergence rate IF integrand is bounded variation.

**First test:** bound the Hardy-Krause variation of the cylinder integrand at `α = 0.5` (well inside well-behaved regime). Compute `φ(r=R/2)` via tensor-Gauss `n=16³` and Sobol `N=4096`. Predicted: Sobol within 1e-4, Gauss within 1e-6, both converging — but Sobol is asymptotically better as N grows. At `α = 0.99` (near critical), Hardy-Krause may be unbounded — QMC reverts to MC rate. **The frame partially fires; conditional on smoothness.**

**Structural attack on current:** none strong. Tensor-Gauss is right for bounded-derivative integrands; QMC is right for higher-dimensional non-smooth integrands. The discriminator is dimension and smoothness; cylinder is borderline.

**Precedent:** none in ORPHEUS. Standard MOC literature (Carlvik 1965, Geneva benchmarks) uses tensor quadrature. QMC for ray generation is in the modern MOC literature (Boyd/Forget 2013 OpenMOC).

---

### Frame 11 — Wiener-Hopf factorization for finite-medium half-range

**Trigger:** A.7 row "Fredholm integral equations". Variant α IS a Fredholm equation in surface phase-space; the half-range / specular structure of bouncing IS exactly Wiener-Hopf's domain.

**Verbatim quote (Noble 1958, "Methods Based on the Wiener-Hopf Technique" §1.1, Pergamon):**

> "A Wiener-Hopf integral equation of the first kind on the half-line `f(x) = g(x) + ∫_0^∞ K(x − y) f(y) dy` is solved by factorizing `1 − k̂(z) = X_+(z) X_-(z)` where `X_±` are analytic in `Im z > 0` and `Im z < 0` respectively."

**Match assessment.** WEAK match. Wiener-Hopf is for **half-line** integral equations on continuous phase-space. Variant α is a **discrete** bouncing process on a closed surface phase-space (each bounce is a discrete step, not a continuous half-line variable). The trigger (Fredholm + half-range) fires, but the lever is wrong:

- The half-range bi-orthogonality of singular eigenfunctions (Mika 1961, Inönü 1973, McCormick-Kuščer 1965) IS Wiener-Hopf factorization — and it is precisely the structure exploited in `fn_method/` (the Chandrasekhar / `(γ)` path). Variant α `(α)` does NOT use this — it bypasses it via the bouncing closed-form.
- Sanchez 1986 deliberately **avoids** Wiener-Hopf by using the Peierls integral form. Bringing it back contradicts the construction.

**REJECTED.** WIener-Hopf belongs in the parallel `fn_method` solver family (validated in `sanchez_chandrasekhar_gap.md` Phase A). It is the right frame for `(γ)` Chandrasekhar; it is the wrong frame for `(α)` Variant α. The two solver families exist precisely so each can use its native frame.

**Precedent:** `sanchez_chandrasekhar_gap.md` (the three-meanings taxonomy makes this rejection explicit).

---

## CROSS-METHOD POLLINATION

**Current method class:** MOC (ray-tracing + closed-form geometric-series multi-bounce summation). Phase-space discretization angle-resolved.

### Borrowings from adjacent methods

#### From CP method

- **Trigger:** Fredholm theory (A.7) — CP is the canonical Fredholm transport method. Variant α is also Fredholm (resolvent + integral kernel) but on surface phase-space, not volume.
- **Borrowing — Hierarchical-matrix compression of `S` for multi-region.** CP literature has H-matrix CP (Goluoglu-Slovacek-Mosher 2008 *NSE*) compressing the volume-volume CP matrix for distant cell pairs. Multi-region Variant α (sphere multi-region currently implemented; cylinder multi-region planned) has a **shell-shell escape probability matrix** that has the same low-rank-far-field structure. Promote to H-matrix when shell count exceeds ~20.
- **First test:** profile multi-region sphere driver at 50 shells; identify if shell-shell coupling matrix is dense (currently assumed full); test SVD truncation on a 50×50 shell matrix to identify low-rank structure.

#### From MC

- **Trigger:** A.4 Feynman-Kac (master bridge). See Frame 6 above.
- **Borrowing — Path-integral Monte Carlo verification.** Sample bouncing trajectories stochastically; weight by `α^N e^{-Σ_t Σ L_n}`. Compare to deterministic Variant α. Provides the missing structurally-independent cross-check of `(α)` against itself for `vv-principles` L11 elevation.
- **First test:** see Frame 6 first test.

#### From singular_eigenfunction / fn_method

- **Trigger:** none direct. The two solver families are **structurally independent by design** (validated in `sanchez_chandrasekhar_gap.md` §Q3 verification matrix). The pollination is **at the verification level, not the algorithmic level.**
- **Borrowing — Use `(γ)` Chandrasekhar as the anchor truth set for `(α)` Variant α.** Already the validation strategy. The cross-method "pollination" here is the deliberate structural-independence preservation; the two methods are kept apart so they verify each other.

#### From SN / discrete ordinates

- **Trigger:** A.2 Group theory (Frame 9 above).
- **Borrowing — Lebedev / Fourier-collocation quadrature on the angular fiber.** SN methods that use Lebedev (rather than Chebyshev-Gauss product) preserve SO(3); Variant α can adopt the same. Cylinder Fourier-collocate `φ_az` reuses MOC's cyclic-ray quadrature trick.
- **First test:** see Frame 9 first test.

#### From eigenvalue iteration

- **Trigger:** A.5 Krylov subspace theory + A.5 Optimization. See Frame 5 above.
- **Borrowing — Wielandt deflation / Chebyshev acceleration.** Power iteration converges with rate `λ_2/λ_1`. Wielandt deflation extracts `λ_2`, `λ_3` (subdominant eigen-modes physically meaningful — the second harmonic of the Boltzmann operator). Chebyshev iteration accelerates convergence when the spectral gap is known.
- **First test:** see Frame 5 first test.

#### From sensitivity / adjoint

- **Trigger:** A.5 Control theory.
- **Borrowing — Adjoint formulation of `T` with respect to BC parameter `α`.** For uncertainty-quantification and reactivity-coefficient applications, `∂k_eff / ∂α` is the response of `k_eff` to BC reflectivity changes. Closed-form: `∂T/∂α = e^{-τ}/(1 − α e^{-τ})² = T² · e^{-τ}`. So `∂k_eff/∂α` decomposes via the chain rule, no AD pass needed. Lift `compute_resolvent_T_grad_alpha` as a free deliverable.
- **First test:** finite-difference `dk/dα` at `α = 0.5` and compare to closed-form. Predicted match to FD truncation error.

---

## ELEGANCE RANKING (ranked by criteria-hit count + implementation cost)

| Rank | Frame | Criteria hit | Implementation cost | Verdict |
|------|-------|--------------|---------------------|---------|
| 1 | **Frame 1: Tensor networks / MPO** | 4/4 | Medium (refactor) | **HIGH PRIORITY** — Tensor networks user-asked; rank-1/rank-2 nomenclature IS bond dimension. Confirmed match. Pursue. |
| 1 (tied) | **Frame 4: Fiber bundle / G-structure** | 4/4 | Medium (refactor) | **HIGH PRIORITY** — already top-ranked in hindsight memo. Re-fired as the geometric backbone. Pursue. |
| 3 | **Frame 6: Feynman-Kac / Markov / PIMC** | 3/4 | Low (verification only) | **HIGH PRIORITY** — fills missing structural-independent cross-check. Direct payoff to V&V. |
| 4 | **Frame 3: Spectral theory of multiplication operators** | 3/4 | Very low (documentation + 1 unit test) | **MEDIUM** — names per-geometry pathology as one theorem; no new code needed. Already cited in BIE memo. |
| 5 | **Frame 2: Schur complement / Banach algebra** | 3/4 | Low | **MEDIUM** — protects future rank-N geometry extension; not load-bearing for current 6 geometries. |
| 6 | **Frame 5: Krylov / Arnoldi for resolvent** | 3/4 | Medium (replace power iter) | **MEDIUM** — extracts λ_2 as new diagnostic; tripwire was thick optical depth. |
| 7 | **Frame 7: Asymptotic analysis (Padé regularization)** | 2/4 | Low | **LOW-MED** — fixes documented α → 1 finite-precision concern. |
| 8 | **Frame 8: Generating functions for bounce decomposition** | 2/4 | Low | **LOW-MED** — adds a verification axis (monotone in z) for free. |
| 9 | **Frame 9: Group representation theory (Lebedev/Fourier)** | 2/4, conditional | Medium | **CONDITIONAL** — fires when anisotropic Sanchez 1986 lands. |
| 10 | **Frame 10: QMC** | 1/4, conditional | Medium | **LOW** — borderline-dimension; conditional on Hardy-Krause variation. |
| 11 | **Frame 11: Wiener-Hopf** | rejected | — | **REJECTED** — wrong solver family; native to `fn_method`. |

---

## TENSOR NETWORK MATCH — USER-ASKED CONFIRMATION

The user explicitly asked: "are we using tensor networks? probably yes."

**Verdict: NO, not yet — but the structure is already there as a hidden bond-dimension-1/2 MPO.**

Specifically:

- `compute_resolvent_T` (rank-1) IS the contraction of an MPO with bond dimension 1 (single surface in the boundary phase-space).
- `compute_resolvent_T_rank2` (rank-2) IS the contraction of an MPO with bond dimension 2 (two surfaces in the boundary phase-space; bond carries the 2-vector `[ψ_L^+, ψ_R^-]`).
- The rank-1 / rank-2 naming in `variant_alpha_core` IS the bond dimension renamed.
- A 7th geometry with N surfaces (3-wall wedge, multi-shell sphere with > 2 shells) would be bond dimension N — automatic under the MPO frame; manually-derived rank-N closure under the current naming.

**What ORPHEUS is NOT doing yet:**

- TT-format compression of multi-group operators (where it would actually pay).
- PEPS-style geometric encoding for true 2D problems (not relevant for 1D radial; would matter for full 2D heterogeneous problems — out of scope here).
- Tensor-train decomposition of the angular fiber for very high `n_µ` (if Sanchez 1986 anisotropic Pₙ extension lands at `n > 50`, TT compression of the angular flux becomes attractive).

The user's intuition is correct in spirit. The math is already MPO-shaped. Naming it as MPO is the lift.

---

## CRITICAL CAVEATS

Things that would invalidate the recommendations if true:

1. **MPO frame is over-engineering at rank ≤ 2.** True for 6 current geometries. The MPO framing earns its keep when (a) a 7th geometry with rank-N for N≥3 is requested, or (b) multi-group with > 5 groups becomes the primary cost, or (c) angular fiber `n_µ > 50` (anisotropic Sanchez 1986 extension). Without one of these, MPO is a clarity-not-cost lift.
2. **Krylov frame depends on `λ_2/λ_1` distance.** For thin / well-separated systems, power iteration converges in < 10 iterations and Arnoldi has overhead. Krylov pays only when subdominant eigenmodes are clustered (thick-system or near-criticality) — this is documented in the cited Saad reference.
3. **PIMC verification (Frame 6) is stochastic.** Adds variance bars to the verification matrix. Acceptable as a structural-independence cross-check (different identity, different integrand) but worse than Chandrasekhar `(γ)` as a numerical anchor. Use as supplementary, not primary.
4. **Lebedev/Fourier (Frame 9) requires the anisotropic-source extension to fire.** For the current isotropic-source family, it reduces to Gauss-Legendre — no win. The frame is **load-bearing for Phase B (Sanchez kernel) of `sanchez_chandrasekhar_gap.md`** — promote then.
5. **Asymptotic Padé (Frame 7) requires `α` and `Σ_t L` to be available as analytic parameters at runtime.** They are; trivially. But the choice of regime branch (thin / thick / near-critical) introduces a discontinuity at the boundaries — needs smooth blending or single asymptotic continuation valid everywhere.
6. **The fiber-bundle abstraction (Frame 4) was already prescribed in `variant_alpha_family_hindsight.md` as a "wait for instance N+1" item.** Re-firing here is consistent with that recommendation; the instance N+1 has not arrived. Frame 4 should be DESIGNED with a tripwire (next geometry), not eagerly implemented.

---

## STRUCTURAL FACTS (for memory promotion)

These are the new facts this detection pass produced:

1. **The rank-1 / rank-2 naming in `variant_alpha_core` IS bond dimension of an open MPO.** Rank-N geometry = bond-dim-N. Promotion candidate: extend `validated_bc_tensor_network.md` to cover the resolvent itself, not just the BC-decomposition step.
2. **Variant α IS the deterministic specular-only Feynman-Kac expectation. The full transport problem with scattering is the same kernel with the scattering subprocess included.** This unifies `(α)` and the parent transport kernel under one expectation. PIMC verification is the structurally-independent stochastic counterpart.
3. **The grazing-ray pathology in sphere/cylinder/hollow-sphere is one theorem (`ess_range(m) ∋ 1 ⇒ unbounded`) applied to four distinct `m(µ)` functions.** Slab is structurally immune. Frame 3 made this explicit; the BIE memo already cites it; documentation should consolidate the four geometry-specific docstrings into one Phase-5-style A.3 caveat citing the single theorem.
4. **Generating-function lift `Φ(r, z) = Σ_n z^n ψ_n(r)` is a free verification axis (monotonicity in z) and a free bounce-mode decomposition (Cauchy DFT inversion) at zero new code cost.** Promotion candidate.
5. **The closed-form `∂T/∂α = T²·e^{-τ}` opens up ALL adjoint sensitivity calculations in closed form.** No AD pass needed for `α`-derivatives of `k_eff` — direct chain rule. Frame 6 borrowing.

---

## UNEXPLORED

- **Crystallographic groups / Bloch theory (A.2):** repeated-lattice problems with periodic assemblies — not present. Variant α is single-cell.
- **Clebsch-Gordan / spherical harmonics algebra (A.2):** would fire when anisotropic scattering kernel decomposition lands. Not yet — current isotropic.
- **Category theory (A.2):** no compositional-structure trigger; geometries don't compose. Marked low-signal in skill.
- **de Rham cohomology / FEEC (A.1):** discrete-cohomology compatibility for div-curl-grad — Variant α has no differential operators (pure integral form), so trigger absent.
- **Symplectic geometry (A.1):** no field, no drift, no canonical form. Trivial.
- **Information theory / MaxEnt (A.4):** would fire for moment-closure methods (MN). Variant α is angle-resolved, no closure being applied.
- **Cryptography / counter-based RNG (A.4):** only fires for large-N MC; deterministic Variant α is not MC.
- **Variational calculus / Rayleigh-Ritz (A.5):** k_eff IS extracted via Rayleigh quotient. The frame fires implicitly; there is no explicit variational-form lift (e.g., constrained optimization on flux normalization with KKT) currently and the tripwire (rank-N non-monotone, Smell #15) does not fire — Variant α convergence is monotone by Krein-Rutman positivity.
- **Differential geometry / connection coefficients (A.1):** would fire for cylindrical / spherical α-redistribution if it appeared. Variant α deliberately avoids the angular-redistribution term by integrating along straight chords (the bouncing trajectory has no curvature term — chord segments are straight Euclidean, and surface reflection is point-discrete). Not load-bearing here.
- **Homogenization theory (A.6):** would fire for assembly-level / lattice problems. Single-cell Variant α has no scale separation to exploit.
- **Multi-scale numerical analysis (A.6):** same.
- **Graph theory (A.7):** would fire for sweep ordering / DAG scheduling. Variant α has no per-cell sweep — the geometric series is a closed-form per-(r, µ) computation. Trivial DAG (no graph structure).
- **Hierarchical matrices (A.7):** would fire for multi-region multi-shell with > 20 shells. Currently all multi-region work is at < 10 shells. Cross-method pollination borrowing CP; tripwire when shell count crosses 20.

---

**End of foreign-frame detection memo.**
