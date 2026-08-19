# DSA restriction/prolongation as an angular frame — structural-frame analysis (issue #2, Phase 3, 3-P0)

**Agent**: cross-domain-attacker. **Task**: is the DSA (R, P) pair an instance of
the landed frame machinery, and which discipline? **Output**: structural detection
with a verdict per question + the algebra shown; NOT an implementation plan.

**Grounding**: branch-verified against the live worktree (read directly, not Nexus —
L-005), files current as of `docs/sn-doc-architecture` @ 2026-07-26. Load-bearing
code facts: `orpheus/numerics/frame.py` (FrameBase/PetrovGalerkinFrame/GalerkinFrame),
`orpheus/numerics/basis/spherical_harmonic_basis.py` (L=0 branch),
`orpheus/numerics/quadrature/directional.py:417` (`angular_frame`),
`orpheus/numerics/iteration.py` (SourceIteration/KrylovAcceleration),
`orpheus/diffusion/operators.py` (LeakageOperator/DiffusionBoundaryOperator),
`orpheus/diffusion/boundary_realizer.py` (Marshak/half-range).

## FIVE-VERDICT SUMMARY (what 3a/3b can consume)

1. **(R, P) is GALERKIN, not Petrov–Galerkin.** It is `angular_frame(0)` — the ℓ=0
   sub-block of the existing spherical-harmonic `GalerkinFrame`. R = `analysis`
   (moment-0 = scalar flux), P = the canonical normalized reconstruction (isotropic
   injection). Π = P∘M is W-self-adjoint under the PLAIN angular quadrature measure —
   there is NO solution weighting, which is exactly what would make it PG. The measure
   that self-adjointizes it is the one `angular_frame` already carries.
2. **Consistent-DSA's low-order operator = a Schur complement of a two-moment (ℓ≤1)
   Galerkin triple product `R₁ A_high P₁`**, with the Fick/P1 closure = the Schur
   elimination of the ℓ=1 (current) block. "Consistent" (Alcouffe) = the triple product
   is taken on the **assembled DD operator** (reduce-the-discrete, not
   discretize-the-reduced — they do not commute). 3a computes exactly this.
3. **The boundary arm is the SAME Galerkin frame family over a DIFFERENT measure** —
   the half-range trace measure `|Ω·n|w`. The partial-current pair `(J⁺, J⁻)` is the
   two-sign ℓ=0 half-range moment (the analysis face); Marshak is the boundary analog
   of Fick (a closure eliminating the incoming current), not a reconstruction.
4. **Foreign frames all fire and all reduce to machinery that already exists**:
   two-grid (angular coarsening), subspace-correction/deflation (the c→1 near-null
   mode), preconditioned Richardson (SI+DSA = the `preconditioner` hook), Fourier
   (the 0.2247c residue). NO new two-grid/deflation class is needed — the
   `KrylovAcceleration.preconditioner` Callable + a `SourceIteration` correction-wrap
   are the insertion points.
5. **INSTANTIATE `angular_frame(0)`, mint NOTHING.** No `IsotropicBasis` (the constant
   IS `SphericalHarmonicBasis(L=0)`), no ad-hoc R/P pair. `Quadrature.angular_frame` is
   already the single source of the angular moment projection (scattering §5.6 kernel +
   in-sweep moment accumulation share it); DSA's R/P must be its ℓ=0 faces, NOT a
   fourth spelling (Smell #16, shape-4 — the pair fires before the code exists).

---

## STRUCTURAL FEATURES

- **Mathematical objects**: the transport residual `r = (A−S−B)ψ − q` (an
  `AngularFlux`, per-ordinate over S²); the scalar flux φ (its ℓ=0 moment); the
  low-order diffusion operator `A_low = −∇·D∇ + Σ_a`; the partial currents `(J⁺, J⁻)`.
- **Symmetry present**: SO(3) rotational invariance of the scattering kernel
  `Σ_s(Ω·Ω')` (zonal ⟹ Funk–Hecke). The ℓ=0 (constant) subspace is the trivial-irrep
  V₀ of SO(3). This is the symmetry that makes the frame Galerkin (L-009).
- **Symmetry absent**: the DD SPATIAL operator is upwinded/non-self-adjoint
  (characteristic-triangular per L-007) — the low-order operator inherits this, so
  "Galerkin" holds on the ANGULAR axis only.
- **Iterative structure**: a fixed-point (source iteration = Richardson on the
  moment-space operator `(I − TΣ_s)`, T = sweep∘moment); DSA is a two-grid correction
  inside it. The eigenvalue outer (`power_iteration`) sits above.
- **Integral structure**: R = `∫·dΩ` against the constant (moment-0); a compact
  angular reduction from N ordinates to 1 moment.
- **Differential structure**: the low-order operator is elliptic (`−∇·D∇`); the P1
  closure is Fick `J = −D∇φ`.
- **Boundary handling**: half-range `|Ω·n|w` trace measure; Marshak (vacuum = 𝒜=0),
  reflective (𝒜=1) — the #290 albedo family `J⁻ = 𝒜 J⁺`.
- **Scale separation**: thick/diffusive (c→1) is exactly where SI stalls and DSA
  deflates — the multiscale trigger DSA exists to answer.
- **Elegance-detector site**: an R and a P about to be hand-built as two independent
  operators, over a projector (moment-0 / isotropic) that ALREADY has ≥2 spellings.

## ELEGANCE DETECTOR HITS

- **Smell #16 (distinct paths/reps to ONE operator) — shape-4 (a path about to be
  written), the dominant hit.** The moment-0 angular projection already exists as:
  (1) `angular_frame(L).analysis` — the canonical frame face; (2) the in-sweep moment
  accumulation (the frame docstring states it shares THIS object); (3) the scattering
  ℓ=0 in-scatter (`full_scatter_kernel = frame.conjugate(Λ+N2n)`, and the
  `add_iso_source`/`apply_p0_in_scatter` fast path). A hand-rolled DSA `R = Σ_n w_n r_n`
  would be a FOURTH spelling. Fix: consume `angular_frame(0).analysis`. Discriminator
  is `array_equal` (0 ULP), not `allclose`: FP non-associativity makes a separate
  einsum differ at ULP from the frame's fused `einsum('n,nlm,n...->lm...')`.
- **Smell #16 — shape-2 (one quantity, two representations bridged by hand).** The
  isotropic prolongation and the moment-0 restriction are the two faces of ONE frame;
  building them as an unrelated `(restrict, inject)` pair invites an ASYMMETRIC
  normalization (the classic DSA inconsistency where R and P disagree by a 4π/2 factor
  and silently break conservation). The frame makes `P = M⁺` (and Π self-adjoint) a
  theorem.
- **Frame-leak naming (L-006), minor.** Naming the pair `dsa_restrict`/`dsa_prolong`
  after the CONSUMER (DSA) rather than the ROLE (ℓ=0 moment / isotropic frame face)
  becomes a lie the moment a second consumer (PN synthesis, CMFD, a P1 detector
  functional) reads the same ℓ=0 frame. Name the role: `angular_frame(0)`.

---

## Q1 — WHICH FRAME IS (R, P)?

**Trigger**: an analysis/reconstruction pair `R∘A∘P` on a change-of-basis onto a
1-dimensional coarse space (the constant on S²), with SO(3)-zonal symmetry in the
operator being accelerated (L-009: symmetry ⟹ eigenbasis-owned Galerkin frame).

**Reformulation** — the pair IS `angular_frame(0)`, the ℓ=0 sub-block of the SH
`GalerkinFrame`. Read off the code (`spherical_harmonic_basis.py`, L=0 branch:
`Y[:,0,0]=1`, `addition_theorem_factor = 2·0+1 = 1`, Gram `4π/(2ℓ+1) = 4π`):

- **Restriction** `M = analysis`: `(Mψ) = Σ_n w_n · Y₀(Ω_n) · ψ_n = Σ_n w_n ψ_n = φ`
  — the moment-0, i.e. the scalar flux (the un-normalized 0th moment `∫ψ dΩ`).
- **Reconstruction** `R`: `(Rφ)_n = (2·0+1)·Y₀·φ = φ` — the isotropic broadcast,
  factor 1.
- **Coefficient Gram** `G = 4π` (`= 4π/(2ℓ+1)` at ℓ=0); `project = G⁻¹M = φ/4π` = the
  isotropic AVERAGE.

**The Galerkin vs PG discrimination, worked (S₂ slab, the smallest case that shows
it).** Two ordinates `μ = ∓1/√3`, weights `w₁=w₂` (embedded on S² so `Σw = 4π`;
`W = diag(w)`). The isotropic projector is `Π = P∘M` with the section `P = R∘G⁻¹`
(so that `M∘P = I`: `∫(δφ/4π)dΩ = δφ`):

```
Π ψ = P(Mψ) = (1/4π)·(w₁ψ₁ + w₂ψ₂)   (broadcast, constant in n)
Π   = (1/4π)·[[w₁, w₂],
              [w₁, w₂]]              ⟹  (WΠ)_nm = w_n·w_m/4π   — SYMMETRIC in (n,m)
```

`WΠ` symmetric ⟹ Π is **W-self-adjoint** ⟹ orthogonal projector ⟹ **GALERKIN**
(`Π* = Π`, and `M* = R` up to the diagonal Gram: `M* = 4π·R` here). There is **no
place a solution weight (φ, φ\*) enters** — the constant is tested against itself under
the PLAIN L²(S²) measure. Contrast spatial homogenization (the headline PG consumer,
`frame.py:343`): there the test is `φ·1_R` (flux-weighted) ⟹ `M* ≠ R` ⟹ PG. DSA's
angular restriction has NO such weight; the "flux-weighted?" branch the task flags
resolves to **no — plain measure, isotropic average**.

**This is L-009 discharged on the ℓ=0 block**: `Σ_s(Ω·Ω')` is SO(3)-zonal ⟹ Funk–Hecke
diagonalizes it in `{Y_ℓ}` with the Legendre moments as eigenvalues ⟹ the SH frame is
scattering's EIGENBASIS ⟹ `GalerkinFrame` BECAUSE `Σ_s` is self-adjoint-zonal. The ℓ=0
block inherits Galerkin as the V₀ (trivial-irrep) sub-block. The measure that
self-adjointizes it is precisely `s2_measure` (`weights = self.weights`) that
`angular_frame` already binds (`directional.py:451`).

**Elegance payoff**: Structure-exposing (the R/P pair IS a spectral-projector face
pair, `P = M⁺`, `Π = PM` orthogonal — not two hand-built operators). Structurally-
simpler (one `angular_frame(0)` object replaces a restrict+inject pair; `P = R.H`-up-
to-Gram is free via the metric-aware `_AdjointOperator`). Expressive (the same object
the scattering ℓ=0 in-scatter and in-sweep moment accumulation already use).

**First test (discriminates Galerkin from PG — can RED)**: build `M =
angular_frame(0).analysis`, form the dense `Π = P∘M` on the ordinate space, assert
`W@Π` is symmetric to machine precision. A flux-weighted (PG) restriction — a naive
implementer weighting the moment by `|μ|` or by the current-flux — gives a
**non-symmetric** `WΠ` and REDs. (An `allclose(P.apply, isotropic_broadcast)` test
would PASS for both and is theatrics.)

**Structural attack on current**: a DSA that hand-builds `R = Σ w_n r_n` and `P =
(1/4π)·broadcast` as two independent operators leaves the adjoint relation `P = M⁺`
(equivalently `M* = R·G`) UN-encoded — the symmetry of `Π` becomes a coincidence the
first normalization typo breaks, and the frame's free `.H` is re-derived by hand.

**VERDICT Q1**: `(R, P) = angular_frame(0)`, a `GalerkinFrame` (Π\* = Π under the plain
quadrature measure). NOT Petrov–Galerkin — no solution weighting on the test side.

---

## Q2 — THE TWO-GRID / CONSISTENCY IDENTITY

**Trigger**: a low-order operator claimed "consistent" with a high-order one — the
signature of a Galerkin coarse-grid operator `A_low = R A_high P` (multigrid), possibly
post-composed with an elimination of retained-but-closed unknowns.

**Reformulation** — express 3a as named frame algebra:

1. **Moment reduction = a Galerkin triple product on the ℓ≤1 frame.** Let `R₁ =
   angular_frame(1).analysis` (restrict to `{Y₀, Y₁}` = the (φ, J) moments) and
   `P₁` its normalized reconstruction (the ℓ≤1 addition-theorem synthesis). The
   P1/moment reduction of the assembled transport operator is the **triple product**
   `Ã = R₁ A_high P₁`, a 2×2 block operator on moment space:

   ```
   Ã = [ A_φφ  A_φJ ]      (the two-moment Galerkin projection of A_high;
       [ A_Jφ  A_JJ ]       test = trial = {Y₀,Y₁} ⟹ GALERKIN on the angular axis)
   ```

   This reduction is Galerkin (test=trial=ℓ≤1 harmonics) — it is the standard
   PN/spherical-harmonics Galerkin projection; the P1 CLOSURE is the *assumption* that
   the trial space is ℓ≤1 (`ψ ≈ ψ₀ + 3 Ω·J`), which IS the truncation `P₁`.

2. **The Fick/P1 closure = the Schur complement over the odd-moment block.** The
   diffusion form eliminates the ℓ=1 current `J`:

   ```
   A_low = A_φφ − A_φJ · A_JJ⁻¹ · A_Jφ      (Schur complement of Ã over the ℓ=1 block)
   ```

   `A_JJ⁻¹ A_Jφ` IS Fick's law made discrete (`J = −D∇φ`); the Schur complement is the
   condensed diffusion operator. Structurally identical to the interior-current
   condensation the landed `LeakageOperator` already performs
   (`operators.py:78`: "interior face currents stay CONDENSED (never trace DOFs)").

3. **"Consistent" (Alcouffe) = reduce-the-DISCRETE, not discretize-the-reduced.** The
   two operations DO NOT COMMUTE: `R₁·(assemble A_high)·P₁` ≠ `assemble(R₁·A_high·P₁
   continuous)`. Consistency demands the FORMER — the triple product on the **assembled
   DD matrix** (which is exactly what Phase 2's assembly enables, and what 3a computes).
   Inconsistent DSA (an independently-discretized diffusion operator) is the
   discretize-then-reduce path — it diverges for optically-thick cells (McCoy–Larsen).

**Elegance payoff**: Structure-exposing (the ad-hoc "P1 moment reduction recipe"
becomes `Schur(R₁ A_high P₁)` — a named triple product + a named complement; the Fick
closure is located precisely as the odd-block Schur elimination). Algorithmic-advantage
(3a's matrix comparison `A_low ≟ LeakageOperator` is a COMPUTED consistency proof —
Alcouffe's hand-derived stencil recovered by object-level algebra; if they differ, the
diff's STRUCTURE (boundary row? interior stencil? scaling?) is what R4 rules on).

**First test (discriminates consistent from inconsistent — can RED)**: assemble the DD
`A_high` (Phase 2 machinery), form `A_low = Schur(R₁ A_high P₁)` as a matrix, compare
**as matrices** to the diffusion family's `LeakageOperator` assembly on the same mesh.
On an optically-thick uniform slab, an INCONSISTENT (continuous-then-discretized)
diffusion operator differs from `Schur(R₁ A_high P₁)` in the interior stencil by O(h)
terms that a hand-imposed operator lacks; the matrix diff REDs the "the landed
LeakageOperator is automatically the consistent partner" claim IF it is not (this is 3a
by construction — R4 is deliberately "decide from the diff").

**Structural attack on current**: naming 3a "the discrete P1 moment reduction" hides
that it is a Galerkin triple product followed by a Schur complement — two named,
independently-testable structures. Without the names, the non-commutativity
(reduce-discrete ≠ discretize-reduce) that DEFINES consistency is invisible, and the
Fick closure's location (the odd-block elimination) is folded into an opaque recipe.

**VERDICT Q2**: consistent-DSA's low-order operator = `Schur_{ℓ=1}(R₁ A_high P₁)` — a
two-moment Galerkin triple product with the ℓ=1 block Schur-eliminated (that
elimination IS the Fick/P1 closure). "Consistent" = the triple product is taken on the
ASSEMBLED DD operator (reduce-the-discrete). 3a's matrix comparison is a computed
consistency proof.

---

## Q3 — THE BOUNDARY ARM

**Trigger**: a boundary restriction on the same axis but under a different (trace)
measure — the measure-is-metric principle (`frame.rst`: the measure carries the axis +
the L² metric). The #290 seam already states it: "the ℓ=0 half-range moment under the
SHARED `|Ω·n|w` metric (partial currents J±)".

**Reformulation** — the boundary arm is the SAME Galerkin frame family, `basis =
SphericalHarmonicBasis(L=0)` (the constant), bound to a DIFFERENT `DiscreteMeasure`:
the **half-range trace measure** `|Ω·n| w` restricted to a half-sphere `Ω·n ≷ 0`.

- **Analysis face (two-sign ℓ=0 half-range moment)** — the partial currents:

  ```
  J⁺ = ∫_{Ω·n>0} |Ω·n| ψ dΩ = Σ_{n: μ_n·n̂>0} |Ω_n·n̂| w_n ψ_n       (outflow)
  J⁻ = ∫_{Ω·n<0} |Ω·n| ψ dΩ = Σ_{n: μ_n·n̂<0} |Ω_n·n̂| w_n ψ_n       (inflow)
  ```

  This is moment-0 (against the constant) under the measure `|Ω·n|w` on each half. The
  weight `|Ω·n|` is the natural boundary/current metric (flow through the face) — a
  DiscreteMeasure swap, NOT a new frame TYPE. The angular carrier already exists:
  `AngularBoundaryFlux.net_current(face)` was minted in Phase 1 as "the angular sibling
  of the scalar `J⁺−J⁻` — the SAME family DSA's restriction consumes later"
  (roadmap:217).

- **The closure is Marshak, the boundary analog of Fick — NOT a reconstruction face.**
  Just as the interior arm (Q2) eliminates the ℓ=1 current via Fick (odd-block Schur
  complement), the boundary eliminates the incoming partial current via the albedo law
  `J⁻ = 𝒜 J⁺` (`boundary_realizer.py`: vacuum 𝒜=0 = Marshak zero-inflow; reflective
  𝒜=1). The 2-D trace change of basis `φ_Γ = 2(J⁺+J⁻)`, `J_net = J⁺−J⁻` (the P1
  dictionary, `operators.py:83`) is a fixed invertible 2×2 on the trace, not a frame
  synthesis. So: **bulk arm = full-range ℓ=0 (1 moment); boundary arm = half-range ℓ=0
  ×2 signs (2 moments) under `|Ω·n|w`; both Galerkin; Fick and Marshak are the same
  structural move (eliminate an odd/incoming moment) at interior vs boundary.**

- **The correction problem's BC is exactly right by construction** (roadmap 3b.2): the
  DSA error problem has zero inflow, so vacuum→Marshak (𝒜=0) is the correct
  correction-equation BC — the #290 vacuum law IS the error problem's boundary.

**Elegance payoff**: Structurally-simpler (bulk and boundary restrictions are ONE frame
construction parametrized by the measure `w` vs `|Ω·n|w` — not two operators).
Structure-exposing (Fick and Marshak unify as odd/incoming-moment Schur eliminations).
Expressive (reuses the landed `AngularBoundaryFlux.net_current` and the
`DiffusionBoundaryOperator` `(J⁺,J⁻)` trace).

**First test (discriminates the shared-metric claim — can RED)**: assert the boundary
restriction consumes the SAME `|Ω·n|w` metric the `DiffusionBoundaryOperator` trace
uses — build `J_net` from the SN residual's `AngularBoundaryFlux.net_current(face)` and
assert it equals the boundary row the diffusion correction consumes to solver tolerance.
A boundary R that drops the `|Ω·n|` weight (uses the plain `w`, i.e. a full-range-style
moment at the boundary) REDs against the current-conserving diffusion trace — the metric
mismatch is a `J`-scale error the plain-measure moment cannot see. (Mode-12 note: a
conservation-only check `⟨1, Rr⟩ = ⟨1, r⟩` is BLIND to the `|Ω·n|` weight; the
discriminator must compare the WEIGHTED partial current, not the unweighted sum.)

**Structural attack on current**: treating the boundary restriction as "moment-0 at the
boundary" without naming the `|Ω·n|w` measure hides that it is a DIFFERENT frame (trace
metric) from the bulk — a bare-`w` moment at the boundary silently computes the wrong
(non-current-conserving) quantity, the exact ERR-067-family hazard (a Euclidean vs
metric-weighted trace).

**VERDICT Q3**: boundary arm = the ℓ=0 `GalerkinFrame` over the half-range trace measure
`|Ω·n|w` (same basis, same discipline, different measure — measure-is-metric). Analysis
face = the partial-current pair `(J⁺,J⁻)`; Marshak = the boundary Fick (an
incoming-moment Schur elimination), not a reconstruction. Reuses
`AngularBoundaryFlux.net_current` + the `DiffusionBoundaryOperator` trace.

---

## Q4 — FOREIGN-FRAME SWEEP

Current method class: **SN + eigenvalue/fixed-point iteration**; DSA itself is an
acceleration/preconditioning technique. Four foreign frames fire; all reduce to
existing machinery.

**Frame: two-grid multigrid (angular coarsening).**
Trigger: a coarse-space correction `x ← x + P A_low⁻¹ R (b − Ax)` with R/P grid
transfers. DSA IS a two-grid V-cycle — the "grids" are angular (fine = N ordinates,
coarse = ℓ=0 isotropic), R/P = the `angular_frame(0)` faces (Q1), the coarse operator =
`A_low` (Q2). Payoff: names the DSA structure as a bona-fide two-grid method, so the
Galerkin-coarse-operator condition `A_low = R A_high P` (Q2) is the multigrid
consistency condition, not a DSA-specific coincidence. First test: the correction→0 gate
(roadmap 3c) IS the multigrid "coarse correction vanishes at the fine fixed point"
property — assert `P A_low⁻¹ R (b − A ψ*) → 0` at the converged transport fixed point;
a wrong R/P normalization (Q1's asymmetry) leaves a residual correction and REDs.

**Frame: subspace correction / deflation.**
Trigger: the correction targets the near-null slow mode of the iteration operator. SI's
iteration operator `(I − TΣ_s)` has its slowest mode as `c→1` (the flat, diffusive
mode = the ℓ=0 near-null space); DSA deflates exactly that mode by solving it in the
diffusion subspace. Payoff: explains WHY ℓ=0 suffices (the slow mode is isotropic) and
WHEN it fails (anisotropic scattering pushes slow modes into ℓ=1 — Morel 1982, the
pollination trigger below). First test: measure the SI iteration operator's dominant
eigenmode on a homogeneous c≈0.99 slab; assert it is (near-)isotropic — if a strongly
anisotropic config puts the slow mode in ℓ=1, ℓ=0-only DSA under-deflates and the
measured ρ exceeds the 0.2247c bound (RED ⟹ arm-2 trigger).

**Frame: preconditioned Richardson (the exact preconditioner form).**
Trigger: SI is Richardson; DSA is its preconditioner. The moment-space fixed point is
`(I − TΣ_s)φ = Tq`, T = sweep∘moment. SI+DSA is Richardson with the **left
preconditioner**

```
M_DSA⁻¹ = I + A_low⁻¹ Σ_s        (the additive Alcouffe correction on the SI increment)
φ^{k+1} = φ^{k+½} + A_low⁻¹ Σ_s (φ^{k+½} − φ^k),   φ^{k+½} = SI(φ^k)
```

so the preconditioned iteration operator is `(I + A_low⁻¹Σ_s)(I − TΣ_s)`, spectral
radius ≤ 0.2247c. Payoff: this is EXACTLY the `KrylovAcceleration.preconditioner`
Callable (`iteration.py:198`, `M ≈ A⁻¹`) — GMRES-preconditioned-by-DSA is folding #200's
intent with ZERO new iteration machinery, and SI+DSA is a `SourceIteration` correction-
wrap. First test: assert the DSA-preconditioned GMRES on a c≈0.999 slab converges in
O(1) iterations independent of mesh, vs unpreconditioned O(1/(1−c)); a preconditioner
that is NOT the consistent `A_low⁻¹Σ_s` (an inconsistent diffusion op) fails the
mesh-independence and REDs.

**Frame: Fourier / spectral (the 0.2247c residue).**
Trigger: constant-coefficient iteration operator ⟹ diagonalized by Fourier modes
`e^{iλx}`. On the infinite homogeneous medium, SI's iteration operator has Fourier
symbol `ρ_SI(λ) = c·(1 − (tan⁻¹Λ)/Λ)`-family peaking at `c` as `λ→0` (the flat mode).
DSA replaces the `λ→0` eigenvalue with ~0 (the diffusion limit is exact there); the
residual peak sits at an INTERMEDIATE λ where the diffusion approximation is least
accurate — `max_λ ρ_DSA(λ) ≈ 0.2247c`. Payoff: the 0.2247c is not a magic number — it
is the non-deflatable Fourier remainder, and it is the gate target (roadmap 3c). First
test: measure ρ on the classic periodic homogeneous problem across c; assert it tracks
0.2247c (not 0, not c) — a scheme that hits ρ≈0 would be over-claiming (the diffusion
op cannot deflate the intermediate mode); a scheme at ρ≈c has a broken correction.

**Existing machinery to reuse (the leverage question)**: NO two-grid or deflation class
exists in `orpheus/numerics/`. The insertion points are (a) `KrylovAcceleration`'s
`preconditioner: Callable[[ndarray], ndarray]` hook (already the M≈A⁻¹ slot — DSA
becomes a preconditioner Callable), and (b) a `SourceIteration` correction-wrap (the
driver already threads `initial_guess` and computes the equation residual `rhs_{n-1} −
rhs_n` = `Σ g_i(ψ_{n-1}−ψ_n)` per step, which is the `Σ_s`-increment DSA corrects). The
roadmap 3b.4 already commits to both — this analysis confirms it is the right layering
(L-007: DSA is a preconditioner of the shared resolvent iteration, not a new engine).

---

## Q5 — THE ANTI-MINT VERDICT

**What 3b should instantiate**: `Quadrature.angular_frame(0)` — the EXISTING
`GalerkinFrame` over the EXISTING `SphericalHarmonicBasis(L=0)`. R = its `.analysis`
face; P = its normalized reconstruction (`.reconstruction ∘ project`-gram, or
equivalently the metric-adjoint `M⁺` — all frame-native via the free `.H`). For the
consistency derivation (3a), `angular_frame(1)` (ℓ≤1) for the triple product, Schur on
the ℓ=1 block (Q2).

**Mint NOTHING new**:
- **No `IsotropicBasis` / `ConstantBasis`.** The isotropic (constant) function IS
  `SphericalHarmonicBasis(L=0)` — `Y₀ = 1`, already the `angular_frame(0)` trial basis.
  Minting a rank-1 `ConstantBasis` FORKS `R∘M` and ULP-drifts vs the fused SH kernel
  (my prior `iso_source_frame_conjugation_unification` verdict: the separate rank-1
  `iso_frame` was the WRONG object for exactly this reason — it is `angular_frame(0)`,
  the ℓ=0/V₀-trivial-irrep sub-block of the harmonic frame, `M₀ = ∫ = ℓ=0 analysis row`,
  `R₀ = broadcast = P₀ = 1`).
- **No ad-hoc R/P pair.** The pair is one frame's two faces (Q1).
- **No new two-grid/deflation class** (Q4) — reuse the `preconditioner` hook + a
  `SourceIteration` wrap.

**The single-source-of-truth resolution (which existing spelling is the shared
primitive)**: `Quadrature.angular_frame` is ALREADY declared "the single source of the
angular frame consumed by `ScatteringOperator` — its §5.6 kernel AND the in-sweep
moment accumulation share THIS object" (`directional.py:434`). The moment-0 projection
therefore already lives in one canonical place. DSA's R/P = its ℓ=0 faces. This is the
Smell-#16 shape-4 fix applied BEFORE the code is written: the pair fires as a
would-be-fourth-spelling, and the fix is to consume `angular_frame(0)`, not to
twin it.

- The scattering P0 fast path (`add_iso_source` / `apply_p0_in_scatter`, the `skip_l0`
  branch) is a SEPARATE reaction-rate spelling of the ℓ=0 in-scatter. Whether IT should
  also be unified onto the frame is scattering's concern (flagged in my iso_source memo
  as a latent fork), NOT DSA's — DSA must not add a spelling regardless. **Do not
  design that refactor here; flag its blast radius only**: unifying the P0 fast path
  onto `angular_frame(0)` touches `scattering.py` (the `add_iso_source`/`skip_l0`
  path), `material_xs_field.py` (`apply_p0_in_scatter`), and the P0-fastpath perf
  diagnostic (`diag_276_scattering_p0_fastpath_perf.py`) — a scattering-side campaign,
  out of Phase 3 scope.

**Elegance payoff**: Structurally-simpler (zero new types; DSA is a consumer of landed
frame + diffusion machinery). Expressive (R = `angular_frame(0).analysis` reads as "the
ℓ=0 moment," the domain term). Algorithmic-advantage (the frame's fused einsum + free
`.H` — no hand-rolled restrict/inject to drift).

**First test (the anti-twin discriminator — `array_equal`, 0 ULP)**: assert
`dsa_restriction.apply(r)` is BITWISE identical to `angular_frame(0).analysis.apply(r)`
AND to the ℓ=0 slice of the in-sweep moment accumulation, on a heterogeneous non-flat
residual. A hand-rolled `np.einsum('n,n...->...', w, r)` (or `(w[:,None]*r).sum(0)`)
differs at ULP by FP non-associativity and REDs — proving there is ONE spelling. (An
`allclose` here is theatrics — it passes for a twin.)

**Structural attack on current**: a DSA R/P built as fresh operators is a fourth
spelling of a projector with an already-canonical home; the codebase's "which of the
four moment-0 spellings is authoritative?" confusion is exactly the Architecture
cardinal-rule violation, and the twin drifts silently (value-correct-by-coincidence)
until a convention change (a weight normalization, an ordinate reorder) moves one
spelling and not the others.

**VERDICT Q5**: instantiate `angular_frame(0)` (bulk) / `angular_frame(1)`+Schur (3a
low-order) — mint nothing. The shared primitive is `Quadrature.angular_frame`; DSA's
R/P are its ℓ=0 faces, pinned bitwise against the existing moment spelling. The
scattering-P0-fastpath unification is a separate, flagged, out-of-scope campaign.

---

## CROSS-METHOD POLLINATION

Current method: **SN acceleration / eigenvalue iteration**.

**From diffusion (the low-order partner — already landed).** Trigger: the DSA low-order
operator IS an elliptic diffusion operator. Borrowing: `LeakageOperator` +
`DiffusionBoundaryOperator` + `DiffusionMesh.from_material_mesh(sn_mesh)` (roadmap 3b.2
— same axes/materials/BCs by construction) supply the coarse solve directly; the P1
face closure (`_boundary_closure`) IS the Marshak arm (Q3). First test: 3a's matrix
comparison `Schur(R₁ A_high P₁) ≟ LeakageOperator` (Q2) — if equal, the diffusion
family is the consistent partner for free.

**From Krylov/eigenvalue (the preconditioner posture — SAILOR).** Trigger: a slow
fixed-point iteration with an SPD-ish near-null mode. Borrowing: Adams & Larsen 2002's
transport-corrected preconditioning — DSA as a GMRES left preconditioner via the
existing `KrylovAcceleration.preconditioner` hook (Q4). First test: DSA-preconditioned
GMRES iteration count is mesh- and c-independent on a c≈0.999 slab; plain-sweep-
preconditioned GMRES is not.

**From CMFD / nonlinear diffusion acceleration (the degradation escape hatch).**
Trigger: consistent linear DSA degrades for heterogeneous / multi-D (Warsa–Wareing–Morel
2004: ρ→0.88 with "no known remedy" for fully-consistent DSA on skewed cells). Borrowing:
CMFD's nonlinear `D̂` correction factor — the low-order operator is re-derived each outer
with a nonlinear current-conservation factor, restoring `M∘(low-order) = M∘(transport)`
exactly at each iterate. First test (deferred, arm-2 trigger): measure ρ on a
high-contrast heterogeneous slab; if consistent-DSA's ρ exceeds ~0.5 (well above
0.2247c), the nonlinear D̂ is the borrowing. NOT in Phase-3 scope (R5: decide from the
measured table).

**From PN / higher-moment synthetic acceleration (Morel 1982).** Trigger: strongly
anisotropic scattering — the slow mode leaves ℓ=0 (Q4 deflation frame). Borrowing:
Morel's ℓ≤1 (or MDSA) synthetic acceleration coarsens to `angular_frame(1)` instead of
`angular_frame(0)` — the SAME frame family at higher L, so it is a one-parameter change,
NOT a new mechanism. First test: with P1 scattering `Σ_s1 ≠ 0` and c→1, ℓ=0-only DSA's
measured ρ exceeds 0.2247c while ℓ≤1 DSA does not — the trigger to bump L. (This is why
Q5 keeps the frame parametrized by L, not hardcoded to 0.)

---

## UNEXPLORED

- **Homology / chain complex** — no `∂²=0`. The R/P adjoint pair is a biproduct
  inject/project (a section/retraction), NOT a differential; two reflections/moments do
  not compose to zero (L-001). No homology payoff.
- **Wiener–Hopf / Chandrasekhar H-function** — wrong solver family (L-001): native to
  the half-space singular-eigenfunction line, structurally incompatible with the
  sweep+diffusion-correction formulation. Keeping them independent is itself a V&V
  requirement.
- **Tensor networks / MPO** — the angular coarsening is rank-1 (ℓ=0) or rank-4 (ℓ≤1,
  the `{Y₀,Y₁^{-1,0,+1}}` moments) FIXED; no bond-dimension truncation knob. Bond-dim-1/4
  degenerate, not a network (L-001).
- **Category theory / operad / PROP** — the compositional structure (`R∘A∘P`, the
  frame-conjugate 2-cell) is already captured concretely by the double-category frame +
  the biproduct (Q1/Q2); no abstract-nonsense lever produces a test (L-001).
- **Differential geometry / Christoffel** — no curvature term in the angular DSA R/P
  (the `(1−μ²)/r ∂_μ` redistribution lives on the curvilinear STREAMING, not the DSA
  coarsening, and curvilinear DSA is OUT — blocked on #282). Not triggered here.
- **Optimal control / Riccati** — the Fick/Marshak closure (Q2/Q3) is a single STATIC
  Schur elimination, not a dynamic-programming Riccati recursion; no time/DP sweep
  structure. The Schur complement is named directly (Q2); no Riccati frame needed.
- **Wavelet / hierarchical (multi-level angular)** — a full angular multigrid (P0→P1→…→SN
  V-cycle) is the multi-LEVEL generalization of the two-GRID DSA (Q4). Not triggered for
  Phase 3 (two-grid suffices for isotropic-dominated slow modes; the multi-level trigger
  is Morel-anisotropy, handled by the PN pollination above, not a wavelet frame).
