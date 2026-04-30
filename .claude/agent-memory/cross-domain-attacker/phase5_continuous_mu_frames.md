---
name: Phase 5 continuous-µ multi-bounce — frame match inventory
description: Six frames matched against the continuous-µ kernel K_bc^mb = ∫ G_in·F_out·f_mb(µ) dµ. Strongest matches: Hilbert-Schmidt separable kernel (the right object — µ is the rank-1 axis, not r), generating-function bounce-resolved expansion, spectral theory of multiplication operators (gives the named theorem behind the per-geometry pathology), and Feynman-Kac surface Markov chain (verification reference). Weak/refuted: Wiener-Hopf, Mittag-Leffler partial fractions, group-theoretic Marshak-basis switch.
type: project
---

# Phase 5 continuous-µ multi-bounce frame attack (2026-04-28)

Branch `feature/peierls-specular-bc`. Question for cross-domain-attacker:
what foreign mathematical structures match the Phase 5 continuous-µ
multi-bounce kernel so the implementation borrows the right machinery
rather than re-deriving from first principles?

## Problem in one paragraph

`K_bc^mb,sph(r_i, r_j) = 2 ∫_0^1 G_in(r_i, µ) · F_out(r_j, µ) · f_mb(µ) dµ`
with `f_mb(µ) = µ/(1 - e^(-σ·2Rµ))`. The integrand has a removable pole
at µ = 0 (cancellation `µ → 0` against `1 - e^(-2σRµ) → 2σRµ + O(µ²)`).
The current matrix-Galerkin form `K_bc = G·R·(I-T·R)^(-1)·P` projects
this onto a rank-N polynomial basis and *loses* the cancellation —
operator norm of `(I-T·R)^(-1)` grows unboundedly with N for sphere.
Slab is structurally immune (chord = L/µ → ∞ at grazing); cyl is
operationally similar to sphere (different mechanism: R-conditioning).

## Strongest matches (use these for Phase 5)

### M1 — Hilbert-Schmidt / separable kernel (A.7)

**Trigger**: kernel is `a(r_i, µ) · b(r_j, µ) · c(µ)` integrated over µ.
This is the textbook degenerate-kernel form. Rank-in-µ is 1 at every µ;
total rank = number of µ-quadrature points Q.

**What it tells us**: the natural finite-dimensional axis of this
operator is µ, not r. Mode-projecting in r is projecting along the
WRONG dimension. The matrix-Galerkin form is structurally inverted.

**Implementation**: `K_bc^mb[i, j] = Σ_q w_q · G_in(r_i, µ_q) · F_out(r_j, µ_q) · f_mb(µ_q)`.
Indices: i, j over spatial cells, q over µ-quadrature. Storage =
O(I·Q + J·Q); cost = O(I·J·Q). NO matrix inverse. NO rank-N basis.
NO `R = (1/2)M^{-1}` conditioning hazard.

**First test**: implement for homogeneous sphere at τ_R=2.5, Q=64
Gauss-Legendre. Compare to MC k_eff = k_inf (verified in
`diag_specular_overshoot_05`). Expected: rel error < 0.01%.
Discriminator: matrix-Galerkin at N=8 gives +1.61%; HS at Q=64 < 0.1%.

**Status**: STRONG. Direct production candidate.

### M2 — Generating function / Bose-Einstein bounce-resolved expansion (A.4)

**Trigger**: `1/(1-e^(-x)) = Σ_(k=0)^∞ e^(-kx)` — the kernel is the
generating function of bounce count k.

**Reformulation**: bounce-by-bounce. `K_bc^mb = Σ_k K_bc^(k)` with
`K_bc^(k) = 2 ∫_0^1 G_in F_out µ e^(-k·σ·2Rµ) dµ`. Each integral is
non-singular. Truncation at K_max gives explicit residual bound.

**Why this is better than M1 in some cases**: gives a physical
meaning to truncation. K_max ≈ 1/(σR) + log(1/tol) bounces. Useful
when the user wants to study "how much do bounces contribute."

**First test**: at τ_R=2.5, plot ‖K_bc^(k)‖ for k=0..10. Geometric
decay with ratio ≈ e^(-σ·2R·µ_eff). K=10 truncated sum matches HS
direct integration to 1e-5.

**Status**: STRONG. Same physics as M1 but reveals bounce-resolved
structure. Use both: M1 for production, M2 for verification.

### M3 — Spectral theory of multiplication operators (A.3)

**Trigger**: continuous-limit `T·R` is multiplication by `φ(µ) =
e^(-σ·2Rµ)` (sphere) or `e^(-σL/µ)` (slab). Spectrum = essential range
of φ. `(I - T·R)` invertible-bounded iff `1 ∉ closure(ess_range(φ))`.

**What it tells us — the per-geometry pathology in one theorem**:

| Geometry | ess_range(φ) | 1 in closure? | (I-TR)^{-1} bounded? |
|----------|--------------|---------------|----------------------|
| Sphere | (e^(-2σR), 1] | YES (at µ→0) | NO |
| Slab | [0, e^(-σL)] | NO | YES |
| Cyl | (4/π·Ki_3(τ_2D)·cos α range) | NO (cos α→0 wins) | YES |

So **slab is bounded** (already known empirically — diag_phase4_03);
**cyl is bounded** (pathology is R-conditioning, not the kernel —
already known empirically); **sphere is unbounded** (the source of
the matrix-Galerkin divergence). The frame gives a NAMED theorem
(matrix Galerkin of unbounded multiplication operators does not
converge — Trefethen-Embree, *Spectra and Pseudospectra*).

**Implementation**: motivates the µ-direction discretization (M1).
The frame is more diagnostic than algorithmic; it's the "why" behind M1.

**First test**: pin λ_max(T·R) → 1 as N → ∞ for sphere; bounded
for slab/cyl. Already in `diag_specular_overshoot_02_TR_spectrum.py`;
formalize as a theorem statement in Sphinx.

**Status**: STRONG diagnostic, MEDIUM algorithmic. Subsumes the
"why does sphere fail and slab succeed" question.

### M4 — Feynman-Kac surface Markov chain (A.4)

**Trigger**: `Σ_k (T·R)^k` is the Neumann series of a Markov chain on
the surface. Each step = one bounce. Geometric distribution at success
probability `1 - e^(-σ·2Rµ)`.

**What it gives us**: an explicit MC verification reference at any
geometry, any multi-region, with no analytic re-derivation. Specular
preserves µ — that's why the chain factorizes by µ and the closed
form `1/(1-e^(-σ·2Rµ))` exists.

**Implementation**: extend `diag_specular_overshoot_05` to multi-region
sphere/cyl. Compare matrix-Galerkin (overshoots) vs HS quadrature M1
(should match MC).

**Status**: STRONG as verification reference. The MC is the ground
truth; HS quadrature M1 is the production reference.

## Cross-method pollination (use selectively)

- **From CP / Hébert white BC**: `1/(1-P_ss)` IS the rank-1 reduction
  of M1's µ-resolved closure. Verify symbolically that
  `K_bc^mb |_(rank-1) = compute_P_esc · 1/(1-P_ss) · compute_P_ss`.
  This pins M1 against an existing closed form at low-rank.

- **From QMC / Koksma-Hlawka**: at multi-region τ with kinks at impact
  parameters, plain GL Q=64 has degraded algebraic convergence. Either
  use **piecewise-GL with knots at impact parameters** (kink-aware
  subdivision), or **Sobol'/scrambled Sobol'** with O(Q^{-1}(log Q)^d)
  convergence. Decision: piecewise-GL is preferred (deterministic,
  cheaper); use Sobol' as cross-check.

- **From cylinder Knyazev**: WARNING — cyl Phase 5 cannot just multiply
  f_mb into the Knyazev expansion. The multi-bounce factor depends on
  `sin θ_p` (3-D chord), introducing `∫ sin θ_p / (1 - e^(-σ·2R cos α/sin θ_p)) dθ_p`
  which is NOT a finite sum of Ki_(2+k) functions. Cyl Phase 5 needs
  either a new "multi-bounce Bickley-Naylor" family, OR adaptive
  2-D (α, θ_p) quadrature with subdivision near the joint grazing
  limit (α=π/2, θ_p=0). LITERATURE PULL needed before implementation.

- **From SN PN closure**: precompute `a_k = ∫ P̃_k(µ) f_mb(µ) dµ` (well-defined
  integrals; no singularity). Then `K_bc^mb = Σ_k a_k G M_k P` with
  fixed convolution matrices M_k. No matrix inverse, but more
  algebra. Equivalent to M1 with µ-quadrature replaced by polynomial
  expansion. Slower convergence than direct quadrature unless f_mb
  expansion converges fast.

## Refuted matches (negative results)

- **Mittag-Leffler partial-fraction expansion of `1/(1-e^(-x))`**: gives
  `1/x + 1/2 + Σ 2x/(x² + 4π²k²)`. Substituting into integrand yields
  a slow-convergent (1/k²) series of integrals. Right function class,
  WRONG decomposition. Use generating function (M2) instead.

- **Wiener-Hopf**: kernel is multiplication, not convolution. WH adds
  nothing to a diagonal-in-µ operator.

- **Group theory / Marshak vs Lambert basis switch**: the µ-Jacobian
  is real, but switching basis doesn't bound an unbounded multiplication
  operator. Insight is correct (one µ inside, one µ outside the
  measure); fix is M1, not basis change.

- **H-matrix / hierarchical low-rank in (i,j)**: M1 already gives
  O(I·J·Q) with rank-Q in (i,j); no need for hierarchical
  compression at typical I, J ~ 10²-10³.

## Concrete Phase 5 design recommendations

1. **Production reference: M1 (Hilbert-Schmidt separable quadrature)**.
   Q=64 GL on µ for homogeneous; Q=32 piecewise-GL with impact-parameter
   knots for multi-region.

2. **Verification cross-check: M2 (bounce-resolved sum)** at K_max=10
   for homogeneous; M4 (MC surface chain) for multi-region truth.

3. **Diagnostic statement: M3 (spectral theorem)** in Sphinx — gives
   the per-geometry "why" in one named theorem.

4. **Cyl is open**. Knyazev structure does NOT compose with f_mb cleanly.
   Literature-researcher needs a pull on Sanchez/Stamm'ler 3-D
   multi-bounce continuous-µ before cyl Phase 5 is designable.

5. **Slab Phase 5 is OPTIONAL** — slab MB already converges with the
   matrix form; M1 reformulation buys little.

## What this changes in the trigger table

Two trigger entries need to be sharpened (do on next revision):

- **A.7 Fredholm** — add "kernel of form a(r_i, µ)·b(r_j, µ)·c(µ)
  integrated over µ → Hilbert-Schmidt separable; rank-in-µ is the
  finite-dim axis, NOT rank in (r_i, r_j)". This is the M1 lever.

- **A.3 Spectral theory** — add "discrete approximations of unbounded
  multiplication operators do not converge in operator norm; check
  ess_range(φ) ∋ 1 before forming `(I - M_φ)^{-1}` matrix Galerkin".
  This is the M3 lever.

## Files referenced

- `.claude/agent-memory/numerics-investigator/specular_mb_overshoot_root_cause.md`
- `.claude/agent-memory/numerics-investigator/specular_mb_phase4_cyl_slab.md`
- `.claude/plans/specular-bc-phase4-multibounce-rollout.md` §10
- `orpheus/derivations/peierls_geometry.py` lines 2032-2400
- `derivations/diagnostics/diag_specular_overshoot_*.py` (12 numbered)
- `derivations/diagnostics/diag_specular_mb_phase4_*.py` (9 numbered)

---

**Promoted to skill on 2026-04-30** — see `cross-domain-frames`
reference.md Parts A.7 (Hilbert-Schmidt / separable kernels)
and A.3 (Spectral theory of multiplication operators)
sharpenings, plus new precedent file
`scripts/validated_hilbert_schmidt_separable.md` capturing
the M1/M2/M4 verification triad. This memory file remains as
evidence / precedent until Phase 5 production lands.
