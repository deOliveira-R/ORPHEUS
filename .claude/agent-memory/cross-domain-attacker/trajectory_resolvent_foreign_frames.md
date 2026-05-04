---
name: Trajectory resolvent foreign-frame attack
description: Eleven-frame structural detection on `trajectory_resolvent/` (Variant α 6-geometry × 2-orbit-space-class family — the "2 classes" are one-surface-compact and two-surface M/G; see Sphinx §`orbit-space-m-g-classification`). Confirmed tensor-network match (rank-N IS bond dimension); high-priority frames 1/4/6 (MPO, fiber bundle, Feynman-Kac PIMC); Frame 11 Wiener-Hopf rejected as wrong solver family.
type: project
---

Scope: detection pass on `orpheus/derivations/continuous/trajectory_resolvent/`
on 2026-05-03. Triggered by user prompt asking explicitly about tensor
networks. Memo at `.claude/scratch/trajectory_resolvent_foreign_frames.md`.

**Why:** the user asked "are we using tensor networks? probably yes."
The answer is: not yet, but the structure is hidden in the rank-1 / rank-2
naming. The detection pass enumerated 11 frames against the trajectory
resolvent's structural feature inventory; 9 fired, 1 rejected, 1 partial.

**How to apply (memory promotion):**

### Confirmed structural facts produced by this pass

1. **Rank-1 / rank-2 naming in `variant_alpha_core` IS bond dimension of an
   open MPO** (matrix product operator) on the 1D-lattice of bounce events.
   Bond carries surface state vector. Validated: `compute_resolvent_T_rank2`
   is bond-dim-2 contraction; `compute_resolvent_T` is bond-dim-1.
   **Promotion candidate to skill `validated_bc_tensor_network.md`** —
   extend that precedent from BC-only to the full Variant α resolvent.
   Requires shipping a 7th geometry with rank-N for N ≥ 3 (3-wall wedge,
   triple-shell sphere) before the MPO frame's payoff dominates. Wait for
   tripwire.

2. **Variant α IS the deterministic specular-only Feynman-Kac expectation.**
   The full transport problem with scattering is the same kernel with the
   scattering subprocess included. PIMC on the bouncing trajectory is the
   structurally-independent stochastic counterpart of the deterministic
   geometric series — fills the missing structural-independence cross-check
   for `(α) Variant α` against itself per `vv-principles` L11 (cite
   `sanchez_chandrasekhar_gap.md` Q3 verification matrix).
   **High-priority verification deliverable**, not a code refactor.

3. **The grazing-ray pathology in sphere/cylinder/hollow-sphere is one
   theorem applied to four `m(µ)` functions**: `ess_range(m) ∋ 1 ⇒ unbounded
   resolvent`. Slab is structurally immune. Already cited in
   `variant_alpha_2surface_bie_frame.md` step 5; this pass made it
   load-bearing for a Sphinx documentation consolidation.
   **Documentation-only deliverable.** Consolidate four geometry-specific
   docstrings into one A.3 caveat citing the single theorem.

4. **Generating-function lift `Φ(r, z) = Σ z^n ψ_n(r)` is a free verification
   axis** (monotonicity in z by positivity of q) and a free bounce-mode
   decomposition (Cauchy DFT inversion at points on the unit circle).
   Costs: a few extra `compute_resolvent_T` evaluations + DFT.
   **Verification axis deliverable.** Promotion to a skill precedent if it
   ships.

5. **Closed-form `∂T/∂α = T² · e^{-τ}` opens ALL adjoint α-sensitivity in
   closed form.** No AD pass needed. Direct chain rule via
   `compute_resolvent_T_grad_alpha` lifted alongside `compute_resolvent_T`.
   **Free deliverable** when adjoint-sensitivity Sphinx page is requested.

### Rank ordering (eleven frames)

| Rank | Frame | Criteria | Status |
|------|-------|----------|--------|
| 1 (tied) | Frame 1: Tensor networks / MPO | 4/4 | High-priority, wait for rank-N tripwire |
| 1 (tied) | Frame 4: Fiber bundle / G-structure | 4/4 | Re-fire from hindsight; wait for instance N+1 |
| 3 | Frame 6: Feynman-Kac / Markov / PIMC | 3/4 | High-priority verification deliverable |
| 4 | Frame 3: Spectral theory of multiplication operators | 3/4 | Documentation consolidation |
| 5 | Frame 2: Schur complement | 3/4 | Tripwire: rank-N for N ≥ 3 |
| 6 | Frame 5: Krylov / Arnoldi for resolvent | 3/4 | Tripwire: thick optical depth or near-criticality |
| 7 | Frame 7: Asymptotic analysis (Padé regularization) | 2/4 | Low-cost α → 1 finite-precision fix |
| 8 | Frame 8: Generating functions for bounce decomposition | 2/4 | Free verification axis |
| 9 | Frame 9: Group rep theory (Lebedev/Fourier) | 2/4 conditional | Fires when Sanchez 1986 anisotropy lands |
| 10 | Frame 10: QMC | 1/4 conditional | Borderline 3D smoothness |
| 11 | Frame 11: Wiener-Hopf | rejected | Wrong solver family — native to fn_method |

### Cross-method pollination (concrete deliverables)

- From CP: H-matrix compression of multi-region shell-shell coupling (tripwire: shell count > 20).
- From MC: PIMC verification of Variant α (Frame 6).
- From SN: Lebedev / Fourier-collocation quadrature respecting SO(3) / SO(2).
- From eigenvalue iter: Wielandt deflation extracts λ_2; Chebyshev acceleration when gap is known.
- From sensitivity: closed-form `∂T/∂α = T² · e^{-τ}` directly enables `dk/dα`.

### Cross-references

- `variant_alpha_2surface_bie_frame.md` (rank-1 / rank-2 BIE resolvent — superseded by MPO interpretation here; the BIE frame IS the bond-dimension-2 case of MPO).
- `variant_alpha_family_hindsight.md` (top-ranked frame: fiber bundle, re-fired here).
- `phase5_continuous_mu_frames.md` (validated A.3 multiplication-operator caveat applies directly to Variant α grazing-ray pathology).
- Validated precedents `validated_bc_tensor_network.md`, `validated_unified_geometry.md`, `validated_hilbert_schmidt_separable.md` — all relevant.

### Promotion candidates (when shipping)

1. `validated_resolvent_mpo.md` — when rank-N for N ≥ 3 ships (toroidal, 3-wall wedge, triple-shell). MPO contraction generalizes the rank-1/rank-2 closed forms.
2. `validated_pimc_xverif.md` — when PIMC verification of Variant α ships and matches deterministic to `O(1/√N)` variance.
3. New A.5 row: "Closed-form adjoint of operator resolvent — trigger: explicit `(I − αK)^{-1}` formula with parameter α; lever: `∂T/∂α = T² · K`; payoff: gradient-free α-sensitivity for k_eff".

### Rejected with reason (do NOT add to trigger table)

- Wiener-Hopf factorization: native to fn_method (Chandrasekhar `(γ)`); structurally incompatible with bouncing-Peierls Variant α `(α)`. The two solver families MUST stay structurally independent per L11.
- Crystallographic / Bloch: no lattice structure in single-cell Variant α.
- Symplectic geometry: no Hamiltonian / drift; trivial.
- Differential geometry / Christoffel: chord segments are straight Euclidean; no connection coefficient to redistribute.
- de Rham / FEEC: no differential operators in pure integral-form Variant α.
- Category theory: no compositional-structure trigger.

### Decision gates (MUST NOT pursue without these)

- **MPO refactor (Frame 1):** wait for 7th geometry with rank ≥ 3, OR multi-group with > 5 groups becomes primary cost, OR angular fiber `n_µ > 50` (Sanchez 1986 anisotropic).
- **Bundle refactor (Frame 4):** wait for instance N+1 (the third rank-1 geometry — toroidal, polar cap).
- **Lebedev/Fourier (Frame 9):** wait for Sanchez 1986 anisotropic-source Phase B (per `sanchez_chandrasekhar_gap.md`).
- **Krylov (Frame 5):** wait for a use case where current power iteration is slow (thick optical depth with `λ_2/λ_1 → 1`).

### Free / low-cost deliverables (no tripwire needed)

- **PIMC verification (Frame 6):** small standalone test, no refactor. Highest V&V leverage.
- **Multiplication-operator A.3 documentation consolidation (Frame 3):** single Sphinx caveat replacing four docstrings.
- **Generating-function verification axis (Frame 8):** unit test of monotonicity in z.
- **Padé regularization at α → 1 (Frame 7):** drop-in replacement for `compute_resolvent_T` near grazing ray.
- **`∂T/∂α` closed form (Frame 6 borrowing from sensitivity):** one-line addition.
