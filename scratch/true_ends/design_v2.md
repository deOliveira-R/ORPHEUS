# True ends: design v2 (2026-11-11)

This revision folds in the rulings of 2026-11-11 (R6, R7, R4 re-ruled), the W5 reviews (`review/xda.md`, `review/elegance.md`) and the R6 census (`census/r6_density_sites.md`). Rulings R1–R7 are in `.claude/plans/true_ends.md`.

## Settled objects

1. **State side.** V and its member classes: the flux classes, plus the importance ψ† as a flux-classed member.
2. **Data side.** V\* = `V.dual()`, stored integrated (R6). The seven `*SourceSink` classes become V\* members beside `Covector`.
   - `[M]` `Covector` is already stored integrated, as G·w (`numerics/space.py:865`, census finding 1). So D1's chart question dissolves: V\* has one stored chart, the integrated one.
   - The density reading is a VIEW (`♯`), queried, never stored.
3. **Operators.**
   - Flux → source atoms are `V → V*`, with no metric legs.
   - The inverses are `V* → V`.
   - The transport adjoint is `.dual()`, metric-free, and `(Aᵀ)ᵀ = A` holds by algebra.
   - `.H`, the Hilbert adjoint, stays for same-side operators (2 of 18 production sites) and refuses on `V → V*`.
   - `[M]` All metric legs are built in one place, `AdjointOperator.__init__` (`operator.py:1672/1678`).
4. **The moment head (R7).**
   - C is the coefficient space, read by the reconstruction R. C\* is the moment space, written by the analysis M. They are linked by the continuum ℓ-diagonal map.
   - `HarmonicMomentFlux` ∈ C\* ⊗ cells.
   - Λ and the fission moment operator are spelled through the tensor algebra on their true ends, with the ends being the product of the factor spaces (a new law on `TensorProductOperator`).
5. **The S_N boundary.** The state is ψ on the non-tangential ordinates (the quotient), with the current view ♭ψ. Reflection and albedo are column-stochastic maps on currents.
6. **Member seam (R5).** `FlattenedOperator`'s template retires outright, and its members come from `domain.member()` / `codomain.member()`.
7. **Retired:**
   - `PredictedKernel.H`'s hand-rolled Riesz legs;
   - the metric legs on flux → source adjoints (15 of 18 `.H` sites move to `.dual()`; that split is `[R]`);
   - the `/W` producer convention at its main home (`angular_lift.py:440`, `[M]` 3.05 M of 6.0 M value-creating mints).

## Defects the census found, fixed in the carve

- **Diffusion k integrates the loss output over volume a second time** (`diffusion/solver.py:431`). Under R6 that either raises or moves k.
- **The adjoint eigenvalue gauge** is a plain coordinate sum of F†ψ\* (`iteration.py:1590`). It is restated as the canonical pairing, which is the global rescale of the importance.
- **`KEigenvalue.compute_production_rate` and `rayleigh(w=1)`** are restated as pairings ⟨covector, flux⟩.

## Open questions for the user

- **Q-A, the angular factor.** R6 states the integrated units as n/(s·sr) per cell, which is integrated in space and still per steradian in angle. That is the discrete per-ordinate balance, since the transport operator is pointwise in angle.
  - Then the data side is dual on the SPATIAL factor only: V\* = A ⊗ X\*. The sweep needs no ÷w_n, so zero-weight and negative-weight ordinates stay harmless (D6, D7).
  - The angular measure enters only where angles couple: in the scattering through M and R, and in the pairing.
  - The alternative is dual on both factors (n/s per cell-ordinate). That puts ÷w_n in the sweep, which is undefined on w_n = 0 and flips sign on negative weights.
  - `[R]` The spatial-only form is native: it is the equation the sweep solves.
- **Q-B, face area.** `[M]` The trace metric is |Ω·n|·w_n with NO face area, while a current carries the area, and `_trace_flux_weights` already includes A_f. On curvilinear and 2-D traces the areas differ per face, so either the trace measure gains the area (a true current) or the metric is a per-unit-area current. `[R]` the area belongs in the measure.
- **Q-C, user inputs.** The external source, the diffusion source and the detector response are specified by users. Physicists write them as densities (n/(cm³·s)). Recommendation: accept densities and convert ONCE at the input boundary (♭ at the parse edge, `coding-standards` "an input boundary"). `phi_fixed/*` in `carry_head.npz` tests whichever convention is chosen.
- **Q-D: RESOLVED, no conflict** (`ld_moments/memo.md`). The stored value is `q_int = M_cell·q_coef` (V·θ^{|o|} on Cartesian cells). The design rows this creates `[R]`:
  - the 1-D `_slope_fold` drops θ from its source half once sources arrive integrated (θ stays in D₂');
  - S, N₂ₙ and F emit mass-weighted moments (`[M]` 0 of 12 files in `transport/operators/` carry θ today);
  - the multi-D UBLD DAG loses its `/mass_diag` in the residual and its transpose, and gains the ×V of census C3, both in one step;
  - the moment-resolved arm of the external-source door (`sn/solver.py:3143`) follows R10 and applies `M_cell` at the input boundary.

## Sequencing `[R]`

1. **Step A:** the data-side classes and chart (R6), and the operator ends `V → V*` / `V* → V`, with `.dual()` as the transport adjoint and the defects above. This is the largest step.
   - `[M]` About 780 gate node ids read source values.
   - Bit-identity on k and fluxes is the acceptance.
   - Source pins re-baseline with the measure factor as their reason (51 snapshot rows in 10 files, plus about 7 npz keys `[R]`).
   - 100 hand-built density references in 33 files are rewritten in integrated form.
   - 314 pairing gates in 51 files are re-posed.
2. **Step B:** the S_N boundary quotient and the current view (R4).
3. **Step C:** the moment head C/C\* (R7), with Λ and fission on their true ends.
4. **Step D:** #705 rebases on top, for the remaining degenerate primal spaces.

Each step goes through the test-architect first (red gates), then the main agent writes, then qa and elegance review.
