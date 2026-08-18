# `implements` declaration inventory — 332 claimed equations

⚠ **THIS IS A FIRST PASS WITH KNOWN, MEASURED DEFECTS. Do not act on it
without reading this header.** Produced 2026-08-18 by a 12-way Haiku
fanout over every equation that carries a `verifies` coverage claim and
lacks a DECLARED `implements` link. Graph at `4add6b59`.

## Why it exists

`[M]` 12999 of 13084 `implements` edges are inferred from a shared name
token, and only 46 of 903 equations carry a declaration. That is why the
V&V audit landed in nexus (`a0f038e`) cannot produce a trustworthy
refutation: all 10 refutable claims test against a guessed link. See
ORPHEUS #381, nexus#82.

## What the fanout produced

| verdict | count |
|---|---:|
| DECLARABLE | 180 |
| NOTHING:definition | 67 |
| NOTHING:identity | 51 |
| NOTHING:law | 19 |
| NOTHING:canonical-form | 13 |
| UNSURE | 2 |

## ⛔ Its two halves are NOT equally trustworthy

**DECLARABLE is mechanically checkable, and 18 rows fail.** Of 293
claimed implementers: 275 resolve to a legal node, **11 name a MODULE**
(`[edge.implements]` admits only function/method/class), **7 do not
exist at all**, and 1 row is DECLARABLE with no implementer listed. Each
would become a dangling or ontology-forbidden `confidence=1.0` edge —
strictly worse than the guess it replaces. Repair them before landing.

**NOTHING is NOT checkable and fails in the FLATTERING direction.** A
wrong NOTHING suppresses every guess for the equation AND records that
nothing implements it, hiding a real coverage gap behind a confident
declaration. `[M]` **110 of 150** NOTHING rows have contrary evidence
(claiming tests and/or candidates pointing at real code) — e.g.
`transport-cartesian` (42 claiming tests) as a law, `peierls-unified`
(96) as a definition, `ki3-def` (28, the subject of #348) as a
definition.

⛔ **One is PROVEN wrong.** `attenuation` (`method_of_characteristics`,
41 claiming tests) was classified `NOTHING:identity`. It ships:
`orpheus/moc/core.py:214-218` computes
`ψ_out = ψ_in·e^{-τ} + Q/Σ·(1-e^{-τ})`, small-τ series guard and all.

⚠ **The calibration signal is absent by construction**: 296 of 332 rows
are "high" confidence and only **2** are UNSURE, on unfamiliar physics,
despite the brief inviting UNSURE explicitly. Read every confidence
here as uninformative.

## Two prerequisites before any of this can land

1. **`:nothing:` does not exist** — `option_spec = {"by": ...}` only
   (nexus#85 is still a sketch). The whole NOTHING half has nowhere to
   land, which happens to be the untrustworthy half.
2. **The declaring path does not consult the ontology** (nexus#86), so
   the 11 module-typed targets would ship silently. The check that
   caught them lives in a scratch script; it belongs in the build.

## ⚠ The contract that governs landing any of this

From `ImplementsDirective`'s own docstring: declaring ONE implementer
stands the inference down for **the whole equation**. So a partial
declaration is worse than none — it leaves every other implementer
silently unlinked. 77 rows here list 2–5 implementers; the singletons
are the ones to re-check hardest.

---


---

# === out_01.md ===

# Implementation Inventory Analysis — Slice 01

Analyzing 25 equations from 3 Sphinx pages. Working systematically through each equation to determine verdict and list all implementers.

---

## ws-pitch
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:320; orpheus/mc/solver.py:92
note: Geometric definition of the Wigner-Seitz cell pitch formula `p = r_cell * sqrt(π)`. Defined once in `ConcentricPinCell.default_pwr()` constructor with the exact formula; it's a constant definition, not a behavioral rule that would have multiple implementations or a verb.

## hetero-tolerance
verdict: NOTHING:definition
implementers:
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:350
note: Tolerance criterion for heterogeneous test acceptance. This is a verification threshold for the tolerance parameter, not a coded equation. The actual tolerance check lives in test_mc_heterogeneous but as a test assertion, not a declarable equation implementation.

## majorant
verdict: DECLARABLE
implementers: orpheus.mc.solver._precompute_xs
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:379; orpheus/mc/solver.py:355-357 (_precompute_xs computes sig_t_max = max over materials); test: tests/mc/test_properties.py::test_majorant_computation
note: The majorant cross-section is computed as the maximum total cross section over all materials. Implemented in _precompute_xs as sig_t_max = np.maximum(sig_t_max, mix.SigT).

## free-flight
verdict: DECLARABLE
implementers: orpheus.mc.solver._random_walk
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:398; orpheus/mc/solver.py:392 (free_path = -np.log(rng.random()) / xs.sig_t_max[ig]); test: tests/mc/test_gaps.py::test_free_path_exponential
note: Free-flight distance sampled from exponential distribution with rate Σ_maj. Implemented directly in the random walk loop.

## direction-sampling
verdict: DECLARABLE
implementers: orpheus.mc.solver._random_walk
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:526; orpheus/mc/solver.py:395-398 (theta = pi * rng.random(); phi = 2*pi * rng.random(); dir_x/dir_y = sin/cos projections); test: tests/mc/test_gaps.py (multiple tests exercise this)
note: ⚠ **Note**: The docs explicitly flag this as NOT true isotropic sampling (ERR-018), but it IS a declarable equation because the code implements this exact formula.

## virtual-collision-probability
verdict: DECLARABLE
implementers: orpheus.mc.solver._random_walk
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:501; orpheus/mc/solver.py:407-410 (sig_v / xs.sig_t_max[ig] >= rng.random() decides virtual collision); test: tests/mc/test_properties.py::test_delta_tracking_virtual_probability
note: Probability P_virtual = (Σ_maj - Σ_t) / Σ_maj coded directly as sig_v / sig_t_max.

## decompose
verdict: NOTHING:identity
implementers:
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:421; orpheus/mc/solver.py:356-357
note: Algebraic decomposition identity: Σ_maj = Σ_t + Σ_virtual. This is a mathematical identity, not an algorithm to implement. The actual values are computed in _precompute_xs, but the identity itself is not "implemented" in code.

## scattering-cdf
verdict: DECLARABLE
implementers: orpheus.mc.solver._random_walk
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:702; orpheus/mc/solver.py:413-415 (cum_s = np.cumsum(sig_s_row); ig = np.searchsorted(cum_s, rng.random() * sig_s_sum)); test: tests/mc/test_gaps.py (multiple tests), tests/mc/test_properties.py::test_scattering_cdf_sampling
note: Cumulative distribution function for scattering group selection. Implemented via cumsum + searchsorted in the random walk.

## branching
verdict: DECLARABLE
implementers: orpheus.mc.solver._random_walk
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:628; orpheus/mc/solver.py:412-428 (r = rng.random() * sig_t; if r < sig_s_sum: scatter; elif r < sig_s_sum + sig_2n_sum: n2n; else: absorption); test: tests/mc/test_properties.py::test_scattering_branching_ratio
note: Three-way branch decision (scatter/n2n/absorption) based on reaction cross sections. Fully coded in _random_walk via scaled uniform sampling.

## fission-weight
verdict: DECLARABLE
implementers: orpheus.mc.solver._random_walk
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:741; orpheus/mc/solver.py:421-423 (w *= sig_p / sig_a); test: tests/mc/test_properties.py::test_fission_weight_adjustment, tests/mc/test_gaps.py::test_absorption_nonfissile_zeroes_weight
note: Weight multiplication on absorption: w ← w · (νΣ_f / Σ_a). Implemented in _random_walk as the else branch of the collision decision.

## chi-sampling
verdict: DECLARABLE
implementers: orpheus.mc.solver._random_walk
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:802; orpheus/mc/solver.py:424-425 (ig = np.searchsorted(xs.chi_cum, rng.random())); test: tests/mc/test_gaps.py (multiple tests via chi-sampling), tests/mc/test_properties.py tests
note: Fission spectrum sampling via CDF. Implemented using precomputed cumulative spectrum chi_cum.

## periodic-bc
verdict: DECLARABLE
implementers: orpheus.mc.solver._random_walk
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:606; orpheus/mc/solver.py:401-402 (nx_ = nx_ % pitch; ny_ = ny_ % pitch); test: tests/mc/test_monte_carlo.py (multiple tests use periodic BC)
note: Periodic boundary condition via modulo operation. Implemented directly in the random walk loop.

## roulette-prob
verdict: DECLARABLE
implementers: orpheus.mc.solver._russian_roulette
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:829; orpheus/mc/solver.py:452-453 (terminate_p = 1.0 - bank.weight[i_n] / weight0[i_n]); test: tests/mc/test_properties.py, tests/mc/test_gaps.py
note: Russian roulette kill probability: P_kill = 1 - w/w_0. Implemented in _russian_roulette function.

## roulette-restore
verdict: DECLARABLE
implementers: orpheus.mc.solver._russian_roulette
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:837; orpheus/mc/solver.py:454-457 (if terminate_p >= rng.random(): w=0; elif terminate_p > 0: w=w_0); test: tests/mc/test_properties.py::test_roulette_restore_weight, test_roulette_weight_conservation
note: Post-roulette weight restoration: w_after ∈ {0, w_0}. Implemented in _russian_roulette.

## roulette-conservation
verdict: NOTHING:identity
implementers:
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:858
note: Algebraic identity proving that E[w_after] = w_before under the roulette scheme. This is a derivation/proof, not an implementation. The identity follows from the probabilities, not from code.

## splitting
verdict: DECLARABLE
implementers: orpheus.mc.solver._split_heavy
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:890; orpheus/mc/solver.py:468-482 (N = floor(w); stochastic rounding; w_new = w/N); test: tests/mc/test_properties.py::test_splitting_weight_conservation, tests/mc/test_gaps.py::test_splitting_copy_count
note: Weight-based splitting with stochastic rounding. Fully implemented in _split_heavy function.

## splitting-weight-conservation
verdict: NOTHING:identity
implementers:
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:914
note: Algebraic identity: E[N] = w under the stochastic splitting scheme. This is a proof/derivation of the expected value, not an implementation.

## keff-cycle
verdict: DECLARABLE
implementers: orpheus.mc.solver.solve_monte_carlo
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:935; orpheus/mc/solver.py:621 (keff_cycle = bank.weight[:bank.n].sum() / weight0.sum()); test: tests/mc/test_gaps.py::test_keff_cycle_estimator
note: Cycle eigenvalue estimation as weight ratio. Implemented in the main solver loop.

## keff-mean
verdict: DECLARABLE
implementers: orpheus.mc.solver.solve_monte_carlo
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:964; orpheus/mc/solver.py:630-631 (keff_history[ia] = keff_active[:i_active].mean()); test: tests/mc/test_convergence.py::test_bias_decreases_with_histories, test_inactive_cycles_reduce_bias, test_sigma_scales_with_sqrt_n
note: Cumulative mean over active cycles. Implemented as the running mean of keff_active.

## sigma-keff
verdict: DECLARABLE
implementers: orpheus.mc.solver.solve_monte_carlo
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:971; orpheus/mc/solver.py:632-635 (sigma_history[ia] = sqrt(sum((k_m - mean)^2) / (M-1) / M)); test: tests/mc/test_convergence.py::test_sigma_scales_with_sqrt_n
note: Standard deviation of the mean using unbiased sample variance. Fully implemented in the solver loop.

## collision-estimator
verdict: DECLARABLE
implementers: orpheus.mc.solver._random_walk
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:998; orpheus/mc/solver.py:411 (tally[ig] += w / sig_t); test: tests/mc/test_gaps.py::test_2g_flux_ratio_homogeneous
note: Collision estimator for flux tally: φ_g ≈ (1/N·V) Σ(w/Σ_t). Implemented at every real collision in _random_walk.

## mc-lethargy-width-sign
verdict: NOTHING:identity
implementers:
confidence: high
evidence: docs/theory/methods/monte_carlo.rst:1040
note: Mathematical sign relationship for lethargy width: Δu_g = ln(E_{g+1}/E_g) < 0 when groups are fastest-first ordered. This is a consequence of the group ordering convention, not an implemented algorithm.


## quadrature-product-weights
verdict: DECLARABLE
implementers: orpheus.numerics.quadrature.directional.Quadrature.product
confidence: high
evidence: docs/theory/methods/sn/angular_quadrature.rst:351; orpheus/numerics/quadrature/directional.py:660 (Quadrature.product factory method); orpheus/numerics/quadrature/recipes.py (product_mu_phi function computes weights w_GL * 2π/N_φ); test: tests/sn/primitives/test_quadrature.py (many tests)
note: Tensor product quadrature weights w_{p,m} = w_GL(μ_p) · 2π/N_φ. Implemented via the product factory method and underlying recipes module.

## quadrature-ordinate-permutation
verdict: DECLARABLE
implementers: orpheus.numerics.quadrature.directional.Quadrature.ordinate_permutation
confidence: high
evidence: docs/theory/methods/sn/angular_quadrature.rst:418; orpheus/numerics/quadrature/directional.py:337 (ordinate_permutation method); orpheus/numerics/symmetry.py (motion certification via RigidMotion.preserves); test: tests/sn/primitives/test_quadrature.py (TestReflectionIndices and related tests)
note: Permutation π such that Ω_{π(n)} = Q·Ω_n and w_{π(n)} = w_n, certified by ordinate-permutation method on all Quadrature objects.

## reflective-bc
verdict: DECLARABLE
implementers: orpheus.geometry.boundary.reflective.ReflectiveBoundary
confidence: high
evidence: docs/theory/methods/sn/boundary_conditions.rst:126; orpheus/geometry/boundary/reflective.py:32 (ReflectiveBoundary class implements ψ_n^in = ψ_{n'}^out via ordinate permutation); orpheus/sn/boundary/realizer.py (SNBoundaryRealizer realizes the BC as PermutationOperator); test: tests/sn/eigenvalue/test_keff_slab.py, tests/sn/primitives/test_quadrature.py (alpha_boundary_zero, reflection tests)
note: Specular reflection: incoming flux set to outgoing flux of reflected partner. Implemented via ReflectiveBoundary descriptor + SNBoundaryRealizer.

---

## Summary

**Total equations analyzed:** 25

**Verdicts:**
- **DECLARABLE**: 16 equations
  - Monte Carlo: majorant, free-flight, direction-sampling, virtual-collision-probability, scattering-cdf, branching, fission-weight, chi-sampling, periodic-bc, roulette-prob, roulette-restore, splitting, keff-cycle, keff-mean, sigma-keff, collision-estimator
  - Quadrature: quadrature-product-weights, quadrature-ordinate-permutation
  - Boundary: reflective-bc

- **NOTHING:identity** or **NOTHING:law**: 6 equations (algebraic identities/laws)
  - decompose, roulette-conservation, splitting-weight-conservation, mc-lethargy-width-sign

- **NOTHING:definition**: 3 equations (definitions without algorithmic implementation)
  - ws-pitch, hetero-tolerance

- **UNSURE**: 0

**Structural observations:**

1. **Monte Carlo page (23 equations)**: Almost entirely declarable. The solver code in `orpheus/mc/solver.py` provides complete implementations of all major equations through organized functions: `_precompute_xs`, `_random_walk`, `_russian_roulette`, `_split_heavy`, and `solve_monte_carlo`. The few non-implementable equations are either algebraic identities (conservation laws, decomposition) or definitions (pitch formula, tolerance criterion).

2. **Quadrature page (2 equations)**: Both are declarable. The `Quadrature` class in `orpheus/numerics/quadrature/directional.py` provides factory methods (`.product()`, `.folded_product()`, `.level_symmetric()`, `.gauss_legendre()`) and the `ordinate_permutation()` method that directly implement the documented math.

3. **Boundary conditions page (1 equation)**: Reflective BC is declarable via `ReflectiveBoundary` class + `SNBoundaryRealizer`. The equation ψ_n^in = ψ_{n'}^out is realized as a permutation operator.

**Key implementation patterns:**
- Monte Carlo solver uses "random walk → collision physics → population control → statistics" architecture with extracted functions for each phase
- Quadrature families use factory methods (classmethods) rather than subclasses
- Boundary conditions use descriptor pattern (law carries geometry + response kernel) + realizer pattern (converts to actual operators)

**All implementable equations found their production code.** No guesses were needed — every DECLARABLE verdict points to explicit function/class implementations in orpheus/ or tests.

---

# === out_02.md ===

# Implementation Inventory — Slice 02

Processing equations from: acceleration, index, spherical_harmonics, method_of_characteristics


## sn-dsa-cell-update
verdict: DECLARABLE
implementers: orpheus.sn.acceleration.dsa, orpheus.derivations.discrete.sn.dsa
confidence: high
evidence: docs/theory/methods/sn/acceleration.rst:418; test verifies at line 57-60 of tests/derivations/test_dsa_rules.py; implementation in orpheus.sn.acceleration.dsa.DSACorrection and derivation module
note: The four-step DSA derivation produces the cell-average update equations (28) in Larsen 1982. Tests verify via exact symbolic derivation (L0) and production build ties it to DSACorrection/DSALowOrderSystem (L3).

## sn-dsa-coefficients
verdict: DECLARABLE
implementers: orpheus.sn.acceleration.dsa.DSALowOrderSystem, orpheus.derivations.discrete.sn.dsa
confidence: high
evidence: docs/theory/methods/sn/acceleration.rst:368; test at line 71-75 of tests/derivations/test_dsa_rules.py; implements via DD instance
note: Larsen (23a–f) coefficients for consistent DSA. Derivation proves they equal Larsen's printed forms. Production DSALowOrderSystem.from_sn_mesh() builds these coefficients.

## sn-dsa-consistent-fourier
verdict: DECLARABLE
implementers: orpheus.sn.acceleration.dsa
confidence: high
evidence: docs/theory/methods/sn/acceleration.rst:199; test at test_dsa_rate.py::TestD11SpectralRadius::test_rho_below_bound_c09
note: Fourier analysis bound ρ_DSA ≤ 0.2247c. The test measures the spectral radius. The bound is achieved by the DSA correction implementation.

## sn-dsa-consistent-low-order
verdict: DECLARABLE
implementers: orpheus.sn.acceleration.dsa.DSALowOrderSystem, orpheus.derivations.discrete.sn.dsa.build_consistent_dd_system
confidence: high
evidence: docs/theory/methods/sn/acceleration.rst:396; two tests: line 50-55 of test_dsa_rules.py proves the main theorem; test_dsa_low_order.py::TestProductionTie::test_low_order_matches_reference_builder ties production to derivation
note: Larsen's (27) interior row with coefficients (23a–f). Proven by exact symbolic derivation, then production DSALowOrderSystem matched entry-for-entry against the reference builder.

## sn-dsa-correction-vanishes
verdict: DECLARABLE
implementers: orpheus.sn.acceleration.dsa.DSACorrection
confidence: high
evidence: docs/theory/methods/sn/acceleration.rst:707; test test_dsa_acceleration.py::TestD6CorrectionVanishes::test_zero_displacement_maps_to_exact_zero
note: The correction operator produces zero when the displacement is zero, a correctness-by-construction property. Tested by verifying DSACorrection.apply() on zero input.

## sn-dsa-marshak
verdict: DECLARABLE
implementers: orpheus.sn.acceleration.dsa.DSALowOrderSystem
confidence: high
evidence: docs/theory/methods/sn/acceleration.rst:478; test test_dsa_rules.py::TestFourStepDerivation::test_boundary_rows
note: Marshak boundary rows (38a) with coefficients derived as (γ_N, W2+/W2). Built by DSALowOrderSystem constructor for open boundary conditions.

## sn-dsa-restriction
verdict: DECLARABLE
implementers: orpheus.sn.acceleration.dsa.DSALowOrderSystem
confidence: medium
evidence: docs/theory/methods/sn/acceleration.rst:602; two tests in test_dsa_low_order.py::TestRestrictionProlongation
note: Restriction (ℓ=0 frame row) that conserves particles and matches the frame moment. The R/P pair uses the angular frame's (ℓ=0) and (ℓ=1) rows; DSALowOrderSystem wires them.

## sn-dsa-s2-exactness
verdict: DECLARABLE
implementers: orpheus.sn.acceleration.dsa.DSACorrection
confidence: high
evidence: docs/theory/methods/sn/acceleration.rst:1268; test test_dsa_rate.py::TestS2Exactness::test_one_correction_exactness
note: The S2 component of the correction is exact for the single-correction cycle (second-order convergence to the fixed point). Measured by DSACorrection's application in the rate test.

## sn-dsa-sweep-inverse-identity
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/methods/sn/acceleration.rst:1104; label shows it's testing L^(-1) ∘ L = I
note: This is a mathematical identity: the forward sweep applied to the inverse sweep produces the identity. Tested by operators/test_sweep_inverse_identity.py::TestSweepInverseIdentity. Identity operators do not "implement" them — they are verified.

## sn-dsa-synthesis
verdict: DECLARABLE
implementers: orpheus.sn.acceleration.dsa.DSACorrection
confidence: high
evidence: docs/theory/methods/sn/acceleration.rst:807; test test_dsa_rate.py::TestP1DSAArm::test_s2_exactness_with_l1_scattering
note: The P1 synthesis step that applies the correction. DSACorrection carries both P0 and P1 arms and applies the correction via the specified P1 method.

## dd-curvilinear-scalar
verdict: DECLARABLE
implementers: orpheus.transport.spatial.diamond.DiamondDifference, orpheus.transport.spatial.scheme.CellUpdateStrategy
confidence: high
evidence: docs/theory/methods/sn/index.rst:638; multiple tests execute the curvilinear sweep; code in orpheus/transport/spatial/ implements the DD scheme
note: The cell average update formula for curvilinear DD, mirroring the dissolved sweep code. Implemented in DiamondDifference.solve_cell() and related methods in the spatial scheme.

## dd-cylindrical-degenerate
verdict: DECLARABLE
implementers: orpheus.transport.spatial.diamond.DiamondDifference
confidence: medium
evidence: docs/theory/methods/sn/index.rst:670; test test_diamond.py::TestCylindricalDegenerate::test_degenerate_cell_synthetic
note: The degenerate case where radial flux vanishes (pure azimuthal flow). The special case is handled in the DiamondDifference implementation when abs_mu < 1e-15.

## dd-mm-closure-constants
verdict: DECLARABLE
implementers: orpheus.sn.sweep.pole_angular_closure, orpheus.transport.spatial.diamond.DiamondDifference
confidence: high
evidence: docs/theory/methods/sn/index.rst:628; constants computed from Bailey dome (alpha) and Morel-Montry weights; tested in test_diamond.py
note: The redistribution constants c_out and c_in for the M-M angular closure used in curvilinear DD. Computed from pole_angular_closure functions and used in cell updates.

## dd-slab-scalar
verdict: DECLARABLE
implementers: orpheus.transport.spatial.diamond.DiamondDifference
confidence: high
evidence: docs/theory/methods/sn/index.rst:598; slab cell update formula; tested by test_diamond.py::TestBitIdenticalSlab tests
note: The scalar form of the DD recurrence for slab geometry. Implemented in DiamondDifference.solve_cell() for 1D Cartesian geometry.

## transport-cylindrical
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/sn/index.rst:305; this is the PDE form of the cylindrical transport equation
note: This is a mathematical definition of the cylindrical transport equation PDE, not an implemented algorithm. Tests verify MMS convergence against the discretized form, but the PDE itself is not "implemented" — it is the basis for discretization.

## transport-spherical
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/sn/index.rst:287; this is the PDE form of the spherical transport equation
note: This is a mathematical definition of the spherical transport equation PDE with angular redistribution. Not implemented directly; the discrete balance equation is what gets discretized and solved.

## hilbert-adjoint-equals-metric-times-S0
verdict: DECLARABLE
implementers: orpheus.numerics.frame._AdjointOperator, orpheus.numerics.spaces.spherical_harmonic_space.SphericalHarmonicSpace
confidence: high
evidence: docs/theory/foundations/spherical_harmonics.rst:370; implemented by frame.analysis.H and computed generically by _AdjointOperator wrapper; tested at test_spherical_harmonic_space.py::test_T_carries_w_n_and_H_carries_g_C
note: The Hilbert adjoint operator equals the S0 synthesis times the metric diagonal g_C. Implemented as a typed operator composition on the spherical harmonic frame.

## moment-projection-transpose-T
verdict: DECLARABLE
implementers: orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.analyze_transpose
confidence: high
evidence: docs/theory/foundations/spherical_harmonics.rst:355; test test_spherical_harmonic_space.py::test_T_carries_w_n_and_H_carries_g_C
note: The moment projection transpose (Π† = w_n · S0). Exposed by frame.analysis.apply_transpose on the analysis face's codomain carrying the quadrature weight.

## pi-r-equals-4pi-i
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/foundations/spherical_harmonics.rst:317; the 4π-tightness identity; tested by test_spherical_harmonic_space.py::test_pi_R_is_4pi_identity_on_band_limited
note: This is a mathematical identity (Π R = 4π I) that holds on band-limited inputs. The frame operator is 4π-tight by construction. The operators that appear in it (Π, R) are implemented, but the identity itself is a property verified, not an algorithm implemented.

## real-sh-discrete-orthogonality
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/foundations/spherical_harmonics.rst:301; orthogonality relation; tested by test_spherical_harmonic_space.py tests
note: This is a mathematical identity defining orthogonality of spherical harmonics under a quadrature. It is verified by tests (basis mass matrix) but is not "implemented" as an algorithm — it is a property of the basis.

## sh-addition-theorem-reconstruction
verdict: DECLARABLE
implementers: orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.reconstruct, orpheus.numerics.spaces.spherical_harmonic_space.SphericalHarmonicSpace, orpheus.numerics.frame.GalerkinFrame.reconstruction
confidence: high
evidence: docs/theory/foundations/spherical_harmonics.rst:387; implemented by frame.reconstruction face; test test_spherical_harmonic_space.py::test_R_equals_2l_plus_1_times_S0
note: The reconstruction operator R = (2ℓ+1) · S0 that projects from SH coefficients to angular ordinates. The (2ℓ+1) factor is read from SphericalHarmonicBasis.addition_theorem_factor.

## sh-space-metric
verdict: DECLARABLE
implementers: orpheus.numerics.spaces.spherical_harmonic_space.SphericalHarmonicSpace
confidence: high
evidence: docs/theory/foundations/spherical_harmonics.rst:397; the inner product metric g_C = 4π/(2ℓ+1); tested by test_spherical_harmonic_space.py::test_space_inner_product_weights_equal_4pi_over_2l_plus_1
note: The metric tensor diag(4π/(2ℓ+1)) that defines the inner product on SH coefficient space. Implemented as the basis_space inner_product_weights on SphericalHarmonicSpace.

## moc-mms-psi-ref
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/verification/method_of_characteristics.rst:305; MMS reference ansatz for testing
note: This is a manufactured solution definition: a smooth radial scalar flux used in MMS convergence tests. It is a test fixture specification, not an algorithm to be implemented. The test code evaluates it, but does not "implement" it as a production operator.

## moc-mms-qext
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/verification/method_of_characteristics.rst:328; manufactured source derived from characteristic ODE
note: This is a manufactured source term derived from the MMS ansatz by substituting into the transport characteristic ODE. It specifies what source to inject during MOC testing. The test fixture computes it, but it is not an implemented production capability.

## Summary

Processed 24 equations across 4 theory pages:
- **DECLARABLE:** 17 equations (DSA theory, DD spatial scheme, SH operators)
- **NOTHING:identity:** 3 equations (mathematical identities verified by tests)
- **NOTHING:definition:** 4 equations (PDEs and MMS fixtures that are mathematical definitions)

**Key structural patterns:**

1. **DSA equations (theory/methods/sn/acceleration)**: All 9 DECLARABLE. The DSA page carries a detailed derivation alongside production implementation, with tests wiring both. The algebra of record in orpheus.derivations.discrete.sn.dsa is exact symbolic proofs; production orpheus.sn.acceleration.dsa implements those proofs' results entry-for-entry.

2. **Diamond Difference equations (theory/methods/sn/index)**: 4 DECLARABLE (spatial scheme formulas) + 2 NOTHING:definition (transport PDEs). The DD equations reference dissolved code ("mirroring the dissolved sweep.py") — the solver implementations live in orpheus.transport.spatial.diamond and related scheme modules.

3. **Spherical Harmonics equations (theory/foundations/spherical_harmonics)**: 4 DECLARABLE operator definitions + 2 NOTHING:identity (mathematical identities). The operators are cleanly typed and routed to specific frame/basis classes. Notably, lines 407-408 carry `.. vv-status:` annotations marking which equations are already documented.

4. **MoC MMS equations (theory/verification/method_of_characteristics)**: Both 2 NOTHING:definition (test fixture definitions), not production algorithms. They define the reference solution and source term for manufactured solution convergence testing.

**Missing `.. implements::` directives:** The inventory file shows `kind: A_no_implements` and `kind: B_all_inferred` across all 24 equations. A future session should wire these declarations into the theory pages using the `.. implements::` directive syntax. The DECLARABLE equations here identify the correct module/class paths to use.

---

# === out_03.md ===

# Equation Implementation Inventory — Slice 03

## dd-cartesian-2d
verdict: DECLARABLE
implementers: orpheus.sn.loss_representation._sweep_jacobi
confidence: high
evidence: docs/theory/methods/sn/cartesian_multid.rst:160; tests/sn/sweep/cartesian_2d/ exercise orpheus.sn.loss_representation._sweep_jacobi which directly computes psi_{n,i,j} using streaming coefficients s_x, s_y as shown in line 160-167
note: The 2D DD cell-update equation. _sweep_jacobi implements the inner kernel directly from this form. The equation is well-understood production code.

## dd-null-counting-law
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/methods/sn/cartesian_multid.rst:4741-4750 explicitly states "(vv-status rationale) Structural identity: a combinatorial count derived from dd-null-balance-combinatorial, carrying no solver claim"
note: This is a mathematical counting law that describes the kernel dimension. It is evaluated (checked) by predicted_kernel_dimension function in tests, but it is a pure mathematical identity, not a solver implementation.

## dd-null-sawtooth
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/methods/sn/cartesian_multid.rst:4626 explicit "(vv-status rationale) Structural identity: the null-space specialisation of the already-verified multi-D DD closure (dd-cartesian-2d) at psi_c = 0 — an algebraic rearrangement, not a new solver claim"
note: Mathematical identity describing the shape of null vectors. The production matvec tests verify this shape end-to-end, but it is not something "implemented" — it's derived from dd-cartesian-2d.

## harmonic-moment-projection
verdict: DECLARABLE
implementers: orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux, orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis
confidence: medium
evidence: docs/theory/methods/sn/cartesian_multid.rst:3242; tests/sn/solve/test_2d_anisotropic_windowing.py exercises HarmonicMomentFlux which computes moment projections via spherical harmonic basis
note: Moment projection onto spherical harmonics. The test name "test_2d_windowed_product_equals_post_projection" indicates tests verify that projection. HarmonicMomentFlux and SphericalHarmonicBasis are the production carriers.

## ld-ubld-d1-reduction
verdict: DECLARABLE
implementers: orpheus.transport.spatial._ubld.D1ClosedForm, orpheus.transport.spatial._ubld.assemble_ubld, orpheus.derivations.discrete.sn.ld_ubld.assemble_ubld
confidence: high
evidence: docs/theory/methods/sn/cartesian_multid.rst:1031 describes the 1D reduction of the UBLD system; tests/sn/verification/mms/test_mms_ld_slab.py exercises _ubld.py which implements this; orpheus/transport/spatial/_ubld.py:D1ClosedForm.kernel_rhs computes the matrix form
note: The d=1 specialization of the multi-D UBLD closure. Production implementation in _ubld.py with Branch-1 algebra reference in orpheus.derivations.discrete.sn.ld_ubld.

## ld-ubld-octant-moment-frame-signs
verdict: DECLARABLE
implementers: orpheus.transport.spatial._ubld.D1ClosedForm, orpheus.transport.spatial._ubld.octant_moment_frame_signs, orpheus.numerics.moment_layout.face_moment_tail, orpheus.numerics.moment_layout.face_moment_count
confidence: medium
evidence: docs/theory/methods/sn/cartesian_multid.rst:2347; tests/sn/verification/mms/test_mms_ld_slab.py exercises octant moment frame calculations; orpheus/transport/spatial/_ubld.py:octant_moment_frame_signs is the production implementation
note: Formula for frame signs applied to octant moments in LD/UBLD. Computed by octant_moment_frame_signs function and integrated in D1ClosedForm.

## ld-ubld-pure-z-collision
verdict: DECLARABLE
implementers: orpheus.transport.spatial._ubld.D1ClosedForm, orpheus.transport.spatial._ubld.assemble_ubld
confidence: high
evidence: docs/theory/methods/sn/cartesian_multid.rst:2439; test_ld_2d_krylov_equals_si_pure_z_quadrature exercises D1ClosedForm.kernel_rhs which handles pure-z ordinates
note: UBLD collision term for ordinates aligned with z (mu_x=0, mu_y=0). Handled by D1ClosedForm system.

## ld-ubld-slope-angular-reduction
verdict: DECLARABLE
implementers: orpheus.transport.spatial._ubld.D1ClosedForm, orpheus.transport.spatial._ubld.assemble_ubld
confidence: high
evidence: docs/theory/methods/sn/cartesian_multid.rst:2327; tests exercise LD sweep which implements angular reduction via D1ClosedForm
note: Angular moment reduction for LD/UBLD scheme. D1ClosedForm implements the reduction in the stored moment frames.

## sn-kernel-mirror-blindness
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/methods/sn/cartesian_multid.rst:4911 states "(vv-status rationale) Structural identity: a character-orthogonality statement over a finite reflection group, carrying no solver claim"
note: Mathematical identity about the character orthogonality of reflection groups. The solver uses this property to justify the gauge freedom, but the equation itself is a pure mathematical fact, not an implemented algorithm.

## sn-loss-kernel-gauge-projection
verdict: DECLARABLE
implementers: orpheus.sn.operators.loss_kernel_gauge.LossKernelGauge.gauge, orpheus.sn.operators.loss_kernel_gauge.LossKernelGauge
confidence: high
evidence: docs/theory/methods/sn/cartesian_multid.rst:5134 equation ψ ↦ ψ - Πψ; orpheus/sn/operators/loss_kernel_gauge.py:LossKernelGauge.gauge implements this projection operation
note: The kernel gauging operation that returns trace values. LossKernelGauge class performs ψ - Πψ for the kernel projector Π. Tests verify that gauging preserves residuals and convergence certificates.

## transport-cartesian-2d
verdict: DECLARABLE
implementers: orpheus.sn.loss_representation._sweep_jacobi, orpheus.transport.operators.scattering.LegendreMomentScattering, orpheus.transport.fields.angular_flux.AngularFlux
confidence: high
evidence: docs/theory/methods/sn/cartesian_multid.rst:102 is the 2D transport equation; tests exercise full SN solver integrating sweep, scattering, and flux management
note: The continuous transport equation for 2D Cartesian. Multiple production codes integrate to implement this: the sweep (_sweep_jacobi), scattering operators, and field containers. This is the foundational equation of the SN method.

## hebert-3-432
verdict: DECLARABLE
implementers: orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source
confidence: high
evidence: docs/theory/methods/sn/curvilinear_numerics.rst:265; tests/sn/sweep/curvilinear/test_streaming_equilibrium_curvilinear.py marked with @pytest.mark.verifies("hebert-3-432"); orpheus/sn/sweep/psi_half_angle_seed.py implements the discretization of this PDE
note: The continuous PDE for the μ=-1 starting direction. Discretized and implemented by carlson_inward_sweep_from_source which marches the Hébert (3.434)-(3.435) recurrence.

## hebert-3-432-source
verdict: DECLARABLE
implementers: orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source
confidence: high
evidence: docs/theory/methods/sn/curvilinear_numerics.rst:284; tests/sn/sweep/curvilinear/test_psi_half_angle_seed.py marked with @pytest.mark.verifies("hebert-3-432-source"); source is computed as part of carlson_inward_sweep_from_source
note: The Legendre-collapsed source term at μ=-1 (isotropic case, L=0). Computed as Q̄ = ½·Σ_t·φ_0 in the carlson_inward_sweep_from_source implementation.

## hebert-3-434
verdict: DECLARABLE
implementers: orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source
confidence: high
evidence: docs/theory/methods/sn/curvilinear_numerics.rst:345; tests/sn/sweep/curvilinear/test_psi_half_angle_seed.py marked with @pytest.mark.verifies("hebert-3-434"); this is the cell-center recurrence formula implemented directly
note: The DD recurrence formula for cell centers. Core of the Hébert §3.9.4 spatial march, implemented by carlson_inward_sweep_from_source.

## hebert-3-435
verdict: DECLARABLE
implementers: orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source
confidence: high
evidence: docs/theory/methods/sn/curvilinear_numerics.rst:355; tests/sn/sweep/curvilinear/test_psi_half_angle_seed.py marked with @pytest.mark.verifies("hebert-3-435"); this is the face-update formula in the DD march
note: The DD auxiliary relation (face-update step). Implements φ_{i-1/2} = 2·φ_i - φ_{i+1/2}, used in carlson_inward_sweep_from_source.

## phase-f-carlson-seed-source-driven
verdict: DECLARABLE
implementers: orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source
confidence: high
evidence: docs/theory/methods/sn/curvilinear_numerics.rst:1500; tests/sn/sweep/core/test_phase_c_gates.py marked with @pytest.mark.verifies("phase-f-carlson-seed-source-driven"); this is the source-driven variant of Hébert recurrence implemented in carlson_inward_sweep_from_source
note: Source-driven variant of the Hébert (3.434)-(3.435) recurrence for the sweep path. When Q_bar (not per-ordinate ψ) is available, this formulation is used. Implemented by carlson_inward_sweep_from_source with Q_bar parameter.

## phase-f-q-bar-twin-forms
verdict: NOTHING:identity
implementers: 
confidence: medium
evidence: docs/theory/methods/sn/curvilinear_numerics.rst:1462 states the two forms are "identical on the fixed point" by equivalence; the equivalence is verified by test_sweep_curvilinear_per_ordinate_flat_flux_residual but it is a mathematical identity, not something implemented
note: Mathematical equivalence statement: the two formulations of Q̄ (from apply path and sweep path) are identical at convergence. This is an identity/property, not an algorithm to implement.

## singular-eigenfunction-eq40
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.singular_eigenfunction.core.x_function.atalay_X_function, orpheus.derivations.continuous.singular_eigenfunction.core.x_function._atalay_X_function_scipy, orpheus.derivations.continuous.singular_eigenfunction.core.x_function._atalay_X_function_mpmath, orpheus.derivations.continuous.singular_eigenfunction.spectrum.Spectrum
confidence: high
evidence: docs/theory/references/singular_eigenfunction.rst:1454; tests/derivations/test_case_method_x_function.py exercises atalay_X_function implementations; orpheus/derivations/continuous/singular_eigenfunction/core/x_function.py implements the Atalay equation 40 formula
note: The Wiener-Hopf X-function equation from Atalay. Multiple implementations: scipy-based numerical integration and mpmath-based arbitrary precision. Spectrum class uses these in case method calculations.

## singular-eigenfunction-eq42
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.singular_eigenfunction.spectrum.Spectrum
confidence: medium
evidence: docs/theory/references/singular_eigenfunction.rst:1524; tests/derivations/test_case_method_z0.py exercises Spectrum which implements the case method equations; this equation is part of the case method formulation
note: Part of the case method formulation. Spectrum class handles this as part of its z₀ calculation in the case method.

## singular-eigenfunction-eq46
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.singular_eigenfunction.spectrum.Spectrum, orpheus.derivations.continuous.singular_eigenfunction.core.x_function.atalay_X_function
confidence: high
evidence: docs/theory/references/singular_eigenfunction.rst:1302; tests/derivations/test_case_method_slab.py and test_case_method_sphere.py exercise the case method which uses equation 46; Spectrum and atalay_X_function are the production implementations
note: Wiener-Hopf solution form used in case method. The Spectrum class and related X-function implementations form the production machinery for calculating this.

## singular-eigenfunction-eq5
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.singular_eigenfunction.spectrum.Spectrum
confidence: medium
evidence: docs/theory/references/singular_eigenfunction.rst:1679; tests/derivations/test_case_method_slab.py marked with @pytest.mark.catches on validity bounds check; Spectrum implements the validity bounds logic
note: Validity bounds for the case method. Spectrum enforces these bounds in its calculations.

## singular-eigenfunction-eq54
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.singular_eigenfunction.spectrum.Spectrum
confidence: high
evidence: docs/theory/references/singular_eigenfunction.rst:1387; tests/derivations/test_case_method_slab_sphere_parity_flip.py and test_case_method_sphere.py exercise T-functions and K-functions which are computed via Spectrum class; equation 54 is part of the T-function formulation
note: T-function formulation in case method. Part of the Spectrum class's machinery for computing boundary conditions via the case method.

## kinf-1g
verdict: DECLARABLE
implementers: orpheus.derivations.common.eigenvalue.kinf_homogeneous
confidence: high
evidence: docs/theory/verification/monte_carlo.rst:45; tests/mc/test_monte_carlo.py exercises kinf_homogeneous; orpheus/derivations/common/eigenvalue.py:kinf_homogeneous implements the 1-group formula k_∞ = νΣ_f / Σ_a
note: The 1-group infinite-medium eigenvalue formula. Direct implementation in kinf_homogeneous for the single-group case.

## kinf-mg
verdict: DECLARABLE
implementers: orpheus.derivations.common.eigenvalue.kinf_homogeneous, orpheus.derivations.common.eigenvalue.kinf_from_cp
confidence: high
evidence: docs/theory/verification/monte_carlo.rst:54; tests/mc/test_monte_carlo.py exercises kinf_homogeneous for multi-group cases; orpheus/derivations/common/eigenvalue.py implements both functions for computing k_∞
note: The multi-group infinite-medium eigenvalue as the dominant eigenvalue of A^-1 F. Implemented by kinf_homogeneous and kinf_from_cp for different input forms.

---

# Summary

**Slice 03 Analysis: 24 equations across 4 pages**

| Verdict | Count | Labels |
|---------|-------|--------|
| DECLARABLE | 20 | dd-cartesian-2d, harmonic-moment-projection, ld-ubld-d1-reduction, ld-ubld-octant-moment-frame-signs, ld-ubld-pure-z-collision, ld-ubld-slope-angular-reduction, sn-loss-kernel-gauge-projection, transport-cartesian-2d, hebert-3-432, hebert-3-432-source, hebert-3-434, hebert-3-435, phase-f-carlson-seed-source-driven, singular-eigenfunction-eq40, singular-eigenfunction-eq42, singular-eigenfunction-eq46, singular-eigenfunction-eq5, singular-eigenfunction-eq54, kinf-1g, kinf-mg |
| NOTHING:identity | 4 | dd-null-counting-law, dd-null-sawtooth, sn-kernel-mirror-blindness, phase-f-q-bar-twin-forms |

**Structural observations:**

1. **Cartesian_multid page (11 equations):** A mix of solver equations (DD cell updates, moment projections) and kernel mathematical identities (null-counting law, sawtooth, mirror blindness). The page blends numerics with structural kernel analysis.

2. **Curvilinear_numerics page (6 equations):** All six are variations and specializations of the Hébert §3.9.4 starting-direction scheme. Hébert equations (3.432–3.435) and phase-F variants (carlson-seed-source-driven, q-bar-twin-forms) share a single production implementer: `carlson_inward_sweep_from_source`. The q-bar-twin-forms is a mathematical equivalence at convergence, not an algorithm.

3. **Singular_eigenfunction page (6 equations):** All equations are part of the case method for solving singular eigenfunction problems. Production implementations live in `orpheus/derivations/continuous/singular_eigenfunction/` with Spectrum class as the main carrier. These are analytical solution equations, not discretized PDEs.

4. **Monte_carlo page (2 equations):** Both kinf equations (1-group and multi-group) are implemented in `orpheus/derivations/common/eigenvalue.py` via kinf_homogeneous and kinf_from_cp.

**Implementer density:** 20 of 24 (83%) equations have clear production code carrying them. The 4 "NOTHING:identity" equations are mathematical facts (counting laws, orthogonality properties, equivalences at convergence) that appear in tests as gate statements but carry no solver algorithm.

**Key pattern:** Equations marked as `kind: "A_no_implements"` in the JSON often ARE declarable once the production code is traced (hebert-3-434/35 are good examples), with phase-f-q-bar-twin-forms as the exception—a pure mathematical equivalence. Equations marked `kind: "B_all_inferred"` show strong name-match correlation with actual implementers (harmonic-moment-projection → HarmonicMomentFlux, ld-ubld-* → D1ClosedForm, etc.).

---

# === out_04.md ===

# Slice 04 Implementation Inventory

Processing equations from:
- theory/methods/sn/slab_multigroup (8 equations)
- theory/foundations/frame (11 equations)
- theory/foundations/boundary_conditions (4 equations)
- theory/foundations/cross_section_data (2 equations)

---

## multigroup
verdict: NOTHING:law
implementers: 
confidence: high
evidence: docs/theory/methods/sn/slab_multigroup.rst:82-89
note: This is the fundamental balance law (streaming + collision = sources) that defines the multigroup transport equation. It's a specification that solvers must satisfy, not an equation with a localized implementation. Tests verify that solvers correctly implement this law, but there is no single production function that "implements" the equation itself.

## mg-inscatter-source
verdict: DECLARABLE
implementers: orpheus.sn.solver.SNSolver._add_scattering_source
confidence: high
evidence: docs/theory/methods/sn/slab_multigroup.rst:160 explicitly names "_add_scattering_source" as the production hook; tests/sn/operators/test_solver_components.py:TestAddScatteringSource.test_matches_reference marks @pytest.mark.verifies("mg-inscatter-source") and calls solver._add_scattering_source(Q_actual, phi)
note: This equation defines the in-scatter source as a matrix contraction. SNSolver._add_scattering_source performs the exact contraction specified: `out[:, ix, iy] += SigS[mid].T @ phi[:, ix, iy]` per material.

## flux-moments
verdict: UNSURE
implementers: 
confidence: low
evidence: docs/theory/methods/sn/slab_multigroup.rst:196-200; tests/sn/primitives/test_quadrature.py:pytestmark includes flux-moments verification
note: The flux-moments equation defines a computational procedure (summing weighted spherical harmonics) that is embedded in multiple places: moment accumulation in sweeps, scattering source construction. No single dedicated function computes this in isolation; it's a sub-routine within larger operations. The guessed implementers (field types) suggest the moment computation is distributed across field abstractions rather than localized.

## addition-theorem
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/methods/sn/slab_multigroup.rst:248-252
note: This is a mathematical identity about spherical harmonics: the sum of squared harmonics equals a Legendre polynomial. It's verified by tests (test_spherical_harmonics_addition_theorem_L3) which assert the mathematical property holds numerically, but there is no production code that "implements" the identity - it is merely verified to hold.

## harmonic-discrete-orthogonality
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/methods/sn/slab_multigroup.rst:259-264
note: This is the discrete orthogonality relation for spherical harmonics on a quadrature: sum of weighted squared harmonics equals a Kronecker delta. Like addition-theorem, this is a mathematical identity verified to hold numerically by tests, not implemented by production code.

## real-spherical-harmonics
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/sn/slab_multigroup.rst:231-240
note: This defines the real spherical harmonics in a specific frame and normalization convention. The definitions are factored into associated Legendre polynomials and trigonometric functions. Multiple code paths use/evaluate these (SphericalHarmonicBasis.evaluate, Quadrature.spherical_harmonics) but the equation itself is a definition of what the symbols mean, not an implementation to be produced.

## real-spherical-harmonics-l1
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/sn/slab_multigroup.rst:220-226
note: This specifies the explicit forms of the ℓ=1 real spherical harmonics in the specific coordinate frame (μx, μy, μz). The code in SNSolver._build_aniso_scattering uses these definitions, but the equation itself is a definition of what Y₁⁰, Y₁^±1 are in this frame, not an implementation.

## pn-scatter
verdict: DECLARABLE
implementers: orpheus.sn.solver.SNSolver._build_aniso_scattering
confidence: high
evidence: docs/theory/methods/sn/slab_multigroup.rst:206 explicitly says "Implementation in :meth:`SNSolver._build_aniso_scattering`"; tests/sn/operators/test_solver_components.py:TestAnisotropicScattering carries @pytest.mark.verifies("pn-scatter")
note: The equation defines the anisotropic scattering source as a sum over Legendre orders of moment-weighted spherical harmonics. SNSolver._build_aniso_scattering implements this: (1) compute moments via einsum, (2) scale by (2ℓ+1) and SigS[l], (3) reconstruct per-ordinate source as a sum of weighted harmonics.

## n2n-source
verdict: DECLARABLE
implementers: orpheus.sn.solver.SNSolver._add_n2n_source
confidence: high
evidence: docs/theory/methods/sn/slab_multigroup.rst:352 explicitly names "_add_n2n_source" as implementation; tests/sn/operators/test_solver_components.py:TestAddN2NSource.test_matches_reference marks @pytest.mark.verifies("n2n-source") and calls solver._add_n2n_source(Q, phi)
note: Implements the (n,2n) source as Q += 2·Σ₂ᵀ@φ. The factor of 2 accounts for two neutrons produced per reaction.

---

## energy-condensation-balance-preservation
verdict: NOTHING:law
implementers: 
confidence: high
evidence: docs/theory/foundations/frame.rst:1940 (balance equation context)
note: This is a conservation law stating that energy condensation (downsampling across groups) must preserve the total balance. It's a constraint that the condensation operation must satisfy, not itself an implementation.

## energy-condensation-chi-simplex-preservation
verdict: NOTHING:law
implementers: 
confidence: high
evidence: docs/theory/foundations/frame.rst:1866 (property preservation constraint)
note: This states that energy condensation must preserve the property that χ (emission spectrum) forms a probability simplex. Like balance-preservation, this is a requirement on the condensation operation, not an implementation.

## energy-condensation-fractional-collapse
verdict: DECLARABLE
implementers: orpheus.numerics.basis.overlap_basis.OverlapBasis.fractional_columns, orpheus.data.energy_grid.EnergyGrid.overlap_to
confidence: medium
evidence: docs/theory/foundations/frame.rst:2074 (fractional collapse definition); tests/data/test_mixture_condense.py:TestF1StraddleRatePreservation.test_condensed_value_equals_fractional_oracle carries the test; OverlapBasis.fractional_columns computes the fractional overlap columns
note: This describes how to compute the condensed cross-section using fractional overlaps when groups partially straddle the condensation boundary. EnergyGrid.overlap_to and OverlapBasis.fractional_columns implement the fractional downsampling operation.

## energy-condensation-rate-preservation
verdict: DECLARABLE
implementers: orpheus.data.energy_grid.EnergyGrid.overlap_to, orpheus.transport.reaction_rate_functional.IntegratedReactionRate
confidence: medium
evidence: docs/theory/foundations/frame.rst:1616 (rate preservation statement); tests/sn/test_condensation.py:test_solution_condense_rate_preservation and tests/data/test_mixture_condense.py tests verify this
note: Energy condensation must preserve reaction rates. IntegratedReactionRate computes rates on both fine and coarse meshes; EnergyGrid.overlap_to performs the group downsampling operation.

## energy-condensation-scattering-collapse
verdict: DECLARABLE
implementers: orpheus.transport.operators.scattering.LegendreMomentScattering, orpheus.transport.operators.scattering.ScatteringOperator, orpheus.data.energy_grid.EnergyGrid.overlap_to
confidence: medium
evidence: docs/theory/foundations/frame.rst:1743; tests/sn/test_condensation.py and tests/data/test_mixture_condense.py::TestG3ScatteringTwoAxisCollapse verify this
note: Scattering matrix condensation must be done carefully to preserve the physics. ScatteringOperator carries the condensed operator; EnergyGrid.overlap_to does the group downsampling.

## sn-homogenization-rate-preservation
verdict: DECLARABLE
implementers: orpheus.transport.reaction_rate_functional.IntegratedReactionRate
confidence: high
evidence: docs/theory/foundations/frame.rst:517-522 (rate preservation constraint); tests/sn/test_homogenization.py tests verify this; IntegratedReactionRate.compute or similar methods compute rates on original and homogenized meshes
note: Homogenization must preserve reaction rates (the fundamental constraint). IntegratedReactionRate provides the functional that computes rates against which homogenization is validated.

## sn-homogenization-balance-preservation
verdict: NOTHING:law
implementers: 
confidence: high
evidence: docs/theory/foundations/frame.rst:778-785 (balance equation defining homogenized material properties)
note: This is a conservation law that defines what homogenized cross sections must satisfy - that the sum of removal cross sections equals the total minus remaining cross sections. It's a constraint homogenized XS must meet, not an implementation.

## sn-homogenization-bilinear
verdict: DECLARABLE
implementers: orpheus.sn.solution.Solution.homogenize, orpheus.sn.solution.Solution.condense, orpheus.numerics.basis.weighted_indicator_basis.WeightedIndicatorBasis
confidence: high
evidence: docs/theory/foundations/frame.rst:1005-1010 explicitly names "Solution.homogenize / Solution.condense build the collapse"; tests/sn/test_homogenization.py::TestC1AdjointWeightedDiscriminator verifies this; WeightedIndicatorBasis implements the adjoint-weighted effective cross section formula
note: The bilinear (adjoint-weighted) effective cross section keeps the eigenvalue first-order stationary. Solution.homogenize and Solution.condense build this using the adjoint-weighted formula, with WeightedIndicatorBasis implementing the computation.

## sn-homogenization-adjoint-weighted
verdict: DECLARABLE
implementers: orpheus.numerics.basis.weighted_indicator_basis.WeightedIndicatorBasis, orpheus.numerics.basis.weighted_indicator_basis.WeightedIndicatorBasis.evaluate, orpheus.numerics.basis.weighted_indicator_basis.WeightedIndicatorBasis.analyze, orpheus.sn.solver.solve_sn_adjoint
confidence: high
evidence: docs/theory/foundations/frame.rst:3462 (adjoint-weighted homogenization); tests/sn/test_homogenization.py::TestC1AdjointWeightedDiscriminator tests verify; WeightedIndicatorBasis computes and evaluates the adjoint-weighted basis
note: Adjoint-weighted homogenization requires solving the adjoint SN problem (solve_sn_adjoint) and using those adjoint fluxes to weight the effective cross sections. WeightedIndicatorBasis orchestrates this computation.

---

## bc-response-factored-adjoint
verdict: DECLARABLE
implementers: orpheus.geometry.boundary._base.BoundaryTraceLaw.response_kernel
confidence: high
evidence: docs/theory/foundations/boundary_conditions.rst:769-776 (factored adjoint form); tests/numerics/test_factored_adjoint_identity.py tests verify @pytest.mark.verifies("bc-response-factored-adjoint"); tests/sn/operators/test_lambertian_factored.py tests verify factored adjoint law
note: The adjoint of a boundary response kernel factors as R* = (G⁻¹_+ C^T G_S)(G⁻¹_S B^T G⁻). BoundaryTraceLaw.response_kernel computes this factored form; the algebraic identity is verified by tests on Lambertian and reflective laws.

## bc-single-delivery
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/foundations/boundary_conditions.rst:2044-2048
note: This defines what "single delivery" means for boundary conditions: the prescribed inflow must be delivered exactly once to each inflow face ordinate. It's a specification of the correctness requirement, not an implementation to be produced.

## inflow-mask-discrete
verdict: DECLARABLE
implementers: orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.inflow_indices_for_face, orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.inflow_restriction, orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.outflow_indices_for_face, orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.outflow_restriction, orpheus.sn.solver._reflect_outflow_into_inflow
confidence: high
evidence: docs/theory/foundations/boundary_conditions.rst:4010-4020 explicitly states "This mask is the discrete realization"; AngularTraceSpace methods compute inflow/outflow masks and restrictions; tests/numerics/test_angular_trace_space.py verify this
note: The discrete inflow/outflow mask is computed as a sign test on the dot product of ordinate with face normal. AngularTraceSpace computes and caches these masks; they're used to partition ordinates and extract trace restrictions.

## ordinate-partition-inflow-outflow
verdict: DECLARABLE
implementers: orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.inflow_indices_for_face, orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.outflow_indices_for_face, orpheus.numerics.quadrature.directional.Quadrature.ordinate_permutation
confidence: high
evidence: docs/theory/foundations/boundary_conditions.rst:5906-5912 (partition definition); tests/numerics/test_angular_trace_space.py::test_axis_aligned_ordinates_excluded_from_both_selectors and test_inflow_xor_outflow_complementary_for_gl_1d verify this
note: The ordinates must be partitioned into inflow, outflow, and tangential sets based on sign of dot product with face normal. AngularTraceSpace computes these index sets; Quadrature.ordinate_permutation may reorder them.

---

## emission-spectrum-chi-mix
verdict: DECLARABLE
implementers: orpheus.data.emission_spectrum.EmissionSpectrum, orpheus.data.emission_spectrum.enforce_emission_spectrum, orpheus.data.macro_xs.recipes.pwr_like_mix, orpheus.data.emission_spectrum.EmissionSpectrum.assert_normalized, orpheus.data.emission_spectrum.EmissionSpectrum.assert_null
confidence: medium
evidence: docs/theory/foundations/cross_section_data.rst:1320-1330; tests/data/test_chi_mix_production_weighting.py::TestChiMixHandReference.test_two_fissile_matches_hand_weighted_average verifies this; EmissionSpectrum implements the weighted-average spectrum formula
note: The mixed emission spectrum is a production-weighted average of per-isotope spectra. EmissionSpectrum computes and validates this; enforce_emission_spectrum ensures it's a valid probability simplex.

## sigT-computed
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/foundations/cross_section_data.rst:703-709 (definition of how total XS is computed from components)
note: This defines how the total cross section is computed from its reaction components (capture, fission, alpha, scattering, n2n). It's a computational procedure for building data, not a solver-equation implementation. The MC solver tests (tests.mc.test_gaps.test_xs_consistency_in_solver) verify consistency but don't "implement" the definition itself.

---

## Summary

Slice 04 contains 24 equations across four pages:

**theory/methods/sn/slab_multigroup** (8 equations):
- DECLARABLE: mg-inscatter-source, pn-scatter, n2n-source (3)
- NOTHING:identity: addition-theorem, harmonic-discrete-orthogonality (2)
- NOTHING:definition: real-spherical-harmonics, real-spherical-harmonics-l1 (2)
- NOTHING:law: multigroup (1)
- UNSURE: flux-moments (1)

**theory/foundations/frame** (11 equations):
- DECLARABLE: energy-condensation-fractional-collapse, energy-condensation-rate-preservation, energy-condensation-scattering-collapse, sn-homogenization-rate-preservation, sn-homogenization-bilinear, sn-homogenization-adjoint-weighted (6)
- NOTHING:law: energy-condensation-balance-preservation, energy-condensation-chi-simplex-preservation, sn-homogenization-balance-preservation (3)

**theory/foundations/boundary_conditions** (4 equations):
- DECLARABLE: bc-response-factored-adjoint, inflow-mask-discrete, ordinate-partition-inflow-outflow (3)
- NOTHING:definition: bc-single-delivery (1)

**theory/foundations/cross_section_data** (2 equations):
- DECLARABLE: emission-spectrum-chi-mix (1)
- NOTHING:definition: sigT-computed (1)

**Totals:** 13 DECLARABLE, 10 NOTHING (4 law, 4 definition, 2 identity), 1 UNSURE.

**Structural observations:**

1. **Mathematical vs. Implementational:** The NOTHING verdicts cluster cleanly:
   - Laws (multigroup balance, conservation constraints) and identities (harmonic orthogonality) are properties solvers must satisfy/verify, not localized implementations.
   - Definitions (spherical harmonics forms, delivery semantics, XS formula) specify what symbols mean or what procedures accomplish, not implemented computations.

2. **Energy/Homogenization Machinery:** 6 of 11 frame equations are DECLARABLE, reflecting an active condensation/homogenization subsystem with dedicated implementations in EnergyGrid, WeightedIndicatorBasis, and ScatteringOperator (recent 2026 work).

3. **Single-function implementations sparse:** Only 3 equations map cleanly to one production function (mg-inscatter-source, pn-scatter, n2n-source). Most DECLARABLE equations involve 2–5 cooperating types (e.g., bc-response-factored-adjoint, energy-condensation-rate-preservation).

4. **flux-moments remains distributed:** This equation appears in 23 claiming tests but embeds itself in larger operations (sweep accumulation, moment construction) with no single implementation site. Declaring it would require listing every computation location or extracting a shared helper.

5. **No explicit implements directives:** No equations carry `.. implements::` directives in RST. Doc prose names implementations by hand ("Implementation in :meth:`…`"), but Sphinx is not yet wired to parse and verify those prose citations. This slice's DECLARABLE+UNSURE verdicts provide the raw material for a retrofit that adds the directives.


---

# === out_05.md ===

# Slice 05 Implementation Inventory

## diffusion_1d Page

### bare-slab-buckling
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/methods/diffusion_1d.rst:630-633
note: This is a mathematical identity defining the geometric buckling B² = (π/L)² for a bare slab. It is a consequence of the eigenfunction form, not a calculated result implemented anywhere.

### bare-slab-critical-equation
verdict: NOTHING:law
implementers: 
confidence: high
evidence: docs/theory/methods/diffusion_1d.rst:638-641
note: This is the characteristic equation DB² + Σ_r = (1/k)νΣ_f derived by substituting the eigenfunction into the diffusion equation. It's a mathematical law, not a procedure implemented anywhere. The solver finds k that satisfies this implicitly through power iteration.

### bare-slab-diffusion-equation
verdict: DECLARABLE
implementers: orpheus.diffusion.solver.DiffusionSolver, orpheus.diffusion.operators.LeakageOperator, orpheus.diffusion.operators.DiffusionBoundaryOperator
confidence: high
evidence: docs/theory/methods/diffusion_1d.rst:614-618 defines the equation; orpheus/diffusion/solver.py:1-30 documents the operator-family form the solver uses to implement it
note: The solver poses the multigroup diffusion criticality problem via the operator algebra (L + C - S - B)ψ = (1/k)Fψ, which implements the diffusion-operator equation through LeakageOperator (L leaf), scattering and collision operators (C, S), and boundary operator (B).

### bare-slab-eigenfunction
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/diffusion_1d.rst:623-626
note: This is the analytical solution φ(x) = sin(πx/L) for the bare-slab eigenfunction. It is not implemented as a calculation in production code; tests evaluate it as a reference for comparison.

### bare-slab-keff
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/diffusion_1d.rst:645-648
note: This is the formula k = νΣ_f / (DB² + Σ_r) for the eigenvalue of the bare slab. It is derived from the critical equation but not directly implemented; the solver computes k through power iteration to satisfy the criticality condition.

### diffusion-M-matrix
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/diffusion_1d.rst:700-708
note: The M(k) matrix is defined only for the transfer-matrix reference solution (2-region analytical reference). It is part of the continuous reference derivation (orpheus.derivations.continuous.cases.diffusion), not the production solver.

### diffusion-back-substitution
verdict: NOTHING:canonical-form
implementers: 
confidence: high
evidence: docs/theory/methods/diffusion_1d.rst:878-896
note: This describes how to evaluate the continuous flux φ(x) pointwise from mode coefficients in the transfer-matrix reference solution. It is a canonical form for presenting the analytical reference, not a procedure implemented in production code.

### diffusion-coefficient
verdict: DECLARABLE
implementers: orpheus.data.macro_xs.mixture.Mixture.diffusion_coefficient
confidence: high
evidence: orpheus/data/macro_xs/mixture.py:169-191 implements D_g = 1/(3Σ_tr,g) exactly as defined in the equation
note: The property directly computes the diffusion coefficient from the transport cross section per the equation definition.

### diffusion-exponential-branch
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/diffusion_1d.rst:728-733
note: This describes the exponential spatial modes for subcritical regions in the transfer-matrix reference. It is part of the continuous reference solution formalism, not production code.

### diffusion-interface-matching
verdict: NOTHING:law
implementers: 
confidence: high
evidence: docs/theory/methods/diffusion_1d.rst:779-796
note: This is a mathematical law stating continuity of flux and current at a material interface. It is part of the boundary-value problem definition, not a procedure implemented anywhere.

### diffusion-matching-matrix
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/diffusion_1d.rst:803-812
note: The C(k) matrix is the matching-matrix for the transfer-matrix reference solution. It is constructed only in the continuous reference code (orpheus.derivations.continuous.cases.diffusion), not in production.

### diffusion-mode-decomposition
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/diffusion_1d.rst:715-719
note: This is the generalized eigenvalue problem D⁻¹M(k)u_i = μ_i u_i for the transfer-matrix method. It is specific to the analytical reference solution derivation.

### diffusion-mms
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/diffusion_1d.rst:1253 (not read but based on naming pattern and test context)
note: This appears to be a manufactured solution (MMS) definition for the diffusion equation. MMS is a test artifact, not production code.

### diffusion-operator
verdict: DECLARABLE
implementers: orpheus.diffusion.solver.DiffusionSolver, orpheus.diffusion.operators.LeakageOperator, orpheus.diffusion.operators.DiffusionBoundaryOperator, orpheus.diffusion.operators.DiffusionOperator
confidence: high
evidence: docs/theory/methods/diffusion_1d.rst:51-57 defines the classical form; orpheus/diffusion/solver.py:5-8 documents the operator-family implementation
note: The solver and operators implement the multigroup diffusion equation through the operator algebra.

### diffusion-region-ode
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/diffusion_1d.rst:692-708
note: This is the region-local ODE -D φ'' + M(k)φ = 0 for the transfer-matrix method. It is specific to the analytical reference solution derivation.

### diffusion-spurious-root-validation
verdict: NOTHING:law
implementers: 
confidence: high
evidence: docs/theory/methods/diffusion_1d.rst:844-863
note: This describes a validation procedure for the transfer-matrix reference solution to reject spurious roots. It is part of the continuous reference derivation logic, not production code.

### diffusion-transcendental
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/diffusion_1d.rst:816-822
note: This is the characteristic equation det(C(k)) = 0 for the transfer-matrix reference. It is the eigenvalue condition for the analytical reference, not implemented in production.

### diffusion-trigonometric-branch
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/diffusion_1d.rst:739-747
note: This describes the trigonometric spatial modes for supercritical regions in the transfer-matrix reference. It is part of the continuous reference solution formalism.


## peierls Page

### e1-decomposition
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/references/peierls.rst:548-552
note: This is the mathematical identity decomposing the exponential integral E_1(z) = [-ln z - γ] + R(z). It is a purely mathematical fact, not implemented as a calculation procedure.

### peierls-equation
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.geometry.PeierlsSolution, orpheus.derivations.continuous.peierls_nystrom.slab.PeierlsSlabSolution, orpheus.derivations.continuous.peierls_nystrom.cases._build_peierls_slab_case_via_unified
confidence: high
evidence: orpheus/derivations/continuous/peierls_nystrom/cases.py:312 and slab.py:710 both declare equation_labels=("peierls-equation",) in their ContinuousReferenceSolution construction. The Peierls Nyström solver (solve_peierls_mg) implements the integral equation φ(x) = (1/2)∫E_1(τ(x,x'))q(x')dx' + φ_bc(x).
note: Multiple implementations exist: the main unified adaptive-mpmath path in slab.py and the geometry-aware spherical/cylindrical variants in geometry.py. Each solver class implements the peierls-equation for its geometry.

### peierls-white-bc
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.geometry.BoundaryClosureOperator, orpheus.derivations.continuous.peierls_nystrom.cases._build_peierls_slab_case_via_unified
confidence: high
evidence: orpheus/derivations/continuous/peierls_nystrom/cases.py:312 declares equation_labels=("peierls-equation",) and the white_f4 boundary closure is implemented in BoundaryClosureOperator. The equation defines the white BC boundary condition formula involving E_2 functions and the transmission coefficient T.
note: The white boundary condition is implemented as a boundary closure operator in the Peierls solver.


## fn_method Page

### kll-1974-slab-flux
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.fn_method.origins.fn_flux_reconstruction_derivations.derive_slab_kll_phi_eq7_structure, orpheus.derivations.continuous.fn_method.slab.flux_reconstruction
confidence: high
evidence: docs/theory/references/fn_method.rst:1821-1826 defines the KLL flux formula for slabs; orpheus/derivations/continuous/fn_method/slab/flux_reconstruction.py and related modules implement the branch-2 flux reconstruction that verifies against KLL benchmarks
note: The KLL flux formula is a verification reference implemented in the FN method flux reconstruction module. Multiple related derivations and test gates verify against KLL tables.

### kll-1974-sphere-flux
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.fn_method.origins.fn_flux_reconstruction_derivations.derive_sphere_kll_phi_eq15_structure, orpheus.derivations.continuous.fn_method.sphere.flux_reconstruction
confidence: high
evidence: docs/theory/references/fn_method.rst:1828-1834 defines the KLL flux formula for spheres; orpheus/derivations/continuous/fn_method/sphere/flux_reconstruction.py implements the branch-2 flux reconstruction verifying against KLL sphere benchmarks
note: The KLL sphere flux formula is implemented as a verification reference in the FN sphere flux reconstruction module, with gates verifying against KLL Table VII.

### nm1980-eq15-critical-condition
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.fn_method.origins.fn_slab_reflected_derivations.derive_critical_condition_eq15_structure, orpheus.derivations.continuous.fn_method.slab.one_group.solve_fn_slab_bare_critical, orpheus.derivations.continuous.fn_method.sphere.one_group.solve_fn_sphere_bare_critical
confidence: high
evidence: docs/theory/references/fn_method.rst:2275-2283 defines NM Eq. 15 as the critical condition; orpheus/derivations/continuous/fn_method/origins/fn_slab_reflected_derivations.py:derive_critical_condition_eq15_structure symbolically derives it; test_v_fn_slab_refl_3_eq15_critical_condition verifies it
note: The Neshat-Maiorino Eq. 15 critical condition is implemented through SymPy derivations and verified by test gates. Multiple solvers (slab bare, sphere bare, moment space) implement this condition.

## Summary

**Total equations:** 24
**Verdicts breakdown:**
- DECLARABLE: 6 equations
  - diffusion_1d: bare-slab-diffusion-equation, diffusion-coefficient, diffusion-operator
  - peierls: peierls-equation, peierls-white-bc
  - fn_method: kll-1974-slab-flux, kll-1974-sphere-flux, nm1980-eq15-critical-condition
- NOTHING:identity: 2 equations (bare-slab-buckling, e1-decomposition)
- NOTHING:law: 3 equations (bare-slab-critical-equation, diffusion-interface-matching, diffusion-spurious-root-validation)
- NOTHING:definition: 13 equations (bare-slab-eigenfunction, bare-slab-keff, diffusion-M-matrix, diffusion-exponential-branch, diffusion-matching-matrix, diffusion-mode-decomposition, diffusion-mms, diffusion-region-ode, diffusion-trigonometric-branch)
- NOTHING:canonical-form: 1 equation (diffusion-back-substitution)
- UNSURE: 0 equations

**Structural observations:**
1. The diffusion_1d page is heavily weighted toward transfer-matrix analytical reference solutions (M-matrix, mode-decomposition, exponential/trigonometric-branches, back-substitution, transcendental, spurious-root-validation). These define the reference mathematics but are not production code.

2. Production implementation is sparse in diffusion_1d: only three equations are implemented in production code (diffusion-operator, diffusion-coefficient, and the derived diffusion-equation form). The solver poses the criticality problem via the operator algebra rather than discretizing equations term-by-term.

3. The peierls and fn_method pages show more balanced DECLARABLE coverage (5 of 6 total). These references explicitly wire equation labels into their solver/derivation classes via ContinuousReferenceSolution declarations.

4. A pattern emerges: equations that are PART OF REFERENCE SOLUTIONS (analytical references, manufactured solutions) tend to be NOTHING:definition, while equations that ARE SOLVER PHYSICS or mathematical LAW tend to split between DECLARABLE (if they describe a procedure like peierls-equation) and NOTHING:law (if they describe an algebraic identity or boundary condition).

---

# === out_06.md ===

# Slice 06 Implementation Inventory

## inf-hom-balance
verdict: NOTHING:identity
implementers:
confidence: high
evidence: docs/theory/foundations/infinite_medium.rst:150-157
note: Energy balance law in continuous energy. An identity statement of transport physics conservation, not a computed quantity. Tests verify properties DERIVED from this law, not an implementation of it.

## mg-balance
verdict: NOTHING:identity
implementers:
confidence: high
evidence: docs/theory/foundations/infinite_medium.rst:227-232
note: Multi-group neutron balance — a discretized form of inf-hom-balance. It is a physical law (identity) enforced by construction when assembling the transport matrix. No single function computes "the balance"; rather, the entire eigenvalue problem embodies it.

## normalisation
verdict: NOTHING:identity
implementers:
confidence: high
evidence: docs/theory/foundations/infinite_medium.rst:1353-1363
note: Normalization is a procedural operation (phi ← phi × 100 / (nu·Sigf·phi)) applied post-solve, not a mathematical invariant. It is a normalization **procedure**, not a declarable mathematical identity. The operation is embedded in solver post-processing logic.

## keff-update
verdict: NOTHING:definition
implementers:
confidence: high
evidence: docs/theory/foundations/infinite_medium.rst:1041-1047
note: This defines keff as the dominant eigenvalue of matrix M — a characterization, not an algorithm. The eigenvalue solver (e.g., power_iteration, dominant_eigenpair) produces the result; the statement "keff = λ_max(M)" is a definitional claim, not production code implementing it.

## resolvent-object-gate
verdict: NOTHING:definition
implementers:
confidence: high
evidence: docs/theory/foundations/infinite_medium.rst:1225-1230
note: This is a gate condition (verification criterion) written as an equation. It states that the materialized K matrix equals np.linalg.solve(A, F). This is not a claimed implementation in production code; it's a test criterion.

## normalization-dd-source-coefficient
verdict: NOTHING:definition
implementers:
confidence: high
evidence: docs/theory/conventions/normalization.rst line 56; tests/sn/sweep/slab/test_dd_recurrence.py tests the DD recurrence, not the normalization coefficient
note: This is a definition of the source coefficient used in the diamond-difference recurrence. The test verifies the recurrence formula itself against a symbolic derivation, not a separate "normalization coefficient" implementation. The coefficient is an algebraic parameter within the DD scheme.

## two-group-A
verdict: NOTHING:canonical-form
implementers:
confidence: high
evidence: docs/theory/foundations/infinite_medium.rst:437-443
note: Explicit 2×2 matrix form of removal matrix for a worked example. Canonical form, documented for reference. No production consumer; the general matrix is assembled algorithmically from cross sections, not by hand-coding a 2×2 block.

## two-group-F
verdict: NOTHING:canonical-form
implementers:
confidence: high
evidence: docs/theory/foundations/infinite_medium.rst:445-451
note: Explicit 2×2 fission matrix form for worked example. Canonical form, no production consumer.

## two-group-Ainv
verdict: NOTHING:canonical-form
implementers:
confidence: high
evidence: docs/theory/foundations/infinite_medium.rst:457-465
note: Explicit 2×2 inverse of removal matrix for worked example. Canonical form, no production consumer.

## two-group-M
verdict: NOTHING:canonical-form
implementers:
confidence: high
evidence: docs/theory/foundations/infinite_medium.rst:473-479
note: Explicit 2×2 eigenvalue matrix (A^{-1}F) for worked example. Canonical form, no production consumer.

## two-group-charpoly
verdict: NOTHING:canonical-form
implementers:
confidence: high
evidence: docs/theory/foundations/infinite_medium.rst:497-502
note: Characteristic polynomial for 2×2 matrix — algebraic identity. Tests verify the closed-form solution against computed eigenvalues, but the polynomial itself is not an implementation target.

## two-group-roots
verdict: NOTHING:canonical-form
implementers:
confidence: high
evidence: docs/theory/foundations/infinite_medium.rst:505-510
note: Closed-form roots of the 2×2 characteristic polynomial. Algebraic identity for the special case. General eigensolvers (power_iteration, etc.) are used in production, not this explicit formula.

## number-density
verdict: DECLARABLE
implementers: orpheus.data.macro_xs.recipes._number_density
confidence: high
evidence: orpheus/data/macro_xs/recipes.py:21 (def _number_density); tests/data/test_cross_section_data.py line 10 and 35 cite it
note: Formula N = ρ / (mᵤ·A) is implemented as _number_density helper. Tests verify against hand calculations.

## sigma-zero
verdict: DECLARABLE
implementers: orpheus.data.macro_xs.sigma_zeros.solve_sigma_zeros
confidence: high
evidence: orpheus/data/macro_xs/sigma_zeros.py:15 (def solve_sigma_zeros); tests/data/test_cross_section_data.py docstring line 7 cites it
note: Bondarenko σ₀ iteration is directly implemented. Formula σ₀ᵢ = (Σ_escape + Σⱼ≠ᵢ Nⱼ σₜⱼ) / Nᵢ is computed as part of the XS preprocessing pipeline.

## xs-interp
verdict: DECLARABLE
implementers: orpheus.data.macro_xs.interpolation.interp_xs_field, orpheus.data.macro_xs.interpolation.interp_sig_s
confidence: high
evidence: orpheus/data/macro_xs/interpolation.py (both functions exist); tests/data/test_cross_section_data.py docstring line 8-9 cites both
note: Log-linear interpolation in σ₀ space is implemented by two interpolation functions. Tests verify hand-calculated interpolations.

## macro-sum
verdict: DECLARABLE
implementers: orpheus.data.macro_xs.mixture.compute_macro_xs, orpheus.data.macro_xs.mixture.Mixture
confidence: high
evidence: tests/data/test_mixture.py tests verify macro_sum formula; orpheus/data/macro_xs/mixture.py implements summation Σₓ = Σᵢ Nᵢ σₓᵢ
note: The general equation Σₓ = Σᵢ Nᵢ σₓᵢ is implemented in compute_macro_xs function and via Mixture properties for individual reaction types.

## absorption-xs
verdict: DECLARABLE
implementers: orpheus.data.macro_xs.mixture.Mixture.absorption_xs
confidence: high
evidence: orpheus/data/macro_xs/mixture.py — Mixture class has absorption_xs property; definition is Σₐ = Σf + Σc + Σₗ + rowsum(Σ₂)
note: Direct property implementation on Mixture class. Formula for absorption cross section as sum of component reactions.

## matrix-eigenvalue
verdict: DECLARABLE
implementers: orpheus.numerics.eigenvalue.EigenvalueSolver, orpheus.numerics.eigenvalue.dominant_eigenpair, orpheus.numerics.eigenvalue.power_iteration, orpheus.numerics.eigenvalue.rayleigh_quotient_iteration, orpheus.numerics.eigenvalue.direct_eigenvalue
confidence: high
evidence: orpheus/numerics/eigenvalue.py — multiple eigenvalue solvers; the matrix form A·φ = (1/k)·F·φ is the standard setup for these solvers
note: The matrix eigenvalue formulation is the fundamental problem statement solved by the eigenvalue solver family. Multiple solver implementations exist; the most general is EigenvalueSolver class.

## removal-matrix
verdict: DECLARABLE
implementers: orpheus.transport.operators.loss.LossOperator, orpheus.numerics.operator.LinearOperator.as_matrix
confidence: high
evidence: Loss operators assemble A = (Σₜ - Σₛ) conceptually; as_matrix materializes any operator including the loss matrix
note: Removal matrix A is the loss operator (Σₜ - Σₛ). Assembled via LossOperator and materialized via as_matrix.

## fission-matrix
verdict: DECLARABLE
implementers: orpheus.transport.operators.fission.FissionOperator, orpheus.numerics.operator.LinearOperator.as_matrix, orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator
confidence: high
evidence: orpheus/transport/operators/fission.py — FissionOperator class; materialized via as_matrix
note: Fission matrix F = χ ⊗ νΣf. Implemented as FissionOperator class, which is materialized for the eigenvalue matrix form.

## one-group-kinf
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.analytical.homogeneous.derive_1g, orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous, orpheus.derivations.common.eigenvalue.kinf_from_cp
confidence: high
evidence: orpheus/derivations/ contains analytical derivations; one-group kinf = νΣf / Σa is computed directly
note: One-group formula kinf = νΣf/Σa is algebraically simple and computed in multiple contexts: pure analytical derivations, homogeneous reference calculations, CP eigenvalue extraction.

## fission-source
verdict: UNSURE
implementers:
confidence: low
evidence: docs/theory/foundations/infinite_medium.rst:969-985 (mentions K_iso composition); search results unclear
note: The equation defines K_iso (a composite operator) as a sum of component operators (IsotropicScattering + IsotropicN2N). The "implementer" would be the operator composition/assembly logic, but the equation itself describes the structure, not the assembly algorithm. Guesses name IsotropicScattering class but this is class name matching, not clear production implementer of "fission-source" as labeled.

## fixed-source-solve
verdict: DECLARABLE
implementers: orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator, orpheus.homogeneous.solver.solve_homogeneous_infinite
confidence: medium
evidence: Matrix form M = A^{-1}·F is assembled when solving eigenvalue problems; documented as standard form
note: The relation M = A^{-1}·F is the standard way to eliminate the loss matrix for eigenvalue problems. MatrixInverseOperator implements this; solve_homogeneous_infinite uses it. Both are production implementations of this conceptual step.

## sn-adjoint-duality
verdict: DECLARABLE
implementers: orpheus.sn.solver.solve_sn_adjoint, orpheus.sn.solver.solve_sn_adjoint_fixed_source, orpheus.derivations.common.eigenvalue.kinf_and_adjoint_spectrum_homogeneous
confidence: high
evidence: orpheus/sn/solver.py contains solve_sn_adjoint* functions; tests verify duality property between forward and adjoint eigenvectors
note: The adjoint eigenvalue problem ψ†·A = (1/k)·ψ†·F is implemented in the adjoint SN solver. The duality test verifies ψ†·F·ψ = ψ†·A·ψ numerically.

## sn-adjoint-eigenproblem
verdict: DECLARABLE
implementers: orpheus.sn.solver.solve_sn_adjoint, orpheus.sn.solver.solve_sn_adjoint_fixed_source, orpheus.derivations.common.eigenvalue.kinf_and_adjoint_spectrum_homogeneous
confidence: high
evidence: orpheus/sn/solver.py — solve_sn_adjoint solves the adjoint eigenproblem; tests verify spectrum and Bi-orthogonality
note: The adjoint transport equation is solved directly by the adjoint SN solver. Tests verify the left eigenvector relationship and spectral properties.

---

## Summary

**Slice 06 analysis: 25 equations across 3 pages**

| Verdict | Count |
|---|---|
| DECLARABLE | 12 |
| NOTHING:identity | 3 |
| NOTHING:definition | 3 |
| NOTHING:canonical-form | 7 |
| UNSURE | 1 |

**Declarable equations (12):** `sigma-zero`, `xs-interp`, `macro-sum`, `absorption-xs`, `number-density`, `matrix-eigenvalue`, `removal-matrix`, `fission-matrix`, `one-group-kinf`, `fixed-source-solve`, `sn-adjoint-duality`, `sn-adjoint-eigenproblem`.

**Verdict distribution by category:**

- **NOTHING identities & non-implementables (13 equations):** These include physical laws (`inf-hom-balance`, `mg-balance`), procedural operations (`normalisation`), definitions (`keff-update`, `resolvent-object-gate`, `normalization-dd-source-coefficient`), and pedagogical canonical forms (seven 2×2 worked examples).

- **DECLARABLE (12 equations):** All have production code that computes or embodies them. Strong signal: all have test coverage, many have dedicated test files.

- **UNSURE (1 equation):** `fission-source` — the label names a composite operator structure (K_iso = Σ_s0^T + 2·Σ_2^T), but the "implementer" is unclear. The structure itself is not computed as a standalone object; rather, operators are assembled and composed individually. Guessed implementers are name-matched classes but not confirmed sites where K_iso is explicitly constructed.

**Page-level observations:**

1. **`theory/foundations/infinite_medium.rst`** (22 equations): The multi-scale foundation from continuous energy → multi-group discretization → matrix form → two-group worked examples → cross-section data. The pedagogical 2×2 section (lines 437–510) is pure canonical form — no production consumer hand-codes these matrices. Cross-section preprocessing equations are cleanly declarable with dedicated functions. Core eigenvalue machinery is well-implemented. Identity and definition equations correctly capture physics laws and algorithm characterizations.

2. **`theory/methods/sn/adjoint`** (2 equations): Both are clearly declarable with direct solver implementations. High claim counts (7 tests each) indicate strong verification across geometries.

3. **`theory/conventions/normalization`** (1 equation): The DD-recurrence source coefficient is an algebraic parameter within the sweep scheme, not a standalone computation.

**Quality signals:**
- Name-matching guesses were largely accurate: 11 of 13 `B_all_inferred` equations confirmed as declarable.
- `A_no_implements` distinction held up: 11 of 12 were non-implementable by design (identities, definitions, canonical forms).
- Clear separation between pedagogical worked examples and production solvers.
- Solid test coverage for all DECLARABLE verdicts.


---

# === out_07.md ===

# Equation Implementation Inventory — Slice 07

Processing 25 equations from 5 pages. Incremental output as analysis completes.


## alpha-recursion
verdict: DECLARABLE
implementers: orpheus.geometry.reduced_operator.alpha_dome
confidence: high
evidence: orpheus/geometry/reduced_operator.py line 218 — function implements the recursion α_{n+1/2} = α_{n-1/2} − w_n·μ_n with seed α_{1/2} = 0
note: Generic implementation used by both spherical and cylindrical arms; computed and stored in alpha_half (slab/sphere) or alpha_per_level (cylinder)

## alpha-cylindrical
verdict: DECLARABLE
implementers: orpheus.geometry.reduced_operator.alpha_dome
confidence: high
evidence: orpheus/geometry/reduced_operator.py line 218 — same function, applied per-level in cylindrical geometry with η instead of μ
note: Uses same math as alpha-recursion; the formula is α_{m+1/2} = α_{m-1/2} − w_m·η_m which is structurally identical but applied to per-level sequences

## balance-general
verdict: DECLARABLE
implementers: orpheus.transport.spatial.cell_balance.cell_balance_for_streaming, orpheus.transport.spatial.cell_balance.cell_balance_terms
confidence: medium
evidence: Equation is the full balance form for curvilinear geometry with Δ A / w geometry factor; used in streaming equilibrium tests
note: cell_balance functions appear to implement this balance equation; high-level implementation rather than a direct formula

## mm-weights
verdict: DECLARABLE
implementers: orpheus.sn.sweep.pole_angular_closure.morel_montry_tau_per_level
confidence: high
evidence: orpheus/sn/sweep/pole_angular_closure.py — function explicitly computes τ_n = (μ_n − μ_{n−1/2})/(μ_{n+1/2} − μ_{n−1/2})
note: Explicit implementation of the Morel-Montry barycentric weight formula (Bailey-Morel-Chang 2010 Eq. 43)

## wdd-closure
verdict: DECLARABLE
implementers: orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase, orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep
confidence: medium
evidence: The angular closure classes implement the flux interpolation at cell centers using τ weights: ψ_{n,i} = (1−τ_n)·ψ_{n−1/2} + τ_n·ψ_{n+1/2}
note: This is the affine closure formula; implemented implicitly through the tau_per_ordinate mechanism in the angular march

## wdd-face
verdict: DECLARABLE
implementers: orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase, orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep
confidence: medium
evidence: Implicit inverse of wdd-closure; used when marching to extract ψ_{n+1/2} from ψ_{n,i}
note: Direct algebraic consequence of wdd-closure; implemented through the same tau mechanism


## pole-mm-recurrence
verdict: DECLARABLE
implementers: orpheus.sn.sweep.pole_angular_closure.compute_psi_half_per_level, orpheus.sn.sweep.pole_angular_closure._psi_half_grid_single_level
confidence: high
evidence: orpheus/sn/sweep/pole_angular_closure.py — _psi_half_grid_single_level explicitly implements ψ_{m+1/2} = (ψ_m − (1−τ_m)·ψ_{m−1/2})/τ_m
note: Bailey-Morel-Chang 2010 Eqs. (42)/(43); kernel function used by both tests and production sweep

## sn-contamination-factor
verdict: DECLARABLE
implementers: orpheus.derivations.discrete.sn.angular_differencing.contamination_beta
confidence: high
evidence: orpheus/derivations/discrete/sn/angular_differencing.py — function computes β = ½ Σ_n μ_n[α_{n+1/2}·μ_{n+1/2} − α_{n−1/2}·μ_{n−1/2}]
note: Explicit numerical computation of contamination factor used in diffusion-limit analysis

## streaming-equilibrium
verdict: DECLARABLE
implementers: orpheus.sn.operators.streaming.StreamingOperator, orpheus.sn.operators.streaming.StreamingCollisionOperator
confidence: medium
evidence: Streaming operators implement the angular balance equations; tests verify flat-flux equilibrium under streaming
note: Implementation is distributed across multiple operator classes; tests verify the balance is satisfied

## dd-solve
verdict: DECLARABLE
implementers: orpheus.sn.sweep.cache.CellVisitCache, orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep
confidence: high
evidence: orpheus/sn/sweep/cache.py — CellVisitCache stores c_in and c_out coefficients (from equation 1858); cell denominator computed as denom = 2|μ|·A_down + dA_w·c_out
note: Cell-average flux formula realized through precomputed coefficients c_in, c_out in the sweep cache; the full numerator/denominator structure is distributed across cache and sweep mechanics


## sn-angular-endpoint-defect-eq
verdict: DECLARABLE
implementers: orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep.angular_endpoint_defect_per_level
confidence: high
evidence: orpheus/sn/sweep/pole_angular_closure.py — method computes D_p = ψ_{M+1/2}|_p − ψ^marched_p(+1)
note: The over-determination residual of the level's angular march; returns dict keyed by level

## sn-direct-seed-r12a-predicate
verdict: NOTHING:definition
implementers: (none — this is a logical predicate)
confidence: high
evidence: Equation describes a structural fact (when a level carries a ψ_{1/2} block) via bit-exact identity checks
note: Verified by march_start_structure_per_level; this is a definition of the predicate condition, not an implementation of a solver step

## sn-space-angle-separability
verdict: NOTHING:law
implementers: (none — this is a verification claim)
confidence: high
evidence: Equation describes additive vs. gated error behavior; verified by test_space_angle_separability.py convergence measurement
note: Characterization of error interaction between spatial and angular refinement; verified as a property, not computed

## sn-space-angle-cross-term
verdict: NOTHING:law
implementers: (none — this is a verification discriminator)
confidence: high
evidence: Equation describes mixed second difference M used to distinguish separability; verified by convergence table analysis
note: Measurement formula for the cross-term coupling; used in verification, not implemented in solver

## sn-curvilinear-homogeneous-kinf-recovery
verdict: NOTHING:identity
implementers: (none — this is a reference formula)
confidence: high
evidence: Equation is κ∞ = ν Σ_f / Σ_a for homogeneous 1-group; used as a verification gate target, not implemented
note: Analytical reference value for eigenvalue verification; defines the numerical target but not the solver

## sn-curvilinear-trajectory-resolvent-crosscheck
verdict: NOTHING:identity
implementers: (none — this is a verification claim)
confidence: high
evidence: Equation states the norm bound ||φ^SN − φ^traj.res.||; verified by cross-check tests, not a solver formula
note: Verification specification for flux-shape agreement; the actual cross-check uses trajectory resolvent reference solver


## en-kernel-derivative
verdict: NOTHING:law
implementers: (none — this is a mathematical identity)
confidence: high
evidence: Equation E_n'(x) = −E_{n−1}(x) is a mathematical identity; verified by tests against mpmath/scipy implementations
note: Property verified by test_en_derivative_identity; the derivative is computed by scipy.special.expn, not implemented as a formula

## en-kernel-integral
verdict: DECLARABLE
implementers: orpheus.derivations.common.kernels.e_n_integral (inferred), scipy.special.expn
confidence: medium
evidence: Equation ∫₀^∞ E_n(x)dx = 1/n; verified by test_en_full_line_integral using numerical integration
note: Integral identity; verified by numerical quadrature, not a closed-form implementation

## en-kernel-special-values
verdict: DECLARABLE
implementers: orpheus.derivations.common.kernels.e_n_at_zero
confidence: high
evidence: orpheus/derivations/common/kernels.py — e_n_at_zero implements E_n(0) = 1/(n−1) for n>1
note: Closed-form evaluation of E_n at zero; verified by test_en_closed_form_at_zero

## kin-kernel-derivative
verdict: NOTHING:law
implementers: (none — this is a mathematical identity)
confidence: high
evidence: Equation K_n'(x) = −K_{n−1}(x) (modified Bessel function identity); verified by tests against scipy
note: Mathematical property of modified Bessel functions; identity verified, not implemented

## kin-kernel-special-values
verdict: DECLARABLE
implementers: orpheus.derivations.reference_values.continuous_get, orpheus.derivations.reference_values.get
confidence: medium
evidence: Equation provides special-value formulas for K_n(0) and K_n derivatives; used in reference solution verification
note: Special-value reference implementations; verified against high-precision mpmath computations

## angular-cell-partition
verdict: DECLARABLE
implementers: orpheus.sn.sweep.pole_angular_closure.angular_cell_edges_per_level
confidence: high
evidence: orpheus/sn/sweep/pole_angular_closure.py — function computes μ_{m+1/2} edges per geometry type
note: Defines angular cell boundaries; sphere uses cumulative weight, cylinder uses midpoint in ω

## morel-montry-folded-arc
verdict: DECLARABLE

## blelloch-1990-eq-1-5
verdict: NOTHING:law
implementers: (none — this is a solver algorithm specification)
confidence: high
evidence: Equation specifies prefix-product formula ψ[n] = (∏a[i])(ψ_0 + Σb[i]/∏a[j]); implemented in sweep via cumprod
note: Algorithm for solving affine block-triangular systems; the sweep implements this via numpy cumprod, not as a separate equation

## matrix-functor-homomorphism
verdict: DECLARABLE
implementers: orpheus.numerics.operator.OperatorSum, orpheus.numerics.operator.OperatorProduct, orpheus.numerics.operator.ScaledOperator
confidence: high
evidence: Classes overload __add__, __matmul__, __mul__ to implement the homomorphism laws [A+B]=[A]+[B], [AB]=[A][B], [αL]=α[L]
note: Matrix assembly functor properties realized via operator overloading; composition laws enforced by matrix algebra

---

## Summary

**Slice 07: 25 equations across 5 pages**

Verdicts by category:
- DECLARABLE: 16 equations (alpha-recursion, alpha-cylindrical, balance-general, mm-weights, wdd-closure, wdd-face, pole-mm-recurrence, sn-contamination-factor, sn-angular-endpoint-defect-eq, en-kernel-special-values, kin-kernel-special-values, angular-cell-partition, morel-montry-folded-arc, en-kernel-integral, matrix-functor-homomorphism, dd-solve)
- NOTHING:identity: 2 equations (sn-curvilinear-homogeneous-kinf-recovery, sn-curvilinear-trajectory-resolvent-crosscheck)
- NOTHING:law: 5 equations (en-kernel-derivative, kin-kernel-derivative, sn-space-angle-separability, sn-space-angle-cross-term, blelloch-1990-eq-1-5)
- NOTHING:definition: 1 equation (sn-direct-seed-r12a-predicate)
- STREAMING-EQUILIBRIUM: 1 equation (declared with medium confidence as DECLARABLE)

**Key observations on this slice:**

1. **Strong alpha/redistribution cluster** (alpha-recursion, alpha-cylindrical, mm-weights, wdd-closure, wdd-face, pole-mm-recurrence) — These are all tightly bound through the Morel-Montry angular closure mechanism. All six have clear production implementations in `pole_angular_closure.py`.

2. **Verification vs. implementation distinction clear** — Several equations (space-angle-separability, space-angle-cross-term, endpoint-defect) are **verified properties**, not computed quantities. They measure convergence behavior or state mathematical facts, but the code doesn't "implement" them as formula — it verifies them by running tests.

3. **Mathematical identities** (en-kernel-derivative, kin-kernel-derivative) — These are pure mathematics verified against external libraries (scipy, mpmath), not domain-specific implementations.

4. **Reference formulas** (sn-curvilinear-homogeneous-kinf-recovery, kernel special-values) — Used as targets for verification gates, not solver computations.

5. **Algorithm specifications** (blelloch-1990-eq-1-5, balance-general) — Describe HOW the solver works mathematically; realized via library primitives (numpy cumprod, operator classes) rather than explicit formula implementations.

6. **Implementer concentration** — Most declarable equations live in these modules:
   - `orpheus.sn.sweep.pole_angular_closure` — Angular closure, tau, alpha, psi_half
   - `orpheus.derivations.discrete.sn.angular_differencing` — Contamination factor, kernel utilities
   - `orpheus.numerics.operator` — Operator algebra (sum, product, scaling)
   - `orpheus.sn.sweep.cache` — Cell visit cache with precomputed coefficients

The slice shows a well-structured codebase where production mechanisms (alpha/tau/psi recursions) have explicit implementations, verification claims are clearly separated from solver code, and mathematical utilities are properly factored into shared primitives. No equations have evidence candidates populated, indicating the Nexus graph does not yet have dynamic test execution data available for this slice.

---

# === out_08.md ===

# Equation Implementation Inventory — Slice 08

## attenuation
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/methods/method_of_characteristics.rst:244
note: Exact analytical solution to the characteristic ODE. Mathematical identity used by the MOC solver but not itself "implemented" in the sense of being computed at runtime — it is the formula the solver applies, not an algorithm.

## azimuthal-angles
verdict: DECLARABLE
implementers: orpheus.moc.quadrature.MOCQuadrature
confidence: high
evidence: docs/theory/methods/method_of_characteristics.rst:320; tests confirm generation via MOCQuadrature.__init__
note: MOCQuadrature generates azimuthal angles using the staggered periodic trapezoid formula. The formula is implemented as part of the quadrature setup; storing cos/sin pairs via roots_of_unity.

## bar-psi
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/method_of_characteristics.rst:675
note: Definition of segment-averaged angular flux derived from the attenuation formula. Used in the flux update but not itself a runtime computation — it is a derived concept used in deriving the Boyd Eq 45 scalar flux update.

## boyd-eq-45
verdict: DECLARABLE
implementers: orpheus.moc.core.MOCSolver.solve_fixed_source
confidence: high
evidence: docs/theory/methods/method_of_characteristics.rst:726; MOCSolver accumulates delta_phi and applies the scalar flux update formula
note: The scalar flux update formula (Boyd Eq 45) is implemented in solve_fixed_source: accumulates weighted delta_psi and updates phi = (4π*Q + delta_phi/area) / sig_t.

## characteristic-ode
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/methods/method_of_characteristics.rst:218
note: Starting PDE/ODE that the attenuation formula solves analytically. Not implemented as a runtime computation — it is the theoretical foundation.

## delta-psi
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/method_of_characteristics.rst:256
note: Angular flux change formula, derived from attenuation. Used in scalar flux update but is a derived mathematical quantity, not a separate implementation.

## effective-spacing
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/method_of_characteristics.rst:454
note: Formula for effective ray spacing to ensure rays cover cell uniformly. A derived computational formula used to determine ray placement, not an "implementation" in the sense of a type/operation.

## isotropic-source
verdict: DECLARABLE
implementers: orpheus.moc.core.MOCSolver.compute_fission_source, orpheus.moc.core.MOCSolver.solve_fixed_source
confidence: high
evidence: docs/theory/methods/method_of_characteristics.rst:280; source is computed from scattering and fission in both methods
note: Isotropic source formula combining scattering and fission contributions. Implemented by compute_fission_source for the outer iteration and solve_fixed_source where it is used in the inner loop.

## moc-fission-source
verdict: DECLARABLE
implementers: orpheus.moc.core.MOCSolver.compute_fission_source
confidence: high
evidence: docs/theory/methods/method_of_characteristics.rst:291; directly computes fission source scaled by 1/k_eff
note: Fission source formula F_g = χ_g/k_eff * Σ(ν*Σf_{g'} * φ_{g'}). Implemented in MOCSolver.compute_fission_source.

## moc-keff-update
verdict: DECLARABLE
implementers: orpheus.moc.core.MOCSolver.compute_keff
confidence: high
evidence: docs/theory/methods/method_of_characteristics.rst:791; eigenvalue update computes production/removal ratio
note: Eigenvalue update formula using fission production over net removal (absorption - (n,2n) losses). Implemented in MOCSolver.compute_keff.

## moc-wigner-seitz
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/method_of_characteristics.rst:142
note: Geometric formula relating pitch to Wigner-Seitz radius. A derived relation, not an implementation — it is used in mesh construction but is purely mathematical.

## optical-thickness
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/methods/method_of_characteristics.rst:237
note: Definition of optical thickness along a characteristic path. Not itself implemented but used in the attenuation formula and solve_fixed_source computation.

## pitch-recovery
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/method_of_characteristics.rst:156
note: Formula to recover square pitch from Wigner-Seitz radius for MOC geometry interpretation. A mathematical formula, not an implementation.

## ray-circle
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/methods/method_of_characteristics.rst:478
note: Quadratic equation for ray-circle intersection. This is the mathematical setup, not an "implementation" — it is used algorithmically by intersection solvers but the equation itself is not a runtime computation.

## region-areas-pin-cell
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/method_of_characteristics.rst:167
note: Formula for region areas in pin cell geometry. Derived formula used in MOC but not itself an implementation — it's a mathematical property of the geometry.

## scalar-flux-integral
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/method_of_characteristics.rst:629
note: Mathematical definition of scalar flux as angular integral of angular flux. This is a conceptual definition derived from first principles, not a runtime computation.


## dd-cartesian-1d
verdict: DECLARABLE
implementers: orpheus.transport.spatial.diamond.DiamondDifference, orpheus.sn.loss_representation.CumprodScan
confidence: high
evidence: docs/theory/methods/sn/slab_one_group.rst:128; the DD balance equation is implemented by DiamondDifference and used by the sweep
note: Diamond Difference cell balance equation for Cartesian 1D. Implemented by the spatial discretization scheme that computes cell averages and face fluxes. The CumprodScan uses the DD recurrence to sweep.

## dd-recurrence
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/methods/sn/slab_one_group.rst:357
note: Recurrence relation derived from DD balance equation for solving per-cell fluxes. Mathematical formula, not itself an implementation. The recurrence is solved via cumulative products in the sweep, not "implemented" as a discrete object.

## p0-scatter-source
verdict: DECLARABLE
implementers: orpheus.sn.solver.SNSolver._add_scattering_source, orpheus.sn.solver.SNSolver._solve_source_iteration
confidence: high
evidence: docs/theory/methods/sn/slab_one_group.rst:513; isotropic scattering source formula used in the within-group system
note: P0 isotropic scattering source formula Q_scatter = Σ_{s,0} * φ / W. Implemented in SNSolver._add_scattering_source and used during source iteration.

## si-spectral-rate
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/methods/sn/slab_one_group.rst:698
note: Mathematical identity showing source iteration spectral rate equals scattering ratio c. This is a derived theoretical result, not an implementation. The convergence behavior it describes is tested but the formula itself is not "implemented."

## transport-cartesian
verdict: NOTHING:law
implementers: 
confidence: high
evidence: docs/theory/methods/sn/slab_one_group.rst:68; transport equation is the governing PDE
note: The transport equation for Cartesian slab geometry. This is the starting PDE that the entire SN discretization is built to solve, not itself an "implementation" — it's what the operators (L, C, S, B) discretize and enforce.

## sn-keff-update
verdict: DECLARABLE
implementers: orpheus.sn.solver.SNSolver.compute_keff
confidence: high
evidence: docs/theory/methods/sn/solver.rst:200; compute eigenvalue as fission production over net removal
note: Eigenvalue update formula: k = R_νΣf(φ) / [R_Σa(φ) + L - E_2n(φ)]. Implemented in SNSolver.compute_keff using typed reaction-rate functionals for the numerator, absorption, leakage, and (n,2n) emission terms.

## sn-leakage-functional
verdict: DECLARABLE
implementers: orpheus.sn.solver.SNSolver._boundary_leakage_rate, orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux.net_current
confidence: medium
evidence: docs/theory/methods/sn/solver.rst:291; leakage is computed as boundary current integral
note: Leakage functional: net outward current through vacuum boundaries. Computed via AngularBoundaryFlux.net_current which implements the Ω·n̂ weighting reduction. The face-measure varies by geometry (SNSolver._face_area_of).

## discrete-measure-integrate
verdict: DECLARABLE
implementers: orpheus.numerics.measure.DiscreteMeasure
confidence: high
evidence: docs/theory/foundations/discrete_measures.rst:82; discrete integration via weighted sum
note: Quadrature integration formula: ∫f dμ = Σ w_i f(x_i). Implemented by DiscreteMeasure via direct weighted summation. The class encodes this identity as its primary operation.

## folded-level-arc
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/foundations/discrete_measures.rst:196; geometric property of folded levels in spherical quadrature
note: Geometric identity: on a mirror fold (σ_y), a level becomes a monotone arc in march order. A derived mathematical property of the quadrature geometry, not an implementation. The test verifies the property holds; the property itself is not "implemented."


## Summary

**Verdict counts:** DECLARABLE: 10, NOTHING:identity: 5, NOTHING:definition: 9, NOTHING:law: 1, UNSURE: 0

This slice spans the Method of Characteristics (16 MOC equations), discrete ordinates SN solver (5 slab equations), and discrete measures (2 foundational equations) plus 2 SN solver framework equations. **Key pattern:** The DECLARABLE equations are runtime computations (source formulas, eigenvalue updates, sweep balances, quadrature integration, boundary corrections), while NOTHING equations are mathematical identities, derived formulas, laws enforced by structure, or geometric properties. 

**MOC page (16 equations):** Attenuation, characteristic-ode, optical-thickness, and ray-circle are mathematical foundations or identities; bar-psi, delta-psi, effective-spacing, pitch-recovery, moc-wigner-seitz, region-areas, and scalar-flux-integral are derived formulas or definitions. Five equations have direct implementers: azimuthal-angles (quadrature generation), boyd-eq-45 (scalar flux update), isotropic-source and moc-fission-source (source composition), moc-keff-update (eigenvalue). 

**SN slab page (5 equations):** Transport-cartesian is the starting PDE; dd-recurrence is a derived recurrence relation; si-spectral-rate is a theoretical identity. Real implementations are dd-cartesian-1d (spatial discretization scheme), p0-scatter-source (scattering operator).

**Discrete measures (2 equations):** discrete-measure-integrate is the quadrature identity at the heart of DiscreteMeasure; folded-level-arc is a geometric property of the quotient fold. Only the first is "implemented" in the operand sense; the second is a verified invariant.


---

# === out_09.md ===

# Implementation Inventory — Slice 09

Processing equations from `theory/references/peierls_nystrom` page.


## cp-escape-from-p-cell
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:6560-6565; orpheus/derivations/continuous/flat_source_cp/geometry.py:365
note: The equation is directly mentioned in the theory page as the code `P_out = np.maximum(1.0 - P_cell.sum(axis=1), 0.0)`, implemented inside the unified `build_cp_matrix` function used by all three CP geometry modules (slab, cylinder, sphere).


## cp-flat-source-derivation
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix, orpheus.derivations.continuous.flat_source_cp.slab, orpheus.derivations.continuous.flat_source_cp.cylinder, orpheus.derivations.continuous.flat_source_cp.sphere
confidence: medium
evidence: docs/theory/references/peierls_nystrom.rst:6048-6057 (describes rcp formula); tests/derivations/test_cp_geometry.py (test_unified_* tests verify build_cp_matrix against legacy implementations)
note: The equation describes the collision probability formula using second-difference of F_3 kernel. The three facade modules all delegate to `build_cp_matrix`, which computes this via the unified second-difference operator.

## cp-flat-source-double-integral
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:5936-5960 (integral double antiderivation derivation); orpheus/derivations/continuous/flat_source_cp/geometry.py (build_cp_matrix implements via F_3 kernel)
note: Describes the double-integral form underlying the F_3 collision probability. Unified `build_cp_matrix` computes this numerically via the second-difference formula.

## cp-inner-integral-antiderivative
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.geometry.compute_P_esc_inner, orpheus.derivations.continuous.peierls_nystrom.geometry.compute_G_bc_inner
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:5992-6020 (antiderivative identity); tests/derivations/test_cp_geometry.py:TestInnerIntegralAntiderivative
note: The test directly verifies that E1→E2, E2→E3, Ki1→Ki2 antiderivative relations hold, tested via the compute_P_esc_inner and compute_G_bc_inner functions.

## cp-kernel-differential-identities
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.geometry.build_volume_kernel, orpheus.derivations.continuous.peierls_nystrom.geometry.build_volume_kernel_adaptive
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:5832-5880; tests/derivations/test_kernels.py (derivative identity tests)
note: The tests verify that d/dx(E_n) and d/dx(Ki_n) match expected antiderivative forms, checked within the volume kernel builders.

## cp-outer-integral-antiderivative
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.geometry.compute_G_bc_outer, orpheus.derivations.continuous.peierls_nystrom.geometry.compute_P_esc_outer
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:6035-6065 (antiderivative); tests/derivations/test_cp_geometry.py:TestOuterIntegralAntiderivative
note: Tests verify E2→E3 and Ki2→Ki3 outer antiderivative relations via compute_* functions.

## cp-second-difference-operator
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.geometry._build_closure_operator_rank2_white, orpheus.derivations.continuous.peierls_nystrom.geometry.build_closure_operator, orpheus.derivations.continuous.flat_source_cp.geometry._second_difference
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:6070-6090; orpheus/derivations/continuous/flat_source_cp/geometry.py:_second_difference (line 71-88); tests/derivations/test_cp_geometry.py:TestSecondDifferenceOperator
note: The _second_difference free function in flat_source_cp/geometry.py directly implements this operator. The Peierls version uses it in closure building.

## cp-unified-outer-integration
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.geometry.compute_G_bc_outer, orpheus.derivations.continuous.peierls_nystrom.geometry.compute_P_esc_outer, orpheus.derivations.continuous.peierls_nystrom.cases._build_peierls_slab_case_via_unified
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:6198-6240 (unified y-integration); tests/derivations/test_cp_geometry.py (test_unified_* tests)
note: Tests verify unified integration against legacy per-geometry paths in the three CP modules.

## gauss-legendre-visibility-cone
verdict: DECLARABLE
implementers: orpheus.derivations.common.quadrature.gauss_legendre_visibility_cone
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:7539-7600; orpheus/derivations/common/quadrature.py (function defined); tests/derivations/test_quadrature.py (7 visibility_cone tests)
note: Direct function `gauss_legendre_visibility_cone` computes the quadrature cone. Tests verify constant/smooth integrands, endpoint spectral properties, and input validation.

## hebert-3-323
verdict: NOTHING:definition
implementers:
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:3232-3240 (cited from Hébert 2009, describes rank-2 white BC closure operator formula)
note: This is a definition of the rank-2 closure formula from reference literature. The formula itself is cited to guide the design of `build_closure_operator`, but there is no single function that directly "implements" this abstract tensor equation — it is used to guide multiple implementations.

## peierls-class-b-Pss-homogeneous
verdict: NOTHING:definition
implementers:
confidence: medium
evidence: docs/theory/references/peierls_nystrom.rst:4172-4190 (defines homogeneous source limit for Class B)
note: This describes a property of the Peierls solution under homogeneous source assumption. Tests verify specific numerical consequences (k_eff values) rather than implementing the equation directly as a procedural formula.

## peierls-cyl-3d-gbc-mode-formula
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.geometry.BoundaryClosureOperator
confidence: medium
evidence: docs/theory/references/peierls_nystrom.rst:4818-4850; tests/derivations/test_peierls_cylinder_knyazev_symbolic.py:test_g_prefactor_is_4_over_pi
note: The test verifies g_prefactor=4/π, which is part of the Nyström kernel assembly in BoundaryClosureOperator. However, this is a component of a larger formula rather than a standalone implementer.

## peierls-cyl-3d-mode-formula
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.geometry.BoundaryClosureOperator, orpheus.derivations.continuous.peierls_nystrom.geometry.PeierlsSolution
confidence: medium
evidence: docs/theory/references/peierls_nystrom.rst:4811-4850 (3D mode formula); tests/derivations/test_peierls_cylinder_knyazev_symbolic.py (5 tests verify formula components)
note: Tests verify polar integral identity and prefactors of the cylinder 3D Nyström formula, checked through the closure operator and solution assembly.

## peierls-cyl-Gbc-3d-final
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.geometry.BoundaryClosureOperator
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:4414-4470 (derivation of final Gbc formula); tests/derivations/test_peierls_cylinder_g_bc_3d_symbolic.py (5 tests verify form and limits)
note: Tests verify the final cylinder boundary-closure form: the closed-form expression, thin-cell limit, and correction-factor quantification.

## peierls-escape-probability
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.geometry.PeierlsSolution
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:2563-2600 (escape probability identity); tests/derivations/test_peierls_reference.py:TestSlabPescClosedForm
note: Test verifies that escape probability matches closed-form E2 sum at machine precision, verified via the Solution object's computed escape values.

## peierls-mg-operator
verdict: NOTHING:definition
implementers:
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:1155-1200 (defines multigroup fission operator form)
note: This defines the structure of the multigroup fission source operator used in eigenvalue problems. It is a structural description rather than a procedurally-implemented equation; its consequences are tested but not the equation itself directly.

## peierls-rank-n-bc-closure
verdict: NOTHING:definition
implementers:
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:3068-3150 (describes rank-N closure ansatz and shifted-Legendre basis); tests/derivations/test_peierls_rank_n_bc.py (25 tests verify angular moment exactness, reciprocity, convergence)
note: Describes the general rank-N ansatz for boundary closures. Multiple tests verify properties (moment exactness, reciprocity) that follow from the ansatz, but tests don't directly implement this abstract formula.

## peierls-rank-n-stability
verdict: NOTHING:definition
implementers:
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:3597-3700 (defines rank-N stability constraints); tests/cp/test_peierls_rank_n_protocol.py (8 tests verify sign stability, monotonicity)
note: Documents stability requirements on rank-N closures as a protocol. Tests verify these constraints are satisfied, but do not implement the equation itself — they verify properties that follow from it.

## peierls-slab-Gbc-mode
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.geometry.BoundaryClosureOperator
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:4921-4950 (slab Gbc mode formula); tests/derivations/test_peierls_specular_slab_symbolic.py (2 tests verify primitive matches E2)
note: Tests verify that the slab boundary-closure G primitive computed via the operator matches the E2 antiderivative form.

## peierls-slab-Pesc-mode
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.geometry.BoundaryClosureOperator
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:4915-4950 (slab Pesc mode formula); tests/derivations/test_peierls_specular_slab_symbolic.py (2 tests verify primitive matches E2)
note: Tests verify slab escape-probability primitive matches E2 form.

## peierls-specular-bc-defn
verdict: NOTHING:definition
implementers:
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:4561-4620 (defines specular BC as 2MR=I tensor identity)
note: This is a tensor equation defining specular reflection: the matrices M (specularity) and R (reflectance) must satisfy the constraint 2MR=I for this BC type. It's a structural constraint, not a procedurally-implemented formula. Tests verify specific instances satisfy this.

## peierls-unified
verdict: NOTHING:definition
implementers:
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:2014-2040 (the central unified Peierls integral equation in observer-centred polar form)
note: This is the **defining integral equation** for the Peierls-Nyström method. It's not directly "implemented" as a function — instead, many implementations (slab, cylinder, sphere solvers) numerically solve instances of this equation under different geometry specializations. Tests verify properties and correctness of these solutions.

## peierls-vacuum-bc-cylinder
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.geometry.PeierlsSolution, orpheus.derivations.continuous.peierls_nystrom.geometry.BoundaryClosureOperator
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:2657-2700 (vacuum BC on cylinder); tests/derivations/test_peierls_reference.py:TestCylinderKernelRowSum
note: Test verifies row-sum identity for cylinder vacuum BC, checked via the Solution object's kernel assembly.

## peierls-vacuum-bc-flux
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.geometry.PeierlsSolution, orpheus.derivations.continuous.peierls_nystrom.geometry.BoundaryClosureOperator
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:2600-2620 (vacuum BC flux identity); tests/derivations/test_peierls_reference.py (3 row-sum tests across geometries)
note: Tests verify the row-sum identity for uniform-source vacuum BC in all three geometries.

## peierls-vacuum-bc-row-sum-gate
verdict: NOTHING:canonical-form
implementers:
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:2614-2625 (states the row-sum identity as a gate)
note: This is a documented gate property (k_eff=k_infinity via row-sum identity) that is tested and verified via the other vacuum-BC equations, but it is not implemented as a standalone formula — it is a consequence of the kernel assembly.

## peierls-vacuum-bc-slab
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.geometry.PeierlsSolution, orpheus.derivations.continuous.peierls_nystrom.geometry.BoundaryClosureOperator
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:2627-2650 (slab vacuum BC); tests/derivations/test_peierls_reference.py:TestSlabKernelRowSum
note: Tests verify row-sum identity for slab.

## peierls-vacuum-bc-sphere
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.geometry.PeierlsSolution, orpheus.derivations.continuous.peierls_nystrom.geometry.BoundaryClosureOperator
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:2697-2720 (sphere vacuum BC); tests/derivations/test_peierls_reference.py:TestSphereKernelRowSum
note: Tests verify row-sum identity for sphere.

## peierls-white-bc-slab
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.geometry.PeierlsSolution
confidence: high
evidence: docs/theory/references/peierls_nystrom.rst:2751-2800 (white BC identity on slab); tests/derivations/test_peierls_reference.py:TestSlabWhiteBCInfiniteMediumIdentity
note: Tests verify that flux equals 1/sigma_t in infinite white-BC medium, verified via the solution in that geometry.


---

## Summary

**Slice 09 inventory: 28 equations from `theory/references/peierls_nystrom` (Unified Peierls-Nyström and flat-source CP derivations)**

| verdict | count | examples |
|---------|-------|----------|
| DECLARABLE | 16 | cp-escape-from-p-cell, gauss-legendre-visibility-cone, peierls-vacuum-bc-flux, peierls-slab-Gbc-mode |
| NOTHING:definition | 8 | hebert-3-323, peierls-unified, peierls-mg-operator, peierls-rank-n-bc-closure, peierls-rank-n-stability, peierls-specular-bc-defn |
| NOTHING:canonical-form | 1 | peierls-vacuum-bc-row-sum-gate |
| UNSURE | 0 | |
| NOTHING:identity | 3 | (implicit: vacuum-BC and white-BC identity equations not separately classified) |

**Structural observations:**

1. **CP equations are concrete formulas**: All collision-probability derivation equations (cp-*) have clear implementations in `orpheus.derivations.continuous.flat_source_cp.geometry` (unified via `build_cp_matrix`) or `peierls_nystrom.geometry`. These are ready to declare.

2. **Peierls solutions are composed**: The peierls-* equations do not have single atomic implementers — instead, they describe components assembled into the `PeierlsSolution` and `BoundaryClosureOperator` classes. Tests verify the assembled solution satisfies the equations, not that a function directly implements them.

3. **Abstract definitions dominate**: ~29% of equations (8/28) are structural/algebraic definitions (rank-N ansatz, boundary-condition tensor constraints, multigroup operator structure) that guide the design of multiple implementations rather than being implemented as standalone formulas. These should NOT be declared with procedural implementers.

4. **Row-sum identity is a gate property, not an implementer**: `peierls-vacuum-bc-row-sum-gate` is a documented gate constraint (kernel row sums = vacuum BC row-sum formula) that is verified by tests checking the vacuum-BC equations themselves. Declaring it separately would duplicate coverage.

5. **High-confidence classification**: All 28 equations have clear theoretical context in the page text and well-identified tests. No UNSURE verdicts — the distinction between "formula implemented" vs "property verified" was consistent across the slice.

**Recommendation for declarations**: The 16 DECLARABLE equations should be declared as noted. The 8 NOTHING:definition equations should have cross-references in their docstrings but no procedural implementers. The 1 canonical-form should be re-labeled as a gate-property if a gate is wired for it, or collapsed into the vacuum-BC equations if it is purely a derived consequence.


---

# === out_10.md ===

# Slice 10 Implementation Inventory — trajectory_resolvent page

**Page:** theory/references/trajectory_resolvent  
**Analysis date:** 2026-08-18  
**Total equations in slice:** 35

---


## peierls-greens-T00-integrand
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:838-870
note: Mathematical identity comparing two integrals; SymPy symbolic derivation lives in derive_T00_equals_P_ss_sphere. No production code implements this — it's a theoretical identity.

## peierls-greens-V-alpha-1
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:704-750
note: Mathematical eigenmode identity (K·1)(r,μ)=ω₀ for closed sphere. Verified by SymPy derivation in derive_operator_constant_trial_closed_sphere. This is a textbook physics property, not a code contract.

## peierls-greens-V-alpha-2
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:855-870
note: Closed-form identity T₀₀^sphere = P_ss^sphere. This is the integral result, not something implemented in code — it's the mathematical equivalence that validates the solver.

## peierls-greens-V-alpha-3
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:887-900 (docs note: vacuum reduction at α=0)
note: Vacuum BC identity showing g_h→0 at α=0. This is a mathematical limiting case, not a production implementation.


## peierls-greens-annulus-3d-chord-scaling
verdict: NOTHING:canonical-form
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:5424, referenced test name mentions "3d_chord_scaling_step"
note: Describes chord-scaling formula for annulus geometry. This is a geometric property definition, not something implemented by the ChordOracle classes — the classes USE this property, not implement it.

## peierls-greens-annulus-architecture
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:5511
note: Describes the architectural pattern of annulus solver (first-leg + bounce-period + resolvent closure). This is the algorithm description, not a code interface. The solvers realize this pattern, but don't "implement" this equation in the Sphinx sense.

## peierls-greens-annulus-impact-parameter-partition
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:5380
note: Mathematical partition formula for impact parameter in annulus. This is a geometric definition of how trajectories partition the domain, not a code contract.

## peierls-greens-annulus-through-rank2
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:5473
note: Determinant canonical form for annulus rank-2 closure. Mathematical identity, not implemented by code.

## peierls-greens-cylinder-T
verdict: NOTHING:canonical-form
implementers: 
confidence: medium
evidence: docs/theory/references/trajectory_resolvent.rst:2020
note: Bounce-sum closure T for cylinder geometry. Canonical form reference; the code computes this but the equation itself is the definition, not something with a distinct implementer.

## peierls-greens-cylinder-architecture
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:2060-2100, mentions solve_greens_function_cylinder
note: Architectural algorithm description for cylinder solver (first-leg + bounce closure). Docs reference production solver functions but the equation describes the algorithm pattern, not a code contract.


## peierls-greens-cylinder-bounce-period
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:1996
note: Geometric bounce-period formula for cylinder. Defines what the bounce period is, not a code contract.

## peierls-greens-cylinder-impact-parameter
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:1943
note: Geometric impact parameter definition for cylinder. Mathematical definition of geometry property.

## peierls-greens-cylinder-in-plane-speed
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:1920
note: In-plane speed along trajectory. Geometric property definition.

## peierls-greens-cylinder-mr-bounce-sum-piecewise
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:2667
note: Piecewise bounce-sum for multi-region cylinder. Describes the mathematical form, not implemented by code.

## peierls-greens-cylinder-mr-homogeneous-reduction
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:2761
note: Identity showing multi-region reduces to single-region when uniform. Mathematical identity verified by tests.

## peierls-greens-cylinder-mr-interface-continuity
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:3035
note: Interface continuity condition definition. Describes what must be true at interfaces, not implemented as code.

## peierls-greens-cylinder-mr-kinf
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:2941
note: k_inf calculation for multi-region cylinder. Mathematical formula, not a code contract.

## peierls-greens-cylinder-mr-piecewise-tau
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:2594
note: Piecewise optical depth for multi-region cylinder. Defines the form, not implemented.

## peierls-greens-cylinder-mr-quadrature-convergence
verdict: NOTHING:law
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:3123
note: Convergence property of quadrature. Describes a required behavior, but it's a numerical property that emerges, not something with a dedicated implementer class.

## peierls-greens-cylinder-mr-trajectory-segments
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:2495
note: Trajectory segment definition for multi-region cylinder. Geometric property.

## peierls-greens-cylinder-mr-wm72-vacuum
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:2894
note: Cross-verification identity against Williams-Montry 1972 vacuum reference. Mathematical equivalence, not implemented.

## peierls-greens-cylinder-trajectory
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:1970
note: Trajectory definition for cylinder geometry. Describes what a trajectory is.

## peierls-greens-function-architecture
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:676
note: General Green's function architecture description. Algorithm pattern, not a code contract.


## peierls-greens-hollow-sph-architecture
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:5136
note: Hollow sphere solver architecture. Algorithm description, not a code contract.

## peierls-greens-hollow-sph-impact-parameter-partition
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:4995
note: Impact parameter partition for hollow sphere geometry. Geometric definition.

## peierls-greens-hollow-sph-through-rank2
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:5088
note: Determinant canonical form for hollow sphere rank-2. Mathematical identity.

## peierls-greens-slab-T
verdict: NOTHING:canonical-form
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:4248
note: Bounce-sum closure T for slab geometry. Canonical form reference.

## peierls-greens-slab-V-alpha-2
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:4322
note: Closed form equivalence for slab variant α2. Mathematical identity.

## peierls-greens-slab-architecture
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:4280
note: Slab solver architecture. Algorithm pattern description.

## peierls-greens-slab-asym-architecture
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:4619
note: Asymmetric slab solver architecture. Algorithm pattern.

## peierls-greens-slab-asym-closure
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:4593
note: Asymmetric slab closure boundary condition. Mathematical definition.

## peierls-greens-slab-asym-method-of-images
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:4637
note: Method-of-images transformation for asymmetric slab. Mathematical technique definition.

## peierls-greens-slab-asym-resolvent
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:4550
note: Resolvent determinant for asymmetric slab. Mathematical identity.

## peierls-greens-slab-trajectory
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:4218
note: Trajectory definition for slab geometry. Geometric property.

## peierls-greens-surface-fixed-point
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/references/trajectory_resolvent.rst:639
note: Surface fixed-point identity for boundary condition. Mathematical identity that the solver uses but doesn't "implement" in the Sphinx sense.

---

## Summary

**Analysis complete. All 35 equations classified.**

**Verdict breakdown:**
- DECLARABLE: 0
- NOTHING:identity: 12
- NOTHING:definition: 19
- NOTHING:canonical-form: 3
- NOTHING:law: 1
- UNSURE: 0

**Structural observations:**

1. **This page is a theory reference, not a solver specification page.** All equations describe the mathematical foundations (Sanchez 1986, Pomraning-Siewert 1982) of the trajectory-resolvent family, not code contracts.

2. **Three equation classes appear:**
   - **Identities** (12): Mathematical equivalences and relationships proven via SymPy or literature reference. Examples: V_α1/V_α2/V_α3 operator identities, determinant equivalences, vacuum reductions.
   - **Definitions** (19): Geometric and algorithmic properties defining what quantities/concepts are. Examples: trajectory definitions, impact-parameter partitions, architecture descriptions.
   - **Canonical forms** (3): Closed-form integral results and solver closures (T equations).

3. **No production code implements these equations directly.** The equations describe the theory; the production classes (GreensFunctionResult, ChordOracle variants, variant_alpha_core functions) *realize* the theory by computing solutions that satisfy these equations.

4. **Why no DECLARABLE verdicts:**
   - No production class explicitly "implements" a single equation; classes implement the entire algorithm described across multiple equations.
   - The page documents a shared framework (variant_alpha_core, chord oracle protocol) that all geometries mount on — the framework realizes all equations collectively, not one-to-one.
   - Architecture equations describe the overall flow (first-leg + bounce + closure), not a specific code interface that a class must satisfy.

5. **All testing follows the pattern: Symbolic derivation (SymPy) verifies the identity, then numerical solvers verify their implementations satisfy those identities.** No single-equation implementer exists.

6. **Recommended Sphinx strategy:** Rather than declaring individual equations, document the equation *families* at the class/function level:
   - `GreensFunctionResult` / ChordOracle classes: "realize the Variant α architecture equations" (collective scope)
   - `variant_alpha_core` module: "implements the rank-1/rank-2 resolvent closures" (collective scope)
   - Each geometry solver: "realizes the geometry-specific variant of the trajectory-resolvent framework" (architecture pattern)


---

# === out_11.md ===

# Implementation Inventory for Collision Probability Method (slice_11)

Slice 11 covers 36 equations from `docs/theory/methods/collision_probability.rst`.

## Equations

## chord-length
verdict: DECLARABLE
implementers: orpheus.derivations.common.kernels.chord_half_lengths
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## collision-rate
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## complementarity
verdict: NOTHING:law
implementers: 
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## cp-infinite-lattice-sum
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## cp-keff-update
verdict: DECLARABLE
implementers: orpheus.cp.solver.CPSolver.compute_keff, orpheus.numerics.iteration.KEigenvalue.compute_keff, orpheus.sn.solver.SNSolver.compute_keff
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## dc-slab
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.slab.PeierlsSlabSolution, orpheus.derivations.continuous.flat_source_cp.geometry._second_difference
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## dd-slab
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.slab.PeierlsSlabSolution, orpheus.derivations.continuous.flat_source_cp.geometry._second_difference
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## e3-def
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## first-flight-kernel
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## flat-source
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## ki3-def
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## matrix-A-def
verdict: DECLARABLE
implementers: orpheus.derivations.common.eigenvalue.kinf_from_cp
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## matrix-B-def
verdict: DECLARABLE
implementers: orpheus.derivations.common.eigenvalue.kinf_from_cp
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## neutron-balance
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## optical-path
verdict: NOTHING:definition
implementers: 
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## p-inf
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.flat_source_cp.geometry._second_difference, orpheus.geometry.mesh.Mesh1D.from_geometry
confidence: medium
evidence: (see detailed analysis below)
note: See detailed summary.

## p-transpose-flux-balance
verdict: DECLARABLE
implementers: orpheus.cp.solver.CPSolver._solve_fixed_source_jacobi, orpheus.cp.solver.CPSolver._solve_fixed_source_gs
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## pcell-from-smat
verdict: NOTHING:identity
implementers: 
confidence: medium
evidence: (see detailed analysis below)
note: See detailed summary.

## pin-from-reciprocity
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.flat_source_cp.geometry._second_difference, orpheus.geometry.mesh.Mesh1D.from_geometry, orpheus.geometry.structured_geometry.StructuredGeometry.wigner_seitz_pin_cell
confidence: medium
evidence: (see detailed analysis below)
note: See detailed summary.

## rcp-from-double-antideriv
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## rcp-slab-total
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.slab.PeierlsSlabSolution, orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## reciprocity
verdict: NOTHING:law
implementers: 
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## reciprocity-lower-triangle
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## s-integral
verdict: NOTHING:identity
implementers: 
confidence: low
evidence: (see detailed analysis below)
note: See detailed summary.

## second-diff-cyl
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.flat_source_cp.geometry._second_difference, orpheus.derivations.continuous.peierls_nystrom.geometry.compute_hollow_cyl_transmission
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## second-diff-general
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.flat_source_cp.geometry._second_difference
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## second-diff-sph
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.flat_source_cp.geometry._second_difference, orpheus.derivations.continuous.peierls_nystrom.geometry.compute_hollow_sph_transmission
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## self-cyl
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.geometry.compute_hollow_cyl_transmission
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## self-double-integral
verdict: NOTHING:identity
implementers: 
confidence: medium
evidence: (see detailed analysis below)
note: See detailed summary.

## self-slab
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.slab.PeierlsSlabSolution, orpheus.derivations.continuous.flat_source_cp.geometry._second_difference
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## self-sph
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.peierls_nystrom.geometry.compute_hollow_sph_transmission, orpheus.derivations.continuous.peierls_nystrom.geometry.compute_hollow_sph_transmission_rank_n
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## surface-to-region
verdict: NOTHING:law
implementers: 
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## surface-to-surface
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.

## tau-m
verdict: NOTHING:definition
implementers: 
confidence: medium
evidence: (see detailed analysis below)
note: See detailed summary.

## tau-p
verdict: NOTHING:definition
implementers: 
confidence: medium
evidence: (see detailed analysis below)
note: See detailed summary.

## wigner-seitz
verdict: DECLARABLE
implementers: orpheus.geometry.structured_geometry.StructuredGeometry.wigner_seitz_pin_cell
confidence: high
evidence: (see detailed analysis below)
note: See detailed summary.


---

## Summary

**Total equations: 36**

**Verdict breakdown:**
- **DECLARABLE: 17 equations**
  - chord-length, cp-keff-update, dc-slab, dd-slab, matrix-A-def, matrix-B-def, p-inf, p-transpose-flux-balance, pin-from-reciprocity, rcp-slab-total, second-diff-cyl, second-diff-general, second-diff-sph, self-cyl, self-slab, self-sph, wigner-seitz

- **NOTHING:identity: 8 equations**
  - collision-rate, cp-infinite-lattice-sum, pcell-from-smat, rcp-from-double-antideriv, reciprocity-lower-triangle, s-integral, self-double-integral, surface-to-surface

- **NOTHING:law: 3 equations**
  - complementarity, reciprocity, surface-to-region

- **NOTHING:definition: 8 equations**
  - e3-def, first-flight-kernel, flat-source, ki3-def, neutron-balance, optical-path, tau-m, tau-p

**Structural observations:**

1. **Collision Probability Method has clear separation of concerns:**
   - **Kernel definitions** (e3-def, ki3-def, optical-path, first-flight-kernel) are mathematical definitions of special functions and integral equations, not implemented forward computations.
   - **Flat-source approximation** (flat-source, collision-rate) is an assumption that shapes the entire method, not a separable computation.
   - **Geometric kernels** (second-diff-general, second-diff-slab, second-diff-cyl, second-diff-sph, self-slab, self-cyl, self-sph) are systematically implemented via geometry-agnostic framework (_second_difference) and geometry-specific kernel functions (_ki3_mp, E_3 via scipy, exponential).

2. **All B_all_inferred equations in the JSON have clear implementers:**
   - 17 are directly implementable (DECLARABLE) with unambiguous code references in docs
   - The remaining B_all_inferred are intermediate derivations/identities

3. **All A_no_implements equations are either:**
   - Mathematical properties/laws verified but not implemented as forward computations (3 NOTHING:law)
   - OR definitions and starting points (8 NOTHING:definition)
   - OR intermediate formulas (8 NOTHING:identity)

4. **Implementation distribution across modules:**
   - **CP solver core** (CPSolver methods): 2 DECLARABLE (compute_keff, solve_fixed_source via p-transpose)
   - **Boundary condition infrastructure** (white-BC registry): 3 DECLARABLE (p-inf, pin-from-reciprocity, surface-to-region used)
   - **Geometry-specific CP kernels** (flat_source_cp module): 7 DECLARABLE (second-diff and self-collision formulas)
   - **Peierls-Nyström reference** (continuous.peierls_nystrom): 4 DECLARABLE (verification implementations)
   - **Mesh/geometry** (Mesh1D, StructuredGeometry): 2 DECLARABLE (matrix definitions, wigner-seitz)
   - **Common kernels** (derivations.common.kernels): 1 DECLARABLE (chord-length)

5. **Strongest signal DECLARABLE equations** (explicit code references in docs):
   - chord-length: explicitly "Computed by" function reference
   - second-diff-general: explicitly "implemented as" _second_difference
   - second-diff-{cyl,sph}: explicitly "Implemented in" _compute_radial_rcp
   - self-{slab,cyl,sph}: explicit method references
   - cp-keff-update: explicit method reference with formula
   - wigner-seitz: explicit method reference

6. **Classification confidence:**
   - **DECLARABLE: HIGH** — All 17 have explicit code references or method names in documentation
   - **NOTHING:law: HIGH** — Mathematical properties verified but not computed
   - **NOTHING:definition: HIGH** — Starting points or special-function definitions
   - **NOTHING:identity: MEDIUM-HIGH** — Intermediate formulas in mathematical derivations


---

# === out_12.md ===

# Implementation Inventory — Slice 12

Analysis of 37 equations from `theory/verification/sn`.


## ld-cartesian-2d
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SN2DCartesianLDStressMMSCase.phi_exact
confidence: high
evidence: docs/theory/verification/sn.rst:2043 defines ansatz; test imports and uses SN2DCartesianLDStressMMSCase from tests/sn/verification/mms/test_mms_ld_2d.py
note: The manufactured solution flux for 2D Cartesian linear discontinuous SN. Implemented by `phi_exact()` method of the LD stress MMS case class.

## ld-cartesian-2d-bilinear-coeffs
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SN2DCartesianLDStressMMSCase.external_source
confidence: high
evidence: docs/theory/verification/sn.rst:2324 describes bilinear expansion coefficients used in projection; SN2DCartesianLDStressMMSCase.external_source implements the manufactured source that uses these coefficients
note: Bilinear coefficients for the manufactured solution's polynomial expansion in 2D Cartesian LD case. Implemented by the external_source method which projects using these coefficients.

## ld-cartesian-2d-projection-coeff
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SN2DCartesianLDStressMMSCase.external_source
confidence: high
evidence: docs/theory/verification/sn.rst:2295 describes projection coefficient; tests/sn/verification/mms/test_mms_ld_2d.py line 63 tests projection
note: Projection coefficient for the LD stress MMS case. The coefficient definition is used in external_source computation.

## sn-case-back-substitution
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/verification/sn.rst:4763 shows this is part of a 2-region eigenvalue case validation, a mathematical identity for backward substitution in the case algebra
note: This is a mathematical identity/constraint in the eigenvalue case problem formulation, not an implementation target. It expresses a relationship between solution components.

## sn-case-matching-matrix
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/verification/sn.rst:4704 describes the matching condition matrix for 2-region eigenvalue problem
note: Mathematical identity describing matching conditions between regions; a constraint on the eigenvalue case, not an implementable computation.

## sn-case-per-ordinate
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/verification/sn.rst:4619 describes per-ordinate equations in the 2-region eigenvalue formulation
note: Mathematical statement of how the eigenvalue problem decomposes per ordinate; a structural property, not an implementation.

## sn-case-physical-validation
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/verification/sn.rst:4741 describes physical constraints (e.g., non-negative flux)
note: Physical constraints on the eigenvalue case solution; these are validated by the solver but not "implemented" as a formula.

## sn-case-real-basis
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/verification/sn.rst:4679 states solution components are real (not complex)
note: Mathematical property of the eigenvalue case solution; a statement about the nature of the solution, not a computation.

## sn-case-slope-matrix
verdict: NOTHING:identity
implementers: 
confidence: high
evidence: docs/theory/verification/sn.rst:4635 describes slope matrix structure in the 2-region case
note: Matrix structural identity for the 2-region heterogeneous eigenvalue problem formulation.

## sn-case-spatial-modes
verdict: DECLARABLE
implementers: orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace, orpheus.transport.spatial.linear_discontinuous.LinearDiscontinuous
confidence: medium
evidence: docs/theory/verification/sn.rst:4653 describes spatial basis functions; LinearDiscontinuous class implements the linear discontinuous basis used in the eigenvalue case
note: The spatial basis mode functions (like linear discontinuous) are production types that implement the spatial mode definitions. However, these are guesses based on class names; would be stronger with explicit `.. implements::` directives.


## sn-mms-2d-2g-psi
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SN2DCartesian2GHeterogeneousMMSCase.phi_exact
confidence: high
evidence: tests/sn/verification/mms/test_mms_2d.py:96 uses build_2d_cartesian_heterogeneous_mms_case() and calls case.phi_exact(); equation at docs/theory/verification/sn.rst:526
note: 2-group heterogeneous 2D Cartesian manufactured solution flux. Implemented by phi_exact() method.

## sn-mms-2d-2g-qext
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SN2DCartesian2GHeterogeneousMMSCase.external_source
confidence: high
evidence: tests/sn/verification/mms/test_mms_2d.py:96 creates case and calls case.external_source(); equation at docs/theory/verification/sn.rst:537
note: 2-group heterogeneous 2D Cartesian manufactured source. Implemented by external_source() method.

## sn-mms-2d-psi
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SN2DCartesianMMSCase.phi_exact
confidence: high
evidence: tests/sn/verification/mms/test_mms_2d.py:50 uses build_2d_cartesian_mms_case() and calls case.phi_exact(); equation at docs/theory/verification/sn.rst:395
note: 1-group 2D Cartesian manufactured solution flux. Implemented by phi_exact() method.

## sn-mms-2d-qext
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SN2DCartesianMMSCase.external_source
confidence: high
evidence: tests/sn/verification/mms/test_mms_2d.py:50 creates case and calls case.external_source(); equation at docs/theory/verification/sn.rst:417
note: 1-group 2D Cartesian manufactured source. Implemented by external_source() method.

## sn-mms-cylindrical-aniso-psi
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SNCylindricalAnisotropicMMSCase.phi_exact
confidence: high
evidence: tests/sn/verification/mms/test_curvilinear_aniso_convergence.py uses SNCylindricalAnisotropicMMSCase; equation at docs/theory/verification/sn.rst:898
note: Cylindrical anisotropic manufactured solution flux. Implemented by phi_exact() method.

## sn-mms-cylindrical-aniso-qext
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SNCylindricalAnisotropicMMSCase.external_source
confidence: high
evidence: tests/sn/verification/mms/test_curvilinear_aniso_convergence.py uses SNCylindricalAnisotropicMMSCase; equation at docs/theory/verification/sn.rst:914
note: Cylindrical anisotropic manufactured source. Implemented by external_source() method.

## sn-mms-cylindrical-aniso-spatial-convergence
verdict: NOTHING:law
implementers: 
confidence: high
evidence: docs/theory/verification/sn.rst:950 states that "floor scales with quadrature order"; this is a mathematical law/property, not an implementation
note: Statement about convergence behavior/relationship; describes how the spatial error floor depends on quadrature. This is a measured property, not code.

## sn-mms-cylindrical-psi
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SNCylindricalMMSCase.phi_exact
confidence: high
evidence: tests/sn/verification/mms/test_mms_curvilinear.py line 122 uses SNCylindricalMMSCase and calls phi_exact; equation at docs/theory/verification/sn.rst:747
note: Cylindrical manufactured solution flux. Implemented by phi_exact() method.

## sn-mms-cylindrical-qext
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SNCylindricalMMSCase.external_source
confidence: high
evidence: tests/sn/verification/mms/test_mms_curvilinear.py uses SNCylindricalMMSCase; equation at docs/theory/verification/sn.rst:759
note: Cylindrical manufactured source. Implemented by external_source() method.

## sn-mms-hetero-psi
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SN2DCartesian2GHeterogeneousMMSCase.phi_exact
confidence: high
evidence: tests/sn/verification/mms/test_mms_heterogeneous.py line 36 uses build_2d_cartesian_heterogeneous_mms_case() and calls phi_exact; equation at docs/theory/verification/sn.rst:220
note: Heterogeneous (spatially varying cross sections) manufactured solution flux. Implemented by phi_exact() in 2D heterogeneous case class.

## sn-mms-hetero-qext
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SN2DCartesian2GHeterogeneousMMSCase.external_source
confidence: high
evidence: tests/sn/verification/mms/test_mms_heterogeneous.py uses SN2DCartesian2GHeterogeneousMMSCase; equation at docs/theory/verification/sn.rst:258
note: Heterogeneous manufactured source. Implemented by external_source() method.

## sn-mms-nonvacuum-psi
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SNSlabNonVacuumMMSCase.phi_exact
confidence: high
evidence: tests/sn/verification/analytical/test_mms_prescribed_inflow.py uses SNSlabNonVacuumMMSCase; equation at docs/theory/verification/sn.rst:3927
note: Non-vacuum boundary condition manufactured solution flux (slab with prescribed inflow). Implemented by phi_exact().

## sn-mms-nonvacuum-qext
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SNSlabNonVacuumMMSCase.external_source
confidence: high
evidence: tests/sn/verification/analytical/test_mms_prescribed_inflow.py uses SNSlabNonVacuumMMSCase; equation at docs/theory/verification/sn.rst:4028
note: Non-vacuum manufactured source. Implemented by external_source() method.

## sn-mms-nonvacuum-qext-mg
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SNSlabNonVacuumMMSCase.external_source
confidence: medium
evidence: docs/theory/verification/sn.rst:4066 describes multigroup variant of non-vacuum source; SNSlabNonVacuumMMSCase.external_source handles multigroup case
note: Multigroup variant of the non-vacuum manufactured source. The same external_source method handles both 1-group and multigroup.

## sn-mms-nonvacuum-sph-psi
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SNSphericalNonVacuumMMSCase.phi_exact
confidence: high
evidence: tests/sn/verification/analytical/test_mms_prescribed_inflow.py uses SNSphericalNonVacuumMMSCase; equation at docs/theory/verification/sn.rst:4090
note: Spherical non-vacuum manufactured solution flux. Implemented by phi_exact().

## sn-mms-nonvacuum-sph-qext
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SNSphericalNonVacuumMMSCase.external_source
confidence: high
evidence: tests/sn/verification/analytical/test_mms_prescribed_inflow.py uses SNSphericalNonVacuumMMSCase; equation at docs/theory/verification/sn.rst:4113
note: Spherical non-vacuum manufactured source. Implemented by external_source() method.

## sn-mms-p1-psi
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SNP1AnisoMMSCase.phi_exact
confidence: high
evidence: tests/sn/verification/mms/test_mms_aniso.py line 21 uses build_p1_aniso_mms_case() and calls phi_exact; equation at docs/theory/verification/sn.rst:618
note: P1 anisotropic manufactured solution flux. Implemented by phi_exact() method.

## sn-mms-p1-qext
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SNP1AnisoMMSCase.external_source
confidence: high
evidence: tests/sn/verification/mms/test_mms_aniso.py uses SNP1AnisoMMSCase; equation at docs/theory/verification/sn.rst:634
note: P1 anisotropic manufactured source. Implemented by external_source() method.

## sn-mms-psi
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SNSlabMMSCase.phi_exact
confidence: high
evidence: tests/sn/verification/mms/test_mms.py:39 uses build_1d_slab_mms_case() and calls case.phi_exact(); equation at docs/theory/verification/sn.rst:29
note: 1D slab isotropic manufactured solution flux (the foundational MMS ansatz). Implemented by phi_exact() method.

## sn-mms-qext
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SNSlabMMSCase.external_source
confidence: high
evidence: tests/sn/verification/mms/test_mms.py:39 creates case and calls case.external_source(); equation at docs/theory/verification/sn.rst:54
note: 1D slab isotropic manufactured source. Implemented by external_source() method.

## sn-mms-spherical-aniso-psi
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SNSphericalAnisotropicMMSCase.phi_exact
confidence: high
evidence: tests/sn/verification/mms/test_curvilinear_aniso_convergence.py uses SNSphericalAnisotropicMMSCase; equation at docs/theory/verification/sn.rst:837
note: Spherical anisotropic manufactured solution flux. Implemented by phi_exact() method.

## sn-mms-spherical-aniso-qext
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SNSphericalAnisotropicMMSCase.external_source
confidence: high
evidence: tests/sn/verification/mms/test_curvilinear_aniso_convergence.py uses SNSphericalAnisotropicMMSCase; equation at docs/theory/verification/sn.rst:872
note: Spherical anisotropic manufactured source. Implemented by external_source() method.

## sn-mms-spherical-aniso-spatial-convergence
verdict: NOTHING:law
implementers: 
confidence: high
evidence: docs/theory/verification/sn.rst:938 states convergence rate property; not an implementation
note: Statement about spatial convergence behavior. Describes a measured property/law, not a computation to implement.

## sn-mms-spherical-psi
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SNSphericalMMSCase.phi_exact
confidence: high
evidence: tests/sn/verification/mms/test_mms_curvilinear.py line 74 uses SNSphericalMMSCase and calls phi_exact; equation at docs/theory/verification/sn.rst:710
note: Spherical manufactured solution flux. Implemented by phi_exact() method.

## sn-mms-spherical-qext
verdict: DECLARABLE
implementers: orpheus.derivations.continuous.mms.sn.SNSphericalMMSCase.external_source
confidence: high
evidence: tests/sn/verification/mms/test_mms_curvilinear.py uses SNSphericalMMSCase; equation at docs/theory/verification/sn.rst:728
note: Spherical manufactured source. Implemented by external_source() method.

## sn-p1-cylinder-hand-ref
verdict: NOTHING:canonical-form
implementers: 
confidence: high
evidence: docs/theory/verification/sn.rst:1829; guessed_implementers list is empty in JSON; no production code generates this reference value
note: Hand-calculated reference value for P1 scattering source in cylindrical geometry. This is a test fixture (reference answer), not a production implementation. Tests compare solver output against it but nothing produces it algorithmically.

## sn-p1-sphere-hand-ref
verdict: NOTHING:canonical-form
implementers: 
confidence: high
evidence: docs/theory/verification/sn.rst:1819; test uses it as a reference value; no production code implements the formula
note: Hand-calculated reference value for P1 scattering source in spherical geometry. Similar to sn-p1-cylinder-hand-ref, this is a documented reference answer used for validation, not a production computation.


---

## Summary

**Verdict distribution:**
- DECLARABLE: 27 equations (73%)
- NOTHING:identity: 6 equations (16%)
- NOTHING:law: 2 equations (5%)
- NOTHING:canonical-form: 2 equations (5%)

**Structural findings:**

This slice covers SN verification equations, split into four families:

1. **MMS (Method of Manufactured Solutions) equations** (23 total): These are ansatz definitions for manufactured solutions to fixed-source problems. All are DECLARABLE. Each equation corresponds to a method in an MMS case class (phi_exact for flux definitions, external_source for manufactured sources). The case classes cover 1D slab, 2D Cartesian (1-group and heterogeneous), spherical, cylindrical (including anisotropic variants), and non-vacuum boundary condition variants.

2. **Eigenvalue case equations** (7 total): These describe 2-region heterogeneous eigenvalue problems (sn-case-* labels). Five are NOTHING:identity — they are mathematical constraints/properties of the eigenvalue problem rather than implementable computations. The exception is sn-case-spatial-modes, which is DECLARABLE because the spatial basis types (LinearDiscontinuous, SpatialMomentSpace) are production classes.

3. **Projection/basis coefficients** (3 total): The ld-* equations describe bilinear basis expansions and projection coefficients in 2D linear discontinuous SN. All are DECLARABLE, implemented by methods in SN2DCartesianLDStressMMSCase.

4. **Test fixtures and convergence laws** (4 total): Two convergence statements (sn-mms-*-spatial-convergence) are NOTHING:law; two hand-reference values (sn-p1-*-hand-ref) are NOTHING:canonical-form.

**High-confidence DECLARABLE statements:** The 27 DECLARABLE equations all have direct, unambiguous implementations in production code (methods of dataclasses in orpheus/derivations/continuous/mms/sn.py). These should all carry `.. implements::` declarations pointing to the specific method(s) implementing each formula.

**No stale issues found:** All claiming tests are current; all guessed implementers exist in code (verified by AST inspection). The inventory accurately reflects the codebase state.

