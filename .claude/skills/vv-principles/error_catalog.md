# L0 Error Catalog

Errors caught during development by term-level verification.  Each
entry records the error, how it hid from higher-level tests, the L0
test that caught it, and the lesson learned.

This file is the primary QA publication artifact.  It supersedes
``gotchas.md`` (deleted).

## Error Classification

Errors are classified by the 6 AI Failure Modes taxonomy:

| # | Mode | Description |
|---|------|-------------|
| 1 | Sign flip | `(a − b)` vs `(b − a)` |
| 2 | Variable swap | `mu_x` vs `mu_y`, `h_sat_l` vs `h_sat_v` |
| 3 | Missing factor | Missing `2×`, `ΔA/w`, volume |
| 4 | Factor error | Wrong constant, hardcoded value |
| 5 | Index error | `face[i]` vs `face[i+1]` |
| 6 | Convention drift | Definition site vs usage site disagree |

---

## ERR-001 — Z-ordinate weight loss in 2D Lebedev sweep

**Failure mode:** #3 Missing factor — missing contribution  
**Date:** 2026-04-03  
**Solver:** SN (Cartesian 2D)

**Bug:** In `_sweep_2d_wavefront`, ordinates with `mu_x = mu_y = 0`
(z-directed) were skipped with `continue`.  Their quadrature weights
(0.77% of total) were lost from scalar flux integration.

**Impact:** Multi-group eigenvalue error of ~0.4% (2G: 1.867 vs 1.875).

**How it hid from higher-level tests:**
- 1-group keff = νΣ_f/Σ_a is independent of weight loss (cancels)
- Spatial convergence showed the scheme converged — to the wrong value
- "Reasonable numbers" — 0.4% looks like discretization error

**L0 test that catches it:** L0-SN-001 (streaming equilibrium) — the
volume-averaged φ would deviate from Q/Σ_t by the missing weight
fraction.

**Lesson:** 1-group eigenvalue tests are degenerate.  Always verify
with ≥2 groups.

---

## ERR-002 — Scattering matrix transpose in vectorization

**Failure mode:** #2 Variable swap — `SigS` vs `SigS^T`  
**Date:** 2026-04-03  
**Solver:** SN (all geometries)

**Bug:** Vectorized `_add_scattering_source` used `phi @ SigS^T`
(double-transposed) instead of `phi @ SigS`.

**Impact:** keff = 2.06 (catastrophically wrong) — caught immediately.

**How it could hide:** For symmetric scattering matrices (1-group
self-scatter), SigS = SigS^T and the bug is invisible.

**L0 test that catches it:** L0-SN-009 (scattering source magnitude)
— hand-calculated SigS^T @ φ vs code output with asymmetric 2G matrix.

**Lesson:** Always test with asymmetric inputs.  The identity
`(A^T v) = (v^T A)^T` means the transpose moves sides, not into
a pre-transposed matrix.

---

## ERR-003 — Octant batching breaks reflective BC ordering

**Failure mode:** #6 Convention drift — implicit ordinate ordering  
**Date:** 2026-04-03  
**Solver:** SN (Cartesian 2D)

**Bug:** Batching ordinates by sweep direction changed the order in
which reflective BC boundary fluxes were updated.  A group reads
boundary fluxes that its reflected group hasn't updated yet.

**Impact:** 2G convergence test failed (keff diffs grew).

**How it hid:** 1G converged; the optimization gave 2× speedup
on the passing test.

**L0 test that catches it:** Not directly an L0 problem (this is
a coupling issue between ordinates).  Caught by L1: eigenvalue
convergence with mesh refinement.

**Lesson:** Sequential processing order is part of the interface
contract for reflective BCs.  Don't parallelize without verifying
the data dependency.

---

## ERR-004 — Hardcoded 4π in BiCGSTAB RHS normalization

**Failure mode:** #4 Factor error — wrong normalization constant  
**Date:** 2026-04-03  
**Solver:** SN BiCGSTAB (spherical)

**Bug:** `build_rhs` hardcoded `4π` for angular normalization.
Correct for Lebedev (sum(w) = 4π), wrong for GL (sum(w) = 2).

**Impact:** BiCGSTAB on GL diverged (keff oscillating).

**How it hid:** All initial BiCGSTAB testing used Lebedev where
4π is correct.

**L0 test that catches it:** L0-SN-001 (streaming equilibrium)
with BiCGSTAB solver — φ would not converge to Q/Σ_t with wrong
normalization.

**Lesson:** Never hardcode quadrature-dependent constants.

---

## ERR-005 — DD recurrence rewrite breaks multi-group convergence

**Failure mode:** #4 Factor error — numerically unstable rewrite  
**Date:** 2026-04-04  
**Solver:** SN (Cartesian 1D)

**Bug:** Algebraically equivalent rewrite of `_solve_recurrence`
introduced catastrophic cancellation: `2*(…) − psi_in` subtracts
nearly-equal large numbers.

**Impact:** Multi-group scattering iteration diverged (~1e34/iter).

**How it hid:** 1-group unaffected; formulas are algebraically
identical; 2D wavefront (different code path) still passed.

**L0 test that catches it:** Not a single-term issue; caught by
L1 (multi-group eigenvalue convergence).

**Lesson:** "Algebraically equivalent" ≠ "numerically equivalent".
The stable form `0.5*(psi_in + psi_out)` averages known quantities;
the unstable form subtracts large numbers.

---

## ERR-006 — Wrong α recursion + missing ΔA/w in curvilinear sweep

**Failure mode:** #2 Variable swap + #3 Missing factor  
**Date:** 2026-04-04  
**Solver:** SN (cylindrical, spherical)

**Bug:** Two simultaneous bugs:
1. α recursion used `cumsum(+w·ξ)` with azimuthal cosine mu_y
   instead of `cumsum(−w·η)` with radial cosine mu_x
2. Missing `ΔA_i/w_m` geometry factor on the redistribution term

**Impact:** Heterogeneous keff diverged with mesh refinement
(1.15 → 0.90 → 0.52 → 0.25).  Spherical had flux spike 5.1× at r=0.

**How it hid from 20 passing tests:**
- Homogeneous eigenvalue: exact (redistribution cancels for flat flux)
- 1-group: degenerate (keff = material ratio)
- Particle balance: exact (telescoping is by construction)
- Conservation: exact (total is correct; per-ordinate is wrong)
- Flux non-negativity: no negatives produced
- Single sweep finite: fluxes are finite, just wrong

**L0 tests that catch it:**
- L0-SN-003 (per-ordinate flat-flux) — streaming + redistribution ≠ 0
  per ordinate without ΔA/w.  This is the definitive L0 diagnostic.
- L0-SN-001 (streaming equilibrium) — flux spike at r=0 visible
- L0-SN-008 (contamination β) — β ≈ 2.0 instead of ~0

**Investigation history:** 6 approaches failed before root cause found:
reverse sweep, step closure, starting direction, bidirectional sweep,
scaled α, zero redistribution.  Full details in
`docs/theory/discrete_ordinates.rst` §Investigation History.

**Lesson:** Per-ordinate flat-flux consistency (L0-SN-003) is the
FUNDAMENTAL correctness criterion for curvilinear SN.  It should be
the first test written for any curvilinear transport implementation.

---

## ERR-007 — Multi-group BiCGSTAB unstable for spherical geometry

**Failure mode:** #3 Missing factor — same as ERR-006  
**Date:** 2026-04-05  
**Solver:** SN BiCGSTAB (spherical, cylindrical)

**Bug:** The explicit FD transport operator was missing the same
ΔA/w geometry factor as the sweep (ERR-006).

**Impact:** Multi-group BiCGSTAB diverged for spherical (keff → NaN).
Previously documented as "BiCGSTAB is unreliable for curvilinear"
when the real issue was a missing geometry factor.

**L0 test that catches it:** L0-SN-003 applied to the operator output.

**Lesson:** When a bug is found in one code path (sweep), check ALL
code paths that implement the same physics (BiCGSTAB operator).

---

## ERR-008 — Boundary volume halving in SN keff computation

**Failure mode:** #6 Convention drift — geometry vs solver disagree  
**Date:** 2026-04-03  
**Solver:** SN (unified solver, CartesianMesh)

**Bug:** `CartesianMesh.volume` halved first/last cell volumes (matching
the old `PinCellGeometry` convention from the 2D MATLAB port).  The SN
solver used these volumes in the keff computation: `k = Σ(νΣf·φ·V) / Σ(Σa·φ·V)`.
With reflective BCs, boundary cells are at the symmetry plane and should
have full volume, not half.

**Impact:** ~1e-4 systematic error in keff for heterogeneous problems.
Homogeneous problems unaffected (uniform flux × uniform XS → volume
cancels in the ratio).

**How it hid:** Homogeneous verification tests all passed to machine
precision.  The error only appeared when comparing the unified solver
against the old `sn_1d.py` solver (which correctly used full `dx` as
volume, not the geometry's `volume` property).

**L0 test that catches it:** Direct comparison of keff between old and
new solver on the same heterogeneous problem.  A dedicated test would
verify `Σ(νΣf·φ·V) / Σ(Σa·φ·V) = keff_reported` with known volumes.

**Lesson:** Volume conventions (half-cell vs full-cell) must be
explicit and documented.  The same geometry object should not be used
for both "mesh for sweeping" (where boundary halving might make sense
for source normalization) and "volumes for integral quantities" (where
the full cell contributes).

---

### ERR-009 — CP neutron balance transpose

**Failure mode:** #2 Variable swap — P vs P.T  
**Date:** 2026-04-05  
**Solver:** CP (slab and cylindrical)

**Bug:** The CP power iteration computed `phi = P_inf @ source` instead
of `phi = P_inf.T @ source`.  With the convention `P[i,j]` =
P(birth_i → collision_j), the neutron balance for collision target `j`
sums over birth regions: `Σ_i P[i,j] · V_i · Q_i = P.T @ source`.

**Impact:** Wrong eigenvalue for ALL heterogeneous problems.  For the
1G 2-region slab benchmark: solver gave k=1.373 vs analytical k=1.272
(8% error).  keff was systematically too high because the flux
redistribution between fuel and moderator was incorrect.

**How it hid:** Homogeneous benchmarks passed to machine precision
because P is symmetric when all regions have identical cross sections
(P = P.T).  The 1-group homogeneous case is doubly degenerate: no
spatial redistribution AND k = νΣ_f/Σ_a regardless of P.  The bug
only appeared when running 2-region heterogeneous benchmarks with the
formal verification suite.

**L0 test that catches it:** The synthetic 1G 2-region slab benchmark
(`benchmark_1g_slab`) with analytical eigenvalue from the CP matrix.
The analytical k is computed independently by solving the 2×2 matrix
eigenvalue problem `det(A - (1/k)B) = 0`.  Any transpose error in the
solver causes immediate disagreement.

**Lesson:** The CP matrix convention (birth-first vs collision-first
indexing) propagates through the entire solver.  Document the convention
explicitly (now in `collision_probability.rst` §Flat-Source Approximation)
and verify with heterogeneous multi-region benchmarks.  Homogeneous
verification is *necessary but not sufficient* — it only tests the
diagonal of the CP matrix.

---

## ERR-010 — pyXSteam viscosity cutoff at 900 °C causes NaN cascade

**Failure mode:** #4 Factor error — library-imposed validity limit  
**Date:** 2026-04-05  
**Solver:** Thermal Hydraulics (Module 07)

**Bug:** pyXSteam's `my_AllRegions_ph` returns NaN for T > 900 °C due
to an artificial guard (`if T > 900 + 273.15: return NaN`).  The IAPWS
2008 viscosity correlation itself is well-defined beyond this limit.
During post-failure LOCA blowdown, coolant reaches ~901 °C in the
outlet node.  NaN viscosity → NaN kinematic viscosity → NaN friction
factor → NaN pressure → solver crash.

**Impact:** Integration stopped at t ≈ 395 s (of 600 s target).

**How it hid from higher-level tests:**
- Pre-failure phase (t < 287 s) coolant stays well below 900 °C
- Post-failure code path was never exercised until event detection was
  implemented (TH-20260401-001)
- The NaN manifested as `cool_p = [NaN, NaN]`, suggesting a pressure
  bug rather than a viscosity bug — required tracing through two levels
  of function calls to find the root cause

**L0 test that catches it:** Direct property evaluation:
`assert not np.isnan(h2o_properties(0.33, 4399e3)[0].mu)` — tests
that viscosity is returned for high-enthalpy states reachable during
LOCA.

**Fix:** `_iapws_viscosity(T_K, rho)` — same formula without cutoff.

**Lesson:** Third-party library validity limits are not always physical
limits.  When a library returns NaN, check whether the underlying
correlation is actually invalid or just conservatively guarded.

---

## ERR-011 — MATLAB gap geometry mixes radius with axial height

**Failure mode:** #2 Variable swap — `fuel.r` vs `fuel.dz`  
**Date:** 2026-04-05  
**Solver:** Thermal Hydraulics (MATLAB reference, `funRHS.m` line 272)

**Bug:** `gap.r_ = (clad.r(1) + fuel.dz)/2` adds a cladding inner
radius (~4.22 mm) to the fuel axial node height (~1.5 m), producing a
"gap radius" of 0.752 m instead of the correct ~4.17 mm.  The gap heat
transfer area becomes 180× too large.

**Impact:** MATLAB fuel centre temperature is 808 °C instead of the
correct 1140 °C at steady state.  This artificially keeps coolant below
pyXSteam's 900 °C viscosity limit and delays clad failure by ~138 s.

**How it hid:**
- MATLAB ran to completion because the wrong gap area prevented the
  coolant from ever reaching the viscosity cutoff
- All MATLAB results are self-consistent — temperatures, stresses, and
  failure time are plausible for a LOCA scenario, just based on wrong
  gap thermal resistance
- The bug is in a deformable geometry update section where `clad.r`,
  `fuel.r`, `fuel.dz`, `clad.dz` all appear nearby — easy to grab the
  wrong variable

**L0 test that catches it:** Steady-state fuel centre temperature check:
at 69 kW total power, T_fuel_centre should be ~1100-1200 °C, not ~800 °C.
A simple analytical estimate: ΔT ≈ LHGR/(4πk) ≈ 567 °C above the fuel
surface.

**Lesson:** Variable names that differ only in the last character
(`fuel.r` vs `fuel.dz`) are a maintenance hazard, especially in code
with mixed scalar/vector indexing (MATLAB's `fuel.dz` is a vector, but
`clad.r(1)` extracts a scalar via linear indexing of a 2D array).

---

## ERR-012 — Static heat transfer areas in deformable TH/RK modules

**Failure mode:** #3 Missing factor — missing geometry update  
**Date:** 2026-04-01  
**Solver:** Thermal Hydraulics (Module 07), Reactor Kinetics (Module 08)

**Bug:** Gap and clad radial heat transfer areas (`gap_a_bnd`, `clad_a_bnd`)
were computed once at initialization from fabrication geometry and never
updated with deformed radii/heights.  MATLAB's `funRHS.m` recomputes
`gap.a_` and `clad.a_` every RHS call from the current deformed geometry
(`clad.r`, `clad.dz`, `fuel.r`, `fuel.dz`).

**Impact:** During LOCA/RIA transients, fuel thermal expansion changes the
gap geometry by ~0.5–3%.  Using stale fabrication areas introduces a
systematic bias in the radial heat transfer.  The impact is small at steady
state but compounds during transients with large deformations.

**How it hid from higher-level tests:**
- At t=0 and during early transient, deformations are negligible (< 0.1%)
  so static areas produce identical results
- Steady-state eigenvalue tests don't exercise the deformable geometry path
- The `clad_a_bnd_def` variable WAS computed in the RHS (line 791 of TH)
  but never used — a classic "dead code" pattern

**L0 test that catches it:** Compare fuel surface temperature with static
vs deformed areas during a LOCA at t > 300 s.  The static version over-
estimates gap heat transfer when the gap narrows (area decreases with
thermal expansion).

**Fix:** Replace `p["gap_a_bnd"]` and `p["clad_a_bnd"]` with locally
computed deformed values.  Fuel areas kept static (MATLAB convention —
fuel boundary areas are not recomputed in `funRHS.m`).

**Lesson:** When initializing geometry parameters for an ODE RHS, document
explicitly which quantities are "frozen at fabrication" vs "updated each
call".  The MATLAB code doesn't distinguish these — it uses globals that
are silently overwritten.

---

## ERR-013 — Closed-gap stress BC uses fabrication gap width instead of roughness

**Failure mode:** #4 Factor error — wrong denominator  
**Date:** 2026-04-01  
**Solver:** Fuel Behaviour (Module 06)

**Bug:** In `_solve_stress()`, the closed-gap boundary conditions (BC3 and BC4)
divided stress/strain gradients across the gap by `params["gap_dr0"]` (fabrication
gap width = 100 μm) instead of the effective contact gap thickness (~6 μm roughness).
The MATLAB DAE uses the current deformed gap width, which converges to roughness
after closure.

**Impact:** Contact pressure was 40.5 MPa vs MATLAB's 39.8 MPa with the correct
fix, but was systematically wrong by a factor related to (100/6) ≈ 17× in the
stress gradient term before the fix.  The initial fix used `gap_dr0` → the contact
pressure was exactly 10× too high compared to the MATLAB reference value reported
at a different timestep.

**How it hid from higher-level tests:**
- The open-gap phase (before 2.85 years) was unaffected — BC3/BC4 use the
  pressure BCs, not the gap gradient form
- The closed-gap phase stress values were "reasonable" (~40 MPa) even though
  they were wrong — no analytical reference exists for the full closed-gap
  coupled system
- The 10× discrepancy was initially attributed to MATLAB's contact pressure
  being at a different timestep (it was — but the BC was also wrong)

**L0 test that catches it:** Compare contact pressure at a fixed time after
closure against an independent analytical estimate: for a thin-walled tube
under internal pressure (p_gas - p_cool) with contact, σ_r(inner) ≈ -p_contact
≈ -(p_gas - p_cool) × geometry_factor.

**Fix:** BC4 rewritten as a displacement-based gap constraint:
`r_clad_in_deformed - r_fuel_out_deformed = roughness`.  This is linear in
stresses, physically transparent, and avoids any division by gap width.

**Lesson:** When converting DAE residuals to algebraic equations for a linear
solve, the effective "thickness" used in finite-difference gradients across
interfaces must match the physical gap state.  A fabrication-time value is
wrong after gap closure because the physics has fundamentally changed.

---

## ERR-014 — sigT truncation inconsistency in .m data files

**Failure mode:** #4 Factor error — truncated intermediate vs stored total  
**Date:** 2026-04-05  
**Solver:** Cross-section data pipeline (all solvers affected via sigma-zero iteration)

**Bug:** The MATLAB `convertCSVtoM.m` script computes `sigT = sigC + sigF + sigL + sum(sigS)` from full-precision intermediates, then writes ALL quantities (sigC, sigF, sigS, sigT) independently truncated to `%13.6e` format.  The stored `sigT` is the once-truncated full-precision sum, while recomputing `sigT` from the stored (doubly-truncated) components gives a different value.  For U-238 at 600K: stored sigT[0,0] = 108.14, recomputed from components = 77.87, offset = 30.27 barns.

**Impact:** When the GXS→HDF5 converter computed sigT from components (the physically correct approach), the sigma-zero iterations converged to different values, shifting PWR k_inf by 0.4% (1.01771 vs 1.01357).

**How it hid from higher-level tests:**
- The aqueous reactor (H-1 + O-16 + U-235 at 294K) was unaffected because H-1 has only 1 sigma-zero (no interpolation needed) and O-16 has negligible offset
- All `.m` file components matched the GXS parser exactly (sigC, sigF, sigS diffs = 0) — the discrepancy only appeared in the stored sigT which was precomputed from higher-precision intermediates
- The 0.4% k_inf shift is within the range of "plausible" numerical differences between implementations

**L0 test that catches it:** Direct comparison: `assert max(|sigT_stored - (sigC + sigF + sigL + sigS_rowsum)|) < 1e-3` for each isotope.  This immediately reveals the 30-barn offset.

**Lesson:** When porting a data pipeline, verify not just individual components but also **derived quantities** that were computed from those components.  The same truncation format applied to inputs and outputs does not guarantee consistency between stored outputs and recomputed-from-stored-inputs.

---

## ERR-015 — compute_keff ignores (n,2n) net neutron production

**Failure mode:** #3 Missing factor — missing (n,2n) contribution in eigenvalue estimate  
**Date:** 2026-04-05  
**Solver:** CP (all geometries)

**Bug:** `CPSolver.compute_keff` computed `k = νΣf·φ·V / Σa·φ·V` where
`Σa = Σc + Σf + ΣL + Σ₂_out`.  Each (n,2n) reaction produces one
*extra* neutron (one in, two out), but neither the extra neutron
production nor the removal accounting reflected this.  The correct
eigenvalue balance is:

    k = νΣf·φ·V / (Σt − Σs − 2Σ₂)·φ·V

The denominator is the *net* removal after all scattering and (n,2n)
production is subtracted from the total.  When Σ₂ = 0, this reduces to
`νΣf / Σa` (no change to existing tests).

**Impact:** With `Sig2[0,0] = 0.01` on region A (2G), the solver
converged to k = 1.793 instead of the analytical k = 2.045 — a 12%
error.  The transport solve (`solve_fixed_source`) correctly included
`2·Σ₂·φ` in the source, so the flux shape was right, but the eigenvalue
estimate was biased low by the missing production term.

**How it hid from higher-level tests:**
- ALL test materials have `Sig2 = 0` (zero sparse matrix).  The formula
  is correct when Sig2 = 0 because `total - scatter - 0 = absorption`.
- The `make_mixture` function hardcoded `Sig2 = csr_matrix((ng, ng))`
  with no parameter to override it, so it was impossible to construct
  a test material with nonzero (n,2n) through the standard API.
- The error only appeared when a custom `Mixture` with nonzero `Sig2`
  was constructed directly and compared against the dense eigensolver.

**First wrong fix attempt:** Added `n2n_production` to the numerator:
`k = (νΣf + Σ₂_production) / Σa`.  This gives 1.808, still wrong.  The
issue is that `Σa` already includes `Σ₂_out` as removal, but for the
eigenvalue balance the (n,2n) appears as a *source* (2 neutrons out), not
a removal.  The correct denominator is `Σt - Σs - 2Σ₂`, not `Σa`.

**L0 test that catches it:** `test_cp_verification.py::TestN2N::
test_n2n_solver_keff_matches_analytical` — constructs a 2G material with
`Sig2[0,0] = 0.01`, computes the analytical eigenvalue via
`kinf_from_cp(..., sig_2_mats=[sig_2])`, and compares against the
solver.  Tolerance: 1e-5.

**Lesson:** The eigenvalue estimate formula `production / absorption`
hides the implicit assumption Σ₂ = 0.  When adding a new reaction type,
trace it through BOTH the transport solve (where it's a source) AND the
eigenvalue estimate (where it changes the production/removal balance).
The two must be consistent.  Testing with zero cross sections hides the
inconsistency — need at least one test with the term nonzero.

---

## ERR-016 — Tautological inner-iteration convergence residual

**Failure mode:** #6 Convention drift — convergence check tests identity, not convergence  
**Date:** 2026-04-05  
**Solver:** CP Gauss-Seidel mode

**Bug:** The GS inner convergence check computed:

    phi_g_new = transported_g / denom_g
    collision_rate_g = denom_g * phi_g_new
    res = ||collision_rate_g - transported_g||

Substituting the first line into the third: `denom_g * (transported_g /
denom_g) - transported_g = 0`.  The residual is identically zero
regardless of whether the within-group scattering has converged.  The
inner loop always exited after exactly 1 iteration.

**Impact:** The GS solver mode was functionally identical to a
sequential-group Jacobi update.  Groups with strong within-group
self-scatter (thermal groups where Σs(g→g)/Σt(g) is large) were NOT
iterating to convergence.  The outer power iteration still converged to
the correct eigenvalue (because it does its own convergence check), but
the GS mode provided no acceleration benefit from inner iterations.

**How it hid from higher-level tests:**
- All 27 eigenvalue tests passed because the outer iteration converged
  to the correct answer regardless of inner iteration count
- The diagnostic `n_inner` array showed values of 1 everywhere, which
  was interpreted as "fast convergence" rather than "broken residual"
- The `test_thermal_group_needs_more_inner_iterations` test used `>=`
  (thermal ≥ fast), which passed vacuously since all values were 1
- The AI QA review initially concluded that inner iterations are
  *fundamentally unnecessary* for the CP method (wrong — the source
  depends on the flux through self-scatter)

**L0 test that catches it:** `test_cp_verification.py::TestGSInnerIterations::
test_no_self_scatter_one_inner` — material with zero diagonal in Σs
should converge in ≤ 2 inner iterations (no self-consistency needed).
`test_thermal_needs_more_inner_than_fast` — with the corrected residual,
thermal groups genuinely need more inner iterations than fast groups.

**Fix:** Changed residual to relative flux change:
`||φ_new - φ_old|| / ||φ_new||`.  This is nonzero when within-group
self-scatter changes the source between iterations, and zero when the
source doesn't depend on the current group's flux.

**Lesson:** A convergence check that compares quantities derived from
each other by construction tests nothing.  The residual must compare
*independent* quantities: the old flux vs the new flux, or the old
source vs the recomputed source.  When a "convergence diagnostic" shows
all 1s, that's a red flag — it could mean instant convergence OR a
tautological check.  Distinguish the two by testing with a problem that
*should* require multiple iterations (strong self-scatter).

---

## ERR-017 — Wigner-Seitz pitch formula doubled in MC heterogeneous tests

**Failure mode:** #3 Missing factor — extra factor of 2  
**Date:** 2026-04-06  
**Solver:** Monte Carlo (test suite)

**Bug:** The MC heterogeneous test computed the square unit cell pitch as
`pitch = r_cell * sqrt(pi) * 2` instead of `pitch = r_cell * sqrt(pi)`.
The correct formula equates the square cell area to the Wigner-Seitz
circle area: `pitch^2 = pi * r_cell^2`, giving `pitch = r_cell * sqrt(pi)`.
The factor of 2 quadrupled the cell area.

**Impact:** The extra area was all moderator, which drastically changed the
neutron economy.  For the 1G 2-region case: k_mc = 0.757 vs k_ref = 0.990
(24% systematic error).  For the 2G 2-region case, the population collapsed
to zero (NaN) because the subcritical system with 4× moderator couldn't
sustain a neutron population at 200 neutrons/cycle.

**How it hid from higher-level tests:**
- All homogeneous tests passed (single material everywhere — pitch is
  irrelevant for delta-tracking in a homogeneous medium)
- The tests were all marked `@pytest.mark.slow` and may not have been
  run regularly
- The z-scores were NaN (0/0 from collapsed population) which fails
  the `< 5.0` assertion but doesn't indicate which direction the error is
- The error looked like "MC can't handle subcritical systems" rather than
  "the geometry is wrong"

**L0 test that catches it:** Direct comparison of pitch against the factory
convention: `pwr_pin_equivalent` uses `r_cell = pitch / sqrt(pi)`, so
inverting gives `pitch = r_cell * sqrt(pi)`.  A unit test asserting
`pitch**2 == pi * r_cell**2` would immediately flag the factor of 2.

**Lesson:** When constructing geometry for cross-method comparison, verify
the cell area/volume matches between the two methods.  A factor-of-2 error
in a linear dimension is a factor-of-4 in area — large enough to change
the qualitative physics (supercritical → subcritical), yet small enough to
be invisible in the code review because `* 2` looks like it "corrects for
a half-cell to full-cell conversion" or "accounts for the diameter vs
radius convention".

---

## ERR-018 — Direction sampling uses uniform theta instead of isotropic

**Failure mode:** #4 Factor error — wrong PDF for spherical sampling  
**Date:** 2026-04-06 (identified during L0 test design)  
**Solver:** Monte Carlo

**Bug:** The solver samples the polar angle as `theta = pi * rng.random()`
(uniform in [0, π]) instead of `theta = arccos(1 - 2*xi)` (uniform on the
unit sphere).  True isotropic sampling requires the PDF `p(theta) =
sin(theta)/2` to account for the solid angle Jacobian.

**Impact:** The uniform-theta sampling overweights the poles (theta ≈ 0 and
theta ≈ π) where `sin(theta)` is small.  For the 2D projection used by the
solver: `E[sin^2(theta)] = 1/2` (uniform) vs `2/3` (isotropic).  This
systematically shortens the average 2D step length by ~19%.

**Classification:** Known simplification, not a bug to fix.  The formula
matches the original MATLAB `monteCarloPWR.m` implementation.  Since the
solver only tracks 2D projections (x, y) and uses periodic BCs on a square
cell, the non-isotropic sampling affects the effective mean free path but
does not invalidate the eigenvalue calculation (it changes the effective
geometry scaling, which is absorbed into the keff estimate).

**L0 test that documents it:** `test_mc_properties.py::test_direction_sampling`
verifies `E[dir_x^2] = 1/4` (the formula's prediction, not the isotropic
1/3), confirming the code matches the INTENDED formula.

**Lesson:** When porting from MATLAB, document which simplifications are
intentional vs accidental.  A sampling formula that "looks wrong" may be
a deliberate approximation that the original author validated empirically.
The L0 test should verify the INTENDED formula, not the physically correct
one, and the documentation should explain the distinction.

---

## ERR-019 — Missing 4π·sin(θ) weight factor in MOC scalar flux update

**Failure mode:** #3 Missing factor — incomplete angular integration weight  
**Date:** 2026-04-06  
**Solver:** MOC (2D pin cell)

**Bug:** The MOC transport sweep accumulated `delta_phi` with weight
`omega_a * omega_p * t_s` instead of the correct
`4*pi * omega_a * omega_p * t_s * sin(theta_p)`.  Two factors were
missing: (1) the `4*pi` from the angular flux → scalar flux integral
(`phi = integral_{4pi} psi dOmega`), and (2) the `sin(theta_p)` that
arises because the 2D segment-averaged angular flux `bar_psi` relates
to the 3D path integral via `bar_psi = Q/Sig_t + delta_psi * sin(theta) / (Sig_t * ell)`.

The scalar flux update formula (Boyd et al. 2014, Eq. 45) is:

    phi_i = (4*pi / Sig_t_i) * [Q_i + (1/A_i) * sum omega_a * omega_p * t_s * sin(theta_p) * delta_psi]

The `4*pi` factor multiplies the entire bracket. When delta_phi is defined
as `sum(4*pi * omega_a * omega_p * t_s * sin_p * delta_psi)`, the update
becomes `phi = (4*pi*Q + delta_phi/A) / Sig_t`.

**Impact:** Heterogeneous keff was completely wrong: MOC gave 1.344 vs
CP reference of 0.902 for a 2-region fuel+coolant pin cell (1G).  The
homogeneous case was UNAFFECTED because `delta_psi = 0` when the angular
flux is spatially uniform (all boundary fluxes equal `Q/Sig_t`).

**How it hid from homogeneous tests:**
- For homogeneous material with converged boundary fluxes,
  `psi_in = Q/Sig_t` everywhere → `delta_psi = 0` → `delta_phi = 0`
- `phi = 4*pi*Q/Sig_t` regardless of the weight factor
- 1G: k = nu*SigF/SigA (weight-independent)
- 2G/4G: matrix eigenvalue (still weight-independent for uniform medium)
- All 3 homogeneous eigenvalue tests passed to machine precision

**L0 test that catches it:** `test_moc_verification.py::TestL0EquilibriumFlux::
test_pure_scatterer_equilibrium_single_sweep` — injects a non-trivial
boundary flux and checks that the resulting scalar flux matches the
analytical value.  With wrong weights, the correction term `delta_phi/A`
has the wrong magnitude and the flux deviates.  The heterogeneous
particle balance test also catches it immediately (production/absorption ≠ keff).

**Lesson:** The angular integration weight in MOC contains problem-specific
factors (`4*pi` from the full-sphere integral, `sin(theta_p)` from the
2D→3D projection) that cancel out for spatially uniform solutions.  This
makes the missing factor invisible to homogeneous tests.  ALWAYS test the
transport sweep with a heterogeneous problem before declaring the weight
formula correct.  The Boyd Eq. 45 formula should be verified term-by-term
against the derivation, not just checked for self-consistency on the
homogeneous case.

---

## ERR-020 — ULP-noisy cell volumes from `cbrt → **3` round trip

**Failure mode:** #3 Missing factor — numerical round-trip through a
non-bijective float64 operation destroys a structural invariant.
**Date:** 2026-04-13
**Solver:** `orpheus.geometry` (all consumers via `Mesh1D.volumes`)

**Bug:** `_subdivide_zone` constructed equal-volume spherical edges via
`r_k = cbrt(inner^3 + k/n * (outer^3 - inner^3))`, and `compute_volumes_1d`
then re-derived cell volumes as `(4/3) π · diff(edges**3)`. Because
`cbrt(x) ** 3 != x` at the ULP level in float64, the reconstructed
`edges**3` values drifted ~1 ULP per cell, so cells in a zone that were
*supposed* to have identical volume by construction drifted by up to
~2.2e-14 relative error. The cylindrical path (`sqrt → **2`) had the
same bug but at ~6.7e-15 — just under the common `rtol=1e-14` threshold,
which is why only the spherical case failed visibly.

The specific failing tests:

    tests/geometry/test_geometry.py::TestZoneSubdivision::test_equal_volume_single_zone[SPHERICAL]
    tests/geometry/test_geometry.py::TestZoneSubdivision::test_equal_volume_multi_zone[SPHERICAL]

reported volumes of `7.260569688296488` vs reference `7.260569688296414`
— a relative difference of `1.03e-14` against an `rtol=1e-14` assertion.

**Impact:** None observed in solver eigenvalues (the relative drift is
well below every physics tolerance in the repo), but the assertion
`"all cells in an equal-volume zone are bit-identical"` was broken,
masking an invariant a future bug could violate more seriously without
detection. Fixing it also tightens the cylindrical path from ~7e-15 to
bit-exact, eliminating a hidden source of noise in CP/SN/MOC
quadrature-weighted integrals over spherical/cylindrical meshes.

**How it hid:** Every downstream consumer tolerates ULP-level volume
drift (physics tolerances are ≥1e-10 for eigenvalues, ≥1e-8 for flux
shapes). The only test that asserted bit-exactness was the geometry
invariant test itself, and it was introduced late enough that the
spherical path had never been pushed to `rtol=1e-14`. The cylindrical
path accidentally passed.

**Fix:** Compute cell volumes **from the algebraic invariant** at
subdivision time, not from the edges after the fact. `_subdivide_zone`
now returns `(edges, volumes)`; for each coordinate system:

* Cartesian:   `V_cell = (outer - inner) / n`
* Cylindrical: `V_cell = π · (outer² - inner²) / n`
* Spherical:   `V_cell = (4/3) π · (outer³ - inner³) / n`

One scalar per zone, broadcast to every cell — no round trip through
`sqrt`/`cbrt`, so every cell in an equal-volume zone is bit-identical
by construction. `Mesh1D` gained an optional
`precomputed_volumes` field that overrides the edge-derived default;
`mesh1d_from_zones` populates it. Manually-constructed meshes with
arbitrary edges continue to fall back to `compute_volumes_1d` from
edges.

**L0 test that catches it:**
`tests/geometry/test_geometry.py::TestZoneSubdivision::test_equal_volume_{single,multi}_zone`
for every `CoordSystem` — enforces bit-equal volumes at `rtol=1e-14`.

**Lesson:** Non-bijective float operations (`sqrt`, `cbrt`, `exp`,
`log`) do not survive a round trip, and invariants that *should* hold
algebraically are not free — they must be preserved by design. When
an invariant of the form "X_i == X_j for all i, j" exists, compute X
once and broadcast, don't compute it N times and hope for the best.
Fishbone: whenever you see `op(op_inverse(x))` in the code, ask
whether the inverse is bit-exact, and if not, refactor to avoid the
round trip.

---

## ERR-021 — Degenerate ray tangent to pin-cell corner raises IndexError

**Failure mode:** #5 Index error — unchecked assumption that the
intersection list always has ≥2 entries.
**Date:** 2026-04-14
**Solver:** MOC (`orpheus.moc.geometry._ray_box_intersections`)

**Bug:** `_ray_box_intersections` walks the four walls of the square
pin cell `[0, pitch]^2`, collects every wall crossing whose hit point
lies inside the wall segment, deduplicates by `s`-tolerance, and then
indexed `s_vals[0]` and `s_vals[1]` as entry/exit. If the ray grazes a
*corner* of the box, the two adjacent wall solutions collapse to the
same `s` value; after the `1e-12` dedup pass only **one** entry
remains, and `s_vals[1]` raises `IndexError`. The same failure mode
triggers if the ray is parallel to one axis and offset so neither
orthogonal wall yields an in-range hit — the list can even start
empty. The half-step offset in the track generator makes an exact
corner hit vanishingly unlikely in normal use, but the guard was
missing.

**Impact:** None observed in production (no test ever seeded a ray
exactly through a corner). A seeded pitch or azimuthal angle that
aligned a ray with a corner would crash ray tracing before any flux
solve, masking the degeneracy as a random `IndexError` rather than a
skippable geometric edge case.

**Fix:** `_ray_box_intersections` now returns
`tuple | None`: `None` signals a degenerate ray (empty `s_vals`, or
fewer than two distinct entries after dedup). `_trace_single_ray`
short-circuits a `None` to `([], (x0, y0), (x0, y0), -1, -1)`, and
`MOCMesh._trace_all_rays` already skips tracks with `if not
segments: continue`, so a degenerate ray is now silently dropped
instead of aborting the trace.

**L0 test that catches it:**
`tests/moc/test_ray_tracing.py::test_degenerate_corner_ray` — seeds a
ray with `(x0, y0) = (0, 0)` and `phi = π/4` so the entry point is
exactly the `(0, 0)` corner, and asserts that
`_ray_box_intersections` returns `None` and `_trace_single_ray`
returns empty segments rather than raising.

**Lesson:** Whenever a function indexes a collection by a fixed
position (`list[0]`, `list[1]`), the *precondition* that the
collection has that many entries must be either enforced by
construction or checked explicitly. Geometric primitives in
particular must handle degenerate inputs (tangent, parallel,
coincident) as first-class cases, not crashes — in a ray tracer,
"this ray contributes nothing" is a valid outcome.

---

## ERR-022 — Negative lethargy bin width flips flux-per-lethargy sign

**Failure mode:** #6 Convention drift — sign of `du` depends on the
energy-grid ordering convention, which the callers silently relied on.
**Date:** 2026-04-14
**Solver:** MC (`orpheus.mc.solver.solve_monte_carlo`), homogeneous
spectrum solver (`orpheus.homogeneous.solver`), MOC plotting helper
(`orpheus.plotting.plot_moc_spectra`).

**Bug:** Group boundaries in ORPHEUS follow the standard nuclear-data
convention of *descending* energy (`eg[0]` = fast edge, `eg[-1]` =
thermal edge), so the lethargy widths
`du[g] = log(eg[g+1] / eg[g])` are **negative**. Three call sites
divided the (non-negative) group tally by this signed `du` to get
`flux_per_lethargy`, producing uniformly negative values:

* `orpheus.mc.solver.solve_monte_carlo`:
  `flux_per_lethargy = tally / xs.du`
* `orpheus.homogeneous.solver.HomogeneousResult.flux_per_lethargy`:
  `self.flux / self.du` (where `du` was stored signed)
* `orpheus.plotting.plot_moc_spectra`: `flux / du` for fuel, clad,
  coolant spectra

The homogeneous solver also stored `de = eg[1:] - eg[:-1]` as a
signed value, so `flux_per_energy` had the same sign flip.

**Impact:** None on eigenvalues — `flux_per_lethargy` is used only
for spectrum visualization, never fed back into a solver. The
visual output of every MC / MOC / homogeneous spectrum plot was
mirror-flipped through `y = 0`, which readers were silently
compensating for by reading the magnitudes and ignoring the sign.

**Fix:** Take the absolute value at the *definition* site of
`du` / `de`, not at the consumer site, so "lethargy bin width" and
"energy bin width" are non-negative by construction regardless of
grid ordering:

    du = np.abs(np.log(eg[1:] / eg[:-1]))
    de = np.abs(eg[1:] - eg[:-1])

Applied in `orpheus.mc.solver` (MCResult assignment),
`orpheus.homogeneous.solver._spectrum_result`, and
`orpheus.plotting.plot_moc_spectra`.

**L0 test that catches it:**
`tests/mc/test_gaps.py::test_flux_per_lethargy_nonnegative` — runs a
small MC with a descending two-group grid, asserts
`result.flux_per_lethargy >= 0` element-wise.

**Lesson:** "Width" quantities should be non-negative by convention,
and sign should be a property of a *direction*, not of a *measure*.
When the same quantity is computed at three different call sites, fix
it at the *definition* site — the code equivalent of normalizing at
the source rather than patching every consumer. Every consumer is an
opportunity for the bug to resurface.

This also reinforces ERR-006 / Meta-Lesson 6: convention-dependent
values (like "is `eg` ascending or descending?") must be pinned down
at a single source of truth, and every helper must be robust to the
convention, not assume it.

**Scope note:** This does *not* address the separate design question
in issue #25 about whether the MC tally should be a scattering
estimator (`w/Σ_s` on scatters, current behavior) or a collision
estimator (`w/Σ_t` on all collisions). That is a choice, not a bug —
both are unbiased estimators of the scalar flux. The sign-flip was
the genuine bug.

---

## ERR-023 — MC solver silently ignores Sig2 (n,2n) reactions

**Failure mode:** #6 Untested code path — every existing MC test
material had `Sig2 = 0`, so the missing branch was invisible.
**Date:** 2026-04-15
**Solver:** MC (`orpheus.mc.solver._random_walk`,
`orpheus.mc.solver._precompute_xs`).

**Bug:** The random walk only computed `sig_t = sig_a + sig_s_sum`
and used a two-way branch between absorption and scatter. The `mat.Sig2`
matrix was never touched — no reaction was sampled and no weight
doubling occurred. At the same time `_precompute_xs` seeded the
majorant with `mix.SigT`, which by the project's convention
(`orpheus.data.macro_xs.mixture._compute_mixture`, line 142) already
includes one copy of `Sig2.sum(axis=1)`. The mismatch meant the Σ_2n
fraction of the majorant was effectively *always* rejected as a
virtual collision: the particle free-flighted past (n,2n) sites
without ever sampling them. Net effect: zero (n,2n) contribution to
the scattering kernel, whereas the CP solver correctly includes
`2·Sig2·φ` as a source (anti-ERR-015).

**Impact:** Bias on `keff` whenever a material has nonzero (n,2n).
For the 2 G Region-A fixture with `Sig2[0,0] = 0.01`, the
analytical `k_inf` is 0.817 vs 0.800 with Sig2 = 0 — a 2 % shift
that the MC was unable to reproduce.

**Fix:** In `_random_walk`, compute `sig_2n_row = Sig2[ig, :]` and
`sig_2n_sum = sig_2n_row.sum()`, set
`sig_t = sig_a + sig_s_sum + sig_2n_sum`, and add a third branch in
the collision decision:

    r = rng.random() * sig_t
    if r < sig_s_sum:            ... # scatter
    elif r < sig_s_sum+sig_2n:   w *= 2.0; sample exit from Sig2 row
    else:                        ... # absorb / fission

The majorant stays as `mix.SigT` (unchanged) — because mixture.SigT
already carries `Σ_2n.sum` once, no `2·` factor is needed at the
majorant level. The weight doubling inside the (n,2n) branch is the
analog-MC convention for "one reaction, two neutrons emitted."

**L1 test that catches it:**
`tests/mc/test_gaps.py::test_mc_n2n_keff_matches_analytical` — builds
the 2 G Region-A mixture with `Sig2[0,0] = 0.01`, solves a scipy
generalised eigenvalue problem with effective loss
`SigT − Σ_s^T − 2·Σ_2n^T`, and checks that the MC keff matches to
`5σ + 5·10⁻³`. Also checks that the MC has moved at least halfway
from the Sig2 = 0 baseline toward the (n,2n) reference.

**Lesson:** Reinforces Meta-Lesson 6 (zero cross sections hide bugs).
A structurally correct-looking `sig_t = sig_a + sig_s_sum` is only
correct if *every* term in the project's total-XS definition is
accounted for. When a mixture field (here `Sig2`) is **never read**
inside a transport kernel, that field is silently dropped on the
floor — and every downstream test that happens to use a zero value
for it gives false confidence.

---

## ERR-024 — MC flux tally: scattering estimator instead of collision estimator

**Failure mode:** #3 Design error rather than term error — the tally
was a well-defined *scattering* estimator but the output field
`flux_per_lethargy` claimed to represent the scalar flux.
**Date:** 2026-04-15
**Solver:** MC (`orpheus.mc.solver._random_walk` tally accumulation).

**Bug:** On each real scattering event the solver accumulated
`tally[ig] += w / sig_s_sum`. This is the textbook *scattering*
estimator for a response-like integral weighted by Σ_s, but the
`MCResult.flux_per_lethargy` field was divided only by `|du|` and
treated as a scalar flux by plotting and by
`tests/mc/test_gaps.py::test_2g_flux_ratio_homogeneous`. Absorption
events contributed nothing, and the per-event weighting was
`1/Σ_s` instead of the `1/Σ_t` required by a collision estimator.
The existing spectral test had to be loosened to `ratio > 0.1` to
accommodate the bias (issue #25).

**Impact:** None on keff (the eigenvalue is computed from the
weight ratio in :eq:`keff-cycle` and never touches the tally). Bias
on every flux-spectrum plot, proportional to the relative shape
difference between Σ_s and Σ_t across groups. For Region A the
shape distortion is visible at the ~10 % level.

**Fix:** Move the tally inside the "real collision" branch
*before* the scatter-vs-(n,2n)-vs-absorb decision and use

    tally[ig] += w / sig_t

where `sig_t` is the real total (not the majorant). This is the
standard collision estimator and is unbiased for any combination of
reactions being sampled.

**L1 test that catches it:**
`tests/mc/test_gaps.py::test_2g_flux_ratio_homogeneous` — now
compares the MC flux shape (per-group, normalised) against the
analytical eigenvector from
`scipy.linalg.eig(F, diag(SigT) - Σ_s^T - 2·Σ_2n^T)` with
`rtol = 10 %`. The previous scattering-estimator bias cannot satisfy
this tolerance because the `Σ_s / Σ_t` ratio differs across groups.

**Lesson:** Unbiasedness alone is not a specification. A scattering
estimator and a collision estimator are *both* unbiased for their
respective integrals, but only one of them estimates the *flux*. If
a result field is labelled `flux_*`, the code that produces it must
be a flux estimator — and that contract should be enforced by a
spectral test, not a marginal `> 0.1` placeholder.

---

## ERR-025 — Diamond-difference cumprod recurrence: missing −Σ_t in numerator and missing 1/W source normalization

**Failure mode:** #3 Missing factor + #4 Factor error (two
compensating factor-of-two errors that cancel for homogeneous
problems)
**Date:** 2026-04-16
**Solver:** SN (`orpheus.sn.sweep._sweep_1d_cumprod`, 1D Cartesian
Gauss-Legendre fast path).

**Bug:** The precomputed face-flux recurrence coefficients were

    a = 2μ / (2μ + Δx·Σ_t)           # WRONG — missing −Σ_t in numerator
    b = 0.5·Δx·Q / (2μ + Δx·Σ_t)     # WRONG — missing 1/W, extra factor 0.5

instead of the canonical diamond-difference (DD) recurrence derived
symbolically in `orpheus.derivations.sn_balance.derive_cumprod_recurrence`:

    a = (2μ − Δx·Σ_t) / (2μ + Δx·Σ_t)
    b = 2·Δx·(Q/W) / (2μ + Δx·Σ_t)

where `W = Σ w_n` is the quadrature weight sum. The `1/W` factor is
needed because `SNSolver._add_scattering_source` produces `Q` in
scalar-flux units while the per-ordinate transport equation sees
`Q/W` on the right-hand side — the same normalization
`_sweep_2d_wavefront` already applied via its `weight_norm = 1/W`
factor (`Q_scaled = Q * weight_norm`). The 1D fast path had been
independently derived without that normalization, and its `a`
formula had an additional sign error in the numerator.

**Why the two errors cancel for homogeneous problems:** the fixed
point of the buggy recurrence is `ψ = Q/(2Σ_t)`, half the correct
`ψ = Q/Σ_t`. The missing `1/W = 1/2` for Gauss-Legendre on `[−1, 1]`
rescales by exactly 2, turning `Q/(2Σ_t)` back into `Q/Σ_t` per
ordinate. The resulting scalar flux is correct up to a uniform
rescaling by `Σ_t(x)`. For eigenvalue problems with a single
material this is invisible, because the Rayleigh quotient
`k = νΣ_f·φ / Σ_a·φ` is invariant under a uniform rescaling of φ.
At a material interface the rescale factor depends on which side
of the interface you are on, so the cancellation breaks and k_eff
shifts.

**Impact:** ~1.48 × 10⁻² error in k_eff on the ORPHEUS Phase 2.1b
2-region A+B reflective slab (fuel Σ_t=1, Σ_s=0.5, νΣ_f=0.75; mod
Σ_t=2, Σ_s=1.9, νΣ_f=0; reflective BCs). Case
singular-eigenfunction reference and CP slab E₃ kernel both give
k ≈ 1.27461, while the buggy solver converged to ≈ 1.25988.

**How it hid from higher-level tests:**
- Homogeneous single-region k_inf tests: exact to machine precision
  (uniform rescaling of φ is eigenvalue-invariant).
- Same-material two-region tests: exact for the same reason.
- Smooth-Σ MMS verification (Phase 2.1a): passed cleanly because
  the MMS consumer test uses `solve_sn_fixed_source`, which goes
  through the `_sweep_2d_wavefront` path with the correct
  `weight_norm = 1/W`.
- Self-referencing Richardson convergence tests on heterogeneous
  problems: saw clean O(h) convergence **to the wrong asymptote**,
  because the Richardson reference was built from the same buggy
  solver. This is exactly the T3 dead-end pattern documented in
  `docs/theory/diffusion_1d.rst` "Investigation history".

**Fix:** `orpheus/sn/sweep.py:119-140` — replaced the wrong
coefficients with the canonical DD recurrence, added a
source-of-truth comment pointing at `derive_cumprod_recurrence`.
One-formula correction; nothing downstream of the coefficients
needed changes.

**Evidence after fix:**

| Method                           | k_eff       |
|----------------------------------|-------------|
| Case singular-eigenfunction (S8) | 1.27461604  |
| CP slab E₃ kernel (converged)    | 1.27442847  |
| solve_sn S8 @ n_per=320 post-fix | 1.27461601  |

Case ↔ solve_sn agreement at matching quadrature order improved
from 1.48 × 10⁻² to 3.4 × 10⁻⁸.

**L1 test that catches it:**
`tests/sn/test_cartesian.py::test_heterogeneous_absolute_keff` — pins
the 2-region A+B reflective slab against the Case singular-eigenfunction
reference to 5 × 10⁻⁴. Without a material interface the bug is
invisible (the Rayleigh quotient's rescale invariance hides it), so
this test is the minimal configuration that exposes it.

---

## ERR-026 — Curvilinear sweep WDD angular closure converges to wrong fixed-source solution

**Failure mode:** #6 Wrong answer accepted — sweep converges to a
stable, balance-satisfying solution that is NOT the correct discrete
transport solution
**Date:** 2026-04-17
**Solver:** SN (`orpheus.sn.sweep._sweep_1d_spherical`,
`_sweep_1d_cylindrical`).

**Bug:** The curvilinear sweeps use a one-directional WDD angular
face-flux closure:

    ψ_{n+1/2} = (ψ_n − (1−τ)·ψ_{n−1/2}) / τ

while the BiCGSTAB transport operator
(`build_transport_linear_operator_spherical`) uses a symmetric closure:

    ψ_{n+1/2} = τ·ψ_{n+1} + (1−τ)·ψ_n

Both are consistent for flat flux analytically. But the sweep's
one-directional WDD, combined with the zero-area face at r=0
(which eliminates spatial coupling at the innermost cell), creates
a system where the iterative sweep converges to a NON-FLAT solution
that still satisfies the discrete balance equation.

**Impact:** `solve_sn_fixed_source` on spherical / cylindrical meshes
produces 35–50% error at cell 0, **growing** with mesh refinement
(divergent, not convergent). MMS verification (Phase 3.3–3.4) is
blocked. Global conservation is exact despite the wrong spatial profile.

**How it hid from higher-level tests:**
- Eigenvalue solver routes to BiCGSTAB for curvilinear geometry — the
  sweep is never exercised by the eigenvalue path
- 1-group k_eff is shape-independent (Rayleigh quotient invariance)
- Multi-group eigenvalue tests use 2% tolerance that absorbs the sweep
  error
- No fixed-source tests existed for curvilinear geometry until the
  Phase 3.3 MMS attempt

**Evidence:**

| Test (constant source, reflective BCs) | Sweep    | BiCGSTAB |
|-----------------------------------------|----------|----------|
| φ at cell 0 (nx=20)                    | 0.64     | 1.000    |
| Error vs refinement                    | Diverges | Zero     |
| Conservation (volume-weighted average)  | Exact    | Exact    |
| Cartesian sweep (same test)             | Exact    | N/A      |

**Fix:** Route `solve_sn_fixed_source` through BiCGSTAB for curvilinear
geometry, adding an external-source slot to `build_rhs_spherical` /
`build_rhs_cylindrical`. GitHub Issue #98 tracks the fix, Issue #99
tracks the blocked MMS verification.

**L1 test that catches it:**
`tests/sn/test_sweep_operator_inconsistency.py::test_spherical_sweep_vs_bicgstab_flat_flux`
— runs constant-source reflective-BC problem via both sweep and BiCGSTAB
and asserts BiCGSTAB is exact while documenting the sweep's deviation.

A cheaper L0 alternative — a direct unit test of the fixed-point
of the 1D cumprod recurrence on a 1-cell uniform-material slab
that would catch the factor-of-2 in milliseconds without needing
any reference solver — is tracked for a future commit as part of
the broader derivation-implementation audit (issue #95).

**Lesson:** When a symbolic derivation module exists
(`orpheus.derivations.sn_balance.derive_cumprod_recurrence`), its
output is the source of truth and the implementation must visibly
match. A comment in the implementation pointing back at the
derivation function would have caught this at review time. Two
opposite-sign factor-of-two errors cancelled exactly for eigenvalue
problems because the only thing that matters to
`k = νΣ_f·φ / Σ_a·φ` is the *shape* of φ, not its scale — and
factor-of-Σ_t(x) cancellations hide everywhere except at material
interfaces. Phase 1.2 learned the same lesson for diffusion
(hardcoded tolerances masking quadratic convergence); Phase 2.1b
is the same pattern repeated for SN. See GitHub issue #95 for the
follow-up audit work checking every solver implementation against
its derivation.

**Status:** **PARTIAL CLOSURE** by Wave E Round 3 (SN reshape
campaign; GitHub Issues #98, #99, #164).

What Round 3 closed (BC infrastructure layer):

The vacuum-BC equation-map gap that blocked Round 2 was closed by
extending :func:`solution_to_angular_flux*` and
:func:`transport_operator_matvec*` to consume the
:class:`~orpheus.geometry.boundary.BoundaryOperator` instances on the
:class:`~orpheus.sn.geometry.SNMesh` (Wave B Issue 7 BC algebra).
Each function now dispatches through ``bc.apply_to_incoming(out, quad)``
so vacuum, reflective, white, periodic, albedo, and mixed BCs are all
honoured uniformly. Bit-identity to the pre-Round-3 reflective-only
fill is preserved for :class:`SpecularBoundaryOperator` (the standard
``BC.reflective`` factory), which is the load-bearing condition for
the 11 frozen regression snapshots to stay green.

This BC plumbing is the load-bearing infrastructure for any future
ERR-026 closure on fixed-source MMS: the
:class:`SNStreamingOperator` now solves the *correct* operator
equation for any BC family, so a Krylov-on-``apply`` solve can give
the right answer if the operator's discretization is otherwise
2nd-order accurate. The constant-source flat-flux test in
:file:`tests/sn/test_sweep_operator_inconsistency.py` confirms this:
``inner_solver="krylov"`` produces the analytical flat flux to
round-off where the sweep produces the documented ERR-026 deviation.

What Round 3 left open (FD operator boundary truncation):

Empirically the symmetric-closure FD operator at the curvilinear
outer face uses cell-center as a face-flux approximation
(``psi_right = fi[:, n, i, 0]`` at ``i = nx-1`` for outgoing μ > 0).
This is exact for *constant* solutions but only first-order accurate
on non-constant solutions like the manufactured ``A(r) = sin(πr/R)``
ansatz used by the curvilinear MMS test suite. Switching the
``solve_sn_fixed_source`` curvilinear default from
``"source_iteration"`` to ``"krylov"`` therefore *regresses* the
MMS convergence rate from the WDD sweep's ~:math:`\mathcal{O}(h^{1.3})`
(ERR-026-affected, but a benign volumetric-error mode for these MMS)
to ~:math:`\mathcal{O}(h^{1})` (FD operator's boundary truncation).
Round 3 therefore *keeps* ``inner_solver="source_iteration"`` as the
default for all geometries; ``"krylov"`` is opt-in and correct for
constant-source problems but not the right default for MMS.

Full ERR-026 closure on MMS depends on a follow-up that extrapolates
the curvilinear outer-face flux at second order (DD diamond relation
at the boundary, or an analogous ghost-cell technique).

The two MMS xfail-strict tripwires
(:file:`tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py`,
both anisotropic spherical and cylindrical cases) therefore stay
xfail through Round 3 with updated reason strings reflecting the
partial closure.

What Wave H Phase A added (commit ``d73ef68``, GH #168 Phase A):

The cell-center / BC-face storage conflation in
:func:`solution_to_angular_flux*` (Defect 2 of GH #168) was closed by
a structural rewrite — :func:`solution_to_angular_flux*` now returns
``(fi, boundary_face_flux)`` with separate storage for cell-centre
values and BC face values, and the matvec uses the new
:class:`BoundaryFaceFlux` Protocol
(:mod:`orpheus.sn.spatial.boundary_face_flux`) for the outgoing outer
face flux. The default strategy
:class:`DDExtrapolation` uses
``psi_face_out = 1.5·psi[N-1] − 0.5·psi[N-2]`` (one-sided
second-order; closes Defect 1 of GH #168).
:class:`CellCenter` reproduces the legacy first-order substitution
for ablation tests. Phase A regenerated the regression snapshot
contract: the 5 Cartesian snapshots stay bit-identical, the 6
curvilinear snapshots were intentionally invalidated and now skip
gracefully via ``if not snapshot_file.exists(): pytest.skip(...)``
pending the Wave H Phase C regeneration.

What Wave H Phase B added (Architecture only, GH #168 Phase B):

The angular-redistribution closure was lifted from inlined matvec
math to a new :class:`PoleAngularClosure` Protocol
(:mod:`orpheus.sn.spatial.pole_angular_closure`) with three concrete
strategies:

- :class:`LegacyTauSymmetricInterpolation` (default) — bit-for-bit
  reproduction of the pre-Phase-B inlined τ-symmetric form on
  arbitrary input. Preserves the regression contract under Phase B.
- :class:`BaileyFlatFluxRedist` — algebraic flat-flux collapse
  ``redist = -μ_n·ΔA_i·ψ_n,i / V_i``. Equivalent to the legacy form
  on flat ψ; differs per-ordinate on angularly-varying ψ
  (the Defect 3 disagreement).
- :class:`MorelMontryAngularSweep` (opt-in) — canonical Hébert §3.9.4
  per-cell M-M weighted DD recurrence (Eqs. 3.428, 3.432–3.439). The
  per-ordinate redistribution is
  ``(α_{n+1/2}·ψ_{n+1/2} − α_{n-1/2}·ψ_{n-1/2}) · ΔS_i / (2·𝒲_n·V_i)``
  with starting condition ``ψ_{1/2}=0`` and DD closure
  ``ψ_{n+1/2} = 2·ψ_n − ψ_{n-1/2}``.

The α-recursion normalisation between ORPHEUS and Hébert was pinned
explicitly: ORPHEUS uses ``α^O[n+1] = α^O[n] − 𝒲_n·μ_n`` with
divisor ``ΔA_i / 𝒲_n``; Hébert uses ``α^H_{n+1/2} = α^H_{n-1/2} −
2·𝒲_n·μ_n`` with divisor ``ΔS_i / (2·𝒲_n)``. The factor of 2 is
absorbed into the recurrence normalisation: ``α^O = α^H / 2``. Both
forms are mathematically equivalent. Documented in
:mod:`orpheus.geometry.reduced_operator`,
:mod:`orpheus.sn.spatial.pole_angular_closure`, and
:doc:`docs/theory/discrete_ordinates`.

Phase B did NOT close ERR-026. The empirical finding pinned in
:func:`tests.sn.test_snstreamingoperator.test_apply_spherical_constant_flux_under_morel_montry_canonical_form`
is that pairing :class:`MorelMontryAngularSweep` with the apply
matvec's existing **spatial** closure (interior arithmetic average
``0.5·(ψ_i + ψ_{i+1})`` + outer DD extrapolation) gives a *worse*
operator than the legacy form on flat ψ:
``test_spherical_sweep_vs_bicgstab_flat_flux`` regresses (φ ranges
0.6–1.004 instead of analytical 1.0); MMS probes diverge to
~10^14 at nc=10. The DD angular recurrence on flat ψ produces
oscillating half-angle face fluxes ``0, 2c, 0, 2c, …``, which
combined with the symmetric spatial average gives garbage.
Mechanically: the apply matvec must use the same WDD spatial closure
the sweep uses (``ψ_face_out = 2·ψ_avg − ψ_face_in``) for the
canonical angular closure to be consistent. That's Phase C scope.

Phase B citation correction:
:mod:`orpheus.geometry.reduced_operator` previously cited Bailey,
Adams, Yang & Zika (2009) JCP 227 — a piecewise-linear FE diffusion
paper unrelated to curvilinear S\ :sub:`N`. Corrected to **Bailey,
Morel & Chang (2010), NSE 165(2):149-169** (LLNL-JRNL-420356, OA at
https://www.osti.gov/servlets/purl/1020346). Hébert (2009) §3.9.4
is the primary source; Bailey-Morel-Chang 2010 is the auxiliary
asymptotic-diffusion-limit τ-clamp justification. Citations updated
in :mod:`orpheus.geometry.reduced_operator`,
:mod:`orpheus.sn.sweep`, :mod:`orpheus.sn.spatial.diamond`.

Phase C scope (full ERR-026 closure):

1. **Spatial-closure alignment**: rewrite the apply matvec's
   interior face-flux closure from arithmetic average
   ``0.5·(ψ_i + ψ_{i+1})`` to the sweep's WDD form
   ``ψ_face_out = 2·ψ_avg − ψ_face_in`` (design memo §6.4 / §11).
   The Phase A :class:`BoundaryFaceFlux` Protocol stays; the
   interior face closure receives a parallel reformulation.
2. **Default flips**: once spatial closures are aligned, flip
   :attr:`SNMesh.pole_angular_closure` default
   :class:`LegacyTauSymmetricInterpolation` →
   :class:`MorelMontryAngularSweep`; flip curvilinear
   ``solve_sn_fixed_source`` default ``"source_iteration"`` →
   ``"krylov"``.
3. **Snapshot regeneration**: regenerate the 6 deleted curvilinear
   regression snapshots with the Phase-C-corrected operator AND
   verify each via FN-method cross-check on the closest Sood
   La13511 case (per :mod:`orpheus.derivations.continuous.sood_registry`).
4. **Marker removal**: the four xfail-strict ERR-026 tripwires
   (the two anisotropic + two isotropic curvilinear MMS tests)
   come off. ERR-026 status: PARTIAL CLOSURE → CLOSED.

Tests (current state through Phase B):

- `tests/sn/test_sweep_operator_inconsistency.py` — pinned as the
  ERR-026 evidence ledger; the sweep still produces the documented
  WDD deviation when explicitly invoked, the krylov path gives the
  correct answer for constant-source problems under the Phase B
  default :class:`LegacyTauSymmetricInterpolation`.
- `tests/sn/spatial/test_pole_angular_closure.py` (NEW) — 28
  foundation tests pinning Protocol contract, α-recursion identity
  (Hébert Eqs. 3.423-3.424), 2-ordinate hand-calcs, and the
  BFF↔MMS Defect-3 disagreement on angularly-varying ψ.
- `tests/sn/l1_analytical/test_pole_closure_flat_flux_identity.py` (NEW) —
  5 L1 tests pinning flat-flux invariants (Legacy↔BFF bit-for-bit
  on flat ψ; MMS preserves angular-integrated invariant by
  α-telescoping).
- `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py` —
  xfail-strict, awaits Phase C spatial-closure alignment.
- `tests/sn/test_mms_curvilinear.py` (legacy isotropic ansatz) —
  fails with order ≈ 0 on the WDD sweep; awaits Phase C.

---

## ERR-027 — Peierls slab K-matrix: naive GL collocation for cross-panel entries

**Failure mode:** #3 Missing factor — missing quadrature resolution
(one-point rule where adaptive is required)
**Date:** 2026-04-19
**Solver:** CP / Peierls (`orpheus.derivations.peierls_slab._build_kernel_matrix`).

**Bug:** For cross-panel entries ``K[i, j]`` (observer node *i* and
source node *j* in different panels) the assembly used the one-point
collocation rule

    K[i, j] = (1/2) * E_1(σ_t * |x_i - x_j|) * w_j

i.e. the integral ∫_panel E_1(σ_t|x_i - x'|) L_j(x') dx' was
approximated by evaluating the kernel at the single node x_j and
multiplying by the GL weight w_j. This is only exact when the
integrand R(σ_t|x_i - x'|) × L_j(x') is polynomial of degree ≤ 2p−1.
E_1 is transcendental with a near-log spike at x' → x_i; the one-point
rule leaves ~1% quadrature error even for panel pairs with modest
optical separation, worst when x_i is within a small optical distance
of the source panel boundary (where the near-log spike sits just
outside the panel).

**Impact:** K[i, j] cross-panel entries wrong by 0.5–1.5% at 2 panels
× p=4; refining panels reduces this only at O(h) rate (non-spectral).
Downstream: slab k-eff tests showed ~0.4% tie-point offsets in the
Sanchez R≈1.98 case.

**How it hid from higher-level tests:**
- Row-sum K·[1, 1, …] is exact to O(1e-15) even with buggy cross-panel
  entries because ∑_j L_j(x') = 1 (partition of unity) kills the kink
  in the summed integrand. Every existing row-sum test was blind.
- k-eff test tolerance was 2% — absorbed the ~0.4% propagated error.
- 1-group and uniform-source tests also cancel the bug via partition
  of unity.

**L1 test that catches it:**
`tests/derivations/test_peierls_reference.py::TestSlabKMatrixElementwiseVsReference::test_cross_panel_boundary_neighbour_elementwise`
— element-wise `K[4, 3]` at n_panels=2, p=4 vs the adaptive
`slab_K_vol_element` reference at 1e-10.

**Fix:** Unified basis-aware Nyström assembly — every ``K[i, j]`` is
``(1/2) ∫_panel E_1(τ(x_i, x')) L_j(x') dx'`` evaluated via adaptive
``mpmath.quad``. Mirrors the adaptive reference
`peierls_reference.slab_K_vol_element`. See issue #113.

**Lesson:** "Row-sum conservation" tests are systematically blind to
basis-individual quadrature errors that happen to sum to zero under
partition-of-unity of the Lagrange basis. Element-wise K[i, j]
verification against an adaptive reference is the only reliable
L0/L1 gate for Nyström kernel assembly.

---

## ERR-028 — Peierls slab K-matrix: GL collocation of remainder R(τ) has unresolved kink at x'=x_i

**Failure mode:** #3 Missing factor — missing subdivision hint
**Date:** 2026-04-19
**Solver:** CP / Peierls (`orpheus.derivations.peierls_slab._build_kernel_matrix`).

**Bug:** The same-panel singularity-subtraction branch split the kernel
into ``E_1(τ) = R(τ) − ln τ − γ`` with ``R(τ) = E_1(τ) + ln τ + γ``
the smooth remainder, then used GL collocation for the ``R`` integral
and exact product-integration weights for the ``−ln`` integral. The
flaw: ``R(σ_t|x_i - x'|)`` is smooth in τ but its argument ``|x_i - x'|``
has a C⁰ kink in x' at x'=x_i. GL cannot integrate a function with a
derivative discontinuity in the interior of the interval — produces
~1% error on the diagonal panel.

**Impact:** Diagonal-panel entries (~40% of all K entries for p=4,
n_panels=8) wrong by ~1% relative. The row-sum identity survives by
partition of unity (see ERR-027), but element-wise K[i, i] deviates
by ~1.2e-2 at n=2, p=4.

**How it hid from higher-level tests:** Same cause as ERR-027 —
partition-of-unity of the Lagrange basis cancels the kink in the
row-sum, hiding from every conservation-style test.

**L1 test that catches it:**
`tests/derivations/test_peierls_reference.py::TestSlabKMatrixElementwiseVsReference::test_small_case_elementwise_agreement`
— element-wise K[i, j] including diagonal entries vs the adaptive
`slab_K_vol_element` reference at 1e-10.

**Fix:** Unified with ERR-027: adaptive `mpmath.quad` with the
subdivision hint ``[panel_a, x_i, panel_b]`` for same-panel entries.
mpmath handles both the log singularity and the derivative kink
natively. See issue #113.

**Lesson:** Singularity subtraction splits a kernel into "smooth" and
"singular" parts — but "smooth in τ" is not the same as "smooth in x'".
A change of variables that folds the singularity into a kink instead
of removing it merely trades one unresolved feature for another. The
adaptive-with-hint approach is more robust and unifies all cases.

---

## ERR-029 — Peierls curvilinear K-matrix: ρ/ω integration does not subdivide at ray-panel crossings or tangent angles

**Failure mode:** #3 Missing factor — missing subdivision hint
**Date:** 2026-04-19
**Solver:** CP / Peierls curvilinear
(`orpheus.derivations.peierls_geometry.build_volume_kernel`).

**Bug:** The unified polar-form Nyström assembly for sphere and cylinder
used fixed-order Gauss–Legendre on both the angular (ω) and radial (ρ)
integration domains with no subdivision. Two independent sources of
non-smoothness went unresolved:

1. **ρ-integration:** along a ray, the source position
   ``r'(ρ) = sqrt(r_i² + 2 r_i ρ cos ω + ρ²)`` crosses spatial panel
   boundaries at specific ρ values; the Lagrange basis L_j(r') has a
   derivative discontinuity at each crossing. Fixed GL cannot
   integrate across these kinks — produces 0.5–5% per-entry error
   that decays only as O(h) in n_ρ under refinement, and
   non-monotonically (GL nodes straddle kinks differently at each
   order).

2. **ω-integration:** for observers strictly outside an interior
   panel boundary ``r_b`` (``r_i > r_b``), there is a critical angle
   ``ω_c = arcsin(r_b/r_i)`` at which the ray is *tangent* to the
   boundary. Just below ω_c the ray crosses r_b twice; just above,
   zero crossings. At ω_c the ρ-crossing breakpoint structure bifurcates
   and the integrand L_j(r'(ρ)) acquires a C¹-discontinuous quadratic
   kink at the tangent ρ. Without subdivision at ω_c the outer GL on
   ω plateaus at ~1e-3 even with n_ang=128.

**Impact:** K[i, j] entries for observers near interior panel boundaries
wrong by 0.5–5%; downstream k_eff affected at the ~0.4% level
(Sanchez cylinder R≈1.98 tie-point). The rank-N investigation's
"K·[1] = 1 error 3%→7% for sphere R=10 at N=2" is consistent with
this bug: the row-sum K·[1] is exact even with buggy entries due to
``∑_j L_j(r') = 1`` (partition of unity) cancelling the kinks, but
non-uniform rank-N closures break this cancellation.

**How it hid from higher-level tests:**
- Row-sum K·[1] conservation passes at 1e-15 because partition of
  unity cancels the L_j kinks in the summed integrand. Every existing
  row-sum test was blind (``test_peierls_sphere_white_bc``, etc.).
- 1-group k_eff ∝ Rayleigh quotient of K — the shape independence
  (ERR-002 lesson) masks ~1% K errors.
- Multi-group eigenvalue tolerances were 2% — absorbed the error.

**L1 test that catches it:**
`tests/derivations/test_peierls_reference.py::TestSphereKMatrixElementwise`
— element-wise K[i, j] vs an independent shell-average reference at
n_ang=64, n_rho=32, p=3 with n_panels=2; all diagonal-entry rel diffs
gate at 1e-5. Companion test asserts monotone non-oscillating
convergence under (n_ang, n_rho) refinement.

**Fix:** Two new geometry methods on ``CurvilinearGeometry`` —
``rho_crossings_for_ray`` (panel-boundary crossings along a ray) and
``omega_tangent_angles`` (tangent critical angles per observer).
``build_volume_kernel`` now subdivides the ρ integration at each
crossing *and* the ω integration at each tangent angle, applying
independent GL rules on each sub-interval. See issue #114.

**Lesson:** Row-sum conservation tests are systematically blind to
basis-individual quadrature errors that cancel under partition of
unity. Element-wise K[i, j] verification with at least two
*mathematically-independent* references (shell-average + polar-form)
is required for curvilinear Peierls assembly. Multi-dimensional
Nyström also needs *each* integration dimension audited for kinks —
ρ-crossings were the obvious candidate, but the ω-tangent bifurcation
is equally load-bearing and easy to miss.

---

## ERR-030 — Peierls rank-N white-BC: mode-0/mode-n≥1 normalization mismatch

**Failure mode:** #1 Wrong formula — inconsistent normalization
between two integrand factors that look like the same partial-current
moment but are not.
**Date:** 2026-04-25
**Solver:** CP / Peierls curvilinear rank-N closure
(`orpheus.derivations.peierls_geometry.build_closure_operator`,
peierls_geometry.py:3618-3642).

**Bug:** The rank-N Marshak white-BC closure routes mode 0 through
the legacy ``compute_P_esc / compute_G_bc`` (no surface Jacobian,
the isotropic-escape-probability form) while modes ``n ≥ 1`` use
``compute_P_esc_mode / compute_G_bc_mode`` (with the canonical
``(ρ_max/R)²`` surface-to-observer Jacobian). The two forms live in
**different normalization spaces** for the partial-current expansion
moment ``J_n = ∫ μ·P̃_n(μ) ψ⁺ dμ`` and therefore cannot be summed
into a single rank-N ``K_bc = Σ_n u_n ⊗ v_n`` consistently.

The mismatch was introduced as a deliberate compromise (documented
in the function docstring as "for bit-exact rank-1 regression") so
that ``n_bc_modes = 1`` would reproduce the existing rank-1 Mark
behavior bit-exactly. In single-region solid cells the mismatch is
hidden because the legacy mode-0 was historically calibrated to make
the rank-1 Mark closure approximately right; adding mode-1 then
adds a small *Marshak* correction that lands the 1R rank-2 within
1 % of ``k_inf``. The published 2026-04-18 1G/1R rank-N table at
``peierls_geometry.py`` lines 3934-3961 records this as a "rank-N
converges" claim — but those numbers were achieved by the
legacy/canonical hybrid that is now demonstrably wrong in MR.

**Impact:** Sphere 1G/2R fuel-A inner / moderator-B outer
(``LAYOUTS=["A","B"]``, ``radii=[0.5, 1.0]``, ``σ_t=[1, 2]``) gives:

- rank-1 Mark: ``k_eff = 0.551`` vs ``k_inf = 0.648`` → -15 %
  (acceptable Mark closure floor)
- **rank-2 Marshak: ``k_eff = 1.015`` → +57 % (sign flip)**
- rank-3..8: ``k_eff`` plateaus at 1.08 → +67 %, NOT converging

Cylinder 1G/2R rank-2 also degrades to +18 % (smaller because Issue
#112 Phase C cyl Ki_{k+2} normalization adds a separate floor that
masks the underlying mode-0 mismatch).

Replacing legacy mode-0 with the canonical ``compute_P_esc_mode(n=0)``
gives a uniformly *worse* result in 1R (rank-2 = -29 %) — the proper
fix requires re-derivation of the rank-N partial-current expansion
end-to-end so mode 0 and mode ``n ≥ 1`` live in the same basis.
Tracked in Issue #132.

**How it hid from higher-level tests:**
- All existing rank-N tests in ``tests/derivations/test_peierls_rank_n_bc.py``
  exercise 1G × 1R only. The 1R control gives rank-2 = -1.10 % on
  sphere — a deceptive "convergence" that masked the mismatch.
- The conservation row-sum test
  (``tests/derivations/test_peierls_rank_n_conservation.py``) uses
  uniform Σ_t = 1, where ``K · 1 = 1`` holds by construction. The
  identity does NOT generalize to piecewise Σ_t (it becomes an
  integrated identity, not pointwise), so the test is structurally
  blind to MR mismatches.
- The 2026-04-22 falsification of the rank-N per-face hollow Class A
  closure (research log L21) closed the rank-N path against the F.4
  reference, but explicitly never audited Class B (solid). The
  N=1 ``boundary="white_f4"`` collapse for solid cells was treated
  as the "shipped" Class B closure — except the rank-N branch was
  still callable through ``boundary="white_rank1_mark"`` with
  ``n_bc_modes ≥ 2`` and silently broken.

**L1 tests that catch it:**
- ``tests/derivations/test_peierls_rank_n_class_b_mr_mg.py::test_class_b_mr_catastrophe_sphere_1g_2r_rank2``
  (XFAIL, strict=True) — pins the +57 % sphere catastrophe
- ``::test_class_b_mr_catastrophe_cylinder_1g_2r_rank2`` (XFAIL,
  strict=True) — pins the +18 % cylinder analog
- ``::test_class_b_mr_routing_invariance_uniform_sigma`` — passes
  today, regression gate against accidental introduction of an
  ``len(radii) > 1`` routing divergence even at uniform σ_t

**Fix:** Pending. Issue #132 tracks the re-derivation work.
Candidate paths: (a) a canonical mode-0 partial-current
normalization that reduces to Mark at ``n_bc_modes = 1``, OR (b)
replacing rank-N with the Sanchez & McCormick 1982 §III.F.1
partial-current-moment basis end-to-end (Eqs. 165-169), OR (c)
the Knyazev 1993 Ki_{2+k} polynomial expansion for cylinder
(Issue #112 Phase C, may share root).

**Lesson:** A "rank-N converges" claim must be verified at MR×MG —
single-region single-group passing rates are degenerate evidence
because two structurally-different normalizations can produce the
same ``k_eff`` by historical calibration. For any partial-current-
moment closure, the rank-1 → rank-2 step must hold across
``radii=[1.0]`` AND ``radii=[0.5, 1.0]`` with non-trivial Σ_t
breakpoints. The Issue #131 anti-pattern audit ("does the multi-
region branch silently differ from the single-region branch?")
must be performed for every closure primitive, not just the volume
kernel.

---

## ERR-031 — Test calls compute_P_ss_cylinder with radii/sig_t arguments swapped

**Failure mode:** #1 Wrong formula — positional argument inversion;
masked by the function not validating its inputs at the time the test
was written.
**Date:** 2026-04-29
**Solver:** CP / Peierls cylinder P_ss
(``orpheus.derivations.peierls_geometry.compute_P_ss_cylinder``,
test ``tests/cp/test_cylinder_pss.py::test_pss_multiregion_layer_order_matters_for_grazing_chords``).

**Bug:** ``compute_P_ss_cylinder`` has signature ``(radii, sig_t, *,
n_quad, dps)``. The test was passing ``compute_P_ss_cylinder(sig_t,
radii, n_quad=64)`` — first arg the cross-section vector, second arg
the radii. Pre-Issue #134, this silently produced a numerically
wrong P_ss because:

- ``radii = np.array([0.1, 2.0])`` (the intended ``sig_t``) — happened
  to be strictly increasing, so neither the "strictly increasing" nor
  "all positive" guards triggered;
- ``sig_t = np.array([0.5, 1.0])`` (the intended ``radii``) — used as
  per-region cross sections in the τ build.

The test compared two such broken calls against each other. Because
both were broken in the same way (same argument inversion), the
**ratio** ``p_thick_inner / p_thin_inner`` was still order-of-
magnitude consistent with the test's qualitative claim ("> 5×
asymmetry under sig_t swap"), so the assertion passed. The test
was numerically meaningless — every value in the ratio was computed
on a permuted-input integrand — but the ratio test held by accident.

**Impact:** None on production code (production callers all pass
keyword args or correct positional order). Test reported PASS for
months while exercising a wrong code path. Caught at Issue #134
when the chord_quadrature recipe migration introduced an upstream
``radii must be strictly increasing`` validation that triggered on
``[2.0, 0.1]`` (the second test call's intended-sig_t now decreasing
from positional inversion).

**How it hid from higher-level tests:**
- The test imports ``compute_P_ss_cylinder`` directly and exercises
  it in isolation; no integration test covered this specific call.
- The function had no upstream validation that ``radii`` was
  strictly increasing — so non-monotone "radii" silently flowed
  through ``chord_half_lengths`` to produce nonsense chord lengths
  rather than raising.
- The qualitative-ratio assertion ``> 5×`` is loose enough that
  even a mis-permuted integrand can satisfy it.

**L0/L1 tests that catch it:**
- After Issue #134 migration: the test now raises
  ``ValueError: radii must be strictly increasing, got [2. 0.1]``
  at the ``chord_quadrature`` contract.
- After test-bug fix (this commit): the test passes with the right
  args and verifies the qualitative ratio claim properly.

**Fix:** ``tests/cp/test_cylinder_pss.py:164-170`` — swap the
positional args to match ``compute_P_ss_cylinder(radii, sig_t, ...)``.

**Lesson:** Contract validation upstream of consumers catches bugs
that loose qualitative tests miss. Recipes that validate inputs at
construction (this is the ``Quadrature1D`` contract pattern) push
errors leftward — into clear failures at the test boundary instead
of silent numerical drift in the integrand. The Issue #134 migration
caught a year-old test bug as a side effect of routing CP through
the recipe.

---

## ERR-032 — Slab white-BC analytical reference: wrong ∫E₂ antiderivative (factor-of-two algebra bug)

**Failure mode:** #1 Wrong formula — wrong exponential-integral
antiderivative used to evaluate the partial-current balance for the
slab white-BC closed form.
**Date:** 2026-04-20 introduced (commit ``2538cfe``); caught
shortly after when the first-order ``K_bc`` row-sum identity was
exercised against the new analytical reference.
**Solver:** Peierls slab, white-BC analytical reference
(``orpheus.derivations.peierls_reference.slab_uniform_source_white_bc_analytical``).

**Bug:** Commit ``2538cfe`` shipped the closed form

.. code-block::

    φ_wrong(x) = (1/(2 Σ_t)) · [2 + (2β − 1) · (E_2(τ) + E_2(τ'))]
    β = (1 − E_3(Σ_t L)) / (1 − 2 · E_3(Σ_t L))

obtained from ``J^+(L)|_vol = (1/(2 Σ_t)) · (1 − E_3(τ_L))`` — which
uses the wrong antiderivative identity ``∫E_2 = 1 − E_3``. The correct
identity is ``∫E_2 = ½ − E_3``, and the algebraic simplification
``(½ − E_3)/(1 − 2 E_3) = ½`` collapses ``J^-`` to the constant
``1/(4 Σ_t)``, giving ``φ ≡ 1/Σ_t`` exactly (the Wigner-Seitz identity
for a uniform absorber cell).

**Why two "independent" derivations agreed at 1e-39:** A fixed-point
diagnostic of the partial-current balance was used as the
"independent" cross-check for the analytical formula. Both the
analytical derivation **and** the fixed-point iteration applied the
same wrong antiderivative when reducing volume integrals of ``E_2`` to
``E_3``. They converged to the same wrong number to machine precision
because they shared a factor-of-two algebra mistake — not because
either was right. **Lesson generalises beyond this bug:** "two
independent derivations agreeing at 1e-39" is worthless evidence if
both derivations share an upstream identity.

**Impact:** No production code consumed ``φ_wrong`` (it was a fresh
analytical reference, not a solver). The bug was caught before any
downstream rank-N or MR code was wired to it.

**How it hid from higher-level tests:** The reference was new — there
were no upstream consumers to drive it through a regression. The
fixed-point cross-check that "should have" caught it was structurally
unable to catch it (same factor-of-two error in both branches). Only
exercising the row-sum identity ``Σ_j K_{ij} · 1 = ∫_0^L ½ E_1(Σ_t |
x_i − x'|) dx' = (1/(2 Σ_t)) [2 − E_2(τ_i) − E_2(τ_i')]`` — a
*different* analytical identity, derived from a *different* integrand
(the kernel itself, not its volume integral) — exposed the
disagreement: factor ``~ 2.2`` between row-sum from ``K`` and
``φ_wrong``.

**L1 test that catches it:**
- ``tests/derivations/test_peierls_reference.py::TestSlabKernelRowSum::test_row_sum_matches_analytical_uniform_source``
  (``@pytest.mark.l1``, ``@pytest.mark.verifies("peierls-unified")``,
  ``@pytest.mark.catches("ERR-032")``). Asserts the row-sum identity
  to ``< 1e-8`` over a small ``(L, Σ_t)`` grid using the *vacuum-BC*
  uniform-source analytical (``slab_uniform_source_analytical``),
  which uses the correct ``∫E_2 = ½ − E_3`` identity. Disagreement of
  ``φ_wrong`` with this row-sum was the first-order fingerprint.

**Fix:** Re-derive ``J^-`` with the correct antiderivative; the result
collapses algebraically to ``φ ≡ 1/Σ_t`` (uniform Wigner-Seitz). The
shipped function in
``orpheus.derivations.peierls_reference.slab_uniform_source_white_bc_analytical``
returns this constant. See ``docs/theory/peierls_nystrom.rst``
``§White-BC analytical flux — slab`` (the canonical derivation; the
"History of the algebra bug" subsection in that page now points at
this catalog entry).

**Lesson:** Cross-checks must be *structurally* independent, not just
*procedurally* independent. A fixed-point iteration written in
parallel with an analytical derivation is procedurally independent
(different code paths) but structurally coupled (both reduce the same
class of integrals using the same antiderivative table from the same
human's working memory). The most reliable cross-checks come from a
different integrand or a different identity — not the same identity
written in two ways. Apply at every analytical-reference shipping
point: pick the cross-check from the *kernel* (row-sum, balance) and
the *closed form* (eigenvalue, asymptotic limit) — not two
derivations of the same closed form.

---

## ERR-033 — Peierls slab single-surface aggregate: finite-N GL fall-through where ½·E₂(τ) closed form applies

**Failure mode:** #3 Missing factor — finite-N Gauss–Legendre branch
left in place where a closed-form exponential-integral identity reduces
the integrand exactly.
**Date:** 2026-04-23 (audit; Issue #131 follow-up).
**Solver:** CP / Peierls slab single-surface aggregate
(``orpheus.derivations.peierls_geometry.compute_P_esc`` /
``compute_G_bc``, peierls_geometry.py around lines 1350 and 1448).

**Bug:** Issue #131 (commit ``3b0b2c9``) replaced finite-N GL with
the closed form ``½·E_2(τ_total)`` in the per-face primitives
``compute_P_esc_outer`` / ``compute_P_esc_inner`` /
``compute_G_bc_outer`` / ``compute_G_bc_inner`` for slab-polar at
``len(radii) > 1``. The legacy single-surface aggregates
``compute_P_esc`` / ``compute_G_bc`` — which are still reachable
from ``build_white_bc_correction`` and from
``build_closure_operator(reflection="marshak")`` mode-0 on the
``boundary="white_rank1_mark"`` legacy rank-1 single-surface layout
— retained the original finite-N GL branch on the
``len(radii) > 1`` path. The slab angular integral has the form
``∫_0^1 f(µ)·exp(-τ/µ) dµ`` with τ piecewise-constant in σ_t (and
therefore µ-independent), so it IS closed-form as ``E_n``
regardless of ``n_regions``. Finite-N GL leaves a
~4 × 10⁻³ to 6 × 10⁻⁵ relative error that decays only as O(N⁻²)
in the GL node count.

**Impact:** None on the shipped 2026-04-23 ``peierls_slab_2eg_2rg``
reference (it routes through ``white_f4`` and bypasses the legacy
aggregates entirely). Any user selecting
``boundary="white_rank1_mark"`` on a multi-region slab hits the
same residual GL artefact that Issue #131 patched on the per-face
primitives.

**How it hid from higher-level tests:**
- The shipped ``peierls_slab_2eg_2rg`` reference uses ``white_f4``,
  so the entire ``boundary="white_rank1_mark"`` × multi-region path
  was untested.
- The per-face Issue #131 fix landed alongside an ``len(radii) == 1``
  guard test, but the audit did not extend to the legacy
  ``compute_P_esc`` / ``compute_G_bc`` callers — they share the
  ``slab-polar`` dispatch but not the same branch.

**L1 test that catches it:** Issue #131 multi-region slab test
(``tests/derivations/test_peierls_slab_multiregion.py``) — extended
to the ``boundary="white_rank1_mark"`` path will surface the same
4 × 10⁻³ → 6 × 10⁻⁵ bit-deficit the per-face fix already eliminates.

**Fix:** Replace the ``len(radii) > 1`` branch in
``compute_P_esc`` / ``compute_G_bc`` with the closed form
``½·(E_2(τ_inner) + E_2(τ_outer))`` (and the corresponding ``2·`` form
for ``compute_G_bc``) using the same ``_slab_tau_to_{inner,outer}_face``
primitive as the per-face fix — identical math, identical precision,
no new quadrature.

**Lesson:** When a closed-form identity has been applied to one of
several call sites that share the same integrand class, treat the
remaining call sites as suspects until each one is audited. The
``compute_P_esc`` / ``compute_G_bc`` aggregates and the
``compute_P_esc_outer`` / ``compute_P_esc_inner`` per-face forms
both dispatch on ``kind == "slab-polar"`` and both compute
``∫_0^1 f(µ)·exp(-τ/µ) dµ`` with µ-independent τ — they should have
been fixed in the same commit. Generalises Probe-Cascade rule
"closed-form detection" (Issue #131 Sphinx ``§theory-peierls-slab-polar-g5-diagnosis``):
**any slab angular integral with µ-independent τ is closed-form as
E_n regardless of n_regions, and finite-N GL on such an integrand
is a bug, not a tolerance choice.** Curvilinear angular integrals
(τ depends on µ through the chord ρ_max) remain irreducible —
this rule applies ONLY to slab-polar.

**Related-pattern audit:** The 19 finite-N GL call sites in
``peierls_geometry.py`` for curvilinear angular integrals
(lines 1378, 1476, 1504, 1528, 1675, 1763, 1838, 1918, 1997, 2021,
2096, 2120, 2188, 2265, 2516, 2582, 2751, 2839, 2867) are
correctly finite-N because τ = σ_t · ρ_max(r_i, cos Ω, R) depends
on Ω. ``build_volume_kernel`` finite-N is explicitly exempted on
the same grounds. W transmission integrals use adaptive
``mpmath.quad`` and are correct by design.

---

## ERR-034 — Slab Variant α first-leg trajectory: missing µ factor in x_traj parametrisation

**Failure mode:** #6 Convention drift — the chord arclength
parameter `s` was used to advance the position with NO direction-
cosine factor (`x_traj = x - s` instead of `x_traj = x - μ·s`),
silent for any uniform-source test.
**Date:** 2026-05-02 (Phase-3B method-of-images dispatch).
**Solver:** Slab Variant α Green's function reference, both
Phase-3A symmetric (`greens_function_slab.py::_apply_operator_slab`)
and the inherited Phase-3B asymmetric-slab module
(`greens_function_slab_asymmetric.py::_apply_operator_slab_asymmetric`).

**Bug:** The first-leg backward chord parametrisation per the
transport equation's formal solution along characteristics is
`x_back(s) = x - μ·s` (1D slab, μ = direction cosine). The Phase-3A
code wrote `x_traj = x - s_pts_first` (no μ factor), incorrectly
treating arclength `s` as x-distance. The exponential factor
`exp(-σ_t · s_pts_first)` IS arclength-correct (σ_t × arclength is
the correct optical depth weighting), so the bug is purely in the
position-lookup of the source `q(x_traj(s))` — INVISIBLE for any
spatially uniform `q`.

**Impact:** For non-uniform source profiles (which arise during
the eigenvalue iteration on a non-uniform flux), the F integral
evaluates `q` at the wrong positions along the chord, producing a
biased operator with a different fixed point than the correct
operator. Phase-3A vacuum α=0 k_eff was off by ~21% relative
(0.130 vs corrected 0.157 at fuel-A τ_L=5). The asymmetric Phase-3B
method-of-images test caught the bug because it requires the
asymmetric eigenmode to peak AT the reflective wall x=0 — the
buggy parametrisation gave a peak at x≈0.16, ~5% relative k_eff
error vs the symmetric vacuum reference, neither of which Phase-3A
exercised.

**How it hid from higher-level tests:**
- Phase-3A V_α1_slab (closed slab α=1) uses constant source `q =
  Σ_s` — `q(x_traj)` is the same value regardless of x_traj
  parametrisation. The bug is invisible.
- Phase-3A vacuum α=0 self-consistency test only checked
  `0.45 < k_eff/k_inf < 0.85` (a wide band). The 21% systematic
  offset was within the band.
- The vacuum slow-convergence floor (5e-4 between (24,40,96) and
  (32,56,128) grids) was REINTERPRETED as a quadrature limitation
  in the Phase-3A closeout memo. With the fix, vacuum converges
  faster (the slow floor was the bug's spatial bias being slowly
  resolved as the Spline interpolation got finer).

**L0 test that catches it:** Method-of-images symmetry test
(`test_method_of_images_reflective_vacuum_equals_double_vacuum`)
in `tests/derivations/test_peierls_greens_function_slab_asymmetric_solver.py`
(tagged `@pytest.mark.catches("ERR-034")` if added). The asymmetric
[0,1] reflective-vacuum eigenvalue must equal the symmetric [0,2]
vacuum-vacuum eigenvalue to ≤ 1e-7; the buggy version differed by
~5% relative.

**Fix:** Replace `x_traj = x - s_pts_first` with `x_traj = x - mu *
s_pts_first` for μ > 0, and `x_traj = x + abs_mu * s_pts_first`
for μ < 0. Applied to both Phase-3A and Phase-3B modules in the
Phase-3B integration commit.

**Lesson:** The transport equation's characteristic parametrisation
along arclength s in direction Ω̂ has Δr = s·Ω̂ (vector). In 1D, this
means Δx = s·μ — NOT Δx = s·sign(μ). The factor of μ matters for
ANY non-uniform source / non-trivial flux profile. **A test
exercising non-trivial spatial profiles (mesh-refinement on a
heterogeneous problem, or a method-of-images consistency check)
catches this immediately; uniform-source-only tests are blind by
construction.**

This is a slab-specific instance of the same "trajectory
parametrisation must respect the direction-cosine chain rule"
pattern as ERR-006 (curvilinear sphere α-recursion + ΔA/w bug).

---

## ERR-035 — Slab Variant α symmetric closure: heuristic α·B_period/(1-α²·e^{-2τ}) is wrong at intermediate α

**Status:** **FIXED 2026-05-02** by delegating Phase-3A symmetric
slab to the Phase-3B rank-2 path at :math:`\alpha_L = \alpha_R =
\alpha`. See "Fix" section below.

**Failure mode:** #1 Sign flip / convention drift — the symmetric
rank-1 closure was constructed by analogy with sphere/cylinder
(applying `alpha_per_period = α²` for the 2-bounce-per-period
geometry) without re-deriving from first principles. The
analogical formula coincides with the correct closure ONLY at
α∈{0, 1} corners.
**Date discovered:** 2026-05-02 (Phase-3B reduce-to-symmetric
consistency check).
**Date fixed:** 2026-05-02 (delegation refactor, same day).
**Solver:** Slab Variant α Green's function reference, Phase-3A
symmetric (`greens_function_slab.py::_apply_operator_slab` calling
`apply_variant_alpha_closure(... alpha_per_period=α²)` with the
full out-and-back bounce-period B integral).

**Bug:** The Phase-3A heuristic closure

    ψ_surf = α · B_period / (1 - α² · e^{-2τ})
    where B_period = ∫_0^{2L/|μ|} q · e^{-σ_t · s} ds (no α factor)

does NOT match the first-principles rank-2-at-symmetric closure

    ψ_surf = α · B / (1 - α · e^{-τ})
    where B = ∫_0^{L/|μ|} q · e^{-σ_t · s} ds (single transit)

— EXCEPT at α=1 (closed slab, where (1-e^{-τ})(1+e^{-τ}) =
1-e^{-2τ} collapses both forms to q/Σ_t for constant source) and
α=0 (vacuum, where both reduce to F). The first-principles
derivation is documented in the Phase-3B closeout memo: starting
from `ψ_L^+ = α_L · ψ(0, -μ̂)` and tracing the trajectory
explicitly through both walls yields the correct rank-2 closure,
which under symmetric BC reduces to `α · B / (1 - α · e^{-τ})`.

**Impact:** For 0 < α < 1 (partial-reflection BCs that are common
in real reactor problems), Phase-3A k_eff is off by ~1.3e-4 relative
at α=0.5 (fuel-A τ_L=5). Phase-3B's Branch-1 SymPy V_α2_slab_asym
proves the discrepancy algebraically: rank-2 at α_L=α_R=α gives
`(T_11 + T_12) = 1/(1 - α·e^{-τ})`, NOT `(1 + e^{-τ})/(1 - α²·e^{-2τ})`
that the Phase-3A formula encodes.

**How it hid from higher-level tests:**
- Phase-3A V_α1_slab uses constant source — both formulas reduce
  to q/Σ_t, agreement at machine precision masks the heuristic.
- Phase-3A vacuum-only and closed-only tests exercise α∈{0,1}
  corners only; intermediate α was never tested.
- Phase-3A's MG-2G-asymmetric test uses α=1 (closed asymmetric
  scattering matrix detector); α<1 path was never multi-group
  tested.
- The 22 Phase-3A tests are all blind to the discrepancy by
  construction — none of them stress the closure formula at
  intermediate α with a spatially non-uniform eigenmode.

**L1 test that captures it (post-fix regression gate):**
`test_rank1_path_now_agrees_with_rank2_via_delegation_after_ERR035_fix`
in `tests/derivations/test_peierls_greens_function_slab_asymmetric_solver.py`,
tagged `@pytest.mark.catches("ERR-035")`. The test was originally
`test_rank2_vs_rank1_at_intermediate_alpha_documented_discrepancy`
with a `[5e-5, 5e-4]` disagreement gate; post-fix it asserts ≤
1e-12 rtol bit-equal agreement at α=0.5, catching any
re-introduction of a heuristic Phase-3A closure.

**Fix (applied 2026-05-02):** Refactored
`solve_greens_function_slab` and `solve_greens_function_slab_mg`
into thin wrappers that delegate to
`solve_greens_function_slab_asymmetric{,_mg}` at α_L=α_R=α. The
heuristic `_apply_operator_slab`, `_first_leg_chord_slab`, and
`_bounce_period_chord_slab` helpers were deleted (dead code after
delegation — keeping them as "alternative paths" would re-create
the very pattern this fix removes). The Phase-3A public API
(`SlabGreensResult`, `SlabGreensMGResult`, the two solver
entrypoints) is preserved via re-wrap of the asymmetric result
into the original dataclass shape. Net effect: all Phase-3A code
paths now use the first-principles single-transit closure
:math:`\alpha\cdot B/(1-\alpha\,e^{-\tau})`.

**Side effects of the fix:**

- The Phase-3A vacuum α=0 convergence floor improved
  dramatically: pre-fix the floor was ~5e-4 between (24,40,96)
  and (32,56,128) grids and was misattributed to slow GL
  quadrature; post-fix the same finest pair self-consistency
  is ~8.85e-6. The improvement comes from ERR-034 (the dominant
  fix at α=0; ERR-035's closure isn't load-bearing at α=0 since
  the closure is bypassed) plus the rank-2 path's single-transit
  B integrals having less integrand variation than the full-
  period integrals. The floor gate in
  `test_alpha_zero_convergence_floor` was retightened from 1e-3
  to 5e-5 to match.
- Phase-3A V_α1_slab k_inf-exactness preserved bit-equal
  (closed slab, constant source — the closure formula difference
  cancels).
- Phase-3A 2G asymmetric scattering at α=1 → k_inf preserved at
  1e-9 rtol.
- Sphere and cylinder solvers untouched (they live in separate
  modules and continue to use the rank-1 `apply_variant_alpha_closure`
  with `alpha_per_period=alpha`).

**Lesson:** **Analogical generalisation of a closure formula
across geometries is an algebraic claim that requires a proof,
not a substitution.** Phase-3A took the sphere/cylinder rank-1
form `α · B_period / (1 - α^N · e^{-N·τ})` and applied it with
`N=2` to slab on the strength of "slab has 2 bounces per period",
WITHOUT re-deriving from first principles for the slab geometry.
The derivation actually requires careful tracking of α factors
inside the bounce-period source integral — a single α factor does
appear inside the second-leg term, NOT just as a denominator α².
**Every Variant α generalisation to a new geometry must include a
first-principles derivation of the closure on a non-uniform
source, with the result cross-checked against the rank-2 closure
at the symmetric corner.**

This is a meta-instance of ERR-032's "structural-independence
illusion": Phase-3A had the V_α1_slab algebraic identity AND the
vacuum α=0 sanity check, both of which agreed with the
implementation by construction (constant source / α=0). The bug
hid in the algebraic interpolation between those two corners that
no test exercised — exactly as ERR-032 hid in the algebraic step
that two procedurally-independent derivations both inherited from
the same identity.

---

## ERR-036 — Slab Peierls Path A.i: log-singular kernel diagonal truncation in plain GL

**Status:** **DOCUMENTED 2026-05-03**, with a structurally-distinct
fix path (Atkinson product-Nyström) implemented and tested. The
legacy plain-GL path (`slab_scalar_flux_fn_projection`) is left in
place at its honest 5–7 % tolerance — it is the *plain* baseline;
the Atkinson path
(`slab_scalar_flux_fn_projection_atkinson`) is the singularity-aware
replacement that hits the F_N moment floor.

**Failure mode:** #3 Missing factor / #4 Factor error — silent
truncation of the divergent log-singular kernel at the diagonal.
**Date discovered:** 2026-05-03 (Wave 2-A Path A.i convergence-rate
investigation).
**Solver:** `orpheus.derivations.continuous.fn_method.slab.flux_reconstruction.slab_scalar_flux_fn_projection`
(Path A.i — power-iterates the Peierls integral operator with plain
Gauss–Legendre on the (z, μ) tensor product).

**Bug:** The Path A.i discrete operator builds

    K[i, j] = (c/2) Σ_k (w_μ_k / μ_k) exp(-|z_i - z_j| / μ_k)

with all-positive `w_μ_k`. After the μ-integral, this should equal
`(c/2) E_1(|z_i - z_j|)`, where `E_1(0+) = +∞` (logarithmic
divergence). Off-diagonal (`|z_i - z_j| > 0.05` mfp) the
μ-quadrature converges to E_1 at machine precision. **At the
diagonal** (`z_i = z_j`), the integrand `(1/μ) exp(0) = 1/μ` is
itself singular at μ → 0+; with finite-N Gauss-Legendre nodes on
(0, 1) the smallest node is `μ_min ~ 1/n_μ²`, and the sum saturates
at the *finite* value `Σ_k w_k / μ_k ≈ −log(μ_min) ≈ 2 log(n_μ)`.
The discrete kernel is therefore qualitatively wrong at the
diagonal: a finite truncation of the true `+∞`. This matches
Atkinson 1972 §III's textbook scenario for log-singular kernels.

**Impact:** On the bare-critical 1G isotropic slab (Wave 2-A cases
Ua-1-0-SL c=1.30, PUa-1-0-SL c=1.50, PUb-1-0-SL c=1.40), the
discrete eigenmode `φ(z)/φ(0)` deviates from the KLL Path B
reference (Wiener–Hopf factorisation, converged to 1e-10) by:

| `z/a`   | rel err at `n_quad_z=128` |
| ------- | ------------------------- |
| 0.0     | 0 (by normalization)      |
| 0.5     | 2.0 %                     |
| 0.9     | 3.4 %                     |
| 0.99    | 4.0 %                     |

The error PEAKS at the slab edges (where the eigenmode's local
spatial structure is most sensitive to short-range kernel
contributions) and is roughly UNIFORM across `c` in 1.3–1.5.

**How it hid from higher-level tests:**

- The existing
  `test_l1_path_ai_vs_path_b_flux_ratios` was authored with an
  "honest" tolerance of 5 × 10⁻², which the legacy plain-GL Path
  A.i passes by ~ 2× margin (~ 7e-2 worst-case at n=128). The
  test design admitted the gap as expected ("singularity-aware
  quadrature out of scope"), so the discrepancy was *known but
  unrepaired* rather than missed.
- The companion `test_l1_path_ai_convergence_under_refinement`
  only checks `ratio > 1.2` for n-doubling (i.e., error decreases
  by any amount), which is consistent with the empirical -0.9
  rate — and would also be consistent with -0.5, -2, or -4
  rates. It cannot distinguish the failure mode.
- Path A.i was originally pitched as "structurally independent"
  of Path B (KLL): the procedural divergence (BTE power
  iteration in `(z, μ)` vs Wiener–Hopf factorisation in
  ν-space) is real, but procedural independence is not the
  same as structural independence (`vv-principles` §three pillars,
  §6 reference contamination). The true issue is that BOTH
  paths share the underlying integral kernel `(c/2) E_1`; what
  Path A.i mis-discretises is exactly what Path B
  reconstructs analytically. Saying they "agree to ~ 5 %" is
  cross-implementation agreement at L4 — it carries no
  *correctness* information beyond "both paths see the same
  eigenmode topology".
- The empirical convergence-rate slope is **−0.9 (close to first
  order with log correction)**, not the literature-predicted **−1/2**
  (Atkinson–Schneider C^{(0,1)} endpoint singularity). The
  `err·n/log(n)` doubling-ratio attribution shows the dominant
  rate-limiting factor is the *kernel* diagonal-truncation bug,
  not the *solution* endpoint singularity. The literature memo
  (`scratch/derivations/peierls_log_singular/atkinson_product_nystrom.md`)
  emphasised the latter; the diagnostic cascade (Phase 1.1–1.5)
  established the former is dominant.

**L1 tests that catch it:**

- `tests/derivations/test_atkinson_product_nystrom.py::test_l1_atkinson_vs_kll_5e_minus_4`
  (parametrised over the three Wave 2-A cases) — asserts
  `sup |err| ≤ 5e-4` at `n_panels=64`, ~100× tighter than the
  legacy 5e-2 threshold.
- `test_l1_atkinson_high_n_hits_fn_moment_floor` — asserts
  `sup |err| ≤ 5e-5` at `n_panels=256`, hitting the F_N moment
  reference floor.
- `test_l1_atkinson_convergence_rate_better_than_first_order` —
  asserts the n-doubling error ratio is `> 2.5` (i.e. faster than
  first order). The legacy plain-GL would give ~ 1.7–2.0.

All three are tagged `@pytest.mark.catches("ERR-036")`.

**Fix (applied 2026-05-03):** New module
`orpheus.derivations.continuous.fn_method.peierls_atkinson_nystrom`
implementing Atkinson 1972 / 1997 §4.2 product-Simpson rule. The
log-singular kernel piece is integrated *analytically* against the
piecewise-quadratic Lagrange basis on each Simpson panel via
closed-form antiderivatives `∫ s^k log|t-s| ds` for `k ∈ {0, 1, 2}`
(verified L0 against scipy.integrate.quad with explicit
singularity points). The `C^∞` smooth remainder `R(τ) = E_1(τ) +
log τ` is integrated with standard Simpson. New entrypoints:
`slab_scalar_flux_fn_projection_atkinson(...)` and its `_ratio`
companion. The legacy `slab_scalar_flux_fn_projection` is
preserved unchanged for backward compatibility with the existing
tolerance-explicit tests.

**Mechanism of the fix:** Product-Nyström replaces the discrete
operator's qualitatively-wrong diagonal entry (a truncation of
+∞) with a *correct* one — namely the integral
`∫_panel (c/2) (-log|z_i - z'|) L_l(z') dz'` evaluated in closed
form, where `L_l` is the Lagrange basis. This is exact for the
log-singular factor against any quadratic trial function. The
discrete operator therefore reproduces the *integrable* log
singularity faithfully, and the eigenmode converges at the
operator's natural Simpson rate (de Hoog–Weiss 1973: `O(h⁴ log h)`
on operator norm; on flux ratios the empirical doubling ratios
are 4×–10× per n-doubling at moderate n, until the F_N moment
floor is hit at ~ 1e-5).

**Lesson:** When a discrete kernel is built from a singular
integrand, the "natural" finite-N quadrature truncation is
silently incorrect at the singularity. The smoking-gun fingerprint
is: error converges at "first order with log correction" — `err *
n / log(n) ≈ const` — instead of either (a) machine precision
(everything's fine), (b) `1/sqrt(n)` (Schneider endpoint
singularity), or (c) `1/n²` (Simpson on a smooth integrand). The
correct fix is product-Nyström — integrate the singular kernel
analytically against the trial basis, never sample it at zero.
This pattern WILL recur for the sphere Peierls kernel
(`r r' / |r-r'|`-style structure also has a log singularity after
the angular reduction); the new
`peierls_atkinson_nystrom` module is structured to be
sphere-extensible.

→ numerical-bug-signatures Signature: "log-singular kernel diagonal
truncation in plain quadrature"; convergence-rate fingerprint
**`err ~ log(n) / n`**, NOT `err ~ 1/sqrt(n)`.

---

## ERR-037 — Atalay Eq 42 z_0 quadrature endpoint at μ=1: bracket pole 1/(1-μ²) cancels algebraically but slowly under plain mp.quad

**Status:** **FIXED 2026-05-03**, with a permanent regression test
(11 cases × Atalay Table 1) that pins z_0(c, f₁=0) to 1e-5 absolute.

**Failure mode:** #4 Factor error — slow numerical cancellation of an
endpoint pole misdiagnosed in the Wave 2-B closeout as a "convention
drift / completeness-sum normalisation discrepancy" with Atalay's
published form. The actual cause is quadrature endpoint convergence,
not convention drift.
**Date discovered:** 2026-05-02 (Wave 2-B parallel waves close;
1.5% gap on Atalay Table 1).
**Date diagnosed and fixed:** 2026-05-03 (numerics-investigator
follow-up cascade).
**Solver:** `orpheus.derivations.continuous.case_method.core.extrapolated_endpoint.atalay_z0`
(Atalay 1997 Prog. Nucl. Energy 31(3), 229-252, Eq 42 — Milne
extrapolated endpoint for linearly-anisotropic scattering).

**Bug:** The Atalay Eq 42 integrand is
$(c/4) \nu_0 g_1(\nu) [d^2(\nu^2)(1+c\nu^2/(1-\nu^2))+3 f_1(1-c)^2 \nu^2 d(-\nu^2)] \ln((\nu_0+\nu)/(\nu_0-\nu))$
on $\nu \in (0, 1)$. The bracket carries a factor $1/(1-\nu^2)$
that becomes singular at $\nu = 1$; this pole is algebraically
cancelled by the $\lambda^2 \sim \ln^2(1-\nu)$ growth in $g_1(\nu)$
(making the integrand integrable), but the cancellation is
**inefficient under direct `mp.quad` over [0, 1]**. Even at dps=35
the integral is under-evaluated by 1.5% relative.

The fix substitutes $\mu = \tanh(t)$, mapping $(0, 1) \to (0, \infty)$
with a Jacobian $\mathrm{sech}^2(t) = 1 - \mu^2$ that **cancels the
pole exactly**. The transformed integrand is exponentially decaying
at $t \to \infty$ ($g_1 \sim 1/t^2$) and `mp.quad` resolves it to
6-7 digits at dps=25 in <1 ms.

**Impact:**

| Output                        | Pre-fix              | Post-fix              | Reduction |
| ----------------------------- | -------------------- | --------------------- | --------- |
| z_0(c=1.30, f₁=0)             | 0.535568             | 0.547144              | 1.5% → 1e-7 |
| z_0(c=1.10..2.00, f₁=0)       | 1.5-2% relative      | 6-7 digit absolute    | 1e5×       |
| Slab 2d_c, R=0, c=1.30, f₁=0  | 1.12% relative       | 0.12% relative        | 10×        |
| Slab 2d_c, R=0, c=1.50, f₁=0  | 1.23% relative       | 0.43% relative        | 3×         |
| Sphere R_c, c=1.30, f₁=0      | 0.48% relative       | 0.001% relative       | 480×       |
| Slab 2d_c, R=0, c=1.30, f₁=0.10 | 1.12%             | 0.14%                 | 8×         |

The reflected (R>0) and high-c anisotropic cases retain a residual
1-4% gap.

**2026-05-03 update — extended-fix hypothesis FALSIFIED.** The Wave 2-B
closeout speculated the same μ=tanh(t) substitution should fix the K_j
moments (`half_range.py::_atalay_K_or_L_moment_value`). A follow-up
investigation tested this directly:

- K_j outer scipy.quad CONVERGES fine (rel drift ~1e-11 between
  baseline epsabs=1e-10 and tight epsabs=1e-13). NOT Signature 7 at
  the K_j level.
- Bottleneck is X(-ν) accuracy: scipy backend has ~1e-3 relative
  error vs mpmath dps=30, propagated as X² into K_j.
- Atalay Eq 40's X-function integrand carries the same `1/(1-ν²)`
  bracket pole and cancelling factor `g_1 ~ 1/log²(1-ν)`. BUT the
  product is **logarithmically divergent** at ν → 1 (primitive
  ∝ -1/log(1-ν), bounded but truncation-dependent).
- mpmath dps 15→60: X(-0.55) drifts 0.866 → 0.860 (drift 6e-3,
  monotone, slow). Both `(0,1)` and `(0,∞)` substituted forms give
  finite values that depend on the algorithm's stopping point near
  the endpoint.
- Implementing μ=tanh(t) in `_atalay_X_function_scipy` AND `_atalay_
  X_function_mpmath` shifted the gap from 4.40% → 4.44% on R=0.75
  (marginal degradation). Reverted.

**Conclusion**: this is NOT Signature 7. The residual gap originates
in the X-function's algorithm-dependent regularisation of a divergent
integrand. Atalay's Eq 40 is either (a) missing a factor lost in
typesetting, (b) interpreted under a specific regularisation
(Cauchy PV or implicit via the Eq 26 definition path), or (c) needs
the closed-form Case-Plazcek-Hofmann 1961 X-function for isotropic.
Tracked as separate investigation, regression-pinned by
`derivations/diagnostics/diag_kj_x_function_divergent_integrand.py`.

**How it hid from higher-level tests:**

1. **L0 tests for z_0 didn't exist.** The Atalay implementation was
   verified at the **derivative level** (V_case.7 SymPy: integrand
   bracket structure matches Eq 42) but not at the **value level**
   against Atalay Table 1.
2. **The Wave 2-B closeout misdiagnosed the cause** as a
   "Case-Zweifel completeness-sum normalisation discrepancy"
   (citing 1% gap in `(c/2)·∫g_1·bracket·dν vs 1` — a separate
   issue, since the Case-Zweifel sum is `discrete + continuum = 1`
   and the continuum-only integral is supposed to be < 1).
3. **The slab/sphere critical-thickness tests were tagged as
   `assert err < 5e-2` "loose tolerances commensurate with
   implementation accuracy"** — accepting the 1.5% gap as a
   feature rather than investigating.
4. **The convergence-with-dps signature was buried:**
   z_0 at dps=15 was 0.529, at dps=25 was 0.535, at dps=35 was 0.539
   (slow monotone increase toward 0.547). This screamed
   "quadrature unconverged" but was attributed to "Atalay's
   published form has a different convention."

**Fingerprint:** Sign-pattern + magnitude scaling: error is
**uniformly negative** across c (output low by ~1.5-2%), with
magnitude weakly dependent on c (1.0-2.5% rel), and **decreasing
with dps but not converging at engineering speed**. This is the
fingerprint of a slow-convergent endpoint, not a structural error.

**L1 tests that catch it:**

- `derivations/diagnostics/diag_05_z0_regression_atalay_table1.py::test_atalay_z0_table1_isotropic`
  (parametrised over Atalay Table 1's 10 c values) — asserts
  `abs(z_0 - z_0_atalay) < 1e-5` at f₁=0. **Tagged
  `@pytest.mark.catches("ERR-037")`** for permanent regression
  binding. Promote to `tests/derivations/test_case_method_z0.py`.
- `test_atalay_z0_consistent_with_milne_form_at_c130` — cross-checks
  z_0 against Davison's Milne form $z_0 = (\pi/2) u_0 - a_c$
  using Atalay's own Table 2 row. This pin is structurally
  independent of Eq 42 — it only requires that Eq 42 reproduce
  the Milne extrapolated endpoint.

**Lesson:** When an integrand carries an endpoint pole that
"algebraically cancels" against another factor, **plain adaptive
quadrature is NOT sufficient**. The substitution $\mu = \tanh(t)$
or analogous (Atkinson product-Nyström, Gauss-Jacobi with weight
matching the singularity) is mandatory. Watch for the fingerprint:
**slow monotone convergence with dps that fails to converge at
engineering rate** is a red flag for endpoint behaviour, NOT a
"different convention." See ERR-036 for the analogous lesson on
log-singular Peierls kernel.

**Anti-pattern caught:**
- Pillar diagnosis: a quadrature accuracy gap masquerading as a
  convention drift makes the convention-bridge investigation a
  red herring. The Wave 2-B closeout asserted "1.5% z_0 gap is
  a normalisation issue, deferred — multi-day fix" — but the
  actual fix was a 1-line variable substitution.

→ numerical-bug-signatures Signature: "Endpoint pole that
algebraically cancels but slow under plain mp.quad";
fingerprint: error fixed-sign, weakly c-dependent, slow-convergent
with dps. Diagnostic probe: substitute the problematic variable
to push the pole to infinity (tanh, log, or domain-mapping) and
compare convergence speed.

---

## ERR-038 — Atalay 1997 Tables 2-5 first-order Fredholm precision floor at small slab thicknesses

**Status:** **PAPER LIMITATION CHARACTERISED 2026-05-03** — NOT a code bug.
This entry documents a *reference* limitation, not a solver defect, and
exists so future investigators don't re-chase the gap as a bug.

**Failure mode:** Reference contamination — an apparent 5% gap to a
published reference (Atalay 1997 Table 2, R=0.99 column) was investigated
as a solver bug across two earlier closeouts (Wave 2-B and the Front-3
follow-up). The actual cause is the published reference's own
documented first-order approximation, which Atalay explicitly states
"we expect some improvement in the accuracy especially for the small
slab thicknesses."
**Date discovered:** 2026-05-03 (Front-3 R=0.99 numerics-investigator
cascade).
**Solver:** `orpheus.derivations.continuous.singular_eigenfunction.slab.solve_case_method_slab_critical`
(Atalay 1997 Eq 46).

**Mechanism:** Atalay 1997 §2 (p.236) derives Eq 46 from the full
Fredholm equation Eq 32 by the explicit step "we here skip the zeroth
order and proceed directly with the first order approximation. This
provides us the required optimum accuracy. **The first order
approximation necessitates that we omit the integral term in Eq.(32)**."
Tables 2-5 are then computed from Eq 46 with first-order accuracy.
Atalay (p.246) further states that "as in the work of Kaper et al.
(1974) for isotropic scattering, one may consider to iterate further
until a better convergence is obtained… **we expect some improvement
in the accuracy especially for the small slab thicknesses**."

Empirical fingerprint (Front-3 cascade, c=1.30 f₁=0):

| 2d_atalay (mfp) | Our 2d (mfp) | Rel error |
| --------------- | ------------ | --------- |
| 20.0            | 19.99961     | 0.002%    |
| 2.0             | 2.00071      | 0.036%    |
| 0.20            | 0.20641      | 3.2%      |
| 0.01456 (R=0.99)| 0.01529      | 5.0%      |

The error scales monotonically with `1/d_crit` — the signature of an
omitted-higher-order term that vanishes at large d (where `T(R,μ)→ -1`
saturates and the omitted Fredholm integral contributes negligibly to
the boundary residue) and dominates at small d.

**Cascade evidence ruling out alternative mechanisms:**

1. **Conditioning / cancellation (mechanism 1)**: ruled out — K_j moments
   are dps-independent at bit-identical level across dps 15→40.
2. **Singular asymptotic R→1 (mechanism 2)**: ruled out — error is smooth
   in R; at R=0.99 c=1.024 with 2d=0.20 (Atalay Table 6) the same 3.2%
   error appears, demonstrating the gap is in `1/d_crit` not in `1/(1-R)`.
3. **Different quadrature pole (mechanism 3)**: ruled out — Phase 1.5
   built a μ=tanh(t) substituted X-function (mpmath) that gives 1.2%
   different X(-0.99) values, but Phase 2.2 confirmed this changes
   K_j by 0.3-0.6% and 2d_crit by <0.1% across all 15 Table 2 cells.
4. **X-function singular branch (mechanism 4)**: ruled out by the same
   evidence as 3.
5. **Atalay paper limitation (mechanism 5)**: confirmed by:
   - Atalay's explicit text on p.236 ("first order approximation
     necessitates that we omit the integral term")
   - Atalay's explicit text on p.246 ("expect some improvement…
     especially for small slab thicknesses")
   - Empirical 1/d_crit error scaling
   - Self-consistency: at moderate d our solver agrees with Atalay
     to machine precision

**X-function endpoint pole — separate, real, but non-load-bearing:**

The X-function integrand (Atalay Eq 40) carries the same `1/(1-ν²)`
bracket pole that ERR-037 fixed for z_0. Phase 1.5 of this cascade
demonstrated a μ=tanh(t) substitution gives an X-function value 1.2%
different from the legacy mpmath path. **However, this 1.2% X-function
shift propagates to 0.3-0.6% in K_j and <0.1% in 2d_crit at all R/c
in Atalay Table 2.** The X-function tanh-fix is a robustness
improvement, NOT a fix for the Table 2 R=0.99 gap. Tracked as
follow-up; not shipped in this session because production results
unchanged.

**Mode-bracketing artefact — separate sub-front, deferred:**

At R=0.99 c=1.30 the bracket-scan in `solve_case_method_slab_critical`
returns mode=1's value (0.01529) when asked for mode=2 (Atalay Table 2:
5.95846). The `sin(diff_wrapped)` residual is π-periodic and produces
spurious sign-changes that contaminate the mode index at high R. This
is a real bug but does not affect the fundamental-mode results (which
all unit tests target); deferred to a separate Issue.

**How it hid from earlier triage:**

1. The Wave 2-B closeout left "case_method R=0.99 perfect reflector
   limit" as Open Follow-up #2 with the description "Atalay's last
   column (R=0.99) is still 10%+ off". The "10%+" was an overestimate
   from a stale n_bracket; with proper bracketing the actual gap is
   2-5%, not 10%+.
2. The Wave 2 cascade strongly oriented around Signature 7 (endpoint
   pole-cancellation) which had been the cause of ERR-037. Phase 1.2
   of this cascade pointed to a similar K_j endpoint, but the
   fingerprint (dps-independent K_j) immediately ruled it out — saving
   a multi-day false trail.
3. The Wave 2-B "extended-fix hypothesis FALSIFIED" subsection of
   ERR-037 already noted that scipy.quad on K_j converges fine and
   that the gap originates in X(-ν) accuracy. This Front-3 cascade
   confirmed and quantified that finding, then identified the residual
   gap as a paper-precision-floor effect.

**Fingerprint:** Sign-pattern + magnitude scaling: error is
**positive** (our 2d_crit > Atalay's, indicating Atalay's first-order
approximation under-predicts critical thickness at small d), **scaling
as 1/d_crit**, **insensitive to dps and quadrature precision tightening**.
This is the fingerprint of a reference-side approximation, not a code-side
quadrature failure. Distinguishing test: verify the solver's value at
moderate d (where the omitted higher-order term is negligible) — if
machine-precision agreement, the small-d gap is paper-side.

**L1 tests that pin it:**

- `tests/derivations/test_case_method_slab.py::test_slab_reflected_isotropic_atalay_table2`
  — tightened from old "10%+" hand-wave to the documented 5e-2 floor
  for R∈[0.25, 0.75] and 7e-2 for R=0.99.
- New regression test: `test_slab_atalay_table2_first_order_floor_at_r099`
  pinning the 5e-2 floor at R=0.99 with the docstring documenting that
  this is the PAPER's first-order approximation floor, not the SOLVER's.
- New consistency test: `test_solver_machine_precision_at_moderate_d`
  pinning that at 2d=20 (where higher-order Fredholm contributions are
  negligible) the solver agrees with Atalay's eigenvalue table to
  1e-4 relative across multiple R values — the structural-independent
  ground for the verdict.

**Lesson:** When investigating a numerical disagreement with a
published reference, **read the paper's stated approximation level
explicitly before assuming the gap is a code bug**. Atalay's text
(p.236, p.246) twice explicitly states the published values are
first-order approximations with degraded precision at small slab
thicknesses. The diagnostic that resolves "code bug vs paper floor"
is **scaling the same physical problem to a regime where the paper's
approximation is exact** (here: large d) and verifying machine
precision there.

**Anti-pattern caught:**

- Reference contamination by under-reading the reference's own
  caveats: the Wave 2-B closeout listed R=0.99 as "still 10%+ off,
  needs careful analysis of the singular limit 2d→0", treating it as
  a **mathematical** singular-limit problem requiring multi-day
  asymptotic analysis. The actual issue is **Atalay's first-order
  approximation precision floor**, fully documented in Atalay's text,
  with no analytic limit needed and nothing to fix in our solver.
- Cross-reference grounding: ERR-037 was the CORRECT diagnosis of
  the z_0 problem (real Signature 7); ERR-038 is the CORRECT
  *non-diagnosis* — same family of solver, similar-shaped error,
  but a fundamentally different cause.

→ This is NOT a numerical-bug-signature entry. The fingerprint
"error scales with 1/d_crit, insensitive to all numerical
precision parameters, paper text states first-order approximation"
is the fingerprint of **a reference precision floor**, which lives
in `vv-principles` reference.md §reference-contamination, not in
`numerical-bug-signatures`.

---

## Meta-Lessons

1. **1-group is degenerate.** k = νΣ_f/Σ_a regardless of flux shape.
   Every solver MUST be verified at ≥2 groups.

2. **Homogeneous is degenerate.** Redistribution cancels for flat flux.
   Every curvilinear solver MUST be verified on heterogeneous problems.

3. **Per-ordinate consistency is fundamental.** L0-SN-003 would have
   caught ERR-006 (the most expensive bug) immediately.

4. **20 passing tests don't mean correct.** The cylindrical sweep
   passed 20 tests including homogeneous exact, particle balance, and
   conservation.  The bug survived because none tested what mattered:
   heterogeneous eigenvalue convergence with mesh refinement.

5. **Test with the pathological case.** Every physics feature has a
   regime where bugs are exposed.  For curvilinear: test near r=0 with
   mesh refinement.  For scattering: test with asymmetric matrices.
   For normalization: test with multiple quadrature types.

6. **Zero cross sections hide bugs.** If a reaction type (n,2n,
   upscatter) is zero in all test materials, every code path touching
   it is untested.  ERR-015 survived because `Sig2 = 0` everywhere.
   For every XS term, there must be at least one test where it's nonzero.

7. **A tautological residual proves nothing.** If the convergence
   check computes `f(x) - g(f(x))` where `g` is the inverse of the
   step that produced `f(x)`, the residual is identically zero.
   ERR-016 survived because "all inner iterations = 1" was mistaken
   for fast convergence.  Always verify with a problem that SHOULD
   require multiple iterations.

8. **Geometry area/volume must match across methods.** When
   comparing solvers with different geometry representations (e.g.,
   MC square cell vs CP Wigner-Seitz cylinder), verify that the
   cell area/volume is equal.  A factor-of-2 in a linear dimension
   is a factor-of-4 in area — enough to change supercritical to
   subcritical.  ERR-017 survived because all homogeneous tests
   passed and the heterogeneous tests were `@pytest.mark.slow`.

9. **Angular integration weights have hidden factors.** In MOC,
   the weight `omega_a * omega_p * t_s` is the spatial/angular
   discretization weight, but the scalar flux integral also requires
   `4*pi` (full-sphere normalization) and `sin(theta_p)` (2D→3D
   projection).  These factors cancel for spatially uniform solutions,
   making them invisible to homogeneous tests.  ERR-019 survived
   three homogeneous tests at machine precision because delta_psi = 0
   for uniform media.  Derive the weight formula from first principles
   and verify against a heterogeneous problem BEFORE trusting the
   homogeneous result.


---

## ERR-039 — HarmonicMomentProjection.apply_transpose claimed Π* = R but used the addition-theorem reconstruction (with (2ℓ+1) factor)

**Status:** **CAUGHT IN QA REVIEW 2026-05-10**, fixed in same commit as
introduction (Wave 0 / C0.5 of SN performance plan).

**Failure mode:** **#6 (convention drift)** — definition site (docstring
claiming "Π* = R, the addition-theorem reconstruction") and the
mathematical adjoint identity ⟨Πψ, c⟩_C = ⟨ψ, Π* c⟩_V disagreed by
the factor (2ℓ+1) summed over ℓ. The author drafted apply_transpose
as a delegation to HarmonicMomentReconstruction.apply, which is the
addition-theorem reconstruction used by the SN scattering source —
not the W-weighted Hilbert adjoint.

**Date discovered:** 2026-05-10 (Wave 0 V&V audit by qa agent).
**Module:** `orpheus.numerics.projection.HarmonicMomentProjection`.

**Mechanism:** Under the V-inner-product
⟨ψ, φ⟩_V = Σ_n w_n ψ_n φ_n (which absorbs the quadrature weights as
part of the inner product, NOT as part of the operator), the adjoint
of Π = Y* W satisfies

    ⟨Π ψ, c⟩_C  =  Σ_n ψ_n w_n (Σ_lm Y_lm(n) c_lm)  =  ⟨ψ, Π* c⟩_V

so (Π* c)_n = Σ_lm Y_lm(n) c_lm — the **naked** reconstruction with
NO (2ℓ+1) factor. The original implementation built a transient
HarmonicMomentReconstruction (with (2ℓ+1) factor for the addition
theorem) and called its apply, returning a vector that is (2ℓ+1)·Π*c
on each ℓ-block. Numerical fingerprint:

| L | quadrature | ⟨Π ψ, c⟩_C | ⟨ψ, R c⟩_V (with factor) | ⟨ψ, Π*c⟩_V (no factor, correct) |
|---|---|---|---|---|
| 3 | Lebedev 13 | 4.7488 | 27.6468 | 4.7488 |

The factor-of-(sum over ℓ of (2ℓ+1)) discrepancy is exactly the
fingerprint of the (2ℓ+1) addition-theorem factor wrongly applied to
the adjoint side.

**How it hid:** The capability set advertised CAP_APPLY_TRANSPOSE, but
the only test of apply_transpose was `test_projection_advertises_apply_and_apply_transpose`
which checked the capability set is non-empty — NOT that the method
works. The Galerkin adjoint-pairing test `TestGalerkinAdjointPairing`
hand-coded the einsum reference but never crossed the boundary into
the production apply_transpose path. The bug was hidden by missing
test coverage of the production code path, not by an algorithmic
near-cancellation.

**How caught:** `qa` agent's V&V review of Wave 0 detected the
capability dishonesty by reading the implementation against the
docstring claim, then numerically checked ⟨Πψ, c⟩_C vs
⟨ψ, R c⟩_V (with factor) on a Lebedev-13 / L=3 instance — found
the (2ℓ+1)-summed discrepancy.

**Fix:** apply_transpose now computes
    Π* c = Σ_lm Y_lm(n) c_lm   (single np.einsum, no factor)
matching the true W-weighted adjoint identity. Direct test
`TestApplyTransposeIsWWeightedAdjoint` added to gate the production
code path against the adjoint identity (lhs == rhs to 1e-12) AND
against the distinction Π* ≠ R (the (2ℓ+1) factor must NOT appear on
the adjoint side).

**Lesson:** **Every capability that an operator advertises MUST be
tested by a check that crosses the production boundary, not by a
hand-coded reference of the same math.** Capability-set advertisement
is a contract; the contract is hollow if the only test is "caps is
non-empty." The Wave 0 V&V audit caught this by reading for
discipline-vs-implementation drift, but only because qa specifically
looked at apply_transpose. Going forward, every new LinearOperator
that advertises CAP_APPLY_TRANSPOSE or CAP_SOLVE MUST ship with a
direct test of that method against an independent reference (the
adjoint identity, or the solve round-trip on a non-trivial input).

**Test reference:** `tests/numerics/test_projection_operators.py::TestApplyTransposeIsWWeightedAdjoint::test_adjoint_identity_matches_production` and `::test_apply_transpose_no_2l_plus_1_factor`. Tagged `@pytest.mark.l1` and `@pytest.mark.catches("ERR-039")`.

